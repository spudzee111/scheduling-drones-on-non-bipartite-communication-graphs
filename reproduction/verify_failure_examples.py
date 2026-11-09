#!/usr/bin/env python3
"""Verify two single-failure examples using exact rational event times.

Standalone: Python's standard library only; no imports from the thesis's
grid engine. Run this file to regenerate failure_examples.json beside it.
Drone and trajectory labels are one-based. All removals happen immediately before the link arrivals at time 1.
Initial schedules are departure states, with no shifts at time zero.
"""

from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json


# Each edge stores the bearing on its lower-labelled trajectory, in turns.
EDGES = ((1, 2, F(0)), (1, 3, F(1, 6)), (2, 3, F(1, 3)))
PARTNER = {}
PORTS = {circle: [] for circle in range(1, 4)}
for left, right, angle in EDGES:
    opposite = (angle + F(1, 2)) % 1
    PARTNER[left, angle] = (right, opposite)
    PARTNER[right, opposite] = (left, angle)
    PORTS[left].append(angle)
    PORTS[right].append(opposite)

# State entries: (persistent drone label, trajectory, angle, orientation).
INITIAL = tuple(
    (2 * (circle - 1) + phase + 1, circle, F(phase, 2), 1)
    for circle in range(1, 4)
    for phase in (0, 1)
)


def state_data(state):
    return [dict(drone=label, trajectory=circle, angle=str(angle),
                 orientation=direction)
            for label, circle, angle, direction in state]


def insert_interval(intervals, circle, start, direction, duration):
    """Store the actual angular interval traversed, splitting at angle zero."""
    if duration >= 1:
        intervals[circle].append((F(0), F(1)))
        return
    left = (start if direction == 1 else start - duration) % 1
    right = left + duration
    intervals[circle].append((left, min(right, F(1))))
    if right > 1:
        intervals[circle].append((F(0), right - 1))


def covered_lengths(intervals):
    answer = []
    for circle in range(1, 4):
        right, length = F(0), F(0)
        for left, end in sorted(intervals[circle]):
            length += max(F(0), end - max(left, right))
            right = max(right, end)
        answer.append(length)
    return answer


def components(labels, edges):
    adjacency = {label: set() for label in labels}
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    unseen, groups = set(labels), []
    while unseen:
        first = min(unseen)
        unseen.remove(first)
        group, stack = {first}, [first]
        while stack:
            for neighbor in adjacency[stack.pop()] & unseen:
                unseen.remove(neighbor)
                group.add(neighbor)
                stack.append(neighbor)
        groups.append(sorted(group))
    return groups


def apply_arrivals(state):
    """Apply the rule to an arrival configuration, preserving drone labels."""
    occupied = {(c, a): label for label, c, a, direction in state}
    if len(occupied) != len(state):
        raise ValueError("Collision before the arrival update")
    post, shifts, pairs = [], [], []
    for label, circle, angle, direction in state:
        other = PARTNER.get((circle, angle))
        partner = occupied.get(other)
        if other is not None and partner is None:
            post.append((label, *other, -direction))
            shifts.append(label)
        else:
            post.append((label, circle, angle, direction))
            if partner is not None and label < partner:
                pairs.append((label, partner))
    return tuple(post), shifts, pairs


def simulate(initial, max_events=10000):
    """Independently simulate a departure state using exact event times."""
    state = first = tuple(initial)
    time = F(0)
    seen, edges = set(), set()
    intervals = {circle: [] for circle in range(1, 4)}
    events, shift_count = [], 0
    for _ in range(max_events):
        if state in seen:
            assert state == first, "Unexpected transient"
            lengths = covered_lengths(intervals)
            labels = [label for label, c, a, s in state]
            groups = components(labels, edges)
            met = {label for pair in edges for label in pair}
            return dict(status="periodic", labelled_period=str(time),
                        collision_free=True, full_coverage=all(x == 1 for x in lengths),
                        covered_lengths=list(map(str, lengths)),
                        synchronization_edges=[list(e) for e in sorted(edges)],
                        synchronization_components=groups,
                        synchronization_connected=len(groups) == 1,
                        starving_drones=[label for label in labels if label not in met],
                        valid=len(groups) == 1 and all(x == 1 for x in lengths),
                        encounter_events=events, shift_count=shift_count)
        seen.add(state)
        occupied = {(c, a): label for label, c, a, s in state}
        assert len(occupied) == len(state), "Co-location"
        # Initial contacts count as communications, without causing shifts.
        for label, c, a, s in state:
            other = occupied.get(PARTNER.get((c, a)))
            if other is not None and label < other:
                edges.add((label, other))
                events.append(dict(time=str(time), drones=[label, other]))
        delay = min(((bearing-a)*s) % 1 or F(1)
                    for label, c, a, s in state for bearing in PORTS[c])
        collisions = []
        for left, right in combinations(state, 2):
            a, c, angle, s = left
            b, other_c, other_angle, other_s = right
            if c == other_c and s != other_s:
                dt = ((other_angle-angle)*s) % 1 / 2
                if dt <= delay:
                    collisions.append((dt, a, b, c, (angle+s*dt) % 1))
        travel = min(x[0] for x in collisions) if collisions else delay
        for label, c, a, s in state:
            insert_interval(intervals, c, a, s, travel)
        if collisions:
            return dict(status="collision", first_collision_time=str(time+travel),
                        first_collisions=[dict(drones=[a,b], trajectory=c, angle=str(angle))
                                          for dt,a,b,c,angle in collisions if dt == travel],
                        collision_free=False, valid=False, labelled_period=None,
                        encounter_events_before_collision=events,
                        shift_count_before_collision=shift_count)
        arrivals = tuple((label,c,(a+s*delay)%1,s) for label,c,a,s in state)
        state, shifts, pairs = apply_arrivals(arrivals)
        shift_count += len(shifts)
        time += delay
    raise RuntimeError("Event limit reached; outcome remains unknown")


def main():
    baseline = simulate(INITIAL)
    assert baseline["valid"] and baseline["labelled_period"] == "1"
    assert baseline["shift_count"] == 0
    assert baseline["covered_lengths"] == ["1"] * 3
    assert len(baseline["synchronization_edges"]) == 6
    # The no-shift baseline reaches INITIAL again at t=1, before processing
    # that arrival. Remove the specified drone then, and process survivors.
    failure_time = F(1)
    cases = []
    for removed in (1, 4):
        survivors = tuple(drone for drone in INITIAL if drone[0] != removed)
        post, shifted, pairs = apply_arrivals(survivors)
        result = simulate(post)
        if result["status"] == "collision":
            result["absolute_collision_time"] = str(failure_time + F(result["first_collision_time"]))
        cases.append(dict(removed_drone=removed, failure_time=str(failure_time),
                          arrivals_after_removal=state_data(survivors),
                          shifts_at_failure=shifted,
                          synchronizations_at_failure=[list(e) for e in pairs],
                          initial_survivors=state_data(post), outcome=result))
    collision, isolation = (c["outcome"] for c in cases)
    assert collision["first_collision_time"] == "1/4"
    assert collision["absolute_collision_time"] == "5/4"
    assert collision["first_collisions"] == [dict(drones=[2,4],trajectory=1,angle="3/4")]
    assert collision["encounter_events_before_collision"] == []
    assert isolation["collision_free"] and isolation["labelled_period"] == "3"
    assert isolation["full_coverage"] and isolation["synchronization_edges"] == []
    assert isolation["starving_drones"] == [1,2,3,5,6]
    assert isolation["shift_count"] == 30 and not isolation["valid"]

    # Cross-check the independent Fraction calculation against the grid engine.
    import sys
    root = Path(__file__).resolve().parents[1] / "Drones"
    sys.path.insert(0, str(root)); sys.path.insert(0, str(root / "scripts"))
    from scs.model import Graph
    from scs.dynamics import Engine
    from run_vse_experiments import grid_properties
    graph = Graph.read(root / "benchmarks/vse/triangle.json")
    engine = Engine(graph)
    for case in cases:
        rows = case["initial_survivors"]
        state = graph.read_state([(v["trajectory"]-1,v["angle"],v["orientation"]) for v in rows])
        grid = engine.simulate(state)
        expected = case["outcome"]
        assert grid.status == expected["status"]
        if grid.status == "collision":
            assert F(grid.collision_tick, graph.Q) == F(expected["first_collision_time"])
            assert sorted(rows[i]["drone"] for i in grid.collision_pair) == [2,4]
        else:
            assert F(grid.period_ticks,graph.Q) == F(expected["labelled_period"])
            assert not grid_properties(graph,grid)[0]
            assert not any(grid.syncs)
    output = dict(
        method="Independent Fraction event simulation, cross-checked against grid simulation",
        state_convention="Departure states; no shifts at time zero",
        failure_convention="Remove immediately before arrivals at absolute time 1; process those arrivals once",
        time_convention="Outcome times are elapsed since failure, with absolute collision time also recorded",
        label_convention="Trajectory labels 1..3; persistent drone labels 1..6",
        geometry=dict(unit_circle_radius="1",communication_radius="2/5",
                      inclusive_link_threshold="12/5",all_center_distances="12/5"),
        initial_state=state_data(INITIAL), baseline=baseline, failures=cases)
    destination = Path(__file__).with_name("failure_examples.json")
    destination.write_text(json.dumps(output,indent=2)+"\n")
    print("Baseline: valid, period 1, no shifts.")
    print("Remove drone 1 before t=1 arrivals: collision at t=5/4 (1/4 after failure).")
    print("Remove drone 4 before t=1 arrivals: residual period 3, full coverage, no communication.")
    print("Both failure outcomes independently verified against the grid engine.")


if __name__ == "__main__":
    main()
