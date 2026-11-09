#!/usr/bin/env python3
"""Complete small-instance searches with full coverage and orientation checks.

The old enumeration core returns collision-free connected orbits. This driver
applies the thesis's stronger validity definition and independently replays all
such orbits with exact rational event times and interval-union coverage.
"""
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from math import comb
import json
import platform
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scs.model import Graph
from scs.dynamics import Engine
from scs.enumeration import enumerate_schedules

SUITE = [("square", range(2, 6)), ("triangle", range(2, 7)),
         ("diamond", range(2, 6))]


def check_geometry(path):
    """Check exact centers, unit-circle separation, bearings, and every link."""
    data = json.loads(path.read_text())
    centers = [[tuple(map(F, component)) for component in point]
               for point in data["exact_centers"]]
    rho, step = F(data["metadata"]["communication_radius"]), F("12/5")
    edges = {(i,j): F(a) for i,j,a in data["edges"]}
    vectors = {F(0): ((F(1),F(0)),(F(0),F(0))),
               F(1,6): ((F(1,2),F(0)),(F(0),F(1,2))),
               F(1,3): ((F(-1,2),F(0)),(F(0),F(1,2))),
               F(1,4): ((F(0),F(0)),(F(1),F(0)))}
    for i in range(len(centers)):
        for j in range(i+1,len(centers)):
            vector = [tuple(y-x for x,y in zip(a,b))
                      for a,b in zip(centers[i],centers[j])]
            squared = (sum(a*a+3*b*b for a,b in vector),
                       sum(2*a*b for a,b in vector))
            assert squared[1] == 0  # All squared distances are rational here.
            assert squared[0] > 4
            linked = squared[0] <= (2+rho)**2
            assert linked == ((i,j) in edges)
            if linked:
                assert squared[0] == step**2 == (2+rho)**2
                expected = [tuple(step*c for c in component)
                            for component in vectors[edges[i,j]]]
                assert vector == expected


def grid_properties(graph, result):
    """Measure angular grid cells traversed from departure states."""
    engine = Engine(graph)
    covered, same = set(), []
    for tick, (state, pairs) in enumerate(zip(result.states, result.syncs)):
        for i, j in pairs:
            if engine.decoded[state[i]][2] == engine.decoded[state[j]][2]:
                same.append((tick, i, j))
        for code in state:
            circle, angle, direction = engine.decoded[code]
            cell = angle if direction == 1 else (angle - 1) % graph.Q
            covered.add((circle, cell))
    missing = sorted(set((c, a) for c in range(graph.n)
                         for a in range(graph.Q)) - covered)
    return missing, same


def add_interval(intervals, circle, start, direction, duration):
    """Insert the actual traversed circular interval, independently of the grid."""
    if duration >= 1:
        intervals[circle].append((F(0), F(1)))
        return
    left = (start if direction == 1 else start - duration) % 1
    right = left + duration
    intervals[circle].append((left, min(right, F(1))))
    if right > 1:
        intervals[circle].append((F(0), right - 1))


def reference_properties(graph, initial, max_events=100000):
    """Independent Fraction replay from an event-aligned departure state.

    Initial paired endpoints communicate; unmatched initial endpoints do not
    shift. Every subsequent update follows a strictly positive travel time.
    """
    state = tuple((c, F(a, graph.Q), s)
                  for c, a, s in map(graph.decode, initial))
    first = state
    ports = [[F(a, graph.Q) for a in p] for p in graph.ports]
    partner = {(c, F(a, graph.Q)):
               (other // graph.Q, F(other % graph.Q, graph.Q))
               for c, p in enumerate(graph.ports) for a in p
               for other in [graph.partner[c * graph.Q + a]]}
    intervals = [[] for _ in range(graph.n)]
    syncs, same, seen, time = set(), False, set(), F(0)
    for _ in range(max_events):
        if state in seen:
            assert state == first, "Unexpected transient"
            break
        seen.add(state)
        occupied = {}
        for i, (c, a, s) in enumerate(state):
            assert (c, a) not in occupied, "Co-located drones"
            occupied[c, a] = i
        for i, (c, a, s) in enumerate(state):
            other = partner.get((c, a))
            if other is None:
                continue
            j = occupied.get(other)
            if j is not None:
                syncs.add(tuple(sorted((i, j))))
                same |= s == state[j][2]
        dt = min(((b-a)*s) % 1 or F(1)
                 for c, a, s in state for b in ports[c])
        for i, (c, a, s) in enumerate(state):
            for cc, aa, ss in state[:i]:
                if c == cc and s != ss:
                    distance = ((aa-a)*s) % 1
                    assert distance != 0 and distance/2 > dt, "Collision"
            add_interval(intervals, c, a, s, dt)
        arrivals = tuple((c, (a+s*dt) % 1, s) for c, a, s in state)
        occupied = {(c, a) for c, a, s in arrivals}
        state = tuple((*partner[c, a], -s)
                      if (c, a) in partner and partner[c, a] not in occupied
                      else (c, a, s) for c, a, s in arrivals)
        time += dt
    else:
        raise AssertionError("Independent replay did not finish")
    covered_lengths = []
    for pieces in intervals:
        right, total = F(0), F(0)
        for left, end in sorted(pieces):
            total += max(F(0), end-max(left, right))
            right = max(right, end)
        covered_lengths.append(total)
    return dict(period=time, syncs=syncs, opposite_only=not same,
                covered_lengths=covered_lengths)


def main():
    output_dir = ROOT / "results" / "vse"
    output_dir.mkdir(parents=True, exist_ok=True)
    cases = []
    replayed = rejected = 0
    for name, sizes in SUITE:
        input_path = ROOT / "benchmarks" / "vse" / (name + ".json")
        check_geometry(input_path)
        graph = Graph.read(input_path)
        for size in sizes:
            result = enumerate_schedules(graph, size)
            assert result["counts"]["complete"], "Cutoff prevents an exact count"
            count_v = count_o = count_no_shift = 0
            for record in result["schedules"]:
                state = graph.read_state(record["assignments"])
                orbit = Engine(graph).simulate(state)
                assert orbit.connected
                missing, same = grid_properties(graph, orbit)
                reference = reference_properties(graph, state)
                assert reference["period"] == F(orbit.period_ticks, graph.Q)
                assert reference["syncs"] == {tuple(p) for p in record["sync_pairs"]}
                assert reference["opposite_only"] == (not same)
                missing_by_circle = [sum(c == i for c, a in missing)
                                     for i in range(graph.n)]
                assert reference["covered_lengths"] == [F(graph.Q-m, graph.Q)
                                                          for m in missing_by_circle]
                record["full_coverage"] = not missing
                record["missing_cells"] = missing
                record["opposite_only"] = not same
                record["first_same_orientation_encounter"] = same[0] if same else None
                record["covered_lengths"] = [str(x) for x in reference["covered_lengths"]]
                record["vse_valid"] = not missing
                record["ovse_valid"] = not missing and not same
                count_v += not missing
                count_o += not missing and not same
                count_no_shift += not missing and orbit.shifts == 0
                replayed += 1
                rejected += bool(missing)
            # Preserve exact output with fields explaining the stricter filter.
            result["vse_counts"] = dict(valid=count_v, opposite_only=count_o,
                                         without_shifts=count_no_shift,
                                         incomplete_coverage=sum(not r["full_coverage"]
                                                                 for r in result["schedules"]))
            (output_dir / f"{name}_d{size}.json").write_text(json.dumps(result, indent=2)+"\n")
            cases.append(dict(graph=name, n=graph.n, rings=len(graph.rings), Q=graph.Q,
                              drones=size, candidates=result["counts"]["candidate_bound"],
                              valid=count_v, opposite_only=count_o,
                              without_shifts=count_no_shift,
                              incomplete_coverage=result["vse_counts"]["incomplete_coverage"],
                              complete=True))
            print(name, size, count_v, count_o, count_no_shift, flush=True)
    sources = [ROOT/"scs"/f for f in ("model.py", "dynamics.py", "enumeration.py")]
    sources += [Path(__file__)]
    sources += [ROOT/"benchmarks"/"vse"/(name+".json") for name, _ in SUITE]
    manifest = dict(state_convention="departure states; no shifts at time zero",
                    fleet_size_semantics="exact counts; decision budget means at most d",
                    python=platform.python_version(), platform=platform.platform(),
                    hashes={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()
                            for p in sources},
                    independent_replays=replayed, incomplete_coverage_rejections=rejected,
                    all_complete=True, exact_geometry_verified=True, cases=cases)
    (output_dir/"summary.json").write_text(json.dumps(manifest, indent=2)+"\n")
    labels = {"square":r"Square", "triangle":r"Triangle", "diamond":r"Diamond"}
    rows = [f"{labels[c['graph']]} & {c['drones']} & {c['candidates']:,} & "
            f"{c['valid']} & {c['opposite_only']} & {c['without_shifts']} \\\\"
            for c in cases]
    (ROOT.parent/"figures"/"vse_enumeration.tex").write_text("\n".join(rows)+"\n")
    print(f"Independent exact replays: {replayed}; incomplete coverage: {rejected}")


if __name__ == "__main__":
    main()
