"""Complete anchored-grid search modulo drone relabelling and time shift."""
from itertools import combinations, product
from math import comb
from time import perf_counter
from .dynamics import Engine


def candidates(graph, m):
    q, n = graph.Q, graph.n
    if not 2 <= m <= n*q:
        raise ValueError("Require 2 <= fleet size <= n*Q")
    all_positions = range(n*q)
    for i, j, a in graph.edges:
        p, r = i*q+a, j*q+(a+q//2) % q
        others = [x for x in all_positions if x not in (p, r)]
        for extra in combinations(others, m-2):
            positions = (p, r)+extra
            for signs in product((0, 1), repeat=m):
                yield tuple(sorted(2*p+s for p, s in zip(positions, signs)))


def enumerate_schedules(graph, m, max_steps=100000, on_found=None):
    start, engine = perf_counter(), Engine(graph)
    # All visited occupancy states are equivalent under labels and grid-time
    # translation; invalid forward paths also certify invalidity from each state.
    visited, records = set(), []
    counts = dict(candidate_bound=len(graph.edges)*comb(graph.n*graph.Q-2, m-2)*2**m,
                  generated=0, simulated=0, collision_runs=0,
                  disconnected_orbits=0, unknown_runs=0, valid_connected_orbits=0)
    for state in candidates(graph, m):
        counts["generated"] += 1
        if state in visited:
            continue
        result = engine.simulate(state, max_steps)
        counts["simulated"] += 1
        if result.status == "unknown":
            # Unknown paths may end in different outcomes. Do not prune them.
            counts["unknown_runs"] += 1
            continue
        visited.update(tuple(sorted(s)) for s in result.states)
        if result.status == "collision":
            counts["collision_runs"] += 1
        elif not result.connected:
            counts["disconnected_orbits"] += 1
        else:
            representative = result.canonical
            result = engine.simulate(representative, max_steps)
            assert result.connected
            record = {"state": list(representative),
                      "assignments": graph.describe_state(representative),
                      "period_ticks": result.period_ticks,
                      "shifts_per_period": result.shifts,
                      "sync_pairs": sorted({p for batch in result.syncs for p in batch})}
            records.append(record)
            counts["valid_connected_orbits"] += 1
            if on_found:
                on_found(record)
    records.sort(key=lambda r: r["state"])
    counts["complete"] = counts["unknown_runs"] == 0
    counts["seconds"] = perf_counter()-start
    counts["visited_states"] = len(visited)
    return {"graph": graph.name, "n": graph.n, "edges": len(graph.edges),
            "Q": graph.Q, "fleet_size": m, "counts": counts, "schedules": records}
