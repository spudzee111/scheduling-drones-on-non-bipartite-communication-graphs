#!/usr/bin/env python3
"""Exact minimum removals destroying residual validity, with a survivor.

At each stage the current schedule is valid. All its grid phases are tried,
with a nonempty proper subset removed immediately before the arrival update.
A removal between ticks can be delayed to the next tick: no surviving drone
changes its motion before that tick. Batches at different times are handled
recursively, with strictly decreasing fleet size. If an intermediate schedule
is invalid, stop: later deletions cannot lower the cost of that witness.

No claim is made about preserving communication or coverage from before a
removal. The surviving execution must satisfy validity on its own.
"""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scs.model import Graph
from scs.dynamics import Engine
from run_vse_experiments import grid_properties, reference_properties


class RemovalSearch:
    def __init__(self, graph, independent=True):
        self.g = graph
        self.engine = Engine(graph)
        self.cache = {}
        self.orbits = {}
        self.periods = {}
        self.answers = {}
        self.independent = independent
        self.replays = 0

    def inspect(self, initial):
        initial = tuple(sorted(initial))
        if initial in self.cache:
            return self.cache[initial]
        result = self.engine.simulate(initial)
        if result.status == 'unknown':
            raise RuntimeError('Simulation cutoff: cannot certify this search')
        collision = result.status == 'collision'
        if collision:
            value = dict(valid=False, collision=True, connected=None,
                         full_coverage=None)
        else:
            missing, _ = grid_properties(self.g, result)
            value = dict(valid=result.connected and not missing, collision=False,
                         connected=result.connected, full_coverage=not missing)
        if self.independent:
            # The event-based reference returns to an event-aligned start.
            # A removal can leave no drone at a link in this grid snapshot;
            # in that case replay the same periodic execution from its next
            # event, keeping the labelled order supplied by the grid engine.
            reference_start = initial
            if not collision:
                reference_start = next(state for state in result.states
                                       if any(code//2 in self.g.partner for code in state))
            try:
                ref = reference_properties(self.g, reference_start)
            except AssertionError as exc:
                if not collision or str(exc) not in ('Collision', 'Co-located drones'):
                    raise
            else:
                assert not collision
                assert ref['period'] == F(result.period_ticks, self.g.Q)
                pairs = {tuple(p) for batch in result.syncs for p in batch}
                assert ref['syncs'] == pairs
                assert ref['covered_lengths'] == [
                    F(self.g.Q-sum(c == i for c, a in missing), self.g.Q)
                    for i in range(self.g.n)]
            self.replays += 1
        if value['valid']:
            phases = {}
            for tick, state in enumerate(result.states):
                phases.setdefault(tuple(sorted(state)), tick)
            key = min(phases)
            value['key'] = key
            origin = phases[key]
            self.orbits[key] = {state: (tick-origin) % result.period_ticks
                                for state, tick in phases.items()}
            self.periods[key] = result.period_ticks
        for state in result.states:
            key = tuple(sorted(state))
            assert key not in self.cache or self.cache[key] == value
            self.cache[key] = value
        return value

    def analyze(self, initial):
        value = self.inspect(initial)
        assert value['valid']
        key = value['key']
        if key in self.answers:
            return self.answers[key]
        size = len(key)
        children = {}
        batch_min = None
        batch_witness = None
        one_counts = Counter()
        scenarios = 0
        # A labelled period may visit an unlabelled state more than once.
        # Subsets at such repeats give the same set of future executions.
        for removed_count in range(1, size):
            found = False
            for departure, tick in self.orbits[key].items():
                # The link update is its own inverse; recover the arrivals.
                arrivals, _, _ = self.engine.link_update(departure)
                for removed in combinations(range(size), removed_count):
                    removed_set = set(removed)
                    survivors = tuple(code for j, code in enumerate(arrivals)
                                      if j not in removed_set)
                    after, _, _ = self.engine.link_update(survivors)
                    outcome = self.inspect(after)
                    step = dict(from_state=list(key), phase_tick=tick,
                                wait_ticks=tick or self.periods[key], Q=self.g.Q,
                                original_departure=list(departure),
                                arrivals=list(arrivals), removed_indices=list(removed),
                                surviving_departure=list(after))
                    scenarios += 1
                    if removed_count == 1:
                        one_counts['scenarios'] += 1
                        one_counts['valid'] += outcome['valid']
                        if outcome['collision']:
                            one_counts['collision'] += 1
                        elif not outcome['valid']:
                            one_counts['collision_free_invalid'] += 1
                            one_counts['disconnected'] += not outcome['connected']
                            one_counts['incomplete_coverage'] += not outcome['full_coverage']
                    if not outcome['valid']:
                        found = True
                        if batch_witness is None:
                            batch_witness = step | {'outcome': outcome}
                    else:
                        child = outcome['key']
                        if child not in children or removed_count < children[child][0]:
                            children[child] = (removed_count, step)
            if found:
                batch_min = removed_count
                break
        best = batch_min
        witness = [batch_witness] if batch_witness else None
        for child, (cost, step) in children.items():
            # At least one more deletion is necessary from a valid child.
            if best is not None and cost + 1 >= best:
                continue
            answer = self.analyze(child)
            extra = answer['minimum_sequence_removals']
            if extra is not None and (best is None or cost + extra < best):
                best = cost + extra
                witness = [step] + answer['sequence_witness']
        # For a threshold of two, also seek an explicitly separated pair of
        # one-drone removals. Wait at least a full surviving period before
        # the second recorded phase if necessary. This witnesses the result
        # even when a deletion sequence is required to have distinct times.
        if best == 2:
            for child, (cost, step) in children.items():
                if cost != 1:
                    continue
                answer = self.analyze(child)
                if answer['minimum_sequence_removals'] == 1:
                    witness = [step] + answer['sequence_witness']
                    break
        answer = dict(initial=list(key), fleet_size=size,
                      minimum_simultaneous_removals=batch_min,
                      minimum_sequence_removals=best,
                      single_removal_counts=dict(one_counts),
                      tested_batch_scenarios=scenarios,
                      simultaneous_witness=batch_witness,
                      sequence_witness=witness)
        self.answers[key] = answer
        return answer

    def verify_sequence(self, root, answer):
        """Check that recorded stages are reachable, delete actual occupants,
        and end in an invalid execution after the reported number of removals.
        State-based removal choices do not require persistent label names.
        """
        current = tuple(root)
        count = 0
        for index, step in enumerate(answer['sequence_witness'] or []):
            orbit = self.engine.simulate(current)
            assert self.inspect(current)['valid']
            assert tuple(step['from_state']) in {tuple(sorted(s)) for s in orbit.states}
            start = tuple(step['from_state'])
            shifted = self.engine.simulate(start)
            departure = shifted.states[step['wait_ticks'] % shifted.period_ticks]
            assert tuple(sorted(departure)) == tuple(step['original_departure'])
            assert step['wait_ticks'] > 0
            # Reconstruct both batches directly from positions and partners.
            def reference_batch(state):
                occupied = {x//2 for x in state}
                result = []
                for code in state:
                    other = self.g.partner.get(code//2)
                    result.append(2*other + 1-code%2
                                  if other is not None and other not in occupied else code)
                return tuple(result)
            arrivals = reference_batch(step['original_departure'])
            assert arrivals == tuple(step['arrivals'])
            survivors = tuple(code for j, code in enumerate(arrivals)
                              if j not in step['removed_indices'])
            assert survivors
            current = reference_batch(survivors)
            assert current == tuple(step['surviving_departure'])
            count += len(step['removed_indices'])
            assert self.inspect(current)['valid'] == (index+1 < len(answer['sequence_witness']))
        if answer['sequence_witness']:
            assert count == answer['minimum_sequence_removals']


def main():
    result_dir = ROOT / 'results/vse'
    cases = []
    replays = 0
    for name in ('square', 'triangle', 'diamond'):
        g = Graph.read(ROOT / f'benchmarks/vse/{name}.json')
        search = RemovalSearch(g)
        for path in sorted(result_dir.glob(name + '_d*.json')):
            source = json.loads(path.read_text())
            records = []
            for index, record in enumerate(source['schedules']):
                if not record['vse_valid']:
                    continue
                answer = search.analyze(record['state'])
                search.verify_sequence(record['state'], answer)
                if answer['minimum_sequence_removals'] == 2:
                    assert len(answer['sequence_witness']) == 2
                    assert all(len(step['removed_indices']) == 1
                               for step in answer['sequence_witness'])
                if record['ovse_valid']:
                    assert not answer['single_removal_counts'].get('collision', 0)
                occupancy = Counter(g.arc_ring[g.arc_at(code)] for code in record['state'])
                ring_occupancy = [occupancy.get(i, 0) for i in range(len(g.rings))]
                if record['ovse_valid']:
                    assert min(ring_occupancy) == answer['minimum_sequence_removals']
                records.append(dict(source_schedule_index=index,
                                    state=record['state'], opposite_only=record['ovse_valid'],
                                    initial_ring_occupancy=ring_occupancy,
                                    **answer))
            dist = Counter(str(r['minimum_sequence_removals']) for r in records)
            case = dict(graph=name, drones=source['fleet_size'], Q=g.Q,
                        valid_classes=len(records), sequence_distribution=dict(dist),
                        simultaneous_distribution=dict(Counter(
                            str(r['minimum_simultaneous_removals']) for r in records)),
                        schedules=records)
            cases.append(case)
            print(name, source['fleet_size'], dict(dist), flush=True)
        replays += search.replays
        print('Independent residual replays:', name, search.replays, flush=True)
    sources = [Path(__file__), ROOT / 'scs/dynamics.py',
               ROOT / 'scripts/run_vse_experiments.py']
    sources += sorted((ROOT / 'benchmarks/vse').glob('*.json'))
    sources += [p for name in ('square', 'triangle', 'diamond')
                for p in sorted(result_dir.glob(name + '_d*.json'))]
    output = dict(complete=True,
                  definition='minimum deletions in a sequence leaving at least one survivor and an invalid residual schedule',
                  timing='all grid phases before arrival actions; arbitrary between-event removals have equivalent futures',
                  none_means='no destructive deletion sequence with a surviving drone',
                  independent_residual_replays=replays,
                  all_sequence_witnesses_replayed=True,
                  all_thresholds_certified_by_single_removals_and_two_removal_witnesses=
                      all(r['minimum_sequence_removals'] in (1, 2)
                          for c in cases for r in c['schedules']),
                  hashes={str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest()
                          for p in sources}, cases=cases)
    (result_dir / 'removal_robustness.json').write_text(json.dumps(output, indent=2)+'\n')
    rows = []
    totals = [0, 0, 0, 0]
    for case in cases:
        if not case['schedules']:
            continue
        counts = [sum(r['opposite_only'] == opposite and
                      r['minimum_sequence_removals'] == threshold
                      for r in case['schedules'])
                  for opposite in (True, False) for threshold in (1, 2)]
        totals = [a+b for a,b in zip(totals, counts)]
        rows.append(f"{case['graph'].capitalize()} & {case['drones']} & " +
                    ' & '.join(map(str, counts)) + r' \\')
    rows += [r'\midrule', r'\multicolumn{2}{l}{Total} & ' +
             ' & '.join(map(str, totals)) + r' \\']
    (ROOT.parent/'figures/removal_robustness.tex').write_text('\n'.join(rows)+'\n')


if __name__ == '__main__':
    main()
