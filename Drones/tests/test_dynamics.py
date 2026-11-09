"""Regression tests for the thesis's departure-state simulation convention."""
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))
from scs.model import Graph
from scs.dynamics import Engine
from scs.enumeration import candidates
from run_vse_experiments import grid_properties, reference_properties, check_geometry


class DynamicsTests(unittest.TestCase):
    def test_unmatched_initial_endpoint_does_not_shift(self):
        g = Graph(2, [(0, 1, '0')])
        e = Engine(g)
        x = (g.encode(0, 0, 1),)
        y, pairs, shifts, collision, _ = e.step(x)
        self.assertEqual(g.describe_state(y), [[0, '1/2', 1]])
        self.assertEqual((pairs, shifts, collision), ((), 0, None))
        z, _, shifts, _, _ = e.step(y)
        self.assertEqual(g.describe_state(z), [[1, '1/2', -1]])
        self.assertEqual(shifts, 1)

    def test_initial_partners_communicate(self):
        g = Graph(2, [(0, 1, '0')])
        e = Engine(g)
        x = (g.encode(0, 0, 1), g.encode(1, F(1, 2), 1))
        result = e.simulate(x)
        self.assertEqual(result.syncs[0], ((0, 1),))
        self.assertTrue(result.connected)
        self.assertEqual(result.shifts, 0)
        self.assertEqual(result.period_ticks, 2)

    def test_collisions_before_any_shift(self):
        g = Graph.read(ROOT / 'benchmarks/vse/square.json')
        e = Engine(g)
        for angle, expected in [(F(1, 4), F(1, 2)), (F(1, 2), F(1))]:
            x = (g.encode(0, 0, 1), g.encode(0, angle, -1))
            result = e.simulate(x)
            self.assertEqual(result.status, 'collision')
            self.assertEqual(result.collision_tick, expected)
            self.assertEqual(result.shifts, 0)

    def test_shift_does_not_trigger_another_shift(self):
        g = Graph(2, [(0, 1, '0')])
        result = Engine(g).simulate((g.encode(0, 0, 1),))
        self.assertEqual(result.status, 'periodic')
        self.assertEqual(result.period_ticks, 4)
        self.assertEqual(result.shifts, 2)

    def test_cutoff_is_unknown_and_exact_bound_is_periodic(self):
        g = Graph(2, [(0, 1, '0')])
        e = Engine(g)
        x = (g.encode(0, 0, 1),)
        self.assertEqual(e.simulate(x, 3).status, 'unknown')
        self.assertEqual(e.simulate(x, 4).status, 'periodic')

    def test_geometry_includes_links_at_threshold(self):
        for name in ('triangle', 'square', 'diamond'):
            check_geometry(ROOT / f'benchmarks/vse/{name}.json')

    def test_independent_replay_of_small_anchored_states(self):
        for name in ('square', 'triangle'):
            g = Graph.read(ROOT / f'benchmarks/vse/{name}.json')
            e = Engine(g)
            for size in (2, 3):
                for x in set(candidates(g, size)):
                    result = e.simulate(x)
                    self.assertNotEqual(result.status, 'unknown')
                    if result.status == 'collision':
                        with self.assertRaisesRegex(AssertionError, 'Collision'):
                            reference_properties(g, x)
                        continue
                    ref = reference_properties(g, x)
                    missing, same = grid_properties(g, result)
                    self.assertEqual(ref['period'], F(result.period_ticks, g.Q))
                    self.assertEqual(ref['syncs'], {p for batch in result.syncs for p in batch})
                    self.assertEqual(ref['opposite_only'], not same)
                    self.assertEqual(ref['covered_lengths'],
                                     [F(g.Q-sum(c == i for c, a in missing), g.Q)
                                      for i in range(g.n)])


if __name__ == '__main__':
    unittest.main()
