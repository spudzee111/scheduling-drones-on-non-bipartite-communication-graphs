"""Independent checks of failure timing, surviving coverage, and witnesses."""
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'scripts'))
from scs.model import Graph
from run_removal_experiments import RemovalSearch


class RemovalTests(unittest.TestCase):
    def test_one_drone_is_retained_and_singleton_communication_is_connected(self):
        g = Graph(2, [(0, 1, '0')])
        search = RemovalSearch(g)
        # With two equally oriented drones, either one alone covers the ring.
        result = search.analyze((g.encode(0, 0, 1), g.encode(1, '1/2', 1)))
        self.assertIsNone(result['minimum_sequence_removals'])
        self.assertIsNone(result['minimum_simultaneous_removals'])

    def test_square_opposite_schedule_needs_two_removals(self):
        g = Graph.read(ROOT/'benchmarks/vse/square.json')
        search = RemovalSearch(g)
        result = search.analyze((0, 13, 17, 28))
        self.assertEqual(result['minimum_sequence_removals'], 2)
        self.assertEqual(result['single_removal_counts']['scenarios'],
                         result['single_removal_counts']['valid'])
        self.assertEqual(len(result['sequence_witness']), 2)
        for index, step in enumerate(result['sequence_witness']):
            self.assertEqual(len(step['removed_indices']), 1)
            after = search.inspect(step['surviving_departure'])
            self.assertEqual(after['valid'], index == 0)

    def test_triangle_antipodal_schedule_can_collide_after_one_removal(self):
        g = Graph.read(ROOT/'benchmarks/vse/triangle.json')
        search = RemovalSearch(g)
        state = tuple(g.encode(c, angle, 1) for c in range(3)
                      for angle in ('0', '1/2'))
        result = search.analyze(state)
        self.assertEqual(result['minimum_sequence_removals'], 1)
        self.assertGreater(result['single_removal_counts']['collision'], 0)

    def test_removal_after_a_tick_matches_before_the_next_update(self):
        g = Graph.read(ROOT/'benchmarks/vse/square.json')
        search = RemovalSearch(g)
        result = search.engine.simulate((0, 13, 17, 28))
        for state in result.states:
            for removed in range(len(state)):
                after_now = tuple(x for i,x in enumerate(state) if i != removed)
                next_arrivals = tuple(search.engine.move[x] for x in after_now)
                after_next, _, _ = search.engine.link_update(next_arrivals)
                self.assertEqual(search.inspect(after_now), search.inspect(after_next))


if __name__ == '__main__':
    unittest.main()
