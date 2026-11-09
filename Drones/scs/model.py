"""Rational link geometry. Floating point is never used for decisions."""
from dataclasses import dataclass
from fractions import Fraction
from math import lcm
from pathlib import Path
import json


@dataclass(frozen=True)
class Arc:
    circle: int
    start: int
    end: int
    length: int


class Graph:
    """Connected simple graph, distinct ports per circle, antipodal partners.

    A state code is 2 * (circle * Q + angle_tick) + (orientation == +1).
    State tuples keep drone identities by tuple index. Saved input snapshots
    are departure states: no shift is applied at time zero.
    """

    def __init__(self, n, edges, name="graph", centers=None, metadata=None):
        if not isinstance(n, int) or n < 2:
            raise ValueError("At least two trajectories are required")
        edges = list(edges)
        if any(not isinstance(i, int) or not isinstance(j, int) or isinstance(a, float)
               for i, j, a in edges):
            raise ValueError("Use integer vertex indices and exact rational angles, not floats")
        raw = [(i, j, Fraction(a) % 1) for i, j, a in edges]
        if any(not 0 <= i < j < n for i, j, _ in raw):
            raise ValueError("Edges must satisfy 0 <= i < j < n")
        if len({(i, j) for i, j, _ in raw}) != len(raw):
            raise ValueError("Parallel edges are not supported")
        self.name, self.n, self.centers = name, n, centers
        self.metadata = metadata or {}
        self.Q = lcm(2, *(a.denominator for _, _, a in raw))
        q = self.Q
        self.edges = tuple((i, j, int(a*q)) for i, j, a in raw)
        self.ports = [dict() for _ in range(n)]
        self.partner = {}
        for e, (i, j, a) in enumerate(self.edges):
            b = (a+q//2) % q
            for c, t, d, u in [(i, a, j, b), (j, b, i, a)]:
                if t in self.ports[c]:
                    raise ValueError(f"Two links use the same port on circle {c}")
                self.ports[c][t] = e
                self.partner[c*q+t] = d*q+u
        seen, todo = {0}, [0]
        while todo:
            c = todo.pop()
            for a in self.ports[c]:
                d = self.partner[c*q+a]//q
                if d not in seen:
                    seen.add(d)
                    todo.append(d)
        if len(seen) != n:
            raise ValueError("The communication graph must be connected")
        self.arcs, self.circle_arcs, self.arc_after = [], [], {}
        for c in range(n):
            angles = sorted(self.ports[c])
            ids = []
            for a, b in zip(angles, angles[1:]+angles[:1]):
                idx = len(self.arcs)
                self.arcs.append(Arc(c, a, b, (b-a) % q or q))
                self.arc_after[c, a, 1] = idx
                self.arc_after[c, b, -1] = idx
                ids.append(idx)
            self.circle_arcs.append(ids)
        self.rings, self.arc_ring, self.arc_sign = self._rings()
        self.ring_lengths = [sum(self.arcs[a].length for a in r) for r in self.rings]

    @classmethod
    def read(cls, path):
        data = json.loads(Path(path).read_text())
        return cls(data["n"], data["edges"], data.get("name", Path(path).stem),
                   data.get("centers"), data.get("metadata"))

    @classmethod
    def read_legacy(cls, path):
        """Drones3/4/6 forward-edge format; blank lines retain circle IDs."""
        edges, largest = [], 0
        lines = Path(path).read_text().splitlines()
        for i, line in enumerate(lines):
            largest = max(largest, i)
            for entry in line.split():
                j, a = entry.split(":")
                j = int(j)
                edges.append((i, j, a))
                largest = max(largest, i, j)
        return cls(largest+1, edges, Path(path).stem)

    def encode(self, circle, angle, direction):
        if direction not in (-1, 1) or not 0 <= circle < self.n:
            raise ValueError("Invalid circle or orientation")
        tick = Fraction(angle)*self.Q
        if tick.denominator != 1:
            raise ValueError("Angle is not on this graph's grid")
        return 2*(circle*self.Q+int(tick) % self.Q)+(direction == 1)

    def decode(self, code):
        p, s = divmod(code, 2)
        c, a = divmod(p, self.Q)
        return c, a, 2*s-1

    def describe_state(self, state):
        return [[c, str(Fraction(a, self.Q)), s] for c, a, s in map(self.decode, state)]

    def read_state(self, assignments):
        return tuple(self.encode(c, a, s) for c, a, s in assignments)

    def arc_at(self, code, outgoing=True):
        c, a, s = self.decode(code)
        if a in self.ports[c]:
            return self.arc_after[c, a, s if outgoing else -s]
        return next(i for i in self.circle_arcs[c]
                    if (a-self.arcs[i].start) % self.Q < self.arcs[i].length)

    def _rings(self):
        rings, arc_ring, signs = [], {}, {}
        for first in range(len(self.arcs)):
            if first in arc_ring:
                continue
            a, s, ring = first, 1, []
            while a not in arc_ring:
                arc_ring[a], signs[a] = len(rings), s
                ring.append(a)
                arc = self.arcs[a]
                endpoint = arc.end if s == 1 else arc.start
                p = self.partner[arc.circle*self.Q+endpoint]
                c, t = divmod(p, self.Q)
                s = -s
                a = self.arc_after[c, t, s]
            if (a, s) != (first, 1):
                raise AssertionError("Ring traversal did not close consistently")
            rings.append(tuple(ring))
        return tuple(rings), arc_ring, signs

    def is_bipartite(self):
        color, todo = {0: 0}, [0]
        while todo:
            c = todo.pop()
            for a in self.ports[c]:
                d = self.partner[c*self.Q+a]//self.Q
                if d in color and color[d] == color[c]:
                    return False
                if d not in color:
                    color[d] = 1-color[c]
                    todo.append(d)
        return True
