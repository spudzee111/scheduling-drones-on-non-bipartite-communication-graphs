"""Exact grid dynamics, including collisions strictly between grid times.

This is an exact finite transition system, not numerical time integration:
ports and initial states lie on 1/Q, and the only possible intervening
collision time is the half tick. All simultaneous arrivals are batched.
"""
from dataclasses import dataclass
from fractions import Fraction


@dataclass
class Result:
    status: str                 # periodic, collision, or unknown
    states: tuple              # labelled departure/post-shift snapshots
    syncs: tuple               # contacts at each snapshot, parallel to states
    shifts: int
    collision_tick: object = None
    collision_pair: object = None
    period_ticks: object = None

    def components(self):
        n = len(self.states[0]) if self.states else 0
        adj = [set() for _ in range(n)]
        for batch in self.syncs:
            for i, j in batch:
                adj[i].add(j)
                adj[j].add(i)
        groups, remaining = [], set(range(n))
        while remaining:
            start = remaining.pop()
            group, todo = {start}, [start]
            while todo:
                for j in adj[todo.pop()] & remaining:
                    remaining.remove(j)
                    group.add(j)
                    todo.append(j)
            groups.append(tuple(sorted(group)))
        return tuple(groups)

    @property
    def connected(self):
        return self.status == "periodic" and len(self.components()) == 1

    @property
    def starving(self):
        if self.status != "periodic":
            return None
        met = {i for batch in self.syncs for pair in batch for i in pair}
        return tuple(i for i in range(len(self.states[0])) if i not in met)

    @property
    def canonical(self):
        if self.status != "periodic":
            return None
        epochs = [tuple(sorted(s)) for s, sync in zip(self.states, self.syncs) if sync]
        return min(epochs) if epochs else min(tuple(sorted(s)) for s in self.states)


class Engine:
    def __init__(self, graph):
        self.g = graph
        q = graph.Q
        self.decoded = [graph.decode(i) for i in range(2*graph.n*q)]
        self.move = [2*(c*q+(a+s) % q)+(s == 1) for c, a, s in self.decoded]
        self.transfer = {2*p+b: 2*other+1-b
                         for p, other in graph.partner.items() for b in (0, 1)}

    def link_update(self, state):
        """Process arrivals once. Input must have distinct occupied positions.

        This is not applied to a schedule's initial state. It is also exposed
        for a removal immediately before a specified positive-time arrival.
        """
        occ = {code // 2: i for i, code in enumerate(state)}
        if len(occ) != len(state):
            raise ValueError("A collision must be rejected before link updates")
        post, sync, shifts = list(state), [], 0
        for i, code in enumerate(state):
            partner = self.g.partner.get(code // 2)
            if partner is None:
                continue
            j = occ.get(partner)
            if j is None:
                post[i] = self.transfer[code]
                shifts += 1
            elif i < j:
                sync.append((i, j))
        return tuple(post), tuple(sync), shifts

    def step(self, state):
        """Move one tick, reject collisions, then process the ending arrivals.

        Return (next departure state, contacts now, ending shifts,
        collision fraction of a tick, colliding pair). Initial link occupants
        communicate with a partner but never shift merely for starting there.
        """
        occ = {}
        for i, code in enumerate(state):
            p = code//2
            if p in occ:
                return None, (), 0, Fraction(0), (occ[p], i)
            occ[p] = i
        sync = []
        for i, code in enumerate(state):
            partner = self.g.partner.get(code//2)
            j = occ.get(partner)
            if j is not None and i < j:
                sync.append((i, j))
        # With no interior link events, relative speed is exactly 2 ticks/tick.
        first_collision, first_pair = None, None
        for i, code in enumerate(state):
            c, a, s = self.decoded[code]
            for j in range(i):
                d, b, r = self.decoded[state[j]]
                if c != d or s == r:
                    continue
                distance = ((b-a) if s == 1 else (a-b)) % self.g.Q
                if 0 < distance <= 2:
                    when = Fraction(distance, 2)
                    if first_collision is None or when < first_collision:
                        first_collision, first_pair = when, (j, i)
        if first_collision is not None:
            return None, tuple(sync), 0, first_collision, first_pair
        arrivals = tuple(self.move[x] for x in state)
        post, _, shifts = self.link_update(arrivals)
        return post, tuple(sync), shifts, None, None

    def simulate(self, initial, max_steps=100000):
        if not isinstance(max_steps, int) or max_steps < 1:
            raise ValueError("max_steps must be a positive integer")
        state = tuple(initial)
        if not state:
            raise ValueError("Simulate at least one surviving drone")
        if any(not isinstance(x, int) or not 0 <= x < len(self.move) for x in state):
            raise ValueError("Invalid state code")
        states, syncs, seen, shifts = [], [], {}, 0
        for t in range(max_steps):
            if state in seen:
                # The collision-free map is injective: a transient is a bug.
                if seen[state] != 0:
                    raise AssertionError("Unexpected transient in reversible dynamics")
                return Result("periodic", tuple(states), tuple(syncs), shifts,
                              period_ticks=t)
            seen[state] = t
            states.append(state)
            state, pairs, count, collision, pair = self.step(state)
            syncs.append(pairs)
            shifts += count
            if collision is not None:
                return Result("collision", tuple(states), tuple(syncs), shifts,
                              Fraction(t)+collision, pair)
        # Returning exactly at the bound is still a proved period.
        if state == tuple(initial):
            return Result("periodic", tuple(states), tuple(syncs), shifts,
                          period_ticks=max_steps)
        return Result("unknown", tuple(states), tuple(syncs), shifts)
