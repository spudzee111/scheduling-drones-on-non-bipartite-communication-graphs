# Scheduling Drones on Non-Bipartite Trajectory Graphs

Exact simulations, enumeration, and drone-removal experiments for Elizabeth
Cockerell's master's thesis, *Scheduling Drones on Non-Bipartite Trajectory Graphs*.

## Get the code

```sh
git clone https://github.com/spudzee111/scheduling-drones-on-non-bipartite-communication-graphs.git
cd scheduling-drones-on-non-bipartite-communication-graphs
```

The repository URL retains the earlier title so that the link in the submitted
thesis continues to work. Saved inputs and results are included; rerun the
commands below to reproduce them.

## Reproduce the results

Use Python 3.10 or later. The simulations, tests, and geometry checks use only
the Python standard library. Run these commands from this directory:

```sh
python3 -m unittest discover -s Drones/tests -v
python3 Drones/scripts/run_vse_experiments.py
python3 Drones/scripts/check_structural_examples.py
python3 Drones/scripts/run_removal_experiments.py
```

The first experiment command performs all 13 complete enumerations and replays
each collision-free connected class with an independent rational event
simulator. The expected totals across the specified graph/fleet-size pairs are:

| Result | Classes |
|---|---:|
| Collision-free with connected communication | 248 |
| Valid (also full coverage) | 246 |
| Opposite-orientation valid | 174 |
| Valid without shifts | 10 |
| Rejected for incomplete coverage | 2 |

The per-instance counts are in `Drones/results/vse/summary.json`; full schedule
records are in the adjacent JSON files. The run also writes the LaTeX table rows
to `figures/vse_enumeration.tex`. A simulation cutoff produces `unknown`, not a
proof of nonexistence. Every reported enumeration has `complete: true`.

The structural checks verify the six-circle obstruction and replay the
alternating-square family for 1 through 16 squares. The infinite-family result
is proved in the thesis, rather than inferred from those finite checks.

The removal search determines the minimum number of removals destroying
validity, with at least one survivor, for all 246 valid classes. It tests
every single removal at every grid phase. For the 90 classes surviving
every such removal, it verifies two removals at distinct times that destroy
validity. The remaining 156 classes have threshold one. Thus the reported
thresholds are exact even for arbitrary deletion sequences. The 1,120 residual
simulation replays agree with the independent rational simulator.

The detailed outcomes and witnesses are in
`Drones/results/vse/removal_robustness.json`; the generated table rows are in
`figures/removal_robustness.tex`. Each witness records a configuration, a
strictly positive wait to the removal phase, removed occupants, and the
resulting departure state. Phase zero is realized after a full period, not
by applying shifts to the schedule's initial configuration.

## Model and input conventions

- Unit circles; normalized angles and unit angular speed.
- Benchmark communication radius `2/5`; all linked center distances equal
  `12/5 = 2 + rho`. The distance threshold is inclusive.
- Input JSON edges are `[i, j, "angle"]`, where `i < j`, indices start at zero,
  and the bearing is measured on circle `i`. The reverse bearing is a half-turn
  away. The common grid denominator includes both bearings.
- Exact benchmark centers use `[a, b]` for `a + b*sqrt(3)`. The decimal `centers`
  field is for plotting. Geometry is checked with exact arithmetic.
- Schedule lists preserve drone identities by list position. Each readable state
  is `[circle, "angle", orientation]`, with `+1` counterclockwise and `-1`
  clockwise.
- Initial configurations are departure states. **No drone shifts at time zero.**
  Drones initially at paired endpoints communicate. Later arrivals apply the
  shifting rule once: communicate and stay if a partner arrives, otherwise
  transfer and reverse orientation.
- Collisions are rejected before shifts, including collisions halfway through
  a grid tick and simultaneous arrivals at one endpoint.
- Coverage includes arrival endpoints even when an immediate shift follows.

The enumerations count classes with **exactly** the specified number of drones.
They identify schedules under time translation and one fixed permutation of
drone identities, while retaining circle labels. The thesis's decision problem
asks for **at most** a given budget. A zero for an exact fleet size does not
imply a NO answer for that budget; smaller fleets must also be considered.

## Removal timing

Removal tests act immediately before arrivals are processed. A removal in
the open interval between two grid times has the same future as one made
just before the latter grid time: survivors cannot reach a link in that
interval. Removing just after an arrival is covered by the next grid time.
At each stage, the surviving fleet must still be nonempty. Communication
and coverage from before the last removal are not credited to the residual
schedule. A one-drone synchronization graph is connected, as in the thesis.

All reported thresholds are one or two. Their lower bounds require only
exhaustive single-removal tests; their upper bounds are supplied by explicit
one- or two-removal witnesses. No claim about larger thresholds is needed.

## Supplementary failure examples

The earlier two individual failure examples remain available separately:

```sh
python3 reproduction/verify_failure_examples.py
```

The baseline antipodal triangle starts at time zero without shifts. At time
`1`, it again reaches its initial angles. Remove the indicated drone immediately
before processing those arrivals. The remaining arrivals are then processed
once, and the residual schedule starts from the resulting departure states.

- Removing drone 1 (`u_1`) causes a collision at absolute time `5/4`, or `1/4`
  time units after failure, on circle 1 at angle `3/4`.
- Removing drone 4 (`v_2`) leaves five collision-free drones with full coverage,
  residual period `3`, and no further communications.

`reproduction/failure_examples.json` gives persistent one-based drone labels,
the arrival and departure states at failure, and both elapsed and absolute
collision times. Ordinary enumeration files use zero-based circle labels.

## Files

| Path | Purpose |
|---|---|
| `Drones/scs/model.py` | Rational graph, angular grid, arcs, and rings |
| `Drones/scs/dynamics.py` | Exact grid simulation and collision detection |
| `Drones/scs/enumeration.py` | Anchored enumeration and orbit equivalence |
| `Drones/scripts/run_vse_experiments.py` | Experiments and independent event replay |
| `Drones/scripts/check_structural_examples.py` | Structural example checks |
| `Drones/scripts/run_removal_experiments.py` | Removal thresholds and verified deletion sequences |
| `Drones/benchmarks/vse/` | Exact benchmark data |
| `Drones/results/vse/` | Recomputed counts and representatives |
| `Drones/tests/` | Simulation-convention and differential tests |
| `reproduction/verify_failure_examples.py` | Independent failure simulation |

For optional figure regeneration, install Matplotlib and NumPy, then run:

```sh
python3 Drones/scripts/make_thesis_figures.py
```

Runtime and platform metadata can differ between reproductions. Counts,
periods, coverage, and collision outcomes should agree exactly. The enumeration
summary contains SHA-256 hashes of its core sources and inputs.

This implementation stores traces and grid cells to reuse work on the small
benchmarks. The thesis's polynomial-space decision procedure instead streams
the simulation and records coverage on the original arcs. The experimental
implementation is not claimed to achieve that space bound.
