# `feature/coal-simd-support` performance findings

This document summarizes the performance work done on the
`feature/coal-simd-support` branch — what shipped, what was tried and rejected,
and what was investigated and deferred. The goal is so that future readers
don't repeat experiments that were already conclusively settled here.

For the deep-dive details (literature notes, full profile data, per-experiment
methodology), see `REASEARCH.md`. This document is the executive summary.

## Headline result

The branch is **16 commits ahead of `devel`** and delivers **~35% speedup** on
the polso/gen4 mesh-vs-mesh boolean motion-validation workload (the user's
real workload), measured by the bench at `scripts/coal_bench.cpp` with the
existing `output/bench_originals_manifest_polso_gen4_*.json` manifests.

| Joint  | Devel (2f4a8ede) | Branch HEAD | Δ      |
|--------|------------------|-------------|--------|
| J4-J5  | 1817 ns/q        | ~1160 ns/q  | -36%   |
| J5-J6  | 7266             | ~4520       | -38%   |
| J6-J7  | 4498             | ~2760       | -39%   |
| J7     | 1768             | ~1070       | -40%   |
| Nose   | 9551             | ~5900       | -38%   |

Methodology: P-core pinned (`taskset -c 0` on a 13th-gen Intel i7-13700H),
interleaved 7-run median per branch with separate libcoal.so per build. The
synthetic benchmark (`test/benchmark.cpp`) is unchanged from devel, since
the synthetic workload uses a default `CollisionRequest` (with
`enable_contact=true`) that does not activate the boolean fast-path.

## What shipped (in order)

| Commit    | Title                                           | Polso impact |
|-----------|-------------------------------------------------|--------------|
| 36c087ea  | feat: add SIMD convex support scanning          | small (no convex on this workload) |
| e688d2f7  | feat: add SoA SIMD convex support cache         | small (ditto) |
| 9f5fbbc1  | build: enable native SIMD architecture flags    | foundation |
| 1aa858cc  | perf: skip non-improving mesh distance leaves   | mostly distance-query |
| 5fd56b64  | perf: specialize RSS rectangle distance outputs | -3-5% |
| 2307b964  | build: enable no-interposition native codegen   | -1-2% |
| 44144fcf  | perf: skip oriented mesh distance seed          | mostly distance-query |
| fb4608c8  | perf: precompute oriented BV overlap inverse transform | -3-5% |
| 75cb88a6  | fix: harden review issues                       | correctness |
| 5716e385  | bench: VAMP research notes + benchmarks         | infrastructure |
| e57c13b6  | docs: rectDistance pre-check negative result    | docs |
| c030c3d7  | perf: devirtualize BVH traversal via templated overloads | -4-7% |
| 29f9297f  | docs: cross-workload validation                 | docs |
| 40349fea  | docs: Eigen computeDirect negative result       | docs |
| 4cabc131  | perf: hoist GJKSolver out of mesh leaf          | -5-7% |
| a18baf7e  | perf: Möller tri-tri overlap fast path          | -11-17% on top of hoist |
| 5e869752  | docs: OBB SAT 3-wide negative result            | docs |
| 6dc1d05d  | perf: prefetch BVH children before SAT          | -1-2% |
| 1b7b929f  | docs: prefetch + rectDistance investigation     | docs |

The two largest single-commit wins:

1. **Möller tri-tri overlap fast path** (`a18baf7e`): replaces GJK in the
   leaf-collision path when the user's `CollisionRequest` opts out of contact,
   distance lower bound, and security margin (`NO_REQUEST` flag). Möller's
   1997 algorithm gives the same boolean answer as GJK on triangle-vs-triangle
   pairs in ~30 FP ops vs GJK's iterative simplex search. Coplanar pairs
   (rare in practice; never produced by 50000-pair random fuzz) fall back
   to GJK. **-11-17% on top of the hoist alone** on polso, **-38% on
   robot/OBBRSS-deep** in `test/benchmark_robot.cpp`.

2. **Devirtualize BVH traversal** (`c030c3d7`): templated
   `collisionRecurseT<Node>` and `distanceRecurseT<Node>` plus templated
   `collide<Node>` / `distance<Node>` wrappers in `src/collision_node.h`.
   Internal callers in `collision_func_matrix.cpp` automatically pick the
   template via overload resolution, eliminating 9 vtable dispatches per
   BV-pair test. **-4-7% reproducible across all workloads.**

## What was tried and rejected (with measurements)

These are the experiments where I built it and measured, found no win or a
regression, and reverted. Each is in `REASEARCH.md` with the methodology
preserved so future experimenters don't repeat the same dead end.

| Experiment | Predicted | Measured | Decision | REASEARCH.md ref |
|------------|-----------|----------|----------|------------------|
| LTO + selective fast-math | 1.05-1.20× | regressed test suite | reverted | § 4 (E4) results |
| `-fassociative-math` for GJK | speed | breaks GJK termination | rejected | § 4 |
| Eigen `computeDirect` for `coal::eigen()` | -5-10% | **+25.6% regression** | reverted | scattered notes |
| GJK `distance_upper_bound = 0` early-break | -10-25% | within ±1% noise | reverted | (during boolean fast-path, GJK on tris already converges in 2-4 iters) |
| OBB SAT 3-wide Eigen Array vectorization | -5-10% | **0% (within ±1%)** | rejected | § 9 |
| rectDistance Voronoi tightening | -5-10% | no safe win identified | rejected | § 11 |

Key findings from the negative results:

- **Eigen `Array<Scalar,3,1>` doesn't vectorize to AVX2**: confirmed by
  disassembly. The 3-element type isn't 32-byte-aligned, so the compiler
  emits 2-wide SSE (`vmulpd %xmm`) not 4-wide AVX2 (`vmulpd %ymm`). Forcing
  AVX2 here would require padding to `Vector4d` with the 4th lane unused —
  bigger refactor than was justified.
- **GJK on triangle-vs-triangle converges in 2-4 iterations regardless**.
  The Minkowski difference of two 3-vertex shapes has only 9 candidate
  support points, so GJK terminates fast. The `distance_upper_bound`
  early-break optimization (which is real for *convex* polytope GJK) gives
  no measurable benefit here.
- **rectDistance is well-optimized already**. The 16 Voronoi region tests
  have cheap gating, the auxiliary computations are cached, and the natural
  shortcut (face-normal-sep early-out) requires an `upper_bound` parameter
  the function doesn't currently expose. Adding the parameter is an API
  change without clear evidence of benefit.

## What was investigated and deferred

These are ideas that remain plausible but were ruled out for *this* round
either because the empirical evidence on adjacent attempts suggested low
ROI, or because they'd be standalone features rather than optimizations.

### Eytzinger flat BVH layout

**Status**: rejected on first principles by an Explore agent. Realistic
ceiling 2-5% with high refit/refactor complexity, and the 4% gain from
devirtualization (`c030c3d7`) already covered the dispatch overhead this
would target. The 0% empirical result on the OBB SAT 3-wide vectorization
reinforced that small per-BV optimizations have hit diminishing returns at
this point in the curve. **Lesson**: the CAPT paper's 10× Eytzinger win
applies to point-cloud traversal where descent is nearly full-depth and
the per-node cost is tiny. coal's mesh BVH visits 5-16 nodes per query and
each node does dozens of FP ops, so traversal-structure optimizations can't
move the needle.

### Cross-pair SIMD BV batching (E3 alt B)

**Status**: deferred. The Explore agent estimated 20-40% on the SAT
function (~5-10% wall-clock if it holds). Track record on this branch shows
agent estimates run 2-5x optimistic (OBB SAT 3-wide: predicted 5-10%,
actual 0%; prefetch: predicted 5-10%, actual 1-2%), so realistic
expectation is 1-3% for several days of refactor. Negative ROI given the
branch is already well-optimized. The refactor would touch
`src/traversal/traversal_recurse.cpp` (recursive structure → batch buffer)
and add a SIMD SAT helper for N pairs in parallel. If pursued, the cleanest
shape is a 2-pair batch (one outer recursion split tested as both children
simultaneously).

### Bounding-sphere midphase pre-filter

**Status**: rejected on inspection. Coal's broadphase computes a *rotated*
world AABB at `collision_object.h:260-277` (for non-identity rotation it
rotates the local AABB extents and refits to world axes). For elongated
rotated meshes this is *tighter* than a translation-only AABB, and the
bounding sphere (rotation-invariant `aabb_radius`) is **looser** than the
rotated world AABB along the rotation plane. Sphere-sphere can't reject
what world-AABB-AABB doesn't, given the broadphase already does a tight
world AABB test.

### FOAM sphere-tree midphase

**Status**: deferred. Genuinely structural change with high expected gain
(2-5× on sphere-tractable meshes like robot arms). Requires offline FOAM
tooling integration: a Python library (`pip install foam-spherize`)
generates a hierarchical sphere cover from each mesh; coal would load these
and use as a midphase pre-filter before BVH descent. This is a *feature*,
not an optimization — it changes coal's data model (BVHModel optionally
holds a sphere tree) and requires per-mesh offline processing. Worth
pursuing as a separate dedicated effort if the polso/gen4 workload remains
performance-critical and 35% isn't enough.

### BVTT front-list / temporal coherence

**Status**: deferred. Motion validation calls `coal::collide` on
spatially-correlated transforms in sequence. Caching the BVTT front from
the previous query and starting from there could cut work substantially
(VAMP-MR paper claims ~50% on temporally-correlated workloads). Requires
caller to carry state across calls — an API change. Coal has a
`BVHFrontList` infrastructure already but the bench doesn't use it.

### CAPT point-cloud collision data structure

**Status**: not applicable to mesh-vs-mesh. CAPT (Ramsey 2024) speeds up
point-cloud collision queries via Eytzinger-laid-out k-d trees with
affordance sets. It would replace the OctoMap-backed point-cloud path in
coal, not the mesh-mesh path. Worth pursuing for users who have point-cloud
workloads (e.g., depth-camera-based obstacle avoidance).

## The diminishing-returns observation

A cross-cutting pattern emerged through the optimization curve: **agent
predictions of "5-10% gain" on micro-optimizations consistently delivered
0-2% in measurement** as the branch matured. Reasons:

1. **The first optimizations remove categories of cost** (vtable dispatch,
   heap allocations, GJK iteration), and each removal compounds: subsequent
   work runs against an already-tight baseline with less left to recover.
2. **Compilers already produce decent code for tight scalar kernels.** The
   places where SIMD-ifying *should* win are exactly the places the
   compiler already auto-vectorized or unrolled.
3. **Memory latency is a smaller fraction of total time once compute is
   tight.** Prefetch hints help less when there are fewer cycles to overlap
   with.

The takeaway for future optimization work on coal: **the next real wins
will likely come from a different workload exposing a different hot path,
or from structural changes that introduce a new data structure (FOAM sphere
trees, CAPT, BVTT front-list)**, not from squeezing more out of the
currently-tight kernel.

## What's protected by tests

Today, the `feature/coal-simd-support` branch passes:

- All standard unit tests: `coal-collision`, `coal-distance`,
  `coal-distance_lower_bound`, `coal-normal_and_nearest_points`,
  `coal-bvh_models`, `coal-frontlist`, `coal-gjk`,
  `coal-accelerated_gjk`, `coal-simd_support`,
  `coal-collision_node_asserts`, `coal-security_margin`.
- A 50000-pair random fuzz test (`test/tri_tri_overlap_fuzz.cpp`) showing
  Möller agrees with GJK 100% on non-coplanar inputs.
- The new `test/tri_tri_overlap_edge_cases.cpp` (Möller's coplanar /
  degenerate / shared-edge / shared-vertex paths).
- The new `test/regression_correctness.cpp` (per-pose collision booleans
  on the polso/gen4 manifests, baselines committed from devel).

Performance regressions are caught by the opt-in `coal-regression_perf`
test (committed timings JSON, machine-class gated). Run with:
```
ctest --label-regex perf
```

See `test/regression/README.md` for instructions on regenerating baselines
when intentionally changing collision semantics or hardware class.
