# VAMP → coal optimization research

This document captures the literature on VAMP and adjacent work in enough depth that we never need to re-read the papers, plus an experiment plan for porting promising techniques into coal. Status as of 2026-04-29.

## Operating constraints (from project owner)

1. **No public-API changes.** All optimizations must be internal — same headers, same function signatures.
2. **Empirical only.** A change ships only if it produces a *measurable* speedup on the existing benchmark suite. No speculative rewrites.
3. **AVX2 only for now.** Skip ARM NEON / WASM SIMD; revisit later if useful.
4. **No robot-specialized geometry.** No FOAM-spherized robot meshes, no specialized sphere covers. Benchmarks must use the existing env/rob meshes from `coal-test-benchmark` and synthetic generic geometry.
5. **`perf` is unavailable** in the dev environment (`/proc/sys/kernel/perf_event_paranoid=4`, no passwordless sudo). Fall back to wall-clock benchmarks. Revisit if `sudo sysctl kernel.perf_event_paranoid=1` becomes possible.

Implications for the experiment list in § 3.2/§ 4:

- **E7 (sphere-tree midphase via FOAM): dropped** — relies on robot-specialized sphere covers per constraint 4.
- **E1 (CAPT): kept but reframed.** CAPT can be implemented as a generic point-cloud collision primitive against existing `BVHModel`-of-spheres or octomap-style inputs. We do *not* generate the input clouds from spherized robots; we use synthetic clouds (uniform random surface samples of the existing benchmark meshes) for measurement.
- **E5 (convex hill-climb), E6 (kDOP SIMD), E9 (AABB SIMD distance), E10 (vectorized tri-tri):** unchanged, all generic.
- **E2 (Eytzinger BVH), E3 (SIMD SAT batching), E4 (LTO + fast-math), E8 (sphere-AABB SIMD prefilter):** unchanged.

## 0. Reading-order TL;DR

VAMP solves a fundamentally different problem than coal:

- **VAMP** = sampling-based motion planner. Batches *N configurations* of one robot, evaluates FK + collision via SIMD lanes, rejects the whole batch on first hit. Robot is approximated as **spheres only**; environment as primitive spheres / capsules / cuboids / cylinders / heightfields / point clouds. No mesh narrowphase. No GJK/EPA. No exact distance, only boolean overlap.
- **coal** = general-purpose exact collision/distance library (HPP-FCL fork). Single query, exact triangle-mesh narrowphase via OBB/RSS/kIOS BVH, GJK+EPA for convex shapes, returns distance and contact information.

This asymmetry means **most VAMP wins do not port directly**. The high-impact wins ("rake" SIMD over configurations; sphere-only robot; tracing-compiled FK) are tied to motion planning and would require coal to grow a planner-flavored API. The portable wins are mostly *machine-sympathetic data layout and codegen patterns*: Eytzinger-layout BVH, branch-free traversal, fast-math flags, AVX2 batching of independent BV tests, and a CAPT-like point-cloud structure to replace OctoMap-backed queries.

## 1. Literature distilled

### 1.1 VAMP (Thomason et al., ICRA 2024) — primary paper

**Headline number:** 35 µs median planning time for 7-DOF Panda over MotionBenchMaker on one core of a Ryzen 9 7950X — > 500× over PyBullet/OMPL, ~100–200× over MoveIt/OMPL. Achieves 25 kHz median planning rate, 10 kHz mean for Panda. ARM Cortex-A76 (Orange Pi 5B) still 20–50× faster than desktop baselines.

**Architecture:**

- AVX2 (8-wide float) is primary; ARM NEON; WASM SIMD.  AVX-512 was *tried and was slower* — attributed to downclocking and lack of 512-bit registers on consumer CPUs. **Implication for us: don't chase AVX-512 on Ryzen-class hardware unconditionally.**
- SoA throughout. Each batch of 8 configurations stored as `FloatVector<8, n_dofs>` (8 lanes × n DOFs). Result of FK: `FloatVector<8, n_spheres>` separately for x, y, z, r — never AoS `Vec3`.
- "Rake": within a single edge-validation, the 8 configurations checked in parallel are *spatially distributed* across the discretized motion (`q[0], q[n/8], q[2n/8], ...`). Probability of hitting any collision early goes up dramatically because nearby states are correlated. After the initial rake passes, "comb through" remaining states by stepping each lane forward by 1.

**Tracing compiler for FK** (key contribution):

- Reads URDF, traces the operations of `forward_kinematics(q) → sphere_positions`, emits *straight-line* C++ with no joint-type polymorphism, no dynamic branching, no spurious data dependencies between link transforms. Constant-folds DH params, sin/cos of fixed angles, etc.
- Per-sphere FK computations are *interleaved* with per-sphere collision checks in the generated code, so an early-out collision skips computing later spheres' positions. Comparable to CRTP-static-dispatch (e.g. Pinocchio's [Carpentier 2019]) but emits even leaner code because precise-operation tracking removes ops that aren't statically detectable.
- Compiler is ~3000 LOC of Python; generated headers are ~600+ floats and 14k operations for Baxter (75 spheres, 14 DOF).

**Collision tests** (all 8-wide AVX2):

- `sphere_sphere_sql2(a, b)`: squared distance `(ax−bx)² + (ay−by)² + (az−bz)²`, then compare against `(ra+rb)²`. **No sqrt anywhere in inner loop.**
- `sphere_capsule`: project sphere center onto capsule axis with `clamp(dot · rdv, 0, 1)`; one squared-distance test. Uses precomputed `rdv = 1 / |axis|²` cached in shape struct.
- `sphere_cuboid`: `dot_3(p, axis_i) - axis_i_r` per axis, `abs()`, `max(., 0)`, sum-of-squares. Branchless via SIMD `_mm256_min/max_ps`. Uses `_mm256_blendv_ps`-style mask blends.
- Sphere-heightfield: grid lookup + interpolation, 8-wide.
- Sphere-pointcloud: delegated to **CAPT** (see § 1.2).

**Mid-phase / broad-phase:**

- VAMP *deliberately abandons BVHs* for robot-vs-environment. Instead, a **hierarchy of sphere models**: e.g., one bounding sphere per link at the coarsest level, refined to ~75 spheres at the finest. Each level is checked branch-free against the obstacle set; pruning happens by skipping refinement. This avoids the recursive, branch-heavy memory access of OBB-trees.
- Environment obstacles sorted at init by `min_distance` to robot base — closest tested first; helps early-out.
- No precomputed obstacle inverse transforms; per-query transforms applied lane-wise.

**Compiler flags** (`vamp/cmake/CompilerSettings.cmake`):

- `-O3 -fno-math-errno -fno-signed-zeros -fno-trapping-math -fno-rounding-math -ffp-contract=fast -fassociative-math` (always, x86 only)
- Clang ≥ 12 also: `-fapprox-func -fno-honor-infinities -fno-honor-nans`
- `-march=native -mavx2` (x86) or `-mcpu=native -mtune=native` (ARM); `-flax-vector-conversions` on GCC ≥ 13.
- LTO enabled by default (`-flto=auto` on GCC, `-flto` on Clang).
- Paper's caveat (footnote 9): "Some issues with `-ffast-math` (handling non-finite values, subnormals) are not particularly relevant for motion planning, where configurations are bounded… we warn practitioners to be wary of issues from reciprocal approximation."

**What VAMP critiques about FCL/HPP-FCL** (relevant to us, paraphrased from § II and § 2.3 of VAMP-MR):

> Recursive BVH traversal introduces significant conditional branching and irregular workloads, complicating efficient task distribution… High memory bandwidth required to access BVH transforms and geometry data often becomes a limiting factor without careful data layout and caching strategies.

This is essentially coal's current shape. The complaints are real — see § 2 below for what coal has and hasn't fixed.

**Key future-work limitations called out by VAMP (relevant to porting):**

- Sphere-only collision is conservative; real geometries differ. Resolved partly by FOAM (Coumar 2025) which generates principled sphere covers from meshes via Adaptive Medial Axis Approximation (Bradshaw & O'Sullivan 2004).
- Mesh-on-mesh collision *not addressed* — explicitly outside VAMP's scope. coal is the natural place to host that.
- Distance / signed-distance / contact info also not addressed. coal already has this.

### 1.2 CAPT (Ramsey et al., RSS 2024) — point-cloud collision data structure

**Headline numbers:** mean query time **9.89 ns** on point clouds up to 50k points; ~10× faster than nanoflann (309 ns), ~1000× faster than OctoMap (10 µs). Mean construction time 4-6 ms on Panda/UR5/Fetch scenes after filtering.

**Data structure:**

- Tuple `(T, A, P)` over `n` points (n padded to power of two with `+∞`).
- `T` is `n−1` median-split test values, **stored in Eytzinger layout** (the array layout used for binary heaps: `T[0]` is root, `T[1]` and `T[2]` are children, recursively `T[2i+1]` left / `T[2i+2]` right). Branch-free traversal: at depth d on dimension `d mod k`, compare `x[d mod k]` against `T[i]`, set `i ← 2i + 1 + (x[d] > T[i])`. Loop runs *exactly* `log₂ n` iterations regardless of input — no termination test.
- `A` is `n` axis-aligned bounding boxes, one per leaf, of the leaf's *affordance set*.
- `P` is `n` ragged arrays — each leaf's affordance set, padded to SIMD-width multiples and stored SoA so a SIMD reduction over distances becomes a `_mm256_loadu_ps + sub + sqr + add` chain.

**Affordance set** (the trick): each leaf cell `c` stores not just its representative point but **all points in the cloud that could be hit by any query sphere whose center lies in `c`** at the maximum query radius `r_max`. This duplicates points across leaves but eliminates the back-tracking phase that kills KD-tree SIMD performance. Lemma V.1 proves correctness: a sphere `s ⊂ c` of radius `r ∈ [r_min, r_max]` collides with the cloud iff it collides with the leaf's affordance set.

**Construction** (Alg. 1): standard recursive median-split (quickselect for O(n log n) expected). At each split, propagate the current "outside-cell-but-near" affordance set into both children, intersect with the new cell's afford-relation. Cells small enough that any in-cell sphere of radius `r_min` strictly contains the rep-point store only that single point (Fig. 2b optimization).

**Query** (Alg. 3, branch-free):
```
i, d = 0
while i < n−1:
    i = 2i + 1 + (x[d mod k] > T[i])
    d += 1
# Now i−n+1 is the leaf index.
if not s.intersects(A[i−n+1]): return false
for p in P[i−n+1]:
    if dist(x, p) ≤ r: return true
return false
```
Two parallelism modes:
- **Across queries** (multiple spheres tested together) for the traversal and AABB rejection.
- **Across affordance points** (within one leaf) for the final exhaustive distance test — different parallelism axis than traversal because the access pattern is contiguous within a leaf.

**Filtering with Z-order curves** (Alg. 2): for dense input clouds, collapse points within `r_filter` of each other by sorting along a Morton (Z-order) curve and walking the sort, keeping only points farther than `r_filter` from the running last-kept. Repeated for each of the 6 axis permutations to handle adjacency artefacts. Lemma V.2 bounds the introduced gap: `r_filter ≤ r_min − δ(O, P_C)` keeps plans valid (or pad query radii by `r_filter` for any value).

**Limitations:**
- Immutable. Add/remove requires rebuild. Fine for sensor-frame data (rebuild per frame); not fine for incremental updates.
- O(k·n²) worst-case memory for very-low-dispersion clouds (every point is afforded by every cell).
- Authors themselves suggest hybrid CAPT + spatial-hash for very large clouds.

### 1.3 VAMP-MR (Huang et al., AAAI 2026) — multi-robot extension

**Headline numbers:** 11–28× speedup on single-config collision check vs FCL with the same spherized geometry; 65–148× on motion validation; 100% RRT-Connect success vs 78–96% with FCL on hard scenes; cache miss rate reduced 59–86% vs FCL.

**Method (Alg. 1, FKCC_MULTI):**

- Per-robot: FK → spheres in robot frame → self-CC → environment-CC against obstacles transformed once into robot base frame (cached because base xforms and obstacles are fixed during planning).
- Inter-robot: pairwise sphere-vs-sphere check between robots, in world frame.
- All steps: 8-wide rake over configurations.

**Insights relevant to coal:**

- "Cache transformed environment obstacles in each robot's base frame" — same idea as our `overlapPrecomputedRTranspose`, but done once per init rather than per BV node.
- "We discard traditional broad-phase collision checking and instead opt for a 'rake' strategy" — confirms the rake's effectiveness as a *replacement* for BVH for this workload.
- Confirms that compiling with `-march=native -mavx2 -O3` + GCC 9 already gives most of the speedup; specialized SIMD intrinsics layered on top.

### 1.4 FOAM (Coumar et al., 2025) — sphere generation

URDF → spheres tool, based on **Adaptive Medial Axis Approximation** (Bradshaw & O'Sullivan 2004). Inputs: any mesh, possibly with defects. Outputs: a hierarchical sphere tree with configurable level-of-detail. This is the upstream tool that makes sphere-only robot collision tractable for arbitrary URDFs. Open-source Python library.

**Implication for coal:** if we build a sphere-tree midphase for `BVHModel`, FOAM is the obvious tool to generate the sphere tree from input meshes — we don't need to invent the sphere-fitter.

## 2. coal's current optimization status

(After the recent SIMD branch including review fixups committed at `75cb88a6`.)

| Component | Status | Detail |
|---|---|---|
| AVX2 max-dot for convex support (linear search) | **Done** | `simd::maxDot` / `maxDotSoA`, 8-wide float / 4-wide double, scalar fallback. |
| SoA cache for convex points | **Done** | `support_points_x/y/z` rebuilt in `initialize`/`deepcopy`/serialization-load. |
| OBB/RSS/OBBRSS/kIOS overlap precomputed transpose | **Done** | One transpose + mat-vec saved per BV node test. |
| RSS rectDistance specialization (no points) | **Done** | `rectDistanceImpl<bool ComputePoints>`. |
| Mesh-distance leaf early-out + seed removal | **Done** | Guarded against `min_distance < 0` and seed `Scalar::max` per review. |
| Fast-math flags | **None** | No `-ffast-math`, `-ffp-contract=fast`, etc. — only `-march=native -mavx2 -fno-semantic-interposition`. |
| LTO | **Off by default** | No CMake plumbing for `-flto`. |
| GJK simplex inner loop | **Scalar** | Eigen-vectorized but no coal-specific SIMD. |
| OBB SAT (15 axes) | **Scalar** | `obbDisjointAndLowerBoundDistance` is straight scalar loops with `Bf` (precomputed `|R|`). |
| RSS rect-distance Voronoi region tests | **Scalar** | Heavy region-conditional branching. |
| BVH traversal | **Pointer-recursive, AoS** | `traversal_recurse.cpp` uses `std::vector<BVPair>` stack, depth-first, no branch-free walk, no SIMD batching of pair tests. |
| Triangle-triangle distance | **Scalar** (delegated to GJK) | `TriangleP` distance via `triangle_triangle.cpp` ⇒ GJK. |
| AABB tests | **Eigen-auto-vectorized** | `BV/AABB.cpp` `overlap` uses `array().max()`; `distance()` is scalar 3-axis loop. |
| kDOP dot products | **Scalar** | Templates for k=5/6/9 each compute scalar combinations; obvious 4× SIMD candidate. |
| Convex hill-climb (`getShapeSupportLog`) | **Scalar** | Neighbor walk evaluates one vertex at a time, then routes the warm-start refresh through `simd::maxDot` once. |
| Memory layout for BVH nodes | **AoS, no align hints** | `BVNode` bundles bv + structural data; default Eigen 16-byte align. |
| Eigen AVX2 macros | **Inherited from `-mavx2`** | `EIGEN_VECTORIZE_AVX2` is auto-detected; coal does not set it explicitly. |
| Runtime CPU dispatch | **Per-call function-pointer-style** | `isAvx2Enabled()` checked once via `static const`. No CRTP dispatch on the hot mesh-mesh path. |

## 3. Gap analysis

### 3.1 VAMP wins that **do not port** (and why)

- **"Rake" 8-wide configuration batching** — coal has no notion of "configuration." It evaluates one (T1, T2) collision query at a time. Adding a batched API (`collide(N transforms vs N transforms)`) is doable but is essentially a new top-level API, not a hot-path optimization. *Out of scope for this round.*
- **Tracing-compiled FK** — coal does not do FK. Pinocchio is the upstream choice in this stack.
- **Sphere-only robot model** — coal supports it via `BVHModel<...>` of sphere shapes, but the speed-up only manifests when the *user* opts in. Not an internal optimization.
- **Discarding broad-phase entirely** — only viable when geometry is trivially fast (sphere-vs-primitive). With triangle meshes, you still want the BVH. So we'd be replacing a good BVH with a cheap-to-traverse-but-bad BVH. Bad trade.
- **`-ffp-contract=fast` + `-fassociative-math` for GJK termination** — risky for coal because GJK convergence is *defined by the cumulative numerical error of the simplex*. Different reductions can flip the termination decision and break tests. Worth investigating but high risk.

### 3.2 VAMP wins that **port well** (ranked by expected ROI)

These are the experiment candidates. Each is fleshed out as a numbered experiment in § 4.

| # | Idea | Expected speedup | Risk | Effort |
|---|---|---|---|---|
| E1 | CAPT-style point-cloud structure replacing OctoMap path | 10-100× on PC scenes | Low (additive) | 2-3 weeks |
| E2 | Eytzinger / branch-free flat BVH traversal | 1.3-1.8× on mesh-mesh BVH | Medium (touches hot core) | 2 weeks |
| E3 | SIMD-batched OBB SAT (8 BV pairs at once) | 1.5-2.0× on mesh-mesh collision | Medium-high (numerical) | 2 weeks |
| E4 | LTO + selective fast-math flags audit | 1.05-1.2× across the board | Low–Medium (test fragility) | 3-4 days |
| E5 | Vectorized hill-climb in `getShapeSupportLog` (warm-start scan and neighbor batch) | 1.2-1.5× on convex-vs-convex | Low | 1 week |
| E6 | kDOP SIMD dot-products (4× per axis) | 1.3-2× on kDOP-using BVHs | Low | 2-3 days |
| E7 | ~~Sphere-tree midphase~~ — **DROPPED** per constraint 4 | — | — | — |
| E8 | AVX2 sphere-vs-AABB primitive used as broadphase prefilter | 1.1-1.3× on AABB-heavy queries | Low | 4-5 days |
| E9 | AABB SIMD distance (replace 3-axis scalar loop) | 1.05-1.10× | Low | 1-2 days |
| E10 | Vectorized triangle-triangle distance (avoid GJK delegation for tris) | 1.5-2× on tri-tri leafs | Medium-high (correctness) | 2-3 weeks |

A note on ordering: **E4 (build flags) is cheap and informative — do it first** because it changes the baseline numbers all subsequent experiments are measured against. **E1 (CAPT) is highest absolute value** but additive, so it can run in parallel with the rest.

## 4. Experiment specs

Format per experiment: hypothesis · method · success criteria · files to touch · benchmark workload · risks.

### E1 — CAPT for point-cloud collision

**Hypothesis.** Replacing coal's current OcTree-backed point-cloud query with a CAPT will yield 10×+ throughput on dense point clouds, and (if the workload exists) bring point-cloud collision down to comparable cost as primitive-vs-primitive.

**Method.**
1. Build a new `CollisionGeometry` subclass `CAPTGeometry` storing `(T[], A[], P[])` in Eytzinger order. Reuse coal's `Vec3s` SoA pattern (we have it for convex).
2. Implement Alg. 1 (construct) using std::nth_element for the median split, the `r_min`/`r_max` bounds passed via constructor.
3. Implement Alg. 2 (Z-order filter) as a free function; use coal's existing `Vec3s` types.
4. Implement Alg. 3 (query), starting with scalar form to verify correctness against brute force, then add the AVX2 batched form.
5. Wire into `collide(...)` dispatch: when `geom2` is `CAPTGeometry` and `geom1` is a sphere/sphere-tree, use the CAPT path; otherwise fall back to existing.
6. Optional: filter algorithm in Python utility (FOAM-style) to populate test scenes.

**Success criteria.**
- Correctness: 100% match vs brute-force collision check across MotionBenchMaker-style scenes.
- Performance: ≥ 10× over coal's current OctoMap path on a 50k-point synthetic cloud.
- Build time: < 10 ms per filtered cloud of 5k points (matches CAPT paper).

**Files.**
- New: `include/coal/collision/capt.h`, `src/collision/capt.cpp`, `test/capt.cpp`, `test/capt_filter.cpp`.
- Modify: `src/collision_func_matrix.cpp` (dispatch), `src/distance_func_matrix.cpp`, possibly `include/coal/collision_object.h` (new `OT_CAPT` object type).

**Workload.** MotionBenchMaker scenes converted to point clouds via uniform surface sampling at 50k points, filtered to ~5k points with `r_filter = 2 cm`. Run 10k random sphere queries.

**Risks.** Incremental updates not supported (matches CAPT paper limitation); we'd document this. Worst-case memory on low-dispersion clouds is O(n²) — must enforce filter step.

---

### E2 — Eytzinger / flat BVH traversal

**Hypothesis.** Converting `traversal_recurse.cpp`'s pointer-walking depth-first traversal to a flat-array Eytzinger-laid-out BVH with branch-free index updates will reduce L1 misses and branch mispredictions enough to give 30-80% speedup on mesh-mesh collision/distance.

**Method.**
1. Add a `flatten()` method to `BVHModel<BV>` that produces a parallel Eytzinger-ordered array of `BVNode<BV>` (interior + leaf). Original tree retained for backward compat.
2. Add `MeshCollisionTraversalNodeFlat<BV, RTIsIdentity>` that traverses the flat array with a small fixed-size stack (no `std::vector` allocations).
3. Use the same `overlapPrecomputedRTranspose` work we did before — should be unchanged.
4. Optional: branchless update on the "which child first" heuristic — current code calls `firstOverSecond()` to pick the larger BV; in flat form this becomes `i = 2i + 1 + maybe_swap_bit`.

**Success criteria.**
- Correctness: identical contact + distance to baseline across all `mesh_mesh`, `mesh_distance`, `distance_lower_bound` test variants.
- Performance: ≥ 30% reduction in `coal-test-benchmark` total time (we already have median 45,956 µs as the post-fix-up baseline).

**Files.**
- New: `include/coal/internal/traversal_node_bvhs_flat.h`, `src/BVH/BVH_flatten.cpp`.
- Modify: `include/coal/BVH/BVH_model.h` (`flatten()` method), `src/traversal/traversal_recurse.cpp` (alt traversal), `src/collision_func_matrix.cpp` (opt-in dispatch via a `BVHModel` flag).

**Workload.** Existing `coal-test-benchmark` (env/rob meshes, all 4 BV types × 3 split methods). Should produce drop-in replacement numbers.

**Risks.** Memory doubles (we keep both tree representations during transition). The "best-first" heuristic may interact badly with strict Eytzinger order — may need to retain a per-node child order bit.

---

### E3 — SIMD-batched OBB SAT (15 axes × 8 BV pairs)

**Hypothesis.** OBB SAT is 15 dot-product-and-compare tests. With AVX2, we can batch 8 BV pairs and run all 15 SAT tests in parallel across them, giving ~5–6× throughput per single-pair-equivalent (after amortizing the gather/scatter).

**Method.**
1. Stage a SoA staging area inside `MeshCollisionTraversalNode` that buffers up to 8 candidate (b1, b2) BV pairs popped from the BVH stack.
2. When the buffer is full (or at end of traversal), run a vectorized 15-axis SAT over the 8 pairs, producing an 8-bit collision mask.
3. For pairs that the SAT marks as overlapping, push their children onto the stack (branchy, scalar — only 0–8 of them).
4. Tail: empty the buffer when no full batch can be formed.

**Success criteria.**
- Correctness: match scalar SAT bit-for-bit on identical fixtures.
- Performance: ≥ 50% reduction in BV-overlap time on collision-heavy scenes (where most BVH work is BV tests, not leaf tests).

**Files.**
- New: `src/BV/OBB_simd.cpp` (or extend `simd_support.cpp`).
- Modify: `include/coal/internal/traversal_node_bvhs.h` (`MeshCollisionTraversalNode::BVTesting` opt-in batched path), `src/BV/OBB.cpp` (declarations).

**Workload.** Same as E2.

**Risks.** Float-precision drift on the absolute-rotation `Bf` matrix. The 15-axis SAT has a published numerical erratum (Gottschalk's 1996 paper has extra epsilon terms); coal's existing scalar code already has these. Re-injecting them into SIMD is doable but error-prone.

---

### E4 — LTO + selective fast-math audit

**Hypothesis.** LTO alone is worth 5–10% on large libraries. `-ffp-contract=fast` (FMA fusion) is safe; `-fno-math-errno`, `-fno-signed-zeros`, `-fno-trapping-math` are very safe; `-fassociative-math` is risky for GJK and SAT termination but cheap to test. Total: maybe 1.05–1.20× across the suite.

**Method.**
1. Add CMake options `COAL_ENABLE_LTO` (default OFF), `COAL_ENABLE_FAST_FP` (default OFF). LTO maps to `-flto=auto`/`-flto`.
2. Decompose the fast-math bag into per-flag CMake options: `COAL_FAST_FP_CONTRACT`, `COAL_FAST_FP_NO_ERRNO`, etc.
3. Build matrix: baseline; LTO only; LTO + safe fast-math; LTO + safe + associative.
4. Run all of: `coal-test-benchmark`, `coal-distance`, `coal-collision`, `coal-distance_lower_bound`, `coal-gjk`, `coal-accelerated_gjk`, `coal-normal_and_nearest_points`, `coal-convex`, `coal-serialization`. Record per-build pass/fail and timing.

**Success criteria.**
- Identify the *largest subset* of fast-math flags that all tests still pass under, and report the speedup.
- LTO must not regress build time by more than 2× (otherwise off by default).

**Files.**
- Modify: top-level `CMakeLists.txt` (options); `src/CMakeLists.txt` (option-conditional `target_compile_options`).
- New (artifact): `bench/fastmath-matrix.md` recording the matrix results.

**Workload.** Whole test suite + benchmark suite, run 3× each per build configuration.

**Risks.** A flag combination might pass tests but introduce silent precision drift in GJK that only manifests at the user-facing distance/contact level. We mitigate by also comparing distance outputs against a known-good baseline within `1e-6` tolerance.

#### E4 Result — 2026-04-29 — **NOT ADOPTED**

Matrix tested on 13th-gen Intel i7-13700H, taskset to one P-core (CPU 4, max 5 GHz), `coal-test-benchmark` 4 interleaved rounds × 6 configs. All flag stacks are *additive on top of* the baseline `-O3 -DNDEBUG -march=native -mavx2 -fno-semantic-interposition`:

| cfg | flags added                         | min (µs) | mean (µs) | vs A min |
|-----|-------------------------------------|---------:|----------:|---------:|
| A   | (baseline)                          | 47486    | 50759     | +0.0%    |
| B   | +`-flto=auto`                       | 53257    | 54532     | **+12.2%** |
| C   | B + `-fno-math-errno -fno-signed-zeros -fno-trapping-math` | 49034 | 52802 | +3.3%   |
| D   | C + `-ffp-contract=fast` (FMA)      | 51255    | 53298     | +7.9%    |
| E   | D + `-fno-rounding-math`            | 51514    | 52990     | +8.5%    |
| F   | E + `-fassociative-math`            | 48726    | 51854     | +2.6%    |

Per-config stdev across rounds is ~1000-3000 µs (noise floor ~6%). Every fast-math configuration is at or below baseline. **LTO actively hurts (+12% mean)** despite producing a smaller `.so` (10.2 MB vs 12.5 MB) — likely because the existing `-fno-semantic-interposition` already enables direct intra-library calls and LTO's additional inlining bloats hot-path code, hurting I-cache. Fast-math additions on top of LTO recover some loss but never beat baseline.

**Decision (per measurement-required policy): no flag adoption.** Skipped Clang-only flags (`-fno-honor-infinities/-nans -fapprox-func`) since the rest of the stack is GCC and they error at configure time.

#### Phase-count triage — 2026-04-29

A new tool, `bin/coal-test-benchmark-phases` (source `test/benchmark_phases.cpp`), enables `enable_statistics` and reports per-query BV-test and leaf-test counts on the same env/rob workload (10000 random transforms, taskset CPU 4):

```
OBB    collide   total 106725 us   bv/q 355289   leaf/q 51270   hits 0/10000
RSS    collide   total 172716 us   bv/q 353751   leaf/q 54093   hits 0/10000
OBBRSS collide   total 137369 us   bv/q 355289   leaf/q 51270   hits 0/10000
kIOS   collide   total 156549 us   bv/q 341928   leaf/q 46609   hits 0/10000
RSS    distance  total   4282 us   bv/q  10039   leaf/q     3
OBBRSS distance  total   4173 us   bv/q  10039   leaf/q     3
kIOS   distance  total   1279 us   bv/q  10189   leaf/q    67
```

Implications:

1. **BV-overlap is the hot path.** For collision: ~7 BV tests per leaf test (350k vs 50k). For distance: 150-3000× (10k vs 3-67). The per-BV cost is ~30 cycles (~7.5 ns at 4 GHz).
2. **RSS / kIOS / OBBRSS are slower per BV test than OBB**: RSS collide takes 1.6× OBB collide for the same BV count. `rectDistance` (16-branch Voronoi) is the inner culprit.
3. **Total time is BV-test bound by a comfortable margin.** Speeding up the BV overlap function 2× would translate near-linearly to total speedup for collision; for distance the gain is even larger because leaf tests are negligible.
4. **Triangle-triangle distance via GJK is NOT the bottleneck.** Even at ~50k leaf tests per collision query, total leaf time is a small fraction of total time (~14% for collision, ~0.03% for distance). E10 (vectorized tri-tri) is therefore deprioritized.

**Implication for E2/E3:** SIMD batching of the BV overlap function (E3) and/or branch-free flat-BVH traversal (E2) might be the well-aimed experiments — but see the *Correction* section below; the original BV-count was wrong by 100,000×.

#### Phase-count triage — Correction & expanded workloads

The above counts had a bug (cumulative counter read as per-query). After fixing and adding rob-vs-rob workloads where transforms are smaller than the mesh extent (1000×1210×850), the corrected picture is:

```
=== env-vs-rob workload: wide [-3000,3000] ===  (broad-phase dominated)
OBB    collide  total 1632 us  bv/q   2.9  leaf/q  0.0  hits 0/10000
RSS    collide  total 2703 us  bv/q   2.9  leaf/q  0.0  hits 0/10000

=== rob-vs-rob workload: overlap [-100,100] ===  (BVH descends to ~depth 6)
OBB    collide  total 2971 us  bv/q   5.9  leaf/q  0.0  hits 0/10000
RSS    collide  total 7737 us  bv/q   6.1  leaf/q  0.0  hits 0/10000
OBBRSS collide  total 2976 us  bv/q   5.9  leaf/q  0.0  hits 0/10000
kIOS   collide  total 4061 us  bv/q   5.9  leaf/q  0.0  hits 0/10000
RSS    distance total 2264 us  bv/q   2.0  leaf/q  0.0
```

Per-BV-test cost (rob-vs-rob, the harder workload):

| BV     | total time | bv/q | per-BV ns |
|--------|-----------:|-----:|----------:|
| OBB    | 2971 µs    | 5.9  | **50 ns** |
| OBBRSS | 2976 µs    | 5.9  | 50 ns     |
| kIOS   | 4061 µs    | 5.9  | 69 ns     |
| RSS    | 7737 µs    | 6.1  | **127 ns** |

**RSS `rectDistance` is the slow-link: 2.5× the per-call cost of OBB SAT.** It is the natural per-call optimization target — the function has 16 Voronoi-region branches, and even a modest scalar tightening would compound across millions of BVH-edge tests.

Note also: even with overlap workloads, **0 leaf-tri-tri tests fire** because the rob mesh has only 216 triangles arranged sparsely — BV-level rejection prunes everything before reaching leaves. This benchmark therefore measures *only the BV-overlap path*, not the leaf path. To measure leaf paths we would need denser, actually-overlapping geometry. For the current optimization round, the BV-overlap path is the right focus given its dominance.

#### BV micro-benchmark — 2026-04-29

A second benchmark, `bin/coal-test-benchmark-bv-micro` (source `test/benchmark_bv_micro.cpp`), times a tight loop over `RSS::overlap` and `OBB::overlap` on 1024 randomly-rotated input pairs across 4 translation extents. This isolates the BV-overlap function from BVH traversal noise (no front-list, no stack pushes).

```
# RSS::overlap (one call per iteration)
trans_extent       calls/sec        ns/call
0.5                 10570253          94.61    # heavy overlap regime
2                   15197267          65.80
5                   17168922          58.24
20                  17766728          56.28    # disjoint regime, early-exit

# OBB::overlap (one call per iteration)
trans_extent       calls/sec        ns/call
0.5                 14362135          69.63    # heavy overlap regime
2                   14897542          67.13
5                   17481596          57.20
20                  18712877          53.44    # disjoint regime, early-exit
```

In the disjoint regime (the dominant case in BVH traversal once descent terminates), OBB and RSS are within 2 ns of each other (~55 ns/call). In the overlap regime, RSS pays a 25 ns penalty over OBB because `rectDistance` walks more Voronoi-region branches before finding the closest pair.

#### Pathological dense-overlap workload — 2026-04-29

To stress the BVH descent path I added a synthetic-soup workload (2000 random triangles per side, 0.4-unit extent each, in a [-1,1] cube) plus a *self-vs-self near-identity* workload (same mesh on both sides, transforms in [-1mm, 1mm]):

```
=== synthetic dense soup [-0.1,0.1] (1000 queries) ===
OBB    collide  total  778 us  bv/q  14.2  leaf/q 0.1  hits 0/1000
RSS    collide  total 1954 us  bv/q  14.2  leaf/q 0.1  hits 0/1000

=== pathological self-vs-self near-identity (100 queries) ===
OBB    collide  total  122 us  bv/q  16.4  leaf/q 0.8  hits 0/100
RSS    collide  total  276 us  bv/q  16.7  leaf/q 0.6  hits 0/100
```

Even *identical* meshes at near-identity transforms only produce 16 BV tests/query and 1 leaf test/query. coal's BVH is extremely effective at pruning. The conclusion: **per-query, a coal mesh-mesh collision is dominated by 50-200 ns of fixed setup + ~10 BV tests at ~50-150 ns each = 500-2000 ns/query.** That's exactly what the wall-clock data shows (1-2 µs/query).

Implications for an optimization budget:

- A **2× per-call BV-test speedup** translates to roughly 1.4-1.7× total query speedup — not the 5-6× we'd hoped for from SIMD batching, because BV tests are not the only cost.
- **SIMD batching of 8 BV pairs** still requires queue-based traversal (E2) to populate the batch. With only ~10 BV tests/query, batching is hard to feed — average batch fill is much less than 8.
- **Virtual function dispatch in `collisionRecurse`** (`isFirstNodeLeaf`, `isSecondNodeLeaf`, `BVDisjoints` are all virtual on `CollisionTraversalNodeBase`) likely adds 10-15 ns per BV test = 20% overhead. CRTP/template devirtualization is a real candidate but is a substantial refactor across the traversal-node class hierarchy. *This is the most promising single-target optimization given current evidence.*
- **`std::vector<BVPair>::reserve(1000)` in `collisionNonRecurse`** is an 8KB allocation per query — but the default uses recursive traversal so this isn't hit.

#### Decision summary as of 2026-04-29

Tested and **not adopted**:
- E4 LTO + fast-math: LTO regresses 12%, fast-math additions don't recover; full table in the E4 section.

Tested but **inconclusive due to workload**:
- The existing benchmark + the new dense/pathological synthetic workloads all keep BV-test count per query under ~16. There is no benchmark in this round where BV tests dominate enough to make per-call BV-test optimization (E3, RSS rectDistance scalar tightening) yield clearly measurable end-to-end gains beyond ~5-10%.

**Best path forward, ranked:**
1. **CRTP/template traversal node (new candidate, was implicit in E2):** eliminate ~10-15 ns of virtual call overhead per BV test. Estimated 15-20% on collision queries with high BV count. Requires non-trivial refactor of `MeshCollisionTraversalNode`/`MeshDistanceTraversalNode` — but still within "no public API changes" since these are internal traversal nodes. *Recommend as next experiment.*
2. **Eytzinger flat BVH (E2):** complementary to (1). Bigger refactor; would also benefit cache behavior.
3. **rectDistance scalar tightening:** small win (5-10% on RSS-dominated queries). Reasonable filler experiment.

Tooling added (committed):
- `test/benchmark_phases.cpp` → `coal-test-benchmark-phases` — operation-count and per-workload timings, with `enable_statistics`.
- `test/benchmark_bv_micro.cpp` → `coal-test-benchmark-bv-micro` — per-call BV-overlap micro-benchmark with controlled input distributions.

---

### E5 — Vectorized convex hill-climb (`getShapeSupportLog`)

**Hypothesis.** The neighbor walk in `getShapeSupportLog` (line 356-376 of `support_functions.cpp`) currently evaluates one neighbor's dot product at a time. Convex polytopes typically have 4-12 neighbors per vertex; we can batch them with AVX2 (4-wide for double, 8-wide for float) and still respect the visited bitmap.

**Method.**
1. New `simd::maxDotIndexed(points, indices, dir, &maxdot)` — gathers via `_mm256_i32gather_ps` (or scalar gather followed by SIMD reduction; the gather instruction is slow on Zen and can lose to scalar) and reduces.
2. In the hill-climb loop, replace the `for (in = 0; in < n.count; ++in)` scalar walk with one SIMD batch + a tail scalar handler when `n.count` is not a multiple of the SIMD width. Keep the visited bitmap update as scalar-after-batch.

**Success criteria.**
- Correctness: identical hint and support point as scalar across a fuzz test of 1000 random convex hulls × 100 directions.
- Performance: ≥ 20% reduction in `coal-gjk` and `coal-accelerated_gjk` on large convex (256+ vertex) inputs.

**Files.**
- Modify: `include/coal/internal/simd_support.h`, `src/narrowphase/simd_support.cpp`, `src/narrowphase/support_functions.cpp` (`getShapeSupportLog`).
- Add: `test/simd_indexed_dot.cpp`.

**Workload.** New benchmark: `coal-test-benchmark-convex` synthesizing convex hulls of 32, 128, 512, 2048 vertices and running 1M GJK queries.

**Risks.** Gather instructions are bandwidth-limited; benefit may be small for small meshes (≤ 32 verts) — keep the existing scalar path as fallback.

---

### E6 — kDOP SIMD dot products

**Hypothesis.** kDOP-N for N ∈ {16, 18, 24} computes N axis-aligned-and-diagonal dot products per BV. With AVX2 these can be done 4-8 wide trivially. Probably 1.5-2× kDOP overlap throughput.

**Method.**
1. In `src/BV/kDOP.cpp`, replace `getDistances<N>()` template specializations with SIMD reductions over the 4-or-8-wide vectorized dot products.
2. Same for the kDOP overlap test.

**Success criteria.** Match scalar bit-for-bit on test fixtures; ≥ 30% reduction in kDOP-using BVH benchmarks.

**Files.** `src/BV/kDOP.cpp`, possibly `include/coal/BV/kDOP.h`.

**Workload.** Add a kDOP variant to `coal-test-benchmark`.

**Risks.** Lowest-risk change in this list — kDOP test data is intrinsically AoS-aligned to dot-products.

---

### E7 — Sphere-tree midphase as opt-in pre-filter

**Hypothesis.** For BVH-vs-BVH queries where one geometry is amenable to FOAM-style sphere approximation (think: robot link meshes), generating a hierarchical sphere tree at construction time and using it as a midphase pre-filter would cut deep BVH descents in obviously-non-overlapping cases. Borrowed directly from VAMP § III-C.

**Method.**
1. Add a `BVHModel<BV>::buildSphereTree(SphereTreeBuilder)` method that, optionally, attaches a multi-level sphere cover to the model.
2. In `MeshCollisionTraversalNode::canStop` / `BVTesting`, check the sphere tree first when both models have one; return non-overlap early if the bounding spheres don't overlap.
3. Provide a Python utility that wraps FOAM (`pip install foam-spherize`) to generate the sphere tree from a `BVHModel`'s points.

**Success criteria.**
- Correctness: cannot return a false negative — sphere tree must conservatively bound the mesh. Tested via random sampling of mesh-vs-mesh queries.
- Performance: 2-5× speedup on the subset of MotionBenchMaker scenes where geometry is sphere-tractable (sphere-cover error < 5% volume).

**Files.**
- New: `include/coal/BVH/sphere_tree.h`, `src/BVH/sphere_tree.cpp`, `python/coal/foam_integration.py`.
- Modify: `BVHModel` to optionally hold a sphere tree.

**Workload.** Existing benchmark; add a "with sphere tree" / "without" toggle and report.

**Risks.** Building the sphere tree is offline work; if used wrong (per-query rebuild) it's a regression.

---

### E8 — Sphere-vs-AABB SIMD prefilter

**Hypothesis.** For mesh queries against a robot represented as spheres, a fast sphere-vs-mesh-AABB pre-pass kills most of the mesh's BVH before traversal even starts. AVX2 batches 8 sphere-AABB tests per instruction.

**Method.** New free function `simd::sphereAabbBatch8(centers_soa, radius, aabbs_soa, &mask_out)` returning an 8-bit overlap mask. Wire into the front-end of mesh-vs-pointcloud-of-spheres queries.

**Files.** `include/coal/internal/simd_support.h`, `src/narrowphase/simd_support.cpp`. New test fixture.

**Risks.** Low. Speedup magnitude depends on sphere-cluster geometry.

---

### E9 — AABB SIMD distance

**Hypothesis.** `BV/AABB.cpp::distance()` uses a scalar 3-axis `for` loop. Replacing with `(min - p).cwiseMax(0).cwiseMax(p - max)` followed by `.norm()` already gives Eigen-AVX2 vectorization, but a hand-coded SSE sequence may shave 5-10%.

**Method.** Replace the loop with explicit SSE/AVX2 minmax-clamp-norm; benchmark vs Eigen-vectorized.

**Files.** `src/BV/AABB.cpp`.

**Risks.** Lowest in the list. Probably falls within Eigen's auto-vectorization window — pick this only if the auto-vectorized version doesn't match hand SIMD.

---

### E10 — Vectorized triangle-triangle distance

**Hypothesis.** `triangle_triangle.cpp` currently delegates to GJK. A direct closed-form algorithm (e.g., the v-clip-style approach in Ericson's *Real-Time Collision Detection* § 5.2.2 or the original Larsen/Gottschalk OBB-tree paper) can be 2-3× faster, and lends itself to AVX2 batching of 8 triangle pairs.

**Method.** Implement Larsen's segment-segment + point-triangle decomposition (15 sub-cases), in scalar first, then in 8-pair AVX2 batch form. Replace GJK delegation when both shapes are `TriangleP`.

**Files.** `src/distance/triangle_triangle.cpp`.

**Risks.** Highest correctness risk in this list — triangle-triangle distance has many degenerate cases (parallel edges, collinear points). Must have a thorough fuzz test before declaring success.

## 5. Order of execution

1. **Sprint 0 (3-4 days).** Run E4 (LTO + fast-math audit). Establishes new baseline. Adopt the safe subset.
2. **Sprint 1 (1 week).** E6 (kDOP SIMD) + E9 (AABB distance). Both small. Validates SIMD test infrastructure.
3. **Sprint 2 (2 weeks).** E5 (convex hill-climb). Builds on the simd_support module we already have.
4. **Sprint 3 (2-3 weeks).** E1 (CAPT). Highest absolute value if pointcloud workloads matter. Independent of other experiments.
5. **Sprint 4 (2 weeks).** E2 (Eytzinger BVH). High value for general mesh-mesh workloads. Touches hot core — needs full test suite + benchmark validation.
6. **Sprint 5 (2 weeks).** E3 (SIMD SAT). Compounds with E2.
7. **Sprint 6+ (as time permits).** E7 (sphere-tree midphase), E8 (sphere-AABB prefilter), E10 (vectorized tri-tri).

## 6. Verification plan (per experiment)

Every experiment must satisfy these gates before merge:

1. **Correctness.** All existing tests under `test/` pass with `-DNDEBUG` and `-DCMAKE_BUILD_TYPE=Debug`. New tests cover the new code paths with both float and double `Scalar`.
2. **Numerical regression.** Output of `coal-distance`, `coal-normal_and_nearest_points`, `coal-distance_lower_bound` agrees with the pre-experiment baseline within `1e-6` absolute (`1e-3` for float `Scalar`).
3. **Performance regression.** `coal-test-benchmark` total time does not regress on any individual cell vs the pre-experiment baseline by more than 5% (some cells may speed up, some may not — but no cell should silently slow down).
4. **Build.** `clean && cmake -DCMAKE_BUILD_TYPE=Release && make -j` from a fresh build dir; `make test` green.
5. **Docs.** Each experiment lands a short note in `REASEARCH.md` in a `## Results` section appended at the bottom: actual speedup vs predicted, surprises, follow-ups.

## 7. Outstanding unknowns / things worth re-checking the literature for

- **AVX-512 on newer CPUs.** VAMP found AVX-512 slower than AVX2 on Ryzen 9 7950X (Zen 4). Intel Sapphire Rapids and beyond may invert this. Do not bake AVX2 in as the only target; keep dispatch flexible.
- **Whether VAMP-MR's 86% cache-miss reduction figure replicates on coal-style mesh workloads.** Their workload is sphere-only and SoA-flat; coal's BVH is fundamentally pointier.
- **CAPT's worst-case memory on coal's typical workload.** CAPT paper warns about low-dispersion clouds. Verify on real depth-camera output (Intel RealSense D455 is what they used).
- **Whether `-fassociative-math` survives GJK termination tests.** This is the single biggest fast-math win candidate — but it can reorder reductions that GJK's epsilon-comparisons depend on.

## 8. References (for re-grounding without re-reading)

Local copies of the extracted PDF text are at:
- `/home/davide/.claude/projects/-home-davide-Asensus-meshopt-coal/60fca005-84ea-4344-a650-62c5b1e56b78/tool-results/webfetch-1777470034984-2pxe7v.txt` — VAMP main paper (Thomason 2024)
- `/home/davide/.claude/projects/-home-davide-Asensus-meshopt-coal/60fca005-84ea-4344-a650-62c5b1e56b78/tool-results/webfetch-1777470042660-l8yi10.txt` — CAPT (Ramsey 2024)
- `/home/davide/.claude/projects/-home-davide-Asensus-meshopt-coal/60fca005-84ea-4344-a650-62c5b1e56b78/tool-results/webfetch-1777470044621-df1blj.txt` — VAMP-MR (Huang 2026)

Online:
- VAMP code: https://github.com/KavrakiLab/vamp (also locally at `/home/davide/Asensus/meshopt/vamp/`)
- VAMP main paper: https://arxiv.org/pdf/2309.14545
- CAPT paper: https://www.roboticsproceedings.org/rss20/p038.pdf
- VAMP-MR paper: https://openreview.net/pdf?id=ePPOoz8KKp
- FOAM paper: https://arxiv.org/pdf/2503.13704
- FOAM repo: https://github.com/CoMMALab/foam

Key VAMP source landmarks (under `/home/davide/Asensus/meshopt/vamp/src/impl/vamp/`):
- `vector/isa/avx.hh` — AVX2 SIMD primitives, `VectorWidth = 8`, `Alignment = 32`.
- `vector/isa/neon.hh` — ARM NEON, 4-wide.
- `vector/isa/wasm.hh` — WebAssembly SIMD, 4-wide.
- `robots/baxter.hh:233+` — Baxter FK (75 spheres, 14 DOF, ~14k generated ops).
- `robots/sphere.hh:103-108` — `fkcc(env, ConfigurationBlock<rake>)` entry point.
- `collision/sphere_sphere.hh:10-40` — sphere-sphere SQL2.
- `collision/sphere_capsule.hh:9-28` — sphere-capsule with cached `rdv`.
- `collision/sphere_cuboid.hh:9-26` — sphere-cuboid via abs+clamp+sum-of-squares.
- `collision/capt.hh:61-400+` — CAPT data structure (the actual implementation referenced by Ramsey 2024).
- `collision/environment.hh:14-69` — typed obstacle collections, sorted by `min_distance`.
- `collision/math.hh:17-42` — `dot_3`, `sql2_3`.
- `cmake/CompilerSettings.cmake:22-58` — fast-math + LTO setup.
