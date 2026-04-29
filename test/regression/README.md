# Regression test infrastructure

This directory holds the committed reference outputs for two regression tests:

- **`baseline_polso_gen4_*_perpose.csv`** — per-pose collision booleans for
  each of the 5 polso/gen4 joints, captured from the `devel` branch.
  Consumed by `test/regression_correctness.cpp`. Any flip on a future
  branch is a correctness regression.
- **`baseline_polso_gen4_timings.json`** — median per-pose query time
  (across 11 runs, P-core pinned) for each joint, captured from
  `feature/coal-simd-support` HEAD on a 13th-gen Intel laptop. Consumed by
  `test/regression_perf.cpp`. A future branch slower by more than the
  per-joint `tolerance_pct` is a performance regression.

## Running the tests

The regression tests need the meshopt dataset on disk (the manifests
embed absolute STL paths). They skip silently when the dataset isn't
found.

```bash
# Default ctest run includes regression_correctness; it auto-discovers
# the dataset in the meshopt workspace layout.
ctest --output-on-failure

# To run the perf gate too (opt-in, machine-class-gated):
COAL_MACHINE_CLASS=13th-gen-Intel-i7-13700H-P-core-pinned \
  taskset -c 0 ctest --label-regex perf --output-on-failure

# Or force-run the perf test on a different machine to see actual numbers
# (won't fail, just prints):
COAL_PERF_FORCE=1 taskset -c 0 ./bin/coal-regression_perf \
  --log_level=message

# To run from a different dataset location:
COAL_REGRESSION_DATASET_DIR=/path/to/output \
  ./bin/coal-regression_correctness --log_level=message
```

## When to regenerate

### Correctness baselines (`*_perpose.csv`)

Regenerate **only** when an intentional change in collision semantics is
made (a known-bug fix that flips some collision answers). The regen MUST
be a separate commit with a clear message explaining what flipped and
why. Reviewer should check the diff of the CSV files manually.

```bash
# 1. Switch to devel (or whichever branch is the "ground truth")
git checkout devel
# 2. Build coal and the bench
cmake --build /path/to/coal/build --target coal -j12
cmake --build /path/to/scripts/build --target coal_bench -j12
# 3. Re-emit each per-pose CSV
for joint in J4-J5 J5-J6 J6-J7 J7 Nose; do
  manifest=/path/to/output/bench_originals_manifest_polso_gen4_${joint}.json
  out_csv=/path/to/coal/test/regression/baseline_polso_gen4_${joint}_perpose.csv
  taskset -c 0 /path/to/scripts/build/coal_bench \
    "$manifest" /tmp/null.csv "$out_csv"
done
# 4. Switch back to your feature branch and commit the regen
git checkout feature/<your-branch>
git add test/regression/baseline_polso_gen4_*_perpose.csv
git commit -m "test: regenerate correctness baselines after <reason>"
```

### Timing baseline (`baseline_polso_gen4_timings.json`)

Regenerate **only** when:

- Hardware changes (different CPU model, different thermal envelope) and
  the perf gate is now permanently failing — replace `machine_class` and
  re-measure.
- A perf-positive change ships and you want the new floor to protect
  future regressions against the new floor (not the old one).

```bash
# 1. Make sure you're on the branch whose perf you want to lock in.
git status
# 2. Build coal and rebuild the bench against it.
cmake --build /path/to/coal/build --target coal -j12
cmake --build /path/to/scripts/build --target coal_bench -j12
# 3. Capture 11 runs per joint.
for joint in J4-J5 J5-J6 J6-J7 J7 Nose; do
  manifest=/path/to/output/bench_originals_manifest_polso_gen4_${joint}.json
  vals=""
  for i in $(seq 1 11); do
    v=$(taskset -c 0 /path/to/scripts/build/coal_bench \
        "$manifest" /tmp/null.csv 2>&1 | grep -oP 'mean=\K[0-9.]+')
    vals="$vals $v"
  done
  median=$(echo $vals | tr ' ' '\n' | sort -n | awk 'NR==6')
  printf "%-8s : median=%s | %s\n" "$joint" "$median" "$vals"
done
# 4. Hand-edit baseline_polso_gen4_timings.json with the new medians.
# 5. Commit.
```

## Tolerance choice

The committed `tolerance_pct` is 5%. With 11-run median and CPU-pinned
measurement, run-to-run noise on the dev hardware is ~1-2%, so a 5% gate
gives a comfortable margin against transient thermal drift while still
catching real regressions ≥ 3-4%.

If you tighten the tolerance, also tighten the run count or pin the
governor more aggressively — the false-positive rate scales with
(noise_floor / tolerance) ratio.
