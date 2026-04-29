// Performance regression test (opt-in, label = "perf").
//
// Runs the polso/gen4 mesh-vs-mesh manifests through coal::collide N times,
// computes the median per-pose query time, and asserts it is within a
// committed tolerance of the median captured from the feature branch when
// the timing JSON was last regenerated. Skips silently when:
//   - the dataset is not present (manifests reference STL files outside coal),
//   - or the machine_class field doesn't match COAL_MACHINE_CLASS env var
//     (override with COAL_PERF_FORCE=1 to run anyway and just print results).
//
// Not in the default ctest run. Invoke with:
//   ctest --label-regex perf
// or directly:
//   taskset -c 0 ./bin/coal-regression_perf

#define BOOST_TEST_MODULE COAL_REGRESSION_PERF
#include <boost/test/included/unit_test.hpp>

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "coal/collision.h"
#include "coal/collision_object.h"
#include "coal/shape/geometric_shapes.h"
#include "coal/BVH/BVH_model.h"
#include "coal/mesh_loader/assimp.h"
#include "coal/broadphase/broadphase_dynamic_AABB_tree.h"
#include "coal/broadphase/default_broadphase_callbacks.h"

namespace pt = boost::property_tree;
namespace fs = boost::filesystem;
using namespace coal;
using clk = std::chrono::high_resolution_clock;

namespace {

struct Pose {
  Vec3s t;
  Matrix3s R;
};

Matrix3s eulerToMatrix(Scalar rx, Scalar ry, Scalar rz) {
  using std::cos;
  using std::sin;
  Matrix3s Rx, Ry, Rz;
  Rx << 1, 0, 0, 0, cos(rx), -sin(rx), 0, sin(rx), cos(rx);
  Ry << cos(ry), 0, sin(ry), 0, 1, 0, -sin(ry), 0, cos(ry);
  Rz << cos(rz), -sin(rz), 0, sin(rz), cos(rz), 0, 0, 0, 1;
  return Rz * Ry * Rx;
}

std::vector<Pose> readPosesCsv(const std::string& path) {
  std::ifstream f(path);
  std::vector<Pose> out;
  std::string line;
  std::getline(f, line);
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    std::stringstream ss(line);
    std::string token;
    std::vector<double> vals;
    while (std::getline(ss, token, ',')) vals.push_back(std::stod(token));
    if (vals.size() < 6) continue;
    Pose p;
    p.t = Vec3s(vals[3], vals[4], vals[5]);
    p.R = eulerToMatrix(vals[0], vals[1], vals[2]);
    out.push_back(p);
  }
  return out;
}

std::shared_ptr<CollisionGeometry> loadMeshBvh(const std::string& path) {
  auto bvh = std::make_shared<BVHModel<OBBRSS>>();
  loadPolyhedronFromResource(path, Vec3s::Ones(), bvh);
  bvh->computeLocalAABB();
  return bvh;
}

struct CollideCallback : CollisionCallBackBase {
  bool any_collision = false;
  CollisionRequest req;
  CollideCallback() : req(CollisionRequestFlag::NO_REQUEST, 1) {}
  bool collide(CollisionObject* o1, CollisionObject* o2) override {
    CollisionResult res;
    coal::collide(o1, o2, req, res);
    if (res.isCollision()) any_collision = true;
    return false;
  }
};

// Returns mean-ns-per-query for one full sweep of a manifest. Setup
// (BVH load, manager construction) is excluded; only the per-pose
// collide loop is timed. Same accounting as scripts/coal_bench.
double measureManifestMeanNs(const std::string& manifest_path) {
  pt::ptree m;
  pt::read_json(manifest_path, m);

  std::string path_a, path_b;
  for (const auto& kv : m.get_child("representations")) {
    const auto& rep = kv.second;
    if (rep.get<std::string>("name") == "original" &&
        rep.get<std::string>("type") == "mesh") {
      path_a = rep.get<std::string>("path_a");
      path_b = rep.get<std::string>("path_b");
      break;
    }
  }
  BOOST_REQUIRE(!path_a.empty());
  auto poses = readPosesCsv(m.get<std::string>("poses_csv"));
  auto geom_a = loadMeshBvh(path_a);
  auto geom_b = loadMeshBvh(path_b);
  // Allocate the two collision objects once; the per-pose hot loop has to
  // mirror scripts/coal_bench.cpp exactly to make the committed timings
  // comparable. Per-pose shared_ptr<CollisionObject> allocation is roughly
  // an order of magnitude per query — easy to forget.
  auto obj_a =
      std::make_shared<CollisionObject>(geom_a, Transform3s::Identity());
  auto obj_b =
      std::make_shared<CollisionObject>(geom_b, Transform3s::Identity());

  DynamicAABBTreeCollisionManager mgr_a, mgr_b;
  mgr_a.registerObject(obj_a.get());
  mgr_a.setup();

  auto t0 = clk::now();
  for (const auto& pose : poses) {
    Transform3s tf;
    tf.setRotation(pose.R);
    tf.setTranslation(pose.t);
    obj_b->setTransform(tf);
    obj_b->computeAABB();
    mgr_b.clear();
    mgr_b.registerObject(obj_b.get());
    mgr_b.setup();
    CollideCallback cb;
    mgr_a.collide(&mgr_b, &cb);
  }
  auto t1 = clk::now();
  return std::chrono::duration<double, std::nano>(t1 - t0).count() /
         static_cast<double>(poses.size());
}

fs::path findDatasetDir() {
  if (const char* env = std::getenv("COAL_REGRESSION_DATASET_DIR")) {
    return fs::path(env);
  }
  fs::path repo = fs::path(__FILE__).parent_path().parent_path();
  return repo.parent_path() / "output";
}

}  // namespace

BOOST_AUTO_TEST_CASE(polso_gen4_perf_within_tolerance) {
  const fs::path baseline_json = fs::path(__FILE__).parent_path() /
                                 "regression" /
                                 "baseline_polso_gen4_timings.json";
  if (!fs::exists(baseline_json)) {
    BOOST_TEST_MESSAGE("Skipping: timings baseline not found at " +
                       baseline_json.string());
    return;
  }

  const fs::path dataset = findDatasetDir();
  if (!fs::exists(dataset)) {
    BOOST_TEST_MESSAGE("Skipping: dataset not found at " + dataset.string());
    return;
  }

  pt::ptree j;
  pt::read_json(baseline_json.string(), j);
  const std::string baseline_machine = j.get<std::string>("machine_class");

  // Machine-class gating. The bench is highly sensitive to CPU model and
  // pinning, so timings captured on a 13th-gen Intel laptop won't hold on
  // (say) a Ryzen workstation. Default behavior: skip silently. Force-run
  // with COAL_PERF_FORCE=1 (useful when establishing a new baseline).
  const char* user_machine = std::getenv("COAL_MACHINE_CLASS");
  const bool force = std::getenv("COAL_PERF_FORCE") != nullptr;
  if (!force && (!user_machine || baseline_machine != user_machine)) {
    BOOST_TEST_MESSAGE(
        "Skipping: COAL_MACHINE_CLASS does not match baseline (" +
        baseline_machine +
        "). Set COAL_MACHINE_CLASS=<that string> on a matching machine, "
        "or set COAL_PERF_FORCE=1 to run anyway.");
    return;
  }

  const int kRuns = 11;  // odd → median is well-defined without averaging
  for (const auto& kv : j.get_child("joints")) {
    const std::string joint = kv.first;
    const double baseline_ns = kv.second.get<double>("median_ns");
    const double tolerance_pct = kv.second.get<double>("tolerance_pct");
    const fs::path manifest =
        dataset /
        ("bench_originals_manifest_polso_gen4_" + joint + ".json");
    if (!fs::exists(manifest)) {
      BOOST_TEST_MESSAGE("Skipping " + joint + ": manifest not found");
      continue;
    }

    std::vector<double> samples;
    samples.reserve(kRuns);
    for (int r = 0; r < kRuns; ++r) {
      samples.push_back(measureManifestMeanNs(manifest.string()));
    }
    std::sort(samples.begin(), samples.end());
    const double median = samples[kRuns / 2];
    const double lo = samples.front();
    const double hi = samples.back();
    const double ratio = median / baseline_ns;

    std::ostringstream report;
    report << std::fixed << std::setprecision(2) << "  " << joint
           << ": median=" << median << " ns (baseline " << baseline_ns
           << ", ratio " << ratio << ", tolerance "
           << (1.0 + tolerance_pct / 100.0) << ", min " << lo << ", max "
           << hi << ")";
    BOOST_TEST_MESSAGE(report.str());

    if (!force) {
      const double allowed = baseline_ns * (1.0 + tolerance_pct / 100.0);
      BOOST_CHECK_LE(median, allowed);
    }
  }
}
