// Correctness regression test.
//
// Runs the polso/gen4 mesh-vs-mesh manifests through coal::collide and checks
// the per-pose collision boolean against committed baselines captured from
// the `devel` branch (see test/regression/README.md for regen instructions).
// Any flip — pose that was collision on devel but not on the current branch,
// or vice versa — is a correctness regression.
//
// The test skips silently if the external dataset is not present (the
// manifests reference STL files outside the coal repo, so this test is
// only meaningful on machines that have the meshopt dataset checked out
// alongside coal). Set COAL_REGRESSION_DATASET_DIR to override the default
// search location.

#define BOOST_TEST_MODULE COAL_REGRESSION_CORRECTNESS
#include <boost/test/included/unit_test.hpp>

#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <chrono>
#include <cstdlib>
#include <fstream>
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
  BOOST_REQUIRE_MESSAGE(f.is_open(), "cannot open poses CSV: " + path);
  std::vector<Pose> out;
  std::string line;
  std::getline(f, line);  // header
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

std::vector<int> readBaselinePerPoseCsv(const std::string& path) {
  std::ifstream f(path);
  BOOST_REQUIRE_MESSAGE(f.is_open(), "cannot open baseline CSV: " + path);
  std::vector<int> out;
  std::string line;
  std::getline(f, line);  // header
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    // Format: pose_idx,original[,...other_reps]
    auto comma1 = line.find(',');
    if (comma1 == std::string::npos) continue;
    auto comma2 = line.find(',', comma1 + 1);
    const std::string col1 =
        line.substr(comma1 + 1, comma2 == std::string::npos
                                    ? std::string::npos
                                    : comma2 - comma1 - 1);
    out.push_back(std::stoi(col1));
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

// Per-pose collision booleans for the "original" mesh-vs-mesh representation
// in a manifest. Mirrors the `original` column emitted by scripts/coal_bench.
std::vector<int> runManifestPerPose(const std::string& manifest_path) {
  pt::ptree m;
  pt::read_json(manifest_path, m);

  // Locate the "original" representation (mesh-vs-mesh path). We don't need
  // to run the other reps for correctness regression — they're for the
  // user's mesh-approximation study, not coal correctness.
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
  BOOST_REQUIRE_MESSAGE(!path_a.empty(),
                        "no 'original' mesh representation in " +
                            manifest_path);

  auto poses_path = m.get<std::string>("poses_csv");
  auto poses = readPosesCsv(poses_path);

  auto geom_a = loadMeshBvh(path_a);
  auto geom_b = loadMeshBvh(path_b);
  // Allocate the two collision objects once and reuse across poses; the
  // scripts/coal_bench.cpp bench does the same to keep timings clean and
  // we need to mirror that here so timings (and per-pose state) match.
  auto obj_a = std::make_shared<CollisionObject>(geom_a, Transform3s::Identity());
  auto obj_b = std::make_shared<CollisionObject>(geom_b, Transform3s::Identity());

  DynamicAABBTreeCollisionManager mgr_a, mgr_b;
  mgr_a.registerObject(obj_a.get());
  mgr_a.setup();

  std::vector<int> per_pose;
  per_pose.reserve(poses.size());
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
    per_pose.push_back(cb.any_collision ? 1 : 0);
  }
  return per_pose;
}

// Where the meshopt dataset (manifests + STLs) lives. The manifests embed
// absolute paths to STL files, so this test is only meaningful when the
// matching dataset is on disk. Default search: parent-of-coal-repo + /output.
fs::path findDatasetDir() {
  if (const char* env = std::getenv("COAL_REGRESSION_DATASET_DIR")) {
    return fs::path(env);
  }
  // Default: <repo>/.. /output (the meshopt workspace layout used here).
  fs::path repo = fs::path(__FILE__).parent_path().parent_path();
  fs::path candidate = repo.parent_path() / "output";
  return candidate;
}

}  // namespace

BOOST_AUTO_TEST_CASE(polso_gen4_per_pose_matches_baseline) {
  const fs::path dataset = findDatasetDir();
  if (!fs::exists(dataset)) {
    BOOST_TEST_MESSAGE(
        "Skipping: dataset directory not found at " + dataset.string() +
        ". Set COAL_REGRESSION_DATASET_DIR to point to the meshopt "
        "output/ directory if you have it locally.");
    return;
  }

  const fs::path baseline_dir =
      fs::path(__FILE__).parent_path() / "regression";

  const std::vector<std::string> joints = {"J4-J5", "J5-J6", "J6-J7", "J7",
                                           "Nose"};
  for (const auto& joint : joints) {
    const fs::path manifest =
        dataset / ("bench_originals_manifest_polso_gen4_" + joint + ".json");
    const fs::path baseline =
        baseline_dir /
        ("baseline_polso_gen4_" + joint + "_perpose.csv");
    if (!fs::exists(manifest)) {
      BOOST_TEST_MESSAGE("Skipping " + joint + ": manifest not found at " +
                         manifest.string());
      continue;
    }
    if (!fs::exists(baseline)) {
      BOOST_TEST_MESSAGE("Skipping " + joint + ": baseline not found at " +
                         baseline.string());
      continue;
    }

    BOOST_TEST_MESSAGE("Running joint " + joint);
    const auto current = runManifestPerPose(manifest.string());
    const auto expected = readBaselinePerPoseCsv(baseline.string());
    BOOST_REQUIRE_EQUAL(current.size(), expected.size());

    std::size_t disagree = 0;
    std::size_t first_disagree_idx = 0;
    for (std::size_t i = 0; i < current.size(); ++i) {
      if (current[i] != expected[i]) {
        if (disagree == 0) first_disagree_idx = i;
        ++disagree;
      }
    }
    BOOST_TEST_MESSAGE("  " + joint + ": " + std::to_string(disagree) +
                       " / " + std::to_string(current.size()) +
                       " poses disagree with baseline");
    if (disagree > 0) {
      BOOST_TEST_MESSAGE("  first disagreement at pose " +
                         std::to_string(first_disagree_idx) +
                         " (baseline=" +
                         std::to_string(expected[first_disagree_idx]) +
                         ", current=" +
                         std::to_string(current[first_disagree_idx]) + ")");
    }
    BOOST_REQUIRE_EQUAL(disagree, 0u);
  }
}
