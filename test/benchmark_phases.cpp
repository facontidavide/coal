// Phase-counter benchmark: reports num_bv_tests and num_leaf_tests density
// per query for mesh-vs-mesh collision and distance, across BV types.
// Used to triage which BVH-side optimization is worth implementing.

#include <boost/filesystem.hpp>
#include <random>

#include "coal/internal/traversal_node_setup.h"
#include "coal/internal/traversal_node_bvhs.h"
#include "../src/collision_node.h"
#include "coal/internal/BV_splitter.h"

#include "utility.h"
#include "fcl_resources/config.h"

using namespace coal;

template <typename BV>
void makeModel(const std::vector<Vec3s>& vertices,
               const std::vector<Triangle32>& triangles, BVHModel<BV>& model) {
  model.bv_splitter.reset(new BVSplitter<BV>(SPLIT_METHOD_MEAN));
  model.beginModel();
  model.addSubModel(vertices, triangles);
  model.endModel();
}

template <typename BV, typename TraversalNode>
void runCollide(const std::vector<Transform3s>& tf, const BVHModel<BV>& m1,
                const BVHModel<BV>& m2, const char* name) {
  Transform3s pose2;
  CollisionRequest request;
  TraversalNode node(request);
  node.enable_statistics = true;

  std::size_t total_bv = 0, total_leaf = 0;
  std::size_t hits = 0;
  CollisionResult local_result;
  BenchTimer t;
  t.start();
  for (std::size_t i = 0; i < tf.size(); ++i) {
    bool ok = initialize(node, m1, tf[i], m2, pose2, local_result);
    (void)ok;
    node.num_bv_tests = 0;
    node.num_leaf_tests = 0;
    CollisionResult result;
    collide(&node, request, result);
    total_bv += static_cast<std::size_t>(node.num_bv_tests);
    total_leaf += static_cast<std::size_t>(node.num_leaf_tests);
    if (result.numContacts() > 0) ++hits;
  }
  double t_us = t.getElapsedTimeInMicroSec();
  std::size_t n = tf.size();
  printf("%-20s collide  total %.0f us  bv/q %6.1f  leaf/q %6.1f  hits %zu/%zu\n",
         name, t_us, double(total_bv) / n, double(total_leaf) / n, hits, n);
}

template <typename BV, typename TraversalNode>
void runDistance(const std::vector<Transform3s>& tf, const BVHModel<BV>& m1,
                 const BVHModel<BV>& m2, const char* name) {
  Transform3s pose2;
  DistanceRequest request(true);
  DistanceResult local_result;
  TraversalNode node;
  node.enable_statistics = true;

  std::size_t total_bv = 0, total_leaf = 0;
  BenchTimer t;
  t.start();
  for (std::size_t i = 0; i < tf.size(); ++i) {
    if (!initialize(node, m1, tf[i], m2, pose2, request, local_result)) {
      std::cout << "init err\n";
    }
    node.num_bv_tests = 0;
    node.num_leaf_tests = 0;
    distance(&node, NULL);
    total_bv += static_cast<std::size_t>(node.num_bv_tests);
    total_leaf += static_cast<std::size_t>(node.num_leaf_tests);
  }
  double t_us = t.getElapsedTimeInMicroSec();
  std::size_t n = tf.size();
  printf("%-20s distance total %.0f us  bv/q %6.1f  leaf/q %6.1f\n",
         name, t_us, double(total_bv) / n, double(total_leaf) / n);
}

int main(int, char*[]) {
  std::vector<Vec3s> p1, p2;
  std::vector<Triangle32> t1, t2;
  boost::filesystem::path path(TEST_RESOURCES_DIR);
  loadOBJFile((path / "env.obj").string().c_str(), p1, t1);
  loadOBJFile((path / "rob.obj").string().c_str(), p2, t2);

  printf("env: %zu verts %zu tris  rob: %zu verts %zu tris\n",
         p1.size(), t1.size(), p2.size(), t2.size());

  // Report bounding boxes so the workload extents make sense.
  auto bbox = [](const std::vector<Vec3s>& v) {
    Vec3s lo = v[0], hi = v[0];
    for (const auto& p : v) {
      lo = lo.cwiseMin(p);
      hi = hi.cwiseMax(p);
    }
    printf("  min %.1f %.1f %.1f  max %.1f %.1f %.1f  size %.1f %.1f %.1f\n",
           lo[0], lo[1], lo[2], hi[0], hi[1], hi[2],
           hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]);
  };
  printf("env bbox:\n"); bbox(p1);
  printf("rob bbox:\n"); bbox(p2);

  BVHModel<OBB> m_obb1, m_obb2;
  BVHModel<RSS> m_rss1, m_rss2;
  BVHModel<OBBRSS> m_obbrss1, m_obbrss2;
  BVHModel<kIOS> m_kios1, m_kios2;
  makeModel(p1, t1, m_obb1);
  makeModel(p2, t2, m_obb2);
  makeModel(p1, t1, m_rss1);
  makeModel(p2, t2, m_rss2);
  makeModel(p1, t1, m_obbrss1);
  makeModel(p2, t2, m_obbrss2);
  makeModel(p1, t1, m_kios1);
  makeModel(p2, t2, m_kios2);

  // The default benchmark spreads transforms uniformly over [-3000, 3000]
  // which makes ~all queries reject at the root BV; that is broad-phase
  // dominated and tells us nothing about narrow-phase BV-test cost. Run two
  // workloads: the wide ("disjoint") one and a tight ("near") one that
  // forces BVH descent.
  struct Workload {
    const char* name;
    Scalar half_extent;
  };
  Workload workloads[] = {
      {"wide [-3000,3000]", 3000},
      {"med  [-300,300]  ", 300},
      {"near [-30,30]    ", 30},
      {"tight [-3,3]     ", 3},
  };

  // env-vs-rob with the original wide extents (broad-phase regime, baseline).
  for (const auto& w : workloads) {
    printf("\n=== env-vs-rob workload: %s ===\n", w.name);
    std::vector<Transform3s> tf;
    Scalar extents[] = {-w.half_extent, -w.half_extent, -w.half_extent,
                         w.half_extent,  w.half_extent,  w.half_extent};
    srand(42);
    generateRandomTransforms(extents, tf, 10000);

    runCollide<OBB, MeshCollisionTraversalNodeOBB>(tf, m_obb1, m_obb2, "OBB");
    runCollide<RSS, MeshCollisionTraversalNodeRSS>(tf, m_rss1, m_rss2, "RSS");
    runCollide<OBBRSS, MeshCollisionTraversalNodeOBBRSS>(tf, m_obbrss1,
                                                        m_obbrss2, "OBBRSS");
    runCollide<kIOS, MeshCollisionTraversalNodekIOS>(tf, m_kios1, m_kios2,
                                                    "kIOS");

    runDistance<RSS, MeshDistanceTraversalNodeRSS>(tf, m_rss1, m_rss2, "RSS");
    runDistance<OBBRSS, MeshDistanceTraversalNodeOBBRSS>(tf, m_obbrss1,
                                                        m_obbrss2, "OBBRSS");
    runDistance<kIOS, MeshDistanceTraversalNodekIOS>(tf, m_kios1, m_kios2,
                                                    "kIOS");
  }

  // rob-vs-rob with translations smaller than the rob bbox (1000x1210x850)
  // so the two copies actually overlap and the BVH descends. This is the
  // narrow-phase regime that exercises the BV-overlap function.
  BVHModel<OBB> m_obb_rob1, m_obb_rob2;
  BVHModel<RSS> m_rss_rob1, m_rss_rob2;
  BVHModel<OBBRSS> m_obbrss_rob1, m_obbrss_rob2;
  BVHModel<kIOS> m_kios_rob1, m_kios_rob2;
  makeModel(p2, t2, m_obb_rob1);
  makeModel(p2, t2, m_obb_rob2);
  makeModel(p2, t2, m_rss_rob1);
  makeModel(p2, t2, m_rss_rob2);
  makeModel(p2, t2, m_obbrss_rob1);
  makeModel(p2, t2, m_obbrss_rob2);
  makeModel(p2, t2, m_kios_rob1);
  makeModel(p2, t2, m_kios_rob2);

  Workload rob_workloads[] = {
      {"overlap [-100,100]", 100},
      {"overlap [-300,300]", 300},
      {"overlap [-600,600]", 600},
  };
  for (const auto& w : rob_workloads) {
    printf("\n=== rob-vs-rob workload: %s ===\n", w.name);
    std::vector<Transform3s> tf;
    Scalar extents[] = {-w.half_extent, -w.half_extent, -w.half_extent,
                         w.half_extent,  w.half_extent,  w.half_extent};
    srand(42);
    generateRandomTransforms(extents, tf, 10000);

    runCollide<OBB, MeshCollisionTraversalNodeOBB>(tf, m_obb_rob1, m_obb_rob2,
                                                  "OBB");
    runCollide<RSS, MeshCollisionTraversalNodeRSS>(tf, m_rss_rob1, m_rss_rob2,
                                                  "RSS");
    runCollide<OBBRSS, MeshCollisionTraversalNodeOBBRSS>(
        tf, m_obbrss_rob1, m_obbrss_rob2, "OBBRSS");
    runCollide<kIOS, MeshCollisionTraversalNodekIOS>(tf, m_kios_rob1,
                                                    m_kios_rob2, "kIOS");
    runDistance<RSS, MeshDistanceTraversalNodeRSS>(tf, m_rss_rob1, m_rss_rob2,
                                                  "RSS");
    runDistance<OBBRSS, MeshDistanceTraversalNodeOBBRSS>(
        tf, m_obbrss_rob1, m_obbrss_rob2, "OBBRSS");
    runDistance<kIOS, MeshDistanceTraversalNodekIOS>(tf, m_kios_rob1,
                                                    m_kios_rob2, "kIOS");
  }

  // Synthetic dense-overlap workload: two random triangle soups in the same
  // [-1,1] cube, ~equal density, so BVHs are guaranteed to descend deeply.
  {
    auto makeSoup = [](unsigned seed, int N, std::vector<Vec3s>& pts,
                       std::vector<Triangle32>& tris) {
      pts.clear();
      tris.clear();
      pts.reserve(static_cast<size_t>(3 * N));
      tris.reserve(static_cast<size_t>(N));
      std::mt19937 rng(seed);
      std::uniform_real_distribution<Scalar> p(Scalar(-1.0), Scalar(1.0));
      std::uniform_real_distribution<Scalar> e(Scalar(-0.4), Scalar(0.4));
      for (int i = 0; i < N; ++i) {
        Vec3s c(p(rng), p(rng), p(rng));
        pts.push_back(Vec3s(c[0] + e(rng), c[1] + e(rng), c[2] + e(rng)));
        pts.push_back(Vec3s(c[0] + e(rng), c[1] + e(rng), c[2] + e(rng)));
        pts.push_back(Vec3s(c[0] + e(rng), c[1] + e(rng), c[2] + e(rng)));
        tris.push_back(Triangle32(static_cast<unsigned>(3 * i),
                                  static_cast<unsigned>(3 * i + 1),
                                  static_cast<unsigned>(3 * i + 2)));
      }
    };

    std::vector<Vec3s> pa, pb;
    std::vector<Triangle32> ta, tb;
    makeSoup(0xA1u, 2000, pa, ta);
    makeSoup(0xB2u, 2000, pb, tb);

    BVHModel<OBB> sa_obb, sb_obb;
    BVHModel<RSS> sa_rss, sb_rss;
    BVHModel<OBBRSS> sa_obbrss, sb_obbrss;
    BVHModel<kIOS> sa_kios, sb_kios;
    makeModel(pa, ta, sa_obb);
    makeModel(pb, tb, sb_obb);
    makeModel(pa, ta, sa_rss);
    makeModel(pb, tb, sb_rss);
    makeModel(pa, ta, sa_obbrss);
    makeModel(pb, tb, sb_obbrss);
    makeModel(pa, ta, sa_kios);
    makeModel(pb, tb, sb_kios);

    std::vector<Transform3s> tf;
    Scalar extents[] = {-Scalar(0.1), -Scalar(0.1), -Scalar(0.1),
                         Scalar(0.1),  Scalar(0.1),  Scalar(0.1)};
    srand(42);
    generateRandomTransforms(extents, tf, 1000);  // fewer queries; each is heavy

    // Pathological: same mesh on both sides, near-identity transforms.
    // Forces every BV pair to overlap and the BVH to descend to leaves.
    BVHModel<OBB> sx_obb;
    BVHModel<RSS> sx_rss;
    BVHModel<OBBRSS> sx_obbrss;
    BVHModel<kIOS> sx_kios;
    makeModel(pa, ta, sx_obb);
    makeModel(pa, ta, sx_rss);
    makeModel(pa, ta, sx_obbrss);
    makeModel(pa, ta, sx_kios);
    std::vector<Transform3s> tf_path;
    Scalar tiny[] = {-Scalar(1e-3), -Scalar(1e-3), -Scalar(1e-3),
                      Scalar(1e-3),  Scalar(1e-3),  Scalar(1e-3)};
    srand(123);
    generateRandomTransforms(tiny, tf_path, 100);

    printf("\n=== synthetic dense soup [-0.1,0.1] (1000 queries) ===\n");
    runCollide<OBB, MeshCollisionTraversalNodeOBB>(tf, sa_obb, sb_obb, "OBB");
    runCollide<RSS, MeshCollisionTraversalNodeRSS>(tf, sa_rss, sb_rss, "RSS");
    runCollide<OBBRSS, MeshCollisionTraversalNodeOBBRSS>(tf, sa_obbrss,
                                                        sb_obbrss, "OBBRSS");
    runCollide<kIOS, MeshCollisionTraversalNodekIOS>(tf, sa_kios, sb_kios,
                                                    "kIOS");
    runDistance<RSS, MeshDistanceTraversalNodeRSS>(tf, sa_rss, sb_rss, "RSS");
    runDistance<OBBRSS, MeshDistanceTraversalNodeOBBRSS>(tf, sa_obbrss,
                                                        sb_obbrss, "OBBRSS");
    runDistance<kIOS, MeshDistanceTraversalNodekIOS>(tf, sa_kios, sb_kios,
                                                    "kIOS");

    printf("\n=== pathological self-vs-self near-identity (100 queries) ===\n");
    runCollide<OBB, MeshCollisionTraversalNodeOBB>(tf_path, sx_obb, sx_obb,
                                                  "OBB");
    runCollide<RSS, MeshCollisionTraversalNodeRSS>(tf_path, sx_rss, sx_rss,
                                                  "RSS");
    runCollide<OBBRSS, MeshCollisionTraversalNodeOBBRSS>(tf_path, sx_obbrss,
                                                        sx_obbrss, "OBBRSS");
    runCollide<kIOS, MeshCollisionTraversalNodekIOS>(tf_path, sx_kios, sx_kios,
                                                    "kIOS");
    runDistance<RSS, MeshDistanceTraversalNodeRSS>(tf_path, sx_rss, sx_rss,
                                                  "RSS");
    runDistance<OBBRSS, MeshDistanceTraversalNodeOBBRSS>(tf_path, sx_obbrss,
                                                        sx_obbrss, "OBBRSS");
    runDistance<kIOS, MeshDistanceTraversalNodekIOS>(tf_path, sx_kios, sx_kios,
                                                    "kIOS");
  }

  return 0;
}
