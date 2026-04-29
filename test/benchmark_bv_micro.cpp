// Per-call micro-benchmark for the inner BV overlap functions.
// Measures rectDistance (used by RSS) and obbDisjointAndLowerBoundDistance
// (used by OBB) across realistic input distributions, in tight loops,
// without BVH-traversal noise.

#include <chrono>
#include <cstdio>
#include <random>

#include "coal/BV/OBB.h"
#include "coal/BV/RSS.h"
#include "coal/collision_data.h"

namespace coal {
// rectDistance is in an anonymous namespace inside RSS.cpp; not callable
// directly. We benchmark by way of the RSS::overlap path, which is the
// single-call form. That includes the wrapping subtraction of radii but
// no transforms, so it isolates rectDistance.
}

using namespace coal;

namespace {

// Small fixture set: precomputed Rab/Tab pairs that exercise different
// regions of rectDistance. We pre-generate to avoid randomness in the timed
// loop.
struct RsCase {
  Matrix3s Rab;
  Vec3s Tab;
  Scalar a[2];
  Scalar b[2];
};

std::vector<RsCase> makeRectCases(int count, Scalar trans_extent) {
  std::mt19937 rng(0xC0A1u);
  std::uniform_real_distribution<Scalar> tdist(-trans_extent, trans_extent);
  std::uniform_real_distribution<Scalar> adist(Scalar(-EIGEN_PI), Scalar(EIGEN_PI));
  std::uniform_real_distribution<Scalar> sdist(Scalar(0.5), Scalar(2.0));
  std::vector<RsCase> v;
  v.reserve(static_cast<size_t>(count));
  for (int i = 0; i < count; ++i) {
    RsCase c;
    Scalar rx = adist(rng), ry = adist(rng), rz = adist(rng);
    Eigen::AngleAxis<Scalar> Rx(rx, Vec3s::UnitX());
    Eigen::AngleAxis<Scalar> Ry(ry, Vec3s::UnitY());
    Eigen::AngleAxis<Scalar> Rz(rz, Vec3s::UnitZ());
    c.Rab = (Rx * Ry * Rz).matrix();
    c.Tab = Vec3s(tdist(rng), tdist(rng), tdist(rng));
    c.a[0] = sdist(rng);
    c.a[1] = sdist(rng);
    c.b[0] = sdist(rng);
    c.b[1] = sdist(rng);
    v.push_back(c);
  }
  return v;
}

struct ObbCase {
  Matrix3s B;
  Vec3s T;
  Vec3s a;
  Vec3s b;
};

std::vector<ObbCase> makeObbCases(int count, Scalar trans_extent) {
  std::mt19937 rng(0xC0BBu);
  std::uniform_real_distribution<Scalar> tdist(-trans_extent, trans_extent);
  std::uniform_real_distribution<Scalar> adist(Scalar(-EIGEN_PI), Scalar(EIGEN_PI));
  std::uniform_real_distribution<Scalar> sdist(Scalar(0.5), Scalar(2.0));
  std::vector<ObbCase> v;
  v.reserve(static_cast<size_t>(count));
  for (int i = 0; i < count; ++i) {
    ObbCase c;
    Scalar rx = adist(rng), ry = adist(rng), rz = adist(rng);
    Eigen::AngleAxis<Scalar> Rx(rx, Vec3s::UnitX());
    Eigen::AngleAxis<Scalar> Ry(ry, Vec3s::UnitY());
    Eigen::AngleAxis<Scalar> Rz(rz, Vec3s::UnitZ());
    c.B = (Rx * Ry * Rz).matrix();
    c.T = Vec3s(tdist(rng), tdist(rng), tdist(rng));
    c.a = Vec3s(sdist(rng), sdist(rng), sdist(rng));
    c.b = Vec3s(sdist(rng), sdist(rng), sdist(rng));
    v.push_back(c);
  }
  return v;
}

}  // namespace

int main(int, char*[]) {
  using namespace std::chrono;

  // ============================================================
  //  RSS::overlap micro-benchmark — exercises rectDistance internally.
  //  Inputs span trans_extent values to cover broad-phase (far) and
  //  narrow-phase (descending) regimes.
  // ============================================================
  printf("# RSS::overlap (one call per iteration)\n");
  printf("%-15s %12s %14s\n", "trans_extent", "calls/sec", "ns/call");
  for (Scalar te : {Scalar(0.5), Scalar(2.0), Scalar(5.0), Scalar(20.0)}) {
    auto cases = makeRectCases(1024, te);

    // Two RSS objects sharing radius=0.1 (so overlap = rectDistance < 0.2)
    RSS rss1, rss2;
    rss1.radius = Scalar(0.1);
    rss2.radius = Scalar(0.1);
    CollisionRequest req;

    // Warm up
    Scalar acc = 0;
    for (int i = 0; i < 16; ++i) {
      const auto& c = cases[i % cases.size()];
      rss1.length[0] = c.a[0];
      rss1.length[1] = c.a[1];
      rss2.length[0] = c.b[0];
      rss2.length[1] = c.b[1];
      // Place rss2 at (Rab,Tab) relative to rss1.
      rss1.axes.setIdentity();
      rss1.Tr.setZero();
      rss2.axes = c.Rab;
      rss2.Tr = c.Tab;
      Scalar lb = 0;
      bool o = overlap(Matrix3s::Identity(), Vec3s::Zero(), rss1, rss2, req, lb);
      acc += o ? Scalar(1) : Scalar(0);
    }

    const int reps = 1000;
    auto t0 = steady_clock::now();
    for (int r = 0; r < reps; ++r) {
      for (const auto& c : cases) {
        rss1.length[0] = c.a[0];
        rss1.length[1] = c.a[1];
        rss2.length[0] = c.b[0];
        rss2.length[1] = c.b[1];
        rss1.axes.setIdentity();
        rss1.Tr.setZero();
        rss2.axes = c.Rab;
        rss2.Tr = c.Tab;
        Scalar lb = 0;
        bool o = overlap(Matrix3s::Identity(), Vec3s::Zero(), rss1, rss2, req, lb);
        acc += o ? Scalar(1) : Scalar(0);
      }
    }
    auto t1 = steady_clock::now();
    double sec = duration<double>(t1 - t0).count();
    double calls = double(reps) * double(cases.size());
    printf("%-15g %12.0f %14.2f   sink=%g\n", te, calls / sec,
           sec * 1e9 / calls, double(acc));
  }

  // ============================================================
  //  OBB SAT micro-benchmark — exercises obbDisjointAndLowerBoundDistance.
  // ============================================================
  printf("\n# OBB::overlap (one call per iteration)\n");
  printf("%-15s %12s %14s\n", "trans_extent", "calls/sec", "ns/call");
  for (Scalar te : {Scalar(0.5), Scalar(2.0), Scalar(5.0), Scalar(20.0)}) {
    auto cases = makeObbCases(1024, te);
    OBB ob1, ob2;
    CollisionRequest req;

    // Warm up
    Scalar acc = 0;
    for (int i = 0; i < 16; ++i) {
      const auto& c = cases[i % cases.size()];
      ob1.axes.setIdentity();
      ob1.To.setZero();
      ob1.extent = c.a;
      ob2.axes = c.B;
      ob2.To = c.T;
      ob2.extent = c.b;
      Scalar lb = 0;
      bool o = overlap(Matrix3s::Identity(), Vec3s::Zero(), ob1, ob2, req, lb);
      acc += o ? Scalar(1) : Scalar(0);
    }

    const int reps = 1000;
    auto t0 = steady_clock::now();
    for (int r = 0; r < reps; ++r) {
      for (const auto& c : cases) {
        ob1.axes.setIdentity();
        ob1.To.setZero();
        ob1.extent = c.a;
        ob2.axes = c.B;
        ob2.To = c.T;
        ob2.extent = c.b;
        Scalar lb = 0;
        bool o = overlap(Matrix3s::Identity(), Vec3s::Zero(), ob1, ob2, req, lb);
        acc += o ? Scalar(1) : Scalar(0);
      }
    }
    auto t1 = steady_clock::now();
    double sec = duration<double>(t1 - t0).count();
    double calls = double(reps) * double(cases.size());
    printf("%-15g %12.0f %14.2f   sink=%g\n", te, calls / sec,
           sec * 1e9 / calls, double(acc));
  }

  return 0;
}
