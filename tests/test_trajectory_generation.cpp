#include <iostream>
#include <iomanip>
#include <vector>
#include <armadillo>

#include <catch2/catch_all.hpp>

#include "modern_robotics/trajectory_generation.hpp"

constexpr double TOLERANCE = 1e-3;

TEST_CASE("Testing cubic time scaling", "[CubicTimeScaling]")
{
  const double Tf = 2.0;
  const double t = 0.6;

  const double st = mr::CubicTimeScaling(Tf, t);
  // std::cout << st << std::endl;

  REQUIRE_THAT(st, Catch::Matchers::WithinAbs(0.216, TOLERANCE));
}

TEST_CASE("Testing quintic time scaling", "[QuinticTimeScaling]")
{
  const double Tf = 2.0;
  const double t = 0.6;

  const double st = mr::QuinticTimeScaling(Tf, t);
  // std::cout << st << std::endl;

  REQUIRE_THAT(st, Catch::Matchers::WithinAbs(0.16308, TOLERANCE));
}

TEST_CASE("Testing joint trajectory", "[JointTrajectory]")
{
  const arma::vec thetastart{1, 0, 0, 1, 1, 0.2, 0, 1};
  const arma::vec thetaend{1.2, 0.5, 0.6, 1.1, 2, 2, 0.9, 1};
  const double Tf = 4.0;
  const size_t N = 6;
  const mr::Method method = mr::Method::Cubic;

  const std::vector<arma::vec> traj = mr::JointTrajectory(thetastart, thetaend, Tf, N, method);
  // for (const auto & t : traj) {
  //   std::cout << t << std::endl;
  // }

  REQUIRE(traj.size() == N);
  REQUIRE_THAT(traj.at(0).at(0), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(3), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(4), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(5), Catch::Matchers::WithinAbs(0.2, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(6), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(7), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(0), Catch::Matchers::WithinAbs(1.0208, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(1), Catch::Matchers::WithinAbs(0.052, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(2), Catch::Matchers::WithinAbs(0.0624, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(3), Catch::Matchers::WithinAbs(1.0104, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(4), Catch::Matchers::WithinAbs(1.104, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(5), Catch::Matchers::WithinAbs(0.3872, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(6), Catch::Matchers::WithinAbs(0.0936, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(7), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(0), Catch::Matchers::WithinAbs(1.0704, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(1), Catch::Matchers::WithinAbs(0.176, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(2), Catch::Matchers::WithinAbs(0.2112, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(3), Catch::Matchers::WithinAbs(1.0352, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(4), Catch::Matchers::WithinAbs(1.352, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(5), Catch::Matchers::WithinAbs(0.8336, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(6), Catch::Matchers::WithinAbs(0.3168, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(7), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(0), Catch::Matchers::WithinAbs(1.1296, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(1), Catch::Matchers::WithinAbs(0.324, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(2), Catch::Matchers::WithinAbs(0.3888, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(3), Catch::Matchers::WithinAbs(1.0648, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(4), Catch::Matchers::WithinAbs(1.648, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(5), Catch::Matchers::WithinAbs(1.3664, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(6), Catch::Matchers::WithinAbs(0.5832, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(7), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(4).at(0), Catch::Matchers::WithinAbs(1.1792, TOLERANCE));
  REQUIRE_THAT(traj.at(4).at(1), Catch::Matchers::WithinAbs(0.448, TOLERANCE));
  REQUIRE_THAT(traj.at(4).at(2), Catch::Matchers::WithinAbs(0.5376, TOLERANCE));
  REQUIRE_THAT(traj.at(4).at(3), Catch::Matchers::WithinAbs(1.0896, TOLERANCE));
  REQUIRE_THAT(traj.at(4).at(4), Catch::Matchers::WithinAbs(1.896, TOLERANCE));
  REQUIRE_THAT(traj.at(4).at(5), Catch::Matchers::WithinAbs(1.8128, TOLERANCE));
  REQUIRE_THAT(traj.at(4).at(6), Catch::Matchers::WithinAbs(0.8064, TOLERANCE));
  REQUIRE_THAT(traj.at(4).at(7), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(5).at(0), Catch::Matchers::WithinAbs(1.2, TOLERANCE));
  REQUIRE_THAT(traj.at(5).at(1), Catch::Matchers::WithinAbs(0.5, TOLERANCE));
  REQUIRE_THAT(traj.at(5).at(2), Catch::Matchers::WithinAbs(0.6, TOLERANCE));
  REQUIRE_THAT(traj.at(5).at(3), Catch::Matchers::WithinAbs(1.1, TOLERANCE));
  REQUIRE_THAT(traj.at(5).at(4), Catch::Matchers::WithinAbs(2, TOLERANCE));
  REQUIRE_THAT(traj.at(5).at(5), Catch::Matchers::WithinAbs(2, TOLERANCE));
  REQUIRE_THAT(traj.at(5).at(6), Catch::Matchers::WithinAbs(0.9, TOLERANCE));
  REQUIRE_THAT(traj.at(5).at(7), Catch::Matchers::WithinAbs(1, TOLERANCE));
}

TEST_CASE("Testing screw trajectory", "[ScrewTrajectory]")
{
  const arma::mat44 Xstart{
    {1, 0, 0, 1},
    {0, 1, 0, 0},
    {0, 0, 1, 1},
    {0, 0, 0, 1}
  };
  const arma::mat44 Xend{
    {0, 0, 1, 0.1},
    {1, 0, 0, 0},
    {0, 1, 0, 4.1},
    {0, 0, 0, 1}
  };
  const double Tf = 5.0;
  const size_t N = 4;
  const mr::Method method = mr::Method::Cubic;

  const std::vector<arma::mat44> traj = mr::ScrewTrajectory(Xstart, Xend, Tf, N, method);
  // for (const auto & t : traj) {
  //   std::cout << t << std::endl;
  // }

  REQUIRE(traj.size() == N);
  REQUIRE_THAT(traj.at(0).at(0, 0), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(0, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(0, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(0, 3), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(1, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(1, 1), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(1, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(1, 3), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(2, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(2, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(2, 2), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(2, 3), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(3, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(3, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(3, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(3, 3), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(0, 0), Catch::Matchers::WithinAbs(0.904, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(0, 1), Catch::Matchers::WithinAbs(-0.25, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(0, 2), Catch::Matchers::WithinAbs(0.346, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(0, 3), Catch::Matchers::WithinAbs(0.441, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(1, 0), Catch::Matchers::WithinAbs(0.346, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(1, 1), Catch::Matchers::WithinAbs(0.904, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(1, 2), Catch::Matchers::WithinAbs(-0.25, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(1, 3), Catch::Matchers::WithinAbs(0.529, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(2, 0), Catch::Matchers::WithinAbs(-0.25, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(2, 1), Catch::Matchers::WithinAbs(0.346, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(2, 2), Catch::Matchers::WithinAbs(0.904, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(2, 3), Catch::Matchers::WithinAbs(1.601, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(3, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(3, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(3, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(3, 3), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(0, 0), Catch::Matchers::WithinAbs(0.346, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(0, 1), Catch::Matchers::WithinAbs(-0.25, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(0, 2), Catch::Matchers::WithinAbs(0.904, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(0, 3), Catch::Matchers::WithinAbs(-0.117, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(1, 0), Catch::Matchers::WithinAbs(0.904, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(1, 1), Catch::Matchers::WithinAbs(0.346, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(1, 2), Catch::Matchers::WithinAbs(-0.25, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(1, 3), Catch::Matchers::WithinAbs(0.473, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(2, 0), Catch::Matchers::WithinAbs(-0.25, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(2, 1), Catch::Matchers::WithinAbs(0.904, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(2, 2), Catch::Matchers::WithinAbs(0.346, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(2, 3), Catch::Matchers::WithinAbs(3.274, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(3, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(3, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(3, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(3, 3), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(0, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(0, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(0, 2), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(0, 3), Catch::Matchers::WithinAbs(0.1, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(1, 0), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(1, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(1, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(1, 3), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(2, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(2, 1), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(2, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(2, 3), Catch::Matchers::WithinAbs(4.1, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(3, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(3, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(3, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(3, 3), Catch::Matchers::WithinAbs(1, TOLERANCE));
}

TEST_CASE("Testing cartesian trajectory", "[CartesianTrajectory]")
{
  const arma::mat44 Xstart{
    {1, 0, 0, 1},
    {0, 1, 0, 0},
    {0, 0, 1, 1},
    {0, 0, 0, 1}
  };
  const arma::mat44 Xend{
    {0, 0, 1, 0.1},
    {1, 0, 0, 0},
    {0, 1, 0, 4.1},
    {0, 0, 0, 1}
  };
  const double Tf = 5.0;
  const size_t N = 4;
  const mr::Method method = mr::Method::Quintic;

  const std::vector<arma::mat44> traj = mr::CartesianTrajectory(Xstart, Xend, Tf, N, method);
  // for (const auto & t : traj) {
  //   std::cout << t << std::endl;
  // }

  REQUIRE(traj.size() == N);
  REQUIRE_THAT(traj.at(0).at(0, 0), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(0, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(0, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(0, 3), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(1, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(1, 1), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(1, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(1, 3), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(2, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(2, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(2, 2), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(2, 3), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(3, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(3, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(3, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(0).at(3, 3), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(0, 0), Catch::Matchers::WithinAbs(0.937, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(0, 1), Catch::Matchers::WithinAbs(-0.214, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(0, 2), Catch::Matchers::WithinAbs(0.277, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(0, 3), Catch::Matchers::WithinAbs(0.811, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(1, 0), Catch::Matchers::WithinAbs(0.277, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(1, 1), Catch::Matchers::WithinAbs(0.937, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(1, 2), Catch::Matchers::WithinAbs(-0.214, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(1, 3), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(2, 0), Catch::Matchers::WithinAbs(-0.214, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(2, 1), Catch::Matchers::WithinAbs(0.277, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(2, 2), Catch::Matchers::WithinAbs(0.937, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(2, 3), Catch::Matchers::WithinAbs(1.651, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(3, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(3, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(3, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(1).at(3, 3), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(0, 0), Catch::Matchers::WithinAbs(0.277, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(0, 1), Catch::Matchers::WithinAbs(-0.214, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(0, 2), Catch::Matchers::WithinAbs(0.937, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(0, 3), Catch::Matchers::WithinAbs(0.289, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(1, 0), Catch::Matchers::WithinAbs(0.937, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(1, 1), Catch::Matchers::WithinAbs(0.277, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(1, 2), Catch::Matchers::WithinAbs(-0.214, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(1, 3), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(2, 0), Catch::Matchers::WithinAbs(-0.214, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(2, 1), Catch::Matchers::WithinAbs(0.937, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(2, 2), Catch::Matchers::WithinAbs(0.277, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(2, 3), Catch::Matchers::WithinAbs(3.449, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(3, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(3, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(3, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(2).at(3, 3), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(0, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(0, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(0, 2), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(0, 3), Catch::Matchers::WithinAbs(0.1, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(1, 0), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(1, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(1, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(1, 3), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(2, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(2, 1), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(2, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(2, 3), Catch::Matchers::WithinAbs(4.1, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(3, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(3, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(3, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(traj.at(3).at(3, 3), Catch::Matchers::WithinAbs(1, TOLERANCE));
}
