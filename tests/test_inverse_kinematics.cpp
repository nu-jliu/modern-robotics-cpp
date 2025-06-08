#include <iostream>
#include <vector>

#include <catch2/catch_all.hpp>

#include "modern_robotics/inverse_kinematics.hpp"

constexpr double TOLERANCE = 1e-3;

TEST_CASE("Test body inverse kinematics", "[IKinBody]")
{
  const std::vector<arma::vec6> Blist{
    {0, 0, -1, 2, 0, 0},
    {0, 0, 0, 0, 1, 0},
    {0, 0, 1, 0, 0, 0.1}
  };
  const arma::mat44 M{
    {-1, 0, 0, 0},
    {0, 1, 0, 6},
    {0, 0, -1, 2},
    {0, 0, 0, 1}
  };
  const arma::mat44 T{
    {0, 1, 0, -5},
    {1, 0, 0, 4},
    {0, 0, -1, 1.6858},
    {0, 0, 0, 1}
  };
  const arma::vec thetalist0{1.5, 2.5, 3};
  const double emog = TOLERANCE;
  const double ev = TOLERANCE;

  const auto &[thetalist, success] = mr::IKinBody(Blist, M, T, thetalist0, emog, ev);
  // std::cout << thetalist << std::endl;
  // std::cout << std::boolalpha << success << std::noboolalpha << std::endl;

  REQUIRE(success);
  REQUIRE(thetalist.size() == Blist.size());
  REQUIRE_THAT(thetalist.at(0), Catch::Matchers::WithinAbs(1.57073819, TOLERANCE));
  REQUIRE_THAT(thetalist.at(1), Catch::Matchers::WithinAbs(2.999667, TOLERANCE));
  REQUIRE_THAT(thetalist.at(2), Catch::Matchers::WithinAbs(3.14153913, TOLERANCE));
}

TEST_CASE("Testing space inverse kinematics", "[IKinSpace]")
{
  const std::vector<arma::vec6> Slist{
    {0, 0, 1, 4, 0, 0},
    {0, 0, 0, 0, 1, 0},
    {0, 0, -1, -6, 0, -0.1}
  };
  const arma::mat44 M{
    {-1, 0, 0, 0},
    {0, 1, 0, 6},
    {0, 0, -1, 2},
    {0, 0, 0, 1}
  };
  const arma::mat44 T{
    {0, 1, 0, -5},
    {1, 0, 0, 4},
    {0, 0, -1, 1.6858},
    {0, 0, 0, 1}
  };
  const arma::vec thetalist0{1.5, 2.5, 3};
  const double emog = TOLERANCE;
  const double ev = TOLERANCE;

  const auto &[thetalist, success] = mr::IKinSpace(Slist, M, T, thetalist0, emog, ev);
  // std::cout << thetalist << std::endl;
  // std::cout << std::boolalpha << success << std::noboolalpha << std::endl;

  REQUIRE(success);
  REQUIRE(thetalist.size() == Slist.size());
  REQUIRE_THAT(thetalist.at(0), Catch::Matchers::WithinAbs(1.57073819, TOLERANCE));
  REQUIRE_THAT(thetalist.at(1), Catch::Matchers::WithinAbs(2.999667, TOLERANCE));
  REQUIRE_THAT(thetalist.at(2), Catch::Matchers::WithinAbs(3.14153913, TOLERANCE));
}
