#include <catch2/catch_all.hpp>

#include "modern_robotics/forward_kinematics.hpp"

constexpr double TOLERANCE = 1e-6;

TEST_CASE("Test forward kinematic body", "[FKinBody]")
{
  const arma::mat44 M{
    {-1, 0, 0, 0},
    {0, 1, 0, 6},
    {0, 0, -1, 2},
    {0, 0, 0, 1}
  };
  const std::vector<arma::vec6> Blist{
    {0, 0, -1, 2, 0, 0},
    {0, 0, 0, 0, 1, 0},
    {0, 0, 1, 0, 0, 0.1}
  };
  const std::vector<double> thetalist{M_PI_2, 3, M_PI};

  const arma::mat44 T = mr::FKinBody(M, Blist, thetalist);
  // std::cout << T << std::endl;

  REQUIRE_THAT(T.at(0, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(0, 1), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(T.at(0, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(0, 3), Catch::Matchers::WithinAbs(-5, TOLERANCE));
  REQUIRE_THAT(T.at(1, 0), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(T.at(1, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(1, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(1, 3), Catch::Matchers::WithinAbs(4, TOLERANCE));
  REQUIRE_THAT(T.at(2, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(2, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(2, 2), Catch::Matchers::WithinAbs(-1, TOLERANCE));
  REQUIRE_THAT(T.at(2, 3), Catch::Matchers::WithinAbs(1.68584073, TOLERANCE));
  REQUIRE_THAT(T.at(3, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(3, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(3, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(3, 3), Catch::Matchers::WithinAbs(1, TOLERANCE));
}

TEST_CASE("Test forward kinematics on space", "[FKinSpace]")
{
  const arma::mat44 M{
    {-1, 0, 0, 0},
    {0, 1, 0, 6},
    {0, 0, -1, 2},
    {0, 0, 0, 1}
  };
  const std::vector<arma::vec6> Slist{
    {0, 0, 1, 4, 0, 0},
    {0, 0, 0, 0, 1, 0},
    {0, 0, -1, -6, 0, -0.1}
  };
  const std::vector<double> thetalist{M_PI_2, 3, M_PI};

  const arma::mat44 T = mr::FKinSpace(M, Slist, thetalist);
  // std::cout << T << std::endl;

  REQUIRE_THAT(T.at(0, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(0, 1), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(T.at(0, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(0, 3), Catch::Matchers::WithinAbs(-5, TOLERANCE));
  REQUIRE_THAT(T.at(1, 0), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(T.at(1, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(1, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(1, 3), Catch::Matchers::WithinAbs(4, TOLERANCE));
  REQUIRE_THAT(T.at(2, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(2, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(2, 2), Catch::Matchers::WithinAbs(-1, TOLERANCE));
  REQUIRE_THAT(T.at(2, 3), Catch::Matchers::WithinAbs(1.68584073, TOLERANCE));
  REQUIRE_THAT(T.at(3, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(3, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(3, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(3, 3), Catch::Matchers::WithinAbs(1, TOLERANCE));
}
