#include <armadillo>
#include <catch2/catch_all.hpp>
#include <tuple>

#include "modern_robotics/rigid_body_motions.hpp"

constexpr double TOLERANCE = 1e-6;

TEST_CASE("Test Inverting Rotational matrix", "[RotInv]")
{
  const arma::mat33 R{
    {0, 0, 1},
    {1, 0, 0},
    {0, 1, 0}
  };
  const arma::mat33 R_inv = modern_robotics::RotInv(R);
  // std::cout << R_inv << std::endl;

  REQUIRE_THAT(R_inv.at(0, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(R_inv.at(0, 1), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(R_inv.at(0, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(R_inv.at(1, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(R_inv.at(1, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(R_inv.at(1, 2), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(R_inv.at(2, 0), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(R_inv.at(2, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(R_inv.at(2, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
}

TEST_CASE("Testing vector3 to so3", "[VecToso3]")
{
  const arma::vec3 omg{1, 2, 3};
  const arma::mat33 so3mat = modern_robotics::VecToso3(omg);
  // std::cout << so3mat << std::endl;

  REQUIRE_THAT(so3mat.at(0, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(so3mat.at(0, 1), Catch::Matchers::WithinAbs(-3, TOLERANCE));
  REQUIRE_THAT(so3mat.at(0, 2), Catch::Matchers::WithinAbs(2, TOLERANCE));
  REQUIRE_THAT(so3mat.at(1, 0), Catch::Matchers::WithinAbs(3, TOLERANCE));
  REQUIRE_THAT(so3mat.at(1, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(so3mat.at(1, 2), Catch::Matchers::WithinAbs(-1, TOLERANCE));
  REQUIRE_THAT(so3mat.at(2, 0), Catch::Matchers::WithinAbs(-2, TOLERANCE));
  REQUIRE_THAT(so3mat.at(2, 1), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(so3mat.at(2, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
}

TEST_CASE("Testing so3 to vec3", "[so3ToVec]")
{
  const arma::mat33 so3mat{
    {0, -3, 2},
    {3, 0, -1},
    {-2, 1, 0}
  };
  const arma::vec3 omg = modern_robotics::so3ToVec(so3mat);
  // std::cout << so3mat << std::endl;

  REQUIRE_THAT(omg.at(0), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(omg.at(1), Catch::Matchers::WithinAbs(2, TOLERANCE));
  REQUIRE_THAT(omg.at(2), Catch::Matchers::WithinAbs(3, TOLERANCE));
}

TEST_CASE("Tesing axis angle", "[AxisAng]")
{
  const arma::vec3 expc3{1, 2, 3};
  const auto &[omghat, theta] = modern_robotics::AxisAng(expc3);

  REQUIRE_THAT(
    omghat.at(0),
    Catch::Matchers::WithinAbs(
      1.0 / std::sqrt(1 * 1 + 2 * 2 + 3 * 3),
      TOLERANCE
    )
  );
  REQUIRE_THAT(
    omghat.at(1),
    Catch::Matchers::WithinAbs(
      2.0 / std::sqrt(1 * 1 + 2 * 2 + 3 * 3),
      TOLERANCE
    )
  );
  REQUIRE_THAT(
    omghat.at(2),
    Catch::Matchers::WithinAbs(
      3.0 / std::sqrt(1 * 1 + 2 * 2 + 3 * 3),
      TOLERANCE
    )
  );
  REQUIRE_THAT(theta, Catch::Matchers::WithinAbs(arma::norm(expc3), TOLERANCE));
}
