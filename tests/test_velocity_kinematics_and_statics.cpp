#include <catch2/catch_all.hpp>
#include "modern_robotics/velocity_kinematics_and_statics.hpp"

constexpr double TOLERANCE = 1e-6;

TEST_CASE("Test body jacobian", "[JacobianBody]")
{
  const std::vector<arma::vec6> Blist{
    {0, 0, 1, 0, 0.2, 0.2},
    {1, 0, 0, 2, 0, 3},
    {0, 1, 0, 0, 2, 1},
    {1, 0, 0, 0.2, 0.3, 0.4}
  };
  const std::vector<double> thetalist{0.2, 1.1, 0.1, 1.2};

  const arma::mat Jb = mr::JacobianBody(Blist, thetalist);
  // std::cout << Jb << std::endl;

  REQUIRE_THAT(Jb.at(0, 0), Catch::Matchers::WithinAbs(-0.04528405, TOLERANCE));
  REQUIRE_THAT(Jb.at(1, 0), Catch::Matchers::WithinAbs(0.74359313, TOLERANCE));
  REQUIRE_THAT(Jb.at(2, 0), Catch::Matchers::WithinAbs(-0.66709716, TOLERANCE));
  REQUIRE_THAT(Jb.at(3, 0), Catch::Matchers::WithinAbs(2.32586047, TOLERANCE));
  REQUIRE_THAT(Jb.at(4, 0), Catch::Matchers::WithinAbs(-1.44321167, TOLERANCE));
  REQUIRE_THAT(Jb.at(5, 0), Catch::Matchers::WithinAbs(-2.06639565, TOLERANCE));
  REQUIRE_THAT(Jb.at(0, 1), Catch::Matchers::WithinAbs(0.99500417, TOLERANCE));
  REQUIRE_THAT(Jb.at(1, 1), Catch::Matchers::WithinAbs(0.09304865, TOLERANCE));
  REQUIRE_THAT(Jb.at(2, 1), Catch::Matchers::WithinAbs(0.03617541, TOLERANCE));
  REQUIRE_THAT(Jb.at(3, 1), Catch::Matchers::WithinAbs(1.66809, TOLERANCE));
  REQUIRE_THAT(Jb.at(4, 1), Catch::Matchers::WithinAbs(2.94561275, TOLERANCE));
  REQUIRE_THAT(Jb.at(5, 1), Catch::Matchers::WithinAbs(1.82881722, TOLERANCE));
  REQUIRE_THAT(Jb.at(0, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(Jb.at(1, 2), Catch::Matchers::WithinAbs(0.36235775, TOLERANCE));
  REQUIRE_THAT(Jb.at(2, 2), Catch::Matchers::WithinAbs(-0.93203909, TOLERANCE));
  REQUIRE_THAT(Jb.at(3, 2), Catch::Matchers::WithinAbs(0.56410831, TOLERANCE));
  REQUIRE_THAT(Jb.at(4, 2), Catch::Matchers::WithinAbs(1.43306521, TOLERANCE));
  REQUIRE_THAT(Jb.at(5, 2), Catch::Matchers::WithinAbs(-1.58868628, TOLERANCE));
  REQUIRE_THAT(Jb.at(0, 3), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(Jb.at(1, 3), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(Jb.at(2, 3), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(Jb.at(3, 3), Catch::Matchers::WithinAbs(0.2, TOLERANCE));
  REQUIRE_THAT(Jb.at(4, 3), Catch::Matchers::WithinAbs(0.3, TOLERANCE));
  REQUIRE_THAT(Jb.at(5, 3), Catch::Matchers::WithinAbs(0.4, TOLERANCE));
}

TEST_CASE("Test jacobian space", "[JacobianSpace]")
{
  const std::vector<arma::vec6> Slist{
    {0, 0, 1, 0, 0.2, 0.2},
    {1, 0, 0, 2, 0, 3},
    {0, 1, 0, 0, 2, 1},
    {1, 0, 0, 0.2, 0.3, 0.4}
  };
  const std::vector<double> thetalist{0.2, 1.1, 0.1, 1.2};

  const arma::mat Js = mr::JacobianSpace(Slist, thetalist);
  // std::cout << Js << std::endl;

  REQUIRE_THAT(Js.at(0, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(Js.at(1, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(Js.at(2, 0), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(Js.at(3, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(Js.at(4, 0), Catch::Matchers::WithinAbs(0.2, TOLERANCE));
  REQUIRE_THAT(Js.at(5, 0), Catch::Matchers::WithinAbs(0.2, TOLERANCE));
  REQUIRE_THAT(Js.at(0, 1), Catch::Matchers::WithinAbs(0.98006658, TOLERANCE));
  REQUIRE_THAT(Js.at(1, 1), Catch::Matchers::WithinAbs(0.19866933, TOLERANCE));
  REQUIRE_THAT(Js.at(2, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(Js.at(3, 1), Catch::Matchers::WithinAbs(1.95218638, TOLERANCE));
  REQUIRE_THAT(Js.at(4, 1), Catch::Matchers::WithinAbs(0.43654132, TOLERANCE));
  REQUIRE_THAT(Js.at(5, 1), Catch::Matchers::WithinAbs(2.96026613, TOLERANCE));
  REQUIRE_THAT(Js.at(0, 2), Catch::Matchers::WithinAbs(-0.09011564, TOLERANCE));
  REQUIRE_THAT(Js.at(1, 2), Catch::Matchers::WithinAbs(0.4445544, TOLERANCE));
  REQUIRE_THAT(Js.at(2, 2), Catch::Matchers::WithinAbs(0.89120736, TOLERANCE));
  REQUIRE_THAT(Js.at(3, 2), Catch::Matchers::WithinAbs(-2.21635216, TOLERANCE));
  REQUIRE_THAT(Js.at(4, 2), Catch::Matchers::WithinAbs(-2.43712573, TOLERANCE));
  REQUIRE_THAT(Js.at(5, 2), Catch::Matchers::WithinAbs(3.23573065, TOLERANCE));
  REQUIRE_THAT(Js.at(0, 3), Catch::Matchers::WithinAbs(0.95749426, TOLERANCE));
  REQUIRE_THAT(Js.at(1, 3), Catch::Matchers::WithinAbs(0.28487557, TOLERANCE));
  REQUIRE_THAT(Js.at(2, 3), Catch::Matchers::WithinAbs(-0.04528405, TOLERANCE));
  REQUIRE_THAT(Js.at(3, 3), Catch::Matchers::WithinAbs(-0.51161537, TOLERANCE));
  REQUIRE_THAT(Js.at(4, 3), Catch::Matchers::WithinAbs(2.77535713, TOLERANCE));
  REQUIRE_THAT(Js.at(5, 3), Catch::Matchers::WithinAbs(2.22512443, TOLERANCE));
}
