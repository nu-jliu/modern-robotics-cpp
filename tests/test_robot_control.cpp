#include <iostream>
#include <catch2/catch_all.hpp>

#include "modern_robotics/robot_control.hpp"

constexpr double TOLERANCE = 1e-6;

TEST_CASE("Testing computing torque", "[ComputeTorque]")
{
  const arma::vec thetalist{0.1, 0.1, 0.1};
  const arma::vec dthetalist{0.1, 0.2, 0.3};
  const arma::vec eint{0.2, 0.2, 0.2};
  const arma::vec g{0, 0, -9.8};

  const arma::mat44 M01{
    {1, 0, 0, 0},
    {0, 1, 0, 0},
    {0, 0, 1, 0.089159},
    {0, 0, 0, 1}
  };
  const arma::mat44 M12{
    {0, 0, 1, 0.28},
    {0, 1, 0, 0.13585},
    {-1, 0, 0, 0},
    {0, 0, 0, 1}
  };
  const arma::mat44 M23{
    {1, 0, 0, 0},
    {0, 1, 0, -0.1197},
    {0, 0, 1, 0.395},
    {0, 0, 0, 1}
  };
  const arma::mat44 M34{
    {1, 0, 0, 0},
    {0, 1, 0, 0},
    {0, 0, 1, 0.14225},
    {0, 0, 0, 1}
  };
  const arma::mat66 G1 = arma::diagmat(arma::vec6{0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7});
  const arma::mat66 G2 = arma::diagmat(
    arma::vec6{0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393}
  );
  const arma::mat66 G3 = arma::diagmat(
    arma::vec6{0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275}
  );

  const std::vector<arma::mat44> Mlist{M01, M12, M23, M34};
  const std::vector<arma::mat66> Glist{G1, G2, G3};
  const std::vector<arma::vec6> Slist{
    {1, 0, 1, 0, 1, 0},
    {0, 1, 0, -0.089, 0, 0},
    {0, 1, 0, -0.089, 0, 0.425}
  };

  const arma::vec thetalistd{1.0, 1.0, 1.0};
  const arma::vec dthetalistd{2, 1.2, 2};
  const arma::vec ddthetalistd{0.1, 0.1, 0.1};

  const double kp = 1.3;
  const double ki = 1.2;
  const double kd = 1.1;

  const arma::vec taulist = mr::ComputeTorque(
    thetalist,
    dthetalist,
    eint,
    g,
    Mlist,
    Glist,
    Slist,
    thetalistd,
    dthetalistd,
    ddthetalistd,
    kp,
    ki,
    kd
  );
  // std::cout << taulist << std::endl;

  REQUIRE_THAT(taulist.at(0), Catch::Matchers::WithinAbs(133.00525246, TOLERANCE));
  REQUIRE_THAT(taulist.at(1), Catch::Matchers::WithinAbs(-29.94223324, TOLERANCE));
  REQUIRE_THAT(taulist.at(2), Catch::Matchers::WithinAbs(-3.03276856, TOLERANCE));
}
