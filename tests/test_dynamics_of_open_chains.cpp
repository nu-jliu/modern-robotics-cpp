#include <iostream>
#include <iomanip>
#include <catch2/catch_all.hpp>

#include "modern_robotics/dynamics_of_open_chains.hpp"

constexpr double TOLERANCE = 1e-6;

TEST_CASE("Test adV", "[ad]")
{
  const arma::vec6 V{1, 2, 3, 4, 5, 6};

  const arma::mat66 adV = mr::ad(V);
  // std::cout << adV << std::endl;

  REQUIRE_THAT(adV.at(0, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(adV.at(0, 1), Catch::Matchers::WithinAbs(-3, TOLERANCE));
  REQUIRE_THAT(adV.at(0, 2), Catch::Matchers::WithinAbs(2, TOLERANCE));
  REQUIRE_THAT(adV.at(0, 3), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(adV.at(0, 4), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(adV.at(0, 5), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(adV.at(1, 0), Catch::Matchers::WithinAbs(3, TOLERANCE));
  REQUIRE_THAT(adV.at(1, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(adV.at(1, 2), Catch::Matchers::WithinAbs(-1, TOLERANCE));
  REQUIRE_THAT(adV.at(1, 3), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(adV.at(1, 4), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(adV.at(1, 5), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(adV.at(2, 0), Catch::Matchers::WithinAbs(-2, TOLERANCE));
  REQUIRE_THAT(adV.at(2, 1), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(adV.at(2, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(adV.at(2, 3), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(adV.at(2, 4), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(adV.at(2, 5), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(adV.at(3, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(adV.at(3, 1), Catch::Matchers::WithinAbs(-6, TOLERANCE));
  REQUIRE_THAT(adV.at(3, 2), Catch::Matchers::WithinAbs(5, TOLERANCE));
  REQUIRE_THAT(adV.at(3, 3), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(adV.at(3, 4), Catch::Matchers::WithinAbs(-3, TOLERANCE));
  REQUIRE_THAT(adV.at(3, 5), Catch::Matchers::WithinAbs(2, TOLERANCE));
  REQUIRE_THAT(adV.at(4, 0), Catch::Matchers::WithinAbs(6, TOLERANCE));
  REQUIRE_THAT(adV.at(4, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(adV.at(4, 2), Catch::Matchers::WithinAbs(-4, TOLERANCE));
  REQUIRE_THAT(adV.at(4, 3), Catch::Matchers::WithinAbs(3, TOLERANCE));
  REQUIRE_THAT(adV.at(4, 4), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(adV.at(4, 5), Catch::Matchers::WithinAbs(-1, TOLERANCE));
  REQUIRE_THAT(adV.at(5, 0), Catch::Matchers::WithinAbs(-5, TOLERANCE));
  REQUIRE_THAT(adV.at(5, 1), Catch::Matchers::WithinAbs(4, TOLERANCE));
  REQUIRE_THAT(adV.at(5, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(adV.at(5, 3), Catch::Matchers::WithinAbs(-2, TOLERANCE));
  REQUIRE_THAT(adV.at(5, 4), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(adV.at(5, 5), Catch::Matchers::WithinAbs(0, TOLERANCE));
}

TEST_CASE("Test inverse dynamics", "[InverseDynamics]")
{
  const arma::vec thetalist{0.1, 0.1, 0.1};
  const arma::vec dthetalist{0.1, 0.2, 0.3};
  const arma::vec ddthetalist{2, 1.5, 1};
  const arma::vec3 g{0, 0, -9.8};
  const arma::vec6 Ftip{1, 1, 1, 1, 1, 1};

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
  const arma::mat44 M34 {
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

  const arma::vec taulist = mr::InverseDynamics(
    thetalist,
    dthetalist,
    ddthetalist,
    g,
    Ftip,
    Mlist,
    Glist,
    Slist
  );
  // std::cout << taulist << std::endl;

  REQUIRE(taulist.size() == 3);
  REQUIRE_THAT(taulist.at(0), Catch::Matchers::WithinAbs(74.69616155, TOLERANCE));
  REQUIRE_THAT(taulist.at(1), Catch::Matchers::WithinAbs(-33.06766016, TOLERANCE));
  REQUIRE_THAT(taulist.at(2), Catch::Matchers::WithinAbs(-3.23057314, TOLERANCE));
}

TEST_CASE("Test constructing mass matrix", "[MassMatrix]")
{
  const arma::vec thetalist{0.1, 0.1, 0.1};
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
  const arma::mat44 M34 {
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

  const arma::mat M = mr::MassMatrix(thetalist, Mlist, Glist, Slist);
  // std::cout << M << std::endl;

  REQUIRE(M.n_rows == 3);
  REQUIRE(M.n_cols == 3);
  REQUIRE_THAT(M.at(0, 0), Catch::Matchers::WithinAbs(2.25433380e+01, TOLERANCE));
  REQUIRE_THAT(M.at(0, 1), Catch::Matchers::WithinAbs(-3.07146754e-01, TOLERANCE));
  REQUIRE_THAT(M.at(0, 2), Catch::Matchers::WithinAbs(-7.18426391e-03, TOLERANCE));
  REQUIRE_THAT(M.at(1, 0), Catch::Matchers::WithinAbs(-3.07146754e-01, TOLERANCE));
  REQUIRE_THAT(M.at(1, 1), Catch::Matchers::WithinAbs(1.96850717e+00, TOLERANCE));
  REQUIRE_THAT(M.at(1, 2), Catch::Matchers::WithinAbs(4.32157368e-01, TOLERANCE));
  REQUIRE_THAT(M.at(2, 0), Catch::Matchers::WithinAbs(-7.18426391e-03, TOLERANCE));
  REQUIRE_THAT(M.at(2, 1), Catch::Matchers::WithinAbs(4.32157368e-01, TOLERANCE));
  REQUIRE_THAT(M.at(2, 2), Catch::Matchers::WithinAbs(1.91630858e-01, TOLERANCE));
}

TEST_CASE("Test crolosis", "[VelQuandraticForces]")
{
  const arma::vec thetalist{0.1, 0.1, 0.1};
  const arma::vec dthetalist{0.1, 0.2, 0.3};

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
  const arma::mat44 M34 {
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

  const arma::vec c = mr::VelQuandraticForces(thetalist, dthetalist, Mlist, Glist, Slist);
  // std::cout << c << std::endl;

  REQUIRE(c.size() == 3);
  REQUIRE_THAT(c.at(0), Catch::Matchers::WithinAbs(0.26453118, TOLERANCE));
  REQUIRE_THAT(c.at(1), Catch::Matchers::WithinAbs(-0.05505157, TOLERANCE));
  REQUIRE_THAT(c.at(2), Catch::Matchers::WithinAbs(-0.00689132, TOLERANCE));
}

TEST_CASE("Test gravitational force", "[GravityForces]")
{
  const arma::vec thetalist{0.1, 0.1, 0.1};
  const arma::vec3 g{0, 0, -9.8};
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
  const arma::mat44 M34 {
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

  const arma::vec grav = mr::GravityForces(thetalist, g, Mlist, Glist, Slist);
  // std::cout << grav << std::endl;

  REQUIRE(grav.size() == 3);
  REQUIRE_THAT(grav.at(0), Catch::Matchers::WithinAbs(28.40331262, TOLERANCE));
  REQUIRE_THAT(grav.at(1), Catch::Matchers::WithinAbs(-37.64094817, TOLERANCE));
  REQUIRE_THAT(grav.at(2), Catch::Matchers::WithinAbs(-5.4415892, TOLERANCE));
}

TEST_CASE("Test end effector force", "[EndEffectorForces]")
{
  const arma::vec thetalist{0.1, 0.1, 0.1};
  const arma::vec6 Ftip{1, 1, 1, 1, 1, 1};

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
  const arma::mat44 M34 {
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

  const arma::vec JTFtip = mr::EndEffectorForces(thetalist, Ftip, Mlist, Glist, Slist);
  // std::cout << JTFtip << std::endl;

  REQUIRE(JTFtip.size() == 3);
  REQUIRE_THAT(JTFtip.at(0), Catch::Matchers::WithinAbs(1.40954608, TOLERANCE));
  REQUIRE_THAT(JTFtip.at(1), Catch::Matchers::WithinAbs(1.85771497, TOLERANCE));
  REQUIRE_THAT(JTFtip.at(2), Catch::Matchers::WithinAbs(1.392409, TOLERANCE));
}

TEST_CASE("Test forward dynamics", "[ForwardDynamics]")
{
  const arma::vec thetalist{0.1, 0.1, 0.1};
  const arma::vec dthetalist{0.1, 0.2, 0.3};
  const arma::vec taulist{0.5, 0.6, 0.7};
  const arma::vec3 g{0, 0, -9.8};
  const arma::vec6 Ftip{1, 1, 1, 1, 1, 1};

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
  const arma::mat44 M34 {
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

  const arma::vec ddthetalist = mr::ForwardDynamics(
    thetalist,
    dthetalist,
    taulist,
    g,
    Ftip,
    Mlist,
    Glist,
    Slist
  );
  // std::cout << ddthetalist << std::endl;

  REQUIRE(ddthetalist.size() == 3);
  REQUIRE_THAT(ddthetalist.at(0), Catch::Matchers::WithinAbs(-0.97392907, TOLERANCE));
  REQUIRE_THAT(ddthetalist.at(1), Catch::Matchers::WithinAbs(25.58466784, TOLERANCE));
  REQUIRE_THAT(ddthetalist.at(2), Catch::Matchers::WithinAbs(-32.91499212, TOLERANCE));
}

TEST_CASE("Test Euler Step", "[EulerStep]")
{
  const arma::vec thetalist{0.1, 0.1, 0.1};
  const arma::vec dthetalist{0.1, 0.2, 0.3};
  const arma::vec ddthetalist{2, 1.5, 1};
  const double dt = 0.1;

  const auto &[thetalistNext, dthetalistNext] = mr::EulerStep(
    thetalist,
    dthetalist,
    ddthetalist,
    dt
  );
  // std::cout << thetalistNext << std::endl;
  // std::cout << dthetalistNext << std::endl;

  REQUIRE_THAT(thetalistNext.at(0), Catch::Matchers::WithinAbs(0.11, TOLERANCE));
  REQUIRE_THAT(thetalistNext.at(1), Catch::Matchers::WithinAbs(0.12, TOLERANCE));
  REQUIRE_THAT(thetalistNext.at(2), Catch::Matchers::WithinAbs(0.13, TOLERANCE));

  REQUIRE_THAT(dthetalistNext.at(0), Catch::Matchers::WithinAbs(0.3, TOLERANCE));
  REQUIRE_THAT(dthetalistNext.at(1), Catch::Matchers::WithinAbs(0.35, TOLERANCE));
  REQUIRE_THAT(dthetalistNext.at(2), Catch::Matchers::WithinAbs(0.4, TOLERANCE));
}

TEST_CASE("Test RK4 Step", "[RK4Step]")
{
  const arma::vec thetalist{0.1, 0.1, 0.1};
  const arma::vec dthetalist{0.1, 0.2, 0.3};
  const arma::vec ddthetalist{2, 1.5, 1};
  const arma::vec3 g{0, 0, -9.8};
  const arma::vec6 Ftip{1, 1, 1, 1, 1, 1};
  const double dt = 0.1;

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
  const arma::mat44 M34 {
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

  const arma::vec taulist = mr::InverseDynamics(
    thetalist,
    dthetalist,
    ddthetalist,
    g,
    Ftip,
    Mlist,
    Glist,
    Slist
  );
  const auto f = [](
    const arma::vec & thetalist,
    const arma::vec & dthetalist,
    const arma::vec & taulist,
    const arma::vec3 & g,
    const arma::vec6 & Ftip,
    const std::vector<arma::mat44> & Mlist,
    const std::vector<arma::mat66> & Glist,
    const std::vector<arma::vec6> & Slist
    ) -> const std::tuple<const arma::vec, const arma::vec> {
      const arma::vec ddthetalist = mr::ForwardDynamics(
        thetalist,
        dthetalist,
        taulist,
        g,
        Ftip,
        Mlist,
        Glist,
        Slist
      );

      return {dthetalist, ddthetalist};
    };

  const auto result = mr::RK4Step(
    thetalist,
    dthetalist,
    std::bind(
      f,
      std::placeholders::_1,
      std::placeholders::_2,
      taulist,
      g,
      Ftip,
      Mlist,
      Glist,
      Slist
    ),
    dt
  );
  const arma::vec thetalistNext = std::get<0>(result);
  const arma::vec dthetalistNext = std::get<1>(result);
  // std::cout << thetalistNext << std::endl;
  // std::cout << dthetalistNext << std::endl;

  REQUIRE_THAT(thetalistNext.at(0), Catch::Matchers::WithinAbs(0.11971397, TOLERANCE));
  REQUIRE_THAT(thetalistNext.at(1), Catch::Matchers::WithinAbs(0.12765397, TOLERANCE));
  REQUIRE_THAT(thetalistNext.at(2), Catch::Matchers::WithinAbs(0.13438573, TOLERANCE));

  REQUIRE_THAT(dthetalistNext.at(0), Catch::Matchers::WithinAbs(0.29073198, TOLERANCE));
  REQUIRE_THAT(dthetalistNext.at(1), Catch::Matchers::WithinAbs(0.35490282, TOLERANCE));
  REQUIRE_THAT(dthetalistNext.at(2), Catch::Matchers::WithinAbs(0.38039816, TOLERANCE));
}
