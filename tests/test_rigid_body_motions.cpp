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

  const arma::mat33 R_inv = mr::RotInv(R);
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

  const arma::mat33 so3mat = mr::VecToso3(omg);
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

  const arma::vec3 omg = mr::so3ToVec(so3mat);
  // std::cout << so3mat << std::endl;

  REQUIRE_THAT(omg.at(0), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(omg.at(1), Catch::Matchers::WithinAbs(2, TOLERANCE));
  REQUIRE_THAT(omg.at(2), Catch::Matchers::WithinAbs(3, TOLERANCE));
}

TEST_CASE("Tesing axis angle", "[AxisAng3]")
{
  const arma::vec3 expc3{1, 2, 3};

  const auto &[omghat, theta] = mr::AxisAng3(expc3);
  // std::cout << omghat << std::endl;
  // std::cout << theta << std::endl;

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

TEST_CASE("Testing Matrix exponential", "[MatrixExp3]")
{
  const arma::mat33 so3mat{
    {0, -3, 2},
    {3, 0, -1},
    {-2, 1, 0}
  };

  const arma::mat33 R = mr::MatrixExp3(so3mat);
  // std::cout << R << std::endl;

  REQUIRE_THAT(R.at(0, 0), Catch::Matchers::WithinAbs(-0.69492056, TOLERANCE));
  REQUIRE_THAT(R.at(0, 1), Catch::Matchers::WithinAbs(+0.71352099, TOLERANCE));
  REQUIRE_THAT(R.at(0, 2), Catch::Matchers::WithinAbs(+0.08929286, TOLERANCE));
  REQUIRE_THAT(R.at(1, 0), Catch::Matchers::WithinAbs(-0.19200697, TOLERANCE));
  REQUIRE_THAT(R.at(1, 1), Catch::Matchers::WithinAbs(-0.30378504, TOLERANCE));
  REQUIRE_THAT(R.at(1, 2), Catch::Matchers::WithinAbs(+0.93319235, TOLERANCE));
  REQUIRE_THAT(R.at(2, 0), Catch::Matchers::WithinAbs(+0.69297817, TOLERANCE));
  REQUIRE_THAT(R.at(2, 1), Catch::Matchers::WithinAbs(+0.63134970, TOLERANCE));
  REQUIRE_THAT(R.at(2, 2), Catch::Matchers::WithinAbs(+0.34810748, TOLERANCE));
}

TEST_CASE("Testing Matrix Logarithm", "[MatrixLog3]")
{
  const arma::mat33 R{
    {0, 0, 1},
    {1, 0, 0},
    {0, 1, 0}
  };

  const arma::mat33 so3mat = mr::MatrixLog3(R);
  // std::cout << so3mat << std::endl;

  REQUIRE_THAT(so3mat.at(0, 0), Catch::Matchers::WithinAbs(0.0, TOLERANCE));
  REQUIRE_THAT(so3mat.at(0, 1), Catch::Matchers::WithinAbs(-1.20919958, TOLERANCE));
  REQUIRE_THAT(so3mat.at(0, 2), Catch::Matchers::WithinAbs(+1.20919958, TOLERANCE));
  REQUIRE_THAT(so3mat.at(1, 0), Catch::Matchers::WithinAbs(+1.20919958, TOLERANCE));
  REQUIRE_THAT(so3mat.at(1, 1), Catch::Matchers::WithinAbs(0.0, TOLERANCE));
  REQUIRE_THAT(so3mat.at(1, 2), Catch::Matchers::WithinAbs(-1.20919958, TOLERANCE));
  REQUIRE_THAT(so3mat.at(2, 0), Catch::Matchers::WithinAbs(-1.20919958, TOLERANCE));
  REQUIRE_THAT(so3mat.at(2, 1), Catch::Matchers::WithinAbs(+1.20919958, TOLERANCE));
  REQUIRE_THAT(so3mat.at(2, 2), Catch::Matchers::WithinAbs(0.0, TOLERANCE));
}

TEST_CASE("Test R p to transformation matrix", "[RpToTrans]")
{
  const arma::mat33 R{
    {1, 0, 0},
    {0, 0, -1},
    {0, 1, 0}
  };
  const arma::vec3 p{1, 2, 5};

  const arma::mat44 T = mr::RpToTrans(R, p);
  // std::cout << T << std::endl;

  REQUIRE_THAT(T.at(0, 0), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(T.at(0, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(0, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(0, 3), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(T.at(1, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(1, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(1, 2), Catch::Matchers::WithinAbs(-1, TOLERANCE));
  REQUIRE_THAT(T.at(1, 3), Catch::Matchers::WithinAbs(2, TOLERANCE));
  REQUIRE_THAT(T.at(2, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(2, 1), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(T.at(2, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(2, 3), Catch::Matchers::WithinAbs(5, TOLERANCE));
  REQUIRE_THAT(T.at(3, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(3, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(3, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(3, 3), Catch::Matchers::WithinAbs(1, TOLERANCE));
}

TEST_CASE("Test transformation matrix to rotation matrix and vector", "TransToRp")
{
  arma::mat44 T{
    {1, 0, 0, 0},
    {0, 0, -1, 0},
    {0, 1, 0, 3},
    {0, 0, 0, 1}
  };

  const auto &[R, p] = mr::TransToRp(T);
  // std::cout << R << std::endl;
  // std::cout << p << std::endl;

  REQUIRE_THAT(R.at(0, 0), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(R.at(0, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(R.at(0, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(R.at(1, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(R.at(1, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(R.at(1, 2), Catch::Matchers::WithinAbs(-1, TOLERANCE));
  REQUIRE_THAT(R.at(2, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(R.at(2, 1), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(R.at(2, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));

  REQUIRE_THAT(p.at(0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(p.at(1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(p.at(2), Catch::Matchers::WithinAbs(3, TOLERANCE));
}

TEST_CASE("Test inverse of T", "[TransInv]")
{
  const arma::mat44 T{
    {1, 0, 0, 0},
    {0, 0, -1, 0},
    {0, 1, 0, 3},
    {0, 0, 0, 1}
  };

  const arma::mat44 Tinv = mr::TransInv(T);
  // std::cout << Tinv << std::endl;

  REQUIRE_THAT(Tinv.at(0, 0), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(Tinv.at(0, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(Tinv.at(0, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(Tinv.at(0, 3), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(Tinv.at(1, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(Tinv.at(1, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(Tinv.at(1, 2), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(Tinv.at(1, 3), Catch::Matchers::WithinAbs(-3, TOLERANCE));
  REQUIRE_THAT(Tinv.at(2, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(Tinv.at(2, 1), Catch::Matchers::WithinAbs(-1, TOLERANCE));
  REQUIRE_THAT(Tinv.at(2, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(Tinv.at(2, 3), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(Tinv.at(3, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(Tinv.at(3, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(Tinv.at(3, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(Tinv.at(3, 3), Catch::Matchers::WithinAbs(1, TOLERANCE));
}

TEST_CASE("Test vector to se3", "[VecTose3]")
{
  const arma::vec6 V{1, 2, 3, 4, 5, 6};

  const arma::mat44 se3mat = mr::VecTose3(V);
  // std::cout << se3mat << std::endl;

  REQUIRE_THAT(se3mat.at(0, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(se3mat.at(0, 1), Catch::Matchers::WithinAbs(-3, TOLERANCE));
  REQUIRE_THAT(se3mat.at(0, 2), Catch::Matchers::WithinAbs(2, TOLERANCE));
  REQUIRE_THAT(se3mat.at(0, 3), Catch::Matchers::WithinAbs(4, TOLERANCE));
  REQUIRE_THAT(se3mat.at(1, 0), Catch::Matchers::WithinAbs(3, TOLERANCE));
  REQUIRE_THAT(se3mat.at(1, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(se3mat.at(1, 2), Catch::Matchers::WithinAbs(-1, TOLERANCE));
  REQUIRE_THAT(se3mat.at(1, 3), Catch::Matchers::WithinAbs(5, TOLERANCE));
  REQUIRE_THAT(se3mat.at(2, 0), Catch::Matchers::WithinAbs(-2, TOLERANCE));
  REQUIRE_THAT(se3mat.at(2, 1), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(se3mat.at(2, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(se3mat.at(2, 3), Catch::Matchers::WithinAbs(6, TOLERANCE));
  REQUIRE_THAT(se3mat.at(3, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(se3mat.at(3, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(se3mat.at(3, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(se3mat.at(3, 3), Catch::Matchers::WithinAbs(0, TOLERANCE));
}

TEST_CASE("Test se3 to vector", "[se3ToVec]")
{
  const arma::mat44 se3mat{
    {0, -3, 2, 4},
    {3, 0, -1, 5},
    {-2, 1, 0, 6},
    {0, 0, 0, 0}
  };

  const arma::vec6 V = mr::se3ToVec(se3mat);
  // std::cout << V << std::endl;

  REQUIRE_THAT(V.at(0), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(V.at(1), Catch::Matchers::WithinAbs(2, TOLERANCE));
  REQUIRE_THAT(V.at(2), Catch::Matchers::WithinAbs(3, TOLERANCE));
  REQUIRE_THAT(V.at(3), Catch::Matchers::WithinAbs(4, TOLERANCE));
  REQUIRE_THAT(V.at(4), Catch::Matchers::WithinAbs(5, TOLERANCE));
  REQUIRE_THAT(V.at(5), Catch::Matchers::WithinAbs(6, TOLERANCE));
}

TEST_CASE("Test adjoint", "[Adjoint]")
{
  const arma::mat44 T{
    {1, 0, 0, 0},
    {0, 0, -1, 0},
    {0, 1, 0, 3},
    {0, 0, 0, 1}
  };

  const arma::mat66 AdT = mr::Adjoint(T);
  // std::cout << AdT << std::endl;

  REQUIRE_THAT(AdT.at(0, 0), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(AdT.at(0, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(0, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(0, 3), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(0, 4), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(0, 5), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(1, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(1, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(1, 2), Catch::Matchers::WithinAbs(-1, TOLERANCE));
  REQUIRE_THAT(AdT.at(1, 3), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(1, 4), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(1, 5), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(2, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(2, 1), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(AdT.at(2, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(2, 3), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(2, 4), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(2, 5), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(3, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(3, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(3, 2), Catch::Matchers::WithinAbs(3, TOLERANCE));
  REQUIRE_THAT(AdT.at(3, 3), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(AdT.at(3, 4), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(3, 5), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(4, 0), Catch::Matchers::WithinAbs(3, TOLERANCE));
  REQUIRE_THAT(AdT.at(4, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(4, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(4, 3), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(4, 4), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(4, 5), Catch::Matchers::WithinAbs(-1, TOLERANCE));
  REQUIRE_THAT(AdT.at(5, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(5, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(5, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(5, 3), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(AdT.at(5, 4), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(AdT.at(5, 5), Catch::Matchers::WithinAbs(0, TOLERANCE));
}

TEST_CASE("Test screw axis", "[ScrewToAxis]")
{
  const arma::vec3 q{3, 0, 0};
  const arma::vec3 s{0, 0, 1};
  const double h = 2;

  const arma::vec6 S = mr::ScrewToAxis(q, s, h);
  // std::cout << S << std::endl;

  REQUIRE_THAT(S.at(0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(S.at(1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(S.at(2), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(S.at(3), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(S.at(4), Catch::Matchers::WithinAbs(-3, TOLERANCE));
  REQUIRE_THAT(S.at(5), Catch::Matchers::WithinAbs(2, TOLERANCE));
}

TEST_CASE("Test Axis angle", "[AxisAng6]")
{
  const arma::vec6 expc6{1, 0, 0, 1, 2, 3};

  const auto &[S, theta] = mr::AxisAng6(expc6);
  // std::cout << S << std::endl;
  // std::cout << theta << std::endl;

  REQUIRE_THAT(theta, Catch::Matchers::WithinAbs(1, TOLERANCE));

  REQUIRE_THAT(S.at(0), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(S.at(1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(S.at(2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(S.at(3), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(S.at(4), Catch::Matchers::WithinAbs(2, TOLERANCE));
  REQUIRE_THAT(S.at(5), Catch::Matchers::WithinAbs(3, TOLERANCE));
}

TEST_CASE("Test Matrix exponetial", "[MatrixExp6]")
{
  const arma::mat44 se3mat{
    {0, 0, 0, 0},
    {0, 0, -1.57079632, 2.35619449},
    {0, 1.57079632, 0, 2.35619449},
    {0, 0, 0, 0}
  };

  const arma::mat44 T = mr::MatrixExp6(se3mat);
  // std::cout << T << std::endl;

  REQUIRE_THAT(T.at(0, 0), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(T.at(0, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(0, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(0, 3), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(1, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(1, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(1, 2), Catch::Matchers::WithinAbs(-1, TOLERANCE));
  REQUIRE_THAT(T.at(1, 3), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(2, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(2, 1), Catch::Matchers::WithinAbs(1, TOLERANCE));
  REQUIRE_THAT(T.at(2, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(2, 3), Catch::Matchers::WithinAbs(3, TOLERANCE));
  REQUIRE_THAT(T.at(3, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(3, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(3, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(T.at(3, 3), Catch::Matchers::WithinAbs(1, TOLERANCE));
}

TEST_CASE("Test matrix log", "[MatrixLog6]")
{
  const arma::mat44 T{
    {1, 0, 0, 0},
    {0, 0, -1, 0},
    {0, 1, 0, 3},
    {0, 0, 0, 1}
  };

  const arma::mat44 se3mat = mr::MatrixLog6(T);
  // std::cout << se3mat << std::endl;

  REQUIRE_THAT(se3mat.at(0, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(se3mat.at(0, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(se3mat.at(0, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(se3mat.at(0, 3), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(se3mat.at(1, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(se3mat.at(1, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(se3mat.at(1, 2), Catch::Matchers::WithinAbs(-1.57079633, TOLERANCE));
  REQUIRE_THAT(se3mat.at(1, 3), Catch::Matchers::WithinAbs(2.35619449, TOLERANCE));
  REQUIRE_THAT(se3mat.at(2, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(se3mat.at(2, 1), Catch::Matchers::WithinAbs(1.57079633, TOLERANCE));
  REQUIRE_THAT(se3mat.at(2, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(se3mat.at(2, 3), Catch::Matchers::WithinAbs(2.35619449, TOLERANCE));
  REQUIRE_THAT(se3mat.at(3, 0), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(se3mat.at(3, 1), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(se3mat.at(3, 2), Catch::Matchers::WithinAbs(0, TOLERANCE));
  REQUIRE_THAT(se3mat.at(3, 3), Catch::Matchers::WithinAbs(0, TOLERANCE));
}
