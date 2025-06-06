#include <limits>
#include <armadillo>

#include "modern_robotics/rigid_body_motions.hpp"

namespace mr
{
const arma::mat33 RotInv(const arma::mat33 & R)
{
  const arma::mat33 invR = R.t();
  return invR;
}

const arma::mat33 VecToso3(const arma::vec3 & omg)
{
  const double tx = omg.at(0);
  const double ty = omg.at(1);
  const double tz = omg.at(2);

  const arma::mat so3mat{
    {0, -tz, ty},
    {tz, 0, -tx},
    {-ty, tx, 0}
  };
  return so3mat;
}

const arma::vec3 so3ToVec(const arma::mat33 & so3mat)
{
  const double tx = so3mat.at(2, 1);
  const double ty = so3mat.at(0, 2);
  const double tz = so3mat.at(1, 0);

  const arma::vec3 omg{tx, ty, tz};
  return omg;
}

const std::tuple<const arma::vec3, double> AxisAng3(const arma::vec3 & expc3)
{
  const arma::vec3 omghat = Normalize(expc3);
  const double theta = arma::norm(expc3, 2);

  return {omghat, theta};
}

const arma::mat33 MatrixExp3(const arma::mat33 & so3mat)
{
  const arma::vec3 omgtheta = so3ToVec(so3mat);
  if (NearZero(arma::norm(omgtheta, 2))) {
    return {arma::fill::eye};
  } else {
    const auto &[omghat, theta] = AxisAng3(omgtheta);
    const arma::mat33 omgmat = so3mat / theta;
    const arma::mat33 I33(arma::fill::eye);

    const arma::mat33 R = I33 + std::sin(theta) * omgmat +
      (1 - std::cos(theta)) * (omgmat * omgmat);
    return R;
  }
}

const arma::mat33 MatrixLog3(const arma::mat33 & R)
{
  const double acosinput = (arma::trace(R) - 1.0) / 2.0;

  if (acosinput >= 1.0) {
    return {arma::fill::zeros};
  } else if (acosinput <= -1.0) {
    const double theta = M_PI;
    arma::vec3 omghat;

    if (!NearZero(1.0 + R.at(2, 2))) {
      const double prefix = 1.0 / std::sqrt(2.0 * (1.0 + R.at(2, 2)));
      omghat = prefix * arma::vec3{R.at(0, 2), R.at(1, 2), 1.0 + R.at(2, 2)};
    } else if (!NearZero(1.0 + R.at(1, 1))) {
      const double prefix = 1.0 / std::sqrt(2.0 * (1.0 + R.at(1, 1)));
      omghat = prefix * arma::vec3{R.at(0, 1), 1.0 + R.at(1, 1), R.at(2, 1)};
    } else {
      const double prefix = 1.0 / std::sqrt(2.0 * (1.0 + R.at(0, 0)));
      omghat = prefix * arma::vec3{1.0 + R.at(0, 0), R.at(1, 0), R.at(2, 0)};
    }

    const arma::mat33 so3mat = VecToso3(theta * omghat);
    return so3mat;
  } else {
    const double theta = std::acos(acosinput);
    const arma::mat33 so3mat = theta / (2.0 * std::sin(theta)) * (R - R.t());
    return so3mat;
  }
}

const arma::mat44 RpToTrans(const arma::mat33 & R, const arma::vec3 & p)
{
  const arma::mat upper = arma::join_horiz(R, p);
  const arma::rowvec4 lower{0, 0, 0, 1};

  const arma::mat44 T = arma::join_vert(upper, lower);
  return T;
}

const std::tuple<const arma::mat33, const arma::vec3> TransToRp(const arma::mat44 & T)
{
  const arma::mat33 R = T.submat(0, 0, 2, 2);
  const arma::vec3 p = T.col(3).subvec(0, 2);

  return {R, p};
}

const arma::mat44 TransInv(const arma::mat44 & T)
{
  const auto &[R, p] = TransToRp(T);
  const arma::mat33 Rt = R.t();
  const arma::vec3 Rtp = Rt * p;

  const arma::mat upper = arma::join_horiz(Rt, -Rtp);
  const arma::rowvec4 lower{0, 0, 0, 1};

  const arma::mat44 invT = arma::join_vert(upper, lower);
  return invT;
}

const arma::mat44 VecTose3(const arma::vec6 & V)
{
  const arma::vec3 omg = V.subvec(0, 2);
  const arma::mat33 so3mat = VecToso3(omg);
  const arma::vec3 v = V.subvec(3, 5);

  const arma::mat upper = arma::join_horiz(so3mat, v);
  const arma::rowvec4 lower(arma::fill::zeros);

  const arma::mat44 se3mat = arma::join_vert(upper, lower);
  return se3mat;
}

const arma::vec6 se3ToVec(const arma::mat44 & se3mat)
{
  const arma::vec3 omg{se3mat.at(2, 1), se3mat.at(0, 2), se3mat.at(1, 0)};
  const arma::vec3 v{se3mat.at(0, 3), se3mat.at(1, 3), se3mat.at(2, 3)};

  const arma::vec6 V = arma::join_cols(omg, v);
  return V;
}

const arma::mat66 Adjoint(const arma::mat44 & T)
{
  const auto &[R, p] = TransToRp(T);
  const arma::mat33 pmat = VecToso3(p);
  const arma::mat33 pR = pmat * R;
  const arma::mat33 Z33(arma::fill::zeros);

  const arma::mat upper = arma::join_horiz(R, Z33);
  const arma::mat lower = arma::join_horiz(pR, R);

  const arma::mat66 AdT = arma::join_vert(upper, lower);
  return AdT;
}

const arma::vec6 ScrewToAxis(const arma::vec3 & q, const arma::vec3 & s, const double h)
{
  const arma::vec3 omg(s);
  const arma::vec3 v = arma::cross(q, s) + h * s;

  const arma::vec6 S = arma::join_cols(omg, v);
  return S;
}

const std::tuple<const arma::vec6, double> AxisAng6(const arma::vec6 & expc6)
{
  const arma::vec3 omg = expc6.subvec(0, 2);
  const arma::vec3 v = expc6.subvec(3, 5);

  double theta;
  theta = arma::norm(omg);

  if (NearZero(theta)) {
    theta = arma::norm(v);
  }

  const arma::vec6 S = expc6 / theta;
  return {S, theta};
}

const arma::mat44 MatrixExp6(const arma::mat44 & se3mat)
{
  const arma::mat33 so3mat = se3mat.submat(0, 0, 2, 2);
  const arma::vec3 vtheta = se3mat.col(3).subvec(0, 2);
  const arma::vec3 omgtheta = so3ToVec(so3mat);
  const double theta = std::get<1>(AxisAng3(omgtheta));
  const arma::mat33 I33(arma::fill::eye);

  // std::cout << "Test" << std::endl;
  arma::mat33 R;
  arma::vec3 p;

  if (NearZero(theta)) {
    R = I33;
    p = vtheta;
  } else {
    const arma::mat33 omgmat = so3mat / theta;
    const arma::vec3 v = vtheta / theta;

    R = MatrixExp3(so3mat);
    p = (I33 * theta + (1.0 - std::cos(theta)) * omgmat +
      (theta - std::sin(theta)) * (omgmat * omgmat)) * v;
  }
  const arma::mat upper = arma::join_horiz(R, p);
  const arma::rowvec4 lower{0, 0, 0, 1};

  const arma::mat44 T = arma::join_vert(upper, lower);
  return T;
}

const arma::mat44 MatrixLog6(const arma::mat44 & T)
{
  const auto &[R, p] = TransToRp(T);
  const arma::mat33 I33{arma::fill::eye};
  const arma::mat33 Z33{arma::fill::zeros};

  arma::mat33 omgmat;
  arma::vec3 v;
  double theta;

  if (NearZero(arma::max(arma::max(arma::abs(R - I33))))) {
    omgmat = Z33;
    v = p / arma::norm(p);
    theta = arma::norm(p);
  } else {
    theta = std::acos((arma::trace(R) - 1.0) / 2.0);
    omgmat = MatrixLog3(R) / theta;

    const arma::mat G = 1.0 / theta * I33 - omgmat / 2.0 +
      (1.0 / theta - 1.0 / std::tan(theta / 2.0) / 2.0) * (omgmat * omgmat);
    v = G * p;
  }

  const arma::mat upper = arma::join_horiz(omgmat * theta, v * theta);
  const arma::rowvec4 lower{arma::fill::zeros};

  const arma::mat44 se3mat = arma::join_vert(upper, lower);
  return se3mat;
}

const arma::mat33 ProjectToSO3(const arma::mat33 & mat)
{
  arma::mat33 U;
  arma::vec3 s;
  arma::mat33 Vh;

  arma::svd(U, s, Vh, mat);
  arma::mat33 R = U * Vh;

  if (arma::det(R) < 0) {
    R.col(2) *= -1;
  }

  return R;
}

const arma::mat44 ProjectToSE3(const arma::mat44 & mat)
{
  const arma::mat33 R = ProjectToSO3(mat.submat(0, 0, 2, 2));
  const arma::vec3 p = mat.col(3).subvec(0, 2);

  const arma::mat44 T = RpToTrans(R, p);
  return T;
}

double DistanceToSO3(const arma::mat33 & mat)
{
  if (arma::det(mat) > 0) {
    const arma::mat33 I33{arma::fill::eye};
    return arma::norm(mat.t() * mat - I33);
  }

  return std::numeric_limits<double>::max();
}

double DistanceToSE3(const arma::mat44 & mat)
{
  const arma::mat33 Rmat = mat.submat(0, 0, 2, 2);
  if (arma::det(Rmat) > 0) {
    const arma::mat44 I44{arma::fill::eye};
    const arma::mat upper = arma::join_horiz(Rmat.t() * Rmat, arma::vec3{arma::fill::zeros});
    const arma::rowvec4 lower = mat.row(3);

    return arma::norm(arma::join_vert(upper, lower) - I44);
  }

  return std::numeric_limits<double>::max();
}
}
