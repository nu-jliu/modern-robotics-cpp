#include <armadillo>

#include "modern_robotics/rigid_body_motions.hpp"

namespace mr
{
const arma::mat33 RotInv(const arma::mat33 & R)
{
  return R.t();
}

const arma::mat33 VecToso3(const arma::vec3 & omg)
{
  const double tx = omg.at(0);
  const double ty = omg.at(1);
  const double tz = omg.at(2);

  return {
    {0, -tz, ty},
    {tz, 0, -tx},
    {-ty, tx, 0}
  };
}

const arma::vec3 so3ToVec(const arma::mat33 & so3mat)
{
  const double tx = so3mat.at(2, 1);
  const double ty = so3mat.at(0, 2);
  const double tz = so3mat.at(1, 0);

  return {tx, ty, tz};
}

const std::tuple<const arma::vec3, double> AxisAng(const arma::vec3 & expc3)
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
    const auto &[omghat, theta] = AxisAng(omgtheta);
    const arma::mat33 omgmat = so3mat / theta;
    const arma::mat33 I33(arma::fill::eye);

    return I33 + std::sin(theta) * omgmat + (1 - std::cos(theta)) * (omgmat * omgmat);
  }
}

const arma::mat33 MatrixLog3(const arma::mat33 & R)
{
  const double acosinput = (arma::trace(R) - 1.0) / 2.0;

  if (acosinput >= 1.0) {
    return {arma::fill::zeros};
  } else if (acosinput <= -1.0) {
    const double theta = M_PI;
    arma::vec3 omg;

    if (!NearZero(1.0 + R.at(2, 2))) {
      const double prefix = 1.0 / std::sqrt(2.0 * (1.0 + R.at(2, 2)));
      omg = prefix * arma::vec3{R.at(0, 2), R.at(1, 2), 1.0 + R.at(2, 2)};
    } else if (!NearZero(1.0 + R.at(1, 1))) {
      const double prefix = 1.0 / std::sqrt(2.0 * (1.0 + R.at(1, 1)));
      omg = prefix * arma::vec3{R.at(0, 1), 1.0 + R.at(1, 1), R.at(2, 1)};
    } else {
      const double prefix = 1.0 / std::sqrt(2.0 * (1.0 + R.at(0, 0)));
      omg = prefix * arma::vec3{1.0 + R.at(0, 0), R.at(1, 0), R.at(2, 0)};
    }

    return VecToso3(theta * omg);
  } else {
    const double theta = std::acos(acosinput);
    return theta / (2.0 * std::sin(theta)) * (R - R.t());
  }
}

const arma::mat44 RpToTrans(const arma::mat33 & R, const arma::vec3 & p)
{
  const arma::mat upper = arma::join_horiz(R, p);
  const arma::rowvec4 lower{0, 0, 0, 1};

  return arma::join_vert(upper, lower);
}

const std::tuple<const arma::mat33, const arma::vec3> TransToRp(const arma::mat44 & T)
{
  const arma::mat33 R = T.submat(0, 0, 2, 2);
  const arma::vec3 p = T.submat(0, 3, 2, 3);

  return {R, p};
}

const arma::mat44 TransInv(const arma::mat44 & T)
{
  const auto &[R, p] = TransToRp(T);
  const arma::mat33 Rt = R.t();
  const arma::vec3 Rtp = Rt * p;

  const arma::mat upper = arma::join_horiz(Rt, -Rtp);
  const arma::rowvec4 lower{0, 0, 0, 1};

  return arma::join_vert(upper, lower);
}
}
