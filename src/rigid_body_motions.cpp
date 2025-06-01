#include <armadillo>

#include "modern_robotics/rigid_body_motions.hpp"

namespace modern_robotics
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
}
