#ifndef MODERN_ROBOTICS__RIGID_BODY_MOTIONS_HPP___
#define MODERN_ROBOTICS__RIGID_BODY_MOTIONS_HPP___

#include <armadillo>
#include "modern_robotics/utils.hpp"

namespace modern_robotics
{
const arma::mat33 RotInv(const arma::mat33 & R);
const arma::mat33 VecToso3(const arma::vec3 & omg);
const arma::vec3 so3ToVec(const arma::mat33 & so3mat);
const std::tuple<const arma::vec3, double> AxisAng(const arma::vec3 & expc3);
const arma::mat33 MatrixExp3(const arma::mat33 & so3mat);
const arma::mat33 MatrixLog3(const arma::mat33 & R);
const arma::mat44 RpToTrans(const arma::mat33 & R, const arma::vec3 & p);
const std::tuple<const arma::mat33, const arma::vec3> TransToRp(const arma::mat44 & T);
const arma::mat33 TransInv(const arma::mat44 & T);
const arma::mat33 VecTose3(const arma::vec3 & V);
const arma::vec3 se3ToVec(const arma::mat33 & se3mat);
const arma::mat66 Adjoint(const arma::mat44 & T);
const arma::vec3 ScrewToAxis(const arma::vec3 & v, const arma::vec3 & p, const double & h);
const std::tuple<const arma::vec6, double> AxisAng6(const arma::vec6 & expc6);
const arma::mat44 MatrixExp6(const arma::mat44 & se3mat);
const arma::mat44 MatrixLog6(const arma::mat44 & T);
} /// namespace modern_robotics

#endif /// MODERN_ROBOTICS__RIGID_BODY_MOTIONS_HPP___
