#ifndef MODERN_ROBOTICS__RIGID_BODY_MOTIONS_HPP___
#define MODERN_ROBOTICS__RIGID_BODY_MOTIONS_HPP___

#include <armadillo>
#include "modern_robotics/utils.hpp"

namespace mr
{
/// \brief Inverts a rotation matrix
/// \param R rotation matrix
/// \return The inverse of R
const arma::mat33 RotInv(const arma::mat33 & R);

/// \brief Converts a 3-vector to an so(3) representation
/// \param omg A 3-vector
/// \return The skew symmetric representation of omg
const arma::mat33 VecToso3(const arma::vec3 & omg);

/// \brief Converts an so(3) representation to a 3-vector
/// \param so3mat A 3x3 skew-symmetric matrix
/// \return The 3-vector corresponding to so3mat
const arma::vec3 so3ToVec(const arma::mat33 & so3mat);

/// \brief Converts a 3-vector of exponential coordinates for rotation into
///        axis-angle form
/// \param expc3 A 3-vector of exponential coordinates for rotation
/// \return omghat: A unit rotation axis
/// \return theta: The corresponding rotation angle
const std::tuple<const arma::vec3, double> AxisAng(const arma::vec3 & expc3);

/// \brief Computes the matrix exponential of a matrix in so(3)
/// \param so3mat A 3x3 skew-symmetric matrix
/// \return The matrix exponential of so3mat
const arma::mat33 MatrixExp3(const arma::mat33 & so3mat);

/// \brief Computes the matrix logarithm of a rotation matrix
/// \param R A 3x3 rotation matrix
/// \return The matrix logarithm of R
const arma::mat33 MatrixLog3(const arma::mat33 & R);

/// \brief Converts a rotation matrix and a position vector into homogeneous
///        transformation matrix
/// \param R A 3x3 rotation matrix
/// \param p A 3-vector
/// \return A homogeneous transformation matrix corresponding to the inputs
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
