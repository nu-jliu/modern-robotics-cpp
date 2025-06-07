#ifndef MODERN_ROBOTICS__RIGID_BODY_MOTIONS_HPP___
#define MODERN_ROBOTICS__RIGID_BODY_MOTIONS_HPP___

#include <armadillo>
#include "modern_robotics/utils.hpp"

namespace mr
{
/// \defgroup rigid_body_motions Chapter 3: Rigid-Body Motions

/// \ingroup rigid_body_motions
/// \brief Inverts a rotation matrix
/// \param R rotation matrix
/// \return The inverse of R
const arma::mat33 RotInv(const arma::mat33 & R);

/// \ingroup rigid_body_motions
/// \brief Converts a 3-vector to an so(3) representation
/// \param omg A 3-vector
/// \return The skew symmetric representation of omg
const arma::mat33 VecToso3(const arma::vec3 & omg);

/// \ingroup rigid_body_motions
/// \brief Converts an so(3) representation to a 3-vector
/// \param so3mat A 3x3 skew-symmetric matrix
/// \return The 3-vector corresponding to so3mat
const arma::vec3 so3ToVec(const arma::mat33 & so3mat);

/// \ingroup rigid_body_motions
/// \brief Converts a 3-vector of exponential coordinates for rotation into
///        axis-angle form
/// \param expc3 A 3-vector of exponential coordinates for rotation
/// \return omghat: A unit rotation axis
/// \return theta: The corresponding rotation angle
const std::tuple<const arma::vec3, double> AxisAng3(const arma::vec3 & expc3);

/// \ingroup rigid_body_motions
/// \brief Computes the matrix exponential of a matrix in so(3)
/// \param so3mat A 3x3 skew-symmetric matrix
/// \return The matrix exponential of so3mat
const arma::mat33 MatrixExp3(const arma::mat33 & so3mat);

/// \ingroup rigid_body_motions
/// \brief Computes the matrix logarithm of a rotation matrix
/// \param R A 3x3 rotation matrix
/// \return The matrix logarithm of R
const arma::mat33 MatrixLog3(const arma::mat33 & R);

/// \ingroup rigid_body_motions
/// \brief Converts a rotation matrix and a position vector into homogeneous
///        transformation matrix
/// \param R A 3x3 rotation matrix
/// \param p A 3-vector
/// \return A homogeneous transformation matrix corresponding to the inputs
const arma::mat44 RpToTrans(const arma::mat33 & R, const arma::vec3 & p);

/// \ingroup rigid_body_motions
/// \brief Converts a homogeneous transformation matrix into a rotation matrix
///        and position vector
/// \param T A homogeneous transformation matrix
/// \return R: The corresponding rotation matrix,
/// \return p: The corresponding position vector.
const std::tuple<const arma::mat33, const arma::vec3> TransToRp(const arma::mat44 & T);

/// \ingroup rigid_body_motions
/// \brief Inverts a homogeneous transformation matrix
/// \param T A homogeneous transformation matrix
/// \return The inverse of T
/// \details Uses the structure of transformation matrices to avoid taking a matrix
///          inverse, for efficiency.
const arma::mat44 TransInv(const arma::mat44 & T);

/// \ingroup rigid_body_motions
/// \brief Converts a spatial velocity vector into a 4x4 matrix in se3
/// \param V A 6-vector representing a spatial velocity
/// \return The 4x4 se3 representation of V
const arma::mat44 VecTose3(const arma::vec6 & V);

/// \ingroup rigid_body_motions
/// \brief Converts an se3 matrix into a spatial velocity vector
/// \param se3mat A 4x4 matrix in se3
/// \return The spatial velocity 6-vector corresponding to se3mat
const arma::vec6 se3ToVec(const arma::mat44 & se3mat);

/// \ingroup rigid_body_motions
/// \brief Computes the adjoint representation of a homogeneous transformation matrix
/// \param T A homogeneous transformation matrix
/// \return The 6x6 adjoint representation [AdT] of T
const arma::mat66 Adjoint(const arma::mat44 & T);

/// \ingroup rigid_body_motions
/// \brief Takes a parametric description of a screw axis and converts it to a
///        normalized screw axis
/// \param q A point lying on the screw axis
/// \param s A unit vector in the direction of the screw axis
/// \param h The pitch of the screw axis
/// \return A normalized screw axis described by the inputs
const arma::vec6 ScrewToAxis(const arma::vec3 & q, const arma::vec3 & s, const double h);

/// \ingroup rigid_body_motions
/// \brief Converts a 6-vector of exponential coordinates into screw axis-angle form
/// \param expc6 A 6-vector of exponential coordinates for rigid-body motion S*theta
/// \return S: The corresponding normalized screw axis
/// \return theta: The distance traveled along/about S
const std::tuple<const arma::vec6, double> AxisAng6(const arma::vec6 & expc6);

/// \ingroup rigid_body_motions
/// \brief Computes the matrix exponential of an se3 representation of
///        exponential coordinates
/// \param se3mat A matrix in se3
/// \return The matrix exponential of se3mat
const arma::mat44 MatrixExp6(const arma::mat44 & se3mat);

/// \ingroup rigid_body_motions
/// \brief Computes the matrix logarithm of a homogeneous transformation matrix
/// \param T A matrix in SE3
/// \return The matrix logarithm of R
const arma::mat44 MatrixLog6(const arma::mat44 & T);

/// \ingroup rigid_body_motions
/// \brief Returns a projection of mat into SO(3)
/// \param mat A matrix near SO(3) to project to SO(3)
/// \return The closest matrix to R that is in SO(3)
/// \details Projects a matrix mat to the closest matrix in SO(3) using singular-value
///          decomposition (see
///          http://hades.mech.northwestern.edu/index.php/Modern_Robotics_Linear_Algebra_Review).
///          This function is only appropriate for matrices close to SO(3).
const arma::mat33 ProjectToSO3(const arma::mat33 & mat);

/// \ingroup rigid_body_motions
/// \brief Returns a projection of mat into SE(3)
/// \param mat A 4x4 matrix to project to SE(3)
/// \return  The closest matrix to T that is in SE(3)
/// \details Projects a matrix mat to the closest matrix in SE(3) using singular-value
///          decomposition (see
///          http://hades.mech.northwestern.edu/index.php/Modern_Robotics_Linear_Algebra_Review).
///          This function is only appropriate for matrices close to SE(3).
const arma::mat44 ProjectToSE3(const arma::mat44 & mat);

/// \ingroup rigid_body_motions
/// \brief Returns the Frobenius norm to describe the distance of mat from the SO(3) manifold
/// \param mat A 3x3 matrix
/// \return A quantity describing the distance of mat from the SO(3) manifold
/// \details Computes the distance from mat to the SO(3) manifold using the following
///          method:
///            If det(mat) <= 0, return a large number.
///            If det(mat) > 0, return norm(mat^T.mat - I).
double DistanceToSO3(const arma::mat33 & mat);

/// \ingroup rigid_body_motions
/// \brief Returns the Frobenius norm to describe the distance of mat from the SE(3) manifold
/// \param mat A 4x4 matrix
/// \return A quantity describing the distance of mat from the SE(3) manifold
/// \details Computes the distance from mat to the SE(3) manifold using the following
///          method:
///          Compute the determinant of matR, the top 3x3 submatrix of mat.
///          If det(matR) <= 0, return a large number.
///          If det(matR) > 0, replace the top 3x3 submatrix of mat with matR^T.matR,
///          and set the first three entries of the fourth column of mat to zero. Then
///          return norm(mat - I).
double DistanceToSE3(const arma::mat44 & mat);

/// \ingroup rigid_body_motions
/// \brief Returns true if mat is close to or on the manifold SO(3)
/// \param mat A 3x3 matrix
/// \return True if mat is very close to or in SO(3), false otherwise
/// \details Computes the distance d from mat to the SO(3) manifold using the
///          following method:
///          If det(mat) <= 0, d = a large number.
///          If det(mat) > 0, d = norm(mat^T.mat - I).
///          If d is close to zero, return true. Otherwise, return false.
inline bool TestIfSO3(const arma::mat33 & mat)
{
  return std::fabs(DistanceToSO3(mat)) < tolerance;
}

/// \ingroup rigid_body_motions
/// \brief Returns true if mat is close to or on the manifold SE(3)
/// \param mat A 4x4 matrix
/// \return True if mat is very close to or in SE(3), false otherwise
/// \details Computes the distance d from mat to the SE(3) manifold using the
///          following method:
///          Compute the determinant of the top 3x3 submatrix of mat.
///          If det(mat) <= 0, d = a large number.
///          If det(mat) > 0, replace the top 3x3 submatrix of mat with mat^T.mat, and
///          set the first three entries of the fourth column of mat to zero.
///          Then d = norm(T - I).
///          If d is close to zero, return true. Otherwise, return false.
inline bool TestIfSE3(const arma::mat44 & mat)
{
  return std::fabs(DistanceToSE3(mat)) < tolerance;
}
} /// namespace modern_robotics

#endif /// MODERN_ROBOTICS__RIGID_BODY_MOTIONS_HPP___
