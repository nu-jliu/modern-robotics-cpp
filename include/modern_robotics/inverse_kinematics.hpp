#ifndef MODERN_ROBOTICS__INVERSE_KINEMATICS_HPP___
#define MODERN_ROBOTICS__INVERSE_KINEMATICS_HPP___

#include <armadillo>
#include <vector>
#include <utility>

namespace mr
{
/// \defgroup inverse_kinematics Chapter 6: Inverse Kinematics

/// \ingroup inverse_kinematics
/// \brief Computes inverse kinematics in the body frame for an open chain robot
/// \param Blist The joint screw axes in the end-effector frame when the
///              manipulator is at the home position, in the format of a
///              matrix with axes as the columns
/// \param M The home configuration of the end-effector
/// \param T The desired end-effector configuration Tsd satisfying Tsd
/// \param thetalist0 An initial guess of joint angles that are close to
/// \param emog A small positive tolerance on the end-effector orientation
///                 error. The returned joint angles must give an end-effector
///                 orientation error less than eomg
/// \param ev A small positive tolerance on the end-effector linear position
///           error. The returned joint angles must give an end-effector
///           position error less than ev
/// \return thetalist: Joint angles that achieve T within the specified tolerances
/// \return success: A logical value where TRUE means that the function found
///                  a solution and FALSE means that it ran through the set
///                  number of maximum iterations without finding a solution
///                  within the tolerances eomg and ev.
/// \details Uses an iterative Newton-Raphson root-finding method.
///          The maximum number of iterations before the algorithm is terminated has
///          been hardcoded in as a variable called maxiterations. It is set to 20 at
///          the start of the function, but can be changed if needed.
const std::pair<const arma::vec, bool> IKinBody(
  const std::vector<arma::vec6> & Blist,
  const arma::mat44 & M,
  const arma::mat44 & T,
  const arma::vec & thetalist0,
  const double emog,
  const double ev
);

const std::pair<const arma::vec, bool> IKinSpace(
  const std::vector<arma::vec6> & Slist,
  const arma::mat44 & M,
  const arma::mat44 & T,
  const arma::vec & thetalist0,
  const double emog,
  const double ev
);
}

#endif
