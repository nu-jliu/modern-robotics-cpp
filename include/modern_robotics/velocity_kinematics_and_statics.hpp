#ifndef MODERN_ROBOTICS__VELOCITY_KINEMATICS_AND_STATICS_HPP___
#define MODERN_ROBOTICS__VELOCITY_KINEMATICS_AND_STATICS_HPP___

#include <armadillo>

namespace mr
{
/// \defgroup velocity_kinematics_and_statics Chapter 5: Velocity Kinematics and Statics

/// \ingroup velocity_kinematics_and_statics
/// \brief Computes the body Jacobian for an open chain robot
/// \param Blist The joint screw axes in the end-effector frame when the
///              manipulator is at the home position, in the format of a
///              matrix with axes as the columns
/// \param thetalist A list of joint coordinates
/// \return The body Jacobian corresponding to the inputs (6xn real numbers)
const arma::mat JacobianBody(
  const std::vector<arma::vec6> & Blist,
  const arma::vec & thetalist
);

/// \ingroup velocity_kinematics_and_statics
/// \brief Computes the space Jacobian for an open chain robot
/// \param Slist The joint screw axes in the space frame when the
///              manipulator is at the home position, in the format of a
///              matrix with axes as the columns
/// \param thetalist A list of joint coordinates
/// \return The space Jacobian corresponding to the inputs (6xn real numbers)
const arma::mat JacobianSpace(
  const std::vector<arma::vec6> & Slist,
  const arma::vec & thetalist
);
}

#endif
