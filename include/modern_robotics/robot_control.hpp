#ifndef MODERN_ROBOTICS__ROBOT_CONTROL_HPP___
#define MODERN_ROBOTICS__ROBOT_CONTROL_HPP___

#include <tuple>
#include <vector>
#include <armadillo>

namespace mr
{
/// \defgroup robot_control Chapter 11. Robot Control

/// \ingroup robot_control
/// \brief Computes the joint control torques at a particular time instant
/// \param thetalist n-vector of joint variables
/// \param dthetalist n-vector of joint rates
/// \param eint n-vector of the time-integral of joint errors
/// \param g Gravity vector g
/// \param Mlist List of link frames {i} relative to {i-1} at the home position
/// \param Glist Spatial inertia matrices Gi of the links
/// \param Slist Screw axes Si of the joints in a space frame, in the format
///              of a matrix with axes as the columns
/// \param thetalistd n-vector of reference joint variables
/// \param dthetalistd n-vector of reference joint velocities
/// \param ddthetalistd n-vector of reference joint accelerations
/// \param kp The feedback proportional gain (identical for each joint)
/// \param ki The feedback integral gain (identical for each joint)
/// \param kd The feedback derivative gain (identical for each joint)
/// \return The vector of joint forces/torques computed by the feedback
///         linearizing controller at the current instant
const arma::vec ComputeTorque(
  const arma::vec & thetalist,
  const arma::vec & dthetalist,
  const arma::vec & eint,
  const arma::vec3 & g,
  const std::vector<arma::mat44> & Mlist,
  const std::vector<arma::mat66> & Glist,
  const std::vector<arma::vec6> & Slist,
  const arma::vec & thetalistd,
  const arma::vec & dthetalistd,
  const arma::vec & ddthetalistd,
  const double kp,
  const double ki,
  const double kd
);

/// \ingroup robot_control
/// \brief Simulates the computed torque controller over a given desired trajectory
/// \param thetalist n-vector of initial joint variables
/// \param dthetalist n-vector of initial joint velocities
/// \param g Actual gravity vector g
/// \param Ftipmat An N x 6 matrix of spatial forces applied by the end-
///                effector (If there are no tip forces the user should
///                input a zero and a zero matrix will be used)
/// \param Mlist Actual list of link frames i relative to i-1 at the home position
/// \param Glist Actual spatial inertia matrices Gi of the links
/// \param Slist Screw axes Si of the joints in a space frame, in the format
///              of a matrix with axes as the columns
/// \param thetamatd An Nxn matrix of desired joint variables from the
///                  reference trajectory
/// \param dthetamatd An Nxn matrix of desired joint velocities
/// \param ddthetamatd An Nxn matrix of desired joint accelerations
/// \param gtilde The gravity vector based on the model of the actual robot
///               (actual values given above)
/// \param Mtildelist The link frame locations based on the model of the
///                   actual robot (actual values given above)
/// \param Gtildelist The link spatial inertias based on the model of the
///                   actual robot (actual values given above)
/// \param kp The feedback proportional gain (identical for each joint)
/// \param ki The feedback integral gain (identical for each joint)
/// \param kd The feedback derivative gain (identical for each joint)
/// \param dt The timestep between points on the reference trajectory
/// \param intRes Integration resolution is the number of times integration
///               (Euler) takes places between each time step. Must be an
///               integer value greater than or equal to 1
/// \return taumat: An Nxn matrix of the controllers commanded joint forces/
///                 torques, where each row of n forces/torques corresponds
///                 to a single time instant
/// \return thetamat: An Nxn matrix of actual joint angles
///                   The end of this function plots all the actual and desired joint angles
///                   using matplotlib and random libraries.
const std::tuple<std::vector<arma::vec>, std::vector<arma::vec>>
SimulateControl(
  const arma::vec & thetalist,
  const arma::vec & dthetalist,
  const arma::vec3 & g,
  const std::vector<arma::vec6> & Ftipmat,
  const std::vector<arma::mat44> & Mlist,
  const std::vector<arma::mat66> & Glist,
  const std::vector<arma::vec6> & Slist,
  const std::vector<arma::vec> & thetamatd,
  const std::vector<arma::vec> & dthetamatd,
  const std::vector<arma::vec> & ddthetamatd,
  const arma::vec3 & gtilde,
  const std::vector<arma::mat44> & Mtildelist,
  const std::vector<arma::mat66> & Gtildelist,
  const double kp,
  const double ki,
  const double kd,
  const double dt,
  const size_t intRes
);
}

#endif
