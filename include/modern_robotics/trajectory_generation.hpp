#ifndef MODERN_ROBOTICS__TRAJECTORY_GENERATION_HPP___
#define MODERN_ROBOTICS__TRAJECTORY_GENERATION_HPP___

#include <armadillo>
#include <vector>

namespace mr
{
/// \defgroup trajectory_generation Chapter 9. Trajectory Generation

/// \brief Trajectory generation method
enum class Method : uint8_t
{
  Cubic, /// cubic time scaling using 3rd-order polynomial
  Quintic /// 5th-order polynomial using 5th-order polynomial
};

/// \ingroup trajectory_generation
/// \brief Computes s(t) for a cubic time scaling
/// \param Tf Total time of the motion in seconds from rest to rest
/// \param t The current time t satisfying 0 < t < Tf
/// \return The path parameter s(t) corresponding to a third-order
///         polynomial motion that begins and ends at zero velocity
double CubicTimeScaling(const double Tf, const double t);

/// \ingroup trajectory_generation
/// \brief Computes s(t) for a quintic time scaling
/// \param Tf Total time of the motion in seconds from rest to rest
/// \param t The current time t satisfying 0 < t < Tf
/// \return The path parameter s(t) corresponding to a fifth-order
///         polynomial motion that begins and ends at zero velocity and zero
///         acceleration
double QuinticTimeScaling(const double Tf, const double t);

/// \ingroup trajectory_generation
/// \brief Computes a straight-line trajectory in joint space
/// \param thetastart The initial joint variables
/// \param thetaend The final joint variables
/// \param Tf Total time of the motion in seconds from rest to rest
/// \param N The number of points N > 1 (Start and stop) in the discrete
///          representation of the trajectory
/// \param method The time-scaling method, where 3 indicates cubic (third-
///               order polynomial) time scaling and 5 indicates quintic
///               (fifth-order polynomial) time scaling
/// \return A trajectory as an N x n matrix, where each row is an n-vector
///         of joint variables at an instant in time. The first row is
///         thetastart and the Nth row is thetaend . The elapsed time
///         between each row is Tf / (N - 1)
const std::vector<arma::vec> JointTrajectory(
  const arma::vec & thetastart,
  const arma::vec & thetaend,
  const double Tf,
  const size_t N,
  const Method & method
);

/// \ingroup trajectory_generation
/// \brief Computes a trajectory as a list of N SE(3) matrices corresponding to
///        the screw motion about a space screw axis
/// \param Xstart The initial end-effector configuration
/// \param Xend The final end-effector configuration
/// \param Tf Total time of the motion in seconds from rest to rest
/// \param N The number of points N > 1 (Start and stop) in the discrete
///           representation of the trajectory
/// \param method The time-scaling method, where 3 indicates cubic (third-
///               order polynomial) time scaling and 5 indicates quintic
///               (fifth-order polynomial) time scaling
/// \return The discretized trajectory as a list of N matrices in SE(3)
///          separated in time by Tf/(N-1). The first in the list is Xstart
///          and the Nth is Xend
const std::vector<arma::mat44> ScrewTrajectory(
  const arma::mat44 & Xstart,
  const arma::mat44 & Xend,
  const double Tf,
  const size_t N,
  const Method & method
);

/// \ingroup trajectory_generation
/// \brief Computes a trajectory as a list of N SE(3) matrices corresponding to
///        the origin of the end-effector frame following a straight line
/// \param Xstart The initial end-effector configuration
/// \param Xend The final end-effector configuration
/// \param Tf Total time of the motion in seconds from rest to rest
/// \param N The number of points N > 1 (Start and stop) in the discrete
///          representation of the trajectory
/// \param method The time-scaling method, where 3 indicates cubic (third-
///               order polynomial) time scaling and 5 indicates quintic
///               (fifth-order polynomial) time scaling
/// \return The discretized trajectory as a list of N matrices in SE(3)
///         separated in time by Tf/(N-1). The first in the list is Xstart
///         and the Nth is Xend
/// \details This function is similar to ScrewTrajectory, except the origin of the
///          end-effector frame follows a straight line, decoupled from the rotational motion.
const std::vector<arma::mat44> CartesianTrajectory(
  const arma::mat44 & Xstart,
  const arma::mat44 & Xend,
  const double Tf,
  const size_t N,
  const Method & method
);
}

#endif
