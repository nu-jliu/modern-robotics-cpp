#ifndef MODERN_ROBOTICS__DYNAMICS_OF_OPEN_CHAINS___
#define MODERN_ROBOTICS__DYNAMICS_OF_OPEN_CHAINS___

#include <armadillo>

namespace mr
{
/// \defgroup dynamics_open_open_chains Chapter 8: Dynamics of Open Chains

/// \brief Calculate the 6x6 matrix [adV] of the given 6-vector
/// \param V A 6-vector spatial velocity
/// \return The corresponding 6x6 matrix [adV]
/// \details Used to calculate the Lie bracket [V1, V2] = [adV1]V2
const arma::mat66 ad(const arma::vec6 & V);

const arma::vec InverseDynamics(
  const arma::vec & thetalist,
  const arma::vec & dthetalist,
  const arma::vec & ddthetalist,
  const arma::vec & g,
  const arma::vec6 & Ftip,
  const std::vector<arma::mat44> & Mlist,
  const std::vector<arma::mat66> & Glist,
  const std::vector<arma::vec6> & Slist
);

/// \ingroup dynamics_open_open_chains
/// \brief Computes the mass matrix of an open chain robot based on the
///        given configuration
/// \param thetalist A list of joint variables
/// \param Mlist List of link frames i relative to i-1 at the home position
/// \param Glist Spatial inertia matrices Gi of the links
/// \param Slist Screw axes Si of the joints in a space frame, in the format
///              of a matrix with axes as the columns
/// \return The numerical inertia matrix M(thetalist) of an n-joint serial
///         chain at the given configuration thetalist
/// \details This function calls InverseDynamics n times, each time passing a
///          ddthetalist vector with a single element equal to one and all other
///          inputs set to zero.
///          Each call of InverseDynamics generates a single column, and these columns
///          are assembled to create the inertia matrix.
const arma::mat MassMatrix(
  const arma::vec & thetalist,
  const std::vector<arma::mat44> & Mlist,
  const std::vector<arma::mat66> & Glist,
  const std::vector<arma::vec6> & Slist
);

/// \ingroup dynamics_open_open_chains
/// \brief Computes the Coriolis and centripetal terms in the inverse dynamics
///        of an open chain robot
/// \param thetalist A list of joint variables
/// \param dthetalist A list of joint rates
/// \param Mlist List of link frames i relative to i-1 at the home position
/// \param Glist Spatial inertia matrices Gi of the links
/// \param Slist Screw axes Si of the joints in a space frame, in the format
///              of a matrix with axes as the columns.
/// \return The vector c(thetalist,dthetalist) of Coriolis and centripetal
///         terms for a given thetalist and dthetalist.
/// \details This function calls InverseDynamics with g = 0, Ftip = 0, and
///          ddthetalist = 0.
const arma::vec VelQuandraticForces(
  const arma::vec & thetalist,
  const arma::vec & dthetalist,
  const std::vector<arma::mat44> Mlist,
  const std::vector<arma::mat66> & Glist,
  const std::vector<arma::vec6> & Slist
);

/// \ingroup dynamics_open_open_chains
/// \brief Computes the joint forces/torques an open chain robot requires to
///        overcome gravity at its configuration
/// \param thetalist A list of joint variables
/// \param g 3-vector for gravitational acceleration
/// \param Mlist List of link frames i relative to i-1 at the home position
/// \param Glist Spatial inertia matrices Gi of the links
/// \param Slist Screw axes Si of the joints in a space frame, in the format
///              of a matrix with axes as the columns
/// \return The joint forces/torques required to overcome gravity at thetalist
/// \details This function calls InverseDynamics with Ftip = 0, dthetalist = 0, and
///          ddthetalist = 0.
const arma::vec GravityForces(
  const arma::vec & thetalist,
  const arma::vec3 & g,
  const std::vector<arma::mat44> & Mlist,
  const std::vector<arma::mat66> & Glist,
  const std::vector<arma::vec6> & Slist
);

/// @brief Computes the joint forces/torques an open chain robot requires only to
///        create the end-effector force Ftip
/// @param thetalist A list of joint variables
/// @param Ftip Spatial force applied by the end-effector expressed in frame {n+1}
/// @param Mlist List of link frames i relative to i-1 at the home position
/// @param Glist Spatial inertia matrices Gi of the links
/// @param Slist Screw axes Si of the joints in a space frame, in the format
///              of a matrix with axes as the columns
/// @return The joint forces and torques required only to create the
///         end-effector force Ftip
const arma::vec EndEffectorForces(
  const arma::vec & thetalist,
  const arma::vec6 & Ftip,
  const std::vector<arma::mat44> & Mlist,
  const std::vector<arma::mat66> & Glist,
  const std::vector<arma::vec6> & Slist
);

const arma::vec ForwardDynamics(
  const arma::vec & thetalist,
  const arma::vec & dthetalist,
  const arma::vec & taulist,
  const arma::vec3 & g,
  const arma::vec6 & Ftip,
  const std::vector<arma::mat44> & Mlist,
  const std::vector<arma::mat66> & Glist,
  const std::vector<arma::vec6> & Slist
);

const std::tuple<const arma::vec, const arma::vec> EulerStep(
  const arma::vec & thetalist,
  const arma::vec & dthetalist,
  const arma::vec & ddthetalist
);

const std::vector<arma::vec> InverseDynamicsTrajectory(
  const std::vector<arma::vec> & thetamat,
  const std::vector<arma::vec> & dthetamat,
  const std::vector<arma::vec> & ddthetamat,
  const arma::vec3 & g,
  const arma::vec6 & Ftip,
  const std::vector<arma::mat44> & Mlist,
  const std::vector<arma::mat66> & Glist,
  const std::vector<arma::vec6> & Slist
);

const std::tuple<const std::vector<arma::vec>, const std::vector<arma::vec>>
ForwardDynamicsTrajectory(
  const arma::vec & thetalist,
  const arma::vec & dthetalist,
  const std::vector<arma::vec> & taumat,
  const arma::vec3 & g,
  const arma::vec6 & Ftip,
  const std::vector<arma::mat44> & Mlist,
  const std::vector<arma::mat66> & Glist,
  const std::vector<arma::vec6> & Slist,
  const double dt,
  const int intRes
);
} /// namespace mr

#endif /// MODERN_ROBOTICS__DYNAMICS_OF_OPEN_CHAINS___
