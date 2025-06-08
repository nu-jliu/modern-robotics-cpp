#ifndef MODERN_ROBOTICS__DYNAMICS_OF_OPEN_CHAINS___
#define MODERN_ROBOTICS__DYNAMICS_OF_OPEN_CHAINS___

#include <armadillo>

namespace mr
{
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

/// \brief Computes the mass matrix of an open chain robot based on the given configuration
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
} /// namespace mr

#endif /// MODERN_ROBOTICS__DYNAMICS_OF_OPEN_CHAINS___
