#include "modern_robotics/rigid_body_motions.hpp"
#include "modern_robotics/forward_kinematics.hpp"

namespace mr
{
const arma::mat44 FKinBody(
  const arma::mat44 & M,
  const std::vector<arma::vec6> & Blist,
  const std::vector<double> & thetalist
)
{
  arma::mat44 T(M);
  const size_t n = thetalist.size();

  for (size_t i = 0; i < n; ++i) {
    const arma::vec6 B = Blist.at(i);
    const double theta = thetalist.at(i);

    const arma::mat44 se3mat = VecTose3(B * theta);
    const arma::mat44 Tij = MatrixExp6(se3mat);
    T = T * Tij;
  }

  return T;
}

const arma::mat44 FKinSpace(
  const arma::mat44 & M,
  const std::vector<arma::vec6> & Slist,
  const std::vector<double> & thetalist
)
{
  arma::mat44 T(M);
  const size_t n = thetalist.size();

  for (size_t j = 0; j < n; ++j) {
    const size_t i = n - 1 - j;
    const arma::vec6 S = Slist.at(i);
    const double theta = thetalist.at(i);

    const arma::mat44 se3mat = VecTose3(S * theta);
    const arma::mat44 Tij = MatrixExp6(se3mat);
    T = Tij * T;
  }

  return T;
}
}
