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
  for (size_t i = 0; i < thetalist.size(); ++i) {
    const arma::vec6 B = Blist.at(i);
    const double theta = thetalist.at(i);
    const arma::mat44 se3mat = VecTose3(B * theta);
    const arma::mat44 Tij = MatrixExp6(se3mat);

    T *= Tij;
  }

  return T;
}
}
