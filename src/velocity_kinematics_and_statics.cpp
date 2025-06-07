#include "modern_robotics/rigid_body_motions.hpp"
#include "modern_robotics/velocity_kinematics_and_statics.hpp"

namespace mr
{
const arma::mat JacobianBody(
  const std::vector<arma::vec6> & Blist,
  const std::vector<double> & thetalist
)
{
  arma::mat Jb{6, Blist.size(), arma::fill::zeros};
  const size_t n = thetalist.size();
  arma::mat44 T{arma::fill::eye};

  for (size_t j = 0; j < n - 1; ++j) {
    const size_t i = n - 2 - j;
    const arma::vec6 B = Blist.at(i + 1);
    const double theta = thetalist.at(i + 1);

    const arma::mat44 se3mat = VecTose3(B * -theta);
    const arma::mat44 Tij = MatrixExp6(se3mat);

    T *= Tij;
    const arma::mat66 AdT = Adjoint(T);

    Jb.col(i) = AdT * Blist.at(i);
  }

  Jb.col(n - 1) = Blist.at(n - 1);
  return Jb;
}

const arma::mat JacobianSpace(
  const std::vector<arma::vec6> & Slist,
  const std::vector<double> & thetalist
)
{
  const size_t n = thetalist.size();
  arma::mat Js{6, Slist.size(), arma::fill::zeros};
  arma::mat44 T{arma::fill::eye};

  for (size_t i = 1; i < n; ++i) {
    const arma::vec6 S = Slist.at(i - 1);
    const double theta = thetalist.at(i - 1);

    const arma::mat44 se3mat = VecTose3(S * theta);
    const arma::mat44 Tij = MatrixExp6(se3mat);

    T *= Tij;
    const arma::mat66 AdT = Adjoint(T);

    Js.col(i) = AdT * Slist.at(i);
  }

  Js.col(0) = Slist.at(0);
  return Js;
}
}
