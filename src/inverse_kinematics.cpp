#include <armadillo>

#include "modern_robotics/rigid_body_motions.hpp"
#include "modern_robotics/forward_kinematics.hpp"
#include "modern_robotics/velocity_kinematics_and_statics.hpp"
#include "modern_robotics/inverse_kinematics.hpp"

namespace mr
{
const std::pair<const arma::vec, bool> IKinBody(
  const std::vector<arma::vec6> & Blist,
  const arma::mat44 & M,
  const arma::mat44 & T,
  const arma::vec & thetalist0,
  const double emog,
  const double ev
)
{
  constexpr int max_iter = 20;
  arma::vec thetalist{thetalist0};
  int i = 0;
  bool err = true;

  do {
    const arma::mat44 Tsb = FKinBody(M, Blist, thetalist);
    const arma::mat44 Tsbinv = TransInv(Tsb);
    const arma::vec6 Vb = se3ToVec(MatrixLog6(Tsbinv * T));

    const arma::mat Jb = JacobianBody(Blist, thetalist);
    const arma::vec dtheta = arma::pinv(Jb) * Vb;
    thetalist += dtheta;

    err = arma::norm(Vb.subvec(0, 2)) > emog || arma::norm(Vb.subvec(3, 5)) > ev;
    ++i;
  } while (err && i < max_iter);

  return {thetalist, !err};
}

const std::pair<const arma::vec, bool> IKinSpace(
  const std::vector<arma::vec6> & Slist,
  const arma::mat44 & M,
  const arma::mat44 & T,
  const arma::vec & thetalist0,
  const double emog,
  const double ev
)
{
  constexpr int max_iter = 20;
  arma::vec thetalist{thetalist0};
  int i = 0;
  bool err = true;

  do {
    const arma::mat44 Tsb = FKinSpace(M, Slist, thetalist);
    const arma::mat44 Tsbinv = TransInv(Tsb);
    const arma::mat66 AdTsb = Adjoint(Tsb);
    const arma::vec6 Vs = AdTsb * se3ToVec(MatrixLog6(Tsbinv * T));

    const arma::mat Js = JacobianSpace(Slist, thetalist);
    const arma::vec dtheta = arma::pinv(Js) * Vs;
    thetalist += dtheta;

    err = arma::norm(Vs.subvec(0, 2)) > emog || arma::norm(Vs.subvec(3, 5)) > ev;
    ++i;
  } while (err && i < max_iter);

  return {thetalist, !err};
}
}
