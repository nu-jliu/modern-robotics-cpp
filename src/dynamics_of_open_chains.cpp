#include <armadillo>

#include "modern_robotics/rigid_body_motions.hpp"
#include "modern_robotics/dynamics_of_open_chains.hpp"

namespace mr
{
const arma::mat66 ad(const arma::vec6 & V)
{
  const arma::mat33 omgmat = VecToso3(V.subvec(0, 2));
  const arma::mat33 vmat = VecToso3(V.subvec(3, 5));
  const arma::mat33 Z33{arma::fill::zeros};

  const arma::mat upper = arma::join_horiz(omgmat, Z33);
  const arma::mat lower = arma::join_horiz(vmat, omgmat);

  const arma::mat66 adV = arma::join_vert(upper, lower);
  return adV;
}

const arma::vec InverseDynamics(
  const arma::vec & thetalist,
  const arma::vec & dthetalist,
  const arma::vec & ddthetalist,
  const arma::vec & g,
  const arma::vec6 & Ftip,
  const std::vector<arma::mat44> & Mlist,
  const std::vector<arma::mat66> & Glist,
  const std::vector<arma::vec6> & Slist
)
{
  const size_t n = thetalist.size();
  arma::mat44 Mi{arma::fill::eye};
  std::vector<arma::vec6> Ai;
  std::vector<arma::mat66> AdTi;
  std::vector<arma::vec6> Vi;
  std::vector<arma::vec6> Vdi;
  arma::vec6 Fi{Ftip};


  Vi.push_back({arma::fill::zeros});
  Vdi.push_back(arma::join_cols(arma::vec3(arma::fill::zeros), -g));

  for (size_t i = 0; i < n; ++i) {
    const arma::mat44 M = Mlist.at(i);
    const arma::vec6 S = Slist.at(i);
    const double theta = thetalist.at(i);
    const double dtheta = dthetalist.at(i);
    const double ddtheta = ddthetalist.at(i);
    const arma::vec6 V_curr = Vi.at(i);
    const arma::vec6 Vd_curr = Vdi.at(i);

    Mi = Mi * M;
    const arma::mat Minv = TransInv(M);
    const arma::vec6 A = Adjoint(TransInv(Mi)) * S;
    const arma::mat66 AdT = Adjoint(MatrixExp6(VecTose3(A * -theta)) * Minv);
    const arma::vec6 V = AdT * V_curr + A * dtheta;
    const arma::vec6 Vd = AdT * Vd_curr + A * ddtheta + ad(V) * A * dtheta;

    Ai.push_back(A);
    AdTi.push_back(AdT);
    Vi.push_back(V);
    Vdi.push_back(Vd);
  }

  AdTi.push_back(Adjoint(TransInv(Mlist.at(n))));

  arma::vec taulist{n, arma::fill::zeros};
  for (size_t j = 0; j < n; ++j) {
    const size_t i = n - 1 - j;
    const arma::mat66 G = Glist.at(i);
    const arma::vec6 A = Ai.at(i);
    const arma::mat66 AdT = AdTi.at(i + 1);
    const arma::vec6 V = Vi.at(i + 1);
    const arma::vec6 Vd = Vdi.at(i + 1);

    Fi = AdT.t() * Fi + G * Vd - ad(V).t() * (G * V);
    const double tau = arma::dot(Fi, A);

    taulist.at(i) = tau;
  }

  return taulist;
}

const arma::mat MassMatrix(
  const arma::vec & thetalist,
  const std::vector<arma::mat44> & Mlist,
  const std::vector<arma::mat66> & Glist,
  const std::vector<arma::vec6> & Slist
)
{
  const size_t n = thetalist.size();
  arma::mat M{n, n, arma::fill::zeros};

  for (size_t i = 0; i < n; ++i) {
    const arma::vec dthetalist{n, arma::fill::zeros};
    arma::vec ddthetalist{n, arma::fill::zeros};
    ddthetalist.at(i) = 1;
    const arma::vec3 g{arma::fill::zeros};
    const arma::vec6 Ftip{arma::fill::zeros};

    const arma::vec taulist = InverseDynamics(
      thetalist,
      dthetalist,
      ddthetalist,
      g,
      Ftip,
      Mlist,
      Glist,
      Slist
    );

    M.col(i) = taulist.as_col();
  }

  return M;
}
}
