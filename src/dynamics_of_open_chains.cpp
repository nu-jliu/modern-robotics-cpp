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

const arma::vec VelQuandraticForces(
  const arma::vec & thetalist,
  const arma::vec & dthetalist,
  const std::vector<arma::mat44> Mlist,
  const std::vector<arma::mat66> & Glist,
  const std::vector<arma::vec6> & Slist
)
{
  const arma::vec ddthetalist{thetalist.size(), arma::fill::zeros};
  const arma::vec3 g{arma::fill::zeros};
  const arma::vec6 Ftip{arma::fill::zeros};

  const arma::vec c = InverseDynamics(
    thetalist,
    dthetalist,
    ddthetalist,
    g,
    Ftip,
    Mlist,
    Glist,
    Slist
  );
  return c;
}

const arma::vec GravityForces(
  const arma::vec & thetalist,
  const arma::vec3 & g,
  const std::vector<arma::mat44> & Mlist,
  const std::vector<arma::mat66> & Glist,
  const std::vector<arma::vec6> & Slist
)
{
  const arma::vec dthetalist{thetalist.size(), arma::fill::zeros};
  const arma::vec ddthetalist{thetalist.size(), arma::fill::zeros};
  const arma::vec6 Ftip{arma::fill::zeros};

  const arma::vec grav = InverseDynamics(
    thetalist,
    dthetalist,
    ddthetalist,
    g,
    Ftip,
    Mlist,
    Glist,
    Slist
  );
  return grav;
}

const arma::vec EndEffectorForces(
  const arma::vec & thetalist,
  const arma::vec6 & Ftip,
  const std::vector<arma::mat44> & Mlist,
  const std::vector<arma::mat66> & Glist,
  const std::vector<arma::vec6> & Slist
)
{
  const arma::vec dthetalist{thetalist.size(), arma::fill::zeros};
  const arma::vec ddthetalist{thetalist.size(), arma::fill::zeros};
  const arma::vec3 g{arma::fill::zeros};

  const arma::vec JTFtip = InverseDynamics(
    thetalist,
    dthetalist,
    ddthetalist,
    g,
    Ftip,
    Mlist,
    Glist,
    Slist
  );
  return JTFtip;
}

const arma::vec ForwardDynamics(
  const arma::vec & thetalist,
  const arma::vec & dthetalist,
  const arma::vec & taulist,
  const arma::vec3 & g,
  const arma::vec6 & Ftip,
  const std::vector<arma::mat44> & Mlist,
  const std::vector<arma::mat66> & Glist,
  const std::vector<arma::vec6> & Slist
)
{
  const arma::mat Mmat = MassMatrix(thetalist, Mlist, Glist, Slist);
  const arma::vec c = VelQuandraticForces(thetalist, dthetalist, Mlist, Glist, Slist);
  const arma::vec grav = GravityForces(thetalist, g, Mlist, Glist, Slist);
  const arma::vec JTFtip = EndEffectorForces(thetalist, Ftip, Mlist, Glist, Slist);

  const arma::mat Minv = Mmat.i();
  const arma::vec rhs = taulist - c - grav - JTFtip;

  return Minv * rhs;
}

const std::tuple<const arma::vec, const arma::vec> EulerStep(
  const arma::vec & thetalist,
  const arma::vec & dthetalist,
  const arma::vec & ddthetalist,
  const double dt
)
{
  const arma::vec thetalistNext = thetalist + dthetalist * dt;
  const arma::vec dthetalistNext = dthetalist + ddthetalist * dt;

  return {thetalistNext, dthetalistNext};
}

const std::tuple<const arma::vec, const arma::vec> RK4Step(
  const arma::vec & thetalist,
  const arma::vec & dthetalist,
  const std::function<
    const std::tuple<
      const arma::vec,
      const arma::vec
    >(
      const arma::vec &,
      const arma::vec &
    )
  > & f,
  const double dt
)
{
  const auto &[k1, dk1] = f(thetalist, dthetalist);
  const auto &[k2, dk2] = f(thetalist + dt / 2. * k1, dthetalist + dt / 2. * dk1);
  const auto &[k3, dk3] = f(thetalist + dt / 2. * k2, dthetalist + dt / 2. * dk2);
  const auto &[k4, dk4] = f(thetalist + dt * k3, dthetalist + dt * dk3);

  const arma::vec thetalistNext = thetalist + dt / 6. * (k1 + 2. * k2 + 2. * k3 + k4);
  const arma::vec dthetalistNext = dthetalist + dt / 6. * (dk1 + 2. * dk2 + 2. * dk3 + dk4);

  return {thetalistNext, dthetalistNext};
}

const std::vector<arma::vec> InverseDynamicsTrajectory(
  const std::vector<arma::vec> & thetamat,
  const std::vector<arma::vec> & dthetamat,
  const std::vector<arma::vec> & ddthetamat,
  const arma::vec3 & g,
  const std::vector<arma::vec6> & Ftipmat,
  const std::vector<arma::mat44> & Mlist,
  const std::vector<arma::mat66> & Glist,
  const std::vector<arma::vec6> & Slist
)
{
  std::vector<arma::vec> taumat;

  for (size_t i = 0; i < thetamat.size(); ++i) {
    const arma::vec thetalist = thetamat.at(i);
    const arma::vec dthetalist = dthetamat.at(i);
    const arma::vec ddthetalist = ddthetamat.at(i);
    const arma::vec6 Ftip = Ftipmat.at(i);

    const arma::vec taulist = mr::InverseDynamics(
      thetalist,
      dthetalist,
      ddthetalist,
      g,
      Ftip,
      Mlist,
      Glist,
      Slist
    );
    taumat.push_back(taulist);
  }

  return taumat;
}

const std::tuple<const std::vector<arma::vec>, const std::vector<arma::vec>>
ForwardDynamicsTrajectory(
  const arma::vec & thetalist,
  const arma::vec & dthetalist,
  const std::vector<arma::vec> & taumat,
  const arma::vec3 & g,
  const std::vector<arma::vec6> & Ftipmat,
  const std::vector<arma::mat44> & Mlist,
  const std::vector<arma::mat66> & Glist,
  const std::vector<arma::vec6> & Slist,
  const double dt,
  const int intRes
)
{
  arma::vec theta{thetalist};
  arma::vec dtheta{dthetalist};

  std::vector<arma::vec> thetamat;
  std::vector<arma::vec> dthetamat;
  thetamat.push_back(theta);
  dthetamat.push_back(dtheta);

  for (size_t i = 0; i < taumat.size() - 1; ++i) {
    const arma::vec taulist = taumat.at(i);
    const arma::vec6 Ftip = Ftipmat.at(i);

    for (int j = 0; j < intRes; ++j) {

      const arma::vec ddtheta = ForwardDynamics(
        theta,
        dtheta,
        taulist,
        g,
        Ftip,
        Mlist,
        Glist,
        Slist
      );

      const auto result = EulerStep(theta, dtheta, ddtheta, dt);
      theta = std::get<0>(result);
      dtheta = std::get<1>(result);
    }

    thetamat.push_back(theta);
    dthetamat.push_back(dtheta);
  }

  return {thetamat, dthetamat};
}
}
