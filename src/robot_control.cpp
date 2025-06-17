#include "modern_robotics/robot_control.hpp"
#include "modern_robotics/dynamics_of_open_chains.hpp"

namespace mr
{
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
)
{
  const size_t n = thetalist.size();
  const arma::mat I{n, n, arma::fill::eye};
  const arma::mat Kp = kp * I;
  const arma::mat Ki = ki * I;
  const arma::mat Kd = kd * I;

  const arma::vec ep = thetalistd - thetalist;
  const arma::vec ei = eint + ep;
  const arma::vec ed = dthetalistd - dthetalist;
  const arma::mat Mmat = MassMatrix(thetalist, Mlist, Glist, Slist);

  const arma::vec taulist1 = Mmat * (Kp * ep + Ki * ei + Kd * ed);
  const arma::vec taulist2 = InverseDynamics(
    thetalist,
    dthetalist,
    ddthetalistd,
    g,
    {arma::fill::zeros},
    Mlist,
    Glist,
    Slist
  );

  return taulist1 + taulist2;
}


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
  const unsigned int intRes
)
{
  const size_t m = thetalist.size();
  const size_t n = thetamatd.size();

  arma::vec thetacurrent(thetalist);
  arma::vec dthetacurrent(dthetalist);
  arma::vec eint{m, arma::fill::zeros};

  std::vector<arma::vec> taumat;
  std::vector<arma::vec> thetamat;

  for (size_t i = 0; i < n; ++i) {
    const arma::vec thetalistd = thetamatd.at(i);
    const arma::vec dthetalistd = dthetamatd.at(i);
    const arma::vec ddthetalistd = ddthetamatd.at(i);
    const arma::vec taulist = ComputeTorque(
      thetacurrent,
      dthetacurrent,
      eint,
      gtilde,
      Mtildelist,
      Gtildelist,
      Slist,
      thetalistd,
      dthetalistd,
      ddthetalistd,
      kp,
      ki,
      kd
    );

    for (size_t j = 0; j < intRes; ++j) {
      const arma::vec6 Ftip = Ftipmat.at(i);
      const arma::vec ddthetalist = ForwardDynamics(
        thetacurrent,
        dthetacurrent,
        taulist,
        g,
        Ftip,
        Mlist,
        Glist,
        Slist
      );

      const auto res = EulerStep(
        thetacurrent,
        dthetacurrent,
        ddthetalist,
        dt / static_cast<double>(intRes)
      );

      thetacurrent = std::get<0>(res);
      dthetacurrent = std::get<1>(res);
    }

    taumat.push_back(taulist);
    thetamat.push_back(thetalist);
    const arma::vec e = thetalistd - thetacurrent;
    eint += e;
  }

  return {taumat, thetamat};
}
}
