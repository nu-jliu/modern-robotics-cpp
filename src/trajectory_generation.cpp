#include "modern_robotics/trajectory_generation.hpp"
#include "modern_robotics/rigid_body_motions.hpp"

namespace mr
{
double CubicTimeScaling(const double Tf, const double t)
{
  const double a0 = 0.0;
  const double a1 = 0.0;
  const double a2 = 3.0 / std::pow(Tf, 2.0);
  const double a3 = -2.0 / std::pow(Tf, 3.0);

  const double st = a0 + a1 * t + a2 * std::pow(t, 2.0) + a3 * std::pow(t, 3.0);
  return st;
}

double QuinticTimeScaling(const double Tf, const double t)
{
  const double a0 = 0.0;
  const double a1 = 0.0;
  const double a2 = 0.0;
  const double a3 = 10.0 / std::pow(Tf, 3.0);
  const double a4 = -15.0 / std::pow(Tf, 4.0);
  const double a5 = 6.0 / std::pow(Tf, 5.0);

  const double st =
    a0 +
    a1 * t +
    a2 * std::pow(t, 2.0) +
    a3 * std::pow(t, 3.0) +
    a4 * std::pow(t, 4.0) +
    a5 * std::pow(t, 5.0);
  return st;
}

const std::vector<arma::vec> JointTrajectory(
  const arma::vec & thetastart,
  const arma::vec & thetaend,
  const double Tf,
  const size_t N,
  const Method & method
)
{
  std::vector<arma::vec> traj;
  const double timestep = Tf / (static_cast<double>(N) - 1.0);
  const auto func = method == Method::Cubic ?
    std::bind(&CubicTimeScaling, std::placeholders::_1, std::placeholders::_2) :
    std::bind(&QuinticTimeScaling, std::placeholders::_1, std::placeholders::_2);

  for (size_t i = 0; i < N; ++i) {
    const double t = timestep * i;
    const double s = func(Tf, t);
    const arma::vec theta = s * thetaend + (1.0 - s) * thetastart;

    traj.push_back(theta);
  }

  return traj;
}

const std::vector<arma::mat44> ScrewTrajectory(
  const arma::mat44 & Xstart,
  const arma::mat44 & Xend,
  const double Tf,
  const size_t N,
  const Method & method
)
{
  std::vector<arma::mat44> traj;
  const double timestep = Tf / (static_cast<double>(N) - 1.0);
  const auto func = method == Method::Cubic ?
    std::bind(&CubicTimeScaling, std::placeholders::_1, std::placeholders::_2) :
    std::bind(&QuinticTimeScaling, std::placeholders::_1, std::placeholders::_2);

  for (size_t i = 0; i < N; ++i) {
    const double t = timestep * i;
    const double s = func(Tf, t);

    const arma::mat44 XstartInv = TransInv(Xstart);
    const arma::mat44 se3mat = MatrixLog6(XstartInv * Xend);

    const arma::mat44 X = Xstart * MatrixExp6(se3mat * s);
    traj.push_back(X);
  }

  return traj;
}

const std::vector<arma::mat44> CartesianTrajectory(
  const arma::mat44 & Xstart,
  const arma::mat44 & Xend,
  const double Tf,
  const size_t N,
  const Method & method
)
{
  std::vector<arma::mat44> traj;
  const double timestep = Tf / (static_cast<double>(N) - 1.0);
  const auto &[Rstart, pstart] = TransToRp(Xstart);
  const auto &[Rend, pend] = TransToRp(Xend);
  const auto func = method == Method::Cubic ?
    std::bind(&CubicTimeScaling, std::placeholders::_1, std::placeholders::_2) :
    std::bind(&QuinticTimeScaling, std::placeholders::_1, std::placeholders::_2);

  for (size_t i = 0; i < N; ++i) {
    const double t = timestep * i;
    const double s = func(Tf, t);

    const arma::mat33 so3mat = MatrixLog3(Rstart.t() * Rend);
    const arma::mat33 R = Rstart * MatrixExp3(so3mat * s);
    const arma::vec3 p = s * pend + (1.0 - s) * pstart;

    const arma::mat44 X = RpToTrans(R, p);
    traj.push_back(X);
  }

  return traj;
}
}
