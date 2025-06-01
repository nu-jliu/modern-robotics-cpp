#ifndef MODERN_ROBOTICS__UTILS_HPP___
#define MODERN_ROBOTICS__UTILS_HPP___

#include <armadillo>
#include <cmath>
#include <cstdlib>

namespace modern_robotics
{
constexpr double tolerance = 1e-6;

/// \brief Check if a number is close to zero
/// \param x The number to check
/// \return True if the number is close to zero, false otherwise
constexpr bool NearZero(const double & x)
{
  return std::fabs(x) < tolerance;
}

/// \brief Normalize a vector to unit length
/// \param vec The vector to normalize
/// \return The normalized vector
const arma::vec Normalize(const arma::vec & vec);
} // namespace modern_robotics

#endif /// MODERN_ROBOTICS__UTILS_HPP___
