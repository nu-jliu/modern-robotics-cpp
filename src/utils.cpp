#include "modern_robotics/utils.hpp"
#include <cmath>

namespace modern_robotics
{
const arma::vec Normalize(const arma::vec & vec)
{
  const arma::vec result = vec / arma::norm(vec);
  return result;
}
} /// namespace modern_robotics
