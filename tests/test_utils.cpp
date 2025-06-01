#include <armadillo>
#include <catch2/catch_all.hpp>

#include "modern_robotics/utils.hpp"

constexpr double TOLERANCE = 1e-6;

TEST_CASE("Test near zero", "[NearZero]") {
  REQUIRE(modern_robotics::NearZero(0.0));
  REQUIRE(!modern_robotics::NearZero(1e-5));
  REQUIRE(!modern_robotics::NearZero(-1e-5));
  REQUIRE(!modern_robotics::NearZero(1e-3));
  REQUIRE(!modern_robotics::NearZero(-1e-3));
  REQUIRE(modern_robotics::NearZero(1e-10));
  REQUIRE(modern_robotics::NearZero(-1e-10));
  REQUIRE(modern_robotics::NearZero(1e-15));
  REQUIRE(modern_robotics::NearZero(-1e-15));
  REQUIRE(modern_robotics::NearZero(1e-20));
}

TEST_CASE("Test normalize vector", "[Normalize]") {
  const arma::vec vec = {1.0, 2.0, 3.0};
  const double norm = std::sqrt(14.0);
  const arma::vec normalized_vec = modern_robotics::Normalize(vec);

  REQUIRE_THAT(arma::norm(normalized_vec), Catch::Matchers::WithinAbs(1.0, TOLERANCE));

  REQUIRE(modern_robotics::NearZero(normalized_vec.at(0) - 1.0 / norm));
  REQUIRE(modern_robotics::NearZero(normalized_vec.at(1) - 2.0 / norm));
  REQUIRE(modern_robotics::NearZero(normalized_vec.at(2) - 3.0 / norm));
}
