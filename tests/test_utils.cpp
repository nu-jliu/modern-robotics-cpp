#include <armadillo>
#include <catch2/catch_all.hpp>

#include "modern_robotics/utils.hpp"

constexpr double TOLERANCE = 1e-6;

TEST_CASE("Test near zero", "[NearZero]") {
  REQUIRE(mr::NearZero(0.0));
  REQUIRE(!mr::NearZero(1e-5));
  REQUIRE(!mr::NearZero(-1e-5));
  REQUIRE(!mr::NearZero(1e-3));
  REQUIRE(!mr::NearZero(-1e-3));
  REQUIRE(mr::NearZero(1e-10));
  REQUIRE(mr::NearZero(-1e-10));
  REQUIRE(mr::NearZero(1e-15));
  REQUIRE(mr::NearZero(-1e-15));
  REQUIRE(mr::NearZero(1e-20));
}

TEST_CASE("Test normalize vector", "[Normalize]") {
  const arma::vec vec{1.0, 2.0, 3.0};
  const double norm = std::sqrt(14.0);
  const arma::vec normalized_vec = mr::Normalize(vec);

  REQUIRE_THAT(arma::norm(normalized_vec), Catch::Matchers::WithinAbs(1.0, TOLERANCE));

  REQUIRE_THAT(normalized_vec.at(0) - 1.0 / norm, Catch::Matchers::WithinAbs(0.0, TOLERANCE));
  REQUIRE_THAT(normalized_vec.at(1) - 2.0 / norm, Catch::Matchers::WithinAbs(0.0, TOLERANCE));
  REQUIRE_THAT(normalized_vec.at(2) - 3.0 / norm, Catch::Matchers::WithinAbs(0.0, TOLERANCE));
}
