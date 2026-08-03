#include <Rcpp.h>
#include <RcppEigen.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

struct NetcropAucObservation {
  double prediction;
  std::uint8_t group;
};

// Memory-conscious rank AUC. Sorting prediction/label pairs avoids allocating
// a second R vector containing every rank and preserves average ranks for ties.
// [[Rcpp::export]]
double netcrop_auc_cpp(Rcpp::NumericVector group,
                       Rcpp::NumericVector predictions) {
  const R_xlen_t n = predictions.size();
  if (group.size() != n) {
    Rcpp::stop("group and predictions must have equal lengths");
  }
  
  std::vector<NetcropAucObservation> observations;
  observations.reserve(static_cast<std::size_t>(n));
  std::uint64_t n1 = 0;
  
  for (R_xlen_t i = 0; i < n; ++i) {
    const double prediction = predictions[i];
    const double label = group[i];
    if (!std::isfinite(prediction) || (label != 0.0 && label != 1.0)) {
      Rcpp::stop("predictions must be finite and group must contain only 0 and 1");
    }
    const std::uint8_t positive = static_cast<std::uint8_t>(label);
    observations.push_back({prediction, positive});
    n1 += positive;
  }
  
  std::sort(observations.begin(), observations.end(),
            [](const NetcropAucObservation& left,
               const NetcropAucObservation& right) {
              return left.prediction < right.prediction;
            });
  
  long double positive_rank_sum = 0.0L;
  std::size_t begin = 0;
  while (begin < observations.size()) {
    std::size_t end = begin + 1;
    std::uint64_t positives_in_tie = observations[begin].group;
    while (end < observations.size() &&
           observations[end].prediction == observations[begin].prediction) {
      positives_in_tie += observations[end].group;
      ++end;
    }
    const long double average_rank =
      (static_cast<long double>(begin + 1) + static_cast<long double>(end)) /
        2.0L;
    positive_rank_sum += average_rank * positives_in_tie;
    begin = end;
  }
  
  const std::uint64_t n0 = static_cast<std::uint64_t>(n) - n1;
  if (n0 == 0 || n1 == 0) {
    return R_NaN;
  }
  const long double u = positive_rank_sum -
    0.5L * static_cast<long double>(n1) * static_cast<long double>(n1 + 1);
  return static_cast<double>(
    u / (static_cast<long double>(n0) * static_cast<long double>(n1))
  );
}

// The original outer-add routine uses y as rows and x as columns.
// no_init_matrix() avoids zero-filling the result before values are set.
// [[Rcpp::export]]
Rcpp::NumericMatrix netcrop_outer_add_cpp(Rcpp::NumericVector x,
                                          Rcpp::NumericVector y) {
  const R_xlen_t nx = x.size();
  const R_xlen_t ny = y.size();
  Rcpp::NumericMatrix result = Rcpp::no_init_matrix(ny, nx);
  double* output = result.begin();
  
  for (R_xlen_t column = 0; column < nx; ++column) {
    const double x_value = x[column];
    const R_xlen_t offset = column * ny;
    for (R_xlen_t row = 0; row < ny; ++row) {
      output[offset + row] = x_value + y[row];
    }
  }
  return result;
}

// Translation-enabled orthogonal Procrustes without materializing an N-by-N
// centering matrix. The rotation-only path remains in R to use BLAS crossprod.
// [[Rcpp::export]]
Rcpp::List netcrop_procrustes_translated_cpp(
    const Eigen::Map<Eigen::MatrixXd>& X,
    const Eigen::Map<Eigen::MatrixXd>& Xstar) {
  if (X.rows() != Xstar.rows() || X.cols() != Xstar.cols()) {
    Rcpp::stop("X and Xstar must have equal dimensions");
  }
  
  const Eigen::Index n = X.rows();
  const Eigen::RowVectorXd mean_x = X.colwise().mean();
  const Eigen::RowVectorXd mean_xstar = Xstar.colwise().mean();
  const Eigen::MatrixXd cross = Xstar.transpose() * X -
    static_cast<double>(n) * mean_xstar.transpose() * mean_x;
  
  Eigen::MatrixXd rotation;
  if (X.cols() == 1) {
    rotation = Eigen::MatrixXd::Ones(1, 1);
  } else {
    const Eigen::JacobiSVD<Eigen::MatrixXd> decomposition(
        cross,
        Eigen::ComputeFullU | Eigen::ComputeFullV
    );
    rotation = decomposition.matrixV() * decomposition.matrixU().transpose();
  }
  
  Eigen::MatrixXd aligned = X * rotation;
  const Eigen::RowVectorXd translation = mean_xstar - mean_x * rotation;
  aligned.rowwise() += translation;
  
  return Rcpp::List::create(
    Rcpp::Named("X.new") = aligned,
    Rcpp::Named("R") = rotation,
    Rcpp::Named("t") = translation
  );
}
