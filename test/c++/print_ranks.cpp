#include <nda/nda.hpp>
#include <cppdlr/cppdlr.hpp>

using namespace cppdlr;
using namespace nda;

/**
* @brief Print ranks for different combinations of DLR cutoff lambda and
* tolerance epsilon
*/
int main() {

  auto lambda    = nda::vector<double>({1e1, 1e2, 1e3, 1e4, 1e5, 1e6});
  double eps     = 1e-10;   // DLR tolerance
  auto statistic = Fermion; // Green's function statistics

  for (int i = 0; i < lambda.size(); ++i) {

    // Get DLR frequencies
    auto dlr_rf     = build_dlr_rf(lambda(i), eps);
    auto dlr_rf_sym = build_dlr_rf(lambda(i), eps, statistic, SYM);

    // Print results
    std::cout << "lambda = " << lambda(i) << std::endl;
    std::cout << "Unsymmetrized DLR rank = " << dlr_rf.size() << std::endl;
    std::cout << "Symmetrized DLR rank = " << dlr_rf_sym.size() << std::endl;
  }
}