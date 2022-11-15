#include <numbers>
#include "utils.hpp"

using std::numbers::pi;
using namespace nda;

  barycheb::barycheb(int n):
    xc(n),wc(n) {
      
    auto c = (2*arange(1,n+1)-1)*pi/(2*n);
    xc = cos(c);
    wc = (1-2*(arange(0,n)%2))*sin(c);
    
  }
