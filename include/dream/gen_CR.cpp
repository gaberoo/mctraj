#include "gen_CR.h"

void gen_CR(const gsl_rng* rng, const vector<double>& pCR, 
    Array2D<int>& CRm, vector<unsigned>& L) 
{
  size_t numChains(CRm.n_x());
  size_t loopSteps(CRm.n_y());
  size_t size(numChains*loopSteps);
  if (L.size() < pCR.size()) {
    cerr << "Bad array size." << endl;
    return;
  }
  // pick candidates for each crossover value
  gsl_ran_multinomial(rng,pCR.size(),size,pCR.data(),L.data());
  size_t j(0);
  for (size_t m(0); m < L.size(); ++m) {
    for (size_t i(0); i < L[m]; ++i) {
      CRm[j] = m+1;
      ++j;
    }
  }
  gsl_ran_shuffle(rng,CRm.pt(0,0),size,sizeof(int));
}


