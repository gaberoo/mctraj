#ifndef __DREAM_GEN_CR_H__
#define __DREAM_GEN_CR_H__

#include <vector>
#include <iostream>
using namespace std;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "array.h"

void gen_CR(const gsl_rng* rng, const vector<double>& pCR, 
            Array2D<int>& CRm, vector<unsigned>& L);

#endif
