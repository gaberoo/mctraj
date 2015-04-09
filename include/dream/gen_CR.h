#ifndef __DREAM_GEN_CR_H__
#define __DREAM_GEN_CR_H__

#include <vector>
#include <iostream>
using namespace std;

#include "dream.h"
#include "array.h"

void gen_CR(rng::RngStream* rng, const vector<double>& pCR, 
            Array2D<int>& CRm, vector<unsigned>& L);

#endif
