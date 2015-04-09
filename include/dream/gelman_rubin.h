#ifndef __GELMAN_RUBIN_H__
#define __GELMAN_RUBIN_H__

#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics.h>

#include "array.h"

void gelman_rubin(Array3D<double>& state, vector<double>& scaleReduction, 
    const vector<bool>& lockVar, int first = 0, int numIter = -1,
    int adjustDF = 0);

#endif // __GELMAN_RUBIN_H__
