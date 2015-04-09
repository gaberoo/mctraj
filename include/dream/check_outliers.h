#ifndef __DREAM_CHECK_OUTLIERS_H__
#define __DREAM_CHECK_OUTLIERS_H__

#include <vector>
#include <cstdlib>
#include <cmath>
#include <algorithm>
using namespace std;

#include <gsl/gsl_statistics.h>

#include "array.h"

void check_outliers(int t, Array2D<double>& lik, vector<double>& meanlik,
                    vector<bool>& outliers);

#endif
