#ifndef __DREAM_H__
#define __DREAM_H__

#include <rapidjson/document.h>
#include <cstdlib>
#include <string>
#include <map>
using namespace std;

double Nmin = 0;           /* minimum total population size */
double Nmax = 1000;        /* maximum total population size */
int root = 0;              /* add a root to the tree */
int SImodel = 1;           /* use non-saturating model */
int n = 0;                 /* number of times */
int rescale = 1;           /* rescale probabilities in expomv */
int extant = 0;            /* number of extant species */
int maxExtant = 0;         /* max number of extant species */
double fixedRatio = -1.0;  /* fix mu-phi ratio */
char jointLikelihood = 'j';   /* forest likelihood type */
int survival = 1;          /* correct for survival probability of the tree */
int modelType = 1;
int reverseTrees = 1;
int saveTraj = 0;

string json_input = "[]";

typedef struct t_dream_pars {
  int vflag;                 /* vebose flag */
  int maxEvals;              /* max number of function evaluations */
  int optimAlg;              /* lock variables */
  int numChains;    
  string fn;                 /* input filename */
  string out_fn;             /* output filename */
  int appendFile;            /* continue from previous state */
  int report_interval;       /* report interval for state */
  int diagnostics;           /* report diagnostics at the end of the run */
  int burnIn;                /* number of steps for which to run an adaptive proposal size */

  // DREAM variables
  int collapseOutliers;
  int gelmanEvals;
  int loopSteps;
  int numPars;
  int deltaMax;
  int pCR_update;
  int nCR;

  double noise;
  double bstar_zero;
  double scaleReductionCrit;
  double reenterBurnin;

  int nfree;
  size_t nvar;
  double* varLo;
  double* varHi;
  double* varInit;
  int* varLock;
  string* varName;
} dream_pars;

void dream_pars_default(dream_pars* p);
void dream_pars_init_vars(dream_pars* p, size_t n);
void dream_pars_free_vars(dream_pars* p);

// jpars.Parse<0>(json_input.c_str());
// assert(jpars.IsObject());
void dream_pars_read_json(dream_pars* p, rapidjson::Document& jpars);

#endif

