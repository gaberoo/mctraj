#include <iostream>
using namespace std;

#include <gsl/gsl_rng.h>

#include "MCTraj.h"
#include "Model.h"
#include "models/Models.h"
using namespace MCTraj;

#ifdef MKLRNG
#include <rng/MKLRng.h>
#else
#include <rng/GSLRng.h>
#endif

#include <tclap/CmdLine.h>

int main(int argc, char** argv) {
  TCLAP::CmdLine cmd("Simulate trees using particle filter models", ' ', "0.9");

  TCLAP::ValueArg<double> N("N","popSize","Total population size",false,1.0,"double",cmd);
  TCLAP::ValueArg<double> beta("b","beta","Transmission rate",true,1.0,"double",cmd);
  TCLAP::ValueArg<double> mu("u","mu","Recovery rate",true,0.1,"double",cmd);
  TCLAP::ValueArg<double> psi("s","psi","Sequential sampling rate",true,0.1,"double",cmd);
  TCLAP::ValueArg<double> rho("o","rho","Homochroneous sampling rate",true,0.5,"double",cmd);
  TCLAP::ValueArg<double> gamma("g","gamma","transition rate (for SEIR)",false,0.1,"double",cmd);

  TCLAP::ValueArg<double> mT("T","maxTime","Maximum simulation time",false,100.0,"double",cmd);
  TCLAP::ValueArg<char> type("t","model","Model type",true,'I',"char",cmd);
  TCLAP::ValueArg<unsigned long> seed("S","seed","Seed",false,time(NULL),"int",cmd);

  try {
    cmd.parse(argc,argv);
  } catch (TCLAP::ArgException &e) { 
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl; 
  }

  EpiState* init;
  Model* model;

  double maxTime = mT.getValue();

  IModel::EpiPars i_pars;
  i_pars.beta  = beta.getValue() / N.getValue();
  i_pars.mu    = mu.getValue();
  i_pars.psi   = psi.getValue();
  i_pars.rho   = rho.getValue();

  SISModel::EpiPars sis_pars;
  sis_pars.N     = N.getValue();
  sis_pars.beta  = beta.getValue();
  sis_pars.mu    = mu.getValue();
  sis_pars.psi   = psi.getValue();
  sis_pars.rho   = rho.getValue();

  SEISModel::EpiPars seis_pars;
  seis_pars.N     = N.getValue();
  seis_pars.beta  = beta.getValue();
  seis_pars.mu    = mu.getValue();
  seis_pars.psi   = psi.getValue();
  seis_pars.rho   = rho.getValue();
  seis_pars.gamma = gamma.getValue();

  switch (type.getValue()) {
    case 'I': // I
      model = new I(&i_pars);
      init = new EpiState(IModel::nstates);
      (*init)[0] = 1;
      (*init)[1] = 0;
      break;

    case 'S': // SIS
      model = new SIS(&sis_pars);
      init = new EpiState(SISModel::nstates);
      (*init)[0] = sis_pars.N-1;
      (*init)[1] = 1;
      (*init)[2] = 0;
      break;

    case 'R': // SIR
      model = new SIR(&sis_pars);
      init = new EpiState(SISModel::nstates);
      (*init)[0] = sis_pars.N-1;
      (*init)[1] = 1;
      (*init)[2] = 0;
      (*init)[3] = 0;
      break;

    case 'E': // SEIR
      model = new SEIS(&seis_pars);
      init = new EpiState(SEISModel::nstates);
      (*init)[0] = seis_pars.N-1;
      (*init)[1] = 0;
      (*init)[2] = 1;
      (*init)[3] = 0;
      (*init)[4] = 0;
      break;

    default:
      cerr << "Not yet implemented!" << endl;
      return 0;
      break;
  }

  Trajectory traj(*init,model);

  // Setup random number generators
  rng::RngStream* rng = NULL;
#ifdef MKLRNG
  rng = new rng::MKLStream;
#else
  rng = new rng::GSLStream;
#endif
  rng->alloc(seed.getValue());

  traj.simulateTrajectory(maxTime,&seis_pars,rng,false);

  vector<TreeNode> tree;
  vector<TreeNode> phylo;

  model->toTree(traj,rng,tree);

  // for (size_t i = 0; i < tree.size(); ++i) cout << tree[i] << endl;

//  string newick = to_newick(tree,0);
//  string sample_newick;
//  only_sampled(tree,sample_newick,0);
//  cout << "begin trees;" << endl;
//  cout << "tree 'full_tree' = " << newick << ";" << endl;
//  cout << "tree 'sampled_tree' = " << sample_newick << ";" << endl;
//  cout << "end;" << endl;

  only_sampled(tree,phylo,0);

  rapidjson::StringBuffer buf;
  rapidjson::Writer<rapidjson::StringBuffer> json_w(buf);
  json_w.StartObject(); {
    json_w.String("model"); {
      switch (type.getValue()) {
        case 'E': seis_pars.json(json_w); break;
        case 'R': 
        case 'S': sis_pars.json(json_w); break;
        case 'I': i_pars.json(json_w); break;
        default: json_w.String(""); break;
      }
    }
    json_w.String("traj"); traj.json(json_w);
    json_w.String("full_tree"); json_w.StartArray(); {
      for (size_t i = 0; i < tree.size(); ++i) tree[i].json(json_w);
    } json_w.EndArray();
    json_w.String("sampled"); json_w.StartArray(); {
      for (size_t i = 0; i < phylo.size(); ++i) phylo[i].json(json_w);
    } json_w.EndArray();
  } json_w.EndObject();

  cout << buf.GetString() << endl;

  delete model;
  delete init;

  return 0;
}
