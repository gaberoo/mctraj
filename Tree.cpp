#include "Tree.h"

// ========================================================================

void Tree::reverse() {
  double mt = maxTime();
  for (size_t i = 0; i < times.size(); ++i) {
    times[i] = mt-times[i];
  }
  std::reverse(times.begin(),times.end());
  std::reverse(ttypes.begin(),ttypes.end());
  std::reverse(ids.begin(),ids.end());
  for (size_t i = 0; i < branches.size(); ++i) {
    branches[i].time = mt - branches[i].time;
  }
  is_rev = ! is_rev;
}

// ========================================================================

void Tree::addRateShift(double t) {
  // find position
  size_t pos(0);
  do {
    if (times[pos] > t) break;
    ++pos;
  } while (pos < times.size());
  vector<double>::iterator p1(times.begin());
  vector<int>::iterator p2(ttypes.begin());
  vector< vector<int> >::iterator p3(ids.begin());
  p1 += pos;
  p2 += pos;
  p3 += pos;
  times.insert(p1,t);
  ttypes.insert(p2,20);
  ids.insert(p3,vector<int>(3,-1));
}

// ========================================================================

int Tree::readFromFile(string fn) {
  ifstream in(fn.c_str());

  if (! in.is_open()) {
    cerr << "Problem opening file '" << fn << "'. Aborting." << endl;
    return -2;
  }

  string input;
  double x1;
  int    x2;
  int    A;

  while (! getline(in,input).eof()) {
    if (input.length() == 0) continue;
    if (input[0] == '#') continue;

    istringstream istr(input);

    if (istr >> x1 >> x2) {
      times.push_back(x1);
      ttypes.push_back(x2);
      ids.push_back(vector<int>());
    } else {
      cerr << "Error reading times/type:" << endl << " > " << fn << endl;
    }

    while (1) {
      if (istr >> A) ids.back().push_back(A);
      else if (istr.eof()) break;
    }
  }
  in.close();

  // make sure times are sorted
  vector<size_t> p(times.size());
  gsl_sort_index(p.data(),times.data(),1,times.size());

  // apply sort
  // gsl_permute(p.data(),times.data(),1,times.size());
  // gsl_permute_int(p.data(),ttypes.data(),1,ttypes.size());

  extant = nroot;
  maxExtant = nroot;
  numBranches = nroot;

  int n(times.size());

  for (int i(n-1); i >= 0; --i) {
    switch (ttypes[i]) {
      case 1:
      case 11:
        ++extant;
        numBranches += 2;
        break;
      case 0:
      case 2:
        --extant;
        break;
      default:
        break;
    }
    if (maxExtant < extant) maxExtant = extant;
  }

  if (extant < 0) {
    fprintf(stderr,"Invalid tree: More samples than infections.\n");
    cout << "extant = " << extant << endl;
    return -1;
  } 

  makeBranches();

  return 0;
}

// ========================================================================

void Tree::makeBranches() {
  branches.resize(max_id()+1);
  int id;
  for (size_t i = 0; i < ids.size(); ++i) {
    if (ids[i].size() > 0) {
      id = ids[i][0];
      if (id >= 0) {
        branches[id].id = id;
        branches[id].type = ttypes[i];
        branches[id].time = times[i];
        branches[id].child1 = ids[i][1];
        branches[id].child2 = ids[i][2];
      } else {
        cerr << "Invalid ID at position " << i << " : " << id << endl;
      }
    } else {
#ifdef DEBUG
      cerr << "No ID information!" << endl;
#endif
    }
  }
}

// ========================================================================

size_t Tree::countTypes(int type) const {
  size_t c = 0;
  for (size_t i = 0; i < ttypes.size(); ++i) {
    if (ttype(i) == type) ++c;
  }
  return c;
}

// ========================================================================

int Tree::max_id() const {
  int max = 0;
  for (size_t i = 0; i < ids.size(); ++i) {
    if (ids[i].size() > 0) {
      if (ids[i][0] > max) max = ids[i][0];
    }
  }
  return max;
}


