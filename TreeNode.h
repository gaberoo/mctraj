#ifndef __TREENODE2_H__
#define __TREENODE2_H__

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
using namespace std;

namespace MCTraj {
  class TreeNode {
    public:
      TreeNode(int c = 0) 
        : code(c), parent(-1), age(0.0), n_off(0), 
          extant_off(0), ancestor_age(0.0), 
          epi_id(0), epi_state(0), off(vector<int>())
      {}

      TreeNode(const TreeNode& t)
        : code(t.code), parent(t.parent), age(t.age), n_off(t.n_off),
          extant_off(t.extant_off), ancestor_age(t.ancestor_age),
          epi_id(t.epi_id), epi_state(t.epi_state), off(t.off)
      {}

      virtual ~TreeNode() {}

      void operator=(const TreeNode& t);
      friend ostream& operator<<(ostream& out, const TreeNode& tn);

      inline void add_off(int id) { off.push_back(id); }

      int code;
      int parent;
      double age;
      int n_off;
      int extant_off;
      double ancestor_age;
      int epi_id;
      int epi_state;
      vector<int> off;
  };

  ostream& operator<<(ostream& out, const TreeNode& tn);

  string to_newick(const vector<TreeNode>& tree, int root);
  int only_sampled(const vector<TreeNode>& tree, vector<TreeNode>& sample, int root);
  int only_sampled(const vector<TreeNode>& tree, string& newick, int root);
}

#endif

