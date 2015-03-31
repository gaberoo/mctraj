#ifndef __TREENODE2_H__
#define __TREENODE2_H__

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
using namespace std;

#include "MCTraj.h"

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

      // =====================================================================

      template<typename T> 
      void json(rapidjson::Writer<T>& json_w) const {
        json_w.StartObject(); {
          json_w.String("parent");     json_w.Int(parent);
          json_w.String("code");       json_w.Int(code);
          json_w.String("age");        json_w.Double(age);
          json_w.String("n_off");      json_w.Int(n_off);
          json_w.String("extant_off"); json_w.Int(extant_off);
          json_w.String("parent_age"); json_w.Double(ancestor_age);
          json_w.String("epi_id");     json_w.Int(epi_id);
          json_w.String("epi_state");  json_w.Int(epi_state);
        } json_w.EndObject();
      }
      string to_json() const;

      // =====================================================================

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

  int add_extant(vector<TreeNode>& tree, int node);
}

#endif

