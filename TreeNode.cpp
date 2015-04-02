#include "TreeNode.h"

/****************************************************************************/

void MCTraj::TreeNode::operator=(const TreeNode& t) {
  code = t.code;
  parent = t.parent;
  age = t.age;
  n_off = t.n_off;
  extant_off = t.extant_off;
  ancestor_age = t.ancestor_age;
  epi_id = t.epi_id;
  epi_state = t.epi_state;
}

/****************************************************************************/

string MCTraj::TreeNode::to_json() const {
  rapidjson::StringBuffer buf;
  rapidjson::Writer<rapidjson::StringBuffer> json_w(buf);
  json(json_w);
  return buf.GetString();
}

/****************************************************************************/

ostream& MCTraj::operator<<(ostream& out, const TreeNode& tn) {
  out << "{"
      << "\"parent\": "     << tn.parent       << ","
      << "\"code\": "       << tn.code         << ","
      << "\"age\": "        << tn.age          << ","
      << "\"n_off\": "      << tn.n_off        << ","
      << "\"extant_off\": " << tn.extant_off   << ","
      << "\"parent_age\": " << tn.ancestor_age << ","
      << "\"epi_id\": "     << tn.epi_id       << ","
      << "\"epi_state\": "  << tn.epi_state
      << "}";
  return out;
}

/****************************************************************************/

int MCTraj::only_sampled(const vector<TreeNode>& tree, vector<TreeNode>& sample, int root) {
  int n = tree.size();
  if (root < 0 || root >= n) {
    return 0;
  } else {
    if (tree[root].off.size() == 0) {     // no offspring, ...
      if (tree[root].extant_off > 0) {    // ... but there are extant offspring ...
        sample.push_back(tree[root]);     // ... then add the node to the tree
        return 1;
      } else {
        return 0;
      }
    } else {
      int n_off = 0;
      int off = 0;
      int t_off = 0;
      for (size_t i = 0; i < tree[root].off.size(); ++i) {
        off = only_sampled(tree,sample,tree[root].off[i]);
        n_off += off;
        if (off > 0) ++t_off;
      }
      if (n_off > 0) {
        sample.push_back(tree[root]);
        sample.back().n_off = t_off;
      }
      return n_off;
    }
  }
}

/****************************************************************************/

int MCTraj::only_sampled(const vector<TreeNode>& tree, string& newick, int root) {
  int n = tree.size();
  if (root < 0 || root >= n) {
    newick = "";
    return 0;
  } else {
    if (tree[root].off.size() == 0) {
      if (tree[root].extant_off > 0) {
        ostringstream out;
        out << tree[root].epi_id << "_" << tree[root].epi_state;
        out << ":" << tree[root].age - tree[root].ancestor_age;
        newick = out.str();
        return 1;
      } else {
        newick = "";
        return 0;
      }
    } else {
      int n_off = 0;
      int off = 0;
      ostringstream out;
      string tmp;
      for (size_t i = 0; i < tree[root].off.size(); ++i) {
        off = only_sampled(tree,tmp,tree[root].off[i]);
        if (off > 0) {
          if (n_off > 0) out << ",";
          out << tmp;
          n_off += off;
        }
      }
      tmp = out.str();
      out.str("");
      out << "(" << tmp << ")";
      out << tree[root].epi_id << "_" << tree[root].epi_state;
      out << ":" << tree[root].age - tree[root].ancestor_age;
      newick = out.str();
      return n_off;
    }
  }
}

/****************************************************************************/

string MCTraj::to_newick(const vector<TreeNode>& tree, int root) {
  int n = tree.size();
  if (root < 0 || root >= n) return "";
  else {
    ostringstream out;
    if (tree[root].off.size() == 0) {
      out << tree[root].epi_id << "_" << tree[root].epi_state;
      if (tree[root].extant_off == 0) out << "+";
      out << ":" << tree[root].age - tree[root].ancestor_age;
      return out.str();
    } else {
      out << "(";
      for (size_t i = 0; i < tree[root].off.size(); ++i) {
        out << to_newick(tree,tree[root].off[i]);
        if (i < tree[root].off.size()-1) out << ",";
      }
      out << ")";
      out << tree[root].epi_id << "_" << tree[root].epi_state;
      out << ":" << tree[root].age - tree[root].ancestor_age;
      return out.str();
    }
  }
}

/****************************************************************************/

int MCTraj::add_extant(vector<TreeNode>& tree, int node) {
  int n = tree.size();
  if (node >= 0 && node < n) {
    tree[node].extant_off++;
    return add_extant(tree,tree[node].parent) + 1;
  } else {
    return 0;
  }
}

