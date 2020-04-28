/* Dstar.h
 * James Neufeld (neufeld@cs.ualberta.ca)
 * Compilation fixed by Arek Sredzki (arek@sredzki.com)
 */

#ifndef DSTAR_H
#define DSTAR_H

#include <math.h>
#include <stack>
#include <queue>
#include <list>
#include <stdio.h>
#include <ext/hash_map>

using namespace std;
using namespace __gnu_cxx;

class state {
 public:
  vector<int> dims;
  // int x;
  // int y;
  pair<double,double> k;

  state(uint size=2){
    dims.clear();
    dims.resize(size);
  }

  state operator = (const state &s2) {
    dims.clear();
    dims.resize(s2.size());
    for (uint i=0; i<size(); i++){
      set(i, s2.get(i));
    }
    return *this;
  }

  bool operator == (const state &s2) const {
    // return ((x == s2.x) && (y == s2.y));
    bool is_so = true;
    for (uint i=0; i<size(); i++){
      is_so &= (dims[i] == s2.dims[i]);
    }
    return is_so;
  }

  bool operator != (const state &s2) const {
    // return ((x != s2.x) || (y != s2.y));
    bool is_so = false;
    for (uint i=0; i<size(); i++){
      is_so |= (dims[i] != s2.dims[i]);
    }
    return is_so;
  }

  bool operator > (const state &s2) const {
    if (k.first-0.00001 > s2.k.first) return true;
    else if (k.first < s2.k.first-0.00001) return false;
    return k.second > s2.k.second;
  }

  bool operator <= (const state &s2) const {
    if (k.first < s2.k.first) return true;
    else if (k.first > s2.k.first) return false;
    return k.second < s2.k.second + 0.00001;
  }

  bool operator < (const state &s2) const {
    if (k.first + 0.000001 < s2.k.first) return true;
    else if (k.first - 0.000001 > s2.k.first) return false;
    return k.second < s2.k.second;
  }

  size_t size() const {
    return dims.size();
  }

  int get(const uint idx) const {
    return dims[idx];
  }
  void set(const uint idx, const int val) {
    dims[idx] = val;
  }

};

// struct ipoint2 {
//   int x,y;
// };

struct cellInfo {

  double g;
  double rhs;
  double cost;

};

class state_hash {
 public:
  size_t operator()(const state &s) const {
    // return s.x + 34245*s.y;
    size_t hash;
    for (uint i=0; i<s.size(); i++) {
      hash += pow(1000,i)*s.get(i);
    }
    return hash;
  }
};


typedef priority_queue<state, vector<state>, greater<state> > ds_pq;
typedef hash_map<state,cellInfo, state_hash, equal_to<state> > ds_ch;
typedef hash_map<state, float, state_hash, equal_to<state> > ds_oh;


class Dstar {

 public:

  Dstar(int);
  void   init(int sX, int sY, int gX, int gY);
  void   init(state s, state g);
  void   updateCell(int x, int y, double val);
  void   updateCell(state s, double val);
  void   updateStart(int x, int y);
  void   updateStart(state s);
  void   updateGoal(int x, int y);
  void   updateGoal(state s);
  bool   replan();
  void   draw();
  void   drawCell(state s,float z);
  int size() { return state_size; }

  list<state> getPath();

 private:

  list<state> path;

  double C1;
  double k_m;
  state s_start, s_goal, s_last;
  int maxSteps;
  uint state_size;

  ds_pq openList;
  ds_ch cellHash;
  ds_oh openHash;

  bool   close(double x, double y);
  void   makeNewCell(state u);
  double getG(state u);
  double getRHS(state u);
  void   setG(state u, double g);
  void   setRHS(state u, double rhs);
  double eightCondist(state a, state b);
  int    computeShortestPath();
  void   updateVertex(state u);
  void   insert(state u);
  void   remove(state u);
  double trueDist(state a, state b);
  double heuristic(state a, state b);
  state  calculateKey(state u);
  void   getSucc(state u, list<state> &s);
  void   getPred(state u, list<state> &s);
  double cost(state a, state b);
  bool   occupied(state u);
  bool   isValid(state u);
  float  keyHashCode(state u);
};

#endif
