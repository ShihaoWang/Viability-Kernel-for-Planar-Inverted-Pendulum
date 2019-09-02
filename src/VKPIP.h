#ifndef VKPIP_H
#define VKPIP_H
#include <iostream>
#include <math.h>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <set>
#include <ctime>
using namespace std;

template <typename T>
std::vector<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N-1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}

template <typename T>
std::vector<int> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  std::vector<int> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](int i1, int i2) {return v[i1] < v[i2];});

  return idx;
}

struct SystemState
{
  // This is the struct that saves the information of the robot state
  SystemState(){ L = 0.0; Ldot = 0.0; Theta = 0.0; Thetadot = 0.0;}
  SystemState(double _L, double _Ldot, double _Theta, double _Thetadot):L(_L), Ldot(_Ldot), Theta(_Theta), Thetadot(_Thetadot){}
  void Update(double _L, double _Ldot, double _Theta, double _Thetadot)
  {
    L = _L;
    Ldot = _Ldot;
    Theta = _Theta;
    Thetadot = _Thetadot;
  }
  double L, Ldot, Theta, Thetadot;
};
struct SystemIndex
{
  // This is the struct of the system index
  SystemIndex(){L_index = 0; Ldot_index = 0; Theta_index = 0; Thetadot_index = 0;}
  SystemIndex(int _L_index, int _Ldot_index, int _Theta_index, int _Thetadot_index):L_index(_L_index), Ldot_index(_Ldot_index), Theta_index(_Theta_index), Thetadot_index(_Thetadot_index){}
  int L_index, Ldot_index, Theta_index, Thetadot_index;
};
typedef struct DBNode *DBNodePtr;    // This is a pointer to the Tree_Node
struct DBNode
{
  // Each data point in the database is considered to be a tree node linked all the way up down to the final node whose obj value is the minimum.
  DBNode()
  {
    SelfIndex = -1;
    NextIndex = -2;               // -2 means has not been initialized, -1 means no next node, >0 gives node index of the next node.
    TransFeasibleFlag = 1;       // default to be feasible

    Viable = false;
    Objective = 0.0;
  }
  DBNode(const int& L_index, const int& Ldot_index, const int& Theta_index, const int& Thetadot_index, const int& L_Grids, const int& Ldot_Grids, const int& Theta_Grids, const int& Thetadot_Grids, const std::vector<double>& L_vector, const std::vector<double>& Ldot_vector, const std::vector<double>& Theta_vector, const std::vector<double>& Thetadot_vector)
  {
    // Node instantiation
    SelfIndex = L_index * Ldot_Grids * Theta_Grids * Thetadot_Grids + Ldot_index * Theta_Grids * Thetadot_Grids + Theta_index * Thetadot_Grids + Thetadot_index;
    NextIndex = -2;
    TransFeasibleFlag = 1;
    Viable = false;
    Objective = 0.0;

    NodeIndex.L_index =         L_index;
    NodeIndex.Ldot_index =      Ldot_index;
    NodeIndex.Theta_index =     Theta_index;
    NodeIndex.Thetadot_index =  Thetadot_index;

    NodeState.L =               L_vector[L_index];
    NodeState.Ldot =            Ldot_vector[Ldot_index];
    NodeState.Theta =           Theta_vector[Theta_index];
    NodeState.Thetadot =        Thetadot_vector[Thetadot_index];
  }

  void InitUpdate(const int& L_index, const int& Ldot_index, const int& Theta_index, const int& Thetadot_index, const int& L_Grids, const int& Ldot_Grids, const int& Theta_Grids, const int& Thetadot_Grids, const std::vector<double>& L_vector, const std::vector<double>& Ldot_vector, const std::vector<double>& Theta_vector, const std::vector<double>& Thetadot_vector)
  {
    SelfIndex = L_index * Ldot_Grids * Theta_Grids * Thetadot_Grids + Ldot_index * Theta_Grids * Thetadot_Grids + Theta_index * Thetadot_Grids + Thetadot_index;
    NextIndex = -2;
    TransFeasibleFlag = 1;
    Viable = false;
    Objective = 0.0;

    NodeIndex.L_index =         L_index;
    NodeIndex.Ldot_index =      Ldot_index;
    NodeIndex.Theta_index =     Theta_index;
    NodeIndex.Thetadot_index =  Thetadot_index;

    NodeState.L =               L_vector[L_index];
    NodeState.Ldot =            Ldot_vector[Ldot_index];
    NodeState.Theta =           Theta_vector[Theta_index];
    NodeState.Thetadot =        Thetadot_vector[Thetadot_index];
  }

  // Each node has the following information
  int SelfIndex;                // The index of its own
  int NextIndex;                // The index of the next node it should go according to the principle of optimality
  int TransFeasibleFlag;       // This flag is used for state transition where the next state is not feasible for dynamical system.
  bool Viable;                  // This boolean variable is used to check whether the state is viable or not.

  float Objective;              // We choose the objective to be the time to reach goal.

  SystemIndex NodeIndex;
  SystemState NodeState;
};
#endif
