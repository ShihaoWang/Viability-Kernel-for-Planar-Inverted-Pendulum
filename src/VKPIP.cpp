#include "VKPIP.h"

// This function computes the failure metric with the forward value iteration method.
// One more change that I will have to make is the control input: I would like to change the control input to be robot's length acceleration.

const double PI = 3.141592653589793238463;

// The lower and upper bound of the pendulum variables
double LLow = 0.35;             double LUpp = 1.05;
double LdotLow = -1.0;          double LdotUpp = 1.0;
double ThetaLow = -PI/6.0;      double ThetaUpp = PI/2.0;
double ThetadotLow = -1.5;      double ThetadotUpp = 1.5;

const int L_Grids = 71;                             // Its unit 0.01 m
const int Ldot_Grids = 21;                          // Its unit 0.1m/s
const int Theta_Grids = 211;                        // Its unit 9.97x10^(-3) rad
const int Thetadot_Grids = 31;                      // Its unit 0.1rad/s

double delta_t = 0.1;
const int F_Grids = 101;                           // This is the discretization of the control points within a range

double g = 9.81;
double gUppCoef = 2.25;

// The corresponding dimension size of the StateMatrix
double L_length =         LUpp - LLow;
double Ldot_length =      LdotUpp - LdotLow;
double Theta_length =     ThetaUpp - ThetaLow;
double Thetadot_length =  ThetadotUpp - ThetadotLow;

double L_unit =           L_length/(1.0 * L_Grids - 1.0);
double Ldot_unit =        Ldot_length/(1.0 * Ldot_Grids - 1.0);
double Theta_unit =       Theta_length/(1.0 * Theta_Grids - 1.0);
double Thetadot_unit =    Thetadot_length/(1.0 * Thetadot_Grids - 1.0);

std::vector<double> L_vector(L_Grids), Ldot_vector(Ldot_Grids), Theta_vector(Theta_Grids), Thetadot_vector(Thetadot_Grids);

static void StateVectorSubs(std::vector<double>& L_vector,std::vector<double>& Ldot_vector,std::vector<double>& Theta_vector, std::vector<double>& Thetadot_vector)
{
  // L_vector filling
  for (int i = 0; i < L_Grids; i++)
  {
    L_vector[i] = LLow + (1.0 * i) * L_unit;
  }
  // Ldot_vector filling
  for (int i = 0; i < Ldot_Grids; i++)
  {
    Ldot_vector[i] = LdotLow + (1.0 * i) * Ldot_unit;
  }
  // Theta_vector filling
  for (int i = 0; i < Theta_Grids; i++)
  {
    Theta_vector[i] = ThetaLow + (1.0 * i) * Theta_unit;
  }
  // Thetadot_vector filling
  for (int i = 0; i < Thetadot_Grids; i++)
  {
    Thetadot_vector[i] = ThetadotLow + (1.0 * i) * Thetadot_unit;
  }
}

static SystemIndex SystemState2Index(const SystemState & x)
{
  // This function is used to transform the System State from double into the corresponding index
  // Here the x must have already been checked with the given pendulum length bound.

  double L_FloatIndex =         (x.L - LLow)/L_unit * 1.0;
  double Ldot_FloatIndex =      (x.Ldot - LdotLow)/Ldot_unit * 1.0;
  double Theta_FloatIndex =     (x.Theta - ThetaLow)/Theta_unit * 1.0;
  double Thetadot_FloatIndex =  (x.Thetadot - ThetadotLow)/Thetadot_unit * 1.0;

  int L_Index =         std::round(L_FloatIndex);
  int Ldot_Index =      std::round(Ldot_FloatIndex);
  int Theta_Index =     std::round(Theta_FloatIndex);
  int Thetadot_Index =  std::round(Thetadot_FloatIndex);

  // The L Index reshift
  if(L_Index<0)
  {
    L_Index = 0;
  }
  if(L_Index>=L_Grids)
  {
    L_Index = L_Grids-1;
  }

  // The Ldot Index reshift
  if(Ldot_Index<0)
  {
    Ldot_Index = 0;
  }
  if(Ldot_Index>=Ldot_Grids)
  {
    Ldot_Index = Ldot_Grids-1;
  }

  // The Theta Index reshift
  if(Theta_Index<0)
  {
    Theta_Index = 0;
  }
  if(Theta_Index>=Theta_Grids)
  {
    Theta_Index = Theta_Grids-1;
  }

  // The Thetadot Index reshift
  if(Thetadot_Index<0)
  {
    Thetadot_Index = 0;
  }
  if(Thetadot_Index>=Thetadot_Grids)
  {
    Thetadot_Index = Thetadot_Grids-1;
  }

  SystemIndex SystemIndex_i(L_Index, Ldot_Index, Theta_Index, Thetadot_Index);
  return SystemIndex_i;
}

SystemIndex SystemState2StateIndex(const SystemState & x)
{
  // This function is used to transform the System State from double into the corresponding index
  // Here the x must have already been checked with the given pendulum length bound.

  double L_FloatIndex =         (x.L - LLow)/L_unit * 1.0;
  double Ldot_FloatIndex =      (x.Ldot - LdotLow)/Ldot_unit * 1.0;
  double Theta_FloatIndex =     (x.Theta - ThetaLow)/Theta_unit * 1.0;
  double Thetadot_FloatIndex =  (x.Thetadot - ThetadotLow)/Thetadot_unit * 1.0;

  int L_Index =         std::round(L_FloatIndex);
  int Ldot_Index =      std::round(Ldot_FloatIndex);
  int Theta_Index =     std::round(Theta_FloatIndex);
  int Thetadot_Index =  std::round(Thetadot_FloatIndex);

  // The L Index reshift
  if(L_Index<0)
  {
    L_Index = 0;
  }
  if(L_Index>=L_Grids)
  {
    L_Index = L_Grids-1;
  }

  // The Ldot Index reshift
  if(Ldot_Index<0)
  {
    Ldot_Index = 0;
  }
  if(Ldot_Index>=Ldot_Grids)
  {
    Ldot_Index = Ldot_Grids-1;
  }

  // The Theta Index reshift
  if(Theta_Index<0)
  {
    Theta_Index = 0;
  }
  if(Theta_Index>=Theta_Grids)
  {
    Theta_Index = Theta_Grids-1;
  }

  // The Thetadot Index reshift
  if(Thetadot_Index<0)
  {
    Thetadot_Index = 0;
  }
  if(Thetadot_Index>=Thetadot_Grids)
  {
    Thetadot_Index = Thetadot_Grids-1;
  }

  SystemIndex SystemIndex_i(L_Index, Ldot_Index, Theta_Index, Thetadot_Index);
  return SystemIndex_i;
}

void ObjectiveUpdate(DBNode& Node_k, const bool & ViableFlag)
{
  // This function is only used in state initialization stage.
  // A positive value indicates that the given state is viable, otherwise the state is not viable.
  switch (ViableFlag)
  {
    case true:
    {
      // This indicates that the current Node_k is at goal
      Node_k.Objective = 0.0;
    }
    break;
    default:
    {
      Node_k.Objective = -1.0;      // A non-viable state will be assigned a negative value.
    }
    break;
  }
}

bool GoalStateTest(const DBNode & Node_k)
{
  // This function is used to tell whether a given state is at goal or not.
  double L = Node_k.NodeState.L;
  double Ldot = Node_k.NodeState.Ldot;
  double Theta = Node_k.NodeState.Theta;
  double Thetadot = Node_k.NodeState.Thetadot;

  // Here the goal state is a state where robot's state lies within safe side while its velocity is non-negative.
  if ((Theta>0)&&(Thetadot>=0)&&(Ldot>=0))
  {
    return true;
  }
  return false;
}

void DBNodeInfoUpdate(DBNode & Node_k)
{
  // This function can only be called after Node_k has been initialized

  double Node_kp1_L = Node_k.NodeState.L + Node_k.NodeState.Ldot * delta_t;
  if((Node_kp1_L<LLow)||(Node_kp1_L>LUpp))
  {
    Node_k.TransFeasibleFlag = 0;
  }
  switch (Node_k.TransFeasibleFlag)
  {
    case 0:
    {
      Node_k.Objective = -1.0;
    }
    break;
    default:
    {
      bool GoalFlag = GoalStateTest(Node_k);
      Node_k.Viable = GoalFlag;
      ObjectiveUpdate(Node_k, GoalFlag);
    }
    break;
  }
}

double Forward_Evaluation(DBNode& Node_k, Eigen::Tensor<DBNode,4>& StateNodeMatrix)
{
  // This function conducts the first-order Euler integration from the kth-layer node to the k+1th layer node.
  /*
    L_{k+1} =         L_{k} + Ldot_{k} * delta_t;
    Ldot_{k+1} =      Ldot_{k} + Lddot_{k} * delta_t;
    Theta_{k+1} =     Theta_{k} + Thetadot_{k} * delta_t;
    Thetadot_{k+1} =  Thetadot_{k} + Thetaddot_{k} * delta_t;
  */

  // If the output is greater than 0, this indicates that this node has found another path to reach goal with a lower cost.

  std::set< std::tuple<int, int, int, int>> Reachable_Set;

  // kth layer
  double L_k=           Node_k.NodeState.L;
  double Ldot_k =       Node_k.NodeState.Ldot;
  double Theta_k =      Node_k.NodeState.Theta;
  double Thetadot_k =   Node_k.NodeState.Thetadot;
  double Thetaddot_k =  g/L_k * sin(Theta_k) - 2.0 * Thetadot_k * Ldot_k/L_k;       // Here the value of g needs to be updated!

  // k+1th layer
  double L_kp1, Ldot_kp1, Theta_kp1, Thetadot_kp1;

  switch (Node_k.TransFeasibleFlag)
  {
    case 0:
    {
      // This means that there is no need to integrate this node.
      return 0.0;
    }
    break;
    default:
    {
      L_kp1 = L_k + Ldot_k * delta_t;     // Update the next node
    }
    break;
  }

  Theta_kp1 = Theta_k + Thetadot_k * delta_t;
  Thetadot_kp1 = Thetadot_k + Thetaddot_k * delta_t;

  // The bounds on the acceleration of Pendulum Length
  double Lddot_Low = -g * cos(Theta_k);
  Lddot_Low = min(Lddot_Low, 0.0);
  double Lddot_Upp = gUppCoef * g * cos(Theta_k)/(LLow - LUpp) * (L_k - LLow) + gUppCoef * g * cos(Theta_k);
  Lddot_Upp = max(Lddot_Upp, 0.0);
  double Lddot_Min = Lddot_Low;
  double Lddot_Max = max(Lddot_Low, Lddot_Upp);
  std::vector<double> Lddot_vec = linspace(Lddot_Min, Lddot_Max, F_Grids);

  // Ldot integration
  for (int i = 0; i < F_Grids; i++)
  {
    Ldot_kp1 = Ldot_k + Lddot_vec[i] * delta_t;     // Now Lddot_vec saves the acceleration of pendulum length at k+1th node
    SystemState x_kp1(L_kp1, Ldot_kp1, Theta_kp1, Thetadot_kp1);
    SystemIndex x_kp1_index = SystemState2StateIndex(x_kp1);
    Reachable_Set.insert(make_tuple(x_kp1_index.L_index, x_kp1_index.Ldot_index, x_kp1_index.Theta_index, x_kp1_index.Thetadot_index));
  }

  float ObjectiveVia = 0.0;

  switch (Reachable_Set.size())
  {
    case 0:
    {
      return ObjectiveVia;
    }
    default:
    {
      const int ReachableNodeNumber = Reachable_Set.size();
      std::vector<double> Reachable_Objective(ReachableNodeNumber);
      std::set<std::tuple<int, int, int, int>>::iterator Reachable_Set_Itr = Reachable_Set.begin();
      for (int i = 0; i < ReachableNodeNumber; i++)
      {
        switch (i)
        {
          case 0:
          break;
          default:
          std::advance(Reachable_Set_Itr, 1);
          break;
        }
        // Here DBNodePtr points to the node at (k+1)th layer
        DBNode* DBNodePtr =&StateNodeMatrix(std::get<0>(*Reachable_Set_Itr),std::get<1>(*Reachable_Set_Itr),std::get<2>(*Reachable_Set_Itr),std::get<3>(*Reachable_Set_Itr));

        float Reachable_Objective_i = 0.0;
        switch (DBNodePtr->Viable)
        {
          case true:
          {
            double UpdatedObj = DBNodePtr->Objective + delta_t;
            switch (Node_k.Viable)
            {
              case true:
              {
                if(Node_k.Objective>UpdatedObj)
                {
                  // This indicates that we have a better path.
                  Reachable_Objective_i = Node_k.Objective - UpdatedObj;
                  Node_k.Objective = UpdatedObj;
                  Node_k.NextIndex = DBNodePtr->SelfIndex;
                }
                else
                {
                  // This means that the current path is not as good as the current one.
                  Reachable_Objective_i = 0.0;
                }
              }
              break;
              default:
              {
                // This means that the current node is not viable, so we should make it to be viable.
                Node_k.Viable = true;
                Reachable_Objective_i = UpdatedObj;
                Node_k.Objective = UpdatedObj;
                Node_k.NextIndex = DBNodePtr->SelfIndex;
              }
              break;
            }
            Reachable_Objective[i] = Reachable_Objective_i;
          }
          break;
          default:
          {
            Reachable_Objective[i] = 0.0;
          }
          break;
        }
      }
      ObjectiveVia = *std::max_element(Reachable_Objective.begin(), Reachable_Objective.end());
    }
    break;
  }
  return ObjectiveVia;
}

void HJBDataWriter(const std::vector<int>& NodeParentIndVector, const std::vector<float>& NodeCostVector, const int & Angle_i)
{
  // *. Node Transition
  FILE * StateTransFile = NULL;
  string StateTransFileName = "PVKNextInd" + std::to_string(Angle_i) + ".bin";
  const char * StateTransFileName_Char = StateTransFileName.c_str();
  StateTransFile = fopen(StateTransFileName_Char, "wb");
  fwrite(&NodeParentIndVector[0], sizeof(int), NodeParentIndVector.size(), StateTransFile);
  fclose(StateTransFile);

  // *. Cost
  FILE * StateCostFile = NULL;
  string StateCostFileName = "PVKFailureMetric" + std::to_string(Angle_i) + ".bin";
  const char * StateCostFileName_Char = StateCostFileName.c_str();
  StateCostFile = fopen(StateCostFileName_Char, "wb");
  fwrite(&NodeCostVector[0], sizeof(float), NodeCostVector.size(), StateCostFile);
  fclose(StateCostFile);
}

void StateNodeMatrixInit(Eigen::Tensor<DBNode,4> &StateNodeMatrix, std::vector<float> &NodeFailureMetricVector)
{
  // This function is used to initialize the State Node Matrix
  DBNode * Node_ptr;
  for (int i = 0; i < L_Grids; i++)
  {
    for (int j = 0; j < Ldot_Grids; j++)
    {
      for (int k = 0; k < Theta_Grids; k++)
      {
        for (int l = 0; l < Thetadot_Grids; l++)
        {
          Node_ptr = &StateNodeMatrix(i,j,k,l);
          Node_ptr->InitUpdate(i, j, k, l, L_Grids, Ldot_Grids, Theta_Grids, Thetadot_Grids, L_vector, Ldot_vector, Theta_vector, Thetadot_vector);
          DBNodeInfoUpdate(*Node_ptr);
          NodeFailureMetricVector[Node_ptr->SelfIndex] = Node_ptr->Objective;
        }
      }
    }
  }
}

int main()
{
  std::vector<double> StateVectorSpecs;

  StateVectorSpecs.push_back(LLow);
  StateVectorSpecs.push_back(LUpp);
  StateVectorSpecs.push_back(LdotLow);
  StateVectorSpecs.push_back(LdotUpp);
  StateVectorSpecs.push_back(ThetaLow);
  StateVectorSpecs.push_back(ThetaUpp);
  StateVectorSpecs.push_back(ThetadotLow);
  StateVectorSpecs.push_back(ThetadotUpp);

  StateVectorSpecs.push_back(L_Grids * 1.0);
  StateVectorSpecs.push_back(Ldot_Grids * 1.0);
  StateVectorSpecs.push_back(Theta_Grids * 1.0);
  StateVectorSpecs.push_back(Thetadot_Grids * 1.0);

  const int AngleLow = 0;
  const int AngleUpp = 80;
  const int AngleDiff = 10;

  // This part needs to be updated with new information.
  StateVectorSpecs.push_back(AngleLow * 1.0);
  StateVectorSpecs.push_back(AngleUpp * 1.0);
  StateVectorSpecs.push_back(AngleDiff * 1.0);
  StateVectorSpecs.push_back(delta_t);        // The last element is the timestep duration.

  // Save the PVkDataSpecs first into PVkDataSpecs.bin
  FILE * StateVectorSpecsFile = NULL;
  StateVectorSpecsFile = fopen("PVKDataSpecs.bin", "wb");
  fwrite(&StateVectorSpecs[0], sizeof(double), StateVectorSpecs.size(), StateVectorSpecsFile);
  fclose(StateVectorSpecsFile);

  /*
    Main computation objects initialization
  */
  StateVectorSubs(L_vector, Ldot_vector, Theta_vector, Thetadot_vector);

  for (int i = 0; i < (AngleUpp - AngleLow)/AngleDiff + 1; i++)
  {
    int Angle_i = AngleLow + AngleDiff * i;
    g = 9.81 * cos(Angle_i*1.0/180 * 3.1415926535897);
    cout<<g<<endl;
    Eigen::Tensor<DBNode,4> StateNodeMatrix(L_Grids, Ldot_Grids, Theta_Grids, Thetadot_Grids);

    std::vector<int> NodeTransIndVector(L_Grids * Ldot_Grids * Theta_Grids * Thetadot_Grids);
    std::vector<float> NodeFailureMetricVector(L_Grids * Ldot_Grids * Theta_Grids * Thetadot_Grids);
    StateNodeMatrixInit(StateNodeMatrix, NodeFailureMetricVector);

    std::clock_t start; double duration; start = std::clock();

    double FailureMetricVia = 0.01;              // Here this value is the evaluation of the failure metric change to determine the convergence.
    double FailureMetricTol = 1e-10;            // Tolerence for convergence completeness
    double FailureMetricVia_i;

    while(FailureMetricVia>FailureMetricTol)
    {
      // Now this inner for loop turns out to be very easy to be solved with value iteration.
      double FailureMetricVia_ref = 0.0;
      for (int i = 0; i < L_Grids; i++)
      {
        for (int j = 0; j < Ldot_Grids; j++)
        {
          for (int k = 0; k < Theta_Grids; k++)
          {
            for (int l = 0; l < Thetadot_Grids; l++)
            {
              double FailureMetricVia_i = Forward_Evaluation(StateNodeMatrix(i,j,k,l), StateNodeMatrix);
              // double FailureMetricVia_i = Forward_Evaluation(StateNodeMatrix(42, 14, 101, 20), StateNodeMatrix);
              if(FailureMetricVia_i>FailureMetricVia_ref)
              {
                FailureMetricVia_ref = FailureMetricVia_i;
              }
            }
          }
        }
      }
      FailureMetricVia = FailureMetricVia_ref;
      printf ("Objective Viation: %f \n", FailureMetricVia);
    }

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    printf ("Total running time: %f s\n", duration);

    int IterIndex = 0;
    for (int i = 0; i < L_Grids; i++)
    {
      for (int j = 0; j < Ldot_Grids; j++)
      {
        for (int k = 0; k < Theta_Grids; k++)
        {
          for (int l = 0; l < Thetadot_Grids; l++)
          {
            NodeTransIndVector[IterIndex] = StateNodeMatrix(i,j,k,l).NextIndex;
            NodeFailureMetricVector[IterIndex] = StateNodeMatrix(i,j,k,l).Objective;
            IterIndex = IterIndex + 1;
          }
        }
      }
    }
    HJBDataWriter(NodeTransIndVector, NodeFailureMetricVector, Angle_i);
  }
  return 0;
}
