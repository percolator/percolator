#ifndef NET_H_
#define NET_H_

/********** Neural Net stuff****************************/
 
class State
{
 public:
  State():len(0), x((double*)0), dx((double*)0) {}
  void resize(int n) {if(len != n) {clear(); x = new double[n]; dx = new double[n]; len = n;}}
  void clear(){ if(len != 0) {delete[] x; delete[] dx; len = 0;}}
  ~State(){clear();}
 
  int len;
  double *x;
  double *dx;
 
};
 

class Sigmoid
{
public:
  Sigmoid():num_neurons(0), a((double*)0){}
  ~Sigmoid() {clear();}
  void clear(){if(num_neurons!=0) delete [] a; num_neurons = 0;}
  void resize(int m);
  Sigmoid& operator=(Sigmoid &S);
 
  inline int getNumNeurons(){return num_neurons;}
 
  void write_to_file(ofstream &outfile);
  void read_from_file(ifstream &infile);
 
  void classify(State& down, State &up);
  void adjust_weights(State &down, State &up);
   
 protected:
  int num_neurons;
  double *a;
};




class Linear
{
 public:
  Linear(): num_neurons(0), len_w(0), weightDecay(0.0) {}
  ~Linear() {clear();}
  void resize(int m, int n);
  void clear();
  Linear& operator=(Linear &L);
 
  inline void set_mu(double m){ mu = m;}
  inline double get_mu(){return mu;}
  inline void set_num_feat(int nf) {num_features = nf;}
  inline void set_len_w(int len){len_w=len;}
  inline int getNumFeatures(){return num_features;}
  inline int getNumNeurons(){return num_neurons;}
  
  inline void set_weightDecay(double lambda) {weightDecay = lambda;}
  inline double get_weightDecay(double lambda) {return weightDecay;}

  void write_to_file(ofstream &outfile);
  void read_from_file(ifstream &infile);
 
  void classify(State& down, State &up);
  void adjust_weights(State &down, State &up);
  void adjust_weights1(State &down, State &up);


  
 protected:
  double mu;
  int num_neurons;
  int num_features;
  int len_w;
  double **w;

  double weightDecay;
};


class LinearLoss
{
 public:
  LinearLoss() {}
  ~LinearLoss() {}

  double classify(State& down, int label);
  void adjust_weights(State &down, int label);
   
};

class SigmoidLoss
{
public:
  SigmoidLoss() {}
  ~SigmoidLoss() {}
  
  double classify(State& down, int label);
  void adjust_weights(State &down, int label);
 
 protected:
  int num_neurons;
  double *a;
};



class SquareLoss
{
 public:
  SquareLoss() {}
  ~SquareLoss() {}
  
  double classify(State& down, int label);
  void adjust_weights(State &down, int label);
   
};


 
class NeuralNet
{
 public:
  NeuralNet() : lin_flag(0){}
  ~NeuralNet() {}
   /**********************************************************************************************
  * nfeatures - number of features in the input vector  
  * num_hidden_units - number of hidden units in the first layer
  * mu - learning rate
  * cost_flag - if 0 - cost1(sigmoid) loss fxn is used, if 1 - cost2(linear) loss function is used
  * lflag - linear flag, if 1 - the net is linear, if 0 - the net is nonlinear, i.e. has hidden units
  * bias - if 1 -  the last linear module has bias, if 0 - last linear module does not have bias
  ************************************************************************************************/
  void initialize(int nfeatures, int num_hidden_units, double mu, int cost_flag, int lflag, int bias);
  void resize_states();
  void remove_bias();
  void clear();

  NeuralNet& operator=(NeuralNet& N);
  void write_to_file(ofstream &outfile);
  void read_from_file(ifstream &infile);
  
  double classify(const double *x);
  double get_err(int label);
  double train(const double *x, int label);
  double train1(const double *x, int label);

  inline void set_cost_flag(int clf){cost_lin_flag = clf;}
  inline double get_mu(){return lin1.get_mu();}
  inline void set_mu(double m) {lin1.set_mu(m); lin2.set_mu(m);}
  
  inline void set_weightDecay(double lambda) {lin1.set_weightDecay(lambda); lin2.set_weightDecay(lambda);}

 protected:
  int lin_flag;
  int cost_lin_flag;

  State start;
  Linear lin1;
  
  State s1;
  Sigmoid sigm1;
  State s2;
  Linear lin2;
  
  State finish;
 
  LinearLoss cost2;
  SigmoidLoss cost1;
    
};
 
#endif /*NET_H_*/
