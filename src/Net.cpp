#include <iostream>
#include <fstream>
#include <cmath>
#include <assert.h>
using namespace std;
#include "Net.h"


/******************** Linear *************************/

void Linear :: resize(int m, int n)
{
  if(num_neurons != m || len_w != n)
    {
      clear();
      num_neurons = m;
      len_w = n;
      w = new double*[num_neurons];
      for(int k = 0; k < num_neurons; k++)
	w[k] = new double[len_w];
    }
  for(int k = 0; k < num_neurons; k++)
    for(int j = 0; j < len_w; j++)
      w[k][j] = ((double)rand()/RAND_MAX - 0.5)/num_features;
  
  
}
 
void Linear :: clear()
{
  if(num_neurons != 0)
    { 
      for(int k = 0; k < num_neurons; k++)
	if(len_w != 0)
	  delete [] w[k];
      delete [] w;
    }
  num_neurons = 0;
  len_w = 0;
}

Linear& Linear :: operator=(Linear &L)
{
  //set the mu
  mu = L.mu;
  num_features = L.num_features;
  resize(L.num_neurons, L.len_w);
  for(int k = 0; k < num_neurons; k++)
    for(int j = 0; j < len_w; j++)
      w[k][j] = L.w[k][j];
  return *this;
}
 

 
void Linear :: write_to_file(ofstream &outfile)
{
  //white out the dimensions of the linear module
  outfile << num_neurons << " " << len_w << " " << num_features << " " << mu << "\n";
  for(int k = 0; k < num_neurons; k++)
    {
      for(int j = 0; j < len_w; j++)
	outfile << w[k][j] << " ";
      outfile << "\n";
    }
 
}
 
void Linear :: read_from_file(ifstream &infile)
{
  int m,n;
  infile >> m;
  infile >> n;
  infile >> num_features;
  infile >> mu;
  resize(m,n);
  for(int k = 0; k < num_neurons; k++)
    {
      for(int j = 0; j < len_w; j++)
	infile >> w[k][j];
    }
}
 
 
void Linear :: classify(State& down, State& up)
{
  double dist = 0;
  for(int k = 0; k < num_neurons; k++)
    {
      int j;
      dist = 0;
      for(j = 0; j < num_features; j++)
	  dist += w[k][j]*down.x[j];
      //if there is a bias
      if(len_w != num_features)
	{
	  assert(len_w == num_features+1);
	  dist += w[k][j];
	}
      up.x[k] = dist;
    }
 }
 
 
void Linear :: adjust_weights(State &down, State &up)
{

  for(int j = 0; j < num_features; j++)
    down.dx[j] = 0.0;
  for(int j = 0; j < num_features; j++)
    for(int k = 0; k < num_neurons; k++)
      down.dx[j] += up.dx[k]*w[k][j];
  
  for(int k = 0; k < num_neurons; k++)
    {
      int j;
      for(j = 0; j < num_features; j++)
	{
	  w[k][j] -= mu*up.dx[k]*down.x[j];
	  if(weightDecay > 0)
	    w[k][j] -= weightDecay*w[k][j];
	}
      
      //if there is a bias
      if(len_w != num_features)
	{
	  assert(len_w == num_features+1);
	  w[k][j] -= mu*up.dx[k];
	  if(weightDecay > 0)
	    w[k][j] -= weightDecay*w[k][j];
	}
    }
  
}



 
void Linear :: adjust_weights1(State &down, State &up)
{

  for(int j = 0; j < num_features; j++)
    down.dx[j] = 0.0;
  for(int j = 0; j < num_features; j++)
    for(int k = 0; k < num_neurons; k++)
      down.dx[j] += up.dx[k]*w[k][j];
  
  for(int k = 0; k < num_neurons; k++)
    {
      int j;
      for(j = 0; j < num_features; j++)
	{
	  w[k][j] -= (mu/(len_w*num_features))*up.dx[k]*down.x[j];
	  if(weightDecay > 0)
	    w[k][j] -= weightDecay*w[k][j];
	}
      
      //if there is a bias
      if(len_w != num_features)
	{
	  assert(len_w == num_features+1);
	  w[k][j] -= (mu/(len_w*num_features))*up.dx[k];
	  if(weightDecay > 0)
	    w[k][j] -= weightDecay*w[k][j];
	}
    }
  
}

/******************************* Sigmoid *******************************/

void Sigmoid :: resize(int m)
{
  if(num_neurons != m)
    {
      clear();
      num_neurons = m;
      a = new double[num_neurons];
    }
  for(int k = 0; k < num_neurons; k++)
    a[k] = 0;
}
 

Sigmoid& Sigmoid :: operator=(Sigmoid &S)
{
  resize(S.num_neurons);
  return *this;
}
 
void Sigmoid :: write_to_file(ofstream &outfile)
{
  outfile << num_neurons << "\n";
}
 

void Sigmoid :: read_from_file(ifstream &infile)
{
  int m;
  infile >> m;
  resize(m);
}


/*
void Sigmoid :: classify(State& down, State &up)
{
  for(int k = 0; k < num_neurons; k++)
    {
      a[k] = exp(down.x[k]);
      up.x[k] = 1/(1+a[k]);
    }
}
 
void Sigmoid :: adjust_weights(State &down, State &up)
{
  for(int k = 0; k < num_neurons; k++)
    down.dx[k] = up.dx[k]*(-a[k]/((1+a[k])*(1+a[k])));
}
*/

void Sigmoid :: classify(State& down, State &up)
{
  for(int k = 0; k < num_neurons; k++)
    {
      a[k] = exp(-down.x[k]);
      up.x[k] = 1/(1+a[k]);
    }
}
 
void Sigmoid :: adjust_weights(State &down, State &up)
{
  for(int k = 0; k < num_neurons; k++)
    down.dx[k] = up.dx[k]*(1-up.x[k])*up.x[k];
}

/******************** Linear Loss***********************************/

double LinearLoss :: classify(State &down, int label)
{
  return -label*down.x[0];
}

void LinearLoss :: adjust_weights(State &down, int label)
{
  down.dx[0] = -label;
}



/***************** Sigmoid Loss *************************/

double SigmoidLoss :: classify(State& down, int label)
{
  double a = exp(label*down.x[0]);
  return  (1/(1+a));
}
 
 

void SigmoidLoss :: adjust_weights(State &down, int label)
{
  double a = exp(label*down.x[0]);
  down.dx[0] = -a/((1+a)*(1+a))*label;
}



/******************** Square Loss***********************************/

double SquareLoss :: classify(State &down, int label)
{
  return (label-down.x[0])*(label-down.x[0])/2;
}

void SquareLoss :: adjust_weights(State &down, int label)
{
  down.dx[0] = -1*(label-down.x[0]);
}



 

 
/****************** NN **************************/
 
void NeuralNet :: resize_states()
{
  int num_features = lin1.getNumFeatures();
  int num_neurons1 = lin1.getNumNeurons();
  int num_neurons2 = 0;
 
  //first state
  start.resize(num_features);
 
  if(lin_flag)
    {
      assert(num_neurons1 == 1);
      finish.resize(num_neurons1);
    }
  else
    {
      assert(num_neurons1 == sigm1.getNumNeurons());
      s1.resize(num_neurons1);
      assert(num_neurons1 == lin2.getNumFeatures());
      s2.resize(num_neurons1);
 
      num_neurons2 = lin2.getNumNeurons();
      assert(num_neurons2 == 1);
      finish.resize(num_neurons2);
    }
}
 


void NeuralNet :: initialize(int nfeatures, int num_hidden_units, double mu, int cost_flag, int lflag, int bias)
{
  lin_flag = lflag;
  cost_lin_flag = cost_flag;

  cerr << "NN " << num_hidden_units <<" " << lin_flag << " " << nfeatures << "\n";

  if(lin_flag == 1)
    {
      int num_neurons1 = 1;
      lin1.set_mu(mu);
      lin1.set_num_feat(nfeatures);
      if(bias)
	lin1.resize(num_neurons1, nfeatures+1);
      else
	lin1.resize(num_neurons1, nfeatures);
    }
  else
    {
      int num_neurons1 = num_hidden_units;
     
      lin1.set_mu(mu);
      lin1.set_num_feat(nfeatures);
      if(bias)
	lin1.resize(num_neurons1, nfeatures+1);
      else
	lin1.resize(num_neurons1, nfeatures);
           
      sigm1.resize(num_neurons1);
         
      int num_neurons2 = 1;
      
      lin2.set_mu(mu);
      lin2.set_num_feat(num_neurons1);
      if(bias)
	lin2.resize(num_neurons2, num_neurons1+1);
      else
	lin2.resize(num_neurons2, num_neurons1);
    }

  resize_states();
}
 

void NeuralNet :: remove_bias()
{
  if(lin_flag)
    lin1.set_len_w(lin1.getNumFeatures());
  else
    lin2.set_len_w(lin2.getNumFeatures());
}


void NeuralNet::clear()
{
  lin_flag = 0;
  cost_lin_flag = 0;

  start.clear();
  lin1.clear();
  
  s1.clear();
  sigm1.clear();
  s2.clear();
  lin2.clear();
  
  finish.clear();
 
}

NeuralNet& NeuralNet :: operator=(NeuralNet& N)
{
  lin_flag = N.lin_flag;
  cost_lin_flag = N.cost_lin_flag;
  if(lin_flag)
    lin1 = N.lin1;
  else
    {
      lin1 = N.lin1;
      sigm1 = N.sigm1;
      lin2 = N.lin2;
    }
  resize_states();
  return *this;
}
 
 
 
void NeuralNet :: write_to_file(ofstream &outfile)
{
 
  //write out the vars
  outfile << "LinFlag\n";
  outfile << lin_flag << "\n";
 
  if(lin_flag)
    {
      outfile << "Linear\n";
      lin1.write_to_file(outfile);
 
      outfile << "Cost\n";
      outfile << cost_lin_flag << "\n";
    }
  else
    {
      outfile << " Linear\n";
      lin1.write_to_file(outfile);
 
      outfile << "Sigmoid\n";
      sigm1.write_to_file(outfile);
 
      outfile << " Linear\n";
      lin2.write_to_file(outfile);
 
      outfile << "Cost\n";
      outfile << cost_lin_flag << "\n";
    }
}
 

void NeuralNet :: read_from_file(ifstream &infile)
{
  string in_str;
 
  infile >> in_str;
  infile >> lin_flag;
 
  if(lin_flag)
    {
      //read the linear module
      infile >> in_str;
      lin1.read_from_file(infile);
 
      infile >> in_str;
      infile >> cost_lin_flag;

    }
  else
    {
      infile >> in_str;
      lin1.read_from_file(infile);
      infile >> in_str;
      sigm1.read_from_file(infile);
 
     
      infile >> in_str;
      lin2.read_from_file(infile);
      infile >> in_str;
      infile >> cost_lin_flag;
    }
  resize_states();
}
 
double NeuralNet::classify(const double *x)
{
  int num_features = lin1.getNumFeatures();
  
  for(int j = 0; j < num_features; j++)
    start.x[j] = x[j];
  
  double r = 0.0;
  if(lin_flag == 1)
    {
      //pass throug the net
      lin1.classify(start,finish);
      
      //get the result of f(x)
      r = finish.x[0];
    }
  else
    {
      //pass throug the net
      lin1.classify(start,s1);
      sigm1.classify(s1,s2);
      lin2.classify(s2,finish);
 
      //get the result
      r = finish.x[0];
    }
 
  return r;
}
 

double NeuralNet :: get_err(int label) 
{
  if(cost_lin_flag)
    return cost2.classify(finish,label);
  else
    return cost1.classify(finish,label);
}


double NeuralNet:: train(const double *x, int label)
{
  int num_features = lin1.getNumFeatures();
  for(int j = 0; j < num_features; j++)
    start.x[j] = x[j];
 
  double r = 0.0;
  if(lin_flag == 1)
    {
      //pass throug the net
      lin1.classify(start,finish);
      r = finish.x[0];
      //pass through the cost and get the result
      if(cost_lin_flag)
	cost2.classify(finish,label);
      else
	cost1.classify(finish,label);
  
      //adjust weights
      if(cost_lin_flag)
	cost2.adjust_weights(finish, label);
      else
	cost1.adjust_weights(finish, label);
      lin1.adjust_weights(start, finish);
    }
  else
    {
      //pass throug the net
      lin1.classify(start,s1);
      sigm1.classify(s1,s2);
      lin2.classify(s2,finish);
      r = finish.x[0];
      //pass through the cost
      if(cost_lin_flag)
	cost2.classify(finish,label);
      else
	cost1.classify(finish,label);
            
      //adjust weights
      if(cost_lin_flag)
	cost2.adjust_weights(finish, label);
      else
	cost1.adjust_weights(finish, label);
      
      lin2.adjust_weights(s2,finish);
      sigm1.adjust_weights(s1,s2);
      lin1.adjust_weights(start, s1);
    }
 
  return r;
}



double NeuralNet:: train1(const double *x, int label)
{
  int num_features = lin1.getNumFeatures();
  for(int j = 0; j < num_features; j++)
    start.x[j] = x[j];
 
  double r = 0.0;
  if(lin_flag == 1)
    {
      //pass throug the net
      lin1.classify(start,finish);
      r = finish.x[0];
      //pass through the cost and get the result
      if(cost_lin_flag)
	cost2.classify(finish,label);
      else
	cost1.classify(finish,label);
  
      //adjust weights
      if(cost_lin_flag)
	cost2.adjust_weights(finish, label);
      else
	cost1.adjust_weights(finish, label);
      lin1.adjust_weights1(start, finish);
    }
  else
    {
      //pass throug the net
      lin1.classify(start,s1);
      sigm1.classify(s1,s2);
      lin2.classify(s2,finish);
      r = finish.x[0];
      //pass through the cost
      if(cost_lin_flag)
	cost2.classify(finish,label);
      else
	cost1.classify(finish,label);
            
      //adjust weights
      if(cost_lin_flag)
	cost2.adjust_weights(finish, label);
      else
	cost1.adjust_weights(finish, label);
      
      lin2.adjust_weights1(s2,finish);
      sigm1.adjust_weights(s1,s2);
      lin1.adjust_weights1(start, s1);
    }
 
  return r;
}







