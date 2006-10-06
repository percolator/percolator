#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "ssl.h"

void exit_with_help()
{
	printf(
	" Copyright: Vikas Sindhwani [vikass@cs.uchicago.edu] \n\n"
	" This software is freely available for academic or commercial use. It must not\n"
        " be modified and distributed without prior permission of the author.\n"
        " The author is not responsible for implications from the use of this software\n"
	"\n Usage\n"
	"     Train: svmlin [options] training_examples training_labels\n"
	"      Test: svmlin -w weights_file test_examples\n"
	"  Evaluate: svmlin -w weights_file test_examples test_labels\n"
	"\n Options:\n"
	" -A algorithm : set algorithm (default 1)\n"
        "        0 -- Regularized Least Squares Classification (RLSC)\n" 
	"	1 -- SVM (L2-SVM-MFN) (default choice)\n"
	"	2 -- Multi-switch Transductive SVM (using L2-SVM-MFN)\n"
        "        3 -- Deterministic Annealing Semi-supervised SVM (using L2-SVM-MFN)\n"
	" -W regularization parameter lambda (default 1)\n" 
        " -U regularization parameter lambda_u (default 1)\n"
	" -S maximum number of switches in TSVM (default 10000)\n"
	" -R positive class fraction of unlabeled data  (default 0.5)\n"
	" -f weights_filename (Test Mode: input_filename is the test file) \n"
	" -Cp relative cost for positive examples (only available with -A 1)\n"
	" -Cn relative cost for positive examples (only available with -A 1)\n"
	" \n References:\n [1] Vikas Sindhwani and S. Sathiya Keerthi.\n     Large Scale Semi-supervised Linear SVMs.\n     Proceedings of ACM SIGIR, 2006\n"
	" [2] V. Sindhwani and S. Sathiya Keerthi.\n     Newton Methods for Fast Solution of Semi-supervised Linear SVMs.\n     Book Chapter in Large Scale Kernel Machines, MIT Press, 2006\n"
	" "
	" \n"
        );
	exit(1);
}

void parse_command_line(int argc, char **argv);
void Read(char *input_file_name, struct data* DATA, const vector_int *Subset);
void Read(char *weights_file_name, struct vector_double *Weights);
void ReadTrainingData();
void run_splits(char *input_file_name, char *labels_file_name, vector_int *Splits);

using namespace std;
 
char weights_file_name[1024];
char outputs_file_name[1024];
char inputs_file_name[1024];
char labels_file_name[1024];
char splits_file_name[1024];

struct options *Options = new options[1];
struct data *Data = new data[1];
struct vector_double *Weights = new vector_double[1];
struct vector_double *Outputs = new vector_double[1];
struct vector_double *Labels = new vector_double[1];

bool train = true; /* mode -- 1 for training and 0 for testing */
bool evaluate = true; /* in test mode, whether to evaluate or not, set to zero if no test labels are given */

int mainy(int argc, char **argv)
{
  cout << endl << " svmlin (v1.0)   \n\n"; 
  parse_command_line(argc, argv);
  if(train) /* train */
    {
	  ReadTrainingData();
	  ssl_train(Data,Options,Weights,Outputs);
	  cout << endl << "Writing Weights to file: " << weights_file_name << endl;
	  Write(weights_file_name,Weights); 
	  Clear(Data);
    }
  else
    {  /* test */
      Read(weights_file_name,Weights);
      ssl_predict(inputs_file_name,Weights,Outputs); 
      if(evaluate) /* evaluate */
	{
	  Read(labels_file_name,Labels);
	  ssl_evaluate(Outputs,Labels);
	}
    }
  
  cout << "Writing Outputs to file: " << outputs_file_name << endl;
  Write(outputs_file_name,Outputs);
  
  Clear(Labels);
  Clear(Weights);
  Clear(Outputs);
  
  return 0;
}

void parse_command_line(int argc, char **argv)
{
  int i;
  Options->algo = 1;
  Options->lambda=1.0;
  Options->lambda_u=1.0;
  Options->S=10000;
  Options->R=0.5;
  Options->epsilon=EPSILON;
  Options->cgitermax=CGITERMAX;
  Options->mfnitermax=MFNITERMAX;
  Options->Cp = 1.0;
  Options->Cn = 1.0;
  // parse options
  for(i=1;i<argc;i++)
    {
      if(argv[i][0] != '-') break;
      if(++i>=argc)
	exit_with_help();
      switch(argv[i-1][1])
	{
	case 'A':
	  Options->algo = atoi(argv[i]);
	  break;
	case 'W':
	  Options->lambda = atof(argv[i]);
	  break;
	case 'U':
	  Options->lambda_u = atof(argv[i]);
	  break;
	case 'S':
	  Options->S = atoi(argv[i]);
	  break;
	case 'R':
	  Options->R = atof(argv[i]); 
	  break;
	case 'C':
	  switch(argv[i-1][2])
	    {
	    case 'p':
	      Options->Cp = atof(argv[i]);
	      break;
	    case 'n':
	      Options->Cn = atof(argv[i]);
	      break;
	    default:
	      fprintf(stderr,"unknown option -- Specify Cp and/or Cn.\n");
	      exit_with_help();
	    }
	  break;
	case 'f':
	  train=0;
	  strcpy(weights_file_name, argv[i]);
	  break;
	default:
	  fprintf(stderr,"unknown option\n");
	  exit_with_help();
	}
    }
  // determine filenames
  
  if(i>=argc)
    exit_with_help();
  
  strcpy(inputs_file_name, argv[i]);
  char *p = strrchr(argv[i],'/');
  if(p==NULL)
    p = argv[i];
  else
    ++p;
  if(train)
    sprintf(weights_file_name,"%s.weights",p);  
  sprintf(outputs_file_name,"%s.outputs",p);  
  ++i;
  if(i>=argc)
    {
      if(train)
	exit_with_help(); // need labels
      else
	evaluate=0; // test mode, no test labels given, cannot evaluate
    }	
  else
    strcpy(labels_file_name, argv[i]);  
}

/* Read a subset of input patterns from input_file_name specified by
 the index Subset and set the val,rowptr,colind fields in DATA */
void Read(char *input_file_name, struct data* DATA, const vector_int *Subset)
{
  cout << "Reading file: " << inputs_file_name  << endl;
  timer tictoc;
  tictoc.restart();
  FILE *fpin;
  fpin = fopen(inputs_file_name, "r");
  if (fpin == NULL)
    {
      std::cerr << "Cannot open input file: " << inputs_file_name << "\n";
      exit(1);
    }  
  int m = 0;
  int n = 0;
  int nz = 0;   
  int t=0;
  int read = 0;
  if(Subset->vec[0]==0) {read=1;t++;}
  
  while (1)
    {
      int c = fgetc(fpin);
      switch(c)
	{
	case '\n':
	  ++m; 
	  if(m==Subset->vec[t]) {read=1; t++;} else read=0;
	  break;
	case ':': 
	  if(read) ++nz;
	  break;
	case EOF:
	  goto out;
	default:
	  ;
	}
    }
 out:  
  rewind(fpin);
  t=0;
  nz += Subset->d ;/* create space to hold bias feature for each example */
  double* VAL = new double[nz];
  int* C = new int[nz];
  int* R = new int[Subset->d+1];
  
  int j=0;
  
  for(int i=0;i<m;i++)
    {
      if(Subset->vec[t]==i)
	{	 
	  read=1;
	  R[t]=j;
	  t++;
	  while(1)
	    {
	      int c;
	      do {
		c = getc(fpin);
		if(c=='\n') goto out2;
	      } while(isspace(c));
	      ungetc(c,fpin);
	      fscanf(fpin,"%d:%lf",&C[j],&VAL[j]);
	      C[j]=C[j]-1;                         
	       if(C[j]>n) n=C[j];                       
	       ++j;
	       if(j>nz) break;
	    }	  	  
	}
      else 
	{
	  read=0;
	  while(1)
	    {
	      int c;
	      do {
		c = getc(fpin);
	      } while(c!='\n');
	      goto out2;
	    }
	}
      
    out2:
      if(read)
	{
	  C[j]=0; /* placeholder for the constant bias feature */
	  VAL[j]=1.0; /* for the constant bias feature */
	  ++j; /* for the constant bias feature */
	}
    }
  R[Subset->d]=nz;
  /* record the column index for the bias feature */
  n++;
  for(int i=0;i<Subset->d;i++)
    C[R[i+1]-1]=n;
  n++;
  Data->n=n;
  Data->m=Subset->d;
  Data->val=VAL;
  Data->rowptr=R;
  Data->colind=C;
  Data->nz = nz;

  cout << "  Input Data Matrix Statistics:" << endl;
  cout << "      Examples: " << Subset->d << endl;
  cout << "      Features: " << n << " (including bias feature)" << endl; 
  cout << "      Non-zeros:  " << nz << " (including bias features)" << endl;
  cout << "      Average sparsity: " << nz*1.0/Subset->d << " non-zero features per example." << endl;
  fclose(fpin);
  tictoc.stop();
  cout <<  "Reading took " << tictoc.time() << " secs. \n" << endl; 
  return;
}      

 
void Read(char *weights_file_name, struct vector_double *Weights)
{
  cout << "Reading file: " << weights_file_name << endl;
  int n=0;
  FILE *fpwts;
  fpwts = fopen(weights_file_name,"r");

  while (1)
    {
      int c = fgetc(fpwts);
      switch(c)
	{
	case '\n':
	  ++n; 
	  break;
	case EOF:
	  goto out;
	default:
	  ;
	}
    }
 out:
  rewind(fpwts);
  Weights->vec = new double[n]; 
  Weights->d=n;
  for(int i=0;i<n;i++)
      fscanf(fpwts,"%lf\n",&Weights->vec[i]);
  return;
}


void ReadTrainingData()
{
  struct vector_int *Subset = new vector_int[1];
  Read(labels_file_name,Labels);
  int pos=0;
  int neg=0;
  Data->l=0;
  Data->u=0;
  for(int i=0;i<Labels->d;i++)
    if(Labels->vec[i]!=0.0 || Options->algo==-1) {
      Data->l++;
      if(Labels->vec[i] > 0) pos++; else neg++;
    }
  if(Options->algo > 1) /* semi-supervised method */
    {
      initialize(Subset,Labels->d);
      Data->m = Labels->d;
      Data->u = Data->m - Data->l;
      Data->Y = new double[Data->m];
      Data->C = new double[Data->m];
      for(int i=0;i<Data->m;i++)
	{
	  Data->C[i] = 1.0;
	  Data->Y[i] = Labels->vec[i];
	}
    }
  else /* supervised method, ignore unlabeled data */
    {
      Data->u=0; 
      Subset->vec = new int[Data->l];
      Subset->d = Data->l;
      Data->Y = new double[Data->l];
      Data->C = new double[Data->l];
      int t=0;
      for(int i=0;i<Labels->d;i++)
	{
	  if(Labels->vec[i]!=0.0 || Options->algo==-1) 
	    {
	      Subset->vec[t]=i; 
	      Data->Y[t]=Labels->vec[i];
	      if(Labels->vec[i]>0) 
		Data->C[t]=Options->Cp; 
	      else 
		Data->C[t]=Options->Cn;
	      t++;
	    }
	}
    }
  
  cout << "  Label Vector Statistics:\n" << "      Labeled = " << Data->l 
       << " (Positive = " << pos << "  Negative = " << neg << ")" << endl 
       << "      Unlabeled = " << Data->u << endl;
  Read(inputs_file_name,Data,Subset);
  Clear(Subset);
}

