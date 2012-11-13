/*******************************************************************************
 Copyright 2006-2010 Lukas Käll <lukas.kall@cbr.su.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/
/*
 * @ Created by Luminita Moruz
 * Sep, 2010
 */
/* This file includes the implementations for the functions defined in LibsvmWrapper.h */
#include <stdlib.h>
#include <string.h>

#include "LibsvmWrapper.h"
#include "Globals.h"
#include "PSMDescription.h"
#include "svm.h"

svm_model* libsvm_wrapper::TrainModel(const std::vector<PSMDescription> &psms, const int &number_features, const svm_parameter &parameter) {
  svm_model *svr_model;
  int number_examples = psms.size();
  svm_problem data;
  data.l = number_examples;
  data.x = new svm_node[number_examples];
  data.y = new double[number_examples];
  for (int i = 0; i < number_examples; i++) {
    data.x[i].values = psms[i].retentionFeatures;
    data.x[i].dim = number_features;
    data.y[i] = psms[i].retentionTime;
  }
  // build a model by training the SVM on the given training set
  char const *error_message = svm_check_parameter(&data, &parameter);
  if (error_message != NULL) {
    delete[] data.x;
    delete[] data.y;
    ostringstream temp;
    temp << "Error : Incorrect parameters for the SVR. Execution aborted. " << endl;
    throw MyException(temp.str());
  }
  svr_model = svm_train(&data, &parameter);
  delete[] data.x;
  delete[] data.y;
  return svr_model;
}

double libsvm_wrapper::PredictRT(const svm_model* svr, const int &number_features, double *features) {
  svm_node node;
  node.values = features;
  node.dim = number_features;

  return svm_predict(svr, &node);
}

int libsvm_wrapper::SaveModel(FILE* fp, const svm_model* model) {
  /*FILE* fp = fopen(model_file_name, "w");
   if (fp == NULL) {
     return -1;
   }*/
   const char* svm_type_table_[] = { "c_svc", "nu_svc", "one_class",
                                   "epsilon_svr", "nu_svr", NULL };

   const char* kernel_type_table_[] = { "linear", "polynomial", "rbf",
                                      "sigmoid", "precomputed", NULL };
   const svm_parameter& param = model->param;
   fprintf(fp, "svm_type %s\n", svm_type_table_[param.svm_type]);
   fprintf(fp, "kernel_type %s\n", kernel_type_table_[param.kernel_type]);
   if (param.kernel_type == POLY) {
     fprintf(fp, "degree %d\n", param.degree);
   }
   if (param.kernel_type == POLY || param.kernel_type == RBF
       || param.kernel_type == SIGMOID) {
     fprintf(fp, "gamma %g\n", param.gamma);
   }
   if (param.kernel_type == POLY || param.kernel_type == SIGMOID) {
     fprintf(fp, "coef0 %g\n", param.coef0);
   }
   int nr_class = model->nr_class;
   int l = model->l;
   fprintf(fp, "nr_class %d\n", nr_class);
   fprintf(fp, "total_sv %d\n", l);
   {
     fprintf(fp, "rho");
     for (int i = 0; i < nr_class * (nr_class - 1) / 2; i++) {
       fprintf(fp, " %g", model->rho[i]);
     }
     fprintf(fp, "\n");
   }
   if (model->label) {
     fprintf(fp, "label");
     for (int i = 0; i < nr_class; i++) {
       fprintf(fp, " %d", model->label[i]);
     }
     fprintf(fp, "\n");
   }
   if (model->probA) { // regression has probA only
     fprintf(fp, "probA");
     for (int i = 0; i < nr_class * (nr_class - 1) / 2; i++) {
       fprintf(fp, " %g", model->probA[i]);
     }
     fprintf(fp, "\n");
   }
   if (model->probB) {
     fprintf(fp, "probB");
     for (int i = 0; i < nr_class * (nr_class - 1) / 2; i++) {
       fprintf(fp, " %g", model->probB[i]);
     }
     fprintf(fp, "\n");
   }
   if (model->nSV) {
     fprintf(fp, "nr_sv");
     for (int i = 0; i < nr_class; i++) {
       fprintf(fp, " %d", model->nSV[i]);
     }
     fprintf(fp, "\n");
   }
   fprintf(fp, "SV\n");
   const double* const * sv_coef = model->sv_coef;
 #ifdef _DENSE_REP
   const svm_node* SV = model->SV;
 #else
   const svm_node* const* SV = model->SV;
 #endif
   for (int i = 0; i < l; i++) {
     for (int j = 0; j < nr_class - 1; j++) {
       fprintf(fp, "%.16g ", sv_coef[j][i]);
     }
 #ifdef _DENSE_REP
     const svm_node* p = (SV + i);
     if (param.kernel_type == PRECOMPUTED) {
       fprintf(fp, "0:%d ", (int)(p->values[0]));
     } else
       for (int j = 0; j < p->dim; j++)
         if (p->values[j] != 0.0) {
           fprintf(fp, "%d:%.8g ", j, p->values[j]);
         }
 #else
     const svm_node* p = SV[i];
     if (param.kernel_type == PRECOMPUTED) {
       fprintf(fp, "0:%d ", (int)(p->value));
     } else
     while (p->index != -1) {
       fprintf(fp, "%d:%.8g ", p->index, p->value);
       p++;
     }
 #endif
     fprintf(fp, "\n");
   }
   if (ferror(fp) != 0) {
     return -1;
   } else {
     return 0;
   }
}

svm_model* libsvm_wrapper::LoadModel(FILE* fp) {
  const char* svm_type_table[] = { "c_svc", "nu_svc", "one_class",
                                   "epsilon_svr", "nu_svr", NULL };

  const char* kernel_type_table[] = { "linear", "polynomial", "rbf",
                                      "sigmoid", "precomputed", NULL };
  int ret;
  // read parameters
  svm_model* model = (svm_model*) malloc(sizeof(svm_model));
  svm_parameter& param = model->param;
  model->rho = NULL;
  model->probA = NULL;
  model->probB = NULL;
  model->label = NULL;
  model->nSV = NULL;
  char cmd[81];
  while (1) {
    ret = fscanf(fp, "%80s", cmd);
    if (strcmp(cmd, "svm_type") == 0) {
      ret = fscanf(fp, "%80s", cmd);
      int i;
      for (i = 0; svm_type_table[i]; i++) {
        if (strcmp(svm_type_table[i], cmd) == 0) {
          param.svm_type = i;
          break;
        }
      }
      if (svm_type_table[i] == NULL) {
        fprintf(stderr, "unknown svm type.\n");
        free(model->rho);
        free(model->label);
        free(model->nSV);
        free(model);
        return NULL;
      }
    } else if (strcmp(cmd, "kernel_type") == 0) {
      ret = fscanf(fp, "%80s", cmd);
      int i;
      for (i = 0; kernel_type_table[i]; i++) {
        if (strcmp(kernel_type_table[i], cmd) == 0) {
          param.kernel_type = i;
          break;
        }
      }
      if (kernel_type_table[i] == NULL) {
        fprintf(stderr, "unknown kernel function.\n");
        free(model->rho);
        free(model->label);
        free(model->nSV);
        free(model);
        return NULL;
      }
    } else if (strcmp(cmd, "degree") == 0) {
      ret = fscanf(fp, "%d", &param.degree);
    } else if (strcmp(cmd, "gamma") == 0) {
      ret = fscanf(fp, "%lf", &param.gamma);
    } else if (strcmp(cmd, "coef0") == 0) {
      ret = fscanf(fp, "%lf", &param.coef0);
    } else if (strcmp(cmd, "nr_class") == 0) {
      ret = fscanf(fp, "%d", &model->nr_class);
    } else if (strcmp(cmd, "total_sv") == 0) {
      ret = fscanf(fp, "%d", &model->l);
    } else if (strcmp(cmd, "rho") == 0) {
      int n = model->nr_class * (model->nr_class - 1) / 2;
      model->rho = (double *)malloc(n*sizeof(double));
      for (int i = 0; i < n; i++) {
        ret = fscanf(fp, "%lf", &model->rho[i]);
      }
    } else if (strcmp(cmd, "label") == 0) {
      int n = model->nr_class;
      model->label = (int *)malloc(n*sizeof(int));
      for (int i = 0; i < n; i++) {
        ret = fscanf(fp, "%d", &model->label[i]);
      }
    } else if (strcmp(cmd, "probA") == 0) {
      int n = model->nr_class * (model->nr_class - 1) / 2;
      model->probA = (double *)malloc(n*sizeof(double));
      for (int i = 0; i < n; i++) {
        ret = fscanf(fp, "%lf", &model->probA[i]);
      }
    } else if (strcmp(cmd, "probB") == 0) {
      int n = model->nr_class * (model->nr_class - 1) / 2;
      model->probB = (double *)malloc(n*sizeof(double));
      for (int i = 0; i < n; i++) {
        ret = fscanf(fp, "%lf", &model->probB[i]);
      }
    } else if (strcmp(cmd, "nr_sv") == 0) {
      int n = model->nr_class;
      model->nSV = (int *)malloc(n*sizeof(int));
      for (int i = 0; i < n; i++) {
        ret = fscanf(fp, "%d", &model->nSV[i]);
      }
    } else if (strcmp(cmd, "SV") == 0) {
      while (1) {
        int c = getc(fp);
        if (c == EOF || c == '\n') {
          break;
        }
      }
      break;
    } else {
      fprintf(stderr, "unknown text in model file: [%s]\n", cmd);
      free(model->rho);
      free(model->label);
      free(model->nSV);
      free(model);
      return NULL;
    }
  }
  // read sv_coef and SV
  int elements = 0;
  long pos = ftell(fp);
#ifdef _DENSE_REP
  char buffer[128], c;
  double value;
  int index = 0;
  // read the max dimension of all vectors
  while ((c = fgetc(fp)) != EOF) {
    if (isspace(c)) {
      index = 0;
    } else if (c == ':') {
      buffer[index] = '\0';
      // variable 'elements' is used to keep max dimension in dense rep.
      elements = max(elements, atoi(buffer) + 1);
      index = 0;
    } else {
      buffer[index++] = c;
    }
  }
#else
  while (1) {
    int c = fgetc(fp);
    switch (c) {
      case '\n':
      // count the '-1' element
      case ':':
      ++elements;
      break;
      case EOF:
      goto out;
      default:
      ;
    }
  }
  out:
#endif
  fseek(fp, pos, SEEK_SET);
  int m = model->nr_class - 1;
  int l = model->l;
  model->sv_coef =  (double **)malloc(m*sizeof(double*));
  int i;
  for (i = 0; i < m; i++) {
    model->sv_coef[i] = (double *)malloc(l*sizeof(double));
  }
#ifdef _DENSE_REP
  model->SV = (svm_node*) malloc(l*sizeof(svm_node));
  for (i = 0; i < l; i++) {
    model->SV[i].values =  (double *)malloc(elements*sizeof(double));
    model->SV[i].dim = 0;
    for (int k = 0; k < m; k++) {
      ret = fscanf(fp, "%lf", &model->sv_coef[k][i]);
    }
    int* d = &(model->SV[i].dim);
    while ((c = getc(fp)) != '\n') {
      if (!isspace(c)) {
        ungetc(c, fp);
        ret = fscanf(fp, "%d:%lf", &index, &value);
        while (*d < index) {
          model->SV[i].values[(*d)++] = 0.0;
        }
        model->SV[i].values[(*d)++] = value;
      }
    }
  }
#else
  model->SV = Malloc(svm_node*, l);
  svm_node* x_space = NULL;
  if (l > 0) {
    x_space = Malloc(svm_node, elements);
  }
  int j = 0;
  for (i = 0; i < l; i++) {
    model->SV[i] = &x_space[j];
    for (int k = 0; k < m; k++) {
      fscanf(fp, "%lf", &model->sv_coef[k][i]);
    }
    while (1) {
      int c;
      do {
        c = getc(fp);
        if (c == '\n') {
          goto out2;
        }
      }while (isspace(c));
      ungetc(c, fp);
      fscanf(fp, "%d:%lf", &(x_space[j].index), &(x_space[j].value));
      ++j;
    }
    out2:
    x_space[j++].index = -1;
  }
#endif
  if (ferror(fp) != 0) {
    return NULL;
  }
  model->free_sv = 1; // XXX
  return model;
}
