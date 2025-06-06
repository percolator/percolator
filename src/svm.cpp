#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <float.h>
#include <string.h>
#include <stdarg.h>
#include "svm.h"
typedef float Qfloat;
typedef signed char schar;
#ifndef min
template<class T> inline T min(T x, T y) {
  return (x < y) ? x : y;
}
#endif
#ifndef max
template<class T> inline T max(T x, T y) {
  return (x > y) ? x : y;
}
#endif
template<class T> inline void swap(T& x, T& y) {
  T t = x;
  x = y;
  y = t;
}
template<class S, class T> inline void clone(T*& dst, S* src, int n) {
  dst = new T[n];
  memcpy((void*)dst, (void*)src, sizeof(T) * static_cast<std::size_t>(n));
}
inline double powi(double base, int times) {
  double tmp = base, ret = 1.0;
  for (int t = times; t > 0; t /= 2) {
    if (t % 2 == 1) {
      ret *= tmp;
    }
    tmp = tmp * tmp;
  }
  return ret;
}
#define INF HUGE_VAL
#define TAU 1e-12
#define Malloc(type,n) (type *)malloc((static_cast<std::size_t>(n))*sizeof(type))
#if 0
void info(const char* fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  vprintf(fmt, ap);
  va_end(ap);
}
void info_flush() {
  fflush(stdout);
}
#else
void info(const char* fmt, ...) {
}
void info_flush() {
}
#endif

//
// Kernel Cache
//
// l is the number of total data items
// size is the cache size limit in bytes
//
class Cache {
  public:
    Cache(std::size_t l, long int size);
    ~Cache();

    // request data [0,len)
    // return some position p where [p,len) need to be filled
    // (p >= len if nothing needs to be filled)
    int get_data(const int index, Qfloat** data, int len);
    void swap_index(int i, int j); // future_option
  private:
    std::size_t l;
    long int size;
    struct head_t {
        head_t* prev, *next; // a cicular list
        Qfloat* data;
        int len; // data[0,len) is cached in this entry
    };

    head_t* head;
    head_t lru_head;
    void lru_delete(head_t* h);
    void lru_insert(head_t* h);
};

Cache::Cache(std::size_t l_, long int size_) :
  l(l_), size(size_) {
  head = (head_t*)calloc(l, sizeof(head_t)); // initialized to 0
  long QfloatSize = static_cast<long>(sizeof(Qfloat));
  size /= QfloatSize;
  size -= static_cast<long>(l * sizeof(head_t) / static_cast<std::size_t>(QfloatSize));
  size = max(size, 2 * (long int)l); // cache must be large enough for two columns
  lru_head.next = lru_head.prev = &lru_head;
}

Cache::~Cache() {
  for (head_t* h = lru_head.next; h != &lru_head; h = h->next) {
    free(h->data);
  }
  free(head);
}

void Cache::lru_delete(head_t* h) {
  // delete from current location
  h->prev->next = h->next;
  h->next->prev = h->prev;
}

void Cache::lru_insert(head_t* h) {
  // insert to last position
  h->next = &lru_head;
  h->prev = lru_head.prev;
  h->prev->next = h;
  h->next->prev = h;
}

int Cache::get_data(const int index, Qfloat** data, int len) {
  head_t* h = &head[index];
  if (h->len) {
    lru_delete(h);
  }
  int more = len - h->len;
  if (more > 0) {
    // free old space
    while (size < more) {
      head_t* old = lru_head.next;
      lru_delete(old);
      free(old->data);
      size += old->len;
      old->data = 0;
      old->len = 0;
    }
    // allocate new space
    h->data = (Qfloat*)realloc(h->data, sizeof(Qfloat) * static_cast<unsigned long>(len));
    size -= more;
    swap(h->len, len);
  }
  lru_insert(h);
  *data = h->data;
  return len;
}

void Cache::swap_index(int i, int j) {
  if (i == j) {
    return;
  }
  if (head[i].len) {
    lru_delete(&head[i]);
  }
  if (head[j].len) {
    lru_delete(&head[j]);
  }
  swap(head[i].data, head[j].data);
  swap(head[i].len, head[j].len);
  if (head[i].len) {
    lru_insert(&head[i]);
  }
  if (head[j].len) {
    lru_insert(&head[j]);
  }
  if (i > j) {
    swap(i, j);
  }
  for (head_t* h = lru_head.next; h != &lru_head; h = h->next) {
    if (h->len > i) {
      if (h->len > j) {
        swap(h->data[i], h->data[j]);
      } else {
        // give up
        lru_delete(h);
        free(h->data);
        size += h->len;
        h->data = 0;
        h->len = 0;
      }
    }
  }
}

//
// Kernel evaluation
//
// the static method k_function is for doing single kernel evaluation
// the constructor of Kernel prepares to calculate the l*l kernel matrix
// the member function get_Q is for getting one column from the Q Matrix
//
class QMatrix {
  public:
    virtual Qfloat* get_Q(int column, int len) const = 0;
    virtual Qfloat* get_QD() const = 0;
    virtual void swap_index(int i, int j) const = 0;
    virtual ~QMatrix() {
    }
};

class Kernel : public QMatrix {
  public:
#ifdef _DENSE_REP
    Kernel(int l, svm_node* x, const svm_parameter& param);
#else
    Kernel(int l, svm_node* const* x, const svm_parameter& param);
#endif
    virtual ~Kernel();

    static double k_function(const svm_node* x, const svm_node* y,
                             const svm_parameter& param);
    virtual Qfloat* get_Q(int column, int len) const = 0;
    virtual Qfloat* get_QD() const = 0;
    virtual void swap_index(int i, int j) const { // no so const...
      swap(x[i], x[j]);
      if (x_square) {
        swap(x_square[i], x_square[j]);
      }
    }
  protected:

    double(Kernel::*kernel_function)(int i, int j) const;

  private:
#ifdef _DENSE_REP
    svm_node* x;
#else
    const svm_node** x;
#endif
    double* x_square;

    // svm_parameter
    const int kernel_type;
    const int degree;
    const double gamma;
    const double coef0;

    static double dot(const svm_node* px, const svm_node* py);
#ifdef _DENSE_REP
    static double dot(const svm_node& px, const svm_node& py);
#endif

    double kernel_linear(int i, int j) const {
      return dot(x[i], x[j]);
    }
    double kernel_poly(int i, int j) const {
      return powi(gamma * dot(x[i], x[j]) + coef0, degree);
    }
    double kernel_rbf(int i, int j) const {
      return exp(-gamma
          * (x_square[i] + x_square[j] - 2 * dot(x[i], x[j])));
    }
    double kernel_sigmoid(int i, int j) const {
      return tanh(gamma * dot(x[i], x[j]) + coef0);
    }
    double kernel_precomputed(int i, int j) const {
#ifdef _DENSE_REP
      return (x + i)->values[(int)((x + j)->values[0])];
#else
      return x[i][(int)(x[j][0].value)].value;
#endif
    }
};

#ifdef _DENSE_REP
Kernel::Kernel(int l, svm_node* x_, const svm_parameter& param)
#else
Kernel::Kernel(int l, svm_node* const* x_, const svm_parameter& param)
#endif
:
  kernel_type(param.kernel_type), degree(param.degree),
      gamma(param.gamma), coef0(param.coef0) {
  switch (kernel_type) {
    case LINEAR:
      kernel_function = &Kernel::kernel_linear;
      break;
    case POLY:
      kernel_function = &Kernel::kernel_poly;
      break;
    case RBF:
      kernel_function = &Kernel::kernel_rbf;
      break;
    case SIGMOID:
      kernel_function = &Kernel::kernel_sigmoid;
      break;
    case PRECOMPUTED:
      kernel_function = &Kernel::kernel_precomputed;
      break;
  }
  clone(x, x_, l);
  if (kernel_type == RBF) {
    x_square = new double[l];
    for (std::size_t i = 0; static_cast<int>(i) < l; i++) {
      x_square[i] = dot(x[i], x[i]);
    }
  } else {
    x_square = 0;
  }
}

Kernel::~Kernel() {
  delete[] x;
  delete[] x_square;
}

#ifdef _DENSE_REP
double Kernel::dot(const svm_node* px, const svm_node* py) {
  double sum = 0;
  int dim = min(px->dim, py->dim);
  for (std::size_t i = 0; static_cast<int>(i) < dim; i++) {
    sum += (px->values)[i] * (py->values)[i];
  }
  return sum;
}

double Kernel::dot(const svm_node& px, const svm_node& py) {
  double sum = 0;
  int dim = min(px.dim, py.dim);
  for (std::size_t i = 0; static_cast<int>(i) < dim; i++) {
    sum += px.values[i] * py.values[i];
  }
  return sum;
}
#else
double Kernel::dot(const svm_node* px, const svm_node* py) {
  double sum = 0;
  while (px->index != -1 && py->index != -1) {
    if (px->index == py->index) {
      sum += px->value * py->value;
      ++px;
      ++py;
    } else {
      if (px->index > py->index) {
        ++py;
      } else {
        ++px;
      }
    }
  }
  return sum;
}
#endif

double Kernel::k_function(const svm_node* x, const svm_node* y,
                          const svm_parameter& param) {
  switch (param.kernel_type) {
    case LINEAR:
      return dot(x, y);
    case POLY:
      return powi(param.gamma * dot(x, y) + param.coef0, param.degree);
    case RBF: {
      double sum = 0;
#ifdef _DENSE_REP
      int dim = min(x->dim, y->dim), i;
      for (i = 0; i < dim; i++) {
        double d = x->values[i] - y->values[i];
        sum += d * d;
      }
      for (; i < x->dim; i++) {
        sum += x->values[i] * x->values[i];
      }
      for (; i < y->dim; i++) {
        sum += y->values[i] * y->values[i];
      }
#else
      while (x->index != -1 && y->index != -1) {
        if (x->index == y->index) {
          double d = x->value - y->value;
          sum += d * d;
          ++x;
          ++y;
        } else {
          if (x->index > y->index) {
            sum += y->value * y->value;
            ++y;
          } else {
            sum += x->value * x->value;
            ++x;
          }
        }
      }
      while (x->index != -1) {
        sum += x->value * x->value;
        ++x;
      }
      while (y->index != -1) {
        sum += y->value * y->value;
        ++y;
      }
#endif
      return exp(-param.gamma * sum);
    }
    case SIGMOID:
      return tanh(param.gamma * dot(x, y) + param.coef0);
    case PRECOMPUTED: //x: test (validation), y: SV
#ifdef _DENSE_REP
      return x->values[(int)(y->values[0])];
#else
      return x[(int)(y->value)].value;
#endif
    default:
      return 0; // Unreachable
  }
}

// An SMO algorithm in Fan et al., JMLR 6(2005), p. 1889--1918
// Solves:
//
//   min 0.5(\alpha^T Q \alpha) + p^T \alpha
//
//        y^T \alpha = \delta
//        y_i = +1 or -1
//        0 <= alpha_i <= Cp for y_i = 1
//        0 <= alpha_i <= Cn for y_i = -1
//
// Given:
//
//   Q, p, y, Cp, Cn, and an initial feasible point \alpha
//   l is the size of vectors and matrices
//   eps is the stopping tolerance
//
// solution will be put in \alpha, objective value will be put in obj
//
class Solver {
  public:
    Solver() {
    }
    ;
    virtual ~Solver() {
    }
    ;

    struct SolutionInfo {
        double obj;
        double rho;
        double upper_bound_p;
        double upper_bound_n;
        double r; // for Solver_NU
    };

    void Solve(int l, const QMatrix& Q, const double* p_, const schar* y_,
               double* alpha_, double Cp, double Cn, double eps,
               SolutionInfo* si, int shrinking);
  protected:
    int active_size;
    schar* y;
    double* G; // gradient of objective function
    enum {
      LOWER_BOUND, UPPER_BOUND, FREE
    };
    char* alpha_status; // LOWER_BOUND, UPPER_BOUND, FREE
    double* alpha;
    const QMatrix* Q;
    const Qfloat* QD;
    double eps;
    double Cp, Cn;
    double* p;
    int* active_set;
    double* G_bar; // gradient, if we treat free variables as 0
    int l;
    bool unshrinked; // XXX

    double get_C(int i) {
      return (y[i] > 0) ? Cp : Cn;
    }
    void update_alpha_status(int i) {
      if (alpha[i] >= get_C(i)) {
        alpha_status[i] = UPPER_BOUND;
      } else if (alpha[i] <= 0) {
        alpha_status[i] = LOWER_BOUND;
      } else {
        alpha_status[i] = FREE;
      }
    }
    bool is_upper_bound(int i) {
      return alpha_status[i] == UPPER_BOUND;
    }
    bool is_lower_bound(int i) {
      return alpha_status[i] == LOWER_BOUND;
    }
    bool is_free(int i) {
      return alpha_status[i] == FREE;
    }
    void swap_index(int i, int j);
    void reconstruct_gradient();
    virtual int select_working_set(int& i, int& j);
    virtual double calculate_rho();
    virtual void do_shrinking();
  private:
    bool be_shrunken(int i, double Gmax1, double Gmax2);
};

void Solver::swap_index(int i, int j) {
  Q->swap_index(i, j);
  swap(y[i], y[j]);
  swap(G[i], G[j]);
  swap(alpha_status[i], alpha_status[j]);
  swap(alpha[i], alpha[j]);
  swap(p[i], p[j]);
  swap(active_set[i], active_set[j]);
  swap(G_bar[i], G_bar[j]);
}

void Solver::reconstruct_gradient() {
  // reconstruct inactive elements of G from G_bar and free variables
  if (active_size == l) {
    return;
  }
  int i;
  for (i = active_size; i < l; i++) {
    G[i] = G_bar[i] + p[i];
  }
  for (i = 0; i < active_size; i++)
    if (is_free(i)) {
      const Qfloat* Q_i = Q->get_Q(i, l);
      double alpha_i = alpha[i];
      for (int j = active_size; j < l; j++) {
        G[j] += alpha_i * Q_i[j];
      }
    }
}

void Solver::Solve(int l, const QMatrix& Q, const double* p_,
                   const schar* y_, double* alpha_, double Cp, double Cn,
                   double eps, SolutionInfo* si, int shrinking) {
  this->l = l;
  this->Q = &Q;
  QD = Q.get_QD();
  clone(p, p_, l);
  clone(y, y_, l);
  clone(alpha, alpha_, l);
  this->Cp = Cp;
  this->Cn = Cn;
  this->eps = eps;
  unshrinked = false;
  // initialize alpha_status
  {
    alpha_status = new char[l];
    for (int i = 0; i < l; i++) {
      update_alpha_status(i);
    }
  }
  // initialize active set (for shrinking)
  {
    active_set = new int[l];
    for (int i = 0; i < l; i++) {
      active_set[i] = i;
    }
    active_size = l;
  }
  // initialize gradient
  {
    G = new double[l];
    G_bar = new double[l];
    int i;
    for (i = 0; i < l; i++) {
      G[i] = p[i];
      G_bar[i] = 0;
    }
    for (i = 0; i < l; i++)
      if (!is_lower_bound(i)) {
        const Qfloat* Q_i = Q.get_Q(i, l);
        double alpha_i = alpha[i];
        int j;
        for (j = 0; j < l; j++) {
          G[j] += alpha_i * Q_i[j];
        }
        if (is_upper_bound(i)) for (j = 0; j < l; j++) {
          G_bar[j] += get_C(i) * Q_i[j];
        }
      }
  }
  // optimization step
  int iter = 0;
  int counter = min(l, 1000) + 1;
  while (1) {
    // show progress and do shrinking
    if (--counter == 0) {
      counter = min(l, 1000);
      if (shrinking) {
        do_shrinking();
      }
      info(".");
      info_flush();
    }
    int i, j;
    if (select_working_set(i, j) != 0) {
      // reconstruct the whole gradient
      reconstruct_gradient();
      // reset active set size and check
      active_size = l;
      info("*");
      info_flush();
      if (select_working_set(i, j) != 0) {
        break;
      } else {
        counter = 1; // do shrinking next iteration
      }
    }
    ++iter;
    // update alpha[i] and alpha[j], handle bounds carefully
    const Qfloat* Q_i = Q.get_Q(i, active_size);
    const Qfloat* Q_j = Q.get_Q(j, active_size);
    double C_i = get_C(i);
    double C_j = get_C(j);
    double old_alpha_i = alpha[i];
    double old_alpha_j = alpha[j];
    if (y[i] != y[j]) {
      double quad_coef = Q_i[i] + Q_j[j] + 2 * Q_i[j];
      if (quad_coef <= 0) {
        quad_coef = TAU;
      }
      double delta = (-G[i] - G[j]) / quad_coef;
      double diff = alpha[i] - alpha[j];
      alpha[i] += delta;
      alpha[j] += delta;
      if (diff > 0) {
        if (alpha[j] < 0) {
          alpha[j] = 0;
          alpha[i] = diff;
        }
      } else {
        if (alpha[i] < 0) {
          alpha[i] = 0;
          alpha[j] = -diff;
        }
      }
      if (diff > C_i - C_j) {
        if (alpha[i] > C_i) {
          alpha[i] = C_i;
          alpha[j] = C_i - diff;
        }
      } else {
        if (alpha[j] > C_j) {
          alpha[j] = C_j;
          alpha[i] = C_j + diff;
        }
      }
    } else {
      double quad_coef = Q_i[i] + Q_j[j] - 2 * Q_i[j];
      if (quad_coef <= 0) {
        quad_coef = TAU;
      }
      double delta = (G[i] - G[j]) / quad_coef;
      double sum = alpha[i] + alpha[j];
      alpha[i] -= delta;
      alpha[j] += delta;
      if (sum > C_i) {
        if (alpha[i] > C_i) {
          alpha[i] = C_i;
          alpha[j] = sum - C_i;
        }
      } else {
        if (alpha[j] < 0) {
          alpha[j] = 0;
          alpha[i] = sum;
        }
      }
      if (sum > C_j) {
        if (alpha[j] > C_j) {
          alpha[j] = C_j;
          alpha[i] = sum - C_j;
        }
      } else {
        if (alpha[i] < 0) {
          alpha[i] = 0;
          alpha[j] = sum;
        }
      }
    }
    // update G
    double delta_alpha_i = alpha[i] - old_alpha_i;
    double delta_alpha_j = alpha[j] - old_alpha_j;
    for (std::size_t k = 0; static_cast<int>(k) < active_size; k++) {
      G[k] += Q_i[k] * delta_alpha_i + Q_j[k] * delta_alpha_j;
    }
    // update alpha_status and G_bar
    {
      bool ui = is_upper_bound(i);
      bool uj = is_upper_bound(j);
      update_alpha_status(i);
      update_alpha_status(j);
      int k;
      if (ui != is_upper_bound(i)) {
        Q_i = Q.get_Q(i, l);
        if (ui)
          for (k = 0; k < l; k++) {
            G_bar[k] -= C_i * Q_i[k];
          }
        else
          for (k = 0; k < l; k++) {
            G_bar[k] += C_i * Q_i[k];
          }
      }
      if (uj != is_upper_bound(j)) {
        Q_j = Q.get_Q(j, l);
        if (uj)
          for (k = 0; k < l; k++) {
            G_bar[k] -= C_j * Q_j[k];
          }
        else
          for (k = 0; k < l; k++) {
            G_bar[k] += C_j * Q_j[k];
          }
      }
    }
  }
  // calculate rho
  si->rho = calculate_rho();
  // calculate objective value
  {
    double v = 0;
    int i;
    for (i = 0; i < l; i++) {
      v += alpha[i] * (G[i] + p[i]);
    }
    si->obj = v / 2;
  }
  // put back the solution
  {
    for (std::size_t i = 0; static_cast<int>(i) < l; i++) {
      alpha_[active_set[i]] = alpha[i];
    }
  }
  si->upper_bound_p = Cp;
  si->upper_bound_n = Cn;
  info("\noptimization finished, #iter = %d\n", iter);
  delete[] p;
  delete[] y;
  delete[] alpha;
  delete[] alpha_status;
  delete[] active_set;
  delete[] G;
  delete[] G_bar;
}

// return 1 if already optimal, return 0 otherwise
int Solver::select_working_set(int& out_i, int& out_j) {
  // return i,j such that
  // i: maximizes -y_i * grad(f)_i, i in I_up(\alpha)
  // j: minimizes the decrease of obj value
  //    (if quadratic coefficeint <= 0, replace it with tau)
  //    -y_j*grad(f)_j < -y_i*grad(f)_i, j in I_low(\alpha)
  double Gmax = -INF;
  double Gmax2 = -INF;
  int Gmax_idx = -1;
  int Gmin_idx = -1;
  double obj_diff_min = INF;
  for (int t = 0; t < active_size; t++)
    if (y[t] == +1) {
      if (!is_upper_bound(t)) if (-G[t] >= Gmax) {
        Gmax = -G[t];
        Gmax_idx = t;
      }
    } else {
      if (!is_lower_bound(t)) if (G[t] >= Gmax) {
        Gmax = G[t];
        Gmax_idx = t;
      }
    }
  int i = Gmax_idx;
  const Qfloat* Q_i = NULL;
  if (i != -1) { // NULL Q_i not accessed: Gmax=-INF if i=-1
    Q_i = Q->get_Q(i, active_size);
  }
  for (int j = 0; j < active_size; j++) {
    if (y[j] == +1) {
      if (!is_lower_bound(j)) {
        double grad_diff = Gmax + G[j];
        if (G[j] >= Gmax2) {
          Gmax2 = G[j];
        }
        if (grad_diff > 0) {
          double obj_diff;
          double quad_coef = Q_i[i] + QD[j] - 2.f * y[i] * Q_i[j];
          if (quad_coef > 0) {
            obj_diff = -(grad_diff * grad_diff) / quad_coef;
          } else {
            obj_diff = -(grad_diff * grad_diff) / TAU;
          }
          if (obj_diff <= obj_diff_min) {
            Gmin_idx = j;
            obj_diff_min = obj_diff;
          }
        }
      }
    } else {
      if (!is_upper_bound(j)) {
        double grad_diff = Gmax - G[j];
        if (-G[j] >= Gmax2) {
          Gmax2 = -G[j];
        }
        if (grad_diff > 0) {
          double obj_diff;
          double quad_coef = Q_i[i] + QD[j] + 2.f * y[i] * Q_i[j];
          if (quad_coef > 0) {
            obj_diff = -(grad_diff * grad_diff) / quad_coef;
          } else {
            obj_diff = -(grad_diff * grad_diff) / TAU;
          }
          if (obj_diff <= obj_diff_min) {
            Gmin_idx = j;
            obj_diff_min = obj_diff;
          }
        }
      }
    }
  }
  if (Gmax + Gmax2 < eps) {
    return 1;
  }
  out_i = Gmax_idx;
  out_j = Gmin_idx;
  return 0;
}

bool Solver::be_shrunken(int i, double Gmax1, double Gmax2) {
  if (is_upper_bound(i)) {
    if (y[i] == +1) {
      return (-G[i] > Gmax1);
    } else {
      return (-G[i] > Gmax2);
    }
  } else if (is_lower_bound(i)) {
    if (y[i] == +1) {
      return (G[i] > Gmax2);
    } else {
      return (G[i] > Gmax1);
    }
  } else {
    return (false);
  }
}

void Solver::do_shrinking() {
  int i;
  double Gmax1 = -INF; // max { -y_i * grad(f)_i | i in I_up(\alpha) }
  double Gmax2 = -INF; // max { y_i * grad(f)_i | i in I_low(\alpha) }
  // find maximal violating pair first
  for (i = 0; i < active_size; i++) {
    if (y[i] == +1) {
      if (!is_upper_bound(i)) {
        if (-G[i] >= Gmax1) {
          Gmax1 = -G[i];
        }
      }
      if (!is_lower_bound(i)) {
        if (G[i] >= Gmax2) {
          Gmax2 = G[i];
        }
      }
    } else {
      if (!is_upper_bound(i)) {
        if (-G[i] >= Gmax2) {
          Gmax2 = -G[i];
        }
      }
      if (!is_lower_bound(i)) {
        if (G[i] >= Gmax1) {
          Gmax1 = G[i];
        }
      }
    }
  }
  // shrink
  for (i = 0; i < active_size; i++)
    if (be_shrunken(i, Gmax1, Gmax2)) {
      active_size--;
      while (active_size > i) {
        if (!be_shrunken(active_size, Gmax1, Gmax2)) {
          swap_index(i, active_size);
          break;
        }
        active_size--;
      }
    }
  // unshrink, check all variables again before final iterations
  if (unshrinked || Gmax1 + Gmax2 > eps * 10) {
    return;
  }
  unshrinked = true;
  reconstruct_gradient();
  for (i = l - 1; i >= active_size; i--)
    if (!be_shrunken(i, Gmax1, Gmax2)) {
      while (active_size < i) {
        if (be_shrunken(active_size, Gmax1, Gmax2)) {
          swap_index(i, active_size);
          break;
        }
        active_size++;
      }
      active_size++;
    }
}

double Solver::calculate_rho() {
  double r;
  int nr_free = 0;
  double ub = INF, lb = -INF, sum_free = 0;
  for (int i = 0; i < active_size; i++) {
    double yG = y[i] * G[i];
    if (is_upper_bound(i)) {
      if (y[i] == -1) {
        ub = min(ub, yG);
      } else {
        lb = max(lb, yG);
      }
    } else if (is_lower_bound(i)) {
      if (y[i] == +1) {
        ub = min(ub, yG);
      } else {
        lb = max(lb, yG);
      }
    } else {
      ++nr_free;
      sum_free += yG;
    }
  }
  if (nr_free > 0) {
    r = sum_free / nr_free;
  } else {
    r = (ub + lb) / 2;
  }
  return r;
}

//
// Solver for nu-svm classification and regression
//
// additional constraint: e^T \alpha = constant
//
class Solver_NU : public Solver {
  public:
    Solver_NU() {
    }
    void Solve(int l, const QMatrix& Q, const double* p, const schar* y,
               double* alpha, double Cp, double Cn, double eps,
               SolutionInfo* si, int shrinking) {
      this->si = si;
      Solver::Solve(l, Q, p, y, alpha, Cp, Cn, eps, si, shrinking);
    }
  private:
    SolutionInfo* si;
    int select_working_set(int& i, int& j);
    double calculate_rho();
    bool be_shrunken(int i, double Gmax1, double Gmax2, double Gmax3,
                     double Gmax4);
    void do_shrinking();
};

// return 1 if already optimal, return 0 otherwise
int Solver_NU::select_working_set(int& out_i, int& out_j) {
  // return i,j such that y_i = y_j and
  // i: maximizes -y_i * grad(f)_i, i in I_up(\alpha)
  // j: minimizes the decrease of obj value
  //    (if quadratic coefficeint <= 0, replace it with tau)
  //    -y_j*grad(f)_j < -y_i*grad(f)_i, j in I_low(\alpha)
  double Gmaxp = -INF;
  double Gmaxp2 = -INF;
  int Gmaxp_idx = -1;
  double Gmaxn = -INF;
  double Gmaxn2 = -INF;
  int Gmaxn_idx = -1;
  int Gmin_idx = -1;
  double obj_diff_min = INF;
  for (int t = 0; t < active_size; t++)
    if (y[t] == +1) {
      if (!is_upper_bound(t)) if (-G[t] >= Gmaxp) {
        Gmaxp = -G[t];
        Gmaxp_idx = t;
      }
    } else {
      if (!is_lower_bound(t)) if (G[t] >= Gmaxn) {
        Gmaxn = G[t];
        Gmaxn_idx = t;
      }
    }
  int ip = Gmaxp_idx;
  int in = Gmaxn_idx;
  const Qfloat* Q_ip = NULL;
  const Qfloat* Q_in = NULL;
  if (ip != -1) { // NULL Q_ip not accessed: Gmaxp=-INF if ip=-1
    Q_ip = Q->get_Q(ip, active_size);
  }
  if (in != -1) {
    Q_in = Q->get_Q(in, active_size);
  }
  for (int j = 0; j < active_size; j++) {
    if (y[j] == +1) {
      if (!is_lower_bound(j)) {
        double grad_diff = Gmaxp + G[j];
        if (G[j] >= Gmaxp2) {
          Gmaxp2 = G[j];
        }
        if (grad_diff > 0) {
          double obj_diff;
          double quad_coef = Q_ip[ip] + QD[j] - 2 * Q_ip[j];
          if (quad_coef > 0) {
            obj_diff = -(grad_diff * grad_diff) / quad_coef;
          } else {
            obj_diff = -(grad_diff * grad_diff) / TAU;
          }
          if (obj_diff <= obj_diff_min) {
            Gmin_idx = j;
            obj_diff_min = obj_diff;
          }
        }
      }
    } else {
      if (!is_upper_bound(j)) {
        double grad_diff = Gmaxn - G[j];
        if (-G[j] >= Gmaxn2) {
          Gmaxn2 = -G[j];
        }
        if (grad_diff > 0) {
          double obj_diff;
          double quad_coef = Q_in[in] + QD[j] - 2 * Q_in[j];
          if (quad_coef > 0) {
            obj_diff = -(grad_diff * grad_diff) / quad_coef;
          } else {
            obj_diff = -(grad_diff * grad_diff) / TAU;
          }
          if (obj_diff <= obj_diff_min) {
            Gmin_idx = j;
            obj_diff_min = obj_diff;
          }
        }
      }
    }
  }
  if (max(Gmaxp + Gmaxp2, Gmaxn + Gmaxn2) < eps) {
    return 1;
  }
  if (y[Gmin_idx] == +1) {
    out_i = Gmaxp_idx;
  } else {
    out_i = Gmaxn_idx;
  }
  out_j = Gmin_idx;
  return 0;
}

bool Solver_NU::be_shrunken(int i, double Gmax1, double Gmax2,
                            double Gmax3, double Gmax4) {
  if (is_upper_bound(i)) {
    if (y[i] == +1) {
      return (-G[i] > Gmax1);
    } else {
      return (-G[i] > Gmax4);
    }
  } else if (is_lower_bound(i)) {
    if (y[i] == +1) {
      return (G[i] > Gmax2);
    } else {
      return (G[i] > Gmax3);
    }
  } else {
    return (false);
  }
}

void Solver_NU::do_shrinking() {
  double Gmax1 = -INF; // max { -y_i * grad(f)_i | y_i = +1, i in I_up(\alpha) }
  double Gmax2 = -INF; // max { y_i * grad(f)_i | y_i = +1, i in I_low(\alpha) }
  double Gmax3 = -INF; // max { -y_i * grad(f)_i | y_i = -1, i in I_up(\alpha) }
  double Gmax4 = -INF; // max { y_i * grad(f)_i | y_i = -1, i in I_low(\alpha) }
  // find maximal violating pair first
  int i;
  for (i = 0; i < active_size; i++) {
    if (!is_upper_bound(i)) {
      if (y[i] == +1) {
        if (-G[i] > Gmax1) {
          Gmax1 = -G[i];
        }
      } else if (-G[i] > Gmax4) {
        Gmax4 = -G[i];
      }
    }
    if (!is_lower_bound(i)) {
      if (y[i] == +1) {
        if (G[i] > Gmax2) {
          Gmax2 = G[i];
        }
      } else if (G[i] > Gmax3) {
        Gmax3 = G[i];
      }
    }
  }
  // shrinking
  for (i = 0; i < active_size; i++)
    if (be_shrunken(i, Gmax1, Gmax2, Gmax3, Gmax4)) {
      active_size--;
      while (active_size > i) {
        if (!be_shrunken(active_size, Gmax1, Gmax2, Gmax3, Gmax4)) {
          swap_index(i, active_size);
          break;
        }
        active_size--;
      }
    }
  // unshrink, check all variables again before final iterations
  if (unshrinked || max(Gmax1 + Gmax2, Gmax3 + Gmax4) > eps * 10) {
    return;
  }
  unshrinked = true;
  reconstruct_gradient();
  for (i = l - 1; i >= active_size; i--)
    if (!be_shrunken(i, Gmax1, Gmax2, Gmax3, Gmax4)) {
      while (active_size < i) {
        if (be_shrunken(active_size, Gmax1, Gmax2, Gmax3, Gmax4)) {
          swap_index(i, active_size);
          break;
        }
        active_size++;
      }
      active_size++;
    }
}

double Solver_NU::calculate_rho() {
  int nr_free1 = 0, nr_free2 = 0;
  double ub1 = INF, ub2 = INF;
  double lb1 = -INF, lb2 = -INF;
  double sum_free1 = 0, sum_free2 = 0;
  for (int i = 0; i < active_size; i++) {
    if (y[i] == +1) {
      if (is_upper_bound(i)) {
        lb1 = max(lb1, G[i]);
      } else if (is_lower_bound(i)) {
        ub1 = min(ub1, G[i]);
      } else {
        ++nr_free1;
        sum_free1 += G[i];
      }
    } else {
      if (is_upper_bound(i)) {
        lb2 = max(lb2, G[i]);
      } else if (is_lower_bound(i)) {
        ub2 = min(ub2, G[i]);
      } else {
        ++nr_free2;
        sum_free2 += G[i];
      }
    }
  }
  double r1, r2;
  if (nr_free1 > 0) {
    r1 = sum_free1 / nr_free1;
  } else {
    r1 = (ub1 + lb1) / 2;
  }
  if (nr_free2 > 0) {
    r2 = sum_free2 / nr_free2;
  } else {
    r2 = (ub2 + lb2) / 2;
  }
  si->r = (r1 + r2) / 2;
  return (r1 - r2) / 2;
}

//
// Q matrices for various formulations
//
class SVC_Q : public Kernel {
  public:
    SVC_Q(const svm_problem& prob, const svm_parameter& param,
          const schar* y_) :
      Kernel(static_cast<int>(prob.l), prob.x, param) {
      clone(y, y_, static_cast<int>(prob.l));
      cache = new Cache(prob.l, (long int)(param.cache_size * (1 << 20)));
      QD = new Qfloat[prob.l];
      for (int i = 0; static_cast<std::size_t>(i) < prob.l; i++) {
        QD[i] = (Qfloat)(this->*kernel_function)(i, i);
      }
    }

    Qfloat* get_Q(int i, int len) const {
      Qfloat* data;
      int start;
      if ((start = cache->get_data(i, &data, len)) < len) {
        for (int j = start; j < len; j++) {
          data[j] = (Qfloat)(y[i] * y[j] * (this->*kernel_function)(i, j));
        }
      }
      return data;
    }

    Qfloat* get_QD() const {
      return QD;
    }

    void swap_index(int i, int j) const {
      cache->swap_index(i, j);
      Kernel::swap_index(i, j);
      swap(y[i], y[j]);
      swap(QD[i], QD[j]);
    }

    ~SVC_Q() {
      delete[] y;
      delete cache;
      delete[] QD;
    }
  private:
    schar* y;
    Cache* cache;
    Qfloat* QD;
};

class ONE_CLASS_Q : public Kernel {
  public:
    ONE_CLASS_Q(const svm_problem& prob, const svm_parameter& param) :
      Kernel(static_cast<int>(prob.l), prob.x, param) {
      cache = new Cache(prob.l, (long int)(param.cache_size * (1 << 20)));
      QD = new Qfloat[prob.l];
      for (int i = 0; static_cast<std::size_t>(i) < prob.l; i++) {
        QD[i] = (Qfloat)(this->*kernel_function)(i, i);
      }
    }

    Qfloat* get_Q(int i, int len) const {
      Qfloat* data;
      int start;
      if ((start = cache->get_data(i, &data, len)) < len) {
        for (int j = start; j < len; j++) {
          data[j] = (Qfloat)(this->*kernel_function)(i, j);
        }
      }
      return data;
    }

    Qfloat* get_QD() const {
      return QD;
    }

    void swap_index(int i, int j) const {
      cache->swap_index(i, j);
      Kernel::swap_index(i, j);
      swap(QD[i], QD[j]);
    }

    ~ONE_CLASS_Q() {
      delete cache;
      delete[] QD;
    }
  private:
    Cache* cache;
    Qfloat* QD;
};

class SVR_Q : public Kernel {
  public:
    SVR_Q(const svm_problem& prob, const svm_parameter& param) :
      Kernel(static_cast<int>(prob.l), prob.x, param) {
      l = prob.l;
      cache = new Cache(l, (long int)(param.cache_size * (1 << 20)));
      QD = new Qfloat[2 * l];
      sign = new schar[2 * l];
      index = new int[2 * l];
      for (std::size_t k = 0; k < l; k++) {
        sign[k] = 1;
        sign[k + l] = -1;
        index[k] = static_cast<int>(k);
        index[k + l] = static_cast<int>(k);
        QD[k] = (Qfloat)(this->*kernel_function)(static_cast<int>(k), static_cast<int>(k));
        QD[k + l] = QD[k];
      }
      buffer[0] = new Qfloat[2 * l];
      buffer[1] = new Qfloat[2 * l];
      next_buffer = 0;
    }

    void swap_index(int i, int j) const {
      swap(sign[i], sign[j]);
      swap(index[i], index[j]);
      swap(QD[i], QD[j]);
    }

    Qfloat* get_Q(int i, int len) const {
      Qfloat* data;
      int real_i = index[i];
      if (cache->get_data(real_i, &data, static_cast<int>(l)) < static_cast<int>(l)) {
        for (int j = 0; static_cast<std::size_t>(j) < l; j++) {
          data[j] = (Qfloat)(this->*kernel_function)(real_i, j);
        }
      }
      // reorder and copy
      Qfloat* buf = buffer[next_buffer];
      next_buffer = 1 - next_buffer;
      schar si = sign[i];
      for (int j = 0; j < len; j++) {
        buf[j] = static_cast<Qfloat>(si * sign[j]) * data[index[j]];
      }
      return buf;
    }

    Qfloat* get_QD() const {
      return QD;
    }

    ~SVR_Q() {
      delete cache;
      delete[] sign;
      delete[] index;
      delete[] buffer[0];
      delete[] buffer[1];
      delete[] QD;
    }
  private:
    std::size_t l;
    Cache* cache;
    schar* sign;
    int* index;
    mutable int next_buffer;
    Qfloat* buffer[2];
    Qfloat* QD;
};

//
// construct and solve various formulations
//
static void solve_c_svc(const svm_problem* prob,
                        const svm_parameter* param, double* alpha,
                        Solver::SolutionInfo* si, double Cp, double Cn) {
  std::size_t l = prob->l;
  double* minus_ones = new double[l];
  schar* y = new schar[l];
  int i;
  for (i = 0; i < static_cast<int>(l); i++) {
    alpha[i] = 0;
    minus_ones[i] = -1;
    if (prob->y[i] > 0) {
      y[i] = +1;
    } else {
      y[i] = -1;
    }
  }
  Solver s;
  s.Solve(static_cast<int>(l),
          SVC_Q(*prob, *param, y),
          minus_ones,
          y,
          alpha,
          Cp,
          Cn,
          param->eps,
          si,
          param->shrinking);
  double sum_alpha = 0;
  for (i = 0; i < static_cast<int>(l); i++) {
    sum_alpha += alpha[i];
  }
  if (Cp == Cn) {
    info("nu = %f\n", sum_alpha / (Cp * static_cast<double>(prob->l)));
  }
  for (i = 0; i < static_cast<int>(l); i++) {
    alpha[i] *= y[i];
  }
  delete[] minus_ones;
  delete[] y;
}

static void solve_nu_svc(const svm_problem* prob,
                         const svm_parameter* param, double* alpha,
                         Solver::SolutionInfo* si) {
  int i;
  std::size_t l = prob->l;
  double nu = param->nu;
  schar* y = new schar[l];
  for (i = 0; i < static_cast<int>(l); i++)
    if (prob->y[i] > 0) {
      y[i] = +1;
    } else {
      y[i] = -1;
    }
  double sum_pos = nu * static_cast<double>(l) / 2.;
  double sum_neg = nu * static_cast<double>(l) / 2.;
  for (i = 0; i < static_cast<int>(l); i++)
    if (y[i] == +1) {
      alpha[i] = min(1.0, sum_pos);
      sum_pos -= alpha[i];
    } else {
      alpha[i] = min(1.0, sum_neg);
      sum_neg -= alpha[i];
    }
  double* zeros = new double[l];
  for (i = 0; i < static_cast<int>(l); i++) {
    zeros[i] = 0;
  }
  Solver_NU s;
  s.Solve(static_cast<int>(l),
          SVC_Q(*prob, *param, y),
          zeros,
          y,
          alpha,
          1.0,
          1.0,
          param->eps,
          si,
          param->shrinking);
  double r = si->r;
  info("C = %f\n", 1 / r);
  for (i = 0; i < static_cast<int>(l); i++) {
    alpha[i] *= y[i] / r;
  }
  si->rho /= r;
  si->obj /= (r * r);
  si->upper_bound_p = 1 / r;
  si->upper_bound_n = 1 / r;
  delete[] y;
  delete[] zeros;
}

static void solve_one_class(const svm_problem* prob,
                            const svm_parameter* param, double* alpha,
                            Solver::SolutionInfo* si) {
  std::size_t l = prob->l;
  double* zeros = new double[l];
  schar* ones = new schar[l];
  int i;
  int n = (int)(param->nu * static_cast<double>(prob->l)); // # of alpha's at upper bound
  for (i = 0; i < n; i++) {
    alpha[i] = 1;
  }
  if (n < static_cast<int>(prob->l)) {
    alpha[n] = param->nu * static_cast<double>(prob->l) - n;
  }
  for (i = n + 1; i < static_cast<int>(l); i++) {
    alpha[i] = 0;
  }
  for (i = 0; i < static_cast<int>(l); i++) {
    zeros[i] = 0;
    ones[i] = 1;
  }
  Solver s;
  s.Solve(static_cast<int>(l),
          ONE_CLASS_Q(*prob, *param),
          zeros,
          ones,
          alpha,
          1.0,
          1.0,
          param->eps,
          si,
          param->shrinking);
  delete[] zeros;
  delete[] ones;
}

static void solve_epsilon_svr(const svm_problem* prob,
                              const svm_parameter* param, double* alpha,
                              Solver::SolutionInfo* si) {
  std::size_t l = prob->l;
  double* alpha2 = new double[2 * l];
  double* linear_term = new double[2 * l];
  schar* y = new schar[2 * l];
  std::size_t i;
  for (i = 0; i < l; i++) {
    alpha2[i] = 0;
    linear_term[i] = param->p - prob->y[i];
    y[i] = 1;
    alpha2[i + l] = 0;
    linear_term[i + l] = param->p + prob->y[i];
    y[i + l] = -1;
  }
  Solver s;
  s.Solve(2 * static_cast<int>(l),
          SVR_Q(*prob, *param),
          linear_term,
          y,
          alpha2,
          param->C,
          param->C,
          param->eps,
          si,
          param->shrinking);
  double sum_alpha = 0;
  for (i = 0; i < l; i++) {
    alpha[i] = alpha2[i] - alpha2[i + l];
    sum_alpha += fabs(alpha[i]);
  }
  info("nu = %f\n", sum_alpha / (param->C * static_cast<double>(l)));
  delete[] alpha2;
  delete[] linear_term;
  delete[] y;
}

static void solve_nu_svr(const svm_problem* prob,
                         const svm_parameter* param, double* alpha,
                         Solver::SolutionInfo* si) {
  std::size_t l = prob->l;
  double C = param->C;
  double* alpha2 = new double[2 * l];
  double* linear_term = new double[2 * l];
  schar* y = new schar[2 * l];
  std::size_t i;
  double sum = C * param->nu * static_cast<double>(l) / 2.;
  for (i = 0; i < l; i++) {
    alpha2[i] = alpha2[i + l] = min(sum, C);
    sum -= alpha2[i];
    linear_term[i] = -prob->y[i];
    y[i] = 1;
    linear_term[i + l] = prob->y[i];
    y[i + l] = -1;
  }
  Solver_NU s;
  s.Solve(2 * static_cast<int>(l),
          SVR_Q(*prob, *param),
          linear_term,
          y,
          alpha2,
          C,
          C,
          param->eps,
          si,
          param->shrinking);
  info("epsilon = %f\n", -si->r);
  for (i = 0; i < l; i++) {
    alpha[i] = alpha2[i] - alpha2[i + l];
  }
  delete[] alpha2;
  delete[] linear_term;
  delete[] y;
}

//
// decision_function
//
struct decision_function {
    double* alpha;
    double rho;
};

decision_function svm_train_one(const svm_problem* prob,
                                const svm_parameter* param, double Cp,
                                double Cn) {
  double* alpha = Malloc(double, prob->l);
  Solver::SolutionInfo si;
  switch (param->svm_type) {
    case C_SVC:
      solve_c_svc(prob, param, alpha, &si, Cp, Cn);
      break;
    case NU_SVC:
      solve_nu_svc(prob, param, alpha, &si);
      break;
    case ONE_CLASS:
      solve_one_class(prob, param, alpha, &si);
      break;
    case EPSILON_SVR:
      solve_epsilon_svr(prob, param, alpha, &si);
      break;
    case NU_SVR:
      solve_nu_svr(prob, param, alpha, &si);
      break;
  }
  info("obj = %f, rho = %f\n", si.obj, si.rho);
  // output SVs
  int nSV = 0;
  int nBSV = 0;
  for (std::size_t i = 0; i < prob->l; i++) {
    if (fabs(alpha[i]) > 0) {
      ++nSV;
      if (prob->y[i] > 0) {
        if (fabs(alpha[i]) >= si.upper_bound_p) {
          ++nBSV;
        }
      } else {
        if (fabs(alpha[i]) >= si.upper_bound_n) {
          ++nBSV;
        }
      }
    }
  }
  info("nSV = %d, nBSV = %d\n", nSV, nBSV);
  decision_function f;
  f.alpha = alpha;
  f.rho = si.rho;
  return f;
}

// Platt's binary SVM Probablistic Output: an improvement from Lin et al.
void sigmoid_train(std::size_t l, const double* dec_values, const double* labels,
                   double& A, double& B) {
  double prior1 = 0, prior0 = 0;
  std::size_t i;
  for (i = 0; i < l; i++)
    if (labels[i] > 0) {
      prior1 += 1;
    } else {
      prior0 += 1;
    }
  int max_iter = 100; // Maximal number of iterations
  double min_step = 1e-10; // Minimal step taken in line search
  double sigma = 1e-12; // For numerically strict PD of Hessian
  double eps = 1e-5;
  double hiTarget = (prior1 + 1.0) / (prior1 + 2.0);
  double loTarget = 1 / (prior0 + 2.0);
  double* t = Malloc(double, l);
  double fApB, p, q, h11, h22, h21, g1, g2, det, dA, dB, gd, stepsize;
  double newA, newB, newf, d1, d2;
  int iter;
  // Initial Point and Initial Fun Value
  A = 0.0;
  B = log((prior0 + 1.0) / (prior1 + 1.0));
  double fval = 0.0;
  for (i = 0; i < l; i++) {
    if (labels[i] > 0) {
      t[i] = hiTarget;
    } else {
      t[i] = loTarget;
    }
    fApB = dec_values[i] * A + B;
    if (fApB >= 0) {
      fval += t[i] * fApB + log(1 + exp(-fApB));
    } else {
      fval += (t[i] - 1) * fApB + log(1 + exp(fApB));
    }
  }
  for (iter = 0; iter < max_iter; iter++) {
    // Update Gradient and Hessian (use H' = H + sigma I)
    h11 = sigma; // numerically ensures strict PD
    h22 = sigma;
    h21 = 0.0;
    g1 = 0.0;
    g2 = 0.0;
    for (i = 0; i < l; i++) {
      fApB = dec_values[i] * A + B;
      if (fApB >= 0) {
        p = exp(-fApB) / (1.0 + exp(-fApB));
        q = 1.0 / (1.0 + exp(-fApB));
      } else {
        p = 1.0 / (1.0 + exp(fApB));
        q = exp(fApB) / (1.0 + exp(fApB));
      }
      d2 = p * q;
      h11 += dec_values[i] * dec_values[i] * d2;
      h22 += d2;
      h21 += dec_values[i] * d2;
      d1 = t[i] - p;
      g1 += dec_values[i] * d1;
      g2 += d1;
    }
    // Stopping Criteria
    if (fabs(g1) < eps && fabs(g2) < eps) {
      break;
    }
    // Finding Newton direction: -inv(H') * g
    det = h11 * h22 - h21 * h21;
    dA = -(h22 * g1 - h21 * g2) / det;
    dB = -(-h21 * g1 + h11 * g2) / det;
    gd = g1 * dA + g2 * dB;
    stepsize = 1; // Line Search
    while (stepsize >= min_step) {
      newA = A + stepsize * dA;
      newB = B + stepsize * dB;
      // New function value
      newf = 0.0;
      for (i = 0; i < l; i++) {
        fApB = dec_values[i] * newA + newB;
        if (fApB >= 0) {
          newf += t[i] * fApB + log(1 + exp(-fApB));
        } else {
          newf += (t[i] - 1) * fApB + log(1 + exp(fApB));
        }
      }
      // Check sufficient decrease
      if (newf < fval + 0.0001 * stepsize * gd) {
        A = newA;
        B = newB;
        fval = newf;
        break;
      } else {
        stepsize = stepsize / 2.0;
      }
    }
    if (stepsize < min_step) {
      info("Line search fails in two-class probability estimates\n");
      break;
    }
  }
  if (iter >= max_iter) {
    info("Reaching maximal iterations in two-class probability estimates\n");
  }
  free(t);
}

double sigmoid_predict(double decision_value, double A, double B) {
  double fApB = decision_value * A + B;
  if (fApB >= 0) {
    return exp(-fApB) / (1.0 + exp(-fApB));
  } else {
    return 1.0 / (1 + exp(fApB));
  }
}

// Method 2 from the multiclass_prob paper by Wu, Lin, and Weng
void multiclass_probability(std::size_t k, double** r, double* p) {
  std::size_t t, j;
  int iter = 0, max_iter = static_cast<int>(max(static_cast<std::size_t>(100), k));
  double** Q = Malloc(double*, k);
  double* Qp = Malloc(double, k);
  double pQp, eps = 0.005 / static_cast<double>(k);
  for (t = 0; t < k; t++) {
    p[t] = 1.0 / static_cast<int>(k); // Valid if k = 1
    Q[t] = Malloc(double, k);
    Q[t][t] = 0;
    for (j = 0; j < t; j++) {
      Q[t][t] += r[j][t] * r[j][t];
      Q[t][j] = Q[j][t];
    }
    for (j = t + 1; j < k; j++) {
      Q[t][t] += r[j][t] * r[j][t];
      Q[t][j] = -r[j][t] * r[t][j];
    }
  }
  for (iter = 0; iter < max_iter; iter++) {
    // stopping condition, recalculate QP,pQP for numerical accuracy
    pQp = 0;
    for (t = 0; t < k; t++) {
      Qp[t] = 0;
      for (j = 0; j < k; j++) {
        Qp[t] += Q[t][j] * p[j];
      }
      pQp += p[t] * Qp[t];
    }
    double max_error = 0;
    for (t = 0; t < k; t++) {
      double error = fabs(Qp[t] - pQp);
      if (error > max_error) {
        max_error = error;
      }
    }
    if (max_error < eps) {
      break;
    }
    for (t = 0; t < k; t++) {
      double diff = (-Qp[t] + pQp) / Q[t][t];
      p[t] += diff;
      pQp = (pQp + diff * (diff * Q[t][t] + 2 * Qp[t])) / (1 + diff) / (1
          + diff);
      for (j = 0; j < k; j++) {
        Qp[j] = (Qp[j] + diff * Q[t][j]) / (1 + diff);
        p[j] /= (1 + diff);
      }
    }
  }
  if (iter >= max_iter) {
    info("Exceeds max_iter in multiclass_prob\n");
  }
  for (t = 0; t < k; t++) {
    free(Q[t]);
  }
  free(Q);
  free(Qp);
}

// Cross-validation decision values for probability estimates
void svm_binary_svc_probability(const svm_problem* prob,
                                const svm_parameter* param, double Cp,
                                double Cn, double& probA, double& probB) {
  int i;
  std::size_t nr_fold = 5;
  int* perm = Malloc(int, prob->l);
  double* dec_values = Malloc(double, prob->l);
  // random shuffle
  for (i = 0; i < static_cast<int>(prob->l); i++) {
    perm[i] = i;
  }
  for (i = 0; i < static_cast<int>(prob->l); i++) {
    int j = i + static_cast<int>(
      PseudoRandom::lcg_rand() % 
      static_cast<long unsigned int>(static_cast<int>(prob->l) - i));
    swap(perm[i], perm[j]);
  }
  for (i = 0; i < static_cast<int>(nr_fold); i++) {
    int begin = static_cast<int>(static_cast<std::size_t>(i) * prob->l / nr_fold);
    int end = static_cast<int>(static_cast<std::size_t>(i + 1) * prob->l / nr_fold);
    int j, k;
    struct svm_problem subprob;
    subprob.l = prob->l - static_cast<std::size_t>(end - begin);
#ifdef _DENSE_REP
    subprob.x = Malloc(struct svm_node, subprob.l);
#else
    subprob.x = Malloc(struct svm_node*, subprob.l);
#endif
    subprob.y = Malloc(double, subprob.l);
    k = 0;
    for (j = 0; j < begin; j++) {
      subprob.x[k] = prob->x[perm[j]];
      subprob.y[k] = prob->y[perm[j]];
      ++k;
    }
    for (j = end; j < static_cast<int>(prob->l); j++) {
      subprob.x[k] = prob->x[perm[j]];
      subprob.y[k] = prob->y[perm[j]];
      ++k;
    }
    int p_count = 0, n_count = 0;
    for (j = 0; j < k; j++)
      if (subprob.y[j] > 0) {
        p_count++;
      } else {
        n_count++;
      }
    if (p_count == 0 && n_count == 0)
      for (j = begin; j < end; j++) {
        dec_values[perm[j]] = 0;
      }
    else if (p_count > 0 && n_count == 0)
      for (j = begin; j < end; j++) {
        dec_values[perm[j]] = 1;
      }
    else if (p_count == 0 && n_count > 0)
      for (j = begin; j < end; j++) {
        dec_values[perm[j]] = -1;
      }
    else {
      svm_parameter subparam = *param;
      subparam.probability = 0;
      subparam.C = 1.0;
      subparam.nr_weight = 2;
      subparam.weight_label = Malloc(int, 2);
      subparam.weight = Malloc(double, 2);
      subparam.weight_label[0] = +1;
      subparam.weight_label[1] = -1;
      subparam.weight[0] = Cp;
      subparam.weight[1] = Cn;
      struct svm_model* submodel = svm_train(&subprob, &subparam);
      for (j = begin; j < end; j++) {
#ifdef _DENSE_REP
        svm_predict_values(submodel,
                           (prob->x + perm[j]),
                           &(dec_values[perm[j]]));
#else
        svm_predict_values(submodel, prob->x[perm[j]], &(dec_values[perm[j]]));
#endif
        // ensure +1 -1 order; reason not using CV subroutine
        dec_values[perm[j]] *= submodel->label[0];
      }
      svm_destroy_model(submodel);
      svm_destroy_param(&subparam);
    }
    free(subprob.x);
    free(subprob.y);
  }
  sigmoid_train(prob->l, dec_values, prob->y, probA, probB);
  free(dec_values);
  free(perm);
}

// Return parameter of a Laplace distribution
double svm_svr_probability(const svm_problem* prob,
                           const svm_parameter* param) {
  int i;
  int nr_fold = 5;
  double* ymv = Malloc(double, prob->l);
  double mae = 0;
  svm_parameter newparam = *param;
  newparam.probability = 0;
  svm_cross_validation(prob, &newparam, nr_fold, ymv);
  for (i = 0; i < static_cast<int>(prob->l); i++) {
    ymv[i] = prob->y[i] - ymv[i];
    mae += fabs(ymv[i]);
  }
  mae /= static_cast<double>(prob->l);
  double std = sqrt(2 * mae * mae);
  int count = 0;
  mae = 0;
  for (i = 0; i < static_cast<int>(prob->l); i++)
    if (fabs(ymv[i]) > 5 * std) {
      count = count + 1;
    } else {
      mae += fabs(ymv[i]);
    }
  mae /= static_cast<double>(static_cast<int>(prob->l) - count);
  info("Prob. model for test data: target value = predicted value + z,\nz: Laplace distribution e^(-|z|/sigma)/(2sigma),sigma= %g\n",
       mae);
  free(ymv);
  return mae;
}

// label: label name, start: begin of each class, count: #data of classes, perm: indices to the original data
// perm, length l, must be allocated before calling this subroutine
void svm_group_classes(const svm_problem* prob, int* nr_class_ret,
                       int** label_ret, int** start_ret, int** count_ret,
                       int* perm) {
  std::size_t l = prob->l;
  std::size_t max_nr_class = 16;
  std::size_t nr_class = 0;
  int* label = Malloc(int, max_nr_class);
  int* count = Malloc(int, max_nr_class);
  int* data_label = Malloc(int, l);
  int i;
  for (i = 0; i < static_cast<int>(l); i++) {
    int this_label = (int)prob->y[i];
    int j;
    for (j = 0; j < static_cast<int>(nr_class); j++) {
      if (this_label == label[j]) {
        ++count[j];
        break;
      }
    }
    data_label[i] = j;
    if (j == static_cast<int>(nr_class)) {
      if (nr_class == max_nr_class) {
        max_nr_class *= 2;
        label = (int*)realloc(label, max_nr_class * sizeof(int));
        count = (int*)realloc(count, max_nr_class * sizeof(int));
      }
      label[nr_class] = this_label;
      count[nr_class] = 1;
      ++nr_class;
    }
  }
  int* start = Malloc(int, nr_class);
  start[0] = 0;
  for (i = 1; i < static_cast<int>(nr_class); i++) {
    start[i] = start[i - 1] + count[i - 1];
  }
  for (i = 0; i < static_cast<int>(l); i++) {
    perm[start[data_label[i]]] = i;
    ++start[data_label[i]];
  }
  start[0] = 0;
  for (i = 1; i < static_cast<int>(nr_class); i++) {
    start[i] = start[i - 1] + count[i - 1];
  }
  *nr_class_ret = static_cast<int>(nr_class);
  *label_ret = label;
  *start_ret = start;
  *count_ret = count;
  free(data_label);
}

//
// Interface functions
//
svm_model* svm_train(const svm_problem* prob, const svm_parameter* param) {
  svm_model* model = Malloc(svm_model, 1);
  model->param = *param;
  model->free_sv = 0; // XXX
  if (param->svm_type == ONE_CLASS || param->svm_type == EPSILON_SVR
      || param->svm_type == NU_SVR) {
    // regression or one-class-svm
    model->nr_class = 2;
    model->label = NULL;
    model->nSV = NULL;
    model->probA = NULL;
    model->probB = NULL;
    model->sv_coef = Malloc(double*, 1);
    if (param->probability && (param->svm_type == EPSILON_SVR
        || param->svm_type == NU_SVR)) {
      model->probA = Malloc(double, 1);
      model->probA[0] = svm_svr_probability(prob, param);
    }
    decision_function f = svm_train_one(prob, param, 0, 0);
    model->rho = Malloc(double, 1);
    model->rho[0] = f.rho;
    std::size_t nSV = 0;
    std::size_t i;
    for (i = 0; i < prob->l; i++)
      if (fabs(f.alpha[i]) > 0) {
        ++nSV;
      }
    model->l = static_cast<int>(nSV);
#ifdef _DENSE_REP
    model->SV = Malloc(svm_node, nSV);
#else
    model->SV = Malloc(svm_node*, nSV);
#endif
    model->sv_coef[0] = Malloc(double, nSV);
    int j = 0;
    for (i = 0; i < prob->l; i++)
      if (fabs(f.alpha[i]) > 0) {
        model->SV[j] = prob->x[i];
        model->sv_coef[0][j] = f.alpha[i];
        ++j;
      }
    free(f.alpha);
  } else {
    // classification
    std::size_t l = prob->l;
    int nr_class;
    int* label = NULL;
    int* start = NULL;
    int* count = NULL;
    int* perm = Malloc(int, l);
    // group training data of the same class
    svm_group_classes(prob, &nr_class, &label, &start, &count, perm);
#ifdef _DENSE_REP
    svm_node* x = Malloc(svm_node, l);
#else
    svm_node** x = Malloc(svm_node*, l);
#endif
    int i;
    for (i = 0; i < static_cast<int>(l); i++) {
      x[i] = prob->x[perm[i]];
    }
    // calculate weighted C
    double* weighted_C = Malloc(double, static_cast<std::size_t>(nr_class));
    for (i = 0; i < nr_class; i++) {
      weighted_C[i] = param->C;
    }
    for (i = 0; i < param->nr_weight; i++) {
      int j;
      for (j = 0; j < nr_class; j++)
        if (param->weight_label[i] == label[j]) {
          break;
        }
      if (j == nr_class) {
        fprintf(stderr,
                "warning: class label %d specified in weight is not found\n",
                param->weight_label[i]);
      } else {
        weighted_C[j] *= param->weight[i];
      }
    }
    // train k*(k-1)/2 models
    bool* nonzero = Malloc(bool, l);
    for (i = 0; i < static_cast<int>(l); i++) {
      nonzero[i] = false;
    }
    std::size_t allocationSize = static_cast<std::size_t>(nr_class * (nr_class - 1) / 2);
    decision_function* f =
        Malloc(decision_function, allocationSize);
    double* probA = NULL, *probB = NULL;
    if (param->probability) {
      probA = Malloc(double, allocationSize);
      probB = Malloc(double, allocationSize);
    }
    int p = 0;
    for (i = 0; i < nr_class; i++)
      for (int j = i + 1; j < nr_class; j++) {
        svm_problem sub_prob;
        int si = start[i], sj = start[j];
        int ci = count[i], cj = count[j];
        sub_prob.l = static_cast<std::size_t>(ci + cj);
#ifdef _DENSE_REP
        sub_prob.x = Malloc(svm_node, sub_prob.l);
#else
        sub_prob.x = Malloc(svm_node*, sub_prob.l);
#endif
        sub_prob.y = Malloc(double, sub_prob.l);
        int k;
        for (k = 0; k < ci; k++) {
          sub_prob.x[k] = x[si + k];
          sub_prob.y[k] = +1;
        }
        for (k = 0; k < cj; k++) {
          sub_prob.x[ci + k] = x[sj + k];
          sub_prob.y[ci + k] = -1;
        }
        if (param->probability) {
          svm_binary_svc_probability(&sub_prob,
                                                           param,
                                                           weighted_C[i],
                                                           weighted_C[j],
                                                           probA[p],
                                                           probB[p]);
        }
        f[p] = svm_train_one(&sub_prob,
                             param,
                             weighted_C[i],
                             weighted_C[j]);
        for (k = 0; k < ci; k++)
          if (!nonzero[si + k] && fabs(f[p].alpha[k]) > 0) {
            nonzero[si + k] = true;
          }
        for (k = 0; k < cj; k++)
          if (!nonzero[sj + k] && fabs(f[p].alpha[ci + k]) > 0) {
            nonzero[sj + k] = true;
          }
        free(sub_prob.x);
        free(sub_prob.y);
        ++p;
      }
    // build output
    model->nr_class = nr_class;
    model->label = Malloc(int, static_cast<std::size_t>(nr_class));
    for (i = 0; i < nr_class; i++) {
      model->label[i] = label[i];
    }
    model->rho = Malloc(double, allocationSize);
    for (i = 0; i < static_cast<int>(allocationSize); i++) {
      model->rho[i] = f[i].rho;
    }
    if (param->probability) {
      model->probA = Malloc(double, allocationSize);
      model->probB = Malloc(double, allocationSize);
      for (i = 0; i < static_cast<int>(allocationSize); i++) {
        model->probA[i] = probA[i];
        model->probB[i] = probB[i];
      }
    } else {
      model->probA = NULL;
      model->probB = NULL;
    }
    int total_sv = 0;
    int* nz_count = Malloc(int, static_cast<std::size_t>(nr_class));
    model->nSV = Malloc(int, static_cast<std::size_t>(nr_class));
    for (i = 0; i < nr_class; i++) {
      int nSV = 0;
      for (int j = 0; j < count[i]; j++)
        if (nonzero[start[i] + j]) {
          ++nSV;
          ++total_sv;
        }
      model->nSV[i] = nSV;
      nz_count[i] = nSV;
    }
    info("Total nSV = %d\n", total_sv);
    model->l = total_sv;
#ifdef _DENSE_REP
    model->SV = Malloc(svm_node, static_cast<std::size_t>(total_sv));
#else
    model->SV = Malloc(svm_node*, static_cast<std::size_t>(total_sv));
#endif
    p = 0;
    for (i = 0; i < static_cast<int>(l); i++)
      if (nonzero[i]) {
        model->SV[p++] = x[i];
      }
    int* nz_start = Malloc(int, static_cast<std::size_t>(nr_class));
    nz_start[0] = 0;
    for (i = 1; i < nr_class; i++) {
      nz_start[i] = nz_start[i - 1] + nz_count[i - 1];
    }
    if (nr_class > 1) {
      model->sv_coef = Malloc(double*, static_cast<std::size_t>(nr_class - 1));
      for (i = 0; i < nr_class - 1; i++) {
        model->sv_coef[i] = Malloc(double, static_cast<std::size_t>(total_sv));
      }
    } else {
      fprintf(stderr, "Error: nr_class (%d) must be at least 2\n", nr_class);
      model->sv_coef =  nullptr;   
    }
    p = 0;
    for (i = 0; i < nr_class; i++)
      for (int j = i + 1; j < nr_class; j++) {
        // classifier (i,j): coefficients with
        // i are in sv_coef[j-1][nz_start[i]...],
        // j are in sv_coef[i][nz_start[j]...]
        int si = start[i];
        int sj = start[j];
        int ci = count[i];
        int cj = count[j];
        int q = nz_start[i];
        int k;
        for (k = 0; k < ci; k++)
          if (nonzero[si + k]) {
            model->sv_coef[j - 1][q++] = f[p].alpha[k];
          }
        q = nz_start[j];
        for (k = 0; k < cj; k++)
          if (nonzero[sj + k]) {
            model->sv_coef[i][q++] = f[p].alpha[ci + k];
          }
        ++p;
      }
    free(label);
    free(probA);
    free(probB);
    free(count);
    free(perm);
    free(start);
    free(x);
    free(weighted_C);
    free(nonzero);
    for (i = 0; i < static_cast<int>(allocationSize); i++) {
      free(f[i].alpha);
    }
    free(f);
    free(nz_count);
    free(nz_start);
  }
  return model;
}

// Stratified cross validation
void svm_cross_validation(const svm_problem* prob,
                          const svm_parameter* param, int nr_fold,
                          double* target) {
  int i;
  int* fold_start = Malloc(int, static_cast<std::size_t>(nr_fold + 1));
  std::size_t l = prob->l;
  int* perm = Malloc(int, l);
  int nr_class;
  // stratified cv may not give leave-one-out rate
  // Each class to l folds -> some folds may have zero elements
  if ((param->svm_type == C_SVC || param->svm_type == NU_SVC) && nr_fold
      < static_cast<int>(l)) {
    int* start = NULL;
    int* label = NULL;
    int* count = NULL;
    svm_group_classes(prob, &nr_class, &label, &start, &count, perm);
    // random shuffle and then data grouped by fold using the array perm
    int* fold_count = Malloc(int, static_cast<std::size_t>(nr_fold));
    int c;
    int* index = Malloc(int, l);
    for (i = 0; i < static_cast<int>(l); i++) {
      index[i] = perm[i];
    }
    for (c = 0; c < nr_class; c++)
      for (i = 0; i < count[c]; i++) {
        int j = i + static_cast<int>(
          PseudoRandom::lcg_rand() % static_cast<unsigned long>(count[c] - i));
        swap(index[start[c] + j], index[start[c] + i]);
      }
    for (i = 0; i < nr_fold; i++) {
      fold_count[i] = 0;
      for (c = 0; c < nr_class; c++) {
        fold_count[i] += (i + 1) * count[c] / nr_fold - i * count[c]
            / nr_fold;
      }
    }
    fold_start[0] = 0;
    for (i = 1; i <= nr_fold; i++) {
      fold_start[i] = fold_start[i - 1] + fold_count[i - 1];
    }
    for (c = 0; c < nr_class; c++)
      for (i = 0; i < nr_fold; i++) {
        int begin = start[c] + i * count[c] / nr_fold;
        int end = start[c] + (i + 1) * count[c] / nr_fold;
        for (int j = begin; j < end; j++) {
          perm[fold_start[i]] = index[j];
          fold_start[i]++;
        }
      }
    fold_start[0] = 0;
    for (i = 1; i <= nr_fold; i++) {
      fold_start[i] = fold_start[i - 1] + fold_count[i - 1];
    }
    free(start);
    free(label);
    free(count);
    free(index);
    free(fold_count);
  } else {
    for (i = 0; i < static_cast<int>(l); i++) {
      perm[i] = i;
    }
    for (i = 0; i < static_cast<int>(l); i++) {
      int j = i + static_cast<int>(
        PseudoRandom::lcg_rand() % static_cast<unsigned long>(static_cast<int>(l) - i));
      swap(perm[i], perm[j]);
    }
    for (i = 0; i <= nr_fold; i++) {
      fold_start[i] = i * static_cast<int>(l) / nr_fold;
    }
  }
  for (i = 0; i < nr_fold; i++) {
    int begin = fold_start[i];
    int end = fold_start[i + 1];
    int j, k;
    struct svm_problem subprob;
    subprob.l = static_cast<std::size_t>(static_cast<int>(l) - (end - begin));
#ifdef _DENSE_REP
    subprob.x = Malloc(struct svm_node, subprob.l);
#else
    subprob.x = Malloc(struct svm_node*, subprob.l);
#endif
    subprob.y = Malloc(double, subprob.l);
    k = 0;
    for (j = 0; j < begin; j++) {
      subprob.x[k] = prob->x[perm[j]];
      subprob.y[k] = prob->y[perm[j]];
      ++k;
    }
    for (j = end; j < static_cast<int>(l); j++) {
      subprob.x[k] = prob->x[perm[j]];
      subprob.y[k] = prob->y[perm[j]];
      ++k;
    }
    struct svm_model* submodel = svm_train(&subprob, param);
    if (param->probability && (param->svm_type == C_SVC || param->svm_type
        == NU_SVC)) {
      double* prob_estimates = Malloc(double, 
        static_cast<std::size_t>(svm_get_nr_class(submodel)));
      for (j = begin; j < end; j++)
#ifdef _DENSE_REP
        target[perm[j]] = svm_predict_probability(submodel, (prob->x
            + perm[j]), prob_estimates);
#else
      target[perm[j]] = svm_predict_probability(submodel, prob->x[perm[j]], prob_estimates);
#endif
      free(prob_estimates);
    } else
      for (j = begin; j < end; j++)
#ifdef _DENSE_REP
        target[perm[j]] = svm_predict(submodel, prob->x + perm[j]);
#else
    target[perm[j]] = svm_predict(submodel, prob->x[perm[j]]);
#endif
    svm_destroy_model(submodel);
    free(subprob.x);
    free(subprob.y);
  }
  free(fold_start);
  free(perm);
}

int svm_get_svm_type(const svm_model* model) {
  return model->param.svm_type;
}

int svm_get_nr_class(const svm_model* model) {
  return model->nr_class;
}

void svm_get_labels(const svm_model* model, int* label) {
  if (model->label != NULL) for (std::size_t i = 0; static_cast<int>(i) < model->nr_class; i++) {
    label[i] = model->label[i];
  }
}

double svm_get_svr_probability(const svm_model* model) {
  if ((model->param.svm_type == EPSILON_SVR || model->param.svm_type
      == NU_SVR) && model->probA != NULL) {
    return model->probA[0];
  } else {
    info("Model doesn't contain information for SVR probability inference\n");
    return 0;
  }
}

void svm_predict_values(const svm_model* model, const svm_node* x,
                        double* dec_values) {
  if (model->param.svm_type == ONE_CLASS || model->param.svm_type
      == EPSILON_SVR || model->param.svm_type == NU_SVR) {
    double* sv_coef = model->sv_coef[0];
    double sum = 0;
    for (std::size_t i = 0; static_cast<int>(i) < model->l; i++)
#ifdef _DENSE_REP
      sum += sv_coef[i] * Kernel::k_function(x,
                                             model->SV + i,
                                             model->param);
#else
    sum += sv_coef[i] * Kernel::k_function(x, model->SV[i], model->param);
#endif
    sum -= model->rho[0];
    *dec_values = sum;
  } else {
    int i;
    int nr_class = model->nr_class;
    std::size_t l = static_cast<std::size_t>(model->l);
    double* kvalue = Malloc(double, l);
    for (i = 0; i < static_cast<int>(l); i++)
#ifdef _DENSE_REP
      kvalue[i] = Kernel::k_function(x, model->SV + i, model->param);
#else
    kvalue[i] = Kernel::k_function(x, model->SV[i], model->param);
#endif
    int* start = Malloc(int, static_cast<std::size_t>(nr_class));
    start[0] = 0;
    for (i = 1; i < nr_class; i++) {
      start[i] = start[i - 1] + model->nSV[i - 1];
    }
    int p = 0;
    for (i = 0; i < nr_class; i++)
      for (int j = i + 1; j < nr_class; j++) {
        double sum = 0;
        int si = start[i];
        int sj = start[j];
        int ci = model->nSV[i];
        int cj = model->nSV[j];
        int k;
        double* coef1 = model->sv_coef[j - 1];
        double* coef2 = model->sv_coef[i];
        for (k = 0; k < ci; k++) {
          sum += coef1[si + k] * kvalue[si + k];
        }
        for (k = 0; k < cj; k++) {
          sum += coef2[sj + k] * kvalue[sj + k];
        }
        sum -= model->rho[p];
        dec_values[p] = sum;
        p++;
      }
    free(kvalue);
    free(start);
  }
}

double svm_predict(const svm_model* model, const svm_node* x) {
  if (model->param.svm_type == ONE_CLASS || model->param.svm_type
      == EPSILON_SVR || model->param.svm_type == NU_SVR) {
    double res;
    svm_predict_values(model, x, &res);
    if (model->param.svm_type == ONE_CLASS) {
      return (res > 0) ? 1 : -1;
    } else {
      return res;
    }
  } else {
    int i;
    std::size_t nr_class = static_cast<std::size_t>(model->nr_class);
    double* dec_values = Malloc(double, nr_class * (nr_class - 1) / 2);
    svm_predict_values(model, x, dec_values);
    int* vote = Malloc(int, nr_class);
    for (i = 0; i < static_cast<int>(nr_class); i++) {
      vote[i] = 0;
    }
    int pos = 0;
    for (i = 0; i < static_cast<int>(nr_class); i++)
      for (int j = i + 1; j < static_cast<int>(nr_class); j++) {
        if (dec_values[pos++] > 0) {
          ++vote[i];
        } else {
          ++vote[j];
        }
      }
    int vote_max_idx = 0;
    for (i = 1; i < static_cast<int>(nr_class); i++)
      if (vote[i] > vote[vote_max_idx]) {
        vote_max_idx = i;
      }
    free(vote);
    free(dec_values);
    return model->label[vote_max_idx];
  }
}

double svm_predict_probability(const svm_model* model, const svm_node* x,
                               double* prob_estimates) {
  if ((model->param.svm_type == C_SVC || model->param.svm_type == NU_SVC)
      && model->probA != NULL && model->probB != NULL) {
    int i;
    std::size_t nr_class = static_cast<std::size_t>(model->nr_class);
    double* dec_values = Malloc(double, nr_class * (nr_class - 1) / 2);
    svm_predict_values(model, x, dec_values);
    double min_prob = 1e-7;
    double** pairwise_prob = Malloc(double*, nr_class);
    for (i = 0; i < static_cast<int>(nr_class); i++) {
      pairwise_prob[i] = Malloc(double, nr_class);
    }
    std::size_t k = 0;
    for (i = 0; i < static_cast<int>(nr_class); i++)
      for (std::size_t j = static_cast<std::size_t>(i) + 1; j < nr_class; j++) {
        pairwise_prob[i][j] = min(max(sigmoid_predict(dec_values[k],
                                                      model->probA[k],
                                                      model->probB[k]),
                                      min_prob), 1 - min_prob);
        pairwise_prob[j][i] = 1 - pairwise_prob[i][j];
        k++;
      }
    multiclass_probability(nr_class, pairwise_prob, prob_estimates);
    int prob_max_idx = 0;
    for (i = 1; i < static_cast<int>(nr_class); i++)
      if (prob_estimates[i] > prob_estimates[prob_max_idx]) {
        prob_max_idx = i;
      }
    for (i = 0; i < static_cast<int>(nr_class); i++) {
      free(pairwise_prob[i]);
    }
    free(dec_values);
    free(pairwise_prob);
    return model->label[prob_max_idx];
  } else {
    return svm_predict(model, x);
  }
}

const char* svm_type_table[] = { "c_svc", "nu_svc", "one_class",
                                 "epsilon_svr", "nu_svr", NULL };

const char* kernel_type_table[] = { "linear", "polynomial", "rbf",
                                    "sigmoid", "precomputed", NULL };

int svm_save_model(const char* model_file_name, const svm_model* model) {
  FILE* fp = fopen(model_file_name, "w");
  if (fp == NULL) {
    return -1;
  }
  const svm_parameter& param = model->param;
  fprintf(fp, "svm_type %s\n", svm_type_table[param.svm_type]);
  fprintf(fp, "kernel_type %s\n", kernel_type_table[param.kernel_type]);
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
  std::size_t l = static_cast<std::size_t>(model->l);
  fprintf(fp, "nr_class %d\n", nr_class);
  fprintf(fp, "total_sv %zu\n", l);
  {
    fprintf(fp, "rho");
    for (std::size_t i = 0; static_cast<int>(i) < nr_class * (nr_class - 1) / 2; i++) {
      fprintf(fp, " %g", model->rho[i]);
    }
    fprintf(fp, "\n");
  }
  if (model->label) {
    fprintf(fp, "label");
    for (std::size_t i = 0; static_cast<int>(i) < nr_class; i++) {
      fprintf(fp, " %d", model->label[i]);
    }
    fprintf(fp, "\n");
  }
  if (model->probA) { // regression has probA only
    fprintf(fp, "probA");
    for (std::size_t i = 0; static_cast<int>(i) < nr_class * (nr_class - 1) / 2; i++) {
      fprintf(fp, " %g", model->probA[i]);
    }
    fprintf(fp, "\n");
  }
  if (model->probB) {
    fprintf(fp, "probB");
    for (std::size_t i = 0; static_cast<int>(i) < nr_class * (nr_class - 1) / 2; i++) {
      fprintf(fp, " %g", model->probB[i]);
    }
    fprintf(fp, "\n");
  }
  if (model->nSV) {
    fprintf(fp, "nr_sv");
    for (std::size_t i = 0; static_cast<int>(i) < nr_class; i++) {
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
  for (std::size_t i = 0; i < l; i++) {
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
  if (ferror(fp) != 0 || fclose(fp) != 0) {
    return -1;
  } else {
    return 0;
  }
}

svm_model* svm_load_model(const char* model_file_name) {
  FILE* fp = fopen(model_file_name, "r");
  if (fp == NULL) {
    return NULL;
  }
  // read parameters
  svm_model* model = Malloc(svm_model, 1);
  svm_parameter& param = model->param;
  model->rho = NULL;
  model->probA = NULL;
  model->probB = NULL;
  model->label = NULL;
  model->nSV = NULL;
  char cmd[81];
  while (1) {
    if(fscanf(fp, "%80s", cmd) == EOF)
    {
      /* do nothing */
    }
    
    if (strcmp(cmd, "svm_type") == 0) {
      if(fscanf(fp, "%80s", cmd) == EOF)
      {
	 /*do nothing*/
      }
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
      
      if(fscanf(fp, "%80s", cmd) == EOF)
      {
	/* do nothing */
      }
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
      if(fscanf(fp, "%d", &param.degree) == EOF)
      {
	/*do nothing*/
      }
    } else if (strcmp(cmd, "gamma") == 0) {
      if(fscanf(fp, "%lf", &param.gamma) == EOF)
      {
	/*do nothing*/
      }
    } else if (strcmp(cmd, "coef0") == 0) {
      if(fscanf(fp, "%lf", &param.coef0) == EOF)
      {
	/*do nothing*/
      }
    } else if (strcmp(cmd, "nr_class") == 0) {
      if(fscanf(fp, "%d", &model->nr_class) == EOF)
      {
	/*do nothing*/
      }
    } else if (strcmp(cmd, "total_sv") == 0) {
      if(fscanf(fp, "%d", &model->l) == EOF)
      {
	/*do nothing*/
      }
    } else if (strcmp(cmd, "rho") == 0) {
      int n = model->nr_class * (model->nr_class - 1) / 2;
      model->rho = Malloc(double, static_cast<std::size_t>(n));
      for (std::size_t i = 0; static_cast<int>(i) < n; i++) {
        if(fscanf(fp, "%lf", &model->rho[i]) == EOF)
	{
	  /*do nothing*/
	}
      }
    } else if (strcmp(cmd, "label") == 0) {
      int n = model->nr_class;
      model->label = Malloc(int, static_cast<std::size_t>(n));
      for (std::size_t i = 0; static_cast<int>(i) < n; i++) {
        if(fscanf(fp, "%d", &model->label[i]) == EOF)
	{
	  /*do nothing*/
	}
      }
    } else if (strcmp(cmd, "probA") == 0) {
      int n = model->nr_class * (model->nr_class - 1) / 2;
      model->probA = Malloc(double, n);
      for (std::size_t i = 0; static_cast<int>(i) < n; i++) {
        if(fscanf(fp, "%lf", &model->probA[i]) == EOF)
	{
	  /*do nothing*/
	}
      }
    } else if (strcmp(cmd, "probB") == 0) {
      int n = model->nr_class * (model->nr_class - 1) / 2;
      model->probB = Malloc(double, n);
      for (std::size_t i = 0; static_cast<int>(i) < n; i++) {
        if(fscanf(fp, "%lf", &model->probB[i]) == EOF)
	{
	  /*do nothing*/
	}
      }
    } else if (strcmp(cmd, "nr_sv") == 0) {
      int n = model->nr_class;
      model->nSV = Malloc(int, n);
      for (std::size_t i = 0; static_cast<int>(i) < n; i++) {
        if(fscanf(fp, "%d", &model->nSV[i]) == EOF)
	{
	  /*do nothing*/
	}
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
  int c_int = 0;
  // read the max dimension of all vectors
  while ((c_int = fgetc(fp)) != EOF) {
    c = (char) c_int;
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
  std::size_t l = static_cast<std::size_t>(model->l);
  model->sv_coef = Malloc(double*, m);
  int i;
  for (i = 0; i < m; i++) {
    model->sv_coef[i] = Malloc(double, l);
  }
#ifdef _DENSE_REP
  model->SV = Malloc(svm_node, l);
  for (i = 0; i < static_cast<int>(l); i++) {
    model->SV[i].values = Malloc(double, elements);
    model->SV[i].dim = 0;
    for (std::size_t k = 0; static_cast<int>(k) < m; k++) {
      if(fscanf(fp, "%lf", &model->sv_coef[k][i]) == EOF)
      {
	/* do nothing */
      }
    }
    int* d = &(model->SV[i].dim);
    while ((c = static_cast<char>(getc(fp))) != '\n') {
      if (!isspace(c)) {
        ungetc(c, fp);
        if(fscanf(fp, "%d:%lf", &index, &value) == EOF)
	{
	  /* do nothing */
	}
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
    for (std::size_t k = 0; k < m; k++) {
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
  if (ferror(fp) != 0 || fclose(fp) != 0) {
    return NULL;
  }
  model->free_sv = 1; // XXX
  return model;
}

void svm_destroy_model(svm_model* model) {
  if (model->free_sv && model->l > 0)
#ifdef _DENSE_REP
    for (std::size_t i = 0; static_cast<int>(i) < model->l; i++) {
      free(model->SV[i].values);
    }
#else
  free((void*)(model->SV[0]));
#endif
  for (std::size_t i = 0; static_cast<int>(i) < model->nr_class - 1; i++) {
    free(model->sv_coef[i]);
  }
  free(model->SV);
  free(model->sv_coef);
  free(model->rho);
  free(model->label);
  free(model->probA);
  free(model->probB);
  free(model->nSV);
  free(model);
}

void svm_destroy_param(svm_parameter* param) {
  free(param->weight_label);
  free(param->weight);
}

const char* svm_check_parameter(const svm_problem* prob,
                                const svm_parameter* param) {
  // svm_type
  int svm_type = param->svm_type;
  if (svm_type != C_SVC && svm_type != NU_SVC && svm_type != ONE_CLASS
      && svm_type != EPSILON_SVR && svm_type != NU_SVR) {
    return "unknown svm type";
  }
  // kernel_type, degree
  int kernel_type = param->kernel_type;
  if (kernel_type != LINEAR && kernel_type != POLY && kernel_type != RBF
      && kernel_type != SIGMOID && kernel_type != PRECOMPUTED) {
    return "unknown kernel type";
  }
  /*if (param->degree < 0) {
    return "degree of polynomial kernel < 0";
  }*/
  // cache_size,eps,C,nu,p,shrinking
  if (param->cache_size <= 0) {
    return "cache_size <= 0";
  }
  if (param->eps <= 0) {
    return "eps <= 0";
  }
  if (svm_type == C_SVC || svm_type == EPSILON_SVR || svm_type == NU_SVR) if (param->C
      <= 0) {
    return "C <= 0";
  }
  if (svm_type == NU_SVC || svm_type == ONE_CLASS || svm_type == NU_SVR) if (param->nu
      <= 0 || param->nu > 1) {
    return "nu <= 0 or nu > 1";
  }
  if (svm_type == EPSILON_SVR) if (param->p < 0) {
    return "p < 0";
  }
  if (param->shrinking != 0 && param->shrinking != 1) {
    return "shrinking != 0 and shrinking != 1";
  }
  if (param->probability != 0 && param->probability != 1) {
    return "probability != 0 and probability != 1";
  }
  if (param->probability == 1 && svm_type == ONE_CLASS) {
    return "one-class SVM probability output not supported yet";
  }
  // check whether nu-svc is feasible
  if (svm_type == NU_SVC) {
    std::size_t l = prob->l;
    std::size_t max_nr_class = 16;
    std::size_t nr_class = 0;
    int* label = Malloc(int, max_nr_class);
    int* count = Malloc(int, max_nr_class);
    int i;
    for (i = 0; static_cast<std::size_t>(i) < l; i++) {
      int this_label = (int)prob->y[i];
      int j;
      for (j = 0; static_cast<std::size_t>(j) < nr_class; j++)
        if (this_label == label[j]) {
          ++count[j];
          break;
        }
      if (static_cast<std::size_t>(j) == nr_class) {
        if (nr_class == max_nr_class) {
          max_nr_class *= 2;
          label = (int*)realloc(label, max_nr_class * sizeof(int));
          count = (int*)realloc(count, max_nr_class * sizeof(int));
        }
        label[nr_class] = this_label;
        count[nr_class] = 1;
        ++nr_class;
      }
    }
    for (i = 0; static_cast<std::size_t>(i) < nr_class; i++) {
      int n1 = count[i];
      for (int j = i + 1; static_cast<std::size_t>(j) < nr_class; j++) {
        int n2 = count[j];
        if (param->nu * (n1 + n2) / 2 > min(n1, n2)) {
          free(label);
          free(count);
          return "specified nu is infeasible";
        }
      }
    }
    free(label);
    free(count);
  }
  return NULL;
}

int svm_check_probability_model(const svm_model* model) {
  return ((model->param.svm_type == C_SVC || model->param.svm_type
      == NU_SVC) && model->probA != NULL && model->probB != NULL)
      || ((model->param.svm_type == EPSILON_SVR || model->param.svm_type
          == NU_SVR) && model->probA != NULL);
}
