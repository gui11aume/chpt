#include <pthread.h>

#ifndef TADBIT_LOADED
#define TADBIT_LOADED

#define TOLERANCE 1e-6
#define MAXITER 10000


typedef struct {
   const int n;
   const int m;
   const double **x;
   const char *skip;
   double *llikmat;
   const int verbose;
} llworker_arg;

typedef struct {
   const int n;
   const double *llikmat;
   double *old_llik;
   double *new_llik;
   int nbreaks;
   int *new_bkpt_list;
   const int *old_bkpt_list;
} dpworker_arg;



// 'chpt' output struct.
typedef struct {
   int maxbreaks;
   int nbreaks_opt;
   //int *passages;
   double *llikmat;
   double *mllik;
   int *bkpts;
} chpt_output;


void
chpt(
  /* input */
  const double **x,
  int n,
  const int m,
  int n_threads,
  const int verbose,
  /* output */
  chpt_output *seg
);
#endif
