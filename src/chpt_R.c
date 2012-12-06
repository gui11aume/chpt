#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "chpt.h"

// Declare and register R/C interface.
SEXP
chpt_R_call(
  SEXP list,
  SEXP n_threads,
  SEXP verbose
);

R_CallMethodDef callMethods[] = {
   {"chpt_R_call", (DL_FUNC) &chpt_R_call, 3},
   {NULL, NULL, 0}
};

void R_init_chpt(DllInfo *info) {
   R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}


SEXP
chpt_R_call(
  SEXP list,
  SEXP n_threads,
  SEXP verbose
){

   R_len_t i, len, m = length(list);
   int first = 1, N, *dim_int;

   SEXP dim;
   PROTECT(dim = allocVector(INTSXP, 2));

   // Convert 'obs_list' to pointer of pointer to double.
   double **x = (double **) malloc(m * sizeof(double **));
   for (i = 0 ; i < m ; i++) {
      // This fails if list element is not numeric.
      x[i] = REAL(coerceVector(VECTOR_ELT(list, i), REALSXP));
      // Check that input is a matrix.
      // Check the length.
      len = length(VECTOR_ELT(list, i));
      if (first) {
         N = (int) len;
         first = 0;
      }
      else {
         if (N != (int) len) {
            error("all vectors must have same length");
         }
      }
   }

   UNPROTECT(1);

   chpt_output *seg = (chpt_output *) malloc(sizeof(chpt_output));
   
   // Call 'chpt'.
   chpt(x, N, m, INTEGER(n_threads)[0], INTEGER(verbose)[0], seg);

   int maxbreaks = seg->maxbreaks;

   // Copy output to R-readable variables.
   SEXP nbreaks_SEXP;
   SEXP llikmat_SEXP;
   SEXP mllik_SEXP;
   SEXP bkpts_SEXP;

   PROTECT(nbreaks_SEXP = allocVector(INTSXP, 1));
   PROTECT(llikmat_SEXP = allocVector(REALSXP, N*N));
   PROTECT(mllik_SEXP = allocVector(REALSXP, maxbreaks));
   PROTECT(bkpts_SEXP = allocVector(INTSXP, N*(maxbreaks-1)));

   int *nbreaks_opt = INTEGER(nbreaks_SEXP); 
   double *llikmat = REAL(llikmat_SEXP);
   double *mllik = REAL(mllik_SEXP);
   int *bkpts = INTEGER(bkpts_SEXP);

   nbreaks_opt[0] = seg->nbreaks_opt;
   for (i = 0 ; i < N*N ; i++) llikmat[i] = seg->llikmat[i];
   for (i = 0 ; i < maxbreaks ; i++) mllik[i] = seg->mllik[i];
   // Remove first column associated with 0 breaks. Itcontains only
   // 0s and shifts the index in R (vectors start at position 1).
   for (i = N ; i < N*(maxbreaks-1) ; i++) bkpts[i-N] = seg->bkpts[i];


   // Set 'dim' attributes.
   SEXP dim_llikmat;
   PROTECT(dim_llikmat = allocVector(INTSXP, 2));
   INTEGER(dim_llikmat)[0] = N;
   INTEGER(dim_llikmat)[1] = N;
   setAttrib(llikmat_SEXP, R_DimSymbol, dim_llikmat);

   SEXP dim_breaks;
   PROTECT(dim_breaks = allocVector(INTSXP, 2));
   INTEGER(dim_breaks)[0] = N;
   INTEGER(dim_breaks)[1] = maxbreaks-1;
   setAttrib(bkpts_SEXP, R_DimSymbol, dim_breaks);

   free(x);
   free(seg->llikmat);
   free(seg->mllik);
   free(seg->bkpts);
   free(seg);

   SEXP list_SEXP;
   PROTECT(list_SEXP = allocVector(VECSXP, 4));
   SET_VECTOR_ELT(list_SEXP, 0, nbreaks_SEXP);
   SET_VECTOR_ELT(list_SEXP, 1, llikmat_SEXP);
   SET_VECTOR_ELT(list_SEXP, 2, mllik_SEXP);
   SET_VECTOR_ELT(list_SEXP, 3, bkpts_SEXP);
   UNPROTECT(7);

   return list_SEXP;

}
