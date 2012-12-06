#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <ctype.h>
#include <unistd.h>
#include <float.h>

#include "chpt.h"


// Global variables. //

int n_processed;              // Number of slices processed so far.
int n_to_process;             // Total number of slices to process.
int taskQ_i;                  // Index used for task queue.
pthread_mutex_t chpt_lock;    // Mutex to access task queue.

double
ll(
  const int    n,
  const double *x
){
// SYNOPSIS:                                                            
//   Find the log-likelihood.
//                                                                      
// ARGUMENTS:                                                           
//   'n': length of the array 'x'.
//   'x': values of the array.
//                                                                      
// RETURN:                                                              
//   The maximum log-likelihood of a segment.                          
//                                                                      

   int i;
   double sum = 0.0;
   double sumsq = 0.0;

   for (i = 0 ; i < n ; i++) {
      sum += x[i];
      sumsq += x[i]*x[i];
   }

   return sum*sum/n - sumsq;

}

void *
fill_DP(
  void *arg
){
// SYNOPSIS:                                                            
//   Thread function to compute the values of the dynamic programming   
//   'DPwalk' for a given number of breaks (value of 'nbreaks').        
//                                                                      
// PARAMETERS:                                                          
//   'arg': arguments for a thread (see header file).                   
//                                                                      
// RETURN:                                                              
//   'void *'                                                           
//                                                                      
// SIDE-EFFECTS:                                                        
//   Update 'old_llik', 'new_llik' and 'new_bkpt_list' in place.        
//                                                                      

   dpworker_arg *myargs = (dpworker_arg *) arg;
   const int n = myargs->n;
   const double *llikmat = (const double *) myargs->llikmat;
   double *old_llik = (double *) myargs->old_llik;
   double *new_llik = (double *) myargs->new_llik;
   const int nbreaks = myargs->nbreaks;
   int *new_bkpt_list = (int *) myargs->new_bkpt_list;
   const int *old_bkpt_list = (const int *) myargs->old_bkpt_list;

   int i;

   while (1) {
      pthread_mutex_lock(&chpt_lock);
      if (taskQ_i > n-1) {
         // Task queue is empty. Exit loop and return
         pthread_mutex_unlock(&chpt_lock);
         break;
      }
      // A task gives an end point 'j'.
      int j = taskQ_i;
      taskQ_i++;
      pthread_mutex_unlock(&chpt_lock);

      new_llik[j] = -INFINITY;
      int new_bkpt = -1;

      // Cycle over start point 'i'.
      for (i = 2 * nbreaks ; i < j ; i++) {

         // If NAN the following condition evaluates to false.
         double tmp = old_llik[i-1] + llikmat[i+j*n];
         if (tmp > new_llik[j]) {
            new_llik[j] = tmp;
            new_bkpt = i-1;
         }
      }

      // Update breakpoint list (skip if log-lik is undefined).
      // No need to use mutex because 'j' is different for every thread.
      if (new_llik[j] > -INFINITY) {
         for (i = 0 ; i < n ; i++) {
            new_bkpt_list[j+i*n] = old_bkpt_list[new_bkpt+i*n];
         }
         new_bkpt_list[j+new_bkpt*n] = 1;
      }
   }

   return NULL;

}

void
DPwalk(
  // input //
  const double *llikmat,
  const int n,
  const int MAXBREAKS,
  int n_threads,
  // output //
  double *mllik,
  int *breakpoints
){
// SYNOPSIS:                                                            
//   Dynamic programming algorithm to compute the most likely position  
//   of breakpoints given a matrix of slice maximum log-likelihood.     
//                                                                      
// PARAMETERS:                                                          
//   '*llikmat': matrix of maximum log-likelihood values.               
//   'n': row/col number of 'llikmat'.                                  
//   'MAXBREAKS': The maximum number of breakpoints.                    
//        -- output arguments --                                        
//   '*mllik': maximum log-likelihood of the segmentations.             
//   '*breakpoints': optimal breakpoints per number of breaks.          
//                                                                      
// RETURN:                                                              
//   'void'                                                             
//                                                                      
// SIDE-EFFECTS:                                                        
//   Update 'breakpoints' in place.                                     
//                                                                      

   int i;
   int j;
   int nbreaks;

   //double new_llik[n];
   //double old_llik[n];
   double *new_llik = (double *) malloc(n * sizeof(double));
   double *old_llik = (double *) malloc(n * sizeof(double));

   // Breakpoint lists. The first index (row) is the end of the slice,
   // the second is 1 if this position is an end (breakpoint).
   // These must be allocated from the heap because 'n' can be large.
   int *new_bkpt_list = (int *) malloc(n*n * sizeof(int));
   int *old_bkpt_list = (int *) malloc(n*n * sizeof(int));

   // Initializations.
   // 'breakpoints' is a 'n' x 'MAXBREAKS' array. The first index (row)
   // is 1 if there is a breakpoint at that location, the second index
   // (column) is the number of breakpoints.
   for (i = 0 ; i < n*MAXBREAKS ; i++) breakpoints[i] = 0;

   for (i = 0 ; i < n*n ; i++) {
      new_bkpt_list[i] = 0;
      old_bkpt_list[i] = 0;
   }

   mllik[0] = llikmat[0+(n-1)*n];
   for (i = 1 ; i < MAXBREAKS ; i++) {
      mllik[i] = NAN;
   }

   // Initialize 'old_llik' to the first line of 'llikmat' containing
   // the log-likelihood of segments starting at index 0.
   for (j = 0 ; j < n ; j++) {
      old_llik[j] = llikmat[0+j*n];
      new_llik[j] = -INFINITY;
   }

   int err = pthread_mutex_init(&chpt_lock, NULL);
   if (err) {
      fprintf(stderr, "error initializing mutex (%d)\n", err);
      free(new_llik);
      free(old_llik);
      free(new_bkpt_list);
      free(old_bkpt_list);
      return;
   }

   dpworker_arg arg = {
      .n = n,
      .llikmat = llikmat,
      .old_llik = old_llik,
      .new_llik = new_llik,
      .nbreaks = 1,
      .new_bkpt_list = new_bkpt_list,
      .old_bkpt_list = old_bkpt_list,
   };

   pthread_t *tid = (pthread_t *) malloc(n_threads * sizeof(pthread_t));

   // Dynamic programming.
   for (nbreaks = 1 ; nbreaks < MAXBREAKS ; nbreaks++) {

      arg.nbreaks = nbreaks;
      // Update breakpoint lists.
      for (i = 0 ; i < n*n ; i++) {
         old_bkpt_list[i] = new_bkpt_list[i];
      }
      // Set 'taskQ_i' to mallest end of segment (variable 'j').
      taskQ_i = 2 * nbreaks + 1;

      for (i = 0 ; i < n_threads ; i++) tid[i] = 0;
      for (i = 0 ; i < n_threads ; i++) {
         err = pthread_create(&(tid[i]), NULL, &fill_DP, &arg);
         if (err) {
            fprintf(stderr, "error creating thread (%d)\n", err);
            free(new_llik);
            free(old_llik);
            free(new_bkpt_list);
            free(old_bkpt_list);
            return;
         }
      }

      // Wait for threads to return.
      for (i = 0 ; i < n_threads ; i++) {
         pthread_join(tid[i], NULL);
      }

      // Update full log-likelihoods.
      mllik[nbreaks] = new_llik[n-1];

      // Record breakpoints.
      for (i = 0 ; i < n ; i++) {
         old_llik[i] = new_llik[i];
         breakpoints[i+nbreaks*n] = new_bkpt_list[n-1+i*n];
      }

   }

   free(new_llik);
   free(old_llik);
   free(new_bkpt_list);
   free(old_bkpt_list);

   return;

}

void *
fill_llikmat(
   void *arg
){
// SYNOPSIS:                                                            
//   Compute the log-likelihood of the fragments.
//                                                                      
// PARAMETERS:                                                          
//   'arg': thread arguments (see header file for definition).          
//                                                                      
// RETURN:                                                              
//   'void'                                                             
//                                                                      
// SIDE-EFFECTS:                                                        
//   Update 'llikmat' in place.                                         
//                                                                      

   llworker_arg *myargs = (llworker_arg *) arg;
   const int n = myargs->n;
   const int m = myargs->m;
   const double **x = myargs->x;
   const char *skip = (const char *) myargs->skip;
   double *llikmat = myargs->llikmat;
   const int verbose = myargs->verbose;

   int i;
   int j;
   int l;
   int job_index;
   
   // Break out of the loop when task queue is empty.
   while (1) {

      pthread_mutex_lock(&chpt_lock);
      while ((taskQ_i < n*n) && (skip[taskQ_i] > 0)) {
         // Fast forward to the next job.
         taskQ_i++;
      }
      if (taskQ_i >= n*n) {
         // Task queue is empty. Exit loop and return
         pthread_mutex_unlock(&chpt_lock);
         break;
      }
      job_index = taskQ_i;
      taskQ_i++;
      pthread_mutex_unlock(&chpt_lock);

      // Compute the log-likelihood of fragment '(i,j)'.
      i = job_index % n;
      j = job_index / n;

      // Distinct parts of the array, no lock needed.
      llikmat[i+j*n] = 0.0;
      for (l = 0 ; l < m ; l++) {
         llikmat[i+j*n] += ll(j-i+1, x[l]+i);
      }
      n_processed++;
      if (verbose) {
         fprintf(stderr, "computing likelihood (%0.f%% done)\r",
            99 * n_processed / (float) n_to_process);
      }
   }

   return NULL;

}

void
chpt(
  // input //
  const double **x,
  const int n,
  const int m,
  int n_threads,
  const int verbose,
  // output //
  chpt_output *seg
){

   // Get thread number if set to 0 (automatic).
   if (n_threads < 1) {
      #ifdef _SC_NPROCESSORS_ONLN
         n_threads = (int) sysconf(_SC_NPROCESSORS_ONLN);
      #else
         n_threads = 1;
      #endif
   }

   int err;       // Used for error checking.

   int i;
   int j;

   const int MAXBREAKS = n/2;

   double *llikmat = (double *) malloc(n*n * sizeof(double));
   for (i = 0 ; i < n*n ; i++) llikmat[i] = NAN;

   // 'skip' will contain only 0 or 1 and can be stored as 'char'.
   char *skip = (char *) malloc(n*n * sizeof(char));

   for (j = 0 ; j < n ; j++)
   for (i = 0 ; i < n ; i++)
      skip[i+j*n] = j > i ? 0 : 1;

   // Allocate 'tid'.
   pthread_t *tid = (pthread_t *) malloc(n_threads * sizeof(pthread_t));

   llworker_arg arg = {
      .n = n,
      .m = m,
      .x = x,
      .skip = skip,
      .llikmat = llikmat,
      .verbose = verbose,
   };

   err = pthread_mutex_init(&chpt_lock, NULL);
   if (err) {
      fprintf(stderr, "error initializing mutex (%d)\n", err);
      return;
   }

   n_to_process = 0;
   for (i = 0 ; i < n*n ; i++) n_to_process += (1-skip[i]);
   n_processed = 0;
   taskQ_i = 0;
   
   // Instantiate threads and start running jobs.
   for (i = 0 ; i < n_threads ; i++) tid[i] = 0;
   for (i = 0 ; i < n_threads ; i++) {
      err = pthread_create(&(tid[i]), NULL, &fill_llikmat, &arg);
      if (err) {
         fprintf(stderr, "error creating thread (%d)\n", err);
         return;
      }
   }

   // Wait for threads to return.
   for (i = 0 ; i < n_threads ; i++) {
      pthread_join(tid[i], NULL);
   }
   if (verbose) {
      fprintf(stderr, "computing likelihood (100%% done)\n");
   }

   // The matrix 'llikmat' now contains the log-likelihood of the
   // segments. The breakpoints are found by dynamic programming.
   double *mllik = (double *) malloc(MAXBREAKS * sizeof(double));
   int *bkpts = (int *) malloc(MAXBREAKS*n * sizeof(int));
   DPwalk(llikmat, n, MAXBREAKS, n_threads, mllik, bkpts);

   // Get optimal number of breaks by AIC.
   int n_params;
   int nbreaks_opt;
   double newAIC = -INFINITY;
   for (nbreaks_opt = 1 ; nbreaks_opt < MAXBREAKS ; nbreaks_opt++) {
      n_params = nbreaks_opt + m*(nbreaks_opt+1);
      if (newAIC > mllik[nbreaks_opt] - n_params) break;
      newAIC = mllik[nbreaks_opt] - n_params;
   }


   pthread_mutex_destroy(&chpt_lock);
   free(skip);
   free(tid);

   /*
   // Compute breakpoint confidence by penalized dynamic progamming.
   double *llikmatcpy = (double *) malloc (n*n * sizeof(double));
   double *mllikcpy = (double *) malloc(MAXBREAKS * sizeof(double));
   int *bkptscpy = (int *) malloc(n*MAXBREAKS * sizeof(int));
   int *passages = (int *) malloc(n * sizeof(int));
   for (i = 0 ; i < n*MAXBREAKS ; i++) bkptscpy[i] = bkpts[i];
   for (i = 0 ; i < n*n ; i++) llikmatcpy[i] = llikmat[i];
   for (i = 0 ; i < n ; i++) passages[i] = 0;

   for (l = 0 ; l < 10 ; l++) {
      i = 0;
      for (j = 0 ; j < n ; j++) {
         if (bkptscpy[j+nbreaks_opt*n]) {
            // Apply a constant penalty every time a TAD is present
            // in the final decomposition. The penalty is set to
            // 'm*6' because it is the expected log-likelihood gain
            // for adding a new TAD around the optimum log-likelihood.
            llikmatcpy[i+j*n] -= m*6;
            passages[j] += bkpts[j+nbreaks_opt*n];
            i = j+1;
         }
      }
      if (i < n) llikmatcpy[i+(n-1)*n] -= m*6;
      DPwalk(llikmatcpy, n, nbreaks_opt+1, n_threads, mllikcpy, bkptscpy);
   }
   free(llikmatcpy);
   free(mllikcpy);
   free(bkptscpy);
   */
   
   // Update output struct.
   seg->maxbreaks = MAXBREAKS;
   seg->nbreaks_opt = nbreaks_opt;
   seg->llikmat = llikmat;
   seg->mllik = mllik;
   seg->bkpts = bkpts;

   return;

}
