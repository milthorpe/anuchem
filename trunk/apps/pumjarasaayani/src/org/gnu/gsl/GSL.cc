#include "GSL.h"

namespace org {
   namespace gnu {   
      namespace gsl {
         int GSLWrapper::eigenSymmv(double *A, double *evec, double *eval, int aSize) {
             gsl_matrix_view aView;
             gsl_matrix_view evecView;
             gsl_vector_view evalView;
             gsl_eigen_symmv_workspace* workspace;

             aView     = gsl_matrix_view_array(A, aSize, aSize);
             evecView  = gsl_matrix_view_array(evec, aSize, aSize);
             evalView  = gsl_vector_view_array(eval, aSize);
             workspace = gsl_eigen_symmv_alloc(aSize);
             
             int returnValue = gsl_eigen_symmv(&(aView.matrix), &(evalView.vector), &(evecView.matrix), workspace);

             gsl_eigen_symmv_sort(&(evalView.vector), &(evecView.matrix), GSL_EIGEN_SORT_VAL_ASC);

             gsl_eigen_symmv_free(workspace);

             return returnValue;
         }

         int GSLWrapper::solve(double *A, double *b, double *x, int aSize) {
             gsl_matrix_view aView = gsl_matrix_view_array(A, aSize, aSize);
             gsl_vector_view bView = gsl_vector_view_array(b, aSize);
             gsl_vector_view xView = gsl_vector_view_array(x, aSize);
             gsl_permutation * p = gsl_permutation_alloc(aSize);

             int info;

             gsl_linalg_LU_decomp(&(aView.matrix), p, &info);
             gsl_linalg_LU_solve(&(aView.matrix), p, &(bView.vector), &(xView.vector));

             gsl_permutation_free(p);

             return info;
         }
      } // gsl
   } // gnu
} // org


