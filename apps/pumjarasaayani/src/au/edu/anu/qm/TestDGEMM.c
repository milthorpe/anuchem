#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

int main(int argc, char* argv[]) {
   int N=700,nOrbital=25,roK=196;
   int i,j;

   double (*K)[N],(*auxIntMat)[N],(*halfAuxMat)[N*roK],(*mos)[N];

   K = (double (*)[N]) malloc(sizeof(double)*N*N);
   auxIntMat = (double (*)[N]) malloc(sizeof(double)*N*N*roK);
   halfAuxMat = (double (*)[N*roK]) malloc(sizeof(double)*N*nOrbital*roK);
   mos = (double (*)[N]) malloc(sizeof(double)*N*nOrbital);

   for (i=0; i<N; i++) for (j=0; j<N; j++)
       K[i][j]=0.;

   for (i=0; i<N*roK; i++) for (j=0; j<N; j++)
       auxIntMat[i][j]=(double)i+j;

// for (i=0; i<nOrbital; i++)  for (j=0; j<N*roK; j++)
//     halfAuxMat[i][j]=(double)i+j;

   for (i=0; i<nOrbital; i++) for (j=0; j<N; j++)
       mos[i][j]=(double)i+j;

   printf("Memalloc and array initialisation completed!\n");
   clock_t start;
   double alpha,beta;
   int LDA,LDB,LDC,m,n,k;

   // printf("mkl_get_max_threads()=%d\n",mkl_get_max_threads());
   // mkl_set_num_threads(1); segmetation fault http://software.intel.com/en-us/forums/topic/299365
   // MKL_Set_Num_Threads(8);
   // printf("mkl_get_max_threads()=%d\n",mkl_get_max_threads());

for (i=1; i<=4; i++) {

   printf("blas_get_cpu_number()=%d\n",blas_get_cpu_number());
   printf("get_num_nodes()=%d\n",get_num_nodes());
   //printf("omp_get_num_procs()=%d\n",omp_get_num_procs());
   goto_set_num_threads(i);
   printf("blas_get_cpu_number()=%d\n",blas_get_cpu_number());
   printf("get_num_nodes()=%d\n",get_num_nodes());


   // Because of Fortran vs C convention DGEMM call will result in C = B x A
   // The LDx variables should be # Columns of corresponding matrix
   // matrix C [n x m]

   // halfAuxMat[nOrbital x (N*roK)] = mos[nOrbital x N] x auxIntMat[(N*roK) x N]^T
   alpha=1.;beta=0.; LDA=N ; LDB=N; LDC=N*roK; m=N*roK; n=nOrbital; k=N;
   start=clock();
   dgemm_("T", "N", &m, &n, &k, &alpha, auxIntMat, &LDA, mos, &LDB, &beta, halfAuxMat, &LDC);
   printf("1st DGEMM completed! %f s\n",(double)(clock()-start)/CLOCKS_PER_SEC);

   double (*halfAuxMat2)[N]= halfAuxMat;

   // K[N x N] += halfAuxMat[(nOrbital*roK) x N]^T x halfAuxMat[(nOrbital*roK) x N]
   alpha=1.;beta=1.; LDA=N ; LDB=N; LDC=N; m=N; n=N; k=nOrbital*roK;
   start=clock();
   dgemm_("N", "T", &m, &n, &k, &alpha, halfAuxMat2, &LDA, halfAuxMat2, &LDB, &beta, K, &LDC);
   printf("2nd DGEMM completed! %f s\n",(double)(clock()-start)/CLOCKS_PER_SEC);
}
   return 0;
}

