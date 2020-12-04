#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <sys/time.h>


int main() {
  // variables to time the different steps of fftw
  struct timeval t1, t2, t3;

  // number of times of execution
  int M0=MRUN;

  double pi = 3.141592653589793;

  // dataset size for 2D
#if LONG >= 1
  int imax=65536;
  int kmax=1024;
#else
  int imax=8;
  int kmax=8;
#endif

  // geometry parameters
  double Lz  = 2.0*pi;
  double dz = Lz/kmax;
  double dzk = 1.0/dz;
  double z;

  // run variables
  int i,k,m;

  // error compared to analytical solution
  double error;

  // fftw stuff
  fftw_plan fftw_f, fftw_b;
    
  // allocate arrays
  double  *in  = (double *) malloc(imax*kmax*sizeof(double));
  double *out  = (double *) malloc(imax*kmax*sizeof(double));

  double *inA  = malloc(kmax*sizeof(double));
  double *zrt  = malloc(kmax*sizeof(double));

  // fftw kind parameter
  fftw_r2r_kind *kind;
  kind = (fftw_r2r_kind*) fftw_malloc(sizeof(fftw_r2r_kind) * 1);


  // definition of fftw_plan_many_r2r
  /* fftw_plan fftw_plan_many_r2r(int rank, const int *n, int howmany, */
  /*                          double *in, const int *inembed, */
  /*                          int istride, int idist, */
  /*                          double *out, const int *onembed, */
  /*                          int ostride, int odist, */
  /*                          const fftw_r2r_kind *kind, unsigned flags); */


  ///
  /// Plan creation
  ///
  gettimeofday(&t1, NULL);
  // which fftew-plan to use
  
  printf("fftw_plan_r2r_2d   ");
  fftw_f = fftw_plan_r2r_2d(kmax,imax,  in, out, FFTW_R2HC,FFTW_R2HC, FFTW_MEASURE);
  fftw_b = fftw_plan_r2r_2d(kmax,imax, out,  in, FFTW_HC2R,FFTW_HC2R, FFTW_MEASURE);

  gettimeofday(&t2, NULL);

  // number of executions
  for(m=0;m<M0;m++) {
    
    // set initial values
    for(i = 0;i <imax;i++){
      for(k = 0;k <kmax;k++){
	z=(0.5+k)*dz;
	in[i*kmax+k]  =  sin(z);
	out[i*kmax+k] =  sin(z);

	// analytical solution
	inA[k] = -sin(z);
      }
    }

    // eigenvalues

 /* http://www.fftw.org/fftw3_doc/Real_002dto_002dReal-Transform-Kinds.html */
/*      FFTW_R2HC computes a real-input DFT with output in “halfcomplex” format, i.e. real and imaginary parts for a transform of size n stored as: */

/* r0, r1, r2, ..., rn/2, i(n+1)/2-1, ..., i2, i1 */
/* (Logical N=n, inverse is FFTW_HC2R.) */
/* FFTW_HC2R computes the reverse of FFTW_R2HC, above. (Logical N=n, inverse is FFTW_R2HC.)  */
    
    zrt[0]      = -2.0*dzk*dzk+2.0*dzk*dzk*cos(0.*2.0*pi/((double)kmax));
    zrt[kmax/2] = -2.0*dzk*dzk+2.0*dzk*dzk*cos((0.+(double)kmax/2)*2.0*pi/((double)kmax));

    for(k = 1;k <kmax/2;k++) {
      zrt[k]=-2.0*dzk*dzk+2.0*dzk*dzk*cos(( ( ((double)k) ) )*2.0*pi/((double)kmax));
      zrt[kmax-1-k+1]  = zrt[k];
      
    }


    // execute fftw-plan  ----> forward
    fftw_execute(fftw_f);

    // to solve poisson equation -> divide by eigenvalues
    for( i=0;i<imax;i++) {
      for( k=0;k<kmax;k++) {
	if(zrt[k] == 0.) {
	  out[i*kmax+k]=0.0;
	}
	else {
	  out[i*kmax+k]=out[i*kmax+k] / zrt[k];
	}
      }
    }
    
    // execute fftw-plan  ----> backward
    fftw_execute(fftw_b);

    //check error compred to analytical solution
    error=0.0;
    for(k = 0;k <kmax;k++){
      for(i = 0;i <imax;i++){
	in[i*kmax+k]/=imax*kmax;
	error=error+fabs(in[i*kmax+k]-inA[i]);
      }
    }

  }
  // time end of execution
  gettimeofday(&t3, NULL);

  // compute & display runtime
  long Dt1 = (((t2.tv_sec - t1.tv_sec) * 1000000) + t2.tv_usec) - (t1.tv_usec);
  long Dt2 = (((t3.tv_sec - t2.tv_sec) * 1000000) + t3.tv_usec) - (t2.tv_usec);

  printf("Dt plan create %10ld mus execute %10ld mus\n",Dt1,Dt2);

#if VERBOSE >=1
  // output  
  for(k = 0;k <kmax;k++){
    z=(0.5+k)*dz;
    printf("%10.4f %10.4f %10.4f %10.4f %10.4f\n",z,in[k+0],in[k+kmax*(imax/2)],in[k+kmax*(imax-1)],inA[k]);
    } 
#endif


  // cleanup    
  fftw_destroy_plan(fftw_f);
  fftw_destroy_plan(fftw_b);
  fftw_cleanup();


  free(inA);
  free(zrt);


  free(in);
  free(out);

  free(kind);

  return 0;
}
