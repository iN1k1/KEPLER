/* file:        whistcY.c
** description: MEX weighted histc function.
**/


/** @file
 ** @brief HISTC MEX function implementation.
 **/

#include"mex.h"
#include<stdlib.h>
#include<math.h>

/** WHISTC(X,W,EDGES) 
 **/
#define min(a,b) ((a<b)?a:b)
#define max(a,b) ((a<b)?b:a)


/** @brief MEX driver.
 ** @param nout number of MATLAB output arguments.
 ** @param out MATLAB output arguments.
 ** @param nin number of MATLAB input arguments.
 ** @param in MATLAB input arguments.
 **/  
void
mexFunction(int nout, mxArray *out[], 
            int nin, const mxArray *in[])
{
  int M, N, NE ;
  double* Xpt ;
  double* Wpt ; 
  double* EDGESpt ;
  double* RESpt ;
  enum {X=0, W, EDGES} ;

  /** -----------------------------------------------------------------
   **                                               Check the arguments
   ** -------------------------------------------------------------- */
  M = mxGetM(in[X]) ;
  N = mxGetN(in[X]) ;

  
  NE = max(mxGetM(in[EDGES]), mxGetN(in[EDGES])) ;

  
  Xpt = mxGetPr(in[X]) ;
  Wpt = mxGetPr(in[W]) ;
  EDGESpt = mxGetPr(in[EDGES]) ;


  /* If the input is a vector, make it a column. */
  if(M == 1) {
    M = N ; 
    N = 1 ;
  }

  /* Alloc the result. */
  out[0] = mxCreateDoubleMatrix(NE, 1, mxREAL) ;
  RESpt = mxGetPr(out[0]) ; 

  /** -----------------------------------------------------------------
   **                                                        Do the job
   ** -------------------------------------------------------------- */
	{  
     int i = 0;
	  for(; i < M; i++){
			int c = 0;
			while(c < NE){
				if((c==NE-1) || (Xpt[i] < EDGESpt[c+1] && Xpt[i] >= EDGESpt[c])){
					RESpt[c] += Wpt[i];
					c = NE+1;
				}
				else
					c++;
			  }
	  }
     
  }
  return ;
}


