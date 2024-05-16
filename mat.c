

#include "include/mat.h"
#include "include/math.h"
#include "include/io.h"
#include "include/random.h"
#include "include/cplx.h"
#include "mat.h"

#include <math.h>
mat filtrage(mat I)
{
     /*
        1 2 1
        2 4 2
        1 2 1
     */
     
     mat X = mat_clone(I);
     int nH = mat_height(I);
     int nW = mat_width(I);
     for (int i=1; i < nH-1; i++){
	     for (int j=1; j < nW-1; j++){
       		X[i][j] = (int)((
               I[i-1][j-1] + 2*I[i-1][j]   + I[i-1][j+1] +
               2*I[i][j-1] + 4*I[i][j]     + 2*I[i][j+2]   +
			   I[i+1][j-1] + 2*I[i+1][j]   + I[i+1][j+1] )/16);	
     	 }
     }
     return X;
        
     
}