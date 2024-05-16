#include "include/math.h"
#include "include/vec.h"
#include "include/io.h"
#include "include/random.h"
#include "vec.h"

#include "include/constants.h"



/* SPLIT : Nicolas */
void split_HFBF(vec all, vec BF, vec HF)
{
     int i, j;
     for (i = 0; i < vec_length( BF ); i++)
     {
         BF[i] = all[i];
     }
     for (j = 0; j < vec_length( HF ); j++, i++)
     {
         HF[j] = all[i];
     }
}


void vec_rand_bin( vec v, int n )
{
  if ( n < vec_length(v) )
  {
       int tmp, cmp;
       cmp = 0;
       while (cmp < n)
       {
         tmp = mt19937_rand_int32() % vec_length(v);
         if (tmp < 0)
         {
           tmp = -tmp;
         }
         if (v[tmp] != 0) 
         {
           cmp++;          
         }
         v[tmp] = 1;
       }
  }
  else
  {
      printf("ERREUR::Dimension des sous-espaces trop grande");
  }
}



void vec_fill_ones(vec v, int dim_SP){
     int i, j;
     int d = (int)sqrt(dim_SP);
     vec_full_ones(v);
     for (i = 0; i < d; i++)
     {
         for (j = 0; j < d; j++)
         {
             int k = (i*dim_SP)/4+j;
             if (k < vec_length(v)) v[(i*dim_SP)/4+j] = 0;
         }
     }
}

void vec_full_ones(vec v)
{
     int i;
     for (i = 0; i < vec_length(v); i++)
     {
         v[i] = 1;
     }
}    

void vec_show(vec v){
     int i, N;
     printf("( ");
     N = vec_length(v);
     for (i = 0; i < N-1; i++)
     {
         printf("%.3f, ", v[i]);
     }
     printf("%.3f )", v[N-1]);
     printf("\nNorme : %.3f", vec_norm(v, 1));
     printf("\nNorme2 : %.3f", vec_norm(v, 2));
     printf("\nSQRT(PS(v,v)) : %.3f", sqrt(vec_inner_product(v, v)));
}

void vec_norm_patern(vec p)
{
     double mean = 0;
     int i;
     for (i = 0; i < vec_length(p); i++)
         mean += p[i];
     mean /= vec_length(p);
     
     for (i = 0; i < vec_length(p); i++)
         p[i] -= mean;
     
     double std = 0;
     for (i = 0; i < vec_length(p); i++)
         std += p[i]*p[i];
     
     std /= vec_length(p);
     std = sqrt(std);
     
     if (std > 0)
     {
        for (i = 0; i < vec_length(p); i++)
            p[i] /= std;
     }
}