#include "../include/mat.h"
#include "../include/vec.h"
#include "../include/parser.h"
#include "../include/separable2D.h"
#include "../include/wavelet.h"
#include "../include/transform2D.h"
#include "../include/wavelet2D.h"
#include "../include/io.h"
#include "../include/random.h"
#include "../include/distance.h"

#include "../include/extract.h"
#include "../include/project.h"
#include "../include/utils.h"
#include "../include/constants.h"

#include <iostream>
using namespace std;

#define LEVELS   4
#define PAS      200

double val_abs(double a);
char* bin2char(int* bin);
int oct2dec(int* oct);
int bin2Lmsg(int* bin);





/* PROGRAMME PRINCPAL ******************************************************* */
int main (int argc, char **argv)
{

  /* Arguments */
  char *inputFile;
  unsigned int key[4];
  unsigned int key1, key2, key3, key4;
  
  key1 = 0; 
  key2 = 0; 
  key3 = 0; 
  key4 = 0;
  
  key[0] = key1; 
  key[1] = key2; 
  key[2] = key3; 
  key[3] = key4;

  mat I_X = NULL;                   /* Matrice image à tatouer */
  mat I_Y = NULL;                   /* Matrice image tatouée */
  unsigned int h_I, w_I, h, w;      /* Dimension image */
  
  /* Parametres d'initialisation */
  int dim_SP, dim_BF, dim_HF, dim;  /* Dimension des espaces */
  vec L1, SP, HF, BF, image;        /* Image decomposée */

      
      
  /*************************************************
   *  Lecture de l'image / Initialisation des porteuses
   *************************************************/
  
   inputFile= argv[1];                              /* Fichier image a tatouer */
   mt19937_srand_by_array(key, 4);                  /* Initialisation des fonctions aléatoires */
    
   I_X = mat_pgm_read(inputFile);                   /* Représentation matricielle dans le domaine spacial */
                                                    // mat_pgm_write ("IMAGE_ORIGINALE.pgm", I_X);
   h_I = mat_height(I_X);                           /* Hauteur de l'image */
   w_I = mat_width(I_X);                            /* Largeur de l'image */

   dim = (h_I * w_I);                               /* Dimension de l'image */
   dim_BF = dim / ((int) pow (2, (2 * 1)));         /* Dimension de l'espace BF */
   dim_HF = dim - dim_BF;                           /* Dimension de l'espace HF */
  
  
  
  /***************************************************
   *  Décomposition dans le domaine des ondelettes   *
   ***************************************************/
   
   mat Wav_X = NULL;
   it_wavelet2D_t *wavelet2D = NULL;
  
   wavelet2D = it_wavelet2D_new (it_wavelet_lifting_97, LEVELS); /* Caractérisation du domaine ondelette */
   Wav_X = it_wavelet2D_transform (wavelet2D, I_X);              /* Décomposition dans le domaine ondelette */ 
  
  
  
  /***************************************************
   *  Décomposition de l'image en 2 sous bandes BF/HF*
   ***************************************************/
  
   L1 = vec_new_zeros (dim_BF);        
   HF = vec_new_zeros (dim_HF);        /* Vecteur du domaine HF : Couche 2 */
   extract (Wav_X, L1, HF, 1);         /* Séparation de l'image en deux domaines BF-HF */
  
  
  
  /***************************************************
   *   Traitement de la couche BF                    * 
   ***************************************************/

   int i, j, k;
   w = w_I / (int)pow(2, LEVELS);
   h = h_I / (int)pow(2, LEVELS);
   
   dim_SP = w * h;
   dim_BF = dim_BF - dim_SP;
   
   SP = vec_new_zeros(dim_SP);               /* Vecteur du domaine SP : Couche 0 */
   BF = vec_new_zeros(dim_BF);               /* Vecteur du domaine BF : Couche 1 */
   
   mat I2 = NULL;
   I2 = vec_to_mat( L1, w_I/2 );
   extract (I2, SP, BF, LEVELS - 1);         /* Séparation de l'image en deux domaines BF-HF */



  /***************************************************
   *   Initialisation des porteuses                  * 
   ***************************************************/

   int nb_bits = 0;
   printf("Nombre de bits a detecter : ");
   scanf("%d", &nb_bits);
   
   int* mot;
   mot = new int[nb_bits];
   
   // Initialisation du mot a inserer
   for (i = 0; i < nb_bits; i++)
   {
       mot[i] = -1;
   }
   
   vec* porteuses_BF;
   porteuses_BF = new vec[nb_bits];
   
   // Initialisation porteuses aleatoire
   for (i = 0; i < nb_bits; i++)
   {
       porteuses_BF[i] = vec_new_zeros(dim_BF);
       vec_randn(porteuses_BF[i]);
   }

   double produit = 0;
   double pas = PAS;
   
   int cellule = 0;
   unsigned int bit = 0;
   
   vec tmp = vec_new_zeros(dim_BF);
   for (i = 0; i < nb_bits; i++)
   {
       vec_normalize(porteuses_BF[i], 2);

       // Produit scalaire avec l'image
       produit = vec_inner_product(BF, porteuses_BF[i]);

       // Cellule dans laquelle je suis
       cellule = (int)floor(produit / pas);
       
       // Bit codant
       bit = (((int)val_abs(cellule)) % 2);
       mot[i] = bit;
  }
  
  cout << endl << "INFORMATION DETECTION" << endl;
  
  char* motInv;
  motInv = new char[bin2Lmsg(mot)];
  motInv = bin2char(mot);
  
  cout << endl << "Message (BF) : ";
  for (i = 0; i < bin2Lmsg(mot); i++)
      cout << motInv[i];
    
  cout << endl; 
  
  
  
    //**************************************************
    //   Initialisation de la clef HF                  * 
    //**************************************************

    // Parametres de generation de la clef
    key[0] = 1 - key1;
    key[1] = 1 - key2; 
    key[2] = 1 - key3; 
    key[3] = 1 - key4;
    
    // Initialisation des fonctions pseudo aleatoires
    mt19937_srand_by_array(key, 4);
        


    //**************************************************
    //   Initialisation des porteuses HF               * 
    //**************************************************
    
    vec* porteuses_HF;
    porteuses_HF = new vec[nb_bits];

    tmp = vec_new_zeros(dim_HF);
    for (i = 0; i < nb_bits; i++)
    {
        porteuses_HF[i] = vec_new_zeros(dim_HF);
        vec_randn(porteuses_HF[i]);
        vec_normalize(porteuses_HF[i], 2);       
    }

   for (i = 0; i < nb_bits; i++)
   {
       mot[i] = -1;
       vec_normalize(porteuses_HF[i], 2);

       // Produit scalaire avec l'image
       produit = vec_inner_product(HF, porteuses_HF[i]);

       // Cellule dans laquelle je suis
       cellule = (int)floor(produit / pas);
       
       // Bit codant
       bit = (((int)val_abs(cellule)) % 2);
       mot[i] = bit;
  }

  motInv = bin2char(mot);
  
  // TODO Ecrire dans un fichier
  cout << "Message (HF) : ";
  for (i = 0; i < bin2Lmsg(mot); i++)
      cout << motInv[i];
    
  cout << endl; 
        
  return (0);
}
/**************************************************************************** */









/* FONCTIONS **************************************************************** */
double val_abs(double a)
{
   if (a < 0)
      a = -a;
   
   return a;
}

char* bin2char(int* bin)
{
     int k = 0;
     int i;
     char* s = new char[bin2Lmsg(bin)-1];
     //cout << "Longueur msg : " << bin2Lmsg(bin) << endl;
     int* tmp = new int[8];
     
     while ( k < bin2Lmsg(bin) )
     {
           for (i = 8*k; i < 8*(k+1); i++)
           {
               tmp[i-8*k] = bin[i];
           }
           s[k] = oct2dec(tmp);
           // cout << k << " : " << s[k] << endl;
           k++;
     }
     
     return s;
}

int oct2dec(int* oct)
{
    int dec = 0;
    int* tmp = new int[8];
    
    int i;
    for (i = 0; i < 8; i++)
    {
        tmp[i] = oct[i];
    }
    
    int b = 1;
    for (i = 0; i < 8; i++)
    {
        dec = dec + b * tmp[7-i];
        b = 2 * b;
    }
    
    return dec;
}

int bin2Lmsg(int* bin)
{
     int k = 0;
     int* tmp = new int[8];
     
     while ( bin[8*k] == 0 || bin[8*k] == 1 )
     {
           k++;
     }
     
     return k;
} 


