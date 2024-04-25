#include "include/mat.h"
#include "include/vec.h"
#include "include/parser.h"
#include "include/separable2D.h"
#include "include/wavelet.h"
#include "include/transform2D.h"
#include "include/wavelet2D.h"
#include "include/io.h"
#include "include/random.h"
#include "include/distance.h"
#include "include/extract.h"
#include "include/project.h"
#include "include/utils.h"
#include "include/constants.h"

#include <iostream>
using namespace std;

#define LEVELS   4
#define PAS      200

double val_abs(double a);          // Valeur absolue

int longueur( char *chaine );      // Longueur d'un char*
int longueurMsgBin(int* bin);      // Longueur du char* a partir du mot binaire

int* dec2bin(int dec);             // Conversion entier -> binaire
int oct2dec(int* oct);             // Conversion octal -> decimal
char* bin2char(int* bin);          // Conversion binaire -> char*





/* PROGRAMME PRINCPAL ******************************************************* */
int main (int argc, char **argv)
{
    int i, j, k;

    //************************************************
    //  Lecture de l'image                           *
    //************************************************

    char *inputFile;
    inputFile = argv[1];                             // Fichier image a tatouer

    mat I_X = NULL;
    I_X = mat_pgm_read(inputFile);                   // Matrice image

    unsigned int h_I, w_I;
    h_I = mat_height(I_X);                           // Hauteur de l'image
    w_I = mat_width(I_X);                            // Largeur de l'image

    unsigned int h, w;
    w = w_I / (int)pow(2, LEVELS);                   // Hauteur de l'imagette
    h = h_I / (int)pow(2, LEVELS);                   // Largeur de l'imagette

    int dim = h_I * w_I;                             // Dimension de l'image
    int dim_L1 = dim / ((int) pow (2, (2 * 1)));
    int dim_HF = dim - dim_L1;                       // Dimension de l'espace HF
    int dim_SP = w * h;                              // Dimension de l'espace SP
    int dim_BF = dim_L1 - dim_SP;                    // Dimension de l'espace BF

    cout << endl;
    cout << "INFORMATIONS IMAGE" << endl;
    cout << "Taille de l'image : " << h_I << "x" << w_I << endl;
    cout << "Dimension Hautes Frequences : " << dim_HF << endl;
    cout << "Dimension Basses Frequences : " << dim_BF << endl;
    cout << "Dimension de l'Imagette : " << dim_SP << endl << endl;



    //**************************************************
    //  Décomposition dans le domaine des ondelettes   *
    //**************************************************

    it_wavelet2D_t *wavelet2D = NULL;
    wavelet2D = it_wavelet2D_new (it_wavelet_lifting_97, LEVELS); // Caractérisation du domaine ondelette

    mat Wav_X = NULL;
    Wav_X = it_wavelet2D_transform (wavelet2D, I_X);              // Décomposition dans le domaine ondelette
    // mat_pgm_write("IMAGE_ONDELETTES.pgm", Wav_X);

    cout << "INFORMATIONS ONDELETTE" << endl;
    cout << "Niveaux de decomposition : " << LEVELS << endl;
    cout << "Type d'ondelettes : Daubechies 9/7" << endl << endl;



    //**************************************************
    //  Décomposition de l'image en 2 sous bandes BF/HF*
    //**************************************************

    vec L1 = vec_new_zeros(dim_L1);
    vec HF = vec_new_zeros(dim_HF);         // Vecteur du domaine HF : Couche 2
    vec BF = vec_new_zeros(dim_BF);         // Vecteur du domaine BF : Couche 1
    vec SP = vec_new_zeros(dim_SP);         // Vecteur du domaine SP : Couche 0
    extract (Wav_X, L1, HF, 1);             // Séparation de l'image en deux domaines BF-HF

    mat I2 = vec_to_mat( L1, w_I/2 );
    extract (I2, SP, BF, LEVELS - 1);       // Séparation de l'image en deux domaines SP-BF



    //**************************************************
    //   Message a inserer dans l'image                *
    //**************************************************

    // TODO Lire un fichier
    char msg[255];
    cout << "Message : ";
    cin >> msg;

    int L = longueur(msg) - 1;
    int nb_bits = 8*L;

    int* tmpt;
    tmpt = new int[8];

    int* mot;
    mot = new int[nb_bits];

    k = 0;
    for (i = 0; i < L; i++)
    {
        tmpt = dec2bin(msg[i]);
        for (j = 0; j < 8; j++)
        {
            mot[k] = tmpt[j];
            k++;
        }
    }



    //**************************************************
    //   Initialisation de la clef BF                  *
    //**************************************************

    // Parametres de generation de la clef
    unsigned int key[4];
    unsigned int key1 = 0;
    unsigned int key2 = 0;
    unsigned int key3 = 0;
    unsigned int key4 = 0;

    key[0] = key1;
    key[1] = key2;
    key[2] = key3;
    key[3] = key4;

    // Initialisation des fonctions pseudo aleatoires
    mt19937_srand_by_array(key, 4);



    //**************************************************
    //   Initialisation des porteuses BF               *
    //**************************************************

    vec* porteuses_BF;
    porteuses_BF = new vec[nb_bits];
    double produit = 0;

    vec tmp = vec_new_zeros(dim_BF);
    for (i = 0; i < nb_bits; i++)
    {
        porteuses_BF[i] = vec_new_zeros(dim_BF);
        vec_randn(porteuses_BF[i]);

        // Principe d'orthonormalisation de Gram-Schmidt
        /*
        for (j = 0; j < i; j++)
        {
            vec_copy(tmp, porteuses_BF[j]);
            produit = vec_inner_product(porteuses_BF[i], tmp);
            produit /= vec_inner_product(tmp, tmp);

            vec_mul_by(tmp, produit);
            vec_sub(porteuses_BF[i], tmp);
        }
        */

        vec_normalize(porteuses_BF[i], 2);
    }




    //**************************************************
    //   Tatouage dans les BF                          *
    //**************************************************

    double pas = PAS;
    double d1, d2, d;

    int cellule = 0;
    unsigned int bit = 0;

    for (i = 0; i < nb_bits; i++)
    {
        // Produit scalaire avec les coeff basses frequences
        produit = vec_inner_product(BF, porteuses_BF[i]);

        // Cellule dans laquelle je suis
        cellule = (int)floor(produit / pas);

        // Bit codant
        bit = (((int)val_abs(cellule)) % 2);

        if (mot[i] == bit)
        {
           // Centre de cette cellule
           d = (cellule + 0.5) * pas - produit;
        }
        else
        {
           // Choix de la cellule la plus proche
           d1 = val_abs((cellule - 0.5) * pas - produit);
           d2 = val_abs((cellule + 1.5) * pas - produit);

           d = (cellule - 0.5) * pas - produit;
           if (d2 < d1)
           {
              d = (cellule + 1.5) * pas  - produit;
           }
        }

        vec_copy(tmp, porteuses_BF[i]);
        vec_mul_by(tmp, d);
        vec_add(BF, tmp);
    }

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

    // Concatenation de l'ensemble des coeff ondelettes

    vec* porteuses_HF;
    porteuses_HF = new vec[nb_bits];

    tmp = vec_new_zeros(dim_HF);
    for (i = 0; i < nb_bits; i++)
    {
        porteuses_HF[i] = vec_new_zeros(dim_HF);
        vec_randn(porteuses_HF[i]);
        vec_normalize(porteuses_HF[i], 2);
    }



    //**************************************************
    //   Tatouage dans les HF                          *
    //**************************************************

    for (i = 0; i < nb_bits; i++)
    {
        // Produit scalaire avec l'image
        produit = vec_inner_product(HF, porteuses_HF[i]);

        // Cellule dans laquelle je suis
        cellule = (int)floor(produit / pas);

        // Bit codant
        bit = (((int)val_abs(cellule)) % 2);

        if (mot[i] == bit)
        {
           // Centre de la cellule
           d = (cellule + 0.5) * pas - produit;
        }
        else
        {
           // Choix de la cellule la plus proche
           d1 = val_abs((cellule - 0.5) * pas - produit);
           d2 = val_abs((cellule + 1.5) * pas - produit);

           d = (cellule - 0.5) * pas - produit;
           if (d2 < d1)
           {
              d = (cellule + 1.5) * pas  - produit;
           }
        }

        vec_copy(tmp, porteuses_HF[i]);
        vec_mul_by(tmp, d);
        vec_add(HF, tmp);
    }



    //**************************************************
    //   Transformations inverses                      *
    //**************************************************

    // Separation de l'ensemble des coefficients modifies
    // split_HFBF(allHFBF, BF, HF);

    mat Wav_Y = NULL;
    Wav_Y = mat_new_zeros(h_I, w_I);                 // Matrice image tatouée dans les ondelettes
    extractInv(SP, BF, LEVELS-1, h_I/2, w_I/2, I2);
    L1 = mat_to_vec(I2);

    // Transformation inverse
    mat I_Y = NULL;
    extractInv(L1, HF, 1, h_I, w_I, Wav_Y);
    I_Y = it_wavelet2D_itransform(wavelet2D, Wav_Y); // Matrice image tatouée
    double psnr = mat_psnr (I_X, I_Y);               // Calcul du PSNR

    mat_pgm_write("IMAGE_TATOUEE.pgm", I_Y);



    //**************************************************
    //   Affichage informations diverses               *
    //**************************************************

    cout << endl << "INFORMATIONS TATOUAGE" << endl;
    cout << "PSNR = " << psnr << " dB" << endl;
    cout << "Message : ";
    for (i = 0; i < nb_bits; i++)
        cout << mot[i];

    cout << endl << "Nb bits " << nb_bits << endl;
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

int longueur( char *chaine )
{
	int index = 0;
	while( chaine[ index ++ ] );
	return index;
}

int* dec2bin(int dec)
{
    int* bin;
    bin = new int[8];

    int i;
    int tmp = dec;
    for (i=7; i>=0; i--)
    {
        bin[i] = tmp%2;
        tmp = (tmp-(tmp%2))/2;
    }
    return bin;
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

char* bin2char(int* bin)
{
     int k = 0;
     int i;
     char* s = new char[longueurMsgBin(bin)-1];
     int* tmp = new int[8];

     while ( k < longueurMsgBin(bin) )
     {
           for (i = 8*k; i < 8*(k+1); i++)
           {
               tmp[i-8*k] = bin[i];
           }
           s[k] = oct2dec(tmp);
           k++;
     }

     return s;
}

int longueurMsgBin(int* bin)
{
     int k = 0;
     int* tmp = new int[8];

     while ( bin[8*k] == 0 || bin[8*k] == 1 )
     {
           k++;
     }

     return k;
}
/**************************************************************************** */
