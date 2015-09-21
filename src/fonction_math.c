// ** ================================================================================================================
// ** R-Package: graphscan 1.1
// ** Fichier : src/fonction_math.c
// ** Description : détection 1D des clusters avec l'indice de Cucala et celui de Kulldorff.
// **               fonction de normalisation, fonction de calcul d'ordre et fonction beta.
// ** License : GPL-2 | GPL-3
// ** Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
// ** ================================================================================================================

#include "fonction_math.h"
//#include <gsl/gsl_rng.h>
#include <time.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#include <Rmath.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>


/** fonction : normalisation_et_distance_entre_stat_dordre-)
 ** normalise puis calcul les distances entre les evenements
 ** entree :
 ** taille_sequence  taille de la sequence
 ** nbEv  nombre d'evenement
 ** vecteur_x  pointeur sur les evenements non normalise venant du R
 ** entree/sortie :
 ** vecteur_D  pointeur sur les distances entre evenements normalise
 ** tab_P  pointeur sur les indices des distances
 **/
void normalisation_et_distance_entre_stat_dordre(int normalisation_debut, int normalisation_fin,int nbEv,double * vecteur_X, long double * vecteur_D, long double * tab_P)
{
	int i;
	
	#pragma omp parallel for
	for(i=0;i<=nbEv;i++)
	{
		vecteur_D[i] = (long double)((vecteur_X[i+1] - vecteur_X[i])/(normalisation_fin-normalisation_debut));
		tab_P[i] = vecteur_X[i];
	}

}

void distance_entre_stat_dordre(int nbEv,double * vecteur_X, long double * vecteur_D)
{
	int i;
	
	#pragma omp parallel for
	for(i=0;i<=nbEv;i++)
	{
		vecteur_D[i] = (long double)(vecteur_X[i+1] - vecteur_X[i]);
	}
}

/** fonction : Beta
 ** appel la fonction pbeta de R
 ** entre :
 ** x  valeur entre 0 et 1
 ** pin  valeur de p
 ** qin  valeur de q
 ** sortie :
 ** valeur de pbeta
 **/
long double Beta(double x, double pin , double qin )
{
	return (long double)pbeta(x,pin,qin,1,0);
}
