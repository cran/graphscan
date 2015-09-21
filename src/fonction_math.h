// ** ================================================================================================================
// ** R-Package: graphscan 1.1
// ** Fichier : src/fonction_math.h
// ** Description : détection 1D des clusters avec l'indice de Cucala et celui de Kulldorff.
// **               Déclarations des fonctions de normalisation, de calcul d'ordre et beta.
// ** License : GPL-2 | GPL-3
// ** Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
// ** ================================================================================================================

#ifndef FONCTION_MATH_H
#define FONCTION_MATH_H

#include "types.h"


void normalisation_et_distance_entre_stat_dordre(int normalisation_debut,int normalisation_fin,int nbEv,double * vecteur_X, long double * vecteur_D, long double * tab_p);

void distance_entre_stat_dordre(int nbEv,double * vecteur_X, long double * vecteur_D);

long double Beta(double x, double pin , double qin );

#endif
