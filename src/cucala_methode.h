// ** ================================================================================================================
// ** R-Package: graphscan
// ** Fichier : src/detection_dagregat.h
// ** Description : détection 1D des clusters avec l'indice de Cucala et celui de Kulldorff. Déclarations des
// **               differentes structures utilisées par les fonctions.
// ** License : GPL-2 | GPL-3
// ** Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
// ** ================================================================================================================

#ifndef CUCALA_METHODE_H
#define CUCALA_METHODE_H

#include "types.h"


agregat_potentiel_indice_cucala_t calcul_agregat_positif_et_indice_cucala(
int nbEv, long double * vecteur_D);

agregat_potentiel_indice_cucala_t calcul_agregat_negatif_et_indice_cucala(
int nbEv, long double * vecteur_D);

double calcul_p_valeur_positif(int nbEv, int nbsim, long double indice_cucala);

double calcul_p_valeur_negatif(int nbEv, int nbsim, long double indice_cucala);

void calcul_p_valeur_negatif_positif(int nbEv, int nbsim, long double indice_cucala_pos,long double indice_cucala_neg,double * p_val_min, double * p_val_max);

int compare_doubles(const void *a, const void *b);


#endif
