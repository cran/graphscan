// ** ================================================================================================================
// ** R-Package: graphscan
// ** Fichier : src/detection_dagregat.h
// ** Description : détection 1D des clusters avec l'indice de Cucala et celui de Kulldorff. Déclarations des fonctions.
// ** License : GPL-2 | GPL-3
// ** Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
// ** ================================================================================================================

#ifndef DETECTION_AGREGAT_H
#define DETECTION_AGREGAT_H

#include "types.h"

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

SEXP detection_multiple_dagregat(SEXP nbEv_r,
	SEXP normalisation_debut_r,SEXP normalisation_fin_r,SEXP vecteur_X_r,
	SEXP alpha_r, SEXP theta_r, SEXP nbSim_r,SEXP choix_detection_r, SEXP choix_type_agregat_r);

void decalage_tableau(long double * vecteur_D, long double * vecteur_P,int indice_debut, int indice_fin, int *nbEv);
#endif
