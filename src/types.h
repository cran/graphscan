// ** ================================================================================================================
// ** R-Package: graphscan
// ** Fichier : src/types.h
// ** Description : contient les structures utilisées par les fonctions.
// ** License : GPL-2 | GPL-3
// ** Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
// ** ================================================================================================================

#ifndef TYPES_H
#define TYPES_H

#include <stdint.h>
#include <float.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

// structure de sortie du programme
typedef struct resultat {
	long double debut; // debut de l'agregat
	long double fin;   // fin de l'agregat
	long double  indice_concentration; // indice de concentration (positif ou negatif)
	double  p_valeur; // p_valeur de l'agregat (positif ou negatif)
	uint8_t  positif; // 0 = negatif , 1 = positif
	int id_agregat; // numero d'ordre
	struct resultat * suiv;
}resultat_t;

typedef struct ListeRepere { 
  resultat_t *debut; 
  resultat_t *fin; 
  int taille; 
}list_t; 

#define MAX_DOUBLE  DBL_MAX


typedef struct agregat_potentiel_indice_cucala {
	int indice_debut;
	int indice_fin;
	long double indice_cucala;
}agregat_potentiel_indice_cucala_t;

#endif
