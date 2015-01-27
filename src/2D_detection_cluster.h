// ** ================================================================================================================
// ** R-Package: graphscan
// ** Fichier : src/2D_detection_cluster.c
// ** Description : détection 2D et 3D des clusters avec l'indice de Cucala et celui de Kulldorff. Déclratations
// **               des fonctions detection_cluster() et calcul_concentration()
// ** License : GPL-2 | GPL-3
// ** Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
// ** ================================================================================================================


#ifndef DETECTION_CLUSTER_2D_H
#define DETECTION_CLUSTER_2D_H

#include "2D_types.h"
#include "2D_fonction_math.h"

SEXP detection_cluster(SEXP nb_point_r,
	SEXP dimension_r,SEXP nb_simulation_r,
	SEXP id_r, SEXP coordonnees_r, SEXP controle_r,SEXP cas_r, SEXP memory_size_r);

int calcul_concentration(Point_element *point_cas, Point_element *point_zero_cas,
	List_cluster_cucala *cluster_cucala, List_cluster_kulldorff *cluster_kulldorff,
	ID_TYPE nb_cas, ID_TYPE nb_controle, ID_TYPE nb_point_cas, ID_TYPE nb_point_zero_cas, ID_TYPE nb_dimension,
	ID_TYPE simulation);
#endif
