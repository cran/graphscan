// ** ================================================================================================================
// ** R-Package: graphscan 1.1
// ** Fichier : src/2D_fonction_math.c
// ** Description : détection 2D et 3D des clusters avec l'indice de Cucala et celui de Kulldorff. Manipulations
// **               des arbres Kd. Déclarations de fonctions.
// ** License : GPL-2 | GPL-3
// ** Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
// ** ================================================================================================================


#ifndef FONCTION_MATH_2D_H
#define FONCTION_MATH_2D_H

#include "2D_types.h"
#include "rbtree.h"
                                                                                             

/*****************************************************/
/********** Manipulation sur les distances ***********/
/*****************************************************/

// Calcul de la distance euclidienne au carré (pas de sqrt pour économiser du temps)
DISTANCE_TYPE distance_euclidienne2(Point_element *point_i, Point_element *point_j, ID_TYPE nb_dimension);

// Fonction qui calcule la distance entre 2 points et ajoute cette distance à l'arbre rouge-noir passé en paramètre
void ajouter_distance(rbtree arbre_distance, Point_element *point1, Point_element *point2, ID_TYPE nb_dimension);

/*****************************************************/
/************* Manipulation des clusters *************/
/*****************************************************/

// Fonction qui crée un nouveau cluster vide
Cluster *creer_cluster(Point_element *point1, Point_element *point2, Point_element *arbre_controle, ID_TYPE nb_dimension);

// Fonction qui ajoute un point à un cluster
void cluster_add_point(Cluster *cluster, Point_element *point, Point_element *arbre_controle, ID_TYPE nb_dimension);

// Fonction qui ajoute un point de controle à un cluster
void cluster_add_point_controle(Cluster *cluster, Point_element *point_controle);

// Fonction qui déplace toutes les distances de l'arbre 2 dans l'arbre 1, puis libère la mémoire de l'arbre 2
void fusion_rbtree_distance(rbtree tree1, rbtree tree2);

// Fonction  qui fusionne l'arbre de controle de l'ancien cluster dans l'arbre de controle du nouveau cluster
// libère l'arbre de controle de l'ancien cluster
void fusion_rbtree_controle(Cluster *new_cluster, Cluster *old_cluster);

// Fonction qui fusionne le plus petit cluster dans le plus grand
Cluster *fusion_cluster(Cluster *cluster1, Cluster *cluster2);

/*****************************************************/
/************ Manipulation des arbres kd *************/
/*****************************************************/

// Fonction récursive qui crée un arbre kd avec tab_point le pointeur vers le premier élément du tableau de point de l'arbre, parent le point pivot précédent
// nb_point le nombre de point à partir du pointeur à considérer, nb_dimension le nombre de dimension de l'arbre et dimension_courante la dimension du niveau courant de l'arbre
// En cas de nombre pair de point, on mettra plus de points dans la partie gauche de l'arbre.
// Du coup pour compenser, quand plusieurs points auront la même valeur que le point pivot, on les mettra dans la partie droite de l'arbre
Point_element *creer_arbre_kd(Point_element *tab_point, Point_element *parent, ID_TYPE nb_point, ID_TYPE nb_dimension, ID_TYPE dimension_courante);

// Fonction récursive qui parcours l'arbre (ou sous-arbre) kd pour renvoyer le plus proche voisin du point passé en paramètre
// arbre_kd est le sous-arbre sur lequel on travaille, point le point dont on cherche le plus proche voisin, is_controle est un booléen pour savoir
// si l'arbre est un arbre de points de controle ou non, nb_dimension est le nombre de dimension
Distance *get_plus_proche_voisin(Point_element *arbre_kd, Point_element *point, int is_controle, ID_TYPE nb_dimension);

// Fonction récursive qui parcours l'arbre (ou sous-arbre) kd pour renvoyer le plus proche voisin du point passé en paramètre
// point_kd est le sous-arbre sur lequel on travaille, point le point dont on cherche le plus proche voisin, meilleure_distance est la meilleure distance actuelle,
// nb_dimension est le nombre de dimension, dimension_courante la dimension courante
void recur_get_plus_proche_voisin_cas(Point_element *point_kd, Point_element *point, Distance *meilleure_distance, ID_TYPE nb_dimension, ID_TYPE dimension_courante);

// Fonction récursive qui parcours l'arbre (ou sous-arbre) kd pour renvoyer le plus proche voisin du point passé en paramètre
// point_kd est le sous-arbre sur lequel on travaille, point le point dont on cherche le plus proche voisin, meilleure_distance est la meilleure distance actuelle,
// nb_dimension est le nombre de dimension, dimension_courante la dimension courante
void recur_get_plus_proche_voisin_controle(Point_element *point_kd, Point_element *point, Distance *meilleure_distance, ID_TYPE nb_dimension, ID_TYPE dimension_courante);

// Fonction récursive qui déplace le noeud node et tous ces fils dans l'arbre tree
void recur_fusion_rbtree_distance(rbtree tree, rbtree_node node);

// Fonction récursive qui fusionne le node (et ses fils) de l'ancien cluster dans l'arbre des controles du nouveau cluster
// node est le noeud courant (de old_cluster->arbre_points_controle) sur lequel on est en train de travailler
void recur_fusion_rbtree_controle(rbtree_node node, Cluster *new_cluster, Cluster *old_cluster);




#endif
