// ** ================================================================================================================
// ** R-Package: graphscan
// ** Fichier : src/2D_types.h
// ** Description : détection 2D et 3D des clusters avec l'indice de Cucala et celui de Kulldorff. Déclarations des
// **               differentes structures utilisées par les fonctions.
// ** License : GPL-2 | GPL-3
// ** Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
// ** ================================================================================================================

#ifndef TYPES_2D_H
#define TYPES_2D_H

#include <stdint.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>



// type de l'identifiant utilisé pour les points, les mesures,...
#define ID_TYPE int
// type de la valeur de la distance (pas d'unsigned, comparaison utilisée)
#define DISTANCE_TYPE long double //long double 
// maximum du type ci-dessus
#define MAX_DISTANCE LDBL_MAX //LDBL_MAX
 // minimum du type ci-dessus
#define MIN_DISTANCE -LDBL_MAX //LDBL_MIN

// Petite macro pour vérifier l'allocation mémoire et arrêter le programme en cas d'erreur d'allocation mémoire
#define VERIFY_ALLOC(MVar) do{ 			\
 	if((MVar)==NULL){ 		\
 		error("\nERROR: unable to allocate the memory for the variable "#MVar" - terminating\n");	\
	 }		\
 }while(0)




// typedef des structures
typedef struct point Point;
typedef struct point_element Point_element;
typedef struct cluster Cluster;
typedef struct distance Distance;
typedef struct list_cluster_cucala List_cluster_cucala;
typedef struct list_cluster_kulldorff List_cluster_kulldorff;
typedef struct rbtree_value_elem_t *rbtree_value_elem;
typedef struct rbtree_node_t *rbtree_node;
typedef struct rbtree_t *rbtree;

enum rbtree_node_color { RED, BLACK };
typedef enum rbtree_node_color rbtree_color;



/*************************************************************/
/******************** Définition du point ********************/
/*************************************************************/

// 	88""Yb  dP"Yb  88 88b 88 888888 
// 	88__dP dP   Yb 88 88Yb88   88   
// 	88"""  Yb   dP 88 88 Y88   88   
// 	88      YbodP  88 88  Y8   88   

//structure representant un point
struct point{
	DISTANCE_TYPE id;                // identifiant du point (DISTANCE_TYPE pour pouvoir l'utiliser comme id dans l'arbre rb)
	double * coordonnees;  // vecteur de coordonnees de dimension n
	ID_TYPE nb_total;          // nombre de cas + nombre de controles
	Point *suiv_cluster_cucala;		// liste chainee pour créer un cluster validé par cucala 
	Point *suiv_cluster_kulldorff;	// liste chainee pour créer un cluster validé par kulldorf.
};



/*************************************************************/
/******** Structures pour créer des listes de points *********/
/*************************************************************/

// 	88""Yb  dP"Yb  88 88b 88 888888            888888 88     888888 8b    d8 888888 88b 88 888888 
// 	88__dP dP   Yb 88 88Yb88   88              88__   88     88__   88b  d88 88__   88Yb88   88   
// 	88"""  Yb   dP 88 88 Y88   88              88""   88  .o 88""   88YbdP88 88""   88 Y88   88   
// 	88      YbodP  88 88  Y8   88   oooooooooo 888888 88ood8 888888 88 YY 88 888888 88  Y8   88   

// structure representant un élément d'un arbre kd, ainsi qu'un élément des cluster
// les éléments de point pouvant être différents selon les simulations, c'est cette structure que l'on
// va utiliser pour stocker les éléments typique uniques aux simulations. Cela permet aussi de paralléliser facilement si besoin.
// En pratique, ces éléments seront déclaré sous forme de tableau cas et controle joints (voir algo)
struct point_element{
	Point *point; 							// pointeur sur le point d'origine
	Point_element *suivant_cluster;			// pointeur sur l'élément suivant dans le cluster
	ID_TYPE nb_cas;        					// nombre de cas
	ID_TYPE nb_controle;    				// nombre de controles
	Cluster *cluster;						// pointeur vers le cluster auquel appartient le point (utilisé que pour les point cas, qui ne peuvent avoir qu'un cluster)
	Cluster *cluster_sous_arbre;			// pointeur vers le cluster auquel appartient le point et tout le sous-arbre du point (utilisé que pour les point cas, qui ne peuvent avoir qu'un cluster)
	rbtree controle_cluster;				// Arbre de pointeur vers les clusters auxquels appartient le point (utilisé que si le point est un controle)
	rbtree controle_cluster_sous_arbre; 	// Arbre de pointeur vers le ou les clusters auxquels appartient le point et l'ensemble de l'arbre sous le point (utilisé que si le point est un controle)
											// (NULL si des points ont plusieurs clusters ou aucun)
											// Obligé d'utiliser un arbre rb à cause des points de controle qui peuvent appartenir à plusieurs clusters
	Point_element *kd_suivant_gauche;	// pointeur sur le Point_element qui correspond à la branche suivante gauche de l'arbre KD
	Point_element *kd_suivant_droite;	// pointeur sur le Point_element qui correspond à la branche suivante droite de l'arbre KD
	Point_element *kd_parent;			// pointeur sur le Point_element qui correspond au point précédant l'arbre KD
};



/*************************************************************/
/************ Structures pour créer des clusters *************/
/*************************************************************/

// 	 dP""b8 88     88   88 .dP"Y8 888888 888888 88""Yb 
// 	dP   `" 88     88   88 `Ybo."   88   88__   88__dP 
// 	Yb      88  .o Y8   8P o.`Y8b   88   88""   88"Yb  
// 	 YboodP 88ood8 `YbodP' 8bodP'   88   888888 88  Yb 

//structure representant la tete d'un cluster, et permettant de faire une liste doublement chainee de clusters
struct cluster{
	DISTANCE_TYPE id_controle;			// identifiant qui va servir à savoir dans l'arbre kd controle à quel cluster appartient le point
										// on utilise un identifiant différent de simplement l'adresse du cluster pour éviter le parcours de tous les points de controle lors de la fusion des arbres controle
										// (on prend le plus petit qu'on fusionne dans le plus grand) (voir fusion)
										// Le problème étant que les clusters étant libérés en mémoire au cours du programme, il est possible que 2 clusters aient la même adresse mémoire (pas en même temps)
	Point_element *premier_point;		// liste chainée des points qui forment le cluster - premier point
	Point_element *dernier_point;		// liste chainée des points qui forment le cluster - dernier point
	ID_TYPE nb_cas;   					// nombre de cas dans le cluster
	ID_TYPE nb_controle;  				// nombre de controles dans le cluster (dans les points cas)
	ID_TYPE nb_controle_pur;			// nombre de controles dans les points de controle de l'entourage du cluster
	ID_TYPE taille;  					// nombre de points dans le cluster
	rbtree arbre_distance_cas_controle;		// arbre des distance de chaque cas avec le plus proche de ses controle
	rbtree arbre_points_controle; 			// arbre des points de controle qui ont déjà été ajouté au cluster	
};



/****************************************************************************/
/** Structures pour gérer les distances dans un arbre binaire de recherche **/
/****************************************************************************/

// 	8888b.  88 .dP"Y8 888888    db    88b 88  dP""b8 888888 
// 	 8I  Yb 88 `Ybo."   88     dPYb   88Yb88 dP   `" 88__   
// 	 8I  dY 88 o.`Y8b   88    dP__Yb  88 Y88 Yb      88""   
// 	8888Y"  88 8bodP'   88   dP""""Yb 88  Y8  YboodP 888888 

// Structure representant une distance entre 2 points
// Dans le cas d'une distance cas-controle le cas est le point1, le controle le point2
struct distance{
	Point_element *point1;           // index du premier point dans le tableau arbre kd correspondant (arbre_kd_cas)
	Point_element *point2;           // index du deuxieme point dans le tableau arbre kd correspondant (arbre_kd_cas ou arbre_kd_control)
	DISTANCE_TYPE valeur;			// valeur de la distance (au carré)
};

/***************************************************************/
/** Structures pour gérer les clusters de Cucala et Kulldorff **/
/***************************************************************/

// 	88     88 .dP"Y8 888888             dP""b8 88     88   88 .dP"Y8 888888 888888 88""Yb             dP""b8 88   88  dP""b8    db    88        db    
// 	88     88 `Ybo."   88              dP   `" 88     88   88 `Ybo."   88   88__   88__dP            dP   `" 88   88 dP   `"   dPYb   88       dPYb   
// 	88  .o 88 o.`Y8b   88              Yb      88  .o Y8   8P o.`Y8b   88   88""   88"Yb             Yb      Y8   8P Yb       dP__Yb  88  .o  dP__Yb  
// 	88ood8 88 8bodP'   88   oooooooooo  YboodP 88ood8 `YbodP' 8bodP'   88   888888 88  Yb oooooooooo  YboodP `YbodP'  YboodP dP""""Yb 88ood8 dP""""Yb 

// structure representant la tete de la liste contenant les points et
// les infos du cluster de cucala
struct list_cluster_cucala{
	Point * debut;
	Point * fin;
	Point_element *tmp_elem_debut;	// Variable temporaire utilisée pour éviter des parcours inutiles de listes de points
	Point_element *tmp_elem_fin;	// Variable temporaire utilisée pour éviter des parcours inutiles de listes de points
	DISTANCE_TYPE dist_max;
	DISTANCE_TYPE concentration_cucala;
	ID_TYPE nb_cas_cucala;
	ID_TYPE nb_controle_cucala;
	ID_TYPE taille;
};

// 	88     88 .dP"Y8 888888             dP""b8 88     88   88 .dP"Y8 888888 888888 88""Yb            88  dP 88   88 88     88     8888b.   dP"Yb  88""Yb 888888 888888 
// 	88     88 `Ybo."   88              dP   `" 88     88   88 `Ybo."   88   88__   88__dP            88odP  88   88 88     88      8I  Yb dP   Yb 88__dP 88__   88__   
// 	88  .o 88 o.`Y8b   88              Yb      88  .o Y8   8P o.`Y8b   88   88""   88"Yb             88"Yb  Y8   8P 88  .o 88  .o  8I  dY Yb   dP 88"Yb  88""   88""   
// 	88ood8 88 8bodP'   88   oooooooooo  YboodP 88ood8 `YbodP' 8bodP'   88   888888 88  Yb oooooooooo 88  Yb `YbodP' 88ood8 88ood8 8888Y"   YbodP  88  Yb 88     88     

// structure representant la tete de la liste contenant les points et
// les infos du cluster de kulldorff
struct list_cluster_kulldorff{
	Point * debut;
	Point * fin;
	Point_element *tmp_elem_debut;	// Variable temporaire utilisée pour éviter des parcours inutiles de listes de points
	Point_element *tmp_elem_fin;	// Variable temporaire utilisée pour éviter des parcours inutiles de listes de points
	DISTANCE_TYPE concentration_kulldorff;
	DISTANCE_TYPE dist_max;
	ID_TYPE nb_cas_kulldorff;
	ID_TYPE nb_controle_kulldorff;
	ID_TYPE taille;
};


/*************************************************/
/** Structures pour gérer les arbres rouge-noir **/
/*************************************************/

// 	88""Yb 88""Yb 888888 88""Yb 888888 888888            Yb    dP    db    88     88   88 888888            888888 88     888888 8b    d8            888888 
// 	88__dP 88__dP   88   88__dP 88__   88__               Yb  dP    dPYb   88     88   88 88__              88__   88     88__   88b  d88              88   
// 	88"Yb  88""Yb   88   88"Yb  88""   88""                YbdP    dP__Yb  88  .o Y8   8P 88""              88""   88  .o 88""   88YbdP88              88   
// 	88  Yb 88oodP   88   88  Yb 888888 888888 oooooooooo    YP    dP""""Yb 88ood8 `YbodP' 888888 oooooooooo 888888 88ood8 888888 88 YY 88 oooooooooo   88   

struct rbtree_value_elem_t{
    void *value;                    
    rbtree_value_elem next;
};

// 	88""Yb 88""Yb 888888 88""Yb 888888 888888            88b 88  dP"Yb  8888b.  888888            888888 
// 	88__dP 88__dP   88   88__dP 88__   88__              88Yb88 dP   Yb  8I  Yb 88__                88   
// 	88"Yb  88""Yb   88   88"Yb  88""   88""              88 Y88 Yb   dP  8I  dY 88""                88   
// 	88  Yb 88oodP   88   88  Yb 888888 888888 oooooooooo 88  Y8  YbodP  8888Y"  888888 oooooooooo   88   

struct rbtree_node_t{
    DISTANCE_TYPE key;
    rbtree_value_elem value_first;
    rbtree_node left;
    rbtree_node right;
    rbtree_node parent;
    rbtree_color color;
};

// 	88""Yb 88""Yb 888888 88""Yb 888888 888888            888888 
// 	88__dP 88__dP   88   88__dP 88__   88__                88   
// 	88"Yb  88""Yb   88   88"Yb  88""   88""                88   
// 	88  Yb 88oodP   88   88  Yb 888888 888888 oooooooooo   88   

struct rbtree_t{
    rbtree_node root;
    rbtree_node smallest;
    ID_TYPE size;    // taille de l'arbre
};

#endif
