// ** ================================================================================================================
// ** R-Package: graphscan 1.1
// ** Fichier : src/2D_fonction_math.c
// ** Description : détection 2D et 3D des clusters avec l'indice de Cucala et celui de Kulldorff. Manipulations
// **               des arbres Kd.
// ** License : GPL-2 | GPL-3
// ** Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
// ** ================================================================================================================

#include "2D_types.h"
#include <math.h>
#include "qsort_elem.h"
#include "2D_fonction_math.h"


/*****************************************************/
/************* Manipulation de l'arbre KD ************/
/*****************************************************/

// 	 dP""b8 88""Yb 888888 888888 88""Yb               db    88""Yb 88""Yb 88""Yb 888888            88  dP 8888b.  
// 	dP   `" 88__dP 88__   88__   88__dP              dPYb   88__dP 88__dP 88__dP 88__              88odP   8I  Yb 
// 	Yb      88"Yb  88""   88""   88"Yb              dP__Yb  88"Yb  88""Yb 88"Yb  88""              88"Yb   8I  dY 
// 	 YboodP 88  Yb 888888 888888 88  Yb oooooooooo dP""""Yb 88  Yb 88oodP 88  Yb 888888 oooooooooo 88  Yb 8888Y"  

// Fonction récursive qui crée un arbre kd avec tab_point le pointeur vers le premier élément du tableau de point de l'arbre, parent le point pivot précédent
// nb_point le nombre de point à partir du pointeur à considérer, nb_dimension le nombre de dimension de l'arbre et dimension_courante la dimension du niveau courant de l'arbre
// En cas de nombre pair de point, on mettra plus de points dans la partie gauche de l'arbre.
// Du coup pour compenser, quand plusieurs points auront la même valeur que le point pivot, on les mettra dans la partie droite de l'arbre
Point_element *creer_arbre_kd(Point_element *tab_point, Point_element *parent, ID_TYPE nb_point, ID_TYPE nb_dimension, ID_TYPE dimension_courante){
	Point_element *pivot;
	ID_TYPE nb_point_gauche, nb_point_droite;

	// On regarde si il y a besoin de classer
	if(nb_point <= 0)
		return NULL;
	if(nb_point == 1)
	{
		tab_point[0].kd_parent = parent;
		return &tab_point[0];
	}
	// Sinon on trie le tableau
	quicksort_elem(tab_point, nb_point, sizeof(Point_element), dimension_courante);	// Amélioration perso basée sur qsort pour la rendre thread-safe
	// On récupère le point pivot (médian) en prenant le milieu du tableau
	pivot = &tab_point[nb_point/2];	// En cas de nombre pair, c'est le point de droite qui est pris
	nb_point_gauche = (nb_point/2);
	// Si des points à gauche ont la même valeur que le pivot, on décale le pivot vers la gauche pour les mettre sur la droite
	while(pivot != tab_point && (pivot-1)->point->coordonnees[dimension_courante] == pivot->point->coordonnees[dimension_courante] )
	{
		pivot--;	// arithmétique de pointeur
		nb_point_gauche--;
	}
	// Pivot est maintenant correctement positionné
	nb_point_droite = nb_point - nb_point_gauche - 1; 	// -1 le pivot
	pivot->kd_parent = parent;
	// Appels récursifs de l'arbre pour créer les parties gauches et droites de l'arbre
	pivot->kd_suivant_droite = creer_arbre_kd(pivot+1, pivot, nb_point_droite, nb_dimension, (dimension_courante+1)%nb_dimension);
	pivot->kd_suivant_gauche = creer_arbre_kd(tab_point, pivot, nb_point_gauche, nb_dimension, (dimension_courante+1)%nb_dimension);

	// On retourne le pivot
	return pivot;
}


// 	88""Yb 888888  dP""b8 88   88 88""Yb             dP""b8 888888 888888            88""Yb 88     88   88 .dP"Y8            88""Yb 88""Yb  dP"Yb   dP""b8 88  88 888888            Yb    dP  dP"Yb  88 .dP"Y8 88 88b 88             dP""b8    db    .dP"Y8 
// 	88__dP 88__   dP   `" 88   88 88__dP            dP   `" 88__     88              88__dP 88     88   88 `Ybo."            88__dP 88__dP dP   Yb dP   `" 88  88 88__               Yb  dP  dP   Yb 88 `Ybo." 88 88Yb88            dP   `"   dPYb   `Ybo." 
// 	88"Yb  88""   Yb      Y8   8P 88"Yb             Yb  "88 88""     88              88"""  88  .o Y8   8P o.`Y8b            88"""  88"Yb  Yb   dP Yb      888888 88""                YbdP   Yb   dP 88 o.`Y8b 88 88 Y88            Yb       dP__Yb  o.`Y8b 
// 	88  Yb 888888  YboodP `YbodP' 88  Yb oooooooooo  YboodP 888888   88   oooooooooo 88     88ood8 `YbodP' 8bodP' oooooooooo 88     88  Yb  YbodP   YboodP 88  88 888888 oooooooooo    YP     YbodP  88 8bodP' 88 88  Y8 oooooooooo  YboodP dP""""Yb 8bodP' 

// Fonction récursive qui parcours l'arbre (ou sous-arbre) kd pour renvoyer le plus proche voisin du point passé en paramètre
// point_kd est le sous-arbre sur lequel on travaille, point le point dont on cherche le plus proche voisin, meilleure_distance est la meilleure distance actuelle,
// nb_dimension est le nombre de dimension, dimension_courante la dimension courante
void recur_get_plus_proche_voisin_cas(Point_element *point_kd, Point_element *point, Distance *meilleure_distance, ID_TYPE nb_dimension, ID_TYPE dimension_courante){
	DISTANCE_TYPE tmp_dist;
	// On ne parcours pas si le sous-arbre appartient entièrement au même cluster que le point
	if( point->cluster!=NULL && point_kd->cluster_sous_arbre==point->cluster)
		return;	// Le sous-arbre de cas est déjà dans le cluster
	// Si on est pas sur le point lui-même, et que le point n'appartient pas déjà au cluster
	if( point_kd != point && (point->cluster==NULL || point->cluster!=point_kd->cluster) )
	{
		// On calcule la distance au carré entre le point cas et le point controle
		tmp_dist = distance_euclidienne2(point, point_kd, nb_dimension);
		// Si elle est meilleure, on met la meilleure distance à jour
		if(tmp_dist < meilleure_distance->valeur)
		{
			meilleure_distance->valeur = tmp_dist;
			meilleure_distance->point2 = point_kd;
		}
	}
	// On regarde la distance entre le point et le plan sécant (pas encore au carré, pour voir si le point se trouve à gauche ou à droite)
	tmp_dist = point->point->coordonnees[dimension_courante] - point_kd->point->coordonnees[dimension_courante];
	// On privilégie le sous-arbre dans lequel se trouve le point
	if(tmp_dist < 0)
	{
		// Le point est à gauche
		// On transforme la différence en carré car on veut comparer avec une distance au carré
		tmp_dist = tmp_dist*tmp_dist;
		// parcours gauche en priorité, il est possible qu'il y est un meilleur point dans le sous-arbre gauche
		if(point_kd->kd_suivant_gauche != NULL)
			recur_get_plus_proche_voisin_cas(point_kd->kd_suivant_gauche, point, meilleure_distance, nb_dimension, (dimension_courante+1)%nb_dimension);
		// puis parcours droite, si c'est possible de trouver un plus proche
		if(point_kd->kd_suivant_droite != NULL && tmp_dist <= meilleure_distance->valeur)
			recur_get_plus_proche_voisin_cas(point_kd->kd_suivant_droite, point, meilleure_distance, nb_dimension, (dimension_courante+1)%nb_dimension);
	}else{
		// Le point est à droite
		// On transforme la différence en carré car on veut comparer avec une distance au carré
		tmp_dist = tmp_dist*tmp_dist;
		// parcours droite en priorité, il est possible qu'il y est un meilleur point dans le sous-arbre droite
		if(point_kd->kd_suivant_droite != NULL)
			recur_get_plus_proche_voisin_cas(point_kd->kd_suivant_droite, point, meilleure_distance, nb_dimension, (dimension_courante+1)%nb_dimension);
		// puis parcours gauche, si c'est possible de trouver un plus proche
		if(point_kd->kd_suivant_gauche != NULL && tmp_dist <= meilleure_distance->valeur)
			recur_get_plus_proche_voisin_cas(point_kd->kd_suivant_gauche, point, meilleure_distance, nb_dimension, (dimension_courante+1)%nb_dimension);
	}
}


// 	88""Yb 888888  dP""b8 88   88 88""Yb             dP""b8 888888 888888            88""Yb 88     88   88 .dP"Y8            88""Yb 88""Yb  dP"Yb   dP""b8 88  88 888888            Yb    dP  dP"Yb  88 .dP"Y8 88 88b 88             dP""b8  dP"Yb  88b 88 888888 88""Yb  dP"Yb  88     888888 
// 	88__dP 88__   dP   `" 88   88 88__dP            dP   `" 88__     88              88__dP 88     88   88 `Ybo."            88__dP 88__dP dP   Yb dP   `" 88  88 88__               Yb  dP  dP   Yb 88 `Ybo." 88 88Yb88            dP   `" dP   Yb 88Yb88   88   88__dP dP   Yb 88     88__   
// 	88"Yb  88""   Yb      Y8   8P 88"Yb             Yb  "88 88""     88              88"""  88  .o Y8   8P o.`Y8b            88"""  88"Yb  Yb   dP Yb      888888 88""                YbdP   Yb   dP 88 o.`Y8b 88 88 Y88            Yb      Yb   dP 88 Y88   88   88"Yb  Yb   dP 88  .o 88""   
// 	88  Yb 888888  YboodP `YbodP' 88  Yb oooooooooo  YboodP 888888   88   oooooooooo 88     88ood8 `YbodP' 8bodP' oooooooooo 88     88  Yb  YbodP   YboodP 88  88 888888 oooooooooo    YP     YbodP  88 8bodP' 88 88  Y8 oooooooooo  YboodP  YbodP  88  Y8   88   88  Yb  YbodP  88ood8 888888 

// Fonction récursive qui parcours l'arbre (ou sous-arbre) kd pour renvoyer le plus proche voisin du point passé en paramètre
// point_kd est le sous-arbre sur lequel on travaille, point le point dont on cherche le plus proche voisin, meilleure_distance est la meilleure distance actuelle,
// nb_dimension est le nombre de dimension, dimension_courante la dimension courante
void recur_get_plus_proche_voisin_controle(Point_element *point_kd, Point_element *point, Distance *meilleure_distance, ID_TYPE nb_dimension, ID_TYPE dimension_courante){
	DISTANCE_TYPE tmp_dist;
	// On ne parcours pas si le sous-arbre appartient entièrement au même cluster que le point
	// Le point appartient forcément à un cluster, vu que l'on cherche son plus proche voisin controle
	if( rbtree_find(point_kd->controle_cluster_sous_arbre, point->cluster->id_controle) != NULL )
		return;
	// Si le point n'appartient pas déjà au cluster
	if( rbtree_find(point_kd->controle_cluster, point->cluster->id_controle) == NULL )
	{
		// On calcule la distance au carré entre le point cas et le point controle
		tmp_dist = distance_euclidienne2(point, point_kd, nb_dimension);
		// Si elle est meilleure, on met la meilleure distance à jour
		if(tmp_dist < meilleure_distance->valeur)
		{
			meilleure_distance->valeur = tmp_dist;
			meilleure_distance->point2 = point_kd;
		}
	}
	// On regarde la distance entre le point et le plan sécant (pas encore au carré, pour voir si le point se trouve à gauche ou à droite)
	tmp_dist = point->point->coordonnees[dimension_courante] - point_kd->point->coordonnees[dimension_courante];
	// On privilégie le sous-arbre dans lequel se trouve le point
	if(tmp_dist < 0)
	{
		// Le point est à gauche
		// On transforme la différence en carré car on veut comparer avec une distance au carré
		tmp_dist = tmp_dist*tmp_dist;
		// parcours gauche en priorité, il est possible qu'il y est un meilleur point dans le sous-arbre gauche
		if(point_kd->kd_suivant_gauche != NULL)
			recur_get_plus_proche_voisin_controle(point_kd->kd_suivant_gauche, point, meilleure_distance, nb_dimension, (dimension_courante+1)%nb_dimension);
		// puis parcours droite, si c'est possible de trouver un plus proche
		if(point_kd->kd_suivant_droite != NULL && tmp_dist <= meilleure_distance->valeur)
			recur_get_plus_proche_voisin_controle(point_kd->kd_suivant_droite, point, meilleure_distance, nb_dimension, (dimension_courante+1)%nb_dimension);
	}else{
		// Le point est à droite
		// On transforme la différence en carré car on veut comparer avec une distance au carré
		tmp_dist = tmp_dist*tmp_dist;
		// parcours droite en priorité, il est possible qu'il y est un meilleur point dans le sous-arbre droite
		if(point_kd->kd_suivant_droite != NULL)
			recur_get_plus_proche_voisin_controle(point_kd->kd_suivant_droite, point, meilleure_distance, nb_dimension, (dimension_courante+1)%nb_dimension);
		// puis parcours gauche, si c'est possible de trouver un plus proche
		if(point_kd->kd_suivant_gauche != NULL && tmp_dist <= meilleure_distance->valeur)
			recur_get_plus_proche_voisin_controle(point_kd->kd_suivant_gauche, point, meilleure_distance, nb_dimension, (dimension_courante+1)%nb_dimension);
	}
}




// 	88""Yb 888888  dP""b8 88   88 88""Yb            888888 88   88 .dP"Y8 88  dP"Yb  88b 88            88""Yb 88""Yb 888888 88""Yb 888888 888888            8888b.  88 .dP"Y8 888888    db    88b 88  dP""b8 888888 
// 	88__dP 88__   dP   `" 88   88 88__dP            88__   88   88 `Ybo." 88 dP   Yb 88Yb88            88__dP 88__dP   88   88__dP 88__   88__               8I  Yb 88 `Ybo."   88     dPYb   88Yb88 dP   `" 88__   
// 	88"Yb  88""   Yb      Y8   8P 88"Yb             88""   Y8   8P o.`Y8b 88 Yb   dP 88 Y88            88"Yb  88""Yb   88   88"Yb  88""   88""               8I  dY 88 o.`Y8b   88    dP__Yb  88 Y88 Yb      88""   
// 	88  Yb 888888  YboodP `YbodP' 88  Yb oooooooooo 88     `YbodP' 8bodP' 88  YbodP  88  Y8 oooooooooo 88  Yb 88oodP   88   88  Yb 888888 888888 oooooooooo 8888Y"  88 8bodP'   88   dP""""Yb 88  Y8  YboodP 888888 

// Fonction récursif qui déplace le noeud node et tous ces fils dans l'arbre tree
void recur_fusion_rbtree_distance(rbtree tree, rbtree_node node){
	// appel récursif 
	if(node->left != NULL)
        recur_fusion_rbtree_distance(tree, node->left);
    if(node->right != NULL)
        recur_fusion_rbtree_distance(tree, node->right);
    // On déplace le noeud
    rbtree_insert(tree, node, 0);
}




// 	88""Yb 888888  dP""b8 88   88 88""Yb            888888 88   88 .dP"Y8 88  dP"Yb  88b 88            88""Yb 88""Yb 888888 88""Yb 888888 888888             dP""b8  dP"Yb  88b 88 888888 88""Yb  dP"Yb  88     888888 
// 	88__dP 88__   dP   `" 88   88 88__dP            88__   88   88 `Ybo." 88 dP   Yb 88Yb88            88__dP 88__dP   88   88__dP 88__   88__              dP   `" dP   Yb 88Yb88   88   88__dP dP   Yb 88     88__   
// 	88"Yb  88""   Yb      Y8   8P 88"Yb             88""   Y8   8P o.`Y8b 88 Yb   dP 88 Y88            88"Yb  88""Yb   88   88"Yb  88""   88""              Yb      Yb   dP 88 Y88   88   88"Yb  Yb   dP 88  .o 88""   
// 	88  Yb 888888  YboodP `YbodP' 88  Yb oooooooooo 88     `YbodP' 8bodP' 88  YbodP  88  Y8 oooooooooo 88  Yb 88oodP   88   88  Yb 888888 888888 oooooooooo  YboodP  YbodP  88  Y8   88   88  Yb  YbodP  88ood8 888888 

// Fonction récursive qui fusionne le node (et ses fils) de l'ancien cluster dans l'arbre des controles du nouveau cluster
// node est le noeud courant (de old_cluster->arbre_points_controle) sur lequel on est en train de travailler
void recur_fusion_rbtree_controle(rbtree_node node, Cluster *new_cluster, Cluster *old_cluster){
	Point_element *point_controle = RBTREE_GET_AS_CONTROLE(node);
	// Cette fonction appelle récursivement les fils du noeud en question
    if(node->left != NULL)
        recur_fusion_rbtree_controle(node->left, new_cluster, old_cluster);
    if(node->right != NULL)
        recur_fusion_rbtree_controle(node->right, new_cluster, old_cluster);

    // On commence par supprimer l'id de l'ancien cluster de l'arbre des clusters du point de controle
    rbtree_delete(point_controle->controle_cluster, old_cluster->id_controle, 0);
    // On supprime aussi l'id de l'ancien cluster de l'arbre des clusters_sous_arbre du point de controle (le cas échéant, si il n'existe pas la fonction ne fait rien)
    rbtree_delete(point_controle->controle_cluster_sous_arbre, old_cluster->id_controle, 0);

    // Ensuite on ajoute le point à l'arbre des controles du cluster (la fonction n'ajoute rien si il existe déjà)
    cluster_add_point_controle(new_cluster, point_controle);

    // Il ne reste plus qu'à libérer la mémoire de l'ancien noeud
    rbtree_free_value_elem(old_cluster->arbre_points_controle, node, 0);    // on libère la mémoire de la liste de valeur
    free(node);
}



/*****************************************************/
/********** Manipulation sur les distances ***********/
/*****************************************************/

// 	8888b.  88 .dP"Y8 888888    db    88b 88  dP""b8 888888            888888 88   88  dP""b8 88     88 8888b.  88 888888 88b 88 88b 88 888888 oP"Yb. 
// 	 8I  Yb 88 `Ybo."   88     dPYb   88Yb88 dP   `" 88__              88__   88   88 dP   `" 88     88  8I  Yb 88 88__   88Yb88 88Yb88 88__   "' dP' 
// 	 8I  dY 88 o.`Y8b   88    dP__Yb  88 Y88 Yb      88""              88""   Y8   8P Yb      88  .o 88  8I  dY 88 88""   88 Y88 88 Y88 88""     dP'  
// 	8888Y"  88 8bodP'   88   dP""""Yb 88  Y8  YboodP 888888 oooooooooo 888888 `YbodP'  YboodP 88ood8 88 8888Y"  88 888888 88  Y8 88  Y8 888888 .d8888 

// Calcul de la distance euclidienne au carré (pas de sqrt pour économiser du temps)
inline DISTANCE_TYPE distance_euclidienne2(Point_element *point_i, Point_element *point_j, ID_TYPE nb_dimension){
	DISTANCE_TYPE tmp, res = 0;
	for(ID_TYPE k = 0; k < nb_dimension; k++)
	{
		tmp = (DISTANCE_TYPE)(point_i->point->coordonnees[k] - point_j->point->coordonnees[k]);
		res += (tmp*tmp);
	}
	return res;
}

// 	   db     88888  dP"Yb  88   88 888888 888888 88""Yb            8888b.  88 .dP"Y8 888888    db    88b 88  dP""b8 888888 
// 	  dPYb       88 dP   Yb 88   88   88   88__   88__dP             8I  Yb 88 `Ybo."   88     dPYb   88Yb88 dP   `" 88__   
// 	 dP__Yb  o.  88 Yb   dP Y8   8P   88   88""   88"Yb              8I  dY 88 o.`Y8b   88    dP__Yb  88 Y88 Yb      88""   
// 	dP""""Yb "bodP'  YbodP  `YbodP'   88   888888 88  Yb oooooooooo 8888Y"  88 8bodP'   88   dP""""Yb 88  Y8  YboodP 888888 

// Fonction qui calcule la distance entre 2 points et ajoute cette distance à l'arbre rouge-noir passé en paramètre
inline void ajouter_distance(rbtree arbre_distance, Point_element *point1, Point_element *point2, ID_TYPE nb_dimension){
	Distance *new_dist = malloc(sizeof(Distance));
	VERIFY_ALLOC(new_dist);
	new_dist->point1 = point1;
	new_dist->point2 = point2;
	new_dist->valeur = distance_euclidienne2(point1, point2, nb_dimension);
	rbtree_insert_new(arbre_distance, new_dist->valeur, (void*)new_dist, 0);	// 0 car on ajoute des distances, pas des points de controle
}


/*****************************************************/
/************* Manipulation des clusters *************/
/*****************************************************/

// 	 dP""b8 88""Yb 888888 888888 88""Yb             dP""b8 88     88   88 .dP"Y8 888888 888888 88""Yb 
// 	dP   `" 88__dP 88__   88__   88__dP            dP   `" 88     88   88 `Ybo."   88   88__   88__dP 
// 	Yb      88"Yb  88""   88""   88"Yb             Yb      88  .o Y8   8P o.`Y8b   88   88""   88"Yb  
// 	 YboodP 88  Yb 888888 888888 88  Yb oooooooooo  YboodP 88ood8 `YbodP' 8bodP'   88   888888 88  Yb 

// Fonction qui crée un nouveau cluster vide
inline Cluster *creer_cluster(Point_element *point1, Point_element *point2, Point_element *arbre_controle, ID_TYPE nb_dimension){
	Cluster *new_cluster = calloc(1,sizeof(Cluster));
	VERIFY_ALLOC(new_cluster);
	new_cluster->arbre_distance_cas_controle = rbtree_create();
	new_cluster->arbre_points_controle = rbtree_create();
	// On utilise comme identifiant du cluster l'adresse de l'identifiant du premier point du cluster
	new_cluster->id_controle = point1->point->id;
	// On ajoute les deux points au cluster
	cluster_add_point(new_cluster, point1, arbre_controle, nb_dimension);
	cluster_add_point(new_cluster, point2, arbre_controle, nb_dimension);
	return new_cluster;
}

// 	 dP""b8 88     88   88 .dP"Y8 888888 888888 88""Yb               db    8888b.  8888b.             88""Yb  dP"Yb  88 88b 88 888888 
// 	dP   `" 88     88   88 `Ybo."   88   88__   88__dP              dPYb    8I  Yb  8I  Yb            88__dP dP   Yb 88 88Yb88   88   
// 	Yb      88  .o Y8   8P o.`Y8b   88   88""   88"Yb              dP__Yb   8I  dY  8I  dY            88"""  Yb   dP 88 88 Y88   88   
// 	 YboodP 88ood8 `YbodP' 8bodP'   88   888888 88  Yb oooooooooo dP""""Yb 8888Y"  8888Y"  oooooooooo 88      YbodP  88 88  Y8   88   

// Fonction qui ajoute un point à un cluster
inline void cluster_add_point(Cluster *cluster, Point_element *point, Point_element *arbre_controle, ID_TYPE nb_dimension){
	Point_element *tmp_point;
	Distance *tmp_distance;
	// On commence par ajouter le point dans la liste du cluster
	if(cluster->premier_point == NULL)
		cluster->premier_point = point;
	else
		cluster->dernier_point->suivant_cluster = point;
	cluster->dernier_point = point;
	cluster->taille++;
	cluster->nb_cas += point->nb_cas;
	cluster->nb_controle += point->nb_controle;
	point->cluster = cluster;
	// On met à jour les informations du point qui vont être utilisées lors du parcours de l'arbre kd
	// On met à jour le cluster_sous_arbre en remontant l'arbre tant que toute la partie inférieure (droite et gauche) au point courant a lui aussi le sous-arbre du même cluster (ou NULL)
	tmp_point = point;
	while( tmp_point != NULL && tmp_point->cluster == cluster &&
		   (tmp_point->kd_suivant_gauche == NULL || tmp_point->kd_suivant_gauche->cluster_sous_arbre == cluster) &&
		   (tmp_point->kd_suivant_droite == NULL || tmp_point->kd_suivant_droite->cluster_sous_arbre == cluster) )
	{
		tmp_point->cluster_sous_arbre = cluster;
		// On remonte l'arbre
		tmp_point = tmp_point->kd_parent;
	}
	// On regarde le plus proche voisin du point parmi les controles et on l'ajoute à l'arbre
	if(arbre_controle != NULL)
	{
		tmp_distance = get_plus_proche_voisin(arbre_controle, point, 1, nb_dimension);
		rbtree_insert_new(cluster->arbre_distance_cas_controle, tmp_distance->valeur, (void*)tmp_distance, 0);	// 0 car on ajoute des distances
	}
}

// 	 dP""b8 88     88   88 .dP"Y8 888888 888888 88""Yb               db    8888b.  8888b.             88""Yb  dP"Yb  88 88b 88 888888             dP""b8  dP"Yb  88b 88 888888 88""Yb  dP"Yb  88     888888 
// 	dP   `" 88     88   88 `Ybo."   88   88__   88__dP              dPYb    8I  Yb  8I  Yb            88__dP dP   Yb 88 88Yb88   88              dP   `" dP   Yb 88Yb88   88   88__dP dP   Yb 88     88__   
// 	Yb      88  .o Y8   8P o.`Y8b   88   88""   88"Yb              dP__Yb   8I  dY  8I  dY            88"""  Yb   dP 88 88 Y88   88              Yb      Yb   dP 88 Y88   88   88"Yb  Yb   dP 88  .o 88""   
// 	 YboodP 88ood8 `YbodP' 8bodP'   88   888888 88  Yb oooooooooo dP""""Yb 8888Y"  8888Y"  oooooooooo 88      YbodP  88 88  Y8   88   oooooooooo  YboodP  YbodP  88  Y8   88   88  Yb  YbodP  88ood8 888888 

// Fonction qui ajoute un point de controle à un cluster
inline void cluster_add_point_controle(Cluster *cluster, Point_element *point_controle){
	Point_element *tmp_point;
	// On insert le point dans l'arbre des points de controle du cluster, et on ajoute le nombre de controle du point (renvoyé par la fonction) à la somme dans le cluster
	// (si le point existait déjà, la fonction retourne 0 et n'ajoute pas le point à l'arbre)
	cluster->nb_controle_pur += rbtree_insert_new(cluster->arbre_points_controle, point_controle->point->id, (void*)point_controle, 1);
	// On met à jour les informations du point qui vont être utilisées lors du parcours de l'arbre kd (on ajoute le cluster dans l'arbre des clusters du point)
	rbtree_insert_new(point_controle->controle_cluster, cluster->id_controle, (void*)point_controle, 1);
	// On met à jour le cluster_sous_arbre en remontant l'arbre tant que toute la partie inférieure (droite et gauche) au point courant a lui aussi le sous-arbre du même cluster (ou NULL)
	tmp_point = point_controle;
	while( 	tmp_point != NULL && rbtree_find(tmp_point->controle_cluster, cluster->id_controle) != NULL &&
			(tmp_point->kd_suivant_gauche == NULL || rbtree_find(tmp_point->kd_suivant_gauche->controle_cluster_sous_arbre, cluster->id_controle)!=NULL ) &&
			(tmp_point->kd_suivant_droite == NULL || rbtree_find(tmp_point->kd_suivant_droite->controle_cluster_sous_arbre, cluster->id_controle)!=NULL ) )
	{
		rbtree_insert_new(tmp_point->controle_cluster_sous_arbre, cluster->id_controle, (void*)tmp_point, 1);
		// On remonte l'arbre
		tmp_point = tmp_point->kd_parent;
	}
}

// 	888888 88   88 .dP"Y8 88  dP"Yb  88b 88            88""Yb 88""Yb 888888 88""Yb 888888 888888            8888b.  88 .dP"Y8 888888    db    88b 88  dP""b8 888888 
// 	88__   88   88 `Ybo." 88 dP   Yb 88Yb88            88__dP 88__dP   88   88__dP 88__   88__               8I  Yb 88 `Ybo."   88     dPYb   88Yb88 dP   `" 88__   
// 	88""   Y8   8P o.`Y8b 88 Yb   dP 88 Y88            88"Yb  88""Yb   88   88"Yb  88""   88""               8I  dY 88 o.`Y8b   88    dP__Yb  88 Y88 Yb      88""   
// 	88     `YbodP' 8bodP' 88  YbodP  88  Y8 oooooooooo 88  Yb 88oodP   88   88  Yb 888888 888888 oooooooooo 8888Y"  88 8bodP'   88   dP""""Yb 88  Y8  YboodP 888888 

// Fonction qui déplace toutes les distances de l'arbre 2 dans l'arbre 1, puis libère la mémoire de l'arbre 2
inline void fusion_rbtree_distance(rbtree tree1, rbtree tree2){
	if(tree2->root != NULL)
		recur_fusion_rbtree_distance(tree1, tree2->root);
	free(tree2);
}

// 	888888 88   88 .dP"Y8 88  dP"Yb  88b 88            88""Yb 88""Yb 888888 88""Yb 888888 888888             dP""b8  dP"Yb  88b 88 888888 88""Yb  dP"Yb  88     888888 
// 	88__   88   88 `Ybo." 88 dP   Yb 88Yb88            88__dP 88__dP   88   88__dP 88__   88__              dP   `" dP   Yb 88Yb88   88   88__dP dP   Yb 88     88__   
// 	88""   Y8   8P o.`Y8b 88 Yb   dP 88 Y88            88"Yb  88""Yb   88   88"Yb  88""   88""              Yb      Yb   dP 88 Y88   88   88"Yb  Yb   dP 88  .o 88""   
// 	88     `YbodP' 8bodP' 88  YbodP  88  Y8 oooooooooo 88  Yb 88oodP   88   88  Yb 888888 888888 oooooooooo  YboodP  YbodP  88  Y8   88   88  Yb  YbodP  88ood8 888888 

// Fonction  qui fusionne l'arbre de controle de l'ancien cluster dans l'arbre de controle du nouveau cluster
// libère l'arbre de controle de l'ancien cluster
inline void fusion_rbtree_controle(Cluster *new_cluster, Cluster *old_cluster){
	if(old_cluster->arbre_points_controle->root != NULL)
		recur_fusion_rbtree_controle(old_cluster->arbre_points_controle->root, new_cluster, old_cluster);
	free(old_cluster->arbre_points_controle);
}

// 	888888 88   88 .dP"Y8 88  dP"Yb  88b 88             dP""b8 88     88   88 .dP"Y8 888888 888888 88""Yb 
// 	88__   88   88 `Ybo." 88 dP   Yb 88Yb88            dP   `" 88     88   88 `Ybo."   88   88__   88__dP 
// 	88""   Y8   8P o.`Y8b 88 Yb   dP 88 Y88            Yb      88  .o Y8   8P o.`Y8b   88   88""   88"Yb  
// 	88     `YbodP' 8bodP' 88  YbodP  88  Y8 oooooooooo  YboodP 88ood8 `YbodP' 8bodP'   88   888888 88  Yb 

// Fonction qui fusionne le plus petit cluster dans le plus grand
inline Cluster *fusion_cluster(Cluster *cluster1, Cluster *cluster2){
	Cluster *tmp_cluster;
	rbtree tmp_tree;
	Point_element *tmp_point;
	DISTANCE_TYPE tmp_id_controle;
	// On s'assure que le cluster1 est bien le plus grand
	if(cluster1->taille < cluster2->taille)
	{
		// On échange les deux clusters
		tmp_cluster = cluster1;
		cluster1 = cluster2;
		cluster2 = tmp_cluster;
	}
	// On ajoute tous les points du 2ème cluster dans le premier
	// Comme on ajoute toujours le plus petit dans le plus grand, on a en tout (n/2)*(log(n)-1) points ajoutés de la sorte dans tout le programme
	tmp_point = cluster2->premier_point;
	while(tmp_point != NULL)
	{
		cluster_add_point(cluster1, tmp_point, NULL, 0);	// On ajoute pas les distances au cluster, vu qu'on va le faire en fusionnant les arbres distance cas-controle
		tmp_point = tmp_point->suivant_cluster;
	}
	// On fusionne les arbres arbre_distance_cas_controle (l'arbre du cluster1 est le plus grand car il a plus de cas)
	fusion_rbtree_distance(cluster1->arbre_distance_cas_controle, cluster2->arbre_distance_cas_controle);
	// On fusionne le plus petit arbre arbre_points_controle dans le plus grand (arbre des points de controle dans l'entourage du cluster)
	if(cluster1->arbre_points_controle->size < cluster2->arbre_points_controle->size)
	{
		// On échange les arbres
		tmp_tree = cluster1->arbre_points_controle;
		cluster1->arbre_points_controle = cluster2->arbre_points_controle;
		cluster2->arbre_points_controle = tmp_tree;
		// On échange les identifiants correspondants
		tmp_id_controle = cluster1->id_controle;
		cluster1->id_controle = cluster2->id_controle;
		cluster2->id_controle = tmp_id_controle;
		// On recopie aussi le nombre de controles purs (pas besoin de les échanger)
		cluster1->nb_controle_pur = cluster2->nb_controle_pur;	
	}
	// La fonction de fusion récupère le nombre de controle effectivement ajoutés (en enlevant les doubles), qu'on ajoute au nombre de controles purs, en passant par la fonction add_point_controle
	fusion_rbtree_controle(cluster1, cluster2);	
	// On libère la mémoire du cluster 2 (les arbres ont déjà été vidés)
	free(cluster2);
	return cluster1;
}


/*****************************************************/
/************ Manipulation des arbres kd *************/
/*****************************************************/

// 	 dP""b8 888888 888888            88""Yb 88     88   88 .dP"Y8            88""Yb 88""Yb  dP"Yb   dP""b8 88  88 888888            Yb    dP  dP"Yb  88 .dP"Y8 88 88b 88 
// 	dP   `" 88__     88              88__dP 88     88   88 `Ybo."            88__dP 88__dP dP   Yb dP   `" 88  88 88__               Yb  dP  dP   Yb 88 `Ybo." 88 88Yb88 
// 	Yb  "88 88""     88              88"""  88  .o Y8   8P o.`Y8b            88"""  88"Yb  Yb   dP Yb      888888 88""                YbdP   Yb   dP 88 o.`Y8b 88 88 Y88 
// 	 YboodP 888888   88   oooooooooo 88     88ood8 `YbodP' 8bodP' oooooooooo 88     88  Yb  YbodP   YboodP 88  88 888888 oooooooooo    YP     YbodP  88 8bodP' 88 88  Y8 

// Fonction récursive qui parcours l'arbre (ou sous-arbre) kd pour renvoyer le plus proche voisin du point passé en paramètre
// arbre_kd est le sous-arbre sur lequel on travaille, point le point dont on cherche le plus proche voisin, is_controle est un booléen pour savoir
// si l'arbre est un arbre de points de controle ou non, nb_dimension est le nombre de dimension
inline Distance *get_plus_proche_voisin(Point_element *arbre_kd, Point_element *point, int is_controle, ID_TYPE nb_dimension){
	Distance *meilleure_distance = calloc(1,sizeof(Distance));
	VERIFY_ALLOC(meilleure_distance);
	meilleure_distance->point1 = point;
	meilleure_distance->valeur = MAX_DISTANCE; 
	if(is_controle == 1)
		recur_get_plus_proche_voisin_controle(arbre_kd, point, meilleure_distance, nb_dimension, 0);
	else
		recur_get_plus_proche_voisin_cas(arbre_kd, point, meilleure_distance, nb_dimension, 0);
	return meilleure_distance;
}
















