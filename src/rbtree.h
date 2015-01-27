// ** ================================================================================================================
// ** R-Package: graphscan
// ** Fichier : src/rbtree.c
// ** Description : Création et manipulations des arbres rouges-noirs.
// **               Modifications, simplifications et optimisation pour du calcul haute performance.
// ** En particulier, ces arbres vont être utilisés dans 2 occasions: des arbres de distance où l'on cherche 
// ** à faciliter la lecture/suppression de la distance la plus petite, tout en conservant le temps d'insertion 
// ** d'une nouvelle distance en log n, et des arbres de pointeurs (casté en distance long double) 
// ** où l'on va juste chercher à savoir si un élément existe déjà où pas.
// **
// ** License : GPL-2 | GPL-3
// ** The authors of this work have released all rights to it and placed it
// ** in the public domain under the Creative Commons CC0 1.0 waiver
// ** (http://creativecommons.org/publicdomain/zero/1.0/).
// **
// **THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// **EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// **MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// **IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
// **CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
// **TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
// **SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
// **
// ** URL : http://en.literateprograms.org/Red-black_tree_(C)?oldid=19567
// **
// **
// ** Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
// ** ================================================================================================================

#ifndef _RBTREE_H_
#define _RBTREE_H_


#include <stdlib.h>
#include "2D_types.h"

//     .88b  d88.  .d8b.   .o88b. d8888b.  .d88b.  .d8888. 
//     88'YbdP`88 d8' `8b d8P  Y8 88  `8D .8P  Y8. 88'  YP 
//     88  88  88 88ooo88 8P      88oobY' 88    88 `8bo.   
//     88  88  88 88~~~88 8b      88`8b   88    88   `Y8b. 
//     88  88  88 88   88 Y8b  d8 88 `88. `8b  d8' db   8D 
//     YP  YP  YP YP   YP  `Y88P' 88   YD  `Y88P'  `8888Y' 
//                                                         
//                                                         

// Macro pour éviter les appels de fonction

// Retourne la première valeur du node avec la plus petite clé (void*)
#define RBTREE_READ_SMALLEST(t) ( (t)->smallest != NULL ? (t)->smallest->value_first->value : NULL )    

// Lire la première valeur d'un noeud comme un pointeur vers un point
#define RBTREE_GET_AS_CONTROLE(n) ( (Point_element*)((n)->value_first->value) )


//     d8888b. d88888b  .o88b. db       .d8b.  d8888b.  .d8b.  d888888b d888888b  .d88b.  d8b   db .d8888. 
//     88  `8D 88'     d8P  Y8 88      d8' `8b 88  `8D d8' `8b `~~88~~'   `88'   .8P  Y8. 888o  88 88'  YP 
//     88   88 88ooooo 8P      88      88ooo88 88oobY' 88ooo88    88       88    88    88 88V8o 88 `8bo.   
//     88   88 88~~~~~ 8b      88      88~~~88 88`8b   88~~~88    88       88    88    88 88 V8o88   `Y8b. 
//     88  .8D 88.     Y8b  d8 88booo. 88   88 88 `88. 88   88    88      .88.   `8b  d8' 88  V888 db   8D 
//     Y8888D' Y88888P  `Y88P' Y88888P YP   YP 88   YD YP   YP    YP    Y888888P  `Y88P'  VP   V8P `8888Y' 
//                                                                                                         
//                                                                                                         

// Fonction à ne pas appeler directement, elles seront appelées le cas échéant par les inline qui suivent
void insert_case1(rbtree t, rbtree_node n);
void delete_case1(rbtree t, rbtree_node n);
// Fonction récursive pour libérer les éventuels enfants des noeuds avant de libérer les noeuds eux-même
// free_value est un booléen qui permet de savoir s'il faut libérer la ou les valeurs stockées dans le noeud (Distance seulement)
void free_node_and_sons(rbtree t, rbtree_node n, int free_value);

// Création d'un arbre
rbtree rbtree_create();
// Fonction qui crée un nouveau node
rbtree_node rbtree_new_node(DISTANCE_TYPE key, void *value);
// Remplace oldn par newn
void rbtree_replace_node(rbtree t, rbtree_node oldn, rbtree_node newn);
// Effectue une rotation vers la gauche
void rbtree_rotate_left(rbtree t, rbtree_node n);
// Effectue une rotation vers la droite
void rbtree_rotate_right(rbtree t, rbtree_node n);
// retourne le noeud si la clé existe, NULL sinon
rbtree_node rbtree_find(rbtree t, DISTANCE_TYPE key);
// Fonction qui libère la mémoire de la liste de valeurs d'un noeud
// free_value est un booléen qui permet de savoir s'il faut libérer la ou les valeurs stockées dans le noeud (Distance seulement)
void rbtree_free_value_elem(rbtree t, rbtree_node n, int free_value);
// Fonction qui retourne le second smallest d'un arbre (à utiliser avant sa supression)
rbtree_node rbtree_get_second_smallest(rbtree t);
// fonction qui supprime le node indexé par key
// si free_value == 1, la distance pointée par value est libérée aussi
// Pas besoin de mettre à jour smallest, le type d'arbre qui appelle cette fonction ne l'utilise pas
void rbtree_delete(rbtree t, DISTANCE_TYPE key, int free_value);
// Fonction qui vide la mémoire de l'arbre
// MFree_value est un booléen qui permet de savoir s'il faut libérer la ou les valeurs stockées dans le noeud
void rbtree_free(rbtree t, int free_value);
// Fonction qui supprime la valeur avec la plus petite clé (supprime le noeud si il n'y a qu'une seule valeur)
// Cette fonction n'est appelée que pour les arbres de distance, et elle retourne le pointeur vers la valeur (distance)
void *rbtree_get_smallest(rbtree t);
// Fonction pour insérer un point dans un arbre
// is_controle est un booléen, si il est à 1 en cas de point qui existe déjà on n'ajoute pas le point, on compte le nombre de controle des points de controle qui ont été ajouté
// La fonction retourne un ID_TYPE qui est la somme des controles qui ont été ajouté (0 si is_controle==0)
ID_TYPE rbtree_insert(rbtree t, rbtree_node new_node, int is_controle);
// Fonction pour insérer un nouveau point dans l'arbre
// is_controle est un booléen, si il est à 1 en cas de point qui existe déjà on n'ajoute pas le point, on compte le nombre de controle des points de controle qui ont été ajouté
// La fonction retourne un ID_TYPE qui est le nombre de controles qui ont été ajoutés (0 si is_controle==0)
ID_TYPE rbtree_insert_new(rbtree t, DISTANCE_TYPE key, void *value, int is_controle);
// Renvoie le frère de node
rbtree_node rbtree_sibling(rbtree_node node);
// Renvoie la couleur de node
rbtree_color rbtree_node_color(rbtree_node node);
// Retourne le grand-parent d'un noeud
rbtree_node rbtree_grand_parent(rbtree_node node);
// Retourne l'oncle d'un noeud
rbtree_node rbtree_uncle(rbtree_node node);





#endif
