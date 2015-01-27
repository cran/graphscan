// ** ================================================================================================================
// ** R-Package: graphscan
// ** Fichier : src/rbtree.c
// ** Description : Création et manipulations des arbres rouges-noirs.
// **               Modifications effectuées pour du calcul haute performance.
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


#include "rbtree.h"

//     d8888b. d88888b  .o88b. db       .d8b.  d8888b.  .d8b.  d888888b d888888b  .d88b.  d8b   db .d8888. 
//     88  `8D 88'     d8P  Y8 88      d8' `8b 88  `8D d8' `8b `~~88~~'   `88'   .8P  Y8. 888o  88 88'  YP 
//     88   88 88ooooo 8P      88      88ooo88 88oobY' 88ooo88    88       88    88    88 88V8o 88 `8bo.   
//     88   88 88~~~~~ 8b      88      88~~~88 88`8b   88~~~88    88       88    88    88 88 V8o88   `Y8b. 
//     88  .8D 88.     Y8b  d8 88booo. 88   88 88 `88. 88   88    88      .88.   `8b  d8' 88  V888 db   8D 
//     Y8888D' Y88888P  `Y88P' Y88888P YP   YP 88   YD YP   YP    YP    Y888888P  `Y88P'  VP   V8P `8888Y' 
//                                                                                                         
//                                                                                                         

// Fonctions qui resteront locales
void insert_case3(rbtree t, rbtree_node n);
void delete_case3(rbtree t, rbtree_node n);

void insert_case2(rbtree t, rbtree_node n);
void insert_case4(rbtree t, rbtree_node n);
void insert_case5(rbtree t, rbtree_node n);

void delete_case2(rbtree t, rbtree_node n);
void delete_case4(rbtree t, rbtree_node n);
void delete_case5(rbtree t, rbtree_node n);
void delete_case6(rbtree t, rbtree_node n);

//     d8888b. d88888b d88888b d888888b d8b   db d888888b d888888b d888888b  .d88b.  d8b   db .d8888. 
//     88  `8D 88'     88'       `88'   888o  88   `88'   `~~88~~'   `88'   .8P  Y8. 888o  88 88'  YP 
//     88   88 88ooooo 88ooo      88    88V8o 88    88       88       88    88    88 88V8o 88 `8bo.   
//     88   88 88~~~~~ 88~~~      88    88 V8o88    88       88       88    88    88 88 V8o88   `Y8b. 
//     88  .8D 88.     88        .88.   88  V888   .88.      88      .88.   `8b  d8' 88  V888 db   8D 
//     Y8888D' Y88888P YP      Y888888P VP   V8P Y888888P    YP    Y888888P  `Y88P'  VP   V8P `8888Y' 
//                                                                                                    
//                                                                                                    

//     8888b.  888888 88     888888 888888 888888             dP""b8    db    .dP"Y8 888888   .d 
//      8I  Yb 88__   88     88__     88   88__              dP   `"   dPYb   `Ybo." 88__   .d88 
//      8I  dY 88""   88  .o 88""     88   88""              Yb       dP__Yb  o.`Y8b 88""     88 
//     8888Y"  888888 88ood8 888888   88   888888 oooooooooo  YboodP dP""""Yb 8bodP' 888888   88 

void delete_case1(rbtree t, rbtree_node n) {
    if (n->parent != NULL && rbtree_sibling(n) != NULL)
    {
        delete_case2(t, n);
    }
}

//     8888b.  888888 88     888888 888888 888888             dP""b8    db    .dP"Y8 888888 oP"Yb. 
//      8I  Yb 88__   88     88__     88   88__              dP   `"   dPYb   `Ybo." 88__   "' dP' 
//      8I  dY 88""   88  .o 88""     88   88""              Yb       dP__Yb  o.`Y8b 88""     dP'  
//     8888Y"  888888 88ood8 888888   88   888888 oooooooooo  YboodP dP""""Yb 8bodP' 888888 .d8888 

void delete_case2(rbtree t, rbtree_node n) {
    if (rbtree_node_color(rbtree_sibling(n)) == RED) {
        n->parent->color = RED;
        rbtree_sibling(n)->color = BLACK;
        if (n == n->parent->left)
            rbtree_rotate_left(t, n->parent);
        else
            rbtree_rotate_right(t, n->parent);
    }
    delete_case3(t, n);
}

//     8888b.  888888 88     888888 888888 888888             dP""b8    db    .dP"Y8 888888 88888 
//      8I  Yb 88__   88     88__     88   88__              dP   `"   dPYb   `Ybo." 88__     .dP 
//      8I  dY 88""   88  .o 88""     88   88""              Yb       dP__Yb  o.`Y8b 88""   o `Yb 
//     8888Y"  888888 88ood8 888888   88   888888 oooooooooo  YboodP dP""""Yb 8bodP' 888888 YbodP 

void delete_case3(rbtree t, rbtree_node n) {
    if (rbtree_node_color(n->parent) == BLACK &&
        rbtree_node_color(rbtree_sibling(n)) == BLACK &&
        rbtree_node_color(rbtree_sibling(n)->left) == BLACK &&
        rbtree_node_color(rbtree_sibling(n)->right) == BLACK)
    {
        rbtree_sibling(n)->color = RED;
        delete_case1(t, n->parent);
    }
    else
        delete_case4(t, n);
}

//     8888b.  888888 88     888888 888888 888888             dP""b8    db    .dP"Y8 888888   dP88  
//      8I  Yb 88__   88     88__     88   88__              dP   `"   dPYb   `Ybo." 88__    dP 88  
//      8I  dY 88""   88  .o 88""     88   88""              Yb       dP__Yb  o.`Y8b 88""   d888888 
//     8888Y"  888888 88ood8 888888   88   888888 oooooooooo  YboodP dP""""Yb 8bodP' 888888     88  

void delete_case4(rbtree t, rbtree_node n) {
    if(rbtree_sibling(n)!=NULL)
    {
        if (rbtree_node_color(n->parent) == RED &&
            rbtree_node_color(rbtree_sibling(n)) == BLACK &&
            rbtree_node_color(rbtree_sibling(n)->left) == BLACK &&
            rbtree_node_color(rbtree_sibling(n)->right) == BLACK)
        {
            rbtree_sibling(n)->color = RED;
            n->parent->color = BLACK;
        }
        else
            delete_case5(t, n);
    }
}

//     8888b.  888888 88     888888 888888 888888             dP""b8    db    .dP"Y8 888888 888888 
//      8I  Yb 88__   88     88__     88   88__              dP   `"   dPYb   `Ybo." 88__   88oo." 
//      8I  dY 88""   88  .o 88""     88   88""              Yb       dP__Yb  o.`Y8b 88""      `8b 
//     8888Y"  888888 88ood8 888888   88   888888 oooooooooo  YboodP dP""""Yb 8bodP' 888888 8888P' 

void delete_case5(rbtree t, rbtree_node n) {
    if (n == n->parent->left &&
        rbtree_node_color(rbtree_sibling(n)) == BLACK &&
        rbtree_node_color(rbtree_sibling(n)->left) == RED &&
        rbtree_node_color(rbtree_sibling(n)->right) == BLACK)
    {
        rbtree_sibling(n)->color = RED;
        rbtree_sibling(n)->left->color = BLACK;
        rbtree_rotate_right(t, rbtree_sibling(n));
    }
    else if (n == n->parent->right &&
             rbtree_node_color(rbtree_sibling(n)) == BLACK &&
             rbtree_node_color(rbtree_sibling(n)->right) == RED &&
             rbtree_node_color(rbtree_sibling(n)->left) == BLACK)
    {
        rbtree_sibling(n)->color = RED;
        rbtree_sibling(n)->right->color = BLACK;
        rbtree_rotate_left(t, rbtree_sibling(n));
    }
    delete_case6(t, n);
}

//     8888b.  888888 88     888888 888888 888888             dP""b8    db    .dP"Y8 888888   dP'   
//      8I  Yb 88__   88     88__     88   88__              dP   `"   dPYb   `Ybo." 88__   .d8'    
//      8I  dY 88""   88  .o 88""     88   88""              Yb       dP__Yb  o.`Y8b 88""   8P"""Yb 
//     8888Y"  888888 88ood8 888888   88   888888 oooooooooo  YboodP dP""""Yb 8bodP' 888888 `YboodP 

void delete_case6(rbtree t, rbtree_node n) {
    rbtree_sibling(n)->color = rbtree_node_color(n->parent);
    n->parent->color = BLACK;
    if (n == n->parent->left) {
        rbtree_sibling(n)->right->color = BLACK;
        rbtree_rotate_left(t, n->parent);
    }
    else
    {
        rbtree_sibling(n)->left->color = BLACK;
        rbtree_rotate_right(t, n->parent);
    }
}

//     88 88b 88 .dP"Y8 888888 88""Yb 888888             dP""b8    db    .dP"Y8 888888   .d 
//     88 88Yb88 `Ybo." 88__   88__dP   88              dP   `"   dPYb   `Ybo." 88__   .d88 
//     88 88 Y88 o.`Y8b 88""   88"Yb    88              Yb       dP__Yb  o.`Y8b 88""     88 
//     88 88  Y8 8bodP' 888888 88  Yb   88   oooooooooo  YboodP dP""""Yb 8bodP' 888888   88 

void insert_case1(rbtree t, rbtree_node n) {
    if (n->parent == NULL)
        n->color = BLACK;
    else
        insert_case2(t, n);
}

//     88 88b 88 .dP"Y8 888888 88""Yb 888888             dP""b8    db    .dP"Y8 888888 oP"Yb. 
//     88 88Yb88 `Ybo." 88__   88__dP   88              dP   `"   dPYb   `Ybo." 88__   "' dP' 
//     88 88 Y88 o.`Y8b 88""   88"Yb    88              Yb       dP__Yb  o.`Y8b 88""     dP'  
//     88 88  Y8 8bodP' 888888 88  Yb   88   oooooooooo  YboodP dP""""Yb 8bodP' 888888 .d8888 

void insert_case2(rbtree t, rbtree_node n) {
    if (rbtree_node_color(n->parent) != BLACK)
        insert_case3(t, n);
}

//     88 88b 88 .dP"Y8 888888 88""Yb 888888             dP""b8    db    .dP"Y8 888888 88888 
//     88 88Yb88 `Ybo." 88__   88__dP   88              dP   `"   dPYb   `Ybo." 88__     .dP 
//     88 88 Y88 o.`Y8b 88""   88"Yb    88              Yb       dP__Yb  o.`Y8b 88""   o `Yb 
//     88 88  Y8 8bodP' 888888 88  Yb   88   oooooooooo  YboodP dP""""Yb 8bodP' 888888 YbodP 

void insert_case3(rbtree t, rbtree_node n) {
    if (rbtree_node_color(rbtree_uncle(n)) == RED) {
        n->parent->color = BLACK;
        rbtree_uncle(n)->color = BLACK;
        rbtree_grand_parent(n)->color = RED;
        insert_case1(t, rbtree_grand_parent(n));
    } else {
        insert_case4(t, n);
    }
}

//     88 88b 88 .dP"Y8 888888 88""Yb 888888             dP""b8    db    .dP"Y8 888888   dP88  
//     88 88Yb88 `Ybo." 88__   88__dP   88              dP   `"   dPYb   `Ybo." 88__    dP 88  
//     88 88 Y88 o.`Y8b 88""   88"Yb    88              Yb       dP__Yb  o.`Y8b 88""   d888888 
//     88 88  Y8 8bodP' 888888 88  Yb   88   oooooooooo  YboodP dP""""Yb 8bodP' 888888     88  

void insert_case4(rbtree t, rbtree_node n) {
    if (n == n->parent->right && n->parent == rbtree_grand_parent(n)->left) {
        rbtree_rotate_left(t, n->parent);
        n = n->left;
    } else if (n == n->parent->left && n->parent == rbtree_grand_parent(n)->right) {
        rbtree_rotate_right(t, n->parent);
        n = n->right;
    }
    insert_case5(t, n);
}

//     88 88b 88 .dP"Y8 888888 88""Yb 888888             dP""b8    db    .dP"Y8 888888 888888 
//     88 88Yb88 `Ybo." 88__   88__dP   88              dP   `"   dPYb   `Ybo." 88__   88oo." 
//     88 88 Y88 o.`Y8b 88""   88"Yb    88              Yb       dP__Yb  o.`Y8b 88""      `8b 
//     88 88  Y8 8bodP' 888888 88  Yb   88   oooooooooo  YboodP dP""""Yb 8bodP' 888888 8888P' 

void insert_case5(rbtree t, rbtree_node n) {
    n->parent->color = BLACK;
    rbtree_grand_parent(n)->color = RED;
    if (n == n->parent->left && n->parent == rbtree_grand_parent(n)->left) {
        rbtree_rotate_right(t, rbtree_grand_parent(n));
    } else {
        rbtree_rotate_left(t, rbtree_grand_parent(n));        
    }
}

//     888888 88""Yb 888888 888888            88b 88  dP"Yb  8888b.  888888               db    88b 88 8888b.             .dP"Y8  dP"Yb  88b 88 .dP"Y8 
//     88__   88__dP 88__   88__              88Yb88 dP   Yb  8I  Yb 88__                dPYb   88Yb88  8I  Yb            `Ybo." dP   Yb 88Yb88 `Ybo." 
//     88""   88"Yb  88""   88""              88 Y88 Yb   dP  8I  dY 88""               dP__Yb  88 Y88  8I  dY            o.`Y8b Yb   dP 88 Y88 o.`Y8b 
//     88     88  Yb 888888 888888 oooooooooo 88  Y8  YbodP  8888Y"  888888 oooooooooo dP""""Yb 88  Y8 8888Y"  oooooooooo 8bodP'  YbodP  88  Y8 8bodP' 

// Fonction récursive pour libérer les éventuels enfants des noeuds avant de libérer les noeuds eux-même
// free_value est un booléen qui permet de savoir s'il faut libérer la ou les valeurs stockées dans le noeud (Distance seulement)
void free_node_and_sons(rbtree t, rbtree_node n, int free_value){
    if(n->right!=NULL)
        free_node_and_sons(t, n->right, free_value);
    if(n->left!=NULL)
        free_node_and_sons(t, n->left, free_value);
    rbtree_free_value_elem(t, n, free_value);
    free(n);    // Libération du noeud
}
                                                                                                

//     88""Yb 88""Yb 888888 88""Yb 888888 888888            .dP"Y8 88 88""Yb 88     88 88b 88  dP""b8 
//     88__dP 88__dP   88   88__dP 88__   88__              `Ybo." 88 88__dP 88     88 88Yb88 dP   `" 
//     88"Yb  88""Yb   88   88"Yb  88""   88""              o.`Y8b 88 88""Yb 88  .o 88 88 Y88 Yb  "88 
//     88  Yb 88oodP   88   88  Yb 888888 888888 oooooooooo 8bodP' 88 88oodP 88ood8 88 88  Y8  YboodP 

// Renvoie le frère de node
rbtree_node rbtree_sibling(rbtree_node node){
    //fprintf(get_debug(), "sibling\n");
    if (node == node->parent->left)
        return node->parent->right;
    else
        return node->parent->left;
}

//     88""Yb 88""Yb 888888 88""Yb 888888 888888            88b 88  dP"Yb  8888b.  888888             dP""b8  dP"Yb  88      dP"Yb  88""Yb 
//     88__dP 88__dP   88   88__dP 88__   88__              88Yb88 dP   Yb  8I  Yb 88__              dP   `" dP   Yb 88     dP   Yb 88__dP 
//     88"Yb  88""Yb   88   88"Yb  88""   88""              88 Y88 Yb   dP  8I  dY 88""              Yb      Yb   dP 88  .o Yb   dP 88"Yb  
//     88  Yb 88oodP   88   88  Yb 888888 888888 oooooooooo 88  Y8  YbodP  8888Y"  888888 oooooooooo  YboodP  YbodP  88ood8  YbodP  88  Yb 

// Renvoie la couleur de node
rbtree_color rbtree_node_color(rbtree_node node){
    //fprintf(get_debug(), "couleur\n");
    return node == NULL ? BLACK : node->color;
}

//     88""Yb 88""Yb 888888 88""Yb 888888 888888             dP""b8 88""Yb    db    88b 88 8888b.             88""Yb    db    88""Yb 888888 88b 88 888888 
//     88__dP 88__dP   88   88__dP 88__   88__              dP   `" 88__dP   dPYb   88Yb88  8I  Yb            88__dP   dPYb   88__dP 88__   88Yb88   88   
//     88"Yb  88""Yb   88   88"Yb  88""   88""              Yb  "88 88"Yb   dP__Yb  88 Y88  8I  dY            88"""   dP__Yb  88"Yb  88""   88 Y88   88   
//     88  Yb 88oodP   88   88  Yb 888888 888888 oooooooooo  YboodP 88  Yb dP""""Yb 88  Y8 8888Y"  oooooooooo 88     dP""""Yb 88  Yb 888888 88  Y8   88   

// Retourne le grand-parent d'un noeud
rbtree_node rbtree_grand_parent(rbtree_node node){
    //fprintf(get_debug(), "gd_parent\n");
    return node->parent->parent;
}

//     88""Yb 88""Yb 888888 88""Yb 888888 888888            88   88 88b 88  dP""b8 88     888888 
//     88__dP 88__dP   88   88__dP 88__   88__              88   88 88Yb88 dP   `" 88     88__   
//     88"Yb  88""Yb   88   88"Yb  88""   88""              Y8   8P 88 Y88 Yb      88  .o 88""   
//     88  Yb 88oodP   88   88  Yb 888888 888888 oooooooooo `YbodP' 88  Y8  YboodP 88ood8 888888 

// Retourne l'oncle d'un noeud
rbtree_node rbtree_uncle(rbtree_node node){
    //fprintf(get_debug(), "uncle\n");
    return rbtree_sibling(node->parent);
}


//     88""Yb 88""Yb 888888 88""Yb 888888 888888             dP""b8 88""Yb 888888    db    888888 888888 
//     88__dP 88__dP   88   88__dP 88__   88__              dP   `" 88__dP 88__     dPYb     88   88__   
//     88"Yb  88""Yb   88   88"Yb  88""   88""              Yb      88"Yb  88""    dP__Yb    88   88""   
//     88  Yb 88oodP   88   88  Yb 888888 888888 oooooooooo  YboodP 88  Yb 888888 dP""""Yb   88   888888 

// Allocation d'un nouvel arbre
rbtree rbtree_create(){
    rbtree new_tree = calloc(1, sizeof(struct rbtree_t));
    VERIFY_ALLOC(new_tree);
    return new_tree;
}

//     88""Yb 88""Yb 888888 88""Yb 888888 888888            88b 88 888888 Yb        dP            88b 88  dP"Yb  8888b.  888888 
//     88__dP 88__dP   88   88__dP 88__   88__              88Yb88 88__    Yb  db  dP             88Yb88 dP   Yb  8I  Yb 88__   
//     88"Yb  88""Yb   88   88"Yb  88""   88""              88 Y88 88""     YbdPYbdP              88 Y88 Yb   dP  8I  dY 88""   
//     88  Yb 88oodP   88   88  Yb 888888 888888 oooooooooo 88  Y8 888888    YP  YP    oooooooooo 88  Y8  YbodP  8888Y"  888888 

// Fonction qui crée un nouveau node
rbtree_node rbtree_new_node(DISTANCE_TYPE key, void *value){
    rbtree_node new_node = calloc(1,sizeof(struct rbtree_node_t));
    VERIFY_ALLOC(new_node);
    new_node->value_first = calloc(1,sizeof(struct rbtree_value_elem_t));
    VERIFY_ALLOC(new_node->value_first);
    new_node->value_first->value = value;
    new_node->key = key;
    new_node->color = RED;
    return new_node;
}

//     88""Yb 88""Yb 888888 88""Yb 888888 888888            88""Yb 888888 88""Yb 88        db     dP""b8 888888            88b 88  dP"Yb  8888b.  888888 
//     88__dP 88__dP   88   88__dP 88__   88__              88__dP 88__   88__dP 88       dPYb   dP   `" 88__              88Yb88 dP   Yb  8I  Yb 88__   
//     88"Yb  88""Yb   88   88"Yb  88""   88""              88"Yb  88""   88"""  88  .o  dP__Yb  Yb      88""              88 Y88 Yb   dP  8I  dY 88""   
//     88  Yb 88oodP   88   88  Yb 888888 888888 oooooooooo 88  Yb 888888 88     88ood8 dP""""Yb  YboodP 888888 oooooooooo 88  Y8  YbodP  8888Y"  888888 

// Remplace oldn par newn
void rbtree_replace_node(rbtree t, rbtree_node oldn, rbtree_node newn){
    //fprintf(get_debug(), "replace\n");
    if (oldn->parent == NULL) {
        t->root = newn;
    } else {
        if (oldn == oldn->parent->left)
            oldn->parent->left = newn;
        else
            oldn->parent->right = newn;
    }
    if (newn != NULL) {
        newn->parent = oldn->parent;
    }
}

//     88""Yb 88""Yb 888888 88""Yb 888888 888888            88""Yb  dP"Yb  888888    db    888888 888888            88     888888 888888 888888 
//     88__dP 88__dP   88   88__dP 88__   88__              88__dP dP   Yb   88     dPYb     88   88__              88     88__   88__     88   
//     88"Yb  88""Yb   88   88"Yb  88""   88""              88"Yb  Yb   dP   88    dP__Yb    88   88""              88  .o 88""   88""     88   
//     88  Yb 88oodP   88   88  Yb 888888 888888 oooooooooo 88  Yb  YbodP    88   dP""""Yb   88   888888 oooooooooo 88ood8 888888 88       88   

// Effectue une rotation vers la gauche
void rbtree_rotate_left(rbtree t, rbtree_node n){
    //fprintf(get_debug(), "rotate left\n");
    rbtree_node r = n->right;
    rbtree_replace_node(t, n, r);
    n->right = r->left;
    if (r->left != NULL)
        r->left->parent = n;
    r->left = n;
    n->parent = r;   
}

//     88""Yb 88""Yb 888888 88""Yb 888888 888888            88""Yb  dP"Yb  888888    db    888888 888888            88""Yb 88  dP""b8 88  88 888888 
//     88__dP 88__dP   88   88__dP 88__   88__              88__dP dP   Yb   88     dPYb     88   88__              88__dP 88 dP   `" 88  88   88   
//     88"Yb  88""Yb   88   88"Yb  88""   88""              88"Yb  Yb   dP   88    dP__Yb    88   88""              88"Yb  88 Yb  "88 888888   88   
//     88  Yb 88oodP   88   88  Yb 888888 888888 oooooooooo 88  Yb  YbodP    88   dP""""Yb   88   888888 oooooooooo 88  Yb 88  YboodP 88  88   88   

// Effectue une rotation vers la droite
void rbtree_rotate_right(rbtree t, rbtree_node n){
    //fprintf(get_debug(), "rotate right\n");
    rbtree_node L = n->left;
    rbtree_replace_node(t, n, L);
    n->left = L->right;
    if (L->right != NULL)
        L->right->parent = n;
    L->right = n;
    n->parent = L;
}

//     88""Yb 88""Yb 888888 88""Yb 888888 888888            888888 88 88b 88 8888b.  
//     88__dP 88__dP   88   88__dP 88__   88__              88__   88 88Yb88  8I  Yb 
//     88"Yb  88""Yb   88   88"Yb  88""   88""              88""   88 88 Y88  8I  dY 
//     88  Yb 88oodP   88   88  Yb 888888 888888 oooooooooo 88     88 88  Y8 8888Y"  

// retourne le noeud si la clé existe, NULL sinon
rbtree_node rbtree_find(rbtree t, DISTANCE_TYPE key){
    rbtree_node n;
    n = t->root;
    while (n != NULL && key != n->key)
    {
        if(key < n->key)
            n = n->left;
        else if(key > n->key)
            n = n->right;
    }
    return n;
}

//     88""Yb 88""Yb 888888 88""Yb 888888 888888            888888 88""Yb 888888 888888            Yb    dP    db    88     88   88 888888            888888 88     888888 8b    d8 
//     88__dP 88__dP   88   88__dP 88__   88__              88__   88__dP 88__   88__               Yb  dP    dPYb   88     88   88 88__              88__   88     88__   88b  d88 
//     88"Yb  88""Yb   88   88"Yb  88""   88""              88""   88"Yb  88""   88""                YbdP    dP__Yb  88  .o Y8   8P 88""              88""   88  .o 88""   88YbdP88 
//     88  Yb 88oodP   88   88  Yb 888888 888888 oooooooooo 88     88  Yb 888888 888888 oooooooooo    YP    dP""""Yb 88ood8 `YbodP' 888888 oooooooooo 888888 88ood8 888888 88 YY 88 

// Fonction qui libère la mémoire de la liste de valeurs d'un noeud
// free_value est un booléen qui permet de savoir s'il faut libérer la ou les valeurs stockées dans le noeud (Distance seulement)
void rbtree_free_value_elem(rbtree t, rbtree_node n, int free_value){
    rbtree_value_elem tmp_elem, tmp_elem_next;
    tmp_elem = n->value_first;
    while(tmp_elem != NULL)
    { 
        tmp_elem_next = tmp_elem->next;
        if(free_value==1)      
            free((Distance*)tmp_elem->value);  // C'est un arbre de distances, donc on supprime les distances
        free(tmp_elem);
        t->size--;
        tmp_elem = tmp_elem_next;
    }
    n->value_first = NULL;
}

//     88""Yb 88""Yb 888888 88""Yb 888888 888888             dP""b8 888888 888888            .dP"Y8 888888  dP""b8  dP"Yb  88b 88 8888b.             .dP"Y8 8b    d8    db    88     88     888888 .dP"Y8 888888 
//     88__dP 88__dP   88   88__dP 88__   88__              dP   `" 88__     88              `Ybo." 88__   dP   `" dP   Yb 88Yb88  8I  Yb            `Ybo." 88b  d88   dPYb   88     88     88__   `Ybo."   88   
//     88"Yb  88""Yb   88   88"Yb  88""   88""              Yb  "88 88""     88              o.`Y8b 88""   Yb      Yb   dP 88 Y88  8I  dY            o.`Y8b 88YbdP88  dP__Yb  88  .o 88  .o 88""   o.`Y8b   88   
//     88  Yb 88oodP   88   88  Yb 888888 888888 oooooooooo  YboodP 888888   88   oooooooooo 8bodP' 888888  YboodP  YbodP  88  Y8 8888Y"  oooooooooo 8bodP' 88 YY 88 dP""""Yb 88ood8 88ood8 888888 8bodP'   88   

// Fonction qui retourne le second smallest d'un arbre (à utiliser avant sa supression)
rbtree_node rbtree_get_second_smallest(rbtree t){
    rbtree_node tmp_node = t->smallest;
    if(tmp_node != NULL)
    {
        // Soit le smallest a un sous-arbre droit
        if(tmp_node->right != NULL)
        {
            tmp_node = tmp_node->right;
            while(tmp_node->left != NULL)
                tmp_node = tmp_node->left;
        }else{  // Soit le second est tout simplement le père
            tmp_node = tmp_node->parent;
        }
    }
    return tmp_node;
}

//     88""Yb 88""Yb 888888 88""Yb 888888 888888            8888b.  888888 88     888888 888888 888888 
//     88__dP 88__dP   88   88__dP 88__   88__               8I  Yb 88__   88     88__     88   88__   
//     88"Yb  88""Yb   88   88"Yb  88""   88""               8I  dY 88""   88  .o 88""     88   88""   
//     88  Yb 88oodP   88   88  Yb 888888 888888 oooooooooo 8888Y"  888888 88ood8 888888   88   888888 

// fonction qui supprime le node indexé par key
// si free_value == 1, la distance pointée par value est libérée aussi
// Pas besoin de mettre à jour smallest, le type d'arbre qui appelle cette fonction ne l'utilise pas
void rbtree_delete(rbtree t, DISTANCE_TYPE key, int free_value){
    rbtree_node child;
    rbtree_node n;
    // On recherche le noeud
    n = rbtree_find(t, key);
    //fprintf(get_debug(), "suppression de la clé : %f\n", key);
    if (n != NULL)
    {
        //fprintf(get_debug(), "suppression effective\n");
        // Si c'est le plus petit point, on met le smallest à jour
        if(n == t->smallest)
            t->smallest = rbtree_get_second_smallest(t);
        // On libère les valeurs
        rbtree_free_value_elem(t, n, free_value);
        if (n->left != NULL && n->right != NULL) {
            // On va chercher le max du sous-arbre gauche et on copie ses valeurs
            rbtree_node pred = n->left;
            // On va chercher le maximum du sous-arbre
            while (pred->right != NULL) {
                pred = pred->right;
            }
            n->key = pred->key;
            n->value_first = pred->value_first;
            // Si pred est le smallest, pred s'est déplacé dans n donc on déplace aussi le smallest
            if(pred == t->smallest)
                t->smallest = n;
            n = pred;
        }
        child = (n->right == NULL) ? n->left  : n->right;
        if (rbtree_node_color(n) == BLACK) {
            n->color = rbtree_node_color(child);
            delete_case1(t, n);
        }
        rbtree_replace_node(t, n, child);
         // Si n est la racine, on la met noire
        if (n->parent == NULL && child != NULL)
            child->color = BLACK;
        free(n);
    }
}

//     88""Yb 88""Yb 888888 88""Yb 888888 888888            888888 88""Yb 888888 888888 
//     88__dP 88__dP   88   88__dP 88__   88__              88__   88__dP 88__   88__   
//     88"Yb  88""Yb   88   88"Yb  88""   88""              88""   88"Yb  88""   88""   
//     88  Yb 88oodP   88   88  Yb 888888 888888 oooooooooo 88     88  Yb 888888 888888 

// Fonction qui vide la mémoire de l'arbre
// MFree_value est un booléen qui permet de savoir s'il faut libérer la ou les valeurs stockées dans le noeud
void rbtree_free(rbtree t, int free_value){
    if(t->root != NULL)
        free_node_and_sons(t, t->root, free_value);
    free(t);
}

//     88""Yb 88""Yb 888888 88""Yb 888888 888888             dP""b8 888888 888888            .dP"Y8 8b    d8    db    88     88     888888 .dP"Y8 888888 
//     88__dP 88__dP   88   88__dP 88__   88__              dP   `" 88__     88              `Ybo." 88b  d88   dPYb   88     88     88__   `Ybo."   88   
//     88"Yb  88""Yb   88   88"Yb  88""   88""              Yb  "88 88""     88              o.`Y8b 88YbdP88  dP__Yb  88  .o 88  .o 88""   o.`Y8b   88   
//     88  Yb 88oodP   88   88  Yb 888888 888888 oooooooooo  YboodP 888888   88   oooooooooo 8bodP' 88 YY 88 dP""""Yb 88ood8 88ood8 888888 8bodP'   88   

// qui supprime la valeur avec la plus petite clé (supprime le noeud si il n'y a qu'une seule valeur)
// Cette fonction n'est appelée que pour les arbres de distance, et elle retourne le pointeur vers la valeur (distance)
void *rbtree_get_smallest(rbtree t){
    rbtree_node child;
    rbtree_node n = t->smallest;
    rbtree_value_elem value_elem;
    void *value = NULL;
    if(n != NULL)
    {
        // On récupère la valeur
        value_elem = n->value_first;
        // On retire la valeur de la liste
        n->value_first = value_elem->next;
        value = value_elem->value;
        // On libère la mémoire du conteneur
        free(value_elem);
        // Un élément de moins dans l'arbre
        t->size--;
        // Si il y a plusieurs éléments value, pas besoin de supprimer le point, sinon:
        if(n->value_first == NULL)
        {
            // Ne pas oublier de mettre à jour la valeur de smallest
            t->smallest = rbtree_get_second_smallest(t);
            // Si le plus petit était la racine, l'arbre est maintenant vide
            if(t->smallest == NULL)
            {
                t->root = NULL;
            }else{
                // Ensuite on met à jour l'arbre (le suivant gauche est forcément == NULL)
                child = n->right ;
                if (rbtree_node_color(n) == BLACK) {
                    n->color = rbtree_node_color(child);
                    delete_case1(t, n);
                }
                rbtree_replace_node(t, n, child);
                 // root should be black
                if (n->parent == NULL && child != NULL)
                    child->color = BLACK;
            }
            free(n);
        }
    }
    return value;
}

//     88""Yb 88""Yb 888888 88""Yb 888888 888888            88 88b 88 .dP"Y8 888888 88""Yb 888888 
//     88__dP 88__dP   88   88__dP 88__   88__              88 88Yb88 `Ybo." 88__   88__dP   88   
//     88"Yb  88""Yb   88   88"Yb  88""   88""              88 88 Y88 o.`Y8b 88""   88"Yb    88   
//     88  Yb 88oodP   88   88  Yb 888888 888888 oooooooooo 88 88  Y8 8bodP' 888888 88  Yb   88   

// Fonction pour insérer un point dans un arbre
// is_controle est un booléen, si il est à 1 en cas de point qui existe déjà on n'ajoute pas le point, on compte le nombre de controle des points de controle qui ont été ajouté
// La fonction retourne un ID_TYPE qui est la somme des controles qui ont été ajouté (0 si is_controle==0)
ID_TYPE rbtree_insert(rbtree t, rbtree_node new_node, int is_controle){
    ID_TYPE sum_controle = 0;
    int exists = -1;
    rbtree_node n;
    rbtree_value_elem tmp_elem;
    // On remet à zéro les suivants et précédants du noeud
    new_node->parent = NULL;
    new_node->left = NULL;
    new_node->right = NULL;
    // Comportement spécifique quand l'arbre est vide
    if(t->root==NULL)
    {
        //fprintf(get_debug(), "insertion root/new smallest : %f\n", (double)(new_node->key));
        new_node->color = BLACK;
        t->root = new_node;
        t->smallest = new_node;
    }else{
        // Si la valeur à insérer est plus petite que la smallest, on met smallest à jour et on insère le point au bon endroit
        if(new_node->key < t->smallest->key)
        {
            //fprintf(get_debug(), "insertion new smallest : %f\n", (double)(new_node->key));
            t->smallest->left = new_node;
            new_node->parent = t->smallest;
            t->smallest = new_node;
            exists = 0;
        }else{
            n = t->root;
            while (exists == -1)
            {
                if (new_node->key == n->key)
                {
                    // Ajout de la valeur du nouveau noeud à la fin de la liste des valeurs (si ce n'est pas des points de controle)
                    if(is_controle!=1)
                    {
                        while(new_node->value_first != NULL)
                        {
                            tmp_elem = new_node->value_first;
                            new_node->value_first = tmp_elem->next;
                            tmp_elem->next = n->value_first;
                            n->value_first = tmp_elem;
                            t->size++;
                        }
                    }else{
                        rbtree_free_value_elem(t, new_node, 0);    // on libère la mémoire de la liste de valeur
                    }
                    free(new_node);
                    exists = 1;            // Indique qu'il ne faut pas vérifier l'intégrité de l'arbre ensuite (inutile)
                } else if (new_node->key < n->key) {
                    if (n->left == NULL)
                    {
                        n->left = new_node;
                        new_node->parent = n;
                        exists = 0;        // Indique qu'il faut vérifier l'intégrité de l'arbre ensuite
                    }else{
                        n = n->left;
                    }
                } else {    // new_node->key > n->key
                    if (n->right == NULL)
                    {
                        n->right = new_node;
                        new_node->parent = n;
                        exists = 0;        // Indique qu'il faut vérifier l'intégrité de l'arbre ensuite
                    }else{
                        n = n->right;
                    }
                }
            }
        }        
    }
    // On met à jour l'arbre pour le garder équilibré
    if(exists != 1)
    {
        if(is_controle==1)
            sum_controle = RBTREE_GET_AS_CONTROLE(new_node)->nb_controle;
        insert_case1(t, new_node);
        // On compte le nombre de valeurs
        tmp_elem = new_node->value_first;
        while(tmp_elem != NULL)
        {
            t->size++;
            tmp_elem = tmp_elem->next;
        }
    }
    return sum_controle;
}

//     88""Yb 88""Yb 888888 88""Yb 888888 888888            88 88b 88 .dP"Y8 888888 88""Yb 888888            88b 88 888888 Yb        dP 
//     88__dP 88__dP   88   88__dP 88__   88__              88 88Yb88 `Ybo." 88__   88__dP   88              88Yb88 88__    Yb  db  dP  
//     88"Yb  88""Yb   88   88"Yb  88""   88""              88 88 Y88 o.`Y8b 88""   88"Yb    88              88 Y88 88""     YbdPYbdP   
//     88  Yb 88oodP   88   88  Yb 888888 888888 oooooooooo 88 88  Y8 8bodP' 888888 88  Yb   88   oooooooooo 88  Y8 888888    YP  YP    

// Fonction pour insérer un nouveau point dans l'arbre
// is_controle est un booléen, si il est à 1 en cas de point qui existe déjà on n'ajoute pas le point, on compte le nombre de controle des points de controle qui ont été ajouté
// La fonction retourne un ID_TYPE qui est le nombre de controles qui ont été ajoutés (0 si is_controle==0)
ID_TYPE rbtree_insert_new(rbtree t, DISTANCE_TYPE key, void *value, int is_controle){
    rbtree_node inserted_node = rbtree_new_node(key, value);
    return rbtree_insert(t, inserted_node, is_controle);
}
