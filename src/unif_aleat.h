// ** ================================================================================================================
// ** R-Package: graphscan 1.1
// ** Fichier : src/unif_aleat.h
// ** Description : Génération uniforme de nombre aléatoire
// **               Modification de la fonction R pour la rendre thread-safe.
// ** D'après R (nmath/standalone/sunif.c), il s'agit d'une version de Marsaglia-MultiCarry.
// **
// ** License : GPL-2 | GPL-3
// **
// ** Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
// ** ================================================================================================================
#include <Rconfig.h>


// Structure pour stocker les graînes de la génération aléatoire
typedef struct seed_t
{
    unsigned int seed_1;
    unsigned int seed_2;
}Seed;

// Déclarations
Seed *unif_aleat_creer_seed();
double unif_aleat_generer(Seed *seed);
