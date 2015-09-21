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

#include <time.h>
#include <stdlib.h>
#include "unif_aleat.h"
#include <Rinternals.h>
#ifdef _OPENMP
    #include <omp.h>
#endif



//Fonction qui génère une nouvelle graîne en utilisant le temps et le numéro de thread
Seed *unif_aleat_creer_seed(){
    Seed *seed = (Seed*)malloc(sizeof(struct seed_t));
    if(seed == NULL)
        error("\nERROR: unable to allocate the memory for the variable Seed - terminating\n");
    // Initialisation avec le temps et le numéro de thread (+1 car le premier ou le master est 0)
    #ifdef _OPENMP
        unsigned int base_seed = ((int)time(NULL))^(omp_get_thread_num()+1);
    #else
        unsigned int base_seed = (int)time(NULL);
    #endif
    seed->seed_1 = base_seed+1234;
    seed->seed_2 = base_seed+5678;
    return seed;
}

// Fonction qui génère un nombre aléatoire entre 0(inclu) et 1(exclu) en mettant à jour la graîne passée en paramètre
double unif_aleat_generer(Seed *seed)
{
    seed->seed_1= 36969*(seed->seed_1 & 0177777) + (seed->seed_1>>16);
    seed->seed_2= 18000*(seed->seed_2 & 0177777) + (seed->seed_2>>16);
    return ((seed->seed_1 << 16)^(seed->seed_2 & 0177777)) * 2.328306437080797e-10; /* in [0,1) */
} 
