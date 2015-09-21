// ** ================================================================================================================
// ** R-Package: graphscan 1.1
// ** Fichier : src/detection_dagregat.h
// ** Description : détection 1D des clusters avec l'indice de Cucala et celui de Kulldorff. 
// ** License : GPL-2 | GPL-3
// ** Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
// ** ================================================================================================================
#ifdef _OPENMP
	#include <omp.h>
#endif
#include <time.h>

#include "cucala_methode.h"
#include "fonction_math.h"
#include "unif_aleat.h"

#include <R.h>
#include <Rinternals.h>

/** fonction calcul_agregat_positif_et_indice_cucala()
 ** donne la position d'un agregat positif et l'indice de cucala associee
 ** entree :
 ** nbEv  nombre d'evenement
 ** vecteur_D  pointeur sur le vecteur de distance
 ** sortie :
 ** resultat  une structure de forme agregat_potentiel_indice_cucala
 **/
agregat_potentiel_indice_cucala_t calcul_agregat_positif_et_indice_cucala(
int nbEv, long double * vecteur_D)
{
	int i,k;
	agregat_potentiel_indice_cucala_t resultat;
	long double distance_min;
	long double beta_min = 1;
	long double * vecteur_T = malloc((nbEv-1)*sizeof(long double));
	
 	if( vecteur_T == NULL )   // si l'allocation échoue
 	{
	  error("Allocation error in cucala_methode.c");     //on affiche un message d'erreur   
 	}
	
	long double beta_tmp;
	int i_max_tmp = 0;

	resultat.indice_debut = 0;
	resultat.indice_fin   = 0;
	resultat.indice_cucala= 0;

	// copie du tableau D dans le tableau T en decalant de 1
	// suppression de la premiere et de la derniere valeur
	for(i=0;i<=(nbEv-2);i++)
	{
		vecteur_T[i] = vecteur_D[i+1];
	}
	for(k=1;k<=(nbEv-2);k++)
	{
		distance_min = 1;
		// remplacer le for par un while avec comme condition supplementaire distance_min = 0
		for(i=0;i<=(nbEv-2-k);i++)
		{
			//calcul des distances entre statistiques d'ordre quelconques en fixant k
			vecteur_T[i] = vecteur_T[i] + vecteur_D[i+1+k];
			// recherche de la distance minimale
			if(vecteur_T[i] < distance_min)
			{
				distance_min = vecteur_T[i];
				i_max_tmp = i;
			}
		}
		beta_tmp = Beta((double)distance_min,(double)k,(double)(nbEv+1-k));
		if(beta_tmp <= beta_min)
		{
			// position du premier élément de l'agregat
			resultat.indice_debut = i_max_tmp;
			// position du dernier élément de l'agregat
			resultat.indice_fin = resultat.indice_debut + k;
			beta_min = beta_tmp;
			if(beta_tmp != 0 && beta_tmp > 1/MAX_DOUBLE)
				resultat.indice_cucala = 1/beta_tmp;
			else
				resultat.indice_cucala = MAX_DOUBLE;
		}
	}
	free(vecteur_T);
	return resultat;
}

/** fonction calcul_agregat_negatif_et_indice_cucala()
 ** donne la position d'un agregat negatif et l'indice de cucala associee
 ** entree :
 ** nbEv  nombre d'evenement
 ** vecteur_D  pointeur sur le vecteur de distance
 ** sortie :
 ** resultat  une structure de forme agregat_potentiel_indice_cucala
 **/
agregat_potentiel_indice_cucala_t calcul_agregat_negatif_et_indice_cucala(
int nbEv, long double * vecteur_D)
{
	int i,k;
	agregat_potentiel_indice_cucala_t resultat;
	long double distance_max;
	long double beta_max = 0;
	long double * vecteur_T= malloc((nbEv-1)*sizeof(long double));
	if( vecteur_T == NULL )                        //si l'allocation échoue
	{
		error("Allocation impossible agregat negatif");             //on affiche un message d'erreur      
	}
	long double beta_tmp;
	int i_max_tmp = 0;

	resultat.indice_debut = 0;
	resultat.indice_fin   = 0;
	resultat.indice_cucala= 0;
	
	// copie du tableau D dans le tableau T en decalant de 1
	// suppression de la premiere et de la derniere valeur
	for(i=0;i<=(nbEv-2);i++)
	{
		vecteur_T[i] = vecteur_D[i+1];
	}

	for(k=1;k<=(nbEv-2);k++)
	{
		distance_max = 0;
		for(i=0;i<=(nbEv-2-k);i++)
		{
			//calcul des distances entre statistiques d'ordre quelconques en fixant k
			vecteur_T[i] = vecteur_T[i] + vecteur_D[i+1+k];
			// recherche de la distance maximale
			if(vecteur_T[i] > distance_max)
			{
				distance_max = vecteur_T[i];
				i_max_tmp = i;
			}
		}
		beta_tmp = Beta((double)distance_max,(double)k,(double)(nbEv+1-k));
		if(beta_tmp >= beta_max)
		{
			// position du premier élément de l'agregat
			resultat.indice_debut = i_max_tmp;
			// position du dernier élément de l'agregat
			resultat.indice_fin = resultat.indice_debut + k;
			beta_max = beta_tmp;
			if(beta_tmp != 0 && beta_tmp > 1/MAX_DOUBLE)
				resultat.indice_cucala = 1/beta_tmp;
			else
				resultat.indice_cucala = MAX_DOUBLE;
		}
	}
	free(vecteur_T);
	return resultat;
}

/** fonction calcul_p_valeur_positif()
 ** calcul la p_valeur de l'agregat positif en comparaison avec des simulations
 ** entree :
 ** nbEv  nombre d'evenements
 ** nbsim  nombre de simulation
 ** indice_cucala  indice de cucala associé a l'agregat
 ** sortie :
 ** p_val_max  valeur de la p_valeur
 **/
double calcul_p_valeur_positif(int nbEv, int nbsim, long double indice_cucala)
{
	
	double p_val_max =1.0;

	#pragma omp parallel shared(p_val_max)
	{
		// Initialisation de la génération aléatoire
		Seed *seed = unif_aleat_creer_seed();
		agregat_potentiel_indice_cucala_t agregat_positif;
		int i,j;

		#pragma omp for 
		for(i=0;i<nbsim;i++)
		{
			double * var_sim = malloc((nbEv+2)*sizeof(double));
			if( var_sim == NULL )                        //si l'allocation échoue
			{
				error("Allocation impossible sim positif");             //on affiche un message d'erreur      
			}
			// Tirage aléatoire d'évènements
			for(j=0;j<nbEv+2;++j)
			{
				var_sim[j] = unif_aleat_generer(seed);	// génère entre [0,1)
				if(var_sim[j]==0.0)
					j--;	// Si on génère pile un 0, on refait le tirage (chance très faible d'arriver)
			}

			//on trie les éléments du vecteur en ordre croissant par le tri par tas(heapsort)
			//gsl_heapsort(var_sim, nbEv+2, sizeof(double), (gsl_comparison_fn_t) &compare_doubles);

			// Trie par qsort(plus rapide)
			qsort(var_sim, nbEv+2, sizeof(double), &compare_doubles);

			var_sim[0] = 0;
			var_sim[nbEv+1] = 1;
			
			long double * vecteur_D = malloc((nbEv+1)*sizeof(long double));
			if( vecteur_D == NULL )                        //si l'allocation échoue
			{
				error("Allocation impossible vecteur_D positif");             //on affiche un message d'erreur      
			}
			
			distance_entre_stat_dordre(nbEv,var_sim,vecteur_D);
			
			agregat_positif = calcul_agregat_positif_et_indice_cucala(nbEv,vecteur_D);
			#pragma omp critical
			{
				if(agregat_positif.indice_cucala >= indice_cucala)
					p_val_max++;
			}
			free(vecteur_D);
			free(var_sim);
		}
		free(seed);
	}
	p_val_max = p_val_max/(nbsim+1);
	return p_val_max;
}

/** fonction calcul_p_valeur_negatif()
 ** calcul la p_valeur de l'agregat negatif en comparaison avec des simulations
 ** entree :
 ** nbEv  nombre d'evenements
 ** nbsim  nombre de simulation
 ** indice_cucala  indice de cucala associé a l'agregat
 ** sortie :
 ** p_val_max  valeur de la p_valeur
 **/
double calcul_p_valeur_negatif(int nbEv, int nbsim, long double indice_cucala)
{
	double p_val_min =1.0;
	#pragma omp parallel shared(p_val_min)
	{
		// Initialisation de la génération aléatoire
		Seed *seed = unif_aleat_creer_seed();
		agregat_potentiel_indice_cucala_t agregat_negatif;
		int i,j;

		#pragma omp for 
		for(i=0;i<nbsim;i++)
		{
			double * var_sim = malloc((nbEv+2)*sizeof(double));
			if( var_sim == NULL )                        //si l'allocation échoue
			{
				error("Allocation impossible sim negatif");             //on affiche un message d'erreur      
			}
			// Tirage aléatoire d'évènements
			for(j=0;j<nbEv+2;++j)
			{
				var_sim[j] = unif_aleat_generer(seed);	// génère entre [0,1)
				if(var_sim[j]==0.0)
					j--;	// Si on génère pile un 0, on refait le tirage (chance très faible d'arriver)
			}

			//on trie les éléments du vecteur en ordre croissant par le tri par tas(heapsort)
			//gsl_heapsort(var_sim, nbEv+2, sizeof(double), (gsl_comparison_fn_t) &compare_doubles);

			// Trie par qsort(plus rapide)
			qsort(var_sim, nbEv+2, sizeof(double), &compare_doubles);

			var_sim[0] = 0;
			var_sim[nbEv+1] = 1;
			
			long double * vecteur_D = malloc((nbEv+1)*sizeof(long double));

			if( vecteur_D == NULL )                        //si l'allocation échoue
			{
				error("Allocation impossible vecteur d negatif");             //on affiche un message d'erreur      
			}

			distance_entre_stat_dordre(nbEv,var_sim,vecteur_D);
			
			agregat_negatif = calcul_agregat_negatif_et_indice_cucala(nbEv,vecteur_D);

			#pragma omp critical
			{
				if(agregat_negatif.indice_cucala <= indice_cucala)
					p_val_min = p_val_min +1 ;
			}
			free(vecteur_D);
			free(var_sim);
		}
		free(seed);
	}
	p_val_min = p_val_min/(nbsim+1);
	return p_val_min;
}

/** fonction calcul_p_valeur_negatif_positif()
 ** calcul la p_valeur d'un agregat negatif et d'un agregat positif en comparaison avec les memes simulations
 ** entree :
 ** nbEv  nombre d'evenements
 ** nbsim  nombre de simulation
 ** indice_cucala_pos  indice de cucala associé a l'agregat positif
 ** indice_cucala_neg  indice de cucala associé a l'agregat negatif
 ** entree/sortie :
 ** p_val_max  valeur de la p_valeur pour l'agregat positif
 ** p_val_min  valeur de la p_valeur pour l'agregat negatif
 **/
void calcul_p_valeur_negatif_positif(int nbEv, int nbsim, long double indice_cucala_pos,long double indice_cucala_neg,double * p_val_min, double * p_val_max)
{
	*p_val_max =1.0;
	*p_val_min =1.0;
		
	#pragma omp parallel shared(p_val_max,p_val_min)
	{
		// Initialisation de la génération aléatoire
		Seed *seed = unif_aleat_creer_seed();
		agregat_potentiel_indice_cucala_t agregat_negatif;
		agregat_potentiel_indice_cucala_t agregat_positif;
		int i,j;
	
		#pragma omp for
		for(i=0;i<nbsim;i++)
		{
			long double * vecteur_D = malloc((nbEv+1)*sizeof(long double));
			if( vecteur_D == NULL )                        //si l'allocation échoue
			{
				error("Allocation impossible vecteur_D both");             //on affiche un message d'erreur      
			}
			double * var_sim = malloc((nbEv+2)*sizeof(double));
			if( var_sim == NULL )                        //si l'allocation échoue
			{
				error("Allocation impossible sim both");             //on affiche un message d'erreur      
			}
			// Tirage aléatoire d'évènements
			for(j=0;j<nbEv+2;++j)
			{
				var_sim[j] = unif_aleat_generer(seed);	// génère entre [0,1)
				if(var_sim[j]==0.0)
					j--;	// Si on génère pile un 0, on refait le tirage (chance très faible d'arriver)
			}

			//on trie les éléments du vecteur en ordre croissant par le tri par tas(heapsort)
			//gsl_heapsort(var_sim, nbEv+2, sizeof(double), (gsl_comparison_fn_t) &compare_doubles);

			// Trie par qsort(plus rapide)
			qsort(var_sim, nbEv+2, sizeof(double), &compare_doubles);

			var_sim[0] = 0;
			var_sim[nbEv+1] = 1;

			distance_entre_stat_dordre(nbEv,var_sim,vecteur_D);
			
			agregat_negatif = calcul_agregat_negatif_et_indice_cucala(nbEv,vecteur_D);
			agregat_positif = calcul_agregat_positif_et_indice_cucala(nbEv,vecteur_D);

			#pragma omp critical
			{
				if(agregat_negatif.indice_cucala <= indice_cucala_neg)
					*p_val_min = *p_val_min +1 ;
				if(agregat_positif.indice_cucala >= indice_cucala_pos)
					*p_val_max = *p_val_max +1 ;
			}	
			free(vecteur_D);
			free(var_sim);				
		}	
		free(seed);
	}
	*p_val_min = *p_val_min/(nbsim+1);
	*p_val_max = *p_val_max/(nbsim+1);	
}

/** fonction compare_doubles()
 ** compare deux long double entre eux
 ** entree :
 ** a  premier long double a comparer
 ** b  deuxieme long double a comparer
 **/
// fonction pour comparer 2 éléments
int compare_doubles(const void *a, const void *b)
{
	if(*(double*)a > *(double*)b)
		return 1;
	if (*(double*)a < *(double*)b)
		return -1;
	return 0;	
}

