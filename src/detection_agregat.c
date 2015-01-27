// ** ================================================================================================================
// ** R-Package: graphscan
// ** Fichier : src/detection_dagregat.c
// ** Description : détection 1D des clusters avec l'indice de Cucala et celui de Kulldorff.
// ** License : GPL-2 | GPL-3
// ** Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
// ** ================================================================================================================
#ifdef _OPENMP
	#include <omp.h>
#endif

#include "detection_agregat.h"
#include "cucala_methode.h"
#include "fonction_math.h"
#include "unif_aleat.h"

/** fonction add_resultat()
 ** ajoute un resultat a la liste chainee
 ** entree :
 ** liste  tete de la liste chainee
 ** result  resultat a ajoute
 **/
void add_resultat(list_t * liste, resultat_t * result)
{
	if(liste->taille == 0)
	{
		liste->debut = result;
		liste->fin   = result;
	}
	else
	{
		liste->fin->suiv = result;
		liste->fin = result;
	}
	result->suiv = liste->debut;
	liste->taille = liste->taille +1;
}

/** fonction detection_multiple_dagregat()
 ** detecte les agregat positif et/ou negatif
 ** entree :
 ** nbEv_r  nombre d'evenement
 ** taille_sequence_r  nombre de site
 ** vecteur_X_r  pointeur sur vecteur de position
 ** alpha_r  seuil significatif
 ** theta_r
 ** nbSim_r  nombre de simulation
 ** choix_detection_r  choix de la detection ( positif,negatif,mixte)
 ** choix_type_agregat_r  choix de l'agregat
 ** sortie :
 ** sortie  
 **/
SEXP detection_multiple_dagregat(SEXP nbEv_r,
	SEXP normalisation_debut_r,SEXP normalisation_fin_r,SEXP vecteur_X_r,
	SEXP alpha_r, SEXP theta_r, SEXP nbSim_r,SEXP choix_detection_r, SEXP choix_type_agregat_r)
{
	long double * tab_P ;
	long double * vecteur_D;
	int i;
	

	SEXP sortie;
	int nb_colonne = 0; // nombre de resultat trouvé 
	int nb_ligne = 6; // nombre d'element à renvoyer

	int nbEv = INTEGER_VALUE(nbEv_r);
	int normalisation_debut =INTEGER_VALUE(normalisation_debut_r);
	int normalisation_fin = INTEGER_VALUE(normalisation_fin_r);
	double * vecteur_X = REAL(vecteur_X_r);
	double alpha = NUMERIC_VALUE(alpha_r);
	int theta = INTEGER_VALUE(theta_r);
	int nbSim = INTEGER_VALUE(nbSim_r);
	int choix_detection = INTEGER_VALUE(choix_detection_r);
	int choix_type_agregat = INTEGER_VALUE(choix_type_agregat_r);

	
	resultat_t * result;
	
	list_t *list_head = malloc(sizeof(list_t));
	list_head->debut = NULL;
	list_head->fin = NULL;
	list_head->taille = 0;
	vecteur_D = malloc((nbEv+1)*sizeof(long double));
	
 	if( vecteur_D == NULL )                        //si l'allocation échoue
 	{
	  error("Allocation error for vecteur_D in detection_agregat.c");     //on affiche un message d'erreur   
 	}
	
	tab_P = malloc((nbEv+1)*sizeof(long double));
	
 	if( tab_P == NULL )                        //si l'allocation échoue
 	{
	  error("Allocation error for tab_P in detection_agregat.c");     //on affiche un message d'erreur
 	}
	
	// utiliser le nombre maximum de thread, soit le nombre de coeurs.
	#ifdef _OPENMP
		omp_set_num_threads(omp_get_num_procs());
	#endif
		
	normalisation_et_distance_entre_stat_dordre(normalisation_debut,normalisation_fin,nbEv,vecteur_X,vecteur_D,tab_P);
	switch(choix_detection)
	{
		case 1 : // detection multiple des agregats positif
			{
				long double position1 =0; // position du premier element de l'agregat
				long double position2 =0; // position du dernier element de l'agregat
				double p_val_max; // p valeur de l'agregat
				agregat_potentiel_indice_cucala_t agregat_positif;
				
				bool fin = false;
				
				while(fin == false && nbEv > 1)
				{
					if(nbEv >1)
					{
						//detection d'un agregat positif potentiel et calcul de son indice de cucala
						agregat_positif = calcul_agregat_positif_et_indice_cucala(nbEv,vecteur_D);
						p_val_max = calcul_p_valeur_positif(nbEv,nbSim,agregat_positif.indice_cucala);
						
						position1= tab_P[agregat_positif.indice_debut];
						position2= tab_P[agregat_positif.indice_fin];
						if(p_val_max <= alpha)
						{
							//decalage des tableaux vecteur_D et tab_P
							decalage_tableau(vecteur_D, tab_P,agregat_positif.indice_debut, agregat_positif.indice_fin, &nbEv);

							result = malloc(sizeof(resultat_t));
							result->debut = position1;
							result->fin   = position2;
							result->indice_concentration = agregat_positif.indice_cucala;
							result->p_valeur = p_val_max;
							result-> positif = 1;
							result->id_agregat = list_head->taille;
							add_resultat(list_head,result);
						}
						else
							fin = true;
					}
					else
						fin = true;
				}
				break;
			}
		case 2 : // detection multiple des agregats negatif
			{
				long double position1 =0; // position du premier element de l'agregat
				long double position2 =0; // position du dernier element de l'agregat
				double p_val_min; // p valeur de l'agregat
				agregat_potentiel_indice_cucala_t agregat_negatif;
				
				bool fin = false;
				while(fin == false && nbEv > 1)
				{
					if(nbEv > 1)
					{
						//detection d'un agregat negatif potentiel et calcul de son indice de cucala
						agregat_negatif = calcul_agregat_negatif_et_indice_cucala(nbEv,vecteur_D);

						p_val_min = calcul_p_valeur_negatif(nbEv,nbSim,agregat_negatif.indice_cucala);

						position1= tab_P[agregat_negatif.indice_debut];
						position2= tab_P[agregat_negatif.indice_fin];
						
						if(p_val_min <= alpha)
						{
							//decalage des tableaux vecteur_D et tab_P
							decalage_tableau(vecteur_D, tab_P,agregat_negatif.indice_debut, agregat_negatif.indice_fin, &nbEv);

							result = malloc(sizeof(resultat_t));
							result->debut = position1;
							result->fin   = position2;
							result->indice_concentration = agregat_negatif.indice_cucala;
							result->p_valeur = p_val_min ;
							result-> positif = 0;
							result->id_agregat = list_head->taille;
							add_resultat(list_head,result);
						}
						else
							fin = true;
					}
					else
						fin = true;
				}
				
				break;
			}
		case 3 : // detection multiple des agregats positif et negatif
			{
				long double position1 =0; // position du premier element pour l'agregat positif
				long double position2 =0; // position du dernier element pour l'agregat positif
				double p_val_max; // p valeur pour l'agregat positif
				agregat_potentiel_indice_cucala_t agregat_positif;
				
				long double position3 =0; // position du premier element pour l'agregat negatif
				long double position4 =0; // position du dernier element pour l'agregat negatif
				double p_val_min; // p valeur pour l'agregat negatif
				agregat_potentiel_indice_cucala_t agregat_negatif;
				
				
				bool fin = false;
				
				while(fin == false && nbEv>1)
				{
					if(nbEv >1)
					{
						//detection d'un agregat positif potentiel et calcul de son indice de cucala
						agregat_positif = calcul_agregat_positif_et_indice_cucala(nbEv,vecteur_D);
						//detection d'un agregat negatif potentiel et calcul de son indice de cucala
						agregat_negatif = calcul_agregat_negatif_et_indice_cucala(nbEv,vecteur_D);
						calcul_p_valeur_negatif_positif(nbEv,nbSim,agregat_positif.indice_cucala,agregat_negatif.indice_cucala,&p_val_min, &p_val_max);
						position1= tab_P[agregat_positif.indice_debut];
						position2= tab_P[agregat_positif.indice_fin];
						
						position3= tab_P[agregat_negatif.indice_debut];
						position4= tab_P[agregat_negatif.indice_fin];
						// les deux agregats sont significatifs
						if(p_val_max <= alpha/2 && p_val_min <= alpha/2)
						{
							
							//l'agregat negatif est plus significatif que l'agregat positif
							if(p_val_max > p_val_min)
							{
								//decalage des tableaux vecteur_D et tab_P
								decalage_tableau(vecteur_D, tab_P,agregat_negatif.indice_debut, agregat_negatif.indice_fin, &nbEv);

								result = malloc(sizeof(resultat_t));
								result->debut = position3;
								result->fin   = position4;
								result->indice_concentration = agregat_negatif.indice_cucala;
								result->p_valeur = p_val_min;
								result-> positif = 0;
								result->id_agregat = list_head->taille;
								add_resultat(list_head,result);
							}
							else if(p_val_min > p_val_max)//l'agregat positif est plus significatif que l'agregat negatif
							{
								//decalage des tableaux vecteur_D et tab_P
								decalage_tableau(vecteur_D, tab_P,agregat_positif.indice_debut, agregat_positif.indice_fin, &nbEv);

								result = malloc(sizeof(resultat_t));
								result->debut = position1;
								result->fin   = position2;
								result->indice_concentration = agregat_positif.indice_cucala;
								result->p_valeur = p_val_max;
								result-> positif = 1;
								result->id_agregat = list_head->taille;
								add_resultat(list_head,result);
							}
							else //les deux agregat sont autant significatif l'un que l'autre
							{
								long double difference_taille_agregat;
								double D_pos=0;
								double D_neg = 0;
								for(i=agregat_positif.indice_debut;i<=agregat_positif.indice_fin;i++)
								{
									D_pos = D_pos + (double)vecteur_D[i];
								}
								for(i=agregat_negatif.indice_debut;i<=agregat_negatif.indice_fin;i++)
								{
									D_neg = D_neg + (double)vecteur_D[i];
								}
								difference_taille_agregat = D_neg - D_pos;

								//la taille de l'agregat negatif est bien plus grande que celle de l'agregat positif
								if(difference_taille_agregat > (theta/nbEv))
								{
									//decalage des tableaux vecteur_D et tab_P
									decalage_tableau(vecteur_D, tab_P,agregat_positif.indice_debut, agregat_positif.indice_fin, &nbEv);

									result = malloc(sizeof(resultat_t));
									result->debut = position1;
									result->fin   = position2;
									result->indice_concentration = agregat_positif.indice_cucala;
									result->p_valeur = p_val_min;
									result-> positif = 1;
									result->id_agregat = list_head->taille;
									add_resultat(list_head,result);
								}
								//la taille de l'agregat negatif est bien plus petite que celle de l'agregat positif
								else if((-difference_taille_agregat) > (theta/nbEv))
								{
									//decalage des tableaux vecteur_D et tab_P
									decalage_tableau(vecteur_D, tab_P,agregat_negatif.indice_debut, agregat_negatif.indice_fin, &nbEv);

									result = malloc(sizeof(resultat_t));
									result->debut = position3;
									result->fin   = position4;
									result->indice_concentration = agregat_negatif.indice_cucala;
									result->p_valeur = p_val_max;
									result-> positif = 0;
									result->id_agregat = list_head->taille;
									add_resultat(list_head,result);
								}
								else
								{
									int choix;
									if(choix_type_agregat == 3)
									{
										Seed *seed = unif_aleat_creer_seed();
										choix = (int)(unif_aleat_generer(seed)*2.0); // choix aleatoire entre entier 0 et 1
										free(seed);
									}
									else
										choix = choix_type_agregat;

									switch(choix)
									{
										case 1 : //choix de l'agregat negatif
												//decalage des tableaux vecteur_D et tab_P
												decalage_tableau(vecteur_D, tab_P,agregat_negatif.indice_debut, agregat_negatif.indice_fin, &nbEv);

												result = malloc(sizeof(resultat_t));
												result->debut = position3;
												result->fin   = position4;
												result->indice_concentration = agregat_negatif.indice_cucala;
												result->p_valeur = p_val_min;
												result-> positif = 0;
												result->id_agregat = list_head->taille;
												add_resultat(list_head,result);
												break;
										case 2 : // choix de l'agregat positif
												//decalage des tableaux vecteur_D et tab_P
												decalage_tableau(vecteur_D, tab_P,agregat_positif.indice_debut, agregat_positif.indice_fin, &nbEv);

												result = malloc(sizeof(resultat_t));
												result->debut = position1;
												result->fin   = position2;
												result->indice_concentration = agregat_positif.indice_cucala;
												result->p_valeur = p_val_max;
												result-> positif = 1;
												result->id_agregat = list_head->taille;
												add_resultat(list_head,result);
												break;
										default :
												break;
									}
								}
							}
						}
						else if(p_val_min <= (alpha/2))
						{
							//decalage des tableaux vecteur_D et tab_P
							decalage_tableau(vecteur_D, tab_P,agregat_negatif.indice_debut, agregat_negatif.indice_fin, &nbEv);

							result = malloc(sizeof(resultat_t));
							result->debut = position3;
							result->fin   = position4;
							result->indice_concentration = agregat_negatif.indice_cucala;
							result->p_valeur = p_val_min;
							result-> positif = 0;
							result->id_agregat = list_head->taille;
							add_resultat(list_head,result);
						}
						else if(p_val_max <= (alpha/2))
						{
							//decalage des tableaux vecteur_D et tab_P
							decalage_tableau(vecteur_D, tab_P,agregat_positif.indice_debut, agregat_positif.indice_fin, &nbEv);

							result = malloc(sizeof(resultat_t));
							result->debut = position1;
							result->fin   = position2;
							result->indice_concentration = agregat_positif.indice_cucala;
							result->p_valeur = p_val_max;
							result-> positif = 1;
							result->id_agregat = list_head->taille;
							add_resultat(list_head,result);
						}
						else
							fin = true;
					}
					else
						fin = true;
				}
				break;
			}
	}
	nb_colonne = list_head->taille;
	//renvoie une matrice de tout les resultat
	PROTECT(sortie = allocMatrix(REALSXP,nb_ligne,nb_colonne));

	i = 0;
	result = list_head->debut;
	while (i<nb_ligne*nb_colonne)
	{

		REAL(sortie)[i] = (double)result->debut;
		REAL(sortie)[i+1] = (double)result->fin;
		REAL(sortie)[i+2] = (double)result->indice_concentration;
		REAL(sortie)[i+3] = (result->p_valeur);
		REAL(sortie)[i+4] = result->positif;
		REAL(sortie)[i+5] = result->id_agregat;

		list_head->debut = result->suiv;
		list_head->taille --;
		free(result);
		result = list_head->debut;
		i =i +6;
	}
	free(tab_P);
	free(vecteur_D);
	UNPROTECT(1);
	return(sortie);
}

/** fonction decalage_tableau()
 ** enleve les distances, les positions de l'agregat trouve et normalise les distances
 ** entree :
 ** vecteur_D  pointeur sur le vecteur de distance
 ** tab_P  pointeur sur le vecteur de position 
 ** indice_debut  position de debut des distances a enlevee
 ** indice_fin  position de fin des distances a enlevee
 ** entree/sortie :
 ** nbEv  nombre d'evenement
 **/
void decalage_tableau(long double * vecteur_D, long double * tab_P,int indice_debut, int indice_fin, int *nbEv)
{
	int i;
	long double somme = 0.0;
	//decalage des tableaux vecteur_D et tab_P	
	for(i=0; i <=(*nbEv -indice_fin);i++)
	{
		vecteur_D[indice_debut+i+1] = vecteur_D[indice_fin+i+1];
		tab_P[indice_debut+i+1] = tab_P[indice_fin+i+1];
	}

	//calcul de la somme des distances D
	for(i=0;i <(*nbEv -indice_fin +indice_debut);i++)
		somme = somme + vecteur_D[i];

	//normalisation des valeurs de vecteur_D
	#pragma omp parallel for
	for(i=0;i <(*nbEv -indice_fin +indice_debut);i++)
		vecteur_D[i] = vecteur_D[i]/somme;

	//calcul du nombre d'evenement restants
	*nbEv =*nbEv + indice_debut - indice_fin ;
}
