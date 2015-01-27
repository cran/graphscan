// ** ================================================================================================================
// ** R-Package: graphscan
// ** Fichier : src/2D_detection_cluster.c
// ** Description : détection 2D et 3D des clusters avec l'indice de Cucala et celui de Kulldorff.
// ** License : GPL-2 | GPL-3
// ** Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
// ** ================================================================================================================

#include "2D_detection_cluster.h"
#include <R_ext/Utils.h> // to allow user interruption in the R interface
#include <Rconfig.h>
#ifdef _OPENMP
	#include <omp.h>
#endif
#include "unif_aleat.h"
#include <time.h>

// L'algorithme est optimisé pour avoir une complexité de O(n log2(n)) en moyenne. 

// 	d8888b. d88888b d888888b d88888b  .o88b. d888888b d888888b  .d88b.  d8b   db          .o88b. db      db    db .d8888. d888888b d88888b d8888b. 
// 	88  `8D 88'     `~~88~~' 88'     d8P  Y8 `~~88~~'   `88'   .8P  Y8. 888o  88         d8P  Y8 88      88    88 88'  YP `~~88~~' 88'     88  `8D 
// 	88   88 88ooooo    88    88ooooo 8P         88       88    88    88 88V8o 88         8P      88      88    88 `8bo.      88    88ooooo 88oobY' 
// 	88   88 88~~~~~    88    88~~~~~ 8b         88       88    88    88 88 V8o88         8b      88      88    88   `Y8b.    88    88~~~~~ 88`8b   
// 	88  .8D 88.        88    88.     Y8b  d8    88      .88.   `8b  d8' 88  V888         Y8b  d8 88booo. 88b  d88 db   8D    88    88.     88 `88. 
// 	Y8888D' Y88888P    YP    Y88888P  `Y88P'    YP    Y888888P  `Y88P'  VP   V8P C88888D  `Y88P' Y88888P ~Y8888P' `8888Y'    YP    Y88888P 88   YD 
// 	                                                                                                                                               
// 	                            
SEXP detection_cluster(SEXP nb_point_r,
	SEXP nb_dimension_r,SEXP nb_simulation_r,
	SEXP id_r, SEXP coordonnees_r, SEXP controle_r,SEXP cas_r, SEXP memory_size_r)
{
	// nombre de thread utilises pour paralleliser 

	// recuperer les donnees de R
	ID_TYPE nb_point = (ID_TYPE)INTEGER_VALUE(nb_point_r);
	ID_TYPE nb_dimension = (ID_TYPE)INTEGER_VALUE(nb_dimension_r);
	ID_TYPE nb_simulation = (ID_TYPE)INTEGER_VALUE(nb_simulation_r);
	ID_TYPE * id = (ID_TYPE*)INTEGER(id_r);
	#ifdef _OPENMP
		ID_TYPE memory_size = (ID_TYPE)INTEGER_VALUE(memory_size_r);
	#endif
	double * coordonnees = REAL(coordonnees_r);
	double * controle = REAL(controle_r);
	double * cas = REAL(cas_r);

	// booléen pour savoir si des valeurs infinies ont étées trouvées lors du calcul
	int valeurs_infinies=0;

	// variables pour mesurer le temps d'exécution du programme
	time_t time_begin, time_end;

	// declarations des variables
	SEXP sortie;

	// variables temporaires utilisées dans les boucles
	ID_TYPE i, j;

	// Pour la création du résultat
	ID_TYPE nb_ligne, nb_colonne;

	ID_TYPE tmp_cas,tmp_controle;

	ID_TYPE nb_cas;
	ID_TYPE nb_controle;
	ID_TYPE nb_point_cas;
	ID_TYPE nb_point_controle;

	Point * point_tmp;

	// affichage de la progression
	Rprintf("in progress ... \n");	

	// temps de calcul
	time_begin = time(NULL);

	// Le tableau points_tab va contenir tous les points (la partie qui ne change pas selon les simulation) sans distinction (tout est initialisé à zéro)
	Point *points_tab = calloc(nb_point,sizeof(Point));
	VERIFY_ALLOC(points_tab);

	// Création du tableaux de points points_cas et points_controle :
	// On commence par créer un tableau de Point_element de la taille de tous les points
	Point_element *points_element_tous = calloc(nb_point,sizeof(Point_element));
	VERIFY_ALLOC(points_element_tous);
	// Ensuite on crée les pointeurs vers les deux tableaux, qui vont se remplir en miroir, le premier en partant du début et le second en partant de la fin
	Point_element *points_cas = points_element_tous;
	Point_element *points_controle = points_element_tous + nb_point;	// arithmétique de pointeur (ici points_controle pointe vers l'élément juste après le dernier élément du tableau,
	
	nb_cas = 0;
	nb_controle = 0;
	nb_point_cas = 0;
	nb_point_controle = 0;

	// 	88""Yb 888888 8b    d8 88""Yb 88     88 .dP"Y8 .dP"Y8    db     dP""b8 888888 
	// 	88__dP 88__   88b  d88 88__dP 88     88 `Ybo." `Ybo."   dPYb   dP   `" 88__   
	// 	88"Yb  88""   88YbdP88 88"""  88  .o 88 o.`Y8b o.`Y8b  dP__Yb  Yb  "88 88""   
	// 	88  Yb 888888 88 YY 88 88     88ood8 88 8bodP' 8bodP' dP""""Yb  YboodP 888888 
	// Remplissage des tableaux de points et d'index
	for(i = 0; i<nb_point; i++)
	{
		// Création du tableau des coordonnées en fonction de la dimension
		points_tab[i].coordonnees = malloc(nb_dimension*sizeof(double));
		VERIFY_ALLOC(points_tab[i].coordonnees);
		// Remplissage des coordonnées
		for(j=0;j< nb_dimension;j++)
		{
			points_tab[i].coordonnees[j] = coordonnees[i*nb_dimension+j];

		}
		
		// Remplissage des données du point
		tmp_cas = (ID_TYPE)cas[i];
		tmp_controle = (ID_TYPE)controle[i];
		points_tab[i].id = (DISTANCE_TYPE)id[i];
		points_tab[i].nb_total = tmp_cas + tmp_controle;

		nb_cas += tmp_cas;
		nb_controle += tmp_controle;

		// Ajout du point dans la liste correspondante
		if(tmp_cas != 0)
		{
			points_cas[nb_point_cas].point = &points_tab[i];	
			points_cas[nb_point_cas].nb_cas = tmp_cas;		// remplissage des données qui vont être utilisées lors du calcul des cluster
			points_cas[nb_point_cas].nb_controle = tmp_controle;
			nb_point_cas++;
		}else{
			// On décale le pointeur le pointeur du tableau points_controle de 1 pour qu'il pointe toujours sur le début du tableau (qui va grandir vers la gauche)
			points_controle--;		// arithmétique des pointeurs
			// Le nouveau point est donc toujours le point 0 (le plus à gauche du tableau)
			points_controle[0].point = &points_tab[i];	
			points_controle[0].nb_controle = tmp_controle;		// remplissage des données qui vont être utilisées lors du calcul des cluster (pas de cas ici)
			// On initialise les arbres qui serviront lors du parcours des arbres kd (que pour les points de controle)
			points_controle[0].controle_cluster = rbtree_create();
			points_controle[0].controle_cluster_sous_arbre = rbtree_create();
			nb_point_controle++;
		}
	}	
	
	// initialisation des listes cluster_cucala et cluster_kulldorff
	// contenants les points format les clusters
	List_cluster_cucala *cluster_cucala = calloc(1,sizeof(List_cluster_cucala));
	VERIFY_ALLOC(cluster_cucala);
	cluster_cucala->concentration_cucala = MIN_DISTANCE;
	List_cluster_kulldorff *cluster_kulldorff = calloc(1,sizeof(List_cluster_kulldorff));
	VERIFY_ALLOC(cluster_kulldorff);
	cluster_kulldorff->concentration_kulldorff = MIN_DISTANCE;
	

	// calcul des indices de concentration de cucala et de kulldorf
	valeurs_infinies += calcul_concentration(points_cas, points_controle, cluster_cucala, cluster_kulldorff, nb_cas, nb_controle, nb_point_cas ,nb_point_controle ,nb_dimension, 0);

	// Affichage

	// Cucala
	Rprintf("\nBest Cucala index: %Lf\n", cluster_cucala->concentration_cucala);
	Rprintf("%d cases and %d control.\n\n", cluster_cucala->nb_cas_cucala, cluster_cucala->nb_controle_cucala);
	
	// Kulldorff
	Rprintf("Best Kulldorff index: %Lf\n", cluster_kulldorff->concentration_kulldorff);
	Rprintf("%d cases and %d control.\n\n", cluster_kulldorff->nb_cas_kulldorff, cluster_kulldorff->nb_controle_kulldorff);


	// 	.d8888. d888888b .88b  d88. db    db db       .d8b.  d888888b d888888b  .d88b.  d8b   db 
	// 	88'  YP   `88'   88'YbdP`88 88    88 88      d8' `8b `~~88~~'   `88'   .8P  Y8. 888o  88 
	// 	`8bo.      88    88  88  88 88    88 88      88ooo88    88       88    88    88 88V8o 88 
	// 	  `Y8b.    88    88  88  88 88    88 88      88~~~88    88       88    88    88 88 V8o88 
	// 	db   8D   .88.   88  88  88 88b  d88 88booo. 88   88    88      .88.   `8b  d8' 88  V888 
	// 	`8888Y' Y888888P YP  YP  YP ~Y8888P' Y88888P YP   YP    YP    Y888888P  `Y88P'  VP   V8P 
	// 	    

	ID_TYPE nb_mesures_total = nb_controle + nb_cas;
	
#ifdef _OPENMP
	// Calcul du nombre total de mesures prélevées, controles et cas confondus
	int nb_simulations_parallelles;
	// Si le nombre de points/tirage aléatoire est faible, inutile de paralléliser, cela ferait au contraire perdre du temps au programme       
	// Les limites choisies ici sont arbitraires                  
	if(nb_point<10000 && nb_mesures_total<1000000)
	{
		nb_simulations_parallelles = 1;
		Rprintf("No parallelization: not enough points.\n");
	}else{
		// Calcul de la taille mémoire (max) nécessaire à l'exécution d'une simulation (estimation)
		size_t taille_mesures = nb_mesures_total*sizeof(ID_TYPE);	// Tableau des mesures (pour tirage aléatoire)
		size_t taille_kd = nb_point*sizeof(Point_element);	// Tableau contenant les points/arbre kd (les 2 arbres en même temps)
		size_t taille_distances = nb_point_cas*2*sizeof(Distance); // *2 -> cas-cas et cas-ppvcontrole. Attention, le nombre de cas dans les simulations peut varier, prévoir une marge d'erreur
		size_t taille_cluster = (nb_point_cas/2)*sizeof(Cluster); // n/2 cluster max (en même temps)
		size_t taille_tree_distance_cluster = (nb_point_cas/2)*sizeof(struct rbtree_t) + nb_point_cas*(sizeof(struct rbtree_node_t)+sizeof(struct rbtree_value_elem_t)); // Pour conserver les distances cas-controle dans les clusters
		size_t taille_tree_controle_cluster = (nb_point_cas/2)*sizeof(struct rbtree_t) + 2*nb_point_controle*(sizeof(struct rbtree_node_t)+sizeof(struct rbtree_value_elem_t)); // Pour conserver les controles dans les clusters. Attention: un controle peut être présent dans plusieurs clusters (1.5 pour gérer mais arbitraire) 
		size_t taille_tree_controle = (nb_point_controle*2)*sizeof(struct rbtree_t) + 2*nb_point_controle*2*(sizeof(struct rbtree_node_t)+sizeof(struct rbtree_value_elem_t)); // arbre dans les points controle. Attention: un controle peut être présent dans plusieurs clusters (1.5 pour gérer mais arbitraire)
		
		// Taille mémoire globale, on rajoute 10mo comme marge d'erreur/variables locales
		size_t taille_simulation = taille_mesures + taille_kd + taille_distances + taille_cluster + taille_tree_distance_cluster + taille_tree_controle_cluster + taille_tree_controle + 10000000;

		// Affichage de la taille mémoire à l'utilisateur
		Rprintf("Memory size estimated by simulation: ");
		Rprintf("%u mega-bytes\n", (unsigned int)taille_simulation/1000000);

		// Récupération du nombre de processeurs de la machine
		int nb_proc = omp_get_num_procs();

		// Calcul du nombre de simulations possibles vu la quantité de mémoire à utiliser
		nb_simulations_parallelles = (int)floor(memory_size*1000000/(double)(taille_simulation)); // memory_size est en mo

		if(nb_simulations_parallelles<=0)
			error("\nERROR: unable to do the simulations: not enough memory (try to increase 'memory_size' parameter, actually to %d) - terminating\n", memory_size);
		if(nb_simulations_parallelles>nb_proc)
			nb_simulations_parallelles = nb_proc;
	}
	Rprintf("Maximum number of threads to do the simulation: %d , using maximum %d mb.\n\n", nb_simulations_parallelles, memory_size);

	// On définit le nombre de threads parallèles à exécuter
	omp_set_num_threads(nb_simulations_parallelles);
#else
	Rprintf("Parallelization not available on this machine, openMP required.\n\n");
#endif	
	double p_valeur_kulldorff = 1.0;
	double p_valeur_cucala = 1.0;
	// Variable pour pouvoir afficher correctement l'avancement du programme
	ID_TYPE suivi_simulation = 0;

	// 	88""Yb    db    88""Yb    db    88     88     888888 88     88 .dP"Y8    db    888888 88  dP"Yb  88b 88 
	// 	88__dP   dPYb   88__dP   dPYb   88     88     88__   88     88 `Ybo."   dPYb     88   88 dP   Yb 88Yb88 
	// 	88"""   dP__Yb  88"Yb   dP__Yb  88  .o 88  .o 88""   88  .o 88 o.`Y8b  dP__Yb    88   88 Yb   dP 88 Y88 
	// 	88     dP""""Yb 88  Yb dP""""Yb 88ood8 88ood8 888888 88ood8 88 8bodP' dP""""Yb   88   88  YbodP  88  Y8 

	#pragma omp parallel private(i,j,nb_point_cas,nb_point_controle) shared(p_valeur_kulldorff,p_valeur_cucala,suivi_simulation,valeurs_infinies)
	{
		ID_TYPE nb_mesures_point, nb_mesures_parcourues, tirage, nb_a_tirer, tmp, tmp_nb_cas_sim;
		Point_element *tmp_point_element;
		List_cluster_cucala *cluster_cucala_sim;
		List_cluster_kulldorff *cluster_kulldorff_sim;
		int tmp_valeurs_infinies=0, valeur_majoritaire, valeur_minoritaire;
		// Création du tableaux de points points_cas_sim et points_controle_sim :
		// On commence par créer un tableau de Point_element de la taille de tous les points
		Point_element *points_element_tous_sim = calloc(nb_point,sizeof(Point_element));
		VERIFY_ALLOC(points_element_tous_sim);
		// Ensuite on crée les pointeurs vers les deux tableaux, qui vont se remplir en miroir, le premier en partant du début et le second en partant de la fin
		Point_element *points_cas_sim = points_element_tous_sim;
		Point_element *points_controle_sim = points_element_tous_sim + nb_point;	// arithmétique de pointeur (ici points_controle_sim pointe vers l'élément juste après le dernier élément du tableau,

		// Initialisation des seeds de la fonction de génération de nombre aléatoire avec la date (en seconde depuis 1970)
		// omp_get_thread_num est utilisé pour être sûr que chaque thread a bien une initialisation du générateur aléatoire différente
		Seed *seed = unif_aleat_creer_seed();		

		// vecteur de taille nb_mesures_total contenant les mesures, et mises à 0 si la mesure est un controle et à 1 si c'est un cas
		// le rapport controle/cas est identique par rapport aux données observées
		ID_TYPE *mesures = malloc(nb_mesures_total*sizeof(ID_TYPE));
		VERIFY_ALLOC(mesures);
		
		// 	888888  dP"Yb  88""Yb 
		// 	88__   dP   Yb 88__dP 
		// 	88""   Yb   dP 88"Yb  
		// 	88      YbodP  88  Yb 

		#pragma omp for
		for(ID_TYPE simulation=0; simulation<nb_simulation; simulation++)
		{				
			//--- simulation des positions des cas ---
			
			// tirage aleatoire de nb_cas positions sur le vecteur position_cas
			// RAND_MAX = 2 147 483 647		-> dépends des machines

			// On va utiliser au choix deux méthodes:
			// - mettre aléatoirement des cases à 1, le reste à zéro
			// - réaliser un tirage aléatoire sans remise, si le nombre de mesures est supérieur à RAND_MAX (et donc la première méthode serait biaisée),
			//	 ou si le nombre de cas est important vis-à-vis du nombre de mesures (>15%) et donc le nombre de chance, lors de tirage aléatoire, de
			//   tomber sur une valeur qui a déjà été tirée devient trop important

			if(nb_mesures_total > RAND_MAX)
			{
				// Tirage aléatoire sans remise au sein des mesures :
				// On commence par définir le vecteur de mesures avec le bon nombre de cas et de controles
				for(i=0; i<nb_cas; i++)
					mesures[i]=1;
				for(; i<nb_mesures_total; i++)
					mesures[i]=0;
				// Puis on tire nos valeurs parmi notre tableau de non tirés, et on met les valeurs ainsi tirées au début du tableau (réduisant la taille du tableau des non tirés)
				for(i=0; i<nb_mesures_total-1; i++)	// -1 car pas besoin de faire un tirage s'il n'y a qu'une seule valeur
				{
					nb_a_tirer = nb_mesures_total-i;
					// Tirage aléatoire au sein des valeurs n'ayant pas encore été tirées (qui sont à la fin du tableau)
					tirage = (ID_TYPE)(unif_aleat_generer(seed) * nb_a_tirer) + i;		// +i car les valeurs à tirer sont à la fin du tableau
					// On échange la valeur tirée avec celle à ajouter à la fin de la liste des tirées (à la position i)
					tmp = mesures[i];
					mesures[i] = mesures[tirage];
					mesures[tirage] = tmp;
				}	
			}else{
				if(nb_cas > nb_controle)
				{
					// On va tirer des contrôles parmi des cas
					valeur_majoritaire = 1;
					valeur_minoritaire = 0;
					nb_a_tirer = nb_controle;
				}else{
					// On va tirer des cas parmi des contrôles
					valeur_majoritaire = 0;
					valeur_minoritaire = 1;
					nb_a_tirer = nb_cas;
				}
				// on met toutes les valeurs à valeur_majoritaire
				for(i=0; i<nb_mesures_total; i++)
					mesures[i] = valeur_majoritaire;
				// on tire aléatoirement dans le tableau des valeurs d'indice et on met la case correspondante à valeur_minoritaire si elle ne l'était pas déjà
				for(i=0; i<nb_a_tirer; )
				{
					tmp = (ID_TYPE)(unif_aleat_generer(seed)*nb_mesures_total); // tmp [0,nb_mesures_total)
					if(mesures[tmp] == valeur_majoritaire)
					{
						mesures[tmp] = valeur_minoritaire;
						i++;
					}
				}
			}

			// On réinitialise le tableau de Point_element
			memset(points_element_tous_sim, 0, nb_point*sizeof(Point_element));

			// On réinitialise le pointeur vers la seconde partie du tableau
			points_controle_sim = points_element_tous_sim + nb_point;

			nb_point_cas = 0;
			nb_point_controle = 0;

			nb_mesures_parcourues = 0;

			// On parcours tous les points
			for(i=0; i<nb_point; i++)
			{
				// nombre de mesures à prendre pour le point
				nb_mesures_point = points_tab[i].nb_total;

				// nombre de cas simulés égale à zéro
				tmp_nb_cas_sim = 0;
				
				// on parcourt une partie du vecteur mesures qui correspond au point i
				for(j=nb_mesures_parcourues; j<nb_mesures_parcourues+nb_mesures_point; j++)
				{
					// on incrémente le nombre de cas simulés
					tmp_nb_cas_sim += mesures[j];
				}

				// Si le nombre de cas est positif, on l'ajoute au tableau des cas
				if(tmp_nb_cas_sim>0)
				{
					tmp_point_element = &points_cas_sim[nb_point_cas];
					nb_point_cas++;
				}else{
					// On décale le pointeur le pointeur du tableau points_controle_sim de 1 pour qu'il pointe toujours sur le début du tableau (qui va grandir vers la gauche)
					points_controle_sim--;		// arithmétique des pointeurs
					// Le nouveau point est donc toujours le point 0 (le plus à gauche du tableau)
					tmp_point_element = &points_controle_sim[0];	
					// On initialise les arbres qui serviront lors du parcours des arbres kd (que pour les points de controle)
					tmp_point_element->controle_cluster = rbtree_create();
					tmp_point_element->controle_cluster_sous_arbre = rbtree_create();
					nb_point_controle++;
				}
				tmp_point_element->point = &points_tab[i];

				tmp_point_element->nb_cas = tmp_nb_cas_sim;
				
				
				// le nombre de controles simulés est égale au nombre total d'individu moins le nombre de cas
				tmp_point_element->nb_controle = points_tab[i].nb_total - tmp_nb_cas_sim;
				
				// on avance sur le vecteur mesures
				nb_mesures_parcourues += nb_mesures_point;	
			}

			// initialiser les listes cluster_cucala_sim et cluster_kulldorff_sim pour stocker les resultats
			cluster_cucala_sim = calloc(1,sizeof(List_cluster_cucala));
			VERIFY_ALLOC(cluster_cucala_sim);
			cluster_cucala_sim->concentration_cucala = MIN_DISTANCE;
			cluster_kulldorff_sim = calloc(1,sizeof(List_cluster_kulldorff));
			VERIFY_ALLOC(cluster_kulldorff_sim);
			cluster_kulldorff_sim->concentration_kulldorff = MIN_DISTANCE;

			// calcul des indices de concentration de cucala et de kulldorf
			tmp_valeurs_infinies += calcul_concentration(points_cas_sim, points_controle_sim, cluster_cucala_sim, cluster_kulldorff_sim, nb_cas, nb_controle, nb_point_cas, nb_point_controle, nb_dimension, 1);

			#pragma omp critical
			{
				R_CheckUserInterrupt(); // To allow user interruption in the R interface
				// affichage de la progression
				Rprintf("simulation %d/%d ->", ++suivi_simulation, nb_simulation);
				// affichage des concentrations
				Rprintf(" Kulldorff: %Lf , Cucala: %Lf\n", cluster_kulldorff_sim->concentration_kulldorff, cluster_cucala_sim->concentration_cucala);

				if( cluster_kulldorff->concentration_kulldorff < cluster_kulldorff_sim->concentration_kulldorff)
					p_valeur_kulldorff++;
					
				if( cluster_cucala->concentration_cucala < cluster_cucala_sim->concentration_cucala)
					p_valeur_cucala++;
			}
			// Libération de la mémoire des listes de cluster cucala et kulldorf
			free(cluster_cucala_sim);
			free(cluster_kulldorff_sim);
		}
		free(points_element_tous_sim);
		free(mesures);
		free(seed);
		#pragma omp critical
		{
			valeurs_infinies += tmp_valeurs_infinies;
		}
	} // Fin de la parallélisation
	p_valeur_cucala = p_valeur_cucala / (nb_simulation+1);
	p_valeur_kulldorff = p_valeur_kulldorff / (nb_simulation+1);

	// temps de calcul
	time_end = time(NULL);

	// 	88""Yb 888888 .dP"Y8 88   88 88     888888    db    888888 
	// 	88__dP 88__   `Ybo." 88   88 88       88     dPYb     88   
	// 	88"Yb  88""   o.`Y8b Y8   8P 88  .o   88    dP__Yb    88   
	// 	88  Yb 888888 8bodP' `YbodP' 88ood8   88   dP""""Yb   88   
	
	// Affichage
	Rprintf("\nCucala pvalue: %f\n", p_valeur_cucala);
	Rprintf("\nKulldorff pvalue: %f\n\n", p_valeur_kulldorff);

	double temp; 
	double hours=0, min=0, sec=0; 
	temp = modf(difftime(time_end, time_begin)/3600., &hours); 
	temp = modf(temp*60., &min); 
	temp = modf(temp*60., &sec); 
	Rprintf("\nCalculation time: ");
	if(hours>0)
		Rprintf("%d h ", (int)hours);
	if(min>0)
		Rprintf("%d min ", (int)min);
	Rprintf("%d sec\n", (int)sec);

	if(valeurs_infinies>0)
		Rprintf("\nWarning: infinite found in the Cucala index, result may not be precise. Try to do the clustering with less points or use the Kulldorff index.\n");

	Rprintf("\n");

	//creation du resultat a retourner
	nb_ligne = 0; 
	nb_colonne = 2 + 2 * nb_dimension; // une pour cucala , une pour kulldorf + pour chaque dimension

	if(cluster_cucala->taille >= cluster_kulldorff->taille)
		nb_ligne = cluster_cucala->taille + 6;
	else
		nb_ligne = cluster_kulldorff->taille + 6;


	PROTECT(sortie = allocMatrix(REALSXP,nb_ligne,nb_colonne));
	/* |concentration cucala	|		|		|concentration kulldorf	|		|		|
	 * |distance				|		|		|distance				|		|		|
	 * |taille					|		|		|taille					|		|		|
	 * |pvaleur					|		|		|pvaleur				|		|		|
	 * |nb_control				|		|		|nb_controle			|		|		|
	 * |nb_cas					|		|		|nb_cas					|		|		|
	 * |id1						|coord_x|coord_y|id1					|coord_x|coord_y|
	 * |id2						|coord_x|coord_y|id2					|coord_x|coord_y|
	 * |id3						|coord_x|coord_y|id3					|coord_x|coord_y|
	 * |						|coord_x|coord_y|id4					|coord_x|coord_y|*/
	//en R on remplie colonne par colonne
	i = 0;
	point_tmp = cluster_cucala->debut;
	REAL(sortie)[0] = (double)cluster_cucala->concentration_cucala;
	REAL(sortie)[1] = sqrt((double)cluster_cucala->dist_max);			// sqrt car on travaille sur les distances au carré, pour économiser une opération
	REAL(sortie)[2] = (double)cluster_cucala->taille;
	REAL(sortie)[3] = (double)p_valeur_cucala;
	REAL(sortie)[4] = (double)cluster_cucala->nb_controle_cucala;
	REAL(sortie)[5] = (double)cluster_cucala->nb_cas_cucala;

	ID_TYPE dim;

	while(i < cluster_cucala->taille)
	{
		REAL(sortie)[i+6] = (double)point_tmp->id;
		for(dim=0;dim<nb_dimension;dim++)
		{
			REAL(sortie)[i+6+(dim+1)*nb_ligne] = (double)point_tmp->coordonnees[dim];
		}
		point_tmp = point_tmp->suiv_cluster_cucala;
		i++;
	}
	point_tmp = cluster_kulldorff->debut;
	//j = 3*nb_ligne;
	j = (nb_dimension+1)*nb_ligne;
	REAL(sortie)[j+0] = (double)cluster_kulldorff->concentration_kulldorff;
	REAL(sortie)[j+1] = sqrt((double)cluster_kulldorff->dist_max);		// sqrt car on travaille sur les distances au carré, pour économiser une opération
	REAL(sortie)[j+2] = (double)cluster_kulldorff->taille;
	REAL(sortie)[j+3] = (double)p_valeur_kulldorff;
	REAL(sortie)[j+4] = (double)cluster_kulldorff->nb_controle_kulldorff;
	REAL(sortie)[j+5] = (double)cluster_kulldorff->nb_cas_kulldorff;
	i = 0;
	while(i < cluster_kulldorff->taille)
	{
		REAL(sortie)[j+i+6] = (double)point_tmp->id;
		for(dim=0;dim<nb_dimension;dim++)
		{
			REAL(sortie)[j+i+6+(dim+1)*nb_ligne] = (double)point_tmp->coordonnees[dim];
		}
		point_tmp = point_tmp->suiv_cluster_kulldorff;
		i++;
	}

	UNPROTECT(1);
	// 	888888 88""Yb 888888 888888 
	// 	88__   88__dP 88__   88__   
	// 	88""   88"Yb  88""   88""   
	// 	88     88  Yb 888888 888888 

	// On libère les clusters
	free(cluster_cucala);
	free(cluster_kulldorff);

	// On libère les points
	for(i = 0; i<nb_point; i++)
		free(points_tab[i].coordonnees);
	free(points_tab);

	// On libère les point_element
	// Pas besoin de libérer les arbres pour les points de controle, ils le sont automatiquement
	// à la fin de la fonction calcul_concentration
	free(points_element_tous);

	return(sortie);
}






                                                                                  

// 	 .o88b.  .d8b.  db       .o88b. db    db db               .o88b.  .d88b.  d8b   db  .o88b. d88888b d8b   db d888888b d8888b.  .d8b.  d888888b d888888b  .d88b.  d8b   db 
// 	d8P  Y8 d8' `8b 88      d8P  Y8 88    88 88              d8P  Y8 .8P  Y8. 888o  88 d8P  Y8 88'     888o  88 `~~88~~' 88  `8D d8' `8b `~~88~~'   `88'   .8P  Y8. 888o  88 
// 	8P      88ooo88 88      8P      88    88 88              8P      88    88 88V8o 88 8P      88ooooo 88V8o 88    88    88oobY' 88ooo88    88       88    88    88 88V8o 88 
// 	8b      88~~~88 88      8b      88    88 88              8b      88    88 88 V8o88 8b      88~~~~~ 88 V8o88    88    88`8b   88~~~88    88       88    88    88 88 V8o88 
// 	Y8b  d8 88   88 88booo. Y8b  d8 88b  d88 88booo.         Y8b  d8 `8b  d8' 88  V888 Y8b  d8 88.     88  V888    88    88 `88. 88   88    88      .88.   `8b  d8' 88  V888 
// 	 `Y88P' YP   YP Y88888P  `Y88P' ~Y8888P' Y88888P C88888D  `Y88P'  `Y88P'  VP   V8P  `Y88P' Y88888P VP   V8P    YP    88   YD YP   YP    YP    Y888888P  `Y88P'  VP   V8P 
// 	                                                                                                                                                                         
// 	                                                                                                                                                                         
// ---------------------------------------------------------------------------
// calcul de l'indice de concentration
// les résultats sont stockés dans cluster_cucala et cluster_kulldorff
// ---------------------------------------------------------------------------
int calcul_concentration(Point_element *point_cas, Point_element *point_zero_cas,
	List_cluster_cucala *cluster_cucala, List_cluster_kulldorff *cluster_kulldorff,
	ID_TYPE nb_cas, ID_TYPE nb_controle, ID_TYPE nb_point_cas, ID_TYPE nb_point_controle, ID_TYPE nb_dimension,
	ID_TYPE simulation)
{
	int retour_infini=0;
	ID_TYPE i;
	// booleen pour determiner si le recalcul de la concentration est necessaire
	int analyse;
	// nombre de cas dans le cluster
	ID_TYPE nb_cas_cluster;
	// nombre de controle dans le cluster
	ID_TYPE nb_controle_cluster;

	// variables de stockage temporaire d'une distance
	Distance *distance_actuelle;
	Distance *distance_ppv, *smallest_distance;
	// variable de stockage temporaire de points
	Point_element *tmp_point, *tmp_point1, *tmp_point2;

	// Variables utilisées lors du calcul des concentrations de cucala et de kulldorf
	DISTANCE_TYPE concentration_cucala;		// Valeur de la concentration de cucala
	DISTANCE_TYPE concentration_kulldorff;	// Valeur de la concentration de kulldorff

	// Création de l'arbre kd des points avec cas
	Point_element *kd_cas = creer_arbre_kd(point_cas, NULL, nb_point_cas, nb_dimension, 0);

	// Création de l'arbre kd des points de controle
	Point_element *kd_controle = creer_arbre_kd(point_zero_cas, NULL, nb_point_controle, nb_dimension, 0);		// retourne NULL si aucun point de controle pur

	// Création d'un arbre binaire de recherche sur les distances cas-cas
	rbtree arbre_cas_cas = rbtree_create();
	// Parcours de la liste des points cas pour remplir l'arbre des distances cas-cas
	for(i=0; i<nb_point_cas; ++i)
	{
		distance_ppv = get_plus_proche_voisin(kd_cas, &point_cas[i], 0, nb_dimension);
		// On n'ajoute que les distances d'un point vers un point supérieur, ce qui évite les doublons de distances
		if(distance_ppv->point2 != NULL)
		{
			rbtree_insert_new(arbre_cas_cas, distance_ppv->valeur, (void*)distance_ppv, 0);	// 0 car on ajoute des distances, pas des points de controle
		}else{
			free(distance_ppv);
		}
	}
	// initialisation de la variable qui va servir à suivre la progression de l'algorithme
	ID_TYPE nb_points_cluster_courant = 0;

	// variable de stockage temporaire de cluster
	Cluster *cluster = NULL;
	Cluster *cluster1;
	Cluster *cluster2;

	// 	88""Yb    db    88""Yb  dP""b8  dP"Yb  88   88 88""Yb .dP"Y8 
	// 	88__dP   dPYb   88__dP dP   `" dP   Yb 88   88 88__dP `Ybo." 
	// 	88"""   dP__Yb  88"Yb  Yb      Yb   dP Y8   8P 88"Yb  o.`Y8b 
	// 	88     dP""""Yb 88  Yb  YboodP  YbodP  `YbodP' 88  Yb 8bodP' 
	// Tant que tous les points ne sont pas dans un seul cluster
	while(nb_points_cluster_courant < nb_point_cas)
	{
		// On retire la plus petite distance
		distance_actuelle = (Distance*)rbtree_get_smallest(arbre_cas_cas);	// retire de l'arbre
		// On prend les 2 points de la distance
		tmp_point1 = distance_actuelle->point1;
		tmp_point2 = distance_actuelle->point2;


		// Aucun des 2 points n'est dans un cluster, on crée un nouveau cluster avec ces deux points
		if(tmp_point1->cluster == NULL && tmp_point2->cluster == NULL)
		{
			// création d'un cluster
			cluster = creer_cluster(tmp_point1, tmp_point2, kd_controle, nb_dimension);
			analyse = 1;	// L'analyse est à refaire
		}else{
			// Un seul des 2 points est dans un cluster, on ajoute l'autre à ce cluster
			if(tmp_point1->cluster == NULL || tmp_point2->cluster == NULL)
			{
				// on ajoute le point qui n'est pas dans le cluster
				if(tmp_point1->cluster == NULL){
					cluster_add_point(tmp_point2->cluster, tmp_point1, kd_controle, nb_dimension);	// Ajout du point 1
					cluster = tmp_point2->cluster;
				}else{
					cluster_add_point(tmp_point1->cluster, tmp_point2, kd_controle, nb_dimension);	// Ajout du point 2
					cluster = tmp_point1->cluster;
				}					
				analyse = 1;	// L'analyse est à refaire
			}else{
				// Les 2 points sont dans un cluster, on récupère les clusters
				cluster1 = tmp_point1->cluster;
				cluster2 = tmp_point2->cluster;
				// Si les clusters sont différents, on les fusionne;
				if(cluster1 != cluster2)
				{
					cluster = fusion_cluster(cluster1, cluster2);		// Fusion des clusters, avec récupération du cluster final
					analyse = 1;	// L'analyse est à refaire
				}else{
					analyse = 0;	// rien à faire, les deux points sont déjà dans le même cluster
				}
			}
		}

		// 	   db    88b 88    db    88     Yb  dP .dP"Y8 888888 
		// 	  dPYb   88Yb88   dPYb   88      YbdP  `Ybo." 88__   
		// 	 dP__Yb  88 Y88  dP__Yb  88  .o   8P   o.`Y8b 88""   
		// 	dP""""Yb 88  Y8 dP""""Yb 88ood8  dP    8bodP' 888888                                         
		// Si l'analyse est nécessaire
		if(analyse == 1)
		{
			nb_cas_cluster = cluster->nb_cas;
			// On parcours l'arbre binaire de recherche qui contient les distances cas-control, et tant qu'il y a des points qui ont des controles dont la distance est
			// inférieure ou égale à la distance actuelle, on l'ajoute à nb_controle_purs du cluster et on l'ajoute dans l'arbre des controles du cluster
			smallest_distance = (Distance*)RBTREE_READ_SMALLEST(cluster->arbre_distance_cas_controle);
			while(smallest_distance != NULL && smallest_distance->valeur < distance_actuelle->valeur)
			{
				rbtree_get_smallest(cluster->arbre_distance_cas_controle);	// On retire la distance de l'arbre
				// On vérifier que le point de controle n'a pas déjà été ajouté au cluster
				// Deux solutions: vérifier le controle dans arbre_points_controle du cluster, ou le cluster dans controle_cluster du point de controle
				// On va utiliser la deuxième solution, en partant du principe qu'il y a plus de chances que cet arbre soit plus petit
				if(rbtree_find(smallest_distance->point2->controle_cluster, cluster->id_controle) == NULL)
				{
					// Le point n'est pas encore dans les controles du cluster, on l'ajoute alors
					cluster_add_point_controle(cluster, smallest_distance->point2);
				}
				// On recalcule le ppv control (hors déjà compté dans cluster) du point cas
				distance_ppv = get_plus_proche_voisin(kd_controle, smallest_distance->point1, 1, nb_dimension);
				// On l'ajoute à l'arbre
				rbtree_insert_new(cluster->arbre_distance_cas_controle, distance_ppv->valeur, (void*)distance_ppv, 0);
				// Libération mémoire
				free(smallest_distance);
				// On récupère la nouvelle plus petite distance
				smallest_distance = (Distance*)RBTREE_READ_SMALLEST(cluster->arbre_distance_cas_controle);
			}
			// Le nombre de controles du cluster est le nombre de controle des points à cas + le nombre de controles purs à portée du cluster
			nb_controle_cluster = cluster->nb_controle + cluster->nb_controle_pur;
			
			// tester si le cluster est positif			
			//if(((double)(nb_cas_cluster+nb_controle_cluster)/(double)(nb_cas+nb_controle))<((double)nb_cas_cluster/(double)nb_cas))
			// -> en fait on cherche à savoir si le ratio nb_controle_cluster/nb_controle est bien inférieur au ratio nb_cas_cluster/nb_cas soit dit autrement
			if(nb_controle_cluster*nb_cas < nb_cas_cluster*nb_controle)
			{
				// 	88  dP 88   88 88     88     8888b.   dP"Yb  88""Yb 888888 888888 
				// 	88odP  88   88 88     88      8I  Yb dP   Yb 88__dP 88__   88__   
				// 	88"Yb  Y8   8P 88  .o 88  .o  8I  dY Yb   dP 88"Yb  88""   88""   
				// 	88  Yb `YbodP' 88ood8 88ood8 8888Y"   YbodP  88  Yb 88     88     
				/***  Calcul de la concentration de Kulldorff  ***/
				if(nb_controle_cluster>0 && nb_controle_cluster< nb_controle && nb_cas_cluster < nb_cas)
				{
					concentration_kulldorff = 0;
					concentration_kulldorff = nb_cas_cluster*log((double)nb_cas_cluster/(double)(nb_cas_cluster+nb_controle_cluster));
					concentration_kulldorff += nb_controle_cluster*log((double)nb_controle_cluster/(double)(nb_cas_cluster+nb_controle_cluster));
					concentration_kulldorff += (nb_cas-nb_cas_cluster)*log((double)(nb_cas-nb_cas_cluster)/(double)(nb_cas+nb_controle-nb_cas_cluster-nb_controle_cluster));
					concentration_kulldorff += (nb_controle-nb_controle_cluster)*log((double)(nb_controle-nb_controle_cluster)/(double)(nb_cas+nb_controle -nb_cas_cluster -nb_controle_cluster));
					
					if( concentration_kulldorff > cluster_kulldorff->concentration_kulldorff || 
						(concentration_kulldorff == cluster_kulldorff->concentration_kulldorff && nb_cas_cluster > cluster_kulldorff->nb_cas_kulldorff) )	//nouvelle concentration max ou même concentration mais plus de points
					{					  
						// sauvegarde du cluster dans cluster_kulldorff
						cluster_kulldorff->dist_max = distance_actuelle->valeur;
						cluster_kulldorff->nb_cas_kulldorff = nb_cas_cluster;
						cluster_kulldorff->nb_controle_kulldorff = nb_controle_cluster;
						cluster_kulldorff->concentration_kulldorff = concentration_kulldorff;
						cluster_kulldorff->taille = cluster->taille;
						cluster_kulldorff->tmp_elem_debut = cluster->premier_point;
						cluster_kulldorff->tmp_elem_fin = cluster->dernier_point;
					}
					
				}
				// 	 dP""b8 88   88  dP""b8    db    88        db    
				// 	dP   `" 88   88 dP   `"   dPYb   88       dPYb   
				// 	Yb      Y8   8P Yb       dP__Yb  88  .o  dP__Yb  
				// 	 YboodP `YbodP'  YboodP dP""""Yb 88ood8 dP""""Yb 
				/***  Calcul de la concentration de Cucala  ***/
				concentration_cucala = - (DISTANCE_TYPE)pbeta(((double)(nb_cas_cluster+nb_controle_cluster)/(double)(nb_cas+nb_controle)), nb_cas_cluster-1, nb_cas-nb_cas_cluster+2, 1, 1);	// log de la beta incomplète
				//concentration_cucala = 1 / (DISTANCE_TYPE)pbeta(((double)(nb_cas_cluster+nb_controle_cluster)/(double)(nb_cas+nb_controle)), nb_cas_cluster-1, nb_cas-nb_cas_cluster+2, 1, 0);	
				//concentration_cucala = 1 / gsl_sf_beta_inc(nb_cas_cluster-1, nb_cas-nb_cas_cluster+2, ((double)(nb_cas_cluster+nb_controle_cluster)/(double)(nb_cas+nb_controle)));

				if(isinf(concentration_cucala))
					retour_infini=1;
				else
					if( concentration_cucala > cluster_cucala->concentration_cucala || 			//nouvelle concentration max ou
						(concentration_cucala == cluster_cucala->concentration_cucala && nb_cas_cluster > cluster_cucala->nb_cas_cucala) )	// même concentration mais plus de points
					{
						// sauvegarde du cluster dans cluster_cucala
						cluster_cucala->dist_max = distance_actuelle->valeur;
						cluster_cucala->nb_cas_cucala = nb_cas_cluster;
						cluster_cucala->nb_controle_cucala = nb_controle_cluster;
						cluster_cucala->concentration_cucala = concentration_cucala;
						cluster_cucala->taille = cluster->taille;
						cluster_cucala->tmp_elem_debut = cluster->premier_point;
						cluster_cucala->tmp_elem_fin = cluster->dernier_point;
					}

			}

		}
		// On libère la mémoire
		free(distance_actuelle);

		// mise à jour de la variable qui sert à suivre la progression de l'algorithme
		nb_points_cluster_courant = tmp_point1->cluster->taille;
		// On recalcule les nouveaux ppv pour les 2 points de la distance qui vient d'être retirée
		// point1
		distance_ppv = get_plus_proche_voisin(kd_cas, tmp_point1, 0, nb_dimension);
		if(distance_ppv->point2 != NULL)
		{
			rbtree_insert_new(arbre_cas_cas, distance_ppv->valeur, (void*)distance_ppv, 0);
		}else{
			free(distance_ppv);
		}
	}
	
	// 	.dP"Y8    db    88   88 Yb    dP 888888  dP""b8    db    88""Yb 8888b.  888888 
	// 	`Ybo."   dPYb   88   88  Yb  dP  88__   dP   `"   dPYb   88__dP  8I  Yb 88__   
	// 	o.`Y8b  dP__Yb  Y8   8P   YbdP   88""   Yb  "88  dP__Yb  88"Yb   8I  dY 88""   
	// 	8bodP' dP""""Yb `YbodP'    YP    888888  YboodP dP""""Yb 88  Yb 8888Y"  888888 
	// Si on n'est pas dans une simulation, on sauvegarde les clusters de Cucala et de Kulldorff
	// On peut se permettre de ne faire ça qu'une seule fois à la fin car l'ordre des clusters n'est jamais
	// modifié, les fusions de clusters se faisant aux extrémités des listes, le début et la fin des clusters restant donc valables.
	if(simulation == 0)
	{
		// sauvegarde du cluster de Cucala
		if(cluster_cucala->taille > 0)
		{
			cluster_cucala->debut = cluster_cucala->tmp_elem_debut->point;
			cluster_cucala->fin = cluster_cucala->tmp_elem_fin->point;
			tmp_point = cluster_cucala->tmp_elem_debut;
			while(tmp_point != cluster_cucala->tmp_elem_fin)
			{
				if(tmp_point->suivant_cluster != NULL)
					tmp_point->point->suiv_cluster_cucala = tmp_point->suivant_cluster->point;	
				else
					tmp_point->point->suiv_cluster_cucala = NULL;	// normalement inutile, car valeur par défaut
				tmp_point = tmp_point->suivant_cluster;
			}
		}
		// sauvegarde du cluster de Kulldorff
		if(cluster_kulldorff->taille > 0)
		{
			cluster_kulldorff->debut = cluster_kulldorff->tmp_elem_debut->point;
			cluster_kulldorff->fin = cluster_kulldorff->tmp_elem_fin->point;
			tmp_point = cluster_kulldorff->tmp_elem_debut;
			while(tmp_point != cluster_kulldorff->tmp_elem_fin)
			{
				if(tmp_point->suivant_cluster != NULL)
					tmp_point->point->suiv_cluster_kulldorff = tmp_point->suivant_cluster->point;	
				else
					tmp_point->point->suiv_cluster_kulldorff = NULL;	// normalement inutile, car valeur par défaut
				tmp_point = tmp_point->suivant_cluster;
			}
		}
	}

	// 	888888 88""Yb 888888 888888 
	// 	88__   88__dP 88__   88__   
	// 	88""   88"Yb  88""   88""   
	// 	88     88  Yb 888888 888888 
	// libération de la memoire
	// Il ne reste plus qu'un seul cluster, on libère ses arbres
	if(cluster != NULL)
	{
		rbtree_free(cluster->arbre_distance_cas_controle, 1);
		rbtree_free(cluster->arbre_points_controle, 0);
		free(cluster);
	}
	// On libère l'arbre de distances
	rbtree_free(arbre_cas_cas, 1);
	// On libère les arbres de tous les points controle (pas d'arbres dans les points cas)
	for(i=0; i<nb_point_controle; i++)
	{
		rbtree_free(point_zero_cas[i].controle_cluster, 0);
		rbtree_free(point_zero_cas[i].controle_cluster_sous_arbre, 0);
	}
	return retour_infini;
}
