// ** ================================================================================================================
// ** R-Package: graphscan
// ** Fichier : src/get_number_proc.c
// ** Description : renvoie le nombre de processeurs (threads) disponibles sur la machine.
// ** License : GPL-2 | GPL-3
// ** Auteurs : Robin Loche, Benoit Giron, David Abrial, Lionel Cucala, Myriam Charras-Garrido, Jocelyn De-Goer
// ** ================================================================================================================
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#ifdef _OPENMP
	#include <omp.h>
#endif


// 	 dP""b8 888888 888888            88b 88 88   88 8b    d8 88""Yb 888888 88""Yb            88""Yb 88""Yb  dP"Yb   dP""b8 
// 	dP   `" 88__     88              88Yb88 88   88 88b  d88 88__dP 88__   88__dP            88__dP 88__dP dP   Yb dP   `" 
// 	Yb  "88 88""     88              88 Y88 Y8   8P 88YbdP88 88""Yb 88""   88"Yb             88"""  88"Yb  Yb   dP Yb      
// 	 YboodP 888888   88   oooooooooo 88  Y8 `YbodP' 88 YY 88 88oodP 888888 88  Yb oooooooooo 88     88  Yb  YbodP   YboodP                                                                                                                                               

void get_number_proc(int *res)
{
	#ifdef _OPENMP
		*res = (int)omp_get_num_procs();
	#else
		*res = (int)1;
	#endif

	return;
}

