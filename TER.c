#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <gmp.h>
#include <time.h>
#include "rationnels.h"
#include "systemes.h"
#include "bareiss.h"
#include "gauss_sys_rat.h"
#include "io.h"
#include "zpz.h"
#include "mod.h"
#include "mod_dets.h"
#include "mod_thrd.h"

//#define GAUSS
#define BAREISS
//#define ZPZ
//#define MOD_OLD
#define MOD
//#define MOD_DETS
//#define MOD_PARA

#define THR 8

int main(){
	// Initialisation
	FILE* f = stdout;
	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state,time(NULL));
	//ecrit_fichier_au_pif("systeme2.txt",5,state,96);
	//ecrit_fichier_au_pif("systeme3.txt",5,state,512);
	//ecrit_fichier_au_pif("systeme4.txt",5,state,2048);
	//ecrit_fichier_au_pif("systeme5.txt",50,state,96);
	//ecrit_fichier_au_pif("systeme6.txt",50,state,512);
	//ecrit_fichier_au_pif("systeme7.txt",50,state,2048);
	//ecrit_fichier_au_pif("systeme8.txt",200,state,96);
	//ecrit_fichier_au_pif("systeme9.txt",200,state,512);
	//ecrit_fichier_au_pif("systeme10.txt",50,state,12);
	//ecrit_fichier_au_pif("systeme11.txt",200,state,12);
	//ecrit_fichier_au_pif("systeme12.txt",700,state,12);
	systeme s,s_ini;
	syst_zpz s_zpz,s_zpzv;
	init_lit_systeme(&s,"systs/systeme-n50c96.txt");
	int n = s.n;		// Nombre de lignes
	//int m = s.m;		// Nombre de colonnes, en comptant le second membre (en pratique : m = n+1)
	init_copie_systeme(&s_ini,&s);	// Copie qui servira à conserver le système initial, pour pouvoir tester notre solution à la fin (sera aussi donnée à l'algo de Gauss sur les rationnels et à zpz_thrd, car ils ne le modifieront pas)
	rationnel* sol_g = malloc(n*sizeof(rationnel));		// Contiendra la solution donnée par l'algorithme de Gauss
	rationnel* sol_b = malloc(n*sizeof(rationnel));		// Contiendra la solution donnée par l'algorithme de Bareiss
	rationnel* sol_mo = malloc(n*sizeof(rationnel));	// Contiendra la solution obtenue par méthode modulaire (ancienne variante) (MOD_OLD)
	rationnel* sol_mp = malloc(n*sizeof(rationnel));	// Contiendra la solution obtenue par méthode modulaire (version en parallèle) (MOD_PARA)
	rationnel* sol_mh = malloc(n*sizeof(rationnel));	// Contiendra la solution obtenue par méthode modulaire (version avec borne de Hadamard) (MOD)
	rationnel* sol_md = malloc(n*sizeof(rationnel));	// Contiendra la solution obtenue par méthode modulaire (version avec déterminants (et borne de Hadamard)) (MOD_DETS)
	long int* sol_zpz = malloc(n*sizeof(long int));		// Contiendra la solution donnée par l'algorithme de Gauss dans Z/pZ
	for(int i = 0;i < n;i++){
		rat_init(&sol_b[i]);
		rat_init(&sol_g[i]);
		rat_init(&sol_mo[i]);
		rat_init(&sol_mp[i]);
		rat_init(&sol_mh[i]);
		rat_init(&sol_md[i]);
	}
	mpz_t p_mpz;
	mpz_init(p_mpz);
	long int p = genere_p(p_mpz,state,62);
	//int p = 701;
	init_copie_syst_zpz(&s_zpz,&s_ini,p);	// Copie du système sur laquelle l'aglo de Gauss dans Z/pZ sera testé
	init_copie_syst_zpz(&s_zpzv,&s_ini,p);	// Copie du système qui servira à conserver le résultat initial modulo p (pour tester le résultat)
	clock_t debut,fin;		// Pour les mesures de temps
#ifdef GAUSS
	double gauss_tps;
#endif
#ifdef BAREISS
	double bareiss_tps;
#endif
#ifdef MOD_OLD
	double mod_old_tps;
#endif
#ifdef MOD
	double mod_tps;
#endif
#ifdef MOD_DETS
	double mod_dets_tps;
#endif
#ifdef MOD_PARA
	double mod_para_tps;
#endif
	
	// Affichage du système de départ
	fprintf(f,"\n\nSYSTÈME DE DÉPART :\n");
	affiche_systeme(&s_ini,f);		// Système de départ
	
#ifdef GAUSS
	// Calcul d'une solution, avec l'algo de Gauss
	// Prend 3 plombes
	fprintf(f,"\nSOLUTION (GAUSS) :\n");
	debut = clock();
	gauss(&s_ini,sol_g);
	fin = clock();
	gauss_tps = ((double)(fin-debut))/CLOCKS_PER_SEC;
	/*for(int i = 0;i < n;i++){
		rat_aff(sol_g[i],f);
		fprintf(f,"\n");
	}*/
#endif
	
#ifdef BAREISS
	// Calcul d'une solution, avec l'algo de Bareiss
	// Prend un temps raisonnable
	fprintf(f,"\n\nSYSTÈME ÉCHELONNÉ (BAREISS) :\n");
	debut = clock();
	bareiss(&s);
	sol_syst_echelonne(&s,sol_b);
	fin = clock();
	bareiss_tps = ((double)(fin-debut))/CLOCKS_PER_SEC;
	//affiche_systeme(&s,f);		// Système échelonné
	fprintf(f,"\n\nSOLUTION (BAREISS) :\n");
	/*for(int i = 0;i < n;i++){
		rat_aff(sol_b[i],f);
		fprintf(f,"\n");
	}*/
#endif
	
#ifdef ZPZ
	// Calcul d'une solution dans Z/pZ (avec p choisi "au hasard" avec à peu près 32 bits)
	// Va très vite
	fprintf(f,"\n\nSYSTÈME DE DÉPART MODULO P :\n(P = %ld)\n",p);
	affiche_syst_zpz(&s_zpzv,f);
	zpz_resol(&s_zpz,sol_zpz);
	fprintf(f,"\nSYSTÈME ÉCHELONNÉ DANS Z/pZ :\n");
	affiche_syst_zpz(&s_zpz,f);
	fprintf(f,"\nSOLUTION dans Z/pZ :\n");
	for(int i = 0;i < n;i++){
		fprintf(f,"%ld\n",sol_zpz[i]);
	}
#endif
	
#ifdef MOD_OLD
	// Calcul d'une solution par méthode modulaire (en prenant des nombres premiers d'au plus 30 bits)
	// Prend trop de temps, car la reconstruction d'une solution rationnelle est trop coûteuse et faite trop souvent
	fprintf(f,"\n\nSOLUTION (MÉTHODE MODULAIRE, ANCIENNE VERSION) :\n");
	debut = clock();
	modulaire_old(&s_ini,sol_mo,state,62);
	fin = clock();
	mod_old_tps = ((double)(fin-debut))/CLOCKS_PER_SEC;
	/*for(int i = 0;i < n;i++){
		rat_aff(sol_mo[i],f);
		fprintf(f,"\n");
	}*/
#endif
	
#ifdef MOD
	// Calcul d'une solution par méthode modulaire (alternative avec borne de Hadamard)
	fprintf(f,"\n\nSOLUTION (MÉTHODE MODULAIRE) :\n");
	debut = clock();
	modulaire(&s_ini,sol_mh,state,62);
	fin = clock();
	mod_tps = ((double)(fin-debut))/CLOCKS_PER_SEC;
	/*for(int i = 0;i < n;i++){
		rat_aff(sol_mh[i],f);
		fprintf(f,"\n");
	}*/
#endif

#ifdef MOD_DETS
	// Calcul d'une solution par méthode modulaire (alternative avec déterminants)
	fprintf(f,"\n\nSOLUTION (MÉTHODE MODULAIRE AVEC DÉTERMINANTS) :\n");
	debut = clock();
	modulaire_dets(&s_ini,sol_md,state,62);
	fin = clock();
	mod_dets_tps = ((double)(fin-debut))/CLOCKS_PER_SEC;
	for(int i = 0;i < n;i++){
		rat_aff(sol_md[i],f);
		fprintf(f,"\n");
	}
#endif

#ifdef MOD_PARA
	// Calcul d'une solution par méthode modulaire en parallèle
	fprintf(f,"\n\nSOLUTION (MÉTHODE MODULAIRE EN PARALLÈLE) :\n");
	debut = clock();
	modulaire_thrd(&s_ini,sol_mp,state,62,THR);
	fin = clock();
	mod_para_tps = ((double)(fin-debut))/CLOCKS_PER_SEC;
	/*for(int i = 0;i < n;i++){
		rat_aff(sol_mp[i],f);
		fprintf(f,"\n");
	}*/
#endif	
	fputc('\n',f);
	
#ifdef GAUSS
	// Temps d'exécution :
	fprintf(f,"Temps (Gauss) : %lf s\n",gauss_tps);
	// Vérifications (Gauss)
	assert(verif_sol(&s_ini,sol_g));	// Vérif sur le système de départ
#endif
	
#ifdef BAREISS
	// Temps d'exécution :
	fprintf(f,"Temps (Bareiss) : %lf s\n",bareiss_tps);
	// Vérifications (Bareiss)
	assert(verif_sol(&s,sol_b));		// Vérif "triviale"
	assert(verif_sol(&s_ini,sol_b));		// Vérif sur le systeme de départ
#endif
	
#ifdef ZPZ
	// Vérifications (Gauss dans Z/pZ)
	assert(verif_sol_zpz(&s_zpz,sol_zpz));	// Vérif "triviale"
	assert(verif_sol_zpz(&s_zpzv,sol_zpz));	// Vérif sur le système de départ
#endif
	
#ifdef MOD_OLD
	// Temps d'exécution :
	fprintf(f,"Temps (Ancienne méthode modulaire) : %lf s\n",mod_old_tps);
	// Vérification (ancienne méthode modulaire)
	assert(verif_sol(&s_ini,sol_mo));		// Vérif sur le système de départ
#endif
	
#ifdef MOD
	// Temps d'exécution :
	fprintf(f,"Temps (Méthode modulaire) : %lf s\n",mod_tps);
	// Vérification (méthode modulaire)
	assert(verif_sol(&s_ini,sol_mh));
#endif

#ifdef MOD_DETS
	// Temps d'exécution :
	fprintf(f,"Temps (Méthode modulaire avec déterminants) : %lf s\n",mod_dets_tps);
	// Vérifications (Gauss)
	assert(verif_sol(&s_ini,sol_md));	// Vérif sur le système de départ
#endif

#ifdef MOD_PARA
	// Temps d'exécution :
	fprintf(f,"Temps (Méthode modulaire en parallèle) : %lf s\n",mod_para_tps);
	// Vérification (méthode modulaire en parallèle)
	assert(verif_sol(&s_ini,sol_mp));		// Vérif sur le système de départ
#endif

	
	// Suppression des objets utilisés
	for(int i = 0;i < n;i++){
		rat_clear(&sol_b[i]);
		rat_clear(&sol_g[i]);
		rat_clear(&sol_mo[i]);
		rat_clear(&sol_mp[i]);
		rat_clear(&sol_mh[i]);
		rat_clear(&sol_md[i]);
	}
	free(sol_b);
	free(sol_g);
	free(sol_mo);
	free(sol_mp);
	free(sol_mh);
	free(sol_md);
	free(sol_zpz);
	mpz_clear(p_mpz);
	detruit_systeme(&s);
	detruit_systeme(&s_ini);
	detruit_syst_zpz(&s_zpz);
	detruit_syst_zpz(&s_zpzv);
	gmp_randclear(state);
	return 0;
}


