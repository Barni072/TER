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
#include "zpz_thrd.h"

int main(){
	// Initialisation
	FILE* f = stdout;
	int thr = 8;
	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state,time(NULL));
	//ecrit_fichier_au_pif("systeme2.txt",25,state,512);		// Les résultats dépassent du terminal
	ecrit_fichier_au_pif("systeme2.txt",10,state,128);
	systeme s,s_ini;
	syst_zpz s_zpz,s_zpzv;
	//init_lit_systeme(&s,"systeme.txt");
	init_lit_systeme(&s,"systeme2.txt");
	int n = s.n;		// Nombre de lignes
	int m = s.m;		// Nombre de colonnes, en comptant le second membre (en pratique : m = n+1)
	init_copie_systeme(&s_ini,&s);	// Copie qui servira à conserver le système initial, pour pouvoir tester notre solution à la fin (sera aussi donnée à l'algo de Gauss sur les rationnels et à zpz_thrd, car ils ne le modifieront pas)
	rationnel* sol_b = malloc(n*sizeof(rationnel));	// Contiendra la solution donnée par l'algorithme de Bareiss
	rationnel* sol_g = malloc(n*sizeof(rationnel));	// Contiendra la solution donnée par l'algorithme de Gauss
	rationnel* sol_m = malloc(n*sizeof(rationnel));	// Contiendra la solution obtenue par méthode modulaire
	rationnel* sol_mp = malloc(n*sizeof(rationnel));	// Contiendra la solution obtenue par méthode modulaire
	int* sol_zpz = malloc(n*sizeof(int));		// Contiendra la solution donnée par l'algorithme de Gauss dans Z/pZ
	for(int i = 0;i < n;i++){
		rat_init(&sol_b[i]);
		rat_init(&sol_g[i]);
		rat_init(&sol_m[i]);
		rat_init(&sol_mp[i]);
	}
	mpz_t p_mpz;
	mpz_init(p_mpz);
	int p = genere_p(p_mpz,state,30);
	init_copie_syst_zpz(&s_zpz,&s_ini,p);	// Copie du système sur laquelle l'aglo de Gauss dans Z/pZ sera testé
	init_copie_syst_zpz(&s_zpzv,&s_ini,p);	// Copie du système qui servira à conserver le résultat initial modulo p (pour tester le résultat)
	
	// Affichage du système de départ
	fprintf(f,"\n\nSYSTÈME DE DÉPART :\n");
	affiche_systeme(&s_ini,f);		// Système de départ
	
	// Calcul d'une solution, avec l'algo de Bareiss
	// Prend un temps raisonnable
	fprintf(f,"\n\nSYSTÈME ÉCHELONNÉ (BAREISS) :\n");
	bareiss(&s);
	sol_syst_echelonne(&s,sol_b);
	affiche_systeme(&s,f);		// Système échelonné
	fprintf(f,"\n\nSOLUTION (BAREISS) :\n");
	for(int i = 0;i < n;i++){
		rat_aff(sol_b[i],f);
		fprintf(f,"\n");
	}
	
	// Calcul d'une solution, avec l'algo de Gauss
	// Prend 3 plombes
	/*fprintf(f,"\nSOLUTION (GAUSS) :\n");
	gauss(&s_ini,sol_g);
	for(int i = 0;i < n;i++){
		rat_aff(sol_g[i],f);
		fprintf(f,"\n");
	}*/
	
	// Calcul d'une solution dans Z/pZ (avec p choisi "au hasard" avec à peu près 32 bits)
	// Va très vite
	/*fprintf(f,"\n\nSYSTÈME DE DÉPART MODULO P :\n(P = %d)\n",p);
	affiche_syst_zpz(&s_zpzv,f);
	zpz_resol(&s_zpz,sol_zpz);
	fprintf(f,"\nSYSTÈME ÉCHELONNÉ DANS Z/pZ :\n");
	affiche_syst_zpz(&s_zpz,f);
	fprintf(f,"\nSOLUTION dans Z/pZ :\n");
	for(int i = 0;i < n;i++){
		fprintf(f,"%d\n",sol_zpz[i]);
	}*/
	
	// Calcul d'une solution par méthode modulaire (en prenant des nombres premiers d'au plus 30 bits)
	// Prend bien plus de temps que Bareiss
	fprintf(f,"\n\nSOLUTION (MÉTHODE MODULAIRE) :\n");
	modulaire(&s_ini,sol_m,state,30);
	for(int i = 0;i < n;i++){
		rat_aff(sol_m[i],f);
		fprintf(f,"\n");
	}
	
	// Calcul d'une solution par méthode modulaire en parallèle
	fprintf(f,"\n\nSOLUTION (MÉTHODE MODULAIRE EN PARALLÈLE) :\n");
	modulaire_thrd(&s_ini,sol_mp,state,30,thr);
	for(int i = 0;i < n;i++){
		rat_aff(sol_mp[i],f);
		fprintf(f,"\n");
	}
	fputc('\n',f);
	
	// Vérifications (Bareiss)
	//assert(verif_sol(&s,sol_b));		// Vérif "triviale"
	assert(verif_sol(&s_ini,sol_b));		// Vérif sur le systeme de départ
	
	// Vérifications (Gauss)
	//assert(verif_sol(&s_ini,sol_g));	// Vérif sur le système de départ
	
	// Vérifications (Gauss dans Z/pZ)
	//assert(verif_sol_zpz(&s_zpz,sol_zpz));	// Vérif "triviale"
	//assert(verif_sol_zpz(&s_zpzv,sol_zpz));	// Vérif sur le système de départ
	
	// Vérification (méthode modulaire)
	assert(verif_sol(&s_ini,sol_m));		// Vérif sur le système de départ
	
	// Vérification (méthode modulaire en parallèle)
	//assert(verif_sol(&s_ini,sol_mp));		// Vérif sur le système de départ
	
	// Essai d'exécution de Gauss sur des Z/pZ en parallèle
	//test_zpz_multi(&s_ini,n,8,state);
	// TEST OK
	
	// Suppression des objets utilisés
	for(int i = 0;i < n;i++){
		rat_clear(&sol_b[i]);
		rat_clear(&sol_g[i]);
		rat_clear(&sol_m[i]);
		rat_clear(&sol_mp[i]);
	}
	free(sol_b);
	free(sol_g);
	free(sol_m);
	free(sol_mp);
	free(sol_zpz);
	mpz_clear(p_mpz);
	detruit_systeme(&s);
	detruit_systeme(&s_ini);
	detruit_syst_zpz(&s_zpz);
	detruit_syst_zpz(&s_zpzv);
	gmp_randclear(state);
	return 0;
}


