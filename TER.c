#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <gmp.h>
#include "rationnels.h"
#include "systemes.h"
#include "bareiss.h"
#include "gauss_sys_rat.h"
#include "io.h"
#include "zpz.h"

int main(){	
	
	// Initialisation
	gmp_randstate_t state;
	gmp_randinit_default(state);	// L'aléatoire semble cassé, il y a de l'amélioration possible par ici...
	ecrit_fichier_au_pif("systeme2.txt",7,state,16);
	systeme s,s_g,s_zpz,s_ini;
	init_lit_systeme(&s,"systeme.txt");
	//init_lit_systeme(&s,"systeme2.txt");
	int n = s.n;
	int m = s.m;
	init_systeme(&s_ini,n,m);
	init_systeme(&s_g,n,m);
	init_systeme(&s_zpz,n,m);
	copie_systeme(&s,&s_ini);	// Copie qui servira à conserver le système initial, pour pouvoir tester notre solution à la fin
	copie_systeme(&s,&s_g);		// Copie qui sera utilisée par l'algo de Gauss (pas vraiment nécessaire, mais c'est plus propre que de lui donner s_ini)	
	copie_systeme(&s,&s_zpz);	// Copie qui sera utilisée pour tester l'algo de Gauss dans Z/pZ
	/*
	// Calcul d'une solution, avec l'algo de Bareiss
	rationnel sol[n];
	for(int i = 0;i < n;i++){
		rat_init(&sol[i]);
	}
	bareiss(&s);
	sol_syst_echelonne(&s,sol);
	
	// Calcul d'une solution, avec l'algo de Gauss
	rationnel sol_g[n];
	for(int i = 0;i < n;i++){
		rat_init(&sol_g[i]);
	}
	gauss(&s_g,sol_g);*/
	
	// Affichage du système et de la solution
	fprintf(stdout,"\n\nSYSTÈME DE DÉPART :\n");
	affiche_systeme(&s_ini);		// Système de départ
	/*
	fprintf(stdout,"\n\nSYSTÈME ÉCHELONNÉ (BAREISS) :\n");
	affiche_systeme(&s);		// Système échelonné
	fprintf(stdout,"\n\nSOLUTION (BAREISS) :\n");
	for(int i = 0;i < n;i++){
		rat_aff(sol[i]);
		fprintf(stdout,"\n");
	}
	fprintf(stdout,"\nSOLUTION (GAUSS) :\n");
	for(int i = 0;i < n;i++){
		rat_aff(sol_g[i]);
		fprintf(stdout,"\n");
	}
	fprintf(stdout,"\n");
	
	// Vérifications (Bareiss)
	//assert(verif_sol(&s,sol));		// Vérif "triviale"
	assert(verif_sol(&s_ini,sol));		// Vérif sur le systeme de départ
	//ecrit_coeff(&s_ini,0,0,t_m[6]);
	//assert(!verif_sol(&s_ini,sol));		// Autre vérif carrément stupide
	
	// Vérifications (Gauss)
	assert(verif_sol(&s_ini,sol_g));*/
	
	// Test rapide de Gauss Z/pZ :
	mpz_t k,p;
	mpz_urandomb(k,state,32);
	mpz_nextprime(p,k);
	zpz_gauss(&s_zpz,p);
	fprintf(stdout,"\n\nSYSTÈME ÉCHELONNÉ (GAUSS Z/pZ) :\n(p = ");
	mpz_out_str(stdout,10,p);
	fprintf(stdout,")\n");
	affiche_systeme(&s_zpz);
	fprintf(stdout,"\n");
	
	// Suppression des objets utilisés
	/*for(int i = 0;i < n;i++){
		rat_clear(&sol[i]);
		rat_clear(&sol_g[i]);
	}*/
	detruit_systeme(&s);
	detruit_systeme(&s_ini);
	detruit_systeme(&s_g);
	detruit_systeme(&s_zpz);
	gmp_randclear(state);
	return 0;
}
