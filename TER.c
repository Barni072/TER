#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <gmp.h>
#include "rationnels.h"
#include "systemes.h"
#include "bareiss.h"
#include "gauss_sys_rat.h"

int main(){	
	// Initialisation
	systeme s,s_g;
	init_systeme(&s,3,4);
	init_systeme(&s_g,3,4);
	int t_i[12] = {
		1,2,-3,		9,
		4,7,-8,		-5,
		-12,0,5,	4};
	mpz_t t_m[12];
	for(int i = 0;i < 12;i++){
		mpz_init_set_si(t_m[i],t_i[i]);
	}
	for(int i = 0;i < 3;i++){
		for(int j = 0;j < 4;j++){
			ecrit_coeff(&s,i,j,t_m[i*4+j]);
		}
	}
	copie_systeme(&s,&s_g);		// Copie qui sera utilisée par l'algo de Gauss (pas vraiment nécessaire, mais c'est plus propre que de lui donner s_ini)
	
	// Création une copie de s, qu'on ne modiefiera pas, et qui va permettre de tester notre solution à la fin
	systeme s_ini;
	init_systeme(&s_ini,3,4);
	copie_systeme(&s,&s_ini);
	
	// Calcul d'une solution, avec l'algo de Bareiss
	rationnel sol[3];
	for(int i = 0;i < 3;i++){
		rat_init(&sol[i]);
	}
	bareiss(&s);
	sol_syst_echelonne(&s,sol);
	
	// Calcul d'une solution, avec l'algo de Gauss
	rationnel sol_g[3];
	for(int i = 0;i < 3;i++){
		rat_init(&sol_g[i]);
	}
	gauss(&s_g,sol_g);
	
	// Affichage du système et de la solution
	fprintf(stdout,"\nSystème de départ :\n");
	affiche_systeme(&s_ini);		// Système de départ
	fprintf(stdout,"\nSystème échelonné (par Bareiss) :\n");
	affiche_systeme(&s);		// Système échelonné
	fprintf(stdout,"\nSolution (Bareiss) :\n");
	for(int i = 0;i < 3;i++){
		rat_aff(sol[i]);
	}
	fprintf(stdout,"\nSolution (Gauss) :\n");
	for(int i = 0;i < 3;i++){
		rat_aff(sol_g[i]);
	}
	fprintf(stdout,"\n");
	
	// Vérifications (Bareiss)
	//assert(verif_sol(&s,sol));		// Vérif "triviale"
	assert(verif_sol(&s_ini,sol));		// Vérif sur le systeme de départ
	//ecrit_coeff(&s_ini,0,0,t_m[6]);
	//assert(!verif_sol(&s_ini,sol));		// Autre vérif carrément stupide
	
	// Vérifications (Gauss)
	assert(verif_sol(&s_ini,sol_g));
	
	// Suppression des objets utilisés
	for(int i = 0;i < 3;i++){
		rat_clear(&sol[i]);
	}
	for(int i = 0;i < 12;i++){
		mpz_clear(t_m[i]);
	}
	detruit_systeme(&s);
	detruit_systeme(&s_g);
	detruit_systeme(&s_ini);
	return 0;
}		// TODO : récupérer le système à tester depuis un fichier texte, comparer les temps d'exécution des algos de Gauss et Bareiss
