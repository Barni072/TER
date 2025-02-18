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
	int n = s.n;		// Nombre de lignes
	int m = s.m;		// Nombre de colonnes, en comptant le second membre (en pratique : m = n+1)
	init_copie_systeme(&s_ini,&s);	// Copie qui servira à conserver le système initial, pour pouvoir tester notre solution à la fin
	init_copie_systeme(&s_g,&s);	// Copie qui sera utilisée par l'algo de Gauss (pas vraiment nécessaire, mais c'est plus propre que de lui donner s_ini)
	init_copie_systeme(&s_zpz,&s);	// Copie qui sera utilisée pour tester l'algo de Gauss dans Z/pZ
	rationnel sol[n],sol_g[n];	// Futur emplacement des solutions
	mpz_t sol_zpz[n];		// Idem
	for(int i = 0;i < n;i++){
		rat_init(&sol[i]);
		rat_init(&sol_g[i]);
		mpz_init(sol_zpz[i]);
	}
	mpz_t k,p;
	mpz_init(k);
	mpz_init(p);
	
	// Calcul d'une solution, avec l'algo de Bareiss
	bareiss(&s);
	sol_syst_echelonne(&s,sol);
	
	// Calcul d'une solution, avec l'algo de Gauss
	gauss(&s_g,sol_g);
	
	// Calcul d'une solution dans Z/pZ (avec p choisi "au hasard")
	mpz_urandomb(k,state,32);
	mpz_nextprime(p,k);
	zpz_sans_copie(sol_zpz,&s_zpz,p);
	
	// Affichage du système et de la solution
	fprintf(stdout,"\n\nSYSTÈME DE DÉPART :\n");
	affiche_systeme(&s_ini);		// Système de départ
	
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
	fprintf(stdout,"\nSYSTÈME ÉCHELONNÉ (GAUSS Z/pZ) :\n(p = ");
	mpz_out_str(stdout,10,p);
	fprintf(stdout,")\n");
	affiche_systeme(&s_zpz);
	fprintf(stdout,"\nSOLUTION dans Z/pZ :\n");
	for(int i = 0;i < n;i++){
		mpz_out_str(stdout,10,sol_zpz[i]);
		fprintf(stdout,"\n");
	}
	fprintf(stdout,"\n");
	
	// Vérifications (Bareiss)
	//assert(verif_sol(&s,sol));		// Vérif "triviale"
	assert(verif_sol(&s_ini,sol));		// Vérif sur le systeme de départ
	//ecrit_coeff(&s_ini,0,0,t_m[6]);
	//assert(!verif_sol(&s_ini,sol));		// Autre vérif carrément stupide
	
	// Vérifications (Gauss)
	assert(verif_sol(&s_ini,sol_g));
	
	// TODO : Vérifications pour Gauss dans Z/pZ...
	// ...
	
	
	
	// Suppression des objets utilisés
	for(int i = 0;i < n;i++){
		rat_clear(&sol[i]);
		rat_clear(&sol_g[i]);
		mpz_clear(sol_zpz[i]);
	}
	mpz_clear(k);
	mpz_clear(p);
	detruit_systeme(&s);
	detruit_systeme(&s_ini);
	detruit_systeme(&s_g);
	detruit_systeme(&s_zpz);
	gmp_randclear(state);
	return 0;
}
