#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <gmp.h>
#include <pthread.h>
#include <time.h>
#include "rationnels.h"
#include "systemes.h"
#include "bareiss.h"
#include "gauss_sys_rat.h"
#include "io.h"
#include "zpz.h"

void test_zpz_multi(systeme* s,int n,int thr,gmp_randstate_t state){
	// Déclarations
	pthread_t* t = malloc(thr * sizeof(pthread_t));
	mpz_t* sol = malloc(n * thr * sizeof(mpz_t));
	zpz_args* a = malloc(thr * sizeof(zpz_args));
	mpz_t k;
	// Initialisations
	mpz_init(k);
	for(int i = 0;i < thr;i++){
		for(int j = 0;j < n;j++){
			mpz_init(sol[i*n + j]);
		}
		mpz_init(a[i].p);
		mpz_urandomb(k,state,32);
		mpz_nextprime(a[i].p,k);
		a[i].s = s;		// zpz_multi va créer sa propre copie de s, donc on peut donner le même s à chaque thread (s ne sera pas modifié)
		a[i].sol = &sol[i*n];
	}
	// Lancement des calculs en parallèle
	for(int i = 0;i < thr;i++){
		pthread_create(&t[i],NULL,zpz_multi,&a[i]);
	}
	// "Fin" des calculs
	for(int i = 0;i < thr;i++){
		pthread_join(t[i],NULL);
	}	
	// Suppression des objets utilisés
	mpz_clear(k);
	for(int i = 0;i < thr;i++){
		for(int j = 0;j < n;j++){
			mpz_clear(sol[i*n + j]);
		}
		mpz_clear(a[i].p);
	}
	free(t);
	free(sol);
	free(a);
	return;
}

int main(){
	// Initialisation
	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state,time(NULL));
	ecrit_fichier_au_pif("systeme2.txt",20,state,16);
	systeme s,s_g,s_zpz,s_ini,s_zpzm;
	init_lit_systeme(&s,"systeme.txt");
	//init_lit_systeme(&s,"systeme2.txt");
	int n = s.n;		// Nombre de lignes
	int m = s.m;		// Nombre de colonnes, en comptant le second membre (en pratique : m = n+1)
	init_copie_systeme(&s_ini,&s);	// Copie qui servira à conserver le système initial, pour pouvoir tester notre solution à la fin
	init_copie_systeme(&s_g,&s);	// Copie qui sera utilisée par l'algo de Gauss (pas vraiment nécessaire, mais c'est plus propre que de lui donner s_ini)
	init_copie_systeme(&s_zpz,&s);	// Copie qui sera utilisée pour tester l'algo de Gauss dans Z/pZ
	init_copie_systeme(&s_zpzm,&s);	// Copie qui sera utilisée pour tester l'algo de Gauss dans plusieurs Z/pZ à la fois
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
	
	// Affichage du système de départ
	fprintf(stdout,"\n\nSYSTÈME DE DÉPART :\n");
	affiche_systeme(&s_ini);		// Système de départ
	
	// Calcul d'une solution, avec l'algo de Bareiss
	// Prend un temps raisonnable
	/*bareiss(&s);
	sol_syst_echelonne(&s,sol);
	fprintf(stdout,"\n\nSYSTÈME ÉCHELONNÉ (BAREISS) :\n");
	affiche_systeme(&s);		// Système échelonné
	fprintf(stdout,"\n\nSOLUTION (BAREISS) :\n");
	for(int i = 0;i < n;i++){
		rat_aff(sol[i]);
		fprintf(stdout,"\n");
	}*/
	
	// Calcul d'une solution, avec l'algo de Gauss
	// Prend 3 plombes
	/*gauss(&s_g,sol_g);
	fprintf(stdout,"\nSOLUTION (GAUSS) :\n");
	for(int i = 0;i < n;i++){
		rat_aff(sol_g[i]);
		fprintf(stdout,"\n");
	}*/
	
	// Calcul d'une solution dans Z/pZ (avec p choisi "au hasard" avec à peu près 32 bits)
	// Va suspicieusement vite
	/*mpz_urandomb(k,state,30);
	mpz_nextprime(p,k);
	//mpz_set_si(p,17);	// TEMPORAIRE
	zpz_sans_copie(sol_zpz,&s_zpz,p);
	fprintf(stdout,"\nSYSTÈME ÉCHELONNÉ (GAUSS Z/pZ) :\n(p = ");
	mpz_out_str(stdout,10,p);
	fprintf(stdout,")\n");
	affiche_systeme(&s_zpz);
	fprintf(stdout,"\nSOLUTION dans Z/pZ :\n");
	for(int i = 0;i < n;i++){
		mpz_out_str(stdout,10,sol_zpz[i]);
		fprintf(stdout,"\n");
	}
	fprintf(stdout,"\n");*/
	
	// Vérifications (Bareiss)
	//assert(verif_sol(&s,sol));		// Vérif "triviale"
	//assert(verif_sol(&s_ini,sol));		// Vérif sur le systeme de départ
	
	// Vérifications (Gauss)
	//assert(verif_sol(&s_ini,sol_g));
	
	// TODO : Vérifications pour Gauss dans Z/pZ...
	//assert(verif_sol_zpz(&s_zpz,sol_zpz,p));
	// ATTENTION, TEST CI-DESSUS PAS OK
	
	// Essai d'exécution de Gauss sur des Z/pZ en parallèle
	test_zpz_multi(&s_zpzm,n,8,state);
	
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
	detruit_systeme(&s_zpzm);
	gmp_randclear(state);
	return 0;
}
