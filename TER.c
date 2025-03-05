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
	int* sol = malloc(n * thr * sizeof(int));
	pthread_t* t = malloc(thr * sizeof(pthread_t));
	zpz_args* a = malloc(thr * sizeof(zpz_args));
	mpz_t k,p_mpz;
	int p;
	// Initialisations
	mpz_init(p_mpz);
	for(int i = 0;i < thr;i++){
		p = genere_p(p_mpz,state,30);
		fprintf(stderr,"%d\n",p);	// Juste pour voir si quelque chose s'exécute (et, tant qu'à faire, s'exécute avec des valeurs de p différentes)
		init_copie_syst_zpz(&(a[i].s),s,p);
		a[i].sol = &sol[i*n];
	}
	// Lancement des calculs en parallèle
	for(int i = 0;i < thr;i++){
		pthread_create(&t[i],NULL,zpz_resol_thrd,&a[i]);
	}
	// "Fin" des calculs
	for(int i = 0;i < thr;i++){
		pthread_join(t[i],NULL);
	}
	// Suppression des objets utilisés
	mpz_clear(p_mpz);
	free(sol);
	free(t);
	free(a);
	return;
}

int main(){
	// Initialisation
	FILE* f = stdout;
	gmp_randstate_t state;
	gmp_randinit_default(state);
	gmp_randseed_ui(state,time(NULL));
	ecrit_fichier_au_pif("systeme2.txt",25,state,512);		// Les résultats dépassent du terminal
	//ecrit_fichier_au_pif("systeme2.txt",8,state,64);
	systeme s,s_ini;
	syst_zpz s_zpz,s_zpzv;
	//init_lit_systeme(&s,"systeme.txt");
	init_lit_systeme(&s,"systeme2.txt");
	int n = s.n;		// Nombre de lignes
	int m = s.m;		// Nombre de colonnes, en comptant le second membre (en pratique : m = n+1)
	init_copie_systeme(&s_ini,&s);	// Copie qui servira à conserver le système initial, pour pouvoir tester notre solution à la fin (sera aussi donnée à l'algo de Gauss sur les rationnels et à zpz_thrd, car ils ne le modifieront pas)
	rationnel* sol_b = malloc(n*sizeof(rationnel));	// Contiendra la solution donnée par l'algorithme de Bareiss
	rationnel* sol_g = malloc(n*sizeof(rationnel));	// Contiendra la solution donnée par l'algorithme de Gauss
	int* sol_zpz = malloc(n*sizeof(int));		// Contiendra la solution donnée par l'algorithme de Gauss dans Z/pZ
	for(int i = 0;i < n;i++){
		rat_init(&sol_b[i]);
		rat_init(&sol_g[i]);
	}
	mpz_t p_mpz;	// On est pour l'instant obligé de déclarer cette variable, même si on ne s'en sert pas directement
	mpz_init(p_mpz);
	int p = genere_p(p_mpz,state,30);	// Est-il possible de faire ça plus simplement avec randint() ou d'autres fonctions standard ?
	init_copie_syst_zpz(&s_zpz,&s_ini,p);	// Copie du système sur laquelle l'aglo de Gauss dans Z/pZ sera testé
	init_copie_syst_zpz(&s_zpzv,&s_ini,p);	// Copie du système qui servira à conserver le résultat initial modulo p (pour tester le résultat)
	
	// Affichage du système de départ
	fprintf(f,"\n\nSYSTÈME DE DÉPART :\n");
	affiche_systeme(&s_ini,f);		// Système de départ
	
	// Calcul d'une solution, avec l'algo de Bareiss
	// Prend un temps raisonnable
	/*bareiss(&s);
	sol_syst_echelonne(&s,sol_b);
	fprintf(f,"\n\nSYSTÈME ÉCHELONNÉ (BAREISS) :\n");
	affiche_systeme(&s,f);		// Système échelonné
	fprintf(f,"\n\nSOLUTION (BAREISS) :\n");
	for(int i = 0;i < n;i++){
		rat_aff(sol_b[i],f);
		fprintf(f,"\n");
	}
	
	// Calcul d'une solution, avec l'algo de Gauss
	// Prend 3 plombes
	gauss(&s_ini,sol_g);
	fprintf(f,"\nSOLUTION (GAUSS) :\n");
	for(int i = 0;i < n;i++){
		rat_aff(sol_g[i],f);
		fprintf(f,"\n");
	}*/
	
	// Calcul d'une solution dans Z/pZ (avec p choisi "au hasard" avec à peu près 32 bits)
	// Va suspicieusement vite
	fprintf(f,"\nP = %d\n\nSYSTÈME DE DÉPART MODULO P :\n",p);
	affiche_syst_zpz(&s_zpzv,f);
	zpz_resol(&s_zpz,sol_zpz);
	fprintf(f,"\nSYSTÈME ÉCHELONNÉ DANS Z/pZ :\n");
	affiche_syst_zpz(&s_zpz,f);
	fprintf(f,"\nSOLUTION dans Z/pZ :\n");
	for(int i = 0;i < n;i++){
		fprintf(f,"%d\n",sol_zpz[i]);
	}
	fputc('\n',f);
	
	// Vérifications (Bareiss)
	//assert(verif_sol(&s,sol_b));		// Vérif "triviale"
	//assert(verif_sol(&s_ini,sol_b));		// Vérif sur le systeme de départ
	
	// Vérifications (Gauss)
	//assert(verif_sol(&s_ini,sol_g));	// Vérif sur le système de départ
	
	// Vérifications (Gauss dans Z/pZ)
	assert(verif_sol_zpz(&s_zpz,sol_zpz));	// Vérif "triviale"
	assert(verif_sol_zpz(&s_zpzv,sol_zpz));	// Vérif sur le système de départ
	
	// Essai d'exécution de Gauss sur des Z/pZ en parallèle
	//test_zpz_multi(&s_ini,n,8,state);
	// TEST OK
	
	// Suppression des objets utilisés
	for(int i = 0;i < n;i++){
		rat_clear(&sol_b[i]);
		rat_clear(&sol_g[i]);
	}
	free(sol_b);
	free(sol_g);
	free(sol_zpz);
	mpz_clear(p_mpz);
	detruit_systeme(&s);
	detruit_systeme(&s_ini);
	detruit_syst_zpz(&s_zpz);
	detruit_syst_zpz(&s_zpzv);
	gmp_randclear(state);
	return 0;
}
