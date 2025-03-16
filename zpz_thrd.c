#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include <pthread.h>
#include "zpz.h"
#include "zpz_thrd.h"

// Résout a->s dans Z/(a->s->p)Z, et écrit le résultat en int dans a->sol et le résultat en mpz_t dans a->sol_m
// Est fait pour être appelé par pthread_create
void* zpz_thrd(void* a_){
	// Conversion des arguments
	zpz_thrd_args* a = (zpz_thrd_args*)a_;
	int n = a->s->n;
	// Calculs
	zpz_resol(a->s,a->sol);
	// Conversion de la solution en mpz_t
	for(int i = 0;i < n;i++){
		mpz_set_si(a->sol_m[i],a->sol[i]);
	}
	return NULL;
}

// Appelle chinois_n, et est fait pour être appelé par pthread_create
/*void* chinois_n_thrd(void* a_){
	// Conversion des arguments
	chinois_n_thrd_args* a = (chinois_n_thrd_args*) a_;
	chinois_n(a->n,a->res,a->x1,a->x2,a->n1,a->n2);
	return NULL;
}*/

void modulaire_thrd(systeme* s,rationnel* sol,gmp_randstate_t state,mp_bitcnt_t b,int thr){
	// DÉCLARATIONS
	int n = s -> n;
	int j = 1; // DEBUG, nombre d'itérations effectuées
	syst_zpz* sz = malloc(sizeof(syst_zpz)*thr);	// On va initialiser/copier et détruire ce système à chaque itération
	int* p = malloc(sizeof(int)*thr);
	mpz_t* p_mpz = malloc(sizeof(mpz_t)*thr);
	mpz_t prod;		// Produit de tous les nombres premiers utilisés jusqu'à présent (y compris ceux de l'étape en cours)
	mpz_t prod_old;	// Produit de tous les nombres premiers utilisés jusqu'à présent (sans ceux de l'étape en cours)
	mpz_t prod_act;	// Produit des nombres premiers utilisés à l'étape en cours (sans ceux d'avant)
	mpz_t r,v;
	int* sol_zpz = malloc(n*sizeof(int)*thr);	// Emplacement des solutions dans Z/pZ que zpz_resol va calculer
	mpz_t* sol_zpz_m = malloc(n*sizeof(mpz_t)*thr);	// De même, pour les solutions converties en mpz_t
	//mpz_t* sol_act = malloc(n*sizeof(mpz_t));	// Emplacement de la solution dans Z/prod_actZ, construite à partir des thr solutions de sol_zpz_m
	mpz_t* sol_tmp = malloc(n*sizeof(mpz_t));	// Contient ce à partir de quoi on va essayer de reconstruire une solution (et que l'on a construit à partir de sol_zpz_m[0] et sol_tmp_old)
	mpz_t* sol_tmp_old = malloc(n*sizeof(mpz_t));	// Contient ce à partir de quoi on a essayé de reconstruire une solution à l'itération précédente
	rationnel* sol_old = malloc(n*sizeof(rationnel));	// Emplacement du candidat de solution de l'itération précédente (nécessaire pour s'arrêter quand on trouve la même solution 2 fois d'affilée)
	pthread_t* thrd = malloc(sizeof(pthread_t)*thr);	// Threads
	zpz_thrd_args* a = malloc(sizeof(zpz_thrd_args)*thr);	// Arguments de zpz_thrd
	
	// INITIALISATIONS
	mpz_init_set_si(prod,1);
	mpz_init_set_si(prod_old,1);
	mpz_init(prod_act);
	mpz_init(r);
	mpz_init(v);
	for(int t = 0;t < thr;t++){
		mpz_init(p_mpz[t]);
	}
	for(int i = 0;i < n;i++){
		//mpz_init(sol_act[i]);
		mpz_init(sol_tmp[i]);
		mpz_init(sol_tmp_old[i]);
		rat_init(&(sol_old[i]));
	}
	for(int k = 0;k < n*thr;k++){
		mpz_init(sol_zpz_m[k]);
	}
	// PREMIÈRE ITÉRATION
	// Initialisation de prod_act
	mpz_set_si(prod_act,1);
	// "Choix" des nombres premiers pour cette itération, MàJ de prod, calcul de prod_act, et création des syst_zpz
	// Ceci n'a pas été mis en parallèle, car de précédents tests semblaient indiquer qu'on ne pouvait pas vraiment lire dans s en parallèle -> À RE-TESTER...
	for(int t = 0;t < thr;t++){
		p[t] = genere_p(p_mpz[t],state,b);
		mpz_mul(prod,prod,p_mpz[t]);
		//mpz_mul(prod_act,p_mpz[t]);	// PLUS TARD
		init_copie_syst_zpz(&(sz[t]),s,p[t]);
	}
	// Initialisation des zpz_args
	for(int t = 0;t < thr;t++){
		a[t].s = &(sz[t]);
		a[t].sol = &(sol_zpz[t*n]);
		a[t].sol_m = &(sol_zpz_m[t*n]);
	}
	// Calcul des solutions en parallèle
	for(int t = 0;t < thr;t++){
		pthread_create(&thrd[t],NULL,zpz_thrd,&a[t]);
	}
	// Fin des calculs
	for(int t = 0;t < thr;t++){
		pthread_join(thrd[t],NULL);
	}
	// Par restes chinois, on bricole une solution dans Z/prod_actZ
	// On écrit directement le résultat dans sol_tmp, puisqu'il n'y aura aucune autre chinoiserie à effectuer
	// VERSION TEMPORAIRE, ON PEUT FAIRE ÇA EN PARALLÈLE
	mpz_set(prod_act,p_mpz[0]);
	for(int t = 1;t < thr;t++){
		chinois_n(n,sol_tmp,&(sol_zpz_m[t]),sol_tmp,p_mpz[t],prod_act);
		mpz_mul(prod_act,prod_act,p_mpz[t]);
	}
	// Construction modulaire d'un candidat solution (pareil que dans zpz.c)
	for(int i = 0;i < n;i++){
		euclide_etendu_borne(r,v,prod,sol_tmp[i]);
		if(mpz_cmp_si(v,0) < 0){	// v < 0
			mpz_neg(r,r);
			mpz_neg(v,v);	// De sorte que v soit positif
			rat_set_pq(sol[i],r,v);
		}else{		// v > 0
			rat_set_pq(sol[i],r,v);
		}
	}
	// MàJ de prod_old
	mpz_set(prod_old,prod);
	// Destruction des syst_zpz, pour pouvoir les réutiliser avec d'autres valeurs de p
	for(int t = 0;t < thr;t++){
		detruit_syst_zpz(&(sz[t]));
	}
	// BOUCLE PRINCIPALE			(QUELQUE CHOSE NE VA PAS LÀ-DEDANS)
	//while(!sol_egales(sol,sol_old,n)){
	while((!sol_egales(sol,sol_old,n)) && (j < 50)){	// 2ème condition temporaire...
		j++;	// DEBUG
		// Copie du précédent candidat solution sur l'emplacement de l'ancien
		for(int i = 0;i < n;i++){
			rat_set(sol_old[i],sol[i]);
			mpz_set(sol_tmp_old[i],sol_tmp[i]);
		}
		// Réinitialisation de prod_act
		mpz_set_si(prod_act,1);
		// "Choix" des nombres premiers pour cette itération, MàJ de prod, calcul de prod_act, et création des syst_zpz
		// Ceci n'a pas été mis en parallèle, car de précédents tests semblaient indiquer qu'on ne pouvait pas vraiment lire dans s en parallèle -> À RE-TESTER...
		for(int t = 0;t < thr;t++){
			p[t] = genere_p(p_mpz[t],state,b);
			mpz_mul(prod,prod,p_mpz[t]);
			//mpz_mul(prod_act,p_mpz[t]);	// PLUS TARD
			init_copie_syst_zpz(&(sz[t]),s,p[t]);
		}
		// Initialisation des zpz_args
		for(int t = 0;t < thr;t++){
			a[t].s = &(sz[t]);
			a[t].sol = &(sol_zpz[t*n]);
			a[t].sol_m = &(sol_zpz_m[t*n]);
		}
		// Calcul des solutions en parallèle
		for(int t = 0;t < thr;t++){
			pthread_create(&thrd[t],NULL,zpz_thrd,&a[t]);
		}
		// Fin des calculs
		for(int t = 0;t < thr;t++){
			pthread_join(thrd[t],NULL);
		}
		// Par restes chinois, on bricole une solution dans Z/prod_actZ
		// On écrit le résultat dans sol_zpz_m[0]
		// VERSION TEMPORAIRE, ON PEUT FAIRE ÇA EN PARALLÈLE
		mpz_set(prod_act,p_mpz[0]);
		for(int t = 1;t < thr;t++){
			chinois_n(n,&(sol_zpz_m[0]),&(sol_zpz_m[t]),&(sol_zpz_m[0]),p_mpz[t],prod_act);
			mpz_mul(prod_act,prod_act,p_mpz[t]);
		}
		// (prod_old * prod_act == prod)
		// On combine la solution dans Z/prod_actZ (calculée juste avant) avec la solution dans Z/prod_oldZ
		chinois_n(n,sol_tmp,&(sol_zpz_m[0]),sol_tmp_old,prod_act,prod_old);
		// Construction modulaire d'un candidat solution (pareil que dans zpz.c)
		for(int i = 0;i < n;i++){
			euclide_etendu_borne(r,v,prod,sol_tmp[i]);
			if(mpz_cmp_si(v,0) < 0){	// v < 0
				mpz_neg(r,r);
				mpz_neg(v,v);	// De sorte que v soit positif
				rat_set_pq(sol[i],r,v);
			}else{		// v > 0
				rat_set_pq(sol[i],r,v);
			}
		}
		// MàJ de prod_old
		mpz_set(prod_old,prod);
		// Destruction des syst_zpz, pour pouvoir les réutiliser avec d'autres valeurs de p
		for(int t = 0;t < thr;t++){
			detruit_syst_zpz(&(sz[t]));
		}
	}
	fprintf(stderr,"Nombre d'itérations : %d * %d\n\n",j,thr);	// DEBUG
	// SUPPRESSIONS
	mpz_clear(prod);
	mpz_clear(prod_old);
	mpz_clear(prod_act);
	mpz_clear(r);
	mpz_clear(v);
	for(int t = 0;t < thr;t++){
		mpz_clear(p_mpz[t]);
	}
	for(int i = 0;i < n;i++){
		//mpz_clear(sol_act[i]);
		mpz_clear(sol_tmp[i]);
		mpz_clear(sol_tmp_old[i]);
		rat_clear(&(sol_old[i]));
	}
	for(int k = 0;k < n*thr;k++){
		mpz_clear(sol_zpz_m[k]);
	}
	free(sz);
	free(p);
	free(p_mpz);
	free(sol_zpz);
	free(sol_zpz_m);
	free(sol_tmp);
	free(sol_tmp_old);
	free(sol_old);
	free(thrd);
	free(a);
	return;
}





// ANCIEN TEST
/*void test_zpz_multi(systeme* s,int n,int thr,gmp_randstate_t state){
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
}*/
