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
		//mpz_out_str(stderr,10,a->sol_m[i]);		// DEBUG
		//fputc('\n',stderr);				// DEBUG
	}
	return NULL;
}

// Appelle chinois_n, et est fait pour être appelé par pthread_create
void* chinois_n_thrd(void* a_){
	// Conversion des arguments
	chinois_n_thrd_args* a = (chinois_n_thrd_args*) a_;
	chinois_n(a->n,a->res,a->x1,a->x2,a->n1,a->n2);
	return NULL;
}

void modulaire_thrd(systeme* s,rationnel* sol,gmp_randstate_t state,mp_bitcnt_t bit,int thr){
	// DÉCLARATIONS
	int n = s -> n;
	int j = 1; // DEBUG, nombre d'itérations effectuées
	syst_zpz* sz = malloc(sizeof(syst_zpz)*thr);	// On va initialiser/copier et détruire ce système à chaque itération
	int* p = malloc(sizeof(int)*thr);
	mpz_t* p_mpz = malloc(sizeof(mpz_t)*thr);
	mpz_t prod;		// Produit de tous les nombres premiers utilisés jusqu'à présent (y compris ceux de l'étape en cours)
	mpz_t prod_old;	// Produit de tous les nombres premiers utilisés jusqu'à présent (sans ceux de l'étape en cours)
	mpz_t prod_act;	// Produit des nombres premiers utilisés à l'étape en cours (sans ceux d'avant)
	mpz_t r,v,hada;
	int* sol_zpz = malloc(n*sizeof(int)*thr);	// Emplacement des solutions dans Z/pZ que zpz_resol va calculer
	mpz_t* sol_zpz_m = malloc(n*sizeof(mpz_t)*thr);	// De même, pour les solutions converties en mpz_t
	//mpz_t* sol_act = malloc(n*sizeof(mpz_t));	// Emplacement de la solution dans Z/prod_actZ, construite à partir des thr solutions de sol_zpz_m
	mpz_t* sol_tmp = malloc(n*sizeof(mpz_t));	// Contient ce à partir de quoi on va essayer de reconstruire une solution (et que l'on a construit à partir de sol_zpz_m[0] et sol_tmp_old)
	mpz_t* sol_tmp_old = malloc(n*sizeof(mpz_t));	// Contient ce à partir de quoi on a essayé de reconstruire une solution à l'itération précédente
	//rationnel* sol_old = malloc(n*sizeof(rationnel));	// Emplacement du candidat de solution de l'itération précédente (nécessaire pour s'arrêter quand on trouve la même solution 2 fois d'affilée)
	pthread_t* thrd = malloc(sizeof(pthread_t)*thr);	// Threads
	zpz_thrd_args* a = malloc(sizeof(zpz_thrd_args)*thr);	// Arguments de zpz_thrd
	chinois_n_thrd_args* b = malloc(sizeof(chinois_n_thrd_args)*thr);	// Arguments de chinois_n_thrd_args (au plus la moitié serviront...)
	
	// INITIALISATIONS
	mpz_init_set_si(prod,1);
	mpz_init_set_si(prod_old,1);
	mpz_init(prod_act);
	mpz_init(r);
	mpz_init(v);
	mpz_init(hada);
	for(int t = 0;t < thr;t++){
		mpz_init(p_mpz[t]);
		mpz_init(b[t].n1);
		mpz_init(b[t].n2);
	}
	for(int i = 0;i < n;i++){
		//mpz_init(sol_act[i]);
		mpz_init(sol_tmp[i]);
		mpz_init(sol_tmp_old[i]);
		//rat_init(&(sol_old[i]));
	}
	for(int k = 0;k < n*thr;k++){
		mpz_init(sol_zpz_m[k]);
	}
	// Calcul de la "borne de Hadamard" (en fait 4 fois le carré de la borne de Hadamard)
	hadamard(hada,s);
	mpz_mul(hada,hada,hada);
	mpz_mul_si(hada,hada,4);
	// PREMIÈRE ITÉRATION (On verra plus tard si on peut la faire rentrer dans la boucle principale)
	// Initialisation de prod_act
	mpz_set_si(prod_act,1);
	// Ce qui suit n'a pas été mis en parallèle, car de précédents tests semblaient indiquer qu'on ne pouvait pas vraiment lire dans s en parallèle -> À RE-TESTER...
	for(int t = 0;t < thr;t++){
		// "Choix" d'un nombre premier
		p[t] = genere_p(p_mpz[t],state,bit);
		// MàJ de prod
		mpz_mul(prod,prod,p_mpz[t]);
		// MàJ de prod_act
		//mpz_mul(prod_act,prod_act,p_mpz[t]);	// PLUS TARD
		// Initialisation des zpz_args
		init_copie_syst_zpz(&(sz[t]),s,p[t]);
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
	
	
	
	// Copie du premier résultat dans sol_tmp (emplacement du futur résultat des restes chinois
	for(int i = 0;i < n;i++){
		mpz_set(sol_tmp[i],sol_zpz_m[i]);
	}
	mpz_set(prod_act,p_mpz[0]);
	// Restes chinois pour obtenir une solution dans Z/prod_actZ
	// Pour la première itération, il n'y a pas d'autres restes chinois à utiliser, et on peut donc mettre directement le résultat dans sol_tmp
	for(int t = 1;t < thr;t++){
		chinois_n(n,sol_tmp,&(sol_zpz_m[t*n]),sol_tmp,p_mpz[t],prod_act);
		mpz_mul(prod_act,prod_act,p_mpz[t]);
	}
	
	// Alternative expérimentale et cassée (restes chinois en parallèle, fera gagner un temps précieux)
	/*int pas = 1;
	// "Pré-initialisation" des arguments
	for(int t = 0;t < thr;t++){
		mpz_set(b[t].n1,p_mpz[t]);
	}
	fprintf(stderr,"Hola\n");
	while(pas/2 < thr){
		pas = pas*2;
		for(int t = 0;t < thr;t += pas){
			b[t].n = n;
			b[t].res = &sol_zpz_m[t*n];
			b[t].x1 = &sol_zpz_m[t*n];
			b[t].x2 = &sol_zpz_m[t*n + (n*pas/2)];
			// b[t].n1 est à ce stade déjà correct		// TODO : le calculer ici ?
			mpz_set(b[t].n2,b[t+pas/2].n1);
		}
		// Lancement des calculs
		for(int t = 0;t < thr;t += pas){
			pthread_create(&thrd[t],NULL,chinois_n_thrd,&b[t]);
		}
		// Fin des calculs
		for(int t = 0;t < thr;t += pas){
			pthread_join(thrd[t],NULL);
		}
		// MàJ des n1 pour la prochaine étape (on en fait 2 fois trop, mais c'est plus simple pour l'instant
		for(int t = 0;t < thr;t += pas){
			mpz_mul(b[t].n1,b[t].n1,b[t].n2);
		}
		fprintf(stderr,"%d\n",pas);
	}
	// Copie du résultat des restes chinois dabs sol_tmp
	for(int i = 0;i < n;i++){
		mpz_set(sol_tmp[i],sol_zpz_m[i]);
	}*/
	
	
	// Construction d'un candidat solution rationnel		// NON
	/*for(int i = 0;i < n;i++){
		euclide_etendu_borne(r,v,prod,sol_tmp[i]);
		if(mpz_cmp_si(v,0) < 0){	// v < 0
			mpz_neg(r,r);
			mpz_neg(v,v);	// De sorte que v soit positif
			rat_set_pq(sol[i],r,v);
		}else{		// v > 0
			rat_set_pq(sol[i],r,v);
		}
	}*/
	// MàJ de prod_old (on pourrait en fait l'initialiser ici, puisqu'on ne s'en sert pas avant)
	mpz_set(prod_old,prod);
	// Destruction des syst_zpz, pour pouvoir les réutiliser
	for(int t = 0;t < thr;t++){
		detruit_syst_zpz(&(sz[t]));
	}
	// BOUCLE PRINCIPALE
	//while(!sol_egales(sol,sol_old,n)){
	while(mpz_cmp(prod,hada) < 0){
		j++;	// DEBUG
		// Copie du précédent candidat solution sur l'emplacement de l'ancien
		for(int i = 0;i < n;i++){
			//rat_set(sol_old[i],sol[i]);
			mpz_set(sol_tmp_old[i],sol_tmp[i]);
		}
		for(int t = 0;t < thr;t++){
			// "Choix" d'un nombre premier pour cette itération
			p[t] = genere_p(p_mpz[t],state,bit);
			// MàJ de prod
			mpz_mul(prod,prod,p_mpz[t]);
			// MàJ de prod_act
			//mpz_mul(prod_act,prod_act,p_mpz[t]);	// PLUS TARD
			// Initialisation des zpz_args
			init_copie_syst_zpz(&(sz[t]),s,p[t]);
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
		// Réinitialisation de prod act
		mpz_set(prod_act,p_mpz[0]);
		// Restes chinois pour obtenir une solution dans Z/prod_actZ
		// Le résultat est écrit dans sol_zpz_m (cases 0 à n-1) (ie sol_zpz_m[0] vu comme un tableau de taille n, bref)
		for(int t = 1;t < thr;t++){
			chinois_n(n,sol_zpz_m,&(sol_zpz_m[t*n]),sol_zpz_m,p_mpz[t],prod_act);
			mpz_mul(prod_act,prod_act,p_mpz[t]);
		}
		// Restes chinois avec les solutions précédentes
		chinois_n(n,sol_tmp,sol_zpz_m,sol_tmp_old,prod_act,prod_old);
		// Construction modulaire d'un candidat solution		// NON
		/*for(int i = 0;i < n;i++){
			euclide_etendu_borne(r,v,prod,sol_tmp[i]);
			if(mpz_cmp_si(v,0) < 0){	// v < 0
				mpz_neg(r,r);
				mpz_neg(v,v);	// De sorte que v soit positif
				rat_set_pq(sol[i],r,v);
			}else{		// v > 0
				rat_set_pq(sol[i],r,v);
			}
		}*/
		// MàJ de prod_old
		mpz_set(prod_old,prod);
		// Destruction des syst_zpz, pour pouvoir les réutiliser
		for(int t = 0;t < thr;t++){
			detruit_syst_zpz(&(sz[t]));
		}
	}
	// Construction de la solution rationelle
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
	fprintf(stderr,"Nombre d'itérations : %d * %d\n\n",j,thr);	// DEBUG
	// SUPPRESSIONS
	mpz_clear(prod);
	mpz_clear(prod_old);
	mpz_clear(prod_act);
	mpz_clear(r);
	mpz_clear(v);
	mpz_clear(hada);
	for(int t = 0;t < thr;t++){
		mpz_clear(p_mpz[t]);
		mpz_clear(b[t].n1);
		mpz_clear(b[t].n2);
	}
	for(int i = 0;i < n;i++){
		//mpz_clear(sol_act[i]);
		mpz_clear(sol_tmp[i]);
		mpz_clear(sol_tmp_old[i]);
		//rat_clear(&(sol_old[i]));
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
	//free(sol_old);
	free(thrd);
	free(a);
	free(b);
	return;
}


