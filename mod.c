#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gmp.h>
#include <time.h>
#include "systemes.h"
#include "zpz.h"
#include "mod.h"

// Initialise le syst_zpz sdest, avec la copie du modulo p du systeme ssrc
void init_copie_syst_zpz(syst_zpz* sdest,systeme* ssrc,long int p){
	int n = ssrc->n;
	int m = ssrc->m;
	mpz_t pm,c;
	mpz_init_set_si(pm,p);
	mpz_init(c);
	sdest -> n = n;
	sdest -> m = m;
	sdest -> p = p;
	sdest -> t = malloc(sizeof(long int)*n*m);
	for(int k = 0;k < n*m;k++){
		mpz_fdiv_r(c,ssrc->t[k],pm);
		sdest->t[k] = mpz_get_ui(c);
	}
	mpz_clear(pm);
	mpz_clear(c);
	return;
}

void detruit_syst_zpz(syst_zpz* s){
	free(s->t);
	return;
}

// Renvoie le coefficient (i,j) de s
long int lit_coeff_zpz(syst_zpz* s,int i,int j){
	return s->t[i*(s->m) + j];
}

// Remplace le coefficient (i,j) de s par n
// n est supposé compris entre 0 et (s->p - 1)
void ecrit_coeff_zpz(syst_zpz* s,int i,int j,long int n){
	s->t[i*(s->m) + j] = n;
	return;
}

void affiche_syst_zpz(syst_zpz* s,FILE* f){
	//fprintf(f,"p = %d\n",s->p);	// Inutile
	int n = s->n;
	for(int i = 0;i < n;i++){
		// Affichage d'une ligne de coeffs
		for(int j = 0;j < n;j++){
			fprintf(f,"%ld ",lit_coeff_zpz(s,i,j));
		}
		// Affichage d'un coeff du second membre
		fprintf(f,"  %ld\n",lit_coeff_zpz(s,i,n));
	}
	return;
}

// Vérifie si sol est une solution de s dans Z/pZ
// Les coefficients de s doivent être compris entre 0 et p-1
bool verif_sol_zpz(syst_zpz* s,long int* sol){
	int n = s->n;
	long int p = s->p;
	bool res = true;
	long int som;
	for(int i = 0; i < n;i++){
		som = 0;
		for(int j = 0;j < n;j++){
			som = zpz_add(som, zpz_mul(sol[j], lit_coeff_zpz(s,i,j), p), p);
		}
		res = res && (som == lit_coeff_zpz(s,i,n));
	}
	return res;
}

// Échelonne le système s dans Z/pZ, avec l'algorithme de Gauss
// Les coefficients de s doivent être compris entre 0 et p-1
void zpz_gauss(syst_zpz* s){
	int n = s->n;
	long int p = s->p;
	//int ivp,a,b,c;
	long int ivp,a;
	// Échelonne le système
	for(int k = 0;k < n;k++){
		ivp = zpz_inv(lit_coeff_zpz(s,k,k),p);	// Inverse du pivot (le pivot est supposé non  nul)
		for(int i = k+1;i < n;i++){
			a = zpz_mul(ivp,lit_coeff_zpz(s,i,k),p);
			// (Ligne i) <- (Ligne i) - a * (Ligne k)
			for(int j = k+1;j < n+1;j++){
				//b = zpz_mul(a, lit_coeff_zpz(s,k,j), p);
				//c = zpz_sub(lit_coeff_zpz(s,i,j),b,p);
				//ecrit_coeff_zpz(s,i,j,c);
				// Version sans intermédiaire, à peu près illisible :
				ecrit_coeff_zpz(s,i,j, zpz_sub(lit_coeff_zpz(s,i,j), zpz_mul(a, lit_coeff_zpz(s,k,j) ,p) ,p));
			}
			ecrit_coeff_zpz(s,i,k,0);
		}
	}
	return;
}

// Copie dans sol la solution du système préalablement échelonné s
// Les coefficients de s doivent être compris entre 0 et p-1
void zpz_sol_syst_ech(long int* sol,syst_zpz* s){
	int n = s->n;
	long int p = s->p;
	for(int k = n-1;k >= 0;k--){
		sol[k] = lit_coeff_zpz(s,k,n);
		for(int l = k+1;l < n;l++){
			sol[k] = zpz_sub(sol[k], zpz_mul(sol[l], lit_coeff_zpz(s,k,l), p), p);
		}
		sol[k] = zpz_mul(sol[k], zpz_inv(lit_coeff_zpz(s,k,k),p), p);
	}
	return;
}

// Écrit dans sol la solution du systeme s dans Z/pZ, en travaillant directement sur s
void zpz_resol(syst_zpz* s,long int* sol){
	zpz_gauss(s);
	zpz_sol_syst_ech(sol,s);
	return;
}

// Génére un nombre premier aléatoire d'au plus b bits, donne sa valeur à p, et le renvoie aussi (en tant qu'int)
long int genere_p(mpz_t p,gmp_randstate_t state,mp_bitcnt_t b){
	mpz_t k;
	mpz_init(k);
	mpz_urandomb(k,state,b);
	mpz_nextprime(p,k);
	long int res = mpz_get_ui(p);
	mpz_clear(k);
	return res;
}

// Calcule l'image réciproque de (x1,x2) (dans Z/(n1*n2)Z) par l'isomorphisme des restes chinois
// Écrit le résultat dans res, et nécessite d'avoir déjà établi une relation de Bézout : n1*u + n2*v == 1
// Fonctionne correctement si l'une des opérandes est aussi l'emplacement du résultat
// Calcul plus rapide dans le cas où n1 >> n2
// Le résultat n'est pas nécessairement entre 0 et n1*n2-1, ni entre -n1*n2/2 et n1*n2/2
void chinois_interne(mpz_t res,mpz_t x1,mpz_t x2,mpz_t n1,mpz_t n2,mpz_t u){
	mpz_sub(res,x2,x1);
	mpz_mul(res,res,u);
	mpz_tdiv_r(res,res,n2);
	mpz_mul(res,res,n1);
	mpz_add(res,x1,res);
	
	return;
}

// Appelle la fonction précédente et "convertit" le résultat pour qu'il soit entre 0 et n1*n1-1
void chinois(mpz_t res,mpz_t x1,mpz_t x2,mpz_t n1,mpz_t n2,mpz_t u){
	chinois_interne(res,x1,x2,n1,n2,u);
	if(mpz_cmp_si(res,0) < 0){	// On veut un résultat entre 0 et n1*n2-1 (est-ce vraiment important ?)
		mpz_t n1n2;
		mpz_init(n1n2);
		mpz_mul(n1n2,n1,n2);
		mpz_add(res,res,n1n2);
		mpz_clear(n1n2);
	}
	return;
}

// Établit une relation de Bézout n1*u +n2*v == 1, et appelle la fonction précédente sur les n premières cases des tableaux res, x1 et x2
void chinois_n(int n,mpz_t* res,mpz_t* x1,mpz_t* x2,mpz_t n1,mpz_t n2){
	mpz_t u;
	mpz_init(u);
	mpz_gcdext(NULL,u,NULL,n1,n2);	// NULL car le PGCD ne nous intéresse pas (il vaudra toujours 1), et v non plus (il ne sera pas utilisé)
	for(int i = 0;i < n;i++){
		chinois(res[i],x1[i],x2[i],n1,n2,u);
	}
	mpz_clear(u);
	return;
}

// Anciennes versions de ce qui précède
void chinois_old(mpz_t res,mpz_t x1,mpz_t x2,mpz_t n1,mpz_t n2,mpz_t u,mpz_t v){
	mpz_t t1,t2,n1n2;
	mpz_init_set(t1,x1);
	mpz_init_set(t2,x2);
	mpz_init(n1n2);
	mpz_mul(n1n2,n1,n2);
	mpz_mul(t1,t1,n2);
	mpz_mul(t1,t1,v);
	mpz_mul(t2,t2,n1);
	mpz_mul(t2,t2,u);
	mpz_add(res,t1,t2);
	mpz_fdiv_r(res,res,n1n2);
	mpz_clear(t1);
	mpz_clear(t2);
	mpz_clear(n1n2);
	return;
}
void chinois_old_n(int n,mpz_t* res,mpz_t* x1,mpz_t* x2,mpz_t n1,mpz_t n2){
	mpz_t u,v;
	mpz_init(u);
	mpz_init(v);
	mpz_gcdext(NULL,u,v,n1,n2);	// NULL car le PGCD ne nous intéresse pas (il vaudra toujours 1)
	for(int i = 0;i < n;i++){
		chinois_old(res[i],x1[i],x2[i],n1,n2,u,v);
	}
	mpz_clear(u);
	mpz_clear(v);
	return;
}

// Une autre implémentation de l'algorithme d'Euclide étendu, utilisée pour la reconstruction rationnelle
// S'arrête dès que le reste est plus petit que sqrt(a), et utilise des mpz_t
// Ne calcule que les rk et les vk (pas les uk)
void euclide_etendu_borne(mpz_t r,mpz_t v,mpz_t a,mpz_t b){
	// Initialisation
	mpz_t r_b,r_new;
	mpz_t v_b,v_new;
	mpz_t q,borne;
	mpz_set(r,a);
	mpz_init_set(r_b,b);
	mpz_set_si(v,0);
	mpz_init_set_si(v_b,1);
	mpz_init(r_new);
	mpz_init(v_new);
	mpz_init(q);
	mpz_init(borne);
	mpz_sqrt(borne,a);
	while(mpz_cmp(r_b,borne) > 0){		// Tant que r_b > sqrt(a)
		// Calculs
		mpz_fdiv_q(q,r,r_b);
		mpz_set(r_new,r);
		mpz_submul(r_new,r_b,q);
		mpz_set(v_new,v);
		mpz_submul(v_new,v_b,q);
		// Passage à l'étape suivante
		mpz_set(r,r_b);
		mpz_set(v,v_b);
		mpz_set(r_b,r_new);
		mpz_set(v_b,v_new);
	}
	// ATTENTION : le résultat souhaité est en fait [r_b et v_b] (et non [r et v])
	// Fix douteux :
	mpz_set(r,r_b);
	mpz_set(v,v_b);
	// TODO : réécrire cette fonction convenablement
	
	// DEBUG 
	/*mpz_out_str(stderr,10,r);
	fputc(' ',stderr);
	mpz_out_str(stderr,10,v);
	fputc('\n',stderr);*/
	// TEST OK sur l'exemple du cours (a = 257, b = 42, r = 5, v = -6)
	
	// Suppression des variables utilisées
	mpz_clear(r_b);
	mpz_clear(r_new);
	mpz_clear(v_b);
	mpz_clear(v_new);
	mpz_clear(q);
	mpz_clear(borne);
	return;
}

// Calcule la borne de Hadamard du système
// En utilisant seulement des entiers de GMP, c'est potentiellement donc un peu approximatif...
void hadamard(mpz_t res,systeme* s){
	int n = s->n;
	int m = s->m;
	mpz_t norm2,tmp,coeff;
	mpz_set_si(res,1);
	mpz_init(norm2);
	mpz_init(tmp);
	mpz_init(coeff);
	for(int j = 0;j < m;j++){	// On prend aussi en compte le second membre
		// Calcul de la norme 2 d'une colonne
		mpz_set_si(tmp,0);
		for(int i = 0;i < n;i++){
			lit_coeff(coeff,s,i,j);
			mpz_addmul(tmp,coeff,coeff);
		}
		mpz_sqrt(norm2,tmp);
		mpz_add_ui(norm2,norm2,1);	// Permet d'avoir un résultat >= (au lieu de <=) à la "vraie" borne de Hadamard, tout en gardant seulement des entiers
		mpz_mul(res,res,norm2);
	}
	mpz_clear(norm2);
	mpz_clear(tmp);
	mpz_clear(coeff);
	return;
}

// Méthode modulaire de résolution du système [...]
void modulaire_old(systeme* s,rationnel* sol,gmp_randstate_t state,mp_bitcnt_t b){
	// INITIALISATION
	int j = 1;		// DEBUG, nombre d'itérations effectuées
	clock_t debut,fin;		// DEBUG
	double zpz_resol_tps = 0.;
	double chinois_tps = 0.;
	double reconstr_tps = 0.;
	int n = s -> n;
	syst_zpz sz;	// On va initialiser/copier et détruire ce système à chaque itération
	long int* sol_zpz = malloc(n*sizeof(long int));	// Emplacement de la solution dans Z/pZ que zpz_resol va calculer
	mpz_t* sol_zpz_m = malloc(n*sizeof(mpz_t));	// De même, pour la solution convertie en mpz
	mpz_t* sol_tmp = malloc(n*sizeof(mpz_t));	// Contient ce à partir de quoi on va essayer de reconstruire une solution
	mpz_t* sol_tmp_old = malloc(n*sizeof(mpz_t));	// Contient ce à partir de quoi on a essayé de reconstruire une solution à l'itération précédente
	rationnel* sol_old = malloc(n*sizeof(rationnel));	// Emplacement du candidat de solution de l'itération précédente (nécessaire pour s'arrêter quand on trouve la même solution 2 fois d'affilée)
	mpz_t p_mpz,prod_old,prod,u,v;
	long int p;		// Nombre premier "en cours d'utilisation", version machine
	mpz_init(p_mpz);		// Nombre premier "en cours d'utilisation", version GMP
	mpz_init(prod_old);		// Produit des nombres premiers précédemment utilisés (sans p_mpz)
	mpz_init(prod);			// Produit des nombres premiers précédemment utilisés (avec p_mpz) -> pendant une itération, prod == p*prod_old
	mpz_init(u);		// Pour la relation de Bézout, et pour la reconstruction modulaire de la solution (où il représente r, ce qui en fait un choix de nom douteux)
	mpz_init(v);		// Pour la relation de Bézout, et pour la reconstruction modulaire de la solution
	for(int i = 0;i < n;i++){
		mpz_init(sol_zpz_m[i]);
		mpz_init(sol_tmp[i]);
		mpz_init(sol_tmp_old[i]);
		rat_init(&(sol_old[i]));
	}
	// PREMIÈRE ITÉRATION (On verra plus tard si on peut la faire rentrer dans la boucle principale)
	// "Choix" d'un nombre premier
	p = genere_p(p_mpz,state,b);
	// Initialisation de prod et prod_old
	mpz_set(prod,p_mpz);
	mpz_set_ui(prod_old,1);
	// Calcul d'une solution dans Z/pZ
	init_copie_syst_zpz(&sz,s,p);
	debut = clock();	// DEBUG
	zpz_resol(&sz,sol_zpz);
	fin = clock();	// DEBUG
	zpz_resol_tps += ((double)(fin-debut))/CLOCKS_PER_SEC;	// DEBUG
	// Conversion de la solution en mpz_t, directement dans sol_tmp puisqu'il n'y a pas de chinoiserie à faire
	for(int i = 0;i < n;i++){
		mpz_set_si(sol_tmp[i],sol_zpz[i]);
	}
	// Construction modulaire d'un candidat solution
	debut = clock();	// DEBUG
	for(int i = 0;i < n;i++){
		euclide_etendu_borne(u,v,prod,sol_tmp[i]);	// Ici, u devrait s'appeler r
		if(mpz_cmp_si(v,0) < 0){	// v < 0
			mpz_neg(u,u);
			mpz_neg(v,v);	// De sorte que v soit positif
			rat_set_pq(sol[i],u,v);
		}else{		// v > 0
			rat_set_pq(sol[i],u,v);
		}
	}
	fin = clock();	// DEBUG
	reconstr_tps += ((double)(fin-debut))/CLOCKS_PER_SEC;	// DEBUG
	// MàJ de prod_old (on pourrait en fait l'initialiser ici, puisqu'on ne s'en sert pas avant)
	mpz_set(prod_old,p_mpz);
	// Destruction du syst_zpz, pour pouvoir le réutiliser (peu propre (?))
	detruit_syst_zpz(&sz);
	// BOUCLE PRINCIPALE
	while(!sol_egales(sol,sol_old,n)){
		j++;	// DEBUG
		// "Choix" d'un nombre premier pour cette itération
		p = genere_p(p_mpz,state,b);
		// MàJ de prod
		mpz_mul(prod,prod,p_mpz);
		// Copie du précédent candidat solution sur l'emplacement de l'ancien
		for(int i = 0;i < n;i++){
			rat_set(sol_old[i],sol[i]);
			mpz_set(sol_tmp_old[i],sol_tmp[i]);
		}
		// Calcul d'une solution dans Z/pZ
		init_copie_syst_zpz(&sz,s,p);
		debut = clock();	//DEBUG
		zpz_resol(&sz,sol_zpz);
		fin = clock();	// DEBUG
		zpz_resol_tps += ((double)(fin-debut))/CLOCKS_PER_SEC;	// DEBUG
		// Conversion de la nouvelle solution en mpz_t
		for(int i = 0;i < n;i++){
			mpz_set_si(sol_zpz_m[i],sol_zpz[i]);
		}
		// Restes chinois avec les solutions précédentes
		debut = clock();	// DEBUG
		chinois_n(n,sol_tmp,sol_tmp_old,sol_zpz_m,prod_old,p_mpz);
		// Construction modulaire d'un candidat solution
		fin = clock();	// DEBUG
		chinois_tps += ((double)(fin-debut))/CLOCKS_PER_SEC;	// DEBUG
		debut = clock();	// DEBUG
		for(int i = 0;i < n;i++){
			euclide_etendu_borne(u,v,prod,sol_tmp[i]);	// Ici, u devrait s'appeler r
			if(mpz_cmp_si(v,0) < 0){	// v < 0
				mpz_neg(u,u);
				mpz_neg(v,v);	// De sorte que v soit positif
				rat_set_pq(sol[i],u,v);
			}else{		// v > 0
				rat_set_pq(sol[i],u,v);
			}
		}
		fin = clock();	// DEBUG
		reconstr_tps += ((double)(fin-debut))/CLOCKS_PER_SEC;	// DEBUG
		// MàJ de prod_old
		//mpz_mul(prod_old,prod_old,p_mpz);	// Pas nécessaire de faire ce calcul, on sait qu'on va se retrouver avec prod_old == prod à ce stade
		mpz_set(prod_old,prod);	// Plus rapide
		// Destruction du syst_zpz, pour pouvoir le réutiliser (peu propre (?))
		detruit_syst_zpz(&sz);
	}
	//fprintf(stderr,"Nombre d'itérations : %d\n\n",j);	// DEBUG
	fprintf(stderr,"Nombre d'itérations : %d\nTemps de résolution dans Z/pZ : %lf s\nTemps de restes chinois : %lf s\nTemps de recontruction rationnelle : %lf s\n\n",j,zpz_resol_tps,chinois_tps,reconstr_tps);	// DEBUG
	// SUPPRESSION DES OBJETS UTILISÉS
	for(int i = 0;i < n;i++){
		mpz_clear(sol_zpz_m[i]);
		mpz_clear(sol_tmp[i]);
		mpz_clear(sol_tmp_old[i]);
		rat_clear(&(sol_old[i]));
	}
	mpz_clear(p_mpz);
	mpz_clear(prod_old);
	mpz_clear(prod);
	mpz_clear(u);
	mpz_clear(v);
	free(sol_zpz);
	free(sol_zpz_m);
	free(sol_tmp);
	free(sol_tmp_old);
	free(sol_old);
	return;
}

void modulaire(systeme* s,rationnel* sol,gmp_randstate_t state,mp_bitcnt_t b){
	// INITIALISATION
	int j = 1;		// DEBUG, nombre d'itérations effectuées
	clock_t debut,fin;		// DEBUG
	double zpz_resol_tps = 0.;
	double chinois_tps = 0.;
	double reconstr_tps = 0.;
	long int n = s -> n;
	syst_zpz sz;	// On va initialiser/copier et détruire ce système à chaque itération
	long int* sol_zpz = malloc(n*sizeof(long int));	// Emplacement de la solution dans Z/pZ que zpz_resol va calculer
	mpz_t* sol_zpz_m = malloc(n*sizeof(mpz_t));	// De même, pour la solution convertie en mpz
	mpz_t* sol_tmp = malloc(n*sizeof(mpz_t));	// Contient ce à partir de quoi on va essayer de reconstruire une solution
	mpz_t* sol_tmp_old = malloc(n*sizeof(mpz_t));	// Contient ce à partir de quoi on a essayé de reconstruire une solution à l'itération précédente		(Ne pas essayer de l'enlever, ou bien les restes chinois foireront !)
	mpz_t p_mpz,prod_old,prod,u,v,hada;
	long int p;		// Nombre premier "en cours d'utilisation", version machine
	mpz_init(p_mpz);		// Nombre premier "en cours d'utilisation", version GMP
	mpz_init(prod_old);		// Produit des nombres premiers précédemment utilisés (sans p_mpz)
	mpz_init(prod);			// Produit des nombres premiers précédemment utilisés (avec p_mpz) -> pendant une itération, prod == p*prod_old
	mpz_init(u);		// Pour la relation de Bézout, et pour la reconstruction modulaire de la solution (où il représente r, ce qui en fait un choix de nom douteux)
	mpz_init(v);		// Pour la relation de Bézout, et pour la reconstruction modulaire de la solution
	mpz_init(hada);		// Contiendra (le carré du double de) la borne de Hadamard du système
	for(int i = 0;i < n;i++){
		mpz_init(sol_zpz_m[i]);
		mpz_init(sol_tmp[i]);
		mpz_init(sol_tmp_old[i]);
	}
	// Calcul de la "borne de Hadamard" (en fait 4 fois le carré de la borne de Hadamard)
	hadamard(hada,s);
	mpz_mul(hada,hada,hada);
	mpz_mul_si(hada,hada,4);
	// PREMIÈRE ITÉRATION (On verra plus tard si on peut la faire rentrer dans la boucle principale)
	// "Choix" d'un nombre premier
	p = genere_p(p_mpz,state,b);
	// Initialisation de prod et prod_old
	mpz_set(prod,p_mpz);
	mpz_set_ui(prod_old,1);
	// Calcul d'une solution dans Z/pZ
	init_copie_syst_zpz(&sz,s,p);
	debut = clock();	// DEBUG
	zpz_resol(&sz,sol_zpz);
	fin = clock();	// DEBUG
	zpz_resol_tps += ((double)(fin-debut))/CLOCKS_PER_SEC;	// DEBUG
	// Conversion de la solution en mpz_t, directement dans sol_tmp puisqu'il n'y a pas de chinoiserie à faire
	for(int i = 0;i < n;i++){
		mpz_set_si(sol_tmp[i],sol_zpz[i]);
	}
	// MàJ de prod_old
	mpz_set(prod_old,p_mpz);
	// Destruction du syst_zpz, pour pouvoir le réutiliser
	detruit_syst_zpz(&sz);
	// BOUCLE PRINCIPALE
	while(mpz_cmp(prod,hada) < 0){
		j++;	// DEBUG
		// "Choix" d'un nombre premier pour cette itération
		p = genere_p(p_mpz,state,b);
		// MàJ de prod
		mpz_mul(prod,prod,p_mpz);
		// Copie du précédent candidat solution sur l'emplacement de l'ancien
		for(int i = 0;i < n;i++){
			mpz_set(sol_tmp_old[i],sol_tmp[i]);
		}
		// Calcul d'une solution dans Z/pZ
		init_copie_syst_zpz(&sz,s,p);
		debut = clock();	//DEBUG
		zpz_resol(&sz,sol_zpz);
		fin = clock();	// DEBUG
		zpz_resol_tps += ((double)(fin-debut))/CLOCKS_PER_SEC;	// DEBUG
		// Conversion de la nouvelle solution en mpz_t
		for(int i = 0;i < n;i++){
			mpz_set_si(sol_zpz_m[i],sol_zpz[i]);
		}
		// Restes chinois avec les solutions précédentes
		debut = clock();	// DEBUG
		chinois_n(n,sol_tmp,sol_tmp_old,sol_zpz_m,prod_old,p_mpz);
		fin = clock();	// DEBUG
		chinois_tps += ((double)(fin-debut))/CLOCKS_PER_SEC;	// DEBUG
		// MàJ de prod_old
		mpz_set(prod_old,prod);	// Plus rapide
		// Destruction du syst_zpz, pour pouvoir le réutiliser (peu propre (?))
		detruit_syst_zpz(&sz);
	}
	// Reconstruction de la solution (rationnelle)
	debut = clock();	// DEBUG
	for(int i = 0;i < n;i++){
		euclide_etendu_borne(u,v,prod,sol_tmp[i]);	// Ici, u devrait s'appeler r
		if(mpz_cmp_si(v,0) < 0){	// v < 0
			mpz_neg(u,u);
			mpz_neg(v,v);	// De sorte que v soit positif
			rat_set_pq(sol[i],u,v);
		}else{		// v > 0
			rat_set_pq(sol[i],u,v);
		}
	}
	fin = clock();	// DEBUG
	reconstr_tps += ((double)(fin-debut))/CLOCKS_PER_SEC;	// DEBUG
	//fprintf(stderr,"Nombre d'itérations : %d\n\n",j);	// DEBUG
	fprintf(stderr,"Nombre d'itérations : %d\nTemps de résolution dans Z/pZ : %lf s\nTemps de restes chinois : %lf s\nTemps de recontruction rationnelle : %lf s\n\n",j,zpz_resol_tps,chinois_tps,reconstr_tps);	// DEBUG
	// SUPPRESSION DES OBJETS UTILISÉS
	for(int i = 0;i < n;i++){
		mpz_clear(sol_zpz_m[i]);
		mpz_clear(sol_tmp[i]);
		mpz_clear(sol_tmp_old[i]);
	}
	mpz_clear(p_mpz);
	mpz_clear(prod_old);
	mpz_clear(prod);
	mpz_clear(u);
	mpz_clear(v);
	mpz_clear(hada);
	free(sol_zpz);
	free(sol_zpz_m);
	free(sol_tmp);
	free(sol_tmp_old);
	return;
}
