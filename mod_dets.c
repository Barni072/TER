#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <gmp.h>
#include "systemes.h"
#include "zpz.h"
#include "mod.h"
#include "mod_dets.h"

// Pour représenter les systèmes à coefficients dans Z/pZ en représentation symétrique, on utilise la même structure que dans zpz.h

// L'initialisation doit être légèrement modifiée
void init_copie_syst_zpzs(syst_zpzs* sdest,systeme* ssrc,int p){
	int n = ssrc->n;
	int m = ssrc->m;
	mpz_t pm,c;
	mpz_init_set_si(pm,p);
	mpz_init(c);
	sdest -> n = n;
	sdest -> m = m;
	sdest -> p = p;
	sdest -> t = malloc(sizeof(int)*n*m);
	for(int k = 0;k < n*m;k++){
		mpz_fdiv_r(c,ssrc->t[k],pm);
		sdest->t[k] = mpz_get_ui(c);
		// "Symétrisation" :
		if(sdest->t[k] > p/2) sdest->t[k] -= p;
	}
	mpz_clear(pm);
	mpz_clear(c);
	return;
}

// Les 4 fonctions suivantes sont simplement des "alias" -> est-il bien utile de les garder ?
void detruit_syst_zpzs(syst_zpzs* s){
	detruit_syst_zpz(s);
	return;
}
int lit_coeff_zpzs(syst_zpzs* s,int i,int j){
	return lit_coeff_zpz(s,i,j);
}
void ecrit_coeff_zpzs(syst_zpzs* s,int i,int j,int n){
	ecrit_coeff_zpz(s,i,j,n);
	return;
}
void affiche_syst_zpzs(syst_zpzs* s,FILE* f){
	affiche_syst_zpz(s,f);
	return;
}



bool verif_sol_zpzs(syst_zpzs* s,int* sol){
	int n = s->n;
	int p = s->p;
	bool res = true;
	int som;
	for(int i = 0; i < n;i++){
		som = 0;
		for(int j = 0;j < n;j++){
			som = zpzs_add(som, zpzs_mul(sol[j], lit_coeff_zpzs(s,i,j), p), p);
		}
		res = res && (som == lit_coeff_zpzs(s,i,n));
	}
	return res;
}

// Échelonne le système avec l'algorithme de Gauss
void zpzs_gauss(syst_zpzs* s){
	int n = s->n;
	int p = s->p;
	int ivp,a;
	// Échelonne le système
	for(int k = 0;k < n;k++){
		ivp = zpzs_inv(lit_coeff_zpzs(s,k,k),p);	// Inverse du pivot (le pivot est supposé non  nul)
		for(int i = k+1;i < n;i++){
			a = zpzs_mul(ivp,lit_coeff_zpzs(s,i,k),p);
			// (Ligne i) <- (Ligne i) - a * (Ligne k)
			for(int j = k+1;j < n+1;j++){
				ecrit_coeff_zpzs(s,i,j, zpzs_sub(lit_coeff_zpzs(s,i,j), zpzs_mul(a, lit_coeff_zpzs(s,k,j) ,p) ,p));
			}
			ecrit_coeff_zpzs(s,i,k,0);
		}
	}
	return;
}

void zpzs_sol_syst_ech(int* sol,syst_zpzs* s){
	int n = s->n;
	int p = s->p;
	for(int k = n-1;k >= 0;k--){
		sol[k] = lit_coeff_zpzs(s,k,n);
		for(int l = k+1;l < n;l++){
			sol[k] = zpzs_sub(sol[k], zpzs_mul(sol[l], lit_coeff_zpzs(s,k,l), p), p);
		}
		sol[k] = zpzs_mul(sol[k], zpzs_inv(lit_coeff_zpzs(s,k,k),p), p);
	}
	return;
}



// estes chinois. Prend le résultat donné par la fonction chinois_interne, et le convertit en représentation symétrique
void chinois_sym(mpz_t res,mpz_t x1,mpz_t x2,mpz_t n1,mpz_t n2,mpz_t u){
	// Ancienne version
	/*mpz_t tmp;
	mpz_init(tmp);
	// Conversion des arguments depuis la représentation symétrique
	if(mpz_cmp_si(x1,0) < 0){
		mpz_add(x1,x1,n1);
	}
	if(mpz_cmp_si(x2,0) < 0){
		mpz_add(x2,x2,n2);
	}
	// Calculs (avec les éléments Z/pZ représenté par des entiers entre 0 et p-1)
	chinois(res,x1,x2,n1,n2,u,v);
	// Conversion du résultat vers la représentation symétrique
	mpz_mul(tmp,n1,n2);
	mpz_fdiv_q_ui(tmp,tmp,2);
	if(mpz_cmp(res,tmp) > 0){
		mpz_mul_ui(tmp,tmp,2);
		mpz_sub(res,res,tmp);
	}
	mpz_clear(tmp);*/
	
	// Calculs
	chinois_interne(res,n1,n2,x1,x2,u);
	// Conversions
	mpz_t n1n2,n1n2s2,mn1n2s2;
	mpz_init(n1n2);
	mpz_init(n1n2s2);
	mpz_init(mn1n2s2);
	mpz_mul(n1n2,n1,n2);
	mpz_fdiv_q_ui(n1n2s2,n1n2,2);
	mpz_neg(mn1n2s2,n1n2s2);
	if(mpz_cmp(res,n1n2s2) > 0){
		mpz_sub(res,res,n1n2);
	}else if(mpz_cmp(res,mn1n2s2) < 0){
		mpz_add(res,res,n1n2);
	}
	//assert((mpz_cmp(res,n1n2s2) <= 0) && (mpz_cmp(res,mn1n2s2) >= 0));	// OK
	mpz_clear(n1n2);
	mpz_clear(n1n2s2);
	mpz_clear(mn1n2s2);
	return;
}

// Semblable à chinois_n et chinois_old_n, mais appelle chinois_sym
void chinois_sym_n(int n,mpz_t* res,mpz_t* x1,mpz_t* x2,mpz_t n1,mpz_t n2){
	mpz_t u;
	mpz_init(u);
	mpz_gcdext(NULL,u,NULL,n1,n2);
	for(int i = 0;i < n;i++){
		chinois_sym(res[i],x1[i],x2[i],n1,n2,u);
	}
	mpz_clear(u);
	return;
}

// Renvoie (ie écrit dans dets, qui est un tableau de taille s->n+1) le déterminant modulo s->p du système (en dernier), et les déterminants modulo s->p du système avec le 2nd membre à la place de chaque colonne (pas très clair)
void zpzs_dets(int* dets,syst_zpzs* s){
	int n = s->n;
	int p = s->p;
	// Échelonnage (ou "échelonnement" ?)
	zpzs_gauss(s);
	// Calcul du déterminant
	dets[n] = 1;
	for(int i = 0;i < n;i++){
		dets[n] = zpzs_mul(dets[n],lit_coeff_zpzs(s,i,i),p);
	}
	// Calcul des déterminants "avec second membre" (en utilisant la règle de Cramer)
	zpzs_sol_syst_ech(dets,s);
	for(int i = 0;i < n;i++){
		dets[i] = zpzs_mul(dets[i],dets[n],p);
	}
	return;
}

// Version avec (borne de Hadamard et) calcul des déterminants
void modulaire_dets(systeme* s,rationnel* sol,gmp_randstate_t state,mp_bitcnt_t b){
	// INITIALISATION
	int j = 1;		// DEBUG, nombre d'itérations effectuées
	clock_t debut,fin;		// DEBUG
	double zpzs_tps = 0.;
	double chinois_tps = 0.;
	int n = s -> n;
	syst_zpzs sz;	// On va initialiser/copier et détruire ce système à chaque itération
	int* dets = malloc((n+1)*sizeof(int));		// Déterminants (dans Z/pZ, du système avec le second membre à la place de la 1ème colonne, puis du système sans modification (en position n)), version entiers machine
	mpz_t* dets_mpz = malloc((n+1)*sizeof(mpz_t));		// Déterminants (dans Z/pZ, du système avec le second membre à la place de la 1ème colonne, puis du système sans modification (en position n)), version GMP
	mpz_t* dets_tot = malloc((n+1)*sizeof(mpz_t));	// De même dans Z/prodZ
	mpz_t* dets_tot_old = malloc((n+1)*sizeof(mpz_t));	// dets_tot de l'itération précédente
	mpz_t p_mpz,prod_old,prod,hada,det;
	int p;		// Nombre premier "en cours d'utilisation", version machine
	mpz_init(p_mpz);		// Nombre premier "en cours d'utilisation", version GMP
	mpz_init(prod_old);		// Produit des nombres premiers précédemment utilisés (sans p_mpz)
	mpz_init(prod);			// Produit des nombres premiers précédemment utilisés (avec p_mpz) -> pendant une itération, prod == p*prod_old
	mpz_init(hada);		// Contiendra (le double de) la borne de Hadamard du système
	mpz_init(det);		// Contiendra le déterminant de sz (dans Z/pZ)
	for(int i = 0;i < n+1;i++){
		mpz_init(dets_mpz[i]);
		mpz_init(dets_tot[i]);
		mpz_init(dets_tot_old[i]);
	}
	// Calcul de la "borne de Hadamard" (en fait 2 fois la borne de Hadamard)
	hadamard(hada,s);
	mpz_mul_si(hada,hada,2);
	//mpz_mul(hada,hada,hada);	// DEBUG pour s'assurer d'un truc...
	// PREMIÈRE ITÉRATION (On verra plus tard si on peut la faire rentrer dans la boucle principale)
	// "Choix" d'un nombre premier
	p = genere_p(p_mpz,state,b);
	//fprintf(stderr,"%d\n\n",p);	// DEBUG
	// Initialisation de prod et prod_old
	mpz_set(prod,p_mpz);
	mpz_set_ui(prod_old,1);
	// Calcul des déterminants dans Z/pZ
	init_copie_syst_zpzs(&sz,s,p);
	debut = clock();	// DEBUG
	zpzs_dets(dets,&sz);
	fin = clock();	// DEBUG
		zpzs_tps += ((double)(fin-debut))/CLOCKS_PER_SEC;	// DEBUG
	// Copie des déterminants dans dets_tot directement (pour cette étape seulement)
	for(int i = 0;i < n+1;i++){
		mpz_set_si(dets_tot[i],dets[i]);
		//fprintf(stderr,"%d\n",dets[i]);			// DEBUG
		//mpz_out_str(stderr,10,dets_tot[i]);			// DEBUG
		//fputc('\n',stderr);			// DEBUG
	}
	// MàJ de prod_old
	mpz_set(prod_old,p_mpz);
	// Destruction du syst_zpzs, pour pouvoir le réutiliser
	detruit_syst_zpzs(&sz);
	// BOUCLE PRINCIPALE
	while(mpz_cmp(prod,hada) < 0){
		j++;	// DEBUG
		// "Choix" d'un nombre premier pour cette itération
		p = genere_p(p_mpz,state,b);
		// MàJ de prod
		mpz_mul(prod,prod,p_mpz);
		// Copie des précédents déterminants sur l'emplacement dédié
		for(int i = 0;i < n+1;i++){
			mpz_set(dets_tot_old[i],dets_tot[i]);
		}
		// Calcul des déterminants dans Z/pZ
		init_copie_syst_zpzs(&sz,s,p);
		debut = clock();	//DEBUG
		zpzs_dets(dets,&sz);
		fin = clock();	// DEBUG
		zpzs_tps += ((double)(fin-debut))/CLOCKS_PER_SEC;	// DEBUG
		// Copie des déterminants dans dets_mpz
		for(int i = 0;i < n+1;i++){
			mpz_set_si(dets_mpz[i],dets[i]);
			/*fprintf(stderr,"%d\n",dets[i]);			// DEBUG
			mpz_out_str(stderr,10,dets_mpz[i]);			// DEBUG
			fputc('\n',stderr);			// DEBUG
			mpz_out_str(stderr,10,dets_tot[i]);			// DEBUG
			fprintf(stderr,"\n\n");			// DEBUG*/
		}
		// Restes chinois avec les déterminants précédents
		debut = clock();	// DEBUG
		chinois_sym_n(n+1,dets_tot,dets_tot_old,dets_mpz,prod_old,p_mpz);
		//chinois_bis_n(n+1,dets_tot,dets_tot,dets_mpz,prod_old,p_mpz);	// Ne fonctionne pas
		fin = clock();	// DEBUG
		chinois_tps += ((double)(fin-debut))/CLOCKS_PER_SEC;	// DEBUG
		// MàJ de prod_old
		mpz_set(prod_old,prod);
		// Destruction du syst_zpzs, pour pouvoir le réutiliser
		detruit_syst_zpzs(&sz);
	}
	// SOLUTION
	for(int i = 0;i < n;i++){
		rat_set_pq(sol[i],dets_tot[i],dets_tot[n]);
	}
	//fprintf(stderr,"Nombre d'itérations : %d\n\n",j);	// DEBUG
	fprintf(stderr,"Nombre d'itérations : %d\nTemps de résolution dans Z/pZ : %lf s\nTemps de restes chinois : %lf s\n\n",j,zpzs_tps,chinois_tps);	// DEBUG
	// SUPPRESSION DES OBJETS UTILISÉS
	for(int i = 0;i < n+1;i++){
		mpz_clear(dets_mpz[i]);
		mpz_clear(dets_tot[i]);
		mpz_clear(dets_tot_old[i]);
	}
	mpz_clear(p_mpz);
	mpz_clear(prod_old);
	mpz_clear(prod);
	mpz_clear(hada);
	mpz_clear(det);
	free(dets);
	free(dets_mpz);
	free(dets_tot);
	free(dets_tot_old);
	return;
}
