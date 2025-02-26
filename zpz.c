#include <assert.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <gmp.h>
#include "systemes.h"
#include "zpz.h"

// On représentera les éléments de Z/pZ par des mpz_t
// TODO : passer à des int (car p s'écrit en au plus 30 bits)

// Opérations basiques sur les éléments de Z/pZ (p premier)
void zpz_add(mpz_t res,mpz_t x,mpz_t y,mpz_t p){
	mpz_add(res,x,y);
	mpz_mod(res,res,p);
	return;
}
void zpz_sub(mpz_t res,mpz_t x,mpz_t y,mpz_t p){
	mpz_sub(res,x,y);
	mpz_mod(res,res,p);
	return;
}
void zpz_mul(mpz_t res,mpz_t x,mpz_t y,mpz_t p){
	mpz_mul(res,x,y);
	mpz_mod(res,res,p);
	return;
}
/*void zpz_inv(mpz_t res,mpz_t x,mpz_t p){
	//assert(mpz_cmp_si(x,0) != 0);
	mpz_t u,v;
	mpz_init(u);
	mpz_init(v);
	mpz_gcdext(g,u,v,x,p);
	mpz_abs(res,u);
	mpz_clear(u);
	mpz_clear(v);
	return;
}*/		// Plutôt utiliser : mpz_invert(res,x,p);

// Initialise le système sdest, et en fait une copie dans Z/pZ du système ssrc
// Les coefficients de sdest sont les restes de la division euclidienne des coefficients de ssrc par p
void init_copie_systeme_zpz(systeme* sdest,systeme* ssrc,mpz_t p){
	int n = ssrc->n;
	int m = ssrc->m;
	sdest -> n = n;
	sdest -> m = m;
	sdest -> t = malloc(sizeof(mpz_t)*n*m);
	for(int k = 0;k < n*m;k++){
		mpz_init(sdest->t[k]);
		mpz_mod(sdest->t[k],ssrc->t[k],p);
	}
	return;
}

// Modifie le système s, en remplaçant tous ses coefficients par les restes de leurs division euclidienne par p
// Il est impératif d'appeler cette fonction sur s avant d'utiliser celles qui suivent
void systeme_mod_p(systeme* s,mpz_t p){
	int n = s->n;
	int m = s->m;
	for(int k = 0;k < n*m;k++){
		mpz_mod(s->t[k],s->t[k],p);
	}
	return;
}

// Vérifie si sol est une solution de s dans Z/pZ
// Les coefficients de s doivent être compris entre 0 et p-1
bool verif_sol_zpz(systeme* s,mpz_t* sol,mpz_t p){
	int n = s->n;
	bool res = true;
	mpz_t ij,som,tmp;
	mpz_init(ij);
	mpz_init(som);
	mpz_init(tmp);
	affiche_systeme(s);		// DEBUG
	for(int i = 0; i < n;i++){
		mpz_set_si(som,0);
		for(int j = 0;j < n;j++){
			lit_coeff(ij,s,i,j);
			zpz_mul(tmp,ij,sol[j],p);
			zpz_add(som,som,tmp,p);
		}
		lit_coeff(tmp,s,i,n);
			{mpz_out_str(stderr,10,som);
			fputc(' ',stderr);
			mpz_out_str(stderr,10,tmp);
			fputc('\n',stderr);}		// DEBUG
		res = res && (mpz_cmp(som,tmp) == 0);
	}
	mpz_clear(ij);
	mpz_clear(som);
	mpz_clear(tmp);
	return res;
}

// Échelonne le système s dans Z/pZ, avec l'algorithme de Gauss
// Les coefficients de s doivent être compris entre 0 et p-1
void zpz_gauss(systeme* s,mpz_t p){
	int n = s->n;
	mpz_t piv,ivp,ij,ik,kj,a,b,c;
	mpz_init(piv);
	mpz_init(ivp);
	mpz_init(ij);
	mpz_init(ik);
	mpz_init(kj);
	mpz_init(a);
	mpz_init(b);
	mpz_init(c);
	// Transforme le système de sorte que tous les coeffs soient entre 0 et p-1
	for(int i = 0;i < n;i++){
		for(int j = 0;j < n+1;j++){
			lit_coeff(ij,s,i,j);
			mpz_mod(a,ij,p);
			ecrit_coeff(s,i,j,a);
		}
	}
	// Échelonne le système
	for(int k = 0;k < n;k++){
		lit_coeff(piv,s,k,k);	// piv : pivot (supposé non nul)
		//zpz_inv(ivp,piv,p);
		mpz_invert(ivp,piv,p);	// ivp : inverse du pivot dans Z/pZ
		for(int i = k+1;i < n;i++){
			lit_coeff(ik,s,i,k);
			zpz_mul(a,ivp,ik,p);	// (Ligne i) <- (Ligne i) - a * (Ligne k)
			for(int j = k+1;j < n+1;j++){
				lit_coeff(ij,s,i,j);
				lit_coeff(kj,s,k,j);
				zpz_mul(b,a,kj,p);
				zpz_sub(c,ij,b,p);
				ecrit_coeff(s,i,j,c);
			}
			mpz_set_si(b,0);
			ecrit_coeff(s,i,k,b);
		}
	}
	mpz_clear(piv);
	mpz_clear(ivp);
	mpz_clear(ij);
	mpz_clear(ik);
	mpz_clear(kj);
	mpz_clear(a);
	mpz_clear(b);
	mpz_clear(c);
	return;
}

// Copie dans sol la solution du système préalablement échelonné s
// Les coefficients de s doivent être compris entre 0 et p-1
void zpz_sol_syst_ech(mpz_t* sol,systeme* s,mpz_t p){
	int n = s -> n;
	mpz_t a,b;
	mpz_init(a);
	mpz_init(b);
	for(int k = n-1;k >= 0;k--){
		lit_coeff(sol[k],s,k,n);
		for(int l = k+1;l < n;l++){
			lit_coeff(a,s,k,l);
			zpz_mul(b,sol[l],a,p);
			zpz_sub(sol[k],sol[k],b,p);	
		}
		lit_coeff(a,s,k,k);
		mpz_invert(b,a,p);
		zpz_mul(sol[k],sol[k],b,p);
	}
	mpz_clear(a);
	mpz_clear(b);
	return;
}

// Écrit dans sol la solution du systeme s dans Z/pZ, en travaillant directement sur s
void zpz_sans_copie(mpz_t* sol,systeme* s,mpz_t p){
	zpz_gauss(s,p);
	zpz_sol_syst_ech(sol,s,p);
	return;
}

// Écrit dans sol la solution du système s dans Z/pZ, sans modifier s
// (sol, s et p sont "contenus" dans *a)
void* zpz_multi(void* a_){
	zpz_args* a = (zpz_args*)a_;
	systeme s_zpz;
	init_copie_systeme_zpz(&s_zpz,a->s,a->p);	// Faire cette étape avant, sinon ça a l'air d'échouer sur le 2ème thread
		{fprintf(stderr,"p = ");
		mpz_out_str(stderr,10,a->p);
		fputc('\n',stderr);}	// DEBUG
	zpz_gauss(&s_zpz,a->p);
	zpz_sol_syst_ech(a->sol,&s_zpz,a->p);
	assert(verif_sol_zpz(&s_zpz,a->sol,a->p));	// Test éclaté car les coeffs de s ne sont pas nécessairement entre 0 et p-1
	assert(verif_sol_zpz(a->s,a->sol,a->p));		// Test mieux, on l'enlèvera plus tard...
	detruit_systeme(&s_zpz);
	return NULL;
}
