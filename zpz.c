#include <gmp.h>
#include <assert.h>
#include "systemes.h"

// On représentera les éléments de Z/pZ par des mpz_t

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
void zpz_inv(mpz_t res,mpz_t x,mpz_t p){
	//assert(mpz_cmp_si(x,0) != 0);
	/*mpz_t u,v;
	mpz_init(u);
	mpz_init(v);
	mpz_gcdext(g,u,v,x,p);
	mpz_abs(res,u);
	mpz_clear(u);
	mpz_clear(v);*/
	
	// Plus simple : (trop simple ?)
	mpz_invert(res,x,p);
	return;
}

// Échelonne le système s dans Z/pZ, avec l'algorithme de Gauss
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
		zpz_inv(ivp,piv,p);		// ivp : inverse du pivot dans Z/pZ
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
