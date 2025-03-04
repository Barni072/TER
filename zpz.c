#include <assert.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <gmp.h>
#include "systemes.h"
#include "zpz.h"

// On représentera les éléments de Z/pZ par des int (car p s'écrit avec au plus 30 bits)

// Opérations basiques sur les éléments de Z/pZ (p premier)
// Les opérandes sont supposées comprises entre 0 et p-1
int zpz_add(int x,int y,int p){
	int res = x + y;
	if(res >= p) res -= p;
	//Alternative : res = res - p*(res >= p);		// On fera ça plus tard
	assert((0 <= res) && (res < p));
	return res;
}
int zpz_sub(int x,int y,int p){
	int res = x - y;
	if(res <= 0) res += p;
	//Alternative : res = res + p*(res <= 0);		// Plus tard
	assert((0 <= res) && (res < p));
	return res;
}
int zpz_mul(int x,int y,int p){
	int res = x * y;
	res = res % p;	// Ce truc n'est pas le modulo...
	if(res <= 0) res += p;	// Au secours
	assert((0 <= res) && (res < p));
	return res;
}
int zpz_inv(int x,int p){ 	// ÉCLATÉ AU SOL (et à implémenter sans GMP...)
	mpz_t xm,pm,resm;
	mpz_init_set_si(xm,x);
	mpz_init_set_si(pm,p);
	mpz_init(resm);
	mpz_invert(resm,xm,pm);		// Problème ici
	//assert(mpz_invert(resm,xm,pm) != 0);
	mpz_out_str(stderr,10,xm);
	fputc(' ',stderr);
	mpz_out_str(stderr,10,pm);
	fputc(' ',stderr);
	mpz_out_str(stderr,10,resm);
	fputc('\n',stderr);
	int res = mpz_get_si(resm);
	mpz_clear(xm);
	mpz_clear(pm);
	mpz_clear(resm);
	fprintf(stderr,"%d %d %d\n",x,p,res);
	assert(zpz_mul(x,res,p) == 1);	// Pas bon du tout
	return res;
}

void init_copie_syst_zpz(syst_zpz* sdest,systeme* ssrc,int p){
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
		// Alternative (devrait suffire, si j'ai bien compris la documentation) (il doit même exister un moyen de faire ça sans déclarer/initialiser c)
		//sdest->t[k] = mpz_fdiv_r(c,ssrc->t[k],p);
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
int lit_coeff_zpz(syst_zpz* s,int i,int j){
	return s->t[i*(s->m) + j];
}

// Remplace le coefficient (i,j) de s par n
// n est supposé compris entre 0 et (s->p - 1)
void ecrit_coeff_zpz(syst_zpz* s,int i,int j,int n){
	s->t[i*(s->m) + j] = n;
	return;
}

void affiche_syst_zpz(syst_zpz* s,FILE* f){
	//fprintf(f,"p = %d\n",s->p);	// Utile ?
	int n = s->n;
	for(int i = 0;i < n;i++){
		// Affichage d'une ligne de coeffs
		for(int j = 0;j < n;j++){
			fprintf(f,"%d ",lit_coeff_zpz(s,i,j));
		}
		// Affichage d'un coeff du second membre
		fprintf(f,"  %d\n",lit_coeff_zpz(s,i,n));
	}
	return;
}

// Vérifie si sol est une solution de s dans Z/pZ
// Les coefficients de s doivent être compris entre 0 et p-1
bool verif_sol_zpz(syst_zpz* s,int* sol,int p){
	int n = s->n;
	bool res = true;
	int ij,som;
	for(int i = 0; i < n;i++){
		//mpz_set_si(som,0);
		som = 0;
		for(int j = 0;j < n;j++){
			//lit_coeff(ij,s,i,j);
			//zpz_mul(tmp,ij,sol[j],p);
			//zpz_add(som,som,tmp,p);
			som = zpz_add(som, zpz_mul(sol[j], lit_coeff_zpz(s,i,j), p), p);
		}
		//lit_coeff(tmp,s,i,n);
		//res = res && (mpz_cmp(som,tmp) == 0);
		res = res && (som == lit_coeff_zpz(s,i,n));
	}
	return res;
}

// Échelonne le système s dans Z/pZ, avec l'algorithme de Gauss
// Les coefficients de s doivent être compris entre 0 et p-1
void zpz_gauss(syst_zpz* s,int p){
	int n = s->n;
	int ivp,a,b,c;
	// Échelonne le système
	for(int k = 0;k < n;k++){
		ivp = zpz_inv(lit_coeff_zpz(s,k,k),p);	// Inverse du pivot (le pivot est supposé non  nul)
		for(int i = k+1;i < n;i++){
			a = zpz_mul(ivp,lit_coeff_zpz(s,i,k),p);
			// (Ligne i) <- (Ligne i) - a * (Ligne k)
			for(int j = k+1;j < n+1;j++){
				b = zpz_mul(a, lit_coeff_zpz(s,k,j), p);
				c = zpz_sub(lit_coeff_zpz(s,i,j),b,p);
				ecrit_coeff_zpz(s,i,j,c);
				// Version sans intermédiaire, à peu près illisible :
				//ecrit_coeff_zpz(s,i,j, zpz_sub(lit_coeff_zpz(s,i,j), zpz_mul(a, lit_coeff_zpz(s,k,j) ,p) ,p));
			}
			ecrit_coeff_zpz(s,i,k,0);
		}
	}
	return;
}

// Copie dans sol la solution du système préalablement échelonné s
// Les coefficients de s doivent être compris entre 0 et p-1
void zpz_sol_syst_ech(int* sol,syst_zpz* s,int p){
	int n = s -> n;
	for(int k = n-1;k >= 0;k--){
		//lit_coeff(sol[k],s,k,n);
		sol[k] = lit_coeff_zpz(s,k,n);
		for(int l = k+1;l < n;l++){
			//lit_coeff(a,s,k,l);
			//zpz_mul(b,sol[l],a,p);
			//zpz_sub(sol[k],sol[k],b,p);
			sol[k] = zpz_sub(sol[k], zpz_mul(sol[l], lit_coeff_zpz(s,k,l), p), p);
		}
		//lit_coeff(a,s,k,k);
		//mpz_invert(b,a,p);
		//zpz_mul(sol[k],sol[k],b,p);
		sol[k] = zpz_mul(sol[k], zpz_inv(lit_coeff_zpz(s,k,k),p), p);
	}
	return;
}

// Écrit dans sol la solution du systeme s dans Z/pZ, en travaillant directement sur s
void zpz(int* sol,syst_zpz* s,int p){
	zpz_gauss(s,p);
	zpz_sol_syst_ech(sol,s,p);
	return;
}

// Écrit dans sol la solution du système s dans Z/pZ, sans modifier s
// (sol, s et p sont "contenus" dans *a)
/*void* zpz_multi(void* a_){
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
}*/
