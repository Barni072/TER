#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <gmp.h>
#include "rationnels.h"
#include "systemes.h"

// L'objectif de ce fichier est d'implémenter (de façon assez rapide et sale) les systèmes à coefficients rationnels et l'algo de Gauss "classique"

// Structure de système à coeffs rationnels, ne servira que pour l'algo de Gauss "classique", donc qu'ici
struct s_sys_rat{
	int n;	// Nombre de lignes (nombre d'équations et de variables)
	int m;	// Nombre de colonnes (pour l'instant = n+1)
	rationnel* t;	// Le tableau contenant tous les coeffients
};
typedef struct s_sys_rat sys_rat;

// Initialise un système à coeffs rationnels, en copiant un système à coeff entiers
void init_sys_rat(sys_rat* sr,systeme* s){
	sr->n = s->n;
	sr->m = s->m;
	sr->t = malloc(sizeof(rationnel) * s->n * s->m);
	for(int k = 0;k < (s->n)*(s->m);k++){
		rat_init(&(sr->t[k]));
		rat_set_ent(sr->t[k],s->t[k]);
	}
	return;
}

void detruit_sys_rat(sys_rat* sr){
	for(int k = 0;k < (sr->n)*(sr->m);k++){
		rat_clear(&(sr->t[k]));
	}
	free(sr->t);
}

void lit_coeff_rat(rationnel res,sys_rat* s,int i,int j){
	rat_set(res,s->t[i*(s->m) + j]);
	return;
}

void ecrit_coeff_rat(sys_rat* s,int i,int j,rationnel n){
	rat_set(s->t[i*(s->m) + j],n);
	return;
}

// Échelonne le système, avec l'algo de Gauss "classique"
// Fait à l'arrache, à bien tester
void gauss_ech(sys_rat* sr){
	int n = sr->n;
	//int m = sr->m;	// Pas besoin, m=n+1
	rationnel p,ij,ik,kj,a,b,c;
	rat_init(&p);
	rat_init(&ij);
	rat_init(&ik);
	rat_init(&kj);
	rat_init(&a);
	rat_init(&b);
	rat_init(&c);
	for(int k = 0;k < n-1;k++){
		lit_coeff_rat(p,sr,k,k);	// Pivot
		for(int i = k+1;i < n;i++){
			lit_coeff_rat(ik,sr,i,k);
			rat_div(a,ik,p);
			for(int j = k+1;j < n+1;j++){
				lit_coeff_rat(ij,sr,i,j);
				lit_coeff_rat(kj,sr,k,j);
				rat_mul(b,kj,a);
				rat_sub(c,ij,b);
				ecrit_coeff_rat(sr,i,j,c);
			}
			rat_set_si(b,0);
			ecrit_coeff_rat(sr,i,k,b);
		}
	}
	rat_clear(&p);
	rat_clear(&ij);
	rat_clear(&ik);
	rat_clear(&kj);
	rat_clear(&a);
	rat_clear(&b);
	rat_clear(&c);
	return;
}

// Donne la solution d'un systeme à coeffs rationnels préalablement échelonné
// Copié assez sauvagement sol_syst_echelonne, à bien tester
void sol_sys_rat_ech(sys_rat* sr,rationnel* sol){
	int n = sr->n;
	//int m = sr->m;	// Pas besoin, m=n+1
	rationnel r,s;
	rat_init(&r);
	rat_init(&s);
	for(int k = n-1;k >= 0;k--){
		lit_coeff_rat(sol[k],sr,k,n);
		for(int l = k+1;l < n;l++){
			lit_coeff_rat(r,sr,k,l);
			rat_mul(s,sol[l],r);
			rat_sub(r,sol[k],s);
			rat_set(sol[k],r);
		}
		lit_coeff_rat(r,sr,k,k);
		rat_div(s,sol[k],r);
		rat_set(sol[k],s);
	}
	rat_clear(&r);
	rat_clear(&s);
	return;
}

void gauss(systeme* s,rationnel* sol){
	sys_rat sr;
	init_sys_rat(&sr,s);
	
	gauss_ech(&sr);
	sol_sys_rat_ech(&sr,sol);
	
	detruit_sys_rat(&sr);
	return;
}
