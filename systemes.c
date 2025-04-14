#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <gmp.h>
#include "rationnels.h"
#include "systemes.h"

void init_systeme(systeme* s,int nb_lign,int nb_col){	// Si cette fonction prend un système (au lieu d'un pointeur vers un système), il y a un problème au moment de libérer la mémoire allouée pour s.t
	s->n = nb_lign;
	s->m = nb_col;
	s->t = malloc(sizeof(mpz_t) * nb_lign * nb_col);	// Il serait bon de vérifier que l'allocation de mémoire fonctionne...
	for(int k = 0;k < nb_lign*nb_col;k++){
		mpz_init(s->t[k]);
	}
	return;
}

void detruit_systeme(systeme* s){
	for(int k = 0;k < (s->n)*(s->m);k++){
		mpz_clear(s->t[k]);
	}
	free(s->t);
}

// Initialise le système sdest, et en fait une copie du système ssrc
void init_copie_systeme(systeme* sdest,systeme* ssrc){
	int n = ssrc->n;
	int m = ssrc->m;
	sdest -> n = n;
	sdest -> m = m;
	sdest -> t = malloc(sizeof(mpz_t)*n*m);
	for(int k = 0;k < n*m;k++){
		mpz_init_set(sdest->t[k],ssrc->t[k]);
	}
	return;
}

// Lit le coefficient (i,j) de s, et donne sa valeur à res
// Suppose évidemment 0 <= i < n et 0 <= j < m
void lit_coeff(mpz_t res,systeme* s,int i,int j){
	mpz_set(res,s->t[i*(s->m) + j]);
	return;
}

// Remplace le coefficient (i,j) de s par n
// Suppose évidemment 0 <= i < n et 0 <= j < m
void ecrit_coeff(systeme* s,int i,int j,mpz_t n){
	mpz_set(s->t[i*(s->m) + j],n);
	return;
}

// Échange les lignes i1 et i2 de s
void echange_lignes(systeme* s,int i1,int i2){
	int m = s->m;
	for(int j = 0;j < m;j++){
		mpz_swap(s->t[i1*m + j],s->t[i2*m + j]);
	}
	return;
}

// Affiche un système dans le terminal, pour l'instant de façon matricielle et très moche
void affiche_systeme(systeme* s,FILE* f){
	int n = s->n;
	//int m = s->m;		// Pas besoin, on s'intéresse ici seulement aux coeffs du "1er membre" donc le nombre de lignes suffit
	mpz_t k;
	mpz_init(k);
	for(int i = 0;i < n;i++){
		// Affichage d'une ligne de coeffs
		for(int j = 0;j < n;j++){
			lit_coeff(k,s,i,j);
			mpz_out_str(f,10,k);
			fputc(' ',f);
		}
		// Affichage d'un coeff du second membre, puis passage à la ligne suivante
		lit_coeff(k,s,i,n);
		fputc(' ',f);
		fputc(' ',f);
		mpz_out_str(f,10,k);
		fputc('\n',f);
	}
	mpz_clear(k);
	return;
}

// Renvoie vrai ssi s est échelonné
bool est_echelonne(systeme* s){
	bool res = true;
	int n = s->n;
	//int m = s->m;		// Pas besoin, on s'intéresse ici seulement aux coeffs du "1er membre" donc le nombre de lignes suffit
	mpz_t k;
	mpz_init(k);
	for(int j = 0;j < n;j++){
		for(int i = j+1;i < n;i++){
			lit_coeff(k,s,i,j);
			res = res && (mpz_cmp_ui(k,0) == 0);
		}
	}
	mpz_clear(k);
	return res;
}

// "Remonte" le système s (supposé échelonné) et calcule sa solution
// Le tableau sol sert d'emplacement où la solution sera écrite (mal dit)
// DANGER SI PIVOTS NULS
void sol_syst_echelonne(systeme* s, rationnel* sol){
	//assert(est_echelonne(s));
	int n = s->n;
	//int m = s->m;		// Pas besoin, vaut n+1 pour l'instant
	mpz_t e;
	mpz_init(e);
	rationnel r,r2;
	rat_init(&r);
	rat_init(&r2);
	for(int k = n-1;k >= 0;k--){	// Calcul de chaque coordonnée de la solution :
		// Initialisation
		lit_coeff(e,s,k,n);
		rat_set_ent(sol[k],e);
		for(int l = k+1;l < n;l++){	// Ajoute les "contributions des lignes suivantes", une par une (pas clair) :
			lit_coeff(e,s,k,l);
			rat_mul_ent(r,sol[l],e);
			//rat_sub(sol[k],sol[k],r);	// A D'EXCELLENTES CHANCES DE FOIRER, DANGER
			rat_sub(r2,sol[k],r);
			rat_set(sol[k],r2);
		}
		// Division par le pivot :
		lit_coeff(e,s,k,k);
		rat_div_ent(r,sol[k],e);		
		rat_set(sol[k],r);
	}
	mpz_clear(e);
	rat_clear(&r);
	rat_clear(&r2);
	return;
}

bool verif_sol(systeme* s,rationnel* sol){
	int n = s->n;
	mpz_t coeff,zero;
	rationnel som,tmp1,tmp2;
	mpz_init(coeff);
	mpz_init_set_si(zero,0);
	rat_init(&som);
	rat_init(&tmp1);
	rat_init(&tmp2);
	// On suppose m = n+1;
	bool res = true;
	for(int i = 0;i < n;i++){
		rat_set_ent(som,zero);
		for(int j = 0;j < n;j++){
			lit_coeff(coeff,s,i,j);
			rat_mul_ent(tmp1,sol[j],coeff);
			rat_add(tmp2,som,tmp1);
			rat_set(som,tmp2);
		}
		// Comparaison avec le coeff du second membre :
		lit_coeff(coeff,s,i,n);
		res = res && rat_comp_ent(som,coeff);
	}
	mpz_clear(coeff);
	mpz_clear(zero);
	rat_clear(&som);
	rat_clear(&tmp1);
	rat_clear(&tmp2);
	return res;
}
