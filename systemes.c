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

// Copie le système ssrc sur le système sdest (qui doit déjà avoir été initialisé)
// En pratique, on utilisera probablement surtout init_copie_systeme
void copie_systeme(systeme* sdest,systeme* ssrc){
	int n = ssrc->n;
	int m = ssrc->m;
	sdest->n = n;
	sdest->m = m;
	for(int k = 0;k < n*m;k++){
		mpz_set(sdest->t[k],ssrc->t[k]);
	}
	return;
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

// Affiche un système dans le terminal, pour l'instant de façon matricielle et très moche
void affiche_systeme(systeme* s){
	int n = s->n;
	//int m = s->m;		// Pas besoin, on s'intéresse ici seulement aux coeffs du "1er membre" donc le nombre de lignes suffit
	mpz_t k;
	mpz_init(k);
	//fprintf(stdout,"Taille : %d\n",n);	 // Cette information n'est pas très utile
	for(int i = 0;i < n;i++){
		// Affichage d'une ligne de coeffs
		for(int j = 0;j < n;j++){
			lit_coeff(k,s,i,j);
			//fprintf(stdout,"%d ",mpz_get_si(k));
			mpz_out_str(stdout,10,k);
			fputc(' ',stdout);
		}
		// Affichage d'un coeff du second membre, puis passage à la ligne suivante
		lit_coeff(k,s,i,n);
		//fprintf(stdout,"  %d\n",mpz_get_si(k));
		fputc(' ',stdout);
		fputc(' ',stdout);
		mpz_out_str(stdout,10,k);
		fputc('\n',stdout);
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
			//rat_addmul(s,coeff,sol[j]);	// On n'a pas implémenté de rat_addmul, c'est ballot
			rat_mul_ent(tmp1,sol[j],coeff);
			rat_add(tmp2,som,tmp1);	// Car "rat_add(s,s,tmp1);" foirerait
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
