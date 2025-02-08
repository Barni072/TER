#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <gmp.h>
#include "rationnels.h"
#include "systemes.h"

// Serait-il nécessaire d'utiliser des mpz_t pour le nombre de lignes et de colonnes ?
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
			//fprintf(stdout,"%d ",s -> a[i][j]);		// ANCIEN
			lit_coeff(k,s,i,j);
			fprintf(stdout,"%d ",mpz_get_si(k));
		}
		// Affichage d'un coeff du second membre, puis passage à la ligne suivante
		//fprintf(stdout,"  %d\n",s -> b[i]);		// ANCIEN
		lit_coeff(k,s,i,n);
		fprintf(stdout,"  %d\n",mpz_get_si(k));
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
			//res = res && (s->a[i][j] == 0);		// ANCIEN
			//fprintf(stdout,"%d %d %d\n",i,j,s->a[i][j]);	// DEBUG, et ancien
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
	assert(est_echelonne(s));
	int n = s->n;
	//int m = s->m;		// Pas besoin, vaut n+1 pour l'instant
	mpz_t e;
	mpz_init(e);
	rationnel r;
	rat_init(&r);
	for(int k = n-1;k >= 0;k--){	// Calcul de chaque coordonnée de la solution :
		// Initialisation
		//sol[k].p = s->b[k];	// ANCIEN
		//sol[k].q = 1;		// ANCIEN
		lit_coeff(e,s,k,n);
		rat_set_ent(sol[k],e);
		for(int l = k+1;l < n;l++){	// Ajoute les "contributions des lignes suivantes", une par une (pas clair) :
			//rationnel x = prod_avec_ent(sol[l],s->a[k][l]);		// ANCIEN
			//sol[k] = somme(sol[k],opp(x));		// ANCIEN
			lit_coeff(e,s,k,l);
			rat_mul_ent(r,sol[l],e);
			rat_sub(sol[k],sol[k],r);	// NE MARCHERA CERTAINEMENT PAS
			// Ce qui précède est ÉCLATÉ, il faudrait peut-être prendre le temps d'implémenter des fonctions rat_addmul et rat_submul
		}
		// Division par le pivot :
		//sol[k] = div_par_ent(sol[k],s->a[k][k]);	// ANCIEN
		lit_coeff(e,s,k,k);
		rat_div_ent(r,sol[k],e);		
		rat_set(sol[k],r);
	}
	mpz_clear(e);
	rat_clear(&r);
	return;
}
