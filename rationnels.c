#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <gmp.h>
#include <stdbool.h>
#include "rationnels.h"

// Initialisation et destruction
void rat_init(rationnel* x){
	*x = malloc(sizeof(struct s_rationnel));
	// Par défaut, x est initialisé à 0 (et surtout, son dénominateur est non nul)
	mpz_init_set_ui((*x)->p,0);
	mpz_init_set_ui((*x)->q,1);
	return;
}
void rat_clear(rationnel* x){
	mpz_clear((*x)->p);
	mpz_clear((*x)->q);
	free(*x);
	return;
}

// "Normalise" le rationnel x de sorte que son numérateur et son dénominateur soient premiers entre eux
void rat_norm(rationnel x){
	mpz_t pgcd;
	mpz_init(pgcd);
	mpz_gcd(pgcd,x->p,x->q);
	mpz_divexact(x->p,x->p,pgcd);
	mpz_divexact(x->q,x->q,pgcd);
	mpz_clear(pgcd);
	return;
}

// Donne à x la valeur de y
void rat_set(rationnel x,rationnel y){
	mpz_set(x->p,y->p);
	mpz_set(x->q,y->q);
	return;
}
// Change le numérateur de x
void rat_set_p(rationnel x,mpz_t p){
	mpz_set(x->p,p);
	rat_norm(x);
	return;
}
// Change le dénominateur de x
void rat_set_q(rationnel x,mpz_t q){
	mpz_set(x->q,q);
	rat_norm(x);
	return;
}
// Change le numérateur et le dénominateur de x
void rat_set_pq(rationnel x,mpz_t p,mpz_t q){
	mpz_set(x->p,p);
	mpz_set(x->q,q);
	rat_norm(x);
	return;
}
// Donne une valeur entière à x
void rat_set_ent(rationnel x,mpz_t n){
	mpz_set(x->p,n);
	mpz_set_ui(x->q,1);
	return;
}
void rat_set_si(rationnel x,long int n){
	mpz_set_si(x->p,n);
	mpz_set_ui(x->q,1);
	return;
}

// Échange deux rationnels "efficacement"
void rat_swap(rationnel x,rationnel y){
	mpz_swap(x->p,y->p);
	mpz_swap(x->q,y->q);
	return;
}

// Affiche le numérateur, le dénominateur et une approximation flottante du rationnel x
void rat_aff(rationnel x,FILE* f){
	//long int p = mpz_get_si(x->p);
	//unsigned long int q = mpz_get_ui(x->q);
	//double esti = ((double)p/(double)q);
	fprintf(f,"Numérateur : ");
	mpz_out_str(f,10,x->p);
	fprintf(f,"\nDénominateur : ");
	mpz_out_str(f,10,x->q);
	//fprintf(f,"\nApproximation : %f\n"),esti;		// Foireux, semble toujours afficher 0 -> recherches à faire
	fputc('\n',f);
	return;
}

// Opérations de base sur les rationnels (somme, produit, opposé, inverse, produit avec un entier, division par un entier) :
// Remplace res par le résultat de l'opération (sur x et éventuellement y) voulue
// Un seul rationnel ne peut pas à la fois servir de résultat et d'opérande
void rat_add(rationnel res,rationnel x,rationnel y){
	mpz_mul(res->p,x->p,y->q);
	mpz_addmul(res->p,y->p,x->q);
	mpz_mul(res->q,x->q,y->q);
	rat_norm(res);
	return;
}
void rat_sub(rationnel res,rationnel x,rationnel y){
	mpz_mul(res->p,x->p,y->q);
	mpz_submul(res->p,y->p,x->q);
	mpz_mul(res->q,x->q,y->q);
	rat_norm(res);
	return;
}
void rat_mul(rationnel res,rationnel x,rationnel y){
	mpz_mul(res->p,x->p,y->p);
	mpz_mul(res->q,x->q,y->q);
	rat_norm(res);
	return;
}
void rat_neg(rationnel res,rationnel x){
	mpz_neg(res->p,x->p);
	mpz_set(res->q,x->q);
	// Pas de normalisation à faire ici
	return;
}
void rat_inv(rationnel res,rationnel x){
	int cmp = mpz_cmp_si(x->p,0);	// Le signe de x est celui de x->p
	//assert(cmp != 0);
	if(cmp > 0){
		mpz_set(res->p,x->q);
		mpz_set(res->q,x->p);
	}else{	// cmp < 0	// On veut que x->q reste positif
		mpz_neg(res->p,x->q);
		mpz_neg(res->q,x->p);
	}
	// Pas de normalisation à faire ici
	return;
}
void rat_div(rationnel res,rationnel x,rationnel y){
	rationnel tmp;
	rat_init(&tmp);
	rat_inv(tmp,y);
	rat_mul(res,x,tmp);
	rat_clear(&tmp);
	return;
}
void rat_add_ent(rationnel res,rationnel x,mpz_t n){
	mpz_set(res->p,x->p);
	mpz_addmul(res->p,x->q,n);
	mpz_set(res->q,x->q);
	rat_norm(res);
	return;
}
void rat_mul_ent(rationnel res,rationnel x,mpz_t n){
	mpz_mul(res->p,x->p,n);
	mpz_set(res->q,x->q);
	rat_norm(res);
	return;
}
void rat_div_ent(rationnel res,rationnel x,mpz_t n){
	int cmp = mpz_sgn(n);
	assert(cmp != 0);
	if(cmp > 0){
		mpz_set(res->p,x->p);
		mpz_mul(res->q,x->q,n);
	}else{	// n < 0	// On veut que x->q reste positif
		mpz_neg(res->p,x->p);
		mpz_mul(res->q,x->q,n);
		mpz_neg(res->q,res->q);
	}
	rat_norm(res);
	return;
}

bool rat_comp(rationnel x,rationnel y){
	return ((mpz_cmp(x->p,y->p) == 0) && (mpz_cmp(x->q,y->q) == 0));
}
bool rat_comp_ent(rationnel x,mpz_t n){
	return ((mpz_cmp(x->p,n) == 0) && (mpz_cmp_si(x->q,1) == 0));
}
bool rat_comp_int(rationnel x,int n){
	return ((mpz_cmp_si(x->p,n) == 0) && (mpz_cmp_si(x->q,1) == 0));
}
