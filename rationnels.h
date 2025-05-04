#ifndef RATIONNELS_H
#define RATIONNELS_H

#include <stdio.h>
#include <stdbool.h>
#include <gmp.h>

// On représente les nombres rationnels par des quotients d'entiers (premiers entre eux, et avec le dénominateur non nul)
// On a besoin de ça, puisque la représentation flottante n'est pas exacte
// On pourrait juste utilser les rationnels mpq_t de GNU MP, mais programmer ceci constitue un bon moyen de prendre en main les entiers de GNU MP
struct s_rationnel{
	mpz_t p;	// Numérateur
	mpz_t q;	// Dénominateur (strictement positif)
	// p et q doivent être premiers entre eux
};
typedef struct s_rationnel* rationnel;

void rat_init(rationnel* x);
void rat_clear(rationnel* x);
//void rat_norm(rationnel x);	// Pas vraiment nécessaire en dehors de rationnels.c
void rat_set(rationnel x,rationnel y);
void rat_set_p(rationnel x,mpz_t p);
void rat_set_q(rationnel x,mpz_t q);
void rat_set_pq(rationnel x,mpz_t p,mpz_t q);
void rat_set_ent(rationnel x,mpz_t n);
void rat_set_si(rationnel x,long int n);
void rat_swap(rationnel x,rationnel y);
void rat_aff(rationnel x,FILE* f);
void rat_add(rationnel res,rationnel x,rationnel y);
void rat_sub(rationnel res,rationnel x,rationnel y);
void rat_mul(rationnel res,rationnel x,rationnel y);
void rat_neg(rationnel res,rationnel x);
void rat_inv(rationnel res,rationnel x);
void rat_div(rationnel res,rationnel x,rationnel y);
void rat_add_ent(rationnel res,rationnel x,mpz_t n);
void rat_mul_ent(rationnel res,rationnel x,mpz_t n);
void rat_div_ent(rationnel res,rationnel x,mpz_t n);
bool rat_comp(rationnel x,rationnel y);
bool rat_comp_ent(rationnel x,mpz_t n);
bool rat_comp_int(rationnel x,int n);

#endif
