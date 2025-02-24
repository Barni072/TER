#ifndef ZPZ_H
#define ZPZ_H

#include <gmp.h>
#include <stdbool.h>
#include "systemes.h"

// Argument de zpz_multi (nécessaire à cause de pthread_create)
struct s_args{
	mpz_t* sol;
	systeme* s;
	mpz_t p;
};
typedef struct s_args zpz_args;

//void zpz_add(mpz_t res,mpz_t x,mpz_t y,mpz_t p);
//void zpz_sub(mpz_t res,mpz_t x,mpz_t y,mpz_t p);
//void zpz_mul(mpz_t res,mpz_t x,mpz_t y,mpz_t p);
//void init_copie_systeme_zpz(systeme* sdest,systeme* ssrc,mpz_t p);
bool verif_sol_zpz(systeme* s,mpz_t* sol,mpz_t p);
//void zpz_gauss(systeme* s,mpz_t p);
//void zpz_sol_syst_ech(systeme* s,mpz_t* sol,mpz_t p);
void zpz_sans_copie(mpz_t* sol,systeme* s,mpz_t p);
void* zpz_multi(zpz_args* a);

#endif
