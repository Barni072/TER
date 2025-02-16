#ifndef ZPZ_H
#define ZPZ_H

#include <gmp.h>
#include "systemes.h"

//void zpz_add(mpz_t res,mpz_t x,mpz_t y,mpz_t p);
//void zpz_sub(mpz_t res,mpz_t x,mpz_t y,mpz_t p);
//void zpz_mul(mpz_t res,mpz_t x,mpz_t y,mpz_t p);
//void zpz_inv(mpz_t res,mpz_t x,mpz_t p);
void zpz_gauss(systeme* s,mpz_t p);

#endif
