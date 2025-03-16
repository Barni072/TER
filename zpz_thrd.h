#ifndef ZPZTH_H
#define ZPZTH_H

#include <gmp.h>

// Argument de zpz_thrd (nécessaire à cause de pthread_create)
struct s_args1{
	syst_zpz* s;
	int* sol;	// Tableau de taille s->n
	mpz_t* sol_m;	// Tableau de taille s->n
};
typedef struct s_args1 zpz_thrd_args;

/*struct s_args2{
	int n;
	mpz_t* res;
	mpz_t* x1;
	mpz_t* x2;
	mpz_t n1;
	mpz_t n2;
};
typedef struct s_args2 chinois_n_thrd_args;*/

//void* zpz_thrd(void* a_);	// a_ est en fait un pointeur vers un zpz_thrd_args
//void* chinois_n_thrd(void* a_);	// a_ est en fait un pointeur vers un chinois_n_thrd_args
void modulaire_thrd(systeme* s,rationnel* sol,gmp_randstate_t state,mp_bitcnt_t b,int thr);

#endif
