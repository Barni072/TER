#ifndef ZPZ_H
#define ZPZ_H

#include <gmp.h>
#include <stdbool.h>
#include "systemes.h"

// Structure de système à coefficients dans Z/pZ
struct s_syst_zpz{
	int n;	// Nombre de lignes
	int m;	// Nombre de colonnes (m = n+1 en pratique)
	int p;	// Nombre premier
	int* t;	// Tableau contenant les coefficients
};
typedef struct s_syst_zpz syst_zpz;

//int zpz_add(int x,int y,int p);
//int zpz_sub(int x,int y,int p);
//int zpz_mul(int x,int y,int p);
//void euclide_etendu(int* pgcd,int* u,int* v,int a,int b);
//int pgcd(int a,int b);
//int zpz_inv(int x,int p);
void init_copie_syst_zpz(syst_zpz* sdest,systeme* ssrc,int p);
void detruit_syst_zpz(syst_zpz* s);
//int lit_coeff_zpz(syst_zpz* s,int i,int j);
//void ecrit_coeff_zpz(syst_zpz* s,int i,int j,int n);
void affiche_syst_zpz(syst_zpz* s,FILE* f);
bool verif_sol_zpz(syst_zpz* s,int* sol);
//void zpz_gauss(syst_zpz* s);
//void zpz_sol_syst_ech(syst_zpz* s,int* sol);
void zpz_resol(syst_zpz* s,int* sol);
int genere_p(mpz_t p,gmp_randstate_t state,mp_bitcnt_t b);
bool sol_egales(rationnel* sol1,rationnel* sol2,int n);
//void chinois(mpz_t res,mpz_t x1,mpz_t x2,mpz_t n1,mpz_t n2,mpz_t u,mpz_t v);
void chinois_n(int n,mpz_t* res,mpz_t* x1,mpz_t* x2,mpz_t n1,mpz_t n2);
void euclide_etendu_borne(mpz_t r,mpz_t v,mpz_t a,mpz_t b);
void modulaire_old(systeme* s,rationnel* sol,gmp_randstate_t state,mp_bitcnt_t b);
void hadamard(mpz_t res,systeme* s);
void modulaire(systeme* s,rationnel* sol,gmp_randstate_t state,mp_bitcnt_t b);
//void zpz_dets(int* dets,syst_zpz* s);
void modulaire_dets(systeme* s,rationnel* sol,gmp_randstate_t state,mp_bitcnt_t b);

#endif
