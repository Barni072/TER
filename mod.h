#ifndef MOD_H
#define MOD_H

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

void init_copie_syst_zpz(syst_zpz* sdest,systeme* ssrc,int p);
void detruit_syst_zpz(syst_zpz* s);
int lit_coeff_zpz(syst_zpz* s,int i,int j);
void ecrit_coeff_zpz(syst_zpz* s,int i,int j,int n);
void affiche_syst_zpz(syst_zpz* s,FILE* f);
bool verif_sol_zpz(syst_zpz* s,int* sol);

//void zpz_gauss(syst_zpz* s);
//void zpz_sol_syst_ech(syst_zpz* s,int* sol);
void zpz_resol(syst_zpz* s,int* sol);

int genere_p(mpz_t p,gmp_randstate_t state,mp_bitcnt_t b);

void chinois_interne(mpz_t res,mpz_t x1,mpz_t x2,mpz_t n1,mpz_t n2,mpz_t u);
//void chinois(mpz_t res,mpz_t x1,mpz_t x2,mpz_t n1,mpz_t n2,mpz_t u);
void chinois_n(int n,mpz_t* res,mpz_t* x1,mpz_t* x2,mpz_t n1,mpz_t n2);
//void chinois_old(mpz_t res,mpz_t x1,mpz_t x2,mpz_t n1,mpz_t n2,mpz_t u,mpz_t v);
//void chinois_old_n(int n,mpz_t* res,mpz_t* x1,mpz_t* x2,mpz_t n1,mpz_t n2);

void euclide_etendu_borne(mpz_t r,mpz_t v,mpz_t a,mpz_t b);
void hadamard(mpz_t res,systeme* s);

void modulaire_old(systeme* s,rationnel* sol,gmp_randstate_t state,mp_bitcnt_t b);
void modulaire(systeme* s,rationnel* sol,gmp_randstate_t state,mp_bitcnt_t b);

#endif
