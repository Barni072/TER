#ifndef MODDETS_H
#define MODDETS	_H

#include <gmp.h>

typedef struct s_syst_zpz syst_zpzs;

void init_copie_syst_zpzs(syst_zpzs* sdest,systeme* ssrc,long int p);
void detruit_syst_zpzs(syst_zpzs* s);
long int lit_coeff_zpzs(syst_zpzs* s,int i,int j);
void ecrit_coeff_zpzs(syst_zpzs* s,int i,int j,long int n);
void affiche_syst_zpzs(syst_zpzs* s,FILE* f);
bool verif_sol_zpzs(syst_zpzs* s,long int* sol);
//zpzs_gauss(syst_zpzs* s);
//void zpzs_sol_syst_ech(long int* sol,syst_zpzs* s);
//void zpz_dets(long int* dets,syst_zpz* s);
//void chinois_sym(mpz_t res,mpz_t x1,mpz_t x2,mpz_t n1,mpz_t n2,mpz_t u);
//void chinois_sym_n(int n,mpz_t res,mpz_t x1,mpz_t x2,mpz_t n1,mpz_t n2,mpz_t u,mpz_t v){
void modulaire_dets(systeme* s,rationnel* sol,gmp_randstate_t state,mp_bitcnt_t b);

#endif
