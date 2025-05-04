#ifndef MODDETS_H
#define MODDETS	_H

#include <gmp.h>

typedef struct s_syst_zpz syst_zpzs;

void init_copie_syst_zpzs(syst_zpzs* sdest,systeme* ssrc,int p);
void detruit_syst_zpzs(syst_zpzs* s);
int lit_coeff_zpzs(syst_zpzs* s,int i,int j);
void ecrit_coeff_zpzs(syst_zpzs* s,int i,int j,int n);
void affiche_syst_zpzs(syst_zpzs* s,FILE* f);
bool verif_sol_zpzs(syst_zpzs* s,int* sol);
//void zpz_dets(int* dets,syst_zpz* s);
void modulaire_dets(systeme* s,rationnel* sol,gmp_randstate_t state,mp_bitcnt_t b);

#endif
