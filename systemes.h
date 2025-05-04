#ifndef SYSTEMES_H
#define SYSTEMES_H

#include <stdbool.h>
#include <gmp.h>
#include "rationnels.h"

// On représente une système linéaire par un tableau de taille n*(n+1), dont la dernière colonne représente le second membre
struct s_systeme{
	int n;	// Nombre de lignes (nombre d'équations et de variables)
	int m;	// Nombre de colonnes (pour l'instant = n+1)
	mpz_t* t;	// Le tableau contenant tous les coeffients
};
typedef struct s_systeme systeme;

void init_systeme(systeme* s,int nb_lign,int nb_col);
void detruit_systeme(systeme* s);
void init_copie_systeme(systeme* sdest,systeme* ssrc);
void lit_coeff(mpz_t res,systeme* s,int i,int j);
void ecrit_coeff(systeme* s,int i,int j,mpz_t n);
void echange_lignes(systeme* s,int i1,int i2);
void affiche_systeme(systeme* s,FILE* f);
//bool est_echelonne(systeme* s);
void sol_syst_echelonne(systeme* s, rationnel* sol);
bool verif_sol(systeme* s,rationnel* sol);
void echange_lignes(systeme* s,int i1,int i2);
bool sol_egales(rationnel* sol1,rationnel* sol2,int n);

#endif
