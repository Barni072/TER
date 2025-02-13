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
void copie_systeme(systeme* sc,systeme* sd);
void lit_coeff(mpz_t res,systeme* s,int i,int j);
void ecrit_coeff(systeme* s,int i,int j,mpz_t n);
void affiche_systeme(systeme* s);
//bool est_echelonne(systeme* s);		// Ne sert que dans sol_syst_echelonne
void sol_syst_echelonne(systeme* s, rationnel* sol);
bool verif_sol(systeme* s,rationnel* sol);

#endif
