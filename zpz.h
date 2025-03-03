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

// Argument de zpz_thrd (nécessaire à cause de pthread_create)
struct s_args{
	syst_zpz s;
	int* sol;
};
typedef struct s_args zpz_args;

//int zpz_add(int x,int y,int p);
//int zpz_sub(int x,int y,int p);
//int zpz_mul(int x,int y,int p);
//void euclide_etendu(int* pgcd,int* u,int* v,int a,int b);
//int pgcd(int a,int b);
//int zpz_inv(int x,int p);
//int chinois(int x1,int x2,int p1,int p2);
//int chinois_uv_connus(int x1,int x2,int p1,int p2,int u,int v);
void init_copie_syst_zpz(syst_zpz* sdest,systeme* ssrc,int p);
void detruit_syst_zpz(syst_zpz* s);
//int lit_coeff_zpz(syst_zpz* s,int i,int j);
//void ecrit_coeff_zpz(syst_zpz* s,int i,int j,int n);
void affiche_syst_zpz(syst_zpz* s,FILE* f);
bool verif_sol_zpz(syst_zpz* s,int* sol);
//void zpz_gauss(syst_zpz* s);
//void zpz_sol_syst_ech(syst_zpz* s,int* sol);
void zpz(int* sol,syst_zpz* s);
void* zpz_thrd(void* a_);	// a_ est en fait un pointeur vers un zpz_args

#endif
