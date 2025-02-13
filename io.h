#ifndef IO_H
#define IO_H

#include <gmp.h>
#include "systemes.h"

void init_lit_systeme(systeme* s,char* nom);
void ecrit_fichier_au_pif(char* nom,int n,gmp_randstate_t state,mp_bitcnt_t b);

#endif
