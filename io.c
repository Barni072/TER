#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <assert.h>
#include <stdbool.h>
#include "systemes.h"

// On suppose que le fichier ne contient que des chiffres, des signes moins, des espaces, et éventuellement des sauts de ligne

// Initialise le système s à partir des valeurs contenues dans un fichier
void init_lit_systeme(systeme* s,char* nom){
	FILE* f = fopen(nom,"r");
	assert(f != NULL);
	//int n = lit_int(f);
	//int m = lit_int(f);
	int n,m;
	long int c;
	fscanf(f,"%d",&n);
	fscanf(f,"%d",&m);
	s->n = n;
	s->m = m;
	s->t = malloc(n*m*sizeof(mpz_t));
	for(int i = 0;i < n;i++){
		for(int j = 0;j < m;j++){
			mpz_init(s->t[i*m + j]);
			//lit_mpz(s->t[i*m + j],f);
			//fscanf(f,"%ld",&c);
			//mpz_set_si(s->t[i*m + j],c);
			mpz_inp_str(s->t[i*m + j],f,10);
		}
	}
	fclose(f);
	return;
}

// Écrit un système dans un fichier (créé si besoin)
// La taille du système (son nombre de lignes) est donnée en argument, et les coefficients sont uniformément tirés entre 0 et (2^b)-1 (pour l'instant)
void ecrit_fichier_au_pif(char* nom,int n,gmp_randstate_t state,mp_bitcnt_t b){
	int m = n+1;
	FILE* f = fopen(nom,"w+");
	assert(f != NULL);
	mpz_t c;
	mpz_init(c);
	fprintf(f,"%d %d ",n,m);
	for(int i = 0;i < n;i++){
		for(int j = 0;j < m;j++){
			mpz_urandomb(c,state,b);
			mpz_out_str(f,10,c);
			fputc(' ',f);
		}
	}
	fclose(f);
	mpz_clear(c);
	return;
}
