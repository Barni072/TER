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
			mpz_inp_str(s->t[i*m + j],f,10);
		}
	}
	fclose(f);
	return;
}

// Écrit un système dans un fichier (créé si besoin), de taille n, et dont les coefficients sont uniformément tirés entre -(2^(b-1)) et (2^(b-1))-1
void ecrit_fichier_au_pif(char* nom,int n,gmp_randstate_t state,mp_bitcnt_t b){
	int m = n+1;
	FILE* f = fopen(nom,"w+");
	assert(f != NULL);
	mpz_t c,tmp,rnd;
	mpz_init(c);	// Représentera un coefficient
	mpz_init(tmp);
	mpz_init(rnd);
	fprintf(f,"%d %d ",n,m);
	for(int i = 0;i < n;i++){
		for(int j = 0;j < m;j++){
			// Donne à c une valeur tirée uniformément entre 0 et (2^(b-1))-1
			mpz_urandomb(c,state,b);
			// Transforme c en nombre négatif (c <- -1-c) avec une probabilité 1/2
			mpz_urandomb(rnd,state,1);
			if(mpz_cmp_ui(rnd,0) == 0){
				mpz_add_ui(tmp,c,1);
				mpz_neg(c,tmp);
			}
			// Écrit c dans le fichier
			mpz_out_str(f,10,c);
			fputc(' ',f);
		}
	}
	fclose(f);
	mpz_clear(c);
	mpz_clear(tmp);
	mpz_clear(rnd);
	return;
}
