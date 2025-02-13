#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <assert.h>
#include <stdbool.h>
#include "systemes.h"

// On suppose que le fichier ne contient que des chiffres, des signes moins, des espaces, et éventuellement des sauts de ligne

/*
// Lit un entier dans le fichier, et le renvoie
int lit_int(FILE* f){	// Le résultat est en fait un unsigned long int...
	char c = fgetc(f);
	int res = 0;
	assert((c != ' ')&&(c != '\n')&&(c != EOF)&&(c != '-'));
	while((c != ' ')&&(c != '\n')&&(c != EOF)){
		int m = (int)c;
		res = 10*res + m;
		c = fgetc(f);
	}
	return res;
}

//Lit un entier dans le fichier, et l'écrit dans n
void lit_mpz(mpz_t n,FILE* f){
	char c = fgetc(f);
	bool neg = false;
	assert((c != ' ')&&(c != '\n')&&(c != EOF));
	mpz_set_ui(n,0);
	if(c == '-'){
		neg = true;
		c = fgetc(f);
	}
	while((c != ' ')&&(c != '\n')&&(c != EOF)){
		int m = (int)c;
		mpz_mul_ui(n,n,10);
		mpz_add_ui(n,n,m);
		c = fgetc(f);
	}
	if(neg) mpz_neg(n,n);
	return;
}*/
// Ce qui précède est bancal et parfaitement inutile, on peut utiliser fscanf (au moins pour des entiers de taille raisonnable)

// Initialise le système s à partir des valeurs contenues dans un fichier
void init_lit_systeme(systeme* s,char* nom){
	FILE* f = fopen(nom,"r");
	assert(f != NULL);	// ?
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
			fscanf(f,"%ld",&c);
			mpz_set_si(s->t[i*m + j],c);
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
