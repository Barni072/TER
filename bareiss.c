#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <gmp.h>
#include "systemes.h"

// Échelonne le système selon l'algorithme de (Gauss-)Bareiss, foire si un pivot nul est rencontré
void bareiss_old(systeme* s){
	int n = s->n;	// Nombre de lignes
	// Le nombre de colonnes (en comptant le second membre) est supposé égal à n+1
	mpz_t d;	// Représentera le pivot de l'étape précédente
	mpz_t p;	// Représentera le pivot de l'étape courante
	mpz_t e,ij,ik,kj;	// Entiers "de secours", on s'en servira quand on devra stocker temporairement des mpz_t (par exemple des coeffs du système...)
	mpz_init(d);
	mpz_init(p);
	mpz_init(e);
	mpz_init(ij);
	mpz_init(ik);
	mpz_init(kj);
	for(int k = 0;k < n-1;k++){
		if(k == 0) mpz_set_si(d,1);
		else lit_coeff(d,s,k-1,k-1);
		assert(mpz_sgn(d) != 0);		// Pas propre, détecte les pivots nuls
		lit_coeff(p,s,k,k);		
		// Pour chaque pivot, on va modifier tous les coeffs "en bas à droite" du pivot, y compris ceux dans le second membre
		for(int i = k+1;i < n;i++){
			for(int j = k+1;j < n+1;j++){	// j va jusqu'à n (inclus) car on souhaite aussi modifier le second membre
				lit_coeff(ij,s,i,j);
				lit_coeff(ik,s,i,k);
				lit_coeff(kj,s,k,j);
				mpz_mul(e,ij,p);
				mpz_submul(e,ik,kj);
				mpz_divexact(ij,e,d);
				ecrit_coeff(s,i,j,ij);
			}
			mpz_set_si(e,0);
			ecrit_coeff(s,i,k,e);
		}
	}
	mpz_clear(d);
	mpz_clear(p);
	mpz_clear(e);
	mpz_clear(ij);
	mpz_clear(ik);
	mpz_clear(kj);
	return;
}

// Échelonne le système selon l'algorithme de (Gauss-)Bareiss
// Pour chaque pivot, échange les lignes de sorte à choisir le pivot (non nul) le plus petit (en valeur absolue)
 void bareiss(systeme* s){
	int n = s->n;	// Nombre de lignes
	// Le nombre de colonnes (en comptant le second membre) est supposé égal à n+1
	mpz_t d;	// Représentera le pivot de l'étape précédente
	mpz_t p;	// Représentera le pivot de l'étape courante
	mpz_t e,ij,ik,kj;	// Entiers "de secours", on s'en servira quand on devra stocker temporairement des mpz_t (par exemple des coeffs du système...)
	mpz_init(d);
	mpz_init(p);
	mpz_init(e);
	mpz_init(ij);
	mpz_init(ik);
	mpz_init(kj);
	for(int k = 0;k < n-1;k++){
		// Valeur du pivot précédent
		if(k == 0) mpz_set_si(d,1);
		else lit_coeff(d,s,k-1,k-1);
		// Recherche du plus petit (en valeur absolue) pivot possible (non nul) :
		int k_new = k;	// Indice de ligne du pivot qui sera choisi
		lit_coeff(p,s,k,k);
		//mpz_out_str(stderr,10,p);fputc('\n',stderr);	// DEBUG
		mpz_abs(p,p);
		for(int l = k+1;l < n;l++){
			lit_coeff(e,s,l,k);
			mpz_abs(e,e);
			if((mpz_cmp(p,e) > 0) && (mpz_cmp_ui(e,0) != 0)){
				k_new = l;
				mpz_set(p,e);
			}
		}
		echange_lignes(s,k,k_new);	// Échange de lignes pour avoir le plus petit pivot possible (non nul)
		lit_coeff(p,s,k,k);	// On doit lire le pivot à nouveau, car ce qui était précédemment stocké dans p ne prenait pas compte du signe
		if(mpz_cmp_ui(p,0) == 0){	// ie si tous les "candidats pivots" sont nuls
			fprintf(stderr,"BAREISS : Pas de %d-ième pivot non nul, le système n'est pas de Cramer, abandon.\n",k+1);
			exit(1);
		}
		//mpz_out_str(stderr,10,p);fputc('\n',stderr);fputc('\n',stderr);// DEBUG
		// Pour chaque pivot, on va modifier tous les coeffs "en bas à droite" du pivot, y compris ceux dans le second membre
		for(int i = k+1;i < n;i++){
			for(int j = k+1;j < n+1;j++){	// j va jusqu'à n (inclus) car on souhaite aussi modifier le second membre
				lit_coeff(ij,s,i,j);
				lit_coeff(ik,s,i,k);
				lit_coeff(kj,s,k,j);
				mpz_mul(e,ij,p);
				mpz_submul(e,ik,kj);
				mpz_divexact(ij,e,d);
				ecrit_coeff(s,i,j,ij);
			}
			mpz_set_si(e,0);
			ecrit_coeff(s,i,k,e);
		}
	}
	mpz_clear(d);
	mpz_clear(p);
	mpz_clear(e);
	mpz_clear(ij);
	mpz_clear(ik);
	mpz_clear(kj);
	return;
}

// Les deux variantes mettent essentiellement le même temps à s'exécuter
// La deuxième est un peu plus "sûre", au sens où elle ne peut pas échouer sur un système de Cramer
// (+ Un test rapide semble prouver que le pivot est bien changé)
