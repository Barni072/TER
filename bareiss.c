#include <stdio.h>
#include <assert.h>
#include <gmp.h>
#include "systemes.h"

// Échelonne le système selon l'algorithme de (Gauss-)Bareiss, NE GÈRE PAS ENCORE PROPREMENT LES PIVOTS NULS
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
		if(k == 0) mpz_set_si(d,1);
		else lit_coeff(d,s,k-1,k-1);
		//assert(mpz_sgn(d) != 0);		// Pas propre, détecte les pivots nuls
		lit_coeff(p,s,k,k);		
		// Pour chaque pivot, on va modifier tous les coeffs "en bas à droite" du pivot, y compris ceux dans le second membre
		for(int i = k+1;i < n;i++){
			for(int j = k+1;j < n+1;j++){	// j va jusqu'à n (inclus) car on souhaite aussi modifier le second membre
				lit_coeff(ij,s,i,j);
				lit_coeff(ik,s,i,k);
				lit_coeff(kj,s,k,j);
				mpz_mul(e,ij,p);
				mpz_submul(e,ik,kj);
				assert(mpz_divisible_p(e,d) != 0);	// DEBUG/VÉRIF
				mpz_divexact(ij,e,d);
				ecrit_coeff(s,i,j,ij);
			}
			mpz_set_si(e,0);
			ecrit_coeff(s,i,k,e);
		}
		// DEBUG, sera enlevé :
		//fprintf(stdout, "k = %d\n", k);
		//affiche_systeme(s);
	}
	mpz_clear(d);
	mpz_clear(p);
	mpz_clear(e);
	mpz_clear(ij);
	mpz_clear(ik);
	mpz_clear(kj);
	return;
}
