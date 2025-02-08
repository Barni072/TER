#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <gmp.h>
#include "rationnels.h"
#include "systemes.h"



// Échelonne le système selon l'algorithme de (Gauss-)Bareiss, NE MODIFIE PAS ENCORE B, NE GÈRE PAS ENCORE PROPREMENT LES PIVOTS NULS
/*void bareiss(systeme* s){
	entier n = s -> taille;
	entier d;	// Représentera le pivot de l'étape précédente
	entier p;	// Représentera le pivot de l'étape courante
	for(entier k = 0;k < n-1;k++){
		if(k == 0) d = 1;
		else d = s->a[k-1][k-1];
		assert(d != 0);		// Pas propre, détecte les pivots nuls (on pourra plutôt permuter les lignes...)
		p = s->a[k][k];		
		// Pour chaque pivot, on va modifier tous les coeffs "en bas à droite" du pivot (MAL DIT)
		for(entier i = k+1;i < n;i++){
			for(entier j = k+1;j < n;j++){
				assert(((s->a[i][j])*p - (s->a[i][k])*(s->a[k][j])) % d == 0); 	// DEBUG/VERIF, a priori ce truc devrait toujours être nul
				s->a[i][j] = ((s->a[i][j])*p - (s->a[i][k])*(s->a[k][j])) / d;	// LA ligne importante, à comparer avec Gauss (faire une version Gauss-PasBareiss ?)
			}
			s->a[i][k] = 0;
			
			// MODIFIER B ICI
		}
		// DEBUG, sera enlevé :
		fprintf(stdout, "k = %d\n", k);
		affiche_systeme(s);
	}
	return;
}*/




int main(){	
	systeme s;
	init_systeme(&s,2,3);
	
	/*for(int i = 0;i<2*3;i++){
		assert(s.t[i] != NULL);
	}*/
	
	// Tests :
	int t_i[6] = {1,2,-3,4,7,8};
	mpz_t t_m[6];
	for(int i = 0;i < 6;i++){
		mpz_init_set_si(t_m[i],t_i[i]);
	}
	ecrit_coeff(&s,0,0,t_m[0]);
	ecrit_coeff(&s,0,1,t_m[1]);
	ecrit_coeff(&s,1,0,t_m[2]);
	ecrit_coeff(&s,1,1,t_m[3]);
	
	ecrit_coeff(&s,0,2,t_m[4]);
	ecrit_coeff(&s,1,2,t_m[5]);
	
	affiche_systeme(&s);
	
	// Test de sol_syst_echelonne -> OK
	/*assert(est_echelonne(s));
	rationnel* sol = malloc(3*sizeof(rationnel));
	sol_syst_echelonne(s,sol);
	for(entier i = 0;i < 3;i++){
		affiche_rationnel(sol[i]);
	}
	free(sol);*/
	
	//assert(!est_echelonne(s));	// Tests OK
	
	//bareiss(s);
	detruit_systeme(&s);
	return 0;
}		// TODO : tests depuis un fichier texte
