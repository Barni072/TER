#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

// Taille maximale des systèmes :
#define N 1000

typedef int entier;

// On représente les nombres rationnels par des quotients d'entiers (pas nécessairement premiers entre eux (pour l'instant ?))
// On a besoin de ça, puisque le type float n'est pas assez précis pour une résolution exacte
struct s_rationnel{
	entier p;	// Numérateur
	entier q;	// Dénominateur (strictement positif)
};
typedef struct s_rationnel rationnel;

// Opérations de base sur les rationnels (somme, produit, opposé, inverse, produit avec un entier, division par un entier) :
// Noms de tous ces trucs à revoir ?
rationnel somme(rationnel x,rationnel y){
	rationnel res;
	res.p = (x.p * y.q) + (y.p * x.q);
	res.q = x.q * y.q;
	return res;
}
rationnel prod(rationnel x,rationnel y){
	rationnel res;
	res.p = x.p * y.p;
	res.q = x.q * y.q;
	return res;
}
rationnel opp(rationnel x){
	rationnel res;
	res.p = -x.p;
	res.q = x.q;
	return res;
}
rationnel inv(rationnel x){
	assert(x.p != 0);
	rationnel res;
	res.p = x.q;
	res.q = x.p;
	return res;
}
rationnel prod_avec_ent(rationnel x,entier n){
	rationnel res;
	res.p = x.p * n;
	res.q = x.q;
	return res;
}
rationnel div_par_ent(rationnel x,entier n){
	assert(n != 0);
	rationnel res;
	if(n > 0){
		res.p = x.p;
		res.q = x.q * n;
	}else{	// On veut que x.q reste positif
		res.p = -x.p;
		res.q = -x.q * n;
	}
	return res;
}

// Affiche le numérateur, le dénominateur et une approcimation flottante du rationnel x
void affiche_rationnel(rationnel x){
	fprintf(stdout,"Numérateur : %d\nDénominateur : %d\nApproximation : %.4f\n",x.p,x.q,((float)(x.p)/(float)(x.q)));
	return;
}

// Pour l'instant, on représente les systèmes linéaires par une matrice et un vecteur, contenus dans des "tableaux assez grands"
struct s_systeme{
	entier taille;	// Supposé <= N
	entier a[N][N];
	entier b[N];
};
typedef struct s_systeme systeme;

// Affiche un système dans le terminal, pour l'instant de façon matricielle et très moche
void affiche_systeme(systeme* s){
	entier n = s -> taille;
	//fprintf(stdout,"Taille : %d\n",n);	 // Cette information n'est pas très utile
	for(entier i = 0;i < n;i++){
		// Affichage d'une ligne de a
		for(entier j = 0;j < n;j++){
			fprintf(stdout,"%d ",s -> a[i][j]);
		}
		// Affichage d'un coeff de b, puis passage à la ligne suivante
		fprintf(stdout,"  %d\n",s -> b[i]);
	}
	return;
}

// Renvoie vrai ssi s est échelonné
bool est_echelonne(systeme* s){
	bool res = true;
	entier n = s -> taille;
	for(entier j = 0;j < n;j++){
		for(entier i = j+1;i < n;i++){
			res = res && (s->a[i][j] == 0);
			//fprintf(stdout,"%d %d %d\n",i,j,s->a[i][j]);	// Debug
		}
	}
	return res;
}

// "Remonte" le système s (supposé échelonné) et calcule sa solution
// Le tableau sol sert d'emplacement où la solution sera écrite (mal dit)
// DANGER SI PIVOTS NULS
void sol_syst_echelonne(systeme* s, rationnel* sol){
	assert(est_echelonne(s));
	entier n = s -> taille;
	for(entier k = n-1;k >= 0;k--){	// Calcul de chaque coordonnée de la solution :
		// Initialisation
		sol[k].p = s->b[k];
		sol[k].q = 1;
		for(entier l = k+1;l < n;l++){	// Ajoute les "contributions des lignes suivantes", une par une (pas clair) :
			rationnel x = prod_avec_ent(sol[l],s->a[k][l]);
			sol[k] = somme(sol[k],opp(x));
		}
		// Division par le pivot :
		sol[k] = div_par_ent(sol[k],s->a[k][k]);
	}
	return;
}

// Échelonne le système selon l'algorithme de (Gauss-)Bareiss, NE MODIFIE PAS ENCORE B, NE GÈRE PAS ENCORE PROPREMENT LES PIVOTS NULS
void bareiss(systeme* s){
	entier n = s -> taille;
	entier d;	// Représentera les pivots
	for(entier k = 0;k < n-1;k++){
		if(k == 0) d = 1;
		else d = s->a[k-1][k-1];
		assert(d != 0);		// Pas propre, détecte les pivots nuls (on pourra plutôt permuter les lignes...)
				
		// Pour chaque pivot, on va modifier tous les coeffs "en bas à droite" du pivot (MAL DIT)
		for(entier i = k+1;i < n;i++){
			for(entier j = k+1;j < n;j++){
				
				s->a[i][j] = ((s->a[i][j])*(s->a[k][k])-(s->a[i][k])*(s->a[k][j]))/d;	// LA ligne importante, à comparer avec Gauss
				assert(((s->a[i][j])*(s->a[k][k])-(s->a[i][k])*(s->a[k][j]))%d == 0); 	// DEBUG/VERIF, a priori ce truc devrait toujours être nul
			}
			s->a[i][k] = 0;
			
			// MODIFIER B ICI
		}
		// DEBUG, sera enlevé :
		fprintf(stdout, "k = %d\n", k);
		affiche_systeme(s);
	}
	return;
}




int main(){
	/*rationnel x,y;
	x.p = 2;
	x.q = 3;
	y.p = 4;
	y.q = 5;
	affiche_rationnel(x);
	affiche_rationnel(y);
	affiche_rationnel(somme(x,y));
	affiche_rationnel(prod(x,y));
	affiche_rationnel(opp(x));
	affiche_rationnel(inv(x));
	affiche_rationnel(inv(y));
	affiche_rationnel(prod_avec_ent(x,0));
	affiche_rationnel(prod_avec_ent(y,-4));
	affiche_rationnel(div_par_ent(x,3));
	affiche_rationnel(div_par_ent(y,-2));*/		// TESTS OK
	
	systeme* s = malloc(sizeof(systeme));
	assert(s != NULL);
	// Tests :
	s -> taille = 2;
	s -> a[0][0] = 1;
	s -> a[0][1] = 2;
	//s -> a[0][2] = 3;
	s -> a[1][0] = 4;
	s -> a[1][1] = 5;
	//s -> a[1][2] = 6;
	//s -> a[2][0] = 7;
	//s -> a[2][1] = 8;
	//s -> a[2][2] = 9;
	s -> b[0] = 10;
	s -> b[1] = 11;
	//s -> b[2] = 12;
	affiche_systeme(s);
	
	// Test de sol_syst_echelonne -> OK
	/*assert(est_echelonne(s));
	rationnel* sol = malloc(3*sizeof(rationnel));
	sol_syst_echelonne(s,sol);
	for(entier i = 0;i < 3;i++){
		affiche_rationnel(sol[i]);
	}
	free(sol);*/
	
	//assert(!est_echelonne(s));	// Tests OK
	
	bareiss(s);
	free(s);
	return 0;
}
