#include <stddef.h>

// On représentera les éléments de Z/pZ par des int (car p s'écrit avec au plus 30 bits)

// Opérations basiques sur les éléments de Z/pZ (p premier)
// Les opérandes sont supposées comprises entre 0 et p-1
long int zpz_add(long int x,long int y,long int p){
	long int res = x + y;
	res = res - p*(res >= p);
	return res;
}
long int zpz_sub(long int x,long int y,long int p){
	long int res = x - y;
	res = res + p*(res < 0);
	return res;
}
long int zpz_mul(long int x,long int y,long int p){
	__int128 lres = (__int128)x * y;	// On a besoin d'un long int ici, sinon les calculs cassent
	long int res = lres % p;	// Attention, % n'est pas exactement le modulo
	res = res + p*(res < 0);
	return res;
}

// Applique l'algorithme d'Euclide étendu à a et b (compris entre 0 et p-1)
// Calcule le PGCD de a et b, et une relation de Bézout : a*u + b*v = pgcd
// On a |u| <= a/(2*pgcd) et |v| <= b/(2*pgcd)
// Si l'on n'est pas intéressé par une partie du résultat, on peut mettre des pointeurs nuls pour pgcd, u ou v
void euclide_etendu(long int* pgcd,long int* u,long int* v,long int a,long int b){
	// Initialisation (avec juste 6 variables plutôt que 3 listes)
	long int r_a = a;
	long int r_b = b;
	long int u_a = 1;
	long int u_b = 0;
	long int v_a = 0;
	long int v_b = 1;
	long int q,r_new,u_new,v_new;
	while(r_b != 0){
		// Calculs
		q = r_a/r_b;
		r_new = r_a - r_b*q;
		u_new = u_a - u_b*q;
		v_new = v_a - v_b*q;
		// Passage à l'étape suivante
		r_a = r_b;
		u_a = u_b;
		v_a = v_b;
		r_b = r_new;
		u_b = u_new;
		v_b = v_new;
	}
	// Résultats
	if(pgcd != NULL) *pgcd = r_a;
	if(u != NULL) *u = u_a;
	if(v != NULL) *v = v_a;
	return;
}

// Calcule l'inverse de x modulo p
// On suppose que x et p sont premiers entre eux (c'est évidemment le cas si p est premier)
long int zpz_inv(long int x,long int p){
	long int u;		// On veut essentiellement renvoyer u, mais on veut quand même qu'il soit positif
	euclide_etendu(NULL,&u,NULL,x,p);
	if(u < 0) u += p;
	return u;
}



// REPRÉSENTATION SYMÉTRIQUE : Variante où les éléments de Z/pZ sont représentés par des int entre -p/2 et p/2
// Nécessite que p soit un nombre premier autre que 2, ce qui ne devrait pas trop poser de problème
long int zpzs_add(long int x,long int y,long int p){
	long int res = x+y;
	return res - p*(res > (p/2)) + p*(res < -(p/2));
}
long int zpzs_sub(long int x,long int y,long int p){
	long int res = x-y;
	return res - p*(res > (p/2)) + p*(res < -(p/2));
}
long int zpzs_mul(long int x,long int y,long int p){
	__int128 lres = (__int128)x * y;
	long int res = lres%p;
	return res - p*(res > (p/2)) + p*(res < -(p/2));
}

long int zpzs_inv(long int x,long int p){
	long int u;
	euclide_etendu(NULL,&u,NULL,x,p);
	// D'après le cours, on a |u| <= p/2, donc il n'y a rien à faire
	return u;
}



// Conversions entre les représentations "naturelle" et "symétrique"
long int zpzs_from_zpz(long int x,long int p){
	return x - (x > p/2)*p;
}
long int zpz_from_zpzs(long int x,long int p){
	return x + (x < 0)*p;
}
