#include <assert.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <gmp.h>
#include "systemes.h"
#include "zpz.h"

// On représentera les éléments de Z/pZ par des int (car p s'écrit avec au plus 30 bits)

// Opérations basiques sur les éléments de Z/pZ (p premier)
// Les opérandes sont supposées comprises entre 0 et p-1
int zpz_add(int x,int y,int p){
	int res = x + y;
	if(res >= p) res -= p;
	//Alternative : res = res - p*(res >= p);		// On fera ça plus tard
	assert((0 <= res) && (res < p));
	return res;
}
int zpz_sub(int x,int y,int p){
	int res = x - y;
	if(res < 0) res += p;
	//Alternative : res = res + p*(res <= 0);		// Plus tard
	assert((0 <= res) && (res < p));
	return res;
}
int zpz_mul(int x,int y,int p){
	long int lres = (long int)x * (long int)y;	// On a besoin d'un long int ici, sinon les calculs cassent
	int res = lres % p;	// Attention, % n'est pas excatement le modulo
	if(res < 0) res += p;	// Au secours
	//fprintf(stderr,"%d %d %d %ld %d\n",x,y,res,lres,p);
	assert((0 <= res) && (res < p));
	return res;
}

// Applique l'algorithme d'Euclide étendu à a et b (compris entre 0 et p-1)
// Calcule le PGCD de a et b, et une relation de Bézout : a*u + b*v = pgcd
// On a |u| <= a/(2*pgcd) et |v| <= b/(2*pgcd)
// Si l'on n'est pas intéressé par une partie du résultat, on peut mettre des pointeurs nuls pour pgcd, u ou v
void euclide_etendu(int* pgcd,int* u,int* v,int a,int b){
	// Initialisation (avec juste 6 variables plutôt que 3 listes)
	int r_a = a;
	int r_b = b;
	int u_a = 1;
	int u_b = 0;
	int v_a = 0;
	int v_b = 1;
	int q,r_new,u_new,v_new;
	while(r_b != 0){
		// Calculs		(Dans Z OK ?)
		q = r_a/r_b;	// Division entière OK ?
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

// Calcule le PGCD de a et b
int pgcd(int a,int b){
	int res;
	euclide_etendu(&res,NULL,NULL,a,b);
	return res;
}

// Calcule l'inverse de x modulo p
int zpz_inv(int x,int p){
	//int pgcd;	// DEBUG, à terme on mettre un NULL à sa place
	int u;		// On veut essentiellement renvoyer u, mais on veut quand même qu'il soit positif
	//euclide_etendu(&pgcd,&u,NULL,x,p);
	euclide_etendu(NULL,&u,NULL,x,p);
	//assert(pgcd == 1);	// DEBUG
	if(u < 0) u += p;
	//assert((u >= 0) && (u < p));	// DEBUG
	//fprintf(stderr,"%d %d %d\n",x,u,p);	// DEBUG
	//assert(zpz_mul(x,u,p) == 1);	// DEBUG
	return u;
}

int chinois(int x1,int x2,int p1,int p2){	// À TESTER
	int pgcd;	// DEBUG (?)
	int u,v;
	euclide_etendu(&pgcd,&u,&v,p1,p2);	// DEBUG
	//euclide_etendu(NULL,&u,&v,x1,x2);
	int res = (x1*p2*v + x2*p1*u)%(p1*p2);
	if(res >= 0) return res;
	else return res+(p1*p2);
}

int chinois_uv_connus(int x1,int x2,int p1,int p2,int u,int v){	// À TESTER
	int res = (x1*p2*v + x2*p1*u)%(p1*p2);
	if(res >= 0) return res;
	else return res+(p1*p2);
}

void init_copie_syst_zpz(syst_zpz* sdest,systeme* ssrc,int p){
	int n = ssrc->n;
	int m = ssrc->m;
	mpz_t pm,c;
	mpz_init_set_si(pm,p);
	mpz_init(c);
	sdest -> n = n;
	sdest -> m = m;
	sdest -> p = p;
	sdest -> t = malloc(sizeof(int)*n*m);
	for(int k = 0;k < n*m;k++){
		mpz_fdiv_r(c,ssrc->t[k],pm);
		sdest->t[k] = mpz_get_ui(c);
		// Alternative (devrait suffire, si j'ai bien compris la documentation) (il doit même exister un moyen de faire ça sans déclarer/initialiser c)
		//sdest->t[k] = mpz_fdiv_r(c,ssrc->t[k],p);
	}
	mpz_clear(pm);
	mpz_clear(c);
	return;
}

void detruit_syst_zpz(syst_zpz* s){
	free(s->t);
	return;
}

// Renvoie le coefficient (i,j) de s
int lit_coeff_zpz(syst_zpz* s,int i,int j){
	return s->t[i*(s->m) + j];
}

// Remplace le coefficient (i,j) de s par n
// n est supposé compris entre 0 et (s->p - 1)
void ecrit_coeff_zpz(syst_zpz* s,int i,int j,int n){
	s->t[i*(s->m) + j] = n;
	return;
}

void affiche_syst_zpz(syst_zpz* s,FILE* f){
	//fprintf(f,"p = %d\n",s->p);	// Inutile
	int n = s->n;
	for(int i = 0;i < n;i++){
		// Affichage d'une ligne de coeffs
		for(int j = 0;j < n;j++){
			fprintf(f,"%d ",lit_coeff_zpz(s,i,j));
		}
		// Affichage d'un coeff du second membre
		fprintf(f,"  %d\n",lit_coeff_zpz(s,i,n));
	}
	return;
}

// Vérifie si sol est une solution de s dans Z/pZ
// Les coefficients de s doivent être compris entre 0 et p-1
bool verif_sol_zpz(syst_zpz* s,int* sol){
	int n = s->n;
	int p = s->p;
	bool res = true;
	int ij,som;
	for(int i = 0; i < n;i++){
		//mpz_set_si(som,0);
		som = 0;
		for(int j = 0;j < n;j++){
			//lit_coeff(ij,s,i,j);
			//zpz_mul(tmp,ij,sol[j],p);
			//zpz_add(som,som,tmp,p);
			som = zpz_add(som, zpz_mul(sol[j], lit_coeff_zpz(s,i,j), p), p);
		}
		//lit_coeff(tmp,s,i,n);
		//res = res && (mpz_cmp(som,tmp) == 0);
		res = res && (som == lit_coeff_zpz(s,i,n));
	}
	return res;
}

// Échelonne le système s dans Z/pZ, avec l'algorithme de Gauss
// Les coefficients de s doivent être compris entre 0 et p-1
void zpz_gauss(syst_zpz* s){
	int n = s->n;
	int p = s->p;
	int ivp,a,b,c;
	// Échelonne le système
	for(int k = 0;k < n;k++){
		ivp = zpz_inv(lit_coeff_zpz(s,k,k),p);	// Inverse du pivot (le pivot est supposé non  nul)
		for(int i = k+1;i < n;i++){
			a = zpz_mul(ivp,lit_coeff_zpz(s,i,k),p);
			// (Ligne i) <- (Ligne i) - a * (Ligne k)
			for(int j = k+1;j < n+1;j++){
				b = zpz_mul(a, lit_coeff_zpz(s,k,j), p);
				c = zpz_sub(lit_coeff_zpz(s,i,j),b,p);
				ecrit_coeff_zpz(s,i,j,c);
				// Version sans intermédiaire, à peu près illisible :
				//ecrit_coeff_zpz(s,i,j, zpz_sub(lit_coeff_zpz(s,i,j), zpz_mul(a, lit_coeff_zpz(s,k,j) ,p) ,p));
			}
			ecrit_coeff_zpz(s,i,k,0);
		}
	}
	return;
}

// Copie dans sol la solution du système préalablement échelonné s
// Les coefficients de s doivent être compris entre 0 et p-1
void zpz_sol_syst_ech(int* sol,syst_zpz* s){
	int n = s->n;
	int p = s->p;
	for(int k = n-1;k >= 0;k--){
		//lit_coeff(sol[k],s,k,n);
		sol[k] = lit_coeff_zpz(s,k,n);
		for(int l = k+1;l < n;l++){
			//lit_coeff(a,s,k,l);
			//zpz_mul(b,sol[l],a,p);
			//zpz_sub(sol[k],sol[k],b,p);
			sol[k] = zpz_sub(sol[k], zpz_mul(sol[l], lit_coeff_zpz(s,k,l), p), p);
		}
		//lit_coeff(a,s,k,k);
		//mpz_invert(b,a,p);
		//zpz_mul(sol[k],sol[k],b,p);
		sol[k] = zpz_mul(sol[k], zpz_inv(lit_coeff_zpz(s,k,k),p), p);
	}
	return;
}

// Écrit dans sol la solution du systeme s dans Z/pZ, en travaillant directement sur s
void zpz(int* sol,syst_zpz* s){
	zpz_gauss(s);
	zpz_sol_syst_ech(sol,s);
	return;
}

// Fait pareil fonction que la précédente (en fait, appelle carrément la fonction prédcédente), mais dans un format que pthread_create accepte
void* zpz_thrd(void* a_){
	zpz_args* a = (zpz_args*)a_;
	zpz((a->sol),&(a->s));	
	return NULL;
}
