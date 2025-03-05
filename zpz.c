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
	//if(res >= p) res -= p;
	res = res - p*(res >= p);
	assert((0 <= res) && (res < p));
	return res;
}
int zpz_sub(int x,int y,int p){
	int res = x - y;
	//if(res < 0) res += p;
	res = res + p*(res < 0);
	assert((0 <= res) && (res < p));
	return res;
}
int zpz_mul(int x,int y,int p){
	long int lres = (long int)x * (long int)y;	// On a besoin d'un long int ici, sinon les calculs cassent
	int res = lres % p;	// Attention, % n'est pas excatement le modulo
	//if(res < 0) res += p;
	res = res + p*(res < 0);
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
// Suppose que x et p sont premiers entre eux (c'est évidemment le cas si p est premier)
int zpz_inv(int x,int p){
	int u;		// On veut essentiellement renvoyer u, mais on veut quand même qu'il soit positif
	euclide_etendu(NULL,&u,NULL,x,p);
	if(u < 0) u += p;
	return u;
}

// Ne serviront pas, on aura besoin de mpz_t pour ça...
/*int chinois(int x1,int x2,int n1,int n2){	// À TESTER
	int pgcd;	// DEBUG (?)
	int u,v;
	euclide_etendu(&pgcd,&u,&v,n1,n2);	// DEBUG
	//euclide_etendu(NULL,&u,&v,n1,n2);
	int res = (x1*n2*v + x2*n1*u)%(n1*n2);
	if(res >= 0) return res;
	else return res+(n1*n2);
	// Alternative : return res*(res >= 0) + (res+(n1*n2))*(res < 0);
}

int chinois_uv_connus(int x1,int x2,int n1,int n2,int u,int v){	// À TESTER
	int res = (x1*n2*v + x2*n1*u)%(n1*n2);
	if(res >= 0) return res;
	else return res+(n1*n2);
}*/

// Initialise le syst_zpz sdest, avec la copie du modulo p du systeme ssrc
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
		som = 0;
		for(int j = 0;j < n;j++){
			som = zpz_add(som, zpz_mul(sol[j], lit_coeff_zpz(s,i,j), p), p);
		}
		res = res && (som == lit_coeff_zpz(s,i,n));
	}
	return res;
}

// Échelonne le système s dans Z/pZ, avec l'algorithme de Gauss
// Les coefficients de s doivent être compris entre 0 et p-1
void zpz_gauss(syst_zpz* s){
	int n = s->n;
	int p = s->p;
	//int ivp,a,b,c;
	int ivp,a;
	// Échelonne le système
	for(int k = 0;k < n;k++){
		ivp = zpz_inv(lit_coeff_zpz(s,k,k),p);	// Inverse du pivot (le pivot est supposé non  nul)
		for(int i = k+1;i < n;i++){
			a = zpz_mul(ivp,lit_coeff_zpz(s,i,k),p);
			// (Ligne i) <- (Ligne i) - a * (Ligne k)
			for(int j = k+1;j < n+1;j++){
				//b = zpz_mul(a, lit_coeff_zpz(s,k,j), p);
				//c = zpz_sub(lit_coeff_zpz(s,i,j),b,p);
				//ecrit_coeff_zpz(s,i,j,c);
				// Version sans intermédiaire, à peu près illisible :
				ecrit_coeff_zpz(s,i,j, zpz_sub(lit_coeff_zpz(s,i,j), zpz_mul(a, lit_coeff_zpz(s,k,j) ,p) ,p));
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
		sol[k] = lit_coeff_zpz(s,k,n);
		for(int l = k+1;l < n;l++){
			sol[k] = zpz_sub(sol[k], zpz_mul(sol[l], lit_coeff_zpz(s,k,l), p), p);
		}
		sol[k] = zpz_mul(sol[k], zpz_inv(lit_coeff_zpz(s,k,k),p), p);
	}
	return;
}

// Écrit dans sol la solution du systeme s dans Z/pZ, en travaillant directement sur s
void zpz_resol(syst_zpz* s,int* sol){
	zpz_gauss(s);
	zpz_sol_syst_ech(sol,s);
	return;
}

// Fait pareil que la fonction précédente (en fait, appelle carrément la fonction prédcédente), mais dans un format que pthread_create accepte
void* zpz_resol_thrd(void* a_){
	zpz_args* a = (zpz_args*)a_;
	zpz_resol(&(a->s),a->sol);	
	return NULL;
}

// Génére un nombre premier aléatoire d'au plus b bits, donne sa valeur à p, et le renvoie aussi (en tant qu'int)
int genere_p(mpz_t p,gmp_randstate_t state,mp_bitcnt_t b){
	mpz_t k;
	mpz_init(k);
	mpz_urandomb(k,state,b);
	mpz_nextprime(p,k);
	int res = mpz_get_si(p);
	mpz_clear(k);
	return res;
}

bool sol_egales(rationnel* sol1,rationnel* sol2,int n){
	bool res = true;
	for(int i = 0;i < n;i++){
		res == res && rat_comp(sol1[i],sol2[i]);	// Douteux, à tester
	}
	return res;
}

// Méthode modulaire de résolution du système [...]		// PAS FINI
void zpz(systeme* s,rationnel* sol,gmp_randstate_t state,mp_bitcnt_t b){
	// Initialisation
	int n = s -> n;
	syst_zpz sz;	// On va initialiser/copier et détruire ce système à chaque itération
	int* sol_zpz = malloc(n*sizeof(int));	// Emplacement de la solution dans Z/pZ que zpz_resol va calculer
	mpz_t* sol_zpz_m = malloc(n*sizeof(mpz_t));	// De même, pour la solution convertie en mpz
	mpz_t* sol_tmp = malloc(n*sizeof(mpz_t));	// Contient ce à partir de quoi on va essayer de reconstruire une solution
	mpz_t* sol_tmp_old = malloc(n*sizeof(mpz_t));	// Contient ce à partir de quoi on a essayé de reconstruire une solution à l'itération précédente
	rationnel* sol_old = malloc(n*sizeof(rationnel));	// Emplacement du candidat de solution de l'itération précédente
	mpz_t p_mpz,prod_old,prod;
	int p;		// Nombre premier "en cours d'utilisation", version machine
	mpz_init(p_mpz);		// Nombre premier "en cours d'utilisation", version GMP
	mpz_init_set_ui(prod_old,1);		// Produit des nombres premiers précédemment utilisés (sans p_mpz)
	mpz_init_set_ui(prod,1);			// Produit des nombres premiers précédemment utilisés (avec p_mpz)
	for(int i = 0;i < n;i++){
		mpz_init(sol_zpz_m[i]);
		mpz_init(sol_tmp[i]);
		mpz_init(sol_tmp_old[i]);
		rat_init(&(sol_old[i]));
	}
	while(!sol_egales(sol,sol_old,n)){
		// MàJ de prod
		mpz_mul(prod,prod,p_mpz);
		// Copie du précédent candidat solution sur l'amplacement de l'ancien
		for(int i = 0;i < n;i++){
			rat_set(sol_old[i],sol[i]);
			mpz_set(sol_tmp_old[i],sol_tmp[i]);
		}
		// Calcul d'une solution dans un Z/pZ
		p = genere_p(p_mpz,state,b);
		init_copie_syst_zpz(&sz,s,p);
		zpz_resol(&sz,sol_zpz);
		// Conversion de la nouvelle solution en mpz_t
		for(int i = 0;i < n;i++){
			mpz_set_si(sol_zpz_m[i],sol_zpz[i]);
		}
		// Restes chinois avec les solutions précédentes
		// CALCULER U ET V (une seule fois)
		for(int i = 0;i < n;i++){
			// CHINOISER AVEC U ET V (calculer sol_tmp[i] à partir de sol_zpz_m[i] et sol_tmp_old[i])
		}
		// Construction modulaire d'un candidat solution
		// APPLIQUER LE COURS (Pour tout i, calculer sol[i] à partir de sol_tmp[i])
		
		// MàJ de prod_old
		mpz_mul(prod_old,prod_old,p_mpz);
	}
	// Suppression des objets utilisés
	for(int i = 0;i < n;i++){
		mpz_clear(sol_zpz_m[i]);
		mpz_clear(sol_tmp[i]);
		mpz_clear(sol_tmp_old[i]);
		rat_clear(&(sol_old[i]));
	}
	mpz_clear(p_mpz);
	mpz_clear(prod_old);
	mpz_clear(prod);
	free(sol_zpz);
	free(sol_zpz_m);
	free(sol_tmp);
	free(sol_tmp_old);
	free(sol_old);
	return;
}
