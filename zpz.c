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
	//assert((0 <= res) && (res < p));
	return res;
}
int zpz_sub(int x,int y,int p){
	int res = x - y;
	//if(res < 0) res += p;
	res = res + p*(res < 0);
	//assert((0 <= res) && (res < p));
	return res;
}
int zpz_mul(int x,int y,int p){
	long int lres = (long int)x * (long int)y;	// On a besoin d'un long int ici, sinon les calculs cassent
	int res = lres % p;	// Attention, % n'est pas exactement le modulo
	//if(res < 0) res += p;
	res = res + p*(res < 0);
	//assert((0 <= res) && (res < p));
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

// Génére un nombre premier aléatoire d'au plus b bits (en pratique, 30), donne sa valeur à p, et le renvoie aussi (en tant qu'int)
int genere_p(mpz_t p,gmp_randstate_t state,mp_bitcnt_t b){
	mpz_t k;
	mpz_init(k);
	mpz_urandomb(k,state,b);
	mpz_nextprime(p,k);
	int res = mpz_get_ui(p);
	mpz_clear(k);
	fprintf(stderr,"p = %d\n",res);		// DEBUG
	return res;
}

bool sol_egales(rationnel* sol1,rationnel* sol2,int n){
	bool res = true;
	for(int i = 0;i < n;i++){
		res = res && rat_comp(sol1[i],sol2[i]);;
		//rat_aff(sol1[i],stderr);		// DEBUG
		//rat_aff(sol2[i],stderr);		// DEBUG
	}
	return res;
}

// Calcule l'image réciproque de (x1,x2) (dans Z/(n1*n2)Z) par l'isomorphisme des restes chinois
// Écrit le résultat dans res, et nécessite d'avoir déjà établi une relation de Bézout : n1*u + n2*v == 1
void chinois(mpz_t res,mpz_t x1,mpz_t x2,mpz_t n1,mpz_t n2,mpz_t u,mpz_t v){
	mpz_t t1,t2,n1n2;
	mpz_init_set(t1,x1);
	mpz_init_set(t2,x2);
	mpz_init(n1n2);
	mpz_mul(n1n2,n1,n2);
	//mpz_addmul(t1,n2,v);	// Non
	//mpz_addmul(t2,n1,u);	// Non plus
	mpz_mul(t1,t1,n2);
	mpz_mul(t1,t1,v);
	mpz_mul(t2,t2,n1);
	mpz_mul(t2,t2,u);
	mpz_add(res,t1,t2);
	//if(mpz_cmp_ui(res,0) < 0) mpz_addmul(res,n1,n2);	// On veut un résultat entre O et n1*n2-1	// Inutile
	mpz_fdiv_r(res,res,n1n2);	// Faire ça en 2 fois, avant et après l'addition ???
	/*{		// TESTS
		// (0 <= x1 < n1) && (0 <= x2 < n2) (semble OK)
		assert(mpz_cmp_si(x1,0) >= 0);
		assert(mpz_cmp_si(x2,0) >= 0);
		assert(mpz_cmp(x1,n1) < 0);
		assert(mpz_cmp(x2,n2) < 0);
		// n1*u + n2*v == 1 (semble OK)
		mpz_t b;
		mpz_init(b);
		mpz_mul(b,n1,u);
		mpz_addmul(b,n2,v);
		assert(mpz_cmp_si(b,1) == 0);
		mpz_clear(b);
		// 0 <= res < n1*n2 (semble OK)
		mpz_mul(n1n2,n1,n2);
		assert((mpz_cmp_si(res,0) >= 0) && (mpz_cmp(res,n1n2) < 0));
		// (res mod n1 == x1) && (res mod n2 == x2) (semble OK)
		mpz_t r1,r2;
		mpz_init(r1);
		mpz_init(r2);
		mpz_fdiv_r(r1,res,n1);
		mpz_fdiv_r(r2,res,n2);
		mpz_out_str(stderr,10,r1);
		fputc('\n',stderr);
		mpz_out_str(stderr,10,x1);
		fputc('\n',stderr);
		mpz_out_str(stderr,10,r2);
		fputc('\n',stderr);
		mpz_out_str(stderr,10,x2);
		fputc('\n',stderr);
		assert(mpz_cmp(r1,x1) == 0);
		assert(mpz_cmp(r2,x2) == 0);
		mpz_clear(r1);
		mpz_clear(r2);
	}*/
	mpz_clear(t1);
	mpz_clear(t2);
	mpz_clear(n1n2);
	return;
}

// Une autre implémentation de l'algorithme d'Euclide étendu
// S'arrête dès que le reste est plus petit que sqrt(a), et utilise des mpz_t
// SEMBLE CASSÉ
void euclide_etendu_borne(mpz_t r,mpz_t v,mpz_t a,mpz_t b){
//void euclide_etendu_borne(mpz_t r,mpz_t u,mpz_t v,mpz_t a,mpz_t b){
	int nb_iter = 0;	// DEBUG
	// Initialisation
	mpz_t r_b,r_new;
	//mpz_t u_b,u_new;
	mpz_t v_b,v_new;
	mpz_t q,borne;
	mpz_set(r,a);
	mpz_init_set(r_b,b);
	//mpz_set_si(u,1);
	//mpz_init_set_si(u_b,0);
	mpz_set_si(v,0);
	mpz_init_set_si(v_b,1);
	mpz_init(r_new);
	//mpz_init(u_new);
	mpz_init(v_new);
	mpz_init(q);
	mpz_init(borne);
	mpz_sqrt(borne,a);
	//mpz_out_str(stderr,10,a); fputc('\n',stderr);		// DEBUG
	//mpz_out_str(stderr,10,borne); fputc('\n',stderr);	// DEBUG	// sqrt semble OK
	while(mpz_cmp(r_b,borne) > 0){		// Tant que r_b > sqrt(a)
		nb_iter++;	// DEBUG
		// Calculs
		mpz_fdiv_q(q,r,r_b);
		mpz_set(r_new,r);
		mpz_submul(r_new,r_b,q);
		//mpz_set(u_new,u);
		//mpz_submul(u_new,u_b,q);
		mpz_set(v_new,v);
		mpz_submul(v_new,v_b,q);
		// Passage à l'étape suivante
		mpz_set(r,r_b);
		//mpz_set(u,u_b);
		mpz_set(v,v_b);
		mpz_set(r_b,r_new);
		//mpz_set(u_b,u_new);
		mpz_set(v_b,v_new);
	}
	// ATTENTION : le résultat souhaité est en fait [r_b et v_b] (et non [r et v])
	// Fix douteux :
	mpz_set(r,r_b);
	mpz_set(v,v_b);
	// TODO : réécrire cette fonction convenablement
	
	// DEBUG 
	/*fprintf(stderr,"%d\n",nb_iter);
	mpz_out_str(stderr,10,r);
	fputc(' ',stderr);
	mpz_out_str(stderr,10,v);
	fputc('\n',stderr);*/
	// TEST OK sur l'exemple du cours (a = 257, b = 42, r = 5, v = -6)
	
	// Suppression des variables utilisées
	mpz_clear(r_b);
	mpz_clear(r_new);
	//mpz_clear(u_b);
	//mpz_clear(u_new);
	mpz_clear(v_b);
	mpz_clear(v_new);
	mpz_clear(q);
	mpz_clear(borne);
	return;
}

// Méthode modulaire de résolution du système [...]
void zpz(systeme* s,rationnel* sol,gmp_randstate_t state,mp_bitcnt_t b){
	// INITIALISATION
	int j = 1;		// DEBUG, nombre d'itérations effectuées
	int n = s -> n;
	syst_zpz sz;	// On va initialiser/copier et détruire ce système à chaque itération
	int* sol_zpz = malloc(n*sizeof(int));	// Emplacement de la solution dans Z/pZ que zpz_resol va calculer
	mpz_t* sol_zpz_m = malloc(n*sizeof(mpz_t));	// De même, pour la solution convertie en mpz
	mpz_t* sol_tmp = malloc(n*sizeof(mpz_t));	// Contient ce à partir de quoi on va essayer de reconstruire une solution
	mpz_t* sol_tmp_old = malloc(n*sizeof(mpz_t));	// Contient ce à partir de quoi on a essayé de reconstruire une solution à l'itération précédente
	rationnel* sol_old = malloc(n*sizeof(rationnel));	// Emplacement du candidat de solution de l'itération précédente
	mpz_t p_mpz,prod_old,prod,u,v;
	int p;		// Nombre premier "en cours d'utilisation", version machine
	mpz_init(p_mpz);		// Nombre premier "en cours d'utilisation", version GMP
	mpz_init(prod_old);		// Produit des nombres premiers précédemment utilisés (sans p_mpz)
	mpz_init(prod);			// Produit des nombres premiers précédemment utilisés (avec p_mpz) -> pendant une itération, prod == p*prod_old
	mpz_init(u);		// Pour la relation de Bézout, et pour la reconstruction modulaire de la solution (où il représente r, ce qui en fait un choix de nom douteux)
	mpz_init(v);		// Pour la relation de Bézout, et pour la reconstruction modulaire de la solution
	for(int i = 0;i < n;i++){
		mpz_init(sol_zpz_m[i]);
		mpz_init(sol_tmp[i]);
		mpz_init(sol_tmp_old[i]);
		rat_init(&(sol_old[i]));
	}
	// PREMIÈRE ITÉRATION (On verra plus tard si on peut la faire rentrer dans la boucle principale)
	// "Choix" d'un nombre premier
	p = genere_p(p_mpz,state,b);
	// Initialisation de prod et prod_old
	mpz_set(prod,p_mpz);
	mpz_set_ui(prod_old,1);
	// Calcul d'une solution dans Z/pZ
	init_copie_syst_zpz(&sz,s,p);
	zpz_resol(&sz,sol_zpz);
	// Conversion de la solution en mpz_t, directement dans sol_tmp puisqu'il n'y a pas de chinoiserie à faire
	for(int i = 0;i < n;i++){
		mpz_set_si(sol_tmp[i],sol_zpz[i]);
	}
	// Construction modulaire d'un candidat solution
	for(int i = 0;i < n;i++){
		euclide_etendu_borne(u,v,prod,sol_tmp[i]);	// Ici, u devrait s'appeler r
		if(mpz_cmp_si(v,0) < 0){	// v < 0
			mpz_neg(u,u);
			mpz_neg(v,v);	// De sorte que v soit positif
			rat_set_pq(sol[i],u,v);
		}else{		// v > 0
			rat_set_pq(sol[i],u,v);
		}
	}
	// MàJ de prod_old (on pourrait en fait l'initialiser ici, puisqu'on ne s'en sert pas avant)
	mpz_set(prod_old,p_mpz);
	// Destruction du syst_zpz, pour pouvoir le réutiliser (peu propre (?))
	detruit_syst_zpz(&sz);
	// BOUCLE PRINCIPALE
	while(!sol_egales(sol,sol_old,n) && (j < 500)){		// Condition en plus temporaire
		j++;	// DEBUG
		// "Choix" d'un nombre premier pour cette itération
		p = genere_p(p_mpz,state,b);
		// MàJ de prod
		mpz_mul(prod,prod,p_mpz);
		// Copie du précédent candidat solution sur l'emplacement de l'ancien
		for(int i = 0;i < n;i++){
			rat_set(sol_old[i],sol[i]);
			mpz_set(sol_tmp_old[i],sol_tmp[i]);
		}
		// Calcul d'une solution dans Z/pZ
		init_copie_syst_zpz(&sz,s,p);
		zpz_resol(&sz,sol_zpz);
		// Conversion de la nouvelle solution en mpz_t
		for(int i = 0;i < n;i++){
			mpz_set_si(sol_zpz_m[i],sol_zpz[i]);
		}
		// Calcul d'une relation de Bézout : u*p_mpz + v*prod_old == 1
		mpz_gcdext(NULL,u,v,p_mpz,prod_old);	// NULL -> on ne veut pas le PGCD (qui vaut 1)
		// Restes chinois avec les solutions précédentes
		for(int i = 0;i < n;i++){
			// sol_tmp[i] est dans Z/prodZ
			// sol_zpz_m[i] est dans Z/pZ
			// sol_tmp_old[i] est dans Z/prod_oldZ
			chinois(sol_tmp[i],sol_zpz_m[i],sol_tmp_old[i],p_mpz,prod_old,u,v);
		}
		// Construction modulaire d'un candidat solution
		for(int i = 0;i < n;i++){
			euclide_etendu_borne(u,v,prod,sol_tmp[i]);	// Ici, u devrait s'appeler r
			if(mpz_cmp_si(v,0) < 0){	// v < 0
				mpz_neg(u,u);
				mpz_neg(v,v);	// De sorte que v soit positif
				rat_set_pq(sol[i],u,v);
			}else{		// v > 0
				rat_set_pq(sol[i],u,v);
			}
		}
		// MàJ de prod_old
		mpz_mul(prod_old,prod_old,p_mpz);
		// Destruction du syst_zpz, pour pouvoir le réutiliser (peu propre (?))
		detruit_syst_zpz(&sz);
	}
	fprintf(stderr,"Nombre d'itérations : %d\n\n",j);	// DEBUG
	// SUPPRESSION DES OBJETS INUTILISÉS
	for(int i = 0;i < n;i++){
		mpz_clear(sol_zpz_m[i]);
		mpz_clear(sol_tmp[i]);
		mpz_clear(sol_tmp_old[i]);
		rat_clear(&(sol_old[i]));
	}
	mpz_clear(p_mpz);
	mpz_clear(prod_old);
	mpz_clear(prod);
	mpz_clear(u);
	mpz_clear(v);
	free(sol_zpz);
	free(sol_zpz_m);
	free(sol_tmp);
	free(sol_tmp_old);
	free(sol_old);
	return;
}
