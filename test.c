#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <gmp.h>
#include <time.h>
#include "rationnels.h"
#include "systemes.h"
#include "bareiss.h"
#include "gauss_sys_rat.h"
#include "io.h"
#include "zpz.h"
#include "mod.h"
#include "mod_dets.h"
#include "mod_thrd.h"

int main(){
	/*int p = 4999;
	int q = 1307;
	// REPRÉSENTATION SYMÉTRIQUE DE Z/pZ
	assert(zpzs_add(0,1,3) == 1);
	assert(zpzs_add(2,2,5) == -1);
	assert(zpzs_mul(2,2,5) == -1);
	for(int x = 0;x < p;x++){
		for(int y = 0;y < p;y++){
			//fprintf(stderr,"%d %d\n",x,y);
			assert(zpz_add(x,y,p) == zpz_from_zpzs(zpzs_add(zpzs_from_zpz(x,p),zpzs_from_zpz(y,p),p),p));
			assert(zpz_sub(x,y,p) == zpz_from_zpzs(zpzs_sub(zpzs_from_zpz(x,p),zpzs_from_zpz(y,p),p),p));
			assert(zpz_mul(x,y,p) == zpz_from_zpzs(zpzs_mul(zpzs_from_zpz(x,p),zpzs_from_zpz(y,p),p),p));
		}
	}
	for(int x = -p/2;x < p/2;x++){
		for(int y = -p/2;y < p/2;y++){
			//fprintf(stderr,"%d %d\n",x,y);
			assert(zpzs_add(x,y,p) == zpzs_from_zpz(zpz_add(zpz_from_zpzs(x,p),zpz_from_zpzs(y,p),p),p));
			assert(zpzs_sub(x,y,p) == zpzs_from_zpz(zpz_sub(zpz_from_zpzs(x,p),zpz_from_zpzs(y,p),p),p));
			assert(zpzs_mul(x,y,p) == zpzs_from_zpz(zpz_mul(zpz_from_zpzs(x,p),zpz_from_zpzs(y,p),p),p));
		}
	}
	// RESTES CHINOIS AMÉLIORÉS
	mpz_t res1,res2,res3,x,y,u,v,m,n,mn;
	mpz_init(res1);
	mpz_init(res2);
	//mpz_init(res3);
	mpz_init(x);
	mpz_init(y);
	mpz_init(u);
	mpz_init(v);
	mpz_init(m);
	mpz_init(n);
	mpz_init(mn);
	mpz_set_si(m,p);
	mpz_set_si(n,q);
	mpz_mul(mn,m,n);
	for(int i = 0;i < p;i++){
		for(int j = 0;j < q;j++){
			mpz_set_si(x,i);
			mpz_set_si(y,j);
			mpz_gcdext(NULL,u,v,m,n);
			chinois_old(res1,x,y,m,n,u,v);
			chinois(res2,x,y,m,n,u,v);
			//chinois_sym(res3,x,y,m,n,u,v);
			assert(mpz_cmp(res1,res2) == 0);
		}
	}
	mpz_clear(res1);
	mpz_clear(res2);
	//mpz_clear(res3);
	mpz_clear(x);
	mpz_clear(y);
	mpz_clear(u);
	mpz_clear(v);
	mpz_clear(m);
	mpz_clear(n);
	mpz_clear(mn);*/
	fprintf(stderr,"%ld\n",CLOCKS_PER_SEC);
	return 0;
}
