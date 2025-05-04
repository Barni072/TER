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
	assert(zpzs_add(0,1,3) == 1);
	assert(zpzs_add(2,2,5) == -1);
	assert(zpzs_mul(2,2,5) == -1);
	int p = 4999;
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
	return 0;
}
