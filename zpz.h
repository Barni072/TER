#ifndef ZPZ_H
#define ZPZ_H

long int zpz_add(long int x,long int y,long int p);
long int zpz_sub(long int x,long int y,long int p);
long int zpz_mul(long int x,long int y,long int p);
//void euclide_etendu(long int* pgcd,long int* u,long int* v,long int a,long int b);
long int zpz_inv(long int x,long int p);

long int zpzs_add(long int x,long int y,long int p);
long int zpzs_sub(long int x,long int y,long int p);
long int zpzs_mul(long int x,long int y,long int p);
long int zpzs_inv(long int x,long int p);

long int zpzs_from_zpz(long int x,long int p);
long int zpz_from_zpzs(long int x,long int p);

#endif
