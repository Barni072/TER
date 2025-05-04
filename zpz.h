#ifndef ZPZ_H
#define ZPZ_H

int zpz_add(int x,int y,int p);
int zpz_sub(int x,int y,int p);
int zpz_mul(int x,int y,int p);
//void euclide_etendu(int* pgcd,int* u,int* v,int a,int b);
int zpz_inv(int x,int p);

int zpzs_add(int x,int y,int p);
int zpzs_sub(int x,int y,int p);
int zpzs_mul(int x,int y,int p);
int zpzs_inv(int x,int p);

int zpzs_from_zpz(int x,int p);
int zpz_from_zpzs(int x,int p);

#endif
