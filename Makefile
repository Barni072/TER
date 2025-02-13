all: TER.c rationnels.c rationnels.h systemes.c systemes.h
	gcc -c rationnels.c -lgmp
	gcc -c systemes.c -lgmp
	gcc -c bareiss.c -lgmp
	gcc -c gauss_sys_rat.c -lgmp
	gcc -c io.c -lgmp
	gcc -c TER.c -lgmp
	gcc rationnels.o systemes.o bareiss.o gauss_sys_rat.o io.o TER.o -o TER -lgmp
clean:
	rm -f *.o
	rm -f TER
