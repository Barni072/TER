OBJS = TER.o rationnels.o systemes.o bareiss.o gauss_sys_rat.o io.o zpz.o mod.o mod_dets.o mod_thrd.o

TESTS = test.o rationnels.o systemes.o bareiss.o gauss_sys_rat.o io.o zpz.o mod.o mod_dets.o mod_thrd.o

CC = gcc

CFLAGS = -c -g -Wall -O2

all: $(OBJS)
	$(CC) $(OBJS) -o TER -lgmp

%.o: %.c
	$(CC) $(CFLAGS) $<

test: $(TESTS)
	$(CC) $(TESTS) -o test -lgmp

clean:
	rm -f *.o
	rm -f TER
	rm -f test
