OBJS = TER.o rationnels.o systemes.o bareiss.o gauss_sys_rat.o io.o zpz.o zpz_thrd.o

CC = gcc

CFLAGS = -c -g

all: $(OBJS)
	gcc $(OBJS) -o TER -lgmp

%.o: %.c
	$(CC) $(CFLAGS) -c $<
clean:
	rm -f *.o
	rm -f TER
