CC=mpicc
CFLAGS+=-std=gnu99 -O3
LDLIBS=-lm
NP=4

.PHONY : heat
.PHONY : data
.PHONY : movie
.PHONY : clean

heat : heat.c
	${CC} -o $@ $< ${LDLIBS}

data :
	mpirun -np ${NP} ./heat

movie : heat data
	${MAKE} -j4 -C data movie

clean :
	rm -f heat && ${MAKE} -C data clean
