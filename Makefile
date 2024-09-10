LIBRARIES=-lm -fopenmp 
OPTIM=-O3 -mavx2 -march=native  -Wall -Wextra -DUSE_NORM
DEBUG=
SRC="src"
VERBOSE=-DVERBOSE

CC=gcc

all: lib 

lib: bin/libadp2d.so

bin/libadp2d.so: bin/adp2d.o 
	${CC} -shared bin/adp2d.o  ${DEBUG} ${OPTIM} ${LIBRARIES} -o bin/libadp2d.so 
bin/adp2d.o: src/adp2d.c 
	${CC} -c src/adp2d.c  ${DEBUG} ${OPTIM} ${LIBRARIES} ${VERBOSE} -fPIC -o bin/adp2d.o 
clean:
	rm bin/*.so bin/*.o *.o driver test
