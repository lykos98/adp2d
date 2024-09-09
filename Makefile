LIBRARIES=-lm -fopenmp 
OPTIM=-O4 -mavx2 -march=native  -Wall -Wextra -DUSE_NORM
DEBUG=
SRC="src"
VERBOSE=-DVERBOSE

CC=gcc

all: lib 

lib: bin/libdadac.so

bin/libdadac.so: bin/dadac.o 
	${CC} -shared bin/dadac.o  ${DEBUG} ${OPTIM} ${LIBRARIES} -o bin/libdadac.so 

bin/dadac.o: src/dadac.c 
	${CC} -c src/dadac.c  ${DEBUG} ${OPTIM} ${LIBRARIES} -fPIC -o bin/dadac.o 
clean:
	rm bin/*.so bin/*.o *.o driver test
