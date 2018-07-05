
CC = mpicc
CFLAGS = -g -Wall -Wextra -O3 -march=native -fopenmp

LD = gcc
LDFLAGS = -lm#-fopenmp 

all: gravitation clean  

run: gravitation
	./gravitation

gravitation: gravitation.o stopwatch.o bodies.o math_helper.o cluster.o eval.o vector.o
	$(CC) $(CFLAGS) gravitation.o stopwatch.o bodies.o math_helper.o cluster.o eval.o vector.o -o gravitation $(LDFLAGS)


gravitation.o: gravitation.c stopwatch.h bodies.h world.h cluster.h

stopwatch.o: stopwatch.c

bodies.o: bodies.c world.h math_helper.h

cluster.o: cluster.c bodies.h world.h math_helper.h #gravitation.h

math_helper.o: math_helper.c

eval.o: eval.h eval.c cluster.h

vector.o: vector.c

#world.o: world.h

clean:
	rm *.o