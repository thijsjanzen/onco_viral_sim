CFLAGS = -std=c++17 -fno-strict-aliasing -ffast-math -march=native -O3 -Wall -I "../" -lpthread
OBJS = main.o simulation.o node.o 

all: ${OBJS}
	g++ ${CFLAGS} ${OBJS} -o testmodel.out

main.o: main.cpp
	g++ ${CFLAGS} -c main.cpp
  
simulation.o: simulation.cpp
	g++ ${CFLAGS} -c simulation.cpp
  
node.o: node.cpp
	g++ ${CFLAGS} -c node.cpp

clean:
	rm -f testmodel.out ${OBJS}
