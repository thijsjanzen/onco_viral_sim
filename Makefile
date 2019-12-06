CFLAGS = -std=c++17 -O3 -Wall -I "../" 
OBJS = simulation.o node.o 

all: ${OBJS}
	g++ ${CFLAGS} ${OBJS} -o testmodel.out

simulation.o: ../Simulation/simulation.cpp
	g++ ${CFLAGS} -c ../Simulation/simulation.cpp
  
node.o: ../Simulation/node.cpp
	g++ ${CFLAGS} -c ../Simulation/node.cpp

clean:
	rm -f testmodel.out ${OBJS}
