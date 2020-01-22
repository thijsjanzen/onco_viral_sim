CFLAGS = -std=c++17 -O3 -Wall -fno-strict-aliasing -ffast-math
OBJS = main.o simulation.o setup.o node.o 

all: ${OBJS}
	g++ ${CFLAGS} ${OBJS} -o testmodel.out

setup.o: ./Simulation/setup.cpp
	g++ ${CFLAGS} -c ./Simulation/setup.cpp

main.o: ./Simulation/main.cpp
	g++ ${CFLAGS} -c ./Simulation/main.cpp

simulation.o: ./Simulation/simulation.cpp 
	g++ ${CFLAGS} -c ./Simulation/simulation.cpp 
  
node.o: ./Simulation/node.cpp
	g++ ${CFLAGS} -c ./Simulation/node.cpp

clean:
	rm -f testmodel.out ${OBJS}
