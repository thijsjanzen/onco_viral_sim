CFLAGS = -std=c++17 -O3 -Wall 
OBJS = main.o simulation.o node.o 

all: ${OBJS}
	g++ ${CFLAGS} ${OBJS} -o testmodel.out

main.o: ./Simulation/main.cpp
	g++ ${CFLAGS} -c ./Simulation/main.cpp

simulation.o: ./Simulation/simulation.cpp
	g++ ${CFLAGS} -c ./Simulation/simulation.cpp
  
node.o: ./Simulation/node.cpp
	g++ ${CFLAGS} -c ./Simulation/node.cpp

clean:
	rm -f testmodel.out ${OBJS}
