CFLAGS = -std=c++17 -O3 -Wall 
OBJS = main.o simulation.o setup.o node.o analysis.o

all: ${OBJS}
	g++ ${CFLAGS} ${OBJS} -o onco_model

setup.o: ./Simulation/setup.cpp
	g++ ${CFLAGS} -c ./Simulation/setup.cpp
	
analysis.o: ./Simulation/analysis.cpp
	g++ ${CFLAGS} -c ./Simulation/analysis.cpp	

main.o: ./console/main.cpp
	g++ ${CFLAGS} -c ./console/main.cpp

simulation.o: ./Simulation/simulation.cpp 
	g++ ${CFLAGS} -c ./Simulation/simulation.cpp 
  
node.o: ./Simulation/node.cpp
	g++ ${CFLAGS} -c ./Simulation/node.cpp

clean:
	rm -f onco_model ${OBJS}
