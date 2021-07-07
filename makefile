all: main parallelEvo backup TestLogNetwork evo entropyEvo testNetwork

main: main.cpp
	g++ -g -O3 main.cpp -o main

TestLogNetwork: TestLogNetwork.cpp
	g++-9 -g -O3 -DNDEBUG -std=c++1y TestLogNetwork.cpp -o TestLogNetwork\
        -I/opt/local/include/ -lgsl -lgslcblas

parallelEvo: parallelEvo.cpp
	g++-9 -g -O3 -DNDEBUG -std=c++1y parallelEvo.cpp -o parallelEvo -fopenmp\
        -I/opt/local/include/ -lgsl -lgslcblas

backup: main.cpp
	g++ -g -O3 -DNDEBUG -std=c++1y main.cpp -o main -lgsl -lgslcblas\
        -I/opt/local/include/

evo: evo.cpp
	g++-9 -g -O3 -DNDEBUG -std=c++1y evo.cpp -o evo\
        -I/opt/local/include/ -lgsl -lgslcblas

entropyEvo: entropyEvo.cpp
	g++-9 -g -O3 -DNDEBUG -std=c++1y entropyEvo.cpp -o entropyEvo\
        -I/opt/local/include/ -lgsl -lgslcblas

testNetwork: testNetwork.cpp
	g++-9 -g -O3 -DNDEBUG -std=c++1y testNetwork.cpp -o testNetwork\
        -I/opt/local/include/ -lgsl -lgslcblas

