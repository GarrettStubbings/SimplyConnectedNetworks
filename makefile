all: main testNetwork clusterMain clusterTestNetwork

main: main.cpp
	g++-9 -g -O3 main.cpp -o main\
        -I/opt/local/include/ -lgsl -lgslcblas

testNetwork: testNetwork.cpp
	g++-9 -g -O3 -DNDEBUG -std=c++1y testNetwork.cpp -o testNetwork\
        -I/opt/local/include/ -lgsl -lgslcblas

clusterMain: main.cpp
	g++ -g -O3 main.cpp -o main\
        -I/opt/local/include/ -lgsl -lgslcblas

clusterTestNetwork: testNetwork.cpp
	g++ -g -O3 -DNDEBUG -std=c++1y testNetwork.cpp -o testNetwork\
        -I/opt/local/include/ -lgsl -lgslcblas

clean:
	rm main testNetwork
