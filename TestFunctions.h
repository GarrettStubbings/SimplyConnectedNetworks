#include <iostream>

//Print contents of a vector function
void PrintVectorDouble(std::vector<double> &myVector){
    for (std::vector<double>::iterator i = myVector.begin();
        i != myVector.end(); ++i)
        std::cout << ' ' << *i;
    std::cout << '\n';
}

void PrintVectorInt(std::vector<int> &myVector){
    for (std::vector<int>::iterator i = myVector.begin();
        i != myVector.end(); ++i)
        std::cout << ' ' << *i;
    std::cout << '\n';
}


//print the the connections in the network

void CheckNeighbours(std::vector<Node> &Network) {

    for(const auto &v: Network){
        std::cout << "Node " << v.id << " has neighbours: ";
        for(auto x: v.Neighbours){
            std::cout <<  Network[x].id << " ";
            }
        std::cout << std::endl;
    }
}


