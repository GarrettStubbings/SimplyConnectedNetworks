#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <assert.h>
#include <algorithm>
#include <unordered_map>
#include <utility>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

std::mt19937 g1(1); //Random number generator, mersenne twister
std::uniform_real_distribution<double> RNG(0.0,1.0); 

#include "Parameters.h"
#include "Files.h"
#include "DeficitType.h"
#include "Node.h"
#include "DataFile.h" 
#include "OutputDataFile.h" 
#include "SaveData.h"
#include "2Dvector.h"
#include "RandomNetwork.h" 
#include "tree.h" 
#include "GillespieTree.h"
#include "UpdateNode.h"
#include "MortalityNodes.h"
#include "ScaleFreeNetworkTree.h" 
#include "MortalityRules.h"
#include "OutputNetwork.h"
#include "AssortativeMixing.h"
#include "SmallWorld.h"
#include "SetUpNetwork.h"

int main(int argc, char *argv[]) {

    using namespace Parameters; //Using the parameters from cmd line
    SetParameters(argc, argv); //Set initial parameters from command line
    SetSimulate();
    //std::cout << "Parameters Loaded" << std::endl;
    
    //Network
    const int OriginalN = N; //some network types remove unconnected nodes, lowering N. Original N will not change
    std::vector<Node> Network;
    std::vector<int> MortalityNodes(nd);

    //set up network, if single topology set up network before the simulation
    g1.seed(Parameters::SingleSeed);
    SetUpNetwork(Network, MortalityNodes, OriginalN);
    //std::cout << "Single Network Created" << std::endl;

}

