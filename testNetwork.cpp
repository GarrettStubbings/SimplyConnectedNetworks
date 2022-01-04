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
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext.hpp>
#include <boost/range/irange.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/exterior_property.hpp>
#include <boost/graph/clustering_coefficient.hpp>
#include <boost/graph/undirected_graph.hpp>

#include <gsl/gsl_histogram.h>

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
#include "SuperHubby.h"
#include "SetUpNetwork.h"
#include "Evaluations.h"
#include "ImportNetwork.h"
#include "TestFunctions.h"
#include "UpdateFI.h"
#include "MortalityCondition.h"
#include "RewireNetwork.h"

int main(int argc, char *argv[]) {
//set up network, if single topology set up network before the simulation

    using namespace Parameters; //Using the parameters from cmd lin
    SetDefaultParameters();
    SetReducedParameters(argc, argv);
    //std::cout << "\n\nFolder:" << Folder << "\n\n";

    //std::cout << "Parameters Loaded" << std::endl;
    g1.seed(SingleSeed);
    
    //Network
    const int OriginalN = N; 

    std::vector<Node> Network;
    
    //Hold all death ages found
    //std::cout << "Number: " << Number << "\n" ;
    std::vector<double> DeathAges;
    DeathAges.reserve(Number);
    ImportEdgeListNetwork(Network,
        Folder + initialDistribution + ".csv", N);

    //std::cout << "Single Network Created" << std::endl;

    int totalNodes = Network.size();

    //Healthy aging vector for each run
    std::vector<double> healthyAgingVector;
    healthyAgingVector.reserve(N);
    
    // Healthy Aging Indices of each individual on a per Network Basis
    std::vector<double> healthyAgingNVector;
    healthyAgingNVector.reserve(N);
   
    std::vector<double> HANormVector;
    HANormVector.reserve(N);
  
    double prefactor = normalizingFactor(Network, mu);
   
    std::vector<std::vector<double>> populationFIs;

    double timeScale = 0.00183;
    int oldest = 1;

    std::vector<DeficitVal> DeficitsValues;

    //std::cout << "\n";
    //loop through runs to average data
    //std::cout << "Running Starts" << std::endl;
    for(int run = 1; run <= Number; run++) {
        
        //Set up network for every seed if not single topology 
        EventTree<double> tree(N, 1.0); 
        double TotalRate = N; //sum of rates
        double mortalityRate = 0.1;
        double Time = 0; //current Time (unscaled age)

        double oldTime = 0.0; // For calculating timestep in HA
        double dt = 0.0; // Timestep for HealthyAging
        double FI = 0.0; // Frailty Of Individual
        double HA = 0; //Health aging
        std::vector<double> FIVector;
        FIVector.emplace_back(0);
        bool dead = false;
        
        int year = 0;

        int numEvents = 0;

        DeficitsValues.emplace_back(0, 0, 0);

        //evolve in time until Mortality occurs	
        int dummy = 0;
        while(not dead) {
            dummy++;
            numEvents++;
            //std::cout << TotalRate << "\n";
            //Find rate to be performed 
            int Index = FindRate(tree, TotalRate - mortalityRate);
            //perform rate and increase time
            UpdateNode(Network[Index], Time, TotalRate, DeficitsValues); 
             //update the local frailty after the node is modified
            UpdateLocalFrailty(Network, Network[Index]);
            
            // record the frailty index yearly
            int yearsPassed = int(Time/timeScale) - year;
            for (int gapYear = 0; gapYear <= yearsPassed; gapYear++){
                FIVector.emplace_back(FI);
                year++;
            }

            //calculate new rates for the node and its neighbours
            CalculateRates(Network, Network[Index], tree, Index, TotalRate); 
            
            dt = oldTime-Time;
            updateFrailty(FI, Network[Index], totalNodes);
            updateHealthyAging(HA, FI, dt);
            // New mortality Condition:
            //std::cout << "Mortality Rate: " << mortalityRate <<
            //    ", Total Rate: " << TotalRate << "\n";
            updateMortalityRate(mortalityRate, TotalRate,
                Network[Index], gammaD, mu, prefactor);


            //std::cout << "Mortality Rate: " << mortalityRate << 
            //    "Total Rate : " << TotalRate << "\n";
            evaluateMortality(dead, TotalRate, mortalityRate);

            //if(dead == true) break;
            //if(Time > 0.3) break;
            //evaluate mortality
            //if(Mortality(Network, MortalityNodes) == 1) break;    
            //if(isnan(TotalRate)) break;

            oldTime = Time;
        }

        //record healthy aging Data
        healthyAgingVector.emplace_back(HA);
        // record healthy aging norm
        HANormVector.emplace_back(HA/Time);
        //record death age data
        DeathAges.emplace_back(Time);
        // Record FI History
        populationFIs.emplace_back(FIVector);
        if (year > oldest){
            oldest = year;
        }
        
        //Reset Rates and fraility if it is single topology
        
        for(auto &x: Network) x.Reset();  
      

        if (run%10000 == 0)
            std::cout << "Run " << run << std::endl;
        
        //DeficitsValues.clear();

        
    }
    
    //std::cout << "Outputing final data" << std::endl;
    OutputMeans(DeathAges, OriginalN,"DeathAges");
    OutputMeans(healthyAgingVector, OriginalN, "QALY");
    OutputMeans(HANormVector, OriginalN, "HealthyAging");
    std::cout << "Oldest: " << oldest << "\n";
    extendFIs(populationFIs, oldest + 1);
    Output2d(populationFIs, OriginalN, "PopulationFI");
    //std::cout << "Average Death Age: " << mean(DeathAges) << "\n";
    //std::cout << "Average Healthy Aging: " << mean(healthyAgingVector) << "\n";
    //std::cout << "Average Normalized Healthy Aging: " << mean(
    //    HANormVector) << "\n";
    //std::cout << "Finished" << std::endl;

   
}

