/*
Code for C++ implementation of bertotti configuration algorithm
*/

// define the graph stuff
typedef boost::adjacency_list <boost::vecS, boost::vecS, 
    boost::undirectedS> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor MyVertex;

// functions for keeping track of airball edges (edges which shouldnt exist
// according to P(k',k)
void removeAirball(std::vector<int> &airballs, int edgeID){
    std::vector<int>::iterator itr = std::find(airballs.begin(),
        airballs.end(), edgeID);
    if (itr != airballs.end()){
        // find the index of that edge in the airballs list
        int index = std::distance(airballs.begin(), itr);
        // remove it
        airballs.erase(airballs.begin() + index);
    }
    else {
        std::cout << "Airball Does not Exist!\n";
    }

}// end remove airball

// Add an edge which shouldnt exist
void addAirball(std::vector<int> &airballs, int edgeID){
    bool skip = false;
    if (edgeID > airballs[airballs.size() - 1]){
        airballs.emplace_back(edgeID);
        skip = true;
    }
    if (not skip){
        for (int i = 0; i < airballs.size(); i++){
            if (edgeID < airballs[i]){
                airballs.insert(airballs.begin() + i, edgeID); 
                break;
            }
        }// end loop over all airball IDs
    }
}// end add airball

// generate a random network with P(k)
Graph randomInitialNetwork(std::vector<double> degreeDistribution, int &N,
    std::vector<int> &nodeDegrees, double avgDeg, bool &initFailure){
    
    //std::cout << "Beginning Network Initialization\n";
    // initialize graph
    using namespace boost;
    int maxDegree = degreeDistribution.size();
    //std::cout << "Maximum degree: " << maxDegree << "\n";


    
    int degreeCount = 0;
    int degreeSum = 0;
    int numNodesSoFar = 0;
    std::vector<int> numNodes (maxDegree, 0);
    int deg = 0;

    //std::cout << "Total number of nodes: " << totalNumberOfNodes << "\n";

    Graph G;
    for (int i = 0; i < N; i++) add_vertex(G);
   
    
    // sample degree Distribution to fill remaining nodes (due to rounding down
    // when calculating N*P(k)

    // generate the sampler
    std::discrete_distribution<int> sampleDist (degreeDistribution.begin(),
        degreeDistribution.end());
    

    //std::vector<int> numNodes (maxDegree, 0);

    // generate node degrees from the sampler
    //std::cout << "Rounding Lead to Only " << numNodesSoFar <<
    //    " Nodes Generated of the " << N << " Total, sample remaining\n";
    //std::cout << "Degree Sum So Far: " << degreeSum << "\n";
    //PrintVectorInt(numNodes);
    bool withinTolerance = false;
    while (not withinTolerance){
        // Reset number of nodes
        std::fill(numNodes.begin(), numNodes.end(), 0); 
        degreeSum = 0;
        numNodesSoFar = 0;
        for (int i=numNodesSoFar; i < N; i++){
            deg = sampleDist(g1);
            //std::cout << "Adding Node of degree " << degree << "\n";
            ++numNodes[deg];
            degreeSum += deg;
        }// end degree assignemnt
        
        // Assert that the average degree is within 5% of desired
        double sampledAvgDeg = double(degreeSum) / double(N);
        double percentError = abs(avgDeg - sampledAvgDeg) / avgDeg;
        //std::cout << "Sampled Avg Deg: " << sampledAvgDeg << ", AvgDeg: " <<
        //   avgDeg << ", Percent error: " << percentError << "\n";
        withinTolerance = percentError < 0.05;

    }// end Avg K within 5% insurance + degree assignment
    
    //std::cout << "Updated degree Sum: " << degreeSum << "\n";
    //PrintVectorInt(numNodes);


    //PrintVectorInt(numNodes);
    //std::cout << "Degree Sum: " << degreeSum << "\n";

    //std::cout << "Distribution Sampled, Generating Stubs\n";

    //std::cout << "Sampled Degree Counts are: ";
    //PrintVectorInt(numNodes);
    //std::cout << "MAx Degree : " << maxDegree << "\n";

    // generate stubs
    int stubsSum = 0;
    for (int i = 0; i < maxDegree; i++){
        // fill degrees in ascending order
        deg = i; // Descending order case: maxDegree - i;
        // number of stubs
        int numStubs = numNodes[deg];
        // add the correct number of stubs for the given degree
        while (numStubs > 0){
            nodeDegrees.emplace_back(deg);
            stubsSum += 1;
            numStubs -= 1;
        }// end stub filling

    }// end stub generation
    //std::cout << "Node Degrees (Pre-Increment etc)\n";
    //PrintVectorInt(nodeDegrees);

    //std::cout << "Stubs Generated, checking for odd degree\n";
    //std::cout << "Stubs Sum: " << stubsSum << "\n";
    //std::cout << "Number of Nodes: " << nodeDegrees.size() << "\n";
    /*
    int realDegreeSum = 0;
    for (auto k: nodeDegrees){
        realDegreeSum += k;
    }
    if (realDegreeSum != degreeSum){
        //std::cout << "Actual Degree Sum is: " << realDegreeSum << "\n";
        //std::cout << "Number of Nodes: " << nodeDegrees.size() << "\n";
        //std::cout << "Node Degrees: "; 
        //PrintVectorInt(nodeDegrees);
    }// PRoblems
    */
    // assert that the sum of node degrees is even
    if (degreeSum %2 == 1) {
        //std::cout << "Odd total Degree: ";
        //std::cout << degreeSum << " Total Degree\n";
        //int doubleCheck = 0;
        //for (auto x: nodeDegrees) doubleCheck += x;
        //std::cout << "Double Check Has total degree: " << doubleCheck << "\n";
        int incrementID = int(RNG(g1) * nodeDegrees.size());
        int incrementDegree = nodeDegrees[incrementID];

        //std::cout << "Node " << incrementID << " Has Degree " <<
        //    incrementDegree << " (the Max degree is " << maxDegree << "\n";
        while (incrementDegree >= maxDegree - 2 ||
            degreeDistribution[incrementDegree + 1] == 0){
            //std::cout << "HELLOOO?\n";
            incrementID = int(RNG(g1) * nodeDegrees.size());
            incrementDegree = nodeDegrees[incrementID];
        }
        std::cout << "Increment degree is " << incrementDegree << "????\n";
        ++nodeDegrees[incrementID];
        //std::cout << "Node " << incrementID << " Incremented\n";
    }// odd total degree case
    //std::cout << "Node Degrees (Post Increment etc)\n";
    //PrintVectorInt(nodeDegrees);
    // now randomly wire those suckers up
    // (Note: can fail so loop over)

    //std::cout << "Odd Degrees Dealt with, wiring\n";
    bool success = false;
    int attempt = 0;
    while (1) {
        if (attempt > 1000){
            std::cout << "Bad Degree Distribution: Failure to Initialize\n";
            initFailure = true;
            std::cout << "From Netwokr Initialization: InitFailure: " <<
                initFailure << ", FFS\n";
            return G;
        }
        bool failure = false;
        std::vector<int> stubs;
        stubs.reserve(N);
        std::vector<int> remainingIDs;
        remainingIDs.reserve(N);

        // fill the stubs using the node degrees
        //std::cout << "About to Fill Stubs\n";
        for (int i=0; i < N; i++){
            stubs.emplace_back(nodeDegrees[i]);
            remainingIDs.emplace_back(i);
        }// end degree assignemnt

        //std::cout << "Vectors Initialized\n";
        //std::cout << "Stubs:";
        //PrintVectorInt(stubs);
        //std::cout << "Remaining IDs:";
        //PrintVectorInt(remainingIDs);
        
        //std::cout << "Stubs filled\n";
        // loop through each stub connecting them fully
        int size = remainingIDs.size();
        while (true) {
            

            if (size == 0){
                failure = false;
                //std::cout << "Success\n";
                break;
            }
            else if (size == 1){
                failure = true;
                //std::cout << "A sole node remaining\n";
                //PrintVectorInt(stubs);
                //PrintVectorInt(remainingIDs);
                break;
            }


            int nodeID = remainingIDs[size - 1]; 
            //std::cout << "nodeID = " << nodeID << "\n";
            //PrintVectorInt(stubs);
            //PrintVectorInt(remainingIDs);
            //remainingIDs.erase(nodeID);
            // The size is reduced because a node cannot be connected to itself
            // (itself is not included in potential neighbours = stubs[:self])
                    // randomly assign neighbours
            int neighboursLeft = stubs[nodeID];

            // Keep track of nodes which can be deleted from the list
            std::vector<int> deleteIndices;
            deleteIndices.reserve(neighboursLeft + 1);
            // the parent node will be deleted since it will be full
            deleteIndices.emplace_back(remainingIDs.size() -1);
            size -= 1;

            // if the node needs more neighbours than there are nodes left
            // it's all fucked (shouldnt happen)
            if (neighboursLeft > size && size > 0){
                failure = true;
                //std::cout << "Needs more neighbours than nodes left\n";
                //PrintVectorInt(stubs);
                //PrintVectorInt(remainingIDs);
                break;
            }
            while (neighboursLeft > 0){
                //std::cout << "neighbours left = " << neighboursLeft << "\n";
                int neighbourIndex = (int)(RNG(g1)*size);
                //std::cout << "neighbourIndex = " << neighbourIndex << "\n";
                int neighbourID = remainingIDs[neighbourIndex];
                //std::cout << "neighbourID = " << neighbourID << "\n";
                
                // if they arent connected connect them
                if (edge(nodeID, neighbourID, G).second == false &&
                    stubs[neighbourID] > 0){
                    //std::cout << "Viable Neighbour\n";
                    add_edge(nodeID, neighbourID, G);
                    //std::cout << "Nodes Connected\n";
                    // update number of stubs for both
                    stubs[nodeID] -= 1;
                    //std::cout << "Parent Stub Decremented\n";
                    stubs[neighbourID] -= 1;
                    //std::cout << "neighbour Stub Decremented\n";
                    // if that stub is now full remove it from queue 
                    if (stubs[neighbourID] == 0){
                        deleteIndices.emplace_back(neighbourIndex);
                        //std::cout << "Neighbour To be deleted: Index: " <<
                        //    neighbourIndex << " ID: " << neighbourID << "\n";
                    }// end full stub case
                    
                    neighboursLeft -= 1;
                }// end successful neighbour attachment
            }// end assigning neighbours to ith node

            // Now remove the stubs from the list which are fully satisfied
            
            // NOTE: have to delete nodes in descending order for ID preserval
            // when using Indices.
            if (deleteIndices.size() > 1){
                std::sort(deleteIndices.begin(), deleteIndices.end(),
                    std::greater<int>());
            }
            for (auto x: deleteIndices){
                //std::cout << "Stub " << remainingIDs[x] << " is full\n";
                remainingIDs.erase(remainingIDs.begin() + x);
                //std::cout << "Erased.\n";
                // decrement size (removal of nodes)
                size -= 1;
            }// end full stub removal
            size += 1;
            
        }// end 1 rewiring attempt
        
        // Failure condition (unsatisfied stubs)
        if (failure){
            // remove existing edges (cant be satisfied)
            //if (attempt%60000 ==0) std::cout << attempt << " attempts\n";
            for (int i = 0; i < N; i++){
                clear_vertex(i, G);
            }// end edge removal
            attempt++;
        }// end unsuccessful case
        else {
            if (initFailure) std::cout << "why are we here?\n";
            initFailure = false;
            success = true;
            break;
        }
    }// end global rewiring (Must be successful)
    //std::cout << "Initial Network Generated\n";
    //std::cout << "Number OF edges on G: " << num_edges(G) << "\n";

    return G;

}//end Network initiation


// Given a pair of edges and the graph they're on, choose whether or not to rewire them (in place)
void conditionalRewiring(Graph &G,
    std::vector<std::array<int, 2>> &edgeVector, std::vector<int> degrees,
    std::vector<std::vector<double>> &correlationMatrix,
    int a, int b, int c, int d, int edge0, int edge1, int &numSuccess,
    std::vector<int> &airballs){
    
    int kMax = correlationMatrix.size() - 1;
    //std::cout << "Edges to rewire are: " << a << " <-> " << b << " and " <<
    //    c << " <-> " << d << "\n";
    // Degrees of each of the vertices
    double ka = degrees[a];
    double kb = degrees[b];
    double kc = degrees[c];
    double kd = degrees[d];
    //std::cout << "With Degrees: " << ka << " <-> " << kb << " and " <<
    //    kc << " <-> " << kd << "\n";


    // The "goodness" evaluation of the edge pairings (original E0 and alternative pairing E1) via Pkk
    double E0 = correlationMatrix[ka][kb] + correlationMatrix[kc][kd];
    double E1 = correlationMatrix[ka][kc] + correlationMatrix[kb][kd];
    //std::cout << "E0: " << E0 << ", E1: " << E1 << ".\n";

    // Whether or not to swap the edges
    bool swap = false;

    // obvious case: one of the links does not exist
    if (E0 == 0){
        swap = true;
        // if this is actual progress (swapping has a net positive impact) count it as a success
        if (E1 > 0){
            //std::cout << "Actual Progress\n";
            numSuccess++;
        }
    }// end obvious case

    // both edges exist
    else {
        // calculate probability of rewiring
        double p = E1/E0;
        // clearly better choice
        if (p >= 1){
            swap = true;
            numSuccess++;
            //std::cout << "P > 1\n";
        }// end strictly favourable case
        // annealing case: less probable but non zero for 'wrong' choice
        else {
            double r = RNG(g1);
            if (r < p) {
                swap = true;
                //std::cout << "P > r\n";
                numSuccess++;
            }
            // Failure
            //else std::cout << "Failed Rewire P < r\n";

        }// end probabilistic case (Simulated Annealing type case)

    }// end non-trivial

    // swap the edges
    if (swap){
        // edge list swap
        edgeVector[edge0][1] = c; 
        edgeVector[edge1][0] = b; 

        // graph edge swap
        remove_edge(a, b, G);
        remove_edge(c, d, G);
        add_edge(a, c, G);
        add_edge(b, d, G);
        
        //std::cout << "Pre-Changes: Airballs vector:\n";
        //PrintVectorInt(airballs);

        // Updating the Airball List

        // check if a->b was an airball before
        bool edge0Airball = std::binary_search(airballs.begin(),
            airballs.end(), edge0);

        // if it was an airball, remove it
        if ((correlationMatrix[ka][kc] > 0) && edge0Airball){
            //std::cout << "Edge0: " << edge0 << " No longer an Airball\n";
            removeAirball(airballs, edge0);
        }// end Was Previously airball Case

        // if the new edge is an airball, add it
        if ((correlationMatrix[ka][kc] == 0) && (not edge0Airball)){
            //std::cout << "Edge0: " << edge0 << " Now an Airball\n";
            addAirball(airballs, edge0);
        }

        // repeat for c->d

        bool edge1Airball = std::binary_search(airballs.begin(),
            airballs.end(), edge1);

        // remove old airball
        if ((correlationMatrix[kb][kd] > 0) && edge1Airball){
            //std::cout << "Edge1: " << edge1 << " No longer an Airball\n";
            removeAirball(airballs, edge1);
        }// end Was Previously airball Case

        // add new
        if ((correlationMatrix[kb][kd] == 0) && (not edge1Airball)){
            //std::cout << "Edge1: " << edge1 << " Now an Airball\n";
            addAirball(airballs, edge1);
        }


    }// end edge swap

}// End conditional rewiring


// Reconfigure a graph using Bertotti rewiring scheme (Newmann rewiring really + modifications)
// Basically take a graph and some target correlations and rewire it until it looks like the target.
Graph reconfigureGraph(Graph &G, std::vector<std::vector<double>>
    correlationMatrix, std::vector<int> nodeDegrees, int rewiresPerNode,
    bool &fullyConnected){
    
    using namespace boost;

    int size = num_vertices(G);

    int numEdges = num_edges(G);
    //std::cout << numEdges << " EDGES, " << size << " Vertices\n";
    if (numEdges == 0) std::cout << "NO EDGES, " << size << " Vertices\n";

    // Generate array of edges
    std::vector<std::array<int, 2>> edgeVector;
    edgeVector.reserve(numEdges + 1);

    // Keep track of the edges which should not exist according to target Pkk
    // These are called airballs because they whiff entirely
    std::vector<int> airballs;
    airballs.reserve(numEdges + 1);
 
    // setup edge index finders with boost:
    typedef property_map<Graph, vertex_index_t>::type IndexMap;
    IndexMap index = get(vertex_index, G);
   
    // iterate through edges and store in edge array
    graph_traits<Graph>::edge_iterator ei, ei_end;

    // index of edge being put into array
    //std::cout << "Begin recording Edges\n";
    int i = 0;
    // dummy node indices
    int out = 0, in = 0;
    // degrees of dummy indices
    int inDegree = 0, outDegree = 0;
    for (tie(ei, ei_end) = edges(G); ei != ei_end; ++ei){
        // nodes on either side of edge
        out = index[source(*ei, G)];
        in = index[target(*ei, G)];
        if (out == in){
            std::cout << "Received a self connected Vertex: " << in << "\n";
        }
        std::array<int, 2> dummyEdge = {out, in};
        edgeVector.emplace_back(dummyEdge);

        inDegree = nodeDegrees[in];
        outDegree = nodeDegrees[out];
        // check if the edge is an airball
        if (correlationMatrix[inDegree][outDegree] == 0) {
            //std::cout << "P(" <<  inDegree << "," << outDegree << ") = " <<
            //    correlationMatrix[inDegree][outDegree] << "\n";
            airballs.emplace_back(i);
        }
        // update node degrees
        // record edges
        i++;
    }// end filling edge array
    //std::cout << "Done recording Edges\n";
    //PrintVectorInt(airballs);

    // ##### Rewiring Phase  ######

    // nodes to be rewired 
    int a = 0, b = 0, c = 0, d = 0;
    // edges
    int edge0 = 0, edge1 = 0;

    // defining stopping point
    int numSuccess = 0;
    int endCondition = size * rewiresPerNode;

    // run it
    int numTries = 0;

    //std::cout << "Begin Rewiring\n";

    int airballIndex = 0;

    while (1){

        // if there are still airball edges, try to satisfy them
        if (airballs.size() > 0){
            //std::cout << "Selecting an airball\n";
            airballIndex = (int)(RNG(g1)*airballs.size());
            edge0 = airballs[airballIndex];
        }
        // if there aren't any airballs left (unlikely) just pick a random edge
        else {
            edge0 = (int)(RNG(g1)*numEdges);
        }
        a = edgeVector[edge0][0];
        b = edgeVector[edge0][1];

        
        edge1 = (int)(RNG(g1)*numEdges);
        c = edgeVector[edge1][0];
        d = edgeVector[edge1][1];

        if (a == b || c == d){
            std::cout << "\n\n\n\n NODE CONNECTED TO ITSELF \n\n\n\n\n";
        }
        if (a == c || b == d) {
            std::cout << "\n\n Swapping to itself \n\n";
        }
        
        // Invalid moves: connect a node to itself or an existing neighbour
        
        int weStuck = 0;
        while (c == a || b == d || edge(a, c, G).second ||
            edge(b, d, G).second){

            edge1 = (int)(RNG(g1)*numEdges);
            c = edgeVector[edge1][0];
            d = edgeVector[edge1][1];
            
            weStuck++;
            if (weStuck == 101){
                //std::cout << "@@@@@@@@@@@@@@@@@@@@WeStuck " << weStuck << "\n";
                //std::cout << "NumEdges = " << numEdges << "\n";
                //giveUp = true;
                break;
            }// Stuck somehow
            
        }// end invalid pairs case
        
        
        //if (numTries%20000 == 0) std::cout << numTries << " Tries\n";
        //if (giveUp) break;

        // Rewire if favourable
        //std::cout << "Enter Conditional Rewiring\n";
        conditionalRewiring(G, edgeVector, nodeDegrees, correlationMatrix,
            a, b, c, d, edge0, edge1, numSuccess, airballs);
        //std::cout << "Done Conditional Rewiring\n";
        numTries++;

        if (numSuccess > endCondition){
            std::cout << "Number of Airballs at Num Success: " <<
                airballs.size() << "\n";
            if (airballs.size() < 100){
                break;
            }
        }
    }// end rewiring
    //std::cout << "End Rewiring, " << numTries << " Total Attempts for " <<
    //    numSuccess << " completed moves.\n";
    //" Percent of moves successful\n";
    //std::cout << "Checking Connectedness\n"; 
    std::cout << "From Bertotti Code:\n";
    for (auto x: airballs){
        std::cout << "Node " << edgeVector[x][0] << " to " <<
            edgeVector[x][1] << " is an airball\n";
    }
    fullyConnected = isGraphConnected(G); 

    numEdges = num_edges(G);
    //std::cout << "Rewired Graph has: " << numEdges << " edges.\n";


    return G;

}// end rewiring algorithm


// generate a network as described in Bertotti paper
std::vector<Node> networkFromCorrelations(
    std::vector<double> degreeDistribution,
    std::vector<std::vector<double>> correlationMatrix, int N,
    int rewiresPerNode, double avgDeg, bool &initFailure){

    Graph G;

    std::vector<int> nodeDegrees;
    nodeDegrees.reserve(N + 1);
    
    // Have to keep track of whether or not this is satisfiable
    bool satisfiable = false;
    std::vector<Node> Network;
    Network.emplace_back(1);

    std::cout << "Generating Initial Graph";
    //std::cout << "initFailure = " << initFailure << "\n";
    G = randomInitialNetwork(degreeDistribution, N, nodeDegrees, avgDeg,
        initFailure);

    //std::cout << "Graph Initialized, initFailure = " << initFailure << "\n";
    if (initFailure){
        std::cout << "\n\n\n\nFailure To Initialize\n\nWhy is this here?\n\n";
        return Network;
    }
    std::cout << "Initial Graph Completed \n";

    bool fullyConnected = true;

    std::cout << "Reconfiguring Graph\n";
    reconfigureGraph(G, correlationMatrix, nodeDegrees,
        rewiresPerNode, fullyConnected);

    std::cout << "Made it out.\n";
    // in the case of a partially connected graph brute force it with original
    // algorithm
    //rewiresPerNode = 1;

    std::string successType = "Full";
    std::vector<Node> BrokenNetwork;
    std::vector<Node> giantComponent;


    int attempts = 0;
    //if (not fullyConnected) std::cout << "Initial Graph not connected\n";
    bool goodEnough = fullyConnected;
    while (not goodEnough){
        std::cout << "Not Connected, Trying again pt. " << attempts << "\n";
        //G.clear();
        //std::cout << "Initializing\n";
        //std::cout << "Generating Graph\n";
        //std::cout << "Initializing\n";
        nodeDegrees.clear();
        G = randomInitialNetwork(degreeDistribution, N, nodeDegrees, avgDeg,
            initFailure);
        if (initFailure){
            std::cout << "already satisfied but somehow broken..\n";
            std::cout << "\n\n\n\n\n\n\nFailure To Initialize\n\n\n\n\n\n\n";
            std::cout << "Init Failure: " << initFailure << "\n";
            return Network;
        }
        //std::cout << "Initialized\n";
        //std::cout << "Generated Graph\n";
        //std::cout << "Initialized, Rewiring\n";
        //std::cout << "Reconfiguring Graph\n";
        reconfigureGraph(G, correlationMatrix, nodeDegrees,
            rewiresPerNode, fullyConnected);
        goodEnough = fullyConnected;
        //std::cout << "Reconfigured Graph\n";
        //std::cout << attempts << "\n";
        //if (attempts % 101 == 0) std::cout << attempts << " Attempts dawg\n";
        // If not fully connected, see if the giant component is good enough
        if (not goodEnough){
            std::cout << "Incomplete network, Checking Giant Component\n";
            BrokenNetwork = networkFromGraph(G, N);
            std::cout << "Calculated network from graph\n";
            giantComponent = CheckConnected(BrokenNetwork);
            std::cout << "Calculated Giant Component Successfully\n";

            // if the giant component is big enough, check if it satisfies
            // average degree
            if (giantComponent.size() >= 0.95*N){

                std::vector<std::vector<double>> giantPkk = 
                    correlationsFromNetwork(giantComponent, avgDeg, N);
                std::vector<double> giantPk = calculateDegreeDistribution(
                    giantPkk, avgDeg);

                double giantAvgDeg = 0.0;
                for (int k = 1; k < giantPk.size(); k++){
                    giantAvgDeg += k*giantPk[k];
                }// end average degree calculation
                std::cout << "Giant Avg Deg: " << giantAvgDeg << "\n";
                if (abs(avgDeg - giantAvgDeg)/avgDeg < 0.05){
                    goodEnough = true;
                    successType = "Giant";
                }
 
            }// End checking avg Degree

        }// Not fully connected, check if good enough giant component
    
        attempts++;
     }

    //std::cout << "Successfully Generated Network!!\n";
    if (successType == "Full"){
        Network = networkFromGraph(G, N);
        for (auto x: Network){
            if (x.k >= 32){
                std::cout << "DISASTER, Degree >= 32\n";
            }
        }
    }
    else {
        Network = giantComponent;
    }
    std::cout << "Size of Network: " << Network.size() << "\n";
    return Network;
    
}// end network generation

std::vector<Node> generateNetwork(
    std::vector<std::vector<double>> &logCorrelationMatrix,
    std::vector<std::vector<double>> &correlationMatrix,
    std::vector<double> &degreeDistribution,
    std::vector<int> kToLConversion, double dP,
    int numChanges, double avgDeg, int N, int rewiresPerNode){

    std::vector<std::vector<double>> tempLogCorrelationMatrix;
    std::vector<std::vector<double>> tempCorrelationMatrix;

    std::vector<Node> TempNetwork;
    int extension = 0;
    std::vector<double> tempDegreeDistribution;
    bool initFailure = true;
    int failureCount = 0;
    while (initFailure){
        initFailure = false;
        // correlation matrix to modify
        if (failureCount > 0){
            std::cout << "Trying again, failure count: " << failureCount <<
            "\n";
        }
        tempLogCorrelationMatrix = logCorrelationMatrix;

        //modify it 
        for (int j = 0; j < numChanges; j++){
            changeLogCorrelationMatrix(tempLogCorrelationMatrix,
                kToLConversion, dP);
        }// End degree correlation modification
        tempCorrelationMatrix = logToLinCorrelationMatrix(
            tempLogCorrelationMatrix, kToLConversion);
        tempDegreeDistribution =
            calculateDegreeDistribution(tempCorrelationMatrix, avgDeg);
        // generate the network
        //std::cout << "After Correlation Changes\n";
        std::cout << "About to generate a new Temp NEtwork\n";
        TempNetwork = networkFromCorrelations(
            tempDegreeDistribution, tempCorrelationMatrix, N, rewiresPerNode,
            avgDeg, initFailure);
        std::cout << " Finished building it\n";
        /*
        if (TempNetwork.size() < N){
            std::cout << "Tiny Baby Network, InitFailure: " << initFailure <<
                "\n";
            std::cout << "FAILURE\n";
        }
        */
        if (initFailure){
            std::cout << "\n\n\n\n\n\n\nFailure To Initialize\n\n\n\n\n\n\n";
        }
        //std::cout << "InitFailure: " << initFailure << "\n";
        failureCount++;

    }// end network generation (asserting solvable Pkk)
    degreeDistribution = tempDegreeDistribution;
    correlationMatrix = tempCorrelationMatrix;
    logCorrelationMatrix = tempLogCorrelationMatrix;
    std::cout << "Finished generating successfully\n";
    return TempNetwork;

}// end network generation


