//
// Remake of Bertotti.h essentially
// That Code is pure fucking dogshit
// Goals: ROBUST CODE: GOOD HANDLING OF FAILURE CONDITIONS
//


// define the graph stuff
typedef boost::adjacency_list <boost::vecS, boost::vecS,
boost::undirectedS> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor MyVertex;

// function to check if a sampled degree distribution can ever be satisfied
bool checkDegreeSequence(std::vector<int> stubs) {
    /*
     * loop through the degree sequence attaching nodes (no self or double
     connections)
     * target node is always the 0th index, attaching first to last should be
     * fastest since
     * nodes are sorted by degree
     * failure condition is when there is a node with degree larger than the
     * remaining number of nodes
     * in the degree sequence
     * code picks this up by having neighbour index equal to target index (0)
    */
    // index of the neighbour being attached to the target node
    int neighbourIndex = 0;
    //std::cout << "Checking The Following Candidate Degree Sequence:\n";
    //PrintVectorInt(stubs);

    // just connect the nodes non-randomly to see if it can be done
    while (stubs.size() > 0){
        // take the first element in the degree sequence and pair it to
        // neighbours on the end
        neighbourIndex = stubs.size() - 1;

        if (stubs[0] > stubs.size()){
            return false;
        }
        // pair it to the last element until it has enough edges
        while (stubs[0] > 0){

            // decrement stubs of neighbour
            stubs[neighbourIndex] -= 1;
            // decrement index of original node
            stubs[0] -= 1;

            // if the neighbour is out of stubs delete it
            if (stubs[neighbourIndex] == 0){
                stubs.erase(stubs.begin() + neighbourIndex);
                neighbourIndex -= 1;
                //std::cout << "Neighbour emptied, neighbour index: " << 
                //    neighbourIndex << "\n";
            }// end neighbour erase condition

            // if that neighbour still has excess stubs, select the next
            // neighbour candidate in line
            else {
                // check that this index is not the same as the target node
                // (avoid self connection)
                if (neighbourIndex == 0){
                    //std::cout << "Only Neighbour is self, Stubs left: " <<
                    //stubs[0] << "\n";
                    return false;
                }// end case where there are no longer enough neighbours to
                //satisfy the needs of parent node

                // decrement neighbour index (select a new neighbour to connect to.
                neighbourIndex -= 1;
            }// end neighbour re-selection case

        }// All edges of original node satisfied

        // remove original node since it is satisfied
        stubs.erase(stubs.begin());
    }// end loop

    // if made it out of the loop the degree sequence can be satisfied
    //std::cout << "Made it out!\n";
    //PrintVectorInt(stubs);
    return true;
}// end degree sequence check

// function to find bin number for arbitrary degree grouping.
int findBinNumber(int degree, std::vector<int> kToLConversion){
    int binNumber = 0;
    int nextBinLowerDegree = kToLConversion[binNumber + 1];
    int size = kToLConversion.size();
    while (degree >= nextBinLowerDegree){
        binNumber++;
        nextBinLowerDegree = kToLConversion[binNumber + 1]; 
    }
    //std::cout << "k = " << degree << ", Bin Number: " << binNumber << 
    //    ", next bin lower degree: " << nextBinLowerDegree << "\n";
    return binNumber;
}// end bin number finding

// Code to sample a degree distribution ENSURING that it satisfies constraints and is
// satisfiable
std::vector<int> sampleDegreeDistribution(
                                        std::vector<double> degreeDistribution,
                                        int N, double avgDeg, bool &countOut){


    std::vector<int> degreeSequence;
    degreeSequence.reserve(N);

    int maxDegree = degreeDistribution.size();
    
    // Check to see if a giant component can exist
    double firstMoment = 0.0;
    double secondMoment = 0.0;
    
    for (int k = 0; k < maxDegree; k++){
        firstMoment += k*degreeDistribution[k];
        secondMoment += k*k*degreeDistribution[k];
    }

    double giantTest = secondMoment - 2*firstMoment;
    if (giantTest <= 0){
        std::cout << "No Giant Component Possible: BREAKING\n";
        countOut = true;
        return degreeSequence;
    }

    // this vector holds the number of times a node of each degree has been sampled
    std::vector<int> degreeCounts (maxDegree, 0);

    // the sum of the degrees in the network is used to calculate average degree (also must be even)
    int degreeSum = 0;

    // placeholders for filling degrees
    int deg = 0;

    // generate the sampler
    std::discrete_distribution<int> sampleDist (degreeDistribution.begin(),
                                                degreeDistribution.end());

    // generate stubs
    int numStubs = 0;
    int stubsSum = 0;
    /*
    for (int i = 0; i < maxDegree; i++){
        // fill degrees in ascending order
        deg = i; // Descending order case: maxDegree - i;
        // number of stubs
        numStubs = degreeCounts[deg];
        // add the correct number of stubs for the given degree
        while (numStubs > 0){
            degreeSequence.emplace_back(deg);
            stubsSum += 1;
            numStubs -= 1;
        }// end stub filling

    }// end stub generation
    PrintVectorInt(degreeSequence);
    */
    // need to be within 5% of the average degree (this is the tolerance threshold)
    bool satisfied = false;
    int numTries = 0;
    double sampledAvgDeg = 0.0;
    while (not satisfied){
        // Reset number of nodes
        
        std::fill(degreeCounts.begin(), degreeCounts.end(), 0);
        degreeSum = 0;
        // sample the degree distribution N times
        for (int i=0; i < N; i++){
            deg = sampleDist(g1);
            //std::cout << "Adding Node of degree " << degree << "\n";
            ++degreeCounts[deg];
            degreeSum += deg;
        }// end degree sampling
        //PrintVectorInt(degreeCounts);

        // Assert that the average degree is within 5% of desired
        sampledAvgDeg = double(degreeSum) / double(N);
        double percentError = abs(avgDeg - sampledAvgDeg) / avgDeg;

        
        degreeSequence.clear();
        // if the average degree is within tolerance and the degree sum is
        // even, check for satisfiability
        //std::cout << "Sampled Average Degree: " << sampledAvgDeg << "\n";
        if (percentError < 0.05 && degreeSum % 2 == 0){
            // generate a degree sequence from the sampled degree counts
            for (int k = 0; k < maxDegree; k++){
                // number of stubs
                numStubs = degreeCounts[k];
                // add the correct number of stubs for the given degree
                while (numStubs > 0){
                    degreeSequence.emplace_back(k);
                    numStubs -= 1;
                }// end stub filling (degree sequence now has the correct
                //number of nodes of that degree)

            }// end stub generation

            // check satisfiability with function above
            satisfied = checkDegreeSequence(degreeSequence);

        }// end satisfiability check
        
        // add a giving up clause
        numTries += 1;
        if (numTries > 1000000000){
            countOut = true;
            std::cout << "Degree Sampling Countout???\n";
            return degreeSequence;
        }//end clause

    }// end Avg K within 5% insurance + degree assignment

    // made it out: the degree sequence satisfies all constraints and can be
    // built into a network.
    //std::cout << "From Degree Sequence Generator: <k> = " << sampledAvgDeg <<
    //            ".\n";
    return degreeSequence;

}// end degree sampling function


// Generate a random initial network from a degree sequence
Graph randomInitialNetwork(std::vector<int> degreeSequence, int N,
    bool &countOut){
    using namespace boost;

    //std::cout << "Just got in\n";
    // graph to keep track of connections and avoid double connections. (size N)
    Graph G(N);

    // exit conditions: fully connected random network
    bool satisfied = false;
    // exit condition: tried 1000 times and still couldn't wire it up
    int numAttempts = 0;

    // stubs is remaining number of connections needed for a given node
    std::vector<int> stubs = degreeSequence;
    stubs.reserve(N);
    // remaining IDs is the IDs of the nodes which still have non-zero needed connections
    std::vector<int> remainingIDs;
    remainingIDs.reserve(N + 1);

    // the node to be connected is always the highest remaining degree (so the
    // last in ascending).
    //std::cout << "Finding First Node ID\n";
    int nodeID;
    // number of neighbours to find for this node
    int neighboursLeft;

    // the IDs of the nodes which will end up satisfied after wiring Up this node
    std::vector<int> deleteIndices;
    deleteIndices.reserve(neighboursLeft + 1);

    //std::cout << "Starting the Wiring!\n";    
    // try to wire it up randomly until it friggin works.
    int numTries = 0;
    while (1){
        // reset the degree sequence
        stubs = degreeSequence;
        // reset the remaining IDs
        remainingIDs.clear();
        for (int i = 0; i < N; i++){
            remainingIDs.emplace_back(i);
        }// end replacing remaining IDs

        // one attempt at rewiring (tried over and over)
        while (1) {

            // Most common Failure Conditions:
            //PrintVectorInt(remainingIDs);
            // if there are no remaining nodes to be wired up it's a victory.
            if (remainingIDs.size() == 0) {
                satisfied = true;
                break;
            }
            // if there's only one node with remaining stubs then we have
            // failed, go again.
            if (remainingIDs.size() == 1) {
                satisfied = false;
                break;
            }

            // select the last node as the parent
            nodeID = remainingIDs[remainingIDs.size() - 1];
            // find number of remaining connections
            neighboursLeft = stubs[nodeID];
            // reset the indices of nodes to be cleared
            deleteIndices.clear();
            // the parent node will necessarily be deleted (satisfied) (remove
            // it now)
            remainingIDs.erase(remainingIDs.end() - 1);

            // if this node cannot possibly be satisfied (needs more neighbours
            // than nodes left)
            if (neighboursLeft > remainingIDs.size()) {
                satisfied = false;
                break;
            }

            // randomly wire the parent node up to neighbours
            while (neighboursLeft > 0) {
                // find a neighbour out of the remaining unsatisfied nodes
                int neighbourIndex = int(RNG(g1) * remainingIDs.size());
                int neighbourID = remainingIDs[neighbourIndex];

                // if these nodes are not already neighbours
                // and the neighbour has to still be unsatisfied
                // shouldn't have to check if neighbour stubs > 0 since
                // To be in remaining IDs going in it must be 0, and it could
                // only connect to the parent node
                if (edge(nodeID, neighbourID, G).second == false) {
                    //wire them up
                    add_edge(nodeID, neighbourID, G);
                    // decrement remaining edges to satisfy parent and target
                    stubs[nodeID] -= 1;
                    stubs[neighbourID] -= 1;
                    // if that neighbour is now satisfied, add it to delete Q
                    if (stubs[neighbourID] == 0) {
                        deleteIndices.emplace_back(neighbourIndex);
                    }
                    // decrement remaining stubs
                    neighboursLeft -= 1;
                }// end viable connection case

            }// end filling neighbours of parent node

            // Now delete all of the nodes in the deleteIndices queue
            // first sort descending ID due to order preserval
            if (deleteIndices.size() > 1){
                std::sort(deleteIndices.begin(),  deleteIndices.end(),
                          std::greater<int>());
            }// end ID sorting
            // Now remove all of those node from the remaining IDs list
            for (auto x: deleteIndices){
                remainingIDs.erase(remainingIDs.begin() + x);
            }

        }// End of 1 rewiring Attempt

        // if the graph is not satisfied, clear it and go again
        if (not satisfied){
            for (int i = 0; i < N; i++){
                clear_vertex(i, G);
            }// end edge removal
            //PrintVectorInt(remainingIDs);
            numAttempts++;
            if (numAttempts > 10000000){
                std::cout << "Num Attempts: " << numAttempts << " :(\n";
                std::cout << "Down for the count, Breaking\n";
                countOut = true;
                return G;
            }
        }// end unsuccessful case
        else {
            //std::cout << "Random Initial Graph Created on attempt: " <<
            //    numAttempts << ".\n";
            return G;
        }// end successful case

    }// end Re-Trying loop.

    //std::cout << "I guess I output nothing?\n";
    return G;

}// end random initial network generation

// function that randomly selects 2 edges in the graph (a-b + c-d) and swaps
// them if (a-c + b-d) is favourable
// Will always find viable edges to swap (IE no self connections etc) so
// failures are due to the degree correlations.
void swapEdges(Graph &G,
               std::vector<std::array<int, 2>> &edgeVector,
               std::vector<int> degrees,
               std::vector<std::vector<double>> logMatrix,
               std::vector<int> kToLConversion,
               int &numSuccess){

    int a, b, c, d;
    int edge0, edge1;
    int numEdges = edgeVector.size();

    // just pick a random edge
    edge0 = int(RNG(g1) * numEdges);

    a = edgeVector[edge0][0];
    b = edgeVector[edge0][1];

    // other edge is randomly selected
    edge1 = int(RNG(g1) * numEdges);
    c = edgeVector[edge1][0];
    d = edgeVector[edge1][1];

    // illegal moves: a and c (or b and d) same node (self wiring case)
    // or a already connected to c (or b to d) (double connected case).
    int retries = 0;
    while (c == a || b == d || edge(a, c, G).second ||
           edge(b, d, G).second) {

        edge1 = (int) (RNG(g1) * numEdges);
        c = edgeVector[edge1][0];
        d = edgeVector[edge1][1];

        retries++;
        // in the case where a viable pair cannot be found (should never happen?)
        if (retries > 2 * numEdges){
            edge0 = (int) (RNG(g1) * numEdges);
            a = edgeVector[edge0][0];
            b = edgeVector[edge0][1];
            retries = 0;
            //std::cout << "Had to find a new Edge-0\n";
        }// end cant find viable pair case

    }// end illegal moves avoiding

    // Now that we have viable candidate edges we can check whether or not to
    // rewire based on correlations
    // need the degrees
    int kA = degrees[a];
    int kB = degrees[b];
    int kC = degrees[c];
    int kD = degrees[d];

    int nA = findBinNumber(kA, kToLConversion);
    int nB = findBinNumber(kB, kToLConversion);
    int nC = findBinNumber(kC, kToLConversion);
    int nD = findBinNumber(kD, kToLConversion);

    // The "goodness" evaluation of the edge pairings (original E0 and
    // alternative pairing E1) via Pkk
    double E0;
    double E1;
    E0 = logMatrix[nA][nB] * logMatrix[nC][nD];
    E1 = logMatrix[nA][nC] * logMatrix[nB][nD];

    // Whether or not to swap the edges
    bool swap = false;

    // obvious case: one of the links does not exist
    if (E0 == 0){
        swap = true;
        // if this is actual progress (swapping has a net positive impact) 
        // count it as a success
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
    if (swap) {
        // edge list swap
        edgeVector[edge0][1] = c;
        edgeVector[edge1][0] = b;

        // graph edge swap
        remove_edge(a, b, G);
        remove_edge(c, d, G);
        add_edge(a, c, G);
        add_edge(b, d, G);

    }// end successful swap case (updating edge lists/graph)

}// end swap Edges


// function that takes a randomly initialized graph and rewires it (according
// to the Newman rewiring scheme) to match a target correlation matrix
Graph imposeCorrelations(Graph &G,
                         std::vector<std::vector<double>> logMatrix,
                         std::vector<int> nodeDegrees, 
                         std::vector<int> kToLConversion,
                         int rewiresPerNode, bool &countOut){

    using namespace boost;

    int size = num_vertices(G);

    int numEdges = num_edges(G);
    //std::cout << numEdges << " EDGES, " << size << " Vertices\n";
    if (numEdges == 0) std::cout << "NO EDGES, " << size << " Vertices\n";

    // Generate array of edges
    std::vector<std::array<int, 2>> edgeVector;
    edgeVector.reserve(numEdges + 1);

    // setup edge index finders with boost:
    typedef property_map<Graph, vertex_index_t>::type IndexMap;
    IndexMap index = get(vertex_index, G);

    // iterate through edges and store in edge array
    graph_traits<Graph>::edge_iterator ei, ei_end;

    // index of edge being put into array
    int i = 0;
    // dummy node indices
    int out = 0, in = 0;

    for (tie(ei, ei_end) = edges(G); ei != ei_end; ++ei){
        // nodes on either side of edge
        out = index[source(*ei, G)];
        in = index[target(*ei, G)];
        if (out == in){
            std::cout << "Received a self connected Vertex: " << in << "\n";
        }
        std::array<int, 2> dummyEdge = {out, in};
        edgeVector.emplace_back(dummyEdge);

        // record edges
        i++;
    }// end filling edge array
    //std::cout << "Done recording Edges\n";

    // Now we rewire the network by swapping edges between pairs of nodes
    int numSuccess = 0;
    int targetRewires = size * rewiresPerNode;
    int numTries = 0;
    while (numSuccess < targetRewires){
        /*
        if (numTries % 1000000 == 0){
            std::cout << numTries << " Tries and only " << numSuccess <<
                " Successes :(\n";
        }
        */
        if (numTries >= 10000000){
            std::cout << "10 Million Tries at rewiring, BREAKING\n";
            countOut = true;
            break;
        }
        swapEdges(G, edgeVector, nodeDegrees, logMatrix, kToLConversion,
                  numSuccess);
        numTries++;
    }// end rewiring phase

    return G;

}// end impose correlations

// a boolean function to check whether a graph is fully connected
// if it isnt fully connected check whether or not the giant component is
// satisfactory ( size +  avg deg within 5% of desired)
bool checkGraph(Graph &G, int N, double avgDeg, int kMax){
    using namespace boost;
    bool fullyConnected = isGraphConnected(G);
    if (fullyConnected){
        return true;
    }
    else{
        std::vector<Node> BrokenNetwork = networkFromGraph(G, N);
        //std::cout << "Calculated network from graph\n";
        std::vector<Node> giantComponent = CheckConnected(BrokenNetwork);
        //std::cout << "Calculated Giant Component Successfully\n";

        // if the giant component is big enough, check if it satisfies
        // average degree
        if (giantComponent.size() >= 0.95*N){

            // calculate correlations and degree distribution (slow way to do it)
            std::vector<std::vector<double>> giantPkk =
                    correlationsFromNetwork(giantComponent, avgDeg, N, kMax);
            std::vector<double> giantPk = calculateDegreeDistribution(
                    giantPkk, avgDeg);
            // calculate the average degree of the giant component
            double giantAvgDeg = 0.0;
            for (int k = 1; k < giantPk.size(); k++){
                giantAvgDeg += k*giantPk[k];
            }// end average degree calculation
            //std::cout << "Giant Avg Deg: " << giantAvgDeg << "\n";

            // check whether this is within tolerance (success exit condition)
            if (abs(avgDeg - giantAvgDeg)/avgDeg < 0.05){
                return true;
            }// checking avg degree tolerance

        }// End checking avg Degree (since the giant component is big enough

    }// end not fully connected case (checking giant component

    // if we made it here the graph is no good
    return false;

}// end checking if graph is satisfactory

// generate network from distributions. Must be successful?
std::vector<Node> generateNetwork(std::vector<double> degreeDistribution,
                    std::vector<std::vector<double>> logMatrix,
                    std::vector<int> kToLConversion,
                    int N, double avgDeg, int rewiresPerNode,
                    bool &giveUp){

    int kMax = kToLConversion[kToLConversion.size()-1];
    // need to generate a graph which satisfies degree distribution/correlation
    Graph G;
    // and is either fully connected or has a satisfactory giant component
    bool satisfiedGraph = false;
    int numGenerationAttempts = 0;
    int maxAttempts = 1000;
    bool countOut = false;
    while (not satisfiedGraph){
        numGenerationAttempts += 1;
        if (numGenerationAttempts > maxAttempts){
            std::cout << "General Failure to Build Full Network, BREAKING\n";
            giveUp = true;
            break;
        }
        //std::cout << "Entered Generation attempt " << numGenerationAttempts <<
        //    "\n";
        // calculate a satisfactory degree sequence by sampling P(k)
        //std::cout << "Sampling Degree Distribution\n";
        std::vector<int> degreeSequence = sampleDegreeDistribution(
                degreeDistribution, N, avgDeg, countOut);
        if (countOut){
            giveUp = true;
            break;
        }
        //std::cout << "Generating initial random Network\n";
        G = randomInitialNetwork(degreeSequence, N, countOut);
        if (countOut){
            giveUp = true;
            break;
        }
        //std::cout << "Imposing Correlations On network\n";
        //PrintVectorInt(degreeSequence);
        imposeCorrelations(G, logMatrix, degreeSequence,
                            kToLConversion, rewiresPerNode, countOut);
        if (countOut){
            giveUp = true;
            break;
        }
        //std::cout << "Checking Satisfiability of network\n";
        satisfiedGraph = checkGraph(G, N, avgDeg, kMax);
    }
    std::vector<Node> Network = networkFromGraph(G, N);
    return Network;

}// end generate network


//NOW DEFUNCT APPROACH: AIRBALL MANAGEMENT 
// functions for keeping track of airball edges (edges which shouldnt exist
// according to P(k',k)
/*
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

*/
