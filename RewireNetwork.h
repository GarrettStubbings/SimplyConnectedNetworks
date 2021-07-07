// exception to throw when found the disconnected node
struct FoundOld : public std::exception {
    const char * what () const throw () {
        return "Disconnected Node Found";
    }
};



// define the graph stuff
typedef boost::adjacency_list <boost::vecS, boost::vecS, 
    boost::undirectedS> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor MyVertex;


// Searching for the disconnected nodes (ensure connectedness)
class MyVisitor : public boost::default_bfs_visitor 
{
    MyVertex match;
public:

    void SetMatch (MyVertex);


    void discover_vertex(MyVertex v, const Graph& g) {
        // if visited vertex matches the one it was disconnected from, or one
        // of the neighbours of the one it was disconnected from the network is
        // whole
        if (v == match || edge(match, v, g).second){
            throw FoundOld();
            }
    
    } //end discover vertex

}; //end custom bfs visitor



// setting the vertex to find to ensure unity
void MyVisitor::SetMatch (MyVertex vert) {
    match = vert;
}

// OLD AND SLOW AND TERRIBLE
// Rewire a Network nRewire Times
std::vector<Node> RewireNetwork(const std::vector<Node> OriginalNetwork,
    const int nRewire) {

    using namespace boost;
    // size of the network
    int size = OriginalNetwork.size();
    //initialization of network
    std::vector<Node> Network = OriginalNetwork;

    // seeding the random function
    // Should be done already in the main code?
    //auto testSeed = RNG(g1);
    
    // initialise the graph
    Graph G;

    // fill the graph
    for(const auto &v: Network) {

	for(auto x: v.Neighbours) {
	    
	    auto res = add_edge(v.id, Network[x].id, G);
	}
    } 
        
    // rewire the network the requisite number of times
    int num = 0;
    int failure = 0;
    bool pass = false;
    
    //Selection of node to be connected
    //disallowed are self, neighbour, disconnected
    //all the nodes
    std::vector<int> connectCandidates;
    connectCandidates.reserve(size);
    push_back(connectCandidates, irange(0, size));

    while (num < nRewire){
        pass = false;
    
        // Initialise the network to be rewired
        std::vector<Node> RewiredNetwork = Network;

        //selection of the node which will have its connection changed
        int switchee = (int)(RNG(g1)*size);

        //selection of the connection to be removed
        //valid candidates must have > 1 connection (unity of network condition)
        std::vector<int> disconnectCandidates;
        // must be more than 0 disconnect candidates
        int discSize = 0;
        while (discSize == 0){
            disconnectCandidates.clear();
            for(auto &v: Network[switchee].Neighbours){
                disconnectCandidates.emplace_back(v);
            }
            discSize = disconnectCandidates.size();
        }
      
        //Select a Random node from candidates and disconnect it
        int discIndex = (int)(RNG(g1)*disconnectCandidates.size());
        int discID = disconnectCandidates[discIndex]; 
        
        // Remove that connection in the graph
        remove_edge(switchee, discID, G);

      
        int conIndex;
        int conID;
        // Pick a node at random ensureing it's not connected or the
        // recently disconnected node
        bool isConnected = true;
        while (isConnected || conID == switchee){
            conIndex = (int)(RNG(g1)*connectCandidates.size());
            conID = connectCandidates[conIndex]; 
            isConnected = edge(switchee, conID, G).second;
        }
        //Add connection in graph
        add_edge(switchee, conID, G);
        
        // Use a search to find out if the network is whole
        // visitor to search the network
        MyVisitor vis;
        // Node to be found
        vis.SetMatch(switchee);
        
        // check if the network is whole
        // exception handling apparently not fantastic for this application
        // dont know how else to terminate the search.
        try{
            // Breadth First Search through the Graph
            breadth_first_search(G, vertex(discID, G), visitor(vis));
        }   catch(FoundOld& e){
            // Successful Rewiring (Disconnected Node found)
            
            //UPDATE THE NETWORK ITSELF
            // Disconnect the Node
            RewiredNetwork[switchee].RemoveConnection(RewiredNetwork[discID]);
            //Connect the new nodes    
            RewiredNetwork[switchee].ConnectNode(RewiredNetwork[conID]);
            // Update Network
            Network = RewiredNetwork;

            num++;
            pass = true;
        }
        if (pass == false){
            // Unsuccessful Rewiring, Reconnect the old edges on the graph
            remove_edge(switchee, conID, G);
            add_edge(switchee, discID, G);
        }
        
    } // End Rewiring loop
    for (auto &x: Network) x.Set_knn(Network);
    return Network;
}// end RewireNode

// guess what the next few Nrewires should be
void newNRewire(int &NRewire, int success, int N){
    double randomExp = - log( RNG(g1) ) / 2;
    double result = 1 + (int)randomExp;
    if (result > N/5) result = N/5;

    NRewire = result;

}// end new N Rewires calculation

Graph graphFromNetwork(const std::vector<Node> Network){
    using namespace boost;
    Graph G;

    // fill the graph
    for(const auto &v: Network) {

	for(auto x: v.Neighbours) {
	    
            auto res = add_edge(v.id, Network[x].id, G);
	}
    } 
    return G; 
}//end graph from network code

std::vector<Node> networkFromGraph(Graph G, int N){
    
    using namespace boost;
    int size = N;
    //std::cout << "Gello\n";
    std::vector<Node> Network;
    Network.reserve(N+1);
    for (int i = 0; i < size; i++) Network.emplace_back(i);

    typedef property_map<Graph, vertex_index_t>::type IndexMap;
    IndexMap index = get(vertex_index, G);
   
    // iterate through edges and store in edge array
    graph_traits<Graph>::edge_iterator ei, ei_end;

    // index of edge being put into array
    int i = 0;
    int out = 0, in = 0;
    for (tie(ei, ei_end) = edges(G); ei != ei_end; ++ei){
        // nodes on either side of edge
        out = index[source(*ei, G)];
        in = index[target(*ei, G)];
        Network[out].ConnectNode(Network[in]);
    }// end filling edge array


    // set up the nearest neighbour stuff.
    for (auto &x: Network) x.Set_knn(Network);
    return Network;
}

bool isGraphConnected(Graph G){
    // check the size of the giant component
    std::vector<int> component(num_vertices(G));
    int numComponents = connected_components(G, &component[0]);
    if (numComponents == 1) return true;
    else return false;
}

// Rewire a Network nRewire Times
std::vector<Node> RewireNetworkFast(const std::vector<Node> OriginalNetwork,
    const int nRewire) {

    using namespace boost;

    // size of the network
    int size = OriginalNetwork.size();
    //initialization of network
    std::vector<Node> Network = OriginalNetwork;

    Graph OriginalG = graphFromNetwork(OriginalNetwork);
    Graph G = OriginalG;
    // seeding the random function
    // Should be done already in the main code?
    //auto testSeed = RNG(g1);

    // Run until theres a completely connect network
    bool isConnected = false;
    
    while (!isConnected){
        // rewire the network the requisite number of times
        int num = 0;
        std::vector<int> viableNeighbours;
        int discIndex = 0;
        int connectIndex = 0;
        int selected = 0;
        int discID;
        Network = OriginalNetwork;
	G = OriginalG;

        while (num < nRewire){
            // Pick a random Node
            selected = (int)(RNG(g1)*size);
            viableNeighbours = Network[selected].Neighbours;
            
            // Pick a neighbour to disconnect   
            discIndex = (int)(RNG(g1)*viableNeighbours.size());
            discID = viableNeighbours[discIndex];
            
            // ensure that neighbour will remain part of the network
            while (Network[discID].k == 1){
                
                // non-viable nodes erased
                viableNeighbours.erase(viableNeighbours.begin() + discIndex);
                
                // run out of viable neighbours: pick a different node
                if (viableNeighbours.size() == 0){
                    selected = (int)(RNG(g1)*size);
                    viableNeighbours = Network[selected].Neighbours;
                }// end k > 1 check / reselect Node

                discIndex = (int)(RNG(g1)*viableNeighbours.size());
                discID = viableNeighbours[discIndex];
                
            }// end check disconnect Candidate k
                
            // disconnect the node
            Network[selected].RemoveConnection(Network[discID]);
            remove_edge(selected, discID, G);

            // pick a node to connect
            connectIndex = (int)(RNG(g1)*size);
            // Cant be self or existing neighbour
            while (connectIndex == selected ||
                edge(selected, connectIndex, G).second){
                connectIndex = (int)(RNG(g1)*size);
            }// end proper connection assertion

            //connect nodes
            Network[selected].ConnectNode(Network[connectIndex]);
            add_edge(selected, connectIndex, G);

            num ++;
        }//End correct number of rewires

        isConnected = isGraphConnected(G);

    }// end connection assertion
    for (auto &x: Network) x.Set_knn(Network);
    return Network;
}// end RewireNode
