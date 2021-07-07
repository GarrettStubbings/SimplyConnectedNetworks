//find node to attach given connectivity
int FindNode(EventTree<double> &tree) {

    double val = tree.TotalSum * RNG(g1);

    return  tree.Search(val);

}



/*
Create a scale free network from the BA algorithm with shifted linear attachment to allow for alpha != 3. 
This uses a tree (https://en.wikipedia.org/wiki/Fenwick_tree) to record partial sums for faster updating and 
calculating of probabilities. See tree.h
*/
std::vector<Node> ScaleFreeNetwork(int N, int m, int m0, double alpha) {
    
    std::vector<Node> Network;
    Network.reserve(N+1); 
    EventTree<double> tree(N, 0);

    //set up initial nodes
    Network.emplace_back(0);
    double lambda = (alpha-3.0) * double(m);
    double SumOfDegrees = 0;

    for(int i = 1; i < m0; i++) {

        Network.emplace_back(i);

	//new node connects to all nodes before it
        for(int j = 0; j < i; j++) {

            Network[i].ConnectNode(Network[j]); //does bidirectional connection, i->j, j->i
            tree.Update(i,Network[i].k + lambda);
            tree.Update(j,Network[j].k + lambda);

        }

    }

    

    //create other nodes
    for(int i = m0; i < N; i++) {

	//store the m connections to make for the new node, then add them at the end
	std::vector<int> ConnectTo;
        for(int j = 0; j < m; j++) {

            int newConnection = FindNode(tree);

	    //find a new node if this connection already exists
	    while(std::find(ConnectTo.begin(), ConnectTo.end(), newConnection) != ConnectTo.end()) {
	        newConnection = FindNode(tree);
	    }
	    
	    ConnectTo.push_back(newConnection);
	    
        }

	//connect all of these to the new node.
        Network.emplace_back(i);
        for(auto x: ConnectTo) {
            Network[i].ConnectNode(Network[x]);
            tree.Update(i,Network[i].k + lambda);
            tree.Update(x,Network[x].k + lambda);
        }

    }

    return Network;

}
