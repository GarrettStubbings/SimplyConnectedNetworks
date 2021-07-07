/* This file builds a network that is the "ideal" longevity network as
 * determined by the previous scale free analysis. Network consists of two
 * enormous hubs, connected only by direct link between hubs */

std::vector<Node> SuperHubby (const int N){

    std::vector<Node> Network;
    Network.reserve(N);

    for(int i = 0; i < N; i++) {

        Network.emplace_back(i);

    }



    Network[0].ConnectNode(Network[1]);
    for (int i = 2; i < N; i++){
        if (i%2 == 0){
            Network[i].ConnectNode(Network[0]);
        }//end even case
        else {
            Network[i].ConnectNode(Network[1]);
        }//end odd case
        
        //fill up extra connections
        if ( i > 3){
            Network[i].ConnectNode(Network[i-2]);
        }//end periferal connections
    }//end Island Formation for loop

    //Finish peripheral connections
    Network[2].ConnectNode(Network[N-2]);
    Network[3].ConnectNode(Network[N-1]);

    //Calculate <k>
    double calcAvgDeg = 0.0;
    int totalK = 0;
    for (int i = 0; i < N; i++){
        totalK += Network[i].k;
    }//end avg k calculation
    calcAvgDeg = (double) totalK / N;
    std::cout << "Average Degree is: " << calcAvgDeg << std::endl;
    std::cout << "Island Degrees are: " << Network[0].k << " And " << 
        Network[1].k << std::endl;
    std::cout << "All others are degree: " << Network[20].k << std::endl;

    std::cout << "Network size is: " << Network.size() << std::endl;
    for (auto &x: Network) x.Set_knn(Network);
    return Network;
}// end Super Hubby
