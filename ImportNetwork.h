/* Builds a network from a csv file containing a list of edges*/
void ImportEdgeListNetwork(std::vector<Node> &Network, std::string fileName,
    int N){
    
    // the source thingy
    std::ifstream f;
    // open the file
    f.open(fileName);
    // check if the file is open
    if (!f.is_open()){
        std::cout << "Couldn't open file.\n";
    }//end file opening check

    // fill the network with nodes (IDs really)
    Network.reserve(N + 1);
    for (int i = 0; i < N; i++){
        Network.emplace_back(i);
    }//end filling the network IDs

    
    // declare line (whole line of file) and value (element of csv)
    std::string line, value;

    // Go through file and connect neighbours 
    while (std::getline (f, line)){

        std::vector<int> nodeIDs;
        // a stream of strings that is the line (value, value, value)
        std::stringstream s (line);
        
        // go through the string stream and pick out values
        // should just be node IDs on either side of the edge
        while (getline(s, value, ',')){
            //save the value to the row data
            nodeIDs.emplace_back(std::stoi (value));
        }//end row of data
        
        // connect the two nodes in the network
        Network[nodeIDs[0]].ConnectNode(Network[nodeIDs[1]]);

    }// end data pulling loop

    // close the file?
    f.close();
    //std::cout << "File Closed\n";


    // do the network stuff
    for (auto &x: Network) x.Set_knn(Network);
	
}//end import network


/* Builds a network from a csv file containing its adjacency matrix*/
void ImportedNetwork(std::vector<Node> &Network, std::string fileName){
    //std::vector<int> &MortalityNodes){

    
    // the source thingy
    std::vector<int> MortalityNodes;
    std::ifstream f;
    // open the file
    f.open(fileName);
    // check if the file is open
    if (!f.is_open()){
        std::cout << "Couldn't open file.\n";
    }//end file opening check
    
    // declare line (whole line of file) and value (element of csv)
    std::string line, value;

    // initialise the adjacency matrix
    std::vector<std::vector<int>> adjacencyMatrix;
    
    //pull the data from the file and fill the matrix
    while (std::getline (f, line)){

        // this row of data initialized
        std::vector<int> rowData;
        // a stream of strings that is the line (value, value, value)
        std::stringstream s (line);
        
        // go through the string stream and pick out values
        while (getline(s, value, ',')){
            //save the value to the row data
            rowData.emplace_back(std::stoi (value));
        }//end row of data
        
        //put the row data in the matrix
        adjacencyMatrix.emplace_back(rowData);

    }// end data pulling loop

    // close the file?
    f.close();
    //std::cout << "File Closed\n";

    //Begin building the network
    Parameters::N = adjacencyMatrix.size();
    //std::cout << "Parameters::N = " << Parameters::N << "\n";;

    Network.reserve(Parameters::N + 1);

    for (int i = 0; i < Parameters::N; i++){
        Network.emplace_back(i);
    }//end filling the network IDs

    /*
    std::cout << "Network Size is " << Network.size() << "\n";
    std::cout << "First 10 Elements:\n";
    for (int j = 0; j < 10; j++){
        //for every sub diagonal element
        for (int i = 0; i < 10; i++){

            std::cout << adjacencyMatrix[j][i] << " ";
        }//end row fill
        std::cout << "\n";
    }// end sample output
    */
    //Add neighbours from matrix
    //for every row
    for (int i = 0; i < Parameters::N; i++){
        //for every sub diagonal element
        for (int j = i + 1; j < Parameters::N; j++){
            //if the element is 1 a connection exists
            if (adjacencyMatrix[i][j] == 1){
                Network[i].ConnectNode(Network[j]);
            }//end conditional connection

        }//end row fill

    }// end filling neighbours
    adjacencyMatrix.clear();


    // do the network stuff
    for (auto &x: Network) x.Set_knn(Network);
	
    std::cout << "Find Mortality" << std::endl;
    MortalityNodes = FindMortalityNodes(Network);
    
    std::cout << "output Degrees" << std::endl;
    OutputNetworkDegrees(Network, Parameters::N);
    double calcAvgDeg = 0.0;
    int totalK = 0;
    for (int i = 0; i < Parameters::N; i++){
        totalK += Network[i].k;
    }//end avg k calculation
    calcAvgDeg = (double) totalK / Parameters::N;
    std::cout << "Average Degree is: " << calcAvgDeg << std::endl;

}//end import network

std::vector<int> importNodes(std::string fileName){

    std::ifstream f;
    // open the file
    f.open(fileName);
    // check if the file is open
    if (!f.is_open()){
        std::cout << "Couldn't open file.\n";
    }//end file opening check
    
    // declare line (whole line of file) and value (element of csv)
    std::string line, value;

    // initialise the column of Node IDs
    std::vector<int> NodeList;
    
    //pull the data from the file and fill the matrix
    while (std::getline (f, line)){

        // a stream of strings that is the line (value, value, value)
        std::stringstream s (line);
        
        // go through the string stream and pick out values
        while (getline(s, value, ',')){
            //save the value to the row data
            NodeList.emplace_back(std::stoi (value));
        }//end row of data
        
    }// end data pulling loop

    return NodeList;

}// end import node IDs

void import1DIntVector(std::vector<int> &vec,
    std::string fileName){

    // the source thingy
    std::ifstream f;
    // open the file
    f.open(fileName);
    // check if the file is open
    if (!f.is_open()){
        std::cout << "Couldn't open file.\n";
    }//end file opening check
    
    // thing to be read from file
    double data;
    f >> data;
    
    //pull the data from the file and fill the matrix
    while (!f.eof()){
        vec.emplace_back(int(data));
        f >> data;

    }// end data pulling loop

    // close the file?
    f.close();

}//end import degree distribution



// Import degree distribution stuff
void importDegreeDistribution(std::vector<double> &degreeDistribution,
    std::string fileName){

    // the source thingy
    std::ifstream f;
    // open the file
    f.open(fileName);
    // check if the file is open
    if (!f.is_open()){
        std::cout << "Couldn't open file.\n";
    }//end file opening check
    
    // thing to be read from file
    double data;
    f >> data;
    
    //pull the data from the file and fill the matrix
    while (!f.eof()){
        degreeDistribution.emplace_back(data);
        f >> data;

    }// end data pulling loop

    // close the file?
    f.close();

}//end import degree distribution

// Import the correlation matrix
void importCorrelationMatrix(std::vector<std::vector<double>>
    &correlationMatrix, std::string fileName){

    // the source thingy
    std::ifstream f;
    // open the file
    f.open(fileName);
    // check if the file is open
    if (!f.is_open()){
        std::cout << "Couldn't open file.\n";
    }//end file opening check
    
    // declare line (whole line of file) and value (element of csv)
    std::string line, value;

    double correlationSum = 0.0;

    //pull the data from the file and fill the matrix
    while (std::getline (f, line)){

        // this row of data initialized
        std::vector<double> rowData;
        // a stream of strings that is the line (value, value, value)
        std::stringstream s (line);
        
        // go through the string stream and pick out values
        while (getline(s, value, ',')){
            //save the value to the row data
            rowData.emplace_back(std::stod (value));
            correlationSum += std::stod(value);
        }//end row of data
        
        //put the row data in the matrix
        correlationMatrix.emplace_back(rowData);

    }// end data pulling loop
    std::cout << "Sum of Imported Matrix: " << correlationSum << "\n";

    // close the file?
    f.close();

}// end import correlation matrix
