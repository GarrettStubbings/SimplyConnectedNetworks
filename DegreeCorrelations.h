/* File for correlation matrix stuff

Mostly for the correlation matrix changing algorithm

*/


double degreeContribution(int k, int kPrime){
    return 1.0/double(k) + 1.0/double(kPrime);
}// end contribution to the normalization constraint

bool inPrevious(int k, int kPrime, std::vector<std::array<int, 2>> previous){

    std::array<int, 2> kArray = {k, kPrime};
    for (auto x: previous){
        if (kArray == x){
            return true;
        }// found kArray in previous k list
    }// end check

    return false;

}// check if the point k, k' is already being modified

std::vector<double> calculateDegreeDistribution(
    std::vector<std::vector<double>> correlationMatrix, double avgDeg){
    
    double degreeDistributionSum = 0.0;

    if (correlationMatrix[1][1] != 0){
        std::cout << "P[1,1] = " << correlationMatrix[1][1] << "\n";
    }

    std::vector<int> degrees;

    int kMax = correlationMatrix.size();
    if (kMax < 1){
        std::cout << "correlation matrix is busted\n";
    }
    //std::cout << "KMAX IS: " << kMax << "\n";
    std::vector<double> degreeDistribution;
    degreeDistribution.emplace_back(0.0);
    degreeDistribution.reserve(kMax);

    bool badValue = false;

    // use P(k) = avg_deg/k sum_k' P(k',k)
    for (int k = 1; k < kMax; k++){
        double pk = 0.0;
        for (auto x: correlationMatrix[k]){
            if (std::isnan(x) || x < 0.0 || x > 1){
                badValue = true;
                //std::cout << "Bad value in Pkk: " << x << "\n";
            }
            pk += avgDeg/k * x;
        }// end specific k value Pk
        degreeDistribution.emplace_back(pk);
        degreeDistributionSum += pk;
        if (pk < 0){
            std::cout << "\n\nNegative Pk\n\n";
            
        }
    }// end degree distribution filling
    if (badValue){
        for (auto x: correlationMatrix){
            PrintVectorDouble(x);
        }
    }
    //std::cout << "Size of Degree Distribution: " << degreeDistribution.size()
    //    << "\n";
    //std::cout << "Degree Distribution Sum: " << degreeDistributionSum << "\n";
    //PrintVectorDouble(degreeDistribution);
    /*
    double averageDegree = 0;
    for (int k = 1; k < degreeDistribution.size(); k++){
        averageDegree += k*degreeDistribution[k];
    }// end average degree calculation
    std::cout << "Average Degree: " << averageDegree << "\n";
    double sumOfCorrelationMatrix = 0;
    for (auto x: correlationMatrix){
        //std::cout << "Length of row: " << x.size() << "\n";
        for (auto y: x){
            sumOfCorrelationMatrix += y;
            //std::cout << y << " ";
        }
        //std::cout << "\n";
    }
    //std::cout << "Sum Of Correlation Matrix: " << sumOfCorrelationMatrix <<
    //    "\n";
    */
    return degreeDistribution;

}// end calculate degree distribution

void changeCorrelationMatrix(
    std::vector<std::vector<double>> &correlationMatrix, double dp,
    int extension){
    
    int kMax = correlationMatrix.size() - 1;

    /*
    for (auto x: correlationMatrix){
        PrintVectorDouble(x);
        std::cout << "Size of Row: " << x.size() << "\n";
    }
    */
    std::vector<std::array<int, 2>> nonZeroCorrelations;
    std::array<int, 2> kVector = {0,0}; 
    std::array<int, 2> kVectorT = {0,0}; 
    //std::cout << "Recording Non Zero Points\n";
    // record the points (k, k') where the degree correlations are nonzero
    int numNonZero = 0;
    for (int i = 0; i <= kMax; i++){

        for (int j = 0; j <= i; j++){
            if (correlationMatrix[i][j] > 0){
                numNonZero++;
                kVector = {i, j};
                nonZeroCorrelations.emplace_back(kVector);
                if (i != j){
                    kVectorT = {j, i};
                    nonZeroCorrelations.emplace_back(kVectorT);
                }// symmetric move
            }// non zero correlation

        }// end column fill

    }// end filling nonzero correlations

    if (numNonZero < 4) std::cout << "LESS THAN 4 NON Zero correlatons\n";
    //std::cout << "Initializing SelectedPoints\n";
    std::vector<std::array<int, 2>> selectedPoints;
    selectedPoints.reserve(8);

    // Select first point To add correlation Juice: anywhere 1 to kMax + 1
    //std::cout << "getting first kup vector\n";
    int kUp = 1 + int(RNG(g1) * (kMax + extension) );
    int kPrimeUp = 1 + int(RNG(g1) * (kMax + extension) );
    while (kUp == 1 && kPrimeUp == 1){
        kUp = 1 + int(RNG(g1) * (kMax + extension) );
        kPrimeUp = 1 + int(RNG(g1) * (kMax + extension) );
    }// dont let it attach degree 1 to degree 1

    //std::cout << "Calculating First kup contribution\n";
    double upContribution = degreeContribution(kUp, kPrimeUp);

    // If either is larger than kMax extend the correlation matrix
    if (kUp > kMax || kPrimeUp > kMax){
        //std::cout << "Extending Matrix\n";
        // extend rows beyond kmax
        for (auto &x: correlationMatrix){
            x.emplace_back(0.0);
        }// end row extension

        //extend columns beyond kmax (add an extra row)
        std::vector<double> newRow(kMax + 2, 0.0);
        correlationMatrix.emplace_back(newRow);
    }// end matrix extension
    
    // record that this point is being changed.
    //std::cout << "Recording kup vector in selected points\n";
    kVector = {kUp, kPrimeUp};
    kVectorT = {kPrimeUp, kUp};
    selectedPoints.emplace_back(kVector);
    selectedPoints.emplace_back(kVectorT);

    // Select First point to steal Correlation: anywhere of k, k' [1, kMax]
    //std::cout << "first k down points\n";
    int kDown = 1 + int(RNG(g1) * (kMax) );
    int kPrimeDown = 1 + int(RNG(g1) * (kMax) );
    
    // check that theres correlation to steal and that its not coming from a
    // point that is already being modified
    while (correlationMatrix[kDown][kPrimeDown] == 0 ||
        inPrevious(kDown, kPrimeDown, selectedPoints) ||
        (kDown == 1 && kPrimeDown == 1)){

        kDown = 1 + int(RNG(g1) * (kMax) );
        kPrimeDown = 1 + int(RNG(g1) * (kMax) );

    }// end good point assertion

    // record that this point is being changed.
    //std::cout << "Record good kdown points\n";
    kVector = {kDown, kPrimeDown};
    kVectorT = {kPrimeDown, kDown};
    selectedPoints.emplace_back(kVector);
    selectedPoints.emplace_back(kVectorT);


    double downContribution = degreeContribution(kDown, kPrimeDown);
    
    // the degree dependent weights of any move deltaP for the first Move
    double firstMoveFactor = upContribution - downContribution;
    // Delta P itself
    double firstMoveDeltaP = correlationMatrix[kDown][kPrimeDown];
    //std::cout << "Degree Correlation of first move: " << firstMoveDeltaP <<
    //    "\n";

    // account for double duty on the diagonal
    if (kDown == kPrimeDown) firstMoveDeltaP /= 2;
    // limit amout of correlation juie that can be moved
    if (firstMoveDeltaP > dp) firstMoveDeltaP = dp;

    //std::cout << "First move is taking " << firstMoveDeltaP <<
    //    " From k,k' with factor " << firstMoveFactor << "\n";

    // Select Second point to Add Correlation: anywhere of k, k' [1, kMax]
    // NOTE: it is written as adding, but the move may need to be done in
    // reverse (to assert that deltaP 1 is postive so we can extend the matrix)
    kUp = 1 + int(RNG(g1) * (kMax) );
    kPrimeUp = 1 + int(RNG(g1) * (kMax) );
    
    // check that theres correlation to steal and that its not coming from a
    // point that is already being modified
    while (correlationMatrix[kUp][kPrimeUp] == 0 ||
        inPrevious(kUp, kPrimeUp, selectedPoints) ||
        (kUp == 1 && kPrimeUp == 1)){

        kUp = 1 + int(RNG(g1) * (kMax) );
        kPrimeUp = 1 + int(RNG(g1) * (kMax) );

    }// end good point assertion

    // record that this point is being changed.
    kVector = {kUp, kPrimeUp};
    kVectorT = {kPrimeUp, kUp};
    selectedPoints.emplace_back(kVector);
    selectedPoints.emplace_back(kVectorT);


    upContribution = degreeContribution(kUp, kPrimeUp);
    
    // The second move's delta P is tentatively set here, The up move and the
    // down move could both be limiting factors since it's unclear which way
    // is which
    double secondMoveDeltaP = correlationMatrix[kUp][kPrimeUp];
    //std::cout << "Degree Correlation of Second move: " << secondMoveDeltaP <<
    //    "\n";


    // account for double duty on the diagonal
    if (kUp == kPrimeUp) secondMoveDeltaP /= 2;
    // limit amout of correlation juice that can be moved
    if (secondMoveDeltaP > dp) secondMoveDeltaP = dp;

    //std::cout << "Second move is taking " << secondMoveDeltaP << "\n"; 


    // Select Second point to steal Correlation: anywhere of k, k' [1, kMax]
    kDown = 1 + int(RNG(g1) * (kMax) );
    kPrimeDown = 1 + int(RNG(g1) * (kMax) );
    
    // check that theres correlation to steal and that its not coming from a
    // point that is already being modified
    while (correlationMatrix[kDown][kPrimeDown] == 0 ||
        inPrevious(kDown, kPrimeDown, selectedPoints) ||
        (kDown == 1 && kPrimeDown == 1)){

        kDown = 1 + int(RNG(g1) * (kMax) );
        kPrimeDown = 1 + int(RNG(g1) * (kMax) );

    }// end good point assertion

    // record that this point is being changed.
    kVector = {kDown, kPrimeDown};
    kVectorT = {kPrimeDown, kDown};
    selectedPoints.emplace_back(kVector);
    selectedPoints.emplace_back(kVectorT);


    downContribution = degreeContribution(kDown, kPrimeDown);
    
    // the degree dependent weights of any move deltaP for the first Move
    double secondMoveFactor = upContribution - downContribution;

    // Delta P 2: Have to check if it's smaller than the up move delta p
    //std::cout << "Degree Correlation of Other Second move: " <<
    //    correlationMatrix[kDown][kPrimeDown] << "\n";


    // Diagonal case
    if (kDown == kPrimeDown && secondMoveDeltaP >
        correlationMatrix[kDown][kPrimeDown]/2){
        //std::cout << "Diagonal and Limitting case\n";
        secondMoveDeltaP = correlationMatrix[kDown][kPrimeDown]/2.0;
    }// end diagonal case
    else if (secondMoveDeltaP > correlationMatrix[kDown][kPrimeDown]){
        //std::cout << "Non Diagonal but still limitting\n";
        secondMoveDeltaP = correlationMatrix[kDown][kPrimeDown];
    }// end non diagonal but still limitting case
    
    // final limitation on deltaP 2
    if (secondMoveDeltaP > dp) secondMoveDeltaP = dp;
    //std::cout << "Second move is taking " << secondMoveDeltaP <<
    //    " From k,k' with factor " << secondMoveFactor << "\n";



    // Ratio between the factors affecting conservation
    double factorRatio = secondMoveFactor/firstMoveFactor;

    // Limitations on deltaPs considering the factor ratio
    if (firstMoveDeltaP < abs(factorRatio * secondMoveDeltaP) ){
        //std::cout << "P1 < |aP2|\n";
        secondMoveDeltaP = -1.0 * firstMoveDeltaP/factorRatio;
    }
    else if (factorRatio > 0){
        //std::cout << "P1 > |aP2| And a > 0\n";
        secondMoveDeltaP *= -1.0;
        firstMoveDeltaP = abs(factorRatio * secondMoveDeltaP);
    }
    else{
        //std::cout << "Else Case\n";
        firstMoveDeltaP = abs(factorRatio * secondMoveDeltaP);
    }

    // Calculate the actual effect of this move (should be 0)
    double moveImpact = firstMoveDeltaP*firstMoveFactor + secondMoveDeltaP * 
        secondMoveFactor;

    //std::cout << "The Move had " << moveImpact << " Effect on Conservation\n";

    double deltaP = firstMoveDeltaP;
    int k = 0;
    int kPrime = 0;
    int sign = 1;
    // Perform all the changes to the correlation Matrix itself
    //std::cout << "Performing Changes\n";
    for (int i = 0; i < selectedPoints.size(); i++){
        if (i == 4) deltaP = secondMoveDeltaP; 
        // Sign goes: ++--++--
        sign = 2 * ( int(i/2 + 1) % 2) - 1;
        k = selectedPoints[i][0];
        kPrime = selectedPoints[i][1];
        /*
        std::cout << "Adding " << sign*deltaP << " to the point " << k <<
            ", " << kPrime << " which only has " << 
            correlationMatrix[k][kPrime] << ".\n";
        */
        
        correlationMatrix[k][kPrime] += sign * deltaP;

    }// end correlation Matrix changes
    //std::cout << "Changes Performed\n";

}// end correlation matrix

std::vector<std::vector<double>> correlationsFromNetwork(
    std::vector<Node> Network, double avgDeg, int N, int kMax){

    double increment = 0.5/(double(N) * avgDeg);
    
    // initialize the correlation matrix to 0s (max deg + 1 by max deg + 1)
    std::vector<std::vector<double>>
        correlationMatrix(kMax, std::vector<double> (kMax, 0.0));

    for (auto x: Network){
        int k = x.k;
        for (auto n: x.Neighbours){
            int kPrime = Network[n].k;
            correlationMatrix[k][kPrime] += increment;
            correlationMatrix[kPrime][k] += increment;
        }
    }
    return correlationMatrix;
}// end calculating correlation matrix from network


// MOVED TO BERTOTTI.H BECAUSE IT DOESNT BELONG HERE
/*
std::vector<Node> generateNetwork(
    std::vector<std::vector<double>> &correlationMatrix,
    std::vector<double> &degreeDistribution, double dP,
    int numChanges, double avgDeg, int N, int rewiresPerNode){

    std::vector<std::vector<double>> tempCorrelationMatrix;

    std::vector<Node> TempNetwork;
    int extension = 0;
    std::vector<double> tempDegreeDistribution;
    bool initFailure = true;
    while (initFailure){
        initFailure = false;
        // correlation matrix to modify
        tempCorrelationMatrix = correlationMatrix;

        //modify it 
        for (int j = 0; j < numChanges; j++){
            if (j < 1){
                extension = 1;
            }
            else {
                extension = 0;
            }
            changeCorrelationMatrix(tempCorrelationMatrix, dP,
                extension);
        }// End degree correlation modification

        tempDegreeDistribution =
            calculateDegreeDistribution(tempCorrelationMatrix, avgDeg);
        // generate the network
        //std::cout << "After Correlation Changes\n";
        TempNetwork = networkFromCorrelations(
            tempDegreeDistribution, tempCorrelationMatrix, N, rewiresPerNode,
            avgDeg, initFailure);
        if (TempNetwork.size() < N){
            std::cout << "Tiny Baby Network, InitFailure: " << initFailure <<
                "\n";
            std::cout << "FAILURE\n";
        }
        if (initFailure){
            std::cout << "\n\n\n\n\n\n\nFailure To Initialize\n\n\n\n\n\n\n";
        }
        //std::cout << "InitFailure: " << initFailure << "\n";

    }// end network generation (asserting solvable Pkk)
    degreeDistribution = tempDegreeDistribution;
    correlationMatrix = tempCorrelationMatrix;

    return TempNetwork;

}// end network generation

*/
