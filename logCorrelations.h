/*
File for calculating stuff in log ccorrelation space
*/

std::vector<double> degreeDistributionFromLogCorrelations(
    std::vector<std::vector<double>> logMatrix,
    std::vector<int> kToLConversion, double avgDeg){
    
    std::vector<double> degreeDistribution;
    degreeDistribution.emplace_back(0.0);
    int kappaRange = kToLConversion.size();
    int maxDegree = kToLConversion[kappaRange - 1];
    int kappa = 1;
    double probability = 0.0;
    double PkSum = 0.0;
    for (int k = 1; k < maxDegree; k++){
        kappa = int(floor(log2(double(k))));
        probability = 0.0; 
        for (int kappaPrime = 0; kappaPrime < kappaRange; kappaPrime++){
            probability += logMatrix[kappa][kappaPrime]/
                                        pow(2.0, double(kappa));
        }// end P(k) calculation
        probability *= avgDeg/double(k);
        degreeDistribution.emplace_back(probability);
        PkSum += probability;
    }// end degree filling loop
    //std::cout << "Sum of Degree Distribution: " << PkSum << "\n";

    return degreeDistribution;
}// end degree distribution calculation

// Area ratio stuff (convert from log to lin correlations)
double areaRatio(std::array<int, 2> lUp, std::array<int, 2> lDown,
    std::vector<int> kToLConversion){
    
    // Primed coordinates are second (index 1)
    int lUpPrimeRange = kToLConversion[lUp[1] + 1] - kToLConversion[lUp[1]];
    int lDownPrimeRange = kToLConversion[lDown[1] + 1] -
        kToLConversion[lDown[1]];


    // Non primed coordinates first (index 0)
    int lUpRange = kToLConversion[lUp[0] + 1] - kToLConversion[lUp[0]];
    int lDownRange = kToLConversion[lDown[0] + 1] - kToLConversion[lDown[0]];
    
    //std::cout << "up*up'/down/down': " << lUpRange << ", " <<
    //    lUpPrimeRange << ", " << lDownRange << ", " << lDownPrimeRange << "\n";

    double ratio = double(lUpPrimeRange*lUpRange)/
        double(lDownPrimeRange*lDownRange);
    
    return ratio;

}// end area ratio

double dpSelection(double dp, double pUp, double pDown,
    std::array<int, 2> lUp, std::array<int, 2> lDown,
    std::vector<int> kToLConversion, bool strictlyUp){

    double aRatio = areaRatio(lUp, lDown, kToLConversion);

    double minDp = 0.0;
    if (strictlyUp){
        minDp = std::min({pDown, dp});
    }
    else {
        minDp = std::min({pDown, pUp, dp});
        if (aRatio > 1){
            minDp /= aRatio;
        }
        else {
            minDp *= aRatio;
        }
    }
    // if it's along the diagonal have to divide by 2
    if (lUp[0] == lUp[1] || lDown[0] == lDown[1]){
        minDp = minDp/2.0;
    }

    return minDp;

}// end dp Selection

double logChangeWeight(std::array<int, 2> lUp, std::array<int, 2> lDown,
    std::vector<int> kToLConversion){
    
    double aRatio = areaRatio(lUp, lDown, kToLConversion);
    double weight = 0.0;

    double upContribution = 0.0;
    double downContribution = 0.0;

    // Loop for the transpose move
    for (int i = 0; i < 2; i++){
    
        int rangeIndex = i;
        int sumIndex = (i+1)%2;

        // Ranges: simplify sum across indices k' since contribution is
        // dependent only on k
        int lUpRange = kToLConversion[lUp[rangeIndex] + 1] -
            kToLConversion[lUp[rangeIndex]];
        int lDownRange = kToLConversion[lDown[rangeIndex] + 1] -
            kToLConversion[lDown[rangeIndex]];
        
        // Do the part dependent on k 
        for (int j = kToLConversion[lUp[sumIndex]];
            j < kToLConversion[lUp[sumIndex] + 1]; j++){
            upContribution += double(lUpRange)/double(j);
        }// end up contribution

        for (int j = kToLConversion[lDown[sumIndex]];
            j < kToLConversion[lDown[sumIndex] + 1]; j++){
            downContribution += double(lDownRange)/double(j);
        }// end down Contribution

    }// end contribution calculation (Looped for transpose as well)

    return upContribution / aRatio - downContribution;

}// end logPkWeight

// THE MEAT AND POTATOES
void changeLogCorrelationMatrix(std::vector<std::vector<double>>
    &logCorrelationMatrix, std::vector<int> kToLConversion, double dp){
    
    int size = kToLConversion.size() - 1;

    std::vector<std::array<int, 2>> nonZeroCorrelations;

    //std::cout << "Initializing SelectedPoints\n";
    std::vector<std::array<int, 2>> selectedPoints;
    selectedPoints.reserve(9);

    // initialize lVector structues (and transpose) (points (l', l))
    std::array<int, 2> lVector = {0,0}; 
    std::array<int, 2> lVectorT = {0,0}; 

    // put (0,0) in points to eliminate it as a legal move
    selectedPoints.emplace_back(lVector);


    //std::cout << "Recording Non Zero Points\n";
    // record the points (k, k') where the degree correlations are nonzero
    int numNonZero = 0;
    for (int i = 0; i < size; i++){

        for (int j = 0; j <= i; j++){
            if (logCorrelationMatrix[i][j] > 0){
                numNonZero++;
                lVector = {i, j};
                nonZeroCorrelations.emplace_back(lVector);
                if (i != j){
                    lVectorT = {j, i};
                    nonZeroCorrelations.emplace_back(lVectorT);
                }// symmetric move
            }// non zero correlation

        }// end column fill

    }// end filling nonzero correlations

    if (numNonZero < 4) std::cout << "LESS THAN 4 NON Zero correlatons\n";
    // Select first point To add correlation Juice: anywhere (besides origin
    // which here is indices l', l = 0,0 (since k[0] = 1).
    //std::cout << "getting first kup vector\n";
    int lUp = int(RNG(g1) * size);
    int lUpPrime = int(RNG(g1) * size);
    while (lUp == 0 && lUpPrime == 0){
        lUp = int(RNG(g1) * size);
        lUpPrime = int(RNG(g1) * size);
    }// dont let it attach degree 1 to degree 1

    double upProbability = logCorrelationMatrix[lUp][lUpPrime];
    std::array<int, 2> lUpVector = {lUp, lUpPrime};

    lVector = {lUp, lUpPrime};
    selectedPoints.emplace_back(lVector);
    lVectorT = {lUpPrime, lUp};
    selectedPoints.emplace_back(lVectorT);

    // Select First point to steal Correlation: anywhere of k, k' [1, kMax]
    //std::cout << "first k down points\n";
    int lDown = int(RNG(g1) * size);
    int lDownPrime = int(RNG(g1) * size);
    double downProbability = logCorrelationMatrix[lDown][lDownPrime];
    // check that theres correlation to steal and that its not coming from a
    // point that is already being modified
    while (downProbability  == 0 ||
        inPrevious(lDown, lDownPrime, selectedPoints)){
        lDown = int(RNG(g1) * size);
        lDownPrime = int(RNG(g1) * size);
        downProbability = logCorrelationMatrix[lDown][lDownPrime];
    }// end good point assertion
    std::array<int, 2> lDownVector = {lDown, lDownPrime};

    // record that this point is being changed.
    //std::cout << "Record good kdown points\n";
    lVector = {lDown, lDownPrime};
    lVectorT = {lDownPrime, lDown};
    selectedPoints.emplace_back(lVector);
    selectedPoints.emplace_back(lVectorT);

    // The key quantities for the first move
    //std::cout << "Calculating important stuff (area etc)\n";
    double firstAreaRatio = areaRatio(lUpVector, lDownVector, kToLConversion);
    double firstMoveWeight = logChangeWeight(lUpVector, lDownVector,
        kToLConversion);
    double firstDeltaP = dpSelection(dp, upProbability, downProbability,
        lUpVector, lDownVector, kToLConversion, true);
    //std::cout << "Up Prob: " << upProbability << ", Down Prob: " <<
    //    downProbability << ", dP1: " << firstDeltaP << "\n";

    // ########## SECOND MOVE TIME ###############


    lUp = int(RNG(g1) * size);
    lUpPrime = int(RNG(g1) * size);

    upProbability = logCorrelationMatrix[lUp][lUpPrime];
    while (upProbability  == 0 ||
        inPrevious(lUp, lUpPrime, selectedPoints)){
        lUp = int(RNG(g1) * size);
        lUpPrime = int(RNG(g1) * size);
        upProbability = logCorrelationMatrix[lUp][lUpPrime];
    }// end good point assertion
 
    lUpVector = {lUp, lUpPrime};

    lVector = {lUp, lUpPrime};
    selectedPoints.emplace_back(lVector);
    lVectorT = {lUpPrime, lUp};
    selectedPoints.emplace_back(lVectorT);

    // Select First point to steal Correlation: anywhere of k, k' [1, kMax]
    //std::cout << "first k down points\n";
    lDown = int(RNG(g1) * size);
    lDownPrime = int(RNG(g1) * size);
    downProbability = logCorrelationMatrix[lDown][lDownPrime];
    // check that theres correlation to steal and that its not coming from a
    // point that is already being modified
    while (downProbability  == 0 ||
        inPrevious(lDown, lDownPrime, selectedPoints)){
        lDown = int(RNG(g1) * size);
        lDownPrime = int(RNG(g1) * size);
        downProbability = logCorrelationMatrix[lDown][lDownPrime];
    }// end good point assertion
    lDownVector = {lDown, lDownPrime};

    // record that this point is being changed.
    //std::cout << "Record good kdown points\n";
    lVector = {lDown, lDownPrime};
    lVectorT = {lDownPrime, lDown};
    selectedPoints.emplace_back(lVector);
    selectedPoints.emplace_back(lVectorT);

    // The key quantities for the first move
    double secondAreaRatio = areaRatio(lUpVector, lDownVector, kToLConversion);
    double secondMoveWeight = logChangeWeight(lUpVector, lDownVector,
        kToLConversion);
    double secondDeltaP = dpSelection(dp, upProbability, downProbability,
        lUpVector, lDownVector, kToLConversion, false);
    //std::cout << "Up Prob: " << upProbability << ", Down Prob: " <<
    //    downProbability << ", dP2: " << secondDeltaP << "\n";


    // Ratio between the factors affecting conservation
    double weightRatio = secondMoveWeight/firstMoveWeight;

    // Limitations on deltaPs considering the weight ratio
    if (firstDeltaP < abs(weightRatio * secondDeltaP) ){
        //std::cout << "P1 < |aP2|\n";
        secondDeltaP = -1.0 * firstDeltaP/weightRatio;
    }
    else if (weightRatio > 0){
        //std::cout << "P1 > |aP2| And a > 0\n";
        secondDeltaP *= -1.0;
        firstDeltaP = abs(weightRatio * secondDeltaP);
    }
    else{
        //std::cout << "Else Case\n";
        firstDeltaP = abs(weightRatio * secondDeltaP);
    }

    // Do the changes to the correlation matrix
    double deltaP = firstDeltaP;
    double aRatio = firstAreaRatio;
    int l = 0;
    int lPrime = 0;
    int sign = 1;
    // loop over the l indices (selected points except first point is (0,0))
    for (int i = 1; i < selectedPoints.size(); i++){
        if (i <= 4){
            deltaP = firstDeltaP; 
            if (deltaP < 0) std::cout << "First Move dp less than 0!\n";
            aRatio = firstAreaRatio;
            if (aRatio < 0) std::cout << "First Move area less than 0!\n";
        }
        else {
            deltaP = secondDeltaP; 
            aRatio = secondAreaRatio;
        }
        // Sign goes: ++--++--
        sign = 2 * ( int((i-1)/2 + 1) % 2) - 1;
        // all + have deltap/area ratio
        if (sign > 0){
            deltaP /= aRatio;
        }
        l = selectedPoints[i][0];
        lPrime = selectedPoints[i][1];
        /*
        std::cout << "Adding " << sign*deltaP << " to the point " << l <<
            ", " << lPrime << " which only has " << 
            logCorrelationMatrix[l][lPrime] << ".\n";
        */
        logCorrelationMatrix[l][lPrime] += sign * deltaP;


    }// end correlation Matrix changes
 
 
}// end change matrix code

std::vector<std::vector<double>> logToLinCorrelationMatrix(
    std::vector<std::vector<double>> logCorrelationMatrix,
    std::vector<int> kToLConversion){
    
    int maxDegree = kToLConversion[kToLConversion.size()-1];
    //2d vector of 0s
    std::vector<std::vector<double>> correlationMatrix (maxDegree,
        std::vector<double> (maxDegree, 0.0));

    int logSize = logCorrelationMatrix.size();
    
    int lPrime = 0;
    int l = 0;
    // loop over the prime coordinates
    for (int kPrime = 1; kPrime < maxDegree; kPrime ++){
        l = 0;
        // if kprime reaches the next block of degrees update lprime
        if (kPrime == kToLConversion[lPrime + 1]){
            lPrime += 1;
        }
        for (int k = 1; k < maxDegree; k++){
            if (k == kToLConversion[l+1]){
                l += 1;
            }

            if ((l >= logSize) || (lPrime > logSize)){
                correlationMatrix[k][kPrime] = 0.0;
            }
            else{
                correlationMatrix[k][kPrime] = logCorrelationMatrix[l][lPrime];
            }
        }// end k loop


    }// end kPrime loop

    return correlationMatrix;

}// end log to lin transform
