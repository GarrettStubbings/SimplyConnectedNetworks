/* File to have statistical thingies for evaluating the goodness
of a specific Network 
Also has various calculations
*/
std::vector<std::vector<double>> denormalizeLogMatrix(
    std::vector<std::vector<double>> normalizedLogMatrix,
    std::vector<int> kToLConversion){

    int size = normalizedLogMatrix.size();
    std::vector<std::vector<double>> logMatrix = normalizedLogMatrix;
    for (int i = 0; i < size; i++){
        int rowWidth = kToLConversion[i + 1] - kToLConversion[i];
        for (int j = 0; j < size; j++){
            int columnHeight = kToLConversion[j+1] - kToLConversion[j];
            double area = double(columnHeight*rowWidth);
            logMatrix[i][j] = normalizedLogMatrix[i][j]/area;
        }// end loop over columns
    }// end loop over rows

    return logMatrix;

}

std::vector<double> compressLogCorrelations(
    std::vector<std::vector<double>> logCorrelationMatrix){

    std::vector<double> compressedLogCorrelations;

    int size = logCorrelationMatrix.size();

    for (int i = 0; i < size; i++){
        for (int j = 0; j <= i; j++){
            compressedLogCorrelations.emplace_back(
                logCorrelationMatrix[i][j]);
        }//end j loop
    }// end i loop
    return compressedLogCorrelations;
}

std::vector<std::vector<double>> normalizeLogMatrix(
    std::vector<std::vector<double>> logMatrix,
    std::vector<int> kToLConversion){

    int size = logMatrix.size();
    std::vector<std::vector<double>> normalizedLogMatrix = logMatrix;
    for (int i = 0; i < size; i++){
        int rowWidth = kToLConversion[i + 1] - kToLConversion[i];
        for (int j = 0; j < size; j++){
            int columnHeight = kToLConversion[j+1] - kToLConversion[j];
            double area = double(columnHeight*rowWidth);
            normalizedLogMatrix[i][j] = area*logMatrix[i][j];
        }// end loop over columns
    }// end loop over rows

    return normalizedLogMatrix;

}// end log matrix normalization

void padVector(std::vector<std::vector<double>> &vec, int desiredSize){
    int size = vec.size();
    int sizeDifference = desiredSize-size;
    for (auto &row: vec){
        for (int i = 0; i < sizeDifference; i++){
            row.emplace_back(0.0);
        }
    }// extend existing rows

    for (int i = 0; i < sizeDifference; i++){
        std::vector<double> zeros(desiredSize, 0.0);
        vec.emplace_back(zeros);
    }// add rows
}

double HellingerDistance2D(
    std::vector<std::vector<double>> p, std::vector<std::vector<double>> q){

    double distance = 0.0;
    double contribution = 0.0;
    int size = p.size();
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            contribution = (sqrt(p[i][j]) - sqrt(q[i][j]));
            distance += pow(contribution, 2.0);
        }// end j loop
    }// i loop
    distance = sqrt(distance)/sqrt(2.0);
    return distance;

}// end hellinger distance calculation for 2d

double KullbackLeiblerDivergence2D(
    std::vector<std::vector<double>> p, std::vector<std::vector<double>> q){
    int size = p.size();
    double divergence = 0.0;
    double contribution = 0.0;
    //std::cout << "Size: " << size << "\n";
    for (int i = 0; i < size; i++){
        //PrintVectorDouble(p[i]);
        //PrintVectorDouble(q[i]);
        for (int j = 0; j < size; j++){
           if (q[i][j] == 0 || p[i][j] == 0){
                continue;
            }
            //std::cout << "Pij: " << p[i][j] << ", Qij: " << q[i][j] << "\n";
            contribution = p[i][j] * log(p[i][j] / q[i][j]);
            //std::cout << "Contribution : " << contribution << "\n";
            divergence += contribution;       
        }// end j loop
    }// end loop over i
    return divergence;

}// end 2d kullback leibler divergence calculation


double KullbackLeiblerDivergence1D(
    std::vector<double> p, std::vector<double> q){
    
    double divergence = 0.0;
    for (int i = 0; i < p.size(); i++){
        if (q[i] == 0 || p[i] == 0) continue;
        divergence += p[i] * log(p[i] / q[i]);
    }// end loop
    return divergence;

}// end 1d kullback leibler divergence calculation

std::vector<std::vector<double>> logFromLinCorrelations(
    std::vector<std::vector<double>> correlationMatrix,
    std::vector<int> kToLConversion, std::string method = "functional"){
    
    int numPoints = kToLConversion.size();

    int kMax = kToLConversion[kToLConversion.size()-1];
    //std::cout << "Max degree: " << kMax << "\n";
    padVector(correlationMatrix, kMax);

    // initialize the logMatrix to be all 0s
    std::vector<std::vector<double>> logMatrix(numPoints,
            std::vector<double>(numPoints, 0.0));

    std::vector<int> kReduced = kToLConversion;
    kReduced.erase(kReduced.begin() + numPoints-1);
    double average = 0.0;
    int logKIndex = 0;
    int logKPrimeIndex = 0;
    int kRange = 0;
    int kPrimeRange = 0;
    double kArea = 0.0;
    for (int k: kReduced){
        if (k >= kMax) break;
        // calculate range of k values (how many equal p values)
        kRange = kToLConversion[logKIndex+1] - k;
        logKPrimeIndex = 0;

        //std::cout << "k range: " << kRange << "\n";
        for (int kPrime: kReduced){
            if (kPrime >= kMax) break;
            //std::cout << "k: " << k << ", k': " << kPrime << "\n";
            
            // how many k' values
            kPrimeRange = kToLConversion[logKPrimeIndex+1] - kPrime;
            //std::cout << "k' range: " << kPrimeRange << "\n";
            // "area" we're averaging over
            kArea = double(kPrimeRange*kRange);
            
            // perform the average over the k and k prime ranges
            for (int i = 0; i < kRange; i++){
                for (int j = 0; j < kPrimeRange; j++){
                    // add to the average
                    int dummyK = k + i;
                    int dummyKPrime = kPrime + j;
                    //std::cout << "Dummy K: " << dummyK << ", Dummy k': " <<
                    //    dummyKPrime << "\n";
                    double contribution =
                                        correlationMatrix[dummyK][dummyKPrime]; 
                    if (method == "normalized"){
                        average += contribution;
                    }
                    else {
                        average += contribution/kArea;
                    }
                }// end loop over k' range
            }// end loop over k range

            // assign the value in the log matrix to the average
            //std::cout << "Average is: " << average << "\n";
            logMatrix[logKIndex][logKPrimeIndex] = average;
            // reset the average
            average = 0.0;
            // move on to the next value
            logKPrimeIndex++;
        }// loop over k' values
        logKIndex++;
    }// loop over k values
    
    return logMatrix;

}// end log correlationa calculation (from correlations in linear space)

std::vector<std::vector<double>> addVector2D(std::vector<std::vector<double>>
    firstVector, std::vector<std::vector<double>> secondVector){
    int N = firstVector.size();
    int M = secondVector.size();
    int size = N;
    //std::cout << "First Vector Size: " << N << ", Second: " << M << "\n";
    // Have to Pad out any size discrepancies
    if (M > N){
        padVector(firstVector, M);
        size = M;
    }
    if (N > M) padVector(secondVector, N);

    std::vector<std::vector<double>> sumVector;
    sumVector.reserve(N);
    for (int i=0; i<size; i++){
        std::vector<double> row;
        row.reserve(size);
        for (int j=0; j<size; j++){
            double componentSum = firstVector[i][j] + secondVector[i][j];
            row.emplace_back(componentSum);
        }// inner loop over j
        sumVector.emplace_back(row);
    }// loop over i
    return sumVector;
}// end 2d vector addition

std::vector<double> addVector(std::vector<double> vec1, std::vector<double>
    vec2){
    
    std::vector<double> vectorSum;
    // Adjust vectors to be the same size (pad with 0s)
    while (vec1.size() > vec2.size()){
        vec2.emplace_back(0.0);
    }
    while (vec2.size() > vec1.size()){
        vec1.emplace_back(0.0);
    }
    int size = vec2.size();
    // Add em up
    for (int i = 0; i < size; i++){
        vectorSum.emplace_back(vec1[i] + vec2[i]);
    }
    return vectorSum;
}

// xlnx function for dealing with 0 case
double xlnx(double x){
    if (x > 0){
        return x * log (x);
    }
    else if (abs(x) < 1.0e-13){
        return 0.0;
    }
    else {
        std::cout << "Value Less Than Zero!\nBad value is: " << x << "\n";
        return 0;
    }
}//end xlnx


std::vector<std::vector<double>> calculatePij(
    std::vector<std::vector<double>> correlationMatrix,
    std::vector<double> degreeDistribution,
    std::vector<int> degreeSequence){
    
    int N = degreeSequence.size();
    //std::cout << "Correlation Matrix size: " << correlationMatrix.size() <<
    //    "\n";
    //PrintVectorInt(degreeSequence);
    std::vector<std::vector<double>> pij;
    pij.reserve(N);
    for (int i = 0; i < N; i++){
        std::vector<double> row;
        row.reserve(N);
        int ki = degreeSequence[i];
        for (int j = 0; j < N; j++){
            int kj = degreeSequence[j];
            double p = correlationMatrix[ki][kj];
            p = p/(degreeDistribution[ki] * degreeDistribution[kj]);
            p = p/double(N*N);
            row.emplace_back(p);
        }// end column loop
        pij.emplace_back(row);        
    }// end row loop
    
    return pij;
}// end Pij calculation

double entropy1D(std::vector<double> probs){
    double s = 0.0;
    for (auto x: probs){
        s -= xlnx(x);
    }
    return s;
}// end 1D shannon entropy calcculation

double entropy2D(std::vector<std::vector<double>> probs){
    double s = 0.0;
    for (auto x: probs){
        for (auto y: x){
            s -= xlnx(y);
        }// end column loop
    }// end row loop
    return s;
}// end 1D shannon entropy calcculation

void calculateDistributionEntropies(
    std::vector<double> &pkEntropies,
    std::vector<double> &pkkEntropies,
    std::vector<double> &logPkkEntropies,
    std::vector<std::vector<double>> correlationMatrix,
    std::vector<double> degreeDistribution,
    std::vector<std::vector<double>> logCorrelationMatrix){

    pkEntropies.emplace_back(entropy1D(degreeDistribution));
    pkkEntropies.emplace_back(entropy2D(correlationMatrix));
    logPkkEntropies.emplace_back(entropy2D(logCorrelationMatrix));
 
}// end distribution entropy claculation (excludes Pij (needs a network))

void calculateEntropies(std::vector<Node> &Network,
    std::vector<double> &pijEntropies,
    std::vector<double> &pkEntropies,
    std::vector<double> &pkkEntropies,
    std::vector<double> &logPkkEntropies,
    std::vector<std::vector<double>> correlationMatrix,
    std::vector<double> degreeDistribution,
    std::vector<std::vector<double>> logCorrelationMatrix){
    //std::cout << "\n\nCalculating Entropies: N = " << Network.size() << "\n\n";
    //std::cout << "\n Calculating Entropies\n\n";
    //std::cout << "Correlation Matrix:\n";
    //for (auto x: correlationMatrix){
    //    PrintVectorDouble(x);
    //}
    //std::cout << "Degree Distributions:\n";
    //PrintVectorDouble(degreeDistribution);

    
    std::vector<int> degreeSequence;
    for (auto x: Network){
        degreeSequence.emplace_back(x.k);
    }

    std::vector<std::vector<double>> pij;
    //std::cout << "Calculating Pij \n";
    pij = calculatePij(correlationMatrix, degreeDistribution, degreeSequence);
    //std::cout << "Finished Calculation\n";
    pijEntropies.emplace_back(entropy2D(pij));
    pkEntropies.emplace_back(entropy1D(degreeDistribution));
    pkkEntropies.emplace_back(entropy2D(correlationMatrix));
    logPkkEntropies.emplace_back(entropy2D(logCorrelationMatrix));
    //std::cout << "Finished Recording\n";
}// end calculate entropies  



// Median Calculation

double vectorMedian(std::vector<double> values){
    
    double median = 0.0;
    size_t n = values.size() / 2;
    // find the middle element in O(n)
    std::nth_element(values.begin(), values.begin() + n, values.end());
    median = values[n];
    // if even number of things have to average for the ''middle''
    if (values.size() %2 == 0){
        std::nth_element(values.begin(), values.begin() + n-1, values.end());
        median = (median + values[n-1])/2.0;
    }
    return median;


}// end median

// Simple mean

double mean (std::vector<double>& values){

    double sum = std::accumulate( values.begin(), values.end(), 0.0);
    
    double average = (double) sum /  values.size();

    return average;
}//end mean

//divide two vectors (for HANorm)
std::vector<double> divideVector(std::vector<double> numerator,
    std::vector<double> denominator){
    
    std::vector<double> result;
    for (int i = 0; i < numerator.size(); i++){
        result.emplace_back(numerator[i]/denominator[i]);
    }
    return result;
}//end divide

double HealthyAging (std::vector<DeficitVal>& Values,  int N){
    
    double FINodes = (double)N; 
    double s = 0.0; //(1.0 - 2.0/(2.0*N))*(Values[0].age);
    std::cout << "Deficit Vals Size: " << Values.size() << "\n";
    for (int i = 1; i < Values.size(); i++){
        
        s += (1.0  + (1 - i + 0.5)/FINodes)
            * (Values[i].age - Values[i-1].age);

    }//end for loop over all deficit values

    return s;

}//end healthy aging


double HealthyAgingFIN (std::vector<DeficitVal>& Values, 
    std::vector<int>& FIIDs){
    
    std::vector<DeficitVal> FIValues;


    for (int i = 0; i < Values.size(); i++){
        if ( std::find( FIIDs.begin(), FIIDs.end(), Values[i].id ) !=
        FIIDs.end() ){
            FIValues.emplace_back(Values[i]);        
        }//end selection if statement
    }// end FI node selection loop

    int size = FIValues.size();
    double FINodes = FIValues.size(); 
    double s = (1 - 2/(2*FINodes))*(FIValues[0].age);
    
    for (int i = 1; i < FIValues.size(); i++){
        
        s += (1.0 - i/FINodes + 1.0/(2.0*FINodes))
            * (FIValues[i].age - FIValues[i-1].age);

    }//end for loop over FI deficit values
    return s;

}// end Healthy Aging with N Most Connected FI Nodes


std::vector<int> FINodes (std::vector<Node>& Network, const int nNodes,
    std::vector<int> mortNodes){
    
    //connectivities
    std::vector<int> ks(Network.size());
    for (int i = 0; i < Network.size(); i++){
        ks[i] = Network[i].k;
    }//end k finding for loop
    
    //sort the connectivities into ascending order
    std::sort( ks.begin(), ks.end());

    
    // minimum k for FI calculations is the nNodesth largest K
    int kMin = ks[ks.size()  - nNodes];
    
    //vector of FI IDs
    std::vector<int> FIIDs;
    //mark the Nodes that are minimum acceptable k, some may habe to be removed
    std::vector<int> tieIndices;
    //find the IDs
    for (int i = 0; i < ks.size(); i++) { 
        if (Network[i].k > kMin){
            FIIDs.emplace_back(Network[i].id);
        }
        if (Network[i].k == kMin){
            tieIndices.emplace_back(FIIDs.size());
            FIIDs.emplace_back(Network[i].id);
        }//end ties
    }// end fi node finding loop
    
    //above loop doesnt handle multiple ks at kmin, if that occurs this handles 
    int index = 0;
    while (FIIDs.size() > nNodes){
        index = (int)(RNG(g1)*FIIDs.size());
        if (mortNodes[0] != FIIDs[index] && mortNodes[1] != FIIDs[index]){
            FIIDs.erase(FIIDs.begin() + index);
        }
    }//end ties for k min
    //std::cout << FIIDs.size();

    
    return FIIDs;
}//end FINodes


//Calculate standard devition of a vector over sqrt(n)

double errorBar(std::vector<double> values){
    //size
    double N = (double) values.size();
    //sum
    double sum = std::accumulate(values.begin(), values.end(), 0.0);
    //mean
    double mean = sum/N;

    double accum = 0.0;
    std::for_each (values.begin(), values.end(), [&](const double d) {
        accum += (d - mean) * (d - mean);
    });
    
    double stdev = sqrt(accum / (N-1));
    
    double error = stdev/sqrt(N);
    
    return error;
}

//Calculate the Reproductive node proc timings
std::vector<double> calculateChildEvents(std::vector<DeficitVal>& Values,
    std::vector<int>& RNodes){
    
    std::vector<double> childTimes;
    childTimes.reserve(RNodes.size());

    for (int i = 0; i < Values.size(); i++){
        if ( std::find( RNodes.begin(), RNodes.end(), Values[i].id ) !=
        RNodes.end() ){
            childTimes.emplace_back(Values[i].age);        
        }//end selection if statement
    }// end time selection loop

    return childTimes;
}//end reproduction node timings

// Directionality Theory Evolutionary Entropy
double evolutionaryEntropy(std::vector<double>& deathAges,
    std::vector<double>& childEvents, double growthRate, double binWidth){
    
    double tot = (double) deathAges.size();

    std::sort( deathAges.begin(), deathAges.end());
    
    double tMax = deathAges[deathAges.size() - 1];
    int numBins = (int) (tMax / binWidth) + 1;

    gsl_histogram * deathDist = gsl_histogram_alloc(numBins);
    gsl_histogram_set_ranges_uniform(deathDist, 0, tMax + binWidth);

    // make a death age distribution
    for(int i = 0; i < tot; i++){
        gsl_histogram_increment(deathDist, deathAges[i]);
    }//end histogram filler

    // use death age distribution to calculate l(x)
    std::vector<double> l;
    l.reserve(numBins);
    l.emplace_back(1.0);
    for(int i = 1; i < numBins; i++){
        l.emplace_back(l[i-1] - (double)gsl_histogram_get(deathDist, i)/tot);
    }//end l calculation

    //calculate m(x) from child events
    gsl_histogram * childDist = gsl_histogram_alloc(numBins);
    gsl_histogram_set_ranges_uniform(childDist, 0, tMax);
    for(int i = 0; i < childEvents.size(); i++){
        gsl_histogram_increment(childDist, childEvents[i]);
    }//end m calculation

    //calculate the time array in this sucker as well
    std::vector<double> times;
    times.reserve(numBins);
    //calculate the product, V(x)
    std::vector<double> v;
    v.reserve(numBins);
    double time = 0.0;
    for(int i = 0; i < numBins; i++){
        time = (double)i*binWidth;
        times.emplace_back(time);
        v.emplace_back(l[i] * (double)gsl_histogram_get(childDist, i)
            * exp (-time*growthRate) );
    }//end V(x) calculation

    //normalize v(x), trapezoidal integration of the function
    double sum = 0.0;
    for(int i = 1; i < numBins - 1; i++){
        sum += binWidth*v[i];
    }//sum for trapezoidal integration
    //endpoints
    sum += binWidth*(v[0] + v[numBins-1])/2.0;

    for(int i = 0; i < numBins; i++){
        v[i] /= sum;
    }//end V normalization

    //calculate numerator for evolutionary entropy (int(plnp))
    double numerator = 0.0;
    for(int i = 1; i < numBins - 1; i++){
        numerator += binWidth*(xlnx(v[i]));
    }//sum for trapezoidal integration
    //endpoints
    numerator += binWidth*(xlnx(v[0]) + xlnx(v[numBins-1]))/2.0;
    
    //calculate denominator
    double denominator = 0.0;
    for(int i = 1; i < numBins - 1; i++){
        denominator+= binWidth*(times[i]*v[i]);
    }//sum for trapezoidal integration
    //endpoints
    denominator+= binWidth*(v[0]*times[0] + v[numBins-1]*times[numBins-1])/2.0;

    double evoEntropy = -numerator/denominator;
    
    //destroy all the vectors and distribution stuff
    times.clear();
    gsl_histogram_free(childDist);
    gsl_histogram_free(deathDist);

    return evoEntropy;
}//end evolutionary Entropy

// Function for adding intermediate data in case of plateaus
bool iterationMatch(int iteration, std::vector<double> &matches){

    for (auto j: matches){
        if (abs(pow(10, (int) log10 (iteration)) * j - iteration) <= 1e-12 ){
            return true;
        }
    }
    return false;

}//end iteration match

// Calculate Clustering Coefficients
std::vector<double> clusteringCoefficients (std::vector<Node> & Network){
    
    std::vector<double> clusteringCoefficients;
    clusteringCoefficients.reserve(Network.size());

    using namespace boost;
    
    typedef adjacency_list <vecS, vecS, undirectedS> Graph;
    typedef exterior_vertex_property<Graph, double> ClusteringProperty;
    typedef ClusteringProperty::container_type ClusteringContainer;
    typedef ClusteringProperty::map_type ClusteringMap;
    Graph G;

    for(const auto &v: Network) {

	for(auto x: v.Neighbours) {
	    
	    //if(G.has_edge(v.id, Network[x].id))
	    auto res = add_edge(v.id, Network[x].id, G);
	    //if (std::get<1>(res)==false) std::cout << "refused" << std::endl;
	}
	    
    }

    ClusteringContainer coefs(num_vertices(G));
    ClusteringMap cm(coefs, G);
    double cc = all_clustering_coefficients(G, cm);
    
    graph_traits<Graph>::vertex_iterator i, end;
    for( tie(i, end) = vertices(G); i!= end; ++i){
        clusteringCoefficients.emplace_back(get(cm, *i));
    }


    return clusteringCoefficients;
}//end Clustering Coefficients calculator

//Calculate shortest paths
std::vector<double> pathLengths(std::vector<Node> &Network){


    using namespace boost;
    
    typedef adjacency_list <vecS, vecS, undirectedS> Graph;
    typedef graph_traits<Graph>::vertex_descriptor Vertex;

    Graph G;

    for(const auto &v: Network) {

	for(auto x: v.Neighbours) {
	    
	    //if(G.has_edge(v.id, Network[x].id))
	    auto res = add_edge(v.id, Network[x].id, G);
	    //if (std::get<1>(res)==false) std::cout << "refused" << std::endl;
	}
	    
    }

    std::vector<Vertex> v(num_vertices(G));
    graph_traits<Graph>::vertices_size_type distances[Network.size()];
    std::fill_n(distances, Network.size(), 0);

    Vertex first = *(vertices(G).first);
    v[first] = first;

    breadth_first_search(G, first, visitor(
        make_bfs_visitor(
            record_distances(distances, on_tree_edge()))));
    
    std::vector<double> distanceVector(distances,
        distances + sizeof distances / sizeof distances[0]);
    

    return distanceVector;
}//end pathLength calculator
