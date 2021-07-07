/*
File to calculate the merits for optimization
*/

bool entropyLagrange(double lambda, std::vector<double> entropy,
    std::vector<double> measure, double &bestEntropy,
    double &bestMeasure, double beta){
    
    double avgEntropy = mean(entropy);
    double avgMeasure = mean(measure);

    double entropyDifference = (bestEntropy - avgEntropy)/bestEntropy;
    double measureDifference = (bestMeasure - avgMeasure)/bestMeasure;
    double loss = lambda*entropyDifference + (1-lambda)*measureDifference;
    double pAccept = exp(-beta*loss);
    //std::cout << "Loss: " << loss << ", Best Loss: " << bestLoss <<
    //    ", pAccept: " << pAccept << ".\n";
    double r = RNG(g1);
    if (r < pAccept){
        //std::cout << "Accepted\n";
        //std::cout << "Entropy Term: " << lambda*entropyDifference <<
        //    ", Measure term: " << (1-lambda)*measureDifference << "," << 
        //    " probability of being accepted: " << pAccept << "\n";;
        bestEntropy = avgEntropy;
        bestMeasure = avgMeasure;
        return true;
    }
    else {
        //std::cout << "Rejected\n";
        return false;
    }
}


// Cooling Schedule Function
double coolingSchedule(int iteration, double beta, double power){
    double scaledIteration = pow(double(iteration), power);
    double inverseTemperature = beta*log(scaledIteration);
    return inverseTemperature;
}
