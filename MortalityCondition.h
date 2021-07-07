/*
File to deal with the new mortality Condition
*/

// Normalizing factor
double normalizingFactor(std::vector<Node> Network, double mu){
    
    double sum = 0.0;
    for (auto x: Network){
        double metric = double(x.k);
        sum += pow(metric, mu);
    }// end loop over network

    return 1.0/sum;
}// end normalizing prefactor

// Update mortality rate
void updateMortalityRate(double &mortalityRate, double &totalRate, Node x,
    double gammaD, double mu, double prefactor){

        double sign = double(x.d) * 2.0 - 1.0;
        
        double hazard = exp(pow(double(x.k), mu) * sign * gammaD * prefactor);
        /*
        if (x.d == 1){
            std::cout << "Mu: " << mu << "\n";
            std::cout << "X.k: " << x.k << "\n";
            std::cout << "x.d: " << x.d << "\n";
            std::cout << "sign: " << sign << "\n";
            std::cout << "gammaD: " << gammaD << "\n";
            std::cout << "prefactor: " << prefactor << "\n";
            std::cout << "mortalityRate: " << mortalityRate << "\n";
            std::cout << "totalRate: " << totalRate << "\n";
        }
        */
        totalRate -= mortalityRate;
        mortalityRate *= hazard;
        totalRate += mortalityRate;
}// end update to mortality rate

// Evaluation of mortality
void evaluateMortality(bool &dead, double totalRate, double mortalityRate){

    dead = bool((RNG(g1) * totalRate) < mortalityRate);

}// end evaluate mortality
