/* Update total frailty of the system */ 

void updateFrailty(double &FI, Node x, int N){
    
    FI += double(2 * x.d - 1)/double(N);

}// end update Frailty

void updateHealthyAging(double &HA, double FI, double dt){
    
    HA -= dt * (1.0 - FI);

}
