#ifndef SAVEDATA_H
#define SAVEDATA_H



void OutputDataFile::SaveData(std::vector<DeficitVal> &DeficitsValues) {

    std::ios::sync_with_stdio(false);
	
    for(const auto &x: DeficitsValues) {
	
        DeficitsOut << x.age << "," << x.id << "," << x.d << "\n";
    		    
    }

    DeficitsValues.clear();

};

void Output2d(const std::vector<std::vector<double>> &Data, int OriginalN,
    std::string dataType) {

    std::string name = TempFolder() + dataType + SetRawName(OriginalN);
    std::ofstream Output;
    Output.open(name.c_str());

    for(std::vector<double> year: Data){
        for (double x: year) Output << x << "\t";
        Output << "\n";
    }
    Output.close();
    
}


void OutputDeathAges(const std::vector<double> &DeathAges, int OriginalN) {

    std::string name = TempFolder() + "RawDeathAgeData" + SetRawName(OriginalN);
    std::ofstream Output;
    Output.open(name.c_str());

    for(double x: DeathAges) Output << x << "\n";

    Output.close();
    
}


void OutputMeans(const std::vector<double> &Means, int OriginalN,
    std::string OutputType) {

    std::string name = TempFolder() + OutputType + SetRawName(OriginalN);
    //std::cout << "Folder: " << TempFolder() << "\n";
    std::ofstream Output;
    //std::cout << "Outputting data to: " << name << "\n";
    Output.open(name.c_str());

    for(double x: Means) Output << x << "\n";

    Output.close();
    
}



void OutputInts(const std::vector<int> &Means, int OriginalN,
    std::string OutputType) {

    std::string name = TempFolder() + OutputType + SetRawName(OriginalN);
    std::ofstream Output;
    Output.open(name.c_str());

    for(double x: Means) Output << x << "\n";

    Output.close();
    
}

void OutputInts2d(const std::vector<std::vector<int>> data, int OriginalN,
                    std::string dataType){
    std::string name = TempFolder() + dataType + SetRawName(OriginalN);
    std::ofstream Output;
    Output.open(name.c_str());

    for(std::vector<int> row: data){
        for (int x: row) Output << x << "\t";
        Output << "\n";
    }
    Output.close();
 
}


#endif
