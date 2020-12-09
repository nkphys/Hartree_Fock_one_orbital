#ifndef Parameters_class
#define Parameters_class
#include "tensor_type.h"

class Parameters{

public:
    int lx, ly, ns, IterMax, RandomSeed;
    double Convergence_Error;
    int lx_cluster, ly_cluster;
    double mus,Total_Particles,pi, mu_old;
    double MuValueFixed;
    bool FixingMu;
    double U_onsite;
    double Disorder_Strength, RandomDisorderSeed;
    Mat_2_doub Onsite_E;
    Mat_2_Complex_doub LongRangeHoppings;
    Mat_2_Complex_doub LongRangeInteractions;

    double dw_dos, eta_dos;
    double w_min, w_max;

    bool Read_OPs;
    bool Create_OPs;
    bool Just_Hartree;
    string File_OPs_in, File_OPs_out;
    string Create_OPs_type;

    bool Simple_Mixing;
    bool Broyden_Mixing;
    bool BroydenSecondMethodMixing;
    double w_minus1,wn;
    int BroydenSecondMethodCounter;
    double alpha_OP;

    double Temperature,beta,Eav,maxmoment;

    bool ReadDisorder;
    string ReadDisorderString, FixingMuString;
    string DisorderSeedFile;

    string file_onsite_energies_;
    string file_hopping_connections_;
    string file_nonlocal_int_connections_;
    string read_onsite_energies;

    char Dflag;

    void Initialize(string inputfile_);
    double matchstring(string file,string match);
    string matchstring2(string file,string match);

};


void Parameters::Initialize(string inputfile_){



    double Simple_Mixing_double, Broyden_Mixing_double, BroydenSecondMethodMixing_double;
    double Read_OPs_double;
    double Create_OPs_double;
    double Just_Hartree_double;

    cout << "____________________________________" << endl;
    cout << "Reading the inputfile: " << inputfile_ << endl;
    cout << "____________________________________" << endl;



    lx = int(matchstring(inputfile_,"Xsite"));
    ly = int(matchstring(inputfile_,"Ysite"));
    BroydenSecondMethodCounter = int(matchstring(inputfile_,"BroydenSecondMethodCounter"));

    ns = lx*ly;
    cout << "TotalNumberOfSites = "<< ns << endl;

    Total_Particles = matchstring(inputfile_,"Total No. of particles");
    cout << "TotalNumberOfParticles = "<< Total_Particles << endl;

    IterMax = int(matchstring(inputfile_,"No_of_SelfConsistency_iters"));
    Convergence_Error=matchstring(inputfile_,"Convergence_Error");
    RandomSeed = matchstring(inputfile_,"RandomSeed");
    FixingMuString = matchstring2(inputfile_, "FixingMu");
    RandomDisorderSeed = matchstring(inputfile_,"RandomDisorderSeed");
    ReadDisorderString = matchstring2(inputfile_,"ReadDisorderConf");
    Disorder_Strength = matchstring(inputfile_,"Disorder_Strength");
    U_onsite = matchstring(inputfile_,"U_Onsite");
    Temperature = matchstring(inputfile_,"Temperature");

    //dw_dos, eta_dos
    dw_dos = matchstring(inputfile_, "dw_dos");
    eta_dos = matchstring(inputfile_, "eta_dos");
    w_min = matchstring(inputfile_, "w_min");
    w_max = matchstring(inputfile_, "w_max");

    alpha_OP = matchstring(inputfile_,"alpha_OP");
    w_minus1 = matchstring(inputfile_,"w_minus1");
    wn = matchstring(inputfile_,"wn");

    Dflag = 'N';

    Simple_Mixing_double=double(matchstring(inputfile_,"Simple_Mixing"));
    Broyden_Mixing_double=double(matchstring(inputfile_,"Broyden_Mixing"));
    BroydenSecondMethodMixing_double=double(matchstring(inputfile_,"Broyden_Second_Method_Mixing"));



    Onsite_E.resize(lx);
    for(int ix=0;ix<lx;ix++){
        Onsite_E[ix].resize(ly);
    }

    file_onsite_energies_=matchstring2(inputfile_,"File_Onsite_Energies");
    ifstream inputfile_Onsite_Energy(file_onsite_energies_.c_str());
    string line_temp;
    int ix_temp, iy_temp;
    double onsite_temp;
    getline(inputfile_Onsite_Energy,line_temp); //#ix iy E[site]
    for(int iy=0;iy<ly;iy++){
        for(int ix=0;ix<lx;ix++){
           inputfile_Onsite_Energy >> ix_temp >> iy_temp >>onsite_temp;
           assert(iy==iy_temp);
           assert(ix==ix_temp);
           Onsite_E[ix][iy]=onsite_temp;
        }
    }


    //mapping b/w site and ix,iy MUST BE CONSISTENT WITH Numbering_lattice() of Coordinate class.
    // basis_ = ix + lx*(iy) + lx*ly*(spin)
    //spin_UP=0, spin_DN=1
    int site_i, site_j, spin_i, spin_j;
    complex<double> Hopp_temp;
    LongRangeHoppings.resize(lx*ly*2);
    for(int basis_=0;basis_<lx*ly*2;basis_++){
        LongRangeHoppings[basis_].resize(lx*ly*2);
        for(int basis2=0;basis2<lx*ly*2;basis2++ ){
          LongRangeHoppings[basis_][basis2]=complex<double> (0.0, 0.0);
        }
    }
    file_hopping_connections_=matchstring2(inputfile_,"File_Hopping_Connections");
    ifstream inputfile_hopping_connections(file_hopping_connections_.c_str());
    getline(inputfile_hopping_connections,line_temp); //#site_i spin_i site_j spin_j Hopping[site_i,spin_i][site_j,spin_j]

    if(inputfile_hopping_connections.is_open())
    {
        while(getline(inputfile_hopping_connections,line_temp))
        {
        stringstream line_temp_ss(line_temp);
        line_temp_ss >> site_i >> spin_i >> site_j >> spin_j >> Hopp_temp;
//        cout<< site_i <<"  "<< spin_i <<"  "<< site_j <<"  "<< spin_j <<"  "<< Hopp_temp<<endl;
        LongRangeHoppings[site_i+(lx*ly*spin_i)][site_j+(lx*ly*spin_j)]=Hopp_temp;
        }
        inputfile_hopping_connections.close();
    }
    else
    {cout<<"Unable to open file = '"<<file_hopping_connections_<<"'"<<endl;}




    //mapping b/w site and ix,iy MUST BE CONSISTENT WITH Numbering_lattice() of Coordinate class.
    // basis_ = ix + lx*(iy) + lx*ly*(spin)
    //spin_UP=0, spin_DN=1
    complex<double> Int_temp;
    LongRangeInteractions.resize(lx*ly);
    for(int basis_=0;basis_<lx*ly;basis_++){
        LongRangeInteractions[basis_].resize(lx*ly);
        for(int basis2=0;basis2<lx*ly;basis2++ ){
          LongRangeInteractions[basis_][basis2]=complex<double> (0.0, 0.0);
        }
    }

    file_nonlocal_int_connections_=matchstring2(inputfile_,"File_NonLocal_Int_Connections");
    ifstream inputfile_nonlocal_int_connections(file_nonlocal_int_connections_.c_str());
    getline(inputfile_nonlocal_int_connections,line_temp); //#site_i site_j Interaction[site_i][site_j]

    if(inputfile_nonlocal_int_connections.is_open())
    {
        while(getline(inputfile_nonlocal_int_connections,line_temp))
        {
        stringstream line_temp_ss2(line_temp);
        line_temp_ss2 >> site_i >> site_j >> Int_temp;
        LongRangeInteractions[site_i][site_j]=Int_temp;
        }
        inputfile_nonlocal_int_connections.close();
    }
    else
    {cout<<"Unable to open file = '"<<file_nonlocal_int_connections_<<"'"<<endl;}




    if(FixingMuString=="true"){
        FixingMu=true;
        MuValueFixed = matchstring(inputfile_,"MuValueFixed");
    }
    else{
        FixingMu=false;
        MuValueFixed=-100000;
    }

    if(ReadDisorderString=="true"){
        ReadDisorder = true;
        DisorderSeedFile = matchstring2(inputfile_,"DisorderSeedFile");
    }
    else{
        ReadDisorder = false;
    }

    if(BroydenSecondMethodMixing_double==1.0){
         BroydenSecondMethodMixing=true;
        Broyden_Mixing=false;
        Simple_Mixing=false;
    }
    else{
         BroydenSecondMethodMixing=false;
        if(Broyden_Mixing_double==1.0){
            Broyden_Mixing=true;
            Simple_Mixing=false;

        }
        else if(Broyden_Mixing_double==0.0){
            Broyden_Mixing=false;
            Simple_Mixing=true;
            cout<<"Broyden_Mixing and  BroydenSecondMethodMixing, both are 0(false). So Simple mixing is used"<<endl;

        }

    }



    Read_OPs_double=double(matchstring(inputfile_,"Read_initial_OPvalues"));
    if(Read_OPs_double==1.0){
        Read_OPs=true;
    }
    else{
        Read_OPs=false;
    }


    Just_Hartree_double=double(matchstring(inputfile_,"Just_Hartree"));
    if(Just_Hartree_double==1.0){
        Just_Hartree=true;
    }
    else{
        Just_Hartree=false;
    }

    Create_OPs_double=double(matchstring(inputfile_,"Create_OPvalues"));
    if(Create_OPs_double==1.0){
        assert(!Read_OPs);
        Create_OPs=true;
    Create_OPs_type=matchstring2(inputfile_,"Create_OPType");
    }
    else{
        Create_OPs=false;
    }


    File_OPs_in=matchstring2(inputfile_,"Read_initial_OPvalues_file");
    File_OPs_out=matchstring2(inputfile_,"Write_Final_OPvalues_file");

    pi=4.00*atan(double(1.0));
    Eav=0.0;


    //beta=(11605.0/Temperature);
     beta=(1.0/Temperature);

    mus=0.25;
    mu_old=0.25;
    cout << "____________________________________" << endl;
}


double Parameters::matchstring(string file,string match) {
    string test;
    string line;
    ifstream readFile(file);
    double amount;
    bool pass=false;
    while (std::getline(readFile, line)) {
        std::istringstream iss(line);
        if (std::getline(iss, test, '=') && pass==false) {
            // ---------------------------------
            if (iss >> amount && test==match) {
                // cout << amount << endl;
                pass=true;
            }
            else {
                pass=false;
            }
            // ---------------------------------
            if(pass) break;
        }
    }
    if (pass==false) {
        string errorout=match;
        errorout+="= argument is missing in the input file!";
        throw std::invalid_argument(errorout);
    }
    cout << match << " = " << amount << endl;
    return amount;
}

string Parameters::matchstring2(string file,string match) {

    string line;
    ifstream readFile(file);
    string amount;
    int offset;

    if(readFile.is_open())
    {
        while(!readFile.eof())
        {
            getline(readFile,line);

            if ((offset = line.find(match, 0)) != string::npos) {
                amount = line.substr (offset+match.length()+1);				}

        }
        readFile.close();
    }
    else
    {cout<<"Unable to open input file while in the Parameters class."<<endl;}




    cout << match << " = " << amount << endl;
    return amount;
}

#endif



