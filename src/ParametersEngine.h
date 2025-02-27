#ifndef Parameters_class
#define Parameters_class
#include "tensor_type.h"

class Parameters{

public:
    int ns, IterMax, RandomSeed;
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
    bool Anderson_Mixing;
    int AM_m; //Anderson_Mixing_m;
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

    bool Restricted_HF;
    string Ansatz;
    string Ansatz_file;
    Mat_1_int CDW_Ansatz_sites;
    double A_charge_modulation;

    void Initialize(string inputfile_);
    double matchstring(string file,string match);
    string matchstring2(string file,string match);

};


void Parameters::Initialize(string inputfile_){


    CDW_Ansatz_sites.clear();

    double Simple_Mixing_double, Broyden_Mixing_double, Anderson_Mixing_double, BroydenSecondMethodMixing_double, Restricted_HF_double;
    double Read_OPs_double;
    double Create_OPs_double;
    double Just_Hartree_double;

    cout << "____________________________________" << endl;
    cout << "Reading the inputfile name: " << inputfile_ << endl;
    cout << "____________________________________" << endl;



    ns = int(matchstring(inputfile_,"Total_sites"));
    BroydenSecondMethodCounter = int(matchstring(inputfile_,"BroydenSecondMethodCounter"));
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

    alpha_OP = matchstring(inputfile_,"alpha_OP");
    w_minus1 = matchstring(inputfile_,"w_minus1");
    wn = matchstring(inputfile_,"wn");

    Dflag = 'N';

    Simple_Mixing_double=double(matchstring(inputfile_,"Simple_Mixing"));
    Broyden_Mixing_double=double(matchstring(inputfile_,"Broyden_Mixing"));
    Anderson_Mixing_double=double(matchstring(inputfile_,"Anderson_Mixing"));
    BroydenSecondMethodMixing_double=double(matchstring(inputfile_,"Broyden_Second_Method_Mixing"));
    Restricted_HF_double=double(matchstring(inputfile_, "Restricted_HF"));
    if(Restricted_HF_double==1.0){
        Restricted_HF=true;
        Ansatz=matchstring2(inputfile_,"HF_Ansatz");

        if(!(Ansatz=="Given_CDW")){
            cout<<"Ansatz can be only in : Given_CDW"<<endl;
            assert(false);
        }

        if(Ansatz=="Given_CDW"){
            CDW_Ansatz_sites.resize(ns);
            for(int i=0;i<ns;i++){
                CDW_Ansatz_sites[i]=-1.0; //ideally empty
            }
            Ansatz_file=matchstring2(inputfile_,"Given_ansatz_file");
            ifstream inputfile_Ansatz_file(Ansatz_file.c_str());
            string line_temp_;
            int i_temp_;
            //getline(inputfile_Ansatz_file,line_temp_);
            while(getline(inputfile_Ansatz_file,line_temp_))
            {
                stringstream line_temp_ss_(line_temp_);
                line_temp_ss_ >> i_temp_;
                CDW_Ansatz_sites[i_temp_]=1.0; //ideally half-filled
            }
        }
    }
    else{
        Ansatz=matchstring2(inputfile_,"HF_Ansatz");
        Ansatz_file=matchstring2(inputfile_,"Given_ansatz_file");
        Restricted_HF=false;
    }


    Onsite_E.resize(ns);
    for(int i=0;i<ns;i++){
        Onsite_E[i].resize(2);
    }

    file_onsite_energies_=matchstring2(inputfile_,"File_Onsite_Energies");
    ifstream inputfile_Onsite_Energy(file_onsite_energies_.c_str());
    string line_temp;
    int i_temp, spin_temp;
    double onsite_temp;

    //cout<<"Here 1"<<endl;
    getline(inputfile_Onsite_Energy,line_temp); //#ix iy spin E[site]
    if(inputfile_Onsite_Energy.is_open())
    {
        while(getline(inputfile_Onsite_Energy,line_temp))
        {
            stringstream line_temp_ss(line_temp);
            cout<<line_temp<<endl;
            line_temp_ss >> i_temp >> spin_temp >> onsite_temp;
            Onsite_E[i_temp][spin_temp]=onsite_temp;
        }
        inputfile_Onsite_Energy.close();
    }
    else
    {cout<<"Unable to open file = '"<<file_onsite_energies_<<"'"<<endl;}




    //mapping b/w site and ix,iy MUST BE CONSISTENT WITH Numbering_lattice() of Coordinate class.
    // basis_ = ix + lx*(iy) + lx*ly*(spin)
    //spin_UP=0, spin_DN=1
    int site_i, site_j, spin_i, spin_j;
    complex<double> Hopp_temp;
    LongRangeHoppings.resize(ns*2);
    for(int basis_=0;basis_<ns*2;basis_++){
        LongRangeHoppings[basis_].resize(ns*2);
        for(int basis2=0;basis2<ns*2;basis2++ ){
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
            LongRangeHoppings[site_i+(ns*spin_i)][site_j+(ns*spin_j)]=Hopp_temp;
        }
        inputfile_hopping_connections.close();
    }
    else
    {cout<<"Unable to open file = '"<<file_hopping_connections_<<"'"<<endl;}




    //mapping b/w site and ix,iy MUST BE CONSISTENT WITH Numbering_lattice() of Coordinate class.
    // basis_ = ix + lx*(iy) + lx*ly*(spin)
    //spin_UP=0, spin_DN=1
    complex<double> Int_temp;
    LongRangeInteractions.resize(ns);
    for(int basis_=0;basis_<ns;basis_++){
        LongRangeInteractions[basis_].resize(ns);
        for(int basis2=0;basis2<ns;basis2++ ){
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

//        for(int i=0;i<LongRangeInteractions.size();i++){
//            for(int j=0;j<LongRangeInteractions.size();j++){
//                cout<<LongRangeInteractions[i][j]<<" ";
//            }
//            cout<<endl;
//        }


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


    BroydenSecondMethodMixing=false;
    Broyden_Mixing=false;
    Simple_Mixing=false;
    Anderson_Mixing=false;
    if(BroydenSecondMethodMixing_double==1.0){
        BroydenSecondMethodMixing=true;
    }
    if(Broyden_Mixing_double==1.0){
        Broyden_Mixing=true;
    }
    if(Simple_Mixing_double==1.0){
        Simple_Mixing=true;
    }
    if(Anderson_Mixing_double==1.0){
        Anderson_Mixing=true;
        AM_m = int(matchstring(inputfile_,"Anderson_Mixing_m"));
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
        cout<<"PERFORMING JUST HARTREE"<<endl;
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



