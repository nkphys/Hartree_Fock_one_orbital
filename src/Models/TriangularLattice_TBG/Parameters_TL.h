#ifndef Parameters_TL_class
#define Parameters_TL_class
#include "../../tensor_type.h"

class Parameters_TL
{

public:
    int lx, ly, ns;
    int n_orbs;
    string ModelType;
    int TBC_mx, TBC_my;
    int TBC_cellsX, TBC_cellsY;
    Mat_2_doub OnSiteE;
    double BoundaryConnection_X, BoundaryConnection_Y;
    bool PBC_X, PBC_Y;

    string File_Onsite_Energies, File_Hoppings, File_LongRange_Ints;

    bool LongRange_interaction;
    double d_screening;
    double a_moire;
    double eps_DE;

    complex<double> t1, t2, t3;
    double U0, U1, U2, U3;
    double AnisotropyZ;
    double hopping_intracell;

    int UnitCellSize_x, UnitCellSize_y;
    bool Just_Hartree;

    bool Using_Initial_Ansatz;
    string Initial_Ansatz_type;

    //For observse
    double beta, mus;
    double eta, domega;

    bool Self_consistency_kspace;
    double Temperature;

    double Total_Particles;
    int IterMax;
    double Convergence_Error;
    int RandomSeed;
    double alpha_OP;
    bool Read_OPs;

    bool Anderson_Mixing;
    int AM_m; //Anderson_Mixing_m;

    string File_OPs_out, File_OPs_in;
    string FockType;

    bool Cooling;
    Mat_1_doub Temperature_points;

    void Initialize(string inputfile_);
    double matchstring(string file, string match);
    string matchstring2(string file, string match);

};

void Parameters_TL::Initialize(string inputfile_)
{


    cout << "____________________________________" << endl;
    cout << "Reading the inputfile for specific model: " << inputfile_ << endl;
    cout << "____________________________________" << endl;

    File_Hoppings=matchstring2(inputfile_,"File_Hopping_Connections");
    File_Onsite_Energies=matchstring2(inputfile_,"File_Onsite_Energies");
    File_LongRange_Ints=matchstring2(inputfile_,"File_NonLocal_Int_Connections");

    string string_t1 = matchstring2(inputfile_, "t1");
    string string_t2 = matchstring2(inputfile_, "t2");
    string string_t3 = matchstring2(inputfile_, "t3");
    stringstream t1_ss(string_t1);t1_ss >> t1;
    stringstream t2_ss(string_t2);t2_ss >> t2;
    stringstream t3_ss(string_t3);t3_ss >> t3;

    cout<<"t1 = "<<t1<<endl;
    cout<<"t2 = "<<t2<<endl;
    cout<<"t3 = "<<t3<<endl;

    LongRange_interaction = int(matchstring(inputfile_, "LongRange_interaction"));
    a_moire=double(matchstring(inputfile_,"a_moire_in_Angstorm"));
    d_screening=double(matchstring(inputfile_,"d_screening_in_Angstorm"));
    eps_DE=double(matchstring(inputfile_,"eps_DE"));
    U0 = double(matchstring(inputfile_, "U0"));
    U1 = double(matchstring(inputfile_, "U1"));
    U2 = double(matchstring(inputfile_, "U2"));
    U3 = double(matchstring(inputfile_, "U3"));
    AnisotropyZ = double(matchstring(inputfile_, "AnisotropyZ"));


    lx = int(matchstring(inputfile_, "Xsite"));
    ly = int(matchstring(inputfile_, "Ysite"));

    UnitCellSize_x =int(matchstring(inputfile_, "UnitCellSize_x"));
    UnitCellSize_y =int(matchstring(inputfile_, "UnitCellSize_y"));

    int Just_Hartree_int;
    Just_Hartree_int = int(matchstring(inputfile_, "Just_Hartree"));
    if(Just_Hartree_int==1){
        Just_Hartree=true;
    }
    else{
        Just_Hartree=false;
        FockType = matchstring2(inputfile_,"Fock_type");
    }

    int Using_Initial_Ansatz_int;
    Using_Initial_Ansatz_int = int(matchstring(inputfile_, "Using_Initial_Ansatz"));
    if(Using_Initial_Ansatz_int==1){
       Using_Initial_Ansatz=true;
       Initial_Ansatz_type=matchstring2(inputfile_,"Initial_Ansatz_type");
    }
    else{
       Using_Initial_Ansatz=false;
    }


    TBC_mx = int(matchstring(inputfile_, "TwistedBoundaryCond_mx"));
    n_orbs = int(matchstring(inputfile_, "N_Orbs"));

    OnSiteE.resize(n_orbs);
    for (int orb=0;orb<n_orbs;orb++){
        OnSiteE[orb].resize(2);
    }

    TBC_my = int(matchstring(inputfile_, "TwistedBoundaryCond_my"));
    TBC_cellsX = int(matchstring(inputfile_, "TBC_cellsX"));
    TBC_cellsY = int(matchstring(inputfile_, "TBC_cellsY"));

    BoundaryConnection_X = double(matchstring(inputfile_, "PBC_X"));
    BoundaryConnection_Y = double(matchstring(inputfile_, "PBC_Y"));
    if(BoundaryConnection_X==1.0){
        PBC_X=true;
    }
    else{
        PBC_X=false;
    }
    if(BoundaryConnection_Y==1.0){
        PBC_Y=true;
    }
    else{
        PBC_Y=false;
    }


    ns = lx * ly;
    cout << "TotalNumberOf Unit cells = " << ns << endl;

    string OnSiteE_up_str = "OnSiteE_up";
    string OnSiteE_dn_str = "OnSiteE_dn";
    string OnSiteE_up_str_temp, OnSiteE_dn_str_temp ;

    for(int n=0;n<n_orbs;n++){
        OnSiteE_up_str_temp = OnSiteE_up_str + to_string(n);
        OnSiteE_dn_str_temp = OnSiteE_dn_str + to_string(n);
        OnSiteE[n][0] = matchstring(inputfile_, OnSiteE_up_str_temp);
        OnSiteE[n][1] = matchstring(inputfile_, OnSiteE_dn_str_temp);
    }



    double Self_consistency_kspace_double;
    Self_consistency_kspace_double=double(matchstring(inputfile_,"Self_consistency_kspace"));
    if(Self_consistency_kspace_double==1){
        Self_consistency_kspace=true;}
    else{
        Self_consistency_kspace=false;
    }

    //New for self-consistency----------------------------

    if(Self_consistency_kspace){
        Total_Particles = matchstring(inputfile_,"Total_particles");
        cout << "TotalNumberOfParticles = "<< Total_Particles << endl;

        IterMax = int(matchstring(inputfile_,"No_of_SelfConsistency_iters"));
        Convergence_Error=matchstring(inputfile_,"Convergence_Error");
        RandomSeed = matchstring(inputfile_,"RandomSeed");
        alpha_OP = matchstring(inputfile_,"alpha_OP");


        double Read_OPs_double=double(matchstring(inputfile_,"Read_initial_OPvalues"));
        if(Read_OPs_double==1.0){
            Read_OPs=true;
        }
        else{
            Read_OPs=false;
        }
    }
    //------------------------------------------------



    double Anderson_Mixing_double=double(matchstring(inputfile_,"Anderson_Mixing"));
    if(Anderson_Mixing_double==1.0){
        Anderson_Mixing=true;
        AM_m = int(matchstring(inputfile_,"Anderson_Mixing_m"));
    }

    double Cooling_double=double(matchstring(inputfile_,"Cooling"));
    if(Cooling_double==1.0){
        Cooling=true;
        string temp_string =  matchstring2(inputfile_,"Temperature");
        stringstream temp_ss(temp_string);
        int T_nos;
        temp_ss>>T_nos;
        Temperature_points.resize(T_nos);
        for(int n=0;n<T_nos;n++){
            temp_ss>>Temperature_points[n];
        }
    }
    else{
        Cooling=false;
        Temperature = matchstring(inputfile_,"Temperature");
        Temperature_points.resize(1);
        Temperature_points[0]=Temperature;
        beta=(1.0/Temperature);
    }

    File_OPs_in=matchstring2(inputfile_,"Read_initial_OPvalues_file");
    File_OPs_out=matchstring2(inputfile_,"Write_Final_OPvalues_file");



    assert(n_orbs==1);
    cout << "____________________________________" << endl;
}

double Parameters_TL::matchstring(string file, string match)
{
    string test;
    string line;
    ifstream readFile(file);
    double amount;
    bool pass = false;
    while (std::getline(readFile, line))
    {
        std::istringstream iss(line);
        if (std::getline(iss, test, '=') && pass == false)
        {
            // ---------------------------------
            if (iss >> amount && test == match)
            {
                // cout << amount << endl;
                pass = true;
            }
            else
            {
                pass = false;
            }
            // ---------------------------------
            if (pass)
                break;
        }
    }
    if (pass == false)
    {
        string errorout = match;
        errorout += "= argument is missing in the input file!";
        throw std::invalid_argument(errorout);
    }
    cout << match << " = " << amount << endl;
    return amount;
}

string Parameters_TL::matchstring2(string file, string match)
{

    string line;
    ifstream readFile(file);
    string amount;
    int offset;

    if (readFile.is_open())
    {
        while (!readFile.eof())
        {
            getline(readFile, line);

            if ((offset = line.find(match, 0)) != string::npos)
            {
                amount = line.substr(offset + match.length() + 1);
            }
        }
        readFile.close();
    }
    else
    {
        cout << "Unable to open input file while in the Parameters_TL class." << endl;
    }

    cout << match << " = " << amount << endl;
    return amount;
}

#endif
