#ifndef Parameters_LL_class
#define Parameters_LL_class
#include "../../tensor_type.h"

class Parameters_LL
{

public:
    int lx, ly, ns;
    int n_orbs;
    string ModelType;
    int TBC_mx, TBC_my;
    int TBC_cellsX, TBC_cellsY;
    Mat_2_doub OnSiteE;
    double lambda_RSOC;
    double BoundaryConnection;

    string File_Onsite_Energies, File_Hoppings, File_LongRange_Ints;

    Matrix<double> hopping_NN_X;
    Matrix<double> hopping_NN_Y;
    Matrix<double> hopping_NNN_PXPY;
    Matrix<double> hopping_NNN_PXMY;
    double hopping_intracell;


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
    double Onsite_U;

    void Initialize(string inputfile_);
    double matchstring(string file, string match);
    string matchstring2(string file, string match);

};

void Parameters_LL::Initialize(string inputfile_)
{


    cout << "____________________________________" << endl;
    cout << "Reading the inputfile for specific model: " << inputfile_ << endl;
    cout << "____________________________________" << endl;

    File_Hoppings=matchstring2(inputfile_,"File_Hopping_Connections");
    File_Onsite_Energies=matchstring2(inputfile_,"File_Onsite_Energies");
    File_LongRange_Ints=matchstring2(inputfile_,"File_NonLocal_Int_Connections");

    lx = int(matchstring(inputfile_, "Xsite"));
    ly = int(matchstring(inputfile_, "Ysite"));

    TBC_mx = int(matchstring(inputfile_, "TwistedBoundaryCond_mx"));
    n_orbs = int(matchstring(inputfile_, "N_Orbs"));

    OnSiteE.resize(n_orbs);
    for (int orb=0;orb<3;orb++){
        OnSiteE[orb].resize(2);
    }

    TBC_my = int(matchstring(inputfile_, "TwistedBoundaryCond_my"));
    TBC_cellsX = int(matchstring(inputfile_, "TBC_cellsX"));
    TBC_cellsY = int(matchstring(inputfile_, "TBC_cellsY"));

    BoundaryConnection = double(matchstring(inputfile_, "PBC"));

    ns = lx * ly;
    cout << "TotalNumberOf Unit cells = " << ns << endl;

    lambda_RSOC =matchstring(inputfile_, "lambda_RSOC");


    string OnSiteE_up_str = "OnSiteE_up";
    string OnSiteE_dn_str = "OnSiteE_dn";
    string OnSiteE_up_str_temp, OnSiteE_dn_str_temp ;

    for(int n=0;n<n_orbs;n++){
        OnSiteE_up_str_temp = OnSiteE_up_str + to_string(n);
        OnSiteE_dn_str_temp = OnSiteE_dn_str + to_string(n);
        OnSiteE[n][0] = matchstring(inputfile_, OnSiteE_up_str_temp);
        OnSiteE[n][1] = matchstring(inputfile_, OnSiteE_dn_str_temp);
    }

    //Hopping matrices -------------------
    hopping_NN_X.resize(n_orbs,n_orbs);
    hopping_NN_Y.resize(n_orbs,n_orbs);
    string Nearest_Neigh_Hopping_basis_X;
    string Nearest_Neigh_Hopping_basis_Y;

    string NN_X_str, NN_Y_str;
    for (int m=0;m<n_orbs;m++){

        NN_X_str = "Nearest_Neigh_Hopping_X_basis_row" + to_string(m);
        NN_Y_str = "Nearest_Neigh_Hopping_Y_basis_row" + to_string(m);
        Nearest_Neigh_Hopping_basis_X=matchstring2(inputfile_, NN_X_str);
        Nearest_Neigh_Hopping_basis_Y=matchstring2(inputfile_, NN_Y_str);

        stringstream stream_X(Nearest_Neigh_Hopping_basis_X);
        stringstream stream_Y(Nearest_Neigh_Hopping_basis_Y);

        for(int n=0;n<n_orbs;n++){
            stream_X >> hopping_NN_X(m,n);
            stream_Y >> hopping_NN_Y(m,n);

        }
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
        Onsite_U =matchstring(inputfile_, "Onsite_U");
        Total_Particles = matchstring(inputfile_,"Total_particles");
        cout << "TotalNumberOfParticles = "<< Total_Particles << endl;

        IterMax = int(matchstring(inputfile_,"No_of_SelfConsistency_iters"));
        Convergence_Error=matchstring(inputfile_,"Convergence_Error");
        RandomSeed = matchstring(inputfile_,"RandomSeed");
        alpha_OP = matchstring(inputfile_,"alpha_OP");
        Temperature = matchstring(inputfile_,"Temperature");
        beta=(1.0/Temperature);

       double Read_OPs_double=double(matchstring(inputfile_,"Read_initial_OPvalues"));
        if(Read_OPs_double==1.0){
            Read_OPs=true;
        }
        else{
            Read_OPs=false;
        }
    }
    //------------------------------------------------



    //Next Nearest hopping------------
    hopping_NNN_PXPY.resize(n_orbs,n_orbs);
    hopping_NNN_PXMY.resize(n_orbs,n_orbs);
    //If needed read from input file

    //Hopping matrices done---------------

    hopping_intracell=matchstring(inputfile_, "Hopping_intracell");
    assert(n_orbs==3);
    cout << "____________________________________" << endl;
}

double Parameters_LL::matchstring(string file, string match)
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

string Parameters_LL::matchstring2(string file, string match)
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
        cout << "Unable to open input file while in the Parameters_LL class." << endl;
    }

    cout << match << " = " << amount << endl;
    return amount;
}

#endif
