#ifndef Parameters_SQL_class
#define Parameters_SQL_class
#include "../../tensor_type.h"

class Parameters_SQL
{

public:
    int lx, ly, ns;
    int n_orbs;
    string ModelType;
    int TBC_mx, TBC_my;
    int TBC_cellsX, TBC_cellsY;
    Mat_2_doub OnSiteE;
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

    void Initialize(string inputfile_);
    double matchstring(string file, string match);
    string matchstring2(string file, string match);

};

void Parameters_SQL::Initialize(string inputfile_)
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
    for (int orb=0;orb<n_orbs;orb++){
        OnSiteE[orb].resize(2);
    }

    TBC_my = int(matchstring(inputfile_, "TwistedBoundaryCond_my"));
    TBC_cellsX = int(matchstring(inputfile_, "TBC_cellsX"));
    TBC_cellsY = int(matchstring(inputfile_, "TBC_cellsY"));

    BoundaryConnection = double(matchstring(inputfile_, "PBC"));

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



    //Next Nearest hopping------------
    hopping_NNN_PXPY.resize(n_orbs,n_orbs);
    hopping_NNN_PXMY.resize(n_orbs,n_orbs);
    //If needed read from input file

    //Hopping matrices done---------------

    assert(n_orbs==1);
    cout << "____________________________________" << endl;
}

double Parameters_SQL::matchstring(string file, string match)
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

string Parameters_SQL::matchstring2(string file, string match)
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
        cout << "Unable to open input file while in the Parameters_SQL class." << endl;
    }

    cout << match << " = " << amount << endl;
    return amount;
}

#endif