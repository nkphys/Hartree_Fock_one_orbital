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
    double BoundaryConnection;

    string File_Onsite_Energies, File_Hoppings, File_LongRange_Ints;

    double t1, t2, t3;
    double U0, U1, U2, U3;
    double hopping_intracell;


    //For observse
    double beta, mus;
    double eta, domega;

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

    t1 = double(matchstring(inputfile_, "t1"));
    t2 = double(matchstring(inputfile_, "t2"));
    t3 = double(matchstring(inputfile_, "t3"));

    U0 = double(matchstring(inputfile_, "U0"));
    U1 = double(matchstring(inputfile_, "U1"));
    U2 = double(matchstring(inputfile_, "U2"));
    U3 = double(matchstring(inputfile_, "U3"));


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
