#ifndef Parameters_HC_class
#define Parameters_HC_class
#include "../../tensor_type.h"

class Parameters_HC
{

public:
    int lx, ly, ns; //ns=lx*ly total number of site, each site has 2 orbitals
    int n_orbs; //n_orbs=2 for Honeycomb
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

    complex<double> t0_upup, t0_dndn;
    Mat_2_Complex_doub t1_plus_a1_upup, t1_minus_a2_upup, t1_plus_a1_dndn, t1_minus_a2_dndn,  t2_upup, t2_dndn, t3;
    double U0;
    double U0_interorb;
    Mat_2_doub U1, U2, U3;
    double U0ByUNN;
    double AnisotropyZ;
    double hopping_intracell;

    int UnitCellSize_x, UnitCellSize_y;
    bool Just_Hartree;

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
    pair_int UnitCellType_intialOPs;
    bool Create_OPs_Ansatz;
    string OP_Ansatz_type;

    bool Fixing_mu;
    double Fixed_mu;

    bool Anderson_Mixing;
    int AM_m; //Anderson_Mixing_m;

    string File_OPs_out, File_OPs_in;
    string FockType;

    bool Cooling;
    Mat_1_doub Temperature_points;

    void Initialize(string inputfile_);
    double matchstring(string file, string match);
    string matchstring2(string file, string match);

    bool OP_only_finite_Int;
    double Truncating_Length_in_am;

};

void Parameters_HC::Initialize(string inputfile_)
{


    cout << "____________________________________" << endl;
    cout << "Reading the inputfile for specific model: " << inputfile_ << endl;
    cout << "____________________________________" << endl;

    File_Hoppings=matchstring2(inputfile_,"File_Hopping_Connections");
    File_Onsite_Energies=matchstring2(inputfile_,"File_Onsite_Energies");
    File_LongRange_Ints=matchstring2(inputfile_,"File_NonLocal_Int_Connections");


    n_orbs = int(matchstring(inputfile_, "N_Orbs"));
    assert(n_orbs==2);

    t1_plus_a1_upup.resize(n_orbs);
    t1_minus_a2_upup.resize(n_orbs);
    t1_plus_a1_dndn.resize(n_orbs);
    t1_minus_a2_dndn.resize(n_orbs);
    t2_upup.resize(n_orbs);
    t2_dndn.resize(n_orbs);
    t3.resize(n_orbs);
    for(int i=0;i<n_orbs;i++){
        t1_plus_a1_upup[i].resize(n_orbs);
        t1_minus_a2_upup[i].resize(n_orbs);
        t1_plus_a1_dndn[i].resize(n_orbs);
        t1_minus_a2_dndn[i].resize(n_orbs);
        t2_upup[i].resize(n_orbs);
        t2_dndn[i].resize(n_orbs);
        t3[i].resize(n_orbs);
    }

    string string_t0_upup = matchstring2(inputfile_, "t0_upup");
    string string_t0_dndn = matchstring2(inputfile_, "t0_dndn");

    string string_t1_plus_a1_upup = matchstring2(inputfile_, "t1_plus_a1_upup");
    string string_t1_minus_a2_upup = matchstring2(inputfile_, "t1_minus_a2_upup");
    string string_t1_plus_a1_dndn= matchstring2(inputfile_, "t1_plus_a1_dndn");
    string string_t1_minus_a2_dndn = matchstring2(inputfile_, "t1_minus_a2_dndn");

    string string_t2_upup = matchstring2(inputfile_, "t2_upup");
    string string_t2_dndn = matchstring2(inputfile_, "t2_dndn");
    string string_t3 = matchstring2(inputfile_, "t3");
    stringstream t0_upup_ss(string_t0_upup);
    stringstream t0_dndn_ss(string_t0_dndn);
    stringstream t1_plus_a1_upup_ss(string_t1_plus_a1_upup);
    stringstream t1_minus_a2_upup_ss(string_t1_minus_a2_upup);
    stringstream t1_plus_a1_dndn_ss(string_t1_plus_a1_dndn);
    stringstream t1_minus_a2_dndn_ss(string_t1_minus_a2_dndn);
    stringstream t2_upup_ss(string_t2_upup);
    stringstream t2_dndn_ss(string_t2_dndn);
    stringstream t3_ss(string_t3);

    t0_upup_ss >> t0_upup;
    t0_dndn_ss >> t0_dndn;

    for(int i=0;i<n_orbs;i++){
        for(int j=0;j<n_orbs;j++){
            t1_plus_a1_upup_ss >> t1_plus_a1_upup[i][j];
            t1_minus_a2_upup_ss >> t1_minus_a2_upup[i][j];
            t1_plus_a1_dndn_ss >> t1_plus_a1_dndn[i][j];
            t1_minus_a2_dndn_ss >> t1_minus_a2_dndn[i][j];
            t2_upup_ss >> t2_upup[i][j];
            t2_dndn_ss >> t2_dndn[i][j];
            t3_ss >> t3[i][j];
        }
    }


    cout<<"t0_upup = "<< t0_upup<<endl;
    cout<<endl;

    cout<<"t0_dndn = "<< t0_dndn<<endl;
    cout<<endl;


    cout<<"t1 plus a1 upup matrix:-----------------"<<endl;
    for(int i=0;i<n_orbs;i++){
        for(int j=0;j<n_orbs;j++){
            cout<< t1_plus_a1_upup[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;

    cout<<"t1 plus a1 dndn matrix:-----------------"<<endl;
    for(int i=0;i<n_orbs;i++){
        for(int j=0;j<n_orbs;j++){
            cout<< t1_plus_a1_dndn[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;



    cout<<"t1 minus a2 upup matrix:-----------------"<<endl;
    for(int i=0;i<n_orbs;i++){
        for(int j=0;j<n_orbs;j++){
            cout<< t1_minus_a2_upup[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;


    cout<<"t1 minus a2 dndn matrix:-----------------"<<endl;
    for(int i=0;i<n_orbs;i++){
        for(int j=0;j<n_orbs;j++){
            cout<< t1_minus_a2_dndn[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;


    cout<<"t2 matrix up:-----------------"<<endl;
    for(int i=0;i<n_orbs;i++){
        for(int j=0;j<n_orbs;j++){
            cout<< t2_upup[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;


    cout<<"t2 matrix dn:-----------------"<<endl;
    for(int i=0;i<n_orbs;i++){
        for(int j=0;j<n_orbs;j++){
            cout<< t2_dndn[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;

    cout<<"t3 matrix:-----------------"<<endl;
    for(int i=0;i<n_orbs;i++){
        for(int j=0;j<n_orbs;j++){
            cout<< t3[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;




    //    cout<<"t1 = "<<t1<<endl;
    //    cout<<"t2 = "<<t2<<endl;
    //    cout<<"t3 = "<<t3<<endl;

    LongRange_interaction = int(matchstring(inputfile_, "LongRange_interaction"));
    a_moire=double(matchstring(inputfile_,"a_moire_in_Angstorm"));
    d_screening=double(matchstring(inputfile_,"d_screening_in_Angstorm"));
    eps_DE=double(matchstring(inputfile_,"eps_DE"));
    //U1ByUNN=double(matchstring(inputfile_,"U0ByUNN"));

    U0 = double(matchstring(inputfile_, "U0"));
    U0_interorb = double(matchstring(inputfile_, "U0_interorb")); //This is U1 too (nearest neighbour)


    //    U1 = double(matchstring(inputfile_, "U1"));
    //    U2 = double(matchstring(inputfile_, "U2"));
    //    U3 = double(matchstring(inputfile_, "U3"));



    U1.resize(n_orbs);
    U2.resize(n_orbs);
    U3.resize(n_orbs);
    for(int i=0;i<n_orbs;i++){
        U1[i].resize(n_orbs);
        U2[i].resize(n_orbs);
        U3[i].resize(n_orbs);
    }

    string string_U1 = matchstring2(inputfile_, "U1");
    string string_U2 = matchstring2(inputfile_, "U2");
    string string_U3 = matchstring2(inputfile_, "U3");
    stringstream U1_ss(string_U1);
    stringstream U2_ss(string_U2);
    stringstream U3_ss(string_U3);

    for(int i=0;i<n_orbs;i++){
        for(int j=0;j<n_orbs;j++){
            U1_ss >> U1[i][j];
            U2_ss >> U2[i][j];
            U3_ss >> U3[i][j];
        }
    }

    cout<<"U1 matrix:-----------------"<<endl;
    for(int i=0;i<n_orbs;i++){
        for(int j=0;j<n_orbs;j++){
            cout<< U1[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;

    cout<<"U2 matrix:-----------------"<<endl;
    for(int i=0;i<n_orbs;i++){
        for(int j=0;j<n_orbs;j++){
            cout<< U2[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;

    cout<<"U3 matrix:-----------------"<<endl;
    for(int i=0;i<n_orbs;i++){
        for(int j=0;j<n_orbs;j++){
            cout<< U3[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;









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



    int Fixing_mu_int;
    Fixing_mu_int = int(matchstring(inputfile_, "Fixing_mu"));
    if(Fixing_mu_int==1){
        Fixing_mu=true;
        Fixed_mu = double(matchstring(inputfile_,"Mu_value"));
    }
    else{
        Fixing_mu=false;
    }




    TBC_mx = int(matchstring(inputfile_, "TwistedBoundaryCond_mx"));


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

            string string_OP_type = matchstring2(inputfile_, "Read_initial_OPvalues_UnitCellsType");
            stringstream OP_type_ss(string_OP_type);
            OP_type_ss >> UnitCellType_intialOPs.first;
            OP_type_ss >> UnitCellType_intialOPs.second;
        }
        else{
            Read_OPs=false;
        }
    }
    //------------------------------------------------



    double OP_only_finite_Int_doub = double(matchstring(inputfile_, "OP_only_finite_Int"));
    if(OP_only_finite_Int_doub==1.0){
        OP_only_finite_Int=true;
    }
    else{
        OP_only_finite_Int=false;
    }

    Truncating_Length_in_am = double(matchstring(inputfile_,"Truncating_Length_in_am"));


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


    cout << "____________________________________" << endl;
}

double Parameters_HC::matchstring(string file, string match)
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

string Parameters_HC::matchstring2(string file, string match)
{

    string line;
    ifstream readFile(file);
    string amount;
    int offset;

    match = match + "=";

    if (readFile.is_open())
    {
        while (!readFile.eof())
        {
            getline(readFile, line);

            if ((offset = line.find(match, 0)) != string::npos)
            {
                amount = line.substr(offset + match.length());
            }
        }
        readFile.close();
    }
    else
    {
        cout << "Unable to open input file while in the Parameters_HC class." << endl;
    }

    cout << match << " " << amount << endl;
    return amount;
}

#endif
