#ifndef Parameters_G2dLatticeNew_class
#define Parameters_G2dLatticeNew_class
#include "../../tensor_type.h"
#include "../../Matrix.h"
#include <fstream>

class Parameters_G2dLatticeNew
{

public:
    int lx, ly, ns; //ns=lx*ly total number of site, each site has n_orbs orbitals and n_atoms unit cell
    int n_orbs; //n_orbs is number of orbitals per atom 
    int n_atoms; //n_atoms in a unit cell

    Mat_2_doub OnSiteE;
    Mat_1_doub PinningFieldZ;
    Mat_1_doub PinningFieldX;


    double BoundaryConnection_X, BoundaryConnection_Y;
    bool PBC_X, PBC_Y;

    string File_Onsite_Energies, File_Hoppings, File_LongRange_Ints;

    bool LongRange_interaction;
    double d_screening;
    double a_moire;
    double eps_DE;

    Mat_2_Complex_doub t0; //index = atom  + orb*2 + spin*4
    Mat_2_Complex_doub t1_plus_a1, t1_minus_a2, t1_plus_a1_minus_a2,t1_plus_a1_plus_a2;


    Mat_3_Complex_doub Susc_OprA, Susc_OprB;
    bool Calculate_Susc;

    Mat_1_doub U0;
    Mat_1_doub JHund;
    Mat_1_doub UPrime;
    Mat_3_doub UInterOrb;


    string File_onsite_U, File_pa1_ma2_U, File_ma2_U;
    double U0ByUNN;
    double AnisotropyZ;
    double hopping_intracell;

    Matrix<int> MUC_Mat;
    bool Just_Hartree;

    bool NoSpinFlipOP;
    //For observse
    double beta, mus;

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

    double PinningFieldValue;
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

    int UnitCellSize_x, UnitCellSize_y;


    double omega_min;
    double omega_max;
    double d_omega;
    double eta;

    int NProcessors;

};

void Parameters_G2dLatticeNew::Initialize(string inputfile_)
{
    cout << "____________________________________" << endl;
    cout << "Reading the inputfile for specific model: " << inputfile_ << endl;
    cout << "____________________________________" << endl;


    MUC_Mat.resize(2,2);


    NProcessors =int(matchstring(inputfile_, "N_Threads"));
    n_atoms = int(matchstring(inputfile_, "N_Atoms"));
    n_orbs = int(matchstring(inputfile_, "N_Orbs"));
    //assert(n_orbs==2);

    t0.resize(n_orbs*n_atoms*2);
    t1_plus_a1.resize(n_orbs*n_atoms*2);
    t1_minus_a2.resize(n_orbs*n_atoms*2);
    t1_plus_a1_minus_a2.resize(n_orbs*n_atoms*2);
    t1_plus_a1_plus_a2.resize(n_orbs*n_atoms*2);
    for(int i=0;i<n_orbs*n_atoms*2;i++){
        t0[i].resize(n_orbs*n_atoms*2);
        t1_plus_a1[i].resize(n_orbs*n_atoms*2);
        t1_minus_a2[i].resize(n_orbs*n_atoms*2);
        t1_plus_a1_minus_a2[i].resize(n_orbs*n_atoms*2);
        t1_plus_a1_plus_a2[i].resize(n_orbs*n_atoms*2);
    }

    File_Hoppings = matchstring2(inputfile_, "File_Hoppings");
    string string_t0 = matchstring2(inputfile_, "t0_mat_file");
    string string_t1_plus_a1 = matchstring2(inputfile_, "t1_plus_a1_file");
    string string_t1_minus_a2 = matchstring2(inputfile_, "t1_minus_a2_file");
    string string_t1_plus_a1_minus_a2 = matchstring2(inputfile_, "t1_plus_a1_minus_a2_mat_file");
    string string_t1_plus_a1_plus_a2 = matchstring2(inputfile_, "t1_plus_a1_plus_a2_mat_file");

    ifstream file_t0(string_t0.c_str());
    ifstream file_t1_plus_a1(string_t1_plus_a1.c_str());
    ifstream file_t1_minus_a2(string_t1_minus_a2.c_str());
    ifstream file_t1_plus_a1_minus_a2(string_t1_plus_a1_minus_a2.c_str());
    ifstream file_t1_plus_a1_plus_a2(string_t1_plus_a1_plus_a2.c_str());


    for(int i=0;i<2*n_atoms*n_orbs;i++){
        for(int j=0;j<2*n_atoms*n_orbs;j++){
            file_t0>>t0[i][j];
            //t0[i][j].imag(0);
            file_t1_plus_a1_minus_a2>>t1_plus_a1_minus_a2[i][j];
            t1_plus_a1_minus_a2[i][j].imag(0);
            file_t1_plus_a1_plus_a2>>t1_plus_a1_plus_a2[i][j];
            t1_plus_a1_plus_a2[i][j].imag(0);
            file_t1_plus_a1>>t1_plus_a1[i][j];
            t1_plus_a1[i][j].imag(0);
            file_t1_minus_a2>>t1_minus_a2[i][j];
            t1_minus_a2[i][j].imag(0);
        }
    }


    string Calculate_Susc_str = matchstring2(inputfile_,"CalculateSusceptibility");
    if(Calculate_Susc_str=="true"){
        Calculate_Susc=true;
    }
    else{
        Calculate_Susc=false;
    }






    if(Calculate_Susc){
    string string_Susc_OprA, string_Susc_OprB;
    string string_Susc_Oprs = matchstring2(inputfile_, "SusceptibilityLocalOprFilesPairWise");
    stringstream Susc_Oprs_stream(string_Susc_Oprs);

    int No_of_ABpairs;
    Susc_Oprs_stream>>No_of_ABpairs;
    Susc_OprA.resize(No_of_ABpairs);
    Susc_OprB.resize(No_of_ABpairs);

    for(int pair_no=0;pair_no<No_of_ABpairs;pair_no++){
    Susc_OprA[pair_no].resize(n_orbs*n_atoms*2);
    Susc_OprB[pair_no].resize(n_orbs*n_atoms*2);
    for(int i=0;i<n_orbs*n_atoms*2;i++){
        Susc_OprA[pair_no][i].resize(n_orbs*n_atoms*2);
        Susc_OprB[pair_no][i].resize(n_orbs*n_atoms*2);
        for(int j=0;j<n_orbs*n_atoms*2;j++){
            Susc_OprA[pair_no][i][j]=0.0;
            Susc_OprB[pair_no][i][j]=0.0;
        }
    }
    }

    for(int pair_no=0;pair_no<No_of_ABpairs;pair_no++){
    Susc_Oprs_stream>>string_Susc_OprA>>string_Susc_OprB;
    ifstream file_Susc_OprA(string_Susc_OprA.c_str());
    ifstream file_Susc_OprB(string_Susc_OprB.c_str());

    for(int i=0;i<2*n_atoms*n_orbs;i++){
        for(int j=0;j<2*n_atoms*n_orbs;j++){
            file_Susc_OprA>>Susc_OprA[pair_no][i][j];
            file_Susc_OprB>>Susc_OprB[pair_no][i][j];
        }
    }

    }

    }



    U0.resize(n_atoms);
    string U0_string = matchstring2(inputfile_, "U0");
    stringstream U0_stream(U0_string);
    for(int atom_no=0;atom_no<n_atoms;atom_no++){
    U0_stream>>U0[atom_no];
        cout<<"U["<<atom_no<<"] = "<<U0[atom_no]<<endl;
    }

    JHund.resize(n_atoms);
    string JHund_string = matchstring2(inputfile_, "JHund");
    stringstream JHund_stream(JHund_string);
    for(int atom_no=0;atom_no<n_atoms;atom_no++){
        JHund_stream>>JHund[atom_no];
    }

    UPrime.resize(n_atoms);
    for(int atom_no=0;atom_no<n_atoms;atom_no++){
        UPrime[atom_no] = U0[atom_no] - 2.0*JHund[atom_no];
    }



    string U_interorb_someatom_str_temp = "U_interorb_atom_";
    string U_interorb_atom_str;
    UInterOrb.resize(n_atoms);
    for(int atom_no=0;atom_no<n_atoms;atom_no++){
        UInterOrb[atom_no].resize(n_orbs);
        for(int orb_no_i=0;orb_no_i<n_orbs;orb_no_i++){
            UInterOrb[atom_no][orb_no_i].resize(n_orbs);
        }
    }

    for(int atom_no=0;atom_no<n_atoms;atom_no++){

        U_interorb_atom_str = U_interorb_someatom_str_temp + to_string(atom_no);
        string UInterOrb_string = matchstring2(inputfile_, U_interorb_atom_str);
        stringstream UInterOrb_stream(UInterOrb_string);

        for(int orb_no_i=0;orb_no_i<n_orbs;orb_no_i++){
            for(int orb_no_j=0;orb_no_j<n_orbs;orb_no_j++){
                UInterOrb_stream>>UInterOrb[atom_no][orb_no_i][orb_no_j];

            }
        }
    }


    lx = int(matchstring(inputfile_, "Xsite"));
    ly = int(matchstring(inputfile_, "Ysite"));

    string MUC_Mat_row0_string =matchstring2(inputfile_, "MUCMat_row0");
    string MUC_Mat_row1_string =matchstring2(inputfile_, "MUCMat_row1");
    stringstream MUC_Mat_row0_stream(MUC_Mat_row0_string);
    stringstream MUC_Mat_row1_stream(MUC_Mat_row1_string);
    MUC_Mat_row0_stream>>MUC_Mat(0,0)>>MUC_Mat(0,1);
    MUC_Mat_row1_stream>>MUC_Mat(1,0)>>MUC_Mat(1,1);


    int Just_Hartree_int;
    Just_Hartree_int = int(matchstring(inputfile_, "Just_Hartree"));
    if(Just_Hartree_int==1){
        Just_Hartree=true;
    }
    else{
        Just_Hartree=false;
        FockType = matchstring2(inputfile_,"Fock_type");
    }

    int NoSpinFlipOP_int;
    NoSpinFlipOP_int=int(matchstring(inputfile_, "NoSpinFlipOP"));
    if(NoSpinFlipOP_int==1){
        NoSpinFlipOP=true;
    }
    else{
        NoSpinFlipOP=false;
    }



    int Create_OPs_Ansatz_int;
    Create_OPs_Ansatz_int = int(matchstring(inputfile_,"Create_OPs_Ansatz"));
    if(Create_OPs_Ansatz_int==1){
        Create_OPs_Ansatz=true;
        OP_Ansatz_type=matchstring2(inputfile_,"OP_Ansatz_type");
    }
    else{
        Create_OPs_Ansatz=false;
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



    OnSiteE.resize(n_atoms*n_orbs);
    for (int orb=0;orb<n_atoms*n_orbs;orb++){
        OnSiteE[orb].resize(2);
    }

    PinningFieldZ.resize(n_atoms*n_orbs);



    string OnSiteE_string = matchstring2(inputfile_, "OnSiteE");
    stringstream OnSiteE_stream(OnSiteE_string);
    for(int atom_no=0;atom_no<n_atoms;atom_no++){
        for(int orb_no=0;orb_no<n_orbs;orb_no++){
            OnSiteE_stream>>OnSiteE[atom_no+n_atoms*orb_no][0];
            OnSiteE[atom_no+n_atoms*orb_no][1]=OnSiteE[atom_no+n_atoms*orb_no][0];
        }}


    string PinningFieldValueZ_string = matchstring2(inputfile_, "PinningFieldValueZ");
    stringstream PinningFieldValueZ_stream(PinningFieldValueZ_string);
    for(int atom_no=0;atom_no<n_atoms;atom_no++){
        for(int orb_no=0;orb_no<n_orbs;orb_no++){
            PinningFieldValueZ_stream>>PinningFieldZ[atom_no+n_atoms*orb_no];
        }}



    PinningFieldX.resize(n_atoms*n_orbs);
    string PinningFieldValueX_string = matchstring2(inputfile_, "PinningFieldValueX");
    stringstream PinningFieldValueX_stream(PinningFieldValueX_string);
    for(int atom_no=0;atom_no<n_atoms;atom_no++){
        for(int orb_no=0;orb_no<n_orbs;orb_no++){
            PinningFieldValueX_stream>>PinningFieldX[atom_no+n_atoms*orb_no];
    }}


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

    //PinningFieldValue=matchstring(inputfile_,"PinningFieldValue");

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




    omega_min=matchstring(inputfile_,"omega_min");
    omega_max=matchstring(inputfile_,"omega_max");
    d_omega=matchstring(inputfile_,"d_omega");
    eta=matchstring(inputfile_,"eta");



    cout << "____________________________________" << endl;
}

double Parameters_G2dLatticeNew::matchstring(string file, string match)
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

string Parameters_G2dLatticeNew::matchstring2(string file, string match)
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
        cout << "Unable to open input file while in the Parameters_G2dLattice class." << endl;
    }

    cout << match << " " << amount << endl;
    return amount;
}

#endif
