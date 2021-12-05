#include <math.h>
#include "tensor_type.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#include "random"
#include <stdlib.h>
#define PI acos(-1.0)

#ifndef MFParams_class
#define MFParams_class

class MFParams{ 
public:

    /* Convention
index = site + spin*N_sites;
*/
    //Define Fields
    Matrix_COO_Complex OParams_;

    Mat_1_doub Disorder;

    // Constructor
    MFParams(Parameters& Parameters__, Coordinates&  Coordinates__, mt19937_64& Generator1__ , mt19937_64& Generator2__)
        :Parameters_(Parameters__),Coordinates_(Coordinates__), Generator1_(Generator1__), Generator2_(Generator2__)
    {
        //setupm_arr();
        initialize();
    }


    double random1();
    double random2();
    void initialize();


    Parameters &Parameters_;
    Coordinates &Coordinates_;
    mt19937_64 &Generator1_; //for random fields
    mt19937_64 &Generator2_; //for random disorder
    int ns_, no_dof_;
    Mat_1_int SI_to_ind;

    uniform_real_distribution<double> dis1_;//for random fields
    uniform_real_distribution<double> dis2_;//for random disorder

    //mt19937_64 mt_rand(Parameters_.RandomSeed);


};



double MFParams::random1(){

    return dis1_(Generator1_);

}

double MFParams::random2(){

    return dis2_(Generator2_);

}


void MFParams::initialize(){


    no_dof_=Coordinates_.no_dof_;
    ns_=Coordinates_.ns_;
    int UP_, DOWN_;
    UP_=0;DOWN_=1;

    // srand(Parameters_.RandomSeed);
    SI_to_ind.resize(2*ns_*2*ns_);
    Disorder.resize(ns_);
    OParams_.nrows=2*ns_;
    OParams_.ncols=2*ns_;
    OParams_.value.clear();
    OParams_.rows.clear();
    OParams_.columns.clear();
    complex<double> comp_temp;
    complex<double> value_temp;
    int row_temp, col_temp;


    ofstream Disorder_conf_file("Disorder_conf_used");
    Disorder_conf_file<<"#seed="<<Parameters_.RandomDisorderSeed<<
                        " for mt19937_64 Generator is used"<<endl;
    Disorder_conf_file<<"#site    Dis[site]"<<endl;

    ofstream Initial_OrderParams_file("Initial_OrderParams_values_generated.txt");

    if(!Parameters_.Read_OPs){

        if(!Parameters_.Create_OPs){

            for(int site_i=0;site_i<ns_;site_i++){
                for(int site_j=0;site_j<ns_;site_j++){

                    if(site_i!=site_j){
                        if(abs(Parameters_.LongRangeInteractions[site_i][site_j])>0.0000001){

                            if(!Parameters_.Just_Hartree){

                                for(int spin_i=0;spin_i<2;spin_i++){
                                    for(int spin_j=spin_i;spin_j<2;spin_j++){
                                        row_temp = Coordinates_.Nc_dof(site_i, spin_i);
                                        col_temp = Coordinates_.Nc_dof(site_j, spin_j);

                                        if(spin_i==spin_j){
                                            if(site_j>site_i){

                                                comp_temp.real(random1());
                                                comp_temp.imag(random1());
                                                OParams_.value.push_back(comp_temp);
                                                OParams_.rows.push_back(row_temp);
                                                OParams_.columns.push_back(col_temp);
                                                SI_to_ind[row_temp + (2*ns_*col_temp)] = OParams_.value.size() - 1;
                                            }
                                        }
                                        else{
                                            comp_temp.real(random1());
                                            comp_temp.imag(random1());
                                            OParams_.value.push_back(comp_temp);
                                            OParams_.rows.push_back(row_temp);
                                            OParams_.columns.push_back(col_temp);
                                            SI_to_ind[row_temp + (2*ns_*col_temp)] = OParams_.value.size() - 1;

                                        }
                                    }
                                }

                            }
                        }
                    }
                    else{
                        for(int spin_i=0;spin_i<2;spin_i++){
                            for(int spin_j=spin_i;spin_j<2;spin_j++){

                                row_temp = Coordinates_.Nc_dof(site_i, spin_i);
                                col_temp = Coordinates_.Nc_dof(site_j, spin_j);

                                comp_temp.real(random1());
                                if(spin_i==spin_j){
                                    comp_temp.imag(0.0);
                                }
                                else{
                                    comp_temp.imag(random1());
                                }

                                if(Parameters_.Just_Hartree){
                                    if(row_temp==col_temp){
                                        OParams_.value.push_back(comp_temp*(Parameters_.Total_Particles/(2*ns_)));
                                        OParams_.rows.push_back(row_temp);
                                        OParams_.columns.push_back(col_temp);
                                        SI_to_ind[row_temp + (2*ns_*col_temp)] = OParams_.value.size() - 1;
                                    }}
                                else{
                                    OParams_.value.push_back(comp_temp);
                                    OParams_.rows.push_back(row_temp);
                                    OParams_.columns.push_back(col_temp);
                                    SI_to_ind[row_temp + (2*ns_*col_temp)] = OParams_.value.size() - 1;
                                }

                            }
                        }
                    }
                }

            }


            Initial_OrderParams_file<<"#seed="<<Parameters_.RandomSeed<<
                                      " for mt19937_64 Generator is used"<<endl;

        }
        else{

            if(Parameters_.Ansatz=="Given_CDW"){
                assert(Parameters_.Restricted_HF);
                cout<<"Creating Anstaz in Restricted Space HF"<<endl;

                Parameters_.A_charge_modulation = random1()*0.5;
                for(int site_i=0;site_i<ns_;site_i++){
                    for(int site_j=0;site_j<ns_;site_j++){

                        if(site_i!=site_j){
                            if(abs(Parameters_.LongRangeInteractions[site_i][site_j])>0.0000001){

                                for(int spin_i=0;spin_i<2;spin_i++){
                                    for(int spin_j=spin_i;spin_j<2;spin_j++){
                                        row_temp = Coordinates_.Nc_dof(site_i, spin_i);
                                        col_temp = Coordinates_.Nc_dof(site_j, spin_j);

                                        if(spin_i==spin_j){
                                            if(site_j>site_i){
                                                comp_temp.real(0.0);
                                                comp_temp.imag(0.0);
                                                OParams_.value.push_back(comp_temp);
                                                OParams_.rows.push_back(row_temp);
                                                OParams_.columns.push_back(col_temp);
                                                SI_to_ind[row_temp + (2*ns_*col_temp)] = OParams_.value.size() - 1;

                                            }
                                        }
                                        else{
                                            comp_temp.real(0.0);
                                            comp_temp.imag(0.0);
                                            OParams_.value.push_back(comp_temp);
                                            OParams_.rows.push_back(row_temp);
                                            OParams_.columns.push_back(col_temp);
                                            SI_to_ind[row_temp + (2*ns_*col_temp)] = OParams_.value.size() - 1;

                                        }
                                    }
                                }
                            }
                        }
                        else{
                            int spin_i_, spin_j_;

                            //DOWN, UP
                            spin_j_=DOWN_;
                            spin_i_=UP_;
                            row_temp = Coordinates_.Nc_dof(site_i, spin_i_);
                            col_temp = Coordinates_.Nc_dof(site_j, spin_j_);
                            comp_temp=complex<double>(random1(),random1());
                            OParams_.value.push_back(comp_temp);
                            OParams_.rows.push_back(row_temp);
                            OParams_.columns.push_back(col_temp);
                            SI_to_ind[row_temp + (2*ns_*col_temp)] = OParams_.value.size() - 1;

                            comp_temp=complex<double>(random1(),0.0);//  sz_i
                            comp_temp = comp_temp*0.5*(0.5 + Parameters_.A_charge_modulation*Parameters_.CDW_Ansatz_sites[site_i]);  //for |sz_i|<=n_i/2
                            //UP, UP
                            value_temp = 0.5*(0.5 + Parameters_.A_charge_modulation*Parameters_.CDW_Ansatz_sites[site_i]) + comp_temp;
                            spin_j_=UP_;
                            spin_i_=UP_;
                            row_temp = Coordinates_.Nc_dof(site_i, spin_i_);
                            col_temp = Coordinates_.Nc_dof(site_j, spin_j_);
                            OParams_.value.push_back(value_temp);
                            OParams_.rows.push_back(row_temp);
                            OParams_.columns.push_back(col_temp);
                            SI_to_ind[row_temp + (2*ns_*col_temp)] = OParams_.value.size() - 1;

                            //DOWN, DOWN
                            value_temp = 0.5*(0.5 + Parameters_.A_charge_modulation*Parameters_.CDW_Ansatz_sites[site_i]) - comp_temp;
                            spin_j_=DOWN_;
                            spin_i_=DOWN_;
                            row_temp = Coordinates_.Nc_dof(site_i, spin_i_);
                            col_temp = Coordinates_.Nc_dof(site_j, spin_j_);
                            OParams_.value.push_back(value_temp);
                            OParams_.rows.push_back(row_temp);
                            OParams_.columns.push_back(col_temp);
                            SI_to_ind[row_temp + (2*ns_*col_temp)] = OParams_.value.size() - 1;

                        }

                    }

                }

            }
            else if(Parameters_.Ansatz=="Given_Tetrahedron"){
                cout<<"Creating Tetrahedron Ansatz using file = '"<<Parameters_.Ansatz_file<<"'"<<endl;
                int uc_x=2;
                int uc_y=2;
                int lx;
                int ly;


                ifstream inputfile_Ansatz_file(Parameters_.Ansatz_file.c_str());
                string line_temp_; //from k-space code
                int dof1, dof2;
                int alpha1, alpha1x, alpha1y, alpha2, alpha2x, alpha2y;
                int spin1, spin2;
                int site_x, site_y;
                Mat_3_Complex_doub OP_TEMP;
                OP_TEMP.resize(4);
                for(int site=0;site<4;site++){
                    OP_TEMP[site].resize(2);
                    for(int spin=0;spin<2;spin++){
                        OP_TEMP[site][spin].resize(2);
                    }
                }

                complex<double> value_OP;
                getline(inputfile_Ansatz_file,line_temp_);
                stringstream line_temp_sizes_(line_temp_);
                line_temp_sizes_ >> lx >> ly;
                int cells_x=lx/2;
                int cells_y=ly/2;
                int nx_cell, ny_cell;
                assert(ns_=lx*ly);


                getline(inputfile_Ansatz_file,line_temp_);
                while(getline(inputfile_Ansatz_file,line_temp_))
                {
                    stringstream line_temp_ss_(line_temp_);
                    line_temp_ss_ >> dof1 >> dof2 >> value_OP;

                    alpha1=dof1%4;
                    spin1=(dof1-alpha1)/4;
                    alpha1x=alpha1%2;
                    alpha1y=(alpha1-alpha1x)/2;

                    alpha2=dof2%4;
                    spin2=(dof2-alpha2)/4;
                    alpha2x=alpha2%2;
                    alpha2y=(alpha2-alpha2x)/2;

                    assert(alpha1==alpha2);

                    OP_TEMP[alpha1][spin1][spin2]=value_OP;
                   // cout<<alpha1<<"  "<<spin1<<"  "<<spin2<<"  "<<value_OP<<endl;
                }



                for(int site_i=0;site_i<ns_;site_i++){
                    for(int site_j=0;site_j<ns_;site_j++){

                        if(site_i!=site_j){
                            if(abs(Parameters_.LongRangeInteractions[site_i][site_j])>0.0000001){

                                for(int spin_i=0;spin_i<2;spin_i++){
                                    for(int spin_j=spin_i;spin_j<2;spin_j++){
                                        row_temp = Coordinates_.Nc_dof(site_i, spin_i);
                                        col_temp = Coordinates_.Nc_dof(site_j, spin_j);

                                        if(spin_i==spin_j){
                                            if(site_j>site_i){
                                                comp_temp.real(0.0);
                                                comp_temp.imag(0.0);
                                                OParams_.value.push_back(comp_temp);
                                                OParams_.rows.push_back(row_temp);
                                                OParams_.columns.push_back(col_temp);
                                                SI_to_ind[row_temp + (2*ns_*col_temp)] = OParams_.value.size() - 1;

                                            }
                                        }
                                        else{
                                            comp_temp.real(0.0);
                                            comp_temp.imag(0.0);
                                            OParams_.value.push_back(comp_temp);
                                            OParams_.rows.push_back(row_temp);
                                            OParams_.columns.push_back(col_temp);
                                            SI_to_ind[row_temp + (2*ns_*col_temp)] = OParams_.value.size() - 1;

                                        }
                                    }
                                }
                            }
                        }
                        else{

                            site_x=site_i%lx;
                            site_y=(site_i - site_x)/lx;
                            nx_cell = site_x/2;
                            ny_cell =site_y/2;
                            alpha1x= site_x - 2*nx_cell;
                            alpha1y= site_y - 2*ny_cell;
                            alpha1 =alpha1x + 2*alpha1y;

                            int spin_i_, spin_j_;

                            //DOWN, UP
                            spin_j_=DOWN_;
                            spin_i_=UP_;
                            row_temp = Coordinates_.Nc_dof(site_i, spin_i_);
                            col_temp = Coordinates_.Nc_dof(site_j, spin_j_);
                            comp_temp=OP_TEMP[alpha1][spin_i_][spin_j_];
                            OParams_.value.push_back(comp_temp);
                            OParams_.rows.push_back(row_temp);
                            OParams_.columns.push_back(col_temp);
                            SI_to_ind[row_temp + (2*ns_*col_temp)] = OParams_.value.size() - 1;

                           // comp_temp=complex<double>(random1(),0.0);//  sz_i
                           // comp_temp = comp_temp*0.5*(0.5 + Parameters_.A_charge_modulation*Parameters_.CDW_Ansatz_sites[site_i]);  //for |sz_i|<=n_i/2
                            //UP, UP
                            spin_j_=UP_;
                            spin_i_=UP_;
                            row_temp = Coordinates_.Nc_dof(site_i, spin_i_);
                            col_temp = Coordinates_.Nc_dof(site_j, spin_j_);
                            comp_temp=OP_TEMP[alpha1][spin_i_][spin_j_];
                            OParams_.value.push_back(comp_temp);
                            OParams_.rows.push_back(row_temp);
                            OParams_.columns.push_back(col_temp);
                            SI_to_ind[row_temp + (2*ns_*col_temp)] = OParams_.value.size() - 1;

                            //DOWN, DOWN
                            spin_j_=DOWN_;
                            spin_i_=DOWN_;
                            row_temp = Coordinates_.Nc_dof(site_i, spin_i_);
                            col_temp = Coordinates_.Nc_dof(site_j, spin_j_);
                            comp_temp=OP_TEMP[alpha1][spin_i_][spin_j_];
                            OParams_.value.push_back(comp_temp);
                            OParams_.rows.push_back(row_temp);
                            OParams_.columns.push_back(col_temp);
                            SI_to_ind[row_temp + (2*ns_*col_temp)] = OParams_.value.size() - 1;

                        }

                    }

                }

            }

        }


    }
    else{

        double small_=0.01;
        //Adding OPs's to  OParams_ to be small value, which in principle can be finite
        for(int site_i=0;site_i<ns_;site_i++){
            for(int site_j=0;site_j<ns_;site_j++){

                if(site_i!=site_j){
                    if(abs(Parameters_.LongRangeInteractions[site_i][site_j])>0.0000001){

                        if(!Parameters_.Just_Hartree){

                            for(int spin_i=0;spin_i<2;spin_i++){
                                for(int spin_j=spin_i;spin_j<2;spin_j++){
                                    row_temp = Coordinates_.Nc_dof(site_i, spin_i);
                                    col_temp = Coordinates_.Nc_dof(site_j, spin_j);

                                    if(spin_i==spin_j){
                                        if(site_j>site_i){

                                            comp_temp.real(small_);
                                            comp_temp.imag(small_);
                                            OParams_.value.push_back(comp_temp);
                                            OParams_.rows.push_back(row_temp);
                                            OParams_.columns.push_back(col_temp);
                                            SI_to_ind[row_temp + (2*ns_*col_temp)] = OParams_.value.size() - 1;
                                        }
                                    }
                                    else{
                                        comp_temp.real(small_);
                                        comp_temp.imag(small_);
                                        OParams_.value.push_back(comp_temp);
                                        OParams_.rows.push_back(row_temp);
                                        OParams_.columns.push_back(col_temp);
                                        SI_to_ind[row_temp + (2*ns_*col_temp)] = OParams_.value.size() - 1;

                                    }
                                }
                            }

                        }
                    }
                }
                else{
                    for(int spin_i=0;spin_i<2;spin_i++){
                        for(int spin_j=spin_i;spin_j<2;spin_j++){

                            row_temp = Coordinates_.Nc_dof(site_i, spin_i);
                            col_temp = Coordinates_.Nc_dof(site_j, spin_j);

                            comp_temp.real(0.0);
                            comp_temp.imag(0.0);

                            if(Parameters_.Just_Hartree){
                                if(row_temp==col_temp){
                                    OParams_.value.push_back(comp_temp);
                                    OParams_.rows.push_back(row_temp);
                                    OParams_.columns.push_back(col_temp);
                                    SI_to_ind[row_temp + (2*ns_*col_temp)] = OParams_.value.size() - 1;
                                }}
                            else{
                                OParams_.value.push_back(comp_temp);
                                OParams_.rows.push_back(row_temp);
                                OParams_.columns.push_back(col_temp);
                                SI_to_ind[row_temp + (2*ns_*col_temp)] = OParams_.value.size() - 1;
                            }

                        }
                    }
                }
            }

        }



        string fl_initial_OP_in = Parameters_.File_OPs_in;
        ifstream file_initial_OP_in(fl_initial_OP_in.c_str());
        string temp1, line_temp;
        int alpha_i, alpha_j;
        complex<double> val_OP;

        //#row col <c^{dagger}_{site_i,spin_i} c_{site_j, spin_j}>
        getline(file_initial_OP_in,temp1);


        if(file_initial_OP_in.is_open())
        {
            while(getline(file_initial_OP_in,line_temp))
            {
                stringstream line_temp_ss(line_temp);
                line_temp_ss >> alpha_i >> alpha_j >> val_OP;
                //cout<<alpha_i<<"  "<<alpha_j<<"  "<<val_OP<<endl;
                if(Parameters_.Just_Hartree){
                    if(alpha_i==alpha_j){
                        OParams_.value[SI_to_ind[alpha_i + (2*ns_*alpha_j)]]=val_OP;
                        assert(OParams_.rows[SI_to_ind[alpha_i + (2*ns_*alpha_j)]]==alpha_i);
                        assert(OParams_.columns[SI_to_ind[alpha_i + (2*ns_*alpha_j)]]==alpha_j);

                    }
                }
                else{
                    OParams_.value[SI_to_ind[alpha_i + (2*ns_*alpha_j)]]=val_OP;
                    //cout<<SI_to_ind[alpha_i + (2*ns_*alpha_j)]<<endl;
                    assert(OParams_.rows[SI_to_ind[alpha_i + (2*ns_*alpha_j)]]==alpha_i);
                    assert(OParams_.columns[SI_to_ind[alpha_i + (2*ns_*alpha_j)]]==alpha_j);
                }

            }
            file_initial_OP_in.close();
        }
        else
        {cout<<"Unable to open file = '"<<fl_initial_OP_in<<"'"<<endl;}


        Initial_OrderParams_file<<"#OParams are read from '"<<Parameters_.File_OPs_in<<"' file"<<endl;

    }


    Initial_OrderParams_file<<"#alpha_i alpha_j OParams_[alpha_i][alpha_j]"<<endl;
    for(int alpha=0;alpha<OParams_.value.size();alpha++){
        Initial_OrderParams_file<<OParams_.rows[alpha]<<setw(15)<<OParams_.columns[alpha]<<setw(15)<<"  "<<OParams_.value[alpha]<<endl;
    }


    cout<<"Total no. of Op's = "<<OParams_.value.size()<<endl;



    int temp_j;
    if(Parameters_.ReadDisorder==false){
        cout <<"Disorder conf is initialized using random seed given in input file"<< endl;
        for(int j=0;j<ns_;j++){
            //RANDOM Disorder
            Disorder[j]=Parameters_.Disorder_Strength*((2.0*random2())-1.0);
            Disorder_conf_file<<j<<"  "<<Disorder[j]<<endl;
        }
    }
    else{
        cout<<"Disorder conf is read from file path given in the input file"<<endl;
        ifstream Disorder_in_file(Parameters_.DisorderSeedFile);
        for(int j=0;j<ns_;j++){
            //RANDOM Disorder
            Disorder_in_file>>temp_j>>Disorder[j];
            assert(j==temp_j);
            Disorder[j]=Parameters_.Disorder_Strength*(Disorder[j]);
            Disorder_conf_file<<j<<"  "<<Disorder[j]<<endl;

        }
    }


} // ----------

#endif
