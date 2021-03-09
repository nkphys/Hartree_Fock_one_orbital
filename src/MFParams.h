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
                                        OParams_.value.push_back(comp_temp*(Parameters_.Total_Particles/2*ns_));
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

            assert(Parameters_.Restricted_HF);

            if(Parameters_.Ansatz=="Given_CDW"){
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


        }


    }
    else{

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

                if(Parameters_.Just_Hartree){
                    if(alpha_i==alpha_j){
                        OParams_.value.push_back(val_OP);
                        OParams_.rows.push_back(alpha_i);
                        OParams_.columns.push_back(alpha_j);
                        SI_to_ind[alpha_i + (2*ns_*alpha_j)] = OParams_.value.size() - 1;
                    }
                }
                else{
                    OParams_.value.push_back(val_OP);
                    OParams_.rows.push_back(alpha_i);
                    OParams_.columns.push_back(alpha_j);
                    SI_to_ind[alpha_i + (2*ns_*alpha_j)] = OParams_.value.size() - 1;
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
