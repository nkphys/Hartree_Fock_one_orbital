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
index = site + spin*Lx*Ly;
*/
    //Define Fields
    Matrix_COO_Complex OParams_;

    Matrix<double> Disorder;

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
    int lx_,ly_,ns_, no_dof_;
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

    lx_=Coordinates_.lx_;
    ly_=Coordinates_.ly_;
    no_dof_=Coordinates_.no_dof_;
    ns_=Coordinates_.ns_;
    int UP_, DOWN_;
    UP_=0;DOWN_=1;

    // srand(Parameters_.RandomSeed);

    SI_to_ind.resize(2*ns_*2*ns_);
    Disorder.resize(lx_,ly_);
    OParams_.nrows=2*ns_;
    OParams_.ncols=2*ns_;
    OParams_.value.clear();
    OParams_.rows.clear();
    OParams_.columns.clear();
    complex<double> comp_temp;
    int row_temp, col_temp;


    ofstream Disorder_conf_file("Disorder_conf_used");
    Disorder_conf_file<<"#seed="<<Parameters_.RandomDisorderSeed<<
                        " for mt19937_64 Generator is used"<<endl;
    Disorder_conf_file<<"#ix   iy    Dis[ix,iy]"<<endl;

    ofstream Initial_OrderParams_file("Initial_OrderParams_values_generated.txt");

    if(!Parameters_.Read_OPs){

        if(!Parameters_.Create_OPs){

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

                                OParams_.value.push_back(comp_temp);
                                OParams_.rows.push_back(row_temp);
                                OParams_.columns.push_back(col_temp);
                                SI_to_ind[row_temp + (2*ns_*col_temp)] = OParams_.value.size() - 1;

                            }
                        }
                    }
                }

            }


            Initial_OrderParams_file<<"#seed="<<Parameters_.RandomSeed<<
                                      " for mt19937_64 Generator is used"<<endl;



        }
        else{

            cout<<"NEED TO ADD ANSATZ STATES"<<endl;


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
            while(!file_initial_OP_in.eof())
            {
                getline(file_initial_OP_in,line_temp);
                stringstream line_temp_ss(line_temp);
                line_temp_ss >> alpha_i >> alpha_j >> val_OP;
                OParams_.value.push_back(val_OP);
                OParams_.rows.push_back(alpha_i);
                OParams_.columns.push_back(alpha_j);
                SI_to_ind[row_temp + (2*ns_*col_temp)] = OParams_.value.size() - 1;

            }
            file_initial_OP_in.close();
        }
        else
        {cout<<"Unable to open file = '"<<fl_initial_OP_in<<"'"<<endl;}


        Initial_OrderParams_file<<"#OParams are read from '"<<Parameters_.File_OPs_in<<"' file"<<endl;

    }


    Initial_OrderParams_file<<"#alpha_i alpha_j OParams_[alpha_i][alpha_j]"<<endl;
    for(int alpha=0;alpha<OParams_.value.size();alpha++){
        Initial_OrderParams_file<<OParams_.rows[alpha]<<setw(15)<<OParams_.columns[alpha]<<setw(15)<<OParams_.value[alpha]<<endl;
    }



    int temp_i, temp_j;
    if(Parameters_.ReadDisorder==false){
        cout <<"Disorder conf is initialized using random seed given in input file"<< endl;
        for(int j=0;j<ly_;j++){
            for(int i=0;i<lx_;i++){
                //RANDOM Disorder
                Disorder(i,j)=Parameters_.Disorder_Strength*((2.0*random2())-1.0);
                Disorder_conf_file<<i<<"  "<<j<<"  "<<Disorder(i,j)<<endl;
            }
            Disorder_conf_file<<endl;
        }
    }
    else{
        cout<<"Disorder conf is read from file path given in the input file"<<endl;
        ifstream Disorder_in_file(Parameters_.DisorderSeedFile);
        for(int j=0;j<ly_;j++){
            for(int i=0;i<lx_;i++){
                //RANDOM Disorder
                Disorder_in_file>>temp_i>>temp_j>>Disorder(i,j);
                assert(i==temp_i);
                assert(j==temp_j);
                Disorder(i,j)=Parameters_.Disorder_Strength*(Disorder(i,j));
                Disorder_conf_file<<i<<"  "<<j<<"  "<<Disorder(i,j)<<endl;
            }
            Disorder_conf_file<<endl;
        }
    }


} // ----------

#endif
