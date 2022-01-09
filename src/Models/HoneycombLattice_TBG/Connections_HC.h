#include <algorithm>
#include <functional>
#include <math.h>
#include "tensor_type.h"
#include "Parameters_HC.h"
#include "Coordinates_HC.h"
#include "../../Matrix.h"
#define PI acos(-1.0)

#ifndef Connections_HC_class
#define Connections_HC_class


class Connections_HC
{
public:
    Connections_HC(Parameters_HC &Parameters__, Coordinates_HC &Coordinates__)
        : Parameters_(Parameters__), Coordinates_(Coordinates__)

    {
        Initialize();
        HTBCreate();
    }

    void Initialize();                                     //::DONE
    void Print_Hopping();                                       //::DONE


    void Print_LongRangeInt();
    void Print_Spin_resolved_OnsiteE();
    void InteractionsCreate();                             //::DONE
    void Interactions_Sorting();
    void Check_Hermiticity();                              //::DONE
    void Check_up_down_symmetry();                         //::DONE
    void HTBCreate();                                      //::DONE
    void Print_Ansatz_LocalDen_CDW();
    void Get_minimum_distance_direction(int l, int gamma_l, int m, int gamma_m, int &r1_, int &r2_, double &dis_min);



    Parameters_HC &Parameters_;
    Coordinates_HC &Coordinates_;
    int lx_, ly_, nsites_, ncells_, n_orbs_;
    Matrix<complex<double>> HTB_;
    Matrix<complex<double>> Hint_;
    Matrix<complex<double>> Ham_;
    Matrix<double> Tx, Ty, Tpxpy, Tpxmy;
    Matrix<double> IntraCell_Hopp, InterCell_px_Hopp, InterCell_py_Hopp, InterCell_pxmy_Hopp;

    int max_l1_dis, max_l2_dis;
};


void Connections_HC::Get_minimum_distance_direction(int l,int gamma_l, int m, int gamma_m, int &r1_, int &r2_, double &dis_min){

    int r1_a,r1_b, r2_a, r2_b;
    int l1_pos, l2_pos;
    int m1_pos, m2_pos;
    double rx_, ry_;
    double dis_x_offset_a, dis_x_offset_b;
    double dis_y_offset_a, dis_y_offset_b;

    l1_pos = Coordinates_.indx_cellwise(l);
    l2_pos = Coordinates_.indy_cellwise(l);
    m1_pos = Coordinates_.indx_cellwise(m);
    m2_pos = Coordinates_.indy_cellwise(m);

    r1_a=m1_pos-l1_pos;
    if(r1_a>0){
        dis_x_offset_a= (gamma_m - gamma_l)*(1.0/2.0);

        r1_b=-1*(lx_-r1_a);
        dis_x_offset_b= (gamma_m - gamma_l)*(1.0/2.0);

    }
    else if (r1_a<0){
        dis_x_offset_a= (gamma_m - gamma_l)*(1.0/2.0);

        r1_b=lx_-abs(r1_a);
        dis_x_offset_b= (gamma_m - gamma_l)*(1.0/2.0);

    }
    else{
        dis_x_offset_a= (gamma_m - gamma_l)*(1.0/2.0);
        dis_x_offset_b= (gamma_m - gamma_l)*(1.0/2.0);
        r1_b=0;
        assert(r1_a==0);
    }

    r2_a=m2_pos-l2_pos;
    if(r2_a>0){
        dis_y_offset_a= (gamma_m - gamma_l)*(1.0/(2.0*sqrt(3.0)));

        r2_b=-1*(ly_-r2_a);
        dis_y_offset_b= (gamma_m - gamma_l)*(1.0/(2.0*sqrt(3.0)));
    }
    else if (r2_a<0){
        dis_y_offset_a= (gamma_m - gamma_l)*(1.0/(2.0*sqrt(3.0)));

        r2_b=ly_-abs(r2_a);
        dis_y_offset_b= (gamma_m - gamma_l)*(1.0/(2.0*sqrt(3.0)));
    }
    else{
        dis_y_offset_a= (gamma_m - gamma_l)*(1.0/(2.0*sqrt(3.0)));
        dis_y_offset_b= (gamma_m - gamma_l)*(1.0/(2.0*sqrt(3.0)));
        r2_b=0;
        assert(r2_a==0);
    }


    double min_dis=1000000000.0;
    double dis;

    //r1a r2a
    rx_ = ((1.0)*(r1_a) +  (1.0/2.0)*(r2_a)) + dis_x_offset_a;
    ry_ =  (0.0*(r1_a) + (sqrt(3.0)/2.0)*(r2_a)) + dis_y_offset_a;
    dis= sqrt(rx_*rx_ + ry_*ry_);
    if(dis<=min_dis){
        r1_=r1_a;
        r2_=r2_a;
        min_dis=dis;
    }

    //r1a r2b
    rx_ = ((1.0)*(r1_a) +  (1.0/2.0)*(r2_b)) + dis_x_offset_a;
    ry_ =  (0.0*(r1_a) + (sqrt(3.0)/2.0)*(r2_b)) + dis_y_offset_b;
    dis= sqrt(rx_*rx_ + ry_*ry_);
    if(dis<=min_dis){
        r1_=r1_a;
        r2_=r2_b;
        min_dis=dis;
    }

    //r1b r2a
    rx_ = ((1.0)*(r1_b) +  (1.0/2.0)*(r2_a)) + dis_x_offset_b;
    ry_ =  (0.0*(r1_b) + (sqrt(3.0)/2.0)*(r2_a)) + dis_y_offset_a;
    dis= sqrt(rx_*rx_ + ry_*ry_);
    if(dis<=min_dis){
        r1_=r1_b;
        r2_=r2_a;
        min_dis=dis;
    }

    //r1b r2b
    rx_ = ((1.0)*(r1_b) +  (1.0/2.0)*(r2_b)) + dis_x_offset_b;
    ry_ =  (0.0*(r1_b) + (sqrt(3.0)/2.0)*(r2_b)) + dis_y_offset_b;
    dis= sqrt(rx_*rx_ + ry_*ry_);
    if(dis<=min_dis){
        r1_=r1_b;
        r2_=r2_b;
        min_dis=dis;
    }


    dis_min = min_dis;
}

void Connections_HC::Initialize()
{


    ly_ = Parameters_.ly;
    lx_ = Parameters_.lx;
    nsites_ = lx_ * ly_;  //total no of sites, each site has 2 orb
    n_orbs_ = Parameters_.n_orbs;
    int space = 2 * n_orbs_*nsites_ * n_orbs_;

    HTB_.resize(space, space);
    Ham_.resize(space, space);
    Hint_.resize(n_orbs_*nsites_,n_orbs_*nsites_);

    if(lx_%2==0){
        max_l1_dis = lx_/2;
    }
    else{
        max_l1_dis = (lx_-1)/2;
    }

    if(ly_%2==0){
        max_l2_dis = ly_/2;
    }
    else{
        max_l2_dis = (ly_-1)/2;
    }


} // ----------


void Connections_HC::Print_Ansatz_LocalDen_CDW(){



    Mat_1_string Ansatz_all;
    Ansatz_all.push_back("Stripe_CDW_NM");
    Ansatz_all.push_back("Stripe_CDW_FM");
    Ansatz_all.push_back("Stripe_CDW_AFM");
    Ansatz_all.push_back("Stripe_CDW_FM2");
    Ansatz_all.push_back("Stripe_CDW_AFM2");
    Ansatz_all.push_back("DiagonalStripe_CDW_AFM");
    Ansatz_all.push_back("Honeycomb_CDW_FM");
    Ansatz_all.push_back("Honeycomb_CDW_AFM");
    Ansatz_all.push_back("size_two_eq_triangles_CDW_FM");
    Ansatz_all.push_back("doubled_unitcell_CDW_FM");


    for(int str_no=0;str_no<Ansatz_all.size();str_no++){
        string ansatz=Ansatz_all[str_no];
        string fileout="Hartree_OP_"+ansatz+".txt";
        ofstream File_Out(fileout.c_str());

        int alpha_i;

        if(ansatz=="Stripe_CDW_NM"){
            //Stripe_CDW_NM
            File_Out<<"#alpha_i   alpha_j   OParams_[alpha_i][alpha_j]   // FOR STRIPE CDW"<<endl;
            for(int ix=0;ix<lx_;ix++){
                for(int iy=0;iy<ly_;iy++){
                    for(int spin_i=0;spin_i<2;spin_i++){
                        alpha_i = Coordinates_.Nbasis(ix,iy,0) + spin_i*(ncells_);

                        if(ix%2==0){
                            File_Out<< alpha_i<<"   "<<alpha_i<<"   "<<one_complex*0.5<<endl;
                        }

                    }
                }
            }
        }

        if(ansatz=="Stripe_CDW_FM"){
            int spin_i;
            //Stripe_CDW_AFM
            File_Out<<"#alpha_i   alpha_j   OParams_[alpha_i][alpha_j]   // FOR STRIPE CDW"<<endl;
            for(int ix=0;ix<lx_;ix++){
                for(int iy=0;iy<ly_;iy++){

                    spin_i=0;

                    alpha_i = Coordinates_.Nbasis(ix,iy,0) + spin_i*(ncells_);
                    if(ix%2==0){
                        File_Out<< alpha_i<<"   "<<alpha_i<<"   "<<one_complex<<endl;
                    }

                }
            }
        }

        if(ansatz=="Stripe_CDW_FM2"){
            int spin_i;

            File_Out<<"#alpha_i   alpha_j   OParams_[alpha_i][alpha_j]   // FOR STRIPE CDW"<<endl;
            for(int ix=0;ix<lx_;ix++){
                for(int iy=0;iy<ly_;iy++){

                    if((ix/2)%2==0){
                        spin_i=0;
                    }
                    else{
                        spin_i=1;
                    }

                    alpha_i = Coordinates_.Nbasis(ix,iy,0) + spin_i*(ncells_);
                    if(ix%2==0){
                        File_Out<< alpha_i<<"   "<<alpha_i<<"   "<<one_complex<<endl;
                    }

                }
            }
        }

        if(ansatz=="Stripe_CDW_AFM2"){
            int spin_i;
            //Stripe_CDW_AFM2
            File_Out<<"#alpha_i   alpha_j   OParams_[alpha_i][alpha_j]   // FOR STRIPE CDW"<<endl;
            for(int ix=0;ix<lx_;ix++){
                for(int iy=0;iy<ly_;iy++){

                    if((ix/2)%2==0){
                        if(iy%2==0){
                            spin_i=0;
                        }
                        else{
                            spin_i=1;
                        }
                    }
                    else{
                        if(iy%2==0){
                            spin_i=1;
                        }
                        else{
                            spin_i=0;
                        }
                    }
                    alpha_i = Coordinates_.Nbasis(ix,iy,0) + spin_i*(ncells_);
                    if(ix%2==0){
                        File_Out<< alpha_i<<"   "<<alpha_i<<"   "<<one_complex<<endl;
                    }

                }
            }
        }

        if(ansatz=="Stripe_CDW_AFM"){
            int spin_i;
            //Stripe_CDW_AFM
            File_Out<<"#alpha_i   alpha_j   OParams_[alpha_i][alpha_j]   // FOR STRIPE CDW"<<endl;
            for(int ix=0;ix<lx_;ix++){
                for(int iy=0;iy<ly_;iy++){

                    if(iy%2==0){
                        spin_i=0;
                    }
                    else{
                        spin_i=1;
                    }
                    alpha_i = Coordinates_.Nbasis(ix,iy,0) + spin_i*(ncells_);
                    if(ix%2==0){
                        File_Out<< alpha_i<<"   "<<alpha_i<<"   "<<one_complex<<endl;
                    }

                }
            }
        }

        if(ansatz=="DiagonalStripe_CDW_AFM"){
            int spin_i;
            //Stripe_CDW_AFM
            File_Out<<"#alpha_i   alpha_j   OParams_[alpha_i][alpha_j]   // FOR STRIPE CDW"<<endl;
            for(int ix=0;ix<lx_;ix++){
                for(int iy=0;iy<ly_;iy++){


                    if(ix%2==0){
                        if(iy%2!=0){
                            spin_i=0;
                            alpha_i = Coordinates_.Nbasis(ix,iy,0) + spin_i*(ncells_);
                            File_Out<< alpha_i<<"   "<<alpha_i<<"   "<<one_complex<<endl;
                        }
                    }
                    else{
                        if(iy%2==0){
                            spin_i=1;
                            alpha_i = Coordinates_.Nbasis(ix,iy,0) + spin_i*(ncells_);
                            File_Out<< alpha_i<<"   "<<alpha_i<<"   "<<one_complex<<endl;
                        }
                    }



                }
            }
        }

        if(ansatz=="Honeycomb_CDW_FM"){
            int spin_i;
            int eff_lx, eff_ly;

            if( lx_%3!=0 || ly_%3!=0){
                cout<<"lx and ly must be multiple of 3 for Honeycomb CDW"<<endl;

                //assert(lx_%3==0);
                //assert(ly_%3==0);
                if(lx_%3!=0  || ly_%3!=0 ){
                    cout<<"WARNING: lx or ly is not multiple of 3"<<endl;
                }
            }
            eff_lx=lx_/3;
            eff_ly=ly_/3;
            int ix, iy;
            File_Out<<"#alpha_i   alpha_j   OParams_[alpha_i][alpha_j]   // FOR STRIPE CDW"<<endl;
            for(int eff_ix=0;eff_ix<eff_lx;eff_ix++){
                for(int eff_iy=0;eff_iy<eff_ly;eff_iy++){

                    //9 points
                    for(int x_=0;x_<3;x_++){
                        for(int y_=0;y_<3;y_++){
                            ix = eff_ix*(3) + x_;
                            iy = eff_iy*(3) + y_;
                            if( (x_==0 && y_==0) || (x_==1 && y_==1) || (x_==2 && y_==2)  ){
                                //nothing
                            }
                            else{
                                spin_i=0;
                                alpha_i = Coordinates_.Nbasis(ix,iy,0) + spin_i*(ncells_);
                                File_Out<< alpha_i<<"   "<<alpha_i<<"   "<<one_complex<<endl;
                            }


                        }
                    }
                }
            }
        }

        if(ansatz=="Honeycomb_CDW_AFM"){
            int spin_i;
            int eff_lx, eff_ly;

            if( lx_%3!=0 || ly_%3!=0){
                cout<<"lx and ly must be multiple of 3 for Honeycomb CDW"<<endl;
                if(lx_%3!=0  || ly_%3!=0 ){
                    cout<<"WARNING: lx or ly is not multiple of 3"<<endl;
                }
            }
            eff_lx=lx_/3;
            eff_ly=ly_/3;
            int ix, iy;
            File_Out<<"#alpha_i   alpha_j   OParams_[alpha_i][alpha_j]   // FOR STRIPE CDW"<<endl;
            for(int eff_ix=0;eff_ix<eff_lx;eff_ix++){
                for(int eff_iy=0;eff_iy<eff_ly;eff_iy++){

                    //9 points
                    for(int x_=0;x_<3;x_++){
                        for(int y_=0;y_<3;y_++){
                            ix = eff_ix*(3) + x_;
                            iy = eff_iy*(3) + y_;
                            if( (x_==0 && y_==0) || (x_==1 && y_==1) || (x_==2 && y_==2)  ){
                                //nothing
                            }
                            else{
                                if( (x_==1 && y_==0)  || (x_==2 && y_==1) || (x_==0 && y_==2) ){
                                    spin_i=0;
                                }
                                else{
                                    spin_i=1;
                                }

                                alpha_i = Coordinates_.Nbasis(ix,iy,0) + spin_i*(ncells_);
                                File_Out<< alpha_i<<"   "<<alpha_i<<"   "<<one_complex<<endl;
                            }


                        }
                    }
                }
            }
        }

        if(ansatz=="size_two_eq_triangles_CDW_FM"){
            int spin_i;
            int eff_lx, eff_ly;

            if( lx_%3!=0 || ly_%3!=0){
                cout<<"lx and ly must be multiple of 3 for Honeycomb CDW"<<endl;
                if(lx_%3!=0  || ly_%3!=0 ){
                    cout<<"WARNING: lx or ly is not multiple of 3"<<endl;
                }
            }
            eff_lx=lx_/3;
            eff_ly=ly_/3;
            int ix, iy;
            File_Out<<"#alpha_i   alpha_j   OParams_[alpha_i][alpha_j]   // FOR STRIPE CDW"<<endl;
            for(int eff_ix=0;eff_ix<eff_lx;eff_ix++){
                for(int eff_iy=0;eff_iy<eff_ly;eff_iy++){

                    //9 points
                    for(int x_=0;x_<3;x_++){
                        for(int y_=0;y_<3;y_++){
                            ix = eff_ix*(3) + x_;
                            iy = eff_iy*(3) + y_;
                            if( (x_==0 && y_==0) || (x_==1 && y_==0) || (x_==0 && y_==1)  ){
                                //nothing
                            }
                            else{
                                spin_i=0;
                                alpha_i = Coordinates_.Nbasis(ix,iy,0) + spin_i*(ncells_);
                                File_Out<< alpha_i<<"   "<<alpha_i<<"   "<<one_complex<<endl;
                            }


                        }
                    }
                }
            }
        }

        if(ansatz=="doubled_unitcell_CDW_FM"){
            int spin_i;
            int eff_lx, eff_ly;

            eff_lx=lx_/2;
            eff_ly=ly_/2;
            int ix, iy;
            File_Out<<"#alpha_i   alpha_j   OParams_[alpha_i][alpha_j]   // FOR STRIPE CDW"<<endl;
            for(int eff_ix=0;eff_ix<eff_lx;eff_ix++){
                for(int eff_iy=0;eff_iy<eff_ly;eff_iy++){

                    //4 points
                    for(int x_=0;x_<2;x_++){
                        for(int y_=0;y_<2;y_++){
                            ix = eff_ix*(2) + x_;
                            iy = eff_iy*(2) + y_;
                            if( (x_==0 && y_==0) ){
                                spin_i=0;
                                alpha_i = Coordinates_.Nbasis(ix,iy,0) + spin_i*(ncells_);
                                File_Out<< alpha_i<<"   "<<alpha_i<<"   "<<one_complex<<endl;
                            }
                            else{
                                //nothing
                            }


                        }
                    }
                }
            }

        }


    }//for str_no

}

void Connections_HC::InteractionsCreate()
{

    Mat_1_int U_neighs;
    Mat_3_doub U_hoppings;
    U_neighs.push_back(0);U_neighs.push_back(2);
    U_hoppings.push_back(Parameters_.U1);U_hoppings.push_back(Parameters_.U1);


    double U_val;


    int l, m, a, b;
    int l1_pos, l2_pos;
    int m1_pos, m2_pos;

    Hint_.fill(0.0);

    if(!Parameters_.LongRange_interaction){

        //Intra-site interactions
        for(l=0;l<nsites_;l++){

            l1_pos = Coordinates_.indx_cellwise(l);
            l2_pos = Coordinates_.indy_cellwise(l);
            int orb0=0;
            int orb1=1;

            a = Coordinates_.Nbasis(l1_pos, l2_pos, orb0);
            b = Coordinates_.Nbasis(l1_pos, l2_pos, orb1);

            if(a!=b){
                Hint_(a,b) = complex<double>(1.0 * Parameters_.U0_interorb, 0.0);
                Hint_(b,a) = Hint_(a,b);
            }

        }


        //inter-sites
        for (l = 0; l < nsites_; l++)
        {
            l1_pos = Coordinates_.indx_cellwise(l);
            l2_pos = Coordinates_.indy_cellwise(l);

            //For U1, U2, U3
            for(int neigh=0;neigh<U_neighs.size();neigh++){
                m = Coordinates_.getneigh(l, U_neighs[neigh]); //+x neighbour cell

                if( (!Coordinates_.HIT_X_BC) || Parameters_.PBC_X){

                    if( (!Coordinates_.HIT_Y_BC) || Parameters_.PBC_Y){

                        m1_pos = Coordinates_.indx_cellwise(m);
                        m2_pos = Coordinates_.indy_cellwise(m);

                        for (int orb1 = 0; orb1 < n_orbs_; orb1++)
                        {
                            for (int orb2 = 0; orb2 < n_orbs_; orb2++)
                            {
                                a = Coordinates_.Nbasis(l1_pos, l2_pos, orb1);
                                b = Coordinates_.Nbasis(m1_pos, m2_pos, orb2);

                                assert(a != b);
                                if (a != b)
                                {
                                    Hint_(a, b) = complex<double>(1.0 * U_hoppings[neigh][orb2][orb1], 0.0);
                                    Hint_(b, a) = conj(Hint_(a, b));
                                }


                            }
                        }
                    }

                }

            }

        }



    }
    else{
        cout <<"WARNING: ALWAYS ASSUMED PBC in Long range int"<<endl;
        int r1_, r2_;
        double rx_, ry_;
        double dis_;
        double dis_min;
        Get_minimum_distance_direction(0 ,0, 0, 1, r1_,r2_, dis_min);
        for(l=0; l < nsites_; l++)
        {
            l1_pos = Coordinates_.indx_cellwise(l);
            l2_pos = Coordinates_.indy_cellwise(l);
            for(int gamma_l=0;gamma_l<n_orbs_;gamma_l++){

                //For U[l][m]
                for(m=0;m<nsites_;m++){
                    m1_pos = Coordinates_.indx_cellwise(m);
                    m2_pos = Coordinates_.indy_cellwise(m);

                    for(int gamma_m=0;gamma_m<n_orbs_;gamma_m++){

                        a = Coordinates_.Nbasis(l1_pos, l2_pos, gamma_l);
                        b = Coordinates_.Nbasis(m1_pos, m2_pos, gamma_m);

                        Get_minimum_distance_direction(l,gamma_l, m,gamma_m, r1_,r2_, dis_);

                        //a1 = (1,0) in (x,y)
                        //a2 = (1/2, (sqrt(3)/2)/2)
                        //                rx_ = ((1.0)*(r1_) +  (1.0/2.0)*(r2_));
                        //                ry_ =  (0.0*(r1_) + (sqrt(3.0)/2.0)*(r2_));
                        //                dis_= sqrt(rx_*rx_ + ry_*ry_);

                        //cout <<l<<"  "<<m<<"  "<<r1_<<"   "<<r2_<<"   "<<dis_<<endl;
                        //(14.3952)*( (1.0/(x*62.6434))  - (1.0/(sqrt( (x*x*62.6434*62.6434)  + d*d)))  )

                        // U_val = ((14.3952*1000)/Parameters_.eps_DE)*( (1.0/(dis_*Parameters_.a_moire))  -
                        //                                               (1.0/(sqrt( (dis_*dis_*Parameters_.a_moire*Parameters_.a_moire)
                        //                                                           + Parameters_.d_screening*Parameters_.d_screening)))  );


                        if(dis_<=Parameters_.Truncating_Length_in_am){
                            U_val = Parameters_.U0_interorb*( (1.0/(dis_*Parameters_.a_moire))  -
                                                              (1.0/(sqrt( (dis_*dis_*Parameters_.a_moire*Parameters_.a_moire)
                                                                          + Parameters_.d_screening*Parameters_.d_screening)))  )*
                                    (1.0 /   (  ( (1.0/(dis_min*Parameters_.a_moire))  -
                                                  (1.0/(sqrt( (dis_min*dis_min*Parameters_.a_moire*Parameters_.a_moire)
                                                              + Parameters_.d_screening*Parameters_.d_screening)))  )  ));
                        }
                        else{
                            U_val=0.0;
                        }


                        // assert(l != m);
                        if (  (l!=m) || (gamma_l!=gamma_m)   )
                        {
                            //cout<<"NOT WORKING"<<endl;
                            //if(dis_<=1.01){
                            Hint_(a, b) = complex<double>(1.0 * U_val, 0.0);
                            //}
                            //Hint_(b, a) = conj(Hint_(a, b));
                        }

                    }

                }

            }
        }
    }


} // ----------



void Connections_HC::Interactions_Sorting()
{

    Mat_1_int U_neighs;
    Mat_3_doub U_hoppings;
    U_neighs.push_back(0);U_neighs.push_back(2);
    U_hoppings.push_back(Parameters_.U1);U_hoppings.push_back(Parameters_.U1);


    double U_val;

    int l, m, a, b;
    int l1_pos, l2_pos;
    int m1_pos, m2_pos;

    Mat_1_doub Int_val;
    Mat_1_doub Dis_val;
    Mat_1_int SiteInd_val;

    if(Parameters_.LongRange_interaction){
        cout <<"WARNING: ALWAYS ASSUMED PBC in Long range int, Sorting is performed for plotting U(r) vs r"<<endl;
        int r1_, r2_;
        double rx_, ry_;
        double dis_;
        double dis_min;
        Get_minimum_distance_direction(0 ,0, 0, 1, r1_,r2_, dis_min);
        l=0;
        l1_pos = Coordinates_.indx_cellwise(l);
        l2_pos = Coordinates_.indy_cellwise(l);
        int gamma_l=0;

        //For U[l][m]
        for(m=0;m<nsites_;m++){
            m1_pos = Coordinates_.indx_cellwise(m);
            m2_pos = Coordinates_.indy_cellwise(m);

            for(int gamma_m=0;gamma_m<n_orbs_;gamma_m++){

                a = Coordinates_.Nbasis(l1_pos, l2_pos, gamma_l);
                b = Coordinates_.Nbasis(m1_pos, m2_pos, gamma_m);

                Get_minimum_distance_direction(l,gamma_l, m,gamma_m, r1_,r2_, dis_);

                //a1 = (1,0) in (x,y)
                //a2 = (1/2, (sqrt(3)/2)/2)
                //                rx_ = ((1.0)*(r1_) +  (1.0/2.0)*(r2_));
                //                ry_ =  (0.0*(r1_) + (sqrt(3.0)/2.0)*(r2_));
                //                dis_= sqrt(rx_*rx_ + ry_*ry_);

                //cout <<l<<"  "<<m<<"  "<<r1_<<"   "<<r2_<<"   "<<dis_<<endl;
                //(14.3952)*( (1.0/(x*62.6434))  - (1.0/(sqrt( (x*x*62.6434*62.6434)  + d*d)))  )

                // U_val = ((14.3952*1000)/Parameters_.eps_DE)*( (1.0/(dis_*Parameters_.a_moire))  -
                //                                               (1.0/(sqrt( (dis_*dis_*Parameters_.a_moire*Parameters_.a_moire)
                //                                                           + Parameters_.d_screening*Parameters_.d_screening)))  );


                if(dis_<=Parameters_.Truncating_Length_in_am){

                    U_val = Parameters_.U0_interorb*( (1.0/(dis_*Parameters_.a_moire))  -
                                                      (1.0/(sqrt( (dis_*dis_*Parameters_.a_moire*Parameters_.a_moire)
                                                                  + Parameters_.d_screening*Parameters_.d_screening)))  )*
                            (1.0 /   (  ( (1.0/(dis_min*Parameters_.a_moire))  -
                                          (1.0/(sqrt( (dis_min*dis_min*Parameters_.a_moire*Parameters_.a_moire)
                                                      + Parameters_.d_screening*Parameters_.d_screening)))  )  ));
                }
                else{
                    U_val=0.0;
                }


                // assert(l != m);
                if (  (l!=m) || (gamma_l!=gamma_m)   )
                {

                    Int_val.push_back(U_val);
                    SiteInd_val.push_back(b);
                    Dis_val.push_back(dis_);

                }

            }

        }
    }




    //Sorting
    double large_val=100000;
    int i_min;
    double min_val;

    Mat_1_doub Dis_val_sorted;
    Mat_1_doub Int_val_sorted;
    Mat_1_int SiteInd_val_sorted;

    for(int j=0;j<Dis_val.size();j++){
        min_val=large_val;
        for(int i=0;i<Dis_val.size();i++){

            if(Dis_val[i]<=min_val){
                min_val=Dis_val[i];
                i_min=i;
            }
        }
        Dis_val_sorted.push_back(Dis_val[i_min]);
        Int_val_sorted.push_back(Int_val[i_min]);
        SiteInd_val_sorted.push_back(SiteInd_val[i_min]);

        Dis_val[i_min]=large_val;

    }


    Dis_val.clear();
    Int_val.clear();
    SiteInd_val.clear();

    Mat_1_int Degeneracy;
    int deg_temp;

    Dis_val.push_back(Dis_val_sorted[0]);
    Int_val.push_back(Int_val_sorted[0]);
    Degeneracy.resize(1);

    deg_temp=1;
    for(int i=1;i<Dis_val_sorted.size();i++){

        if(abs(Int_val[Int_val.size()-1] - Int_val_sorted[i])>0.0000001){
        Dis_val.push_back(Dis_val_sorted[i]);
        Int_val.push_back(Int_val_sorted[i]);
        Degeneracy[Degeneracy.size()-1]=deg_temp;
        Degeneracy.resize(Degeneracy.size()+1);
        deg_temp=1;
        }
        else{
            deg_temp++;
        }
    }

    Degeneracy[Degeneracy.size()-1]=deg_temp;

    string fileout="Sorted_Ur_" + Parameters_.File_LongRange_Ints;
    ofstream file_out(fileout.c_str());
    file_out<<"#i r(i) U(r)"<<endl;
    for(int i=0;i<Dis_val_sorted.size();i++){
        file_out<<SiteInd_val_sorted[i]<<"   "<<Dis_val_sorted[i]<<"   "<<Int_val_sorted[i]<<endl;
    }


    string fileout2="Sorted2_Ur_" + Parameters_.File_LongRange_Ints;
    ofstream file_out2(fileout2.c_str());
    file_out2<<"#i    degeneracy(i)   r(i) U(r)"<<endl;
    file_out2<<0<<"   "<< 1<<"     "<<0<<"   "<<Parameters_.U0<<endl;

    for(int i=0;i<Dis_val.size();i++){
        file_out2<<i+1<<"   "<< Degeneracy[i]<<"     "<<Dis_val[i]<<"   "<<Int_val[i]<<endl;
    }




} // ----------


void Connections_HC::Check_up_down_symmetry()

{
    complex<double> temp(0, 0);
    complex<double> temp2;

    for (int i = 0; i < ncells_ * n_orbs_; i++)
    {
        for (int j = 0; j < ncells_ * n_orbs_; j++)
        {
            temp2 = Ham_(i, j) - Ham_(i + ncells_ * n_orbs_, j + ncells_ * n_orbs_); //+ Ham_(i+orbs_*ns_,j) + Ham_(i,j+orbs_*ns_);
            temp += temp2 * conj(temp2);
        }
    }

    cout << "Assymetry in up-down sector: " << temp << endl;
}

void Connections_HC::Check_Hermiticity()

{


    for (int i = 0; i < Ham_.n_row(); i++)
    {
        for (int j = 0; j < Ham_.n_row(); j++)
        {
            if (Ham_(i, j) != conj(Ham_(j, i)))
            {
                cout << i << "," << j << endl;
                cout << "i,j = " << Ham_(i, j) << ", j,i=" << conj(Ham_(j, i)) << endl;
            }
            assert(Ham_(i, j) == conj(Ham_(j, i))); //+ Ham_(i+orbs_*ns_,j) + Ham_(i,j+orbs_*ns_);
            //temp +=temp2*conj(temp2);
        }
    }

    // cout<<"Hermiticity: "<<temp<<endl;
}


void Connections_HC::HTBCreate()
{

    //Convention used
    //orb=0=Blue  //See Cartoon in notes
    //orb=1=Red

    Mat_1_int t_neighs;
    Mat_3_Complex_doub t_hoppings;
    //t1 hoppings
    t_neighs.push_back(0);t_neighs.push_back(2);
    //t_neighs.push_back(5);
    t_hoppings.push_back(Parameters_.t1);t_hoppings.push_back(Parameters_.t1);
    //t_hoppings.push_back(Parameters_.t1);

    t_neighs.push_back(0);t_neighs.push_back(2);t_neighs.push_back(5);
    t_hoppings.push_back(Parameters_.t2);t_hoppings.push_back(Parameters_.t2);t_hoppings.push_back(Parameters_.t2);

    //    t_neighs.push_back(10);t_neighs.push_back(15);t_neighs.push_back(16);
    //    t_hoppings.push_back(Parameters_.t3);t_hoppings.push_back(Parameters_.t3);t_hoppings.push_back(Parameters_.t3);




    Matrix<complex <double>> sigma_x, sigma_y, sigma_z;
    sigma_x.resize(2, 2);
    sigma_y.resize(2, 2);
    sigma_z.resize(2, 2);

    // X
    sigma_x(0, 0) = 0.0;
    sigma_x(0, 1) = 1.0;
    sigma_x(1, 0) = 1.0;
    sigma_x(1, 1) = 0.0;

    // y
    sigma_y(0, 0) = 0.0;
    sigma_y(0, 1) = -1.0*iota_complex;
    sigma_y(1, 0) = iota_complex;
    sigma_y(1, 1) = 0.0;

    // Z
    sigma_z(0, 0) = 1.0;
    sigma_z(0, 1) = 0.0;
    sigma_z(1, 0) = 0.0;
    sigma_z(1, 1) = -1.0;


    int l, m, a, b;
    int lx_pos, ly_pos;
    int mx_pos, my_pos;

    HTB_.fill(0.0);

    for (l = 0; l < nsites_; l++)
    {
        lx_pos = Coordinates_.indx_cellwise(l);
        ly_pos = Coordinates_.indy_cellwise(l);


        //Intra-site hoppings
        for (int spin = 0; spin < 2; spin++)
        {
            int orb1 = 0;
            int orb2 = 1;

            a = Coordinates_.Nbasis(lx_pos, ly_pos, orb1) + nsites_ * n_orbs_ * spin;
            b = Coordinates_.Nbasis(lx_pos, ly_pos, orb2) + nsites_ * n_orbs_ * spin;
            assert(a != b);
            if (a != b)
            {
                if(spin==0){
                    HTB_(b, a) = -1.0*Parameters_.t0;
                    HTB_(a, b) = conj(HTB_(b, a));
                }
                else{
                    HTB_(b, a) = -1.0*Parameters_.t0;
                    HTB_(a, b) = conj(HTB_(b, a));
                }
            }

        }





        //For t1,t2,t3 hoppings
        for(int neigh=0;neigh<t_neighs.size();neigh++){
            m = Coordinates_.getneigh(l, t_neighs[neigh]);

            if( (!Coordinates_.HIT_X_BC) || Parameters_.PBC_X){

                if( (!Coordinates_.HIT_Y_BC) || Parameters_.PBC_Y){

                    mx_pos = Coordinates_.indx_cellwise(m);
                    my_pos = Coordinates_.indy_cellwise(m);

                    for (int spin = 0; spin < 2; spin++)
                    {
                        for (int orb1 = 0; orb1 < n_orbs_; orb1++)
                        {
                            for (int orb2 = 0; orb2 < n_orbs_; orb2++)
                            {
                                a = Coordinates_.Nbasis(lx_pos, ly_pos, orb1) + nsites_ * n_orbs_ * spin;
                                b = Coordinates_.Nbasis(mx_pos, my_pos, orb2) + nsites_ * n_orbs_ * spin;
                                assert(a != b);
                                if (a != b)
                                {
                                    if(spin==0){
                                        HTB_(b, a) += -1.0*t_hoppings[neigh][orb2][orb1];
                                        HTB_(a, b) = conj(HTB_(b, a));
                                    }
                                    else{
                                        HTB_(b, a) += -1.0*conj(t_hoppings[neigh][orb2][orb1]);
                                        HTB_(a, b) = conj(HTB_(b, a));
                                    }
                                }

                            }
                        }
                    }
                }

            }

        }

    }


} // ----------


void Connections_HC::Print_Hopping(){

    string fileout=Parameters_.File_Hoppings;
    ofstream file_Hopping_out(fileout.c_str());

    file_Hopping_out<<"#site_i spin_i site_j spin_j Hopping[site_i,spin_i][site_j,spin_j]"<<endl;

    int index_i, index_j;
    for(int i=0;i<nsites_*n_orbs_;i++){
        for(int spin_i=0;spin_i<2;spin_i++){
            index_i=i+spin_i*nsites_*n_orbs_;

            for(int j=0;j<nsites_*n_orbs_;j++){
                for(int spin_j=0;spin_j<2;spin_j++){
                    index_j=j+spin_j*nsites_*n_orbs_;

                    if(abs(HTB_(index_i,index_j))>0.0000001){
                        file_Hopping_out<<i<<"  "<<spin_i<<"  "<<j<<"  "<<spin_j<<"  "<<HTB_(index_i,index_j)<<endl;
                    }
                }
            }
        }
    }

}



void Connections_HC::Print_LongRangeInt(){
    string fileout=Parameters_.File_LongRange_Ints;
    ofstream file_out(fileout.c_str());

    file_out<<"#i j Uij"<<endl;

    for(int i=0;i<nsites_*n_orbs_;i++){

        for(int j=0;j<nsites_*n_orbs_;j++){

            if(abs(Hint_(i,j))>0.0000001){
                file_out<<i<<"   "<<j<<"   "<<Hint_(i,j)<<endl;
            }

        }
    }

}

void Connections_HC::Print_Spin_resolved_OnsiteE(){

    string fileout=Parameters_.File_Onsite_Energies;
    ofstream file_out(fileout.c_str());
    file_out<<"#Site Spin Onsite_E"<<endl;
    int index, i_posx, i_posy;

    for (int i = 0; i < nsites_; i++)
    { // For each cell
        i_posx = Coordinates_.indx_cellwise(i);
        i_posy = Coordinates_.indy_cellwise(i);

        for (int orb = 0; orb < n_orbs_; orb++)
        {
            index = Coordinates_.Nbasis(i_posx, i_posy, orb);

            // On-Site potential
            for (int spin = 0; spin < 2; spin++)
            {
                file_out<<index<<"  "<<spin<<"  "<<Parameters_.OnSiteE[orb][spin]<<endl;
            }
        }
    }



}


#endif
