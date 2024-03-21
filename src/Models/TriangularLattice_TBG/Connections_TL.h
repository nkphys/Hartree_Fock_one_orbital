#include <algorithm>
#include <functional>
#include <math.h>
#include "tensor_type.h"
#include "Parameters_TL.h"
#include "Coordinates_TL.h"
#include "../../Matrix.h"
#define PI acos(-1.0)

#ifndef Connections_TL_class
#define Connections_TL_class


class Connections_TL
{
public:
    Connections_TL(Parameters_TL &Parameters__, Coordinates_TL &Coordinates__)
        : Parameters_(Parameters__), Coordinates_(Coordinates__)

    {
        Initialize();
        HTBCreate();
    }

    void Initialize();                                     //::DONE
    void Print_Hopping();                                       //::DONE
    void Print_Hopping2();                                       //::DONE
    void Print_Hopping3();                                       //::DONE
    void Print_Hopping_Snake3();
    void Print_Hopping_YCn();
    void Print_LongRangeInt_SNAKE3();
    void Print_RingExchange();
    void Print_RingExchange_YCn();

    void Print_LongRangeInt();
    void Print_Spin_resolved_OnsiteE();
    void InteractionsCreate();                             //::DONE
    void Check_Hermiticity();                              //::DONE
    void Check_up_down_symmetry();                         //::DONE
    void HTBCreate();                                      //::DONE
    void Print_Ansatz_LocalDen_CDW();
    void Get_minimum_distance_direction(int l,int m,int &r1_, int &r2_);

    int Convert_to_snake3_coordinates_YCn(int i);


    Parameters_TL &Parameters_;
    Coordinates_TL &Coordinates_;
    int lx_, ly_, ncells_, n_orbs_;
    Matrix<complex<double>> HTB_;
    Matrix<complex<double>> Hint_;
    Matrix<complex<double>> Ham_;
    Matrix<double> Tx, Ty, Tpxpy, Tpxmy;
    Matrix<double> IntraCell_Hopp, InterCell_px_Hopp, InterCell_py_Hopp, InterCell_pxmy_Hopp;

    int max_l1_dis, max_l2_dis;
};


void Connections_TL::Get_minimum_distance_direction(int l,int m,int &r1_, int &r2_){

    int r1_a,r1_b, r2_a, r2_b;
    int l1_pos, l2_pos;
    int m1_pos, m2_pos;
    double rx_, ry_;

    l1_pos = Coordinates_.indx_cellwise(l);
    l2_pos = Coordinates_.indy_cellwise(l);
    m1_pos = Coordinates_.indx_cellwise(m);
    m2_pos = Coordinates_.indy_cellwise(m);

    r1_a=m1_pos-l1_pos;
    if(r1_a>0){
        r1_b=-1*(lx_-r1_a);
    }
    else if (r1_a<0){
        r1_b=lx_-abs(r1_a);
    }
    else{
        r1_b=0;
        assert(r1_a==0);
    }

    r2_a=m2_pos-l2_pos;
    if(r2_a>0){
        r2_b=-1*(ly_-r2_a);
    }
    else if (r2_a<0){
        r2_b=ly_-abs(r2_a);
    }
    else{
        r2_b=0;
        assert(r2_a==0);
    }


    double min_dis=1000000.0;
    double dis;

    //r1a r2a
    rx_ = ((1.0)*(r1_a) +  (1.0/2.0)*(r2_a));
    ry_ =  (0.0*(r1_a) + (sqrt(3.0)/2.0)*(r2_a));
    dis= sqrt(rx_*rx_ + ry_*ry_);
    if(dis<=min_dis){
        r1_=r1_a;
        r2_=r2_a;
        min_dis=dis;
    }

    //r1a r2b
    rx_ = ((1.0)*(r1_a) +  (1.0/2.0)*(r2_b));
    ry_ =  (0.0*(r1_a) + (sqrt(3.0)/2.0)*(r2_b));
    dis= sqrt(rx_*rx_ + ry_*ry_);
    if(dis<=min_dis){
        r1_=r1_a;
        r2_=r2_b;
        min_dis=dis;
    }

    //r1b r2a
    rx_ = ((1.0)*(r1_b) +  (1.0/2.0)*(r2_a));
    ry_ =  (0.0*(r1_b) + (sqrt(3.0)/2.0)*(r2_a));
    dis= sqrt(rx_*rx_ + ry_*ry_);
    if(dis<=min_dis){
        r1_=r1_b;
        r2_=r2_a;
        min_dis=dis;
    }

    //r1b r2b
    rx_ = ((1.0)*(r1_b) +  (1.0/2.0)*(r2_b));
    ry_ =  (0.0*(r1_b) + (sqrt(3.0)/2.0)*(r2_b));
    dis= sqrt(rx_*rx_ + ry_*ry_);
    if(dis<=min_dis){
        r1_=r1_b;
        r2_=r2_b;
        min_dis=dis;
    }

}

void Connections_TL::Initialize()
{


    ly_ = Parameters_.ly;
    lx_ = Parameters_.lx;
    ncells_ = lx_ * ly_;
    n_orbs_ = Parameters_.n_orbs;
    int space = 2 * ncells_ * n_orbs_;

    HTB_.resize(space, space);
    Ham_.resize(space, space);
    Hint_.resize(ncells_,ncells_);

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


void Connections_TL::Print_Ansatz_LocalDen_CDW(){



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

void Connections_TL::InteractionsCreate()
{

    Mat_1_int U_neighs;
    Mat_1_doub U_hoppings;
    U_neighs.push_back(0);U_neighs.push_back(3);U_neighs.push_back(5);
    U_hoppings.push_back(Parameters_.U1);U_hoppings.push_back(Parameters_.U1);U_hoppings.push_back(Parameters_.U1);

    U_neighs.push_back(4);U_neighs.push_back(13);U_neighs.push_back(14);
    U_hoppings.push_back(Parameters_.U2);U_hoppings.push_back(Parameters_.U2);U_hoppings.push_back(Parameters_.U2);

    U_neighs.push_back(10);U_neighs.push_back(15);U_neighs.push_back(16);
    U_hoppings.push_back(Parameters_.U3);U_hoppings.push_back(Parameters_.U3);U_hoppings.push_back(Parameters_.U3);


    double U_val;




    int l, m, a, b;
    int l1_pos, l2_pos;
    int m1_pos, m2_pos;

    Hint_.fill(0.0);

    if(!Parameters_.LongRange_interaction){
        for (l = 0; l < ncells_; l++)
        {
            l1_pos = Coordinates_.indx_cellwise(l);
            l2_pos = Coordinates_.indy_cellwise(l);

            //For U1, U2, U3
            for(int neigh=0;neigh<9;neigh++){
                m = Coordinates_.getneigh(l, U_neighs[neigh]); //+x neighbour cell

                if( (!Coordinates_.HIT_X_BC) || Parameters_.PBC_X){

                    if( (!Coordinates_.HIT_Y_BC) || Parameters_.PBC_Y){

                        m1_pos = Coordinates_.indx_cellwise(m);
                        m2_pos = Coordinates_.indy_cellwise(m);

                        assert(l != m);
                        if (l != m)
                        {
                            Hint_(l, m) = complex<double>(1.0 * U_hoppings[neigh], 0.0);
                            Hint_(m, l) = conj(Hint_(l, m));
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
        for(l=0; l < ncells_; l++)
        {
            l1_pos = Coordinates_.indx_cellwise(l);
            l2_pos = Coordinates_.indy_cellwise(l);

            //For U[l][m]
            for(m=0;m<ncells_;m++){
                m1_pos = Coordinates_.indx_cellwise(m);
                m2_pos = Coordinates_.indy_cellwise(m);

                Get_minimum_distance_direction(l,m,r1_,r2_);

                //a1 = (1,0) in (x,y)
                //a2 = (1/2, (sqrt(3)/2)/2)
                rx_ = ((1.0)*(r1_) +  (1.0/2.0)*(r2_));
                ry_ =  (0.0*(r1_) + (sqrt(3.0)/2.0)*(r2_));
                dis_= sqrt(rx_*rx_ + ry_*ry_);

                //cout <<l<<"  "<<m<<"  "<<r1_<<"   "<<r2_<<"   "<<dis_<<endl;
                //(14.3952)*( (1.0/(x*62.6434))  - (1.0/(sqrt( (x*x*62.6434*62.6434)  + d*d)))  )

                U_val = ((14.3952*1000)/Parameters_.eps_DE)*( (1.0/(dis_*Parameters_.a_moire))  -
                                                              (1.0/(sqrt( (dis_*dis_*Parameters_.a_moire*Parameters_.a_moire)
                                                                          + Parameters_.d_screening*Parameters_.d_screening)))  );
                // assert(l != m);
                if (l != m)
                {
                    Hint_(l, m) = complex<double>(1.0 * U_val, 0.0);
                    //Hint_(m, l) += conj(Hint_(l, m));
                }

            }

        }

    }


} // ----------


void Connections_TL::Check_up_down_symmetry()

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

void Connections_TL::Check_Hermiticity()

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


void Connections_TL::HTBCreate()
{

    //Convention used
    //orb=0=A
    //orb=1=B
    //orb=2=C

    Mat_1_int t_neighs;
    Mat_1_Complex_doub t_hoppings;
    //t1 hoppings
    t_neighs.push_back(0);t_neighs.push_back(3);t_neighs.push_back(5);
    t_hoppings.push_back(Parameters_.t1);t_hoppings.push_back(Parameters_.t1);t_hoppings.push_back(Parameters_.t1);

    t_neighs.push_back(4);t_neighs.push_back(13);t_neighs.push_back(14);
    t_hoppings.push_back(Parameters_.t2);t_hoppings.push_back(Parameters_.t2);t_hoppings.push_back(Parameters_.t2);

    t_neighs.push_back(10);t_neighs.push_back(15);t_neighs.push_back(16);
    t_hoppings.push_back(Parameters_.t3);t_hoppings.push_back(Parameters_.t3);t_hoppings.push_back(Parameters_.t3);




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

    for (l = 0; l < ncells_; l++)
    {
        lx_pos = Coordinates_.indx_cellwise(l);
        ly_pos = Coordinates_.indy_cellwise(l);

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
                                a = Coordinates_.Nbasis(lx_pos, ly_pos, orb1) + ncells_ * n_orbs_ * spin;
                                b = Coordinates_.Nbasis(mx_pos, my_pos, orb2) + ncells_ * n_orbs_ * spin;
                                assert(a != b);
                                if (a != b)
                                {
                                    if(spin==0){
                                        HTB_(b, a) = 1.0*t_hoppings[neigh];
                                        HTB_(a, b) = conj(HTB_(b, a));
                                    }
                                    else{
                                        HTB_(b, a) = 1.0*conj(t_hoppings[neigh]);
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


void Connections_TL::Print_Hopping(){

    string fileout=Parameters_.File_Hoppings;
    ofstream file_Hopping_out(fileout.c_str());

    file_Hopping_out<<"#site_i spin_i site_j spin_j Hopping[site_i,spin_i][site_j,spin_j]"<<endl;

    int index_i, index_j;
    for(int i=0;i<ncells_*n_orbs_;i++){
        for(int spin_i=0;spin_i<2;spin_i++){
            index_i=i+spin_i*ncells_*n_orbs_;

            for(int j=0;j<ncells_*n_orbs_;j++){
                for(int spin_j=0;spin_j<2;spin_j++){
                    index_j=j+spin_j*ncells_*n_orbs_;

                    if(abs(HTB_(index_i,index_j))>0.0000001){
                        file_Hopping_out<<i<<"  "<<spin_i<<"  "<<j<<"  "<<spin_j<<"  "<<HTB_(index_i,index_j)<<endl;
                    }
                }
            }
        }
    }

}


void Connections_TL::Print_Hopping2(){

    string fileout= "XCn_ZZ_REALMatrixUPPERTRIANGLE_form_SPINUP" + Parameters_.File_Hoppings;
    ofstream file_Hopping_out(fileout.c_str());

    //file_Hopping_out<<"#site_i spin_i site_j spin_j Hopping[site_i,spin_i][site_j,spin_j]"<<endl;

    int spin_i=0;
    int index_i, index_j;
    for(int i=0;i<ncells_*n_orbs_;i++){
        index_i=i+spin_i*ncells_*n_orbs_;

        for(int j=0;j<ncells_*n_orbs_;j++){
            index_j=j+spin_i*ncells_*n_orbs_;
            if(index_i<index_j){
                file_Hopping_out<<HTB_(index_i,index_j).real()<<"   ";
            }
            else{
                file_Hopping_out<<0<<"   ";
            }
        }
        file_Hopping_out<<endl;
    }


    string fileout2= "REALMatrix_form_SPINUP" + Parameters_.File_Hoppings;
    ofstream file_Hopping_out2(fileout2.c_str());

    //file_Hopping_out<<"#site_i spin_i site_j spin_j Hopping[site_i,spin_i][site_j,spin_j]"<<endl;
    for(int i=0;i<ncells_*n_orbs_;i++){
        index_i=i+spin_i*ncells_*n_orbs_;

        for(int j=0;j<ncells_*n_orbs_;j++){
            index_j=j+spin_i*ncells_*n_orbs_;
            file_Hopping_out2<<HTB_(index_i,index_j).real()<<"   ";

        }
        file_Hopping_out2<<endl;
    }



}


void Connections_TL::Print_RingExchange(){



    int ix,iy, i1, i2, i3;
    string fileout2= "Ring_Exchange_" + Parameters_.File_Hoppings;
    ofstream file_Hopping_out2(fileout2.c_str());
    bool Hit_boundary_x, Hit_boundary_y;

    //file_Hopping_out<<"#site_i spin_i site_j spin_j Hopping[site_i,spin_i][site_j,spin_j]"<<endl;
    for(int i=0;i<ncells_;i++){
        ix=Coordinates_.indx_cellwise(i);
        iy=Coordinates_.indy_cellwise(i);


        //ring-1
        i1=Coordinates_.getneigh(i,0);
        Hit_boundary_x=Coordinates_.HIT_X_BC;
        Hit_boundary_y=Coordinates_.HIT_Y_BC;

        i2=Coordinates_.getneigh(i,8);
        Hit_boundary_x= (Hit_boundary_x || Coordinates_.HIT_X_BC);
        Hit_boundary_y= (Hit_boundary_y || Coordinates_.HIT_Y_BC);

        i3=Coordinates_.getneigh(i,7);
        Hit_boundary_x= (Hit_boundary_x || Coordinates_.HIT_X_BC);
        Hit_boundary_y= (Hit_boundary_y || Coordinates_.HIT_Y_BC);

        if( (!(Hit_boundary_x) || Parameters_.PBC_X)
                &&
                (!(Hit_boundary_y) || Parameters_.PBC_Y)
                ){
            file_Hopping_out2<<i<<" "<<i1<<" "<<i2<<" "<<i3<<"   VAL1"<<endl;
            file_Hopping_out2<<i<<" "<<i3<<" "<<i1<<" "<<i2<<"   VAL1"<<endl;
            file_Hopping_out2<<i<<" "<<i2<<" "<<i1<<" "<<i3<<"   VAL2"<<endl;


        }


        //ring-2
        i1=Coordinates_.getneigh(i,2);
        Hit_boundary_x=Coordinates_.HIT_X_BC;
        Hit_boundary_y=Coordinates_.HIT_Y_BC;

        i2=Coordinates_.getneigh(i,4);
        Hit_boundary_x= (Hit_boundary_x || Coordinates_.HIT_X_BC);
        Hit_boundary_y= (Hit_boundary_y || Coordinates_.HIT_Y_BC);

        i3=Coordinates_.getneigh(i,0);
        Hit_boundary_x= (Hit_boundary_x || Coordinates_.HIT_X_BC);
        Hit_boundary_y= (Hit_boundary_y || Coordinates_.HIT_Y_BC);

        if( (!(Hit_boundary_x) || Parameters_.PBC_X)
                &&
                (!(Hit_boundary_y) || Parameters_.PBC_Y)
                ){
            file_Hopping_out2<<i<<" "<<i1<<" "<<i2<<" "<<i3<<"   VAL1"<<endl;
            file_Hopping_out2<<i<<" "<<i3<<" "<<i1<<" "<<i2<<"   VAL1"<<endl;
            file_Hopping_out2<<i<<" "<<i2<<" "<<i1<<" "<<i3<<"   VAL2"<<endl;

        }



        //ring-3
        i1=Coordinates_.getneigh(i,5);
        Hit_boundary_x=Coordinates_.HIT_X_BC;
        Hit_boundary_y=Coordinates_.HIT_Y_BC;

        i2=Coordinates_.getneigh(i,9);
        Hit_boundary_x= (Hit_boundary_x || Coordinates_.HIT_X_BC);
        Hit_boundary_y= (Hit_boundary_y || Coordinates_.HIT_Y_BC);

        i3=Coordinates_.getneigh(i,2);
        Hit_boundary_x= (Hit_boundary_x || Coordinates_.HIT_X_BC);
        Hit_boundary_y= (Hit_boundary_y || Coordinates_.HIT_Y_BC);

        if( (!(Hit_boundary_x) || Parameters_.PBC_X)
                &&
                (!(Hit_boundary_y) || Parameters_.PBC_Y)
                ){
            file_Hopping_out2<<i<<" "<<i1<<" "<<i2<<" "<<i3<<"   VAL1"<<endl;
            file_Hopping_out2<<i<<" "<<i3<<" "<<i1<<" "<<i2<<"   VAL1"<<endl;
            file_Hopping_out2<<i<<" "<<i2<<" "<<i1<<" "<<i3<<"   VAL2"<<endl;

        }



    }



}



void Connections_TL::Print_Hopping3(){

    // assert (Parameters_.PBC_Y==false);
    //assert (Parameters_.PBC_Y==true);

    string fileout= "REALMatrix_Snake3_form_SPINUP" + Parameters_.File_Hoppings;
    ofstream file_Hopping_out(fileout.c_str());

    //file_Hopping_out<<"#site_i spin_i site_j spin_j Hopping[site_i,spin_i][site_j,spin_j]"<<endl;

    int spin_i=0;
    int index_i, index_j;
    int ix, iy, ix_new, iy_new, i_new;
    int jx, jy, jx_new, jy_new, j_new;
    for(int i=0;i<ncells_*n_orbs_;i++){
        ix = Coordinates_.indx_cellwise(i);
        iy = Coordinates_.indy_cellwise(i);
        if(iy%2==0){
            ix_new=ix;
        }
        else{
            ix_new = (lx_-1 - ix);
        }
        i_new = Coordinates_.Ncell(ix_new, iy);
        index_i=i_new+spin_i*ncells_*n_orbs_;

        for(int j=0;j<ncells_*n_orbs_;j++){
            jx = Coordinates_.indx_cellwise(j);
            jy = Coordinates_.indy_cellwise(j);

            if(jy%2==0){
                jx_new=jx;
            }
            else{
                jx_new = (lx_-1 - jx);
            }

            j_new = Coordinates_.Ncell(jx_new, jy);
            index_j=j_new+spin_i*ncells_*n_orbs_;

            if(j>=i){
                file_Hopping_out<<HTB_(index_i,index_j).real()<<"   ";}
            else{
                file_Hopping_out<<0.0<<"   ";
            }
        }
        file_Hopping_out<<endl;
    }

}



int Connections_TL::Convert_to_snake3_coordinates_YCn(int i){

    int ix_new, iy_new;
    int i_new;
    int ix, iy;
    ix = Coordinates_.indx_cellwise(i);
    iy = Coordinates_.indy_cellwise(i);

    if(iy%2==0){
        ix_new=ix;
    }
    else{
        ix_new = (lx_-1 - ix);
    }
    i_new = Coordinates_.Ncell(ix_new, iy);

    return i_new;

}

void Connections_TL::Print_Hopping_YCn(){


    string fileout_HF="YCn_Lattice_Hoppings.txt";
    ofstream file_Hopping_out_HF(fileout_HF.c_str());

    file_Hopping_out_HF<<"#site_i spin_i site_j spin_j Hopping[site_i,spin_i][site_j,spin_j]"<<endl;



    Matrix<double> Hopp_mat;
    Hopp_mat.resize(ncells_*n_orbs_, ncells_*n_orbs_);
    string fileout= "REALMatrix_YCn_Snake3_form_SPINUP" + Parameters_.File_Hoppings;
    ofstream file_Hopping_out(fileout.c_str());

    string fileout3= "YCn_ZZ_REALMatrixUPPERTRIANGLE_form_SPINUP" + Parameters_.File_Hoppings;
    ofstream file_Hopping_out3(fileout3.c_str());


    int index_i, index_j;
    int ix, iy, ix_new, iy_new, i_new;
    int jx, jy, jx_new, jy_new, j_new;
    int ix_neigh_,iy_neigh_;


    //Hopping mat in conventional notation------

    for(int i=0;i<ncells_*n_orbs_;i++){
        ix = Coordinates_.indx_cellwise(i);
        iy = Coordinates_.indy_cellwise(i);


        Mat_1_int ix_neigh, iy_neigh;
        int i_neigh;

        //+x
        if(ix==(lx_-1)){
            if(Parameters_.PBC_X){
                ix_neigh_ = (ix + 1)%lx_;
                iy_neigh_ = iy;
                ix_neigh.push_back(ix_neigh_);
                iy_neigh.push_back(iy_neigh_);
            }
        }
        else{
            ix_neigh_ = (ix + 1);
            iy_neigh_ = iy;
            ix_neigh.push_back(ix_neigh_);
            iy_neigh.push_back(iy_neigh_);
        }


        //+y
        if(iy==(ly_-1)){
            if(Parameters_.PBC_Y){
                ix_neigh_ = ix;
                iy_neigh_ = (iy+1)%ly_;
                ix_neigh.push_back(ix_neigh_);
                iy_neigh.push_back(iy_neigh_);
            }
        }
        else{
            ix_neigh_ = ix;
            iy_neigh_ = (iy+1)%ly_;
            ix_neigh.push_back(ix_neigh_);
            iy_neigh.push_back(iy_neigh_);
        }

        //diagonal
        bool allowed_=true;
        if(iy%2==0){
            ix_neigh_ = (ix +1)%lx_;
            iy_neigh_ = (iy+1)%ly_;

            if( ix==(lx_-1)){
                if(Parameters_.PBC_X){
            allowed_= allowed_ && true;}
                else{
            allowed_=allowed_ && false;
                }
            }

            if( iy==(ly_-1)){
                if(Parameters_.PBC_Y){
            allowed_=allowed_ && true;}
                else{
            allowed_=allowed_ && false;
                }
            }


        }
        else{ //iy%2 !=0
            ix_neigh_ = (ix-1+lx_)%lx_;
            iy_neigh_ = (iy+1)%ly_;

            if( ix==0){
                if(Parameters_.PBC_X){
            allowed_= allowed_ && true;}
                else{
            allowed_=allowed_ && false;
                }
            }

            if( iy==(ly_-1)){
                if(Parameters_.PBC_Y){
            allowed_=allowed_ && true;}
                else{
            allowed_=allowed_ && false;
                }
            }
        }
        if(allowed_){
        ix_neigh.push_back(ix_neigh_);
        iy_neigh.push_back(iy_neigh_);
        }

       
        // cout<<"--------- i = "<<i<<"("<<ix <<","<<iy<<")"<<"-----------"<<endl;
        // for(int neigh_no=0;neigh_no<ix_neigh.size();neigh_no++){
        //     cout<< ix_neigh[neigh_no]<<"   "<<iy_neigh[neigh_no]<<"   "<<Coordinates_.Ncell(ix_neigh[neigh_no], iy_neigh[neigh_no])<<endl;
        // }
        // cout<<"================================"<<endl;

        for(int neigh_no=0;neigh_no<ix_neigh.size();neigh_no++){
            i_neigh = Coordinates_.Ncell(ix_neigh[neigh_no], iy_neigh[neigh_no]);
            Hopp_mat(i_neigh,i)=-1.0;
            Hopp_mat(i,i_neigh)=-1.0;
        }

    }
    //-------------------------------------------


    double val_temp;
     for(int i=0;i<ncells_*n_orbs_;i++){
        for(int j=0;j<ncells_*n_orbs_;j++){
         if(j>=i){
            val_temp=Hopp_mat(i,j);
         }
         else{
            val_temp=0;
         }
        file_Hopping_out3<<val_temp<<"   ";
        }
        file_Hopping_out3<<endl;
    }



    //file_Hopping_out<<"#site_i spin_i site_j spin_j Hopping[site_i,spin_i][site_j,spin_j]"<<endl;


    //Hopping for Snake-3------used later in DMRG----
    for(int i=0;i<ncells_*n_orbs_;i++){
        ix = Coordinates_.indx_cellwise(i);
        iy = Coordinates_.indy_cellwise(i);
        if(iy%2==0){
            ix_new=ix;
        }
        else{
            ix_new = (lx_-1 - ix);
        }
        i_new = Coordinates_.Ncell(ix_new, iy);
        //index_i=i_new+spin_i*ncells_*n_orbs_;

        for(int j=0;j<ncells_*n_orbs_;j++){
            jx = Coordinates_.indx_cellwise(j);
            jy = Coordinates_.indy_cellwise(j);

            if(jy%2==0){
                jx_new=jx;
            }
            else{
                jx_new = (lx_-1 - jx);
            }

            j_new = Coordinates_.Ncell(jx_new, jy);
            //index_j=j_new+spin_i*ncells_*n_orbs_;

            file_Hopping_out<<Hopp_mat(i_new,j_new)<<"   ";
        }
        file_Hopping_out<<endl;
    }
    //----------------------------










    //----------------------------------------------
    //-------------------------------------------------------
    int i1, i2, i3;
    int i1_new, i2_new, i3_new;
    string fileout2= "Ring_Exchange_YCn_Snake3_" + Parameters_.File_Hoppings;
    ofstream file_Hopping_out2(fileout2.c_str());
    bool Hit_boundary_x, Hit_boundary_y;

    //file_Hopping_out<<"#site_i spin_i site_j spin_j Hopping[site_i,spin_i][site_j,spin_j]"<<endl;
    for(int i=0;i<ncells_;i++){
        ix=Coordinates_.indx_cellwise(i);
        iy=Coordinates_.indy_cellwise(i);

        //ring-1
        if(iy%2==0){
            i1=Coordinates_.getneigh(i,4);}
        else{
            i1=Coordinates_.getneigh(i,2);
        }
        Hit_boundary_x=Coordinates_.HIT_X_BC;
        Hit_boundary_y=Coordinates_.HIT_Y_BC;

        if(iy%2==0){
            i2=Coordinates_.getneigh(i,11);}
        else{
            i2=Coordinates_.getneigh(i,11);
        }
        Hit_boundary_x= (Hit_boundary_x || Coordinates_.HIT_X_BC);
        Hit_boundary_y= (Hit_boundary_y || Coordinates_.HIT_Y_BC);

        if(iy%2==0){
            i3=Coordinates_.getneigh(i,2);}
        else{
            i3=Coordinates_.getneigh(i,5);
        }
        Hit_boundary_x= (Hit_boundary_x || Coordinates_.HIT_X_BC);
        Hit_boundary_y= (Hit_boundary_y || Coordinates_.HIT_Y_BC);

        if( (!(Hit_boundary_x) || Parameters_.PBC_X)
                &&
                (!(Hit_boundary_y) || Parameters_.PBC_Y)
                ){

            i_new = Convert_to_snake3_coordinates_YCn(i);
            i1_new = Convert_to_snake3_coordinates_YCn(i1);
            i2_new = Convert_to_snake3_coordinates_YCn(i2);
            i3_new = Convert_to_snake3_coordinates_YCn(i3);

            file_Hopping_out2<<i_new<<" "<<i1_new<<" "<<i2_new<<" "<<i3_new<<"   VAL1"<<endl;
            file_Hopping_out2<<i_new<<" "<<i3_new<<" "<<i1_new<<" "<<i2_new<<"   VAL1"<<endl;
            file_Hopping_out2<<i_new<<" "<<i2_new<<" "<<i1_new<<" "<<i3_new<<"   VAL2"<<endl;


        }


        //ring-2
        if(iy%2==0){
            i1=Coordinates_.getneigh(i,0);}
        else{
            i1=Coordinates_.getneigh(i,0);
        }
        Hit_boundary_x=Coordinates_.HIT_X_BC;
        Hit_boundary_y=Coordinates_.HIT_Y_BC;

        if(iy%2==0){
            i2=Coordinates_.getneigh(i,17);}
        else{
            i2=Coordinates_.getneigh(i,4);
        }
        Hit_boundary_x= (Hit_boundary_x || Coordinates_.HIT_X_BC);
        Hit_boundary_y= (Hit_boundary_y || Coordinates_.HIT_Y_BC);

        if(iy%2==0){
            i3=Coordinates_.getneigh(i,4);}
        else{
            i3=Coordinates_.getneigh(i,2);
        }
        Hit_boundary_x= (Hit_boundary_x || Coordinates_.HIT_X_BC);
        Hit_boundary_y= (Hit_boundary_y || Coordinates_.HIT_Y_BC);

        if( (!(Hit_boundary_x) || Parameters_.PBC_X)
                &&
                (!(Hit_boundary_y) || Parameters_.PBC_Y)
                ){
            i_new = Convert_to_snake3_coordinates_YCn(i);
            i1_new = Convert_to_snake3_coordinates_YCn(i1);
            i2_new = Convert_to_snake3_coordinates_YCn(i2);
            i3_new = Convert_to_snake3_coordinates_YCn(i3);

            file_Hopping_out2<<i_new<<" "<<i1_new<<" "<<i2_new<<" "<<i3_new<<"   VAL1"<<endl;
            file_Hopping_out2<<i_new<<" "<<i3_new<<" "<<i1_new<<" "<<i2_new<<"   VAL1"<<endl;
            file_Hopping_out2<<i_new<<" "<<i2_new<<" "<<i1_new<<" "<<i3_new<<"   VAL2"<<endl;
        }



        //ring-3
        if(iy%2==0){
            i1=Coordinates_.getneigh(i,7);}
        else{
            i1=Coordinates_.getneigh(i,3);
        }
        Hit_boundary_x=Coordinates_.HIT_X_BC;
        Hit_boundary_y=Coordinates_.HIT_Y_BC;

        if(iy%2==0){
            i2=Coordinates_.getneigh(i,8);}
        else{
            i2=Coordinates_.getneigh(i,7);
        }
        Hit_boundary_x= (Hit_boundary_x || Coordinates_.HIT_X_BC);
        Hit_boundary_y= (Hit_boundary_y || Coordinates_.HIT_Y_BC);

        if(iy%2==0){
            i3=Coordinates_.getneigh(i,0);}
        else{
            i3=Coordinates_.getneigh(i,0);
        }
        Hit_boundary_x= (Hit_boundary_x || Coordinates_.HIT_X_BC);
        Hit_boundary_y= (Hit_boundary_y || Coordinates_.HIT_Y_BC);

        if( (!(Hit_boundary_x) || Parameters_.PBC_X)
                &&
                (!(Hit_boundary_y) || Parameters_.PBC_Y)
                ){
            i_new = Convert_to_snake3_coordinates_YCn(i);
            i1_new = Convert_to_snake3_coordinates_YCn(i1);
            i2_new = Convert_to_snake3_coordinates_YCn(i2);
            i3_new = Convert_to_snake3_coordinates_YCn(i3);

            file_Hopping_out2<<i_new<<" "<<i1_new<<" "<<i2_new<<" "<<i3_new<<"   VAL1"<<endl;
            file_Hopping_out2<<i_new<<" "<<i3_new<<" "<<i1_new<<" "<<i2_new<<"   VAL1"<<endl;
            file_Hopping_out2<<i_new<<" "<<i2_new<<" "<<i1_new<<" "<<i3_new<<"   VAL2"<<endl;

        }



    }


    //-----------------------------------------


    //Hopping file for HF Run-----
    file_Hopping_out_HF<<"#site_i spin_i site_j spin_j Hopping[site_i,spin_i][site_j,spin_j]"<<endl;
    double hopp_val;
    for(int i=0;i<ncells_*n_orbs_;i++){
        ix = Coordinates_.indx_cellwise(i);
        iy = Coordinates_.indy_cellwise(i);
        if(iy%2==0){
            ix_new=ix;
        }
        else{
            ix_new = (lx_-1 - ix);
        }
        i_new = Coordinates_.Ncell(ix_new, iy);


        for(int spin_i=0;spin_i<2;spin_i++){
            index_i=i+spin_i*ncells_*n_orbs_;

            for(int j=0;j<ncells_*n_orbs_;j++){
                jx = Coordinates_.indx_cellwise(j);
                jy = Coordinates_.indy_cellwise(j);

                if(jy%2==0){
                    jx_new=jx;
                }
                else{
                    jx_new = (lx_-1 - jx);
                }

                j_new = Coordinates_.Ncell(jx_new, jy);

                for(int spin_j=0;spin_j<2;spin_j++){
                    index_j=j+spin_j*ncells_*n_orbs_;


                    if(i>j){
                        hopp_val = Hopp_mat(i_new,j_new);
                    }
                    else{
                        hopp_val = Hopp_mat(j_new,i_new);
                    }

                    if(spin_j!=spin_i){
                        hopp_val=0.0;
                    }

                    file_Hopping_out_HF<<i<<"  "<<spin_i<<"  "<<j<<"  "<<spin_j<<"  "<<hopp_val<<endl;

                }
            }
        }
    }
    //---------------------------------



}



void Connections_TL::Print_RingExchange_YCn(){


    Matrix<double> Hopp_mat;
    Hopp_mat.resize(ncells_*n_orbs_, ncells_*n_orbs_);
    string fileout= "HoppMatrix_YCn_form_SPINUP" + Parameters_.File_Hoppings;
    ofstream file_Hopping_out(fileout.c_str());

    int index_i, index_j;
    int ix, iy, ix_new, iy_new, i_new;
    int jx, jy, jx_new, jy_new, j_new;
    int ix_neigh_,iy_neigh_;


    //Hopping mat in conventional notation------

    for(int i=0;i<ncells_*n_orbs_;i++){
        ix = Coordinates_.indx_cellwise(i);
        iy = Coordinates_.indy_cellwise(i);


        Mat_1_int ix_neigh, iy_neigh;
        int i_neigh;

        //+x
        if(ix==(lx_-1)){
            if(Parameters_.PBC_X){
                ix_neigh_ = (ix + 1)%lx_;
                iy_neigh_ = iy;
                ix_neigh.push_back(ix_neigh_);
                iy_neigh.push_back(iy_neigh_);
            }
        }
        else{
            ix_neigh_ = (ix + 1);
            iy_neigh_ = iy;
            ix_neigh.push_back(ix_neigh_);
            iy_neigh.push_back(iy_neigh_);
        }


        //+y
        if(iy==(ly_-1)){
            if(Parameters_.PBC_Y){
                ix_neigh_ = ix;
                iy_neigh_ = (iy+1)%ly_;
                ix_neigh.push_back(ix_neigh_);
                iy_neigh.push_back(iy_neigh_);
            }
        }
        else{
            ix_neigh_ = ix;
            iy_neigh_ = (iy+1)%ly_;
            ix_neigh.push_back(ix_neigh_);
            iy_neigh.push_back(iy_neigh_);
        }

        //diagonal
        if(ix%2==0){
            ix_neigh_ = (ix +1)%lx_;
            iy_neigh_ = (iy+1)%ly_;
        }
        else{
            ix_neigh_ = (ix+1)%lx_;
            iy_neigh_ = (iy-1+ly_)%ly_;
        }
        if( (iy==(ly_-1)) && (ix==(lx_-1)) ){
            if(Parameters_.PBC_Y && Parameters_.PBC_X){
                ix_neigh.push_back(ix_neigh_);
                iy_neigh.push_back(iy_neigh_);
            }

        }
        if((iy==(ly_-1)) && (ix!=(lx_-1)) ){
            if(Parameters_.PBC_Y){
                ix_neigh.push_back(ix_neigh_);
                iy_neigh.push_back(iy_neigh_);
            }

        }
        if((iy!=(ly_-1)) && (ix==(lx_-1))  ){
            if(Parameters_.PBC_X){
                ix_neigh.push_back(ix_neigh_);
                iy_neigh.push_back(iy_neigh_);
            }

        }
        if((iy!=(ly_-1)) && (ix!=(lx_-1)) ){
            ix_neigh.push_back(ix_neigh_);
            iy_neigh.push_back(iy_neigh_);
        }



        cout<<"--------- i = "<<i<<"("<<ix <<","<<iy<<")"<<"-----------"<<endl;
        for(int neigh_no=0;neigh_no<ix_neigh.size();neigh_no++){
            cout<< ix_neigh[neigh_no]<<"   "<<iy_neigh[neigh_no]<<"   "<<Coordinates_.Ncell(ix_neigh[neigh_no], iy_neigh[neigh_no])<<endl;
        }
        cout<<"================================"<<endl;

        for(int neigh_no=0;neigh_no<ix_neigh.size();neigh_no++){
            i_neigh = Coordinates_.Ncell(ix_neigh[neigh_no], iy_neigh[neigh_no]);

            Hopp_mat(i_neigh,i)=-1.0;
            Hopp_mat(i,i_neigh)=-1.0;
        }




    }
    //-------------------------------------------



    //file_Hopping_out<<"#site_i spin_i site_j spin_j Hopping[site_i,spin_i][site_j,spin_j]"<<endl;


    //Hopping for YC-n------used later in Lanczos----
    for(int i=0;i<ncells_*n_orbs_;i++){
        ix = Coordinates_.indx_cellwise(i);
        iy = Coordinates_.indy_cellwise(i);
        if(iy%2==0){
            ix_new=ix;
        }
        else{
            ix_new = (ix);
        }
        i_new = Coordinates_.Ncell(ix_new, iy);
        //index_i=i_new+spin_i*ncells_*n_orbs_;

        for(int j=0;j<ncells_*n_orbs_;j++){
            jx = Coordinates_.indx_cellwise(j);
            jy = Coordinates_.indy_cellwise(j);

            if(jy%2==0){
                jx_new=jx;
            }
            else{
                jx_new = (jx);
            }

            j_new = Coordinates_.Ncell(jx_new, jy);
            //index_j=j_new+spin_i*ncells_*n_orbs_;

            file_Hopping_out<<Hopp_mat(i_new,j_new)<<"   ";
        }
        file_Hopping_out<<endl;
    }
    //----------------------------


    //-------------------------------------------------------
    int i1, i2, i3;
    string fileout2= "Ring_Exchange_YCn_" + Parameters_.File_Hoppings;
    ofstream file_Hopping_out2(fileout2.c_str());
    bool Hit_boundary_x, Hit_boundary_y;

    //file_Hopping_out<<"#site_i spin_i site_j spin_j Hopping[site_i,spin_i][site_j,spin_j]"<<endl;
    for(int i=0;i<ncells_;i++){
        ix=Coordinates_.indx_cellwise(i);
        iy=Coordinates_.indy_cellwise(i);


        //ring-1
        if(ix%2==0){
            i1=Coordinates_.getneigh(i,4);}
        else{
            i1=Coordinates_.getneigh(i,0);
        }
        Hit_boundary_x=Coordinates_.HIT_X_BC;
        Hit_boundary_y=Coordinates_.HIT_Y_BC;

        if(ix%2==0){
            i2=Coordinates_.getneigh(i,0);}
        else{
            i2=Coordinates_.getneigh(i,7);
        }
        Hit_boundary_x= (Hit_boundary_x || Coordinates_.HIT_X_BC);
        Hit_boundary_y= (Hit_boundary_y || Coordinates_.HIT_Y_BC);

        if(ix%2==0){
            i3=Coordinates_.getneigh(i,3);}
        else{
            i3=Coordinates_.getneigh(i,3);
        }
        Hit_boundary_x= (Hit_boundary_x || Coordinates_.HIT_X_BC);
        Hit_boundary_y= (Hit_boundary_y || Coordinates_.HIT_Y_BC);

        if( (!(Hit_boundary_x) || Parameters_.PBC_X)
                &&
                (!(Hit_boundary_y) || Parameters_.PBC_Y)
                ){
            file_Hopping_out2<<i<<" "<<i1<<" "<<i2<<" "<<i3<<"   VAL1"<<endl;
            file_Hopping_out2<<i<<" "<<i3<<" "<<i1<<" "<<i2<<"   VAL1"<<endl;
            file_Hopping_out2<<i<<" "<<i2<<" "<<i1<<" "<<i3<<"   VAL2"<<endl;


        }


        //ring-2
        if(ix%2==0){
            i1=Coordinates_.getneigh(i,2);}
        else{
            i1=Coordinates_.getneigh(i,2);
        }
        Hit_boundary_x=Coordinates_.HIT_X_BC;
        Hit_boundary_y=Coordinates_.HIT_Y_BC;

        if(ix%2==0){
            i2=Coordinates_.getneigh(i,4);}
        else{
            i2=Coordinates_.getneigh(i,0);
        }
        Hit_boundary_x= (Hit_boundary_x || Coordinates_.HIT_X_BC);
        Hit_boundary_y= (Hit_boundary_y || Coordinates_.HIT_Y_BC);

        if(ix%2==0){
            i3=Coordinates_.getneigh(i,0);}
        else{
            i3=Coordinates_.getneigh(i,7);
        }
        Hit_boundary_x= (Hit_boundary_x || Coordinates_.HIT_X_BC);
        Hit_boundary_y= (Hit_boundary_y || Coordinates_.HIT_Y_BC);

        if( (!(Hit_boundary_x) || Parameters_.PBC_X)
                &&
                (!(Hit_boundary_y) || Parameters_.PBC_Y)
                ){
            file_Hopping_out2<<i<<" "<<i1<<" "<<i2<<" "<<i3<<"   VAL1"<<endl;
            file_Hopping_out2<<i<<" "<<i3<<" "<<i1<<" "<<i2<<"   VAL1"<<endl;
            file_Hopping_out2<<i<<" "<<i2<<" "<<i1<<" "<<i3<<"   VAL2"<<endl;

        }



        //ring-3
        if(ix%2==0){
            i1=Coordinates_.getneigh(i,5);}
        else{
            i1=Coordinates_.getneigh(i,1);
        }
        Hit_boundary_x=Coordinates_.HIT_X_BC;
        Hit_boundary_y=Coordinates_.HIT_Y_BC;

        if(ix%2==0){
            i2=Coordinates_.getneigh(i,2);}
        else{
            i2=Coordinates_.getneigh(i,2);
        }
        Hit_boundary_x= (Hit_boundary_x || Coordinates_.HIT_X_BC);
        Hit_boundary_y= (Hit_boundary_y || Coordinates_.HIT_Y_BC);

        if(ix%2==0){
            i3=Coordinates_.getneigh(i,4);}
        else{
            i3=Coordinates_.getneigh(i,0);
        }
        Hit_boundary_x= (Hit_boundary_x || Coordinates_.HIT_X_BC);
        Hit_boundary_y= (Hit_boundary_y || Coordinates_.HIT_Y_BC);

        if( (!(Hit_boundary_x) || Parameters_.PBC_X)
                &&
                (!(Hit_boundary_y) || Parameters_.PBC_Y)
                ){
            file_Hopping_out2<<i<<" "<<i1<<" "<<i2<<" "<<i3<<"   VAL1"<<endl;
            file_Hopping_out2<<i<<" "<<i3<<" "<<i1<<" "<<i2<<"   VAL1"<<endl;
            file_Hopping_out2<<i<<" "<<i2<<" "<<i1<<" "<<i3<<"   VAL2"<<endl;

        }



    }



}



void Connections_TL::Print_Hopping_Snake3(){

    string fileout="Snake3_" + Parameters_.File_Hoppings;
    ofstream file_Hopping_out(fileout.c_str());

    file_Hopping_out<<"#site_i spin_i site_j spin_j Hopping[site_i,spin_i][site_j,spin_j]"<<endl;

    int index_i, index_j;
    int ix, iy, ix_new, iy_new, i_new;
    int jx, jy, jx_new, jy_new, j_new;

    for(int i=0;i<ncells_*n_orbs_;i++){
        ix = Coordinates_.indx_cellwise(i);
        iy = Coordinates_.indy_cellwise(i);
        if(iy%2==0){
            ix_new=ix;
        }
        else{
            ix_new = (lx_-1 - ix);
        }
        i_new = Coordinates_.Ncell(ix_new, iy);

        for(int spin_i=0;spin_i<2;spin_i++){
            index_i=i_new+spin_i*ncells_*n_orbs_;


            for(int j=0;j<ncells_*n_orbs_;j++){
                jx = Coordinates_.indx_cellwise(j);
                jy = Coordinates_.indy_cellwise(j);

                if(jy%2==0){
                    jx_new=jx;
                }
                else{
                    jx_new = (lx_-1 - jx);
                }

                j_new = Coordinates_.Ncell(jx_new, jy);
                for(int spin_j=0;spin_j<2;spin_j++){

                    index_j=j_new+spin_j*ncells_*n_orbs_;

                    if(abs(HTB_(index_i,index_j))>0.0000001){
                        file_Hopping_out<<i<<"  "<<spin_i<<"  "<<j<<"  "<<spin_j<<"  "<<HTB_(index_i,index_j)<<endl;
                    }
                }
            }
        }
    }

}

void Connections_TL::Print_LongRangeInt(){
    string fileout=Parameters_.File_LongRange_Ints;
    ofstream file_out(fileout.c_str());

    file_out<<"#i j Uij"<<endl;

    for(int i=0;i<ncells_;i++){

        for(int j=0;j<ncells_;j++){

            if(abs(Hint_(i,j))>0.0000001){
                file_out<<i<<"   "<<j<<"   "<<Hint_(i,j)<<endl;
            }

        }
    }

}

void Connections_TL::Print_LongRangeInt_SNAKE3(){
    string fileout="REALMatrix_Snake3_LONG_RANGE_INT" + Parameters_.File_LongRange_Ints;
    ofstream file_out(fileout.c_str());

    int index_i, index_j;
    int ix, iy, ix_new, iy_new, i_new;
    int jx, jy, jx_new, jy_new, j_new;

    file_out<<"#i j Uij"<<endl;

    for(int i=0;i<ncells_;i++){
        ix = Coordinates_.indx_cellwise(i);
        iy = Coordinates_.indy_cellwise(i);
        if(iy%2==0){
            ix_new=ix;
        }
        else{
            ix_new = (lx_-1 - ix);
        }
        i_new = Coordinates_.Ncell(ix_new, iy);

        for(int j=0;j<ncells_;j++){

            jx = Coordinates_.indx_cellwise(j);
            jy = Coordinates_.indy_cellwise(j);

            if(jy%2==0){
                jx_new=jx;
            }
            else{
                jx_new = (lx_-1 - jx);
            }

            j_new = Coordinates_.Ncell(jx_new, jy);

            file_out<<Hint_(i_new,j_new).real()<<"   ";


        }
        file_out<<endl;
    }

}

void Connections_TL::Print_Spin_resolved_OnsiteE(){

    string fileout=Parameters_.File_Onsite_Energies;
    ofstream file_out(fileout.c_str());
    file_out<<"#Site Spin Onsite_E"<<endl;
    int index, i_posx, i_posy;

    for (int i = 0; i < ncells_; i++)
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
