#include <algorithm>
#include <functional>
#include <math.h>
#include "Parameters_G2dLatticeNew.h"
#include "Coordinates_G2dLatticeNew.h"
#define PI acos(-1.0)

#ifndef Connections_G2dLatticeNew_class
#define Connections_G2dLatticeNew_class


class Connections_G2dLatticeNew
{
public:
Connections_G2dLatticeNew(Parameters_G2dLatticeNew &Parameters__, Coordinates_G2dLatticeNew &Coordinates__)
: Parameters_(Parameters__), Coordinates_(Coordinates__)

{
Initialize();
HTBCreateNew();
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
void HTBCreateNew();                                      //::DONE

void Print_Ansatz_LocalDen_CDW();
void Get_minimum_distance_direction(int l, int gamma_l, int m, int gamma_m, int &r1_, int &r2_, double &dis_min);
void Get_MUC_details();


Parameters_G2dLatticeNew &Parameters_;
Coordinates_G2dLatticeNew &Coordinates_;
int lx_, ly_, nsites_, ncells_, n_orbs_, n_atoms_;
Matrix<complex<double>> HTB_;
Matrix<complex<double>> Hint_;
Matrix<complex<double>> Ham_;
Matrix<double> Tx, Ty, Tpxpy, Tpxmy;
Matrix<double> IntraCell_Hopp, InterCell_px_Hopp, InterCell_py_Hopp, InterCell_pxmy_Hopp;

int max_l1_dis, max_l2_dis;
Matrix<int> Mat_MUC;
int NSites_in_MUC;
int lx_cells, ly_cells;
int ncells;

Mat_1_intpair Intra_MUC_positions;
Mat_2_int Intra_MUC_index;

};

void Connections_G2dLatticeNew::Get_MUC_details(){


    NSites_in_MUC = (Mat_MUC(0,0)*Mat_MUC(1,1) - Mat_MUC(0,1)*Mat_MUC(1,0));

    if(NSites_in_MUC<0){
        cout<<"Change the MUC matrix so that Determinant is >0"<<endl;
        assert(false);
    }

    if((NSites_in_MUC) == 0){
        cout<<"MUC given is empty"<<endl;
        assert((NSites_in_MUC) != 0);
    }

    if((lx_*ly_)%(NSites_in_MUC) !=0){
        cout<<"MUC given is not commensurate with the Lattice"<<endl;
        assert((lx_*ly_)%(NSites_in_MUC) ==0);
    }


    ncells_ = (lx_*ly_)/(NSites_in_MUC);


    int max_i1p=0;
    int max_i2p=0;
    int i1p=0;
    int i2p=0;
    Mat_1_int Site_included;
    Site_included.resize(lx_*ly_);
    for(int i=0;i<lx_*ly_;i++){
        Site_included[i]=0;
    }

    int i1_temp, i2_temp;
    while(i1p<lx_*ly_){

        i2p=0;
        while(i2p<lx_*ly_){

            i1_temp = ((i1p*Mat_MUC(0,0) + i2p*Mat_MUC(1,0)) + lx_*ly_)%lx_;
            i2_temp = ((i1p*Mat_MUC(0,1) + i2p*Mat_MUC(1,1)) + lx_*ly_)%ly_;

                if(Site_included[i1_temp + lx_*i2_temp]==0){
                    max_i1p = max(max_i1p,i1p);
                    max_i2p = max(max_i2p,i2p);
                    Site_included[i1_temp + lx_*i2_temp]=1;
                }
            i2p += 1;
        }
        i1p += 1;
    }

    lx_cells = max_i1p+1;
    ly_cells = max_i2p+1;



    //Intra_MUC_positions

    Intra_MUC_index.resize(lx_);
    for(int ix=0;ix<lx_;ix++){
        Intra_MUC_index[ix].resize(ly_);
        for(int iy=0;iy<ly_;iy++){
            Intra_MUC_index[ix][iy]=-1000;
        }
    }
    bool inside_MUC;
    for(int h1=0;h1<lx_*ly_;h1++){
        for(int h2=0;h2<lx_*ly_;h2++){
            inside_MUC = ( (h1*Mat_MUC(1,1) - h2*Mat_MUC(1,0))<NSites_in_MUC);
            inside_MUC = inside_MUC &&  ( (h1*Mat_MUC(1,1) - h2*Mat_MUC(1,0))>=0);
            inside_MUC = inside_MUC && ((-h1*Mat_MUC(0,1) + h2*Mat_MUC(0,0))<NSites_in_MUC);
            inside_MUC = inside_MUC && ((-h1*Mat_MUC(0,1) + h2*Mat_MUC(0,0))>=0);

            if(inside_MUC){
                pair_int pos_;
                pos_.first  = h1;
                pos_.second = h2;
                Intra_MUC_positions.push_back(pos_);
                Intra_MUC_index[h1][h2]=Intra_MUC_positions.size()-1;
            }
        }
    }

    assert(Intra_MUC_positions.size() == NSites_in_MUC);



}



void Connections_G2dLatticeNew::Get_minimum_distance_direction(int l,int gamma_l, int m, int gamma_m, int &r1_, int &r2_, double &dis_min){


//in this routine, gamma_l(m) is atom index
cout<<"Connections_G2dLatticeNew::Get_minimum_distan... is not generic yet"<<endl;
assert(false);

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

void Connections_G2dLatticeNew::Initialize()
{

ly_ = Parameters_.ly;
lx_ = Parameters_.lx;
nsites_ = lx_ * ly_;  //total no of sites, each site has 2 orb
n_orbs_ = Parameters_.n_orbs;
n_atoms_ = Parameters_.n_atoms;
int space = 2 * n_orbs_*n_atoms_*nsites_ ;

//Hint_.resize(n_atoms_*n_orbs_*nsites_,n_atoms_*n_orbs_*nsites_);

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


Mat_MUC = Parameters_.MUC_Mat;
Get_MUC_details();
assert(ncells_ == lx_cells * ly_cells);


int space_col = 2 * n_orbs_*n_atoms_*NSites_in_MUC;

HTB_.resize(space, space_col);
//Ham_.resize(space, space);

} // ----------


void Connections_G2dLatticeNew::Print_Ansatz_LocalDen_CDW(){



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

void Connections_G2dLatticeNew::InteractionsCreate()
{


} // ----------



void Connections_G2dLatticeNew::Interactions_Sorting()
{

} // ----------


void Connections_G2dLatticeNew::Check_up_down_symmetry()

{

}

void Connections_G2dLatticeNew::Check_Hermiticity()

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




void Connections_G2dLatticeNew::HTBCreateNew()
{

    double EPS_doub=0.000001;

    Mat_1_int t_neighs;
    Mat_3_Complex_doub t_hoppings;
    t_neighs.clear();
    t_hoppings.clear();


    Mat_1_int t_neighs_conj;
    t_neighs_conj.clear();




    //PX
    t_neighs.push_back(0);t_hoppings.push_back(Parameters_.t1_plus_a1);

    //MY
    t_neighs.push_back(3);t_hoppings.push_back(Parameters_.t1_minus_a2);

    //PXMY
    t_neighs.push_back(7);t_hoppings.push_back(Parameters_.t1_plus_a1_minus_a2);

    //PXPY
    t_neighs.push_back(4);t_hoppings.push_back(Parameters_.t1_plus_a1_plus_a2);

    //For conj
    //MX
    t_neighs_conj.push_back(1);t_hoppings.push_back(Parameters_.t1_plus_a1);

    //PY
    t_neighs_conj.push_back(2);t_hoppings.push_back(Parameters_.t1_minus_a2);

    //MXPY
    t_neighs_conj.push_back(5);t_hoppings.push_back(Parameters_.t1_plus_a1_minus_a2);

    //MXMY
    t_neighs_conj.push_back(6);t_hoppings.push_back(Parameters_.t1_plus_a1_plus_a2);



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

    //HERE
    for (l = 0; l < Intra_MUC_positions.size(); l++)
    {

        lx_pos = Intra_MUC_positions[l].first; //cell x-position
        ly_pos = Intra_MUC_positions[l].second;


        //Intra-cell  hoppings
        for (int spin1 = 0; spin1 < 2; spin1++)
        {

            for (int spin2 = 0; spin2 < 2; spin2++)
            {

                for(int atom1=0;atom1<n_atoms_;atom1++)
                {
                    for(int atom2=0;atom2<n_atoms_;atom2++)
                    {

                        for(int orb1=0;orb1<n_orbs_;orb1++)
                        {
                            for(int orb2=0;orb2<n_orbs_;orb2++)
                            {

                                a = Coordinates_.Nbasis(lx_pos, ly_pos, atom1 + n_atoms_*orb1) + nsites_ * n_orbs_* n_atoms_ * spin1;
                                b = Coordinates_.Nbasis(lx_pos, ly_pos, atom2 + n_atoms_*orb2) + nsites_ * n_orbs_* n_atoms_ * spin2;

                                int a_col = (atom1 + n_atoms_*orb1) +l*(n_atoms_*n_orbs_) + NSites_in_MUC*n_orbs_*n_atoms_*spin1;
                                int b_col = (atom2 + n_atoms_*orb2) +l*(n_atoms_*n_orbs_) + NSites_in_MUC*n_orbs_*n_atoms_*spin2;


                                if (a != b && (abs(Parameters_.t0[atom2 + n_atoms_*orb2 + n_atoms_*n_orbs_*spin2][atom1 + n_atoms_*orb1 + n_atoms_*n_orbs_*spin1])>EPS_doub) )
                                {
                                    HTB_(b, a_col) += 1.0*Parameters_.t0[atom2 + n_atoms_*orb2 + n_atoms_*n_orbs_*spin2][atom1 + n_atoms_*orb1 + n_atoms_*n_orbs_*spin1];
                                    HTB_(a, b_col) += conj(1.0*Parameters_.t0[atom2 + n_atoms_*orb2 + n_atoms_*n_orbs_*spin2][atom1 + n_atoms_*orb1 + n_atoms_*n_orbs_*spin1]);

                                    //HTB_(a, b) = conj(HTB_(b, a));

                                }

                            }
                        }
                    }
                }
            }
        }





        //For t1,t2,t3 hoppings
        for (int spin1 = 0; spin1 < 2; spin1++)
        {
            for (int spin2 = 0; spin2 < 2; spin2++)
            {
                for(int neigh=0;neigh<t_neighs.size();neigh++){

                    int l_ind_sys = Coordinates_.Ncell(lx_pos,ly_pos);
                    m = Coordinates_.getneigh(l_ind_sys, t_neighs[neigh]);

                    if( (!Coordinates_.HIT_X_BC) || Parameters_.PBC_X){

                        if( (!Coordinates_.HIT_Y_BC) || Parameters_.PBC_Y){


                            mx_pos = Coordinates_.indx_cellwise(m);
                            my_pos = Coordinates_.indy_cellwise(m);

                            for (int atom1 = 0; atom1 < n_atoms_; atom1++)
                            {
                                for (int atom2 = 0; atom2 < n_atoms_; atom2++)
                                {

                                    for (int orb1 = 0; orb1 < n_orbs_; orb1++)
                                    {
                                        for (int orb2 = 0; orb2 < n_orbs_; orb2++)
                                        {
                                            a = Coordinates_.Nbasis(lx_pos, ly_pos, atom1 + n_atoms_*orb1) + nsites_ * n_orbs_* n_atoms_*spin1;
                                            b = Coordinates_.Nbasis(mx_pos, my_pos, atom2 + n_atoms_*orb2) + nsites_ * n_orbs_* n_atoms_*spin2;
                                            assert(a != b);

                                            int a_col = (atom1 + n_atoms_*orb1) +l*(n_atoms_*n_orbs_) + NSites_in_MUC*n_orbs_*n_atoms_*spin1;


                                            if (a != b && (abs(t_hoppings[neigh][atom2 + n_atoms_*orb2 + n_atoms_*n_orbs_*spin2][atom1 + n_atoms_*orb1 + n_atoms_*n_orbs_*spin1])>EPS_doub))
                                            {

                                                HTB_(b, a_col) += 1.0*t_hoppings[neigh][atom2+n_atoms_*orb2+n_atoms_*n_orbs_*spin2][atom1+n_atoms_*orb1+n_atoms_*n_orbs_*spin1];
                                                //HTB_(a, b) = conj(HTB_(b, a));


                                            }

                                        }
                                    }
                                }
                            }
                        }
                    }

                }

            }
        }




        //conj
        //For t1,t2,t3 hoppings
        for (int spin1 = 0; spin1 < 2; spin1++)
        {
            for (int spin2 = 0; spin2 < 2; spin2++)
            {
                for(int neigh=0;neigh<t_neighs_conj.size();neigh++){

                    int l_ind_sys = Coordinates_.Ncell(lx_pos,ly_pos);
                    m = Coordinates_.getneigh(l_ind_sys, t_neighs_conj[neigh]);

                    if( (!Coordinates_.HIT_X_BC) || Parameters_.PBC_X){

                        if( (!Coordinates_.HIT_Y_BC) || Parameters_.PBC_Y){


                            mx_pos = Coordinates_.indx_cellwise(m);
                            my_pos = Coordinates_.indy_cellwise(m);

                            for (int atom1 = 0; atom1 < n_atoms_; atom1++)
                            {
                                for (int atom2 = 0; atom2 < n_atoms_; atom2++)
                                {

                                    for (int orb1 = 0; orb1 < n_orbs_; orb1++)
                                    {
                                        for (int orb2 = 0; orb2 < n_orbs_; orb2++)
                                        {
                                            a = Coordinates_.Nbasis(lx_pos, ly_pos, atom2 + n_atoms_*orb2) + nsites_ * n_orbs_* n_atoms_*spin2;
                                            b = Coordinates_.Nbasis(mx_pos, my_pos, atom1 + n_atoms_*orb1) + nsites_ * n_orbs_* n_atoms_*spin1;
                                            assert(a != b);

                                            int a_col = (atom2 + n_atoms_*orb2) +l*(n_atoms_*n_orbs_) + NSites_in_MUC*n_orbs_*n_atoms_*spin2;


                                            if (a != b && (abs(t_hoppings[neigh][atom2 + n_atoms_*orb2 + n_atoms_*n_orbs_*spin2][atom1 + n_atoms_*orb1 + n_atoms_*n_orbs_*spin1])>EPS_doub))
                                            {

                                                HTB_(b, a_col) += 1.0*conj(t_hoppings[neigh][atom2+n_atoms_*orb2+n_atoms_*n_orbs_*spin2][atom1+n_atoms_*orb1+n_atoms_*n_orbs_*spin1]);
                                                //HTB_(a, b) = conj(HTB_(b, a));


                                            }

                                        }
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



void Connections_G2dLatticeNew::HTBCreate()
{

double EPS_doub=0.000001;

Mat_1_int t_neighs;
Mat_3_Complex_doub t_hoppings;
t_neighs.clear();
t_hoppings.clear();


//PX
t_neighs.push_back(0);t_hoppings.push_back(Parameters_.t1_plus_a1);

//MY
t_neighs.push_back(3);t_hoppings.push_back(Parameters_.t1_minus_a2);

//PXMY
t_neighs.push_back(7);t_hoppings.push_back(Parameters_.t1_plus_a1_minus_a2);

//PXPY
t_neighs.push_back(4);t_hoppings.push_back(Parameters_.t1_plus_a1_plus_a2);



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


//Intra-unitcell  hoppings
for (int spin1 = 0; spin1 < 2; spin1++)
{

    for (int spin2 = 0; spin2 < 2; spin2++)
{

for(int atom1=0;atom1<n_atoms_;atom1++)
{
        for(int atom2=0;atom2<n_atoms_;atom2++)
        {

for(int orb1=0;orb1<n_orbs_;orb1++)
{
        for(int orb2=0;orb2<n_orbs_;orb2++)
        {

    a = Coordinates_.Nbasis(lx_pos, ly_pos, atom1 + n_atoms_*orb1) + nsites_ * n_orbs_* n_atoms_ * spin1;
    b = Coordinates_.Nbasis(lx_pos, ly_pos, atom2 + n_atoms_*orb2) + nsites_ * n_orbs_* n_atoms_ * spin2;

    if (a != b && (abs(Parameters_.t0[atom2 + n_atoms_*orb2 + n_atoms_*n_orbs_*spin2][atom1 + n_atoms_*orb1 + n_atoms_*n_orbs_*spin1])>EPS_doub) )
    {
            HTB_(b, a) = 1.0*Parameters_.t0[atom2 + n_atoms_*orb2 + n_atoms_*n_orbs_*spin2][atom1 + n_atoms_*orb1 + n_atoms_*n_orbs_*spin1];
            HTB_(a, b) = conj(HTB_(b, a));
        
    }

    }
}
        }
        }
}
}





//For t1,t2,t3 hoppings
for (int spin1 = 0; spin1 < 2; spin1++)
{
    for (int spin2 = 0; spin2 < 2; spin2++)
{
for(int neigh=0;neigh<t_neighs.size();neigh++){

    m = Coordinates_.getneigh(l, t_neighs[neigh]);

    if( (!Coordinates_.HIT_X_BC) || Parameters_.PBC_X){

        if( (!Coordinates_.HIT_Y_BC) || Parameters_.PBC_Y){

 
            mx_pos = Coordinates_.indx_cellwise(m);
            my_pos = Coordinates_.indy_cellwise(m);

            for (int atom1 = 0; atom1 < n_atoms_; atom1++)
                {
                    for (int atom2 = 0; atom2 < n_atoms_; atom2++)
                    {

                for (int orb1 = 0; orb1 < n_orbs_; orb1++)
                {
                    for (int orb2 = 0; orb2 < n_orbs_; orb2++)
                    {
                        a = Coordinates_.Nbasis(lx_pos, ly_pos, atom1 + n_atoms_*orb1) + nsites_ * n_orbs_* n_atoms_*spin1;
                        b = Coordinates_.Nbasis(mx_pos, my_pos, atom2 + n_atoms_*orb2) + nsites_ * n_orbs_* n_atoms_*spin2;
                        assert(a != b);

                        if (a != b && (abs(t_hoppings[neigh][atom2 + n_atoms_*orb2 + n_atoms_*n_orbs_*spin2][atom1 + n_atoms_*orb1 + n_atoms_*n_orbs_*spin1])>EPS_doub))
                        {

                            HTB_(b, a) = 1.0*t_hoppings[neigh][atom2+n_atoms_*orb2+n_atoms_*n_orbs_*spin2][atom1+n_atoms_*orb1+n_atoms_*n_orbs_*spin1];
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
}

}




} // ----------


void Connections_G2dLatticeNew::Print_Hopping(){

string fileout=Parameters_.File_Hoppings;
ofstream file_Hopping_out(fileout.c_str());

file_Hopping_out<<"#site_i atom_i orb_i spin_i site_j atom_j orb_j spin_j Hopping[i][j]"<<endl;

int index_i, index_j;
int ix_pos, iy_pos, jx_pos, jy_pos; 
for(int site_i=0;site_i<nsites_;site_i++){
for(int atom_i=0;atom_i<n_atoms_;atom_i++){
for(int orb_i=0;orb_i<n_orbs_;orb_i++){
for(int spin_i=0;spin_i<1;spin_i++){
    ix_pos = Coordinates_.indx_cellwise(site_i);
    iy_pos = Coordinates_.indy_cellwise(site_i);
    index_i= Coordinates_.Nbasis(ix_pos, iy_pos, atom_i + n_atoms_*orb_i) + nsites_ * n_orbs_* n_atoms_*spin_i;
    

    for(int site_j=0;site_j<nsites_;site_j++){
    for(int atom_j=0;atom_j<n_atoms_;atom_j++){
    for(int orb_j=0;orb_j<n_orbs_;orb_j++){
        for(int spin_j=0;spin_j<1;spin_j++){
            jx_pos = Coordinates_.indx_cellwise(site_j);
            jy_pos = Coordinates_.indy_cellwise(site_j);
            index_j= Coordinates_.Nbasis(jx_pos, jy_pos, atom_j + n_atoms_*orb_j) + nsites_ * n_orbs_* n_atoms_*spin_j;

            if(abs(HTB_(index_i,index_j))>0.0000001){
                file_Hopping_out<<site_i<<"  "<<atom_i<<"  "<<orb_i<<"  "<<spin_i<<"  "<<site_j<<"  "<<atom_j<<"  "<<orb_j<<"  "<<spin_j<<"  "<<HTB_(index_i,index_j)<<endl;
            }
        }
    }
}
}
}
}
}
}
}



void Connections_G2dLatticeNew::Print_LongRangeInt(){
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



#endif
