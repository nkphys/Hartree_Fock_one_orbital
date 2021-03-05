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
    void Print_LongRangeInt();
    void Print_Spin_resolved_OnsiteE();
    void InteractionsCreate();                             //::DONE
    void Check_Hermiticity();                              //::DONE
    void Check_up_down_symmetry();                         //::DONE
    void HTBCreate();                                      //::DONE



    Parameters_TL &Parameters_;
    Coordinates_TL &Coordinates_;
    int lx_, ly_, ncells_, n_orbs_;
    Matrix<complex<double>> HTB_;
    Matrix<complex<double>> Hint_;
    Matrix<complex<double>> Ham_;
    Matrix<double> Tx, Ty, Tpxpy, Tpxmy;
    Matrix<double> IntraCell_Hopp, InterCell_px_Hopp, InterCell_py_Hopp, InterCell_pxmy_Hopp;

};


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

} // ----------


void Connections_TL::InteractionsCreate()
{

    Mat_1_int U_neighs;
    Mat_1_doub U_hoppings;
    U_neighs.push_back(0);U_neighs.push_back(2);U_neighs.push_back(7);
    U_hoppings.push_back(Parameters_.U1);U_hoppings.push_back(Parameters_.U1);U_hoppings.push_back(Parameters_.U1);

    U_neighs.push_back(4);U_neighs.push_back(8);U_neighs.push_back(9);
    U_hoppings.push_back(Parameters_.U2);U_hoppings.push_back(Parameters_.U2);U_hoppings.push_back(Parameters_.U2);

    U_neighs.push_back(10);U_neighs.push_back(11);U_neighs.push_back(12);
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
            m = Coordinates_.neigh(l, U_neighs[neigh]); //+x neighbour cell
            m1_pos = Coordinates_.indx_cellwise(m);
            m2_pos = Coordinates_.indy_cellwise(m);

            assert(l != m);
            if (l != m)
            {
                Hint_(l, m) += complex<double>(1.0 * U_hoppings[neigh], 0.0);
                Hint_(m, l) += conj(Hint_(l, m));
            }

        }

    }
    }
    else{
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

//                if(m1_pos>l1_pos){
//                r1_ = min(abs(m1_pos-l1_pos),abs(l1_pos+lx_-m1_pos)) ;
//                }
                r1_= min(abs(m1_pos-l1_pos),lx_-abs(m1_pos-l1_pos));
                r2_= min(abs(m2_pos-l2_pos),ly_-abs(m2_pos-l2_pos));

                rx_ = ((sqrt(3.0)/2.0)*(r1_) +  (sqrt(3.0)/2.0)*(r2_));
                ry_ =  (-0.5*(r1_) + 0.5*(r2_));
                dis_= sqrt(rx_*rx_ + ry_*ry_);
//(14.3952)*( (1.0/(x*62.6434))  - (1.0/(sqrt( (x*x*62.6434*62.6434)  + d*d)))  )

                U_val = ((14.3952*1000)/Parameters_.eps_DE)*( (1.0/(dis_*Parameters_.a_moire))  -
                                    (1.0/(sqrt( (dis_*dis_*Parameters_.a_moire*Parameters_.a_moire)
                                    + Parameters_.d_screening*Parameters_.d_screening)))  );
               // assert(l != m);
                if (l != m)
                {
                    Hint_(l, m) += complex<double>(1.0 * U_val, 0.0);
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
    Mat_1_doub t_hoppings;
    //t1 hoppings
    t_neighs.push_back(0);t_neighs.push_back(2);t_neighs.push_back(7);
    t_hoppings.push_back(Parameters_.t1);t_hoppings.push_back(Parameters_.t1);t_hoppings.push_back(Parameters_.t1);

    t_neighs.push_back(4);t_neighs.push_back(8);t_neighs.push_back(9);
    t_hoppings.push_back(Parameters_.t2);t_hoppings.push_back(Parameters_.t2);t_hoppings.push_back(Parameters_.t2);

    t_neighs.push_back(10);t_neighs.push_back(11);t_neighs.push_back(12);
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
        for(int neigh=0;neigh<9;neigh++){
            m = Coordinates_.neigh(l, t_neighs[neigh]); //+x neighbour cell
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
                            HTB_(b, a) = complex<double>(1.0 * t_hoppings[neigh], 0.0);
                            HTB_(a, b) = conj(HTB_(b, a));
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
