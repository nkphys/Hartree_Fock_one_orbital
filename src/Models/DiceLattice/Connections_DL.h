#include <algorithm>
#include <functional>
#include <math.h>
#include "tensor_type.h"
#include "Parameters_DL.h"
#include "Coordinates_DL.h"
#include "../../Matrix.h"
#define PI acos(-1.0)

#ifndef Connections_DL_class
#define Connections_DL_class


class Connections_DL
{
public:
    Connections_DL(Parameters_DL &Parameters__, Coordinates_DL &Coordinates__)
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



    Parameters_DL &Parameters_;
    Coordinates_DL &Coordinates_;
    int lx_, ly_, ncells_, n_orbs_;
    Matrix<complex<double>> HTB_;
    Matrix<complex<double>> Ham_;
    Matrix<double> Tx, Ty, Tpxpy, Tpxmy;
    Matrix<double> IntraCell_Hopp, InterCell_px_Hopp, InterCell_py_Hopp, InterCell_pxmy_Hopp;

};


void Connections_DL::Initialize()
{


    ly_ = Parameters_.ly;
    lx_ = Parameters_.lx;
    ncells_ = lx_ * ly_;
    n_orbs_ = Parameters_.n_orbs;
    int space = 2 * ncells_ * n_orbs_;

    HTB_.resize(space, space);
    Ham_.resize(space, space);

} // ----------


void Connections_DL::InteractionsCreate()
{

    Ham_ = HTB_;
    // Ham_.print();

} // ----------


void Connections_DL::Check_up_down_symmetry()

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

void Connections_DL::Check_Hermiticity()

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


void Connections_DL::HTBCreate()
{

    //Convention used
    //orb=0=A
    //orb=1=B
    //orb=2=C


    Matrix<complex <double>> sigma_x, sigma_y, sigma_z, Value_mat;
    sigma_x.resize(2, 2);
    sigma_y.resize(2, 2);
    sigma_z.resize(2, 2);
    Value_mat.resize(2,2);

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

        // * +x direction Neighbor
        m = Coordinates_.neigh(l, 0); //+x neighbour cell
        mx_pos = Coordinates_.indx_cellwise(m);
        my_pos = Coordinates_.indy_cellwise(m);

        for (int spin = 0; spin < 2; spin++)
        {
            for (int orb1 = 0; orb1 < n_orbs_; orb1++)
            {
                for (int orb2 = 0; orb2 < n_orbs_; orb2++)
                {
                    if (Parameters_.hopping_NN_X(orb1, orb2) != 0.0)
                    {
                        a = Coordinates_.Nbasis(lx_pos, ly_pos, orb1) + ncells_ * n_orbs_ * spin;
                        b = Coordinates_.Nbasis(mx_pos, my_pos, orb2) + ncells_ * n_orbs_ * spin;
                        assert(a != b);
                        if (a != b)
                        {
                            HTB_(b, a) = complex<double>(1.0 * Parameters_.hopping_NN_X(orb1, orb2), 0.0);
                            HTB_(a, b) = conj(HTB_(b, a));
                        }
                    }
                }
            }
        }

        // * +y direction Neighbor
        m = Coordinates_.neigh(l, 2); //+y neighbour cell
        mx_pos = Coordinates_.indx_cellwise(m);
        my_pos = Coordinates_.indy_cellwise(m);

        for (int spin = 0; spin < 2; spin++)
        {
            for (int orb1 = 0; orb1 < n_orbs_; orb1++)
            {
                for (int orb2 = 0; orb2 < n_orbs_; orb2++)
                {
                    if (Parameters_.hopping_NN_Y(orb1, orb2) != 0.0)
                    {

                        a = Coordinates_.Nbasis(lx_pos, ly_pos, orb1) + ncells_ * n_orbs_ * spin;
                        b = Coordinates_.Nbasis(mx_pos, my_pos, orb2) + ncells_ * n_orbs_ * spin;
                        assert(a != b);
                        if (a != b)
                        {
                            HTB_(b, a) = complex<double>(1.0 * Parameters_.hopping_NN_Y(orb1, orb2), 0.0);
                            HTB_(a, b) = conj(HTB_(b, a));
                        }
                    }
                }
            }
        }

        //intra unit cell A(0)<-->B(1)
        for (int spin = 0; spin < 2; spin++)
        {

            //A(0)<-->B(1)
            a = Coordinates_.Nbasis(lx_pos, ly_pos, 0) + ncells_ * n_orbs_ * spin;
            b = Coordinates_.Nbasis(lx_pos, ly_pos, 1) + ncells_ * n_orbs_ * spin;
            assert(a != b);
            HTB_(b, a) = complex<double>(1.0 * Parameters_.hopping_intracell, 0.0);
            HTB_(a, b) = conj(HTB_(b, a));

            //A(1)<-->C(2)
            a = Coordinates_.Nbasis(lx_pos, ly_pos, 1) + ncells_ * n_orbs_ * spin;
            b = Coordinates_.Nbasis(lx_pos, ly_pos, 2) + ncells_ * n_orbs_ * spin;
            assert(a != b);
            HTB_(b, a) = complex<double>(1.0 * Parameters_.hopping_intracell, 0.0);
            HTB_(a, b) = conj(HTB_(b, a));
        }
    }





    //RASHBA SOC
    // Rashba SOC (Strictly for Dice lattice)
    int X_,Y_;
    X_=0;
    Y_=1;
    int lx_new, ly_new;

    //6 bond vectors
    Mat_2_doub D_;
    D_.resize(6);
    for(int i=0;i<6;i++){
        D_[i].resize(2); //z component is never used :)
    }


    D_[0][0]=-0.5; D_[0][1]=(sqrt(3.0)/2.0);
    D_[3][0]=-0.5; D_[3][1]=(sqrt(3.0)/2.0);

    D_[1][0]=1.0; D_[1][1]=0.0;
    D_[4][0]=1.0; D_[4][1]=0.0;

    D_[2][0]=-0.5; D_[2][1]=-1.0*(sqrt(3.0)/2.0);
    D_[5][0]=-0.5; D_[5][1]=-1.0*(sqrt(3.0)/2.0);




    for(int n=0;n<ncells_;n++){
        //cout<<"n = "<<n<<endl;
        lx_pos = Coordinates_.indx_cellwise(n);
        ly_pos = Coordinates_.indy_cellwise(n);
        for(int spin1=0;spin1<2;spin1++){
            a = Coordinates_.Nbasis(lx_pos, ly_pos, 1) + ncells_ * n_orbs_ * spin1;

            for(int spin2=0;spin2<2;spin2++){

                //0 : 1,r to 0,r+e1
                lx_new= (lx_pos+1+lx_)%lx_;
                ly_new= (ly_pos+0+ly_)%ly_;
                b = Coordinates_.Nbasis(lx_new, ly_new, 0) + ncells_ * n_orbs_ * spin2;
                for(int r=0;r<2;r++){
                    for(int c=0;c<2;c++){
                        Value_mat(r,c) = D_[0][X_]*sigma_x(r,c) + D_[0][Y_]*sigma_y(r,c);
                    }
                }

                assert(a != b);
                HTB_(b, a) += -1.0 * Parameters_.lambda_RSOC * Value_mat(spin2,spin1)*iota_complex;
                HTB_(a, b) = conj(HTB_(b, a));

                //1 : 1,r to 2,r
                lx_new= (lx_pos+0+lx_)%lx_;
                ly_new= (ly_pos+0+ly_)%ly_;
                b = Coordinates_.Nbasis(lx_new, ly_new, 2) + ncells_ * n_orbs_ * spin2;
                for(int r=0;r<2;r++){
                    for(int c=0;c<2;c++){
                        Value_mat(r,c) = D_[1][X_]*sigma_x(r,c) + D_[1][Y_]*sigma_y(r,c);
                    }
                }

                assert(a != b);
                HTB_(b, a) += -1.0 * Parameters_.lambda_RSOC * Value_mat(spin2,spin1)*iota_complex;
                HTB_(a, b) = conj(HTB_(b, a));


                //2 : 1,r to 0,r+e2
                lx_new= (lx_pos+0+lx_)%lx_;
                ly_new= (ly_pos+1+ly_)%ly_;
                b = Coordinates_.Nbasis(lx_new, ly_new, 0) + ncells_ * n_orbs_ * spin2;
                for(int r=0;r<2;r++){
                    for(int c=0;c<2;c++){
                        Value_mat(r,c) = D_[2][X_]*sigma_x(r,c) + D_[2][Y_]*sigma_y(r,c);
                    }
                }

                assert(a != b);
                HTB_(b, a) += -1.0 * Parameters_.lambda_RSOC * Value_mat(spin2,spin1)*iota_complex;
                HTB_(a, b) = conj(HTB_(b, a));


                //3 : 1,r to 2,r-e1
                lx_new= (lx_pos-1+lx_)%lx_;
                ly_new= (ly_pos+0+ly_)%ly_;
                b = Coordinates_.Nbasis(lx_new, ly_new, 2) + ncells_ * n_orbs_ * spin2;
                for(int r=0;r<2;r++){
                    for(int c=0;c<2;c++){
                        Value_mat(r,c) = D_[3][X_]*sigma_x(r,c) + D_[3][Y_]*sigma_y(r,c);
                    }
                }

                assert(a != b);
                HTB_(b, a) += -1.0 * Parameters_.lambda_RSOC * Value_mat(spin2,spin1)*iota_complex;
                HTB_(a, b) = conj(HTB_(b, a));


                //4 : 1,r to 0,r
                lx_new= (lx_pos+0+lx_)%lx_;
                ly_new= (ly_pos+0+ly_)%ly_;
                b = Coordinates_.Nbasis(lx_new, ly_new, 0) + ncells_ * n_orbs_ * spin2;
                for(int r=0;r<2;r++){
                    for(int c=0;c<2;c++){
                        Value_mat(r,c) = D_[4][X_]*sigma_x(r,c) + D_[4][Y_]*sigma_y(r,c);
                    }
                }

                assert(a != b);
                HTB_(b, a) += -1.0 * Parameters_.lambda_RSOC * Value_mat(spin2,spin1)*iota_complex;
                HTB_(a, b) = conj(HTB_(b, a));


                //5 : 1,r to 2,r-e2
                lx_new= (lx_pos+0+lx_)%lx_;
                ly_new= (ly_pos-1+ly_)%ly_;
                b = Coordinates_.Nbasis(lx_new, ly_new, 2) + ncells_ * n_orbs_ * spin2;
                for(int r=0;r<2;r++){
                    for(int c=0;c<2;c++){
                        Value_mat(r,c) = D_[5][X_]*sigma_x(r,c) + D_[5][Y_]*sigma_y(r,c);
                    }
                }

                assert(a != b);
                HTB_(b, a) += -1.0 * Parameters_.lambda_RSOC * Value_mat(spin2,spin1)*iota_complex;
                HTB_(a, b) = conj(HTB_(b, a));

            }
        }
    }



} // ----------


void Connections_DL::Print_Hopping(){

    string fileout=Parameters_.File_Hoppings;
    ofstream file_Hopping_out(fileout.c_str());

    file_Hopping_out<<"#site_i spin_i site_j spin_j Hopping[site_i,spin_i][site_j,spin_j]"<<endl;

    int index_i, index_j;
    for(int i=0;i<ncells_*3;i++){
        for(int spin_i=0;spin_i<2;spin_i++){
            index_i=i+spin_i*ncells_*3;

            for(int j=0;j<ncells_*3;j++){
                for(int spin_j=0;spin_j<2;spin_j++){
                    index_j=j+spin_j*ncells_*3;

                    file_Hopping_out<<i<<"  "<<spin_i<<"  "<<j<<"  "<<spin_j<<"  "<<HTB_(index_i,index_j)<<endl;

                }
            }
        }
    }

}

void Connections_DL::Print_LongRangeInt(){
    string fileout=Parameters_.File_LongRange_Ints;
    ofstream file_out(fileout.c_str());

    file_out<<"#i j Uij"<<endl;

}

void Connections_DL::Print_Spin_resolved_OnsiteE(){

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
