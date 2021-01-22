#include <algorithm>
#include <functional>
#include <math.h>
#include "tensor_type.h"
#include "Parameters_LL.h"
#include "Coordinates_LL.h"
#include "../../Matrix.h"
#define PI acos(-1.0)

#ifndef Connections_LL_class
#define Connections_LL_class


class Connections_LL
{
public:
    Connections_LL(Parameters_LL &Parameters__, Coordinates_LL &Coordinates__)
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



    Parameters_LL &Parameters_;
    Coordinates_LL &Coordinates_;
    int lx_, ly_, ncells_, n_orbs_;
    Matrix<complex<double>> HTB_;
    Matrix<complex<double>> Ham_;
    Matrix<double> Tx, Ty, Tpxpy, Tpxmy;
    Matrix<double> IntraCell_Hopp, InterCell_px_Hopp, InterCell_py_Hopp, InterCell_pxmy_Hopp;

};


void Connections_LL::Initialize()
{


    ly_ = Parameters_.ly;
    lx_ = Parameters_.lx;
    ncells_ = lx_ * ly_;
    n_orbs_ = Parameters_.n_orbs;
    int space = 2 * ncells_ * n_orbs_;

    HTB_.resize(space, space);
    Ham_.resize(space, space);

} // ----------


void Connections_LL::InteractionsCreate()
{

    Ham_ = HTB_;
    // Ham_.print();

} // ----------


void Connections_LL::Check_up_down_symmetry()

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

void Connections_LL::Check_Hermiticity()

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


void Connections_LL::HTBCreate()
{

    //Convention used
    //orb=0=A
    //orb=1=B
    //orb=2=C


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

            //A(0)<-->C(2)
            a = Coordinates_.Nbasis(lx_pos, ly_pos, 0) + ncells_ * n_orbs_ * spin;
            b = Coordinates_.Nbasis(lx_pos, ly_pos, 2) + ncells_ * n_orbs_ * spin;
            assert(a != b);
            HTB_(b, a) = complex<double>(1.0 * Parameters_.hopping_intracell, 0.0);
            HTB_(a, b) = conj(HTB_(b, a));
        }
    }





    //RASHBA SOC
    // Rashba SOC (Strictly for Lieb lattice)
    vector<double> bond_vector;
    bond_vector.resize(3);
    Matrix<complex <double>> rashba_mat; //CONVENTION OF X,Y is Fig-1(a) of Physics Letters A 381 (2017) 944-948
    rashba_mat.resize(2, 2);
    int i, col_index, row_index, ix_new, iy_new, i_new;
    for (int ix = 0; ix < lx_; ix++)
    {
        for (int iy = 0; iy < ly_; iy++)
        {
            i = Coordinates_.Ncell(ix, iy);

            for (int alpha = 0; alpha < 3; alpha++)
            {

                //intra unit cell
                for (int beta = 0; beta < 3; beta++)
                {
                    if (alpha == 0 && beta == 1)
                    {
                        bond_vector[0] = 1.0;
                        bond_vector[1] = 0.0;
                        bond_vector[2] = 0.0;
                    }
                    else if (alpha == 1 && beta == 0)
                    {
                        bond_vector[0] = -1.0;
                        bond_vector[1] = 0.0;
                        bond_vector[2] = 0.0;
                    }
                    else if (alpha == 0 && beta == 2)
                    {
                        bond_vector[0] = 0.0;
                        bond_vector[1] = -1.0;
                        bond_vector[2] = 0.0;
                    }
                    else if (alpha == 2 && beta == 0)
                    {
                        bond_vector[0] = 0.0;
                        bond_vector[1] = 1.0;
                        bond_vector[2] = 0.0;
                    }
                    else
                    {
                        bond_vector[0] = 0.0;
                        bond_vector[1] = 0.0;
                        bond_vector[2] = 0.0;
                    }

                    // sigma_x d_y - d_x sigma_y
                    for (int p = 0; p < 2; p++)
                    {
                        for (int q = 0; q < 2; q++)
                        {
                            rashba_mat(p, q) =
                                    sigma_x(p, q) * bond_vector[1] -
                                    sigma_y(p, q) * bond_vector[0];
                        }
                    }

                    for (int spin_i = 0; spin_i < 2; spin_i++)
                    {
                        for (int spin_j = 0; spin_j < 2; spin_j++)
                        {

                            col_index = Coordinates_.Nbasis(ix, iy, alpha) + spin_i*lx_*ly_*3;
                            row_index = Coordinates_.Nbasis(ix, iy, beta) + spin_j*lx_*ly_*3;
                            HTB_(row_index, col_index) += iota_complex * Parameters_.lambda_RSOC * rashba_mat(spin_j, spin_i);
                            //cout<<col_index<<"  "<<row_index<<"  "<< iota_complex*Parameters_.lambda_RSOC * rashba_mat(spin_j, spin_i)<<endl;
                        }
                    }
                }

                // inter +X
                ix_new = (ix + 1) % lx_;
                iy_new = iy;
                i_new = iy_new + ly_ * ix_new;
                for (int beta = 0; beta < 3; beta++)
                {

                    if (alpha == 1 && beta == 0)
                    {
                        bond_vector[0] = 1.0;
                        bond_vector[1] = 0.0;
                        bond_vector[2] = 0.0;
                    }
                    else
                    {
                        bond_vector[0] = 0.0;
                        bond_vector[1] = 0.0;
                        bond_vector[2] = 0.0;
                    }
                    for (int spin_i = 0; spin_i < 2; spin_i++)
                    {
                        for (int spin_j = 0; spin_j < 2; spin_j++)
                        {

                            // sigma_x d_y - d_x sigma_y
                            for (int p = 0; p < 2; p++)
                            {
                                for (int q = 0; q < 2; q++)
                                {
                                    rashba_mat(p, q) =
                                            sigma_x(p, q) * bond_vector[1] -
                                            sigma_y(p, q) * bond_vector[0];
                                }
                            }

                            col_index = Coordinates_.Nbasis(ix, iy, alpha) + spin_i*lx_*ly_*3;

                            row_index = Coordinates_.Nbasis(ix_new, iy_new, beta) + spin_j*lx_*ly_*3;

                            HTB_(row_index, col_index) +=
                                    iota_complex *
                                    Parameters_.lambda_RSOC * rashba_mat(spin_j, spin_i);
                            HTB_(col_index, row_index) +=
                                    conj(iota_complex *
                                         Parameters_.lambda_RSOC * rashba_mat(spin_j, spin_i));
                        }
                    }
                }

                // inter +Y
                ix_new = ix;
                iy_new = (iy + 1) % ly_;
                i_new = iy_new + ly_ * ix_new;
                for (int beta = 0; beta < 3; beta++)
                {

                    if (alpha == 2 && beta == 0)
                    {
                        bond_vector[0] = 0.0;
                        bond_vector[1] = -1.0;
                        bond_vector[2] = 0.0;
                    }
                    else
                    {
                        bond_vector[0] = 0.0;
                        bond_vector[1] = 0.0;
                        bond_vector[2] = 0.0;
                    }
                    for (int spin_i = 0; spin_i < 2; spin_i++)
                    {
                        for (int spin_j = 0; spin_j < 2; spin_j++)
                        {

                            // sigma_x d_y - d_x sigma_y
                            for (int p = 0; p < 2; p++)
                            {
                                for (int q = 0; q < 2; q++)
                                {
                                    rashba_mat(p, q) =
                                            sigma_x(p, q) * bond_vector[1] -
                                            sigma_y(p, q) * bond_vector[0];
                                }
                            }

                            col_index = Coordinates_.Nbasis(ix, iy, alpha) + spin_i*lx_*ly_*3;


                            row_index = Coordinates_.Nbasis(ix_new, iy_new, beta) + spin_j*lx_*ly_*3;

                            HTB_(row_index, col_index) +=
                                    iota_complex *
                                    Parameters_.lambda_RSOC * rashba_mat(spin_j, spin_i);
                            HTB_(col_index, row_index) +=
                                    conj(iota_complex *
                                         Parameters_.lambda_RSOC * rashba_mat(spin_j, spin_i));
                            //cout<<"here 1"<<iota_complex*Parameters_.lambda_RSOC * rashba_mat(spin_j, spin_i)<<endl;
                        }
                    }
                }
            }
        }
    }

} // ----------


void Connections_LL::Print_Hopping(){

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

void Connections_LL::Print_LongRangeInt(){
    string fileout=Parameters_.File_LongRange_Ints;
    ofstream file_out(fileout.c_str());

    file_out<<"#i j Uij"<<endl;

}

void Connections_LL::Print_Spin_resolved_OnsiteE(){

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
