#include "Parameters_LL.h"
#include "Coordinates_LL.h"
#include "Connections_LL.h"

#ifndef Observables_LL_H
#define Observables_LL_H

#define PI acos(-1.0)

class Observables_LL
{
public:
    Observables_LL(Parameters_LL &Parameters__, Coordinates_LL &Coordinates__,
                Connections_LL &Connections__)
        : Parameters_(Parameters__), Coordinates_(Coordinates__),
          Connections_(Connections__),
          lx_(Parameters_.lx), ly_(Parameters_.ly), ncells_(Coordinates_.ncells_), nbasis_(Coordinates_.nbasis_), n_orbs_(Coordinates_.n_orbs_)
    {
        Initialize();
    }

    void Initialize();
    void Calculate_Akw();
    void Calculate_Nw();
    void RealSpaceLocal_ChernNumber();
    double Lorentzian(double x, double brd);
    void Calculate_OrbResolved_Nw();
    double fermi_function(int n);
    double fermi_function(int n, double mu_);
    void calculate_quantum_SiSj();
    void calculate_local_density();
    void Hall_conductance();
    void Create_Current_Oprs();


    Matrix<complex<double>> Ham_;
    vector<double> eigs_;

    double BandWidth;
    Matrix<complex<double>> SiSjQ_;
    Parameters_LL &Parameters_;
    Coordinates_LL &Coordinates_;
    Connections_LL &Connections_;
    int lx_, ly_, ncells_, nbasis_;
    int n_orbs_;
    Matrix<double> SiSj_;
    vector<double> sx_, sy_, sz_;

    Mat_2_doub local_density; 

    Matrix<complex<double>> quantum_SiSjQ_;
    Matrix<complex<double>> quantum_SiSj_;
    Matrix<complex<double>> J_KE_e1, J_KE_e2, J_RSOC_e1, J_RSOC_e2;

};
/*
 * ***********
 *  Functions in Class Observables_LL ------
 *  ***********
*/




void Observables_LL::Hall_conductance(){

    double hall_cond=0.0;
    double eps_temp =0.00000001;

    string fileout_="sigmaxy_vs_mu.txt";
    ofstream fileout(fileout_.c_str());

    for(int m=0;m<6*ncells_;m++){
        for(int n=0;n<6*ncells_;n++){
              if(abs(eigs_[m]-eigs_[n])>=eps_temp){
                hall_cond += ((2.0*PI)/(3.0*ncells_))*(fermi_function(m)-fermi_function(n))*
                        (1.0/( (Parameters_.eta*Parameters_.eta) + ((eigs_[n]-eigs_[m])*(eigs_[n]-eigs_[m]))   ))*
                        ((1.0*J_KE_e1(m,n) + 1.0*J_RSOC_e1(m,n))*(1.0*J_KE_e2(n,m) + 1.0*J_RSOC_e2(n,m))).imag();
              }
        }
    }

    cout<<"Hall Conductance = "<<hall_cond<<endl;


    double mu=eigs_[0]-5.0;
    while(mu<eigs_[eigs_.size()-1]+5.0){
    hall_cond=0.0;

    for(int m=0;m<6*ncells_;m++){
        for(int n=0;n<6*ncells_;n++){
              if(abs(eigs_[m]-eigs_[n])>=eps_temp){
                hall_cond += ((2.0*PI)/(3.0*ncells_))*(fermi_function(m,mu)-fermi_function(n,mu))*
                        (1.0/( (Parameters_.eta*Parameters_.eta) + ((eigs_[n]-eigs_[m])*(eigs_[n]-eigs_[m]))   ))*
                        ((1.0*J_KE_e1(m,n) + 1.0*J_RSOC_e1(m,n))*(1.0*J_KE_e2(n,m) + 1.0*J_RSOC_e2(n,m))).imag();
              }
        }
    }
    fileout<<mu<<"  "<<hall_cond<<endl;
    mu=mu+0.01;
    }

}


void Observables_LL::Create_Current_Oprs(){

    J_KE_e1.resize(ncells_*6, ncells_*6);
    J_RSOC_e1.resize(ncells_*6, ncells_*6);
    J_KE_e2.resize(ncells_*6, ncells_*6);
    J_RSOC_e2.resize(ncells_*6, ncells_*6);

    double Y_conv=1.0;

    //Convention used
    //orb=0=A
    //orb=1=B
    //orb=2=C


    Matrix<complex <double>> sigma_x, sigma_y, sigma_z, Value_mat;
    sigma_x.resize(2, 2);
    sigma_y.resize(2, 2);
    sigma_z.resize(2, 2);
    Value_mat.resize(2,2);

    //X
    sigma_x(0, 0) = 0.0;
    sigma_x(0, 1) = 1.0;
    sigma_x(1, 0) = 1.0;
    sigma_x(1, 1) = 0.0;

    //y
    sigma_y(0, 0) = 0.0;
    sigma_y(0, 1) = -1.0*iota_complex;
    sigma_y(1, 0) = iota_complex;
    sigma_y(1, 1) = 0.0;

    //Z
    sigma_z(0, 0) = 1.0;
    sigma_z(0, 1) = 0.0;
    sigma_z(1, 0) = 0.0;
    sigma_z(1, 1) = -1.0;


    int l, m, a, b;
    int lx_pos, ly_pos;
    int mx_pos, my_pos;

   // HTB_.fill(0.0);

    for (int p = 0; p < ncells_*6; p++){
        for (int n = 0; n < ncells_*6; n++){

            J_KE_e1(p,n)=0.0;
            J_RSOC_e1(p,n)=0.0;
            J_KE_e2(p,n)=0.0;
            J_RSOC_e2(p,n)=0.0;


            //K.E.
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
                                    J_KE_e1(p,n) += iota_complex*complex<double>(1.0 * Parameters_.hopping_NN_X(orb1, orb2), 0.0)*Ham_(a,n)*conj(Ham_(b,p));
                                    J_KE_e1(p,n) -= iota_complex*complex<double>(1.0 * Parameters_.hopping_NN_X(orb1, orb2), 0.0)*Ham_(b,n)*conj(Ham_(a,p));
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
                                    J_KE_e2(p,n) += Y_conv*iota_complex*complex<double>(1.0 * Parameters_.hopping_NN_Y(orb1, orb2), 0.0)*Ham_(a,n)*conj(Ham_(b,p));
                                    J_KE_e2(p,n) -= Y_conv*iota_complex*complex<double>(1.0 * Parameters_.hopping_NN_Y(orb1, orb2), 0.0)*Ham_(b,n)*conj(Ham_(a,p));
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
                    J_KE_e1(p,n) += iota_complex*complex<double>(1.0 * Parameters_.hopping_intracell, 0.0)*Ham_(a,n)*conj(Ham_(b,p));
                    J_KE_e1(p,n) -= iota_complex*complex<double>(1.0 * Parameters_.hopping_intracell, 0.0)*Ham_(b,n)*conj(Ham_(a,p));


                    //A(0)<-->C(2)
                    a = Coordinates_.Nbasis(lx_pos, ly_pos, 0) + ncells_ * n_orbs_ * spin;
                    b = Coordinates_.Nbasis(lx_pos, ly_pos, 2) + ncells_ * n_orbs_ * spin;
                    assert(a != b);
                    J_KE_e2(p,n) += Y_conv*iota_complex*complex<double>(1.0 * Parameters_.hopping_intracell, 0.0)*Ham_(a,n)*conj(Ham_(b,p));
                    J_KE_e2(p,n) -= Y_conv*iota_complex*complex<double>(1.0 * Parameters_.hopping_intracell, 0.0)*Ham_(b,n)*conj(Ham_(a,p));
                }
            }




            //R-SOC
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
                            for (int qp = 0; qp < 2; qp++)
                            {
                                for (int q = 0; q < 2; q++)
                                {
                                    rashba_mat(qp, q) =
                                            sigma_x(qp, q) * bond_vector[1] -
                                            sigma_y(qp, q) * bond_vector[0];
                                }
                            }

                            for (int spin_i = 0; spin_i < 2; spin_i++)
                            {
                                for (int spin_j = 0; spin_j < 2; spin_j++)
                                {

                                    col_index = Coordinates_.Nbasis(ix, iy, alpha) + spin_i*lx_*ly_*3;
                                    row_index = Coordinates_.Nbasis(ix, iy, beta) + spin_j*lx_*ly_*3;

                                    if(alpha==0 && beta==1){
                                        J_RSOC_e1(p,n) += -one_complex*Parameters_.lambda_RSOC * rashba_mat(spin_j, spin_i)*Ham_(col_index,n)*conj(Ham_(row_index,p));
                                        J_RSOC_e1(p,n) += -one_complex*Parameters_.lambda_RSOC * (rashba_mat(spin_i, spin_j))*Ham_(row_index,n)*conj(Ham_(col_index,p));
                                    }
                                    if(alpha==0 && beta==2){
                                        J_RSOC_e2(p,n) += +one_complex*Parameters_.lambda_RSOC * rashba_mat(spin_j, spin_i)*Ham_(col_index,n)*conj(Ham_(row_index,p));
                                        J_RSOC_e2(p,n) += +one_complex*Parameters_.lambda_RSOC * (rashba_mat(spin_i, spin_j))*Ham_(row_index,n)*conj(Ham_(col_index,p));
                                    }
                                    //HTB_(row_index, col_index) += iota_complex * Parameters_.lambda_RSOC * rashba_mat(spin_j, spin_i);
                                    //cout<<col_index<<"  "<<row_index<<"  "<< iota_complex*Parameters_.lambda_RSOC * rashba_mat(spin_j, spin_i)<<endl;
                                }
                            }
                        }

                        // inter +X
                        ix_new = (ix + 1 + lx_) % lx_;
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
                                    for (int qp = 0; qp < 2; qp++)
                                    {
                                        for (int q = 0; q < 2; q++)
                                        {
                                            rashba_mat(qp, q) =
                                                    sigma_x(qp, q) * bond_vector[1] -
                                                    sigma_y(qp, q) * bond_vector[0];
                                        }
                                    }

                                    col_index = Coordinates_.Nbasis(ix, iy, alpha) + spin_i*lx_*ly_*3;
                                    row_index = Coordinates_.Nbasis(ix_new, iy_new, beta) + spin_j*lx_*ly_*3;

                                    J_RSOC_e1(p,n) += -one_complex*Parameters_.lambda_RSOC * rashba_mat(spin_j, spin_i)*Ham_(col_index,n)*conj(Ham_(row_index,p));
                                    J_RSOC_e1(p,n) += -one_complex*Parameters_.lambda_RSOC * rashba_mat(spin_i, spin_j)*Ham_(row_index,n)*conj(Ham_(col_index,p));

                                }
                            }
                        }

                        // inter +Y
                        ix_new = ix;
                        iy_new = (iy + 1 + ly_) % ly_;
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
                                    for (int qp = 0; qp < 2; qp++)
                                    {
                                        for (int q = 0; q < 2; q++)
                                        {
                                            rashba_mat(qp, q) =
                                                    sigma_x(qp, q) * bond_vector[1] -
                                                    sigma_y(qp, q) * bond_vector[0];
                                        }
                                    }

                                    col_index = Coordinates_.Nbasis(ix, iy, alpha) + spin_i*lx_*ly_*3;
                                    row_index = Coordinates_.Nbasis(ix_new, iy_new, beta) + spin_j*lx_*ly_*3;

                                    J_RSOC_e2(p,n) += +one_complex*Parameters_.lambda_RSOC * rashba_mat(spin_j, spin_i)*Ham_(col_index,n)*conj(Ham_(row_index,p));
                                    J_RSOC_e2(p,n) += +one_complex*Parameters_.lambda_RSOC * rashba_mat(spin_i, spin_j)*Ham_(row_index,n)*conj(Ham_(col_index,p));

                                    //cout<<"here 1"<<iota_complex*Parameters_.lambda_RSOC * rashba_mat(spin_j, spin_i)<<endl;
                                }
                            }
                        }
                    }
                }
            }

        }
        cout<<"Current opr row p  = "<<p <<" done"<<endl;
    }



}

void Observables_LL::RealSpaceLocal_ChernNumber(){

    Mat_2_Complex_doub Q, P, XQ, YP, Chern;

    Q.resize(ncells_*6);P.resize(ncells_*6);XQ.resize(ncells_*6);YP.resize(ncells_*6);
    for(int i=0;i<ncells_*6;i++){
        Q[i].resize(ncells_*6);
        P[i].resize(ncells_*6);
        XQ[i].resize(ncells_*6);
        YP[i].resize(ncells_*6);
    }


    Chern.resize(6); //alpha, sigma
    for(int i =0;i<6;i++){
        Chern[i].resize(ncells_);
    }


    for(int bi=0;bi<ncells_*6;bi++){
        for(int bj=0;bj<ncells_*6;bj++){
            Q[bi][bj]=0.0;
            P[bi][bj]=0.0;
            for(int n=0;n<eigs_.size();n++){
                if(eigs_[n]>Parameters_.mus){
                    Q[bi][bj] += Ham_(bi,n)*conj(Ham_(bj,n));
                }
                if(eigs_[n]<Parameters_.mus){
                    P[bi][bj] += Ham_(bi,n)*conj(Ham_(bj,n));
                }
            }
        }
    }


    int bk_;
    for(int bi=0;bi<ncells_*6;bi++){
        for(int bj=0;bj<ncells_*6;bj++){

            YP[bj][bi]=0.0;
            XQ[bi][bj]=0.0;
            for(int bk=0;bk<ncells_*3;bk++){
                for(int spin=0;spin<2;spin++){
                    bk_=bk + spin*ncells_*3;
                    YP[bj][bi] += P[bj][bk_]*(1.0*Coordinates_.indy_basiswise(bk)*Q[bk_][bi]);
                    XQ[bi][bj] += Q[bi][bk_]*(1.0*Coordinates_.indx_basiswise(bk)*P[bk_][bj]);
                }

            }
        }
    }


    int dof, i_, bi_;
    for(int orb=0;orb<3;orb++){
        for (int spin=0;spin<2;spin++){
            dof = orb + spin*3;

            string fileout = "Chern_" + to_string(orb) + "_" + to_string(spin) + ".txt" ;
            ofstream file_out(fileout.c_str());



            for(int iy=0;iy<ly_;iy++){
                for(int ix=0;ix<lx_;ix++){


                    i_=Coordinates_.Ncell(ix,iy);
                    bi_ =Coordinates_.Nbasis(ix, iy, orb) + spin*ncells_*3;

                    Chern[dof][i_] =0.0;
                    for(int bj=0;bj<ncells_*6;bj++){
                        Chern[dof][i_] += XQ[bi_][bj]*YP[bj][bi_];
                    }

                    file_out<<ix<<"    "<<iy<<"    "<<Chern[dof][i_].real()<<"    "<<Chern[dof][i_].imag()<<endl;

                }
                file_out<<endl;
            }
        }
    }





}

void Observables_LL::Calculate_Akw()
{

    //---------Read from input file-----------------------//
    string fileout = "Akw.txt" ;
    double omega_min, omega_max, d_omega;
    double eta = Parameters_.eta;
    omega_min = eigs_[0] - 1.0;
    omega_max = eigs_[eigs_.size() - 1] + 1.0;
    d_omega = Parameters_.domega;
    //---------------------------------------------------//

    int omega_index_max = int((omega_max - omega_min) / (d_omega));

    ofstream file_Akw_out(fileout.c_str());

    int c1;

    Mat_3_Complex_doub A_nk;
    A_nk.resize(6);
    for(int i=0;i<6;i++){
        A_nk[i].resize(eigs_.size());
        for(int n=0;n<eigs_.size();n++){
            A_nk[i][n].resize(Coordinates_.ncells_);

        }
    }


    complex<double> temp_doub;
    int local_dof;
    double kx_val, ky_val;
    //----- A_\alpha\sigma(n,k) = sum_{i} Psi_{i,alpha, sigma ;n}exp(iota*K_vec \cdot r_vec) --------
    for(int n=0;n<eigs_.size();n++){
        for(int kx=0;kx<lx_;kx++){
            kx_val = (2.0*PI*kx)/(1.0*lx_) ;
            for(int ky=0;ky<ly_;ky++){
                ky_val = (2.0*PI*ky)/(1.0*ly_) ;
                for(int orb=0;orb<3;orb++){
                    for(int spin=0;spin<2;spin++){
                        local_dof=orb + spin*(3);

                        temp_doub=0.0;
                        for(int ix=0;ix<lx_;ix++){
                            for(int iy=0;iy<ly_;iy++){
                                c1 = Coordinates_.Nbasis(ix, iy, orb) + spin*(lx_*ly_*n_orbs_);
                                temp_doub += Ham_(c1, n)*exp(-1.0*iota_complex* (kx_val*ix  + ky_val*iy));
                            }
                        }
                        A_nk[local_dof][n][Coordinates_.Ncell(kx,ky)]=temp_doub;
                    }
                }
            }
        }
        if(n%10==0){
            cout << "n = "<<n<<" done"<<endl;}

    }


    cout<< "A_{alpha,sigma}(n,k) is done"<<endl;


    Mat_3_Complex_doub A_kw;

    A_kw.resize(6);
    for(int i=0;i<6;i++){
        A_kw[i].resize(Coordinates_.ncells_);
        for(int n=0;n<Coordinates_.ncells_;n++){
            A_kw[i][n].resize(omega_index_max);
        }
    }



    complex<double> Nup_check(0, 0);
    complex<double> Ndn_check(0, 0);



    // A_alpha_sigma(k,w) = 1/N^2 \sum_{n} |A_alpha, sigma(n,k)|^2 delta(w-eps_{n})
    int k;
    for(int orb=0;orb<3;orb++){
        for(int spin=0;spin<2;spin++){
            local_dof=orb + spin*(3);

            for(int kx=0;kx<lx_;kx++){
                for(int ky=0;ky<ly_;ky++){
                    k=Coordinates_.Ncell(kx,ky);
                    for(int w_no=0;w_no<omega_index_max;w_no++){

                        A_kw[local_dof][k][w_no]=0.0;
                        for(int n=0;n<eigs_.size();n++){
                            A_kw[local_dof][k][w_no] += abs(A_nk[local_dof][n][k])*abs(A_nk[local_dof][n][k])*
                                    Lorentzian(omega_min + (w_no * d_omega) -eigs_[n], eta);

                        }
                        A_kw[local_dof][k][w_no] = A_kw[local_dof][k][w_no]*(1.0/(ncells_*ncells_));

                    }
                    if(k%10==0){
                        cout << "k = "<<k<<" done"<<endl;
                    }
                }
            }
        }
    }




    cout << "Nup_check = " << Nup_check << endl;
    cout << "Ndn_check = " << Ndn_check << endl;


    double kx, ky;
    int kx_i, ky_i;

    Mat_1_intpair k_path;
    k_path.clear();
    Mat_1_intpair k_path2;
    k_path2.clear();
    pair_int temp_pair;


    // ---k_path---------

    //--------\Gamma to X-----------------
    ky_i = 0;
    for (kx_i = 0; kx_i <= (Parameters_.lx / 2); kx_i++)
    {
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------

    //--------X to M-----------------
    kx_i = (Parameters_.lx / 2);
    for (ky_i = 1; ky_i <= (Parameters_.ly / 2); ky_i++)
    {
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------

    //--------M to \Gamma[with one extra point,
    //                  because in gnuplot use "set pm3d corners2color c1"
    //                  ]-----------------
    kx_i = (Parameters_.lx / 2) - 1;
    ky_i = (Parameters_.ly / 2) - 1;
    for (kx_i = (Parameters_.lx / 2) - 1; kx_i >= 0; kx_i--)
    {
        temp_pair.first = kx_i;
        temp_pair.second = kx_i;
        k_path.push_back(temp_pair);
    }

    temp_pair.first = 0;
    temp_pair.second = 0;
    k_path.push_back(temp_pair);

    //----------------------------------

    //----k_path done-------



    // ---k_path2---------

    //--------\Gamma to 2X-----------------
    ky_i = 0;
    for (kx_i = 0; kx_i <= (Parameters_.lx)-1; kx_i++)
    {
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path2.push_back(temp_pair);
    }
    //----------------------------------

    //--------2X to 2M-----------------
    kx_i = (Parameters_.lx)-1;
    for (ky_i = 1; ky_i <= (Parameters_.ly)-1; ky_i++)
    {
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path2.push_back(temp_pair);
    }
    //----------------------------------

    //--------2M to \Gamma[with one extra point,
    //                  because in gnuplor use "set pm3d corners2color c1"
    //                  ]-----------------
    kx_i = (Parameters_.lx) - 1;
    ky_i = (Parameters_.ly) - 1;
    for (kx_i = (Parameters_.lx) - 1; kx_i >= -1; kx_i--)
    {
        temp_pair.first = kx_i;
        temp_pair.second = kx_i;
        k_path2.push_back(temp_pair);
    }
    //----------------------------------

    //----k_path2 done-------



    file_Akw_out<<"#k_point     kx      ky    k_index      omega_val    omega_ind        Akw[orb=0,spin=0]     Akw[orb=0,spin=1]      Akw[orb=1,spin=0]     Akw[orb=1,spin=1]       Akw[orb=2,spin=0]     Akw[orb=2,spin=1]"<<endl;
    for (int k_point = 0; k_point < k_path.size(); k_point++)
    {

        kx_i = k_path[k_point].first;
        ky_i = k_path[k_point].second;
        kx = (2.0 * PI * kx_i) / (1.0 * Parameters_.lx);
        ky = (2.0 * PI * ky_i) / (1.0 * Parameters_.ly);

        for (int omega_ind = 0; omega_ind < omega_index_max; omega_ind++)
        {

            //Use 1:6:7----for gnuplot
            file_Akw_out << k_point << "   " << kx_i << "   " << ky_i << "   " << Coordinates_.Ncell(kx_i, ky_i) << "    " << omega_min + (d_omega * omega_ind) << "   " << omega_ind << "    ";


            for(int orb=0;orb<3;orb++){
                for(int spin=0;spin<2;spin++){
                    local_dof = orb + spin*3;
                    file_Akw_out << A_kw[local_dof][Coordinates_.Ncell(kx_i, ky_i)][omega_ind].real()<<"     ";
                }
            }
            file_Akw_out << endl;
        }
        file_Akw_out << endl;
    }

}

void Observables_LL::Calculate_Nw()
{

    //---------Read from input file-----------------------//
    string fileout = "Nw_total_LL.txt" ;
    double omega_min, omega_max, d_omega;
    double eta = Parameters_.eta;
    omega_min = eigs_[0] - 5.0;
    omega_max = eigs_[eigs_.size() - 1] + 5.0;
    d_omega = Parameters_.domega;
    //---------------------------------------------------//

    int omega_index_max = int((omega_max - omega_min) / (d_omega));

    ofstream file_Nw_out(fileout.c_str());

    double temp_val;

    for (int omega_ind = 0; omega_ind < omega_index_max; omega_ind++)
    {

        temp_val = 0.0;

        for (int n = 0; n < Ham_.n_row(); n++)
        {

            temp_val +=
                    Lorentzian(omega_min + (omega_ind * d_omega) -eigs_[n], eta);
        }

        file_Nw_out << omega_min + (omega_ind * d_omega) << "     "
                    << (1.0 / eigs_.size()) * temp_val
                    << "     " << endl;
    }

    file_Nw_out << "#mu = " << Parameters_.mus << endl;

    ofstream spectrum;
    spectrum.open("spectrum_LL.txt");

    for (int i=0; i < eigs_.size(); i++){
        spectrum << i << "\t" << eigs_[i] << endl;
    }
    spectrum.close();
}

void Observables_LL::Calculate_OrbResolved_Nw()
{

    //---------Read from input file-----------------------//
    string fileout = "NwOrbResolved_total.txt";
    double omega_min, omega_max, d_omega;

    double eta = Parameters_.eta;
    omega_min = eigs_[0] - 5.0;
    omega_max = eigs_[eigs_.size() - 1] + 5.0;
    d_omega = Parameters_.domega;
    //---------------------------------------------------//

    int omega_index_max = int((omega_max - omega_min) / (d_omega));

    ofstream file_Nw_out(fileout.c_str());

    double temp_val;

    for (int omega_ind = 0; omega_ind < omega_index_max; omega_ind++)
    {

        temp_val = 0.0;

        for (int n = 0; n < Ham_.n_row(); n++)
        {

            temp_val += Lorentzian(omega_min + (omega_ind * d_omega) - eigs_[n], eta);
        }

        file_Nw_out << omega_min + (omega_ind * d_omega) << "     " << temp_val << "     " << endl;
    }

    file_Nw_out << "#mu = " << Parameters_.mus << endl;
}


double Observables_LL::Lorentzian(double x, double brd)
{
    double temp;

    temp = (1.0 / PI) * ((brd / 2.0) / ((x * x) + ((brd * brd) / 4.0)));

    return temp;
}


double Observables_LL::fermi_function(int n)
{
    double value;
    value = 1.0 / (exp(Parameters_.beta * (eigs_[n] - Parameters_.mus)) + 1.0);
    return value;

}

double Observables_LL::fermi_function(int n, double mu_)
{
    double value;
    value = 1.0 / (exp(Parameters_.beta * (eigs_[n] - mu_)) + 1.0);
    return value;

}


void Observables_LL::calculate_quantum_SiSj()
{
    /*
    Matrix<complex<double>> F_u_u;
    Matrix<complex<double>> F_d_d;
    Matrix<complex<double>> F_u_d;
    Matrix<complex<double>> F_d_u;
    Matrix<complex<double>> omF_u_u;
    Matrix<complex<double>> omF_d_d;
    Matrix<complex<double>> omF_u_d;
    Matrix<complex<double>> omF_d_u;
    int nx, ny;
    int jx, jy;
    F_u_u.resize(ns_, ns_);
    F_d_d.resize(ns_, ns_);
    F_u_d.resize(ns_, ns_);
    F_d_u.resize(ns_, ns_);
    omF_u_u.resize(ns_, ns_);
    omF_d_d.resize(ns_, ns_);
    omF_u_d.resize(ns_, ns_);
    omF_d_u.resize(ns_, ns_);
    for (int i = 0; i < ns_; i++)
    {
        for (int j = 0; j < ns_; j++)
        {
            for (int n = 0; n < Hamiltonian_.eigs_.size(); n++)
            {
                F_u_u(i, j) += (conj(Hamiltonian_.Ham_(i, n)) * Hamiltonian_.Ham_(j, n) * fermi_function(n));
                F_d_d(i, j) += (conj(Hamiltonian_.Ham_(i + ns_, n)) * Hamiltonian_.Ham_(j + ns_, n) * fermi_function(n));
                F_u_d(i, j) += (conj(Hamiltonian_.Ham_(i, n)) * Hamiltonian_.Ham_(j + ns_, n) * fermi_function(n));
                F_d_u(i, j) += (conj(Hamiltonian_.Ham_(i + ns_, n)) * Hamiltonian_.Ham_(j, n) * fermi_function(n));
                omF_u_u(i, j) += (conj(Hamiltonian_.Ham_(i, n)) * Hamiltonian_.Ham_(j, n) * (1.0 - fermi_function(n)));
                omF_d_d(i, j) += (conj(Hamiltonian_.Ham_(i + ns_, n)) * Hamiltonian_.Ham_(j + ns_, n) * (1.0 - fermi_function(n)));
                omF_u_d(i, j) += (conj(Hamiltonian_.Ham_(i, n)) * Hamiltonian_.Ham_(j + ns_, n) * (1.0 - fermi_function(n)));
                omF_d_u(i, j) += (conj(Hamiltonian_.Ham_(i + ns_, n)) * Hamiltonian_.Ham_(j, n) * (1.0 - fermi_function(n)));
            }
        }
    }

    int i_;
    for (int ix = 0; ix < lx_; ix++)
    {
        for (int iy = 0; iy < ly_; iy++)
        {
            // i = Coordinates_.Nc(ix, iy);

            quantum_SiSj_(ix, iy) = 0.0;
            for (int j = 0; j < ns_; j++)
            {
                jx = Coordinates_.indx(j);
                jy = Coordinates_.indy(j);
                nx = (jx + ix) % lx_;
                ny = (jy + iy) % ly_;
                i_ = Coordinates_.Nc(nx, ny);
                quantum_SiSj_(ix, iy) += (
                0.25*(F_u_u(i_, i_) * F_u_u(j, j) + F_u_u(i_, j) * omF_u_u(j, i_)
                    - ( F_u_u(i_, i_) * F_d_d(j, j) + F_u_d(i_, j) * omF_d_u(j, i_) )
                    - ( F_d_d(i_, i_) * F_u_u(j, j) + F_d_u(i_, j) * omF_u_d(j, i_) )
                    + F_d_d(i_, i_) * F_d_d(j, j) + F_d_d(i_, j) * omF_d_d(j, i_))
                    + 0.5 * (F_u_d(i_, i_) * F_d_u(j, j) + F_u_u(i_, j) * omF_d_d(j, i_))
                    + 0.5 * (F_d_u(i_, i_) * F_u_d(j, j) + F_d_d(i_, j) * omF_u_u(j, i_))
                    ).real();
            }
            quantum_SiSj_(ix, iy) /= (ns_ * 1.0);
        }
    }

    //Fourier transform
    double phase, Cos_ij, Sin_ij;
    for (int qx = 0; qx < lx_; qx++)
    {
        for (int qy = 0; qy < ly_; qy++)
        {
            quantum_SiSjQ_(qx, qy) = complex<double>(0.0, 0.0);
            for (int xr = 0; xr < lx_; xr++)
            {
                for (int yr = 0; yr < ly_; yr++)
                {
                    phase = 2.0 * Parameters_.pi * (double(qx * xr) / double(lx_) + double(qy * yr) / double(ly_));
                    Cos_ij = cos(phase);
                    Sin_ij = sin(phase);
                    quantum_SiSjQ_(qx, qy) += quantum_SiSj_(xr, yr) * complex<double>(Cos_ij, Sin_ij);
                }
            }
            quantum_SiSjQ_(qx, qy) *= double(1.0 / (lx_ * ly_));
            //cout << qx << " "<< qy<< " "<<  SiSjQ_(qx,qy) << endl;
        }
    }

    */
}

void Observables_LL::calculate_local_density()
{
    int c1;
    complex<double> value = zero_complex;
    // cout <<"Parameter mus="<< Parameters_.mus<<endl;
    for (int i = 0; i < nbasis_; i++)
    {
        for (int sigma = 0; sigma < 2; sigma++)
        {
            local_density[i][sigma] = 0.0;
            c1 = i + (sigma * nbasis_);
            for (int n = 0; n < eigs_.size(); n++)
            {
                local_density[i][sigma] += (conj(Ham_(c1, n)) * Ham_(c1, n) * fermi_function(n)).real();
            }

            // value += (conj(Ham_(c1, 1)) * Ham_(c1, 1));
        }
    }
}


void Observables_LL::Initialize()
{

    complex<double> zero(0.0, 0.0);
    int space = 2 * ncells_ * n_orbs_;
    sx_.resize(ncells_ * n_orbs_);
    sy_.resize(ncells_ * n_orbs_);
    sz_.resize(ncells_ * n_orbs_);

    local_density.resize(nbasis_);

    for (int i = 0; i < local_density.size(); i++)
    {
        local_density[i].resize(2);      
    }

    quantum_SiSj_.resize(lx_, ly_);
    quantum_SiSjQ_.resize(lx_, ly_);

} // ----------


#endif // Observables_LL_H
