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
    void calculate_quantum_SiSj();
    void calculate_local_density();


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

};
/*
 * ***********
 *  Functions in Class Observables_LL ------
 *  ***********
*/


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
