#include <algorithm>
#include <functional>
#include <math.h>
#include "Parameters_G2dLatticeNew.h"
#include "Coordinates_G2dLatticeNew.h"
#include "Connections_G2dLatticeNew.h"
#include "random"
#define PI acos(-1.0)
#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef Kspace_calculation_G2dLatticeNew_class
#define Kspace_calculation_G2dLatticeNew_class

extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
                         std::complex<double> *,int *, double *, int *);

extern "C" void dgesdd_ (char *, int *, int *, double *, int *, double *, double *, int *, double *, int*,
                         double *, int *, int *, int *);


extern "C" void zgesdd_ (char *, int *, int *, std::complex<double> *, int *, double *, std::complex<double> *, int *, std::complex<double> *, int*,
                         std::complex<double> *, int *, double * , int *, int *);

extern "C" void zgetri_(int *, std::complex<double> *, int *, int *, std::complex<double> *, int *, int * );
extern "C" void zgetrf_(int *, int* , std::complex<double> *, int* , int *, int *);



class Kspace_calculation_G2dLatticeNew
{
public:
    Kspace_calculation_G2dLatticeNew(Parameters_G2dLatticeNew &Parameters__, Coordinates_G2dLatticeNew &Coordinates__, Connections_G2dLatticeNew &Connections__, mt19937_64& Generator1__ )
        : Parameters_(Parameters__), Coordinates_(Coordinates__), Connections_(Connections__), Generator1_(Generator1__)

    {
        Initialize();
    }

    void Initialize();                                     //::DONE
    void Diagonalize(char option);
    void DiagonalizeGivenMat(char option, Matrix<complex<double>> & HamGiven_, vector<double> & EigsOut_);
    double random1();
    void SelfConsistency();
    void Create_Kspace_Spectrum();
    void Arranging_spectrum();
    double chemicalpotential(double Particles);
    void Get_new_OPs_and_error();
    void Get_Energies();
    void Get_Energies_new();
    void Get_Bands();
    void Get_minimum_distance_direction(int l,int m,int &r1_, int &r2_);
    complex<double> h_KE(int alpha, int gamma, int sigma, int alpha_p, int gamma_p, int sigma_p, int k1, int k2);
    complex<double> IntraCell_U(int alpha, int gamma, int alpha_p, int gamma_p);

    void Perform_SVD(Matrix<double> & A_, Matrix<double> & VT_, Matrix<double> & U_, vector<double> & Sigma_);
    void Perform_SVD_complex(Matrix<complex<double>> & A_, Matrix<complex<double>> & VT_, Matrix<complex<double>> & U_, vector<double> & Sigma_);
    void Update_OrderParameters_AndersonMixing(int iter);
    void Get_spin_resolved_local_densities();
    void Get_local_spins();
    void Get_Tau_Pseudospins();
    void Calculate_Nw();
    double Lorentzian(double x, double brd);
    void Calculate_ChernNumbers();
    void Create_Current_Oprs();
    void Create_Current_Oprs_Faster();
    void Hall_conductance();
    void Update_Total_Density();
    void Calculate_Akw();
    void Create_K_Path(string PathType, Mat_1_intpair & k_path);
    void Create_K_Path_NBZ(string PathType, Mat_1_intpair & k_path);
    void Get_Bare_Susceptibility_New();
    void Get_Bare_Susceptibility();
    void Get_RPA_Susceptibility();
    void Get_RPA_Susceptibility_old();
    void Calculate_RPA_Susc(Matrix<complex<double>> Chi_bare, Matrix<complex<double>> & Chi_RPA);
    void Inverse(Matrix<complex<double>> & A);
    int KroneckerDelta(int i, int j);
    void Connections_bw_MUC();
    void Create_InverseMapping();
    void Get_Bare_Susceptibility_in_NBZ();
    void Get_Bare_Susceptibility_in_NBZ_old();
    void Calculate_InteractionKernel();
    void Get_RPA_Susceptibility_Matrix();
    void Get_RPA_Susceptibility_in_NBZ();
    void Get_Spin_resolved_dc_conductivity();
    void Create_HK(int k1_, int k2_, Matrix<complex<double>> & Hk_);
    void Get_Spin_resolved_dc_conductivity_way2();
    void Get_Spin_resolved_dc_conductivity_way3();
    double FermiDis(double E, double mu, double InvKbT);
    double d_FermiDis_dE(double E, double mu, double InvKbT);
    void Create_Kspace_Spectrum_spin_resolved();

    //::DONE



    mt19937_64 &Generator1_;
    uniform_real_distribution<double> dis1_;
    Parameters_G2dLatticeNew &Parameters_;
    Coordinates_G2dLatticeNew &Coordinates_; //this in cell wise representation, n_orbs=no. of atoms in unitcell
    Connections_G2dLatticeNew &Connections_;
    int lx_, ly_, ncells_, n_orbs_, n_atoms_, lx_cells, ly_cells;
    Mat_1_intpair IntraMUC_sites_array;
    Matrix<complex<double>> Ham_;
    vector<double> eigs_;
    Mat_2_Complex_doub Eigvectors_;
    Mat_1_doub Eigenvalues_;
    Mat_2_Complex_doub Eigvectors_saved;
    Mat_1_doub Eigenvalues_saved;
    Mat_1_doub Kx_values, Ky_values;
    //Mat_1_Complex_doub V_mat;

    Mat_1_doub Eigenvalues_spin_up, Eigenvalues_spin_down;

    Mat_2_Complex_doub Susc_OprA, Susc_OprB;

    Matrix<double> Tx, Ty, Tpxpy, Tpxmy;
    Matrix<double> IntraCell_Hopp, InterCell_px_Hopp, InterCell_py_Hopp, InterCell_pxmy_Hopp;


    Matrix_COO_doub_tensor V_int;
    Matrix<double> V_int_OP_check;

    Matrix< Matrix<complex<double>> > P_mat, W_mat;
    Matrix<complex<double>> M_mat_up, M_mat_dn;
    Matrix<complex<double>> V_bar, A_mat, B_mat, C_mat;



    Matrix_COO_Complex OPs_, OPs_new_;
    double OPs_total_den;
    Mat_1_int SI_to_ind;

    double mu_;
    double OP_error_;
    double E_quant, E_class;

    double E_class_onsite_U0_Hartree, E_class_longrange_Hartree;
    double E_class_onsite_U0_Fock, E_class_longrange_Fock;

    Mat_1_intpair OP_Dcell;
    Mat_1_intpair OP_alpha;
    Mat_1_intpair OP_sigma;

    Mat_2_intpair cell_no;
    Mat_2_int magnetic_atom_no;

    //Declarations for Anderson Mixing
    Mat_1_doub x_km1_, x_k_, Del_x_km1;
    Mat_1_doub f_k_, f_km1_, Del_f_km1;
    Mat_1_doub xbar_k_, fbar_k_, gamma_k_, x_kp1_;
    Matrix<double> X_mat, F_mat;


    vector< Matrix<complex<double>> > J_KE_e1, J_KE_e2;
    vector< Matrix<complex<double>> > J_KE_X, J_KE_Y;

    double kick_while_cooling;

    double Global_Eps;
    bool OP_only_finite_Int;

    int NOT_AVAIL_INT;

    Matrix<int> Mat_MUC;
    int NSites_in_MUC;
    Mat_1_intpair Intra_MUC_positions;

    Mat_4_Complex_doub ChiBare_AB_q_omega;
    Mat_5_Complex_doub IntKer;

    Mat_6_Complex_doub ChiBareMat;
    Mat_6_Complex_doub ChiBareTimesI;
    vector < vector < Matrix <complex<double>> > > OneMinusChiBareTimesI_Inv;
    Mat_6_Complex_doub ChiRPAMat;

    Mat_4_Complex_doub AkwMat;


    //Mat_2_Complex_doub Chi0_;


};



double Kspace_calculation_G2dLatticeNew::FermiDis(double E, double mu, double InvKbT){

    double val;
    val  = 1.0/(exp((E-mu)*InvKbT) + 1.0);
    return val;
}

double Kspace_calculation_G2dLatticeNew::d_FermiDis_dE(double E, double mu, double InvKbT){


    double val;

    if( ((E-mu)*InvKbT)>30.0 ){
        val=0.0;
    }
    else{
    val = -InvKbT*exp((E-mu)*InvKbT)*
            (1.0/(exp((E-mu)*InvKbT) + 1.0))*
            (1.0/(exp((E-mu)*InvKbT) + 1.0));
    }
    return val;
}

void Kspace_calculation_G2dLatticeNew::Get_Spin_resolved_dc_conductivity(){

    int S_=NSites_in_MUC;
    int NBands_ =2*n_orbs_*n_atoms_*S_;
    Mat_1_doub theta_array;
    int NThetaSlices=100;
    double d_theta=(2.0*PI)/(NThetaSlices);

    theta_array.resize(NThetaSlices);
    for(int n=0;n<NThetaSlices;n++){
        theta_array[n] = d_theta*n;
    }


    double Gamma_brdn=Parameters_.eta;

    Mat_5_Complex_doub C_mat_vel_p1; //[k1][k2][m][n][s]
    Mat_5_Complex_doub C_mat_vel_p2; //[k1][k2][m][n][s]
    Mat_5_Complex_doub C_mat_vel_p3; //[k1][k2][m][n][s]



    C_mat_vel_p1.resize(lx_cells);
    C_mat_vel_p2.resize(lx_cells);
    C_mat_vel_p3.resize(lx_cells);
    for(int n1=0;n1<lx_cells;n1++){
        C_mat_vel_p1[n1].resize(ly_cells);
        C_mat_vel_p2[n1].resize(ly_cells);
        C_mat_vel_p3[n1].resize(ly_cells);
        for(int n2=0;n2<ly_cells;n2++){
            C_mat_vel_p1[n1][n2].resize(NBands_);
            C_mat_vel_p2[n1][n2].resize(NBands_);
            C_mat_vel_p3[n1][n2].resize(NBands_);
            for(int m=0;m<NBands_;m++){
                C_mat_vel_p1[n1][n2][m].resize(NBands_);
                C_mat_vel_p2[n1][n2][m].resize(NBands_);
                C_mat_vel_p3[n1][n2][m].resize(NBands_);
                for(int n=0;n<NBands_;n++){
                    C_mat_vel_p1[n1][n2][m][n].resize(2);
                    C_mat_vel_p2[n1][n2][m][n].resize(2);
                    C_mat_vel_p3[n1][n2][m][n].resize(2);
                    for(int s=0;s<2;s++){
                        C_mat_vel_p1[n1][n2][m][n][s]=0.0;
                        C_mat_vel_p2[n1][n2][m][n][s]=0.0;
                        C_mat_vel_p3[n1][n2][m][n][s]=0.0;
                    }
                }
            }
        }
    }

    Mat_2_Complex_doub Cond_dc;//[s][theta]
    Cond_dc.resize(2);
    for(int spin=0;spin<2;spin++){
        Cond_dc[spin].resize(NThetaSlices);
    }


    Matrix<complex<double>> Hk_plus_dkx, Hk_plus_dky, Hk_;
    complex<double> part1_, part2_, part3_;
    int k_index, state_n_k, state_m_k;


    for(int s=0;s<2;s++){
        for(int m=0;m<NBands_;m++){
            for(int n=0;n<NBands_;n++){

                for(int k1_=0;k1_<lx_cells;k1_++){
                    for(int k2_=0;k2_<ly_cells;k2_++){

                        cout<<s<<"  "<<m<<"  "<<n<<"  "<<k1_<<"  "<<k2_<<endl;
                        k_index = Coordinates_.Ncell(k1_,k2_);
                        state_n_k = NBands_*k_index + n;
                        state_m_k = NBands_*k_index + m;

                        Create_HK(k1_, k2_, Hk_);
                        Create_HK( (k1_+1)%lx_cells, k2_, Hk_plus_dkx);
                        Create_HK( k1_, (k2_+1)%ly_cells, Hk_plus_dky);


                        C_mat_vel_p1[k1_][k2_][m][n][s]=0.0;
                        C_mat_vel_p2[k1_][k2_][m][n][s]=0.0;
                        C_mat_vel_p3[k1_][k2_][m][n][s]=0.0;
                        for(int ci_=0;ci_<n_orbs_*n_atoms_*S_;ci_++){ //without spin
                            for(int cj_=0;cj_<n_orbs_*n_atoms_*S_;cj_++){
                                int ci= ci_ + s*(n_orbs_*n_atoms_*S_);
                                int cj= cj_ + s*(n_orbs_*n_atoms_*S_);

                                C_mat_vel_p1[k1_][k2_][m][n][s] += conj(Eigvectors_[state_m_k][ci])*Hk_plus_dkx(ci,cj)*Eigvectors_[state_n_k][cj];
                                C_mat_vel_p2[k1_][k2_][m][n][s] += conj(Eigvectors_[state_m_k][ci])*Hk_plus_dky(ci,cj)*Eigvectors_[state_n_k][cj];
                                C_mat_vel_p3[k1_][k2_][m][n][s] += conj(Eigvectors_[state_m_k][ci])*Hk_(ci,cj)*Eigvectors_[state_n_k][cj];

                            }
                        }

                    }
                }

            }
        }
    }


    //sigma[s][theta_slice]
    complex<double> C_elm;
    for(int spin=0;spin<2;spin++){
        for(int theta_slice=0;theta_slice<NThetaSlices;theta_slice++){
            Cond_dc[spin][theta_slice]=0.0;

            for(int m=0;m<NBands_;m++){
                for(int n=0;n<NBands_;n++){

                    for(int k1_=0;k1_<lx_cells;k1_++){
                        for(int k2_=0;k2_<ly_cells;k2_++){
                            k_index = Coordinates_.Ncell(k1_,k2_);
                            state_n_k = NBands_*k_index + n;
                            state_m_k = NBands_*k_index + m;

                            C_elm = C_mat_vel_p1[k1_][k2_][m][n][spin]*cos(theta_array[theta_slice])
                                    + C_mat_vel_p2[k1_][k2_][m][n][spin]*sin(theta_array[theta_slice])
                                    - C_mat_vel_p3[k1_][k2_][m][n][spin]*(cos(theta_array[theta_slice]) + sin(theta_array[theta_slice]) );


                            Cond_dc[spin][theta_slice] += (PI/(1.0*lx_*ly_))*
                                                         conj(C_elm)*C_elm*
                                    Lorentzian(mu_-Eigenvalues_saved[state_n_k], 2.0*Gamma_brdn)*
                                    Lorentzian(mu_-Eigenvalues_saved[state_m_k], 2.0*Gamma_brdn);
                        }}
                }}
        }
    }




    string file_out="dc_cond_spin_resolved.txt";
    ofstream file_out_strm(file_out.c_str());
    file_out_strm<<"#theta  sigma_up_up.real()  sigam_dn_dn.real()    imag    imag | imag should be zero"<<endl;
    for(int theta_slice=0;theta_slice<NThetaSlices;theta_slice++){
        file_out_strm<<theta_array[theta_slice]<<"   "<<Cond_dc[0][theta_slice].real()<<"   "<<Cond_dc[1][theta_slice].real()<<"    "<<Cond_dc[0][theta_slice].imag()<<"   "<<Cond_dc[1][theta_slice].imag()<<endl;
    }





}





void Kspace_calculation_G2dLatticeNew::Get_Spin_resolved_dc_conductivity_way2(){

    int S_=NSites_in_MUC;
    int NBands_ =2*n_orbs_*n_atoms_*S_;
    Mat_1_doub theta_array;
    int NThetaSlices=100;
    double d_theta=(2.0*PI)/(NThetaSlices);

    theta_array.resize(NThetaSlices);
    for(int n=0;n<NThetaSlices;n++){
        theta_array[n] = d_theta*n;
    }


    double Gamma_brdn=Parameters_.eta;

    Mat_1_doub cond_xx, cond_yy, cond_xy, cond_yx;
    cond_xx.resize(2);cond_yy.resize(2);cond_xy.resize(2);cond_yx.resize(2);
    double spin_contr;

    for(int spin=0;spin<2;spin++){
    cond_xx[spin]=0.0;cond_yy[spin]=0.0;
    cond_xy[spin]=0.0;cond_yx[spin]=0.0;

    for(int n=0;n<NBands_;n++){
        for(int k1_=0;k1_<lx_cells;k1_++){
            for(int k2_=0;k2_<ly_cells;k2_++){

                //(k1_+1)%lx_cells, k2_
                //k_index = Coordinates_.Ncell(k1_,k2_);
                //state_n_k = NBands_*k_index + n;
               //Eigenvalues_saved[state_n_k]

                int k1_p_x = (k1_+1)%lx_cells;
                int k2_p_y = (k2_+1)%ly_cells;

                int k1_m_x = (k1_-1+lx_cells)%lx_cells;
                int k2_m_y = (k2_-1+ly_cells)%ly_cells;

                int k_index = Coordinates_.Ncell(k1_,k2_);
                int state_n_k = NBands_*k_index + n;

                spin_contr=0.0;
                for(int alpha=0;alpha<n_orbs_*n_atoms_*S_;alpha++){
                    int ci= alpha + spin*(n_orbs_*n_atoms_*S_);
                    spin_contr += abs(Eigvectors_[state_n_k][ci])*abs(Eigvectors_[state_n_k][ci]);
                }

                int k_index_p_x = Coordinates_.Ncell(k1_p_x,k2_);
                int k_index_p_y = Coordinates_.Ncell(k1_,k2_p_y);

                int k_index_m_x = Coordinates_.Ncell(k1_m_x,k2_);
                int k_index_m_y = Coordinates_.Ncell(k1_,k2_m_y);


                int state_n_k_p_x = NBands_*k_index_p_x + n;
                int state_n_k_p_y = NBands_*k_index_p_y + n;

                int state_n_k_m_x = NBands_*k_index_m_x + n;
                int state_n_k_m_y = NBands_*k_index_m_y + n;

                cond_xx[spin] += spin_contr*0.25*
                        (Eigenvalues_saved[state_n_k_p_x]-Eigenvalues_saved[state_n_k_m_x])*
                        (Eigenvalues_saved[state_n_k_p_x]-Eigenvalues_saved[state_n_k_m_x])*
                        d_FermiDis_dE(Eigenvalues_saved[state_n_k], mu_, Parameters_.beta);
                       // Lorentzian(mu_-Eigenvalues_saved[state_n_k], 2.0*Gamma_brdn);

                cond_yy[spin] += spin_contr*0.25*
                        (Eigenvalues_saved[state_n_k_p_y]-Eigenvalues_saved[state_n_k_m_y])*
                        (Eigenvalues_saved[state_n_k_p_y]-Eigenvalues_saved[state_n_k_m_y])*
                         d_FermiDis_dE(Eigenvalues_saved[state_n_k], mu_, Parameters_.beta);

                cond_xy[spin] += spin_contr*0.25*
                        (Eigenvalues_saved[state_n_k_p_x]-Eigenvalues_saved[state_n_k_m_x])*
                        (Eigenvalues_saved[state_n_k_p_y]-Eigenvalues_saved[state_n_k_m_y])*
                         d_FermiDis_dE(Eigenvalues_saved[state_n_k], mu_, Parameters_.beta);


    }
    }
    }
    }

    cout<<"xx : "<<cond_xx[0]<<"  "<<cond_xx[1]<<endl;
    cout<<"yy : "<<cond_yy[0]<<"  "<<cond_yy[1]<<endl;
    cout<<"xy : "<<cond_xy[0]<<"  "<<cond_xy[1]<<endl;

    string file_out="dc_cond_spin_resolved_way2.txt";
    ofstream file_out_strm(file_out.c_str());
    file_out_strm<<"#theta  sigma_up_up.real()  sigam_dn_dn.real()"<<endl;
    for(int theta_slice=0;theta_slice<NThetaSlices;theta_slice++){
        file_out_strm<<theta_array[theta_slice];
        for(int spin=0;spin<2;spin++){
        file_out_strm<<"   "<<
            (cond_xx[spin]*cos(theta_array[theta_slice])*cos(theta_array[theta_slice]))
         +  (cond_yy[spin]*sin(theta_array[theta_slice])*sin(theta_array[theta_slice]))
         +  1.0*(cond_xy[spin]*2.0*sin(theta_array[theta_slice])*cos(theta_array[theta_slice]));
        }
        file_out_strm<<endl;
    }





}




void Kspace_calculation_G2dLatticeNew::Get_Spin_resolved_dc_conductivity_way3(){



    cout<<"dc cond [way 3]"<<endl;

    Create_Kspace_Spectrum_spin_resolved();


    int S_=NSites_in_MUC;
    int NBands_ = n_orbs_*n_atoms_*S_; //spin is separate here
    Mat_1_doub theta_array;
    int NThetaSlices=100;
    double d_theta=(2.0*PI)/(NThetaSlices);

    theta_array.resize(NThetaSlices);
    for(int n=0;n<NThetaSlices;n++){
        theta_array[n] = d_theta*n;
    }


    double Gamma_brdn=Parameters_.eta;

    Mat_1_doub cond_xx, cond_yy, cond_xy, cond_yx;
    cond_xx.resize(2);cond_yy.resize(2);cond_xy.resize(2);cond_yx.resize(2);
    double spin_contr;

    for(int spin=0;spin<2;spin++){
    cond_xx[spin]=0.0;cond_yy[spin]=0.0;
    cond_xy[spin]=0.0;cond_yx[spin]=0.0;

    for(int n=0;n<NBands_;n++){
        for(int k1_=0;k1_<lx_cells;k1_++){
            for(int k2_=0;k2_<ly_cells;k2_++){

                //(k1_+1)%lx_cells, k2_
                //k_index = Coordinates_.Ncell(k1_,k2_);
                //state_n_k = NBands_*k_index + n;
               //Eigenvalues_saved[state_n_k]

                int k1_p_x = (k1_+1)%lx_cells;
                int k2_p_y = (k2_+1)%ly_cells;

                int k1_p_2x = (k1_+2)%lx_cells;
                int k2_p_2y = (k2_+2)%ly_cells;

                int k1_m_x = (k1_-1+lx_cells)%lx_cells;
                int k2_m_y = (k2_-1+ly_cells)%ly_cells;

                int k1_m_2x = (k1_-2+lx_cells)%lx_cells;
                int k2_m_2y = (k2_-2+ly_cells)%ly_cells;

                int k_index = Coordinates_.Ncell(k1_,k2_);
                int state_n_k = NBands_*k_index + n;


                int k_index_p_x = Coordinates_.Ncell(k1_p_x,k2_);
                int k_index_p_y = Coordinates_.Ncell(k1_,k2_p_y);

                int k_index_p_2x = Coordinates_.Ncell(k1_p_2x,k2_);
                int k_index_p_2y = Coordinates_.Ncell(k1_,k2_p_2y);

                int k_index_m_x = Coordinates_.Ncell(k1_m_x,k2_);
                int k_index_m_y = Coordinates_.Ncell(k1_,k2_m_y);

                int k_index_m_2x = Coordinates_.Ncell(k1_m_2x,k2_);
                int k_index_m_2y = Coordinates_.Ncell(k1_,k2_m_2y);

                int state_n_k_p_x = NBands_*k_index_p_x + n;
                int state_n_k_p_y = NBands_*k_index_p_y + n;

                int state_n_k_p_2x = NBands_*k_index_p_2x + n;
                int state_n_k_p_2y = NBands_*k_index_p_2y + n;

                int state_n_k_m_x = NBands_*k_index_m_x + n;
                int state_n_k_m_y = NBands_*k_index_m_y + n;

                int state_n_k_m_2x = NBands_*k_index_m_2x + n;
                int state_n_k_m_2y = NBands_*k_index_m_2y + n;


                //Using Five-point Midpoint Formula
                if(spin==0){
                cond_xx[spin] += (1.0/144.0)*
                        (Eigenvalues_spin_up[state_n_k_m_2x] - 8.0*Eigenvalues_spin_up[state_n_k_m_x] + 8.0*Eigenvalues_spin_up[state_n_k_p_x]-Eigenvalues_spin_up[state_n_k_p_2x])*
                        (Eigenvalues_spin_up[state_n_k_m_2x] - 8.0*Eigenvalues_spin_up[state_n_k_m_x] + 8.0*Eigenvalues_spin_up[state_n_k_p_x]-Eigenvalues_spin_up[state_n_k_p_2x])*
                        d_FermiDis_dE(Eigenvalues_spin_up[state_n_k], mu_, Parameters_.beta);
                       // Lorentzian(mu_-Eigenvalues_saved[state_n_k], 2.0*Gamma_brdn);

                cond_yy[spin] +=  (1.0/144.0)*
                        (Eigenvalues_spin_up[state_n_k_m_2y] - 8.0*Eigenvalues_spin_up[state_n_k_m_y] + 8.0*Eigenvalues_spin_up[state_n_k_p_y]-Eigenvalues_spin_up[state_n_k_p_2y])*
                        (Eigenvalues_spin_up[state_n_k_m_2y] - 8.0*Eigenvalues_spin_up[state_n_k_m_y] + 8.0*Eigenvalues_spin_up[state_n_k_p_y]-Eigenvalues_spin_up[state_n_k_p_2y])*
                         d_FermiDis_dE(Eigenvalues_spin_up[state_n_k], mu_, Parameters_.beta);

                cond_xy[spin] +=  (1.0/144.0)*
                        (Eigenvalues_spin_up[state_n_k_m_2x] - 8.0*Eigenvalues_spin_up[state_n_k_m_x] + 8.0*Eigenvalues_spin_up[state_n_k_p_x]-Eigenvalues_spin_up[state_n_k_p_2x])*
                        (Eigenvalues_spin_up[state_n_k_m_2y] - 8.0*Eigenvalues_spin_up[state_n_k_m_y] + 8.0*Eigenvalues_spin_up[state_n_k_p_y]-Eigenvalues_spin_up[state_n_k_p_2y])*
                         d_FermiDis_dE(Eigenvalues_spin_up[state_n_k], mu_, Parameters_.beta);
                }
                if(spin==1){
                    cond_xx[spin] += (1.0/144.0)*
                            (Eigenvalues_spin_down[state_n_k_m_2x] - 8.0*Eigenvalues_spin_down[state_n_k_m_x] + 8.0*Eigenvalues_spin_down[state_n_k_p_x]-Eigenvalues_spin_down[state_n_k_p_2x])*
                            (Eigenvalues_spin_down[state_n_k_m_2x] - 8.0*Eigenvalues_spin_down[state_n_k_m_x] + 8.0*Eigenvalues_spin_down[state_n_k_p_x]-Eigenvalues_spin_down[state_n_k_p_2x])*
                            d_FermiDis_dE(Eigenvalues_spin_down[state_n_k], mu_, Parameters_.beta);
                           // Lorentzian(mu_-Eigenvalues_saved[state_n_k], 2.0*Gamma_brdn);

                    cond_yy[spin] +=  (1.0/144.0)*
                            (Eigenvalues_spin_down[state_n_k_m_2y] - 8.0*Eigenvalues_spin_down[state_n_k_m_y] + 8.0*Eigenvalues_spin_down[state_n_k_p_y]-Eigenvalues_spin_down[state_n_k_p_2y])*
                            (Eigenvalues_spin_down[state_n_k_m_2y] - 8.0*Eigenvalues_spin_down[state_n_k_m_y] + 8.0*Eigenvalues_spin_down[state_n_k_p_y]-Eigenvalues_spin_down[state_n_k_p_2y])*
                             d_FermiDis_dE(Eigenvalues_spin_down[state_n_k], mu_, Parameters_.beta);

                    cond_xy[spin] +=  (1.0/144.0)*
                            (Eigenvalues_spin_down[state_n_k_m_2x] - 8.0*Eigenvalues_spin_down[state_n_k_m_x] + 8.0*Eigenvalues_spin_down[state_n_k_p_x]-Eigenvalues_spin_down[state_n_k_p_2x])*
                            (Eigenvalues_spin_down[state_n_k_m_2y] - 8.0*Eigenvalues_spin_down[state_n_k_m_y] + 8.0*Eigenvalues_spin_down[state_n_k_p_y]-Eigenvalues_spin_down[state_n_k_p_2y])*
                             d_FermiDis_dE(Eigenvalues_spin_down[state_n_k], mu_, Parameters_.beta);
                }


    }
    }
    }
    }

    cout<<"xx : "<<cond_xx[0]<<"  "<<cond_xx[1]<<endl;
    cout<<"yy : "<<cond_yy[0]<<"  "<<cond_yy[1]<<endl;
    cout<<"xy : "<<cond_xy[0]<<"  "<<cond_xy[1]<<endl;

    string file_out="dc_cond_spin_resolved_way3.txt";
    ofstream file_out_strm(file_out.c_str());
    file_out_strm<<"#theta  sigma_up_up.real()  sigam_dn_dn.real()"<<endl;
    for(int theta_slice=0;theta_slice<NThetaSlices;theta_slice++){
        file_out_strm<<theta_array[theta_slice];
        for(int spin=0;spin<2;spin++){
        file_out_strm<<"   "<<
            (cond_xx[spin]*cos(theta_array[theta_slice])*cos(theta_array[theta_slice]))
         +  (cond_yy[spin]*sin(theta_array[theta_slice])*sin(theta_array[theta_slice]))
         +  1.0*(cond_xy[spin]*2.0*sin(theta_array[theta_slice])*cos(theta_array[theta_slice]));
        }
        file_out_strm<<endl;
    }





}







void Kspace_calculation_G2dLatticeNew::Create_Kspace_Spectrum_spin_resolved(){


    assert(Parameters_.NoSpinFlipOP);
    assert(Parameters_.Just_Hartree);
    bool NA_spinflip=Parameters_.NoSpinFlipOP;

    // bool PinningField=true;
    // double PinningFieldValue=Parameters_.PinningFieldValue;


    int S_= NSites_in_MUC;
    //Calculating density using OP's will be going inside Hamiltonian
    OPs_total_den=0.0;
    int index_OP;
    int row_OP, col_OP;


    Matrix<complex<double>> Hamup_, Hamdn_;
    Hamup_.resize(n_orbs_*n_atoms_*S_,n_orbs_*n_atoms_*S_);
    Hamdn_.resize(n_orbs_*n_atoms_*S_,n_orbs_*n_atoms_*S_);


    int row_, col_;
    int UP_, DOWN_;
    UP_=0;DOWN_=1;
    int spin_bar, gamma_bar;
    int k_index;
    double k1x_val, k2x_val, k1y_val, k2y_val;
    double fac_;
    complex<double> OP_value;

    //cout <<"Here 1"<<endl;



    // Create_M_mat();
    // Create_P_mat();

    for(int k1=0;k1<lx_cells;k1++){
        k1x_val = (2.0*PI*k1)/(1.0*lx_cells);
        k1y_val = ((2.0*PI*k1)/(1.0*lx_cells))*(-1.0/sqrt(3));

        for(int k2=0;k2<ly_cells;k2++){
            k2x_val = 0.0;
            k2y_val = ((2.0*PI*k2)/(1.0*ly_cells))*(2.0/sqrt(3));

            k_index = Coordinates_.Ncell(k1,k2);
            Kx_values[k_index]=k1x_val + k2x_val;
            Ky_values[k_index]=k1y_val + k2y_val;

            //cout<<"k1, k2 :" <<k1<<" "<<k2<<endl;

            Hamup_.fill(0.0);
            Hamdn_.fill(0.0);

            //Hoppings
            for(int alpha_p=0;alpha_p<S_;alpha_p++){
                for(int gamma_p=0;gamma_p<n_orbs_*n_atoms_;gamma_p++){
                    //for(int sigma_p=0;sigma_p<2;sigma_p++){
                        col_ = alpha_p + gamma_p*(S_);

                        for(int alpha=0;alpha<S_;alpha++){
                            for(int gamma=0;gamma<n_orbs_*n_atoms_;gamma++){
                               // for(int sigma=0;sigma<2;sigma++){
                               // int sigma=sigma_p;
                                    row_ = alpha + gamma*(S_);

                                    //Ham_(row_, col_) += h_KE(alpha, gamma, sigma, alpha_p, gamma_p, sigma_p, k1, k2);
                                    Hamup_(row_,col_) += h_KE(alpha, gamma, UP_, alpha_p, gamma_p, UP_, k1, k2);
                                    Hamdn_(row_,col_) += h_KE(alpha, gamma, DOWN_, alpha_p, gamma_p, DOWN_, k1, k2);
                                //}
                            }
                        }
                    //}
                }
            }





            /*
            //Anisotropy
            for(int alpha=0;alpha<S_;alpha++){
                for(int spin1=0;spin1<2;spin1++){
                    for(int spin2=0;spin2<2;spin2++){

                        row_ = alpha + spin1*(S_);
                        col_=row_;

                        row_OP = 0*(2*S_) + alpha + spin2*(S_);
                        col_OP = 0*(2*S_) + alpha + spin2*(S_);
                        index_OP = SI_to_ind[col_OP + row_OP*(2*ncells_*S_)];

                        fac_ = 1.0 - (2.0*abs(spin1-spin2)); //i.e 0 --> 1, 1 --> -1

                        Ham_(row_,col_) += (-1.0/2.0)*Parameters_.AnisotropyZ*fac_*OPs_.value[index_OP];
                    }
                }
            }
            */


            //Onsite Energies
            //OnSiteE_up
            for(int alpha=0;alpha<S_;alpha++){
                for(int gamma=0;gamma<n_atoms_*n_orbs_;gamma++){
                    //for(int spin1=0;spin1<2;spin1++){
                        row_ = alpha + gamma*(S_);
                        col_=row_;
                        Hamup_(row_,col_) += Parameters_.OnSiteE[gamma][UP_];
                        Hamdn_(row_,col_) += Parameters_.OnSiteE[gamma][DOWN_];
                    //}
                }
            }


            //magnetic field along z
            for(int alpha=0;alpha<S_;alpha++){
                for(int gamma=0;gamma<n_atoms_*n_orbs_;gamma++){
                   // for(int spin1=0;spin1<2;spin1++){
                        row_ = alpha + gamma*(S_);
                        col_=row_;
                        Hamup_(row_,col_) += Parameters_.PinningFieldZ[gamma]*(0.5 - 1.0*UP_);
                        Hamdn_(row_,col_) += Parameters_.PinningFieldZ[gamma]*(0.5 - 1.0*DOWN_);
                    //}
                }}




            complex<double> OP_val;
            //Interaction local intra-orbital U0:
            for(int j=0;j<S_;j++){
                for(int b_beta=0;b_beta<n_atoms_;b_beta++){
                    for(int gamma=0;gamma<n_orbs_;gamma++){
                        int atom_plus_orb=b_beta +  n_atoms_*gamma;
                        //for(int spin_=0;spin_<2;spin_++){
                            col_ = j + atom_plus_orb*(S_);

                           // for(int spin_p=0;spin_p<2;spin_p++){

                                row_ = j + atom_plus_orb*(S_);

                                //if(spin_==spin_p){
                                    row_OP = j + atom_plus_orb*(S_)
                                            + (1-UP_)*(n_atoms_*n_orbs_*S_);
                                    col_OP = row_OP;
                                    index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                    Hamup_(row_,col_) += Parameters_.U0[b_beta]*OPs_.value[index_OP];

                                    row_OP = j + atom_plus_orb*(S_)
                                            + (1-DOWN_)*(n_atoms_*n_orbs_*S_);
                                    col_OP = row_OP;
                                    index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                    Hamdn_(row_,col_) += Parameters_.U0[b_beta]*OPs_.value[index_OP];

                                //}

                            //}
                       // }
                    }
                }
            }


            //     if(k1==0 && k2==0){
            //      cout<<"---------- Printing matrix kx=0=ky -------------"<<endl;
            //      Ham_.print();
            //      cout<<"-----------------"<<endl;
            //  }

            //            cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
            //            Ham_.print();
            //            cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;




            char Dflag='V';
            Mat_1_doub Eigval_up, Eigval_dn;
            DiagonalizeGivenMat(Dflag, Hamup_, Eigval_up);
            DiagonalizeGivenMat(Dflag, Hamdn_, Eigval_dn);

            for(int row=0;row<n_orbs_*n_atoms_*S_;row++){
                Eigenvalues_spin_up[n_orbs_*n_atoms_*S_*k_index + row]=Eigval_up[row];
                Eigenvalues_spin_down[n_orbs_*n_atoms_*S_*k_index + row]=Eigval_dn[row];
            }


        }
    }




//----------------------------------------------------------
    //Bands
    char temp_char[50];
    sprintf(temp_char, "%.10f", Parameters_.Temperature);

    string File_Out_Bands;
    File_Out_Bands = "BandsSpinResolved_" + string(temp_char)+ ".txt";
    ofstream file_out_bands(File_Out_Bands.c_str());


    int kx_i, ky_i;

    Mat_1_intpair k_path;
    k_path.clear();
    pair_int temp_pair;


    // ---k_path---------

    //--------\Gamma to M-----------------
    ky_i = 0;
    for (kx_i = 0; kx_i <= (lx_cells/ 2); kx_i++)
    {
        temp_pair.first = kx_i;
        temp_pair.second = kx_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------


    //--------M to X-----------------
    kx_i = (lx_cells / 2);
    for (ky_i = (ly_cells / 2)-1; ky_i >=0; ky_i--)
    {
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------


    //--------X to Gamma-----------------
    ky_i = 0;
    for (kx_i = (lx_cells / 2)-1; kx_i >=0; kx_i--)
    {
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------


    //--------Gamma to Y-----------------
    kx_i = 0;
    ky_i = 0;
    for (ky_i = 1; ky_i <=(ly_cells/2); ky_i++)
    {
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------

    //--------Y to M-----------------
    kx_i = 0;
    ky_i = ly_cells/2;
    for (kx_i = 1; kx_i <=(lx_cells/2); kx_i++)
    {
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------


    //--------M to \Gamma[with one extra point,
    //                  because in gnuplot use "set pm3d corners2color c1"
    //                  ]-----------------
    kx_i = (lx_cells / 2) - 1;
    ky_i = (ly_cells / 2) - 1;
    for (kx_i = (lx_cells / 2) - 1; kx_i >= 0; kx_i--)
    {
        ky_i=kx_i;
        if(ky_i < ly_cells){
            temp_pair.first = kx_i;
            temp_pair.second = kx_i;
            k_path.push_back(temp_pair);
        }
    }

    temp_pair.first = 0;
    temp_pair.second = 0;
    k_path.push_back(temp_pair);

    file_out_bands<<"#kindex  E_up(k,1)  E_up(k,2) ... DOWN  ....."<<endl;
    for (int k_point = 0; k_point < k_path.size(); k_point++)
    {
        k_index=Coordinates_.Ncell(k_path[k_point].first, k_path[k_point].second);

        //Eigvectors_[state][c1]
        //k_index = Coordinates_.Ncell(k1,k2);
        //state = 2*n_atoms_*n_orbs_*S_*k_index + n;
        //c1 = alpha_row + gamma_row*(S_) +  sigma_row*(n_atoms_*n_orbs_*S_);
        file_out_bands<<k_point<<"   ";
        for(int band=0;band<n_orbs_*n_atoms_*S_;band++){
            int state = n_atoms_*n_orbs_*S_*k_index + band;
            file_out_bands<<Eigenvalues_spin_up[state]<<"   ";
        }
        for(int band=0;band<n_orbs_*n_atoms_*S_;band++){
            int state = n_atoms_*n_orbs_*S_*k_index + band;
            file_out_bands<<Eigenvalues_spin_down[state]<<"   ";
        }
        file_out_bands<<endl;
    }







}






void Kspace_calculation_G2dLatticeNew::Create_HK(int k1_, int k2_, Matrix<complex<double>> & Hk_){


    int S_= NSites_in_MUC;
    int NBands_ =2*n_orbs_*n_atoms_*S_;

    Hk_.clear();
    Hk_.resize(NBands_,NBands_);

    //Calculating density using OP's will be going inside Hamiltonian
    int index_OP;
    int row_OP, col_OP;


    int row_, col_;
    int UP_, DOWN_;
    UP_=0;DOWN_=1;
    int spin_bar, gamma_bar;
    int k_index;
    double k1x_val, k2x_val, k1y_val, k2y_val;
    double fac_;
    complex<double> OP_value;

    //cout <<"Here 1"<<endl;



    // Create_M_mat();
    // Create_P_mat();

    int k1=k1_;
    int k2=k2_;
    k1x_val = (2.0*PI*k1)/(1.0*lx_cells);
    k1y_val = ((2.0*PI*k1)/(1.0*lx_cells))*(-1.0/sqrt(3));

    k2x_val = 0.0;
    k2y_val = ((2.0*PI*k2)/(1.0*ly_cells))*(2.0/sqrt(3));

    k_index = Coordinates_.Ncell(k1,k2);
    Kx_values[k_index]=k1x_val + k2x_val;
    Ky_values[k_index]=k1y_val + k2y_val;

    Hk_.fill(0.0);

    //cout<<"k1, k2 :" <<k1<<" "<<k2<<endl;

    //Hoppings
    for(int alpha_p=0;alpha_p<S_;alpha_p++){
        for(int gamma_p=0;gamma_p<n_orbs_*n_atoms_;gamma_p++){
            for(int sigma_p=0;sigma_p<2;sigma_p++){
                col_ = alpha_p + gamma_p*(S_)
                        + sigma_p*(n_orbs_*n_atoms_*S_);

                for(int alpha=0;alpha<S_;alpha++){
                    for(int gamma=0;gamma<n_orbs_*n_atoms_;gamma++){
                        for(int sigma=0;sigma<2;sigma++){
                            row_ = alpha + gamma*(S_)
                                    + sigma*(n_orbs_*n_atoms_*S_);

                            Hk_(row_, col_) += h_KE(alpha, gamma, sigma, alpha_p, gamma_p, sigma_p, k1, k2);
                        }
                    }
                }
            }
        }
    }





    /*
            //Anisotropy
            for(int alpha=0;alpha<S_;alpha++){
                for(int spin1=0;spin1<2;spin1++){
                    for(int spin2=0;spin2<2;spin2++){

                        row_ = alpha + spin1*(S_);
                        col_=row_;

                        row_OP = 0*(2*S_) + alpha + spin2*(S_);
                        col_OP = 0*(2*S_) + alpha + spin2*(S_);
                        index_OP = SI_to_ind[col_OP + row_OP*(2*ncells_*S_)];

                        fac_ = 1.0 - (2.0*abs(spin1-spin2)); //i.e 0 --> 1, 1 --> -1

                        Hk_(row_,col_) += (-1.0/2.0)*Parameters_.AnisotropyZ*fac_*OPs_.value[index_OP];
                    }
                }
            }
            */


    //Onsite Energies
    //OnSiteE_up
    for(int alpha=0;alpha<S_;alpha++){
        for(int gamma=0;gamma<n_atoms_*n_orbs_;gamma++){
            for(int spin1=0;spin1<2;spin1++){
                row_ = alpha + gamma*(S_)
                        + spin1*(n_atoms_*n_orbs_*S_);
                col_=row_;
                Hk_(row_,col_) += Parameters_.OnSiteE[gamma][spin1];
            }
        }
    }


    //magnetic field along z
    for(int alpha=0;alpha<S_;alpha++){
        for(int gamma=0;gamma<n_atoms_*n_orbs_;gamma++){
            for(int spin1=0;spin1<2;spin1++){
                row_ = alpha + gamma*(S_)
                        + spin1*(n_atoms_*n_orbs_*S_);
                col_=row_;
                Hk_(row_,col_) += Parameters_.PinningFieldZ[gamma]*(0.5 - 1.0*spin1);
            }
        }}

    //magnetic field along x
    for(int alpha=0;alpha<S_;alpha++){
        for(int gamma=0;gamma<n_atoms_*n_orbs_;gamma++){
            for(int spin1=0;spin1<2;spin1++){
                row_ = alpha + gamma*(S_)
                        + spin1*(n_atoms_*n_orbs_*S_);
                col_= alpha + gamma*(S_)
                        + (1-spin1)*(n_atoms_*n_orbs_*S_);
                Hk_(row_,col_) += Parameters_.PinningFieldX[gamma]*0.5;
            }
        }}


    bool NA_spinflip=Parameters_.NoSpinFlipOP;
    complex<double> OP_val;
    //Interaction local intra-orbital U0:
    for(int j=0;j<S_;j++){
        for(int b_beta=0;b_beta<n_atoms_;b_beta++){
            for(int gamma=0;gamma<n_orbs_;gamma++){
                int atom_plus_orb=b_beta +  n_atoms_*gamma;
                for(int spin_=0;spin_<2;spin_++){
                    col_ = j + atom_plus_orb*(S_) +
                            spin_*(n_atoms_*n_orbs_*S_);

                    for(int spin_p=0;spin_p<2;spin_p++){

                        row_ = j + atom_plus_orb*(S_) +
                                spin_p*(n_atoms_*n_orbs_*S_);

                        if(spin_==spin_p){
                            row_OP = j + atom_plus_orb*(S_)
                                    + (1-spin_)*(n_atoms_*n_orbs_*S_);
                            col_OP = row_OP;
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];

                            Hk_(row_,col_) += Parameters_.U0[b_beta]*OPs_.value[index_OP];
                        }
                        else{
                            if(!Parameters_.Just_Hartree){
                                if(!NA_spinflip){
                                    row_OP = j + atom_plus_orb*(S_)
                                            + (spin_)*(n_atoms_*n_orbs_*S_);
                                    col_OP = j + atom_plus_orb*(S_)
                                            + (spin_p)*(n_atoms_*n_orbs_*S_);

                                    if(col_OP>row_OP){
                                        index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                        OP_val = OPs_.value[index_OP];
                                    }
                                    else{
                                        index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                        OP_val = conj(OPs_.value[index_OP]);
                                    }

                                    Hk_(row_,col_) -=Parameters_.U0[b_beta]*OP_val;
                                }
                            }
                        }
                    }
                }
            }
        }
    }


    //Inter-orbital nn
    //Hartree
    for(int j=0;j<S_;j++){
        for(int atom_no=0;atom_no<n_atoms_;atom_no++){
            for(int alpha=0;alpha<n_orbs_;alpha++){
                for(int spin_=0;spin_<2;spin_++){
                    for(int beta=(alpha+1);beta<n_orbs_;beta++){
                        for(int spin_p=0;spin_p<2;spin_p++){
                            col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                    spin_*(n_atoms_*n_orbs_*S_);
                            row_=col_;
                            row_OP = j + (atom_no + n_atoms_*beta)*(S_)
                                    + (spin_p)*(n_atoms_*n_orbs_*S_);
                            col_OP = row_OP;
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            Hk_(row_,col_) += (Parameters_.UInterOrb[atom_no][alpha][beta])*OPs_.value[index_OP];

                            col_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                    spin_*(n_atoms_*n_orbs_*S_);
                            row_=col_;
                            row_OP = j + (atom_no + n_atoms_*alpha)*(S_)
                                    + (spin_p)*(n_atoms_*n_orbs_*S_);
                            col_OP = row_OP;
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            Hk_(row_,col_) += (Parameters_.UInterOrb[atom_no][alpha][beta])*OPs_.value[index_OP];

                        }
                    }
                }
            }}}



    //Inter-orbital nn
    //Fock
    for(int j=0;j<S_;j++){
        for(int atom_no=0;atom_no<n_atoms_;atom_no++){
            for(int alpha=0;alpha<n_orbs_;alpha++){
                for(int beta=(alpha+1);beta<n_orbs_;beta++){

                    for(int spin_=0;spin_<2;spin_++){

                        //1 and 7 [2 and 8 are h.c. of 1 and 7 ]
                        col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        row_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);

                        col_OP = row_;
                        row_OP = col_;

                        if(col_OP>row_OP){
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = OPs_.value[index_OP];
                        }
                        else{
                            index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = conj(OPs_.value[index_OP]);
                        }
                        //1 and 7
                        Hk_(row_,col_) -= (Parameters_.UInterOrb[atom_no][alpha][beta])*OP_val;
                        //2 and 8
                        Hk_(col_,row_) -= conj((Parameters_.UInterOrb[atom_no][alpha][beta])*OP_val);


                        //3 and 5 [4 and 6 are h.c. of 3 and 5 ]
                        col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        row_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);

                        col_OP = row_;
                        row_OP = col_;


                        if(!NA_spinflip){
                            if(col_OP>row_OP){
                                index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val = OPs_.value[index_OP];
                            }
                            else{
                                index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val = conj(OPs_.value[index_OP]);
                            }
                            //3 and 5
                            Hk_(row_,col_) -= (Parameters_.UInterOrb[atom_no][alpha][beta])*OP_val;
                            //4 and 6
                            Hk_(col_,row_) -= conj((Parameters_.UInterOrb[atom_no][alpha][beta])*OP_val);
                        }
                    }

                }}}}


    //Hunds Coupling
    //Hartree
    for(int j=0;j<S_;j++){
        for(int atom_no=0;atom_no<n_atoms_;atom_no++){
            for(int alpha=0;alpha<n_orbs_;alpha++){
                for(int beta=(alpha+1);beta<n_orbs_;beta++){

                    for(int spin_=0;spin_<2;spin_++){

                        //1 and 3
                        col_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        row_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);

                        col_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);
                        row_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);

                        //NA_spinflip = Parameters_.NoSpinFlipOP && (spin_p!=spin_);
                        if(!NA_spinflip){
                            if(col_OP>row_OP){
                                index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val = OPs_.value[index_OP];
                            }
                            else{
                                index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val = conj(OPs_.value[index_OP]);
                            }
                            //1 and 3
                            Hk_(row_,col_) -= (0.5*Parameters_.JHund[atom_no])*2.0*OP_val;
                        }


                        //2 and 4
                        col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        row_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);

                        col_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);
                        row_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);

                        //NA_spinflip = Parameters_.NoSpinFlipOP && (spin_p!=spin_);
                        if(!NA_spinflip){
                            if(col_OP>row_OP){
                                index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val = OPs_.value[index_OP];
                            }
                            else{
                                index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val = conj(OPs_.value[index_OP]);
                            }
                            //2 and 4
                            Hk_(row_,col_) -= (0.5*Parameters_.JHund[atom_no])*2.0*OP_val;
                        }

                        //5 and 11
                        col_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        row_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);

                        col_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);
                        row_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);

                        if(col_OP>row_OP){
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = OPs_.value[index_OP];
                        }
                        else{
                            index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = conj(OPs_.value[index_OP]);
                        }
                        //5 and 11
                        Hk_(row_,col_) -= (0.5*Parameters_.JHund[atom_no])*OP_val;


                        //6 and 12
                        col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        row_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);

                        col_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);
                        row_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);

                        if(col_OP>row_OP){
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = OPs_.value[index_OP];
                        }
                        else{
                            index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = conj(OPs_.value[index_OP]);
                        }
                        //6 and 12
                        Hk_(row_,col_) -= (0.5*Parameters_.JHund[atom_no])*OP_val;



                        //7 and 9
                        col_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        row_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);

                        col_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);
                        row_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);

                        if(col_OP>row_OP){
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = OPs_.value[index_OP];
                        }
                        else{
                            index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = conj(OPs_.value[index_OP]);
                        }
                        //7 and 9
                        Hk_(row_,col_) += (0.5*Parameters_.JHund[atom_no])*OP_val;


                        //8 and 10
                        col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        row_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);

                        col_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);
                        row_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);

                        if(col_OP>row_OP){
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = OPs_.value[index_OP];
                        }
                        else{
                            index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = conj(OPs_.value[index_OP]);
                        }
                        //8 and 10
                        Hk_(row_,col_) += (0.5*Parameters_.JHund[atom_no])*OP_val;

                    }

                }}}}





    //Hunds Coupling
    //Fock
    for(int j=0;j<S_;j++){
        for(int atom_no=0;atom_no<n_atoms_;atom_no++){
            for(int alpha=0;alpha<n_orbs_;alpha++){
                for(int beta=(alpha+1);beta<n_orbs_;beta++){

                    for(int spin_=0;spin_<2;spin_++){

                        //1 and 3 [hc : 2 and 4]
                        col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        row_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);

                        col_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);
                        row_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);

                        if(col_OP>row_OP){
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = OPs_.value[index_OP];
                        }
                        else{
                            index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = conj(OPs_.value[index_OP]);
                        }
                        //1 and 3
                        Hk_(row_,col_) += (0.5*Parameters_.JHund[atom_no])*2.0*OP_val;
                        //2 and 4
                        Hk_(col_,row_) += conj((0.5*Parameters_.JHund[atom_no])*2.0*OP_val);



                        //5 and 11 [hc : 6 and 12]
                        col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        row_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);

                        col_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);
                        row_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);

                        if(col_OP>row_OP){
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = OPs_.value[index_OP];
                        }
                        else{
                            index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = conj(OPs_.value[index_OP]);
                        }
                        //5 and 11
                        Hk_(row_,col_) += (0.5*Parameters_.JHund[atom_no])*OP_val;
                        //6 and 12
                        Hk_(col_,row_) += conj((0.5*Parameters_.JHund[atom_no])*OP_val);


                        //7 and 9 [hc : 8 and 10]
                        col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        row_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);

                        col_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);
                        row_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);

                        //NA_spinflip = Parameters_.NoSpinFlipOP && (spin_p!=spin_);
                        if(!NA_spinflip){
                            if(col_OP>row_OP){
                                index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val = OPs_.value[index_OP];
                            }
                            else{
                                index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val = conj(OPs_.value[index_OP]);
                            }
                            //7 and 9
                            Hk_(row_,col_) -= (0.5*Parameters_.JHund[atom_no])*OP_val;
                            //8 and 10
                            Hk_(col_,row_) -= conj((0.5*Parameters_.JHund[atom_no])*OP_val);
                        }

                    }

                }}}}




    //Pair Hopping
    //Hartree+Fock
    for(int j=0;j<S_;j++){
        for(int atom_no=0;atom_no<n_atoms_;atom_no++){
            for(int alpha=0;alpha<n_orbs_;alpha++){
                for(int beta=(alpha+1);beta<n_orbs_;beta++){

                    for(int spin_=0;spin_<2;spin_++){

                        //1 and 2 [hc : 3 and 4]
                        col_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        row_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);

                        col_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);
                        row_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);

                        //NA_spinflip = Parameters_.NoSpinFlipOP && (spin_p!=spin_);
                        if(!NA_spinflip){
                            if(col_OP>row_OP){
                                index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val = OPs_.value[index_OP];
                            }
                            else{
                                index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val = conj(OPs_.value[index_OP]);
                            }
                            //1 and 2
                            Hk_(row_,col_) -= (1.0*Parameters_.JHund[atom_no])*OP_val;
                            //3 and 4
                            Hk_(col_,row_) -= conj((1.0*Parameters_.JHund[atom_no])*OP_val);
                        }

                        //5 and 6 [hc : 7 and 8]
                        col_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        row_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);

                        col_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);
                        row_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);

                        if(col_OP>row_OP){
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = OPs_.value[index_OP];
                        }
                        else{
                            index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = conj(OPs_.value[index_OP]);
                        }
                        //5 and 6
                        Hk_(row_,col_) += (1.0*Parameters_.JHund[atom_no])*OP_val;
                        //7 and 8
                        Hk_(col_,row_) += conj((1.0*Parameters_.JHund[atom_no])*OP_val);

                    }
                }}}}




    //     if(k1==0 && k2==0){
    //      cout<<"---------- Printing matrix kx=0=ky -------------"<<endl;
    //      Ham_.print();
    //      cout<<"-----------------"<<endl;
    //  }

    //            cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
    //            Ham_.print();
    //            cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;



}




void Kspace_calculation_G2dLatticeNew::Inverse(Matrix<complex<double>> & A){

    int M=A.n_col();;
    vector<complex<double>> work(M);
    vector<int> ipiv(M);
    int n,info;

    n=M;

    zgetrf_(&n, &n, &(A(0,0)), &n, &(ipiv[0]), &info);

    if(info != 0){
        std::cerr<<"zgetrf failed"<<endl;
    }

    int lwork=64*M;
    work.resize(lwork);
    //cout<<lwork<<endl;
    zgetri_(&n, &(A(0,0)), &n,  &(ipiv[0]), &(work[0]), &lwork, &info);

    if(info != 0){
        std::cerr<<"zgetri failed"<<endl;
    }



}


void Kspace_calculation_G2dLatticeNew::Hall_conductance(){

    Mat_2_doub hall_cond, hall_cond_xy;
    hall_cond.resize(2); hall_cond_xy.resize(2);
    for(int spin_=0;spin_<2;spin_++){
        hall_cond[spin_].resize(2);
        hall_cond_xy[spin_].resize(2);
    }

    double cond_11, cond_22;
    double eps_temp =0.00;
    Parameters_.eta=0.001;

    double FACTOR_;
    FACTOR_ = 1.0;

    string fileout_="sigmaxy_vs_mu.txt";
    ofstream fileout(fileout_.c_str());


    for(int spin_1=0;spin_1<2;spin_1++){
        for(int spin_2=0;spin_2<2;spin_2++){
            for(int m=0;m<lx_*ly_*2*n_atoms_*n_orbs_;m++){
                for(int n=0;n<lx_*ly_*2*n_atoms_*n_orbs_;n++){
                    if(abs(Eigenvalues_saved[m]-Eigenvalues_saved[n])>=eps_temp){
                        hall_cond[spin_1][spin_2] += FACTOR_*((2.0*PI)/(lx_*ly_))*( (1.0/( exp((Eigenvalues_saved[m]-mu_)*Parameters_.beta ) + 1.0))-  (1.0/( exp((Eigenvalues_saved[n]-mu_)*Parameters_.beta ) + 1.0)) )*
                                (1.0/( (Parameters_.eta*Parameters_.eta) + ((Eigenvalues_saved[n]-Eigenvalues_saved[m])*(Eigenvalues_saved[n]-Eigenvalues_saved[m]))   ))*
                                ((1.0*J_KE_e1[spin_1](m,n))*(1.0*J_KE_e2[spin_2](n,m))).imag();
                    }
                }
            }
        }
    }

    for(int spin_1=0;spin_1<2;spin_1++){
        for(int spin_2=0;spin_2<2;spin_2++){
            for(int m=0;m<lx_*ly_*2*n_atoms_*n_orbs_;m++){
                for(int n=0;n<lx_*ly_*2*n_atoms_*n_orbs_;n++){
                    if(abs(Eigenvalues_saved[m]-Eigenvalues_saved[n])>=eps_temp){
                        hall_cond_xy[spin_1][spin_2] += FACTOR_*((2.0*PI)/(lx_*ly_))*( (1.0/( exp((Eigenvalues_saved[m]-mu_)*Parameters_.beta ) + 1.0))-  (1.0/( exp((Eigenvalues_saved[n]-mu_)*Parameters_.beta ) + 1.0)) )*
                                (1.0/( (Parameters_.eta*Parameters_.eta) + ((Eigenvalues_saved[n]-Eigenvalues_saved[m])*(Eigenvalues_saved[n]-Eigenvalues_saved[m]))   ))*
                                ((1.0*J_KE_X[spin_1](m,n))*(1.0*J_KE_Y[spin_2](n,m))).imag();
                    }
                }
            }
        }
    }


    cout<<"========= Hall conductance_12 [s1][s2] =============="<<endl;
    for(int spin_1=0;spin_1<2;spin_1++){
        for(int spin_2=0;spin_2<2;spin_2++){
            cout<<hall_cond[spin_1][spin_2]<<"   ";
        }
        cout<<endl;
    }
    cout<<endl;


    cout<<"========= Hall conductance_XY[s1][s2] =============="<<endl;
    for(int spin_1=0;spin_1<2;spin_1++){
        for(int spin_2=0;spin_2<2;spin_2++){
            cout<<hall_cond_xy[spin_1][spin_2]<<"   ";
        }
        cout<<endl;
    }
    cout<<endl;


    //    cout<<"Hall Conductance_12 upup = "<<hall_cond[0]<<endl;
    //    cout<<"Hall Conductance_12 dndn = "<<hall_cond[1]<<endl;
    //    cout<<"Hall Conductance_XY upup = "<<hall_cond_xy[0]<<endl;
    //    cout<<"Hall Conductance_XY dndn = "<<hall_cond_xy[1]<<endl;

    fileout<<"#mu  hall_cond_upup    hall_cond_xy_upup     cond_11_upup    cond_22_upup"<<endl;
    double mu=Eigenvalues_[0]-1.0;

    /*
    while(mu<Eigenvalues_[Eigenvalues_.size()-1]+1.0){
        hall_cond=0.0;
        hall_cond_xy=0.0;
        cond_11=0.0;
        cond_22=0.0;

        for(int m=0;m<lx_*ly_*2*n_orbs_;m++){
            for(int n=0;n<lx_*ly_*2*n_orbs_;n++){
                if(abs(Eigenvalues_saved[m]-Eigenvalues_saved[n])>=eps_temp){
                    hall_cond += FACTOR_*((2.0*PI)/(lx_*ly_))*( (1.0/( exp((Eigenvalues_saved[m]-mu)*Parameters_.beta ) + 1.0))-  (1.0/( exp((Eigenvalues_saved[n]-mu)*Parameters_.beta ) + 1.0)) )*
                            (1.0/( (Parameters_.eta*Parameters_.eta) + ((Eigenvalues_saved[n]-Eigenvalues_saved[m])*(Eigenvalues_saved[n]-Eigenvalues_saved[m]))   ))*
                            ((1.0*J_KE_e1(m,n))*(1.0*J_KE_e2(n,m))).imag();
                    hall_cond_xy += FACTOR_*((2.0*PI)/(lx_*ly_))*( (1.0/( exp((Eigenvalues_saved[m]-mu)*Parameters_.beta ) + 1.0))-  (1.0/( exp((Eigenvalues_saved[n]-mu)*Parameters_.beta ) + 1.0)) )*
                            (1.0/( (Parameters_.eta*Parameters_.eta) + ((Eigenvalues_saved[n]-Eigenvalues_saved[m])*(Eigenvalues_saved[n]-Eigenvalues_saved[m]))   ))*
                            ((1.0*J_KE_X(m,n))*(1.0*J_KE_Y(n,m))).imag();
                    cond_11 += FACTOR_*((2.0*PI)/(lx_*ly_))*( (1.0/( exp((Eigenvalues_saved[m]-mu)*Parameters_.beta ) + 1.0))-  (1.0/( exp((Eigenvalues_saved[n]-mu)*Parameters_.beta ) + 1.0)) )*
                            (1.0/( (Parameters_.eta*Parameters_.eta) + ((Eigenvalues_saved[n]-Eigenvalues_saved[m])*(Eigenvalues_saved[n]-Eigenvalues_saved[m]))   ))*
                            ((1.0*J_KE_e1(m,n))*(1.0*J_KE_e1(n,m))).imag();
                    cond_22 += FACTOR_*((2.0*PI)/(lx_*ly_))*( (1.0/( exp((Eigenvalues_saved[m]-mu)*Parameters_.beta ) + 1.0))-  (1.0/( exp((Eigenvalues_saved[n]-mu)*Parameters_.beta ) + 1.0)) )*
                            (1.0/( (Parameters_.eta*Parameters_.eta) + ((Eigenvalues_saved[n]-Eigenvalues_saved[m])*(Eigenvalues_saved[n]-Eigenvalues_saved[m]))   ))*
                            ((1.0*J_KE_e2(m,n))*(1.0*J_KE_e2(n,m))).imag();

                }
            }
        }
        fileout<<mu<<"  "<<hall_cond<<"  "<<hall_cond_xy<<"  "<<cond_11<<"  "<<cond_22<<endl;
        mu=mu+0.05;
    }
*/

}


void Kspace_calculation_G2dLatticeNew::Create_Current_Oprs(){
    /*

    double spin_factor;


    J_KE_e1[0].resize(lx_*ly_*2*n_orbs_, lx_*ly_*2*n_orbs_);
    J_KE_e2[0].resize(lx_*ly_*2*n_orbs_, lx_*ly_*2*n_orbs_);

    J_KE_X.resize(lx_*ly_*2*n_orbs_, lx_*ly_*2*n_orbs_);
    J_KE_Y.resize(lx_*ly_*2*n_orbs_, lx_*ly_*2*n_orbs_);


    int row_, col_;
    int d1_, d2_;
    int d1_org, d2_org;
    int alpha_1, alpha_2, alpha_p_1, alpha_p_2;
    int c1, c2;



    int S_ = Nsites_MUC;
    int k_index;
    int state_p,state_n;
    int sigma, sigma_p;
    sigma=0;sigma_p=0;
    int d1_net, d2_net;
    double fac1, fac2, facx, facy, dx_, dy_;

    //2*n_orbs_*S_*k_index + row
    for(int l=0;l<2*n_orbs_*S_;l++){
        for(int m=0;m<2*n_orbs_*S_;m++){

            for(int k1=0;k1<lx_cells;k1++){
                for(int k2=0;k2<ly_cells;k2++){
                    k_index = Coordinates_.Ncell(k1,k2);
                    state_p = 2*n_orbs_*S_*k_index + l;
                    state_n = 2*n_orbs_*S_*k_index + m;

                    for(int sigma=0;sigma<2;sigma++){


                        for(int alpha=0;alpha<S_;alpha++){
                            for(int gamma=0;gamma<2;gamma++){

                                for(int alpha_p=0;alpha_p<S_;alpha_p++){
                                    for(int gamma_p=0;gamma_p<2;gamma_p++){

                                        //alpha + gamma*(S_)
                                        // + spin1*(n_orbs_*S_);
                                        sigma_p=sigma;
                                        c1 = alpha + gamma*(S_) + sigma*(2*S_);
                                        c2 = alpha_p + gamma_p*(S_) + sigma*(2*S_);

                                        alpha_1 = alpha % UnitCellSize_x; //a = a1 +a2*lx
                                        alpha_p_1 = alpha_p % UnitCellSize_x;
                                        alpha_2 = (alpha - alpha_1)/UnitCellSize_x;
                                        alpha_p_2 = (alpha_p - alpha_p_1)/UnitCellSize_x;

                                        for(int cell_no=0;cell_no<ncells_;cell_no++){

                                            //------------------------


                                            d1_org = Coordinates_.indx_cellwise(cell_no);
                                            d2_org = Coordinates_.indy_cellwise(cell_no);

                                            row_ = gamma  + ( (d1_org*UnitCellSize_x) + alpha_1)*n_orbs_
                                                    + (( (d2_org*UnitCellSize_y) + alpha_2)*lx_*n_orbs_)
                                                    + sigma*(lx_*ly_*n_orbs_);
                                            col_ = ( gamma_p + (0 + alpha_p_1)*n_orbs_ + ((0 + alpha_p_2)*lx_*n_orbs_)) + sigma_p*(lx_*ly_*n_orbs_);

                                            Get_minimum_distance_direction(0, cell_no, d1_, d2_);

                                            d1_net = ( (d1_*UnitCellSize_x) + alpha_1) - (0 + alpha_p_1);
                                            d2_net = ((d2_*UnitCellSize_y) + alpha_2) - (0 + alpha_p_2);


                                            dx_ = 1.0*d1_net + 0.5*(d2_net);
                                            dy_ = (sqrt(3.0)/2.0)*d2_net;
                                            if(((dx_*dx_) + (dy_*dy_))>0.00001){
                                                facx = (dx_/sqrt((dx_*dx_) + (dy_*dy_)));
                                                facy = (dy_/sqrt((dx_*dx_) + (dy_*dy_)));
                                            }
                                            else{
                                                facx=0;
                                                facy=0;
                                            }
                                            //                                            if(dx_==0){
                                            //                                                facx=0;
                                            //                                            }
                                            //                                            else{
                                            //                                                facx=dx_/abs(dx_);
                                            //                                            }
                                            //                                            if(dy_==0){
                                            //                                                facy=0;
                                            //                                            }
                                            //                                            else{
                                            //                                                facy=dy_/abs(dy_);
                                            //                                            }



                                            //                                    fac1 = (d1_net/sqrt((d1_net*d1_net) + (d2_net*d2_net)));
                                            //                                    fac2 = (d2_net/sqrt((d1_net*d1_net) + (d2_net*d2_net)));

                                            if(d1_net==0){
                                                fac1=0;
                                            }
                                            else{
                                                fac1=d1_net/abs(d1_net);
                                            }
                                            if(d2_net==0){
                                                fac2=0;
                                            }
                                            else{
                                                fac2=d2_net/abs(d2_net);
                                            }

                                            if(abs(Connections_.HTB_(row_,col_)) > 0.0000001){
                                                J_KE_X(state_p,state_n) += (1.0)*facx*iota_complex*(
                                                            (Connections_.HTB_(row_,col_)*conj(Eigvectors_[state_p][c2])*Eigvectors_[state_n][c1]*exp(iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_ )))
                                                            -(Connections_.HTB_(col_,row_)*conj(Eigvectors_[state_p][c1])*Eigvectors_[state_n][c2]*exp(-1.0*iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_ )))
                                                            );

                                                J_KE_Y(state_p,state_n) += (1.0)*facy*iota_complex*(
                                                            (Connections_.HTB_(row_,col_)*conj(Eigvectors_[state_p][c2])*Eigvectors_[state_n][c1]*exp(iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_ )))
                                                            -(Connections_.HTB_(col_,row_)*conj(Eigvectors_[state_p][c1])*Eigvectors_[state_n][c2]*exp(-1.0*iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_ )))
                                                            );

                                                J_KE_e1(state_p,state_n) += (1.0)*fac1*iota_complex*(
                                                            (Connections_.HTB_(row_,col_)*conj(Eigvectors_[state_p][c2])*Eigvectors_[state_n][c1]*exp(iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_ )))
                                                            -(Connections_.HTB_(col_,row_)*conj(Eigvectors_[state_p][c1])*Eigvectors_[state_n][c2]*exp(-1.0*iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_ )))
                                                            );

                                                J_KE_e2(state_p,state_n) += (1.0)*fac2*iota_complex*(
                                                            (Connections_.HTB_(row_,col_)*conj(Eigvectors_[state_p][c2])*Eigvectors_[state_n][c1]*exp(iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_ )))
                                                            -(Connections_.HTB_(col_,row_)*conj(Eigvectors_[state_p][c1])*Eigvectors_[state_n][c2]*exp(-1.0*iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_ )))
                                                            );

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

        cout<<"Current Ops loop (" <<2*n_orbs_*S_<<"): l="<<l<<" done"<<endl;
    }



*/
}


void Kspace_calculation_G2dLatticeNew::Create_Current_Oprs_Faster(){}
void Kspace_calculation_G2dLatticeNew::Calculate_ChernNumbers(){}

double Kspace_calculation_G2dLatticeNew::Lorentzian(double x, double brd)
{
    double temp;

    temp = (1.0 / PI) * ((brd / 2.0) / ((x * x) + ((brd * brd) / 4.0)));

    return temp;
}


void Kspace_calculation_G2dLatticeNew::Calculate_Akw(){
    double omega_min=Parameters_.omega_min;
    double omega_max=Parameters_.omega_max;
    double d_omega=Parameters_.d_omega;
    double eta=Parameters_.eta;
    int N_omega = int(((omega_max-omega_min)/d_omega) + 0.5) + 1;


    int S_ = NSites_in_MUC;

    AkwMat.resize(2*n_atoms_*n_orbs_*S_);
    for(int c1=0;c1<2*n_atoms_*n_orbs_*S_;c1++){
        AkwMat[c1].resize(2*n_atoms_*n_orbs_*S_);
        for(int c2=0;c2<2*n_atoms_*n_orbs_*S_;c2++){
            AkwMat[c1][c2].resize(lx_cells*ly_cells);
            for(int qind=0;qind<lx_cells*ly_cells;qind++){
                AkwMat[c1][c2][qind].resize(N_omega);
            }
        }
    }


    // int k_index_m = Coordinates_.Ncell(k1_ind_m,k2_ind_m);
    // int state_k_m = 2*n_atoms_*n_orbs_*S_*k_index_m + m;
    // ABmat[IntraMUC_indexA][IntraMUC_indexB][state_k_n][state_k_m] += A_elmt*B_elmt
    //*conj(Eigvectors_[state_k_n][alpha_comp])


    for(int c1=0;c1<2*n_atoms_*n_orbs_*S_;c1++){
        for(int c2=0;c2<2*n_atoms_*n_orbs_*S_;c2++){
            for(int qind=0;qind<lx_cells*ly_cells;qind++){
                for(int omega_ind=0;omega_ind<N_omega;omega_ind++){
                    double omega_val = omega_min + d_omega*omega_ind;
                    AkwMat[c1][c2][qind][omega_ind]=0.0;

                    for(int n=0;n<2*n_atoms_*n_orbs_*S_;n++){
                        int state_q_n = 2*n_atoms_*n_orbs_*S_*qind + n;
                        AkwMat[c1][c2][qind][omega_ind] += conj(Eigvectors_[state_q_n][c1])*(Eigvectors_[state_q_n][c2])*
                                Lorentzian(omega_val - Eigenvalues_saved[state_q_n],eta);
                    }
                }
            }}}




    //Transforming to Normal Brillioun Zone
    Mat_1_intpair q_path;
    //Create_K_Path_NBZ("MGKMpY_HXGBZ", q_path);
    Create_K_Path_NBZ("GMXGYMG", q_path);
    //Create_K_Path_NBZ("FullBZ", q_path);




    int q1_ind, q2_ind;

    Mat_2_Complex_doub ChiRPA_AB_q_omega_NBZ;

    int n1_p, n2_p;

    for(int atom_no_alpha=0;atom_no_alpha<n_atoms_;atom_no_alpha++){
        for(int alpha_orb=0;alpha_orb<n_orbs_;alpha_orb++){
            for(int spin_alpha=0;spin_alpha<2;spin_alpha++){

                for(int atom_no_beta=0;atom_no_beta<n_atoms_;atom_no_beta++){
                    for(int beta_orb=0;beta_orb<n_orbs_;beta_orb++){
                        for(int spin_beta=0;spin_beta<2;spin_beta++){

                            complex<double> Akw_temp=0.0;
                            string File_Out_Akw_str = "Akw_NBZ_atom_orb_spin"
                                    + to_string(atom_no_alpha)+"_"
                                    + to_string(alpha_orb)+"_"
                                    + to_string(spin_alpha)+"_"
                                    + to_string(atom_no_beta)+"_"
                                    + to_string(beta_orb)+"_"
                                    + to_string(spin_beta)+"_"
                                    +".txt";

                            ofstream file_out_Akw(File_Out_Akw_str.c_str());
                            file_out_Akw<<"#q omega Akw(q,omega)"<<endl;

                            for(int q_ind_temp=0;q_ind_temp<q_path.size();q_ind_temp++){
                                //cout<<"ChiRPA_q_omega (NBZ) "<<q_ind_temp<<"("<<q_path.size()<<")"<<endl;

                                q1_ind=q_path[q_ind_temp].first;
                                q2_ind=q_path[q_ind_temp].second;
                                int q_index = q1_ind + q2_ind*lx_;

                                assert( (Mat_MUC(0,0)*lx_cells)%(lx_)==0);
                                assert( (Mat_MUC(0,1)*lx_cells)%(ly_)==0);
                                assert( (Mat_MUC(1,0)*ly_cells)%(lx_)==0);
                                assert( (Mat_MUC(1,1)*ly_cells)%(ly_)==0);


                                n1_p = int(((q1_ind*Mat_MUC(0,0)*lx_cells)/(lx_)) + 0.5)
                                        + int(((q2_ind*Mat_MUC(0,1)*lx_cells)/(ly_)) + 0.5);
                                n1_p = (n1_p + lx_cells)%lx_cells;

                                n2_p = int(((q1_ind*Mat_MUC(1,0)*ly_cells)/(lx_)) + 0.5)
                                        + int(((q2_ind*Mat_MUC(1,1)*ly_cells)/(ly_)) + 0.5);
                                n2_p = (n2_p + ly_cells)%ly_cells;

                                int qp_ind = Coordinates_.Ncell(n1_p, n2_p);

                                for(int omega_ind=0;omega_ind<N_omega;omega_ind++){
                                    double omega_val = omega_min + d_omega*omega_ind;

                                    // ChiRPA_AB_q_omega_NBZ[q_index][omega_ind]=0.0;

                                    Akw_temp=0.0;
                                    for(int g=0;g<NSites_in_MUC;g++){
                                        int g1=Intra_MUC_positions[g].first;
                                        int g2=Intra_MUC_positions[g].second;
                                        for(int h=0;h<NSites_in_MUC;h++){
                                            int h1=Intra_MUC_positions[h].first;
                                            int h2=Intra_MUC_positions[h].second;


                                            //alpha
                                            int ind1 = h + (atom_no_alpha + n_atoms_*alpha_orb)*(S_) +
                                                    spin_alpha*(n_atoms_*n_orbs_*S_);
                                            int ind2 = g + (atom_no_beta + n_atoms_*beta_orb)*(S_) +
                                                    spin_beta*(n_atoms_*n_orbs_*S_);

                                            Akw_temp += (1.0/(NSites_in_MUC))*exp(iota_complex*2.0*PI*( (1.0*q1_ind*(g1-h1))/(lx_) + (1.0*q2_ind*(g2-h2))/(ly_) ))*
                                                    AkwMat[ind1][ind2][qp_ind][omega_ind];
                                        }
                                    }

                                    file_out_Akw<<q_ind_temp<<"  "<<q1_ind<<"  "<<q2_ind<<"  "<<omega_val<<"   "<<Akw_temp.real()<<"  "<<Akw_temp.imag()<<endl;

                                    //file_out_BarreSusc<<q_ind_temp<<"  "<<q1_ind<<"  "<<q2_ind<<"  "<<omega<<"  ";
                                    //file_out_BarreSusc<<ChiRPA_AB_q_omega_NBZ[q_index][omega_ind].real()<<"  "<<ChiRPA_AB_q_omega_NBZ[q_index][omega_ind].imag()<<endl;


                                }
                                file_out_Akw<<endl;
                                //file_out_BarreSusc<<endl;

                            }


                        }}}
            }}}





}

void Kspace_calculation_G2dLatticeNew::Calculate_Nw()
{

    char temp_char[50];
    sprintf(temp_char, "%.10f", Parameters_.Temperature);

    int ind_max;
    ind_max=Eigenvalues_.size();

    //---------Read from input file-----------------------//
    string fileout = "Nw" + string(temp_char)+ ".txt";
    double omega_min, omega_max, d_omega;
    double eta = 0.1;
    omega_min = Eigenvalues_[0]-10.0;
    omega_max = Eigenvalues_[ind_max-1] + 10.0;
    d_omega = 0.001;
    //---------------------------------------------------//

    int omega_index_max = int((omega_max - omega_min) / (d_omega));

    ofstream file_Nw_out(fileout.c_str());

    double temp_val;
    double filling_counter;
    filling_counter=0.0;
    for (int omega_ind = 0; omega_ind < omega_index_max; omega_ind++)
    {

        temp_val = 0.0;

        for (int n = 0; n < Eigenvalues_.size(); n++)
        {

            temp_val += Lorentzian(omega_min + (omega_ind * d_omega) - Eigenvalues_[n], eta);
        }

        filling_counter +=temp_val*(d_omega)*(1.0/Eigenvalues_.size());

        file_Nw_out << omega_min + (omega_ind * d_omega) << "     " << filling_counter<<"    "<<temp_val << "     " << endl;
    }

    file_Nw_out << "#mu = " << mu_ << endl;



    //     string fileout_eigs = "Spectrum.txt";
    //     ofstream file_eigs_out(fileout_eigs.c_str());
    //     for (int n = 0; n < Hamiltonian_.eigs_.size(); n++)
    //     {
    //         file_eigs_out<<n<<"  "<<Hamiltonian_.eigs_[n]<<endl;
    //     }



}


void Kspace_calculation_G2dLatticeNew::Get_minimum_distance_direction(int l,int m,int &r1_, int &r2_){

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
        r1_b=-1*(lx_cells-r1_a);
    }
    else if (r1_a<0){
        r1_b=lx_cells-abs(r1_a);
    }
    else{
        r1_b=0;
        assert(r1_a==0);
    }

    r2_a=m2_pos-l2_pos;
    if(r2_a>0){
        r2_b=-1*(ly_cells-r2_a);
    }
    else if (r2_a<0){
        r2_b=ly_cells-abs(r2_a);
    }
    else{
        r2_b=0;
        assert(r2_a==0);
    }


    double min_dis=1000000.0;
    double dis;

    //r1a r2a
    rx_ = r1_a;//((1.0)*(r1_a) +  (1.0/2.0)*(r2_a));
    ry_ =  r2_a;//(0.0*(r1_a) + (sqrt(3.0)/2.0)*(r2_a));
    dis= sqrt(rx_*rx_ + ry_*ry_);
    if(dis<=min_dis){
        r1_=r1_a;
        r2_=r2_a;
        min_dis=dis;
    }

    //r1a r2b
    rx_ = r1_a;//((1.0)*(r1_a) +  (1.0/2.0)*(r2_b));
    ry_ =  r2_b;//(0.0*(r1_a) + (sqrt(3.0)/2.0)*(r2_b));
    dis= sqrt(rx_*rx_ + ry_*ry_);
    if(dis<=min_dis){
        r1_=r1_a;
        r2_=r2_b;
        min_dis=dis;
    }

    //r1b r2a
    rx_ = r1_b;//((1.0)*(r1_b) +  (1.0/2.0)*(r2_a));
    ry_ =  r2_a;//(0.0*(r1_b) + (sqrt(3.0)/2.0)*(r2_a));
    dis= sqrt(rx_*rx_ + ry_*ry_);
    if(dis<=min_dis){
        r1_=r1_b;
        r2_=r2_a;
        min_dis=dis;
    }

    //r1b r2b
    rx_ = r1_b;//((1.0)*(r1_b) +  (1.0/2.0)*(r2_b));
    ry_ =  r2_b;//(0.0*(r1_b) + (sqrt(3.0)/2.0)*(r2_b));
    dis= sqrt(rx_*rx_ + ry_*ry_);
    if(dis<=min_dis){
        r1_=r1_b;
        r2_=r2_b;
        min_dis=dis;
    }

}


void Kspace_calculation_G2dLatticeNew::Get_Bands(){


    int S_=NSites_in_MUC;

    char temp_char[50];
    sprintf(temp_char, "%.10f", Parameters_.Temperature);

    string File_Out_Bands;
    File_Out_Bands = "Bands" + string(temp_char)+ ".txt";
    ofstream file_out_bands(File_Out_Bands.c_str());


    string File_Out_Bands3;
    File_Out_Bands3 = "Bands_Grid" + string(temp_char)+ ".txt";
    ofstream file_out_bands3(File_Out_Bands3.c_str());


    string File_Out_Bands2;
    File_Out_Bands2 = "Bands_Path2_" + string(temp_char) +".txt";
    ofstream file_out_bands2(File_Out_Bands2.c_str());


    string File_Out_Bands4;
    File_Out_Bands4 = "Bands_Path4_" + string(temp_char) +".txt";
    ofstream file_out_bands4(File_Out_Bands4.c_str());

    string File_Out_Bands5;
    File_Out_Bands5 = "Bands_Path5_" + string(temp_char) +".txt";
    ofstream file_out_bands5(File_Out_Bands5.c_str());

    int k_index;
    int kx_i, ky_i;

    Mat_1_intpair k_path;
    k_path.clear();
    pair_int temp_pair;


    // ---k_path---------

    //--------\Gamma to M-----------------
    ky_i = 0;
    for (kx_i = 0; kx_i <= (lx_cells/ 2); kx_i++)
    {
        temp_pair.first = kx_i;
        temp_pair.second = kx_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------


    //--------M to X-----------------
    kx_i = (lx_cells / 2);
    for (ky_i = (ly_cells / 2)-1; ky_i >=0; ky_i--)
    {
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------


    //--------X to Gamma-----------------
    ky_i = 0;
    for (kx_i = (lx_cells / 2)-1; kx_i >=0; kx_i--)
    {
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------


    //--------Gamma to Y-----------------
    kx_i = 0;
    ky_i = 0;
    for (ky_i = 1; ky_i <=(ly_cells/2); ky_i++)
    {
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------

    //--------Y to M-----------------
    kx_i = 0;
    ky_i = ly_cells/2;
    for (kx_i = 1; kx_i <=(lx_cells/2); kx_i++)
    {
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------


    //--------M to \Gamma[with one extra point,
    //                  because in gnuplot use "set pm3d corners2color c1"
    //                  ]-----------------
    kx_i = (lx_cells / 2) - 1;
    ky_i = (ly_cells / 2) - 1;
    for (kx_i = (lx_cells / 2) - 1; kx_i >= 0; kx_i--)
    {
        ky_i=kx_i;
        if(ky_i < ly_cells){
            temp_pair.first = kx_i;
            temp_pair.second = kx_i;
            k_path.push_back(temp_pair);
        }
    }

    temp_pair.first = 0;
    temp_pair.second = 0;
    k_path.push_back(temp_pair);





    //Path5------------------------
    Mat_1_intpair k_path5;
    k_path5.clear();
    
    //--------\Gamma to K-----------------
    kx_i = 0;
    ky_i = 0;
    for (ky_i = 0; ky_i <= ((ly_cells)/3); ky_i++)
    {
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        kx_i = kx_i +2;
        k_path5.push_back(temp_pair);
    }


    //K(2l/3, l/3) to M(l/2, l/2)
    kx_i = ((2*(ly_cells)/3))-1;
    ky_i = (((ly_cells)/3))+1;
    for (ky_i = (((ly_cells)/3) +1); ky_i<=(ly_cells/2); ky_i++)
    {
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        kx_i=kx_i-1;
        k_path5.push_back(temp_pair);
    }


    //M to Gamma
    kx_i = (lx_cells/2)-1;
    ky_i = (ly_cells/2)-1;
    for (kx_i = (((lx_cells)/2) -1); kx_i>=0; kx_i--)
    {
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        ky_i=ky_i-1;
        k_path5.push_back(temp_pair);
    }
    


    //-------------------------------------



    //----------------------------
    int counter_Gamma1, counter_K,counter_Kprime, counter_Gamma2, counter_M;
    Mat_1_intpair k_path4;
    k_path4.clear();


    counter_Gamma1=0;
    //Gamma to K
    kx_i = 0;
    ky_i = 0;
    for (ky_i = 0; ky_i <= ((ly_cells)/3); ky_i++)
    {
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        kx_i = kx_i +2;
        k_path4.push_back(temp_pair);
    }
    counter_K= k_path4.size()-1;


    //K to K'
    kx_i = ((2*lx_cells)/3) -1;
    for (ky_i = ((ly_cells)/3)+1; ky_i <= ((2*ly_cells)/3); ky_i++)
    {
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        kx_i = kx_i -1;
        k_path4.push_back(temp_pair);
    }
    counter_Kprime = k_path4.size()-1;

    //K' to Gamma
    kx_i = ((lx_cells)/3) -1;
    ky_i = ((2*ly_cells)/3) -2;
    for (kx_i = ((lx_cells)/3)-1; kx_i >=0; kx_i--)
    {
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        ky_i = ky_i -2;
        k_path4.push_back(temp_pair);
    }
    counter_Gamma2= k_path4.size()-1;


    //Gamma to mu(M)
    kx_i = 1;
    ky_i=1;
    for (ky_i = 1; ky_i <= ((ly_cells)/2); ky_i++)
    {
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        kx_i = kx_i +1;
        k_path4.push_back(temp_pair);
    }

    counter_M= k_path4.size()-1;

    //------------------------------






    //k_path2
    //----------------------------------

    int n1, n2, n1_, n2_;
    Mat_1_intpair k_path2;
    k_path2.clear();
    int L1_,L2_;
    L1_=lx_cells;
    L2_=ly_cells;

    //K+' to K-
    n1=int(2*L1_/3);
    n2=int(L2_/3);
    while(n2>=int(-L2_/3)){
        n1_ = (n1 + L1_)%L1_;
        n2_ = (n2 + L2_)%L2_;
        temp_pair.first = n1_;
        temp_pair.second = n2_;
        k_path2.push_back(temp_pair);
        n2--;
        n1=int(2*L1_/L2_)*n2;
    }

    //K- to K+
    n1=int(-2*L1_/3);
    n2=int(-L2_/3);
    n2--;n1++;
    while(n1<=int(-L1_/3)){
        n1_ = (n1 + L1_)%L1_;
        n2_ = (n2 + L2_)%L2_;
        temp_pair.first = n1_;
        temp_pair.second = n2_;
        k_path2.push_back(temp_pair);
        n2--;
        n1++;
    }

    //K+ to K+'
    n1=int(-L1_/3);
    n2=int(-2*L2_/3);
    n2=n2+2; //in principle it should be n2=n2+1, n1=n1+1
    n1=n1+2;
    while(n1<=int(2*L1_/3)){
        n1_ = (n1 + L1_)%L1_;
        n2_ = (n2 + L2_)%L2_;
        temp_pair.first = n1_;
        temp_pair.second = n2_;
        k_path2.push_back(temp_pair);
        n2=n2+2;  //in principle it should be n2=n2+1, n1=n1+1
        n1=n1+2;
    }

    //----------------------------------




    cout<<"PRINTING PATH"<<endl;
    for (int k_point = 0; k_point < k_path.size(); k_point++)
    {
        cout<<k_path[k_point].first<< "   "<<k_path[k_point].second<<endl;
    }
    //----k_path done-------


    file_out_bands<<"#kindex  E(k,1)  orb=0,UP  orb=0,DOWN orb=1,UP  orb=1,DOWN ... E(k,2)  UP DOWN  ....."<<endl;
    for (int k_point = 0; k_point < k_path.size(); k_point++)
    {
        k_index=Coordinates_.Ncell(k_path[k_point].first, k_path[k_point].second);

        //Eigvectors_[state][c1]
        //k_index = Coordinates_.Ncell(k1,k2);
        //state = 2*n_atoms_*n_orbs_*S_*k_index + n;
        //c1 = alpha_row + gamma_row*(S_) +  sigma_row*(n_atoms_*n_orbs_*S_);
        file_out_bands<<k_point<<"   ";
        for(int band=0;band<2*n_orbs_*n_atoms_*S_;band++){

            int state = 2*n_atoms_*n_orbs_*S_*k_index + band;
            file_out_bands<<Eigenvalues_[state];
            double up_contr =0;
            double dn_contr =0;
            int c_up, c_dn;
            for(int gamma_temp=0;gamma_temp<n_atoms_*n_orbs_;gamma_temp++){
                up_contr=0.0;
                dn_contr=0.0;
                for(int alpha_temp=0;alpha_temp<S_;alpha_temp++){
                    c_up = alpha_temp + gamma_temp*(S_) +  0*(n_atoms_*n_orbs_*S_);
                    c_dn = alpha_temp + gamma_temp*(S_) +  1*(n_atoms_*n_orbs_*S_);

                    up_contr += abs(Eigvectors_[state][c_up])*abs(Eigvectors_[state][c_up]);
                    dn_contr += abs(Eigvectors_[state][c_dn])*abs(Eigvectors_[state][c_dn]);

                }
                file_out_bands<<"   "<<up_contr<<"   "<<dn_contr;
            }

            file_out_bands<<"   ";
        }
        file_out_bands<<endl;
    }


    // for (int k_point = 0; k_point < k_path2.size(); k_point++)
    // {
    //     k_index=Coordinates_.Ncell(k_path2[k_point].first, k_path2[k_point].second);

    //     file_out_bands2<<k_point<<"   ";
    //     for(int band=0;band<2*n_orbs_*n_atoms_*S_;band++){
    //         file_out_bands2<<Eigenvalues_[2*n_orbs_*n_atoms_*S_*k_index + band]<<"   ";
    //     }
    //     file_out_bands2<<endl;
    // }



    // for (int k_point = 0; k_point < k_path4.size(); k_point++)
    // {
    //     k_index=Coordinates_.Ncell(k_path4[k_point].first, k_path4[k_point].second);

    //     file_out_bands4<<k_point<<"   ";
    //     for(int band=0;band<2*n_orbs_*n_atoms_*S_;band++){
    //         file_out_bands4<<Eigenvalues_[2*n_orbs_*n_atoms_*S_*k_index + band]<<"   ";
    //     }
    //     file_out_bands4<<endl;
    // }

    //file_out_bands4<<"# Gamma(" <<counter_Gamma1<<")"<<"-->"<<"K(" <<counter_K<<")"<<"-->"<<"K_prime(" <<counter_Kprime<<")"<<"-->"<<"Gamma(" <<counter_Gamma2<<")"<<"-->"<<"M(" <<counter_M<<")"<<endl;




    // for (int k_point = 0; k_point < k_path5.size(); k_point++)
    // {
    //     k_index=Coordinates_.Ncell(k_path5[k_point].first, k_path5[k_point].second);

    //     file_out_bands5<<k_point<<"   "<<k_path5[k_point].first<<"   "<<k_path5[k_point].second<<"   ";
    //     for(int band=0;band<2*n_orbs_*n_atoms_*S_;band++){
    //         file_out_bands5<<Eigenvalues_[2*n_orbs_*n_atoms_*S_*k_index + band]<<"   ";
    //     }
    //     file_out_bands5<<endl;
    // }

    //file_out_bands5<<"# Gamma"<<"-->"<<"K-->"<<"M-->"<<"Gamma"<<endl;



    for (int ky_ind_=0;ky_ind_<ly_cells;ky_ind_++)
    {
        for (int kx_ind_=0;kx_ind_<lx_cells;kx_ind_++)
        {
            k_index=Coordinates_.Ncell(kx_ind_, ky_ind_);

            file_out_bands3<<k_index<<"   "<<kx_ind_<<"   "<<ky_ind_<<"   ";
            for(int band=0;band<2*n_orbs_*n_atoms_*S_;band++){
                file_out_bands3<<Eigenvalues_[2*n_orbs_*n_atoms_*S_*k_index + band]<<"   ";
            }
            file_out_bands3<<endl;
        }
        file_out_bands3<<endl;
    }



}


void Kspace_calculation_G2dLatticeNew::Get_Energies_new(){



    int S_=NSites_in_MUC;
    complex<double> E_class_temp=0.0;
    E_class= 0.0;

    E_class_onsite_U0_Hartree=0.0; E_class_longrange_Hartree=0.0;
    E_class_onsite_U0_Fock=0.0; E_class_longrange_Fock=0.0;


    int row_, col_, index, rowM_, colM_;
    int row_OP, col_OP, index_OP, index_OP_p;

    int UP_, DOWN_;
    UP_=0;DOWN_=1;
    int spin_bar, gamma_bar;
    int k_index;
    double k1x_val, k2x_val, k1y_val, k2y_val;
    double fac_;
    complex<double> OP_value, OP_value2;
    int org_site_nx, org_site_ny, org_site_npx, org_site_npy;
    int rel_org_site_npjpx, rel_org_site_npjpy;
    int rel_np1, rel_np2, rel_jpx, rel_jpy;
    int rel_cell_np1, rel_cell_np2;

    int col_temp, row_temp, cellx_new, celly_new, col_temp_new, row_temp_new, comp_temp ;
    int col_temp_up, row_temp_up, col_temp_dn, row_temp_dn;
    int comp_temp_up, comp_temp_dn;
    complex<double> value_new;

    complex<double> OP_val, OP_val_p;

    bool NA_spinflip=Parameters_.NoSpinFlipOP;

    //Interaction local intra-orbital U0:
    for(int j=0;j<S_;j++){
        for(int b_beta=0;b_beta<n_atoms_;b_beta++){
            for(int gamma=0;gamma<n_orbs_;gamma++){
                int atom_plus_orb=b_beta +  n_atoms_*gamma;
                for(int spin_=0;spin_<2;spin_++){
                    col_ = j + atom_plus_orb*(S_) +
                            spin_*(n_atoms_*n_orbs_*S_);

                    for(int spin_p=0;spin_p<2;spin_p++){

                        row_ = j + atom_plus_orb*(S_) +
                                spin_p*(n_atoms_*n_orbs_*S_);

                        if(col_>row_){
                            index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val_p = OPs_.value[index_OP_p];
                        }
                        else{
                            index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val_p = conj(OPs_.value[index_OP_p]);
                        }


                        if(spin_==spin_p){
                            row_OP = j + atom_plus_orb*(S_)
                                    + (1-spin_)*(n_atoms_*n_orbs_*S_);
                            col_OP = row_OP;
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];

                            E_class_temp += Parameters_.U0[b_beta]*OPs_.value[index_OP]*OP_val_p;
                        }
                        else{
                            if(!Parameters_.Just_Hartree){
                                row_OP = j + atom_plus_orb*(S_)
                                        + (spin_)*(n_atoms_*n_orbs_*S_);
                                col_OP = j + atom_plus_orb*(S_)
                                        + (spin_p)*(n_atoms_*n_orbs_*S_);

                                if(!NA_spinflip){
                                    if(col_OP>row_OP){
                                        index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                        OP_val = OPs_.value[index_OP];
                                    }
                                    else{
                                        index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                        OP_val = conj(OPs_.value[index_OP]);
                                    }

                                    E_class_temp -=Parameters_.U0[b_beta]*OP_val*OP_val_p;
                                }
                            }
                        }
                    }
                }
            }
        }
    }



    //Inter-orbital nn
    //Hartree
    for(int j=0;j<S_;j++){
        for(int atom_no=0;atom_no<n_atoms_;atom_no++){
            for(int alpha=0;alpha<n_orbs_;alpha++){
                for(int spin_=0;spin_<2;spin_++){
                    for(int beta=(alpha+1);beta<n_orbs_;beta++){
                        for(int spin_p=0;spin_p<2;spin_p++){
                            col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                    spin_*(n_atoms_*n_orbs_*S_);
                            row_=col_;
                            index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val_p = OPs_.value[index_OP_p];

                            row_OP = j + (atom_no + n_atoms_*beta)*(S_)
                                    + (spin_p)*(n_atoms_*n_orbs_*S_);
                            col_OP = row_OP;
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            E_class_temp += (Parameters_.UInterOrb[atom_no][alpha][beta])*OPs_.value[index_OP]*OP_val_p;

                            col_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                    spin_*(n_atoms_*n_orbs_*S_);
                            row_=col_;
                            index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val_p = OPs_.value[index_OP_p];

                            row_OP = j + (atom_no + n_atoms_*alpha)*(S_)
                                    + (spin_p)*(n_atoms_*n_orbs_*S_);
                            col_OP = row_OP;
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            E_class_temp += (Parameters_.UInterOrb[atom_no][alpha][beta])*OPs_.value[index_OP]*OP_val_p;

                        }
                    }
                }
            }}}



    //Inter-orbital nn
    //Fock
    for(int j=0;j<S_;j++){
        for(int atom_no=0;atom_no<n_atoms_;atom_no++){
            for(int alpha=0;alpha<n_orbs_;alpha++){
                for(int beta=(alpha+1);beta<n_orbs_;beta++){

                    for(int spin_=0;spin_<2;spin_++){

                        //1 and 7 [2 and 8 are h.c. of 1 and 7 ]
                        col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        row_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        if(col_>row_){
                            index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val_p = OPs_.value[index_OP_p];
                        }
                        else{
                            index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val_p = conj(OPs_.value[index_OP_p]);
                        }

                        col_OP = row_;
                        row_OP = col_;

                        if(col_OP>row_OP){
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = OPs_.value[index_OP];
                        }
                        else{
                            index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = conj(OPs_.value[index_OP]);
                        }
                        //1 and 7
                        E_class_temp -= (Parameters_.UInterOrb[atom_no][alpha][beta])*OP_val*OP_val_p;
                        //2 and 8
                        E_class_temp -= conj((Parameters_.UInterOrb[atom_no][alpha][beta])*OP_val*OP_val_p);


                        //3 and 5 [4 and 6 are h.c. of 3 and 5 ]
                        col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        row_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);

                        if(!NA_spinflip){
                            if(col_>row_){
                                index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val_p = OPs_.value[index_OP_p];
                            }
                            else{
                                index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val_p = conj(OPs_.value[index_OP_p]);
                            }

                            col_OP = row_;
                            row_OP = col_;

                            if(col_OP>row_OP){
                                index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val = OPs_.value[index_OP];
                            }
                            else{
                                index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val = conj(OPs_.value[index_OP]);
                            }
                            //3 and 5
                            E_class_temp -= (Parameters_.UInterOrb[atom_no][alpha][beta])*OP_val*OP_val_p;
                            //4 and 6
                            E_class_temp -= conj((Parameters_.UInterOrb[atom_no][alpha][beta])*OP_val*OP_val_p);
                        }
                    }

                }}}}





    //Hunds Coupling
    //Hartree
    for(int j=0;j<S_;j++){
        for(int atom_no=0;atom_no<n_atoms_;atom_no++){
            for(int alpha=0;alpha<n_orbs_;alpha++){
                for(int beta=(alpha+1);beta<n_orbs_;beta++){

                    for(int spin_=0;spin_<2;spin_++){

                        //1 and 3
                        col_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        row_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);

                        if(!NA_spinflip){
                            if(col_>row_){
                                index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val_p = OPs_.value[index_OP_p];
                            }
                            else{
                                index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val_p = conj(OPs_.value[index_OP_p]);
                            }


                            col_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                    (1-spin_)*(n_atoms_*n_orbs_*S_);
                            row_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                    (spin_)*(n_atoms_*n_orbs_*S_);

                            if(col_OP>row_OP){
                                index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val = OPs_.value[index_OP];
                            }
                            else{
                                index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val = conj(OPs_.value[index_OP]);
                            }
                            //1 and 3
                            E_class_temp -= (0.5*Parameters_.JHund[atom_no])*2.0*OP_val*OP_val_p;
                        }

                        //2 and 4
                        col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        row_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);

                        if(!NA_spinflip){
                            if(col_>row_){
                                index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val_p = OPs_.value[index_OP_p];
                            }
                            else{
                                index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val_p = conj(OPs_.value[index_OP_p]);
                            }

                            col_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                    (1-spin_)*(n_atoms_*n_orbs_*S_);
                            row_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                    (spin_)*(n_atoms_*n_orbs_*S_);

                            if(col_OP>row_OP){
                                index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val = OPs_.value[index_OP];
                            }
                            else{
                                index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val = conj(OPs_.value[index_OP]);
                            }
                            //2 and 4
                            E_class_temp -= (0.5*Parameters_.JHund[atom_no])*2.0*OP_val*OP_val_p;
                        }

                        //5 and 11
                        col_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        row_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);
                        if(col_>row_){
                            index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val_p = OPs_.value[index_OP_p];
                        }
                        else{
                            index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val_p = conj(OPs_.value[index_OP_p]);
                        }

                        col_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);
                        row_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);

                        if(col_OP>row_OP){
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = OPs_.value[index_OP];
                        }
                        else{
                            index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = conj(OPs_.value[index_OP]);
                        }
                        //5 and 11
                        E_class_temp -= (0.5*Parameters_.JHund[atom_no])*OP_val*OP_val_p;


                        //6 and 12
                        col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        row_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);
                        if(col_>row_){
                            index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val_p = OPs_.value[index_OP_p];
                        }
                        else{
                            index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val_p = conj(OPs_.value[index_OP_p]);
                        }

                        col_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);
                        row_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);

                        if(col_OP>row_OP){
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = OPs_.value[index_OP];
                        }
                        else{
                            index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = conj(OPs_.value[index_OP]);
                        }
                        //6 and 12
                        E_class_temp -= (0.5*Parameters_.JHund[atom_no])*OP_val*OP_val_p;



                        //7 and 9
                        col_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        row_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);
                        if(col_>row_){
                            index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val_p = OPs_.value[index_OP_p];
                        }
                        else{
                            index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val_p = conj(OPs_.value[index_OP_p]);
                        }

                        col_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);
                        row_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);

                        if(col_OP>row_OP){
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = OPs_.value[index_OP];
                        }
                        else{
                            index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = conj(OPs_.value[index_OP]);
                        }
                        //7 and 9
                        E_class_temp += (0.5*Parameters_.JHund[atom_no])*OP_val*OP_val_p;


                        //8 and 10
                        col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        row_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);
                        if(col_>row_){
                            index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val_p = OPs_.value[index_OP_p];
                        }
                        else{
                            index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val_p = conj(OPs_.value[index_OP_p]);
                        }

                        col_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);
                        row_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);

                        if(col_OP>row_OP){
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = OPs_.value[index_OP];
                        }
                        else{
                            index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = conj(OPs_.value[index_OP]);
                        }
                        //8 and 10
                        E_class_temp += (0.5*Parameters_.JHund[atom_no])*OP_val*OP_val_p;

                    }

                }}}}





    //Hunds Coupling
    //Fock
    for(int j=0;j<S_;j++){
        for(int atom_no=0;atom_no<n_atoms_;atom_no++){
            for(int alpha=0;alpha<n_orbs_;alpha++){
                for(int beta=(alpha+1);beta<n_orbs_;beta++){

                    for(int spin_=0;spin_<2;spin_++){

                        //1 and 3 [hc : 2 and 4]
                        col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        row_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);
                        if(col_>row_){
                            index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val_p = OPs_.value[index_OP_p];
                        }
                        else{
                            index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val_p = conj(OPs_.value[index_OP_p]);
                        }

                        col_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);
                        row_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);

                        if(col_OP>row_OP){
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = OPs_.value[index_OP];
                        }
                        else{
                            index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = conj(OPs_.value[index_OP]);
                        }
                        //1 and 3
                        E_class_temp += (0.5*Parameters_.JHund[atom_no])*2.0*OP_val*OP_val_p;
                        //2 and 4
                        E_class_temp += conj((0.5*Parameters_.JHund[atom_no])*2.0*OP_val*OP_val_p);



                        //5 and 11 [hc : 6 and 12]
                        col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        row_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);
                        if(col_>row_){
                            index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val_p = OPs_.value[index_OP_p];
                        }
                        else{
                            index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val_p = conj(OPs_.value[index_OP_p]);
                        }

                        col_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);
                        row_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);

                        if(col_OP>row_OP){
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = OPs_.value[index_OP];
                        }
                        else{
                            index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = conj(OPs_.value[index_OP]);
                        }
                        //5 and 11
                        E_class_temp += (0.5*Parameters_.JHund[atom_no])*OP_val*OP_val_p;
                        //6 and 12
                        E_class_temp += conj((0.5*Parameters_.JHund[atom_no])*OP_val*OP_val_p);


                        //7 and 9 [hc : 8 and 10]
                        col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        row_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);
                        if(!NA_spinflip){
                            if(col_>row_){
                                index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val_p = OPs_.value[index_OP_p];
                            }
                            else{
                                index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val_p = conj(OPs_.value[index_OP_p]);
                            }

                            col_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                    (1-spin_)*(n_atoms_*n_orbs_*S_);
                            row_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                    (spin_)*(n_atoms_*n_orbs_*S_);

                            if(col_OP>row_OP){
                                index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val = OPs_.value[index_OP];
                            }
                            else{
                                index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val = conj(OPs_.value[index_OP]);
                            }
                            //7 and 9
                            E_class_temp -= (0.5*Parameters_.JHund[atom_no])*OP_val*OP_val_p;
                            //8 and 10
                            E_class_temp -= conj((0.5*Parameters_.JHund[atom_no])*OP_val*OP_val_p);
                        }

                    }

                }}}}





    //Pair Hopping
    //Hartree+Fock
    for(int j=0;j<S_;j++){
        for(int atom_no=0;atom_no<n_atoms_;atom_no++){
            for(int alpha=0;alpha<n_orbs_;alpha++){
                for(int beta=(alpha+1);beta<n_orbs_;beta++){

                    for(int spin_=0;spin_<2;spin_++){

                        //1 and 2 [hc : 3 and 4]
                        col_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        row_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);
                        if(!NA_spinflip){
                            if(col_>row_){
                                index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val_p = OPs_.value[index_OP_p];
                            }
                            else{
                                index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val_p = conj(OPs_.value[index_OP_p]);
                            }

                            col_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                    (1-spin_)*(n_atoms_*n_orbs_*S_);
                            row_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                    (spin_)*(n_atoms_*n_orbs_*S_);

                            if(col_OP>row_OP){
                                index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val = OPs_.value[index_OP];
                            }
                            else{
                                index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                OP_val = conj(OPs_.value[index_OP]);
                            }
                            //1 and 2
                            E_class_temp -= (1.0*Parameters_.JHund[atom_no])*OP_val*OP_val_p;
                            //3 and 4
                            E_class_temp -= conj((1.0*Parameters_.JHund[atom_no])*OP_val*OP_val_p);
                        }

                        //5 and 6 [hc : 7 and 8]
                        col_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                spin_*(n_atoms_*n_orbs_*S_);
                        row_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (spin_)*(n_atoms_*n_orbs_*S_);
                        if(col_>row_){
                            index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val_p = OPs_.value[index_OP_p];
                        }
                        else{
                            index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val_p = conj(OPs_.value[index_OP_p]);
                        }

                        col_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);
                        row_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                (1-spin_)*(n_atoms_*n_orbs_*S_);

                        if(col_OP>row_OP){
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = OPs_.value[index_OP];
                        }
                        else{
                            index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                            OP_val = conj(OPs_.value[index_OP]);
                        }
                        //5 and 6
                        E_class_temp += (1.0*Parameters_.JHund[atom_no])*OP_val*OP_val_p;
                        //7 and 8
                        E_class_temp += conj((1.0*Parameters_.JHund[atom_no])*OP_val*OP_val_p);

                    }
                }}}}




    E_class = -0.5*ncells_*E_class_temp.real();


    //----------Eclass done


    E_quant=0.0;
    for(int n=0;n<Eigenvalues_.size();n++){
        E_quant += (Eigenvalues_[n]*
                    (1.0/( exp((Eigenvalues_[n]-mu_)*Parameters_.beta ) + 1.0))
                    );
    }

}

void Kspace_calculation_G2dLatticeNew::Get_Energies(){}

void Kspace_calculation_G2dLatticeNew::Get_new_OPs_and_error(){

    //For <c_{c1}* c_{c2}>
    int c1;
    int c2;
    int S_=NSites_in_MUC;
    int cell_row, alpha_row, gamma_row, sigma_row, cell_col, gamma_col, alpha_col, sigma_col;
    int alpha_plus_gamma, alpha_plus_gamma_sigma;
    int row_temp, col_temp;

    int k_index;
    int d1_, d2_;
    int state;

    for(int OP_no=0;OP_no<OPs_.value.size();OP_no++){
        row_temp=OPs_.rows[OP_no];
        col_temp=OPs_.columns[OP_no];


        //alpha + gamma*(S_) +  sigma*(n_atoms_*n_orbs_*S_) + cell_*(2*n_atoms_*n_orbs_*S_)
        cell_row = 0;
        alpha_plus_gamma = row_temp%(n_atoms_*n_orbs_*S_);
        alpha_row = alpha_plus_gamma%S_;
        gamma_row = (alpha_plus_gamma -alpha_row)/S_;
        sigma_row = (row_temp - alpha_plus_gamma)/(n_atoms_*n_orbs_*S_);
        c1 = alpha_row + gamma_row*(S_) +  sigma_row*(n_atoms_*n_orbs_*S_);
        assert(c1==row_temp);

        alpha_plus_gamma_sigma = col_temp%(2*n_atoms_*n_orbs_*S_);
        alpha_plus_gamma = alpha_plus_gamma_sigma%(n_atoms_*n_orbs_*S_);
        alpha_col = alpha_plus_gamma%S_;
        gamma_col = (alpha_plus_gamma -alpha_col)/S_;
        sigma_col = (alpha_plus_gamma_sigma - alpha_plus_gamma)/(n_atoms_*n_orbs_*S_);
        cell_col = (col_temp - alpha_plus_gamma_sigma)/(2*n_atoms_*n_orbs_*S_);
        Get_minimum_distance_direction(0, cell_col, d1_, d2_);
        c2 = alpha_col + gamma_col*(S_) +  sigma_col*(n_atoms_*n_orbs_*S_);

        OPs_new_.value[OP_no]=0.0;
        for(int n=0;n<2*n_atoms_*n_orbs_*S_;n++){ //band_index
            for(int k1=0;k1<lx_cells;k1++){
                for(int k2=0;k2<ly_cells;k2++){
                    k_index = Coordinates_.Ncell(k1,k2);
                    state = 2*n_atoms_*n_orbs_*S_*k_index + n;

                    // cout<<(1.0/ncells_)*(
                    //                               (conj(Eigvectors_[state][c1])*Eigvectors_[state][c2])*
                    //                               (exp(-1.0*iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_   )  ))
                    //                               *(1.0/( exp((Eigenvalues_saved[state]-mu_)*Parameters_.beta ) + 1.0))
                    //                            )<<endl;
                    OPs_new_.value[OP_no] += (1.0/ncells_)*(
                                (conj(Eigvectors_[state][c1])*Eigvectors_[state][c2])*
                                (exp(-1.0*iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_   )  ))
                                *(1.0/( exp((Eigenvalues_saved[state]-mu_)*Parameters_.beta ) + 1.0))
                                );
                }
            }
        }


    }




    OP_error_=0.0;
    for(int OP_no=0;OP_no<OPs_.value.size();OP_no++){
        OP_error_ += abs((OPs_.value[OP_no] - OPs_new_.value[OP_no])*conj(OPs_.value[OP_no] - OPs_new_.value[OP_no]));
    }
    OP_error_ = sqrt(OP_error_);


}


void Kspace_calculation_G2dLatticeNew::Get_spin_resolved_local_densities(){

    //For <c_{c1}* c_{c2}>
    int c1;
    int c2;
    int S_=NSites_in_MUC;

    int k_index;
    int state;

    double Total_den_up, Total_den_dn;

    char temp_char[50];
    sprintf(temp_char, "%.10f", Parameters_.Temperature);

    string File_Out_Local_orb_densities = "Local_spin_resolved_densities" + string(temp_char)+ ".txt";
    ofstream file_out_Local_orb_densities(File_Out_Local_orb_densities.c_str());
    file_out_Local_orb_densities<<"#alpha(site in real lattice)   0up   0dn  1up     1dn"<<endl;

    complex<double> val;

    Total_den_up=0.0;
    Total_den_dn=0.0;
    int site_x, site_y,site;
    for(int alpha=0;alpha<NSites_in_MUC;alpha++){

        for(int cell_1=0;cell_1<lx_cells;cell_1++){

            // site_x = (cell_1*) + alpha_1;
            for(int cell_2=0;cell_2<ly_cells;cell_2++){
                //       site_y = (cell_2*UnitCellSize_y) + alpha_2;

                //     alpha = alpha_1 + alpha_2*(UnitCellSize_x);
                //    site = site_x + site_y*(lx_);
                int alpha_1 = Intra_MUC_positions[alpha].first;
                int alpha_2 = Intra_MUC_positions[alpha].second;

                int i1_ = (((cell_1)*Mat_MUC(0,0) + (cell_2)*Mat_MUC(1,0)) + alpha_1 + lx_*ly_)%lx_;
                int i2_ = (((cell_1)*Mat_MUC(0,1) + (cell_2)*Mat_MUC(1,1)) + alpha_2 + lx_*ly_)%ly_;
                site = i1_ + i2_*lx_;

                file_out_Local_orb_densities<<site<<"    ";

                for(int gamma=0;gamma<n_orbs_*n_atoms_;gamma++){
                    for(int sigma=0;sigma<2;sigma++){

                        c1 = alpha + gamma*S_ + sigma*(n_orbs_*n_atoms_*S_);
                        c2=c1;
                        val=0.0;
                        for(int n=0;n<2*n_orbs_*n_atoms_*S_;n++){ //band_index
                            for(int k1=0;k1<lx_cells;k1++){
                                for(int k2=0;k2<ly_cells;k2++){
                                    k_index = Coordinates_.Ncell(k1,k2);
                                    state = 2*n_atoms_*n_orbs_*S_*k_index + n;

                                    val += (1.0/ncells_)*(
                                                (conj(Eigvectors_[state][c1])*Eigvectors_[state][c2])
                                                *(1.0/( exp((Eigenvalues_saved[state]-mu_)*Parameters_.beta ) + 1.0))
                                                );
                                }
                            }
                        }

                        file_out_Local_orb_densities<<val.real()<<"    ";

                        if(sigma==0){
                            Total_den_up +=val.real();
                        }
                        else{
                            Total_den_dn +=val.real();
                        }

                    }
                }
                file_out_Local_orb_densities<<endl;
            }

        }
    }



    cout<< "Total electrons in Lattice, UP = "<<Total_den_up<<endl;
    cout<< "Total electrons in Lattice, DOWN = "<<Total_den_dn<<endl;


}




void Kspace_calculation_G2dLatticeNew::Update_Total_Density(){}

void Kspace_calculation_G2dLatticeNew::Get_Tau_Pseudospins(){}



void Kspace_calculation_G2dLatticeNew::Create_InverseMapping(){

    cell_no.resize(lx_);
    magnetic_atom_no.resize(lx_);
    for(int ix=0;ix<lx_;ix++){
        cell_no[ix].resize(ly_);
        magnetic_atom_no[ix].resize(ly_);
    }


    pair_int cell_temp, alpha_temp;
    for(int alpha=0;alpha<NSites_in_MUC;alpha++){
        for(int cell_1=0;cell_1<lx_cells;cell_1++){
            for(int cell_2=0;cell_2<ly_cells;cell_2++){

                int alpha_1 = Intra_MUC_positions[alpha].first;
                int alpha_2 = Intra_MUC_positions[alpha].second;

                int i1_ = (((cell_1)*Mat_MUC(0,0) + (cell_2)*Mat_MUC(1,0)) + alpha_1 + lx_*ly_)%lx_;
                int i2_ = (((cell_1)*Mat_MUC(0,1) + (cell_2)*Mat_MUC(1,1)) + alpha_2 + lx_*ly_)%ly_;

                cell_temp.first = cell_1;
                cell_temp.second = cell_2;

                cell_no[i1_][i2_] = cell_temp;
                magnetic_atom_no[i1_][i2_] = alpha;
            }}}

}


void Kspace_calculation_G2dLatticeNew::Get_local_spins(){

    //For <c_{c1}* c_{c2}>


    double a_NN;
    a_NN=1.0;
    int c1;
    int c2;
    int S_=NSites_in_MUC;
    string atomorb_name;

    int k_index;
    int state;

    double Total_Sz, Total_Sx, Total_Sy, Total_den;
    int UP_, DOWN_;
    UP_=0;DOWN_=1;

    char temp_char[50];
    sprintf(temp_char, "%.10f", Parameters_.Temperature);

    string File_Out_Local_orb_densities = "Local_spins"  + string(temp_char) +  ".txt";
    ofstream file_out_Local_orb_densities(File_Out_Local_orb_densities.c_str());
    file_out_Local_orb_densities<<"#site (in real lattice) site_1 site_2  rx  ry   Sz   Sx  Sy  |S|  <n>"<<endl;

    complex<double> val;

    Total_Sz=0.0;
    Total_Sx=0.0;
    Total_Sy=0.0;
    Total_den=0.0;
    double rx_, ry_;
    complex<double> splus_val;
    double sz_val, sx_val, sy_val, den_val;
    int site_x, site_y, alpha, site;



    //site_x = (cell_1*UnitCellSize_x) + alpha_1;

    for(int site_x=0;site_x<lx_;site_x++){
        for(int site_y=0;site_y<ly_;site_y++){
            alpha = magnetic_atom_no[site_x][site_y];
            // int alpha_1 = Intra_MUC_positions[alpha].first;
            // int alpha_2 =  Intra_MUC_positions[alpha].second;
            // int cell_1 = cell_no[site_x][site_y].first;
            // int cell_2 = cell_no[site_x][site_y].second;

            for(int gamma=0;gamma<n_orbs_*n_atoms_;gamma++){

                int orb_, subltc_;
                subltc_= gamma%n_atoms_;
                orb_ = int((gamma - subltc_+0.5)/n_atoms_);

                atomorb_name="subltc_"+to_string(subltc_)+"_orb_"+to_string(orb_);

                // rx_ = a_NN*((sqrt(3.0))*(site_x) +  ((sqrt(3.0))/2.0)*(site_y)  + ((sqrt(3.0))/2.0)*(subltc_) );
                // ry_ = (0.0*(site_x) + (3.0/2.0)*(site_y)   +   (1.0/2.0)*subltc_)*a_NN;

                rx_ = 2*a_NN*(site_x);
                ry_ = 2*a_NN*(site_y);
                if(subltc_==1){
                    rx_ += a_NN;
                }
                if(subltc_==2){
                    ry_ +=a_NN;
                }

                //alpha = alpha_1 + alpha_2*(UnitCellSize_x);
                site = site_x + site_y*(lx_);

                file_out_Local_orb_densities<<site<<"    "<<site_x<<"    "<<site_y<<"    "<<site<<atomorb_name<<"   "<<rx_<<"     "<<ry_<<"     ";


                //Local Splus
                val=0.0;
                c1 = alpha + gamma*(S_) + UP_*(n_atoms_*n_orbs_*S_);
                c2 = alpha + gamma*(S_)  + DOWN_*(n_atoms_*n_orbs_*S_);
                for(int n=0;n<2*n_atoms_*n_orbs_*S_;n++){ //band_index
                    for(int k1=0;k1<lx_cells;k1++){
                        for(int k2=0;k2<ly_cells;k2++){
                            k_index = Coordinates_.Ncell(k1,k2);
                            state = 2*n_atoms_*n_orbs_*S_*k_index + n;
                            val += (1.0/ncells_)*(
                                        (conj(Eigvectors_[state][c1])*Eigvectors_[state][c2])
                                        *(1.0/( exp((Eigenvalues_saved[state]-mu_)*Parameters_.beta ) + 1.0))
                                        );
                        }
                    }
                }
                splus_val = val;
                sx_val = splus_val.real();
                sy_val = splus_val.imag();

                //Local Sz //HERE
                val=0.0;
                double fac_;
                for(int sigma=0;sigma<2;sigma++){
                    fac_=1.0 - (2.0*sigma);
                    c1 = alpha + gamma*(S_) + sigma*(n_orbs_*n_atoms_*S_);
                    c2=c1;
                    for(int n=0;n<2*n_atoms_*n_orbs_*S_;n++){ //band_index
                        for(int k1=0;k1<lx_cells;k1++){
                            for(int k2=0;k2<ly_cells;k2++){
                                k_index = Coordinates_.Ncell(k1,k2);
                                state = 2*n_atoms_*n_orbs_*S_*k_index + n;
                                val += 0.5*(fac_/ncells_)*(
                                            (conj(Eigvectors_[state][c1])*Eigvectors_[state][c2])
                                            *(1.0/( exp((Eigenvalues_saved[state]-mu_)*Parameters_.beta ) + 1.0))
                                            );
                            }
                        }
                    }
                }
                sz_val=val.real();



                //Local den //HERE
                val=0.0;
                for(int sigma=0;sigma<2;sigma++){
                    fac_=1.0;
                    c1 = alpha + gamma*(S_) + sigma*(n_atoms_*n_orbs_*S_);
                    c2=c1;
                    for(int n=0;n<2*n_atoms_*n_orbs_*S_;n++){ //band_index
                        for(int k1=0;k1<lx_cells;k1++){
                            for(int k2=0;k2<ly_cells;k2++){
                                k_index = Coordinates_.Ncell(k1,k2);
                                state = 2*n_atoms_*n_orbs_*S_*k_index + n;
                                val += (fac_/ncells_)*(
                                            (conj(Eigvectors_[state][c1])*Eigvectors_[state][c2])
                                            *(1.0/( exp((Eigenvalues_saved[state]-mu_)*Parameters_.beta ) + 1.0))
                                            );
                            }
                        }
                    }
                }
                den_val=val.real();


                file_out_Local_orb_densities<<sz_val<<"    "<<sx_val<<"   "<<sy_val<<"    "<<sqrt(sz_val*sz_val + sy_val*sy_val + sx_val*sx_val)<<"   "<<den_val<<"     ";


                Total_Sz +=sz_val;
                Total_Sy +=sy_val;
                Total_Sx +=sx_val;



                file_out_Local_orb_densities<<endl;
            }
        }
        file_out_Local_orb_densities<<endl;

    }



    cout<< "Total Sz = "<<Total_Sz<<endl;
    cout<< "Total Sx = "<<Total_Sx<<endl;
    cout<< "Total Sy = "<<Total_Sy<<endl;
    cout<< "|Total_S| = "<<sqrt(Total_Sz*Total_Sz + Total_Sx*Total_Sx + Total_Sy*Total_Sy)<<endl;



}


double Kspace_calculation_G2dLatticeNew::chemicalpotential(double Particles){


    double mu_out;
    double n1,N;
    double dMubydN;
    double muin;
    N=Particles;
    int N_ = int(N);
    muin = 0.5*(Eigenvalues_[N_-1] - Eigenvalues_[N_]);
    double nstate = Eigenvalues_.size();
    dMubydN = (0.0005*(Eigenvalues_[nstate-1] - Eigenvalues_[0])/nstate);

    //temp=Parameters_.temp;
    mu_out = muin;
    bool converged=false;
    int final_i;


    if(1==2){
        for(int i=0;i<50000;i++){
            n1=0.0;
            for(int j=0;j<nstate;j++){
                n1+=double(1.0/( exp( (Eigenvalues_[j]-mu_out)*Parameters_.beta ) + 1.0));
            }
            //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
            if(abs(N-n1)<double(0.0001)){
                //cout<<abs(N-n1)<<endl;
                converged=true;
                final_i=i;
                break;
            }
            else {
                mu_out += (N-n1)*dMubydN;
                cout<<i<<"    "<<n1<<"    "<<mu_out<<endl;

            }
        }

        if(!converged){
            cout<<"mu_not_converged, N = "<<n1<<endl;
        }
        else{
            //cout<<"mu converged, N = "<<n1<<" in "<<final_i<<" iters"<<endl;
        }

    }


    double mu1, mu2;
    double mu_temp = muin;
    //cout<<"mu_input = "<<mu_temp<<endl;
    if(1==1){
        mu1=Eigenvalues_[0]- (10000*(1.0/Parameters_.beta));
        mu2=Eigenvalues_[nstate-1]+ (10000*(1.0/Parameters_.beta));
        for(int i=0;i<40000;i++){
            n1=0.0;
            for(int j=0;j<nstate;j++){
                n1+=double(1.0/( exp( (Eigenvalues_[j]-mu_temp)*Parameters_.beta ) + 1.0));
            }
            //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_temp<< endl;
            if(abs(N-n1)<double(0.0001)){
                //cout<<abs(N-n1)<<endl;
                converged=true;
                break;
            }
            else {
                if(n1 >N){
                    mu2=mu_temp;
                    mu_temp=0.5*(mu1 + mu_temp);
                }
                else{
                    mu1=mu_temp;
                    mu_temp=0.5*(mu2 + mu_temp);
                }

            }
            //cout<<"mu_temp = "<<mu_temp<<"   "<<n1<<endl;
        }

        if(!converged){
            cout<<"mu_not_converged, N,mu = "<<n1<<", "<<mu_temp<<endl;
        }
        else{
            //cout<<"mu converged, N = "<<n1<<endl;
        }

        mu_out = mu_temp;
    }



    return mu_out;
} // ----------

double Kspace_calculation_G2dLatticeNew::random1(){

    return dis1_(Generator1_);

}



void Kspace_calculation_G2dLatticeNew::Initialize()
{
    kick_while_cooling=0.0;
    
    NOT_AVAIL_INT=-1000;

    OP_only_finite_Int=true;
    Global_Eps=0.0000001;

    cout<<"Parameters_.FockType = '" << Parameters_.FockType<<"'"<<endl;
    ly_ = Parameters_.ly;
    lx_ = Parameters_.lx;

    Mat_MUC = Parameters_.MUC_Mat;
    lx_cells=Connections_.lx_cells;
    ly_cells=Connections_.ly_cells;
    ncells_ = Connections_.ncells_;
    NSites_in_MUC = Connections_.NSites_in_MUC;
    Intra_MUC_positions = Connections_.Intra_MUC_positions;

    //Connections_bw_MUC();
    Create_InverseMapping();

    //Get_MUC_details();
    assert(ncells_ == lx_cells * ly_cells);

    n_orbs_ = Parameters_.n_orbs;
    n_atoms_ = Parameters_.n_atoms;


    //assert(n_orbs_==1);
    //assert(n_atoms_==3);

    Ham_.resize(2*n_orbs_*n_atoms_*NSites_in_MUC, 2*n_orbs_*n_atoms_*NSites_in_MUC);
    Eigvectors_.resize(2*n_orbs_*n_atoms_*lx_*ly_);
    Eigenvalues_.resize(2*n_orbs_*n_atoms_*lx_*ly_);
    for(int i=0;i<2*n_atoms_*n_orbs_*lx_*ly_;i++){
        Eigvectors_[i].resize(2*n_orbs_*n_atoms_*NSites_in_MUC); //Eigenvector number
    }


    Eigenvalues_spin_down.resize(n_orbs_*n_atoms_*lx_*ly_);
    Eigenvalues_spin_up.resize(n_orbs_*n_atoms_*lx_*ly_);

    Kx_values.resize(ncells_);
    Ky_values.resize(ncells_);


    // Create_V_int();
    // Create_B_and_C_mat();


    //Order Parameters:
    /*
    Inter-unit cell OP's ( intra-unit cell when d_vec=0)
    <c_{alpha, gamma, sigma}^{dagger} c _{ d_vec, alpha', gamma', sigma' }> ; for all d_vec,
    and alpha+ gamma*(S_) + sigma*(n_orbs_*S_)
       <= alpha'+ gamma'*(S_) + sigma'*(n_orbs_*S_)

    For <c_{alpha, gamma, sigma}^{dagger} c _{ d_vec, alpha',gamma', sigma' }> ,
   when alpha+ gamma*(S_) + sigma*(n_orbs_*S_)
       > alpha'+ gamma'*(S_) + sigma'*(n_orbs_*S_)
    use
    <c_{alpha, gamma, sigma}^{dagger} c _{ d_vec, alpha',gamma', sigma' }>
    = conj(<c_{alpha', gamma', sigma'}^{dagger} c _{ -d_vec, alpha,gamma, sigma}> )

    P=2*n_orbs_*S_

    Total no of Order parameters (Hartree + all Fock):
    ncells_*(P(P+1))/2

    Total no of Order parameters (Hartree + only intracell Fock):
    (P(P+1))/2

    Total no of Order parameters (Hartree):
    P
    */

    int S_=NSites_in_MUC; //no. of sites in single unit cell
    int r1, r2;
    OPs_.value.clear();
    OPs_.rows.clear();
    OPs_.columns.clear();
    OPs_new_.value.clear();
    OPs_new_.rows.clear();
    OPs_new_.columns.clear();
    SI_to_ind.resize(2*n_orbs_*n_atoms_*S_*ncells_*2*n_orbs_*n_atoms_*S_);

    for(int i=0;i<SI_to_ind.size();i++){
        SI_to_ind[i]=NOT_AVAIL_INT; //Not available unless given by following routine.
    }

    int row_temp, col_temp; // col_ = alpha + gamma*(S_) + sigma*(n_atoms_*n_orbs_*S_) + cell_no*(2*n_atoms_*n_orbs_*S_) ;

    double temp_den;


    //Parameters_.Create_OPs_Ansatz=false;
    //Parameters_.Create_OPs_Ansatz=true;
    //Parameters_.OP_Ansatz_type="TL_WC_120AFM";
    //Parameters_.OP_Ansatz_type="AFM_Wigner";
    //Parameters_.OP_Ansatz_type="NM_Wigner";
    //Parameters_.OP_Ansatz_type="Tetrahedron_n3by4";
    //Parameters_.OP_Ansatz_type="Tri_Prism_Wigner";
    //Parameters_.OP_Ansatz_type="Vertices_12_NonCopl_State_Wigner";
    //Parameters_.OP_Ansatz_type = "Vertices_4_NonCopl_Tetrahedron_Wigner";
    //Parameters_.OP_Ansatz_type ="Vertices_10_NonCopl_State_Wigner_n2by3";


    if(!Parameters_.Read_OPs){

        if(!Parameters_.Create_OPs_Ansatz){
            //Hartree Terms
            for(int alpha=0;alpha<NSites_in_MUC;alpha++){
                for(int gamma=0;gamma<n_atoms_*n_orbs_;gamma++){
                    for(int sigma=0;sigma<2;sigma++){
                        // OPs_.value.push_back(complex<double> (random1(),0.0));
                        //OPs_.value.push_back(complex<double> (0.2,0.0));
                        //OPs_new_.value.push_back(0.0);

                        row_temp=  alpha + gamma*(S_) +  sigma*(n_atoms_*n_orbs_*S_) + 0*(2*n_atoms_*n_orbs_*S_);
                        col_temp=row_temp;

                        OPs_.rows.push_back(row_temp);OPs_new_.rows.push_back(row_temp);
                        OPs_.columns.push_back(col_temp);OPs_new_.columns.push_back(col_temp);

                        OPs_new_.value.push_back(0.0);
                        if(Parameters_.Fixing_mu){
                            //assuming near half-filling
                            temp_den = 0.5;//half-filling
                            temp_den += (random1()-0.5)*0.0;
                            OPs_.value.push_back(complex<double> (temp_den,0.0));
                        }
                        else{
                            temp_den = random1();
                            OPs_.value.push_back(complex<double> (temp_den,0.0));
                        }

                        SI_to_ind[col_temp + row_temp*(ncells_*2*n_atoms_*n_orbs_*S_)] = OPs_.value.size()-1;
                    }
                }
            }

            //Fock
            bool check_;
            bool NA_spinflip;
            int alpha_1, alpha_2, alpha_p_1, alpha_p_2, d1_org, d2_org;
            int row_,col_;
            int atom_orb_no, atom_orb_no_p;
            if(!Parameters_.Just_Hartree){
                for(int alpha=0;alpha<NSites_in_MUC;alpha++){
                    // alpha_1 = alpha % UnitCellSize_x;
                    // alpha_2 = (alpha - alpha_1)/UnitCellSize_x;

                    for(int atom_no=0;atom_no<n_atoms_;atom_no++){
                        for(int orb_no=0;orb_no<n_orbs_;orb_no++){
                            atom_orb_no = atom_no + n_atoms_*orb_no;
                            for(int sigma=0;sigma<2;sigma++){

                                //for(int atom_no_p=0;atom_no_p<n_atoms_;atom_no_p++){
                                int atom_no_p=atom_no;
                                for(int orb_no_p=0;orb_no_p<n_orbs_;orb_no_p++){
                                    atom_orb_no_p = atom_no_p + n_atoms_*orb_no_p;
                                    for(int sigma_p=0;sigma_p<2;sigma_p++){


                                        //check_= (sigma_p > sigma);
                                        row_temp = alpha + (atom_orb_no)*(S_) +  sigma*(n_atoms_*n_orbs_*S_) + 0*(2*n_atoms_*n_orbs_*S_);
                                        col_temp = alpha + (atom_orb_no_p)*(S_) +  sigma_p*(n_atoms_*n_orbs_*S_) + 0*(2*n_atoms_*n_orbs_*S_);
                                        check_ = (col_temp > row_temp);

                                        NA_spinflip = Parameters_.NoSpinFlipOP && (sigma_p!=sigma);

                                        if( check_ && (!NA_spinflip)){

                                            OPs_.value.push_back(complex<double> (random1(),random1()));
                                            OPs_new_.value.push_back(0.0);

                                            OPs_.rows.push_back(row_temp);OPs_new_.rows.push_back(row_temp);
                                            OPs_.columns.push_back(col_temp);OPs_new_.columns.push_back(col_temp);

                                            SI_to_ind[col_temp + row_temp*(ncells_*2*n_atoms_*n_orbs_*S_)] = OPs_.value.size()-1;
                                        }
                                    }
                                }
                                // }
                            }
                        }
                    }
                }
            }

            
        }
        else{
            //cout<<"Ansatz State not working right now"<<endl;
            //assert(false);

            if(Parameters_.OP_Ansatz_type=="AFM_AFO_3orb"){
                assert(n_orbs_==3);
                assert(n_atoms_==1);

                int XZ_=0;
                int YZ_=1;
                int XY_=2;
                int SPINUP_=0;
                int SPINDN_=1;

                //Hartree Terms
                for(int alpha=0;alpha<NSites_in_MUC;alpha++){
                    for(int gamma=0;gamma<n_orbs_;gamma++){
                        for(int sigma=0;sigma<2;sigma++){


                            row_temp=  alpha + gamma*(S_) +  sigma*(n_atoms_*n_orbs_*S_) + 0*(2*n_atoms_*n_orbs_*S_);
                            col_temp=row_temp;

                            OPs_.rows.push_back(row_temp);OPs_new_.rows.push_back(row_temp);
                            OPs_.columns.push_back(col_temp);OPs_new_.columns.push_back(col_temp);
                            OPs_new_.value.push_back(0.0);

                            int alpha_1 = Intra_MUC_positions[alpha].first;
                            int alpha_2 = Intra_MUC_positions[alpha].second;
                            int site_x_val =(((0)*Mat_MUC(0,0) + (0)*Mat_MUC(1,0)) + alpha_1 + lx_*ly_)%lx_;
                            int site_y_val =(((0)*Mat_MUC(0,1) + (0)*Mat_MUC(1,1)) + alpha_2 + lx_*ly_)%ly_;
                            int sublattice_type = pow(-1,site_x_val)*pow(-1,site_y_val);

                            if(sublattice_type==1){


                                if(gamma==XZ_ && sigma==SPINUP_){
                                    temp_den=1;
                                }
                                if(gamma==YZ_ && sigma==SPINUP_){
                                    temp_den=0;
                                }
                                if(gamma==XY_ && sigma==SPINUP_){
                                    temp_den=1;
                                }
                                if(sigma==SPINDN_){
                                    temp_den=0;
                                }

                            }
                            else{
                                assert(sublattice_type==-1);
                                if(gamma==XZ_ && sigma==SPINDN_){
                                    temp_den=0;
                                }
                                if(gamma==YZ_ && sigma==SPINDN_){
                                    temp_den=1;
                                }
                                if(gamma==XY_ && sigma==SPINDN_){
                                    temp_den=1;
                                }
                                if(sigma==SPINUP_){
                                    temp_den=0;
                                }

                            }

                            OPs_.value.push_back(complex<double> (temp_den,0.0));
                            SI_to_ind[col_temp + row_temp*(ncells_*2*n_atoms_*n_orbs_*S_)] = OPs_.value.size()-1;
                        }
                    }
                }




                //Fock
                bool check_;
                bool NA_spinflip;
                int atom_orb_no, atom_orb_no_p;
                if(!Parameters_.Just_Hartree){
                    for(int alpha=0;alpha<NSites_in_MUC;alpha++){
                        for(int atom_no=0;atom_no<n_atoms_;atom_no++){
                            for(int orb_no=0;orb_no<n_orbs_;orb_no++){
                                atom_orb_no = atom_no + n_atoms_*orb_no;
                                for(int sigma=0;sigma<2;sigma++){

                                    int atom_no_p=atom_no;
                                    for(int orb_no_p=0;orb_no_p<n_orbs_;orb_no_p++){
                                        atom_orb_no_p = atom_no_p + n_atoms_*orb_no_p;
                                        for(int sigma_p=0;sigma_p<2;sigma_p++){

                                            row_temp = alpha + (atom_orb_no)*(S_) +  sigma*(n_atoms_*n_orbs_*S_) + 0*(2*n_atoms_*n_orbs_*S_);
                                            col_temp = alpha + (atom_orb_no_p)*(S_) +  sigma_p*(n_atoms_*n_orbs_*S_) + 0*(2*n_atoms_*n_orbs_*S_);
                                            check_ = (col_temp > row_temp);

                                            NA_spinflip = Parameters_.NoSpinFlipOP && (sigma_p!=sigma);

                                            if( check_ && (!NA_spinflip)){

                                                OPs_.value.push_back(0.0);
                                                OPs_new_.value.push_back(0.0);

                                                OPs_.rows.push_back(row_temp);OPs_new_.rows.push_back(row_temp);
                                                OPs_.columns.push_back(col_temp);OPs_new_.columns.push_back(col_temp);

                                                SI_to_ind[col_temp + row_temp*(ncells_*2*n_atoms_*n_orbs_*S_)] = OPs_.value.size()-1;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }




            }
            else{
                cout<<"ANSATZ NOT ADDED YET"<<endl;
                assert(false);
            }


        }

    }

    else{


        //Hartree Terms
        for(int alpha=0;alpha<NSites_in_MUC;alpha++){
            for(int gamma=0;gamma<n_atoms_*n_orbs_;gamma++){
                for(int sigma=0;sigma<2;sigma++){
                    // OPs_.value.push_back(complex<double> (random1(),0.0));
                    //OPs_.value.push_back(complex<double> (0.2,0.0));
                    //OPs_new_.value.push_back(0.0);

                    row_temp=  alpha + gamma*(S_) +  sigma*(n_atoms_*n_orbs_*S_) + 0*(2*n_atoms_*n_orbs_*S_);
                    col_temp=row_temp;

                    OPs_.rows.push_back(row_temp);OPs_new_.rows.push_back(row_temp);
                    OPs_.columns.push_back(col_temp);OPs_new_.columns.push_back(col_temp);

                    OPs_new_.value.push_back(0.0);
                    if(Parameters_.Fixing_mu){
                        OPs_.value.push_back(complex<double> (0.0,0.0));
                    }
                    else{
                        OPs_.value.push_back(complex<double> (0.0,0.0));
                    }

                    SI_to_ind[col_temp + row_temp*(ncells_*2*n_atoms_*n_orbs_*S_)] = OPs_.value.size()-1;
                }
            }
        }

        //Fock
        bool check_;
        bool NA_spinflip;
        int alpha_1, alpha_2, alpha_p_1, alpha_p_2, d1_org, d2_org;
        int row_,col_;
        int atom_orb_no, atom_orb_no_p;
        if(!Parameters_.Just_Hartree){
            for(int alpha=0;alpha<NSites_in_MUC;alpha++){
                // alpha_1 = alpha % UnitCellSize_x;
                // alpha_2 = (alpha - alpha_1)/UnitCellSize_x;

                for(int atom_no=0;atom_no<n_atoms_;atom_no++){
                    for(int orb_no=0;orb_no<n_orbs_;orb_no++){
                        atom_orb_no = atom_no + n_atoms_*orb_no;
                        for(int sigma=0;sigma<2;sigma++){

                            //for(int atom_no_p=0;atom_no_p<n_atoms_;atom_no_p++){
                            int atom_no_p=atom_no;
                            for(int orb_no_p=0;orb_no_p<n_orbs_;orb_no_p++){
                                atom_orb_no_p = atom_no_p + n_atoms_*orb_no_p;
                                for(int sigma_p=0;sigma_p<2;sigma_p++){


                                    //check_= (sigma_p > sigma);
                                    row_temp = alpha + (atom_orb_no)*(S_) +  sigma*(n_atoms_*n_orbs_*S_) + 0*(2*n_atoms_*n_orbs_*S_);
                                    col_temp = alpha + (atom_orb_no_p)*(S_) +  sigma_p*(n_atoms_*n_orbs_*S_) + 0*(2*n_atoms_*n_orbs_*S_);
                                    check_ = (col_temp > row_temp);

                                    NA_spinflip = Parameters_.NoSpinFlipOP && (sigma_p!=sigma);

                                    if( check_ && (!NA_spinflip)){

                                        OPs_.value.push_back(complex<double> (0.0,0.0));
                                        OPs_new_.value.push_back(0.0);

                                        OPs_.rows.push_back(row_temp);OPs_new_.rows.push_back(row_temp);
                                        OPs_.columns.push_back(col_temp);OPs_new_.columns.push_back(col_temp);

                                        SI_to_ind[col_temp + row_temp*(ncells_*2*n_atoms_*n_orbs_*S_)] = OPs_.value.size()-1;
                                    }
                                }
                            }
                            // }
                        }
                    }
                }
            }
        }



        //Now Reading
        string File_In_Local_OP = Parameters_.File_OPs_in;
        ifstream file_in_Local_OP(File_In_Local_OP.c_str());
        string temp1, line_temp;
        int lx_old_, ly_old_;
        double temp_doub;
        complex<double> temp_complex;
        getline(file_in_Local_OP,temp1);
        getline(file_in_Local_OP,line_temp);
        stringstream line_temp_ss1(line_temp);
        line_temp_ss1>>lx_old_>>ly_old_;

        int S_= NSites_in_MUC;
        int row_temp, col_temp;
        int cell_row, alpha_plus_gamma, alpha_row,gamma_row, sigma_row, c1;
        int alpha_plus_gamma_sigma, alpha_col, gamma_col, sigma_col, cell_col;
        int alpha_row_read, atom_no_row_read, orb_no_row_read, sigma_row_read;
        int alpha_col_read, atom_no_col_read, orb_no_col_read, sigma_col_read;
        int atom_no_row, atom_no_col, orb_no_row, orb_no_col;
        for(int OP_no=0;OP_no<OPs_.value.size();OP_no++){
            row_temp=OPs_.rows[OP_no];
            col_temp=OPs_.columns[OP_no];


            //alpha + gamma*(S_) +  sigma*(n_atoms_*n_orbs_*S_) + cell_*(2*n_atoms_*n_orbs_*S_)
            cell_row = 0;
            alpha_plus_gamma = row_temp%(n_atoms_*n_orbs_*S_);
            alpha_row = alpha_plus_gamma%S_;
            gamma_row = (alpha_plus_gamma -alpha_row)/S_;
            //gamma = atom_no + n_Atoms*orb_no
            atom_no_row= gamma_row%n_atoms_;
            orb_no_row= (gamma_row - atom_no_row)/n_atoms_;
            sigma_row = (row_temp - alpha_plus_gamma)/(n_atoms_*n_orbs_*S_);
            c1 = alpha_row + gamma_row*(S_) +  sigma_row*(n_atoms_*n_orbs_*S_);
            assert(c1==row_temp);

            alpha_plus_gamma_sigma = col_temp%(2*n_atoms_*n_orbs_*S_);
            alpha_plus_gamma = alpha_plus_gamma_sigma%(n_atoms_*n_orbs_*S_);
            alpha_col = alpha_plus_gamma%S_;
            gamma_col = (alpha_plus_gamma -alpha_col)/S_;
            atom_no_col= gamma_col%n_atoms_;
            orb_no_col= (gamma_col - atom_no_col)/n_atoms_;
            sigma_col = (alpha_plus_gamma_sigma - alpha_plus_gamma)/(n_atoms_*n_orbs_*S_);
            cell_col = (col_temp - alpha_plus_gamma_sigma)/(2*n_atoms_*n_orbs_*S_);

            file_in_Local_OP>>alpha_row_read>>atom_no_row_read>>orb_no_row_read>>sigma_row_read
                    >>alpha_col_read>>atom_no_col_read>>orb_no_col_read>>sigma_col_read
                    >>temp_complex>>temp_doub;
            assert(alpha_row_read==alpha_row) ; assert(atom_no_row_read==atom_no_row) ;assert(orb_no_row_read==orb_no_row) ;assert(sigma_row_read==sigma_row);
            assert(alpha_col_read==alpha_col) ;assert(atom_no_col_read==atom_no_col) ;assert(orb_no_col_read==orb_no_col) ;assert(sigma_col_read==sigma_col);
            OPs_.value[OP_no] =temp_complex;
        }

    }


    string fl_initial_OP_out = "OP_initial.txt";
    ofstream file_initial_OP_out(fl_initial_OP_out.c_str());
    for(int i=0;i<OPs_.value.size();i++){
        file_initial_OP_out<<i<<"  "<<OPs_.rows[i]<<"  "<<OPs_.columns[i]<<"  "<<OPs_.value[i]<<endl;
    }




    cout<<"Total no. of OP's = "<<OPs_.value.size()<<endl;



}



void Kspace_calculation_G2dLatticeNew::Arranging_spectrum(){

    // Eigvectors_saved=Eigvectors_;
    Eigenvalues_saved = Eigenvalues_;

    //    Eigvectors_.resize(ncells_*6);
    //    Eigenvalues_.resize(ncells_*6);
    //    for(int i=0;i<ncells_*6;i++){
    //        Eigvectors_[i].resize(6); //Eigenvector number
    //    }

    double value_;
    //    int S_=2*S_;
    //    Mat_1_Complex_doub Vec_temp;
    //    Vec_temp.resize(S_);

    for(int i=0;i<Eigenvalues_.size();i++){

        for(int j=i+1;j<Eigenvalues_.size();j++){
            if(Eigenvalues_[j]<Eigenvalues_[i]){

                value_=Eigenvalues_[i];
                //                for(int comp=0;comp<S_;comp++){
                //                    Vec_temp[comp]=Eigvectors_[i][comp];
                //                }

                Eigenvalues_[i]=Eigenvalues_[j];
                //                for(int comp=0;comp<S_;comp++){
                //                    Eigvectors_[i][comp]=Eigvectors_[j][comp];
                //                }


                Eigenvalues_[j]=value_;
                //                for(int comp=0;comp<S_;comp++){
                //                    Eigvectors_[j][comp]=Vec_temp[comp];
                //                }
            }
        }

    }


}


void Kspace_calculation_G2dLatticeNew::Connections_bw_MUC(){

    int S_=NSites_in_MUC;
    Mat_4_Complex_doub Connections_MUC;

    Connections_MUC.resize(ncells_);
    for(int c1=0;c1<ncells_;c1++){
        Connections_MUC[c1].resize(ncells_);
        for(int c2=0;c2<ncells_;c2++){
            Connections_MUC[c1][c2].resize(S_*n_atoms_*n_orbs_*2);
            for(int comp1=0;comp1<S_*n_atoms_*n_orbs_*2;comp1++){
                Connections_MUC[c1][c2][comp1].resize(S_*n_atoms_*n_orbs_*2);
                for(int comp2=0;comp2<S_*n_atoms_*n_orbs_*2;comp2++){
                    Connections_MUC[c1][c2][comp1][comp2]=0.0;
                }
            }
        }
    }

    int d1_org, d2_org, d1_org_p, d2_org_p;
    int alpha_1, alpha_2, alpha_p_1, alpha_p_2;

    for(int cell=0;cell<ncells_;cell++){
        d1_org = Coordinates_.indx_cellwise(cell);
        d2_org = Coordinates_.indy_cellwise(cell);
        for(int cell_p=0;cell_p<ncells_;cell_p++){
            d1_org_p = Coordinates_.indx_cellwise(cell_p);
            d2_org_p = Coordinates_.indy_cellwise(cell_p);
            for(int alpha=0;alpha<S_;alpha++){
                alpha_1 = Intra_MUC_positions[alpha].first;
                alpha_2 = Intra_MUC_positions[alpha].second;
                for(int alpha_p=0;alpha_p<S_;alpha_p++){
                    alpha_p_1 = Intra_MUC_positions[alpha_p].first;
                    alpha_p_2 = Intra_MUC_positions[alpha_p].second;
                    int i1_ = (((d1_org)*Mat_MUC(0,0) + (d2_org)*Mat_MUC(1,0)) + alpha_1 + lx_*ly_)%lx_;
                    int i2_ = (((d1_org)*Mat_MUC(0,1) + (d2_org)*Mat_MUC(1,1)) + alpha_2 + lx_*ly_)%ly_;
                    int i1_p = (((d1_org_p)*Mat_MUC(0,0) + (d2_org_p)*Mat_MUC(1,0)) + alpha_p_1 + lx_*ly_)%lx_;
                    int i2_p = (((d1_org_p)*Mat_MUC(0,1) + (d2_org_p)*Mat_MUC(1,1)) + alpha_p_2 + lx_*ly_)%ly_;
                    for(int gamma=0;gamma<n_atoms_*n_orbs_;gamma++){
                        for(int gamma_p=0;gamma_p<n_atoms_*n_orbs_;gamma_p++){

                            for(int sigma=0;sigma<2;sigma++){
                                for(int sigma_p=0;sigma_p<S_;sigma_p++){

                                    int comp=alpha + S_*(gamma) + S_*n_atoms_*n_orbs_*sigma;
                                    int comp_p=alpha_p + S_*(gamma_p) + S_*n_atoms_*n_orbs_*sigma_p;

                                    int row_ = gamma  + (i1_)*n_orbs_*n_atoms_
                                            + (i2_)*lx_*n_orbs_*n_atoms_
                                            + sigma*(lx_*ly_*n_orbs_*n_atoms_);

                                    int col_ = ( gamma_p + (i1_p)*n_orbs_*n_atoms_
                                                 + ((i2_p)*lx_*n_orbs_*n_atoms_))
                                            + sigma_p*(lx_*ly_*n_orbs_*n_atoms_);

                                    Connections_MUC[cell][cell_p][comp][comp_p] += Connections_.HTB_(row_,col_);
                                }}
                        }}
                }}
        }}



    int d1_org_p_m_d1_org, d2_org_p_m_d2_org, cell_p_m_cell;
    //Checking PBC in MUC connections
    for(int cell=0;cell<ncells_;cell++){
        d1_org = Coordinates_.indx_cellwise(cell);
        d2_org = Coordinates_.indy_cellwise(cell);
        for(int cell_p=0;cell_p<ncells_;cell_p++){
            d1_org_p = Coordinates_.indx_cellwise(cell_p);
            d2_org_p = Coordinates_.indy_cellwise(cell_p);

            d1_org_p_m_d1_org = (d1_org_p - d1_org + lx_cells*ly_cells)%lx_cells;
            d2_org_p_m_d2_org = (d2_org_p - d2_org + lx_cells*ly_cells)%ly_cells;

            cell_p_m_cell = Coordinates_.Ncell(d1_org_p_m_d1_org,d2_org_p_m_d2_org);

            for(int alpha=0;alpha<S_;alpha++){
                alpha_1 = Intra_MUC_positions[alpha].first;
                alpha_2 = Intra_MUC_positions[alpha].second;
                for(int alpha_p=0;alpha_p<S_;alpha_p++){
                    alpha_p_1 = Intra_MUC_positions[alpha_p].first;
                    alpha_p_2 = Intra_MUC_positions[alpha_p].second;

                    int i1_ = (((d1_org)*Mat_MUC(0,0) + (d2_org)*Mat_MUC(1,0)) + alpha_1 + lx_*ly_)%lx_;
                    int i2_ = (((d1_org)*Mat_MUC(0,1) + (d2_org)*Mat_MUC(1,1)) + alpha_2 + lx_*ly_)%ly_;
                    int i1_p = (((d1_org_p)*Mat_MUC(0,0) + (d2_org_p)*Mat_MUC(1,0)) + alpha_p_1 + lx_*ly_)%lx_;
                    int i2_p = (((d1_org_p)*Mat_MUC(0,1) + (d2_org_p)*Mat_MUC(1,1)) + alpha_p_2 + lx_*ly_)%ly_;

                    for(int gamma=0;gamma<n_atoms_*n_orbs_;gamma++){
                        for(int gamma_p=0;gamma_p<n_atoms_*n_orbs_;gamma_p++){

                            for(int sigma=0;sigma<2;sigma++){
                                for(int sigma_p=0;sigma_p<S_;sigma_p++){

                                    int comp=alpha + S_*(gamma) + S_*n_atoms_*n_orbs_*sigma;
                                    int comp_p=alpha_p + S_*(gamma_p) + S_*n_atoms_*n_orbs_*sigma_p;


                                    if(Connections_MUC[cell][cell_p][comp][comp_p] !=  Connections_MUC[0][cell_p_m_cell][comp][comp_p]){
                                        cout<<"MUC does not respect the PBC of orginal lattice"<<endl;
                                        assert(Connections_MUC[cell][cell_p][comp][comp_p] ==  Connections_MUC[0][cell_p_m_cell][comp][comp_p]);
                                    }

                                }}
                        }}
                }}
        }}


}


complex<double> Kspace_calculation_G2dLatticeNew::h_KE(int alpha, int gamma, int sigma, int alpha_p, int gamma_p, int sigma_p, int k1, int k2){


    complex<double> temp_val;
    int row_, col_;
    int d1_, d2_;
    int d1_org, d2_org;
    int alpha_1, alpha_2, alpha_p_1, alpha_p_2;
    alpha_1 = Intra_MUC_positions[alpha].first;
    alpha_2 = Intra_MUC_positions[alpha].second;

    alpha_p_1 = Intra_MUC_positions[alpha_p].first;
    alpha_p_2 = Intra_MUC_positions[alpha_p].second;

    temp_val=0.0;
    for(int cell_no=0;cell_no<ncells_;cell_no++){
        d1_org = Coordinates_.indx_cellwise(cell_no);
        d2_org = Coordinates_.indy_cellwise(cell_no);

        int i1_ = (((d1_org)*Mat_MUC(0,0) + (d2_org)*Mat_MUC(1,0)) + alpha_1 + lx_*ly_)%lx_;
        int i2_ = (((d1_org)*Mat_MUC(0,1) + (d2_org)*Mat_MUC(1,1)) + alpha_2 + lx_*ly_)%ly_;

       // int i1_cell0_ = (0 + alpha_p_1 + lx_*ly_)%lx_;
       // int i2_cell0_ = (0 + alpha_p_2 + lx_*ly_)%ly_;

        if(i1_>=lx_ || i2_>=ly_ || i1_<0 || i2_<0){
            assert(false);
        }

        row_ = gamma  + (i1_)*n_orbs_*n_atoms_
                + (i2_)*lx_*n_orbs_*n_atoms_
                + sigma*(lx_*ly_*n_orbs_*n_atoms_);

        //col_ = ( gamma_p + (i1_cell0_)*n_orbs_*n_atoms_ + ((i2_cell0_)*lx_*n_orbs_*n_atoms_)) + sigma_p*(lx_*ly_*n_orbs_*n_atoms_);

        col_ = (gamma_p) +alpha_p*(n_atoms_*n_orbs_) + NSites_in_MUC*n_orbs_*n_atoms_*sigma_p;

        //Coordinates_.Nbasis(lx_pos, ly_pos, atom1 + n_atoms_*orb1) + nsites_ * n_orbs_* n_atoms_ * spin1;


        Get_minimum_distance_direction(0, cell_no, d1_, d2_);

        temp_val += Connections_.HTB_(row_,col_)*exp(iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_   ));

    }


    return temp_val;
}



void Kspace_calculation_G2dLatticeNew::Create_Kspace_Spectrum(){


    // bool PinningField=true;
    // double PinningFieldValue=Parameters_.PinningFieldValue;



    int S_= NSites_in_MUC;
    //Calculating density using OP's will be going inside Hamiltonian
    OPs_total_den=0.0;
    int index_OP;
    int row_OP, col_OP;
    
    for(int k1=0;k1<lx_cells;k1++){
        for(int k2=0;k2<ly_cells;k2++){
            for(int alpha=0;alpha<S_;alpha++){
                for(int gamma=0;gamma<n_atoms_*n_orbs_;gamma++){
                    for(int spin=0;spin<2;spin++){
                        row_OP = alpha + gamma*(S_)
                                + spin*(n_atoms_*n_orbs_*S_);
                        col_OP = alpha + gamma*(S_)
                                + spin*(n_atoms_*n_orbs_*S_);
                        index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                        OPs_total_den += OPs_.value[index_OP].real();
                    }
                }
            }
        }
    }






    int row_, col_;
    int UP_, DOWN_;
    UP_=0;DOWN_=1;
    int spin_bar, gamma_bar;
    int k_index;
    double k1x_val, k2x_val, k1y_val, k2y_val;
    double fac_;
    complex<double> OP_value;

    //cout <<"Here 1"<<endl;



    // Create_M_mat();
    // Create_P_mat();
    
    for(int k1=0;k1<lx_cells;k1++){
        k1x_val = (2.0*PI*k1)/(1.0*lx_cells);
        k1y_val = ((2.0*PI*k1)/(1.0*lx_cells))*(-1.0/sqrt(3));

        for(int k2=0;k2<ly_cells;k2++){
            k2x_val = 0.0;
            k2y_val = ((2.0*PI*k2)/(1.0*ly_cells))*(2.0/sqrt(3));

            k_index = Coordinates_.Ncell(k1,k2);
            Kx_values[k_index]=k1x_val + k2x_val;
            Ky_values[k_index]=k1y_val + k2y_val;


            Ham_.fill(0.0);

            //cout<<"k1, k2 :" <<k1<<" "<<k2<<endl;

            //Hoppings
            for(int alpha_p=0;alpha_p<S_;alpha_p++){
                for(int gamma_p=0;gamma_p<n_orbs_*n_atoms_;gamma_p++){
                    for(int sigma_p=0;sigma_p<2;sigma_p++){
                        col_ = alpha_p + gamma_p*(S_)
                                + sigma_p*(n_orbs_*n_atoms_*S_);

                        for(int alpha=0;alpha<S_;alpha++){
                            for(int gamma=0;gamma<n_orbs_*n_atoms_;gamma++){
                                for(int sigma=0;sigma<2;sigma++){
                                    row_ = alpha + gamma*(S_)
                                            + sigma*(n_orbs_*n_atoms_*S_);

                                    Ham_(row_, col_) += h_KE(alpha, gamma, sigma, alpha_p, gamma_p, sigma_p, k1, k2);
                                }
                            }
                        }
                    }
                }
            }





            /*
            //Anisotropy
            for(int alpha=0;alpha<S_;alpha++){
                for(int spin1=0;spin1<2;spin1++){
                    for(int spin2=0;spin2<2;spin2++){

                        row_ = alpha + spin1*(S_);
                        col_=row_;

                        row_OP = 0*(2*S_) + alpha + spin2*(S_);
                        col_OP = 0*(2*S_) + alpha + spin2*(S_);
                        index_OP = SI_to_ind[col_OP + row_OP*(2*ncells_*S_)];

                        fac_ = 1.0 - (2.0*abs(spin1-spin2)); //i.e 0 --> 1, 1 --> -1

                        Ham_(row_,col_) += (-1.0/2.0)*Parameters_.AnisotropyZ*fac_*OPs_.value[index_OP];
                    }
                }
            }
            */


            //Onsite Energies
            //OnSiteE_up
            for(int alpha=0;alpha<S_;alpha++){
                for(int gamma=0;gamma<n_atoms_*n_orbs_;gamma++){
                    for(int spin1=0;spin1<2;spin1++){
                        row_ = alpha + gamma*(S_)
                                + spin1*(n_atoms_*n_orbs_*S_);
                        col_=row_;
                        Ham_(row_,col_) += Parameters_.OnSiteE[gamma][spin1];
                    }
                }
            }


            //magnetic field along z
            for(int alpha=0;alpha<S_;alpha++){
                for(int gamma=0;gamma<n_atoms_*n_orbs_;gamma++){
                    for(int spin1=0;spin1<2;spin1++){
                        row_ = alpha + gamma*(S_)
                                + spin1*(n_atoms_*n_orbs_*S_);
                        col_=row_;
                        Ham_(row_,col_) += Parameters_.PinningFieldZ[gamma]*(0.5 - 1.0*spin1);
                    }
                }}

            //magnetic field along x
            for(int alpha=0;alpha<S_;alpha++){
                for(int gamma=0;gamma<n_atoms_*n_orbs_;gamma++){
                    for(int spin1=0;spin1<2;spin1++){
                        row_ = alpha + gamma*(S_)
                                + spin1*(n_atoms_*n_orbs_*S_);
                        col_= alpha + gamma*(S_)
                                + (1-spin1)*(n_atoms_*n_orbs_*S_);
                        Ham_(row_,col_) += Parameters_.PinningFieldX[gamma]*0.5;
                    }
                }}


            bool NA_spinflip=Parameters_.NoSpinFlipOP;
            complex<double> OP_val;
            //Interaction local intra-orbital U0:
            for(int j=0;j<S_;j++){
                for(int b_beta=0;b_beta<n_atoms_;b_beta++){
                    for(int gamma=0;gamma<n_orbs_;gamma++){
                        int atom_plus_orb=b_beta +  n_atoms_*gamma;
                        for(int spin_=0;spin_<2;spin_++){
                            col_ = j + atom_plus_orb*(S_) +
                                    spin_*(n_atoms_*n_orbs_*S_);

                            for(int spin_p=0;spin_p<2;spin_p++){

                                row_ = j + atom_plus_orb*(S_) +
                                        spin_p*(n_atoms_*n_orbs_*S_);

                                if(spin_==spin_p){
                                    row_OP = j + atom_plus_orb*(S_)
                                            + (1-spin_)*(n_atoms_*n_orbs_*S_);
                                    col_OP = row_OP;
                                    index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];

                                    Ham_(row_,col_) += Parameters_.U0[b_beta]*OPs_.value[index_OP];
                                }
                                else{
                                    if(!Parameters_.Just_Hartree){
                                        if(!NA_spinflip){
                                            row_OP = j + atom_plus_orb*(S_)
                                                    + (spin_)*(n_atoms_*n_orbs_*S_);
                                            col_OP = j + atom_plus_orb*(S_)
                                                    + (spin_p)*(n_atoms_*n_orbs_*S_);

                                            if(col_OP>row_OP){
                                                index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                                OP_val = OPs_.value[index_OP];
                                            }
                                            else{
                                                index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                                OP_val = conj(OPs_.value[index_OP]);
                                            }

                                            Ham_(row_,col_) -=Parameters_.U0[b_beta]*OP_val;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            

            //Inter-orbital nn
            //Hartree
            for(int j=0;j<S_;j++){
                for(int atom_no=0;atom_no<n_atoms_;atom_no++){
                    for(int alpha=0;alpha<n_orbs_;alpha++){
                        for(int spin_=0;spin_<2;spin_++){
                            for(int beta=(alpha+1);beta<n_orbs_;beta++){
                                for(int spin_p=0;spin_p<2;spin_p++){
                                    col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                            spin_*(n_atoms_*n_orbs_*S_);
                                    row_=col_;
                                    row_OP = j + (atom_no + n_atoms_*beta)*(S_)
                                            + (spin_p)*(n_atoms_*n_orbs_*S_);
                                    col_OP = row_OP;
                                    index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                    Ham_(row_,col_) += (Parameters_.UInterOrb[atom_no][alpha][beta])*OPs_.value[index_OP];

                                    col_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                            spin_*(n_atoms_*n_orbs_*S_);
                                    row_=col_;
                                    row_OP = j + (atom_no + n_atoms_*alpha)*(S_)
                                            + (spin_p)*(n_atoms_*n_orbs_*S_);
                                    col_OP = row_OP;
                                    index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                    Ham_(row_,col_) += (Parameters_.UInterOrb[atom_no][alpha][beta])*OPs_.value[index_OP];

                                }
                            }
                        }
                    }}}



            //Inter-orbital nn
            //Fock
            for(int j=0;j<S_;j++){
                for(int atom_no=0;atom_no<n_atoms_;atom_no++){
                    for(int alpha=0;alpha<n_orbs_;alpha++){
                        for(int beta=(alpha+1);beta<n_orbs_;beta++){

                            for(int spin_=0;spin_<2;spin_++){

                                //1 and 7 [2 and 8 are h.c. of 1 and 7 ]
                                col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                        spin_*(n_atoms_*n_orbs_*S_);
                                row_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                        spin_*(n_atoms_*n_orbs_*S_);

                                col_OP = row_;
                                row_OP = col_;

                                if(col_OP>row_OP){
                                    index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                    OP_val = OPs_.value[index_OP];
                                }
                                else{
                                    index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                    OP_val = conj(OPs_.value[index_OP]);
                                }
                                //1 and 7
                                Ham_(row_,col_) -= (Parameters_.UInterOrb[atom_no][alpha][beta])*OP_val;
                                //2 and 8
                                Ham_(col_,row_) -= conj((Parameters_.UInterOrb[atom_no][alpha][beta])*OP_val);


                                //3 and 5 [4 and 6 are h.c. of 3 and 5 ]
                                col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                        spin_*(n_atoms_*n_orbs_*S_);
                                row_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                        (1-spin_)*(n_atoms_*n_orbs_*S_);

                                col_OP = row_;
                                row_OP = col_;


                                if(!NA_spinflip){
                                    if(col_OP>row_OP){
                                        index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                        OP_val = OPs_.value[index_OP];
                                    }
                                    else{
                                        index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                        OP_val = conj(OPs_.value[index_OP]);
                                    }
                                    //3 and 5
                                    Ham_(row_,col_) -= (Parameters_.UInterOrb[atom_no][alpha][beta])*OP_val;
                                    //4 and 6
                                    Ham_(col_,row_) -= conj((Parameters_.UInterOrb[atom_no][alpha][beta])*OP_val);
                                }
                            }

                        }}}}


            //Hunds Coupling
            //Hartree
            for(int j=0;j<S_;j++){
                for(int atom_no=0;atom_no<n_atoms_;atom_no++){
                    for(int alpha=0;alpha<n_orbs_;alpha++){
                        for(int beta=(alpha+1);beta<n_orbs_;beta++){

                            for(int spin_=0;spin_<2;spin_++){

                                //1 and 3
                                col_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                        spin_*(n_atoms_*n_orbs_*S_);
                                row_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                        (1-spin_)*(n_atoms_*n_orbs_*S_);

                                col_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                        (1-spin_)*(n_atoms_*n_orbs_*S_);
                                row_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                        (spin_)*(n_atoms_*n_orbs_*S_);

                                //NA_spinflip = Parameters_.NoSpinFlipOP && (spin_p!=spin_);
                                if(!NA_spinflip){
                                    if(col_OP>row_OP){
                                        index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                        OP_val = OPs_.value[index_OP];
                                    }
                                    else{
                                        index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                        OP_val = conj(OPs_.value[index_OP]);
                                    }
                                    //1 and 3
                                    Ham_(row_,col_) -= (0.5*Parameters_.JHund[atom_no])*2.0*OP_val;
                                }


                                //2 and 4
                                col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                        spin_*(n_atoms_*n_orbs_*S_);
                                row_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                        (1-spin_)*(n_atoms_*n_orbs_*S_);

                                col_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                        (1-spin_)*(n_atoms_*n_orbs_*S_);
                                row_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                        (spin_)*(n_atoms_*n_orbs_*S_);

                                //NA_spinflip = Parameters_.NoSpinFlipOP && (spin_p!=spin_);
                                if(!NA_spinflip){
                                    if(col_OP>row_OP){
                                        index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                        OP_val = OPs_.value[index_OP];
                                    }
                                    else{
                                        index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                        OP_val = conj(OPs_.value[index_OP]);
                                    }
                                    //2 and 4
                                    Ham_(row_,col_) -= (0.5*Parameters_.JHund[atom_no])*2.0*OP_val;
                                }

                                //5 and 11
                                col_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                        spin_*(n_atoms_*n_orbs_*S_);
                                row_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                        (spin_)*(n_atoms_*n_orbs_*S_);

                                col_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                        (spin_)*(n_atoms_*n_orbs_*S_);
                                row_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                        (spin_)*(n_atoms_*n_orbs_*S_);

                                if(col_OP>row_OP){
                                    index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                    OP_val = OPs_.value[index_OP];
                                }
                                else{
                                    index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                    OP_val = conj(OPs_.value[index_OP]);
                                }
                                //5 and 11
                                Ham_(row_,col_) -= (0.5*Parameters_.JHund[atom_no])*OP_val;


                                //6 and 12
                                col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                        spin_*(n_atoms_*n_orbs_*S_);
                                row_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                        (spin_)*(n_atoms_*n_orbs_*S_);

                                col_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                        (spin_)*(n_atoms_*n_orbs_*S_);
                                row_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                        (spin_)*(n_atoms_*n_orbs_*S_);

                                if(col_OP>row_OP){
                                    index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                    OP_val = OPs_.value[index_OP];
                                }
                                else{
                                    index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                    OP_val = conj(OPs_.value[index_OP]);
                                }
                                //6 and 12
                                Ham_(row_,col_) -= (0.5*Parameters_.JHund[atom_no])*OP_val;



                                //7 and 9
                                col_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                        spin_*(n_atoms_*n_orbs_*S_);
                                row_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                        (spin_)*(n_atoms_*n_orbs_*S_);

                                col_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                        (1-spin_)*(n_atoms_*n_orbs_*S_);
                                row_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                        (1-spin_)*(n_atoms_*n_orbs_*S_);

                                if(col_OP>row_OP){
                                    index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                    OP_val = OPs_.value[index_OP];
                                }
                                else{
                                    index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                    OP_val = conj(OPs_.value[index_OP]);
                                }
                                //7 and 9
                                Ham_(row_,col_) += (0.5*Parameters_.JHund[atom_no])*OP_val;


                                //8 and 10
                                col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                        spin_*(n_atoms_*n_orbs_*S_);
                                row_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                        (spin_)*(n_atoms_*n_orbs_*S_);

                                col_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                        (1-spin_)*(n_atoms_*n_orbs_*S_);
                                row_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                        (1-spin_)*(n_atoms_*n_orbs_*S_);

                                if(col_OP>row_OP){
                                    index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                    OP_val = OPs_.value[index_OP];
                                }
                                else{
                                    index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                    OP_val = conj(OPs_.value[index_OP]);
                                }
                                //8 and 10
                                Ham_(row_,col_) += (0.5*Parameters_.JHund[atom_no])*OP_val;

                            }

                        }}}}





            //Hunds Coupling
            //Fock
            for(int j=0;j<S_;j++){
                for(int atom_no=0;atom_no<n_atoms_;atom_no++){
                    for(int alpha=0;alpha<n_orbs_;alpha++){
                        for(int beta=(alpha+1);beta<n_orbs_;beta++){

                            for(int spin_=0;spin_<2;spin_++){

                                //1 and 3 [hc : 2 and 4]
                                col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                        spin_*(n_atoms_*n_orbs_*S_);
                                row_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                        (spin_)*(n_atoms_*n_orbs_*S_);

                                col_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                        (1-spin_)*(n_atoms_*n_orbs_*S_);
                                row_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                        (1-spin_)*(n_atoms_*n_orbs_*S_);

                                if(col_OP>row_OP){
                                    index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                    OP_val = OPs_.value[index_OP];
                                }
                                else{
                                    index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                    OP_val = conj(OPs_.value[index_OP]);
                                }
                                //1 and 3
                                Ham_(row_,col_) += (0.5*Parameters_.JHund[atom_no])*2.0*OP_val;
                                //2 and 4
                                Ham_(col_,row_) += conj((0.5*Parameters_.JHund[atom_no])*2.0*OP_val);



                                //5 and 11 [hc : 6 and 12]
                                col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                        spin_*(n_atoms_*n_orbs_*S_);
                                row_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                        (spin_)*(n_atoms_*n_orbs_*S_);

                                col_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                        (spin_)*(n_atoms_*n_orbs_*S_);
                                row_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                        (spin_)*(n_atoms_*n_orbs_*S_);

                                if(col_OP>row_OP){
                                    index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                    OP_val = OPs_.value[index_OP];
                                }
                                else{
                                    index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                    OP_val = conj(OPs_.value[index_OP]);
                                }
                                //5 and 11
                                Ham_(row_,col_) += (0.5*Parameters_.JHund[atom_no])*OP_val;
                                //6 and 12
                                Ham_(col_,row_) += conj((0.5*Parameters_.JHund[atom_no])*OP_val);


                                //7 and 9 [hc : 8 and 10]
                                col_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                        spin_*(n_atoms_*n_orbs_*S_);
                                row_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                        (1-spin_)*(n_atoms_*n_orbs_*S_);

                                col_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                        (1-spin_)*(n_atoms_*n_orbs_*S_);
                                row_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                        (spin_)*(n_atoms_*n_orbs_*S_);

                                //NA_spinflip = Parameters_.NoSpinFlipOP && (spin_p!=spin_);
                                if(!NA_spinflip){
                                    if(col_OP>row_OP){
                                        index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                        OP_val = OPs_.value[index_OP];
                                    }
                                    else{
                                        index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                        OP_val = conj(OPs_.value[index_OP]);
                                    }
                                    //7 and 9
                                    Ham_(row_,col_) -= (0.5*Parameters_.JHund[atom_no])*OP_val;
                                    //8 and 10
                                    Ham_(col_,row_) -= conj((0.5*Parameters_.JHund[atom_no])*OP_val);
                                }

                            }

                        }}}}




            //Pair Hopping
            //Hartree+Fock
            for(int j=0;j<S_;j++){
                for(int atom_no=0;atom_no<n_atoms_;atom_no++){
                    for(int alpha=0;alpha<n_orbs_;alpha++){
                        for(int beta=(alpha+1);beta<n_orbs_;beta++){

                            for(int spin_=0;spin_<2;spin_++){

                                //1 and 2 [hc : 3 and 4]
                                col_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                        spin_*(n_atoms_*n_orbs_*S_);
                                row_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                        (1-spin_)*(n_atoms_*n_orbs_*S_);

                                col_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                        (1-spin_)*(n_atoms_*n_orbs_*S_);
                                row_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                        (spin_)*(n_atoms_*n_orbs_*S_);

                                //NA_spinflip = Parameters_.NoSpinFlipOP && (spin_p!=spin_);
                                if(!NA_spinflip){
                                    if(col_OP>row_OP){
                                        index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                        OP_val = OPs_.value[index_OP];
                                    }
                                    else{
                                        index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                        OP_val = conj(OPs_.value[index_OP]);
                                    }
                                    //1 and 2
                                    Ham_(row_,col_) -= (1.0*Parameters_.JHund[atom_no])*OP_val;
                                    //3 and 4
                                    Ham_(col_,row_) -= conj((1.0*Parameters_.JHund[atom_no])*OP_val);
                                }

                                //5 and 6 [hc : 7 and 8]
                                col_ = j + (atom_no + n_atoms_*beta)*(S_) +
                                        spin_*(n_atoms_*n_orbs_*S_);
                                row_ = j + (atom_no + n_atoms_*alpha)*(S_) +
                                        (spin_)*(n_atoms_*n_orbs_*S_);

                                col_OP = j + (atom_no + n_atoms_*beta)*(S_) +
                                        (1-spin_)*(n_atoms_*n_orbs_*S_);
                                row_OP = j + (atom_no + n_atoms_*alpha)*(S_) +
                                        (1-spin_)*(n_atoms_*n_orbs_*S_);

                                if(col_OP>row_OP){
                                    index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                    OP_val = OPs_.value[index_OP];
                                }
                                else{
                                    index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*S_)];
                                    OP_val = conj(OPs_.value[index_OP]);
                                }
                                //5 and 6
                                Ham_(row_,col_) += (1.0*Parameters_.JHund[atom_no])*OP_val;
                                //7 and 8
                                Ham_(col_,row_) += conj((1.0*Parameters_.JHund[atom_no])*OP_val);

                            }
                        }}}}




            //     if(k1==0 && k2==0){
            //      cout<<"---------- Printing matrix kx=0=ky -------------"<<endl;
            //      Ham_.print();
            //      cout<<"-----------------"<<endl;
            //  }

            //            cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
            //            Ham_.print();
            //            cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;


            char Dflag='V';
            Diagonalize(Dflag);


            for(int row=0;row<2*n_orbs_*n_atoms_*S_;row++){
                Eigenvalues_[2*n_orbs_*n_atoms_*S_*k_index + row]=eigs_[row];
                for(int col=0;col<2*n_orbs_*n_atoms_*S_;col++){
                    Eigvectors_[2*n_orbs_*n_atoms_*S_*k_index + col][row]=Ham_(row,col);
                    //(conj(Ham_(1,0))/abs(Ham_(1,0))) ; //fixing phase w.r.t 0th comp

                    //cout<<k_index<<"  "<<col<<"  "<<row<<"  "<<Eigvectors_[2*n_orbs_*n_atoms_*S_*k_index + col][row]<<endl;
                }
            }


        }
    }

    

}


void Kspace_calculation_G2dLatticeNew::Calculate_InteractionKernel(){


    //index = (orb1 + spin1*n_orb);

    IntKer.resize(n_atoms_);

    for(int atom_ind=0;atom_ind<n_atoms_;atom_ind++){
        IntKer[atom_ind].resize((n_orbs_*2));
        for(int i=0;i<(n_orbs_*2);i++){
            IntKer[atom_ind][i].resize((n_orbs_*2));
            for(int j=0;j<(n_orbs_*2);j++){
                IntKer[atom_ind][i][j].resize((n_orbs_*2));
                for(int k=0;k<(n_orbs_*2);k++){
                    IntKer[atom_ind][i][j][k].resize(n_orbs_*2);
                }
            }
        }
    }


    //V[Vrow_,Vcol_] c{Vrow_}^{dag} c{Vcol_}
    int Vcol_, Vrow_;

    //<c{row_OP}^{dag} c{col_OP}>
    int row_OP, col_OP;

    //Interaction local intra-orbital U0:
    for(int atom_ind=0;atom_ind<n_atoms_;atom_ind++){
        for(int gamma=0;gamma<n_orbs_;gamma++){
            // int atomplusorb = beta +  n_atoms_*gamma;
            for(int spin_=0;spin_<2;spin_++){
                Vcol_ = gamma + spin_*(n_orbs_);

                for(int spin_p=0;spin_p<2;spin_p++){
                    Vrow_ = gamma + spin_p*(n_orbs_);

                    if(spin_==spin_p){
                        row_OP = gamma + (1-spin_)*(n_orbs_);
                        col_OP = row_OP;
                        IntKer[atom_ind][Vrow_][Vcol_][row_OP][col_OP] += Parameters_.U0[atom_ind];
                    }
                    else{
                        row_OP = gamma + (spin_)*(n_orbs_);
                        col_OP = gamma + (spin_p)*(n_orbs_);

                        IntKer[atom_ind][Vrow_][Vcol_][row_OP][col_OP] -=Parameters_.U0[atom_ind];
                    }
                }
            }
        }
    }




    //Inter-orbital nn
    //Hartree
    for(int atom_ind=0;atom_ind<n_atoms_;atom_ind++){
        for(int alpha=0;alpha<n_orbs_;alpha++){
            for(int spin_=0;spin_<2;spin_++){
                for(int beta=(alpha+1);beta<n_orbs_;beta++){
                    for(int spin_p=0;spin_p<2;spin_p++){
                        Vcol_ = alpha + spin_*(n_orbs_);
                        Vrow_ = Vcol_;

                        row_OP = beta + (spin_p)*(n_orbs_);
                        col_OP = row_OP;
                        IntKer[atom_ind][Vrow_][Vcol_][row_OP][col_OP] += Parameters_.UInterOrb[atom_ind][alpha][beta];

                        Vcol_ = beta + spin_*(n_orbs_);
                        Vrow_=Vcol_;
                        row_OP = alpha + (spin_p)*(n_orbs_);
                        col_OP = row_OP;
                        IntKer[atom_ind][Vrow_][Vcol_][row_OP][col_OP] += Parameters_.UInterOrb[atom_ind][alpha][beta];

                    }
                }
            }
        }
    }



    //Inter-orbital nn
    //Fock
    for(int atom_no=0;atom_no<n_atoms_;atom_no++){
        for(int alpha=0;alpha<n_orbs_;alpha++){
            for(int beta=(alpha+1);beta<n_orbs_;beta++){
                for(int spin_=0;spin_<2;spin_++){

                    //1 and 7 [2 and 8]
                    Vrow_ = beta + spin_*(n_orbs_);
                    Vcol_ = alpha + spin_*(n_orbs_);
                    row_OP = Vcol_;
                    col_OP = Vrow_;
                    IntKer[atom_no][Vrow_][Vcol_][row_OP][col_OP] -=Parameters_.UInterOrb[atom_no][alpha][beta];
                    IntKer[atom_no][Vcol_][Vrow_][col_OP][row_OP] -=Parameters_.UInterOrb[atom_no][alpha][beta];



                    //3 and 5 [4 and 6]
                    Vrow_ = beta + spin_*(n_orbs_);
                    Vcol_ = alpha + (1-spin_)*(n_orbs_);
                    row_OP = Vcol_;
                    col_OP = Vrow_;
                    IntKer[atom_no][Vrow_][Vcol_][row_OP][col_OP] -=Parameters_.UInterOrb[atom_no][alpha][beta];
                    IntKer[atom_no][Vcol_][Vrow_][col_OP][row_OP] -=Parameters_.UInterOrb[atom_no][alpha][beta];

                }

            }}}



    //Hunds Coupling
    //Hartree
    for(int atom_no=0;atom_no<n_atoms_;atom_no++){
        for(int alpha=0;alpha<n_orbs_;alpha++){
            for(int beta=(alpha+1);beta<n_orbs_;beta++){

                for(int spin_=0;spin_<2;spin_++){

                    //1 and 3
                    Vrow_ = beta + (1-spin_)*(n_orbs_);
                    Vcol_ = beta + spin_*(n_orbs_);
                    row_OP = alpha + (spin_)*(n_orbs_);
                    col_OP = alpha + (1-spin_)*(n_orbs_);
                    IntKer[atom_no][Vrow_][Vcol_][row_OP][col_OP] -=(0.5*Parameters_.JHund[atom_no])*2.0;


                    //2 and 4
                    Vrow_ = alpha + (1-spin_)*(n_orbs_);
                    Vcol_ = alpha + spin_*(n_orbs_);
                    row_OP = beta + (spin_)*(n_orbs_);
                    col_OP = beta + (1-spin_)*(n_orbs_);
                    IntKer[atom_no][Vrow_][Vcol_][row_OP][col_OP] -=(0.5*Parameters_.JHund[atom_no])*2.0;



                    //5 and 11
                    Vrow_ = beta + (spin_)*(n_orbs_);
                    Vcol_ = beta + (spin_)*(n_orbs_);
                    row_OP = alpha + (spin_)*(n_orbs_);
                    col_OP = alpha + (spin_)*(n_orbs_);
                    IntKer[atom_no][Vrow_][Vcol_][row_OP][col_OP] -=(0.5*Parameters_.JHund[atom_no]);


                    //6 and 12
                    Vrow_ = alpha + (spin_)*(n_orbs_);
                    Vcol_ = alpha + (spin_)*(n_orbs_);
                    row_OP = beta + (spin_)*(n_orbs_);
                    col_OP = beta + (spin_)*(n_orbs_);
                    IntKer[atom_no][Vrow_][Vcol_][row_OP][col_OP] -=(0.5*Parameters_.JHund[atom_no]);


                    //7 and 9
                    Vrow_ = beta + (spin_)*(n_orbs_);
                    Vcol_ = beta + (spin_)*(n_orbs_);
                    row_OP = alpha + (1-spin_)*(n_orbs_);
                    col_OP = alpha + (1-spin_)*(n_orbs_);
                    IntKer[atom_no][Vrow_][Vcol_][row_OP][col_OP] +=(0.5*Parameters_.JHund[atom_no]);


                    //8 and 10
                    Vrow_ = alpha + (spin_)*(n_orbs_);
                    Vcol_ = alpha + (spin_)*(n_orbs_);
                    row_OP = beta + (1-spin_)*(n_orbs_);
                    col_OP = beta + (1-spin_)*(n_orbs_);
                    IntKer[atom_no][Vrow_][Vcol_][row_OP][col_OP] +=(0.5*Parameters_.JHund[atom_no]);

                }

            }}}



    //Hunds Coupling
    //Fock
    for(int atom_no=0;atom_no<n_atoms_;atom_no++){
        for(int alpha=0;alpha<n_orbs_;alpha++){
            for(int beta=(alpha+1);beta<n_orbs_;beta++){

                for(int spin_=0;spin_<2;spin_++){

                    //1 and 3
                    Vrow_ = beta + (spin_)*(n_orbs_);
                    Vcol_ = alpha + spin_*(n_orbs_);
                    row_OP = alpha + (1-spin_)*(n_orbs_);
                    col_OP = beta + (1-spin_)*(n_orbs_);
                    IntKer[atom_no][Vrow_][Vcol_][row_OP][col_OP] +=(0.5*Parameters_.JHund[atom_no])*2.0;


                    //2 and 4
                    Vrow_ = alpha + (spin_)*(n_orbs_);
                    Vcol_ = beta + spin_*(n_orbs_);
                    row_OP = beta + (1-spin_)*(n_orbs_);
                    col_OP = alpha + (1-spin_)*(n_orbs_);
                    IntKer[atom_no][Vrow_][Vcol_][row_OP][col_OP] +=(0.5*Parameters_.JHund[atom_no])*2.0;



                    //5 and 11
                    Vrow_ = beta + (spin_)*(n_orbs_);
                    Vcol_ = alpha + (spin_)*(n_orbs_);
                    row_OP = alpha + (spin_)*(n_orbs_);
                    col_OP = beta + (spin_)*(n_orbs_);
                    IntKer[atom_no][Vrow_][Vcol_][row_OP][col_OP] +=(0.5*Parameters_.JHund[atom_no]);


                    //6 and 12
                    Vrow_ = alpha + (spin_)*(n_orbs_);
                    Vcol_ = beta + (spin_)*(n_orbs_);
                    row_OP = beta + (spin_)*(n_orbs_);
                    col_OP = alpha + (spin_)*(n_orbs_);
                    IntKer[atom_no][Vrow_][Vcol_][row_OP][col_OP] +=(0.5*Parameters_.JHund[atom_no]);


                    //7 and 9
                    Vrow_ = beta + (spin_)*(n_orbs_);
                    Vcol_ = alpha + (1-spin_)*(n_orbs_);
                    row_OP = alpha + (1-spin_)*(n_orbs_);
                    col_OP = beta + (spin_)*(n_orbs_);
                    IntKer[atom_no][Vrow_][Vcol_][row_OP][col_OP] -=(0.5*Parameters_.JHund[atom_no]);


                    //8 and 10
                    Vrow_ = alpha + (spin_)*(n_orbs_);
                    Vcol_ = beta + (1-spin_)*(n_orbs_);
                    row_OP = beta + (1-spin_)*(n_orbs_);
                    col_OP = alpha + (spin_)*(n_orbs_);
                    IntKer[atom_no][Vrow_][Vcol_][row_OP][col_OP] -=(0.5*Parameters_.JHund[atom_no]);

                }

            }}}




    //Pair Hopping
    //Hartree+Fock
    for(int atom_no=0;atom_no<n_atoms_;atom_no++){
        for(int alpha=0;alpha<n_orbs_;alpha++){
            for(int beta=(alpha+1);beta<n_orbs_;beta++){

                for(int spin_=0;spin_<2;spin_++){

                    //1 and 2
                    Vrow_ = alpha + (1-spin_)*(n_orbs_);
                    Vcol_ = beta + spin_*(n_orbs_);
                    row_OP = alpha + (spin_)*(n_orbs_);
                    col_OP = beta + (1-spin_)*(n_orbs_);
                    IntKer[atom_no][Vrow_][Vcol_][row_OP][col_OP] -=(1.0*Parameters_.JHund[atom_no])*1.0;

                    //3 and 4
                    Vrow_ = beta + (1-spin_)*(n_orbs_);
                    Vcol_ = alpha + spin_*(n_orbs_);
                    row_OP = beta + (spin_)*(n_orbs_);
                    col_OP = alpha + (1-spin_)*(n_orbs_);
                    IntKer[atom_no][Vrow_][Vcol_][row_OP][col_OP] -=(1.0*Parameters_.JHund[atom_no])*1.0;


                    //5 and 6
                    Vrow_ = alpha + (spin_)*(n_orbs_);
                    Vcol_ = beta + spin_*(n_orbs_);
                    row_OP = alpha + (1-spin_)*(n_orbs_);
                    col_OP = beta + (1-spin_)*(n_orbs_);
                    IntKer[atom_no][Vrow_][Vcol_][row_OP][col_OP] +=(1.0*Parameters_.JHund[atom_no])*1.0;


                    //7 and 8
                    Vrow_ = beta + (spin_)*(n_orbs_);
                    Vcol_ = alpha + spin_*(n_orbs_);
                    row_OP = beta + (1-spin_)*(n_orbs_);
                    col_OP = alpha + (1-spin_)*(n_orbs_);
                    IntKer[atom_no][Vrow_][Vcol_][row_OP][col_OP] +=(1.0*Parameters_.JHund[atom_no])*1.0;

                }
            }}}





    for(int atom_ind=0;atom_ind<n_atoms_;atom_ind++){
        cout<<"XXXXXXXXX Interaction Kernel XXXXXXXXX"<<endl;
        for(int i=0;i<(n_orbs_*2);i++){
            for(int j=0;j<(n_orbs_*2);j++){

                for(int k=0;k<(n_orbs_*2);k++){
                    for(int l=0;l<(n_orbs_*2);l++){
                        cout<<IntKer[atom_ind][i][j][k][l]<<" ";
                    }
                }
                cout<<endl;
            }
        }
        cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
    }


}

void Kspace_calculation_G2dLatticeNew::Create_K_Path(string PathType, Mat_1_intpair & k_path){

    int kx_i, ky_i;

    k_path.clear();
    pair_int temp_pair;

    if(PathType=="GMXGYMG"){
        // ---k_path---------

        //--------\Gamma to M-----------------
        ky_i = 0;
        for (kx_i = 0; kx_i <= (lx_cells/ 2); kx_i++)
        {
            temp_pair.first = kx_i;
            temp_pair.second = kx_i;
            k_path.push_back(temp_pair);
        }
        //----------------------------------


        //--------M to X-----------------
        kx_i = (lx_cells / 2);
        for (ky_i = (ly_cells / 2)-1; ky_i >=0; ky_i--)
        {
            temp_pair.first = kx_i;
            temp_pair.second = ky_i;
            k_path.push_back(temp_pair);
        }
        //----------------------------------


        //--------X to Gamma-----------------
        ky_i = 0;
        for (kx_i = (lx_cells / 2)-1; kx_i >=0; kx_i--)
        {
            temp_pair.first = kx_i;
            temp_pair.second = ky_i;
            k_path.push_back(temp_pair);
        }
        //----------------------------------


        //--------Gamma to Y-----------------
        kx_i = 0;
        ky_i = 0;
        for (ky_i = 1; ky_i <=(ly_cells/2); ky_i++)
        {
            temp_pair.first = kx_i;
            temp_pair.second = ky_i;
            k_path.push_back(temp_pair);
        }
        //----------------------------------

        //--------Y to M-----------------
        kx_i = 0;
        ky_i = ly_cells/2;
        for (kx_i = 1; kx_i <=(lx_cells/2); kx_i++)
        {
            temp_pair.first = kx_i;
            temp_pair.second = ky_i;
            k_path.push_back(temp_pair);
        }
        //----------------------------------


        //--------M to \Gamma[with one extra point,
        //                  because in gnuplot use "set pm3d corners2color c1"
        //                  ]-----------------
        kx_i = (lx_cells / 2) - 1;
        ky_i = (ly_cells / 2) - 1;
        for (kx_i = (lx_cells / 2) - 1; kx_i >= 0; kx_i--)
        {
            ky_i=kx_i;
            if(ky_i < ly_cells){
                temp_pair.first = kx_i;
                temp_pair.second = kx_i;
                k_path.push_back(temp_pair);
            }
        }

        temp_pair.first = 0;
        temp_pair.second = 0;
        k_path.push_back(temp_pair);

    }


    if(PathType=="FullBZ"){
        // ---k_path---------

        for (ky_i = 0; ky_i < ly_cells; ky_i++)
        {
            for (kx_i = 0; kx_i < lx_cells; kx_i++)
            {
                temp_pair.first = kx_i;
                temp_pair.second = ky_i;
                k_path.push_back(temp_pair);
            }
        }

    }




}


void Kspace_calculation_G2dLatticeNew::Create_K_Path_NBZ(string PathType, Mat_1_intpair & k_path){

    int kx_i, ky_i;

    k_path.clear();
    pair_int temp_pair;



    if(PathType=="MGKMpY_HXGBZ"){
        // ---k_path---------

        //--------M to Gamma-----------------
        kx_i = lx_/2;
        ky_i = 0;

        for (kx_i = lx_/2; kx_i >=0; kx_i--)
        {
            temp_pair.first = kx_i;
            temp_pair.second = ky_i;
            k_path.push_back(temp_pair);
        }
        //----------------------------------


        //--------G to K-----------------
        kx_i = 2;
        ky_i = 1;
        for (ky_i = 1; ky_i <=ly_/3; ky_i++)
        {
            temp_pair.first = kx_i;
            temp_pair.second = ky_i;
            k_path.push_back(temp_pair);
            kx_i +=2;
        }
        //----------------------------------


        //--------K to Mp-----------------
        kx_i = 2*(lx_/3) -1;
        ky_i = (ly_/3) + 1;
        for (ky_i = (ly_/3) + 1 ; ky_i <= ly_/2; ky_i++)
        {
            temp_pair.first = kx_i;
            temp_pair.second = ky_i;
            k_path.push_back(temp_pair);
            kx_i -= 1;
        }
        //----------------------------------


        // //--------Mp to Y-----------------
        // kx_i = (lx_/2) - 1;
        // ky_i = (ly_/2) - 2;
        // for (kx_i = (lx_/2) - 1; kx_i >=(lx_/4); kx_i--)
        // {
        //     temp_pair.first = kx_i;
        //     temp_pair.second = ky_i;
        //     k_path.push_back(temp_pair);
        //     ky_i -=2;
        // }
        // //----------------------------------

        // assert(ky_i==-2);
        // assert(kx_i==(lx_/4));

        //----------------------------------


        temp_pair.first = lx_/2;
        temp_pair.second = ly_/2;
        k_path.push_back(temp_pair);

    }





    if(PathType=="GMXGYMG"){
        // ---k_path---------

        //--------\Gamma to M-----------------
        ky_i = 0;
        for (kx_i = 0; kx_i <= (lx_/ 2); kx_i++)
        {
            temp_pair.first = kx_i;
            temp_pair.second = kx_i;
            k_path.push_back(temp_pair);
        }
        //----------------------------------


        //--------M to X-----------------
        kx_i = (lx_ / 2);
        for (ky_i = (ly_ / 2)-1; ky_i >=0; ky_i--)
        {
            temp_pair.first = kx_i;
            temp_pair.second = ky_i;
            k_path.push_back(temp_pair);
        }
        //----------------------------------


        //--------X to Gamma-----------------
        ky_i = 0;
        for (kx_i = (lx_ / 2)-1; kx_i >=0; kx_i--)
        {
            temp_pair.first = kx_i;
            temp_pair.second = ky_i;
            k_path.push_back(temp_pair);
        }
        //----------------------------------


        //--------Gamma to Y-----------------
        kx_i = 0;
        ky_i = 0;
        for (ky_i = 1; ky_i <=(ly_/2); ky_i++)
        {
            temp_pair.first = kx_i;
            temp_pair.second = ky_i;
            k_path.push_back(temp_pair);
        }
        //----------------------------------

        //--------Y to M-----------------
        kx_i = 0;
        ky_i = ly_/2;
        for (kx_i = 1; kx_i <=(lx_/2); kx_i++)
        {
            temp_pair.first = kx_i;
            temp_pair.second = ky_i;
            k_path.push_back(temp_pair);
        }
        //----------------------------------


        //--------M to \Gamma[with one extra point,
        //                  because in gnuplot use "set pm3d corners2color c1"
        //                  ]-----------------
        kx_i = (lx_ / 2) - 1;
        ky_i = (ly_ / 2) - 1;
        for (kx_i = (lx_ / 2) - 1; kx_i >= 0; kx_i--)
        {
            ky_i=kx_i;
            if(ky_i < ly_){
                temp_pair.first = kx_i;
                temp_pair.second = kx_i;
                k_path.push_back(temp_pair);
            }
        }

        temp_pair.first = 0;
        temp_pair.second = 0;
        k_path.push_back(temp_pair);

    }


    if(PathType=="FullBZ"){
        // ---k_path---------

        for (ky_i = 0; ky_i < ly_; ky_i++)
        {
            for (kx_i = 0; kx_i < lx_; kx_i++)
            {
                temp_pair.first = kx_i;
                temp_pair.second = ky_i;
                k_path.push_back(temp_pair);
            }
        }

    }
}



int Kspace_calculation_G2dLatticeNew::KroneckerDelta(int i, int j){
    int val;
    if(i==j){
        val=1;
    }
    else{
        val=0;
    }
    return val;
}

void Kspace_calculation_G2dLatticeNew::Calculate_RPA_Susc(Matrix<complex<double>> Chi_bare, Matrix<complex<double>> & Chi_RPA){


    int S_=NSites_in_MUC;
    //assert(Chi_bare.n_col()==2);

    double U_temp;//=Parameters_.U0[0];
    Matrix<complex<double>> PoleMat;
    PoleMat.resize(Chi_bare.n_row(), Chi_bare.n_col());
    Chi_RPA.resize(Chi_bare.n_row(), Chi_bare.n_col());

    complex<double> Temp_PoleMat;
    Matrix<complex<double>> U_mat;
    U_mat.resize(4*S_*S_*n_atoms_*n_orbs_, 4*S_*S_*n_atoms_*n_orbs_);
    U_mat.fill(0.0);

    for(int sublattice_row1=0;sublattice_row1<S_;sublattice_row1++){
        for(int sublattice_row2=0;sublattice_row2<S_;sublattice_row2++){
            int sublattice_row = sublattice_row1 + S_*sublattice_row2;

            for(int atom_no_row=0;atom_no_row<n_atoms_;atom_no_row++){
                for(int alpha_row=0;alpha_row<n_orbs_;alpha_row++){
                    for(int spin_row1=0;spin_row1<2;spin_row1++){
                        for(int spin_row2=0;spin_row2<2;spin_row2++){
                            int row = (spin_row1+2*spin_row2) + 4*sublattice_row + 4*S_*S_*atom_no_row +
                                    4*S_*S_*n_atoms_*alpha_row;

                            for(int sublattice_col1=0;sublattice_col1<S_;sublattice_col1++){
                                for(int sublattice_col2=0;sublattice_col2<S_;sublattice_col2++){
                                    int sublattice_col = sublattice_col1 + S_*sublattice_col2;
                                    for(int atom_no_col=0;atom_no_col<n_atoms_;atom_no_col++){
                                        for(int alpha_col=0;alpha_col<n_orbs_;alpha_col++){
                                            for(int spin_col1=0;spin_col1<2;spin_col1++){
                                                for(int spin_col2=0;spin_col2<2;spin_col2++){
                                                    int col = (spin_col1+spin_col2*2) + 4*sublattice_col + 4*S_*S_*atom_no_col +
                                                            4*S_*S_*n_atoms_*alpha_col;

                                                    U_mat(row,col) =  Parameters_.U0[0]*
                                                            KroneckerDelta(sublattice_row1,sublattice_row2)*
                                                            KroneckerDelta(sublattice_row1,sublattice_col1)*
                                                            KroneckerDelta(sublattice_row1,sublattice_col2)*
                                                            //(KroneckerDelta(spin_row1, 1-spin_col1)*KroneckerDelta(spin_row1, spin_row2)*KroneckerDelta(spin_row1, 1-spin_col2))
                                                            //KroneckerDelta(spin_row1+2*spin_row2, spin_col1+2*spin_col2)
                                                            //KroneckerDelta(spin_row1, spin_row2)*KroneckerDelta(spin_row1, spin_col1)*KroneckerDelta(spin_row1, spin_col2)
                                                            (KroneckerDelta(spin_row1, 1-spin_col1)*KroneckerDelta(spin_row1, spin_row2)*KroneckerDelta(spin_row1, 1-spin_col2)
                                                             -
                                                             KroneckerDelta(spin_row1, 1-spin_col1)*KroneckerDelta(spin_row1, 1-spin_row2)*KroneckerDelta(spin_row1, spin_col2)
                                                             )
                                                            /*(KroneckerDelta(spin_row1, 1-spin_row2)*KroneckerDelta(spin_row1, spin_col1)*KroneckerDelta(spin_row1, 1-spin_col2)
                                                                                                                                                                                             -
                                                                                                                                                                                             KroneckerDelta(spin_row1, spin_row2)*KroneckerDelta(spin_row1, 1-spin_col1)*KroneckerDelta(spin_row1, 1-spin_col2)
                                                                                                                                                                                             )*/
                                                            ;

                                                }}}}}}
                        }}}}}}

    for(int sublattice_row=0;sublattice_row<S_*S_;sublattice_row++){
        for(int atom_no_row=0;atom_no_row<n_atoms_;atom_no_row++){
            for(int alpha_row=0;alpha_row<n_orbs_;alpha_row++){
                for(int spin_row1=0;spin_row1<2;spin_row1++){
                    for(int spin_row2=0;spin_row2<2;spin_row2++){
                        int row = (spin_row1+2*spin_row2) + 4*sublattice_row + 4*S_*S_*atom_no_row +
                                4*S_*S_*n_atoms_*alpha_row;


                        for(int sublattice_col=0;sublattice_col<S_*S_;sublattice_col++){
                            for(int atom_no_col=0;atom_no_col<n_atoms_;atom_no_col++){
                                for(int alpha_col=0;alpha_col<n_orbs_;alpha_col++){
                                    for(int spin_col1=0;spin_col1<2;spin_col1++){
                                        for(int spin_col2=0;spin_col2<2;spin_col2++){
                                            int col = (spin_col1+spin_col2*2) + 4*sublattice_col + 4*S_*S_*atom_no_col +
                                                    4*S_*S_*n_atoms_*alpha_col;

                                            Temp_PoleMat=0.0;
                                            for(int i=0;i<U_mat.n_col();i++){
                                                //Temp_PoleMat +=  1.0*U_mat(row,i)*Chi_bare(i,col);
                                                Temp_PoleMat -=  1.0*Chi_bare(row,i)*U_mat(i,col);
                                            }

                                            PoleMat(row,col) = (1.0*KroneckerDelta(row,col)) - Temp_PoleMat;
                                            //PoleMat(row,col) = (1.0*KroneckerDelta(row,col)) - Parameters_.U0[0]*Chi_bare(row,col);
                                        }}}}}
                    }}}}}




    Inverse(PoleMat);

    for(int row=0;row<Chi_bare.n_row();row++){
        for(int col=0;col<Chi_bare.n_col();col++){
            Chi_RPA(row,col)=0.0;
            for(int i=0;i<Chi_bare.n_row();i++){
                //Chi_RPA(row,col) += PoleMat(row,i)*Chi_bare(i,col);
                Chi_RPA(row,col) += Chi_bare(row,i)*PoleMat(i,col);
            }
        }}

}


void Kspace_calculation_G2dLatticeNew::Get_RPA_Susceptibility(){


    int S_=NSites_in_MUC;

    string File_Out_BarreSusc_str = "All_Bare_Susc_pm.txt";
    ofstream file_out_BarreSusc(File_Out_BarreSusc_str.c_str());
    file_out_BarreSusc<<"#q omega Chi_pm(q,omega)"<<endl;


    string File_Out_RPASusc_str = "RPA_Susc_pm.txt";
    ofstream file_out_RPASusc(File_Out_RPASusc_str.c_str());
    file_out_RPASusc<<"#q omega Chi_pm(q,omega)"<<endl;

    vector <Matrix<complex<double>>> PauliMatrices;
    PauliMatrices.resize(3); //X,Y,Z
    for(int i=0;i<3;i++){
        PauliMatrices[i].resize(2,2);
        PauliMatrices[i].fill(0.0);
    }

    //X
    PauliMatrices[0](0,1)=1.0;PauliMatrices[0](1,0)=1.0;

    //Y
    PauliMatrices[1](0,1)=-1.0*iota_complex;PauliMatrices[1](1,0)=iota_complex;

    //Z
    PauliMatrices[2](0,0)=1.0;PauliMatrices[2](1,1)=-1.0;


    int SPIN_UP=0;
    int SPIN_DN=1;
    Mat_1_intpair q_path;
    Create_K_Path("GMXGYMG", q_path);

    int q1_ind, q2_ind;

    Matrix<complex<double>> Chi_bare;
    Chi_bare.resize(4*S_*S_*n_atoms_*n_orbs_,4*S_*S_*n_atoms_*n_orbs_);

    Matrix<complex<double>> Chi_RPA;
    Chi_RPA.resize(4*S_*S_*n_atoms_*n_orbs_,4*S_*S_*n_atoms_*n_orbs_);




    double omega_max=10.0;
    double d_omega=0.005;
    double eta=0.01;

    complex<double> V_row, V_col, V_val;
    for(int q_ind_temp=0;q_ind_temp<q_path.size();q_ind_temp++){
        q1_ind=q_path[q_ind_temp].first;
        q2_ind=q_path[q_ind_temp].second;
        int q_index = Coordinates_.Ncell(q1_ind, q2_ind);

        double omega=-0.0;
        while(omega<omega_max){
            // file_out_BarreSusc<<q_ind_temp<<"  "<<q1_ind<<"  "<<q2_ind<<"  "<<omega<<"  ";
            file_out_RPASusc<<q_ind_temp<<"  "<<q1_ind<<"  "<<q2_ind<<"  "<<omega<<"  ";


            for(int sublattice_row1=0;sublattice_row1<S_;sublattice_row1++){
                for(int sublattice_row2=0;sublattice_row2<S_;sublattice_row2++){
                    int sublattice_row = sublattice_row1 + sublattice_row2*(S_);
                    for(int atom_no_row=0;atom_no_row<n_atoms_;atom_no_row++){
                        for(int alpha_row=0;alpha_row<n_orbs_;alpha_row++){
                            for(int spin_row1=0;spin_row1<2;spin_row1++){
                                for(int spin_row2=0;spin_row2<2;spin_row2++){
                                    int bare_row = (spin_row1+2*spin_row2) + 4*sublattice_row + 4*S_*S_*atom_no_row +
                                            4*S_*S_*n_atoms_*alpha_row;

                                    for(int sublattice_col1=0;sublattice_col1<S_;sublattice_col1++){
                                        for(int sublattice_col2=0;sublattice_col2<S_;sublattice_col2++){
                                            int sublattice_col = sublattice_col1 + sublattice_col2*(S_);
                                            for(int atom_no_col=0;atom_no_col<n_atoms_;atom_no_col++){
                                                for(int alpha_col=0;alpha_col<n_orbs_;alpha_col++){
                                                    for(int spin_col1=0;spin_col1<2;spin_col1++){
                                                        for(int spin_col2=0;spin_col2<2;spin_col2++){
                                                            int bare_col = (spin_col1+spin_col2*2) + 4*sublattice_col + 4*S_*S_*atom_no_col +
                                                                    4*S_*S_*n_atoms_*alpha_col;

                                                            Chi_bare(bare_row,bare_col)=0.0;
                                                            for(int m=0;m<2*n_atoms_*n_orbs_*S_;m++){
                                                                for(int l=0;l<2*n_atoms_*n_orbs_*S_;l++){

                                                                    for(int k1_ind=0;k1_ind<lx_cells;k1_ind++){
                                                                        for(int k2_ind=0;k2_ind<ly_cells;k2_ind++){
                                                                            int k_index = Coordinates_.Ncell(k1_ind,k2_ind);
                                                                            int state_k_m = 2*n_atoms_*n_orbs_*S_*k_index + m;

                                                                            int kpq1_ind = (k1_ind + q1_ind + lx_cells)%lx_cells;
                                                                            int kpq2_ind = (k2_ind + q2_ind + ly_cells)%ly_cells;
                                                                            int kpq_index = Coordinates_.Ncell(kpq1_ind,kpq2_ind);
                                                                            int state_kpq_l = 2*n_atoms_*n_orbs_*S_*kpq_index + l;


                                                                            /*
                                V_row = 0.5*((conj(Eigvectors_[state_q_m][sublattice_row + (atom_no_row + n_atoms_*alpha_row)*(S_) +
                                                                     SPIN_UP*(n_atoms_*n_orbs_*S_)])
                                         *
                                         Eigvectors_[state_q_l][sublattice_row + (atom_no_row + n_atoms_*alpha_row)*(S_) +
                                                                  SPIN_UP*(n_atoms_*n_orbs_*S_)]) )
                                        -
                                        0.5*((conj(Eigvectors_[state_q_m][sublattice_row + (atom_no_row + n_atoms_*alpha_row)*(S_) +
                                                                      SPIN_DN*(n_atoms_*n_orbs_*S_)])
                                          *
                                          Eigvectors_[state_q_l][sublattice_row + (atom_no_row + n_atoms_*alpha_row)*(S_) +
                                                                 SPIN_DN*(n_atoms_*n_orbs_*S_)]) );

                                V_col = 0.5*((conj(Eigvectors_[state_kpq_l][sublattice_col + (atom_no_col + n_atoms_*alpha_col)*(S_) +
                                                                            SPIN_UP*(n_atoms_*n_orbs_*S_)])
                                                *
                                                Eigvectors_[state_kpq_m][sublattice_col + (atom_no_col + n_atoms_*alpha_col)*(S_) +
                                                                       SPIN_UP*(n_atoms_*n_orbs_*S_)]) )
                                        -
                                        0.5*((conj(Eigvectors_[state_kpq_l][sublattice_col + (atom_no_col + n_atoms_*alpha_col)*(S_) +
                                                                              SPIN_DN*(n_atoms_*n_orbs_*S_)])
                                                *
                                                Eigvectors_[state_kpq_m][sublattice_col + (atom_no_col + n_atoms_*alpha_col)*(S_) +
                                                                         SPIN_DN*(n_atoms_*n_orbs_*S_)]) );

*/




                                                                            V_val = ((conj(Eigvectors_[state_k_m][sublattice_row1 + (atom_no_row + n_atoms_*alpha_row)*(S_) +
                                                                                           spin_row1*(n_atoms_*n_orbs_*S_)])*
                                                                                      Eigvectors_[state_kpq_l][sublattice_row2 + (atom_no_row + n_atoms_*alpha_row)*(S_) +
                                                                                      spin_row2*(n_atoms_*n_orbs_*S_)]) )*
                                                                                    Eigvectors_[state_k_m][sublattice_col1 + (atom_no_col + n_atoms_*alpha_col)*(S_) +
                                                                                    spin_col1*(n_atoms_*n_orbs_*S_)]*
                                                                                    conj(Eigvectors_[state_kpq_l][sublattice_col2 + (atom_no_col + n_atoms_*alpha_col)*(S_) +
                                                                                         spin_col2*(n_atoms_*n_orbs_*S_)])

                                                                                    ;




                                                                            double a_temp=(omega + Eigenvalues_saved[state_kpq_l] - Eigenvalues_saved[state_k_m]);
                                                                            double absolute_val_sqr= a_temp*a_temp + eta*eta;
                                                                            double real_part = -1.0*a_temp/absolute_val_sqr;
                                                                            double imag_part = 1.0*eta/absolute_val_sqr;

                                                                            complex<double> inverse_pole = complex<double>(real_part, imag_part);
                                                                            Chi_bare(bare_row,bare_col) += (1.0/(lx_cells*ly_cells))*
                                                                                    (-(1.0/( exp((Eigenvalues_saved[state_k_m]-mu_)*Parameters_.beta ) + 1.0)) + (1.0/( exp((Eigenvalues_saved[state_kpq_l]-mu_)*Parameters_.beta ) + 1.0)))
                                                                                    *(V_val)*inverse_pole;


                                                                            /*                    OPs_new_.value[OP_no] += (1.0/ncells_)*(
                                                 (conj(Eigvectors_[state][c1])*Eigvectors_[state][c2])*
                                                 (exp(-1.0*iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_   )  ))
                                                 *(1.0/( exp((Eigenvalues_saved[state]-mu_)*Parameters_.beta ) + 1.0))
                                                 );
*/

                                                                        }
                                                                    }
                                                                }
                                                            }
                                                            //too big file
                                                            // file_out_BarreSusc<<Chi_bare(bare_row,bare_col).real()<<"  "<<Chi_bare(bare_row,bare_col).imag()<<"  ";
                                                            //file_out_RPASusc<<Chi_RPA(bare_row,bare_col).real()<<"  "<<Chi_RPA(bare_row,bare_col).imag()<<"  ";
                                                        }}}}}}
                                }}}}}}

            Calculate_RPA_Susc(Chi_bare, Chi_RPA);

            /*
            for(int sublattice_row=0;sublattice_row<S_;sublattice_row++){
                for(int atom_no_row=0;atom_no_row<n_atoms_;atom_no_row++){
                    for(int alpha_row=0;alpha_row<n_orbs_;alpha_row++){
                        for(int spin_row=0;spin_row<3;spin_row++){
                        int bare_row = spin_row + 3*sublattice_row + 3*S_*atom_no_row +
                                           3*S_*n_atoms_*alpha_row;

                        for(int sublattice_col=0;sublattice_col<S_;sublattice_col++){
                            for(int atom_no_col=0;atom_no_col<n_atoms_;atom_no_col++){
                                for(int alpha_col=0;alpha_col<n_orbs_;alpha_col++){
                                    for(int spin_col=0;spin_col<3;spin_col++){
                                        int bare_col = spin_col + 3*sublattice_col + 3*S_*atom_no_col +
                                                       3*S_*n_atoms_*alpha_col;
                            file_out_RPASusc<<Chi_RPA(bare_row,bare_col).real()<<"  "<<Chi_RPA(bare_row,bare_col).imag()<<"  ";

                                    }}}}
                        }}}}
            */

            complex<double> temp_val;
            for(int sublattice_row1=0;sublattice_row1<S_;sublattice_row1++){
                int sublattice_row = sublattice_row1 + sublattice_row1*S_;

                for(int atom_no_row=0;atom_no_row<n_atoms_;atom_no_row++){
                    for(int alpha_row=0;alpha_row<n_orbs_;alpha_row++){

                        for(int sublattice_col1=0;sublattice_col1<S_;sublattice_col1++){
                            int sublattice_col = sublattice_col1 + sublattice_col1*S_;
                            for(int atom_no_col=0;atom_no_col<n_atoms_;atom_no_col++){
                                for(int alpha_col=0;alpha_col<n_orbs_;alpha_col++){

                                    temp_val=0.0;
                                    for(int spin_row1=0;spin_row1<2;spin_row1++){
                                        for(int spin_row2=0;spin_row2<2;spin_row2++){
                                            for(int spin_col1=0;spin_col1<2;spin_col1++){
                                                for(int spin_col2=0;spin_col2<2;spin_col2++){
                                                    int bare_row = (spin_row1 + 2*spin_row2) + 4*sublattice_row + 4*S_*S_*atom_no_row +
                                                            4*S_*S_*n_atoms_*alpha_row;

                                                    int bare_col = (spin_col1 + 2*spin_col2) + 4*sublattice_col + 4*S_*S_*atom_no_col +
                                                            4*S_*S_*n_atoms_*alpha_col;

                                                    for(int S_comp=0;S_comp<3;S_comp++){
                                                        temp_val += Chi_RPA(bare_row,bare_col)*PauliMatrices[S_comp](spin_row1,spin_row2)
                                                                *PauliMatrices[S_comp](spin_col1,spin_col2);
                                                    }

                                                }}}}
                                    file_out_RPASusc<<temp_val.real()<<"  "<<temp_val.imag()<<"  ";
                                }}}
                    }}}


            //file_out_BarreSusc<<endl;
            file_out_RPASusc<<endl;
            omega=omega+d_omega;
        }
        // file_out_BarreSusc<<endl;
        file_out_RPASusc<<endl;

    }

}





void Kspace_calculation_G2dLatticeNew::Get_RPA_Susceptibility_old(){


    int S_ = NSites_in_MUC;
    string File_Out_BarreSusc_str = "All_Bare_Susc_pm.txt";
    ofstream file_out_BarreSusc(File_Out_BarreSusc_str.c_str());
    file_out_BarreSusc<<"#q omega Chi_pm(q,omega)"<<endl;


    string File_Out_RPASusc_str = "RPA_Susc_pm.txt";
    ofstream file_out_RPASusc(File_Out_RPASusc_str.c_str());
    file_out_RPASusc<<"#q omega Chi_pm(q,omega)"<<endl;

    vector <Matrix<complex<double>>> PauliMatrices;
    PauliMatrices.resize(3); //X,Y,Z
    for(int i=0;i<3;i++){
        PauliMatrices[i].resize(2,2);
        PauliMatrices[i].fill(0.0);
    }

    //X
    PauliMatrices[0](0,1)=1.0;PauliMatrices[0](1,0)=1.0;

    //Y
    PauliMatrices[1](0,1)=-1.0*iota_complex;PauliMatrices[1](1,0)=iota_complex;

    //Z
    PauliMatrices[2](0,0)=1.0;PauliMatrices[2](1,1)=-1.0;


    int SPIN_UP=0;
    int SPIN_DN=1;
    Mat_1_intpair q_path;
    Create_K_Path("GMXGYMG", q_path);

    int q1_ind, q2_ind;

    Matrix<complex<double>> Chi_bare;
    Chi_bare.resize(3*S_*n_atoms_*n_orbs_,3*S_*n_atoms_*n_orbs_);
    Matrix<complex<double>> Chi_RPA;
    Chi_RPA.resize(3*S_*n_atoms_*n_orbs_,3*S_*n_atoms_*n_orbs_);




    double omega_max=10.0;
    double d_omega=0.002;
    double eta=0.02;

    complex<double> V_row, V_col;
    for(int q_ind_temp=0;q_ind_temp<q_path.size();q_ind_temp++){
        q1_ind=q_path[q_ind_temp].first;
        q2_ind=q_path[q_ind_temp].second;
        int q_index = Coordinates_.Ncell(q1_ind, q2_ind);

        double omega=0.0;
        while(omega<omega_max){
            file_out_BarreSusc<<q_ind_temp<<"  "<<q1_ind<<"  "<<q2_ind<<"  "<<omega<<"  ";
            file_out_RPASusc<<q_ind_temp<<"  "<<q1_ind<<"  "<<q2_ind<<"  "<<omega<<"  ";

            for(int sublattice_row=0;sublattice_row<S_;sublattice_row++){
                for(int atom_no_row=0;atom_no_row<n_atoms_;atom_no_row++){
                    for(int alpha_row=0;alpha_row<n_orbs_;alpha_row++){
                        for(int spin_row=0;spin_row<3;spin_row++){
                            int bare_row = spin_row + 3*sublattice_row + 3*S_*atom_no_row +
                                    3*S_*n_atoms_*alpha_row;

                            for(int sublattice_col=0;sublattice_col<S_;sublattice_col++){
                                for(int atom_no_col=0;atom_no_col<n_atoms_;atom_no_col++){
                                    for(int alpha_col=0;alpha_col<n_orbs_;alpha_col++){
                                        for(int spin_col=0;spin_col<3;spin_col++){
                                            int bare_col = spin_col + 3*sublattice_col + 3*S_*atom_no_col +
                                                    3*S_*n_atoms_*alpha_col;

                                            Chi_bare(bare_row,bare_col)=0.0;
                                            for(int m=0;m<2*n_atoms_*n_orbs_*S_;m++){
                                                for(int l=0;l<2*n_atoms_*n_orbs_*S_;l++){

                                                    for(int k1_ind=0;k1_ind<lx_cells;k1_ind++){
                                                        for(int k2_ind=0;k2_ind<ly_cells;k2_ind++){
                                                            int k_index = Coordinates_.Ncell(k1_ind,k2_ind);
                                                            int state_k_m = 2*n_atoms_*n_orbs_*S_*k_index + m;

                                                            int kpq1_ind = (k1_ind + q1_ind + lx_cells)%lx_cells;
                                                            int kpq2_ind = (k2_ind + q2_ind + ly_cells)%ly_cells;
                                                            int kpq_index = Coordinates_.Ncell(kpq1_ind,kpq2_ind);
                                                            int state_kpq_l = 2*n_atoms_*n_orbs_*S_*kpq_index + l;
                                                            int state_kpq_m = 2*n_atoms_*n_orbs_*S_*kpq_index + m;

                                                            int state_q_m = 2*n_atoms_*n_orbs_*S_*q_index + m;
                                                            int state_q_l = 2*n_atoms_*n_orbs_*S_*q_index + l;

                                                            /*
                                V_row = 0.5*((conj(Eigvectors_[state_q_m][sublattice_row + (atom_no_row + n_atoms_*alpha_row)*(S_) +
                                                                     SPIN_UP*(n_atoms_*n_orbs_*S_)])
                                         *
                                         Eigvectors_[state_q_l][sublattice_row + (atom_no_row + n_atoms_*alpha_row)*(S_) +
                                                                  SPIN_UP*(n_atoms_*n_orbs_*S_)]) )
                                        -
                                        0.5*((conj(Eigvectors_[state_q_m][sublattice_row + (atom_no_row + n_atoms_*alpha_row)*(S_) +
                                                                      SPIN_DN*(n_atoms_*n_orbs_*S_)])
                                          *
                                          Eigvectors_[state_q_l][sublattice_row + (atom_no_row + n_atoms_*alpha_row)*(S_) +
                                                                 SPIN_DN*(n_atoms_*n_orbs_*S_)]) );

                                V_col = 0.5*((conj(Eigvectors_[state_kpq_l][sublattice_col + (atom_no_col + n_atoms_*alpha_col)*(S_) +
                                                                            SPIN_UP*(n_atoms_*n_orbs_*S_)])
                                                *
                                                Eigvectors_[state_kpq_m][sublattice_col + (atom_no_col + n_atoms_*alpha_col)*(S_) +
                                                                       SPIN_UP*(n_atoms_*n_orbs_*S_)]) )
                                        -
                                        0.5*((conj(Eigvectors_[state_kpq_l][sublattice_col + (atom_no_col + n_atoms_*alpha_col)*(S_) +
                                                                              SPIN_DN*(n_atoms_*n_orbs_*S_)])
                                                *
                                                Eigvectors_[state_kpq_m][sublattice_col + (atom_no_col + n_atoms_*alpha_col)*(S_) +
                                                                         SPIN_DN*(n_atoms_*n_orbs_*S_)]) );

*/




                                                            V_col = PauliMatrices[spin_col](SPIN_UP,SPIN_UP)*((conj(Eigvectors_[state_kpq_l][sublattice_col + (atom_no_col + n_atoms_*alpha_col)*(S_) +
                                                                                                                    SPIN_UP*(n_atoms_*n_orbs_*S_)])*
                                                                                                               Eigvectors_[state_q_m][sublattice_col + (atom_no_col + n_atoms_*alpha_col)*(S_) +
                                                                                                               SPIN_UP*(n_atoms_*n_orbs_*S_)]) )
                                                                    +
                                                                    PauliMatrices[spin_col](SPIN_UP,SPIN_DN)*((conj(Eigvectors_[state_kpq_l][sublattice_col + (atom_no_col + n_atoms_*alpha_col)*(S_) +
                                                                                                                    SPIN_UP*(n_atoms_*n_orbs_*S_)])*
                                                                                                               Eigvectors_[state_q_m][sublattice_col + (atom_no_col + n_atoms_*alpha_col)*(S_) +
                                                                                                               SPIN_DN*(n_atoms_*n_orbs_*S_)]) )
                                                                    +
                                                                    PauliMatrices[spin_col](SPIN_DN,SPIN_UP)*((conj(Eigvectors_[state_kpq_l][sublattice_col + (atom_no_col + n_atoms_*alpha_col)*(S_) +
                                                                                                                    SPIN_DN*(n_atoms_*n_orbs_*S_)])*
                                                                                                               Eigvectors_[state_q_m][sublattice_col + (atom_no_col + n_atoms_*alpha_col)*(S_) +
                                                                                                               SPIN_UP*(n_atoms_*n_orbs_*S_)]) )
                                                                    +
                                                                    PauliMatrices[spin_col](SPIN_DN,SPIN_DN)*((conj(Eigvectors_[state_kpq_l][sublattice_col + (atom_no_col + n_atoms_*alpha_col)*(S_) +
                                                                                                                    SPIN_DN*(n_atoms_*n_orbs_*S_)])*
                                                                                                               Eigvectors_[state_q_m][sublattice_col + (atom_no_col + n_atoms_*alpha_col)*(S_) +
                                                                                                               SPIN_DN*(n_atoms_*n_orbs_*S_)]) )
                                                                    ;


                                                            V_row = PauliMatrices[spin_row](SPIN_UP,SPIN_UP)*((conj(Eigvectors_[state_kpq_l][sublattice_row + (atom_no_row + n_atoms_*alpha_row)*(S_) +
                                                                                                                    SPIN_UP*(n_atoms_*n_orbs_*S_)])*
                                                                                                               Eigvectors_[state_q_m][sublattice_row + (atom_no_row + n_atoms_*alpha_row)*(S_) +
                                                                                                               SPIN_UP*(n_atoms_*n_orbs_*S_)]) )
                                                                    +
                                                                    PauliMatrices[spin_row](SPIN_UP,SPIN_DN)*((conj(Eigvectors_[state_kpq_l][sublattice_row + (atom_no_row + n_atoms_*alpha_row)*(S_) +
                                                                                                                    SPIN_UP*(n_atoms_*n_orbs_*S_)])*
                                                                                                               Eigvectors_[state_q_m][sublattice_row + (atom_no_row + n_atoms_*alpha_row)*(S_) +
                                                                                                               SPIN_DN*(n_atoms_*n_orbs_*S_)]) )
                                                                    +
                                                                    PauliMatrices[spin_row](SPIN_DN,SPIN_UP)*(conj(Eigvectors_[state_kpq_l][sublattice_row + (atom_no_row + n_atoms_*alpha_row)*(S_) +
                                                                                                                   SPIN_DN*(n_atoms_*n_orbs_*S_)]))*
                                                                    ((Eigvectors_[state_q_m][sublattice_row + (atom_no_row + n_atoms_*alpha_row)*(S_) +
                                                                      SPIN_UP*(n_atoms_*n_orbs_*S_)]))
                                                                    +
                                                                    PauliMatrices[spin_row](SPIN_DN,SPIN_DN)*(conj(Eigvectors_[state_kpq_l][sublattice_row + (atom_no_row + n_atoms_*alpha_row)*(S_) +
                                                                                                                   SPIN_DN*(n_atoms_*n_orbs_*S_)]))*
                                                                    ((Eigvectors_[state_q_m][sublattice_row + (atom_no_row + n_atoms_*alpha_row)*(S_) +
                                                                      SPIN_DN*(n_atoms_*n_orbs_*S_)]))

                                                                    ;






                                                            double a_temp=(omega + Eigenvalues_saved[state_k_m] - Eigenvalues_saved[state_kpq_l]);
                                                            double absolute_val_sqr= a_temp*a_temp + eta*eta;
                                                            double real_part = a_temp/absolute_val_sqr;
                                                            double imag_part = eta/absolute_val_sqr;

                                                            complex<double> inverse_pole = complex<double>(real_part, imag_part);
                                                            Chi_bare(bare_row,bare_col) += (1.0/(lx_cells*ly_cells))*
                                                                    ((1.0/( exp((Eigenvalues_saved[state_k_m]-mu_)*Parameters_.beta ) + 1.0)) - (1.0/( exp((Eigenvalues_saved[state_kpq_l]-mu_)*Parameters_.beta ) + 1.0)))*
                                                                    V_row*conj(V_col)*inverse_pole;


                                                            /*                    OPs_new_.value[OP_no] += (1.0/ncells_)*(
                                                 (conj(Eigvectors_[state][c1])*Eigvectors_[state][c2])*
                                                 (exp(-1.0*iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_   )  ))
                                                 *(1.0/( exp((Eigenvalues_saved[state]-mu_)*Parameters_.beta ) + 1.0))
                                                 );
*/

                                                        }
                                                    }
                                                }
                                            }

                                            file_out_BarreSusc<<Chi_bare(bare_row,bare_col).real()<<"  "<<Chi_bare(bare_row,bare_col).imag()<<"  ";
                                            //file_out_RPASusc<<Chi_RPA(bare_row,bare_col).real()<<"  "<<Chi_RPA(bare_row,bare_col).imag()<<"  ";
                                        }}}}
                        }}}}

            Calculate_RPA_Susc(Chi_bare, Chi_RPA);

            /*
            for(int sublattice_row=0;sublattice_row<S_;sublattice_row++){
                for(int atom_no_row=0;atom_no_row<n_atoms_;atom_no_row++){
                    for(int alpha_row=0;alpha_row<n_orbs_;alpha_row++){
                        for(int spin_row=0;spin_row<3;spin_row++){
                        int bare_row = spin_row + 3*sublattice_row + 3*S_*atom_no_row +
                                           3*S_*n_atoms_*alpha_row;

                        for(int sublattice_col=0;sublattice_col<S_;sublattice_col++){
                            for(int atom_no_col=0;atom_no_col<n_atoms_;atom_no_col++){
                                for(int alpha_col=0;alpha_col<n_orbs_;alpha_col++){
                                    for(int spin_col=0;spin_col<3;spin_col++){
                                        int bare_col = spin_col + 3*sublattice_col + 3*S_*atom_no_col +
                                                       3*S_*n_atoms_*alpha_col;
                            file_out_RPASusc<<Chi_RPA(bare_row,bare_col).real()<<"  "<<Chi_RPA(bare_row,bare_col).imag()<<"  ";

                                    }}}}
                        }}}}
            */
            for(int sublattice_row=0;sublattice_row<S_;sublattice_row++){
                for(int atom_no_row=0;atom_no_row<n_atoms_;atom_no_row++){
                    for(int alpha_row=0;alpha_row<n_orbs_;alpha_row++){

                        for(int sublattice_col=0;sublattice_col<S_;sublattice_col++){
                            for(int atom_no_col=0;atom_no_col<n_atoms_;atom_no_col++){
                                for(int alpha_col=0;alpha_col<n_orbs_;alpha_col++){

                                    for(int spin_row=0;spin_row<3;spin_row++){
                                        int bare_row = spin_row + 3*sublattice_row + 3*S_*atom_no_row +
                                                3*S_*n_atoms_*alpha_row;

                                        int bare_col = spin_row + 3*sublattice_col + 3*S_*atom_no_col +
                                                3*S_*n_atoms_*alpha_col;
                                        file_out_RPASusc<<Chi_RPA(bare_row,bare_col).real()<<"  "<<Chi_RPA(bare_row,bare_col).imag()<<"  ";

                                    }
                                }}}
                    }}}


            file_out_BarreSusc<<endl;
            file_out_RPASusc<<endl;
            omega=omega+d_omega;
        }
        file_out_BarreSusc<<endl;
        file_out_RPASusc<<endl;

    }

}



void Kspace_calculation_G2dLatticeNew::Get_Bare_Susceptibility(){


    cout <<"Doing old Bare Susc. calculation"<<endl;
    string File_Out_BarreSusc_str = "Bare_Susc_pm.txt";
    ofstream file_out_BarreSusc(File_Out_BarreSusc_str.c_str());
    file_out_BarreSusc<<"#q omega Chi_pm(q,omega)"<<endl;


    int S_ = NSites_in_MUC;

    int SPIN_UP=0;
    int SPIN_DN=1;
    Mat_1_intpair q_path;
    Create_K_Path("GMXGYMG", q_path);

    int q1_ind, q2_ind;

    complex<double> Chi_;


    double omega_min=Parameters_.omega_min;
    double omega_max=Parameters_.omega_max;
    double d_omega=Parameters_.d_omega;
    double eta=Parameters_.eta;

    complex<double> V_row, V_col;
    for(int q_ind_temp=0;q_ind_temp<q_path.size();q_ind_temp++){
        q1_ind=q_path[q_ind_temp].first;
        q2_ind=q_path[q_ind_temp].second;
        int q_index = Coordinates_.Ncell(q1_ind, q2_ind);

        double omega=0.0;
        while(omega<omega_max){
            file_out_BarreSusc<<q_ind_temp<<"  "<<q1_ind<<"  "<<q2_ind<<"  "<<omega<<"  ";

            for(int sublattice_row=0;sublattice_row<S_;sublattice_row++){
                for(int atom_no_row=0;atom_no_row<n_atoms_;atom_no_row++){
                    for(int alpha_row=0;alpha_row<n_orbs_;alpha_row++){

                        for(int sublattice_col=0;sublattice_col<S_;sublattice_col++){
                            for(int atom_no_col=0;atom_no_col<n_atoms_;atom_no_col++){
                                for(int alpha_col=0;alpha_col<n_orbs_;alpha_col++){


                                    Chi_=0.0;
                                    for(int m=0;m<2*n_atoms_*n_orbs_*S_;m++){
                                        for(int l=0;l<2*n_atoms_*n_orbs_*S_;l++){

                                            for(int k1_ind=0;k1_ind<lx_cells;k1_ind++){
                                                for(int k2_ind=0;k2_ind<ly_cells;k2_ind++){
                                                    int k_index = Coordinates_.Ncell(k1_ind,k2_ind);
                                                    int state_k_m = 2*n_atoms_*n_orbs_*S_*k_index + m;

                                                    int kpq1_ind = (k1_ind + q1_ind + lx_cells)%lx_cells;
                                                    int kpq2_ind = (k2_ind + q2_ind + ly_cells)%ly_cells;
                                                    int kpq_index = Coordinates_.Ncell(kpq1_ind,kpq2_ind);
                                                    int state_kpq_l = 2*n_atoms_*n_orbs_*S_*kpq_index + l;
                                                    int state_kpq_m = 2*n_atoms_*n_orbs_*S_*kpq_index + m;

                                                    int state_q_m = 2*n_atoms_*n_orbs_*S_*q_index + m;
                                                    int state_q_l = 2*n_atoms_*n_orbs_*S_*q_index + l;

                                                    /*
                                V_row = 0.5*((conj(Eigvectors_[state_q_m][sublattice_row + (atom_no_row + n_atoms_*alpha_row)*(S_) +
                                                                     SPIN_UP*(n_atoms_*n_orbs_*S_)])
                                         *
                                         Eigvectors_[state_q_l][sublattice_row + (atom_no_row + n_atoms_*alpha_row)*(S_) +
                                                                  SPIN_UP*(n_atoms_*n_orbs_*S_)]) )
                                        -
                                        0.5*((conj(Eigvectors_[state_q_m][sublattice_row + (atom_no_row + n_atoms_*alpha_row)*(S_) +
                                                                      SPIN_DN*(n_atoms_*n_orbs_*S_)])
                                          *
                                          Eigvectors_[state_q_l][sublattice_row + (atom_no_row + n_atoms_*alpha_row)*(S_) +
                                                                 SPIN_DN*(n_atoms_*n_orbs_*S_)]) );

                                V_col = 0.5*((conj(Eigvectors_[state_kpq_l][sublattice_col + (atom_no_col + n_atoms_*alpha_col)*(S_) +
                                                                            SPIN_UP*(n_atoms_*n_orbs_*S_)])
                                                *
                                                Eigvectors_[state_kpq_m][sublattice_col + (atom_no_col + n_atoms_*alpha_col)*(S_) +
                                                                       SPIN_UP*(n_atoms_*n_orbs_*S_)]) )
                                        -
                                        0.5*((conj(Eigvectors_[state_kpq_l][sublattice_col + (atom_no_col + n_atoms_*alpha_col)*(S_) +
                                                                              SPIN_DN*(n_atoms_*n_orbs_*S_)])
                                                *
                                                Eigvectors_[state_kpq_m][sublattice_col + (atom_no_col + n_atoms_*alpha_col)*(S_) +
                                                                         SPIN_DN*(n_atoms_*n_orbs_*S_)]) );

*/


                                                    V_row = (conj(Eigvectors_[state_q_m][sublattice_row + (atom_no_row + n_atoms_*alpha_row)*(S_) +
                                                                  SPIN_UP*(n_atoms_*n_orbs_*S_)]))*
                                                            ((Eigvectors_[state_kpq_l][sublattice_row + (atom_no_row + n_atoms_*alpha_row)*(S_) +
                                                              SPIN_DN*(n_atoms_*n_orbs_*S_)]));
                                                    V_col = ((conj(Eigvectors_[state_kpq_l][sublattice_col + (atom_no_col + n_atoms_*alpha_col)*(S_) +
                                                                   SPIN_DN*(n_atoms_*n_orbs_*S_)])*
                                                              Eigvectors_[state_q_m][sublattice_col + (atom_no_col + n_atoms_*alpha_col)*(S_) +
                                                              SPIN_UP*(n_atoms_*n_orbs_*S_)]) );



                                                    double a_temp=(omega + Eigenvalues_saved[state_k_m] - Eigenvalues_saved[state_kpq_l]);
                                                    double absolute_val_sqr= a_temp*a_temp + eta*eta;
                                                    double real_part = a_temp/absolute_val_sqr;
                                                    double imag_part = eta/absolute_val_sqr;

                                                    complex<double> inverse_pole = complex<double>(real_part, imag_part);
                                                    Chi_ += (1.0/(lx_cells*ly_cells))*
                                                            ((1.0/( exp((Eigenvalues_saved[state_k_m]-mu_)*Parameters_.beta ) + 1.0)) - (1.0/( exp((Eigenvalues_saved[state_kpq_l]-mu_)*Parameters_.beta ) + 1.0)))*
                                                            V_row*(V_col)*inverse_pole;


                                                    /*                    OPs_new_.value[OP_no] += (1.0/ncells_)*(
                                                 (conj(Eigvectors_[state][c1])*Eigvectors_[state][c2])*
                                                 (exp(-1.0*iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_   )  ))
                                                 *(1.0/( exp((Eigenvalues_saved[state]-mu_)*Parameters_.beta ) + 1.0))
                                                 );
*/

                                                }
                                            }
                                        }
                                    }

                                    file_out_BarreSusc<<Chi_.real()<<"  "<<Chi_.imag()<<"  ";
                                }}}
                    }}}

            file_out_BarreSusc<<endl;
            omega=omega+d_omega;
        }
        file_out_BarreSusc<<endl;

    }

}





void Kspace_calculation_G2dLatticeNew::Get_Bare_Susceptibility_in_NBZ_old(){


    Mat_1_intpair q_path;
    Create_K_Path_NBZ("GMXGYMG", q_path);
    //Create_K_Path_NBZ("FullBZ", q_path);

    int q1_ind, q2_ind;


    double omega_min=Parameters_.omega_min;
    double omega_max=Parameters_.omega_max;
    double d_omega=Parameters_.d_omega;
    double eta=Parameters_.eta;

    int N_omega = int(((omega_max-omega_min)/d_omega) + 0.5) + 1;


    Mat_2_Complex_doub ChiBare_AB_q_omega_NBZ;

    ChiBare_AB_q_omega_NBZ.resize(lx_*ly_);
    for(int q_ind=0;q_ind<lx_*ly_;q_ind++){
        ChiBare_AB_q_omega_NBZ[q_ind].resize(N_omega);
    }


    string File_Out_BarreSusc_str = "Bare_Susc_NBZ.txt";
    ofstream file_out_BarreSusc(File_Out_BarreSusc_str.c_str());
    file_out_BarreSusc<<"#q omega Chi_bare_gh(q,omega)"<<endl;

    int n1_p, n2_p;

    for(int q_ind_temp=0;q_ind_temp<q_path.size();q_ind_temp++){
        cout<<"Chi_q_omega (NBZ) "<<q_ind_temp<<"("<<q_path.size()<<")"<<endl;

        q1_ind=q_path[q_ind_temp].first;
        q2_ind=q_path[q_ind_temp].second;
        int q_index = q1_ind + q2_ind*lx_;

        assert( (Mat_MUC(0,0)*lx_cells)%(lx_)==0);
        assert( (Mat_MUC(0,1)*lx_cells)%(ly_)==0);
        assert( (Mat_MUC(1,0)*ly_cells)%(lx_)==0);
        assert( (Mat_MUC(1,1)*ly_cells)%(ly_)==0);


        n1_p = int(((q1_ind*Mat_MUC(0,0)*lx_cells)/(lx_)) + 0.5)
                + int(((q2_ind*Mat_MUC(0,1)*lx_cells)/(ly_)) + 0.5);
        n1_p = (n1_p + lx_cells)%lx_cells;

        n2_p = int(((q1_ind*Mat_MUC(1,0)*ly_cells)/(lx_)) + 0.5)
                + int(((q2_ind*Mat_MUC(1,1)*ly_cells)/(ly_)) + 0.5);
        n2_p = (n2_p + ly_cells)%ly_cells;

        int qp_ind = Coordinates_.Ncell(n1_p, n2_p);

        for(int omega_ind=0;omega_ind<N_omega;omega_ind++){
            double omega = omega_ind*d_omega + omega_min;
            ChiBare_AB_q_omega_NBZ[q_index][omega_ind]=0.0;

            for(int g=0;g<NSites_in_MUC;g++){
                int g1=Intra_MUC_positions[g].first;
                int g2=Intra_MUC_positions[g].second;
                for(int h=0;h<NSites_in_MUC;h++){
                    int h1=Intra_MUC_positions[h].first;
                    int h2=Intra_MUC_positions[h].second;
                    ChiBare_AB_q_omega_NBZ[q_index][omega_ind] += (1.0/(NSites_in_MUC))*exp(iota_complex*2.0*PI*( (1.0*q1_ind*(g1-h1))/(lx_) + (1.0*q2_ind*(g2-h2))/(ly_) ))*
                            ChiBare_AB_q_omega[g][h][qp_ind][omega_ind];
                }
            }


            file_out_BarreSusc<<q_ind_temp<<"  "<<q1_ind<<"  "<<q2_ind<<"  "<<omega<<"  ";
            file_out_BarreSusc<<ChiBare_AB_q_omega_NBZ[q_index][omega_ind].real()<<"  "<<ChiBare_AB_q_omega_NBZ[q_index][omega_ind].imag()<<endl;


        }
        file_out_BarreSusc<<endl;

    }





}




void Kspace_calculation_G2dLatticeNew::Get_Bare_Susceptibility_in_NBZ(){


    Mat_1_intpair q_path;
    //Create_K_Path_NBZ("MGKMpY_HXGBZ", q_path);
    Create_K_Path_NBZ("GMXGYMG", q_path);
    //Create_K_Path_NBZ("FullBZ", q_path);

    int S_=NSites_in_MUC;
    double EPS_ZERO=0.0001;

    double omega_min=Parameters_.omega_min;
    double omega_max=Parameters_.omega_max;
    double d_omega=Parameters_.d_omega;

    int N_omega = int(((omega_max-omega_min)/d_omega) + 0.5) + 1;



    //Raw Bare Susc
    Mat_1_intpair q_path_MBZ;
    Create_K_Path("FullBZ", q_path_MBZ);

    for(int pair_no=0;pair_no<Parameters_.Susc_OprA.size();pair_no++){
        for(int g=0;g<NSites_in_MUC;g++){
            for(int h=0;h<NSites_in_MUC;h++){

                string File_Out_BarreSusc_str = "Bare_RawSusc_MBZ_PairNo" + to_string(pair_no)+"MUCs_" +to_string(g)+"_"+to_string(h) + ".txt";
                ofstream file_out_BarreSusc(File_Out_BarreSusc_str.c_str());
                file_out_BarreSusc<<"#q omega Chi_Bare(q,omega)"<<endl;


                Susc_OprA=Parameters_.Susc_OprA[pair_no];
                Susc_OprB=Parameters_.Susc_OprB[pair_no];

                int q1_ind, q2_ind;

                Mat_2_Complex_doub ChiRPA_AB_q_omega_NBZ; //it is bare in here :)

                ChiRPA_AB_q_omega_NBZ.resize(lx_cells*ly_cells);
                for(int q_ind=0;q_ind<lx_cells*ly_cells;q_ind++){
                    ChiRPA_AB_q_omega_NBZ[q_ind].resize(N_omega);
                }


                int n1_p, n2_p;

                for(int atom_no_alpha=0;atom_no_alpha<n_atoms_;atom_no_alpha++){
                    for(int alpha_orb=0;alpha_orb<n_orbs_;alpha_orb++){
                        for(int spin_alpha=0;spin_alpha<2;spin_alpha++){

                            for(int atom_no_beta=0;atom_no_beta<n_atoms_;atom_no_beta++){
                                for(int beta_orb=0;beta_orb<n_orbs_;beta_orb++){
                                    for(int spin_beta=0;spin_beta<2;spin_beta++){

                                        for(int atom_no_alpha_p=0;atom_no_alpha_p<n_atoms_;atom_no_alpha_p++){
                                            for(int alpha_orb_p=0;alpha_orb_p<n_orbs_;alpha_orb_p++){
                                                for(int spin_alpha_p=0;spin_alpha_p<2;spin_alpha_p++){


                                                    for(int atom_no_beta_p=0;atom_no_beta_p<n_atoms_;atom_no_beta_p++){
                                                        for(int beta_orb_p=0;beta_orb_p<n_orbs_;beta_orb_p++){
                                                            for(int spin_beta_p=0;spin_beta_p<2;spin_beta_p++){


                                                                complex<double> A_elmt = Susc_OprA[atom_no_alpha + n_atoms_*alpha_orb + n_atoms_*n_orbs_*spin_alpha][atom_no_beta + n_atoms_*beta_orb + n_atoms_*n_orbs_*spin_beta];
                                                                complex<double> B_elmt = Susc_OprB[atom_no_alpha_p + n_atoms_*alpha_orb_p + n_atoms_*n_orbs_*spin_alpha_p][atom_no_beta_p + n_atoms_*beta_orb_p + n_atoms_*n_orbs_*spin_beta_p];


                                                                if(abs(A_elmt*B_elmt)>EPS_ZERO){



                                                                    for(int q_ind_temp=0;q_ind_temp<lx_cells*ly_cells;q_ind_temp++){
                                                                        //cout<<"ChiRPA_q_omega (NBZ) "<<q_ind_temp<<"("<<q_path.size()<<")"<<endl;

                                                                        q1_ind=q_path_MBZ[q_ind_temp].first;
                                                                        q2_ind=q_path_MBZ[q_ind_temp].second;
                                                                        int q_index = q1_ind + q2_ind*lx_cells;

                                                                        for(int omega_ind=0;omega_ind<N_omega;omega_ind++){

                                                                            //ChiRPA_AB_q_omega_NBZ[q_index][omega_ind]=0.0;


                                                                            //alpha
                                                                            int ind1 = g + (atom_no_alpha + n_atoms_*alpha_orb)*(S_) +
                                                                                    spin_alpha*(n_atoms_*n_orbs_*S_);
                                                                            int ind2 = g + (atom_no_beta + n_atoms_*beta_orb)*(S_) +
                                                                                    spin_beta*(n_atoms_*n_orbs_*S_);
                                                                            int ind3 = h + (atom_no_alpha_p + n_atoms_*alpha_orb_p)*(S_) +
                                                                                    spin_alpha_p*(n_atoms_*n_orbs_*S_);
                                                                            int ind4 = h + (atom_no_beta_p + n_atoms_*beta_orb_p)*(S_) +
                                                                                    spin_beta_p*(n_atoms_*n_orbs_*S_);

                                                                            ChiRPA_AB_q_omega_NBZ[q_index][omega_ind] += A_elmt*B_elmt*ChiBareMat[q_index][omega_ind][ind1][ind2][ind3][ind4];

                                                                            // }
                                                                            //}


                                                                            //file_out_BarreSusc<<q_ind_temp<<"  "<<q1_ind<<"  "<<q2_ind<<"  "<<omega<<"  ";
                                                                            //file_out_BarreSusc<<ChiRPA_AB_q_omega_NBZ[q_index][omega_ind].real()<<"  "<<ChiRPA_AB_q_omega_NBZ[q_index][omega_ind].imag()<<endl;


                                                                        }
                                                                        //file_out_BarreSusc<<endl;

                                                                    }







                                                                }

                                                            }}}
                                                }}}
                                    }}}
                        }}}



                for(int q_ind_temp=0;q_ind_temp<lx_cells*ly_cells;q_ind_temp++){
                    //cout<<"ChiRPA_q_omega (NBZ) "<<q_ind_temp<<"("<<q_path.size()<<")"<<endl;
                    // q1_ind=q_path[q_ind_temp].first;
                    // q2_ind=q_path[q_ind_temp].second;
                    // int q_index = q1_ind + q2_ind*lx_;

                    for(int omega_ind=0;omega_ind<N_omega;omega_ind++){
                        double omega = omega_ind*d_omega + omega_min;

                        file_out_BarreSusc<<q_ind_temp<<"  "<<q1_ind<<"  "<<q2_ind<<"  "<<omega<<"  ";
                        file_out_BarreSusc<<ChiRPA_AB_q_omega_NBZ[q_ind_temp][omega_ind].real()<<"  "<<ChiRPA_AB_q_omega_NBZ[q_ind_temp][omega_ind].imag()<<endl;
                    }
                    file_out_BarreSusc<<endl;
                }



            }
        }
    }



    for(int pair_no=0;pair_no<Parameters_.Susc_OprA.size();pair_no++){

        string File_Out_BarreSusc_str = "Bare_Susc_NBZ_PairNo" + to_string(pair_no)+".txt";
        ofstream file_out_BarreSusc(File_Out_BarreSusc_str.c_str());
        file_out_BarreSusc<<"#q omega Chi_Bare(q,omega)"<<endl;


        Susc_OprA=Parameters_.Susc_OprA[pair_no];
        Susc_OprB=Parameters_.Susc_OprB[pair_no];

        int q1_ind, q2_ind;

        Mat_2_Complex_doub ChiRPA_AB_q_omega_NBZ; //it is bare in here :)

        ChiRPA_AB_q_omega_NBZ.resize(lx_*ly_);
        for(int q_ind=0;q_ind<lx_*ly_;q_ind++){
            ChiRPA_AB_q_omega_NBZ[q_ind].resize(N_omega);
        }


        int n1_p, n2_p;

        for(int atom_no_alpha=0;atom_no_alpha<n_atoms_;atom_no_alpha++){
            for(int alpha_orb=0;alpha_orb<n_orbs_;alpha_orb++){
                for(int spin_alpha=0;spin_alpha<2;spin_alpha++){

                    for(int atom_no_beta=0;atom_no_beta<n_atoms_;atom_no_beta++){
                        for(int beta_orb=0;beta_orb<n_orbs_;beta_orb++){
                            for(int spin_beta=0;spin_beta<2;spin_beta++){

                                for(int atom_no_alpha_p=0;atom_no_alpha_p<n_atoms_;atom_no_alpha_p++){
                                    for(int alpha_orb_p=0;alpha_orb_p<n_orbs_;alpha_orb_p++){
                                        for(int spin_alpha_p=0;spin_alpha_p<2;spin_alpha_p++){


                                            for(int atom_no_beta_p=0;atom_no_beta_p<n_atoms_;atom_no_beta_p++){
                                                for(int beta_orb_p=0;beta_orb_p<n_orbs_;beta_orb_p++){
                                                    for(int spin_beta_p=0;spin_beta_p<2;spin_beta_p++){


                                                        complex<double> A_elmt = Susc_OprA[atom_no_alpha + n_atoms_*alpha_orb + n_atoms_*n_orbs_*spin_alpha][atom_no_beta + n_atoms_*beta_orb + n_atoms_*n_orbs_*spin_beta];
                                                        complex<double> B_elmt = Susc_OprB[atom_no_alpha_p + n_atoms_*alpha_orb_p + n_atoms_*n_orbs_*spin_alpha_p][atom_no_beta_p + n_atoms_*beta_orb_p + n_atoms_*n_orbs_*spin_beta_p];


                                                        if(abs(A_elmt*B_elmt)>EPS_ZERO){



                                                            for(int q_ind_temp=0;q_ind_temp<q_path.size();q_ind_temp++){
                                                                //cout<<"ChiRPA_q_omega (NBZ) "<<q_ind_temp<<"("<<q_path.size()<<")"<<endl;

                                                                q1_ind=q_path[q_ind_temp].first;
                                                                q2_ind=q_path[q_ind_temp].second;
                                                                int q_index = q1_ind + q2_ind*lx_;

                                                                assert( (Mat_MUC(0,0)*lx_cells)%(lx_)==0);
                                                                assert( (Mat_MUC(0,1)*lx_cells)%(ly_)==0);
                                                                assert( (Mat_MUC(1,0)*ly_cells)%(lx_)==0);
                                                                assert( (Mat_MUC(1,1)*ly_cells)%(ly_)==0);


                                                                n1_p = int(((q1_ind*Mat_MUC(0,0)*lx_cells)/(lx_)) + 0.5)
                                                                        + int(((q2_ind*Mat_MUC(0,1)*lx_cells)/(ly_)) + 0.5);
                                                                n1_p = (n1_p + lx_cells)%lx_cells;

                                                                n2_p = int(((q1_ind*Mat_MUC(1,0)*ly_cells)/(lx_)) + 0.5)
                                                                        + int(((q2_ind*Mat_MUC(1,1)*ly_cells)/(ly_)) + 0.5);
                                                                n2_p = (n2_p + ly_cells)%ly_cells;

                                                                int qp_ind = Coordinates_.Ncell(n1_p, n2_p);

                                                                for(int omega_ind=0;omega_ind<N_omega;omega_ind++){

                                                                    ChiRPA_AB_q_omega_NBZ[q_index][omega_ind]=0.0;

                                                                    for(int g=0;g<NSites_in_MUC;g++){
                                                                        int g1=Intra_MUC_positions[g].first;
                                                                        int g2=Intra_MUC_positions[g].second;
                                                                        for(int h=0;h<NSites_in_MUC;h++){
                                                                            int h1=Intra_MUC_positions[h].first;
                                                                            int h2=Intra_MUC_positions[h].second;


                                                                            //alpha
                                                                            int ind1 = g + (atom_no_alpha + n_atoms_*alpha_orb)*(S_) +
                                                                                    spin_alpha*(n_atoms_*n_orbs_*S_);
                                                                            int ind2 = g + (atom_no_beta + n_atoms_*beta_orb)*(S_) +
                                                                                    spin_beta*(n_atoms_*n_orbs_*S_);
                                                                            int ind3 = h + (atom_no_alpha_p + n_atoms_*alpha_orb_p)*(S_) +
                                                                                    spin_alpha_p*(n_atoms_*n_orbs_*S_);
                                                                            int ind4 = h + (atom_no_beta_p + n_atoms_*beta_orb_p)*(S_) +
                                                                                    spin_beta_p*(n_atoms_*n_orbs_*S_);

                                                                            ChiRPA_AB_q_omega_NBZ[q_index][omega_ind] += A_elmt*B_elmt*(1.0/(NSites_in_MUC))*exp(iota_complex*2.0*PI*( (1.0*q1_ind*(g1-h1))/(lx_) + (1.0*q2_ind*(g2-h2))/(ly_) ))*
                                                                                    ChiBareMat[qp_ind][omega_ind][ind1][ind2][ind3][ind4];

                                                                        }
                                                                    }


                                                                    //file_out_BarreSusc<<q_ind_temp<<"  "<<q1_ind<<"  "<<q2_ind<<"  "<<omega<<"  ";
                                                                    //file_out_BarreSusc<<ChiRPA_AB_q_omega_NBZ[q_index][omega_ind].real()<<"  "<<ChiRPA_AB_q_omega_NBZ[q_index][omega_ind].imag()<<endl;


                                                                }
                                                                //file_out_BarreSusc<<endl;

                                                            }







                                                        }

                                                    }}}
                                        }}}
                            }}}
                }}}



        for(int q_ind_temp=0;q_ind_temp<q_path.size();q_ind_temp++){
            //cout<<"ChiRPA_q_omega (NBZ) "<<q_ind_temp<<"("<<q_path.size()<<")"<<endl;
            q1_ind=q_path[q_ind_temp].first;
            q2_ind=q_path[q_ind_temp].second;
            int q_index = q1_ind + q2_ind*lx_;

            for(int omega_ind=0;omega_ind<N_omega;omega_ind++){
                double omega = omega_ind*d_omega + omega_min;

                file_out_BarreSusc<<q_ind_temp<<"  "<<q1_ind<<"  "<<q2_ind<<"  "<<omega<<"  ";
                file_out_BarreSusc<<ChiRPA_AB_q_omega_NBZ[q_index][omega_ind].real()<<"  "<<ChiRPA_AB_q_omega_NBZ[q_index][omega_ind].imag()<<endl;
            }
            file_out_BarreSusc<<endl;
        }



    }

}


/*
void Kspace_calculation_G2dLatticeNew::Get_BareOrbital_Susceptibility(){


    double EPS_ZERO=0.0001;

    int S_ = NSites_in_MUC;

    Mat_1_intpair q_path;

    int q1_ind, q2_ind;

    double omega_min=Parameters_.omega_min;
    double omega_max=Parameters_.omega_max;
    double d_omega=Parameters_.d_omega;
    double eta=Parameters_.eta;

    int N_omega = int(((omega_max-omega_min)/d_omega) + 0.5) + 1;

    //ind = uc_no + (atom + orb*n_atoms)*(S_) +  sigma*(n_atoms_*n_orbs_*S_)



    int ChiMatsize = S_*n_atoms_*n_orbs_*2;





#ifdef _OPENMP
        int N_p = omp_get_max_threads();
        omp_set_num_threads(Parameters_.NProcessors);

        cout<<"Max threads which can be used parallely = "<<N_p<<endl;
        cout<<"No. of threads used parallely = "<<Parameters_.NProcessors<<endl;
#endif



    cout<<"CALCULATING BARE CHI"<<endl;
    //Calculate ChiBare

int thread_id=0;
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int q_ind_temp=0;q_ind_temp<lx_cells*ly_cells;q_ind_temp++){
#ifdef _OPENMP
            thread_id = omp_get_thread_num();
#endif
    cout<<"q_ind = "<<q_ind_temp<<" ("<<q_path.size()<<") in thread "<<thread_id<<endl;



    int q1_ind_local=q_path[q_ind_temp].first;
    int q2_ind_local=q_path[q_ind_temp].second;
    int q_index = Coordinates_.Ncell(q1_ind_local, q2_ind_local);

    for(int g=0;g<S_;g++){
        for(int m=0;m<n_atoms_;m++){
            for(int alpha=0;alpha<n_orbs_;alpha++){
                for(int sigma=0;sigma<2;sigma++){
                    int ind1 = g + (m + alpha*n_atoms_)*S_ + sigma*(n_atoms_*n_orbs_*S_);

                    int h=g;
                    int n=m;
                  //  for(int h=0;h<S_;h++){??
                    //    for(int n=0;n<n_atoms_;n++){
                            for(int beta=0;beta<n_orbs_;beta++){
                                for(int sigma_til=0;sigma_til<2;sigma_til++){
                                    int ind2 = h + (n + beta*n_atoms_)*S_ + sigma_til*(n_atoms_*n_orbs_*S_);


                                    for(int hp=0;hp<S_;hp++){
                                        for(int np=0;np<n_atoms_;np++){
                                            for(int alphap=0;alphap<n_orbs_;alphap++){
                                                for(int sigmap=0;sigmap<2;sigmap++){
                                                    int ind3 = hp + (np + alphap*n_atoms_)*S_ + sigmap*(n_atoms_*n_orbs_*S_);



                                                    for(int betap=0;betap<n_orbs_;betap++){
                                                        for(int sigma_tilp=0;sigma_tilp<2;sigma_tilp++){
                                                            int ind4 = hp + (np + betap*n_atoms_)*S_ + sigma_tilp*(n_atoms_*n_orbs_*S_);

                                                                cout<<"BARE CHI : "<<ind1<<" ("<<S_*n_atoms_*n_orbs_*2<<") "<<ind2<<" ("<<S_*n_atoms_*n_orbs_*2<<") "<<ind3<<" ("<<S_*n_atoms_*n_orbs_*2<<") "<<betap + sigma_tilp*n_orbs_<<" ("<<n_orbs_*2<<") "<<q_ind_temp<<" ("<<q_path.size()<<") "<<endl;


                                                                for(int omega_ind=0;omega_ind<N_omega;omega_ind++){
                                                                    double omega = omega_ind*d_omega + omega_min;
                                                                    //   file_out_BarreSusc<<q_ind_temp<<"  "<<q1_ind<<"  "<<q2_ind<<"  "<<omega<<"  ";

                                                                    ChiBareMat[q_index][omega_ind][ind1][ind2][ind3][ind4]=0.0;

                                                                    for(int k1_ind=0;k1_ind<lx_cells;k1_ind++){
                                                                        for(int k2_ind=0;k2_ind<ly_cells;k2_ind++){

                                                                            int kpq1_ind = (k1_ind + q1_ind_local + lx_cells)%lx_cells;
                                                                            int kpq2_ind = (k2_ind + q2_ind_local + ly_cells)%ly_cells;

                                                                            int k_index = Coordinates_.Ncell(k1_ind,k2_ind);
                                                                            int kpq_index = Coordinates_.Ncell(kpq1_ind,kpq2_ind);





                                                                            for(int lambda=0;lambda<2*n_atoms_*n_orbs_*S_;lambda++){
                                                                                for(int lambdap=0;lambdap<2*n_atoms_*n_orbs_*S_;lambdap++){
                                                                                    int state_k_lambda = 2*n_atoms_*n_orbs_*S_*k_index + lambda;
                                                                                    int state_kpq_lambdap = 2*n_atoms_*n_orbs_*S_*kpq_index + lambdap;

                                                                                    double a_temp=(omega + Eigenvalues_saved[state_k_lambda] - Eigenvalues_saved[state_kpq_lambdap]);
                                                                                    double absolute_val_sqr= a_temp*a_temp + eta*eta;
                                                                                    double real_part = 1.0*a_temp/absolute_val_sqr;
                                                                                    double imag_part = -1.0*eta/absolute_val_sqr;

                                                                                    complex<double> inverse_pole = complex<double>(real_part, imag_part);

                                                                                    ChiBareMat[q_index][omega_ind][ind1][ind2][ind3][ind4] += ((1.0/(lx_cells*ly_cells)))*((1.0/( exp((Eigenvalues_saved[state_k_lambda]-mu_)*Parameters_.beta ) + 1.0)) - (1.0/( exp((Eigenvalues_saved[state_kpq_lambdap]-mu_)*Parameters_.beta ) + 1.0)))
                                                                                            *inverse_pole
                                                                                            *conj(Eigvectors_[state_k_lambda][ind1])*Eigvectors_[state_kpq_lambdap][ind2]
                                                                                            *conj(Eigvectors_[state_kpq_lambdap][ind3])*(Eigvectors_[state_k_lambda][ind4]);



                                                                                    //ChiBare_AB_q_omega[IntraMUC_indexA][IntraMUC_indexB][q_index][omega_ind] += ((1.0/(lx_cells*ly_cells)))*((1.0/( exp((Eigenvalues_saved[state_k_n]-mu_)*Parameters_.beta ) + 1.0)) - (1.0/( exp((Eigenvalues_saved[state_kpq_m]-mu_)*Parameters_.beta ) + 1.0)))*
                                                                                    //        ABmat[IntraMUC_indexA][IntraMUC_indexB][state_k_n][state_kpq_m]*inverse_pole;

                                                                                }}
                                                                        }}
                                                                }
                                                            //}
                                                        }
                                                    }

                                                }
                                            }
                                        }
                                    }

                                }
                            }
                      //  }
                    //}



                }
            }
        }
    }
}


    cout<<"CALCULATING OneMinusChiBareTimesI_Inv"<<endl;
    //Calculate ChiBareTimesI and OneMinusChiBareTimesI_Inv
    for(int q_ind_temp=0;q_ind_temp<q_path.size();q_ind_temp++){
        cout<<"OneMinusChiBareTimesI_Inv   "<<q_ind_temp<<" ("<<q_path.size()<<") "<<endl;

        q1_ind=q_path[q_ind_temp].first;
        q2_ind=q_path[q_ind_temp].second;
        int q_index = Coordinates_.Ncell(q1_ind, q2_ind);

        for(int omega_ind=0;omega_ind<N_omega;omega_ind++){

            for(int g=0;g<S_;g++){
                for(int m=0;m<n_atoms_;m++){
                    for(int alpha=0;alpha<n_orbs_;alpha++){
                        for(int sigma=0;sigma<2;sigma++){
                            int ind1 = g + (m + alpha*n_atoms_)*S_ + sigma*(n_atoms_*n_orbs_*S_);


                            int h=g;
                            int n=m;
                            // for(int h=0;h<S_;h++){
                            //     for(int n=0;n<n_atoms_;n++){
                                    for(int beta=0;beta<n_orbs_;beta++){
                                        for(int sigma_til=0;sigma_til<2;sigma_til++){
                                            int ind2 = h + (n + beta*n_atoms_)*S_ + sigma_til*(n_atoms_*n_orbs_*S_);

                                            for(int hp=0;hp<S_;hp++){
                                                for(int np=0;np<n_atoms_;np++){
                                                    for(int alphap=0;alphap<n_orbs_;alphap++){
                                                        for(int sigmap=0;sigmap<2;sigmap++){
                                                            int ind3 = hp + (np + alphap*n_atoms_)*S_ + sigmap*(n_atoms_*n_orbs_*S_);


                                                            for(int betap=0;betap<n_orbs_;betap++){
                                                                for(int sigma_tilp=0;sigma_tilp<2;sigma_tilp++){
                                                                    int ind4 = hp + (np + betap*n_atoms_)*S_ + sigma_tilp*(n_atoms_*n_orbs_*S_);

                                                                    ChiBareTimesI[q_index][omega_ind][ind1][ind2][ind3][ind4]=0.0;

                                                                    for(int alphapp=0;alphapp<n_orbs_;alphapp++){
                                                                        for(int sigmapp=0;sigmapp<2;sigmapp++){
                                                                            int ind3_p = hp + (np + alphapp*n_atoms_)*S_ + sigmapp*(n_atoms_*n_orbs_*S_);

                                                                            for(int betapp=0;betapp<n_orbs_;betapp++){
                                                                                for(int sigma_tilpp=0;sigma_tilpp<2;sigma_tilpp++){
                                                                                    int ind4_p = hp + (np + betapp*n_atoms_)*S_ + sigma_tilpp*(n_atoms_*n_orbs_*S_);

                                                                                    ChiBareTimesI[q_index][omega_ind][ind1][ind2][ind3][ind4] += ChiBareMat[q_index][omega_ind][ind1][ind2][ind3_p][ind4_p]
                                                                                            *IntKer[np][alphapp+sigmapp*n_orbs_][betapp+sigma_tilpp*n_orbs_][alphap+sigmap*n_orbs_][betap+sigma_tilp*n_orbs_];

                                                                                }
                                                                            }
                                                                        }
                                                                    }

                                                                    OneMinusChiBareTimesI_Inv[q_index][omega_ind](ind1 + ind2*ChiMatsize,ind3 + ind4*ChiMatsize) -= ChiBareTimesI[q_index][omega_ind][ind1][ind2][ind3][ind4];

                                                                }
                                                            }


                                                        }
                                                    }
                                                }
                                            }


                                        }
                                    }
                               // }
                           // }



                        }
                    }
                }
            }


            Inverse(OneMinusChiBareTimesI_Inv[q_index][omega_ind]);

        }}





    //Chi_RPA_Mat
    for(int q_ind_temp=0;q_ind_temp<q_path.size();q_ind_temp++){
        q1_ind=q_path[q_ind_temp].first;
        q2_ind=q_path[q_ind_temp].second;
        int q_index = Coordinates_.Ncell(q1_ind, q2_ind);

        for(int omega_ind=0;omega_ind<N_omega;omega_ind++){

            for(int g=0;g<S_;g++){
                for(int m=0;m<n_atoms_;m++){
                    for(int alpha=0;alpha<n_orbs_;alpha++){
                        for(int sigma=0;sigma<2;sigma++){
                            int ind1 = g + (m + alpha*n_atoms_)*S_ + sigma*(n_atoms_*n_orbs_*S_);

                            for(int beta=0;beta<n_orbs_;beta++){
                                for(int sigma_til=0;sigma_til<2;sigma_til++){
                                    int ind2 = g + (m + beta*n_atoms_)*S_ + sigma_til*(n_atoms_*n_orbs_*S_);

                                    for(int h=0;h<S_;h++){
                                        for(int n=0;n<n_atoms_;n++){
                                            for(int alphap=0;alphap<n_orbs_;alphap++){
                                                for(int sigmap=0;sigmap<2;sigmap++){
                                                    int ind3 = h + (n + alphap*n_atoms_)*S_ + sigmap*(n_atoms_*n_orbs_*S_);


                                                    for(int betap=0;betap<n_orbs_;betap++){
                                                        for(int sigma_tilp=0;sigma_tilp<2;sigma_tilp++){
                                                            int ind4 = h + (n + betap*n_atoms_)*S_ + sigma_tilp*(n_atoms_*n_orbs_*S_);

                                                            ChiRPAMat[q_index][omega_ind][ind1][ind2][ind3][ind4]=0.0;


                                                            for(int hp=0;hp<S_;hp++){
                                                                for(int np=0;np<n_atoms_;np++){
                                                                    for(int alphapp=0;alphapp<n_orbs_;alphapp++){
                                                                        for(int sigmapp=0;sigmapp<2;sigmapp++){
                                                                            int ind3_p = hp + (np + alphapp*n_atoms_)*S_ + sigmapp*(n_atoms_*n_orbs_*S_);

                                                                            for(int gp=0;gp<S_;gp++){
                                                                                for(int mp=0;mp<n_atoms_;mp++){
                                                                                    for(int betapp=0;betapp<n_orbs_;betapp++){
                                                                                        for(int sigma_tilpp=0;sigma_tilpp<2;sigma_tilpp++){
                                                                                            int ind4_p = gp + (mp + betapp*n_atoms_)*S_ + sigma_tilpp*(n_atoms_*n_orbs_*S_);


                                                                                            ChiRPAMat[q_index][omega_ind][ind1][ind2][ind3][ind4] += OneMinusChiBareTimesI_Inv[q_index][omega_ind](ind1 + ChiMatsize*ind2, ind3_p + ChiMatsize*ind4_p)
                                                                                                    *ChiBareMat[q_index][omega_ind][ind3_p][ind4_p][ind3][ind4];

                                                                                        }}}}
                                                                        }}}}
                                                        }}
                                                }}}}
                                }}
                        }}}}
        }}



}


*/


void Kspace_calculation_G2dLatticeNew::Get_RPA_Susceptibility_in_NBZ(){


    Mat_1_intpair q_path;
    //Create_K_Path_NBZ("MGKMpY_HXGBZ", q_path);
    Create_K_Path_NBZ("GMXGYMG", q_path);
    //Create_K_Path_NBZ("FullBZ", q_path);

    int S_=NSites_in_MUC;
    double EPS_ZERO=0.0001;

    double omega_min=Parameters_.omega_min;
    double omega_max=Parameters_.omega_max;
    double d_omega=Parameters_.d_omega;

    int N_omega = int(((omega_max-omega_min)/d_omega) + 0.5) + 1;



    for(int pair_no=0;pair_no<Parameters_.Susc_OprA.size();pair_no++){

        string File_Out_BarreSusc_str = "RPA_Susc_NBZ_PairNo" + to_string(pair_no)+".txt";
        ofstream file_out_BarreSusc(File_Out_BarreSusc_str.c_str());
        file_out_BarreSusc<<"#q omega Chi_RPA_gh(q,omega)"<<endl;


        Susc_OprA=Parameters_.Susc_OprA[pair_no];
        Susc_OprB=Parameters_.Susc_OprB[pair_no];

        int q1_ind, q2_ind;

        Mat_2_Complex_doub ChiRPA_AB_q_omega_NBZ;

        ChiRPA_AB_q_omega_NBZ.resize(lx_*ly_);
        for(int q_ind=0;q_ind<lx_*ly_;q_ind++){
            ChiRPA_AB_q_omega_NBZ[q_ind].resize(N_omega);
        }


        int n1_p, n2_p;

        for(int atom_no_alpha=0;atom_no_alpha<n_atoms_;atom_no_alpha++){
            for(int alpha_orb=0;alpha_orb<n_orbs_;alpha_orb++){
                for(int spin_alpha=0;spin_alpha<2;spin_alpha++){

                    for(int atom_no_beta=0;atom_no_beta<n_atoms_;atom_no_beta++){
                        for(int beta_orb=0;beta_orb<n_orbs_;beta_orb++){
                            for(int spin_beta=0;spin_beta<2;spin_beta++){

                                for(int atom_no_alpha_p=0;atom_no_alpha_p<n_atoms_;atom_no_alpha_p++){
                                    for(int alpha_orb_p=0;alpha_orb_p<n_orbs_;alpha_orb_p++){
                                        for(int spin_alpha_p=0;spin_alpha_p<2;spin_alpha_p++){


                                            for(int atom_no_beta_p=0;atom_no_beta_p<n_atoms_;atom_no_beta_p++){
                                                for(int beta_orb_p=0;beta_orb_p<n_orbs_;beta_orb_p++){
                                                    for(int spin_beta_p=0;spin_beta_p<2;spin_beta_p++){


                                                        complex<double> A_elmt = Susc_OprA[atom_no_alpha + n_atoms_*alpha_orb + n_atoms_*n_orbs_*spin_alpha][atom_no_beta + n_atoms_*beta_orb + n_atoms_*n_orbs_*spin_beta];
                                                        complex<double> B_elmt = Susc_OprB[atom_no_alpha_p + n_atoms_*alpha_orb_p + n_atoms_*n_orbs_*spin_alpha_p][atom_no_beta_p + n_atoms_*beta_orb_p + n_atoms_*n_orbs_*spin_beta_p];


                                                        if(abs(A_elmt*B_elmt)>EPS_ZERO){



                                                            for(int q_ind_temp=0;q_ind_temp<q_path.size();q_ind_temp++){
                                                                //cout<<"ChiRPA_q_omega (NBZ) "<<q_ind_temp<<"("<<q_path.size()<<")"<<endl;

                                                                q1_ind=q_path[q_ind_temp].first;
                                                                q2_ind=q_path[q_ind_temp].second;
                                                                int q_index = q1_ind + q2_ind*lx_;

                                                                assert( (Mat_MUC(0,0)*lx_cells)%(lx_)==0);
                                                                assert( (Mat_MUC(0,1)*lx_cells)%(ly_)==0);
                                                                assert( (Mat_MUC(1,0)*ly_cells)%(lx_)==0);
                                                                assert( (Mat_MUC(1,1)*ly_cells)%(ly_)==0);


                                                                n1_p = int(((q1_ind*Mat_MUC(0,0)*lx_cells)/(lx_)) + 0.5)
                                                                        + int(((q2_ind*Mat_MUC(0,1)*lx_cells)/(ly_)) + 0.5);
                                                                n1_p = (n1_p + lx_cells)%lx_cells;

                                                                n2_p = int(((q1_ind*Mat_MUC(1,0)*ly_cells)/(lx_)) + 0.5)
                                                                        + int(((q2_ind*Mat_MUC(1,1)*ly_cells)/(ly_)) + 0.5);
                                                                n2_p = (n2_p + ly_cells)%ly_cells;

                                                                int qp_ind = Coordinates_.Ncell(n1_p, n2_p);

                                                                for(int omega_ind=0;omega_ind<N_omega;omega_ind++){

                                                                    // ChiRPA_AB_q_omega_NBZ[q_index][omega_ind]=0.0;

                                                                    for(int g=0;g<NSites_in_MUC;g++){
                                                                        int g1=Intra_MUC_positions[g].first;
                                                                        int g2=Intra_MUC_positions[g].second;
                                                                        for(int h=0;h<NSites_in_MUC;h++){
                                                                            int h1=Intra_MUC_positions[h].first;
                                                                            int h2=Intra_MUC_positions[h].second;


                                                                            //alpha
                                                                            int ind1 = g + (atom_no_alpha + n_atoms_*alpha_orb)*(S_) +
                                                                                    spin_alpha*(n_atoms_*n_orbs_*S_);
                                                                            int ind2 = g + (atom_no_beta + n_atoms_*beta_orb)*(S_) +
                                                                                    spin_beta*(n_atoms_*n_orbs_*S_);
                                                                            int ind3 = h + (atom_no_alpha_p + n_atoms_*alpha_orb_p)*(S_) +
                                                                                    spin_alpha_p*(n_atoms_*n_orbs_*S_);
                                                                            int ind4 = h + (atom_no_beta_p + n_atoms_*beta_orb_p)*(S_) +
                                                                                    spin_beta_p*(n_atoms_*n_orbs_*S_);

                                                                            ChiRPA_AB_q_omega_NBZ[q_index][omega_ind] += A_elmt*B_elmt*(1.0/(NSites_in_MUC))*exp(iota_complex*2.0*PI*( (1.0*q1_ind*(g1-h1))/(lx_) + (1.0*q2_ind*(g2-h2))/(ly_) ))*
                                                                                    ChiRPAMat[qp_ind][omega_ind][ind1][ind2][ind3][ind4];

                                                                        }
                                                                    }


                                                                    //file_out_BarreSusc<<q_ind_temp<<"  "<<q1_ind<<"  "<<q2_ind<<"  "<<omega<<"  ";
                                                                    //file_out_BarreSusc<<ChiRPA_AB_q_omega_NBZ[q_index][omega_ind].real()<<"  "<<ChiRPA_AB_q_omega_NBZ[q_index][omega_ind].imag()<<endl;


                                                                }
                                                                //file_out_BarreSusc<<endl;

                                                            }







                                                        }

                                                    }}}
                                        }}}
                            }}}
                }}}



        for(int q_ind_temp=0;q_ind_temp<q_path.size();q_ind_temp++){
            //cout<<"ChiRPA_q_omega (NBZ) "<<q_ind_temp<<"("<<q_path.size()<<")"<<endl;
            q1_ind=q_path[q_ind_temp].first;
            q2_ind=q_path[q_ind_temp].second;
            int q_index = q1_ind + q2_ind*lx_;

            for(int omega_ind=0;omega_ind<N_omega;omega_ind++){
                double omega = omega_ind*d_omega + omega_min;

                file_out_BarreSusc<<q_ind_temp<<"  "<<q1_ind<<"  "<<q2_ind<<"  "<<omega<<"  ";
                file_out_BarreSusc<<ChiRPA_AB_q_omega_NBZ[q_index][omega_ind].real()<<"  "<<ChiRPA_AB_q_omega_NBZ[q_index][omega_ind].imag()<<endl;
            }
            file_out_BarreSusc<<endl;
        }



    }

}



void Kspace_calculation_G2dLatticeNew::Get_RPA_Susceptibility_Matrix(){


    double EPS_ZERO=0.0001;

    int S_ = NSites_in_MUC;

    Mat_1_intpair q_path;
    //Create_K_Path("FullBZ", q_path_plot);
    Create_K_Path("FullBZ", q_path);

    int q1_ind, q2_ind;

    double omega_min=Parameters_.omega_min;
    double omega_max=Parameters_.omega_max;
    double d_omega=Parameters_.d_omega;
    double eta=Parameters_.eta;

    int N_omega = int(((omega_max-omega_min)/d_omega) + 0.5) + 1;

    //ind = uc_no + (atom + orb*n_atoms)*(S_) +  sigma*(n_atoms_*n_orbs_*S_)
    ChiBareMat.resize(lx_cells*ly_cells);
    for(int q_ind=0;q_ind<lx_cells*ly_cells;q_ind++){
        ChiBareMat[q_ind].resize(N_omega);
        for(int omega_ind=0;omega_ind<N_omega;omega_ind++){
            ChiBareMat[q_ind][omega_ind].resize(S_*n_atoms_*n_orbs_*2);
            for(int ind1=0;ind1<S_*n_atoms_*n_orbs_*2;ind1++){
                ChiBareMat[q_ind][omega_ind][ind1].resize(S_*n_atoms_*n_orbs_*2);
                for(int ind2=0;ind2<S_*n_atoms_*n_orbs_*2;ind2++){
                    ChiBareMat[q_ind][omega_ind][ind1][ind2].resize(S_*n_atoms_*n_orbs_*2);
                    for(int ind3=0;ind3<S_*n_atoms_*n_orbs_*2;ind3++){
                        ChiBareMat[q_ind][omega_ind][ind1][ind2][ind3].resize(S_*n_atoms_*n_orbs_*2);
                    }
                }
            }
        }
    }


    ChiRPAMat.resize(lx_cells*ly_cells);
    for(int q_ind=0;q_ind<lx_cells*ly_cells;q_ind++){
        ChiRPAMat[q_ind].resize(N_omega);
        for(int omega_ind=0;omega_ind<N_omega;omega_ind++){
            ChiRPAMat[q_ind][omega_ind].resize(S_*n_atoms_*n_orbs_*2);
            for(int ind1=0;ind1<S_*n_atoms_*n_orbs_*2;ind1++){
                ChiRPAMat[q_ind][omega_ind][ind1].resize(S_*n_atoms_*n_orbs_*2);
                for(int ind2=0;ind2<S_*n_atoms_*n_orbs_*2;ind2++){
                    ChiRPAMat[q_ind][omega_ind][ind1][ind2].resize(S_*n_atoms_*n_orbs_*2);
                    for(int ind3=0;ind3<S_*n_atoms_*n_orbs_*2;ind3++){
                        ChiRPAMat[q_ind][omega_ind][ind1][ind2][ind3].resize(S_*n_atoms_*n_orbs_*2);
                    }
                }
            }
        }
    }



    Matrix<complex<double>> IdentityMat;
    IdentityMat.resize(S_*n_atoms_*n_orbs_*2*S_*n_atoms_*n_orbs_*2, S_*n_atoms_*n_orbs_*2*S_*n_atoms_*n_orbs_*2);
    for(int i=0;i<S_*n_atoms_*n_orbs_*2*S_*n_atoms_*n_orbs_*2;i++){
        IdentityMat(i,i)=1.0;
    }


    int ChiMatsize = S_*n_atoms_*n_orbs_*2;

    ChiBareTimesI.resize(lx_cells*ly_cells);
    OneMinusChiBareTimesI_Inv.resize(lx_cells*ly_cells);
    for(int q_ind=0;q_ind<lx_cells*ly_cells;q_ind++){
        ChiBareTimesI[q_ind].resize(N_omega);
        OneMinusChiBareTimesI_Inv[q_ind].resize(N_omega);
        for(int omega_ind=0;omega_ind<N_omega;omega_ind++){
            OneMinusChiBareTimesI_Inv[q_ind][omega_ind].resize(S_*n_atoms_*n_orbs_*2*S_*n_atoms_*n_orbs_*2, S_*n_atoms_*n_orbs_*2*S_*n_atoms_*n_orbs_*2 );
            OneMinusChiBareTimesI_Inv[q_ind][omega_ind] = IdentityMat;
            ChiBareTimesI[q_ind][omega_ind].resize(S_*n_atoms_*n_orbs_*2);
            for(int ind1=0;ind1<S_*n_atoms_*n_orbs_*2;ind1++){
                ChiBareTimesI[q_ind][omega_ind][ind1].resize(S_*n_atoms_*n_orbs_*2);
                for(int ind2=0;ind2<S_*n_atoms_*n_orbs_*2;ind2++){
                    ChiBareTimesI[q_ind][omega_ind][ind1][ind2].resize(S_*n_atoms_*n_orbs_*2);
                    for(int ind3=0;ind3<S_*n_atoms_*n_orbs_*2;ind3++){
                        ChiBareTimesI[q_ind][omega_ind][ind1][ind2][ind3].resize(S_*n_atoms_*n_orbs_*2);
                    }
                }
            }
        }
    }



#ifdef _OPENMP
    int N_p = omp_get_max_threads();
    omp_set_num_threads(Parameters_.NProcessors);

    cout<<"Max threads which can be used parallely = "<<N_p<<endl;
    cout<<"No. of threads used parallely = "<<Parameters_.NProcessors<<endl;
#endif



    cout<<"CALCULATING BARE CHI"<<endl;
    //Calculate ChiBare

    int thread_id=0;
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
    for(int q_ind_temp=0;q_ind_temp<q_path.size();q_ind_temp++){
#ifdef _OPENMP
        thread_id = omp_get_thread_num();
#endif
        cout<<"q_ind = "<<q_ind_temp<<" ("<<q_path.size()<<") in thread "<<thread_id<<endl;

        int q1_ind_local=q_path[q_ind_temp].first;
        int q2_ind_local=q_path[q_ind_temp].second;
        int q_index = Coordinates_.Ncell(q1_ind_local, q2_ind_local);

        for(int g=0;g<S_;g++){
            for(int m=0;m<n_atoms_;m++){
                for(int alpha=0;alpha<n_orbs_;alpha++){
                    for(int sigma=0;sigma<2;sigma++){
                        int ind1 = g + (m + alpha*n_atoms_)*S_ + sigma*(n_atoms_*n_orbs_*S_);

                        int h=g;
                        int n=m;
                        //  for(int h=0;h<S_;h++){??
                        //    for(int n=0;n<n_atoms_;n++){
                        for(int beta=0;beta<n_orbs_;beta++){
                            for(int sigma_til=0;sigma_til<2;sigma_til++){
                                int ind2 = h + (n + beta*n_atoms_)*S_ + sigma_til*(n_atoms_*n_orbs_*S_);


                                for(int hp=0;hp<S_;hp++){
                                    for(int np=0;np<n_atoms_;np++){
                                        for(int alphap=0;alphap<n_orbs_;alphap++){
                                            for(int sigmap=0;sigmap<2;sigmap++){
                                                int ind3 = hp + (np + alphap*n_atoms_)*S_ + sigmap*(n_atoms_*n_orbs_*S_);



                                                for(int betap=0;betap<n_orbs_;betap++){
                                                    for(int sigma_tilp=0;sigma_tilp<2;sigma_tilp++){
                                                        int ind4 = hp + (np + betap*n_atoms_)*S_ + sigma_tilp*(n_atoms_*n_orbs_*S_);

                                                        cout<<"BARE CHI : "<<ind1<<" ("<<S_*n_atoms_*n_orbs_*2<<") "<<ind2<<" ("<<S_*n_atoms_*n_orbs_*2<<") "<<ind3<<" ("<<S_*n_atoms_*n_orbs_*2<<") "<<betap + sigma_tilp*n_orbs_<<" ("<<n_orbs_*2<<") "<<q_ind_temp<<" ("<<q_path.size()<<") "<<endl;


                                                        for(int omega_ind=0;omega_ind<N_omega;omega_ind++){
                                                            double omega = omega_ind*d_omega + omega_min;
                                                            //   file_out_BarreSusc<<q_ind_temp<<"  "<<q1_ind<<"  "<<q2_ind<<"  "<<omega<<"  ";

                                                            ChiBareMat[q_index][omega_ind][ind1][ind2][ind3][ind4]=0.0;

                                                            for(int k1_ind=0;k1_ind<lx_cells;k1_ind++){
                                                                for(int k2_ind=0;k2_ind<ly_cells;k2_ind++){

                                                                    int kpq1_ind = (k1_ind + q1_ind_local + lx_cells)%lx_cells;
                                                                    int kpq2_ind = (k2_ind + q2_ind_local + ly_cells)%ly_cells;

                                                                    int k_index = Coordinates_.Ncell(k1_ind,k2_ind);
                                                                    int kpq_index = Coordinates_.Ncell(kpq1_ind,kpq2_ind);





                                                                    for(int lambda=0;lambda<2*n_atoms_*n_orbs_*S_;lambda++){
                                                                        for(int lambdap=0;lambdap<2*n_atoms_*n_orbs_*S_;lambdap++){
                                                                            int state_k_lambda = 2*n_atoms_*n_orbs_*S_*k_index + lambda;
                                                                            int state_kpq_lambdap = 2*n_atoms_*n_orbs_*S_*kpq_index + lambdap;

                                                                            double a_temp=(omega + Eigenvalues_saved[state_k_lambda] - Eigenvalues_saved[state_kpq_lambdap]);
                                                                            double absolute_val_sqr= a_temp*a_temp + eta*eta;
                                                                            double real_part = 1.0*a_temp/absolute_val_sqr;
                                                                            double imag_part = -1.0*eta/absolute_val_sqr;

                                                                            complex<double> inverse_pole = complex<double>(real_part, imag_part);

                                                                            ChiBareMat[q_index][omega_ind][ind1][ind2][ind3][ind4] += ((1.0/(lx_cells*ly_cells)))*((1.0/( exp((Eigenvalues_saved[state_k_lambda]-mu_)*Parameters_.beta ) + 1.0)) - (1.0/( exp((Eigenvalues_saved[state_kpq_lambdap]-mu_)*Parameters_.beta ) + 1.0)))
                                                                                    *inverse_pole
                                                                                    *conj(Eigvectors_[state_k_lambda][ind1])*Eigvectors_[state_kpq_lambdap][ind2]
                                                                                    *conj(Eigvectors_[state_kpq_lambdap][ind3])*(Eigvectors_[state_k_lambda][ind4]);



                                                                            //ChiBare_AB_q_omega[IntraMUC_indexA][IntraMUC_indexB][q_index][omega_ind] += ((1.0/(lx_cells*ly_cells)))*((1.0/( exp((Eigenvalues_saved[state_k_n]-mu_)*Parameters_.beta ) + 1.0)) - (1.0/( exp((Eigenvalues_saved[state_kpq_m]-mu_)*Parameters_.beta ) + 1.0)))*
                                                                            //        ABmat[IntraMUC_indexA][IntraMUC_indexB][state_k_n][state_kpq_m]*inverse_pole;

                                                                        }}
                                                                }}
                                                        }
                                                        //}
                                                    }
                                                }

                                            }
                                        }
                                    }
                                }

                            }
                        }
                        //  }
                        //}



                    }
                }
            }
        }
    }


    cout<<"CALCULATING OneMinusChiBareTimesI_Inv"<<endl;
    //Calculate ChiBareTimesI and OneMinusChiBareTimesI_Inv
    for(int q_ind_temp=0;q_ind_temp<q_path.size();q_ind_temp++){
        cout<<"OneMinusChiBareTimesI_Inv   "<<q_ind_temp<<" ("<<q_path.size()<<") "<<endl;

        q1_ind=q_path[q_ind_temp].first;
        q2_ind=q_path[q_ind_temp].second;
        int q_index = Coordinates_.Ncell(q1_ind, q2_ind);

        for(int omega_ind=0;omega_ind<N_omega;omega_ind++){

            for(int g=0;g<S_;g++){
                for(int m=0;m<n_atoms_;m++){
                    for(int alpha=0;alpha<n_orbs_;alpha++){
                        for(int sigma=0;sigma<2;sigma++){
                            int ind1 = g + (m + alpha*n_atoms_)*S_ + sigma*(n_atoms_*n_orbs_*S_);


                            int h=g;
                            int n=m;
                            // for(int h=0;h<S_;h++){
                            //     for(int n=0;n<n_atoms_;n++){
                            for(int beta=0;beta<n_orbs_;beta++){
                                for(int sigma_til=0;sigma_til<2;sigma_til++){
                                    int ind2 = h + (n + beta*n_atoms_)*S_ + sigma_til*(n_atoms_*n_orbs_*S_);

                                    for(int hp=0;hp<S_;hp++){
                                        for(int np=0;np<n_atoms_;np++){
                                            for(int alphap=0;alphap<n_orbs_;alphap++){
                                                for(int sigmap=0;sigmap<2;sigmap++){
                                                    int ind3 = hp + (np + alphap*n_atoms_)*S_ + sigmap*(n_atoms_*n_orbs_*S_);


                                                    for(int betap=0;betap<n_orbs_;betap++){
                                                        for(int sigma_tilp=0;sigma_tilp<2;sigma_tilp++){
                                                            int ind4 = hp + (np + betap*n_atoms_)*S_ + sigma_tilp*(n_atoms_*n_orbs_*S_);

                                                            ChiBareTimesI[q_index][omega_ind][ind1][ind2][ind3][ind4]=0.0;

                                                            for(int alphapp=0;alphapp<n_orbs_;alphapp++){
                                                                for(int sigmapp=0;sigmapp<2;sigmapp++){
                                                                    int ind3_p = hp + (np + alphapp*n_atoms_)*S_ + sigmapp*(n_atoms_*n_orbs_*S_);

                                                                    for(int betapp=0;betapp<n_orbs_;betapp++){
                                                                        for(int sigma_tilpp=0;sigma_tilpp<2;sigma_tilpp++){
                                                                            int ind4_p = hp + (np + betapp*n_atoms_)*S_ + sigma_tilpp*(n_atoms_*n_orbs_*S_);

                                                                            ChiBareTimesI[q_index][omega_ind][ind1][ind2][ind3][ind4] += ChiBareMat[q_index][omega_ind][ind1][ind2][ind3_p][ind4_p]
                                                                                    *IntKer[np][alphapp+sigmapp*n_orbs_][betapp+sigma_tilpp*n_orbs_][alphap+sigmap*n_orbs_][betap+sigma_tilp*n_orbs_];

                                                                        }
                                                                    }
                                                                }
                                                            }

                                                            OneMinusChiBareTimesI_Inv[q_index][omega_ind](ind1 + ind2*ChiMatsize,ind3 + ind4*ChiMatsize) -= ChiBareTimesI[q_index][omega_ind][ind1][ind2][ind3][ind4];

                                                        }
                                                    }


                                                }
                                            }
                                        }
                                    }


                                }
                            }
                            // }
                            // }



                        }
                    }
                }
            }


            Inverse(OneMinusChiBareTimesI_Inv[q_index][omega_ind]);

        }}





    //Chi_RPA_Mat
    for(int q_ind_temp=0;q_ind_temp<q_path.size();q_ind_temp++){
        q1_ind=q_path[q_ind_temp].first;
        q2_ind=q_path[q_ind_temp].second;
        int q_index = Coordinates_.Ncell(q1_ind, q2_ind);

        for(int omega_ind=0;omega_ind<N_omega;omega_ind++){

            for(int g=0;g<S_;g++){
                for(int m=0;m<n_atoms_;m++){
                    for(int alpha=0;alpha<n_orbs_;alpha++){
                        for(int sigma=0;sigma<2;sigma++){
                            int ind1 = g + (m + alpha*n_atoms_)*S_ + sigma*(n_atoms_*n_orbs_*S_);

                            for(int beta=0;beta<n_orbs_;beta++){
                                for(int sigma_til=0;sigma_til<2;sigma_til++){
                                    int ind2 = g + (m + beta*n_atoms_)*S_ + sigma_til*(n_atoms_*n_orbs_*S_);

                                    for(int h=0;h<S_;h++){
                                        for(int n=0;n<n_atoms_;n++){
                                            for(int alphap=0;alphap<n_orbs_;alphap++){
                                                for(int sigmap=0;sigmap<2;sigmap++){
                                                    int ind3 = h + (n + alphap*n_atoms_)*S_ + sigmap*(n_atoms_*n_orbs_*S_);


                                                    for(int betap=0;betap<n_orbs_;betap++){
                                                        for(int sigma_tilp=0;sigma_tilp<2;sigma_tilp++){
                                                            int ind4 = h + (n + betap*n_atoms_)*S_ + sigma_tilp*(n_atoms_*n_orbs_*S_);

                                                            ChiRPAMat[q_index][omega_ind][ind1][ind2][ind3][ind4]=0.0;


                                                            for(int hp=0;hp<S_;hp++){
                                                                for(int np=0;np<n_atoms_;np++){
                                                                    for(int alphapp=0;alphapp<n_orbs_;alphapp++){
                                                                        for(int sigmapp=0;sigmapp<2;sigmapp++){
                                                                            int ind3_p = hp + (np + alphapp*n_atoms_)*S_ + sigmapp*(n_atoms_*n_orbs_*S_);

                                                                            for(int gp=0;gp<S_;gp++){
                                                                                for(int mp=0;mp<n_atoms_;mp++){
                                                                                    for(int betapp=0;betapp<n_orbs_;betapp++){
                                                                                        for(int sigma_tilpp=0;sigma_tilpp<2;sigma_tilpp++){
                                                                                            int ind4_p = gp + (mp + betapp*n_atoms_)*S_ + sigma_tilpp*(n_atoms_*n_orbs_*S_);


                                                                                            ChiRPAMat[q_index][omega_ind][ind1][ind2][ind3][ind4] += OneMinusChiBareTimesI_Inv[q_index][omega_ind](ind1 + ChiMatsize*ind2, ind3_p + ChiMatsize*ind4_p)
                                                                                                    *ChiBareMat[q_index][omega_ind][ind3_p][ind4_p][ind3][ind4];

                                                                                        }}}}
                                                                        }}}}
                                                        }}
                                                }}}}
                                }}
                        }}}}
        }}



}




void Kspace_calculation_G2dLatticeNew::Get_Bare_Susceptibility_New(){



    //Only for pair_no=0
    Susc_OprA=Parameters_.Susc_OprA[0];
    Susc_OprB=Parameters_.Susc_OprB[0];


    double EPS_ZERO=0.0001;

    int S_ = NSites_in_MUC;

    Mat_1_intpair q_path, q_path_plot;
    Create_K_Path("GMXGYMG", q_path_plot);

    //Create_K_Path("FullBZ", q_path_plot);
    Create_K_Path("FullBZ", q_path);

    int q1_ind, q2_ind;


    double omega_min=Parameters_.omega_min;
    double omega_max=Parameters_.omega_max;
    double d_omega=Parameters_.d_omega;
    double eta=Parameters_.eta;

    int N_omega = int(((omega_max-omega_min)/d_omega) + 0.5) + 1;



    ChiBare_AB_q_omega.resize(NSites_in_MUC);
    for(int g=0;g<NSites_in_MUC;g++){
        ChiBare_AB_q_omega[g].resize(NSites_in_MUC);
        for(int h=0;h<NSites_in_MUC;h++){
            ChiBare_AB_q_omega[g][h].resize(lx_cells*ly_cells);
            for(int q_ind=0;q_ind<lx_cells*ly_cells;q_ind++){
                ChiBare_AB_q_omega[g][h][q_ind].resize(N_omega);
            }
        }
    }



    Mat_4_Complex_doub ABmat;
    ABmat.resize(NSites_in_MUC);
    for(int g=0;g<NSites_in_MUC;g++){
        ABmat[g].resize(NSites_in_MUC);
        for(int h=0;h<NSites_in_MUC;h++){
            ABmat[g][h].resize(Eigenvalues_saved.size());
            for(int state_n=0;state_n<Eigenvalues_saved.size();state_n++){
                ABmat[g][h][state_n].resize(Eigenvalues_saved.size());
            }
        }
    }



    //Calculate ABmat
    for(int IntraMUC_indexA=0;IntraMUC_indexA<S_;IntraMUC_indexA++){
        for(int IntraMUC_indexB=0;IntraMUC_indexB<S_;IntraMUC_indexB++){

            for(int atom_no_alpha=0;atom_no_alpha<n_atoms_;atom_no_alpha++){
                for(int alpha_orb=0;alpha_orb<n_orbs_;alpha_orb++){
                    for(int spin_alpha=0;spin_alpha<2;spin_alpha++){

                        //alpha
                        int alpha_comp = IntraMUC_indexA + (atom_no_alpha + n_atoms_*alpha_orb)*(S_) +
                                spin_alpha*(n_atoms_*n_orbs_*S_);


                        for(int atom_no_beta=0;atom_no_beta<n_atoms_;atom_no_beta++){
                            for(int beta_orb=0;beta_orb<n_orbs_;beta_orb++){
                                for(int spin_beta=0;spin_beta<2;spin_beta++){

                                    //beta
                                    int beta_comp = IntraMUC_indexA + (atom_no_beta + n_atoms_*beta_orb)*(S_) +
                                            spin_beta*(n_atoms_*n_orbs_*S_);


                                    for(int atom_no_alpha_p=0;atom_no_alpha_p<n_atoms_;atom_no_alpha_p++){
                                        for(int alpha_orb_p=0;alpha_orb_p<n_orbs_;alpha_orb_p++){
                                            for(int spin_alpha_p=0;spin_alpha_p<2;spin_alpha_p++){

                                                //alpha_p
                                                int alpha_comp_p = IntraMUC_indexB + (atom_no_alpha_p + n_atoms_*alpha_orb_p)*(S_) +
                                                        spin_alpha_p*(n_atoms_*n_orbs_*S_);


                                                for(int atom_no_beta_p=0;atom_no_beta_p<n_atoms_;atom_no_beta_p++){
                                                    for(int beta_orb_p=0;beta_orb_p<n_orbs_;beta_orb_p++){
                                                        for(int spin_beta_p=0;spin_beta_p<2;spin_beta_p++){

                                                            //beta_p
                                                            int beta_comp_p = IntraMUC_indexB + (atom_no_beta_p + n_atoms_*beta_orb_p)*(S_) +
                                                                    spin_beta_p*(n_atoms_*n_orbs_*S_);

                                                            complex<double> A_elmt = Susc_OprA[atom_no_alpha + n_atoms_*alpha_orb + n_atoms_*n_orbs_*spin_alpha][atom_no_beta + n_atoms_*beta_orb + n_atoms_*n_orbs_*spin_beta];
                                                            complex<double> B_elmt = Susc_OprB[atom_no_alpha_p + n_atoms_*alpha_orb_p + n_atoms_*n_orbs_*spin_alpha_p][atom_no_beta_p + n_atoms_*beta_orb_p + n_atoms_*n_orbs_*spin_beta_p];


                                                            if(abs(A_elmt*B_elmt)>EPS_ZERO){

                                                                for(int n=0;n<2*n_atoms_*n_orbs_*S_;n++){
                                                                    for(int k1_ind_n=0;k1_ind_n<lx_cells;k1_ind_n++){
                                                                        for(int k2_ind_n=0;k2_ind_n<ly_cells;k2_ind_n++){
                                                                            int k_index_n = Coordinates_.Ncell(k1_ind_n,k2_ind_n);
                                                                            int state_k_n = 2*n_atoms_*n_orbs_*S_*k_index_n + n;

                                                                            for(int m=0;m<2*n_atoms_*n_orbs_*S_;m++){
                                                                                for(int k1_ind_m=0;k1_ind_m<lx_cells;k1_ind_m++){
                                                                                    for(int k2_ind_m=0;k2_ind_m<ly_cells;k2_ind_m++){
                                                                                        int k_index_m = Coordinates_.Ncell(k1_ind_m,k2_ind_m);
                                                                                        int state_k_m = 2*n_atoms_*n_orbs_*S_*k_index_m + m;

                                                                                        ABmat[IntraMUC_indexA][IntraMUC_indexB][state_k_n][state_k_m] += A_elmt*B_elmt
                                                                                                *conj(Eigvectors_[state_k_n][alpha_comp])*Eigvectors_[state_k_m][beta_comp]
                                                                                                *conj(Eigvectors_[state_k_m][alpha_comp_p])*(Eigvectors_[state_k_n][beta_comp_p]);
                                                                                        //*conj(Eigvectors_[state_k_n][beta_comp])*Eigvectors_[state_k_n][alpha_comp_p];

                                                                                    }}}

                                                                        }}}

                                                            }

                                                        }}}
                                            }}}
                                }}}
                    }}}

        }}



    //calculating Chi_bare_q_w
    for(int IntraMUC_indexA=0;IntraMUC_indexA<S_;IntraMUC_indexA++){
        for(int IntraMUC_indexB=0;IntraMUC_indexB<S_;IntraMUC_indexB++){

            string File_Out_BarreSusc_str = "Bare_Susc_gh_" +to_string(IntraMUC_indexA) + "_" + to_string(IntraMUC_indexB) + "_pair_no" +to_string(0)+".txt";
            ofstream file_out_BarreSusc(File_Out_BarreSusc_str.c_str());
            file_out_BarreSusc<<"#q omega Chi_bare_gh(q,omega)"<<endl;



            cout<<"(g,h) = "<<IntraMUC_indexA<<" , "<<IntraMUC_indexB<<" ----------------"<<endl;
            for(int q_ind_temp=0;q_ind_temp<q_path.size();q_ind_temp++){
                cout<<"Chi_q_omega "<<q_ind_temp<<"("<<q_path.size()<<")"<<endl;

                q1_ind=q_path[q_ind_temp].first;
                q2_ind=q_path[q_ind_temp].second;
                int q_index = Coordinates_.Ncell(q1_ind, q2_ind);

                for(int omega_ind=0;omega_ind<N_omega;omega_ind++){
                    double omega = omega_ind*d_omega + omega_min;
                    //   file_out_BarreSusc<<q_ind_temp<<"  "<<q1_ind<<"  "<<q2_ind<<"  "<<omega<<"  ";

                    ChiBare_AB_q_omega[IntraMUC_indexA][IntraMUC_indexB][q_index][omega_ind]=0.0;
                    for(int k1_ind=0;k1_ind<lx_cells;k1_ind++){
                        for(int k2_ind=0;k2_ind<ly_cells;k2_ind++){

                            int kpq1_ind = (k1_ind + q1_ind + lx_cells)%lx_cells;
                            int kpq2_ind = (k2_ind + q2_ind + ly_cells)%ly_cells;

                            int k_index = Coordinates_.Ncell(k1_ind,k2_ind);
                            int kpq_index = Coordinates_.Ncell(kpq1_ind,kpq2_ind);

                            for(int n=0;n<2*n_atoms_*n_orbs_*S_;n++){
                                for(int m=0;m<2*n_atoms_*n_orbs_*S_;m++){
                                    int state_k_n = 2*n_atoms_*n_orbs_*S_*k_index + n;
                                    int state_kpq_m = 2*n_atoms_*n_orbs_*S_*kpq_index + m;

                                    double a_temp=(omega + Eigenvalues_saved[state_k_n] - Eigenvalues_saved[state_kpq_m]);
                                    double absolute_val_sqr= a_temp*a_temp + eta*eta;
                                    double real_part = 1.0*a_temp/absolute_val_sqr;
                                    double imag_part = -1.0*eta/absolute_val_sqr;

                                    complex<double> inverse_pole = complex<double>(real_part, imag_part);

                                    ChiBare_AB_q_omega[IntraMUC_indexA][IntraMUC_indexB][q_index][omega_ind] += ((1.0/(lx_cells*ly_cells)))*((1.0/( exp((Eigenvalues_saved[state_k_n]-mu_)*Parameters_.beta ) + 1.0)) - (1.0/( exp((Eigenvalues_saved[state_kpq_m]-mu_)*Parameters_.beta ) + 1.0)))*
                                            ABmat[IntraMUC_indexA][IntraMUC_indexB][state_k_n][state_kpq_m]*inverse_pole;

                                }}
                        }}

                    //file_out_BarreSusc<<ChiBare_AB_q_omega[IntraMUC_indexA][IntraMUC_indexB][q_index][omega_ind].real()<<"  "<<ChiBare_AB_q_omega[IntraMUC_indexA][IntraMUC_indexB][q_index][omega_ind].imag()<<endl;
                }
                //file_out_BarreSusc<<endl;

            }


            for(int q_ind_temp=0;q_ind_temp<q_path_plot.size();q_ind_temp++){
                int q1_ind_plot=q_path_plot[q_ind_temp].first;
                int q2_ind_plot=q_path_plot[q_ind_temp].second;
                int q_index = Coordinates_.Ncell(q1_ind_plot, q2_ind_plot);

                for(int omega_ind=0;omega_ind<N_omega;omega_ind++){
                    double omega = omega_ind*d_omega + omega_min;
                    file_out_BarreSusc<<q_ind_temp<<"  "<<q1_ind_plot<<"  "<<q2_ind_plot<<"  "<<omega<<"  ";
                    file_out_BarreSusc<<ChiBare_AB_q_omega[IntraMUC_indexA][IntraMUC_indexB][q_index][omega_ind].real()<<"  "<<ChiBare_AB_q_omega[IntraMUC_indexA][IntraMUC_indexB][q_index][omega_ind].imag()<<endl;
                }
                file_out_BarreSusc<<endl;
            }


        }}





}




void Kspace_calculation_G2dLatticeNew::SelfConsistency(){


    string File_Out_Local_OP_Initial = "Initial_" + Parameters_.File_OPs_out + ".txt";
    ofstream file_out_Local_OP_Initial(File_Out_Local_OP_Initial.c_str());
    file_out_Local_OP_Initial<<"#row col OParams_[row][col]"<<endl;

    for(int alpha=0;alpha<OPs_.value.size();alpha++){
        file_out_Local_OP_Initial<<OPs_.rows[alpha]<<setw(15)<<OPs_.columns[alpha]<<setw(15)<<"  "<<OPs_.value[alpha]<<endl;
    }


    //assert(false);

    kick_while_cooling=0.0;
    for(int Temp_no=0;Temp_no<Parameters_.Temperature_points.size();Temp_no++){

        Parameters_.Temperature = Parameters_.Temperature_points[Temp_no];
        Parameters_.beta = 1.0/Parameters_.Temperature;

        cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
        cout<<"Temperature = "<<Parameters_.Temperature<<endl;

        char temp_char[50];
        sprintf(temp_char, "%.10f", Parameters_.Temperature);


        string File_Out_progress;
        File_Out_progress = "output_Kspace_SelfConsistency_Temp" + string(temp_char) + ".txt";
        ofstream file_out_progress(File_Out_progress.c_str());

        file_out_progress<<"# iter   OP_error_   E_class   E_quant    mu_   Total_Particles   OPs_.total_Den" <<endl;

        cout<<"error targetted = "<<Parameters_.Convergence_Error<<endl;
        cout<<"Max iterations = "<<Parameters_.IterMax<<endl;

        OP_error_=10.0;
        int iter=0;
        while( (OP_error_>=Parameters_.Convergence_Error) && (iter<Parameters_.IterMax)){

            Create_Kspace_Spectrum();
            Arranging_spectrum();

            //            for(int ei=0;ei<Eigenvalues_.size();ei++){
            //                cout<<ei<<"  "<<Eigenvalues_[ei]<<endl;
            //            }

            if(Parameters_.Fixing_mu){
                mu_ = Parameters_.Fixed_mu;
                Update_Total_Density();
            }
            else{
                mu_=chemicalpotential(Parameters_.Total_Particles);
            }

            Get_new_OPs_and_error();
            // Get_Energies();
            Get_Energies_new();


            file_out_progress<<setprecision(15)<<iter<<"   "<<OP_error_<<"   "<<E_class<<"   "<<E_quant<<"    "<<mu_<<"   "<<Parameters_.Total_Particles <<"   "<< OPs_total_den <<endl;
            //        for(int OP_no=0;OP_no<6;OP_no++){
            //            file_out_progress<<OPs_[OP_no].real()<<"    "<<OPs_[OP_no].imag()<<"    ";
            //        }
            //        file_out_progress<<endl;

            if(Parameters_.Anderson_Mixing){
                Update_OrderParameters_AndersonMixing(iter);
                //cout<<"using Anderson mixing"<<endl;
            }
            else{ //SimpleMixing
                for(int i=0;i<OPs_.value.size();i++){
                    OPs_.value[i]=Parameters_.alpha_OP*OPs_.value[i] + (1.0-Parameters_.alpha_OP)*OPs_new_.value[i];
                }
            }


            iter++;

        }




        Create_Kspace_Spectrum();
        if(lx_cells==ly_cells){
            Get_Bands();
        }
        //Calculate_ChernNumbers();

        Arranging_spectrum();
        if(Parameters_.Fixing_mu){
            mu_ = Parameters_.Fixed_mu;
            Update_Total_Density();
        }
        else{
            mu_=chemicalpotential(Parameters_.Total_Particles);
        }

        //Create_Current_Oprs_Faster();
        //        J_KE_e1.print();
        //        cout<<"============ PRINT EIGVECS ============================="<<endl;
        //        for(int r_=0;r_<Eigvectors_.size();r_++){
        //            for(int c_=0;c_<Eigvectors_[r_].size();c_++){
        //                cout<<Eigvectors_[r_][c_]<<" ";
        //            }
        //            cout<<endl;
        //        }
        //        cout<<"=========================================================="<<endl;
        //Create_Current_Oprs();
        //Hall_conductance();

        cout<<"mu = "<<mu_<<endl;
        cout<<"energies shown below , Eclass and Equant:"<<endl;
        Get_new_OPs_and_error();
        //Get_Energies();
        Get_Energies_new();
        cout<<"Eclass = "<<E_class<<" , Equant = "<<E_quant<<"  "<<endl;
        cout<<"ETotal = "<<E_class+E_quant<<endl;


        // cout<<"================== Classical Energies by parts ====================="<<endl;
        // cout<<"E_class_onsite_U0_Hartree = "<<E_class_onsite_U0_Hartree <<endl;
        // cout<<"E_class_longrange_Hartree = "<<E_class_longrange_Hartree<<endl;
        // cout<<"E_class_onsite_U0_Fock = "<<E_class_onsite_U0_Fock<<endl;
        // cout<<"E_class_longrange_Fock = "<<E_class_longrange_Fock<<endl;
        // cout<<"===================================================================="<<endl;

        Get_spin_resolved_local_densities();
        Get_local_spins();

        if(n_orbs_>1){
            Get_Tau_Pseudospins();
        }
        Calculate_Nw();
        //Calculate_Akw();

        string Eigenvalues_fl_out = "Eigenvalues.txt";
        ofstream Eigenvalues_file_out(Eigenvalues_fl_out.c_str());
        for(int ie=0;ie<Eigenvalues_.size();ie++){
            Eigenvalues_file_out<<ie<<"  "<<Eigenvalues_[ie]<<endl;
        }







        string File_Out_Local_OP = Parameters_.File_OPs_out + "_Temp"+string(temp_char) + ".txt";
        ofstream file_out_Local_OP(File_Out_Local_OP.c_str());
        file_out_Local_OP<<"#unitcell_row  atom_row orb sigma row_val unitcella_col tom_col orb sigma col_val OParams_[row][col]"<<endl;
        file_out_Local_OP<<lx_<<  "  "<<ly_<<endl;

        int S_= NSites_in_MUC;
        int row_temp, col_temp;
        int cell_row, alpha_plus_gamma, alpha_row,gamma_row, sigma_row, c1;
        int alpha_plus_gamma_sigma, alpha_col, gamma_col, sigma_col, cell_col;
        int atom_no_row, atom_no_col, orb_no_row, orb_no_col;
        for(int OP_no=0;OP_no<OPs_.value.size();OP_no++){
            row_temp=OPs_.rows[OP_no];
            col_temp=OPs_.columns[OP_no];


            //alpha + gamma*(S_) +  sigma*(n_atoms_*n_orbs_*S_) + cell_*(2*n_atoms_*n_orbs_*S_)
            cell_row = 0;
            alpha_plus_gamma = row_temp%(n_atoms_*n_orbs_*S_);
            alpha_row = alpha_plus_gamma%S_;
            gamma_row = (alpha_plus_gamma -alpha_row)/S_;
            //gamma = atom_no + n_Atoms*orb_no
            atom_no_row= gamma_row%n_atoms_;
            orb_no_row= (gamma_row - atom_no_row)/n_atoms_;
            sigma_row = (row_temp - alpha_plus_gamma)/(n_atoms_*n_orbs_*S_);
            c1 = alpha_row + gamma_row*(S_) +  sigma_row*(n_atoms_*n_orbs_*S_);
            assert(c1==row_temp);

            alpha_plus_gamma_sigma = col_temp%(2*n_atoms_*n_orbs_*S_);
            alpha_plus_gamma = alpha_plus_gamma_sigma%(n_atoms_*n_orbs_*S_);
            alpha_col = alpha_plus_gamma%S_;
            gamma_col = (alpha_plus_gamma -alpha_col)/S_;
            atom_no_col= gamma_col%n_atoms_;
            orb_no_col= (gamma_col - atom_no_col)/n_atoms_;
            sigma_col = (alpha_plus_gamma_sigma - alpha_plus_gamma)/(n_atoms_*n_orbs_*S_);
            cell_col = (col_temp - alpha_plus_gamma_sigma)/(2*n_atoms_*n_orbs_*S_);

            file_out_Local_OP<<alpha_row<<"  "<<atom_no_row<<"  "<<orb_no_row<<"  "<<sigma_row<<"  "
                            <<setw(15)<<alpha_col<<"  "<<atom_no_col<<"  "<<orb_no_col<<"  "<<sigma_col<<"  "
                           <<setw(15)<<"  "<<OPs_new_.value[OP_no]<<"   "<<abs(OPs_new_.value[OP_no])<<endl;
        }




        //Getting ready for next Temperature
        kick_while_cooling=0.0;
        Parameters_.Read_OPs=true;
        Parameters_.File_OPs_in = Parameters_.File_OPs_out + "_Temp"+string(temp_char) + ".txt";
        // Initialize();
    }


}

void Kspace_calculation_G2dLatticeNew::Diagonalize(char option){

    //extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
    //                       std::complex<double> *,int *, double *, int *);


    char jobz=option;
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=Ham_.n_row();
    int lda=Ham_.n_col();
    vector<complex<double>> work(3);
    vector<double> rwork(3*n -2);
    int info;
    int lwork= -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(),eigs_.end(),0);
    // query:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0].real()));
    work.resize(lwork);
    // real work:
    zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}


}




void Kspace_calculation_G2dLatticeNew::DiagonalizeGivenMat(char option, Matrix<complex<double>> & HamGiven_, vector<double> & EigsOut_){


    //extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
    //                       std::complex<double> *,int *, double *, int *);


    char jobz=option;
    char uplo='L'; //WHY ONLY 'L' WORKS?
    int n=HamGiven_.n_row();
    int lda=HamGiven_.n_col();
    vector<complex<double>> work(3);
    vector<double> rwork(3*n -2);
    int info;
    int lwork= -1;

    EigsOut_.resize(HamGiven_.n_row());
    fill(EigsOut_.begin(),EigsOut_.end(),0);
    // query:
    zheev_(&jobz,&uplo,&n,&(HamGiven_(0,0)),&lda,&(EigsOut_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0].real()));
    work.resize(lwork);
    // real work:
    zheev_(&jobz,&uplo,&n,&(HamGiven_(0,0)),&lda,&(EigsOut_[0]),&(work[0]),&lwork,&(rwork[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}



}

void Kspace_calculation_G2dLatticeNew::Update_OrderParameters_AndersonMixing(int iter){

    int S_ = NSites_in_MUC;
    bool with_SVD=true;
    int m_;
    int row_, col_;
    int OP_size;
    Mat_1_int NewInd_to_OldInd;
    NewInd_to_OldInd.clear();
    for(int ind=0;ind<OPs_.value.size();ind++){
        NewInd_to_OldInd.push_back(ind);
    }

    for(int ind=0;ind<OPs_.value.size();ind++){
        row_=OPs_.rows[ind];
        col_=OPs_.columns[ind];
        if(row_!=col_){
            NewInd_to_OldInd.push_back(ind);
        }
    }

    assert(NewInd_to_OldInd.size()==OPs_.value.size() + (OPs_.value.size() - 2*n_orbs_*n_atoms_*S_));
    OP_size=NewInd_to_OldInd.size();


    if(iter==0){
        //        cout<<"Anderson mixing for iter "<<iter<<endl;

        x_k_.clear();x_k_.resize(OP_size);
        for(int i=0;i<OP_size;i++){
            if(i<OPs_.value.size()){
                x_k_[i] = OPs_.value[NewInd_to_OldInd[i]].real();
            }
            else{
                x_k_[i] = OPs_.value[NewInd_to_OldInd[i]].imag();
            }
        }

        f_k_.clear();f_k_.resize(OP_size);
        for(int i=0;i<OP_size;i++){
            if(i<OPs_new_.value.size()){
                f_k_[i] = OPs_new_.value[NewInd_to_OldInd[i]].real() - OPs_.value[NewInd_to_OldInd[i]].real();
            }
            else{
                f_k_[i] = OPs_new_.value[NewInd_to_OldInd[i]].imag()- OPs_.value[NewInd_to_OldInd[i]].imag();
            }
        }
        assert(OPs_new_.value.size() == OPs_.value.size());

        //f_k = OPs_new_.value;
        //x_k = OPs_.value;

        for(int ind=0;ind<OPs_.value.size();ind++){
            OPs_.value[ind] = (1-Parameters_.alpha_OP)*OPs_.value[ind]
                    + Parameters_.alpha_OP*OPs_new_.value[ind];
        }


        x_km1_=x_k_;
        X_mat.resize(0,0);

        f_km1_=f_k_;
        F_mat.resize(0,0);
        //f_km1=

    }
    else{
        //      cout<<"Anderson mixing for iter "<<iter<<endl;
        x_k_.clear();x_k_.resize(OP_size);
        for(int i=0;i<OP_size;i++){
            if(i<OPs_.value.size()){
                x_k_[i] = OPs_.value[NewInd_to_OldInd[i]].real();
            }
            else{
                x_k_[i] = OPs_.value[NewInd_to_OldInd[i]].imag();
            }
        }

        f_k_.clear();f_k_.resize(OP_size);
        for(int i=0;i<OP_size;i++){
            if(i<OPs_new_.value.size()){
                f_k_[i] = OPs_new_.value[NewInd_to_OldInd[i]].real() - OPs_.value[NewInd_to_OldInd[i]].real();
            }
            else{
                f_k_[i] = OPs_new_.value[NewInd_to_OldInd[i]].imag() - OPs_.value[NewInd_to_OldInd[i]].imag();
            }
        }

        Del_x_km1.clear();Del_x_km1.resize(OP_size);
        Del_f_km1.clear();Del_f_km1.resize(OP_size);
        for(int i=0;i<OP_size;i++){
            Del_x_km1[i] = x_k_[i] - x_km1_[i];
            Del_f_km1[i] = f_k_[i] - f_km1_[i];
        }


        m_=min(Parameters_.AM_m,iter);
        gamma_k_.clear();
        gamma_k_.resize(m_);

        //updating X_mat_k
        Matrix <double> Xmat_temp;
        Xmat_temp.resize(X_mat.n_row(),X_mat.n_col());
        Xmat_temp=X_mat;
        X_mat.resize(OP_size,m_);

        if(iter<=Parameters_.AM_m){
            for(col_=0;col_<Xmat_temp.n_col();col_++){
                for(row_=0;row_<Xmat_temp.n_row();row_++){
                    X_mat(row_,col_) = Xmat_temp(row_,col_);
                }
            }

            for(col_=m_-1;col_<m_;col_++){
                for(row_=0;row_<OP_size;row_++){
                    X_mat(row_,col_) = Del_x_km1[row_];
                }
            }
        }
        else{
            for(col_=1;col_<Xmat_temp.n_col();col_++){
                for(row_=0;row_<Xmat_temp.n_row();row_++){
                    X_mat(row_,col_-1) = Xmat_temp(row_,col_);
                }
            }
            for(row_=0;row_<OP_size;row_++){
                X_mat(row_,m_-1) = Del_x_km1[row_];
            }
        }



        //updating F_mat_k
        Matrix <double> Fmat_temp;
        Fmat_temp.resize(F_mat.n_row(),F_mat.n_col());
        Fmat_temp=F_mat;
        F_mat.resize(OP_size,m_);

        if(iter<=Parameters_.AM_m){
            for(col_=0;col_<Fmat_temp.n_col();col_++){
                for(row_=0;row_<Fmat_temp.n_row();row_++){
                    F_mat(row_,col_) = Fmat_temp(row_,col_);
                }
            }

            for(col_=m_-1;col_<m_;col_++){
                for(row_=0;row_<OP_size;row_++){
                    F_mat(row_,col_) = Del_f_km1[row_];
                }
            }
        }
        else{
            for(col_=1;col_<Fmat_temp.n_col();col_++){
                for(row_=0;row_<Fmat_temp.n_row();row_++){
                    F_mat(row_,col_-1) = Fmat_temp(row_,col_);
                }
            }
            for(row_=0;row_<OP_size;row_++){
                F_mat(row_,m_-1) = Del_f_km1[row_];
            }
        }

        //cout<<"here 1"<<endl;

        //Update gamma_k using Total least sqaure minimaztion (using SVD of F_mat)
        if(with_SVD==false){
            for(int i=0;i<m_;i++){
                gamma_k_[i] = 1.0/(1.0*m_);
            }
        }
        else{
            int r_;
            r_=min(OP_size,m_);
            Matrix<double> A_;  //nxm; n=OP_size
            Matrix<double> VT_; //mxm
            Matrix<double> U_;  //nxn
            vector<double> Sigma_; //atmost non-zero min(n,m) values
            A_.resize(F_mat.n_row(), F_mat.n_col());
            A_=F_mat;

            Perform_SVD(A_,VT_,U_,Sigma_);

            Matrix<double> UT_f;
            Matrix<double> Sinv_UT_f;

            UT_f.resize(r_,1);
            for(int i=0;i<r_;i++){
                for(int j=0;j<OP_size;j++){
                    UT_f(i,0) += U_(j,i)*f_k_[j];
                }
            }

            Sinv_UT_f.resize(r_,1);//S-inv in Pseudoinverse of Sigma_
            for(int i=0;i<r_;i++){
                if(abs(Sigma_[i])>=0.001){
                    Sinv_UT_f(i,0) = (1.0/Sigma_[i])*UT_f(i,0);
                }
                else{
                    Sinv_UT_f(i,0)=0.0;
                }
            }

            double sum_gamma=0.0;
            for(int i=0;i<m_;i++){
                gamma_k_[i]=0.0;
                for(int j=0;j<r_;j++){
                    gamma_k_[i] += VT_(j,i)*Sinv_UT_f(j,0);
                }
                sum_gamma += abs(gamma_k_[i]);

            }

            if(sum_gamma>1){
                for(int i=0;i<m_;i++){
                    gamma_k_[i] = gamma_k_[i]*(1.0/sum_gamma);
                }
            }


            //            A_.clear();
            //            VT_.clear();
            //            U_.clear();
            //            Sigma_.clear();
            //            UT_f.clear();
            //            Sinv_UT_f.clear();

        }


        //cout<<"here 3"<<endl;


        //Mat_1_doub Temp_F_gamma_k, Temp_X_gamma_k;
        xbar_k_.clear();fbar_k_.clear();
        xbar_k_.resize(OP_size);
        fbar_k_.resize(OP_size);
        double temp_f, temp_x;
        for(int i=0;i<OP_size;i++){
            temp_f=0.0;
            temp_x=0.0;
            for(int j=0;j<m_;j++){
                temp_f +=F_mat(i,j)*gamma_k_[j];
                temp_x +=X_mat(i,j)*gamma_k_[j];
            }
            xbar_k_[i] = x_k_[i] - 1.0*temp_x;
            fbar_k_[i] = f_k_[i] - 1.0*temp_f;
        }


        x_kp1_.clear();
        x_kp1_.resize(OP_size);
        for(int i=0;i<OP_size;i++){
            x_kp1_[i] = (1.0 - 0.0*Parameters_.alpha_OP)*xbar_k_[i]  +
                    Parameters_.alpha_OP*fbar_k_[i];
        }

        for(int i=0;i<OP_size;i++){
            if(i<OPs_.value.size()){
                OPs_.value[NewInd_to_OldInd[i]].real( x_kp1_[i]);
            }
            else{
                OPs_.value[NewInd_to_OldInd[i]].imag( x_kp1_[i]);
            }
        }


        //---saving arrays for next iteration-----
        x_km1_=x_k_;
        f_km1_=f_k_;

    }





    //    double OPs_total_den_=0.0;
    //    int index_OP;
    //    int row_OP, col_OP;
    //    for(int k1=0;k1<lx_cells;k1++){
    //        for(int k2=0;k2<ly_cells;k2++){
    //            for(int alpha=0;alpha<S_;alpha++){
    //                for(int gamma=0;gamma<n_orbs_;gamma++){
    //                    for(int spin=0;spin<2;spin++){
    //                        row_OP = alpha + gamma*(S_)
    //                                + spin*(n_orbs_*S_);
    //                        col_OP = alpha + gamma*(S_)
    //                                + spin*(n_orbs_*S_);
    //                        index_OP = SI_to_ind[col_OP + row_OP*(2*n_orbs_*ncells_*S_)];
    //                        OPs_total_den_ += OPs_.value[index_OP].real();
    //                    }
    //                }
    //            }
    //        }
    //    }


    //    double ratio = Parameters_.Total_Particles/(OPs_total_den_+0.0001);

    //    for(int alpha=0;alpha<S_;alpha++){
    //        for(int gamma=0;gamma<n_orbs_;gamma++){
    //            for(int spin=0;spin<2;spin++){
    //                row_OP = alpha + gamma*(S_)
    //                        + spin*(n_orbs_*S_);
    //                col_OP = alpha + gamma*(S_)
    //                        + spin*(n_orbs_*S_);
    //                index_OP = SI_to_ind[col_OP + row_OP*(2*n_orbs_*ncells_*S_)];
    //                OPs_.value[index_OP] += ratio/(1.0+ratio);
    //            }
    //        }
    //    }






}

void Kspace_calculation_G2dLatticeNew::Perform_SVD(Matrix<double> & A_, Matrix<double> & VT_, Matrix<double> & U_, vector<double> & Sigma_){


    char jobz='S'; //A,S,O,N

    // m>n willb ethe case in general
    int m=A_.n_row();
    int n=A_.n_col();
    int lda=m;
    int ldu=m;
    int ldvt=n;
    int min_mn=min(m,n);

    Sigma_.clear();
    Sigma_.resize(min(m,n));

    U_.resize(ldu,min_mn);

    VT_.resize(ldvt,n);


    //cout<<"here 1"<<endl;

    vector<double> work(3);
    int info;
    int lwork= -1;
    vector<int> iwork(8*min(m,n));

    // query:
    dgesdd_(&jobz, &m, &n, &(A_(0,0)),&lda, &(Sigma_[0]),&(U_(0,0)), &ldu, &(VT_(0,0)), &ldvt,
            &(work[0]), &lwork, &(iwork[0]), &info);
    //cout<<"here 1.1"<<endl;
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0]));
    work.resize(lwork);

    //cout<<"here 1.2"<<endl;

    // real work:
    dgesdd_(&jobz, &m, &n, &(A_(0,0)),&lda, &(Sigma_[0]),&(U_(0,0)), &ldu, &(VT_(0,0)), &ldvt,
            &(work[0]), &lwork, &(iwork[0]), &info);
    if (info!=0) {
        if(info>0){
            std::cerr<<"info="<<info<<"\n";
            perror("diag: zheev: failed with info>0.\n");}
        if(info<0){
            std::cerr<<"info="<<info<<"\n";
            perror("diag: zheev: failed with info<0.\n");
        }
    }

    //cout<<"here 2"<<endl;

}


void Kspace_calculation_G2dLatticeNew::Perform_SVD_complex(Matrix<complex<double>> & A_, Matrix<complex<double>> & VT_, Matrix<complex<double>> & U_, vector<double> & Sigma_){


    char jobz='A'; //A,S,O,N

    int m=A_.n_row();
    int n=A_.n_col();
    int lda=A_.n_row();
    int ldu=A_.n_row();
    int ldvt=n;

    Sigma_.clear();
    Sigma_.resize(min(m,n));

    U_.resize(ldu,m);

    VT_.resize(ldvt,n);


    vector<complex<double>> work(3);
    int info;
    int lwork= -1;
    vector<int> iwork(8*min(m,n));
    int lrwork = max( (5*min(m,n)*min(m,n)) + 5*min(m,n), (2*max(m,n)*min(m,n)) + (2*min(m,n)*min(m,n)) + min(m,n) );
    vector<double> rwork(lrwork);

    // query:
    zgesdd_(&jobz, &m, &n, &(A_(0,0)),&lda, &(Sigma_[0]),&(U_(0,0)), &ldu, &(VT_(0,0)), &ldvt,
            &(work[0]), &lwork, &(rwork[0]), &(iwork[0]), &info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0]).real());
    work.resize(lwork);
    // real work:
    zgesdd_(&jobz, &m, &n, &(A_(0,0)),&lda, &(Sigma_[0]),&(U_(0,0)), &ldu, &(VT_(0,0)), &ldvt,
            &(work[0]), &lwork, &(rwork[0]), &(iwork[0]), &info);
    if (info!=0) {
        if(info>0){
            std::cerr<<"info="<<info<<"\n";
            perror("diag: zheev: failed with info>0.\n");}
        if(info<0){
            std::cerr<<"info="<<info<<"\n";
            perror("diag: zheev: failed with info<0.\n");
        }
    }

    // Ham_.print();



}




#endif
