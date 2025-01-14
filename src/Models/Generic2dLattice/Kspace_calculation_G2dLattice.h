#include <algorithm>
#include <functional>
#include <math.h>
#include "tensor_type.h"
#include "Parameters_G2dLattice.h"
#include "Coordinates_G2dLattice.h"
#include "Connections_G2dLattice.h"
#include "random"
#include "../../Matrix.h"
#define PI acos(-1.0)

#ifndef Kspace_calculation_G2dLattice_class
#define Kspace_calculation_G2dLattice_class

extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
                       std::complex<double> *,int *, double *, int *);

extern "C" void dgesdd_ (char *, int *, int *, double *, int *, double *, double *, int *, double *, int*,
                        double *, int *, int *, int *);


extern "C" void zgesdd_ (char *, int *, int *, std::complex<double> *, int *, double *, std::complex<double> *, int *, std::complex<double> *, int*,
                        std::complex<double> *, int *, double * , int *, int *);

class Kspace_calculation_G2dLattice
{
public:
    Kspace_calculation_G2dLattice(Parameters_G2dLattice &Parameters__, Coordinates_G2dLattice &Coordinates__, Connections_G2dLattice &Connections__, mt19937_64& Generator1__ )
        : Parameters_(Parameters__), Coordinates_(Coordinates__), Connections_(Connections__), Generator1_(Generator1__)

    {
        Initialize();
    }

    void Initialize();                                     //::DONE
    void Diagonalize(char option);
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
    void Calculate_Akw(double omega_wrt_mu);
    //::DONE



    mt19937_64 &Generator1_;
    uniform_real_distribution<double> dis1_;
    Parameters_G2dLattice &Parameters_;
    Coordinates_G2dLattice &Coordinates_; //this in cell wise representation, n_orbs=no. of atoms in unitcell
    Connections_G2dLattice &Connections_;
    int lx_, ly_, ncells_, n_orbs_, n_atoms_, UnitCellSize_x, UnitCellSize_y, lx_cells, ly_cells;
    Matrix<complex<double>> Ham_;
    vector<double> eigs_;
    Mat_2_Complex_doub Eigvectors_;
    Mat_1_doub Eigenvalues_;
    Mat_2_Complex_doub Eigvectors_saved;
    Mat_1_doub Eigenvalues_saved;
    Mat_1_doub Kx_values, Ky_values;
    //Mat_1_Complex_doub V_mat;

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


};

void Kspace_calculation_G2dLattice::Hall_conductance(){

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
    FACTOR_ = 1.0;//2.0/(sqrt(UnitCellSize_x*UnitCellSize_y*1.0));

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


void Kspace_calculation_G2dLattice::Create_Current_Oprs(){
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



    int S_ = UnitCellSize_x*UnitCellSize_y;
    int k_index;
    int state_p,state_n;
    int sigma, sigma_p;
    sigma=0;sigma_p=0;
    int d1_net, d2_net;
    double fac1, fac2, facx, facy, dx_, dy_;

    //2*n_orbs_*UnitCellSize_x*UnitCellSize_y*k_index + row
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

                                        //alpha + gamma*(UnitCellSize_x*UnitCellSize_y)
                                        // + spin1*(n_orbs_*UnitCellSize_x*UnitCellSize_y);
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


void Kspace_calculation_G2dLattice::Create_Current_Oprs_Faster(){


    J_KE_e1.resize(2);J_KE_e2.resize(2);
    J_KE_X.resize(2); J_KE_Y.resize(2);


    for(int sigma=0;sigma<2;sigma++){
        J_KE_e1[sigma].resize(lx_*ly_*2*n_atoms_*n_orbs_, lx_*ly_*2*n_atoms_*n_orbs_);
        J_KE_e2[sigma].resize(lx_*ly_*2*n_atoms_*n_orbs_, lx_*ly_*2*n_atoms_*n_orbs_);

        J_KE_X[sigma].resize(lx_*ly_*2*n_atoms_*n_orbs_, lx_*ly_*2*n_atoms_*n_orbs_);
        J_KE_Y[sigma].resize(lx_*ly_*2*n_atoms_*n_orbs_, lx_*ly_*2*n_atoms_*n_orbs_);
    }

    //    c1 = alpha + UP_*(S_);
    //    c2 = alpha + DOWN_*(S_);;
    //    for(int n=0;n<2*S_;n++){ //band_index
    //        for(int k1=0;k1<lx_cells;k1++){
    //            for(int k2=0;k2<ly_cells;k2++){
    //                k_index = Coordinates_.Ncell(k1,k2);
    //                state = 2*S_*k_index + n;
    //                val += (1.0/ncells_)*(
    //                            (conj(Eigvectors_[state][c1])*Eigvectors_[state][c2])
    //                            *(1.0/( exp((Eigenvalues_saved[state]-mu_)*Parameters_.beta ) + 1.0))
    //                            );
    //            }
    //        }
    //    }



    int row_, col_;
    int d1_, d2_;
    int d1_org, d2_org;
    int alpha_1, alpha_2, alpha_p_1, alpha_p_2;
    int c1, c2;


    double spin_factor;
    //    temp_val=0.0;
    //    for(int cell_no=0;cell_no<ncells_;cell_no++){
    //        d1_org = Coordinates_.indx_cellwise(cell_no);
    //        d2_org = Coordinates_.indy_cellwise(cell_no);

    //        row_ = ( (d1_org*UnitCellSize_x) + alpha_1) + ((d2_org*UnitCellSize_y) + alpha_2)*lx_ + sigma*(lx_*ly_);
    //        col_ = (0 + alpha_p_1) + (0 + alpha_p_2)*lx_ + sigma_p*(lx_*ly_);

    //        Get_minimum_distance_direction(0, cell_no, d1_, d2_);

    //        temp_val += Connections_.HTB_(row_,col_)*exp(iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_   ));

    //    }


    int S_ = UnitCellSize_x*UnitCellSize_y;
    int k_index;
    int state_p,state_n;
    int sigma, sigma_p;
    sigma=0;sigma_p=0;
    int d1_net, d2_net;
    double fac1, fac2, facx, facy, dx_, dy_;

    //2*n_orbs_*UnitCellSize_x*UnitCellSize_y*k_index + row


    for(int sigma=0;sigma<2;sigma++){
        //spin_factor=((2.0*sigma) - 1.0)*1.0;
        //spin_factor=(1-sigma);
        spin_factor=1.0;
        //cout<<"spin_factor for spin "<<  sigma <<   " = "<<spin_factor<<endl;
        //     spin_factor=1.0;
        //   spin_factor=abs(1-sigma);

        for(int alpha=0;alpha<S_;alpha++){
            for(int gamma=0;gamma<n_atoms_*2;gamma++){

                for(int alpha_p=0;alpha_p<S_;alpha_p++){
                    for(int gamma_p=0;gamma_p<n_atoms_*2;gamma_p++){

                        //alpha + gamma*(UnitCellSize_x*UnitCellSize_y)
                        // + spin1*(n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        sigma_p=sigma;
                        c1 = alpha + gamma*(S_) + sigma*(2*n_atoms_*S_);
                        c2 = alpha_p + gamma_p*(S_) + sigma*(2*n_atoms_*S_);

                        alpha_1 = alpha % UnitCellSize_x; //a = a1 +a2*lx
                        alpha_p_1 = alpha_p % UnitCellSize_x;
                        alpha_2 = (alpha - alpha_1)/UnitCellSize_x;
                        alpha_p_2 = (alpha_p - alpha_p_1)/UnitCellSize_x;

                        for(int cell_no=0;cell_no<ncells_;cell_no++){

                            //------------------------


                            d1_org = Coordinates_.indx_cellwise(cell_no);
                            d2_org = Coordinates_.indy_cellwise(cell_no);

                            row_ = gamma  + ( (d1_org*UnitCellSize_x) + alpha_1)*n_atoms_*n_orbs_
                                   + (( (d2_org*UnitCellSize_y) + alpha_2)*lx_*n_atoms_*n_orbs_)
                                   + sigma*(lx_*ly_*n_atoms_*n_orbs_);
                            col_ = ( gamma_p + (0 + alpha_p_1)*n_atoms_*n_orbs_ + ((0 + alpha_p_2)*lx_*n_atoms_*n_orbs_)) + sigma_p*(lx_*ly_*n_atoms_*n_orbs_);


                            if(abs(Connections_.HTB_(row_,col_)) > 0.0000001){


                                Get_minimum_distance_direction(0, cell_no, d1_, d2_);

                                d1_net = ( (d1_*UnitCellSize_x) + alpha_1) - (0 + alpha_p_1);
                                d2_net = ((d2_*UnitCellSize_y) + alpha_2) - (0 + alpha_p_2);


                                dx_ = d1_net + 0.5*d2_net + ((1.0/sqrt(3.0))*(gamma-gamma_p)*(sqrt(3.0)/2.0));
                                dy_ = (sqrt(3.0)/2.0)*d2_net   + ((1.0/sqrt(3.0))*(gamma-gamma_p)*(1.0/2.0));

                                if( ((dx_*dx_) + (dy_*dy_))>0.000000001){
                                    facx = (dx_/sqrt((dx_*dx_) + (dy_*dy_)));
                                    facy = (dy_/sqrt((dx_*dx_) + (dy_*dy_)));
                                }
                                else{
                                    facx=0;
                                    facy=0;
                                }
                                //                                if(dx_==0){
                                //                                    facx=0;
                                //                                }
                                //                                else{
                                //                                    facx=dx_/abs(dx_);
                                //                                }
                                //                                if(dy_==0){
                                //                                    facy=0;
                                //                                }
                                //                                else{
                                //                                    facy=dy_/abs(dy_);
                                //                                }



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


                                for(int l=0;l<2*n_orbs_*n_atoms_*S_;l++){
                                    for(int m=0;m<2*n_orbs_*n_atoms_*S_;m++){

                                        for(int k1=0;k1<lx_cells;k1++){
                                            for(int k2=0;k2<ly_cells;k2++){
                                                k_index = Coordinates_.Ncell(k1,k2);
                                                state_p = 2*n_orbs_*S_*k_index + l;
                                                state_n = 2*n_orbs_*S_*k_index + m;


                                                J_KE_X[sigma](state_p,state_n) += spin_factor*(0.5)*facx*iota_complex*(
                                                                                       (Connections_.HTB_(row_,col_)*conj(Eigvectors_[state_p][c2])*Eigvectors_[state_n][c1]*exp(iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_ )))
                                                                                       -(Connections_.HTB_(col_,row_)*conj(Eigvectors_[state_p][c1])*Eigvectors_[state_n][c2]*exp(-1.0*iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_ )))
                                                                                       );



                                                J_KE_Y[sigma](state_p,state_n) += (0.5)*facy*iota_complex*(
                                                                                       (Connections_.HTB_(row_,col_)*conj(Eigvectors_[state_p][c2])*Eigvectors_[state_n][c1]*exp(iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_ )))
                                                                                       -(Connections_.HTB_(col_,row_)*conj(Eigvectors_[state_p][c1])*Eigvectors_[state_n][c2]*exp(-1.0*iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_ )))
                                                                                       );



                                                J_KE_e1[sigma](state_p,state_n) += spin_factor*(0.5)*fac1*iota_complex*(
                                                                                        (Connections_.HTB_(row_,col_)*conj(Eigvectors_[state_p][c2])*Eigvectors_[state_n][c1]*exp(iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_ )))
                                                                                        -(Connections_.HTB_(col_,row_)*conj(Eigvectors_[state_p][c1])*Eigvectors_[state_n][c2]*exp(-1.0*iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_ )))
                                                                                        );



                                                J_KE_e2[sigma](state_p,state_n) += (0.5)*fac2*iota_complex*(
                                                                                        (Connections_.HTB_(row_,col_)*conj(Eigvectors_[state_p][c2])*Eigvectors_[state_n][c1]*exp(iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_ )))
                                                                                        -(Connections_.HTB_(col_,row_)*conj(Eigvectors_[state_p][c1])*Eigvectors_[state_n][c2]*exp(-1.0*iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_ )))
                                                                                        );


                                                //                                                J_KE_e1(state_p,state_n) += 0.0;



                                                //                                                J_KE_e2(state_p,state_n) += 0.0;







                                            }




                                        }

                                    }
                                }
                            }
                        }
                    }

                }
                cout<<"sigma(2), alpha(" << S_ << "), gamma(4) = "<<sigma<<" "<<alpha<<" "<<gamma<<endl;
            }
        }
    }




}



void Kspace_calculation_G2dLattice::Calculate_ChernNumbers(){



    Matrix<complex<double>> F_mat; //F1, F2, F3, F4, F5;
    F_mat.resize(n_atoms_*4*UnitCellSize_x*UnitCellSize_y, lx_cells*ly_cells);

    complex<double> Ux_k, Uy_k, Ux_kpy, Uy_kpx;
    vector<complex<double>> F_bands;
    F_bands.resize(n_atoms_*4*UnitCellSize_x*UnitCellSize_y);
    vector<complex<double>> Chern_num;
    Chern_num.resize(n_atoms_*4*UnitCellSize_x*UnitCellSize_y);

    vector<complex<double>> F_bands_orgnl;
    F_bands_orgnl.resize(n_atoms_*4*UnitCellSize_x*UnitCellSize_y);
    vector<complex<double>> Chern_num_orgnl;
    Chern_num_orgnl.resize(n_atoms_*4*UnitCellSize_x*UnitCellSize_y);
    for (int band = 0; band < n_atoms_*4*UnitCellSize_x*UnitCellSize_y; band++)
    {
        string file_Fk="Fk_band"+to_string(band)+".txt";
        ofstream fl_Fk_out(file_Fk.c_str());
        fl_Fk_out<<"#nx  ny  tilde_F(nx,ny).real()  tilde_F(nx,ny).imag()  ArgofLog.real()  ArgofLog.imag()"<<endl;
        fl_Fk_out<<"#Extra momentum point for pm3d corners2color c1"<<endl;

        //        string file_Fk_orgnl="Fk_original_band"+to_string(band)+".txt";
        //        ofstream fl_Fk_orgnl_out(file_Fk_orgnl.c_str());
        //        fl_Fk_orgnl_out<<"#nx  ny  F(nx,ny).real()*(2pi/Lx)*(2pi/Ly)  F(nx,ny).imag()*(2pi/Lx)*(2pi/Ly)"<<endl;

        F_bands[band] = 0.0;
        F_bands_orgnl[band] = 0.0;
        for (int nx = 0; nx < lx_cells; nx++)
        {
            for (int ny = 0; ny < ly_cells; ny++)
            {
                int n = Coordinates_.Ncell(nx,ny);
                int n_left, n_right, nx_left, ny_left, nx_right, ny_right;

                //U1_k
                Ux_k = 0;
                n_left = n;
                nx_right = (nx + 1) % lx_cells;
                ny_right = ny;
                n_right = Coordinates_.Ncell(nx_right, ny_right);
                for (int comp = 0; comp < 4*n_atoms_*UnitCellSize_x*UnitCellSize_y; comp++)
                {
                    Ux_k +=
                        conj(Eigvectors_[4*n_atoms_*UnitCellSize_x*UnitCellSize_y*n_left +  band][comp]) *
                        Eigvectors_[4*n_atoms_*UnitCellSize_x*UnitCellSize_y*n_right + band][comp];
                }
                Ux_k = Ux_k * (1.0 / abs(Ux_k));

                //U2_kpx
                Uy_kpx = 0;
                nx_left = (nx + 1) % lx_cells;
                ny_left = ny;
                n_left = Coordinates_.Ncell(nx_left, ny_left);
                nx_right = nx_left;
                ny_right = (ny_left + 1) % ly_cells;
                n_right = Coordinates_.Ncell(nx_right, ny_right);
                for (int comp = 0; comp < 4*n_atoms_*UnitCellSize_x*UnitCellSize_y; comp++)
                {
                    Uy_kpx +=
                        conj(Eigvectors_[4*n_atoms_*UnitCellSize_x*UnitCellSize_y*n_left +  band][comp]) *
                        Eigvectors_[4*n_atoms_*UnitCellSize_x*UnitCellSize_y*n_right + band][comp];
                }
                Uy_kpx = Uy_kpx * (1.0 / abs(Uy_kpx));

                //U1_kpy
                Ux_kpy = 0;
                nx_left = nx;
                ny_left = (ny + 1) % ly_cells;
                n_left = Coordinates_.Ncell(nx_left, ny_left);
                nx_right = (nx_left + 1) % lx_cells;
                ny_right = ny_left;
                n_right = Coordinates_.Ncell(nx_right, ny_right);
                for (int comp = 0; comp < 4*n_atoms_*UnitCellSize_x*UnitCellSize_y; comp++)
                {
                    Ux_kpy +=
                        conj(Eigvectors_[4*n_atoms_*UnitCellSize_x*UnitCellSize_y*n_left +  band][comp]) *
                        Eigvectors_[4*n_atoms_*UnitCellSize_x*UnitCellSize_y*n_right + band][comp];
                }
                Ux_kpy = Ux_kpy * (1.0 / abs(Ux_kpy));

                //U2_k
                Uy_k = 0;
                nx_left = nx;
                ny_left = ny;
                n_left = Coordinates_.Ncell(nx_left, ny_left);
                nx_right = nx_left;
                ny_right = (ny_left + 1) % ly_cells;
                n_right = Coordinates_.Ncell(nx_right, ny_right);
                for (int comp = 0; comp < 4*n_atoms_*UnitCellSize_x*UnitCellSize_y; comp++)
                {
                    Uy_k +=
                        conj(Eigvectors_[4*n_atoms_*UnitCellSize_x*UnitCellSize_y*n_left +  band][comp]) *
                        Eigvectors_[4*n_atoms_*UnitCellSize_x*UnitCellSize_y*n_right + band][comp];
                }
                Uy_k = Uy_k * (1.0 / abs(Uy_k));

                // Calculating tilde F12
                F_mat(band, n) = log(Ux_k *
                                     Uy_kpx *
                                     conj(Ux_kpy) * conj(Uy_k));

                F_bands[band] += F_mat(band, n);


                fl_Fk_out.precision(10);

                fl_Fk_out<<nx<<"  "<<ny<<"  "<<F_mat(band, n).real()<<"  "<<F_mat(band, n).imag()<<
                    "  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).real()<<
                    "  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).imag()<<endl;




                if(abs((abs(F_mat(band, n).imag()) - PI))<0.0000001){
                    cout<<ny<<"  "<<nx<<"  gives Pi for band"<< band <<endl;
                    // assert (abs((abs(F_mat(band, n).imag()) - M_PI))>0.0000001);
                }


                if(ny==ly_cells-1){//For pm3d corners2color c1
                    fl_Fk_out<<nx<<"  "<<ny<<"  "<<F_mat(band, n).real()<<"  "<<F_mat(band, n).imag()<<
                        "  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).real()<<
                        "  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).imag()<<endl;
                }
            }
            if(nx==lx_cells-1){//For pm3d corners2color c1
                fl_Fk_out<<endl;
                for(int ny_=0;ny_<ly_cells;ny_++){
                    int n_ = nx + lx_cells*ny_;
                    fl_Fk_out<<nx<<"  "<<ny_<<"  "<<F_mat(band, n_).real()<<"  "<<F_mat(band, n_).imag()<<
                        "  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).real()<<
                        "  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).imag()<<endl;
                }
            }
            fl_Fk_out<<endl;
        }


        Chern_num[band] = (-1.0 * iota_complex / (2 * PI)) * F_bands[band];
        // Chern_num_orgnl[band] = (-1.0 * iota_complex / (2 * PI)) * F_bands_orgnl[band];
        fl_Fk_out<<"#Chern no*2pi*Iota= "<<F_bands[band].real()<<"  "<<F_bands[band].imag()<<endl;
        cout << "tilde Chern number [" << band << "] = " << Chern_num[band].real() << "        " << Chern_num[band].imag() << endl;
        //  cout << "Chern number [" << band << "] = " << Chern_num_orgnl[band].real() << " " << Chern_num_orgnl[band].imag() << endl;

    }

}

double Kspace_calculation_G2dLattice::Lorentzian(double x, double brd)
{
    double temp;

    temp = (1.0 / PI) * ((brd / 2.0) / ((x * x) + ((brd * brd) / 4.0)));

    return temp;
}


void Kspace_calculation_G2dLattice::Calculate_Akw(double omega_wrt_mu)
{

    assert(Parameters_.UnitCellSize_x==1);
    assert(Parameters_.UnitCellSize_y==1);

    //---------Read from input file-----------------------//
    string fileout = "AkwCut.txt";
    double omega_min, omega_max, d_omega;

    double eta = 0.01;
    double omega_target = omega_wrt_mu + mu_;

    //---------------------------------------------------//

    ofstream file_Akw_out(fileout.c_str());
    file_Akw_out<<"#k1  k2  Akw_up_orb0(w_target)  Akw_up_orb1(w_target) Akw_up_orb2(w_target)  Akw_dn_orbj(w_target) ...."<<endl;
    double Akw_spin=0.0;
    int S_ = UnitCellSize_x*UnitCellSize_y;
    assert(S_==1);
    int c1;

    int k_index, state;
    for(int k1=0;k1<(lx_cells+1);k1++){
        int k1_temp=k1%lx_cells;

        for(int k2=0;k2<(ly_cells+1);k2++){
            int k2_temp=k2%ly_cells;
            k_index = Coordinates_.Ncell(k1_temp,k2_temp);
            file_Akw_out<<k1<<"  "<<k2;

            for(int spin=0;spin<2;spin++){
                for(int gamma_=0;gamma_<n_atoms_*n_orbs_;gamma_++){
                    Akw_spin=0.0;

                    c1 = 0 + gamma_*(S_) +  spin*(n_atoms_*n_orbs_*S_);

                    for(int n=0;n<2*n_atoms_*n_orbs_*S_;n++){ //band_index

                        state = 2*n_atoms_*n_orbs_*S_*k_index + n;
                        Akw_spin += abs(Eigvectors_[state][c1])*abs(Eigvectors_[state][c1])*
                                    Lorentzian( omega_target-Eigenvalues_saved[state], eta);
                    }

                    file_Akw_out<<"  "<<Akw_spin;
                }
            }
            file_Akw_out<<endl;
        }
        file_Akw_out<<endl;
    }

}



void Kspace_calculation_G2dLattice::Calculate_Nw()
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


void Kspace_calculation_G2dLattice::Get_minimum_distance_direction(int l,int m,int &r1_, int &r2_){

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


void Kspace_calculation_G2dLattice::Get_Bands(){

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




    //    cout<<"PRINTING PATH"<<endl;
    //    for (int k_point = 0; k_point < k_path.size(); k_point++)
    //    {
    //        cout<<k_path[k_point].first<< "   "<<k_path[k_point].second<<endl;
    //    }
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
        for(int band=0;band<2*n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y;band++){

            int state = 2*n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y*k_index + band;
            file_out_bands<<Eigenvalues_[state];
            double up_contr =0;
            double dn_contr =0;
            int c_up, c_dn;
            for(int gamma_temp=0;gamma_temp<n_atoms_*n_orbs_;gamma_temp++){
                up_contr=0.0;
                dn_contr=0.0;
                for(int alpha_temp=0;alpha_temp<UnitCellSize_x*UnitCellSize_y;alpha_temp++){
                    c_up = alpha_temp + gamma_temp*(UnitCellSize_x*UnitCellSize_y) +  0*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                    c_dn = alpha_temp + gamma_temp*(UnitCellSize_x*UnitCellSize_y) +  1*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                    up_contr += abs(Eigvectors_[state][c_up])*abs(Eigvectors_[state][c_up]);
                    dn_contr += abs(Eigvectors_[state][c_dn])*abs(Eigvectors_[state][c_dn]);

                }
                file_out_bands<<"   "<<up_contr<<"   "<<dn_contr;
            }

            file_out_bands<<"   ";
        }
        file_out_bands<<endl;
    }


    for (int k_point = 0; k_point < k_path2.size(); k_point++)
    {
        k_index=Coordinates_.Ncell(k_path2[k_point].first, k_path2[k_point].second);

        file_out_bands2<<k_point<<"   ";
        for(int band=0;band<2*n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y;band++){
            file_out_bands2<<Eigenvalues_[2*n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y*k_index + band]<<"   ";
        }
        file_out_bands2<<endl;
    }



    for (int k_point = 0; k_point < k_path4.size(); k_point++)
    {
        k_index=Coordinates_.Ncell(k_path4[k_point].first, k_path4[k_point].second);

        file_out_bands4<<k_point<<"   ";
        for(int band=0;band<2*n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y;band++){
            file_out_bands4<<Eigenvalues_[2*n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y*k_index + band]<<"   ";
        }
        file_out_bands4<<endl;
    }

    file_out_bands4<<"# Gamma(" <<counter_Gamma1<<")"<<"-->"<<"K(" <<counter_K<<")"<<"-->"<<"K_prime(" <<counter_Kprime<<")"<<"-->"<<"Gamma(" <<counter_Gamma2<<")"<<"-->"<<"M(" <<counter_M<<")"<<endl;




    for (int k_point = 0; k_point < k_path5.size(); k_point++)
    {
        k_index=Coordinates_.Ncell(k_path5[k_point].first, k_path5[k_point].second);

        file_out_bands5<<k_point<<"   "<<k_path5[k_point].first<<"   "<<k_path5[k_point].second<<"   ";
        for(int band=0;band<2*n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y;band++){
            file_out_bands5<<Eigenvalues_[2*n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y*k_index + band]<<"   ";
        }
        file_out_bands5<<endl;
    }

    file_out_bands5<<"# Gamma"<<"-->"<<"K-->"<<"M-->"<<"Gamma"<<endl;



    for (int ky_ind_=0;ky_ind_<ly_cells;ky_ind_++)
    {
        for (int kx_ind_=0;kx_ind_<lx_cells;kx_ind_++)
        {
            k_index=Coordinates_.Ncell(kx_ind_, ky_ind_);

            file_out_bands3<<k_index<<"   "<<kx_ind_<<"   "<<ky_ind_<<"   ";
            for(int band=0;band<2*n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y;band++){
                file_out_bands3<<Eigenvalues_[2*n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y*k_index + band]<<"   ";
            }
            file_out_bands3<<endl;
        }
        file_out_bands3<<endl;
    }



}


void Kspace_calculation_G2dLattice::Get_Energies_new(){



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
    for(int j=0;j<UnitCellSize_x*UnitCellSize_y;j++){
        for(int b_beta=0;b_beta<n_atoms_;b_beta++){
            for(int gamma=0;gamma<n_orbs_;gamma++){
                int atom_plus_orb=b_beta +  n_atoms_*gamma;
                for(int spin_=0;spin_<2;spin_++){
                    col_ = j + atom_plus_orb*(UnitCellSize_x*UnitCellSize_y) +
                           spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                    for(int spin_p=0;spin_p<2;spin_p++){

                        row_ = j + atom_plus_orb*(UnitCellSize_x*UnitCellSize_y) +
                               spin_p*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                        if(col_>row_){
                            index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val_p = OPs_.value[index_OP_p];
                        }
                        else{
                            index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val_p = conj(OPs_.value[index_OP_p]);
                        }


                        if(spin_==spin_p){
                            row_OP = j + atom_plus_orb*(UnitCellSize_x*UnitCellSize_y)
                            + (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                            col_OP = row_OP;
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];

                            E_class_temp += Parameters_.U0[b_beta]*OPs_.value[index_OP]*OP_val_p;
                        }
                        else{
                            if(!Parameters_.Just_Hartree){
                                row_OP = j + atom_plus_orb*(UnitCellSize_x*UnitCellSize_y)
                                + (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                col_OP = j + atom_plus_orb*(UnitCellSize_x*UnitCellSize_y)
                                         + (spin_p)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                if(!NA_spinflip){
                                    if(col_OP>row_OP){
                                        index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                        OP_val = OPs_.value[index_OP];
                                    }
                                    else{
                                        index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
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
    for(int j=0;j<UnitCellSize_x*UnitCellSize_y;j++){
        for(int atom_no=0;atom_no<n_atoms_;atom_no++){
            for(int alpha=0;alpha<n_orbs_;alpha++){
                for(int spin_=0;spin_<2;spin_++){
                    for(int beta=(alpha+1);beta<n_orbs_;beta++){
                        for(int spin_p=0;spin_p<2;spin_p++){
                            col_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                   spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                            row_=col_;
                            index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val_p = OPs_.value[index_OP_p];

                            row_OP = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y)
                                     + (spin_p)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                            col_OP = row_OP;
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            E_class_temp += (Parameters_.UPrime[atom_no] - 0.5*Parameters_.JHund[atom_no])*OPs_.value[index_OP]*OP_val_p;

                            col_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                   spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                            row_=col_;
                            index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val_p = OPs_.value[index_OP_p];

                            row_OP = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y)
                                     + (spin_p)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                            col_OP = row_OP;
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            E_class_temp += (Parameters_.UPrime[atom_no] - 0.5*Parameters_.JHund[atom_no])*OPs_.value[index_OP]*OP_val_p;

                        }
                    }
                }
            }}}



    //Inter-orbital nn
    //Fock
    for(int j=0;j<UnitCellSize_x*UnitCellSize_y;j++){
        for(int atom_no=0;atom_no<n_atoms_;atom_no++){
            for(int alpha=0;alpha<n_orbs_;alpha++){
                for(int beta=(alpha+1);beta<n_orbs_;beta++){

                    for(int spin_=0;spin_<2;spin_++){

                        //1 and 7 [2 and 8 are h.c. of 1 and 7 ]
                        col_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                               spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        row_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                               spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        if(col_>row_){
                            index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val_p = OPs_.value[index_OP_p];
                        }
                        else{
                            index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val_p = conj(OPs_.value[index_OP_p]);
                        }

                        col_OP = row_;
                        row_OP = col_;

                        if(col_OP>row_OP){
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val = OPs_.value[index_OP];
                        }
                        else{
                            index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val = conj(OPs_.value[index_OP]);
                        }
                        //1 and 7
                        E_class_temp -= (Parameters_.UPrime[atom_no] - 0.5*Parameters_.JHund[atom_no])*OP_val*OP_val_p;
                        //2 and 8
                        E_class_temp -= conj((Parameters_.UPrime[atom_no] - 0.5*Parameters_.JHund[atom_no])*OP_val*OP_val_p);


                        //3 and 5 [4 and 6 are h.c. of 3 and 5 ]
                        col_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                               spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        row_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                               (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                        if(!NA_spinflip){
                            if(col_>row_){
                                index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_val_p = OPs_.value[index_OP_p];
                            }
                            else{
                                index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_val_p = conj(OPs_.value[index_OP_p]);
                            }

                            col_OP = row_;
                            row_OP = col_;

                            if(col_OP>row_OP){
                                index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_val = OPs_.value[index_OP];
                            }
                            else{
                                index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_val = conj(OPs_.value[index_OP]);
                            }
                            //3 and 5
                            E_class_temp -= (Parameters_.UPrime[atom_no] - 0.5*Parameters_.JHund[atom_no])*OP_val*OP_val_p;
                            //4 and 6
                            E_class_temp -= conj((Parameters_.UPrime[atom_no] - 0.5*Parameters_.JHund[atom_no])*OP_val*OP_val_p);
                        }
                    }

                }}}}





    //Hunds Coupling
    //Hartree
    for(int j=0;j<UnitCellSize_x*UnitCellSize_y;j++){
        for(int atom_no=0;atom_no<n_atoms_;atom_no++){
            for(int alpha=0;alpha<n_orbs_;alpha++){
                for(int beta=(alpha+1);beta<n_orbs_;beta++){

                    for(int spin_=0;spin_<2;spin_++){

                            //1 and 3
                        col_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                               spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        row_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                               (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                        if(!NA_spinflip){
                            if(col_>row_){
                                index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_val_p = OPs_.value[index_OP_p];
                            }
                            else{
                                index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_val_p = conj(OPs_.value[index_OP_p]);
                            }


                            col_OP = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                     (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                            row_OP = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                     (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                            if(col_OP>row_OP){
                                index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_val = OPs_.value[index_OP];
                            }
                            else{
                                index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_val = conj(OPs_.value[index_OP]);
                            }
                            //1 and 3
                            E_class_temp -= (0.5*Parameters_.JHund[atom_no])*2.0*OP_val*OP_val_p;
                        }

                        //2 and 4
                        col_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                               spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        row_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                               (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                        if(!NA_spinflip){
                            if(col_>row_){
                                index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_val_p = OPs_.value[index_OP_p];
                            }
                            else{
                                index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_val_p = conj(OPs_.value[index_OP_p]);
                            }

                            col_OP = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                     (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                            row_OP = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                     (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                            if(col_OP>row_OP){
                                index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_val = OPs_.value[index_OP];
                            }
                            else{
                                index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_val = conj(OPs_.value[index_OP]);
                            }
                            //2 and 4
                            E_class_temp -= (0.5*Parameters_.JHund[atom_no])*2.0*OP_val*OP_val_p;
                        }

                        //5 and 11
                        col_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                               spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        row_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                               (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        if(col_>row_){
                            index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val_p = OPs_.value[index_OP_p];
                        }
                        else{
                            index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val_p = conj(OPs_.value[index_OP_p]);
                        }

                        col_OP = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                 (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        row_OP = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                 (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                        if(col_OP>row_OP){
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val = OPs_.value[index_OP];
                        }
                        else{
                            index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val = conj(OPs_.value[index_OP]);
                        }
                        //5 and 11
                        E_class_temp -= (0.5*Parameters_.JHund[atom_no])*OP_val*OP_val_p;


                        //6 and 12
                        col_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                               spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        row_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                               (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        if(col_>row_){
                            index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val_p = OPs_.value[index_OP_p];
                        }
                        else{
                            index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val_p = conj(OPs_.value[index_OP_p]);
                        }

                        col_OP = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                 (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        row_OP = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                 (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                        if(col_OP>row_OP){
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val = OPs_.value[index_OP];
                        }
                        else{
                            index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val = conj(OPs_.value[index_OP]);
                        }
                        //6 and 12
                        E_class_temp -= (0.5*Parameters_.JHund[atom_no])*OP_val*OP_val_p;



                        //7 and 9
                        col_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                               spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        row_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                               (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        if(col_>row_){
                            index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val_p = OPs_.value[index_OP_p];
                        }
                        else{
                            index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val_p = conj(OPs_.value[index_OP_p]);
                        }

                        col_OP = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                 (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        row_OP = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                 (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                        if(col_OP>row_OP){
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val = OPs_.value[index_OP];
                        }
                        else{
                            index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val = conj(OPs_.value[index_OP]);
                        }
                        //7 and 9
                        E_class_temp += (0.5*Parameters_.JHund[atom_no])*OP_val*OP_val_p;


                        //8 and 10
                        col_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                               spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        row_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                               (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        if(col_>row_){
                            index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val_p = OPs_.value[index_OP_p];
                        }
                        else{
                            index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val_p = conj(OPs_.value[index_OP_p]);
                        }

                        col_OP = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                 (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        row_OP = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                 (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                        if(col_OP>row_OP){
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val = OPs_.value[index_OP];
                        }
                        else{
                            index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val = conj(OPs_.value[index_OP]);
                        }
                        //8 and 10
                        E_class_temp += (0.5*Parameters_.JHund[atom_no])*OP_val*OP_val_p;

                    }

                }}}}





    //Hunds Coupling
    //Fock
    for(int j=0;j<UnitCellSize_x*UnitCellSize_y;j++){
        for(int atom_no=0;atom_no<n_atoms_;atom_no++){
            for(int alpha=0;alpha<n_orbs_;alpha++){
                for(int beta=(alpha+1);beta<n_orbs_;beta++){

                    for(int spin_=0;spin_<2;spin_++){

                        //1 and 3 [hc : 2 and 4]
                        col_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                               spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        row_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                               (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        if(col_>row_){
                            index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val_p = OPs_.value[index_OP_p];
                        }
                        else{
                            index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val_p = conj(OPs_.value[index_OP_p]);
                        }

                        col_OP = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                 (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        row_OP = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                 (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                        if(col_OP>row_OP){
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val = OPs_.value[index_OP];
                        }
                        else{
                            index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val = conj(OPs_.value[index_OP]);
                        }
                        //1 and 3
                        E_class_temp += (0.5*Parameters_.JHund[atom_no])*2.0*OP_val*OP_val_p;
                        //2 and 4
                        E_class_temp += conj((0.5*Parameters_.JHund[atom_no])*2.0*OP_val*OP_val_p);



                        //5 and 11 [hc : 6 and 12]
                        col_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                               spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        row_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                               (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        if(col_>row_){
                            index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val_p = OPs_.value[index_OP_p];
                        }
                        else{
                            index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val_p = conj(OPs_.value[index_OP_p]);
                        }

                        col_OP = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                 (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        row_OP = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                 (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                        if(col_OP>row_OP){
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val = OPs_.value[index_OP];
                        }
                        else{
                            index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val = conj(OPs_.value[index_OP]);
                        }
                        //5 and 11
                        E_class_temp += (0.5*Parameters_.JHund[atom_no])*OP_val*OP_val_p;
                        //6 and 12
                        E_class_temp += conj((0.5*Parameters_.JHund[atom_no])*OP_val*OP_val_p);


                        //7 and 9 [hc : 8 and 10]
                        col_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                               spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        row_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                               (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        if(!NA_spinflip){
                            if(col_>row_){
                                index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_val_p = OPs_.value[index_OP_p];
                            }
                            else{
                                index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_val_p = conj(OPs_.value[index_OP_p]);
                            }

                            col_OP = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                     (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                            row_OP = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                     (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                            if(col_OP>row_OP){
                                index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_val = OPs_.value[index_OP];
                            }
                            else{
                                index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
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
    for(int j=0;j<UnitCellSize_x*UnitCellSize_y;j++){
        for(int atom_no=0;atom_no<n_atoms_;atom_no++){
            for(int alpha=0;alpha<n_orbs_;alpha++){
                for(int beta=(alpha+1);beta<n_orbs_;beta++){

                    for(int spin_=0;spin_<2;spin_++){

                        //1 and 2 [hc : 3 and 4]
                        col_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                               spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        row_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                               (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        if(!NA_spinflip){
                            if(col_>row_){
                                index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_val_p = OPs_.value[index_OP_p];
                            }
                            else{
                                index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_val_p = conj(OPs_.value[index_OP_p]);
                            }

                            col_OP = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                     (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                            row_OP = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                     (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                            if(col_OP>row_OP){
                                index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_val = OPs_.value[index_OP];
                            }
                            else{
                                index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_val = conj(OPs_.value[index_OP]);
                            }
                            //1 and 2
                            E_class_temp -= (1.0*Parameters_.JHund[atom_no])*OP_val*OP_val_p;
                            //3 and 4
                            E_class_temp -= conj((1.0*Parameters_.JHund[atom_no])*OP_val*OP_val_p);
                        }

                        //5 and 6 [hc : 7 and 8]
                        col_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                               spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        row_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                               (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        if(col_>row_){
                            index_OP_p = SI_to_ind[col_ + row_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val_p = OPs_.value[index_OP_p];
                        }
                        else{
                            index_OP_p = SI_to_ind[row_ + col_*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val_p = conj(OPs_.value[index_OP_p]);
                        }

                        col_OP = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                 (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        row_OP = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                 (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                        if(col_OP>row_OP){
                            index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            OP_val = OPs_.value[index_OP];
                        }
                        else{
                            index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
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

void Kspace_calculation_G2dLattice::Get_Energies(){


    complex<double> E_class_temp=0.0;
    E_class= 0.0;

    E_class_onsite_U0_Hartree=0.0; E_class_longrange_Hartree=0.0;
    E_class_onsite_U0_Fock=0.0; E_class_longrange_Fock=0.0;


    int row_, col_, index, rowM_, colM_;
    int row_OP, col_OP, index_OP;

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

    for(int jpx=0;jpx<UnitCellSize_x;jpx++){
        for(int jpy=0;jpy<UnitCellSize_y;jpy++){
            for(int b_beta_p=0;b_beta_p<n_atoms_*n_orbs_;b_beta_p++){

                org_site_npx=0*UnitCellSize_x + jpx;
                org_site_npy=0*UnitCellSize_y + jpy;
                colM_= b_beta_p + org_site_npx*(n_atoms_*n_orbs_) + org_site_npy*(lx_*n_atoms_*n_orbs_);


                //Hartree
                col_temp_up = (jpx + UnitCellSize_x*jpy) + b_beta_p*(UnitCellSize_x*UnitCellSize_y) +
                              UP_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y) +
                              0*(2*n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                row_temp_up = col_temp_up;
                col_temp_dn = (jpx + UnitCellSize_x*jpy) + b_beta_p*(UnitCellSize_x*UnitCellSize_y) +
                              DOWN_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y) +
                              0*(2*n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                row_temp_dn =col_temp_dn;
                comp_temp_up = SI_to_ind[col_temp_up + row_temp_up*(ncells_*2*n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y)];
                comp_temp_dn = SI_to_ind[col_temp_dn + row_temp_dn*(ncells_*2*n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y)];

                E_class_temp -= 1.0*(ncells_)*Parameters_.U0[b_beta_p]*(OPs_.value[comp_temp_up]*OPs_.value[comp_temp_dn]);


                //Fock
                if(!Parameters_.Just_Hartree){
                    int comp_Splus = SI_to_ind[col_temp_dn + row_temp_up*(ncells_*2*n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y)];
                    E_class_temp += 1.0*(ncells_)*Parameters_.U0[b_beta_p]*(OPs_.value[comp_Splus]*conj(OPs_.value[comp_Splus]));
                }
            }
        }
    }










    E_class = E_class_temp.real();

    //----------Eclass done


    E_quant=0.0;
    for(int n=0;n<Eigenvalues_.size();n++){
        E_quant += (Eigenvalues_[n]*
                    (1.0/( exp((Eigenvalues_[n]-mu_)*Parameters_.beta ) + 1.0))
                    );
    }




}

void Kspace_calculation_G2dLattice::Get_new_OPs_and_error(){

    //For <c_{c1}* c_{c2}>
    int c1;
    int c2;
    int S_=UnitCellSize_x*UnitCellSize_y;
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

void Kspace_calculation_G2dLattice::Get_spin_resolved_local_densities(){

    //For <c_{c1}* c_{c2}>
    int c1;
    int c2;
    int S_=UnitCellSize_x*UnitCellSize_y;

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
    int site_x, site_y, alpha, site;
    for(int alpha_1=0;alpha_1<UnitCellSize_x;alpha_1++){
        for(int cell_1=0;cell_1<lx_cells;cell_1++){
            site_x = (cell_1*UnitCellSize_x) + alpha_1;

            for(int alpha_2=0;alpha_2<UnitCellSize_y;alpha_2++){
                for(int cell_2=0;cell_2<ly_cells;cell_2++){
                    site_y = (cell_2*UnitCellSize_y) + alpha_2;

                    alpha = alpha_1 + alpha_2*(UnitCellSize_x);
                    site = site_x + site_y*(lx_);

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
    }



    cout<< "Total electrons in Lattice, UP = "<<Total_den_up<<endl;
    cout<< "Total electrons in Lattice, DOWN = "<<Total_den_dn<<endl;


}



void Kspace_calculation_G2dLattice::Update_Total_Density(){

    //For sum_c1 <c_{c1}* c_{c1}>

    double val;

    val=0.0;

    for(int state=0;state<Eigenvalues_.size();state++){
        val += (1.0/( exp((Eigenvalues_[state]-mu_)*Parameters_.beta ) + 1.0));
    }

    Parameters_.Total_Particles = val;


}


void Kspace_calculation_G2dLattice::Get_Tau_Pseudospins(){

    //For <c_{c1}* c_{c2}>
    assert(n_orbs_>=2);

    double a_NN;
    a_NN=1.0;
    int c1;
    int c2;
    int S_=UnitCellSize_x*UnitCellSize_y;
    string atomorb_name;

    int k_index;
    int state;

    double Total_Sz, Total_Sx, Total_Sy, Total_den;
    int UP_, DOWN_;
    UP_=0;DOWN_=1;



    char temp_char[50];
    sprintf(temp_char, "%.10f", Parameters_.Temperature);

    string File_Out_Local_orb_densities = "Local_TauPseudospins"  + string(temp_char) +  ".txt";
    ofstream file_out_Local_orb_densities(File_Out_Local_orb_densities.c_str());
    file_out_Local_orb_densities<<"#site (in real lattice) site_1 site_2  rx  ry   Tz   Tx  Ty  |T|"<<endl;

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
    for(int site_x=0;site_x<lx_cells*UnitCellSize_x;site_x++){
        int alpha_1 = site_x%UnitCellSize_x;
        int cell_1 = int((site_x - alpha_1 +0.5)/UnitCellSize_x);

        // for(int alpha_2=0;alpha_2<UnitCellSize_y;alpha_2++){
        //     for(int cell_2=0;cell_2<ly_cells;cell_2++){
        //         site_y = (cell_2*UnitCellSize_y) + alpha_2;
        for(int site_y=0;site_y<lx_cells*UnitCellSize_x;site_y++){
            int alpha_2 = site_y%UnitCellSize_y;
            int cell_2 = int((site_y - alpha_2 +0.5)/UnitCellSize_y);

            for(int gamma=0;gamma<n_atoms_;gamma++){

                int subltc_;
                subltc_= gamma;


                atomorb_name="subltc_"+to_string(subltc_);

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

                alpha = alpha_1 + alpha_2*(UnitCellSize_x);
                site = site_x + site_y*(lx_);

                file_out_Local_orb_densities<<site<<"    "<<site_x<<"    "<<site_y<<"    "<<site<<atomorb_name<<"   "<<rx_<<"     "<<ry_<<"     ";


                //Local Tplus
                val=0.0;
                for(int spin_=0;spin_<2;spin_++){
                    c1 = alpha + (gamma+0*n_atoms_)*(S_) + spin_*(n_atoms_*n_orbs_*S_);
                    c2 = alpha + (gamma+1*n_atoms_)*(S_)  + spin_*(n_atoms_*n_orbs_*S_);
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
                }
                splus_val = val;
                sx_val = splus_val.real();
                sy_val = splus_val.imag();

                //Local Tauz
                val=0.0;
                double fac_;
                for(int sigma=0;sigma<2;sigma++){
                    for(int orb=0;orb<2;orb++){
                        fac_=1.0 - (2.0*orb);
                        c1 = alpha + (gamma+orb*n_atoms_)*(S_) + sigma*(n_orbs_*n_atoms_*S_);
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
                }
                sz_val=val.real();



                file_out_Local_orb_densities<<sz_val<<"    "<<sx_val<<"   "<<sy_val<<"    "<<sqrt(sz_val*sz_val + sy_val*sy_val + sx_val*sx_val)<<"     ";


                Total_Sz +=sz_val;
                Total_Sy +=sy_val;
                Total_Sx +=sx_val;



                file_out_Local_orb_densities<<endl;
            }
        }
        file_out_Local_orb_densities<<endl;

    }



    cout<< "Total Tauz(orb=0,1) = "<<Total_Sz<<endl;
    cout<< "Total Tx = "<<Total_Sx<<endl;
    cout<< "Total Ty = "<<Total_Sy<<endl;
    cout<< "|Total_T| = "<<sqrt(Total_Sz*Total_Sz + Total_Sx*Total_Sx + Total_Sy*Total_Sy)<<endl;




}


void Kspace_calculation_G2dLattice::Get_local_spins(){

    //For <c_{c1}* c_{c2}>


    double a_NN;
    a_NN=1.0;
    int c1;
    int c2;
    int S_=UnitCellSize_x*UnitCellSize_y;
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
    for(int site_x=0;site_x<lx_cells*UnitCellSize_x;site_x++){
        int alpha_1 = site_x%UnitCellSize_x;
        int cell_1 = int((site_x - alpha_1 +0.5)/UnitCellSize_x);

        // for(int alpha_2=0;alpha_2<UnitCellSize_y;alpha_2++){
        //     for(int cell_2=0;cell_2<ly_cells;cell_2++){
        //         site_y = (cell_2*UnitCellSize_y) + alpha_2;
        for(int site_y=0;site_y<lx_cells*UnitCellSize_x;site_y++){
            int alpha_2 = site_y%UnitCellSize_y;
            int cell_2 = int((site_y - alpha_2 +0.5)/UnitCellSize_y);

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

                alpha = alpha_1 + alpha_2*(UnitCellSize_x);
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

double Kspace_calculation_G2dLattice::chemicalpotential(double Particles){


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

double Kspace_calculation_G2dLattice::random1(){

    return dis1_(Generator1_);

}


void Kspace_calculation_G2dLattice::Initialize()
{


    kick_while_cooling=0.0;
    
    NOT_AVAIL_INT=-1000;

    OP_only_finite_Int=true;
    Global_Eps=0.0000001;

    cout<<"Parameters_.FockType = '" << Parameters_.FockType<<"'"<<endl;
    ly_ = Parameters_.ly;
    lx_ = Parameters_.lx;
    UnitCellSize_x = Parameters_.UnitCellSize_x;
    UnitCellSize_y = Parameters_.UnitCellSize_y;

    lx_cells = lx_ / UnitCellSize_x;
    ly_cells = ly_ / UnitCellSize_y;

    ncells_ = lx_cells * ly_cells;
    n_orbs_ = Parameters_.n_orbs;
    n_atoms_ = Parameters_.n_atoms;


    //assert(n_orbs_==1);
    //assert(n_atoms_==3);

    Ham_.resize(2*n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y, 2*n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y);
    Eigvectors_.resize(2*n_orbs_*n_atoms_*lx_*ly_);
    Eigenvalues_.resize(2*n_orbs_*n_atoms_*lx_*ly_);
    for(int i=0;i<2*n_atoms_*n_orbs_*lx_*ly_;i++){
        Eigvectors_[i].resize(2*n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y); //Eigenvector number
    }

    Kx_values.resize(ncells_);
    Ky_values.resize(ncells_);


    // Create_V_int();
    // Create_B_and_C_mat();


    //Order Parameters:
    /*
    Inter-unit cell OP's ( intra-unit cell when d_vec=0)
    <c_{alpha, gamma, sigma}^{dagger} c _{ d_vec, alpha', gamma', sigma' }> ; for all d_vec,
    and alpha+ gamma*(UnitCellSize_x*UnitCellSize_y) + sigma*(n_orbs_*UnitCellSize_x*UnitCellSize_y)
       <= alpha'+ gamma'*(UnitCellSize_x*UnitCellSize_y) + sigma'*(n_orbs_*UnitCellSize_x*UnitCellSize_y)

    For <c_{alpha, gamma, sigma}^{dagger} c _{ d_vec, alpha',gamma', sigma' }> ,
   when alpha+ gamma*(UnitCellSize_x*UnitCellSize_y) + sigma*(n_orbs_*UnitCellSize_x*UnitCellSize_y)
       > alpha'+ gamma'*(UnitCellSize_x*UnitCellSize_y) + sigma'*(n_orbs_*UnitCellSize_x*UnitCellSize_y)
    use
    <c_{alpha, gamma, sigma}^{dagger} c _{ d_vec, alpha',gamma', sigma' }>
    = conj(<c_{alpha', gamma', sigma'}^{dagger} c _{ -d_vec, alpha,gamma, sigma}> )

    P=2*n_orbs_*UnitCellSize_x*UnitCellSize_y

    Total no of Order parameters (Hartree + all Fock):
    ncells_*(P(P+1))/2

    Total no of Order parameters (Hartree + only intracell Fock):
    (P(P+1))/2

    Total no of Order parameters (Hartree):
    P
    */

    int S_=UnitCellSize_x*UnitCellSize_y; //no. of sites in single unit cell
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
            for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
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
                for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                    alpha_1 = alpha % UnitCellSize_x;
                    alpha_2 = (alpha - alpha_1)/UnitCellSize_x;

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

        }

    }

    else{
    }



    string fl_initial_OP_out = "OP_initial.txt";
    ofstream file_initial_OP_out(fl_initial_OP_out.c_str());
    for(int i=0;i<OPs_.value.size();i++){
        file_initial_OP_out<<i<<"  "<<OPs_.rows[i]<<"  "<<OPs_.columns[i]<<"  "<<OPs_.value[i]<<endl;
    }




    cout<<"Total no. of OP's = "<<OPs_.value.size()<<endl;



}


void Kspace_calculation_G2dLattice::Arranging_spectrum(){

    // Eigvectors_saved=Eigvectors_;
    Eigenvalues_saved = Eigenvalues_;

    //    Eigvectors_.resize(ncells_*6);
    //    Eigenvalues_.resize(ncells_*6);
    //    for(int i=0;i<ncells_*6;i++){
    //        Eigvectors_[i].resize(6); //Eigenvector number
    //    }

    double value_;
    //    int S_=2*UnitCellSize_x*UnitCellSize_y;
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

complex<double> Kspace_calculation_G2dLattice::h_KE(int alpha, int gamma, int sigma, int alpha_p, int gamma_p, int sigma_p, int k1, int k2){


    complex<double> temp_val;
    int row_, col_;
    int d1_, d2_;
    int d1_org, d2_org;
    int alpha_1, alpha_2, alpha_p_1, alpha_p_2;
    alpha_1 = alpha % UnitCellSize_x; //a = a1 +a2*lx
    alpha_p_1 = alpha_p % UnitCellSize_x;

    alpha_2 = (alpha - alpha_1)/UnitCellSize_x;
    alpha_p_2 = (alpha_p - alpha_p_1)/UnitCellSize_x;

    temp_val=0.0;
    for(int cell_no=0;cell_no<ncells_;cell_no++){
        d1_org = Coordinates_.indx_cellwise(cell_no);
        d2_org = Coordinates_.indy_cellwise(cell_no);

        row_ = gamma  + ( (d1_org*UnitCellSize_x) + alpha_1)*n_orbs_*n_atoms_
               + (( (d2_org*UnitCellSize_y) + alpha_2)*lx_*n_orbs_*n_atoms_)
               + sigma*(lx_*ly_*n_orbs_*n_atoms_);
        col_ = ( gamma_p + (0 + alpha_p_1)*n_orbs_*n_atoms_ + ((0 + alpha_p_2)*lx_*n_orbs_*n_atoms_)) + sigma_p*(lx_*ly_*n_orbs_*n_atoms_);

        Get_minimum_distance_direction(0, cell_no, d1_, d2_);

        temp_val += Connections_.HTB_(row_,col_)*exp(iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_   ));

    }


    return temp_val;
}



void Kspace_calculation_G2dLattice::Create_Kspace_Spectrum(){


    bool PinningField=true;
    double PinningFieldValue=Parameters_.PinningFieldValue;


    //Calculating density using OP's will be going inside Hamiltonian
    OPs_total_den=0.0;
    int index_OP;
    int row_OP, col_OP;
    
    for(int k1=0;k1<lx_cells;k1++){
        for(int k2=0;k2<ly_cells;k2++){
            for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                for(int gamma=0;gamma<n_atoms_*n_orbs_;gamma++){
                    for(int spin=0;spin<2;spin++){
                        row_OP = alpha + gamma*(UnitCellSize_x*UnitCellSize_y)
                        + spin*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        col_OP = alpha + gamma*(UnitCellSize_x*UnitCellSize_y)
                                 + spin*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
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
            for(int alpha_p=0;alpha_p<UnitCellSize_x*UnitCellSize_y;alpha_p++){
                for(int gamma_p=0;gamma_p<n_orbs_*n_atoms_;gamma_p++){
                    for(int sigma_p=0;sigma_p<2;sigma_p++){
                        col_ = alpha_p + gamma_p*(UnitCellSize_x*UnitCellSize_y)
                        + sigma_p*(n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y);

                        for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                            for(int gamma=0;gamma<n_orbs_*n_atoms_;gamma++){
                                for(int sigma=0;sigma<2;sigma++){
                                    row_ = alpha + gamma*(UnitCellSize_x*UnitCellSize_y)
                                    + sigma*(n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y);

                                    Ham_(row_, col_) += h_KE(alpha, gamma, sigma, alpha_p, gamma_p, sigma_p, k1, k2);
                                }
                            }
                        }
                    }
                }
            }





            /*
            //Anisotropy
            for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                for(int spin1=0;spin1<2;spin1++){
                    for(int spin2=0;spin2<2;spin2++){

                        row_ = alpha + spin1*(UnitCellSize_x*UnitCellSize_y);
                        col_=row_;

                        row_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin2*(UnitCellSize_x*UnitCellSize_y);
                        col_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin2*(UnitCellSize_x*UnitCellSize_y);
                        index_OP = SI_to_ind[col_OP + row_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];

                        fac_ = 1.0 - (2.0*abs(spin1-spin2)); //i.e 0 --> 1, 1 --> -1

                        Ham_(row_,col_) += (-1.0/2.0)*Parameters_.AnisotropyZ*fac_*OPs_.value[index_OP];
                    }
                }
            }
            */


            //Onsite Energies
            //OnSiteE_up
            for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                for(int gamma=0;gamma<n_atoms_*n_orbs_;gamma++){
                    for(int spin1=0;spin1<2;spin1++){
                        row_ = alpha + gamma*(UnitCellSize_x*UnitCellSize_y)
                        + spin1*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                        col_=row_;
                        Ham_(row_,col_) += Parameters_.OnSiteE[gamma][spin1];
                    }
                }
            }


            //small pinning field
            if(PinningField){
                int alpha=0;
                int gamma=0;
                for(int spin1=0;spin1<2;spin1++){
                    row_ = alpha + gamma*(UnitCellSize_x*UnitCellSize_y)
                    + spin1*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                    col_=row_;
                    Ham_(row_,col_) += PinningFieldValue*(spin1);
                }

            }

            bool NA_spinflip=Parameters_.NoSpinFlipOP;
            complex<double> OP_val;
            //Interaction local intra-orbital U0:
            for(int j=0;j<UnitCellSize_x*UnitCellSize_y;j++){
                for(int b_beta=0;b_beta<n_atoms_;b_beta++){
                    for(int gamma=0;gamma<n_orbs_;gamma++){
                        int atom_plus_orb=b_beta +  n_atoms_*gamma;
                        for(int spin_=0;spin_<2;spin_++){
                            col_ = j + atom_plus_orb*(UnitCellSize_x*UnitCellSize_y) +
                                   spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                            for(int spin_p=0;spin_p<2;spin_p++){

                                row_ = j + atom_plus_orb*(UnitCellSize_x*UnitCellSize_y) +
                                       spin_p*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                if(spin_==spin_p){
                                    row_OP = j + atom_plus_orb*(UnitCellSize_x*UnitCellSize_y)
                                    + (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                    col_OP = row_OP;
                                    index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];

                                    Ham_(row_,col_) += Parameters_.U0[b_beta]*OPs_.value[index_OP];
                                }
                                else{
                                    if(!Parameters_.Just_Hartree){
                                        if(!NA_spinflip){
                                            row_OP = j + atom_plus_orb*(UnitCellSize_x*UnitCellSize_y)
                                            + (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                            col_OP = j + atom_plus_orb*(UnitCellSize_x*UnitCellSize_y)
                                                     + (spin_p)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                            if(col_OP>row_OP){
                                                index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                                OP_val = OPs_.value[index_OP];
                                            }
                                            else{
                                                index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
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
            for(int j=0;j<UnitCellSize_x*UnitCellSize_y;j++){
                for(int atom_no=0;atom_no<n_atoms_;atom_no++){
                    for(int alpha=0;alpha<n_orbs_;alpha++){
                        for(int spin_=0;spin_<2;spin_++){
                            for(int beta=(alpha+1);beta<n_orbs_;beta++){
                                for(int spin_p=0;spin_p<2;spin_p++){
                                    col_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                           spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                    row_=col_;
                                    row_OP = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y)
                                             + (spin_p)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                    col_OP = row_OP;
                                    index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                    Ham_(row_,col_) += (Parameters_.UPrime[atom_no] - 0.5*Parameters_.JHund[atom_no])*OPs_.value[index_OP];

                                    col_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                           spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                    row_=col_;
                                    row_OP = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y)
                                             + (spin_p)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                    col_OP = row_OP;
                                    index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                    Ham_(row_,col_) += (Parameters_.UPrime[atom_no] - 0.5*Parameters_.JHund[atom_no])*OPs_.value[index_OP];

                                }
                            }
                        }
                    }}}



            //Inter-orbital nn
            //Fock
            for(int j=0;j<UnitCellSize_x*UnitCellSize_y;j++){
                for(int atom_no=0;atom_no<n_atoms_;atom_no++){
                    for(int alpha=0;alpha<n_orbs_;alpha++){
                        for(int beta=(alpha+1);beta<n_orbs_;beta++){

                            for(int spin_=0;spin_<2;spin_++){

                                //1 and 7 [2 and 8 are h.c. of 1 and 7 ]
                                col_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                       spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                row_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                       spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                col_OP = row_;
                                row_OP = col_;

                                if(col_OP>row_OP){
                                    index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                    OP_val = OPs_.value[index_OP];
                                }
                                else{
                                    index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                    OP_val = conj(OPs_.value[index_OP]);
                                }
                                //1 and 7
                                Ham_(row_,col_) -= (Parameters_.UPrime[atom_no] - 0.5*Parameters_.JHund[atom_no])*OP_val;
                                //2 and 8
                                Ham_(col_,row_) -= conj((Parameters_.UPrime[atom_no] - 0.5*Parameters_.JHund[atom_no])*OP_val);


                                //3 and 5 [4 and 6 are h.c. of 3 and 5 ]
                                col_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                       spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                row_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                       (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                col_OP = row_;
                                row_OP = col_;


                                if(!NA_spinflip){
                                    if(col_OP>row_OP){
                                        index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                        OP_val = OPs_.value[index_OP];
                                    }
                                    else{
                                        index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                        OP_val = conj(OPs_.value[index_OP]);
                                    }
                                    //3 and 5
                                    Ham_(row_,col_) -= (Parameters_.UPrime[atom_no] - 0.5*Parameters_.JHund[atom_no])*OP_val;
                                    //4 and 6
                                    Ham_(col_,row_) -= conj((Parameters_.UPrime[atom_no] - 0.5*Parameters_.JHund[atom_no])*OP_val);
                                }
                            }

                        }}}}


            //Hunds Coupling
            //Hartree
            for(int j=0;j<UnitCellSize_x*UnitCellSize_y;j++){
                for(int atom_no=0;atom_no<n_atoms_;atom_no++){
                    for(int alpha=0;alpha<n_orbs_;alpha++){
                        for(int beta=(alpha+1);beta<n_orbs_;beta++){

                            for(int spin_=0;spin_<2;spin_++){

                                    //1 and 3
                                col_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                       spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                row_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                       (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                col_OP = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                         (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                row_OP = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                         (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                //NA_spinflip = Parameters_.NoSpinFlipOP && (spin_p!=spin_);
                                if(!NA_spinflip){
                                    if(col_OP>row_OP){
                                        index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                        OP_val = OPs_.value[index_OP];
                                    }
                                    else{
                                        index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                        OP_val = conj(OPs_.value[index_OP]);
                                    }
                                    //1 and 3
                                    Ham_(row_,col_) -= (0.5*Parameters_.JHund[atom_no])*2.0*OP_val;
                                }


                                //2 and 4
                                col_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                       spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                row_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                       (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                col_OP = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                         (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                row_OP = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                         (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                //NA_spinflip = Parameters_.NoSpinFlipOP && (spin_p!=spin_);
                                if(!NA_spinflip){
                                    if(col_OP>row_OP){
                                        index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                        OP_val = OPs_.value[index_OP];
                                    }
                                    else{
                                        index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                        OP_val = conj(OPs_.value[index_OP]);
                                    }
                                    //2 and 4
                                    Ham_(row_,col_) -= (0.5*Parameters_.JHund[atom_no])*2.0*OP_val;
                                }

                                //5 and 11
                                col_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                       spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                row_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                       (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                col_OP = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                         (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                row_OP = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                         (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                if(col_OP>row_OP){
                                    index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                    OP_val = OPs_.value[index_OP];
                                }
                                else{
                                    index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                    OP_val = conj(OPs_.value[index_OP]);
                                }
                                //5 and 11
                                Ham_(row_,col_) -= (0.5*Parameters_.JHund[atom_no])*OP_val;


                                //6 and 12
                                col_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                       spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                row_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                       (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                col_OP = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                         (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                row_OP = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                         (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                if(col_OP>row_OP){
                                    index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                    OP_val = OPs_.value[index_OP];
                                }
                                else{
                                    index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                    OP_val = conj(OPs_.value[index_OP]);
                                }
                                //6 and 12
                                Ham_(row_,col_) -= (0.5*Parameters_.JHund[atom_no])*OP_val;



                                //7 and 9
                                col_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                       spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                row_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                       (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                col_OP = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                         (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                row_OP = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                         (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                if(col_OP>row_OP){
                                    index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                    OP_val = OPs_.value[index_OP];
                                }
                                else{
                                    index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                    OP_val = conj(OPs_.value[index_OP]);
                                }
                                //7 and 9
                                Ham_(row_,col_) += (0.5*Parameters_.JHund[atom_no])*OP_val;


                                //8 and 10
                                col_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                       spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                row_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                       (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                col_OP = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                         (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                row_OP = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                         (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                if(col_OP>row_OP){
                                    index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                    OP_val = OPs_.value[index_OP];
                                }
                                else{
                                    index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                    OP_val = conj(OPs_.value[index_OP]);
                                }
                                //8 and 10
                                Ham_(row_,col_) += (0.5*Parameters_.JHund[atom_no])*OP_val;

                            }

                        }}}}





            //Hunds Coupling
            //Fock
            for(int j=0;j<UnitCellSize_x*UnitCellSize_y;j++){
                for(int atom_no=0;atom_no<n_atoms_;atom_no++){
                    for(int alpha=0;alpha<n_orbs_;alpha++){
                        for(int beta=(alpha+1);beta<n_orbs_;beta++){

                            for(int spin_=0;spin_<2;spin_++){

                                //1 and 3 [hc : 2 and 4]
                                col_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                       spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                row_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                       (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                col_OP = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                         (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                row_OP = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                         (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                if(col_OP>row_OP){
                                    index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                    OP_val = OPs_.value[index_OP];
                                }
                                else{
                                    index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                    OP_val = conj(OPs_.value[index_OP]);
                                }
                                //1 and 3
                                Ham_(row_,col_) += (0.5*Parameters_.JHund[atom_no])*2.0*OP_val;
                                //2 and 4
                                Ham_(col_,row_) += conj((0.5*Parameters_.JHund[atom_no])*2.0*OP_val);



                                //5 and 11 [hc : 6 and 12]
                                col_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                       spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                row_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                       (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                col_OP = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                         (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                row_OP = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                         (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                if(col_OP>row_OP){
                                    index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                    OP_val = OPs_.value[index_OP];
                                }
                                else{
                                    index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                    OP_val = conj(OPs_.value[index_OP]);
                                }
                                //5 and 11
                                Ham_(row_,col_) += (0.5*Parameters_.JHund[atom_no])*OP_val;
                                //6 and 12
                                Ham_(col_,row_) += conj((0.5*Parameters_.JHund[atom_no])*OP_val);


                                //7 and 9 [hc : 8 and 10]
                                col_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                       spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                row_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                       (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                col_OP = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                         (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                row_OP = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                         (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                //NA_spinflip = Parameters_.NoSpinFlipOP && (spin_p!=spin_);
                                if(!NA_spinflip){
                                    if(col_OP>row_OP){
                                        index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                        OP_val = OPs_.value[index_OP];
                                    }
                                    else{
                                        index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
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
            for(int j=0;j<UnitCellSize_x*UnitCellSize_y;j++){
                for(int atom_no=0;atom_no<n_atoms_;atom_no++){
                    for(int alpha=0;alpha<n_orbs_;alpha++){
                        for(int beta=(alpha+1);beta<n_orbs_;beta++){

                            for(int spin_=0;spin_<2;spin_++){

                                //1 and 2 [hc : 3 and 4]
                                col_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                       spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                row_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                       (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                col_OP = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                         (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                row_OP = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                         (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                //NA_spinflip = Parameters_.NoSpinFlipOP && (spin_p!=spin_);
                                if(!NA_spinflip){
                                    if(col_OP>row_OP){
                                        index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                        OP_val = OPs_.value[index_OP];
                                    }
                                    else{
                                        index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                        OP_val = conj(OPs_.value[index_OP]);
                                    }
                                    //1 and 2
                                    Ham_(row_,col_) -= (1.0*Parameters_.JHund[atom_no])*OP_val;
                                    //3 and 4
                                    Ham_(col_,row_) -= conj((1.0*Parameters_.JHund[atom_no])*OP_val);
                                }

                                //5 and 6 [hc : 7 and 8]
                                col_ = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                       spin_*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                row_ = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                       (spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                col_OP = j + (atom_no + n_atoms_*beta)*(UnitCellSize_x*UnitCellSize_y) +
                                         (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);
                                row_OP = j + (atom_no + n_atoms_*alpha)*(UnitCellSize_x*UnitCellSize_y) +
                                         (1-spin_)*(n_atoms_*n_orbs_*UnitCellSize_x*UnitCellSize_y);

                                if(col_OP>row_OP){
                                    index_OP = SI_to_ind[col_OP + row_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                    OP_val = OPs_.value[index_OP];
                                }
                                else{
                                    index_OP = SI_to_ind[row_OP + col_OP*(2*n_atoms_*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
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


            for(int row=0;row<2*n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y;row++){
                Eigenvalues_[2*n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y*k_index + row]=eigs_[row];
                for(int col=0;col<2*n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y;col++){
                    Eigvectors_[2*n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y*k_index + col][row]=Ham_(row,col);
                }
            }
        }
    }

    

}


void Kspace_calculation_G2dLattice::SelfConsistency(){


    string File_Out_Local_OP_Initial = "Initial_" + Parameters_.File_OPs_out + ".txt";
    ofstream file_out_Local_OP_Initial(File_Out_Local_OP_Initial.c_str());
    file_out_Local_OP_Initial<<"#row col OParams_[row][col]"<<endl;

    for(int alpha=0;alpha<OPs_.value.size();alpha++){
        file_out_Local_OP_Initial<<OPs_.rows[alpha]<<setw(15)<<OPs_.columns[alpha]<<setw(15)<<"  "<<OPs_.value[alpha]<<endl;
    }


    //assert(false);

    if(UnitCellSize_x==lx_ && UnitCellSize_y==ly_){
        //---------For Real space--------------------------
        string File_Out_Local_OP_Initial_RS = "Initial_" + Parameters_.File_OPs_out + "_for_RealSpace.txt";
        ofstream file_out_Local_OP_Initial_RS(File_Out_Local_OP_Initial_RS.c_str());
        file_out_Local_OP_Initial_RS<<"#row col OParams_[row][col]"<<endl;

        int col_,row_,col_OP, row_OP, col_new,row_new, row_int,col_int;
        complex<double> value_new;
        int alpha, alpha_p, sigma, sigma_p, gamma, gamma_p;
        int alpha_plus_gamma, alphap_plus_gammap;
        int alpha_1, alpha_2, alpha_p_1, alpha_p_2;

        for(int index_=0;index_<OPs_.value.size();index_++){

            //        row_ = gamma  + ( (d1_org*UnitCellSize_x) + alpha_1)*n_orbs_
            //   + (( (d2_org*UnitCellSize_y) + alpha_2)*lx_*n_orbs_)
            //  + sigma*(lx_*ly_*n_orbs_);

            //        //alpha_p + gamma_p*(S_) +  sigma_p*(n_orbs_*S_) + cell_*(2*n_orbs_*S_);
            //        row_OP = 0*(2*n_orbs_*UnitCellSize_x*UnitCellSize_y) + alpha + gamma*(UnitCellSize_x*UnitCellSize_y) + sigma*(n_orbs_*UnitCellSize_x*UnitCellSize_y);
            //        col_OP = cell_no*(2*n_orbs_*UnitCellSize_x*UnitCellSize_y) + alpha_p + gamma_p*(UnitCellSize_x*UnitCellSize_y) + sigma_p*(n_orbs_*UnitCellSize_x*UnitCellSize_y);
            //        cell_no =0
            row_OP = OPs_.rows[index_];
            col_OP = OPs_.columns[index_];
            alpha_plus_gamma = row_OP%(n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y);
            alpha = alpha_plus_gamma%(UnitCellSize_x*UnitCellSize_y);
            gamma = (alpha_plus_gamma -alpha)/(UnitCellSize_x*UnitCellSize_y);
            sigma = (row_OP - alpha_plus_gamma)/(n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y);

            alphap_plus_gammap = col_OP%(n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y);
            alpha_p = alphap_plus_gammap%(UnitCellSize_x*UnitCellSize_y);
            gamma_p = (alphap_plus_gammap -alpha_p)/(UnitCellSize_x*UnitCellSize_y);
            sigma_p = (col_OP - alphap_plus_gammap)/(n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y);

            alpha_1 = alpha % UnitCellSize_x; //a = a1 +a2*lx
            alpha_p_1 = alpha_p % UnitCellSize_x;

            alpha_2 = (alpha - alpha_1)/UnitCellSize_x;
            alpha_p_2 = (alpha_p - alpha_p_1)/UnitCellSize_x;


            row_ = gamma  + (alpha_1)*n_orbs_*n_atoms_
                   + ((alpha_2)*lx_*n_orbs_*n_atoms_)
                   + sigma*(lx_*ly_*n_orbs_*n_atoms_);

            col_ = gamma_p  + (alpha_p_1)*n_orbs_*n_atoms_
                   + ((alpha_p_2)*lx_*n_orbs_*n_atoms_)
                   + sigma_p*(lx_*ly_*n_orbs_*n_atoms_);

            row_int = gamma  + (alpha_1)*n_orbs_*n_atoms_
                      + ((alpha_2)*lx_*n_orbs_*n_atoms_);

            col_int = gamma_p  + (alpha_p_1)*n_orbs_*n_atoms_
                      + ((alpha_p_2)*lx_*n_orbs_*n_atoms_);

            if(col_<row_){
                col_new=row_;
                row_new=col_;
                value_new=conj(OPs_.value[index_]);
            }
            else{
                col_new=col_;
                row_new=row_;
                value_new=OPs_.value[index_];
            }
            if((abs(Connections_.Hint_(row_int, col_int))>0.00000001) || (row_int==col_int)){
                file_out_Local_OP_Initial_RS<<row_new<<setw(15)<<col_new<<setw(15)<<"  "<<value_new<<endl;
            }
        }
        //-------------------------------------------------
    }


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
        Get_Bands();

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
        cout<<E_class<<"   "<<E_quant<<"    "<<endl;


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
        //Calculate_Akw(0.0);

        string Eigenvalues_fl_out = "Eigenvalues.txt";
        ofstream Eigenvalues_file_out(Eigenvalues_fl_out.c_str());
        for(int ie=0;ie<Eigenvalues_.size();ie++){
            Eigenvalues_file_out<<ie<<"  "<<Eigenvalues_[ie]<<endl;
        }







        string File_Out_Local_OP = Parameters_.File_OPs_out + "_Temp"+string(temp_char) + ".txt";
        ofstream file_out_Local_OP(File_Out_Local_OP.c_str());
        file_out_Local_OP<<"#row col OParams_[row][col]"<<endl;
        file_out_Local_OP<<lx_<<  "  "<<ly_<<endl;

        for(int alpha=0;alpha<OPs_.value.size();alpha++){
            file_out_Local_OP<<OPs_new_.rows[alpha]<<setw(15)<<OPs_new_.columns[alpha]<<setw(15)<<"  "<<OPs_new_.value[alpha]<<endl;
        }




        if(UnitCellSize_x==lx_ && UnitCellSize_y==ly_){
            //---------For Real space--------------------------
            string File_Out_Local_OP_RS = Parameters_.File_OPs_out + "_Temp"+string(temp_char) + "_for_RealSpace.txt";
            ofstream file_out_Local_OP_RS(File_Out_Local_OP_RS.c_str());
            file_out_Local_OP_RS<<"#row col OParams_[row][col]"<<endl;

            int col_,row_,col_OP, row_OP, col_new, row_new, row_int, col_int;
            complex<double> value_new;
            int alpha, alpha_p, sigma, sigma_p, gamma, gamma_p;
            int alpha_plus_gamma, alphap_plus_gammap;
            int alpha_1, alpha_2, alpha_p_1, alpha_p_2;

            for(int index_=0;index_<OPs_.value.size();index_++){

                //        row_ = gamma  + ( (d1_org*UnitCellSize_x) + alpha_1)*n_orbs_
                //   + (( (d2_org*UnitCellSize_y) + alpha_2)*lx_*n_orbs_)
                //  + sigma*(lx_*ly_*n_orbs_);

                //        //alpha_p + gamma_p*(S_) +  sigma_p*(n_orbs_*S_) + cell_*(2*n_orbs_*S_);
                //        row_OP = 0*(2*n_orbs_*UnitCellSize_x*UnitCellSize_y) + alpha + gamma*(UnitCellSize_x*UnitCellSize_y) + sigma*(n_orbs_*UnitCellSize_x*UnitCellSize_y);
                //        col_OP = cell_no*(2*n_orbs_*UnitCellSize_x*UnitCellSize_y) + alpha_p + gamma_p*(UnitCellSize_x*UnitCellSize_y) + sigma_p*(n_orbs_*UnitCellSize_x*UnitCellSize_y);
                //        cell_no =0
                row_OP = OPs_.rows[index_];
                col_OP = OPs_.columns[index_];
                alpha_plus_gamma = row_OP%(n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y);
                alpha = alpha_plus_gamma%(UnitCellSize_x*UnitCellSize_y);
                gamma = (alpha_plus_gamma -alpha)/(UnitCellSize_x*UnitCellSize_y);
                sigma = (row_OP - alpha_plus_gamma)/(n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y);

                alphap_plus_gammap = col_OP%(n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y);
                alpha_p = alphap_plus_gammap%(UnitCellSize_x*UnitCellSize_y);
                gamma_p = (alphap_plus_gammap -alpha_p)/(UnitCellSize_x*UnitCellSize_y);
                sigma_p = (col_OP - alphap_plus_gammap)/(n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y);

                alpha_1 = alpha % UnitCellSize_x; //a = a1 +a2*lx
                alpha_p_1 = alpha_p % UnitCellSize_x;

                alpha_2 = (alpha - alpha_1)/UnitCellSize_x;
                alpha_p_2 = (alpha_p - alpha_p_1)/UnitCellSize_x;


                row_ = gamma  + (alpha_1)*n_orbs_*n_atoms_
                       + ((alpha_2)*lx_*n_orbs_*n_atoms_)
                       + sigma*(lx_*ly_*n_orbs_*n_atoms_);

                col_ = gamma_p  + (alpha_p_1)*n_orbs_*n_atoms_
                       + ((alpha_p_2)*lx_*n_orbs_*n_atoms_)
                       + sigma_p*(lx_*ly_*n_orbs_*n_atoms_);


                row_int = gamma  + (alpha_1)*n_orbs_*n_atoms_
                          + ((alpha_2)*lx_*n_orbs_*n_atoms_);

                col_int = gamma_p  + (alpha_p_1)*n_orbs_*n_atoms_
                          + ((alpha_p_2)*lx_*n_orbs_*n_atoms_);

                if(col_<row_){
                    col_new=row_;
                    row_new=col_;
                    value_new=conj(OPs_new_.value[index_]);
                }
                else{
                    col_new=col_;
                    row_new=row_;
                    value_new=OPs_new_.value[index_];
                }

                if((abs(Connections_.Hint_(row_int, col_int))>0.00000001) || (row_int==col_int)){
                    file_out_Local_OP_RS<<row_new<<setw(15)<<col_new<<setw(15)<<"  "<<value_new<<endl;
                }

            }
            //-------------------------------------------------
        }




        //Getting ready for next Temperature
        kick_while_cooling=0.0;
        Parameters_.Read_OPs=true;
        Parameters_.UnitCellType_intialOPs.first=Parameters_.UnitCellSize_x;
        Parameters_.UnitCellType_intialOPs.second=Parameters_.UnitCellSize_y;
        Parameters_.File_OPs_in = Parameters_.File_OPs_out + "_Temp"+string(temp_char) + ".txt";
        Initialize();
    }


}

void Kspace_calculation_G2dLattice::Diagonalize(char option){

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


void Kspace_calculation_G2dLattice::Update_OrderParameters_AndersonMixing(int iter){

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

    assert(NewInd_to_OldInd.size()==OPs_.value.size() + (OPs_.value.size() - 2*n_orbs_*n_atoms_*UnitCellSize_x*UnitCellSize_y));
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
    //            for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
    //                for(int gamma=0;gamma<n_orbs_;gamma++){
    //                    for(int spin=0;spin<2;spin++){
    //                        row_OP = alpha + gamma*(UnitCellSize_x*UnitCellSize_y)
    //                                + spin*(n_orbs_*UnitCellSize_x*UnitCellSize_y);
    //                        col_OP = alpha + gamma*(UnitCellSize_x*UnitCellSize_y)
    //                                + spin*(n_orbs_*UnitCellSize_x*UnitCellSize_y);
    //                        index_OP = SI_to_ind[col_OP + row_OP*(2*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
    //                        OPs_total_den_ += OPs_.value[index_OP].real();
    //                    }
    //                }
    //            }
    //        }
    //    }


    //    double ratio = Parameters_.Total_Particles/(OPs_total_den_+0.0001);

    //    for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
    //        for(int gamma=0;gamma<n_orbs_;gamma++){
    //            for(int spin=0;spin<2;spin++){
    //                row_OP = alpha + gamma*(UnitCellSize_x*UnitCellSize_y)
    //                        + spin*(n_orbs_*UnitCellSize_x*UnitCellSize_y);
    //                col_OP = alpha + gamma*(UnitCellSize_x*UnitCellSize_y)
    //                        + spin*(n_orbs_*UnitCellSize_x*UnitCellSize_y);
    //                index_OP = SI_to_ind[col_OP + row_OP*(2*n_orbs_*ncells_*UnitCellSize_x*UnitCellSize_y)];
    //                OPs_.value[index_OP] += ratio/(1.0+ratio);
    //            }
    //        }
    //    }






}

void Kspace_calculation_G2dLattice::Perform_SVD(Matrix<double> & A_, Matrix<double> & VT_, Matrix<double> & U_, vector<double> & Sigma_){


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


void Kspace_calculation_G2dLattice::Perform_SVD_complex(Matrix<complex<double>> & A_, Matrix<complex<double>> & VT_, Matrix<complex<double>> & U_, vector<double> & Sigma_){


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
