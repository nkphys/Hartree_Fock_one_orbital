#include <algorithm>
#include <functional>
#include <math.h>
#include "tensor_type.h"
#include "Parameters_TL.h"
#include "Coordinates_TL.h"
#include "Connections_TL.h"
#include "random"
#include "../../Matrix.h"
#define PI acos(-1.0)

#ifndef Kspace_calculation_TL_class
#define Kspace_calculation_TL_class

extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
                         std::complex<double> *,int *, double *, int *);

extern "C" void dgesdd_ (char *, int *, int *, double *, int *, double *, double *, int *, double *, int*,
                         double *, int *, int *, int *);


extern "C" void zgesdd_ (char *, int *, int *, std::complex<double> *, int *, double *, std::complex<double> *, int *, std::complex<double> *, int*,
                         std::complex<double> *, int *, double * , int *, int *);

class Kspace_calculation_TL
{
public:
    Kspace_calculation_TL(Parameters_TL &Parameters__, Coordinates_TL &Coordinates__, Connections_TL &Connections__, mt19937_64& Generator1__ )
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
    void Get_Bands();
    void Get_minimum_distance_direction(int l,int m,int &r1_, int &r2_);
    complex<double> h_KE(int alpha, int sigma, int alpha_p, int sigma_p, int k1, int k2);
    complex<double> IntraCell_U(int alpha, int alpha_p);
    complex<double> U_Bar(int alpha, int alpha_p);
    complex<double> V_fock(int alpha, int sigma, int alpha_p, int sigma_p, int k1, int k2);

    void Perform_SVD(Matrix<double> & A_, Matrix<double> & VT_, Matrix<double> & U_, vector<double> & Sigma_);
    void Perform_SVD_complex(Matrix<complex<double>> & A_, Matrix<complex<double>> & VT_, Matrix<complex<double>> & U_, vector<double> & Sigma_);
    void Update_OrderParameters_AndersonMixing(int iter);
    void Get_spin_resolved_local_densities();
    void Get_local_spins();
    void Calculate_Nw();
    double Lorentzian(double x, double brd);
    void Calculate_ChernNumbers();
    void Calculate_ChernNumbers_in_ntimes_brillouin_zone(int N_);
    void Create_Kspace_Spectrum_in_ntimes_brillouin_zone(int N_);
    void Create_Current_Oprs();
    void Hall_conductance();
    //::DONE



    mt19937_64 &Generator1_;
    uniform_real_distribution<double> dis1_;
    Parameters_TL &Parameters_;
    Coordinates_TL &Coordinates_; //this in cell wise representation, n_orbs=no. of atoms in unitcell
    Connections_TL &Connections_;
    int lx_, ly_, ncells_, n_orbs_, UnitCellSize_x, UnitCellSize_y, lx_cells, ly_cells;
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

    Matrix_COO_Complex OPs_, OPs_new_;
    Mat_1_int SI_to_ind;

    double mu_;
    double OP_error_;
    double E_quant, E_class;

    Mat_1_intpair OP_Dcell;
    Mat_1_intpair OP_alpha;
    Mat_1_intpair OP_sigma;

    //Declarations for Anderson Mixing
    Mat_1_doub x_km1_, x_k_, Del_x_km1;
    Mat_1_doub f_k_, f_km1_, Del_f_km1;
    Mat_1_doub xbar_k_, fbar_k_, gamma_k_, x_kp1_;
    Matrix<double> X_mat, F_mat;


    vector<Matrix<complex<double>>> J_KE_e1, J_KE_e2;
    vector<Matrix<complex<double>>> J_KE_X, J_KE_Y;


    double Global_Eps;
    bool OP_only_finite_Int;
    double Temperature_global;

};

void Kspace_calculation_TL::Hall_conductance(){

    Mat_2_doub hall_cond, hall_cond_xy;
    hall_cond.resize(2); hall_cond_xy.resize(2);
    for(int spin_=0;spin_<2;spin_++){
        hall_cond[spin_].resize(2);
        hall_cond_xy[spin_].resize(2);
    }


    double cond_11, cond_22;
    double eps_temp =0.00000001;
    Parameters_.eta=0.01;

    string fileout_="sigmaxy_vs_mu.txt";
    ofstream fileout(fileout_.c_str());


    for(int spin_1=0;spin_1<2;spin_1++){
        for(int spin_2=0;spin_2<2;spin_2++){
            hall_cond[spin_1][spin_2]=0.0;
    for(int m=0;m<lx_*ly_*2;m++){
        for(int n=0;n<lx_*ly_*2;n++){
            //1if(abs(Eigenvalues_saved[m]-Eigenvalues_saved[n])>=-eps_temp){
                hall_cond[spin_1][spin_2] += ((2.0*PI)/(lx_*ly_))*( (1.0/( exp((Eigenvalues_saved[m]-mu_)*Parameters_.beta ) + 1.0))-  (1.0/( exp((Eigenvalues_saved[n]-mu_)*Parameters_.beta ) + 1.0)) )*
                        (1.0/( (Parameters_.eta*Parameters_.eta) + ((Eigenvalues_saved[n]-Eigenvalues_saved[m])*(Eigenvalues_saved[n]-Eigenvalues_saved[m]))   ))*
                        ((1.0*J_KE_e1[spin_1](m,n))*(1.0*J_KE_e2[spin_2](n,m))).imag();
            //}
        }
    }
        }
    }

    for(int spin_1=0;spin_1<2;spin_1++){
        for(int spin_2=0;spin_2<2;spin_2++){
            hall_cond_xy[spin_1][spin_2]=0.0;
    for(int m=0;m<lx_*ly_*2;m++){
        for(int n=0;n<lx_*ly_*2;n++){
            if(abs(Eigenvalues_saved[m]-Eigenvalues_saved[n])>=-eps_temp){
                hall_cond_xy[spin_1][spin_2] += ((2.0*PI)/(lx_*ly_))*( (1.0/( exp((Eigenvalues_saved[m]-mu_)*Parameters_.beta ) + 1.0))-  (1.0/( exp((Eigenvalues_saved[n]-mu_)*Parameters_.beta ) + 1.0)) )*
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

    /*
    double mu=Eigenvalues_saved[0]-1.0;
    while(mu<Eigenvalues_saved[Eigenvalues_saved.size()-1]+1.0){
        hall_cond=0.0;
        hall_cond_xy=0.0;
        cond_11=0.0;
        cond_22=0.0;

        for(int m=0;m<lx_*ly_*2;m++){
            for(int n=0;n<lx_*ly_*2;n++){
                if(abs(Eigenvalues_saved[m]-Eigenvalues_saved[n])>=eps_temp){
                    hall_cond += ((2.0*PI)/(lx_*ly_))*( (1.0/( exp((Eigenvalues_saved[m]-mu)*Parameters_.beta ) + 1.0))-  (1.0/( exp((Eigenvalues_saved[n]-mu)*Parameters_.beta ) + 1.0)) )*
                            (1.0/( (Parameters_.eta*Parameters_.eta) + ((Eigenvalues_saved[n]-Eigenvalues_saved[m])*(Eigenvalues_saved[n]-Eigenvalues_saved[m]))   ))*
                            ((1.0*J_KE_e1(m,n))*(1.0*J_KE_e2(n,m))).imag();
                    hall_cond_xy += ((2.0*PI)/(lx_*ly_))*( (1.0/( exp((Eigenvalues_saved[m]-mu)*Parameters_.beta ) + 1.0))-  (1.0/( exp((Eigenvalues_saved[n]-mu)*Parameters_.beta ) + 1.0)) )*
                            (1.0/( (Parameters_.eta*Parameters_.eta) + ((Eigenvalues_saved[n]-Eigenvalues_saved[m])*(Eigenvalues_saved[n]-Eigenvalues_saved[m]))   ))*
                            ((1.0*J_KE_X(m,n))*(1.0*J_KE_Y(n,m))).imag();
                    cond_11 += ((2.0*PI)/(lx_*ly_))*( (1.0/( exp((Eigenvalues_saved[m]-mu)*Parameters_.beta ) + 1.0))-  (1.0/( exp((Eigenvalues_saved[n]-mu)*Parameters_.beta ) + 1.0)) )*
                            (1.0/( (Parameters_.eta*Parameters_.eta) + ((Eigenvalues_saved[n]-Eigenvalues_saved[m])*(Eigenvalues_saved[n]-Eigenvalues_saved[m]))   ))*
                            ((1.0*J_KE_e1(m,n))*(1.0*J_KE_e1(n,m))).imag();
                    cond_22 += ((2.0*PI)/(lx_*ly_))*( (1.0/( exp((Eigenvalues_saved[m]-mu)*Parameters_.beta ) + 1.0))-  (1.0/( exp((Eigenvalues_saved[n]-mu)*Parameters_.beta ) + 1.0)) )*
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


void Kspace_calculation_TL::Create_Current_Oprs(){



    J_KE_e1.resize(2);J_KE_e2.resize(2);
    J_KE_X.resize(2); J_KE_Y.resize(2);


    for(int sigma=0;sigma<2;sigma++){
        J_KE_e1[sigma].resize(lx_*ly_*2, lx_*ly_*2);
        J_KE_e2[sigma].resize(lx_*ly_*2, lx_*ly_*2);

        J_KE_X[sigma].resize(lx_*ly_*2, lx_*ly_*2);
        J_KE_Y[sigma].resize(lx_*ly_*2, lx_*ly_*2);
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


     for(int sigma=0;sigma<2;sigma++){

    for(int l=0;l<2*S_;l++){
        for(int m=0;m<2*S_;m++){

            for(int k1=0;k1<lx_cells;k1++){
                for(int k2=0;k2<ly_cells;k2++){
                    k_index = Coordinates_.Ncell(k1,k2);
                    state_p = 2*S_*k_index + l;
                    state_n = 2*S_*k_index + m;



                        for(int alpha=0;alpha<S_;alpha++){
                            for(int alpha_p=0;alpha_p<S_;alpha_p++){

                                sigma_p=sigma;
                                c1 = alpha + sigma*(S_);
                                c2 = alpha_p + sigma*(S_);

                                alpha_1 = alpha % UnitCellSize_x; //a = a1 +a2*lx
                                alpha_p_1 = alpha_p % UnitCellSize_x;
                                alpha_2 = (alpha - alpha_1)/UnitCellSize_x;
                                alpha_p_2 = (alpha_p - alpha_p_1)/UnitCellSize_x;

                                for(int cell_no=0;cell_no<ncells_;cell_no++){

                                    d1_org = Coordinates_.indx_cellwise(cell_no);
                                    d2_org = Coordinates_.indy_cellwise(cell_no);

                                    row_ = ( (d1_org*UnitCellSize_x) + alpha_1) + ((d2_org*UnitCellSize_y) + alpha_2)*lx_ + sigma*(lx_*ly_);
                                    col_ = (0 + alpha_p_1) + (0 + alpha_p_2)*lx_ + sigma_p*(lx_*ly_);

                                    Get_minimum_distance_direction(0, cell_no, d1_, d2_);

                                    d1_net = ( (d1_*UnitCellSize_x) + alpha_1) - (0 + alpha_p_1);
                                    d2_net = ((d2_*UnitCellSize_y) + alpha_2) - (0 + alpha_p_2);


                                    dx_ = 1.0*d1_net + 0.5*(d2_net);
                                    dy_ = (sqrt(3.0)/2.0)*d2_net;
                                    //                                    facx = (dx_/sqrt((dx_*dx_) + (dy_*dy_)));
                                    //                                    facy = (dy_/sqrt((dx_*dx_) + (dy_*dy_)));
                                    if(dx_==0){
                                        facx=0;
                                    }
                                    else{
                                        facx=dx_/abs(dx_);
                                    }
                                    if(dy_==0){
                                        facy=0;
                                    }
                                    else{
                                        facy=dy_/abs(dy_);
                                    }



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
                                        J_KE_X[sigma](state_p,state_n) += (1.0)*facx*iota_complex*(
                                                    (Connections_.HTB_(row_,col_)*conj(Eigvectors_[state_p][c2])*Eigvectors_[state_n][c1]*exp(iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_ )))
                                                    -(Connections_.HTB_(col_,row_)*conj(Eigvectors_[state_p][c1])*Eigvectors_[state_n][c2]*exp(-1.0*iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_ )))
                                                    );

                                        J_KE_Y[sigma](state_p,state_n) += (1.0)*facy*iota_complex*(
                                                    (Connections_.HTB_(row_,col_)*conj(Eigvectors_[state_p][c2])*Eigvectors_[state_n][c1]*exp(iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_ )))
                                                    -(Connections_.HTB_(col_,row_)*conj(Eigvectors_[state_p][c1])*Eigvectors_[state_n][c2]*exp(-1.0*iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_ )))
                                                    );

                                        J_KE_e1[sigma](state_p,state_n) += (1.0)*fac1*iota_complex*(
                                                    (Connections_.HTB_(row_,col_)*conj(Eigvectors_[state_p][c2])*Eigvectors_[state_n][c1]*exp(iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_ )))
                                                    -(Connections_.HTB_(col_,row_)*conj(Eigvectors_[state_p][c1])*Eigvectors_[state_n][c2]*exp(-1.0*iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_ )))
                                                    );

                                        J_KE_e2[sigma](state_p,state_n) += (1.0)*fac2*iota_complex*(
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




void Kspace_calculation_TL::Calculate_ChernNumbers_in_ntimes_brillouin_zone(int N_){


    cout<<"Chern numbers in "<<N_<<" times Brilloun Zone XXXXXXXXXXXXXXXX"<<endl;
    //cout<<"here"<<endl;
    Matrix<complex<double>> F_mat; //F1, F2, F3, F4, F5;
    F_mat.resize(2*UnitCellSize_x*UnitCellSize_y, N_*N_*lx_cells*ly_cells);
    //cout<<"here"<<endl;

    int no_of_bands=2*UnitCellSize_x*UnitCellSize_y;
    complex<double> Ux_k, Uy_k, Ux_kpy, Uy_kpx;
    vector<complex<double>> F_bands;
    F_bands.resize(2*UnitCellSize_x*UnitCellSize_y);
    vector<complex<double>> Chern_num;
    Chern_num.resize(2*UnitCellSize_x*UnitCellSize_y);


    for (int band = 0; band < 2*UnitCellSize_x*UnitCellSize_y; band++)
    {
        //cout<<"here"<<endl;
        string file_Fk="Fk_band"+to_string(band)+"_" + to_string(N_) +"times_brillouin_zone.txt";
        ofstream fl_Fk_out(file_Fk.c_str());
        fl_Fk_out<<"#nx  ny  tilde_F(nx,ny).real()  tilde_F(nx,ny).imag()  ArgofLog.real()  ArgofLog.imag()"<<endl;
        fl_Fk_out<<"#Extra momentum point for pm3d corners2color c1"<<endl;

        F_bands[band] = 0.0;
        for (int nx = 0; nx < N_*lx_cells; nx++)
        {
            for (int ny = 0; ny < N_*ly_cells; ny++)
            {

                int n = nx + N_*lx_cells*ny;
                int n_left, n_right, nx_left, ny_left, nx_right, ny_right;

                //U1_k
                Ux_k = 0;
                n_left = n;
                nx_right = (nx + 1) % (N_*lx_cells);
                ny_right = ny;
                n_right = nx_right + N_*lx_cells*ny_right;
                for (int comp = 0; comp < no_of_bands; comp++)
                {
                    Ux_k +=
                            conj(Eigvectors_[no_of_bands*n_left +  band][comp]) *
                            Eigvectors_[no_of_bands*n_right + band][comp];
                }
                Ux_k = Ux_k * (1.0 / abs(Ux_k));

                //U2_kpx
                Uy_kpx = 0;
                nx_left = (nx + 1) % (N_*lx_cells);
                ny_left = ny;
                n_left = nx_left + N_*lx_cells*ny_left;
                nx_right = nx_left;
                ny_right = (ny_left + 1) % (N_*ly_cells);
                n_right = nx_right + N_*lx_cells*ny_right;
                for (int comp = 0; comp < no_of_bands; comp++)
                {
                    Uy_kpx +=
                            conj(Eigvectors_[no_of_bands*n_left +  band][comp]) *
                            Eigvectors_[no_of_bands*n_right + band][comp];
                }
                Uy_kpx = Uy_kpx * (1.0 / abs(Uy_kpx));

                //U1_kpy
                Ux_kpy = 0;
                nx_left = nx;
                ny_left = (ny + 1) % (N_*ly_cells);
                n_left = nx_left + N_*lx_cells*ny_left;
                nx_right = (nx_left + 1) % (N_*lx_cells);
                ny_right = ny_left;
                n_right = nx_right + N_*lx_cells*ny_right;
                for (int comp = 0; comp < no_of_bands; comp++)
                {
                    Ux_kpy +=
                            conj(Eigvectors_[no_of_bands*n_left +  band][comp]) *
                            Eigvectors_[no_of_bands*n_right + band][comp];
                }
                Ux_kpy = Ux_kpy * (1.0 / abs(Ux_kpy));

                //U2_k
                Uy_k = 0;
                nx_left = nx;
                ny_left = ny;
                n_left = nx_left + N_*lx_cells*ny_left;
                nx_right = nx_left;
                ny_right = (ny_left + 1) % (N_*ly_cells);
                n_right = nx_right + N_*lx_cells*ny_right;
                for (int comp = 0; comp < no_of_bands; comp++)
                {
                    Uy_k +=
                            conj(Eigvectors_[no_of_bands*n_left +  band][comp]) *
                            Eigvectors_[no_of_bands*n_right + band][comp];
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


                if(ny==N_*ly_cells-1){//For pm3d corners2color c1
                    fl_Fk_out<<nx<<"  "<<ny<<"  "<<F_mat(band, n).real()<<"  "<<F_mat(band, n).imag()<<
                               "  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).real()<<
                               "  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).imag()<<endl;
                }
            }
            if(nx==N_*lx_cells-1){//For pm3d corners2color c1
                fl_Fk_out<<endl;
                for(int ny_=0;ny_<N_*ly_cells;ny_++){
                    int n_ = nx + N_*lx_cells* ny_;
                    fl_Fk_out<<nx<<"  "<<ny_<<"  "<<F_mat(band, n_).real()<<"  "<<F_mat(band, n_).imag()<<
                               "  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).real()<<
                               "  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).imag()<<endl;
                }
            }
            fl_Fk_out<<endl;
        }


        Chern_num[band] = (-1.0 * iota_complex / (2 * PI*N_*N_)) * F_bands[band];
        // Chern_num_orgnl[band] = (-1.0 * iota_complex / (2 * PI)) * F_bands_orgnl[band];
        fl_Fk_out<<"#Chern no*2pi*Iota= "<<F_bands[band].real()<<"  "<<F_bands[band].imag()<<endl;
        cout << "tilde Chern number [" << band << "] = " << Chern_num[band].real() << "        " << Chern_num[band].imag() << endl;
        //  cout << "Chern number [" << band << "] = " << Chern_num_orgnl[band].real() << " " << Chern_num_orgnl[band].imag() << endl;

    }

}


void Kspace_calculation_TL::Calculate_ChernNumbers(){



    Matrix<complex<double>> F_mat; //F1, F2, F3, F4, F5;
    F_mat.resize(2*UnitCellSize_x*UnitCellSize_y, lx_cells*ly_cells);

    complex<double> Ux_k, Uy_k, Ux_kpy, Uy_kpx;
    vector<complex<double>> F_bands;
    F_bands.resize(2*UnitCellSize_x*UnitCellSize_y);
    vector<complex<double>> Chern_num;
    Chern_num.resize(2*UnitCellSize_x*UnitCellSize_y);

    vector<complex<double>> F_bands_orgnl;
    F_bands_orgnl.resize(2*UnitCellSize_x*UnitCellSize_y);
    vector<complex<double>> Chern_num_orgnl;
    Chern_num_orgnl.resize(2*UnitCellSize_x*UnitCellSize_y);
    for (int band = 0; band < 2*UnitCellSize_x*UnitCellSize_y; band++)
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
                for (int comp = 0; comp < 2*UnitCellSize_x*UnitCellSize_y; comp++)
                {
                    Ux_k +=
                            conj(Eigvectors_[2*UnitCellSize_x*UnitCellSize_y*n_left +  band][comp]) *
                            Eigvectors_[2*UnitCellSize_x*UnitCellSize_y*n_right + band][comp];
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
                for (int comp = 0; comp < 2*UnitCellSize_x*UnitCellSize_y; comp++)
                {
                    Uy_kpx +=
                            conj(Eigvectors_[2*UnitCellSize_x*UnitCellSize_y*n_left +  band][comp]) *
                            Eigvectors_[2*UnitCellSize_x*UnitCellSize_y*n_right + band][comp];
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
                for (int comp = 0; comp < 2*UnitCellSize_x*UnitCellSize_y; comp++)
                {
                    Ux_kpy +=
                            conj(Eigvectors_[2*UnitCellSize_x*UnitCellSize_y*n_left +  band][comp]) *
                            Eigvectors_[2*UnitCellSize_x*UnitCellSize_y*n_right + band][comp];
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
                for (int comp = 0; comp < 2*UnitCellSize_x*UnitCellSize_y; comp++)
                {
                    Uy_k +=
                            conj(Eigvectors_[2*UnitCellSize_x*UnitCellSize_y*n_left +  band][comp]) *
                            Eigvectors_[2*UnitCellSize_x*UnitCellSize_y*n_right + band][comp];
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

double Kspace_calculation_TL::Lorentzian(double x, double brd)
{
    double temp;

    temp = (1.0 / PI) * ((brd / 2.0) / ((x * x) + ((brd * brd) / 4.0)));

    return temp;
}

void Kspace_calculation_TL::Calculate_Nw()
{

    char temp_char[50];
    sprintf(temp_char, "%.10f", Parameters_.Temperature);

    int ind_max;
    ind_max=Eigenvalues_.size();

    //---------Read from input file-----------------------//
    string fileout = "Nw" + string(temp_char)+ ".txt";
    double omega_min, omega_max, d_omega;
    double eta = 0.2;
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


void Kspace_calculation_TL::Get_minimum_distance_direction(int l,int m,int &r1_, int &r2_){

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

void Kspace_calculation_TL::Get_Bands(){

    char temp_char[50];
    sprintf(temp_char, "%.10f", Parameters_.Temperature);

    string File_Out_Bands;
    File_Out_Bands = "Bands" + string(temp_char)+ ".txt";
    ofstream file_out_bands(File_Out_Bands.c_str());

    string File_Out_Bands2;
    File_Out_Bands2 = "Bands_Path2_" + string(temp_char) +".txt";
    ofstream file_out_bands2(File_Out_Bands2.c_str());

    string File_Out_Bands3;
    File_Out_Bands3 = "Bands_All_k_" + string(temp_char) +".txt";
    ofstream file_out_bands3(File_Out_Bands3.c_str());

    int k_index;
    int kx_i, ky_i;

    Mat_1_intpair k_path;
    k_path.clear();
    pair_int temp_pair;


    // ---k_path---------

    //--------\Gamma to X-----------------
    ky_i = 0;
    for (kx_i = 0; kx_i <= (lx_cells/ 2); kx_i++)
    {
        temp_pair.first = kx_i;
        temp_pair.second = ky_i;
        k_path.push_back(temp_pair);
    }
    //----------------------------------

    //--------X to M-----------------
    kx_i = (lx_cells / 2);
    for (ky_i = 1; ky_i <= (ly_cells / 2); ky_i++)
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
        temp_pair.first = kx_i;
        temp_pair.second = kx_i;
        k_path.push_back(temp_pair);
    }

    temp_pair.first = 0;
    temp_pair.second = 0;
    k_path.push_back(temp_pair);



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


    for (int k_point = 0; k_point < k_path.size(); k_point++)
    {
        k_index=Coordinates_.Ncell(k_path[k_point].first, k_path[k_point].second);

        file_out_bands<<k_point<<"   ";
        for(int band=0;band<2*UnitCellSize_x*UnitCellSize_y;band++){
            file_out_bands<<Eigenvalues_[2*UnitCellSize_x*UnitCellSize_y*k_index + band]<<"   ";
        }
        file_out_bands<<endl;
    }


    for (int k_point = 0; k_point < k_path2.size(); k_point++)
    {
        k_index=Coordinates_.Ncell(k_path2[k_point].first, k_path2[k_point].second);

        file_out_bands2<<k_point<<"   ";
        for(int band=0;band<2*UnitCellSize_x*UnitCellSize_y;band++){
            file_out_bands2<<Eigenvalues_[2*UnitCellSize_x*UnitCellSize_y*k_index + band]<<"   ";
        }
        file_out_bands2<<endl;
    }



    for(int k1_=0;k1_<=lx_cells;k1_++){
        for(int k2_=0;k2_<=ly_cells;k2_++){

            int k1_temp = k1_%lx_cells;
            int k2_temp = k2_%ly_cells;
            k_index=Coordinates_.Ncell(k1_temp, k2_temp);

            file_out_bands3<<k_index<<"  "<<k1_<<"  "<<k2_<<"  ";
            for(int band=0;band<2*UnitCellSize_x*UnitCellSize_y;band++){
                file_out_bands3<<Eigenvalues_[2*UnitCellSize_x*UnitCellSize_y*k_index + band]<<"   ";
            }
            file_out_bands3<<endl;

        }
        file_out_bands3<<endl;

    }

}

void Kspace_calculation_TL::Get_Energies(){


    complex<double> E_class_temp=0.0;
    E_class= 0.0;


    int row_, col_, index;
    int row_OP, col_OP, index_OP;

    int UP_, DOWN_;
    UP_=0;DOWN_=1;
    int spin_bar;
    int k_index;
    double k1x_val, k2x_val, k1y_val, k2y_val;
    double fac_;
    complex<double> OP_value, OP_value2;
    for(int k1=0;k1<lx_cells;k1++){
        k1x_val = (2.0*PI*k1)/(1.0*lx_cells);
        k1y_val = ((2.0*PI*k1)/(1.0*lx_cells))*(-1.0/sqrt(3));

        for(int k2=0;k2<ly_cells;k2++){
            k2x_val = 0.0;
            k2y_val = ((2.0*PI*k2)/(1.0*ly_cells))*(2.0/sqrt(3));

            k_index = Coordinates_.Ncell(k1,k2);
            Kx_values[k_index]=k1x_val + k2x_val;
            Ky_values[k_index]=k1y_val + k2y_val;


            //Anisotropy
            for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                for(int spin1=0;spin1<2;spin1++){
                    for(int spin2=0;spin2<2;spin2++){

                        row_ = alpha + spin1*(UnitCellSize_x*UnitCellSize_y);
                        col_=row_;
                        index = SI_to_ind[col_ + row_*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];

                        row_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin2*(UnitCellSize_x*UnitCellSize_y);
                        col_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin2*(UnitCellSize_x*UnitCellSize_y);
                        index_OP = SI_to_ind[col_OP + row_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];

                        fac_ = 1.0 - (2.0*abs(spin1-spin2)); //i.e 0 --> 1, 1 --> -1
                        E_class_temp += (1.0/4.0)*fac_*Parameters_.AnisotropyZ*OPs_.value[index_OP]*OPs_.value[index];

                    }
                }
            }


            //Interaction:Hartree
            //On-site Hartree
            for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                for(int spin=0;spin<2;spin++){
                    spin_bar = (spin +1 )%2;
                    row_ = alpha + spin*(UnitCellSize_x*UnitCellSize_y);
                    col_=row_;
                    index = SI_to_ind[col_ + row_*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];

                    row_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin_bar*(UnitCellSize_x*UnitCellSize_y);
                    col_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin_bar*(UnitCellSize_x*UnitCellSize_y);
                    index_OP = SI_to_ind[col_OP + row_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                    E_class_temp += -0.5*Parameters_.U0*OPs_.value[index_OP]*OPs_.value[index];
                }
            }

            //Intra-Cell Hartree
            for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                for(int spin=0;spin<2;spin++){
                    for(int alpha_p=0;alpha_p<UnitCellSize_x*UnitCellSize_y;alpha_p++){
                        if(alpha != alpha_p){
                            for(int spin_p=0;spin_p<2;spin_p++){

                                row_ = alpha_p + spin_p*(UnitCellSize_x*UnitCellSize_y);
                                col_=row_;
                                index = SI_to_ind[col_ + row_*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];

                                row_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin*(UnitCellSize_x*UnitCellSize_y);
                                col_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin*(UnitCellSize_x*UnitCellSize_y);
                                index_OP = SI_to_ind[col_OP + row_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                E_class_temp += -0.5*IntraCell_U(alpha, alpha_p)*OPs_.value[index_OP]*OPs_.value[index];
                            }
                        }
                    }
                }
            }

            //Inter-Cell Hartree
            for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                for(int spin=0;spin<2;spin++){
                    for(int alpha_p=0;alpha_p<UnitCellSize_x*UnitCellSize_y;alpha_p++){
                        for(int spin_p=0;spin_p<2;spin_p++){

                            row_ = alpha_p + spin_p*(UnitCellSize_x*UnitCellSize_y);
                            col_=row_;
                            index = SI_to_ind[col_ + row_*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];

                            row_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin*(UnitCellSize_x*UnitCellSize_y);
                            col_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin*(UnitCellSize_x*UnitCellSize_y);
                            index_OP = SI_to_ind[col_OP + row_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            E_class_temp += -0.5*U_Bar(alpha, alpha_p)*OPs_.value[index_OP]*OPs_.value[index];
                        }
                    }
                }
            }



            //Interaction:Fock
            if(!Parameters_.Just_Hartree){

                //Onsite-Fock
                if( (Parameters_.FockType=="Onsite" || Parameters_.FockType=="Onsite_Intra") || Parameters_.FockType=="Onsite_Intra_Inter" ){
                    for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                        for(int spin=0;spin<2;spin++){
                            spin_bar = (spin +1 )%2;
                            row_ = alpha + spin*(UnitCellSize_x*UnitCellSize_y);
                            col_= alpha + spin_bar*(UnitCellSize_x*UnitCellSize_y);
                            if(col_>=row_){
                                index = SI_to_ind[col_ + row_*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_value2=OPs_.value[index];
                            }
                            else{
                                index = SI_to_ind[row_ + col_*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_value2= conj(OPs_.value[index]);
                            }

                            row_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin_bar*(UnitCellSize_x*UnitCellSize_y);
                            col_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin*(UnitCellSize_x*UnitCellSize_y);

                            if(col_OP>=row_OP){
                                index_OP = SI_to_ind[col_OP + row_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_value=OPs_.value[index_OP];
                            }
                            else{
                                index_OP = SI_to_ind[row_OP + col_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_value= conj(OPs_.value[index_OP]);
                            }

                            E_class_temp += -1.0*Parameters_.U0*OP_value*(-0.5*OP_value2);
                        }
                    }
                }

                //Intra-Fock
                if( Parameters_.FockType=="Onsite_Intra" || Parameters_.FockType=="Onsite_Intra_Inter" ){

                    int alpha_1, alpha_2, alpha_p_1, alpha_p_2, d1_org, d2_org;
                    int row_U,col_U;
                    //----
                    for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                        for(int spin=0;spin<2;spin++){
                            for(int alpha_p=0;alpha_p<UnitCellSize_x*UnitCellSize_y;alpha_p++){
                                if(alpha != alpha_p){
                                    for(int spin_p=0;spin_p<2;spin_p++){

                                        alpha_1 = alpha % UnitCellSize_x;
                                        alpha_p_1 = alpha_p % UnitCellSize_x;
                                        alpha_2 = (alpha - alpha_1)/UnitCellSize_x;
                                        alpha_p_2 = (alpha_p - alpha_p_1)/UnitCellSize_x;
                                        d1_org = Coordinates_.indx_cellwise(0);
                                        d2_org = Coordinates_.indy_cellwise(0);
                                        row_U = (0 + alpha_1) + (0 + alpha_2)*lx_;
                                        col_U = ((d1_org*UnitCellSize_x) + alpha_p_1) + ((d2_org*UnitCellSize_y) + alpha_p_2)*lx_;


                                        row_ = alpha_p + spin_p*(UnitCellSize_x*UnitCellSize_y);
                                        col_= alpha + spin*(UnitCellSize_x*UnitCellSize_y);
                                        if(col_>=row_){
                                            index = SI_to_ind[col_ + row_*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                            OP_value2=OPs_.value[index];
                                        }
                                        else{
                                            index = SI_to_ind[row_ + col_*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                            OP_value2= conj(OPs_.value[index]);
                                        }


                                        col_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha_p + spin_p*(UnitCellSize_x*UnitCellSize_y);
                                        row_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin*(UnitCellSize_x*UnitCellSize_y);

                                        if(col_OP>=row_OP){
                                            index_OP = SI_to_ind[col_OP + row_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                            OP_value=OPs_.value[index_OP];
                                        }
                                        else{
                                            index_OP = SI_to_ind[row_OP + col_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                            OP_value= conj(OPs_.value[index_OP]);
                                        }
                                        if( (!OP_only_finite_Int) ||  (abs(Connections_.Hint_(row_U, col_U)) > Global_Eps )){
                                            E_class_temp += -1.0*IntraCell_U(alpha, alpha_p)*OP_value*(-0.5*OP_value2);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    //-----
                }

                //Inter UnitCell-Fock
                if(Parameters_.FockType=="Onsite_Intra_Inter" ){

                    //                    for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                    //                        for(int spin=0;spin<2;spin++){
                    //                            for(int alpha_p=0;alpha_p<UnitCellSize_x*UnitCellSize_y;alpha_p++){
                    //                                for(int spin_p=0;spin_p<2;spin_p++){

                    //                                    row_ = alpha_p + spin_p*(UnitCellSize_x*UnitCellSize_y);
                    //                                    col_= alpha + spin*(UnitCellSize_x*UnitCellSize_y);
                    //                                    if(col_>=row_){
                    //                                        index = SI_to_ind[col_ + row_*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                    //                                        OP_value2=OPs_.value[index];
                    //                                    }
                    //                                    else{
                    //                                        index = SI_to_ind[row_ + col_*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                    //                                        OP_value2= conj(OPs_.value[index]);
                    //                                    }

                    //                                    E_class_temp += -1.0*V_fock(alpha,spin,alpha_p,spin_p,k1,k2)*(-0.5*OP_value2);
                    //                                }
                    //                            }
                    //                        }
                    //                    }


                    //....................

                    for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                        for(int spin=0;spin<2;spin++){
                            for(int alpha_p=0;alpha_p<UnitCellSize_x*UnitCellSize_y;alpha_p++){
                                for(int spin_p=0;spin_p<2;spin_p++){
                                    for(int cell_no=1;cell_no<lx_cells*ly_cells;cell_no++){
                                        int alpha_1, alpha_2, alpha_p_1, alpha_p_2;
                                        int d1_org = Coordinates_.indx_cellwise(cell_no);
                                        int d2_org = Coordinates_.indy_cellwise(cell_no);

                                        alpha_1 = alpha % UnitCellSize_x;
                                        alpha_p_1 = alpha_p % UnitCellSize_x;
                                        alpha_2 = (alpha - alpha_1)/UnitCellSize_x;
                                        alpha_p_2 = (alpha_p - alpha_p_1)/UnitCellSize_x;
                                        row_ = (0 + alpha_1)+ (0 + alpha_2)*lx_;
                                        col_ = ((d1_org*UnitCellSize_x) + alpha_p_1) + ((d2_org*UnitCellSize_y) + alpha_p_2)*(lx_);


                                        col_OP = cell_no*(2*UnitCellSize_x*UnitCellSize_y) + alpha_p + spin_p*(UnitCellSize_x*UnitCellSize_y);
                                        row_OP= alpha + spin*(UnitCellSize_x*UnitCellSize_y);
                                        if((alpha_p + spin_p*(UnitCellSize_x*UnitCellSize_y))  >= (alpha + spin*(UnitCellSize_x*UnitCellSize_y))){
                                            index_OP = SI_to_ind[col_OP + row_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                            OP_value2=OPs_.value[index_OP];
                                        }
                                        else{
                                            int d1_new = (lx_cells - d1_org)%lx_cells;
                                            int d2_new = (ly_cells - d2_org)%ly_cells;
                                            int cell_new = Coordinates_.Ncell(d1_new, d2_new);
                                            int row_new = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha_p + spin_p*(UnitCellSize_x*UnitCellSize_y);
                                            int col_new = cell_new*(2*UnitCellSize_x*UnitCellSize_y) + alpha +spin*(UnitCellSize_x*UnitCellSize_y);
                                            index_OP = SI_to_ind[col_new + row_new*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                            OP_value2= conj(OPs_.value[index_OP]);
                                        }

                                        if( (!OP_only_finite_Int) ||  (abs(Connections_.Hint_(row_, col_)) > Global_Eps )){

                                            E_class_temp +=0.5*Connections_.Hint_(row_,col_)*OP_value2*conj(OP_value2);
                                        }

                                    }}
                            }
                        }
                    }


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

void Kspace_calculation_TL::Get_new_OPs_and_error(){

    //For <c_{c1}* c_{c2}>
    int c1;
    int c2;
    int S_=UnitCellSize_x*UnitCellSize_y;
    int cell_row, alpha_row, sigma_row, cell_col, alpha_col, sigma_col;
    int alpha_plus_sigma;
    int row_temp, col_temp;

    int k_index;
    int d1_, d2_;
    int state;

    for(int OP_no=0;OP_no<OPs_.value.size();OP_no++){
        row_temp=OPs_.rows[OP_no];
        col_temp=OPs_.columns[OP_no];

        cell_row = 0;
        alpha_row = row_temp%S_;
        sigma_row = (row_temp - alpha_row)/S_;
        c1 = alpha_row + sigma_row*(S_);

        alpha_plus_sigma = col_temp%(2*S_);
        alpha_col = alpha_plus_sigma%S_;
        sigma_col = (alpha_plus_sigma - alpha_col)/S_;
        cell_col = (col_temp - alpha_plus_sigma)/(2*S_);
        Get_minimum_distance_direction(0, cell_col, d1_, d2_);
        c2 = alpha_col + sigma_col*(S_);

        OPs_new_.value[OP_no]=0.0;
        for(int n=0;n<2*S_;n++){ //band_index
            for(int k1=0;k1<lx_cells;k1++){
                for(int k2=0;k2<ly_cells;k2++){
                    k_index = Coordinates_.Ncell(k1,k2);
                    state = 2*S_*k_index + n;

                    OPs_new_.value[OP_no] += (1.0/ncells_)*(
                                (conj(Eigvectors_[state][c1])*Eigvectors_[state][c2])*(exp(-1.0*iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_   )  ))
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

void Kspace_calculation_TL::Get_spin_resolved_local_densities(){

    //For <c_{c1}* c_{c2}>
    int c1;
    int c2;
    int S_=UnitCellSize_x*UnitCellSize_y;
    int cell_row, alpha_row, sigma_row, cell_col, alpha_col, sigma_col;
    int alpha_plus_sigma;
    int row_temp, col_temp;

    int k_index;
    int d1_, d2_;
    int state;

    double Total_den_up, Total_den_dn;

    char temp_char[50];
    sprintf(temp_char, "%.10f", Parameters_.Temperature);

    string File_Out_Local_orb_densities = "Local_spin_resolved_densities" + string(temp_char)+ ".txt";
    ofstream file_out_Local_orb_densities(File_Out_Local_orb_densities.c_str());
    file_out_Local_orb_densities<<"#alpha(site in real lattice)   up   dn"<<endl;

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
                    for(int sigma=0;sigma<2;sigma++){

                        c1 = alpha + sigma*(S_);
                        c2=c1;
                        val=0.0;
                        for(int n=0;n<2*S_;n++){ //band_index
                            for(int k1=0;k1<lx_cells;k1++){
                                for(int k2=0;k2<ly_cells;k2++){
                                    k_index = Coordinates_.Ncell(k1,k2);
                                    state = 2*S_*k_index + n;

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
                    file_out_Local_orb_densities<<endl;
                }
            }
        }
    }


    cout<< "Total electrons in Lattice, UP = "<<Total_den_up<<endl;
    cout<< "Total electrons in Lattice, DOWN = "<<Total_den_dn<<endl;


}

void Kspace_calculation_TL::Get_local_spins(){

    //For <c_{c1}* c_{c2}>
    int c1;
    int c2;
    int S_=UnitCellSize_x*UnitCellSize_y;
    int cell_row, alpha_row, sigma_row, cell_col, alpha_col, sigma_col;
    int alpha_plus_sigma;
    int row_temp, col_temp;

    int k_index;
    int d1_, d2_;
    int state;

    double Total_Sz, Total_Sx, Total_Sy;
    int UP_, DOWN_;
    UP_=0;DOWN_=1;

    char temp_char[50];
    sprintf(temp_char, "%.10f", Parameters_.Temperature);

    string File_Out_Local_orb_densities = "Local_spins"  + string(temp_char) +  ".txt";
    ofstream file_out_Local_orb_densities(File_Out_Local_orb_densities.c_str());
    file_out_Local_orb_densities<<"#site (in real lattice) site_1 site_2  rx     ry   Sz   Sx  Sy"<<endl;

    complex<double> val;

    Total_Sz=0.0;
    Total_Sx=0.0;
    Total_Sy=0.0;
    double rx_, ry_;
    complex<double> splus_val;
    double sz_val, sx_val, sy_val;
    int site_x, site_y, alpha, site;
    for(int alpha_1=0;alpha_1<UnitCellSize_x;alpha_1++){
        for(int cell_1=0;cell_1<lx_cells;cell_1++){
            site_x = (cell_1*UnitCellSize_x) + alpha_1;

            for(int alpha_2=0;alpha_2<UnitCellSize_y;alpha_2++){
                for(int cell_2=0;cell_2<ly_cells;cell_2++){
                    site_y = (cell_2*UnitCellSize_y) + alpha_2;

                    rx_ = ((1.0)*(site_x) +  (1.0/2.0)*(site_y));
                    ry_ =  (0.0*(site_x) + (sqrt(3.0)/2.0)*(site_y));

                    alpha = alpha_1 + alpha_2*(UnitCellSize_x);
                    site = site_x + site_y*(lx_);

                    file_out_Local_orb_densities<<site<<"    "<<site_x<<"    "<<site_y<<"    "<<rx_<<"     "<<ry_<<"     ";


                    //Local Splus
                    val=0.0;
                    c1 = alpha + UP_*(S_);
                    c2 = alpha + DOWN_*(S_);;
                    for(int n=0;n<2*S_;n++){ //band_index
                        for(int k1=0;k1<lx_cells;k1++){
                            for(int k2=0;k2<ly_cells;k2++){
                                k_index = Coordinates_.Ncell(k1,k2);
                                state = 2*S_*k_index + n;
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

                    //Local Sz
                    val=0.0;
                    double fac_;
                    for(int sigma=0;sigma<2;sigma++){
                        fac_=1.0 - (2.0*sigma);
                        c1 = alpha + sigma*(S_);
                        c2=c1;
                        for(int n=0;n<2*S_;n++){ //band_index
                            for(int k1=0;k1<lx_cells;k1++){
                                for(int k2=0;k2<ly_cells;k2++){
                                    k_index = Coordinates_.Ncell(k1,k2);
                                    state = 2*S_*k_index + n;

                                    val += 0.5*(fac_/ncells_)*(
                                                (conj(Eigvectors_[state][c1])*Eigvectors_[state][c2])
                                                *(1.0/( exp((Eigenvalues_saved[state]-mu_)*Parameters_.beta ) + 1.0))
                                                );
                                }
                            }
                        }
                    }
                    sz_val=val.real();


                    file_out_Local_orb_densities<<sz_val<<"    "<<sx_val<<"   "<<sy_val<<"    ";


                    Total_Sz +=sz_val;
                    Total_Sy +=sy_val;
                    Total_Sx +=sx_val;



                    file_out_Local_orb_densities<<endl;
                }
            }
        }
    }


    cout<< "Total Sz = "<<Total_Sz<<endl;
    cout<< "Total Sx = "<<Total_Sx<<endl;
    cout<< "Total Sy = "<<Total_Sy<<endl;



}

double Kspace_calculation_TL::chemicalpotential(double Particles){


    double mu_out;
    double n1,N;
    double dMubydN;
    double muin;
    N=Particles;
    int N_ = int(N);
    muin = 0.5*(Eigenvalues_[N_-1] - Eigenvalues_[N_]);
    int nstate = Eigenvalues_.size();
    dMubydN = 0.0005*(Eigenvalues_[nstate-1] - Eigenvalues_[0])/nstate;

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
                //cout<<i<<"    "<<n1<<"    "<<N-n1<<endl;

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
        mu1=Eigenvalues_[0]- (5.0/Parameters_.beta);;
        mu2=Eigenvalues_[nstate-1]+ (5.0/Parameters_.beta);;
        for(int i=0;i<40000;i++){
            n1=0.0;
            for(int j=0;j<nstate;j++){
                n1+=double(1.0/( exp( (Eigenvalues_[j]-mu_temp)*Parameters_.beta ) + 1.0));
            }
            //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
            if(abs(N-n1)<double(0.00001)){
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
            cout<<"mu_not_converged, N = "<<n1<<endl;
        }
        else{
            //cout<<"mu converged, N = "<<n1<<endl;
        }

        mu_out = mu_temp;
    }



    return mu_out;
} // ----------

double Kspace_calculation_TL::random1(){

    return dis1_(Generator1_);

}

void Kspace_calculation_TL::Initialize()
{

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


    Ham_.resize(2*UnitCellSize_x*UnitCellSize_y, 2*UnitCellSize_x*UnitCellSize_y);
    Eigvectors_.resize(2*lx_*ly_);
    Eigenvalues_.resize(2*lx_*ly_);
    for(int i=0;i<2*lx_*ly_;i++){
        Eigvectors_[i].resize(2*UnitCellSize_x*UnitCellSize_y); //Eigenvector number
    }

    Kx_values.resize(ncells_);
    Ky_values.resize(ncells_);


    string  str_Initial_OrderParams_file = "Initial_OPs_ToUse_In_RealSpace.txt";
    ofstream Initial_OrderParams_file(str_Initial_OrderParams_file.c_str());
    Initial_OrderParams_file<<"#alpha_i     alpha_j        OP"<<endl;

    //Order Parameters :
    /*
    Inter-unit cell OP's ( intra-unit cell when d_vec=0)
    <c_{alpha, sigma}^{dagger} c _{ d_vec, alpha', sigma' }> ; for all d_vec, and alpha+sigma*(UnitCellSize_x*UnitCellSize_y) <= alpha'+sigma'*(UnitCellSize_x*UnitCellSize_y)

    For <c_{alpha, sigma}^{dagger} c _{ d_vec, alpha', sigma' }> , when alpha+sigma*(UnitCellSize_x*UnitCellSize_y) > alpha'+sigma'*(UnitCellSize_x*UnitCellSize_y)
    use
    <c_{alpha, sigma}^{dagger} c _{ d_vec, alpha', sigma' }>  =  conj(<c_{alpha', sigma'}^{dagger} c _{-d_vec, alpha, sigma}>)


    S=2*UnitCellSize_x*UnitCellSize_y

    Total no of Order parameters (Hartree + all Fock):
    ncells_*(S(S+1))/2

    Total no of Order parameters (Hartree + only intracell Fock):
    (S(S+1))/2

    Total no of Order parameters (Hartree):
    S

    */

    int S_=UnitCellSize_x*UnitCellSize_y;
    int r1, r2;
    OPs_.value.clear();
    OPs_.rows.clear();
    OPs_.columns.clear();
    OPs_new_.value.clear();
    OPs_new_.rows.clear();
    OPs_new_.columns.clear();
    SI_to_ind.resize(2*S_ *  ncells_*2*S_);

    int row_temp, col_temp; // col_ = (cell_*2*S_) + alpha + sigma*(S_);




    if(!Parameters_.Read_OPs){


        if(!Parameters_.Using_Initial_Ansatz){
            //Hartree Terms
            for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                for(int sigma=0;sigma<2;sigma++){
                    OPs_.value.push_back(complex<double> (random1(),0.0));
                    //OPs_.value.push_back(complex<double> (0.2,0.0));

                    OPs_new_.value.push_back(0.0);

                    row_temp= 0*(2*S_) + alpha + sigma*S_;
                    col_temp=row_temp;

                    OPs_.rows.push_back(row_temp);OPs_new_.rows.push_back(row_temp);
                    OPs_.columns.push_back(col_temp);OPs_new_.columns.push_back(col_temp);

                    SI_to_ind[col_temp + row_temp*(ncells_*2*S_)] = OPs_.value.size()-1;

                }
            }


            //Fock Terms
            if(!Parameters_.Just_Hartree){

                bool check_;
                int alpha_1, alpha_2, alpha_p_1, alpha_p_2, d1_org, d2_org;
                int row_,col_;
                for(int cell_=0;cell_<ncells_;cell_++){
                    for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                        for(int sigma=0;sigma<2;sigma++){
                            for(int alpha_p=0;alpha_p<UnitCellSize_x*UnitCellSize_y;alpha_p++){
                                for(int sigma_p=0;sigma_p<2;sigma_p++){

                                    alpha_1 = alpha % UnitCellSize_x;
                                    alpha_p_1 = alpha_p % UnitCellSize_x;
                                    alpha_2 = (alpha - alpha_1)/UnitCellSize_x;
                                    alpha_p_2 = (alpha_p - alpha_p_1)/UnitCellSize_x;

                                    d1_org = Coordinates_.indx_cellwise(cell_);
                                    d2_org = Coordinates_.indy_cellwise(cell_);

                                    row_ = (0 + alpha_1) + (0 + alpha_2)*lx_;
                                    col_ = ((d1_org*UnitCellSize_x) + alpha_p_1) + ((d2_org*UnitCellSize_y) + alpha_p_2)*lx_;

                                    if(Parameters_.FockType=="Onsite_Intra_Inter"){
                                        if(cell_==0){
                                            check_= (alpha_p + sigma_p*(UnitCellSize_x*UnitCellSize_y)) > (alpha + sigma*(UnitCellSize_x*UnitCellSize_y));
                                        }
                                        else{
                                            check_ =(alpha_p + sigma_p*(UnitCellSize_x*UnitCellSize_y)) >= (alpha + sigma*(UnitCellSize_x*UnitCellSize_y));
                                        }
                                    }


                                    if(Parameters_.FockType=="Onsite_Intra"){
                                        if(cell_==0){
                                            check_ = (alpha_p + sigma_p*(UnitCellSize_x*UnitCellSize_y)) > (alpha + sigma*(UnitCellSize_x*UnitCellSize_y));
                                        }
                                        else{ //cell_!=0
                                            check_=false;
                                        }
                                    }

                                    if(Parameters_.FockType=="Onsite"){
                                        if(cell_==0){
                                            if(alpha_p==alpha){
                                                check_ = (alpha_p + sigma_p*(UnitCellSize_x*UnitCellSize_y)) > (alpha + sigma*(UnitCellSize_x*UnitCellSize_y));
                                            }
                                            else{ //alpha !=alpha_p
                                                check_=false;
                                            }
                                        }
                                        else{ //cell_!=0
                                            check_=false;
                                        }
                                    }

                                    if(row_!=col_){
                                        if( (OP_only_finite_Int) && (abs(Connections_.Hint_(row_, col_)) < Global_Eps )){
                                            check_=false;
                                        }
                                    }

                                    if( check_ ){

                                        row_temp = 0*(2*S_) + alpha + sigma*S_;
                                        col_temp = cell_*(2*S_) + alpha_p + sigma_p*S_;

                                        OPs_.value.push_back(complex<double> (random1(),random1()));
                                        OPs_new_.value.push_back(0.0);

                                        OPs_.rows.push_back(row_temp);OPs_new_.rows.push_back(row_temp);
                                        OPs_.columns.push_back(col_temp);OPs_new_.columns.push_back(col_temp);

                                        SI_to_ind[col_temp + row_temp*(ncells_*2*S_)] = OPs_.value.size()-1;


                                    }

                                }
                            }
                        }
                    }
                }

            }


        }

        else{

            if(Parameters_.Initial_Ansatz_type=="Tetrahedron"){
            cout<<"Using Tetrahedron Ansatz"<<endl;

                /*
             1: (1.0/sqrt(3))*(1,1,1)=(Sx,Sy,Sz)
             2: (1.0/sqrt(3))*(1,-1,-1)
             3: (1.0/sqrt(3))*(-1,1,-1)
             4: (1.0/sqrt(3))*(-1,-1,1)
             <n>=0.5

             <n_up>=0.5*<n> + Sz
             <n_up>=0.5*<n> - Sz
             <c_up^dag c_dn> = Sx + i*Sy
            */

                double den_,Sz_,Sx_,Sy_;
                int alpha_1, alpha_2;
                double value_temp;
                // int Site_type;
                den_=0.5;
                //            Mat_2_doub den_array, Sx_array, Sy_array, Sz_array;
                //            den_array.push_back();

                //Hartree Terms
                for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){

                    alpha_1 = alpha % UnitCellSize_x; //a = a1 +a2*lx
                    alpha_2 = (alpha - alpha_1)/UnitCellSize_x;

                    if( (alpha_1%2==0) && (alpha_2%2==0)){
                        Sx_=(1.0/sqrt(3.0));Sy_=(1.0/sqrt(3.0));Sz_=(1.0/sqrt(3.0));
                    }
                    if( (alpha_1%2==1) && (alpha_2%2==0)){
                        Sx_=(1.0/sqrt(3.0));Sy_=(-1.0/sqrt(3.0));Sz_=(-1.0/sqrt(3.0));
                    }
                    if( (alpha_1%2==0) && (alpha_2%2==1)){
                        Sx_=(-1.0/sqrt(3.0));Sy_=(1.0/sqrt(3.0));Sz_=(-1.0/sqrt(3.0));
                    }
                    if( (alpha_1%2==1) && (alpha_2%2==1)){
                        Sx_=(-1.0/sqrt(3.0));Sy_=(-1.0/sqrt(3.0));Sz_=(1.0/sqrt(3.0));
                    }


                    for(int sigma=0;sigma<2;sigma++){

                        value_temp = 0.5*den_ + (1.0 - (2.0*sigma))*Sz_;
                        OPs_.value.push_back(complex<double> (value_temp,0.0));
                        OPs_new_.value.push_back(0.0);

                        row_temp= 0*(2*S_) + alpha + sigma*S_;
                        col_temp=row_temp;

                        OPs_.rows.push_back(row_temp);OPs_new_.rows.push_back(row_temp);
                        OPs_.columns.push_back(col_temp);OPs_new_.columns.push_back(col_temp);

                        SI_to_ind[col_temp + row_temp*(ncells_*2*S_)] = OPs_.value.size()-1;

                    }

                }


                //Fock Terms
                if(!Parameters_.Just_Hartree){

                    bool check_;
                    for(int cell_=0;cell_<ncells_;cell_++){
                        for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){

                            alpha_1 = alpha % UnitCellSize_x; //a = a1 +a2*lx
                            alpha_2 = (alpha - alpha_1)/UnitCellSize_x;

                            if( (alpha_1%2==0) && (alpha_2%2==0)){
                                Sx_=(1.0/sqrt(3.0));Sy_=(1.0/sqrt(3.0));Sz_=(1.0/sqrt(3.0));
                            }
                            if( (alpha_1%2==1) && (alpha_2%2==0)){
                                Sx_=(1.0/sqrt(3.0));Sy_=(-1.0/sqrt(3.0));Sz_=(-1.0/sqrt(3.0));
                            }
                            if( (alpha_1%2==0) && (alpha_2%2==1)){
                                Sx_=(-1.0/sqrt(3.0));Sy_=(1.0/sqrt(3.0));Sz_=(-1.0/sqrt(3.0));
                            }
                            if( (alpha_1%2==1) && (alpha_2%2==1)){
                                Sx_=(-1.0/sqrt(3.0));Sy_=(-1.0/sqrt(3.0));Sz_=(1.0/sqrt(3.0));
                            }


                            for(int sigma=0;sigma<2;sigma++){
                                for(int alpha_p=0;alpha_p<UnitCellSize_x*UnitCellSize_y;alpha_p++){
                                    for(int sigma_p=0;sigma_p<2;sigma_p++){

                                        if(Parameters_.FockType=="Onsite_Intra_Inter"){
                                            if(cell_==0){
                                                check_= (alpha_p + sigma_p*(UnitCellSize_x*UnitCellSize_y)) > (alpha + sigma*(UnitCellSize_x*UnitCellSize_y));
                                            }
                                            else{
                                                check_ =(alpha_p + sigma_p*(UnitCellSize_x*UnitCellSize_y)) >= (alpha + sigma*(UnitCellSize_x*UnitCellSize_y));
                                            }
                                        }


                                        if(Parameters_.FockType=="Onsite_Intra"){
                                            if(cell_==0){
                                                check_ = (alpha_p + sigma_p*(UnitCellSize_x*UnitCellSize_y)) > (alpha + sigma*(UnitCellSize_x*UnitCellSize_y));
                                            }
                                            else{ //cell_!=0
                                                check_=false;
                                            }
                                        }

                                        if(Parameters_.FockType=="Onsite"){
                                            if(cell_==0){
                                                if(alpha_p==alpha){
                                                    check_ = (alpha_p + sigma_p*(UnitCellSize_x*UnitCellSize_y)) > (alpha + sigma*(UnitCellSize_x*UnitCellSize_y));
                                                }
                                                else{ //alpha !=alpha_p
                                                    check_=false;
                                                }
                                            }
                                            else{ //cell_!=0
                                                check_=false;
                                            }
                                        }

                                        if( check_ ){

                                            row_temp = 0*(2*S_) + alpha + sigma*S_;
                                            col_temp = cell_*(2*S_) + alpha_p + sigma_p*S_;

                                            if(!(Parameters_.FockType=="Onsite")){
                                                OPs_.value.push_back(complex<double> (random1(),random1()));
                                            }
                                            else{
                                                OPs_.value.push_back(complex<double> (Sx_,Sy_));
                                            }

                                            OPs_new_.value.push_back(0.0);

                                            OPs_.rows.push_back(row_temp);OPs_new_.rows.push_back(row_temp);
                                            OPs_.columns.push_back(col_temp);OPs_new_.columns.push_back(col_temp);

                                            SI_to_ind[col_temp + row_temp*(ncells_*2*S_)] = OPs_.value.size()-1;


                                        }

                                    }
                                }
                            }
                        }
                    }

                }


            }


            //-----------
            if(Parameters_.Initial_Ansatz_type=="Tri_Prism"){

                // Using Trigonal Prism  with height=d and eq. triangle edge=a

                /* n=1, s_max=n/2=0.5 for all 6 spins
                   d/a is free parameter.
                   "a" is fixed that |S_{i}|=s_max


                1, A: (0,-a/sqrt(3), d/2) = (Sx,Sy,Sz)
                2, B: (a/2,a/(2*sqrt(3)), d/2)
                3, C: (-a/2, a/(2*sqrt(3)), d/2)
                4, D: (0,-a/sqrt(3), -d/2)
                5, E: (a/2,a/(2*sqrt(3)), -d/2)
                6, F: (-a/2, a/(2*sqrt(3)), -d/2)
                */

                Mat_1_doub Sz_vals, Sx_vals, Sy_vals;
                Sz_vals.resize(18);Sx_vals.resize(18);Sy_vals.resize(18);

                double den_;
                double s_max=0.5;
                double d_by_a=0.5;
                double a_;
                double d_;
                double Sz_, Sx_, Sy_;
                int alpha_1, alpha_2;
                int alpha_eff;
                double value_temp;

                a_= s_max/(sqrt( (1.0/3.0)  +  (d_by_a*d_by_a*0.25)   ));
                d_ = a_*d_by_a;

                for(int i_=0;i_<18;i_++){
                    if(i_==1 || i_==6 || i_==14 ){ //A
                        Sx_vals[i_]=0.0; Sy_vals[i_]=-1.0*a_*(1.0/sqrt(3.0));Sz_vals[i_]=0.5*d_;
                    }
                    if(i_==4 || i_==9 || i_==17 ){ //B
                        Sx_vals[i_]=a_*0.5; Sy_vals[i_]=0.5*a_*(1.0/sqrt(3.0)); Sz_vals[i_]=0.5*d_;
                    }
                    if(i_==3 || i_==11 || i_==16 ){ //C
                        Sx_vals[i_]=-0.5*a_; Sy_vals[i_]=0.5*a_*(1.0/sqrt(3.0));Sz_vals[i_]=0.5*d_;
                    }

                    if(i_==5 || i_==10 || i_==15 ){ //D
                        Sx_vals[i_]=0.0; Sy_vals[i_]=-1.0*a_*(1.0/sqrt(3.0));Sz_vals[i_]=-0.5*d_;
                    }
                    if(i_==0 || i_==8 || i_==13 ){ //E
                        Sx_vals[i_]=a_*0.5; Sy_vals[i_]=0.5*a_*(1.0/sqrt(3.0)); Sz_vals[i_]=-0.5*d_;
                    }
                    if(i_==2 || i_==7 || i_==12 ){ //F
                        Sx_vals[i_]=-0.5*a_; Sy_vals[i_]=0.5*a_*(1.0/sqrt(3.0));Sz_vals[i_]=-0.5*d_;
                    }

                }


                //Hartree Terms
                for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){

                    alpha_1 = alpha % UnitCellSize_x; //a = a1 +a2*lx
                    alpha_2 = (alpha - alpha_1)/UnitCellSize_x;

                    alpha_eff = (alpha_1%3) +  3*(alpha_2%6);

                    for(int sigma=0;sigma<2;sigma++){
                        for(int gamma=0;gamma<n_orbs_;gamma++){
                            if(gamma==0){
                                den_=1.0;
                                Sz_=Sz_vals[alpha_eff];
                                value_temp = 0.5*den_ + (1.0 - (2.0*sigma))*Sz_;
                                OPs_.value.push_back(complex<double> (value_temp,0.0));
                            }
                            else{
                                den_=0.0;
                                OPs_.value.push_back(complex<double> (0.0,0.0));
                            }


                            OPs_new_.value.push_back(0.0);

                            row_temp=  alpha + gamma*(S_) +  sigma*(n_orbs_*S_) + 0*(2*n_orbs_*S_);
                            col_temp=row_temp;

                            OPs_.rows.push_back(row_temp);OPs_new_.rows.push_back(row_temp);
                            OPs_.columns.push_back(col_temp);OPs_new_.columns.push_back(col_temp);
                            SI_to_ind[col_temp + row_temp*(ncells_*2*n_orbs_*S_)] = OPs_.value.size()-1;

                        }
                    }
                }


                //Fock Terms
                if(!Parameters_.Just_Hartree){
                    bool check_;

                    int alpha_p_1, alpha_p_2, d1_org, d2_org;
                    int row_,col_;
                    for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                        alpha_1 = alpha % UnitCellSize_x; //a = a1 +a2*lx
                        alpha_2 = (alpha - alpha_1)/UnitCellSize_x;

                        alpha_eff = (alpha_1%3) +  3*(alpha_2%6);

                        for(int gamma=0;gamma<n_orbs_;gamma++){
                            if(gamma==1){
                                Sx_=0.0;Sy_=0.0;
                            }
                            else{
                                Sx_=Sx_vals[alpha_eff];
                                Sy_=Sy_vals[alpha_eff];
                            }

                            for(int sigma=0;sigma<2;sigma++){

                                for(int cell_=0;cell_<ncells_;cell_++){
                                    for(int alpha_p=0;alpha_p<UnitCellSize_x*UnitCellSize_y;alpha_p++){
                                        for(int gamma_p=0;gamma_p<n_orbs_;gamma_p++){



                                            alpha_p_1 = alpha_p % UnitCellSize_x;
                                            alpha_p_2 = (alpha_p - alpha_p_1)/UnitCellSize_x;
                                            d1_org = Coordinates_.indx_cellwise(cell_);
                                            d2_org = Coordinates_.indy_cellwise(cell_);
                                            row_ = gamma + (0 + alpha_1)*n_orbs_ + (0 + alpha_2)*(n_orbs_*lx_);
                                            col_ = gamma_p + ((d1_org*UnitCellSize_x) + alpha_p_1)*n_orbs_ + ((d2_org*UnitCellSize_y) + alpha_p_2)*(n_orbs_*lx_);




                                            for(int sigma_p=0;sigma_p<2;sigma_p++){

                                                if(Parameters_.FockType=="Onsite_Intra_Inter"){
                                                    if(cell_==0){
                                                        check_= ((alpha_p + gamma_p*(S_) +  sigma_p*(n_orbs_*S_)) > (alpha + gamma*(S_) +  sigma*(n_orbs_*S_)));
                                                    }
                                                    else{
                                                        check_ = ((alpha_p + gamma_p*(S_) +  sigma_p*(n_orbs_*S_)) >= (alpha + gamma*(S_) +  sigma*(n_orbs_*S_)));
                                                    }
                                                }


                                                if(Parameters_.FockType=="Onsite_Intra"){
                                                    if(cell_==0){
                                                        check_ = ((alpha_p + gamma_p*(S_) +  sigma_p*(n_orbs_*S_)) > (alpha + gamma*(S_) +  sigma*(n_orbs_*S_)));
                                                    }
                                                    else{ //cell_!=0
                                                        check_=false;
                                                    }
                                                }

                                                if(Parameters_.FockType=="Onsite"){
                                                    if(cell_==0){
                                                        if(alpha_p==alpha){
                                                            check_ = ((alpha_p + gamma_p*(S_) +  sigma_p*(n_orbs_*S_)) > (alpha + gamma*(S_) +  sigma*(n_orbs_*S_)));
                                                        }
                                                        else{ //alpha !=alpha_p
                                                            check_=false;
                                                        }
                                                    }
                                                    else{ //cell_!=0
                                                        check_=false;
                                                    }
                                                }

                                                if(row_!=col_){
                                                    if( (OP_only_finite_Int) && (abs(Connections_.Hint_(row_, col_)) < Global_Eps )){
                                                        check_=false;
                                                    }
                                                }

                                                if( check_ ){



                                                    //--------------------
                                                    row_temp = alpha + gamma*(S_) +  sigma*(n_orbs_*S_) + 0*(2*n_orbs_*S_);
                                                    col_temp = alpha_p + gamma_p*(S_) +  sigma_p*(n_orbs_*S_) + cell_*(2*n_orbs_*S_);

                                                    if(row_!=col_){
                                                        OPs_.value.push_back(complex<double> (0.05*random1(),0.05*random1()));
                                                    }
                                                    else{
                                                        OPs_.value.push_back(complex<double> (1*Sx_,1*Sy_));
                                                    }

                                                    OPs_new_.value.push_back(0.0);

                                                    OPs_.rows.push_back(row_temp);OPs_new_.rows.push_back(row_temp);
                                                    OPs_.columns.push_back(col_temp);OPs_new_.columns.push_back(col_temp);

                                                    SI_to_ind[col_temp + row_temp*(ncells_*2*n_orbs_*S_)] = OPs_.value.size()-1;
                                                    //-------------------


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


            //----------




        }
    }

    else{
        //Hartree Terms
        for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
            for(int sigma=0;sigma<2;sigma++){
                OPs_.value.push_back(complex<double> (0.0,0.0));
                OPs_new_.value.push_back(0.0);
                row_temp= 0*(2*S_) + alpha + sigma*S_;
                col_temp=row_temp;
                OPs_.rows.push_back(row_temp);OPs_new_.rows.push_back(row_temp);
                OPs_.columns.push_back(col_temp);OPs_new_.columns.push_back(col_temp);
                SI_to_ind[col_temp + row_temp*(ncells_*2*S_)] = OPs_.value.size()-1;
            }
        }
        //Fock Terms
        if(!Parameters_.Just_Hartree){
            bool check_;
            for(int cell_=0;cell_<ncells_;cell_++){
                for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                    for(int sigma=0;sigma<2;sigma++){
                        for(int alpha_p=0;alpha_p<UnitCellSize_x*UnitCellSize_y;alpha_p++){
                            for(int sigma_p=0;sigma_p<2;sigma_p++){

                                if(Parameters_.FockType=="Onsite_Intra_Inter"){
                                    if(cell_==0){
                                        check_= (alpha_p + sigma_p*(UnitCellSize_x*UnitCellSize_y)) > (alpha + sigma*(UnitCellSize_x*UnitCellSize_y));
                                    }
                                    else{
                                        check_ =(alpha_p + sigma_p*(UnitCellSize_x*UnitCellSize_y)) >= (alpha + sigma*(UnitCellSize_x*UnitCellSize_y));
                                    }
                                }


                                if(Parameters_.FockType=="Onsite_Intra"){
                                    if(cell_==0){
                                        check_ = (alpha_p + sigma_p*(UnitCellSize_x*UnitCellSize_y)) > (alpha + sigma*(UnitCellSize_x*UnitCellSize_y));
                                    }
                                    else{ //cell_!=0
                                        check_=false;
                                    }
                                }

                                if(Parameters_.FockType=="Onsite"){
                                    if(cell_==0){
                                        if(alpha_p==alpha){
                                            check_ = (alpha_p + sigma_p*(UnitCellSize_x*UnitCellSize_y)) > (alpha + sigma*(UnitCellSize_x*UnitCellSize_y));
                                        }
                                        else{ //alpha !=alpha_p
                                            check_=false;
                                        }
                                    }
                                    else{ //cell_!=0
                                        check_=false;
                                    }
                                }

                                if( check_ ){

                                    row_temp = 0*(2*S_) + alpha + sigma*S_;
                                    col_temp = cell_*(2*S_) + alpha_p + sigma_p*S_;

                                    //cout<<"Here"<<endl;
                                    OPs_.value.push_back(complex<double> (0.01,0.0));
                                    OPs_new_.value.push_back(0.0);

                                    OPs_.rows.push_back(row_temp);OPs_new_.rows.push_back(row_temp);
                                    OPs_.columns.push_back(col_temp);OPs_new_.columns.push_back(col_temp);

                                    SI_to_ind[col_temp + row_temp*(ncells_*2*S_)] = OPs_.value.size()-1;


                                }

                            }
                        }
                    }
                }
            }
        }



        string fl_initial_OP_in = Parameters_.File_OPs_in;
        ifstream file_initial_OP_in(fl_initial_OP_in.c_str());
        string temp1, line_temp;
        int row_, col_;
        int index_;
        complex<double> val_OP;

        //#row col <c^{dagger}_{site_i,spin_i} c_{site_j, spin_j}>
        getline(file_initial_OP_in,temp1);


        if(file_initial_OP_in.is_open())
        {
            while(getline(file_initial_OP_in,line_temp))
            {
                stringstream line_temp_ss(line_temp);
                line_temp_ss >> row_ >> col_ >> val_OP;

                if(Parameters_.Just_Hartree){
                    if(row_==col_){
                        index_ =  SI_to_ind[col_ + row_*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                        OPs_.value[index_]=val_OP;
                        assert(OPs_.rows[index_]==row_);
                        assert(OPs_.columns[index_]==col_);

                    }
                }
                else{
                    index_ =  SI_to_ind[col_ + row_*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                    OPs_.value[index_]=val_OP;
                    assert(OPs_.rows[index_]==row_);
                    assert(OPs_.columns[index_]==col_);
                }

            }
            file_initial_OP_in.close();
        }
        else
        {cout<<"Unable to open file = '"<<fl_initial_OP_in<<"'"<<endl;}

    }



    cout<<"Total no. of OP's = "<<OPs_.value.size()<<endl;


    for(int i=0;i<OPs_.value.size();i++){
        Initial_OrderParams_file<<OPs_.rows[i]<<"  "<<OPs_.columns[i]<<"  "<<OPs_.value[i]<<endl;
    }


    //    int alpha_1, alpha_p_1 ,alpha_2 ,alpha_p_2;
    //    int d1_org, d2_org;
    //    int row_, col_;
    //    int row_OP, col_OP;
    //    int index_OP;
    //    for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
    //        for(int sigma=0;sigma<2;sigma++){
    //            for(int alpha_p=0;alpha_p<UnitCellSize_x*UnitCellSize_y;alpha_p++){
    //                for(int sigma_p=0;sigma_p<2;sigma_p++){
    //                    alpha_1 = alpha % UnitCellSize_x;
    //                    alpha_p_1 = alpha_p % UnitCellSize_x;
    //                    alpha_2 = (alpha - alpha_1)/UnitCellSize_x;
    //                    alpha_p_2 = (alpha_p - alpha_p_1)/UnitCellSize_x;
    //                    for(int cell_no=0;cell_no<ncells_;cell_no++){
    //                        d1_org = Coordinates_.indx_cellwise(cell_no);
    //                        d2_org = Coordinates_.indy_cellwise(cell_no);
    //                        row_ = (0 + alpha_1) + (0 + alpha_2)*lx_;
    //                        col_ = ((d1_org*UnitCellSize_x) + alpha_p_1) + ((d2_org*UnitCellSize_y) + alpha_p_2)*lx_;

    //                        row_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + sigma*(UnitCellSize_x*UnitCellSize_y);
    //                        col_OP = cell_no*(2*UnitCellSize_x*UnitCellSize_y) + alpha_p + sigma_p*(UnitCellSize_x*UnitCellSize_y);

    //                        index_OP = SI_to_ind[col_OP + row_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];

    //                    Initial_OrderParams_file<<row_ + sigma*(UnitCellSize_x*UnitCellSize_y)<<"  "<<
    //                                              col_ + sigma_p*(UnitCellSize_x*UnitCellSize_y)<<"  "<<
    //                                              OPs_.value[index_OP]<<endl;
    //                    }
    //                }
    //            }
    //        }
    //    }

}


void Kspace_calculation_TL::Arranging_spectrum(){

    // Eigvectors_saved=Eigvectors_;
    Eigenvalues_saved = Eigenvalues_;

    //    Eigvectors_.resize(ncells_*6);
    //    Eigenvalues_.resize(ncells_*6);
    //    for(int i=0;i<ncells_*6;i++){
    //        Eigvectors_[i].resize(6); //Eigenvector number
    //    }

    double value_;
    int S_=2*UnitCellSize_x*UnitCellSize_y;
    Mat_1_Complex_doub Vec_temp;
    Vec_temp.resize(S_);

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

complex<double> Kspace_calculation_TL::h_KE(int alpha, int sigma, int alpha_p, int sigma_p, int k1, int k2){


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

        row_ = ( (d1_org*UnitCellSize_x) + alpha_1) + ((d2_org*UnitCellSize_y) + alpha_2)*lx_ + sigma*(lx_*ly_);
        col_ = (0 + alpha_p_1) + (0 + alpha_p_2)*lx_ + sigma_p*(lx_*ly_);

        Get_minimum_distance_direction(0, cell_no, d1_, d2_);

        temp_val += Connections_.HTB_(row_,col_)*exp(iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_   ));

    }


    return temp_val;
}


complex<double> Kspace_calculation_TL::IntraCell_U(int alpha, int alpha_p){

    complex<double> temp_val;

    int row_, col_;
    int alpha_1, alpha_2, alpha_p_1, alpha_p_2;
    alpha_1 = alpha % UnitCellSize_x; //a = a1 +a2*lx
    alpha_p_1 = alpha_p % UnitCellSize_x;

    alpha_2 = (alpha - alpha_1)/UnitCellSize_x;
    alpha_p_2 = (alpha_p - alpha_p_1)/UnitCellSize_x;

    temp_val=0.0;

    row_ = (0 + alpha_1) + (0 + alpha_2)*lx_;
    col_ = (0 + alpha_p_1) + (0 + alpha_p_2)*lx_;

    assert(alpha_1<=lx_);assert(alpha_2<=ly_);
    assert(alpha_p_1<=lx_);assert(alpha_p_2<=ly_);

    if(row_!=col_){
        temp_val = Connections_.Hint_(row_, col_);
    }
    else{
        temp_val = Parameters_.U0;
    }

    return temp_val;
}


complex<double> Kspace_calculation_TL::U_Bar(int alpha, int alpha_p){

    complex<double> temp_val;
    int d1_org, d2_org;
    int row_, col_;
    int alpha_1, alpha_2, alpha_p_1, alpha_p_2;
    alpha_1 = alpha % UnitCellSize_x;
    alpha_p_1 = alpha_p % UnitCellSize_x;

    alpha_2 = (alpha - alpha_1)/UnitCellSize_x;
    alpha_p_2 = (alpha_p - alpha_p_1)/UnitCellSize_x;

    temp_val=0.0;


    for(int cell_no=1;cell_no<ncells_;cell_no++){
        d1_org = Coordinates_.indx_cellwise(cell_no);
        d2_org = Coordinates_.indy_cellwise(cell_no);

        row_ = (0 + alpha_1) + (0 + alpha_2)*lx_;
        col_ = ((d1_org*UnitCellSize_x) + alpha_p_1) + ((d2_org*UnitCellSize_y) + alpha_p_2)*lx_;

        temp_val += Connections_.Hint_(row_, col_);

    }

    return temp_val;
}


complex<double> Kspace_calculation_TL::V_fock(int alpha, int sigma, int alpha_p, int sigma_p, int k1, int k2){

    complex<double> temp_val;
    int d1_org, d2_org, d1_new, d2_new;
    int d1_, d2_;
    int row_, col_;
    int row_OP, col_OP, index_OP;
    int row_new, col_new, cell_new;
    int alpha_1, alpha_2, alpha_p_1, alpha_p_2;
    complex<double> OP_value;

    alpha_1 = alpha % UnitCellSize_x;
    alpha_p_1 = alpha_p % UnitCellSize_x;

    alpha_2 = (alpha - alpha_1)/UnitCellSize_x;
    alpha_p_2 = (alpha_p - alpha_p_1)/UnitCellSize_x;

    temp_val=0.0;

    for(int cell_no=1;cell_no<ncells_;cell_no++){
        d1_org = Coordinates_.indx_cellwise(cell_no);
        d2_org = Coordinates_.indy_cellwise(cell_no);
        Get_minimum_distance_direction(0, cell_no, d1_, d2_);

        row_ = (0 + alpha_1) + (0 + alpha_2)*lx_;
        col_ = ((d1_org*UnitCellSize_x) + alpha_p_1) + ((d2_org*UnitCellSize_y) + alpha_p_2)*lx_;

        row_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + sigma*(UnitCellSize_x*UnitCellSize_y);
        col_OP = cell_no*(2*UnitCellSize_x*UnitCellSize_y) + alpha_p + sigma_p*(UnitCellSize_x*UnitCellSize_y);

        if( (alpha_p + sigma_p*(UnitCellSize_x*UnitCellSize_y)) >= (alpha + sigma*(UnitCellSize_x*UnitCellSize_y))){
            index_OP = SI_to_ind[col_OP + row_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
            OP_value=OPs_.value[index_OP];
        }
        else{
            d1_new = (lx_cells - d1_org)%lx_cells;
            d2_new = (ly_cells - d2_org)%ly_cells;
            cell_new = Coordinates_.Ncell(d1_new, d2_new);
            row_new = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha_p + sigma_p*(UnitCellSize_x*UnitCellSize_y);
            col_new = cell_new*(2*UnitCellSize_x*UnitCellSize_y) + alpha + sigma*(UnitCellSize_x*UnitCellSize_y);
            index_OP = SI_to_ind[col_new + row_new*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
            OP_value= conj(OPs_.value[index_OP]);
        }

        index_OP = SI_to_ind[col_OP + row_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];

        if( (!OP_only_finite_Int) ||  (abs(Connections_.Hint_(row_, col_)) > Global_Eps )){
            temp_val += Connections_.Hint_(row_, col_)*OPs_.value[index_OP]*
                    exp(iota_complex* ( ((2.0*PI*k1)/(1.0*lx_cells))*d1_   +   ((2.0*PI*k2)/(1.0*ly_cells))*d2_   ));
        }

    }

    return temp_val;
}

void Kspace_calculation_TL::Create_Kspace_Spectrum(){

    int row_, col_;
    int row_OP, col_OP, index_OP;
    int UP_, DOWN_;
    UP_=0;DOWN_=1;
    int spin_bar;
    int k_index;
    double k1x_val, k2x_val, k1y_val, k2y_val;
    double fac_;
    complex<double> OP_value;
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


            //Hoppings
            for(int alpha_p=0;alpha_p<UnitCellSize_x*UnitCellSize_y;alpha_p++){
                for(int sigma_p=0;sigma_p<2;sigma_p++){
                    row_ = alpha_p + sigma_p*(UnitCellSize_x*UnitCellSize_y);

                    for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                        for(int sigma=0;sigma<2;sigma++){
                            col_ = alpha + sigma*(UnitCellSize_x*UnitCellSize_y);

                            Ham_(row_, col_) += h_KE(alpha_p, sigma_p, alpha, sigma, k1, k2);
                        }
                    }
                }
            }


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


            //Onsite Energy: Tetrahedron, temporary
//                        for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
//                            int alpha_1 = alpha % UnitCellSize_x; //a = a1 +a2*lx
//                            int alpha_2 = (alpha - alpha_1)/UnitCellSize_x;
//                            double val_temp;
//                            if( (alpha_1%2==0) && (alpha_2%2==0)){
//                                val_temp=0.05;
//                            }
//                            if( (alpha_1%2==1) && (alpha_2%2==0)){
//                                val_temp=0.02;
//                            }
//                            if( (alpha_1%2==0) && (alpha_2%2==1)){
//                                val_temp=0.03;
//                            }
//                            if( (alpha_1%2==1) && (alpha_2%2==1)){
//                                val_temp=0.04;
//                            }

//                            for(int spin1=0;spin1<2;spin1++){
//                                row_ = alpha + spin1*(UnitCellSize_x*UnitCellSize_y);
//                                col_=row_;

//                                row_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin1*(UnitCellSize_x*UnitCellSize_y);
//                                col_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin1*(UnitCellSize_x*UnitCellSize_y);
//                                index_OP = SI_to_ind[col_OP + row_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];

//                                Ham_(row_,col_) += val_temp;

//                            }
//                        }




            //Magnetic field: Tetrahedron, temporary
//                        for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){

//                            for(int spin1=0;spin1<2;spin1++){
//                                    row_ = alpha + spin1*(UnitCellSize_x*UnitCellSize_y);
//                                    col_=row_;

//                                    Ham_(row_,col_) += 0.01*(1.0 - (2.0*spin1));

//                            }
//                        }




            //Interaction:Hartree
            //On-site Hartree
            for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                for(int spin=0;spin<2;spin++){
                    spin_bar = (spin +1 )%2;
                    row_ = alpha + spin*(UnitCellSize_x*UnitCellSize_y);
                    col_=row_;

                    row_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin_bar*(UnitCellSize_x*UnitCellSize_y);
                    col_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin_bar*(UnitCellSize_x*UnitCellSize_y);
                    index_OP = SI_to_ind[col_OP + row_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                    Ham_(row_,col_) += Parameters_.U0*OPs_.value[index_OP];
                }
            }

            //Intra-Cell Hartree
            for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                for(int spin=0;spin<2;spin++){
                    for(int alpha_p=0;alpha_p<UnitCellSize_x*UnitCellSize_y;alpha_p++){
                        if(alpha != alpha_p){
                            for(int spin_p=0;spin_p<2;spin_p++){

                                row_ = alpha_p + spin_p*(UnitCellSize_x*UnitCellSize_y);
                                col_=row_;

                                row_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin*(UnitCellSize_x*UnitCellSize_y);
                                col_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin*(UnitCellSize_x*UnitCellSize_y);
                                index_OP = SI_to_ind[col_OP + row_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                Ham_(row_,col_) += IntraCell_U(alpha, alpha_p)*OPs_.value[index_OP];
                            }
                        }
                    }
                }
            }

            //Inter-Cell Hartree
            for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                for(int spin=0;spin<2;spin++){
                    for(int alpha_p=0;alpha_p<UnitCellSize_x*UnitCellSize_y;alpha_p++){
                        for(int spin_p=0;spin_p<2;spin_p++){

                            row_ = alpha_p + spin_p*(UnitCellSize_x*UnitCellSize_y);
                            col_=row_;

                            row_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin*(UnitCellSize_x*UnitCellSize_y);
                            col_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin*(UnitCellSize_x*UnitCellSize_y);
                            index_OP = SI_to_ind[col_OP + row_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            Ham_(row_,col_) += U_Bar(alpha, alpha_p)*OPs_.value[index_OP];
                        }
                    }
                }
            }



            //Interaction:Fock
            if(!Parameters_.Just_Hartree){

                //Onsite-Fock
                if( (Parameters_.FockType=="Onsite" || Parameters_.FockType=="Onsite_Intra") || Parameters_.FockType=="Onsite_Intra_Inter" ){
                    for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                        for(int spin=0;spin<2;spin++){
                            spin_bar = (spin +1 )%2;
                            row_ = alpha + spin*(UnitCellSize_x*UnitCellSize_y);
                            col_= alpha + spin_bar*(UnitCellSize_x*UnitCellSize_y);

                            row_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin_bar*(UnitCellSize_x*UnitCellSize_y);
                            col_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin*(UnitCellSize_x*UnitCellSize_y);
                            if(col_OP>=row_OP){
                                index_OP = SI_to_ind[col_OP + row_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_value=OPs_.value[index_OP];
                            }
                            else{
                                index_OP = SI_to_ind[row_OP + col_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_value= conj(OPs_.value[index_OP]);
                            }

                            Ham_(row_,col_) += -1.0*Parameters_.U0*OP_value;

                        }
                    }
                }

                //Intra-Fock
                if( Parameters_.FockType=="Onsite_Intra" || Parameters_.FockType=="Onsite_Intra_Inter" ){

                    int alpha_1, alpha_2, alpha_p_1, alpha_p_2, d1_org, d2_org;
                    int row_U,col_U;

                    for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                        for(int spin=0;spin<2;spin++){
                            for(int alpha_p=0;alpha_p<UnitCellSize_x*UnitCellSize_y;alpha_p++){
                                if(alpha != alpha_p){
                                    for(int spin_p=0;spin_p<2;spin_p++){


                                        alpha_1 = alpha % UnitCellSize_x;
                                        alpha_p_1 = alpha_p % UnitCellSize_x;
                                        alpha_2 = (alpha - alpha_1)/UnitCellSize_x;
                                        alpha_p_2 = (alpha_p - alpha_p_1)/UnitCellSize_x;
                                        d1_org = Coordinates_.indx_cellwise(0);
                                        d2_org = Coordinates_.indy_cellwise(0);
                                        row_U = (0 + alpha_1) + (0 + alpha_2)*(lx_);
                                        col_U = ((d1_org*UnitCellSize_x) + alpha_p_1) + ((d2_org*UnitCellSize_y) + alpha_p_2)*(lx_);


                                        row_ = alpha_p + spin_p*(UnitCellSize_x*UnitCellSize_y);
                                        col_= alpha + spin*(UnitCellSize_x*UnitCellSize_y);

                                        col_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha_p + spin_p*(UnitCellSize_x*UnitCellSize_y);
                                        row_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin*(UnitCellSize_x*UnitCellSize_y);

                                        if(col_OP>=row_OP){
                                            index_OP = SI_to_ind[col_OP + row_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                            OP_value=OPs_.value[index_OP];
                                        }
                                        else{
                                            index_OP = SI_to_ind[row_OP + col_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                            OP_value= conj(OPs_.value[index_OP]);
                                        }

                                        if( (!OP_only_finite_Int) ||  (abs(Connections_.Hint_(row_U, col_U)) > Global_Eps )){
                                            Ham_(row_,col_) += -1.0*IntraCell_U(alpha, alpha_p)*OP_value;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                if(Parameters_.FockType=="Onsite_Intra_Inter" ){

                    for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                        for(int spin=0;spin<2;spin++){
                            for(int alpha_p=0;alpha_p<UnitCellSize_x*UnitCellSize_y;alpha_p++){
                                for(int spin_p=0;spin_p<2;spin_p++){

                                    row_ = alpha_p + spin_p*(UnitCellSize_x*UnitCellSize_y);
                                    col_= alpha + spin*(UnitCellSize_x*UnitCellSize_y);

                                    Ham_(row_,col_) += -1.0*V_fock(alpha,spin,alpha_p,spin_p,k1,k2);
                                }
                            }
                        }
                    }
                }




            }


            //cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
            //Ham_.print();
            //cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;

            char Dflag='V';
            Diagonalize(Dflag);


            for(int row=0;row<2*UnitCellSize_x*UnitCellSize_y;row++){
                Eigenvalues_[2*UnitCellSize_x*UnitCellSize_y*k_index + row]=eigs_[row];
                for(int col=0;col<2*UnitCellSize_x*UnitCellSize_y;col++){
                    Eigvectors_[2*UnitCellSize_x*UnitCellSize_y*k_index + col][row]=Ham_(row,col);
                }
            }
        }
    }

}




void Kspace_calculation_TL::Create_Kspace_Spectrum_in_ntimes_brillouin_zone(int N_){


    int no_of_bands = 2*UnitCellSize_x*UnitCellSize_y;;
    Eigvectors_.resize(ncells_*no_of_bands*N_*N_);
    Eigenvalues_.resize(ncells_*no_of_bands*N_*N_);
    for(int i=0;i<ncells_*no_of_bands*N_*N_;i++){
        Eigvectors_[i].resize(no_of_bands); //Eigenvector number
    }

    Kx_values.resize(ncells_*N_*N_);
    Ky_values.resize(ncells_*N_*N_);



    int row_, col_;
    int row_OP, col_OP, index_OP;
    int UP_, DOWN_;
    UP_=0;DOWN_=1;
    int spin_bar;
    int k_index;
    double k1x_val, k2x_val, k1y_val, k2y_val;
    double fac_;
    complex<double> OP_value;


    for(int k1=0;k1<N_*lx_cells;k1++){
        k1x_val = (2.0*PI*k1)/(1.0*lx_cells);
        k1y_val = ((2.0*PI*k1)/(1.0*lx_cells))*(-1.0/sqrt(3));

        for(int k2=0;k2<N_*ly_cells;k2++){
            k2x_val = 0.0;
            k2y_val = ((2.0*PI*k2)/(1.0*ly_cells))*(2.0/sqrt(3));

            k_index = k1 + N_*lx_cells*k2; //Coordinates_.Ncell(k1,k2);
            Kx_values[k_index]=k1x_val + k2x_val;
            Ky_values[k_index]=k1y_val + k2y_val;


            Ham_.fill(0.0);


            //Hoppings
            for(int alpha_p=0;alpha_p<UnitCellSize_x*UnitCellSize_y;alpha_p++){
                for(int sigma_p=0;sigma_p<2;sigma_p++){
                    row_ = alpha_p + sigma_p*(UnitCellSize_x*UnitCellSize_y);

                    for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                        for(int sigma=0;sigma<2;sigma++){
                            col_ = alpha + sigma*(UnitCellSize_x*UnitCellSize_y);

                            Ham_(row_, col_) += h_KE(alpha_p, sigma_p, alpha, sigma, k1, k2);
                        }
                    }
                }
            }


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


            //Onsite Energy: Tetrahedron, temporary
            //            for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
            //                int alpha_1 = alpha % UnitCellSize_x; //a = a1 +a2*lx
            //                int alpha_2 = (alpha - alpha_1)/UnitCellSize_x;
            //                double val_temp;
            //                if( (alpha_1%2==0) && (alpha_2%2==0)){
            //                    val_temp=0.01;
            //                }
            //                if( (alpha_1%2==1) && (alpha_2%2==0)){
            //                    val_temp=0.02;
            //                }
            //                if( (alpha_1%2==0) && (alpha_2%2==1)){
            //                    val_temp=0.03;
            //                }
            //                if( (alpha_1%2==1) && (alpha_2%2==1)){
            //                    val_temp=0.04;
            //                }

            //                for(int spin1=0;spin1<2;spin1++){
            //                    row_ = alpha + spin1*(UnitCellSize_x*UnitCellSize_y);
            //                    col_=row_;

            //                    row_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin1*(UnitCellSize_x*UnitCellSize_y);
            //                    col_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin1*(UnitCellSize_x*UnitCellSize_y);
            //                    index_OP = SI_to_ind[col_OP + row_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];

            //                    Ham_(row_,col_) += val_temp;

            //                }
            //            }



            //Magnetic field: Tetrahedron, temporary
            //            for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){

            //                for(int spin1=0;spin1<2;spin1++){
            //                        row_ = alpha + spin1*(UnitCellSize_x*UnitCellSize_y);
            //                        col_=row_;

            //                        Ham_(row_,col_) += 0.05*(1.0 - (2.0*spin1));

            //                }
            //            }


            //Interaction:Hartree
            //On-site Hartree
            for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                for(int spin=0;spin<2;spin++){
                    spin_bar = (spin +1 )%2;
                    row_ = alpha + spin*(UnitCellSize_x*UnitCellSize_y);
                    col_=row_;

                    row_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin_bar*(UnitCellSize_x*UnitCellSize_y);
                    col_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin_bar*(UnitCellSize_x*UnitCellSize_y);
                    index_OP = SI_to_ind[col_OP + row_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                    Ham_(row_,col_) += Parameters_.U0*OPs_.value[index_OP];
                }
            }

            //Intra-Cell Hartree
            for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                for(int spin=0;spin<2;spin++){
                    for(int alpha_p=0;alpha_p<UnitCellSize_x*UnitCellSize_y;alpha_p++){
                        if(alpha != alpha_p){
                            for(int spin_p=0;spin_p<2;spin_p++){

                                row_ = alpha_p + spin_p*(UnitCellSize_x*UnitCellSize_y);
                                col_=row_;

                                row_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin*(UnitCellSize_x*UnitCellSize_y);
                                col_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin*(UnitCellSize_x*UnitCellSize_y);
                                index_OP = SI_to_ind[col_OP + row_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                Ham_(row_,col_) += IntraCell_U(alpha, alpha_p)*OPs_.value[index_OP];
                            }
                        }
                    }
                }
            }

            //Inter-Cell Hartree
            for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                for(int spin=0;spin<2;spin++){
                    for(int alpha_p=0;alpha_p<UnitCellSize_x*UnitCellSize_y;alpha_p++){
                        for(int spin_p=0;spin_p<2;spin_p++){

                            row_ = alpha_p + spin_p*(UnitCellSize_x*UnitCellSize_y);
                            col_=row_;

                            row_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin*(UnitCellSize_x*UnitCellSize_y);
                            col_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin*(UnitCellSize_x*UnitCellSize_y);
                            index_OP = SI_to_ind[col_OP + row_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                            Ham_(row_,col_) += U_Bar(alpha, alpha_p)*OPs_.value[index_OP];
                        }
                    }
                }
            }



            //Interaction:Fock
            if(!Parameters_.Just_Hartree){

                //Onsite-Fock
                if( (Parameters_.FockType=="Onsite" || Parameters_.FockType=="Onsite_Intra") || Parameters_.FockType=="Onsite_Intra_Inter" ){
                    for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                        for(int spin=0;spin<2;spin++){
                            spin_bar = (spin +1 )%2;
                            row_ = alpha + spin*(UnitCellSize_x*UnitCellSize_y);
                            col_= alpha + spin_bar*(UnitCellSize_x*UnitCellSize_y);

                            row_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin_bar*(UnitCellSize_x*UnitCellSize_y);
                            col_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin*(UnitCellSize_x*UnitCellSize_y);
                            if(col_OP>=row_OP){
                                index_OP = SI_to_ind[col_OP + row_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_value=OPs_.value[index_OP];
                            }
                            else{
                                index_OP = SI_to_ind[row_OP + col_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                OP_value= conj(OPs_.value[index_OP]);
                            }

                            Ham_(row_,col_) += -1.0*Parameters_.U0*OP_value;

                        }
                    }
                }

                //Intra-Fock
                if( Parameters_.FockType=="Onsite_Intra" || Parameters_.FockType=="Onsite_Intra_Inter" ){

                    for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                        for(int spin=0;spin<2;spin++){
                            for(int alpha_p=0;alpha_p<UnitCellSize_x*UnitCellSize_y;alpha_p++){
                                if(alpha != alpha_p){
                                    for(int spin_p=0;spin_p<2;spin_p++){

                                        row_ = alpha_p + spin_p*(UnitCellSize_x*UnitCellSize_y);
                                        col_= alpha + spin*(UnitCellSize_x*UnitCellSize_y);

                                        col_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha_p + spin_p*(UnitCellSize_x*UnitCellSize_y);
                                        row_OP = 0*(2*UnitCellSize_x*UnitCellSize_y) + alpha + spin*(UnitCellSize_x*UnitCellSize_y);

                                        if(col_OP>=row_OP){
                                            index_OP = SI_to_ind[col_OP + row_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                            OP_value=OPs_.value[index_OP];
                                        }
                                        else{
                                            index_OP = SI_to_ind[row_OP + col_OP*(2*ncells_*UnitCellSize_x*UnitCellSize_y)];
                                            OP_value= conj(OPs_.value[index_OP]);
                                        }

                                        Ham_(row_,col_) += -1.0*IntraCell_U(alpha, alpha_p)*OP_value;
                                    }
                                }
                            }
                        }
                    }
                }

                if(Parameters_.FockType=="Onsite_Intra_Inter" ){

                    for(int alpha=0;alpha<UnitCellSize_x*UnitCellSize_y;alpha++){
                        for(int spin=0;spin<2;spin++){
                            for(int alpha_p=0;alpha_p<UnitCellSize_x*UnitCellSize_y;alpha_p++){
                                for(int spin_p=0;spin_p<2;spin_p++){

                                    row_ = alpha_p + spin_p*(UnitCellSize_x*UnitCellSize_y);
                                    col_= alpha + spin*(UnitCellSize_x*UnitCellSize_y);

                                    Ham_(row_,col_) += -1.0*V_fock(alpha,spin,alpha_p,spin_p,k1,k2);
                                }
                            }
                        }
                    }
                }




            }


            char Dflag='V';
            Diagonalize(Dflag);


            for(int row=0;row<2*UnitCellSize_x*UnitCellSize_y;row++){
                Eigenvalues_[2*UnitCellSize_x*UnitCellSize_y*k_index + row]=eigs_[row];
                for(int col=0;col<2*UnitCellSize_x*UnitCellSize_y;col++){
                    Eigvectors_[2*UnitCellSize_x*UnitCellSize_y*k_index + col][row]=Ham_(row,col);
                }
            }
        }
    }

}






void Kspace_calculation_TL::SelfConsistency(){


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

            mu_=chemicalpotential(Parameters_.Total_Particles);
            Get_new_OPs_and_error();
            Get_Energies();


            file_out_progress<<setprecision(15)<<iter<<"   "<<OP_error_<<"   "<<mu_<<"     "<<E_class<<"   "<<E_quant<<"    "<<endl;
            //        for(int OP_no=0;OP_no<6;OP_no++){
            //            file_out_progress<<OPs_[OP_no].real()<<"    "<<OPs_[OP_no].imag()<<"    ";
            //        }
            //        file_out_progress<<endl;

            if(Parameters_.Anderson_Mixing){
                Update_OrderParameters_AndersonMixing(iter);
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
        Calculate_ChernNumbers();
        Create_Current_Oprs();


        //J_KE_e1[0].print();
        //J_KE_e1[1].print();
        Hall_conductance();
        Arranging_spectrum();
        mu_=chemicalpotential(Parameters_.Total_Particles);
        cout<<"mu = "<<mu_<<endl;
        cout<<"energies shown below , Eclass and Equant:"<<endl;
        Get_new_OPs_and_error();
        Get_Energies();
        cout<<E_class<<"   "<<E_quant<<"    "<<endl;

        Get_spin_resolved_local_densities();
        Get_local_spins();
        Calculate_Nw();

        string Eigenvalues_fl_out = "Eigenvalues.txt";
        ofstream Eigenvalues_file_out(Eigenvalues_fl_out.c_str());
        for(int ie=0;ie<Eigenvalues_.size();ie++){
            Eigenvalues_file_out<<ie<<"  "<<Eigenvalues_[ie]<<endl;
        }


        //        for(int Ntimes=1;Ntimes<=4;Ntimes++){
        //            cout<<"Ntimes = "<<Ntimes<<" Spectrum is being created"<<endl;
        //         Create_Kspace_Spectrum_in_ntimes_brillouin_zone(Ntimes);
        //            cout<<"Ntimes ="<<Ntimes<<" Spectrum done"<<endl;
        //         Calculate_ChernNumbers_in_ntimes_brillouin_zone(Ntimes);
        //         cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
        //        }


        string File_Out_Local_OP = Parameters_.File_OPs_out + "_Temp"+string(temp_char) + ".txt";
        ofstream file_out_Local_OP(File_Out_Local_OP.c_str());
        file_out_Local_OP<<"#row col OParams_[row][col]"<<endl;

        for(int alpha=0;alpha<OPs_.value.size();alpha++){
            file_out_Local_OP<<OPs_new_.rows[alpha]<<setw(15)<<OPs_new_.columns[alpha]<<setw(15)<<"  "<<OPs_new_.value[alpha]<<endl;
        }



        //Getting ready for next Temperature
        //        Parameters_.Read_OPs=true;
        //        Parameters_.File_OPs_in = Parameters_.File_OPs_out + "_Temp"+string(temp_char) + ".txt";
        //        Initialize();
    }


}

void Kspace_calculation_TL::Diagonalize(char option){

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


void Kspace_calculation_TL::Update_OrderParameters_AndersonMixing(int iter){
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

    assert(NewInd_to_OldInd.size()==OPs_.value.size() + (OPs_.value.size() - 2*n_orbs_*UnitCellSize_x*UnitCellSize_y));
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

void Kspace_calculation_TL::Perform_SVD(Matrix<double> & A_, Matrix<double> & VT_, Matrix<double> & U_, vector<double> & Sigma_){


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


void Kspace_calculation_TL::Perform_SVD_complex(Matrix<complex<double>> & A_, Matrix<complex<double>> & VT_, Matrix<complex<double>> & U_, vector<double> & Sigma_){


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
