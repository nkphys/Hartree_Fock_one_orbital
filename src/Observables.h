#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#include "Hamiltonian.h"
#include "tensor_type.h"
#include "functions.h"
#include <algorithm>

#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#define PI acos(-1.0)

//n, a, lda, ipvt, work, lwork, info
extern "C" void   zgetri_(int *,std::complex<double> *, int *, int *, std::complex<double> *, int *, int *);

//zgetrf_ (&n, &n, &(_TEMP(0,0)), &n, &(ipvt[0]), &info);
extern "C" void   zgetrf_(int *,int *, std::complex<double> *, int *, int *, int *);


//zhetri (character UPLO, integer N, complex*16 dimension( lda, * ) A, integer LDA,
//integer IPIV, complex*16 dimension( * ) WORK, integer INFO)
extern "C" void   zhetri_(char *, int *, std::complex<double> *, int *, int *, std::complex<double> *, int *);


//zhetrf(uplo,n,a,lda,ipiv,work,lwork,info)
extern "C" void   zhetrf_(char *, int *, std::complex<double> *, int *, int *, std::complex<double> *, int *, int *);

//dsptrf(character UPLO, integer N, double precision dimension( * ) AP
//       , integer dimension( * ) IPIV, integer INFO)
extern "C" void dsptrf_(char *, int *, double *, int *, int *);

//dsptri	(character 	UPLO, integer N, double precision dimension( * ) AP,
//            integer dimension( * ) IPIV, double precision dimension( * ) WORK,
//            integer INFO)
extern "C" void dsptri_(char *, int *, double *, int *, double *, int *);

//dgesdd	(character JOBZ, integer M, integer N, double precision  dimension( lda, * ) A,
//             integer LDA, double precision dimension( * ) S, double precision dimension( ldu, * ) U,
//             integer LDU, double precision dimension( ldvt, * ) VT, integer LDVT,
//             double precision dimension( * ) WORK, integer LWORK, integer dimension( * ) IWORK,
//             integer INFO)
extern "C" void dgesdd_ (char *, int *, int *, double *, int *, double *, double *, int *, double *, int*,
                         double *, int *, int *, int *);


extern "C" void zgesdd_ (char *, int *, int *, std::complex<double> *, int *, double *, std::complex<double> *, int *, std::complex<double> *, int*,
                         std::complex<double> *, int *, double * , int *, int *);


class Observables{
public:

    Observables(Parameters& Parameters__, Coordinates& Coordinates__,
                MFParams& MFParams__, Hamiltonian& Hamiltonian__)
        : Parameters_(Parameters__), Coordinates_(Coordinates__), MFParams_(MFParams__),
          Hamiltonian_(Hamiltonian__),
          ns_(Parameters_.ns)
    {
        Initialize();
    }

    void Initialize();

    void Calculate_Nw();
    double Lorentzian(double x, double brd);
    complex<double> DOT_P(Mat_1_Complex_doub left, Mat_1_Complex_doub right);
    double DOT_P(Mat_1_doub left, Mat_1_doub right);

    void Total_Energy_Average(double Curr_QuantE, double CurrE);

    void Calculate_Local_n_orb_resolved();
    void Calculate_Local_spins_resolved();
    void Calculate_Order_Params();
    void Get_OrderParameters_diffs();
    void Update_OrderParameters(int iter);
    void Hartree_Filter();
    void Update_OrderParameters_Second_Broyden(int iter_count);
    void Update_OrderParameters_AndersonMixing(int iter);
    void Invert_Beta();
    void Invert_Beta_double();
    void Calculate_Single_Particle_Density_Matrix();
    void Calculate_SpinSpincorrelations_Smartly();
    void Calculate_DenDencorrelations_Smartly();
    void Perform_SVD(Matrix<double> & A_, Matrix<double> & VT_, Matrix<double> & U_, vector<double> & Sigma_);
    void Perform_SVD_complex(Matrix<complex<double>> & A_, Matrix<complex<double>> & VT_, Matrix<complex<double>> & U_, vector<double> & Sigma_);


    complex<double> Two_particle_Den_Mat(int _alpha, int _beta, int _gamma, int _delta);


    double Omega(int i);

    double Avg_local_Nup, Avg_local_Ndn;
    Matrix<complex<double>> Transformation;
    Matrix<double> SiSj_Mean_, SiSj_square_Mean_;
    double Nematic_order_mean_, Nematic_order_square_mean_;
    Parameters& Parameters_;
    Coordinates& Coordinates_;
    MFParams& MFParams_;
    Hamiltonian& Hamiltonian_;
    int lx_,ly_,ns_;
    double dosincr_,tpi_;
    vector<double> nia_,nib_,nic_;
    Matrix<double> SiSj_,dos;
    vector<double> sx_,sy_,sz_;
    Matrix_COO_Complex OParams;
    Mat_1_doub Local_n_orb_resolved;
    Mat_4_Complex_doub F_Exciton;

    Mat_2_Complex_doub SP_Density_Matrix;

    // Declare Fields
    Matrix<double> Sz_obs, Sx_obs, Sy_obs;
    Matrix<double> Local_density_obs;
    double Error_OP_;
    double AVG_Total_Energy, AVG_Total_Energy_sqr;


    //Declare Broyden_Mixing vectors, if Mapped_ro_real==true
    vector<double> F_n; //F_n=x_n_out - x_n_in []
    vector<double> F_nm1;
    vector<double> DeltaF_n; //DeltaF_n=F_n - F-nm1;
    vector<double> Delta_x_n; //Delta_x_n= x_n_in - x_nm1_in;
    Mat_2_doub Jinv_n;
    Mat_2_doub Jinv_np1;



    //Declare Broyden_Mixing vectors, if Mapped_ro_real==false
    vector<complex<double>> F_n_; //F_n=x_n_out - x_n_in []
    vector<complex<double>> F_nm1_;
    vector<complex<double>> DeltaF_n_; //DeltaF_n=F_n - F-nm1;
    vector<complex<double>> Delta_x_n_; //Delta_x_n= x_n_in - x_nm1_in;
    Mat_2_Complex_doub Jinv_n_;
    Mat_2_Complex_doub Jinv_np1_;

    //Declarations for Second Broyden Method
    Mat_2_doub _Delta_F;
    Mat_2_doub _u;
    Mat_2_doub _A;
    Mat_1_doub _Fm, _Fm_minus1;
    Mat_1_doub _Delta_OPm, _Delta_OPm_minus1;
    Matrix<double> _Beta;
    Mat_1_doub _cm, _gammam;
    double w_minus1;
    Mat_1_doub w;


    //Declarations for Anderson Mixing
    Mat_1_doub x_km1_, x_k_, Del_x_km1;
    Mat_1_doub f_k_, f_km1_, Del_f_km1;
    Mat_1_doub xbar_k_, fbar_k_, gamma_k_, x_kp1_;
    Matrix<double> X_mat, F_mat;

};
/*
 * ***********
 *  Functions in Class Observables ------
 *  ***********
*/


void Observables::Update_OrderParameters_AndersonMixing(int iter){

    bool with_SVD=false;
    int Offset_;
    int m_;
    int row_, col_;
    int old_ind;
    int OP_size;
    Mat_1_int NewInd_to_OldInd;
    NewInd_to_OldInd.clear();
    for(int ind=0;ind<MFParams_.OParams_.value.size();ind++){
        NewInd_to_OldInd.push_back(ind);
    }

    for(int ind=0;ind<MFParams_.OParams_.value.size();ind++){
        row_=MFParams_.OParams_.rows[ind];
        col_=MFParams_.OParams_.columns[ind];
        if(row_!=col_){
            NewInd_to_OldInd.push_back(ind);
        }
    }

    assert(NewInd_to_OldInd.size()==MFParams_.OParams_.value.size() + (MFParams_.OParams_.value.size() - 2*ns_));
    OP_size=NewInd_to_OldInd.size();


    if(iter==0){
        //        cout<<"Anderson mixing for iter "<<iter<<endl;

        x_k_.clear();x_k_.resize(OP_size);
        for(int i=0;i<OP_size;i++){
            if(i<MFParams_.OParams_.value.size()){
                x_k_[i] = MFParams_.OParams_.value[NewInd_to_OldInd[i]].real();
            }
            else{
                x_k_[i] = MFParams_.OParams_.value[NewInd_to_OldInd[i]].imag();
            }
        }

        f_k_.clear();f_k_.resize(OP_size);
        for(int i=0;i<OP_size;i++){
            if(i<OParams.value.size()){
                f_k_[i] = OParams.value[NewInd_to_OldInd[i]].real() - MFParams_.OParams_.value[NewInd_to_OldInd[i]].real();
            }
            else{
                f_k_[i] = OParams.value[NewInd_to_OldInd[i]].imag()- MFParams_.OParams_.value[NewInd_to_OldInd[i]].imag();
            }
        }
        assert(OParams.value.size() == MFParams_.OParams_.value.size());

        //f_k = OParams.value;
        //x_k = MFParams_.OParams_.value;

        for(int ind=0;ind<MFParams_.OParams_.value.size();ind++){
            MFParams_.OParams_.value[ind] = (1-Parameters_.alpha_OP)*MFParams_.OParams_.value[ind]
                    + Parameters_.alpha_OP*OParams.value[ind];
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
            if(i<MFParams_.OParams_.value.size()){
                x_k_[i] = MFParams_.OParams_.value[NewInd_to_OldInd[i]].real();
            }
            else{
                x_k_[i] = MFParams_.OParams_.value[NewInd_to_OldInd[i]].imag();
            }
        }

        f_k_.clear();f_k_.resize(OP_size);
        for(int i=0;i<OP_size;i++){
            if(i<OParams.value.size()){
                f_k_[i] = OParams.value[NewInd_to_OldInd[i]].real() - MFParams_.OParams_.value[NewInd_to_OldInd[i]].real();
            }
            else{
                f_k_[i] = OParams.value[NewInd_to_OldInd[i]].imag() - MFParams_.OParams_.value[NewInd_to_OldInd[i]].imag();
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

            //cout<<"here 2"<<endl;
            Perform_SVD(A_,VT_,U_,Sigma_);
            //cout<<"here 2.5"<<endl;

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
            if(i<MFParams_.OParams_.value.size()){
                MFParams_.OParams_.value[NewInd_to_OldInd[i]].real( x_kp1_[i]);
            }
            else{
                MFParams_.OParams_.value[NewInd_to_OldInd[i]].imag( x_kp1_[i]);
            }
        }


        //---saving arrays for next iteration-----
        x_km1_=x_k_;
        f_km1_=f_k_;

    }


}


void Observables::Perform_SVD(Matrix<double> & A_, Matrix<double> & VT_, Matrix<double> & U_, vector<double> & Sigma_){


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


    vector<double> work(3);
    int info;
    int lwork= -1;
    vector<int> iwork(8*min(m,n));

    // query:
    dgesdd_(&jobz, &m, &n, &(A_(0,0)),&lda, &(Sigma_[0]),&(U_(0,0)), &ldu, &(VT_(0,0)), &ldvt,
            &(work[0]), &lwork, &(iwork[0]), &info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0]));
    work.resize(lwork);
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

    // Ham_.print();



}


void Observables::Perform_SVD_complex(Matrix<complex<double>> & A_, Matrix<complex<double>> & VT_, Matrix<complex<double>> & U_, vector<double> & Sigma_){


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

void Observables::Invert_Beta(){
    /*
    int n=_Beta.n_row();
    int lda=_Beta.n_col();
    int info;
    vector<int> ipvt;
    ipvt.resize(_Beta.n_col());
    vector<complex<double>> work(n);
    int lwork= n;


    char uplo='U';
    zhetrf_(&uplo, &n, &(_Beta(0,0)),&lda,&(ipvt[0]),&(work[0]),&lwork,&info);
    //cout<<"FACTORIZATION OF MATRIX:"<<endl;
    //_TEMP.print();

    zhetri_(&uplo, &n, &(_Beta(0,0)),&lda,&(ipvt[0]),&(work[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("Inverse: zgetri: failed with info!=0.\n");
    }

    //cout<<"INVERSE OF MATRIX:"<<endl;
    //_TEMP.print();
*/
}

void Observables::Invert_Beta_double(){

    int n=_Beta.n_row();
    int info;
    vector<int> ipvt;
    ipvt.resize(_Beta.n_col());
    vector<double> work(n);

    char uplo='U';
    dsptrf_(&uplo, &n, &(_Beta(0,0)),&(ipvt[0]),&info);
    //cout<<"FACTORIZATION OF MATRIX:"<<endl;
    //_TEMP.print();

    dsptri_(&uplo, &n, &(_Beta(0,0)),&(ipvt[0]),&(work[0]),&info);
    if (info!=0) {
        std::cerr<<"info="<<info<<"\n";
        perror("Inverse: dsptri: failed with info!=0.\n");
    }

    //cout<<"INVERSE OF MATRIX:"<<endl;
    //_TEMP.print();

}



complex<double> Observables::DOT_P(Mat_1_Complex_doub left, Mat_1_Complex_doub right){
    complex<double> temp_;
    temp_=zero_complex;

    assert(left.size()==right.size());

    for(int i=0;i<left.size();i++){
        temp_ += conj(left[i])*right[i];
    }
    return temp_;

}

double Observables::DOT_P(Mat_1_doub left, Mat_1_doub right){
    double temp_;
    temp_=0.0;

    assert(left.size()==right.size());

    for(int i=0;i<left.size();i++){
        temp_ += left[i]*right[i];
    }
    return temp_;
}

void Observables::Get_OrderParameters_diffs(){

    Error_OP_=0.0;

    for(int ind=0;ind<MFParams_.OParams_.value.size();ind++){
        Error_OP_ += abs(OParams.value[ind] - MFParams_.OParams_.value[ind])*
                abs(OParams.value[ind] - MFParams_.OParams_.value[ind]);
    }

    Error_OP_=sqrt(Error_OP_);

}


void Observables::Hartree_Filter(){

    int row_,col_;
    for(int ind=0;ind<MFParams_.OParams_.value.size();ind++){
        row_=MFParams_.OParams_.rows[ind];
        col_=MFParams_.OParams_.columns[ind];
        if(row_!=col_){
            MFParams_.OParams_.value[ind] = zero_complex;
        }
    }


}

void Observables::Update_OrderParameters(int iter){


    int UP_, DOWN_;
    UP_=0;DOWN_=1;
    bool Mapped_to_real=true;
    bool Good_Broyden=true; //If Good is false it runs Bad broyden

    //Simple mixing
    double alpha_OP=Parameters_.alpha_OP;

    if(Parameters_.Simple_Mixing==true){

        if(iter==0){
            cout<<"Using Simple Mixing to gain Self-Consistency"<<endl;
        }

        for(int ind=0;ind<MFParams_.OParams_.value.size();ind++){
            MFParams_.OParams_.value[ind] = (1-alpha_OP)*MFParams_.OParams_.value[ind]
                    + alpha_OP*OParams.value[ind];
        }

        //      Parameters_.mu_old = (1-alpha_OP)*Parameters_.mu_old + alpha_OP*Parameters_.mus;

    }



    /*
    //Declared in initialize function();
    vector<double> F_n; //F_n=x_n_out - x_n_in []
    vector<double> F_nm1;
    vector<double> DeltaF_n; //DeltaF_n=F_n - F-nm1;
    vector<double> Delta_x_n; //Delta_x_n= x_n_in - x_nm1_in;
    Mat_2_doub Jinv_n;
    Mat_2_doub Jinv_np1;
     */

    //cout<<MFParams_.OParams_.value.size()<<endl;
    int Offset_;
    int row_, col_;
    int old_ind;
    Mat_1_int NewInd_to_OldInd;
    NewInd_to_OldInd.clear();
    for(int ind=0;ind<MFParams_.OParams_.value.size();ind++){
        NewInd_to_OldInd.push_back(ind);
    }

    for(int ind=0;ind<MFParams_.OParams_.value.size();ind++){
        row_=MFParams_.OParams_.rows[ind];
        col_=MFParams_.OParams_.columns[ind];
        if(row_!=col_){
            NewInd_to_OldInd.push_back(ind);
        }
    }

    //    int no_local_ops;
    //    if(Parameters_.Restricted_HF){
    //        no_local_ops=Parameters_.CDW_Ansatz_sites.size();
    //    }
    //    else{
    //         no_local_ops=ns_;
    //    }

    assert(NewInd_to_OldInd.size()==MFParams_.OParams_.value.size() + (MFParams_.OParams_.value.size() - 2*ns_));






    if(Parameters_.Broyden_Mixing==true){
        //assert(false);

        if(Mapped_to_real){
            vector<double> vec_V, vec_U;
            double Denominator_;
            vector<double> vec_L;
            Offset_=MFParams_.OParams_.value.size();
            int size_vecs;
            size_vecs=NewInd_to_OldInd.size();
            if(iter==0){

                //Initializing:
                F_n.resize(size_vecs); //F_n=x_n_out - x_n_in []
                F_nm1.resize(size_vecs);
                DeltaF_n.resize(size_vecs); //DeltaF_n=F_n - F-nm1;
                Delta_x_n.resize(size_vecs); //Delta_x_n= x_n_in - x_nm1_in;
                Jinv_n.resize(size_vecs);
                Jinv_np1.resize(size_vecs);

                _Fm.resize((size_vecs));
                _Delta_OPm.resize(size_vecs);
                _Fm_minus1.resize(size_vecs);
                _Delta_OPm_minus1.resize(size_vecs);

                for(int i=0;i<size_vecs;i++){
                    Jinv_n[i]. resize(size_vecs);
                    Jinv_np1[i]. resize(size_vecs);
                }
                //----------------------------------
                cout<<"Using Broyden Mixing to gain Self-Consistency"<<endl;


                //Get Jinv_np1
                for(int i=0;i<size_vecs;i++){
                    for(int j=0;j<size_vecs;j++){
                        if(i==j){
                            Jinv_np1[i][j]=-1.0*alpha_OP;
                        }
                        else{
                            Jinv_np1[i][j]=0.0;
                        }
                    }}

                //Get F_n
                for(int new_ind=0;new_ind<size_vecs;new_ind++){
                    old_ind=NewInd_to_OldInd[new_ind];
                    if(new_ind<OParams.value.size()){
                        F_n[new_ind]= (OParams.value[old_ind] - MFParams_.OParams_.value[old_ind]).real();
                    }
                    else{
                        F_n[new_ind] =  (OParams.value[old_ind] - MFParams_.OParams_.value[old_ind]).imag();
                    }
                }



                for(int i=0;i<size_vecs;i++){
                    Delta_x_n[i] =0.0;
                    for(int j=0;j<size_vecs;j++){
                        Delta_x_n[i] +=  -1.0*Jinv_np1[i][j]*F_n[j];
                    }
                }


                for(int new_ind=0;new_ind<size_vecs;new_ind++){
                    old_ind=NewInd_to_OldInd[new_ind];
                    if(new_ind<OParams.value.size()){
                        MFParams_.OParams_.value[old_ind].real( MFParams_.OParams_.value[old_ind].real() + Delta_x_n[new_ind]);
                    }
                    else{
                        MFParams_.OParams_.value[old_ind].imag( MFParams_.OParams_.value[old_ind].imag() + Delta_x_n[new_ind]);
                    }
                }

                //Copy Jinv_np1 to Jinv_n
                Jinv_n = Jinv_np1;

                //Copy F_n to F_nm1
                F_nm1=F_n;
            }
            else{
                //Get F_n
                for(int new_ind=0;new_ind<size_vecs;new_ind++){
                    old_ind=NewInd_to_OldInd[new_ind];
                    if(new_ind<OParams.value.size()){
                        F_n[new_ind]= (OParams.value[old_ind] - MFParams_.OParams_.value[old_ind]).real();
                    }
                    else{
                        F_n[new_ind] =  (OParams.value[old_ind] - MFParams_.OParams_.value[old_ind]).imag();
                    }
                }

                //Get DeltaF_n
                for (int i=0;i<size_vecs;i++){
                    DeltaF_n[i] = F_n[i] - F_nm1[i];
                }

                //Get vec_V = Jinv_n*DeltaF_n
                vec_V.clear();
                vec_V.resize(size_vecs);

                for(int i=0;i<size_vecs;i++){
                    vec_V[i] =0.0;
                    for(int j=0;j<size_vecs;j++){
                        vec_V[i] += Jinv_n[i][j]*DeltaF_n[j];
                    }
                }

                //Get vec_U = Delta_x_n^dagg*Jinv_n
                vec_U.clear();
                vec_U.resize(size_vecs);

                for(int i=0;i<size_vecs;i++){
                    vec_U[i] =0.0;
                    for(int j=0;j<size_vecs;j++){
                        vec_U[i] += Delta_x_n[j]*Jinv_n[j][i];
                    }
                }

                // Get Denominator_=<Delta_x_n|vec_V>
                Denominator_=0.0;
                for(int i=0;i<size_vecs;i++){
                    if(Good_Broyden){
                        Denominator_ +=Delta_x_n[i]*vec_V[i]; //Good
                    }
                    else{
                        Denominator_ +=DeltaF_n[i]*DeltaF_n[i];  //"Bad broyden"}
                    }
                }

                //Get vec_L=  Delta_x_n - vec_V;
                vec_L.clear();
                vec_L.resize(size_vecs);
                for(int i=0;i<size_vecs;i++){
                    vec_L[i] = Delta_x_n[i] - vec_V[i];
                }


                //Get Mat_Temp [Remember to clear later on];
                Mat_2_doub Mat_Temp;
                Mat_Temp.resize(size_vecs);
                for(int i=0;i<size_vecs;i++){
                    Mat_Temp[i].resize(size_vecs);
                    for(int j=0;j<size_vecs;j++){
                        if(Good_Broyden){
                            Mat_Temp[i][j] = (vec_L[i]*vec_U[j])/(Denominator_); //Good
                        }
                        else{
                            Mat_Temp[i][j] = (vec_L[i]*DeltaF_n[j])/(Denominator_); //Bad
                        }
                    }
                }


                //Get Jinv_np1
                for(int i=0;i<size_vecs;i++){
                    for(int j=0;j<size_vecs;j++){
                        Jinv_np1[i][j]  = Jinv_n[i][j]  + Mat_Temp[i][j];
                    }
                }

                for(int i=0;i<size_vecs;i++){
                    Delta_x_n[i] =0.0;
                    for(int j=0;j<size_vecs;j++){
                        Delta_x_n[i] +=  -1.0*Jinv_np1[i][j]*F_n[j];
                    }
                }


                for(int new_ind=0;new_ind<size_vecs;new_ind++){
                    old_ind=NewInd_to_OldInd[new_ind];
                    if(new_ind<OParams.value.size()){
                        MFParams_.OParams_.value[old_ind].real( MFParams_.OParams_.value[old_ind].real() + Delta_x_n[new_ind]);
                    }
                    else{
                        MFParams_.OParams_.value[old_ind].imag( MFParams_.OParams_.value[old_ind].imag() + Delta_x_n[new_ind]);
                    }
                }


                //Copy Jinv_np1 to Jinv_n
                Jinv_n = Jinv_np1;

                //Copy F_n to F_nm1
                F_nm1=F_n;


                //Clear Mat_Temp
                for(int i=0;i<size_vecs;i++){
                    Mat_Temp[i].clear();
                }
                Mat_Temp.clear();
            }

        }

        else{
            vector<complex<double>> vec_V, vec_U;
            complex<double> Denominator_;
            vector<complex<double>> vec_L;
            int size_vecs;
            size_vecs=MFParams_.OParams_.value.size();
            if(iter==0){

                //Initializing:
                F_n_.resize(size_vecs); //F_n=x_n_out - x_n_in []
                F_nm1_.resize(size_vecs);
                DeltaF_n_.resize(size_vecs); //DeltaF_n=F_n - F-nm1;
                Delta_x_n_.resize(size_vecs); //Delta_x_n= x_n_in - x_nm1_in;
                Jinv_n_.resize(size_vecs);
                Jinv_np1_.resize(size_vecs);

                //                _Fm.resize((size_vecs));
                //                _Delta_OPm.resize(size_vecs);
                //                _Fm_minus1.resize(size_vecs);
                //                _Delta_OPm_minus1.resize(size_vecs);

                for(int i=0;i<size_vecs;i++){
                    Jinv_n_[i]. resize(size_vecs);
                    Jinv_np1_[i]. resize(size_vecs);
                }
                //----------------------------------
                cout<<"Using Broyden Mixing to gain Self-Consistency"<<endl;


                //Get Jinv_np1
                for(int i=0;i<size_vecs;i++){
                    for(int j=0;j<size_vecs;j++){
                        if(i==j){
                            Jinv_np1_[i][j]=-1.0*alpha_OP;
                        }
                        else{
                            Jinv_np1_[i][j]=0.0;
                        }
                    }}

                //Get F_n
                for(int new_ind=0;new_ind<size_vecs;new_ind++){
                    F_n_[new_ind]= OParams.value[new_ind] - MFParams_.OParams_.value[new_ind];
                }



                for(int i=0;i<size_vecs;i++){
                    Delta_x_n_[i] =0.0;
                    for(int j=0;j<size_vecs;j++){
                        Delta_x_n_[i] +=  -1.0*Jinv_np1_[i][j]*F_n_[j];
                    }
                }


                for(int new_ind=0;new_ind<size_vecs;new_ind++){
                    MFParams_.OParams_.value[new_ind] +=  MFParams_.OParams_.value[new_ind] + Delta_x_n_[new_ind];
                }

                //Copy Jinv_np1 to Jinv_n
                Jinv_n_ = Jinv_np1_;

                //Copy F_n to F_nm1
                F_nm1_=F_n_;
            }
            else{
                //Get F_n
                for(int new_ind=0;new_ind<size_vecs;new_ind++){
                    F_n_[new_ind]= (OParams.value[new_ind] - MFParams_.OParams_.value[new_ind]);
                }

                //Get DeltaF_n
                for (int i=0;i<size_vecs;i++){
                    DeltaF_n_[i] = F_n_[i] - F_nm1_[i];
                }

                //Get vec_V = Jinv_n*|DeltaF_n>
                vec_V.clear();
                vec_V.resize(size_vecs);

                for(int i=0;i<size_vecs;i++){
                    vec_V[i] =0.0;
                    for(int j=0;j<size_vecs;j++){
                        vec_V[i] += Jinv_n_[i][j]*DeltaF_n_[j];
                    }
                }

                //Get |vec_U> = (J_inv_n)^dagger|Delta_x_n>
                vec_U.clear();
                vec_U.resize(size_vecs);

                for(int i=0;i<size_vecs;i++){
                    vec_U[i] =0.0;
                    for(int j=0;j<size_vecs;j++){
                        vec_U[i] += (conj(Jinv_n_[j][i]))*Delta_x_n_[j];
                    }
                }

                // Get Denominator_=<Delta_x_n|vec_V>
                Denominator_=0.0;
                for(int i=0;i<size_vecs;i++){
                    if(Good_Broyden){
                        Denominator_ +=conj(Delta_x_n_[i])*vec_V[i]; //Good Broyden
                    }
                    else{
                        Denominator_ +=conj(DeltaF_n_[i])*DeltaF_n_[i];  //"Bad broyden"
                    }
                }


                //Get vec_L=  Delta_x_n - vec_V;
                vec_L.clear();
                vec_L.resize(size_vecs);
                for(int i=0;i<size_vecs;i++){
                    vec_L[i] = Delta_x_n_[i] - vec_V[i];
                }


                //Get Mat_Temp [Remember to clear later on];
                Mat_2_Complex_doub Mat_Temp;
                Mat_Temp.resize(size_vecs);
                for(int i=0;i<size_vecs;i++){
                    Mat_Temp[i].resize(size_vecs);
                    for(int j=0;j<size_vecs;j++){
                        if(Good_Broyden){
                            Mat_Temp[i][j] = (vec_L[i]*conj(vec_U[j]))/(Denominator_); //Good Broyden
                        }
                        else{
                            Mat_Temp[i][j] = (vec_L[i]*conj(DeltaF_n_[j]))/(Denominator_); //Bad Broyden
                        }
                    }
                }


                //Get Jinv_np1
                for(int i=0;i<size_vecs;i++){
                    for(int j=0;j<size_vecs;j++){
                        Jinv_np1_[i][j]  = Jinv_n_[i][j]  + Mat_Temp[i][j];
                    }
                }

                for(int i=0;i<size_vecs;i++){
                    Delta_x_n_[i] =0.0;
                    for(int j=0;j<size_vecs;j++){
                        Delta_x_n_[i] +=  -1.0*Jinv_np1_[i][j]*F_n_[j];
                    }
                }


                for(int new_ind=0;new_ind<size_vecs;new_ind++){
                    MFParams_.OParams_.value[new_ind] +=  MFParams_.OParams_.value[new_ind] + Delta_x_n_[new_ind];
                }


                //Copy Jinv_np1 to Jinv_n
                Jinv_n_ = Jinv_np1_;

                //Copy F_n to F_nm1
                F_nm1_=F_n_;


                //Clear Mat_Temp
                for(int i=0;i<size_vecs;i++){
                    Mat_Temp[i].clear();
                }
                Mat_Temp.clear();
            }

            for(int new_ind=0;new_ind<size_vecs;new_ind++){
                if(MFParams_.OParams_.rows[new_ind]==MFParams_.OParams_.columns[new_ind]){
                    MFParams_.OParams_.value[new_ind].imag(0.0);
                }
            }

        }

    }
    else if (Parameters_.Anderson_Mixing){
        Update_OrderParameters_AndersonMixing(iter);
    }



    //-----------------Ansatz_Filter

    if(Parameters_.Restricted_HF && false){
        int row_temp, col_temp, ind;
        complex<double> comp_temp, value_temp, ni_up, ni_dn, sz_i;

        if(Parameters_.Ansatz=="Given_CDW"){

            Parameters_.A_charge_modulation = 0.0;
            for(int site_i=0;site_i<ns_;site_i++){

                value_temp=0.0;
                for(int spin_i=0;spin_i<2;spin_i++){
                    row_temp = Coordinates_.Nc_dof(site_i, spin_i);
                    ind =MFParams_.SI_to_ind[row_temp + (2*ns_*row_temp)];
                    value_temp +=MFParams_.OParams_.value[ind];
                }
                Parameters_.A_charge_modulation += abs((value_temp -0.5)*(1.0/Parameters_.CDW_Ansatz_sites[site_i]));
            }
            Parameters_.A_charge_modulation = (1.0/ns_)*Parameters_.A_charge_modulation;

            if(Parameters_.A_charge_modulation > 0.5){
                //cout <<"Parameters_.A_charge_modulation = "<<Parameters_.A_charge_modulation<<endl;
                assert(Parameters_.A_charge_modulation <= 0.5);
                Parameters_.A_charge_modulation=0.5;
            }

            for(int site_i=0;site_i<ns_;site_i++){
                for(int site_j=0;site_j<ns_;site_j++){

                    if(site_i!=site_j){
                        if(abs(Parameters_.LongRangeInteractions[site_i][site_j])>0.0000001){

                            for(int spin_i=0;spin_i<2;spin_i++){
                                for(int spin_j=spin_i;spin_j<2;spin_j++){
                                    row_temp = Coordinates_.Nc_dof(site_i, spin_i);
                                    col_temp = Coordinates_.Nc_dof(site_j, spin_j);

                                    if(spin_i==spin_j){
                                        if(site_j>site_i){
                                            comp_temp.real(0.0);
                                            comp_temp.imag(0.0);
                                            MFParams_.OParams_.value[MFParams_.SI_to_ind[row_temp + (2*ns_*col_temp)]]=comp_temp;
                                            assert(MFParams_.OParams_.rows[MFParams_.SI_to_ind[row_temp + (2*ns_*col_temp)]]==row_temp);
                                            assert(MFParams_.OParams_.columns[MFParams_.SI_to_ind[row_temp + (2*ns_*col_temp)]]==col_temp);


                                        }
                                    }
                                    else{
                                        comp_temp.real(0.0);
                                        comp_temp.imag(0.0);
                                        MFParams_.OParams_.value[MFParams_.SI_to_ind[row_temp + (2*ns_*col_temp)]]=comp_temp;
                                        assert(MFParams_.OParams_.rows[MFParams_.SI_to_ind[row_temp + (2*ns_*col_temp)]]==row_temp);
                                        assert(MFParams_.OParams_.columns[MFParams_.SI_to_ind[row_temp + (2*ns_*col_temp)]]==col_temp);

                                    }
                                }
                            }
                        }
                    }
                    else{
                        int spin_i_, spin_j_;

                        //Get ni_up
                        spin_j_=UP_;
                        spin_i_=UP_;
                        row_temp = Coordinates_.Nc_dof(site_i, spin_i_);
                        col_temp = Coordinates_.Nc_dof(site_j, spin_j_);
                        ni_up = MFParams_.OParams_.value[MFParams_.SI_to_ind[row_temp + (2*ns_*col_temp)]];

                        //Get ni_dn
                        spin_j_=DOWN_;
                        spin_i_=DOWN_;
                        row_temp = Coordinates_.Nc_dof(site_i, spin_i_);
                        col_temp = Coordinates_.Nc_dof(site_j, spin_j_);
                        ni_dn = MFParams_.OParams_.value[MFParams_.SI_to_ind[row_temp + (2*ns_*col_temp)]];

                        //Get sz_i
                        sz_i = (ni_up -ni_dn)*0.5;
                        if(abs(sz_i.real()) > 0.5*(0.5 + Parameters_.A_charge_modulation*Parameters_.CDW_Ansatz_sites[site_i])){
                            //sz_i= sz_i*(0.5/abs(sz_i));
                            //assert(sz_i.real() <= 0.5*(0.5 + Parameters_.A_charge_modulation*Parameters_.CDW_Ansatz_sites[site_i]));
                        }
                        //

                        //UP, UP
                        value_temp = 0.5*(0.5 + Parameters_.A_charge_modulation*Parameters_.CDW_Ansatz_sites[site_i]) + sz_i;
                        spin_j_=UP_;
                        spin_i_=UP_;
                        row_temp = Coordinates_.Nc_dof(site_i, spin_i_);
                        col_temp = Coordinates_.Nc_dof(site_j, spin_j_);
                        MFParams_.OParams_.value[MFParams_.SI_to_ind[row_temp + (2*ns_*col_temp)]]=value_temp;
                        assert(MFParams_.OParams_.rows[MFParams_.SI_to_ind[row_temp + (2*ns_*col_temp)]]==row_temp);
                        assert(MFParams_.OParams_.columns[MFParams_.SI_to_ind[row_temp + (2*ns_*col_temp)]]==col_temp);

                        //DOWN, DOWN
                        value_temp = 0.5*(0.5 + Parameters_.A_charge_modulation*Parameters_.CDW_Ansatz_sites[site_i]) - sz_i;
                        spin_j_=DOWN_;
                        spin_i_=DOWN_;
                        row_temp = Coordinates_.Nc_dof(site_i, spin_i_);
                        col_temp = Coordinates_.Nc_dof(site_j, spin_j_);
                        MFParams_.OParams_.value[MFParams_.SI_to_ind[row_temp + (2*ns_*col_temp)]]=value_temp;
                        assert(MFParams_.OParams_.rows[MFParams_.SI_to_ind[row_temp + (2*ns_*col_temp)]]==row_temp);
                        assert(MFParams_.OParams_.columns[MFParams_.SI_to_ind[row_temp + (2*ns_*col_temp)]]==col_temp);
                    }
                }
            }
        }
    }
    //--------------Ansatz filter




}

void Observables::Update_OrderParameters_Second_Broyden(int iter_count){


    //DONT USE IT, NEED TO BE CHANGED for THIS 1-orbital CODE

    //For details see your own notes at
    //"https://github.com/nkphys/3_ORB_SOC_Hartree_Fock/tree/master/Notes/Modified_Broyden_Method"

    //XXXXXXXXXXLiteratureXXXXXXXXXXXXX
    //Second Broyden is used from "D. D. Johnson, Phys. Rev. B 38, 12807, 1988".
    //Look into this "https://arxiv.org/pdf/0805.4446.pdf" as well.
    //XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    //    double alpha_OP=Parameters_.alpha_OP;
    //    double normalization_;
    //    int site;
    //    int iter;
    //    Mat_1_int Offsets_;
    //    Offsets_.resize(5);
    //    Offsets_[0]=5;Offsets_[1]=10;Offsets_[2]=14;Offsets_[3]=17;Offsets_[4]=19;


    //    iter=iter_count%Parameters_.BroydenSecondMethodCounter;


    //    if(iter==0){
    //        //****Getting ready for iters>0*********
    //        _Delta_F.clear();
    //        _u.clear();
    //        _A.clear();
    //        _cm.clear();
    //        _gammam.clear();
    //        w_minus1=Parameters_.w_minus1;
    //        w.clear();
    //        //************************************

    //        cout<<"Using Modified Broyden Mixing to gain Self-Consistency"<<endl;


    //        //Get Fm
    //        for(int site=0;site<ns_;site++){
    //            for(int state1=0;state1<6;state1++){
    //                for(int state2=state1;state2<6;state2++){
    //                    if(state1==state2){
    //                        _Fm[36*(site) + state2]=(OParams[site][state1][state2] -  MFParams_.OParams_[site][state1][state2]).real();
    //                    }
    //                    else{
    //                        _Fm[36*(site) + Offsets_[state1] + (state2-state1)]=(OParams[site][state1][state2] -  MFParams_.OParams_[site][state1][state2]).real();
    //                        _Fm[36*(site) + Offsets_[state1] + (state2-state1) + (21-6)]=(OParams[site][state1][state2] -  MFParams_.OParams_[site][state1][state2]).imag();
    //                    }
    //                }
    //            }
    //        }


    //        for(int j=0;j<36*ns_;j++){
    //            _Delta_OPm[j] = alpha_OP*_Fm[j];
    //        }



    //        for(int site=0;site<ns_;site++){
    //            for(int state1=0;state1<6;state1++){
    //                for(int state2=state1;state2<6;state2++){
    //                    if(state1==state2){
    //                        MFParams_.OParams_[site][state1][state2] += one_complex*(_Delta_OPm[36*(site) + state2]);
    //                    }
    //                    else{
    //                        MFParams_.OParams_[site][state1][state2] += complex<double>(_Delta_OPm[36*(site) + Offsets_[state1] + (state2-state1)],
    //                                _Delta_OPm[36*(site) + Offsets_[state1] + (state2-state1) + (21-6)]);
    //                    }
    //                }
    //            }
    //        }

    //        //Copy Jinv_np1 to Jinv_n
    //        _Delta_OPm_minus1 = _Delta_OPm;

    //        //Copy F_n to F_nm1
    //        _Fm_minus1=_Fm;

    //    }

    //    else{

    //        w.resize(iter);
    //        w[iter-1]=Parameters_.wn;


    //        for(int site=0;site<ns_;site++){
    //            for(int state1=0;state1<6;state1++){
    //                for(int state2=state1;state2<6;state2++){
    //                    if(state1==state2){
    //                        _Fm[36*(site) + state2]=(OParams[site][state1][state2] -  MFParams_.OParams_[site][state1][state2]).real();
    //                    }
    //                    else{
    //                        _Fm[36*(site) + Offsets_[state1] + (state2-state1)]=(OParams[site][state1][state2] -  MFParams_.OParams_[site][state1][state2]).real();
    //                        _Fm[36*(site) + Offsets_[state1] + (state2-state1) + (21-6)]=(OParams[site][state1][state2] -  MFParams_.OParams_[site][state1][state2]).imag();
    //                    }
    //                }
    //            }
    //        }
    //        //******************************


    //        //Get DeltaFm/|DeltaFm|-------------------------//
    //        _Delta_F.resize(iter);
    //        _Delta_F[iter-1].resize(36*ns_);

    //        for(int i_=0;i_<36*ns_;i_++){
    //            _Delta_F[iter-1][i_]=_Fm[i_] - _Fm_minus1[i_];
    //        }


    //        normalization_=0.0;
    //        for(int i_=0;i_<36*ns_;i_++){
    //            normalization_ += ((_Delta_F[iter-1][i_])*_Delta_F[iter-1][i_]);
    //        }
    //        normalization_=sqrt(normalization_);



    //        for(int i_=0;i_<36*ns_;i_++){
    //            _Delta_F[iter-1][i_]=_Delta_F[iter-1][i_]*(1.0/normalization_);
    //        }
    //        //--------------------------------------------//


    //        //Getting Delta_n/|DeltaFm|-------------------//
    //        for(int i=0;i<36*ns_;i++){
    //            _Delta_OPm_minus1[i]=_Delta_OPm_minus1[i]*(1.0/normalization_);
    //        }
    //        //-----------------------------------------------//


    //        //Get u[iter-1]------------------------------//
    //        _u.resize(iter);
    //        _u[iter-1].resize(36*ns_);
    //        for(int i=0;i<36*ns_;i++){
    //            _u[iter-1][i] = (alpha_OP*_Delta_F[iter-1][i])  +  _Delta_OPm_minus1[i];
    //        }
    //        //-------------------------------------------------//


    //        //UPDATE _A----------------------------//

    //        Mat_1_doub temp_vec;
    //        temp_vec.resize(iter);
    //        _A.push_back(temp_vec);
    //        temp_vec.clear();

    //        for(int i=0;i<_A.size();i++){
    //            _A[i].resize(iter);
    //        }


    //        for(int i=0;i<iter;i++){
    //            if(i==(iter-1)){
    //                _A[i][i]=w[i]*w[i]*DOT_P(_Delta_F[i],_Delta_F[i]);
    //            }
    //            else{
    //                _A[iter-1][i]=w[iter-1]*w[i]*DOT_P(_Delta_F[i],_Delta_F[iter-1]);
    //                _A[i][iter-1]=w[iter-1]*w[i]*DOT_P(_Delta_F[iter-1],_Delta_F[i]);
    //            }
    //        }
    //        //---------------------------------------------//

    //        //Get Beta--------------------------------------//
    //        _Beta.resize(iter,iter);
    //        for(int i=0;i<iter;i++){
    //            for(int j=0;j<iter;j++){
    //                if(i==j){
    //                    _Beta(i,j) = (w_minus1*w_minus1) + _A[i][j];
    //                }
    //                else{
    //                    _Beta(i,j) = _A[i][j];
    //                }
    //            }

    //        }


    //        Invert_Beta_double();
    //        for(int i=0;i<_Beta.n_col();i++){
    //            for(int j=0;j<i;j++){
    //                _Beta(i,j)=_Beta(j,i);
    //            }
    //        }



    //        //-----------------------------------------------//

    //        //Get _cm-------------------------------------------//
    //        _cm.clear();
    //        _cm.resize(iter);
    //        for(int i=0;i<iter;i++){
    //            _cm[i]=w[i]*DOT_P(_Delta_F[i],_Fm);
    //        }
    //        //---------------------------------------------------//

    //        //Get _gammam------------------------------------------//
    //        _gammam.clear();
    //        _gammam.resize(iter);
    //        for(int l=0;l<iter;l++){
    //            _gammam[l]=0.0;
    //            for(int k=0;k<iter;k++){
    //                _gammam[l] += _cm[k]*_Beta(k,l);
    //            }
    //        }
    //        //--------------------------------------------------//


    //        //Get _Delta_OPm-----------------------------------------//
    //        for(int i=0;i<36*ns_;i++){
    //            _Delta_OPm[i]=0.0;
    //            for(int n=0;n<iter;n++){
    //                _Delta_OPm[i] += (-1.0)*w[n]*_gammam[n]*_u[iter-1][i];
    //            }
    //        }

    //        for(int i=0;i<36*ns_;i++){
    //            _Delta_OPm[i] += alpha_OP*_Fm[i];
    //        }
    //        //---------------------------------------------//


    //        for(int site=0;site<ns_;site++){
    //            for(int state1=0;state1<6;state1++){
    //                for(int state2=state1;state2<6;state2++){
    //                    if(state1==state2){
    //                        MFParams_.OParams_[site][state1][state2] += one_complex*(_Delta_OPm[36*(site) + state2]);
    //                    }
    //                    else{
    //                        MFParams_.OParams_[site][state1][state2] += complex<double>(_Delta_OPm[36*(site) + Offsets_[state1] + (state2-state1)],
    //                                _Delta_OPm[36*(site) + Offsets_[state1] + (state2-state1) + (21-6)]);
    //                    }
    //                }
    //            }
    //        }

    //        //Copy Jinv_np1 to Jinv_n
    //        _Delta_OPm_minus1 = _Delta_OPm;

    //        //Copy F_n to F_nm1
    //        _Fm_minus1=_Fm;

    //    }

    //    for(int site=0;site<ns_;site++){
    //        for(int state1=0;state1<6;state1++){
    //            for(int state2=state1;state2<6;state2++){
    //                if(state2!=state1){
    //                    MFParams_.OParams_[site][state2][state1] = conj(MFParams_.OParams_[site][state1][state2]);
    //                }
    //            }
    //        }
    //    }



}

void Observables::Calculate_Local_n_orb_resolved(){

    double N_total_temp=0.0;

    Local_n_orb_resolved.resize(ns_*2);

    int c1;
    for(int site=0;site<ns_;site++){
        for(int spin=0;spin<2;spin++){

            c1=Coordinates_.Nc_dof(site,spin);
            Local_n_orb_resolved[c1]=0;

            for(int n=0;n<Hamiltonian_.eigs_.size();n++){
                Local_n_orb_resolved[c1] += (
                            ( conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c1,n))*
                            (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))
                            ).real();
            }

            N_total_temp +=Local_n_orb_resolved[c1];

        }

    }

    cout<<"Total number of particles calculated = "<<N_total_temp<<endl;
    cout<<"Mu (Chemical Potential) = "<<Parameters_.mus<<endl;

    if(Parameters_.FixingMu){
        Parameters_.Total_Particles = N_total_temp;}


    Mat_1_doub Spin_den_Avg;
    Spin_den_Avg.resize(2);
    string File_Out_Local_orb_densities = "Local_spin_resolved_densities.txt";
    ofstream file_out_Local_orb_densities(File_Out_Local_orb_densities.c_str());
    file_out_Local_orb_densities<<"#site   up   dn"<<endl;
    for(int site_i=0;site_i<ns_;site_i++){

        file_out_Local_orb_densities<<site_i;
        for(int spin=0;spin<2;spin++){
            c1=Coordinates_.Nc_dof(site_i,spin);
            file_out_Local_orb_densities<<setw(15)<<Local_n_orb_resolved[c1];
            Spin_den_Avg[spin] +=Local_n_orb_resolved[c1];
        }
        file_out_Local_orb_densities<<endl;
    }

    for(int spin=0;spin<2;spin++){
        Spin_den_Avg[spin] = Spin_den_Avg[spin]*(1.0/(1.0*ns_));
        cout<<"Spin "<<spin<<" avg. den = "<<Spin_den_Avg[spin]<<endl;
    }
    Avg_local_Nup = Spin_den_Avg[0];
    Avg_local_Ndn = Spin_den_Avg[1];




}

void Observables::Calculate_Local_spins_resolved(){

    string fileoutLocalS="LocalS_HF_Snake3.txt";
    ofstream LocalS(fileoutLocalS.c_str());
    double sz, sx, sy;
    int UP_=0;
    int DOWN_=1;

    int LX_=3;


    int c1, c2;
    int site_x, site_y;
    double rx_, ry_;
    int site_;
    int ix, iy, ix_new, iy_new, i_new;
    int jx, jy, jx_new, jy_new, j_new;

    for(int site=0;site<ns_;site++){

        ix = site%LX_;
        iy =  (site-ix)/LX_;

        if(iy%2==0){
            ix_new=ix;
        }
        else{
            ix_new = (LX_-1 - ix);
        }

        site_ = ix_new + iy*(LX_);


        site_x = site_%LX_;
        site_y = (site_-site_x)/LX_;
        rx_ = ((1.0)*(site_x) +  (1.0/2.0)*(site_y));
        ry_ =  (0.0*(site_x) + (sqrt(3.0)/2.0)*(site_y));

        //sz
        sz=0.0;
        c1=Coordinates_.Nc_dof(site,UP_);
        c2=Coordinates_.Nc_dof(site,DOWN_);
        for(int n=0;n<Hamiltonian_.eigs_.size();n++){
            sz += 0.5*(
                        ( (conj(Hamiltonian_.Ham_(c1,n))*Hamiltonian_.Ham_(c1,n))
                          - (conj(Hamiltonian_.Ham_(c2,n))*Hamiltonian_.Ham_(c2,n))
                          )*
                        (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))
                        ).real();
        }

        //sx
        sx=0.0;
        c1=Coordinates_.Nc_dof(site,UP_);
        c2=Coordinates_.Nc_dof(site,DOWN_);
        for(int n=0;n<Hamiltonian_.eigs_.size();n++){
            sx += (( conj(Hamiltonian_.Ham_(c2,n))*Hamiltonian_.Ham_(c1,n)
                     )*
                   (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))
                   ).real();
        }

        //sy
        sy=0.0;
        c1=Coordinates_.Nc_dof(site,UP_);
        c2=Coordinates_.Nc_dof(site,DOWN_);
        for(int n=0;n<Hamiltonian_.eigs_.size();n++){
            sy += (( conj(Hamiltonian_.Ham_(c2,n))*Hamiltonian_.Ham_(c1,n)
                     )*
                   (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))
                   ).imag();
        }


        // LocalS<<site<<setw(15)<<sz<<setw(15)<<sx<<setw(15)<<sy<<endl;
        //LocalS<<site<<"  "<<sz<<"  "<<sx<<"   "<<sy<<endl;
        LocalS<<site<<"  "<<site_x<<"   "<< site_y <<"   "<<rx_<<"   "<<ry_<<"   "<<sz<<"  "<<sx<<"   "<<sy<<endl;

    }

}
void Observables::Calculate_Order_Params(){


    OParams.value.resize(MFParams_.OParams_.value.size());
    OParams.rows.resize(MFParams_.OParams_.rows.size());
    OParams.columns.resize(MFParams_.OParams_.columns.size());

    int row_, col_;
    complex<double> val_;
    for(int ind=0;ind<OParams.value.size();ind++){
        row_=MFParams_.OParams_.rows[ind];
        col_=MFParams_.OParams_.columns[ind];
        OParams.rows[ind]=row_;
        OParams.columns[ind]=col_;

        val_=zero_complex;
        for(int n=0;n<Hamiltonian_.eigs_.size();n++){
            val_+= (
                        (( conj(Hamiltonian_.Ham_(row_,n))*Hamiltonian_.Ham_(col_,n)))*
                        (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0))
                        );

        }

        OParams.value[ind]=val_;

    }

    //checking Hermiticity of order parameters


}


void Observables::Calculate_Nw()
{

    //---------Read from input file-----------------------//
    string fileout = "Nw_total.txt";
    double omega_min, omega_max, d_omega;
    double eta = Parameters_.eta_dos;
    omega_min = Hamiltonian_.eigs_[0] - 5.0;
    omega_max = Hamiltonian_.eigs_[Hamiltonian_.eigs_.size() - 1] + 5.0;
    d_omega = Parameters_.dw_dos;
    //---------------------------------------------------//

    int omega_index_max = int((omega_max - omega_min) / (d_omega));

    ofstream file_Nw_out(fileout.c_str());

    double temp_val;

    for (int omega_ind = 0; omega_ind < omega_index_max; omega_ind++)
    {

        temp_val = 0.0;

        for (int n = 0; n < Hamiltonian_.Ham_.n_row(); n++)
        {

            temp_val +=
                    Lorentzian(omega_min + (omega_ind * d_omega) -Hamiltonian_.eigs_[n], eta);
        }

        file_Nw_out << omega_min + (omega_ind * d_omega) << "     "
                    << (1.0 / Hamiltonian_.eigs_.size()) * temp_val
                    << "     " << endl;
    }

    file_Nw_out << "#mu = " << Parameters_.mus << endl;

    ofstream spectrum;
    spectrum.open("spectrum.txt");

    for (int i=0; i < Hamiltonian_.eigs_.size(); i++){
        spectrum << i << "\t" << Hamiltonian_.eigs_[i] << endl;
    }
    spectrum.close();
}



void Observables::Calculate_Single_Particle_Density_Matrix(){

    /*
      NOTE:
      SP_Density_Matrix[alpha][beta] = <c_{alpha^{daggger}} c_{beta}>
     */
    SP_Density_Matrix.resize(ns_*2);
    for(int i=0;i<ns_*2;i++){
        SP_Density_Matrix[i].resize(ns_*2);
    }

    for(int alpha_=0;alpha_<ns_*2;alpha_++){
        for(int beta_=0;beta_<ns_*2;beta_++){
            SP_Density_Matrix[alpha_][beta_] = zero_complex;
            for(int n=0;n<ns_*2;n++){
                SP_Density_Matrix[alpha_][beta_] += conj(Hamiltonian_.Ham_(alpha_,n))*Hamiltonian_.Ham_(beta_,n)*
                        (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0));
            }
        }
    }

}


complex<double> Observables::Two_particle_Den_Mat(int _alpha, int _beta, int _gamma, int _delta){

    complex<double> temp;
    complex<double> delta_gamma_beta;

    if(_gamma == _beta){
        delta_gamma_beta=one_complex;
    }
    else{
        assert(_gamma != _beta);
        delta_gamma_beta=zero_complex;
    }

    temp = (SP_Density_Matrix[_alpha][_beta]*SP_Density_Matrix[_gamma][_delta])
            +
            (SP_Density_Matrix[_alpha][_delta]*(delta_gamma_beta - SP_Density_Matrix[_gamma][_beta]));
    return temp;
}


double Observables::Lorentzian(double x, double brd){
    double temp;

    temp = (1.0/PI)*( (brd/2.0)/ ( (x*x) + ((brd*brd)/4.0) ) );

    return temp;

}


void Observables::Total_Energy_Average(double Curr_QuantE, double CurrE){

    AVG_Total_Energy += Curr_QuantE + CurrE;
    AVG_Total_Energy_sqr += (Curr_QuantE + CurrE)*(Curr_QuantE + CurrE);
}

void Observables::Initialize(){

} // ----------
void Observables::Calculate_SpinSpincorrelations_Smartly(){


    string S2_out = "Local_S2.txt";
    ofstream file_S2_out(S2_out.c_str());
    file_S2_out<<"#site_i  S^2[site_i]"<<endl;


    string SSr_out = "SSr.txt";
    ofstream file_SSr_out(SSr_out.c_str());
    file_SSr_out<<"#site_i  site_j  SS[site_i][site_j]"<<endl;

    int ci,cj, spin_i,spin_j;

    Mat_4_Complex_doub UP_UP_Fermi, DOWN_DOWN_Fermi, UP_DOWN_Fermi, DOWN_UP_Fermi ;
    Mat_4_Complex_doub UP_UP_1mFermi, DOWN_DOWN_1mFermi, UP_DOWN_1mFermi, DOWN_UP_1mFermi ;

    UP_UP_Fermi.resize(ns_);DOWN_DOWN_Fermi.resize(ns_);UP_DOWN_Fermi.resize(ns_);DOWN_UP_Fermi.resize(ns_);
    UP_UP_1mFermi.resize(ns_);DOWN_DOWN_1mFermi.resize(ns_);UP_DOWN_1mFermi.resize(ns_);DOWN_UP_1mFermi.resize(ns_);

    for(int n=0;n<ns_;n++){
        UP_UP_Fermi[n].resize(ns_);DOWN_DOWN_Fermi[n].resize(ns_);
        UP_DOWN_Fermi[n].resize(ns_);DOWN_UP_Fermi[n].resize(ns_);
        UP_UP_1mFermi[n].resize(ns_);DOWN_DOWN_1mFermi[n].resize(ns_);
        UP_DOWN_1mFermi[n].resize(ns_);DOWN_UP_1mFermi[n].resize(ns_);

        for(int m=0;m<ns_;m++){
            UP_UP_Fermi[n][m].resize(1);DOWN_DOWN_Fermi[n][m].resize(1);
            UP_DOWN_Fermi[n][m].resize(1);DOWN_UP_Fermi[n][m].resize(1);
            UP_UP_1mFermi[n][m].resize(1);DOWN_DOWN_1mFermi[n][m].resize(1);
            UP_DOWN_1mFermi[n][m].resize(1);DOWN_UP_1mFermi[n][m].resize(1);
            for(int orb=0;orb<1;orb++){
                UP_UP_Fermi[n][m][orb].resize(1);DOWN_DOWN_Fermi[n][m][orb].resize(1);
                UP_DOWN_Fermi[n][m][orb].resize(1);DOWN_UP_Fermi[n][m][orb].resize(1);
                UP_UP_1mFermi[n][m][orb].resize(1);DOWN_DOWN_1mFermi[n][m][orb].resize(1);
                UP_DOWN_1mFermi[n][m][orb].resize(1);DOWN_UP_1mFermi[n][m][orb].resize(1);
            }
        }
    }


    for(int i=0;i<ns_;i++){

        for(int j=0;j<ns_;j++){

            for(int orbi=0;orbi<1;orbi++){
                for(int orbj=0;orbj<1;orbj++){

                    UP_UP_Fermi[i][j][orbi][orbj] = zero_complex;
                    DOWN_DOWN_Fermi[i][j][orbi][orbj] = zero_complex;
                    UP_DOWN_Fermi[i][j][orbi][orbj] = zero_complex;
                    DOWN_UP_Fermi[i][j][orbi][orbj] = zero_complex;

                    UP_UP_1mFermi[i][j][orbi][orbj] = zero_complex;
                    DOWN_DOWN_1mFermi[i][j][orbi][orbj] = zero_complex;
                    UP_DOWN_1mFermi[i][j][orbi][orbj] = zero_complex;
                    DOWN_UP_1mFermi[i][j][orbi][orbj] = zero_complex;


                    for(int n=0;n<Hamiltonian_.eigs_.size();n++){

                        spin_i=0;spin_j=0;
                        ci=Coordinates_.Nc_dof(i,orbi + 1*spin_i);
                        cj=Coordinates_.Nc_dof(j,orbj + 1*spin_j);
                        UP_UP_Fermi[i][j][orbi][orbj] +=  conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0));

                        spin_i=1;spin_j=1;
                        ci=Coordinates_.Nc_dof(i,orbi + 1*spin_i);
                        cj=Coordinates_.Nc_dof(j,orbj + 1*spin_j);
                        DOWN_DOWN_Fermi[i][j][orbi][orbj] +=  conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0));

                        spin_i=0;spin_j=1;
                        ci=Coordinates_.Nc_dof(i,orbi + 1*spin_i);
                        cj=Coordinates_.Nc_dof(j,orbj + 1*spin_j);
                        UP_DOWN_Fermi[i][j][orbi][orbj] +=  conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0));

                        spin_i=1;spin_j=0;
                        ci=Coordinates_.Nc_dof(i,orbi + 1*spin_i);
                        cj=Coordinates_.Nc_dof(j,orbj + 1*spin_j);
                        DOWN_UP_Fermi[i][j][orbi][orbj] +=  conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0));

                        spin_i=0;spin_j=0;
                        ci=Coordinates_.Nc_dof(i,orbi + 1*spin_i);
                        cj=Coordinates_.Nc_dof(j,orbj + 1*spin_j);
                        UP_UP_1mFermi[i][j][orbi][orbj] +=  conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                ( 1.0 - (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)));

                        spin_i=1;spin_j=1;
                        ci=Coordinates_.Nc_dof(i,orbi + 1*spin_i);
                        cj=Coordinates_.Nc_dof(j,orbj + 1*spin_j);
                        DOWN_DOWN_1mFermi[i][j][orbi][orbj] +=  conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                ( 1.0 - (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)));

                        spin_i=0;spin_j=1;
                        ci=Coordinates_.Nc_dof(i,orbi + 1*spin_i);
                        cj=Coordinates_.Nc_dof(j,orbj + 1*spin_j);
                        UP_DOWN_1mFermi[i][j][orbi][orbj] +=  conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                ( 1.0 - (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)));

                        spin_i=1;spin_j=0;
                        ci=Coordinates_.Nc_dof(i,orbi + 1*spin_i);
                        cj=Coordinates_.Nc_dof(j,orbj + 1*spin_j);
                        DOWN_UP_1mFermi[i][j][orbi][orbj] +=  conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                ( 1.0 - (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)));

                    }
                }
            }
        }
    }

    Mat_4_Complex_doub SS_ri_rj;
    SS_ri_rj.resize(ns_);
    for(int site_i=0;site_i<ns_;site_i++){
        SS_ri_rj[site_i].resize(ns_);
        for(int site_j=0;site_j<ns_;site_j++){
            SS_ri_rj[site_i][site_j].resize(1);
            for(int orb=0;orb<1;orb++){
                SS_ri_rj[site_i][site_j][orb].resize(1);
            }
        }
    }


    for(int i=0;i<ns_;i++){
        for(int j=0;j<ns_;j++){
            for(int orbi=0;orbi<1;orbi++){
                for(int orbj=0;orbj<1;orbj++){

                    SS_ri_rj[i][j][orbi][orbj]=zero_complex;

                    //SzSz..
                    SS_ri_rj[i][j][orbi][orbj] +=(0.25*one_complex)*(
                                (UP_UP_Fermi[i][j][orbi][orbj]*UP_UP_1mFermi[j][i][orbj][orbi])
                                - (UP_DOWN_Fermi[i][j][orbi][orbj]*DOWN_UP_1mFermi[j][i][orbj][orbi])
                                + (DOWN_DOWN_Fermi[i][j][orbi][orbj]*DOWN_DOWN_1mFermi[j][i][orbj][orbi])
                                - (DOWN_UP_Fermi[i][j][orbi][orbj]*UP_DOWN_1mFermi[j][i][orbj][orbi])
                                + (UP_UP_Fermi[i][i][orbi][orbi]*UP_UP_Fermi[j][j][orbj][orbj])
                                - (UP_UP_Fermi[i][i][orbi][orbi]*DOWN_DOWN_Fermi[j][j][orbj][orbj])
                                + (DOWN_DOWN_Fermi[i][i][orbi][orbi]*DOWN_DOWN_Fermi[j][j][orbj][orbj])
                                - (DOWN_DOWN_Fermi[i][i][orbi][orbi]*UP_UP_Fermi[j][j][orbj][orbj])
                                );


                    //0.5*(S+S- + S-S+)
                    SS_ri_rj[i][j][orbi][orbj] +=(0.5*one_complex)*(

                                (UP_UP_Fermi[i][j][orbi][orbj]*DOWN_DOWN_1mFermi[j][i][orbj][orbi])
                                + (DOWN_DOWN_Fermi[i][j][orbi][orbj]*UP_UP_1mFermi[j][i][orbj][orbi])
                                + (UP_DOWN_Fermi[i][i][orbi][orbi]*DOWN_UP_Fermi[j][j][orbj][orbj])
                                + (DOWN_UP_Fermi[i][i][orbi][orbi]*UP_DOWN_Fermi[j][j][orbj][orbj])

                                );


                    //cout<<i<<"\t"<<j<<" done"<<endl;
                }
            }
        }

    }

    Mat_2_Complex_doub SS_ij;
    SS_ij.resize(ns_);
    for(int site_i=0;site_i<ns_;site_i++){
        SS_ij[site_i].resize(ns_);
    }

    for(int site_i=0;site_i<ns_;site_i++){
        for(int site_j=0;site_j<ns_;site_j++){
            SS_ij[site_i][site_j]=zero_complex;
            for(int orbi=0;orbi<1;orbi++){
                for(int orbj=0;orbj<1;orbj++){
                    SS_ij[site_i][site_j] += SS_ri_rj[site_i][site_j][orbi][orbj];
                }
            }
        }
    }




    for(int site_i=0;site_i<ns_;site_i++){
        for(int site_j=0;site_j<ns_;site_j++){
            file_SSr_out<<site_i<<"\t"<<site_j<<"\t"<<
                          real(SS_ij[site_i][site_j])<<"\t"<<
                          imag(SS_ij[site_i][site_j])<<endl;
        }
    }



    complex<double> Avg_S2=0.0;

    for(int site_i=0;site_i<ns_;site_i++){
        file_S2_out<<site_i<<"\t"<<real(SS_ij[site_i][site_i])
                  <<"\t"<<imag(SS_ij[site_i][site_i])<<endl;

        Avg_S2 += one_complex*(SS_ij[site_i][site_i]);

    }

    cout<<"Avg Local Moment (S^2) = "<<real(Avg_S2)/(1.0*ns_)<<"\t"<<imag(Avg_S2)/(1.0*ns_)<<endl;



    //Assuming square lattice
    int LY_=3; //This are what geometry looks like, flipped from model_TL_..inp file
    int LX_=3;
    int mLX_, mLY_;
    mLX_=-1*LX_;mLY_=-1*LY_;

    int ix,iy,ixp,iyp,iy_new, iyp_new,site_, site_p;
    string SSq_out = "SSq_Lx"   + to_string(LX_) + "_Ly" + to_string(LY_) + "_OUR_XC.txt";
    ofstream file_SSq_out(SSq_out.c_str());
    file_SSq_out<<"#These results are not true in general, use cautiously"<<endl;
    file_SSq_out<<"#qx qy  SS_q"<<endl;

    double x_dis, y_dis, xp_dis, yp_dis;

    complex<double> temp_doub;
    for(int qx_i=mLX_*5;qx_i<LX_*5;qx_i++){
        for(int qy_i=mLY_*5;qy_i<LY_*5;qy_i++){
            temp_doub=0.0;


            //------------------------------------------
            for(int ix=0;ix<LX_;ix++){

                for(int iy=0;iy<LY_;iy++){

//                    x_dis = (1.0*ix) ;
//                    y_dis = (1.0)*iy;

                    x_dis = (1.0*ix) + (0.5*iy);
                    y_dis = (sqrt(3.0)/2.0)*iy;

                    if(ix%2==0){
                        iy_new=iy;
                    }
                    else{
                        iy_new = (LY_-1 - iy);
                    }
                    site_ = iy_new + ix*(LY_);


                    for(int ixp=0;ixp<LX_;ixp++){
                        for(int iyp=0;iyp<LY_;iyp++){

//                            xp_dis = (1.0*ixp);
//                            yp_dis = (1.0)*iyp;

                            xp_dis = (1.0*ixp) + (0.5*iyp);
                            yp_dis = (sqrt(3.0)/2.0)*iyp;

                            if(ixp%2==0){
                                iyp_new=iyp;
                            }
                            else{
                                iyp_new = (LY_-1 - iyp);
                            }
                            site_p = iyp_new + ixp*(LY_);


                            temp_doub += exp(iota_complex*( ((yp_dis-y_dis)*qy_i*(PI/(1.0*LY_)))  + ((xp_dis-x_dis)*qx_i*(PI/(1.0*LX_))) ) )*
                                         (SS_ij[site_][site_p]);


                        }
                    }


                }
            }

            //--------------------------------------------


            file_SSq_out<<qx_i*(PI/(1.0*LX_))<<"    "<<qy_i*(PI/(1.0*LY_))<<"   "<<temp_doub.real()<<"   "<<temp_doub.imag()<<endl;




        }
        file_SSq_out<<endl;
    }



}


void Observables::Calculate_DenDencorrelations_Smartly(){

    string SSr_out = "NNr.txt";
    ofstream file_SSr_out(SSr_out.c_str());
    file_SSr_out<<"#site_i    site_j   NN[site_i][site_j]"<<endl;

    int ci,cj, spin_i,spin_j;

    Mat_4_Complex_doub UP_UP_Fermi, DOWN_DOWN_Fermi, UP_DOWN_Fermi, DOWN_UP_Fermi ;
    Mat_4_Complex_doub UP_UP_1mFermi, DOWN_DOWN_1mFermi, UP_DOWN_1mFermi, DOWN_UP_1mFermi ;

    UP_UP_Fermi.resize(ns_);DOWN_DOWN_Fermi.resize(ns_);UP_DOWN_Fermi.resize(ns_);DOWN_UP_Fermi.resize(ns_);
    UP_UP_1mFermi.resize(ns_);DOWN_DOWN_1mFermi.resize(ns_);UP_DOWN_1mFermi.resize(ns_);DOWN_UP_1mFermi.resize(ns_);

    for(int n=0;n<ns_;n++){
        UP_UP_Fermi[n].resize(ns_);DOWN_DOWN_Fermi[n].resize(ns_);
        UP_DOWN_Fermi[n].resize(ns_);DOWN_UP_Fermi[n].resize(ns_);
        UP_UP_1mFermi[n].resize(ns_);DOWN_DOWN_1mFermi[n].resize(ns_);
        UP_DOWN_1mFermi[n].resize(ns_);DOWN_UP_1mFermi[n].resize(ns_);

        for(int m=0;m<ns_;m++){
            UP_UP_Fermi[n][m].resize(1);DOWN_DOWN_Fermi[n][m].resize(1);
            UP_DOWN_Fermi[n][m].resize(1);DOWN_UP_Fermi[n][m].resize(1);
            UP_UP_1mFermi[n][m].resize(1);DOWN_DOWN_1mFermi[n][m].resize(1);
            UP_DOWN_1mFermi[n][m].resize(1);DOWN_UP_1mFermi[n][m].resize(1);
            for(int orb=0;orb<1;orb++){
                UP_UP_Fermi[n][m][orb].resize(1);DOWN_DOWN_Fermi[n][m][orb].resize(1);
                UP_DOWN_Fermi[n][m][orb].resize(1);DOWN_UP_Fermi[n][m][orb].resize(1);
                UP_UP_1mFermi[n][m][orb].resize(1);DOWN_DOWN_1mFermi[n][m][orb].resize(1);
                UP_DOWN_1mFermi[n][m][orb].resize(1);DOWN_UP_1mFermi[n][m][orb].resize(1);
            }
        }
    }


    for(int i=0;i<ns_;i++){
        for(int j=0;j<ns_;j++){
            for(int orbi=0;orbi<1;orbi++){
                for(int orbj=0;orbj<1;orbj++){

                    UP_UP_Fermi[i][j][orbi][orbj] = zero_complex;
                    DOWN_DOWN_Fermi[i][j][orbi][orbj] = zero_complex;
                    UP_DOWN_Fermi[i][j][orbi][orbj] = zero_complex;
                    DOWN_UP_Fermi[i][j][orbi][orbj] = zero_complex;

                    UP_UP_1mFermi[i][j][orbi][orbj] = zero_complex;
                    DOWN_DOWN_1mFermi[i][j][orbi][orbj] = zero_complex;
                    UP_DOWN_1mFermi[i][j][orbi][orbj] = zero_complex;
                    DOWN_UP_1mFermi[i][j][orbi][orbj] = zero_complex;


                    for(int n=0;n<Hamiltonian_.eigs_.size();n++){

                        spin_i=0;spin_j=0;
                        ci=Coordinates_.Nc_dof(i,orbi + 1*spin_i);
                        cj=Coordinates_.Nc_dof(j,orbj + 1*spin_j);
                        UP_UP_Fermi[i][j][orbi][orbj] +=  conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0));

                        spin_i=1;spin_j=1;
                        ci=Coordinates_.Nc_dof(i,orbi + 1*spin_i);
                        cj=Coordinates_.Nc_dof(j,orbj + 1*spin_j);
                        DOWN_DOWN_Fermi[i][j][orbi][orbj] +=  conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0));

                        spin_i=0;spin_j=1;
                        ci=Coordinates_.Nc_dof(i,orbi + 1*spin_i);
                        cj=Coordinates_.Nc_dof(j,orbj + 1*spin_j);
                        UP_DOWN_Fermi[i][j][orbi][orbj] +=  conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0));

                        spin_i=1;spin_j=0;
                        ci=Coordinates_.Nc_dof(i,orbi + 1*spin_i);
                        cj=Coordinates_.Nc_dof(j,orbj + 1*spin_j);
                        DOWN_UP_Fermi[i][j][orbi][orbj] +=  conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0));

                        spin_i=0;spin_j=0;
                        ci=Coordinates_.Nc_dof(i,orbi + 1*spin_i);
                        cj=Coordinates_.Nc_dof(j,orbj + 1*spin_j);
                        UP_UP_1mFermi[i][j][orbi][orbj] +=  conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                ( 1.0 - (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)));

                        spin_i=1;spin_j=1;
                        ci=Coordinates_.Nc_dof(i,orbi + 1*spin_i);
                        cj=Coordinates_.Nc_dof(j,orbj + 1*spin_j);
                        DOWN_DOWN_1mFermi[i][j][orbi][orbj] +=  conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                ( 1.0 - (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)));

                        spin_i=0;spin_j=1;
                        ci=Coordinates_.Nc_dof(i,orbi + 1*spin_i);
                        cj=Coordinates_.Nc_dof(j,orbj + 1*spin_j);
                        UP_DOWN_1mFermi[i][j][orbi][orbj] +=  conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                ( 1.0 - (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)));

                        spin_i=1;spin_j=0;
                        ci=Coordinates_.Nc_dof(i,orbi + 1*spin_i);
                        cj=Coordinates_.Nc_dof(j,orbj + 1*spin_j);
                        DOWN_UP_1mFermi[i][j][orbi][orbj] +=  conj(Hamiltonian_.Ham_(ci,n))*Hamiltonian_.Ham_(cj,n)*
                                ( 1.0 - (1.0/( exp((Hamiltonian_.eigs_[n]-Parameters_.mus)*Parameters_.beta ) + 1.0)));

                    }
                }
            }
        }
    }

    Mat_4_Complex_doub SS_ri_rj;
    SS_ri_rj.resize(ns_);
    for(int site_i=0;site_i<ns_;site_i++){
        SS_ri_rj[site_i].resize(ns_);
        for(int site_j=0;site_j<ns_;site_j++){
            SS_ri_rj[site_i][site_j].resize(1);
            for(int orb=0;orb<1;orb++){
                SS_ri_rj[site_i][site_j][orb].resize(1);
            }
        }
    }


    for(int i=0;i<ns_;i++){
        for(int j=0;j<ns_;j++){
            for(int orbi=0;orbi<1;orbi++){
                for(int orbj=0;orbj<1;orbj++){

                    SS_ri_rj[i][j][orbi][orbj]=zero_complex;

                    //denden..
                    SS_ri_rj[i][j][orbi][orbj] +=(1.0*one_complex)*(
                                (UP_UP_Fermi[i][j][orbi][orbj]*UP_UP_1mFermi[j][i][orbj][orbi])
                                + (UP_DOWN_Fermi[i][j][orbi][orbj]*DOWN_UP_1mFermi[j][i][orbj][orbi])
                                + (DOWN_DOWN_Fermi[i][j][orbi][orbj]*DOWN_DOWN_1mFermi[j][i][orbj][orbi])
                                + (DOWN_UP_Fermi[i][j][orbi][orbj]*UP_DOWN_1mFermi[j][i][orbj][orbi])
                                + (UP_UP_Fermi[i][i][orbi][orbi]*UP_UP_Fermi[j][j][orbj][orbj])
                                + (UP_UP_Fermi[i][i][orbi][orbi]*DOWN_DOWN_Fermi[j][j][orbj][orbj])
                                + (DOWN_DOWN_Fermi[i][i][orbi][orbi]*DOWN_DOWN_Fermi[j][j][orbj][orbj])
                                + (DOWN_DOWN_Fermi[i][i][orbi][orbi]*UP_UP_Fermi[j][j][orbj][orbj])
                                );
                    //-<den_i><den_j>
                    /*den_i_up=Local_n_orb_resolved[Coordinates_.Nc_dof(i,0)];
                    den_i_dn=Local_n_orb_resolved[Coordinates_.Nc_dof(i,1)];
                    den_j_up=Local_n_orb_resolved[Coordinates_.Nc_dof(j,0)];
                    den_j_dn=Local_n_orb_resolved[Coordinates_.Nc_dof(j,1)];
                    SS_ri_rj[i][j][orbi][orbj] +=(-1.0*one_complex)*(den_i_up+den_i_dn)*
                            (den_j_up+den_j_dn)*/;

                    SS_ri_rj[i][j][orbi][orbj] +=(-1.0*one_complex)*(Avg_local_Nup+Avg_local_Ndn)*
                            (Avg_local_Nup+Avg_local_Ndn);


                    //cout<<i<<"\t"<<j<<" done"<<endl;
                }
            }
        }

    }



    Mat_2_Complex_doub SS_ij;
    SS_ij.resize(ns_);
    for(int site_i=0;site_i<ns_;site_i++){
        SS_ij[site_i].resize(ns_);
    }

    for(int site_i=0;site_i<ns_;site_i++){
        for(int site_j=0;site_j<ns_;site_j++){
            SS_ij[site_i][site_j]=zero_complex;
            for(int orbi=0;orbi<1;orbi++){
                for(int orbj=0;orbj<1;orbj++){
                    SS_ij[site_i][site_j] += SS_ri_rj[site_i][site_j][orbi][orbj];
                }
            }
        }
    }



    for(int site_i=0;site_i<ns_;site_i++){
        for(int site_j=0;site_j<ns_;site_j++){
            file_SSr_out<<site_i<<"\t"<<site_j<<"\t"<<
                          real(SS_ij[site_i][site_j])<<"\t"<<
                          imag(SS_ij[site_i][site_j])<<endl;

        }
    }



}




#endif // OBSERVABLES_H
