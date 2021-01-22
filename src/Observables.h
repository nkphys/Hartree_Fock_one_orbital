#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#include "Hamiltonian.h"
#include "tensor_type.h"
#include "functions.h"

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
    void Calculate_Order_Params();
    void Get_OrderParameters_diffs();
    void Update_OrderParameters(int iter);
    void Hartree_Filter();
    void Update_OrderParameters_Second_Broyden(int iter_count);
    void Invert_Beta();
    void Invert_Beta_double();
    void Calculate_Single_Particle_Density_Matrix();
    void Calculate_SpinSpincorrelations_Smartly();
    void Calculate_DenDencorrelations_Smartly();

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


    //Declare Broyden_Mixing vectors
    vector<double> F_n; //F_n=x_n_out - x_n_in []
    vector<double> F_nm1;
    vector<double> DeltaF_n; //DeltaF_n=F_n - F-nm1;
    vector<double> Delta_x_n; //Delta_x_n= x_n_in - x_nm1_in;
    Mat_2_doub Jinv_n;
    Mat_2_doub Jinv_np1;


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




};
/*
 * ***********
 *  Functions in Class Observables ------
 *  ***********
*/



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

    assert(NewInd_to_OldInd.size()==MFParams_.OParams_.value.size() + (MFParams_.OParams_.value.size() - 2*ns_));

    Offset_=MFParams_.OParams_.value.size();
    int size_vecs;
    size_vecs=NewInd_to_OldInd.size();
    vector<double> vec_V, vec_U;
    double Denominator_;
    vector<double> vec_L;



    if(Parameters_.Broyden_Mixing==true){
        //assert(false);

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
                Denominator_ +=Delta_x_n[i]*vec_V[i];
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
                    Mat_Temp[i][j] = (vec_L[i]*vec_U[j])/(Denominator_);
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



}


void Observables::Calculate_DenDencorrelations_Smartly(){



    string SSr_out = "NNr.txt";
    ofstream file_SSr_out(SSr_out.c_str());
    file_SSr_out<<"#site_i   site_i(x)    site_i(y)    site_j   site_j(x)    site_j(y)  NN[site_i][site_j]"<<endl;

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

                   // SS_ri_rj[i][j][orbi][orbj] +=(-1.0*one_complex)*(Avg_local_Nup+Avg_local_Ndn)*
                    //        (Avg_local_Nup+Avg_local_Ndn);


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
