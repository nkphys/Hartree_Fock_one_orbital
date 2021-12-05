#include <algorithm>
#include <functional>
#include <math.h>
#include "tensor_type.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#define PI acos(-1.0)

#ifndef Hamiltonian_class
#define Hamiltonian_class

extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
                         std::complex<double> *,int *, double *, int *);
//zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);

class Hamiltonian {
public:

    Hamiltonian(Parameters& Parameters__, Coordinates&  Coordinates__, MFParams& MFParams__ )
        :Parameters_(Parameters__),Coordinates_(Coordinates__),MFParams_(MFParams__)

    {
        Initialize();
        Hoppings();
        HTBCreate();
    }


    void Initialize();    //::DONE
    void Hoppings();        //::DONE
    double GetCLEnergy();    //::DONE
    void InteractionsCreate();   //::DONE
    void Check_Hermiticity();  //::DONE
    void HTBCreate();   //::DONE
    double chemicalpotential(double muin,double Particles);    //::DONE

    double TotalDensity();   //::DONE
    double E_QM();   //::DONE

    void Diagonalize(char option);   //::DONE
    void copy_eigs(int i);  //::DONE

    int convert_jm_to_int(string jm_val);

    Parameters &Parameters_;
    Coordinates &Coordinates_;
    MFParams &MFParams_;
    int ns_;
    Matrix<complex<double>> HTB_;
    Matrix<complex<double>> Ham_;
    Matrix<double> Tx,Ty,Tpxpy,Tpxmy;
    vector<double> eigs_,eigs_saved_,sx_,sy_,sz_;

};



double Hamiltonian::chemicalpotential(double muin,double Particles){


    double mu_out;

    if(!Parameters_.FixingMu)
    {
        double n1,N;
        double dMubydN;
        int nstate = eigs_.size();
        dMubydN = 0.0005*(eigs_[nstate-1] - eigs_[0])/nstate;
        N=Particles;
        //temp=Parameters_.temp;
        mu_out = muin;
        bool converged=false;
        int final_i;


        if(1==2){
            for(int i=0;i<50000;i++){
                n1=0.0;
                for(int j=0;j<nstate;j++){
                    n1+=double(1.0/( exp( (eigs_[j]-mu_out)*Parameters_.beta ) + 1.0));
                }
                //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
                if(abs(N-n1)<double(0.000001)){
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
            mu1=eigs_[0]- (5.0/Parameters_.beta);
            mu2=eigs_[nstate-1] + (5.0/Parameters_.beta);
            for(int i=0;i<40000;i++){
                n1=0.0;
                for(int j=0;j<nstate;j++){
                    n1+=double(1.0/( exp( (eigs_[j]-mu_temp)*Parameters_.beta ) + 1.0));
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
                //cout<<"mu_temp = "<<mu_temp<<"   "<<mu1<<"    "<<mu2<<"   "<<eigs_[nstate-1]<<"  "<<eigs_[0]<<"  "<<n1<<endl;
            }

            if(!converged){
                cout<<"mu_not_converged, N = "<<n1<<endl;
            }
            else{
                //cout<<"mu converged, N = "<<n1<<endl;
            }

            mu_out = mu_temp;
        }

    }

    else{
        mu_out = Parameters_.MuValueFixed;
    }

    return mu_out;
} // ----------


void Hamiltonian::Initialize(){

    ns_=Parameters_.ns;

    int space=Coordinates_.no_dof_ ;

    HTB_.resize(space,space);
    Ham_.resize(space,space);
    eigs_.resize(space);
    eigs_saved_.resize(space);

} // ----------

double Hamiltonian::TotalDensity(){

    double n1=0.0;
    /*
    for(int j=0;j<eigs_.size();j++){
        n1 +=  1.0f/( exp(Parameters_.beta*(eigs_[j]-Parameters_.mus) ) + 1.0);
    }
    */
    return n1;

} // ----------



double Hamiltonian::E_QM(){

    double E=0.0;
    for(int j=0;j<eigs_.size();j++){
        //E +=  (eigs_[j]-Parameters_.mus)/( exp(Parameters_.beta*(eigs_[j]-Parameters_.mus) ) + 1.0);
        E +=  (eigs_[j])/( exp(Parameters_.beta*(eigs_[j]-Parameters_.mus) ) + 1.0);
    }
    return E;

} // ----------



double Hamiltonian::GetCLEnergy(){

    complex<double> EClassical=zero_complex;
    double EClassical_doub;

    int UP_, DOWN_;
    UP_=0;DOWN_=1;

    int spin_row, spin_col, site_row, site_col; //For Hamiltonian

    int alpha, beta;
    int site_alpha, site_beta;
    int spin_alpha, spin_beta;
    int alpha_new, beta_new, ind_new, SI_new;
    for(int ind=0;ind<MFParams_.OParams_.value.size();ind++){

        alpha=MFParams_.OParams_.rows[ind];
        beta=MFParams_.OParams_.columns[ind];

        site_alpha = alpha%ns_;
        site_beta = beta%ns_;
        spin_alpha =  (alpha - site_alpha)/ns_;
        spin_beta =  (beta - site_beta)/ns_;

        //<c^{dagger}_{alpha}c_{beta}>; beta>=alpha

        if(alpha==beta){

            //Onsite Hartree
            //spin_row = spin_col = BAR(spin_alpha)
            spin_row = abs(1-spin_alpha);
            spin_col = spin_row;
            site_row = site_alpha;
            site_col = site_alpha;
            alpha_new = site_row + spin_row*ns_;
            beta_new = site_col + spin_col*ns_;
            SI_new = alpha_new + (2*ns_*beta_new);
            ind_new = MFParams_.SI_to_ind[SI_new];

            EClassical += (-0.5)*Parameters_.U_onsite*MFParams_.OParams_.value[ind]*
                    MFParams_.OParams_.value[ind_new];

            //Offsite Hartree
            for(int site_j=0;site_j<ns_;site_j++){
                if(abs(Parameters_.LongRangeInteractions[site_alpha][site_j])>0.000001){
                    for(int spin_j=0;spin_j<2;spin_j++){
                        spin_row=spin_j;
                        spin_col=spin_j;
                        site_row=site_j;
                        site_col=site_j;
                        alpha_new = site_row + spin_row*ns_;
                        beta_new = site_col + spin_col*ns_;
                        SI_new = alpha_new + (2*ns_*beta_new);
                        ind_new = MFParams_.SI_to_ind[SI_new];

                        EClassical += (-0.5)*Parameters_.LongRangeInteractions[site_alpha][site_j]*
                                MFParams_.OParams_.value[ind]*MFParams_.OParams_.value[ind_new];

                    }
                }
            }

        }
        else{ //Fock Terms

            if(!Parameters_.Just_Hartree){
            //Onsite Fock
            if(site_alpha==site_beta){
                assert(spin_alpha==UP_);
                assert(spin_beta==DOWN_);
                site_row=site_alpha;
                site_col=site_alpha;
                spin_col=DOWN_;
                spin_row=UP_;
                alpha_new = site_row + spin_row*ns_;
                beta_new = site_col + spin_col*ns_;
                SI_new = alpha_new + (2*ns_*beta_new);
                ind_new = MFParams_.SI_to_ind[SI_new];

                EClassical += 1.0*Parameters_.U_onsite *
                        MFParams_.OParams_.value[ind]*conj(MFParams_.OParams_.value[ind_new]);

            }
            else{ //OFFSITE Fock Terms


                //up-up or dn-dn
                if(spin_alpha==spin_beta){
                    assert(site_beta > site_alpha);
                    if(abs(Parameters_.LongRangeInteractions[site_alpha][site_beta])>0.000001){

                        site_row = site_alpha;
                        site_col = site_beta;
                        spin_row = spin_alpha;
                        spin_col = spin_alpha;
                        alpha_new = site_row + spin_row*ns_;
                        beta_new = site_col + spin_col*ns_;
                        SI_new = alpha_new + (2*ns_*beta_new);
                        ind_new = MFParams_.SI_to_ind[SI_new];

                        EClassical +=  (1.0)*Parameters_.LongRangeInteractions[site_alpha][site_beta]*
                                MFParams_.OParams_.value[ind]*conj(MFParams_.OParams_.value[ind_new]);
//                        EClassical +=  conj((1.0*0.5)*Parameters_.LongRangeInteractions[site_alpha][site_beta]*
//                                MFParams_.OParams_.value[ind]*MFParams_.OParams_.value[ind_new]);


                    }

                }
                else{
                    assert(spin_beta==DOWN_);
                    assert(spin_alpha==UP_);
                    if(abs(Parameters_.LongRangeInteractions[site_alpha][site_beta])>0.000001){
                        site_row = site_alpha;
                        site_col = site_beta;
                        spin_row = spin_alpha;
                        spin_col = spin_beta;
                        alpha_new = site_row + spin_row*ns_;
                        beta_new = site_col + spin_col*ns_;
                        SI_new = alpha_new + (2*ns_*beta_new);
                        ind_new = MFParams_.SI_to_ind[SI_new];

                        EClassical +=  (1.0)*Parameters_.LongRangeInteractions[site_alpha][site_beta]*
                                MFParams_.OParams_.value[ind]*conj(MFParams_.OParams_.value[ind_new]);

                    }

                }

            }
        }
        }

    }



    if(abs(EClassical.imag())>=0.1){
        cout<<"EClassical.imag() = "<< EClassical.imag()<<endl;
        assert(abs(EClassical.imag())<0.0000001);
    }
    EClassical_doub = EClassical.real();



    return EClassical_doub;

} // ----------



void Hamiltonian::InteractionsCreate(){

    int space=Coordinates_.no_dof_;

    int UP_, DOWN_;
    UP_=0;DOWN_=1;

    int spin_row, spin_col, site_row, site_col; //For Hamiltonian
    for(int i=0;i<space;i++) {
        for(int j=0;j<space;j++) {
            Ham_(i,j)=HTB_(i,j);
        }
    }

    int alpha, beta;
    int site_alpha, site_beta;
    int spin_alpha, spin_beta;
    for(int ind=0;ind<MFParams_.OParams_.value.size();ind++){

        alpha=MFParams_.OParams_.rows[ind];
        beta=MFParams_.OParams_.columns[ind];
        site_alpha = alpha%ns_;
        site_beta = beta%ns_;
        spin_alpha =  (alpha - site_alpha)/ns_;
        spin_beta =  (beta - site_beta)/ns_;
        //<c^{dagger}_{alpha}c_{beta}>; beta>=alpha

        if(alpha==beta){

            //Onsite Hartree
            //spin_row = spin_col = BAR(spin_alpha)
            spin_row = abs(1-spin_alpha);
            spin_col = spin_row;
            site_row = site_alpha;
            site_col = site_alpha;
            Ham_(site_row + spin_row*ns_, site_col + spin_col*ns_) += Parameters_.U_onsite *
                    MFParams_.OParams_.value[ind];

            //Offsite Hartree
            for(int site_j=0;site_j<ns_;site_j++){
                if(abs(Parameters_.LongRangeInteractions[site_alpha][site_j])>0.0000001){
                    for(int spin_j=0;spin_j<2;spin_j++){
                        spin_row=spin_j;
                        spin_col=spin_j;
                        site_row=site_j;
                        site_col=site_j;
                        Ham_(site_row + spin_row*ns_, site_col + spin_col*ns_) += Parameters_.LongRangeInteractions[site_alpha][site_j]*
                                MFParams_.OParams_.value[ind];

                    }
                }
            }
        }
        else{ //Fock Terms

            if(!Parameters_.Just_Hartree){
            //Onsite Fock
            if(site_alpha==site_beta){
                assert(spin_alpha==UP_);
                assert(spin_beta==DOWN_);
                site_row=site_alpha;
                site_col=site_alpha;
                spin_col=UP_;
                spin_row=DOWN_;
                Ham_(site_row + spin_row*ns_, site_col + spin_col*ns_) += -1.0*Parameters_.U_onsite *
                        MFParams_.OParams_.value[ind];
                //HERMITIAN CONJUGATE
                Ham_(site_col + spin_col*ns_, site_row + spin_row*ns_) += -1.0*Parameters_.U_onsite *
                        conj(MFParams_.OParams_.value[ind]);
            }
            else{ //OFFSITE Fock Terms

                //up-up or dn-dn
                if(spin_alpha==spin_beta){
                    assert(site_beta > site_alpha);
                    if(abs(Parameters_.LongRangeInteractions[site_alpha][site_beta])>0.0000001){

                        site_row = site_beta;
                        site_col = site_alpha;
                        spin_row = spin_alpha;
                        spin_col = spin_alpha;
                        Ham_(site_row + spin_row*ns_, site_col + spin_col*ns_) += (-1.0)*Parameters_.LongRangeInteractions[site_alpha][site_beta]*
                                MFParams_.OParams_.value[ind]; //0.5 is not needed here
                        Ham_(site_col + spin_col*ns_, site_row + spin_row*ns_) += (-1.0)*Parameters_.LongRangeInteractions[site_beta][site_alpha]*
                                conj(MFParams_.OParams_.value[ind]);


                    }

                }
                else{
                    assert(spin_beta==DOWN_);
                    assert(spin_alpha==UP_);
                    if(abs(Parameters_.LongRangeInteractions[site_alpha][site_beta])>0.0000001){
                        site_row = site_beta;
                        site_col = site_alpha;
                        spin_row = spin_beta;
                        spin_col = spin_alpha;
                        Ham_(site_row + spin_row*ns_, site_col + spin_col*ns_) += (-1.0)*Parameters_.LongRangeInteractions[site_alpha][site_beta]*
                                MFParams_.OParams_.value[ind];
                        Ham_(site_col + spin_col*ns_, site_row + spin_row*ns_) += (-1.0)*Parameters_.LongRangeInteractions[site_beta][site_alpha]*
                                conj(MFParams_.OParams_.value[ind]);

                    }

                }

            }
        }
        }

    }


} // ----------


int Hamiltonian::convert_jm_to_int(string jm_val){

    int val;
    if(jm_val=="3by2_m3by2"){val=0;}
    if(jm_val=="3by2_3by2"){val=1;}
    if(jm_val=="3by2_m1by2"){val=2;}
    if(jm_val=="3by2_1by2"){val=3;}
    if(jm_val=="1by2_m1by2"){val=4;}
    if(jm_val=="1by2_1by2"){val=5;}
    return val;
}

void Hamiltonian::Check_Hermiticity()

{
    complex<double> temp(0,0);
    complex<double>temp2;

    for(int i=0;i<Ham_.n_row();i++) {
        for(int j=0;j<Ham_.n_row();j++) {
            if(
                    abs(Ham_(i,j) - conj(Ham_(j,i)))>0.00001
                    ) {
                cout<<Ham_(i,j)<<endl;
                cout<<conj(Ham_(j,i))<<endl;

            }
            assert(
                        abs(Ham_(i,j) - conj(Ham_(j,i)))<0.00001
                        ); //+ Ham_(i+orbs_*ns_,j) + Ham_(i,j+orbs_*ns_);
            //temp +=temp2*conj(temp2);
        }
    }

    // cout<<"Hermiticity: "<<temp<<endl;
}





void Hamiltonian::Diagonalize(char option){

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


void Hamiltonian::HTBCreate(){

    int space=Coordinates_.no_dof_;

    int UP_, DOWN_;
    UP_=0;DOWN_=1;
    int ind_temp;

    for(int i=0;i<space;i++) {
        for(int j=0;j<space;j++) {
            HTB_(i,j)= Parameters_.LongRangeHoppings[i][j];
        }
    }

    for(int i=0;i<ns_;i++){
            for(int spin=0;spin<2;spin++){
                ind_temp = Coordinates_.Nc_dof(i,spin);
                HTB_(ind_temp,ind_temp) +=Parameters_.Onsite_E[i][spin];
            }
    }



} // ----------





void Hamiltonian::Hoppings(){

} // ----------

void Hamiltonian::copy_eigs(int i){

    int space=2*ns_;

    if (i == 0) {
        for(int j=0;j<space;j++) {
            eigs_[j] = eigs_saved_[j];
        }
    }
    else {
        for(int j=0;j<space;j++) {
            eigs_saved_[j] = eigs_[j];
        }
    }

}


#endif
