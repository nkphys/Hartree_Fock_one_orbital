#include "tensor_type.h"
#include "ParametersEngine.h"


#ifndef Coordinates_class
#define Coordinates_class
class Coordinates { 

public:
    Coordinates(int ns)
        :ns_(ns)
    {
        Numbering_lattice();
        Numbering_DOF();
    }

    /* Convention
 index = site + spin*Lx*Ly;
*/
    enum {PX=0,MX,PY,MY,PXPY,MXPY,MXMY,PXMY};	// not needed - but keep as a "key"

    void Numbering_lattice();
    int indx(int i);
    int indy(int i);

    void Numbering_DOF();
    int indx_dof(int i);
    int indy_dof(int i);
    int jm_state_dof(int i);
    int Nc_dof(int site, int dof);

    int NNc(int x, int y);
    int neigh(int site, int wneigh);
    int getneigh(int site,int wneigh);

    int ns_;
    Mat_1_int indx_,indy_;
    int no_dof_;
    Mat_1_int indx_dof_,indy_dof_,jm_state_dof_;
    Matrix<int> Nc_,neigh_;
    Matrix<int> Nc_dof_, neigh_dof_;
};

/*  
 * ***********
 *  Functions in Class Coordinates ------
 *  ***********
*/  

int Coordinates::indx(int i){
    if(i>ns_-1){perror("Coordinates.h:x-coordinate of lattice excede limit");}
    return indx_[i];
    // ----------
}


int Coordinates::indy(int i){
    if(i>ns_-1){perror("Coordinates.h:y-coordinate of lattice excede limit");}
    return indy_[i];
    // ----------
}


int Coordinates::indx_dof(int i){
    if(i>2*ns_-1){perror("Coordinates.h:x-coordinate of lattice excede limit");}
    return indx_dof_[i];
    // ----------
}


int Coordinates::indy_dof(int i){
    if(i>2*ns_-1){perror("Coordinates.h:y-coordinate of lattice excede limit");}
    return indy_dof_[i];
    // ----------
}


int Coordinates::Nc_dof(int site, int dof){
    if(site>ns_ || dof>1){
        cout<<"site ="<<site<<endl;
        cout<<"dof = "<<dof<<endl;
        perror("Coordinates.h:ith-sitelabel of lattice excede limit (site>ns_ || dof>1)");}
    return Nc_dof_(site,dof);
    // ----------
}


int Coordinates::neigh(int site, int wneigh){
    if(site>ns_-1 || wneigh>7){
        cout<<site<<endl;
        cout<<wneigh<<endl;
        perror("Coordinates.h:getneigh -> ifstatement-2");}
    return neigh_(site,wneigh);
} // ----------


void Coordinates::Numbering_lattice(){


} // ----------


void Coordinates::Numbering_DOF(){

    no_dof_=ns_*2; // 2 is for spin
    Nc_dof_.resize(ns_,2);

    // dof labeling
    int dofcount=0;

    for(int dof=0;dof<2;dof++){ //It is just SPIN
        for(int site=0;site<ns_;site++){
            Nc_dof_(site,dof)=dofcount;
            dofcount++;
        }}

} // ----------


int Coordinates::getneigh(int site,int wneigh){

    return  -10000;
} // ----------

#endif
