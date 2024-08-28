#include "Parameters_G2dLattice.h"


#ifndef Coordinates_G2dLattice_class
#define Coordinates_G2dLattice_class
class Coordinates_G2dLattice {

public:
    Coordinates_G2dLattice(int lx, int ly, int n_atoms, int n_orbs)
        :lx_(lx),ly_(ly),n_atoms_(n_atoms),n_orbs_(n_orbs)
    {
        Numbering(); 
    }

    enum {PX=0,MX,PY,MY,PXPY,MXPY,MXMY,PXMY};	// not needed - but keep as a "key"
    void Numbering();
    int indx_basiswise(int i);
    int indy_basiswise(int i);
    int indorb_basiswise(int i);
    int Nbasis(int x, int y, int gamma);

    int indx_cellwise(int i);
    int indy_cellwise(int i);
    int Ncell(int x, int y);

    int neigh(int cell, int wneigh);
    int getneigh(int cell,int wneigh);

    int lx_,ly_,n_orbs_,n_atoms_, nbasis_, ncells_;
    Mat_1_int indx_basiswise_,indy_basiswise_, indorb_basiswise_, indatom_basiswise_;
    Mat_3_int Nbasis_;
    Matrix<int> Ncell_, neigh_;
    Mat_1_int indx_cellwise_, indy_cellwise_;

    bool HIT_X_BC, HIT_Y_BC;

};

/*  
 * ***********
 *  Functions in Class Coordinates_G2dLattice ------
 *  ***********
*/  

int Coordinates_G2dLattice::indx_basiswise(int i){
    if(i>nbasis_-1){perror("Coordinates_G2dLattice.h:x-coordinate of lattice excede limit");}
    return indx_basiswise_[i];
    // ----------
}


int Coordinates_G2dLattice::indy_basiswise(int i){
    if(i>nbasis_-1){perror("Coordinates_G2dLattice.h:y-coordinate of lattice excede limit");}
    return indy_basiswise_[i];
    // ----------
}

int Coordinates_G2dLattice::indx_cellwise(int i){
    if(i>ncells_-1){perror("Coordinates_G2dLattice.h:x-coordinate of lattice excede limit");}
    return indx_cellwise_[i];
    // ----------
}


int Coordinates_G2dLattice::indy_cellwise(int i){
    if(i>ncells_-1){perror("Coordinates_G2dLattice.h:y-coordinate of lattice excede limit");}
    return indy_cellwise_[i];
    // ----------
}

int Coordinates_G2dLattice::Nbasis(int x, int y, int gamma){
    if(!( (x<lx_&& y<ly_) &&  gamma<n_atoms_*n_orbs_)){perror("Coordinates_G2dLattice.h:ith-sitelabel of lattice excede limit");}
    return Nbasis_[x][y][gamma];
    // ----------
}


int Coordinates_G2dLattice::Ncell(int x, int y){
    if(!(x<lx_&& y<ly_)){
        cout<<x<<"   "<<y<<endl;
        perror("Coordinates_G2dLattice.h:ith-sitelabel of lattice excede limit : Ncell(int x, int y) ");
    }
    return Ncell_(x,y);
    // ----------
}


int Coordinates_G2dLattice::neigh(int cell, int wneigh){
    if(cell> (lx_*ly_)-1 || wneigh>=17){
        cout<<cell<<"  "<<wneigh<<endl;
        perror("Coordinates_G2dLattice.h:getneigh -> ifstatement-10");}
    return neigh_(cell,wneigh);
} // ----------


void Coordinates_G2dLattice::Numbering(){

    nbasis_=lx_*ly_*n_orbs_*n_atoms_;
    ncells_= lx_*ly_;

    indx_basiswise_.clear(); 	indx_basiswise_.resize(nbasis_);
    indy_basiswise_.clear();	indy_basiswise_.resize(nbasis_);
    indorb_basiswise_.clear();  indorb_basiswise_.resize(nbasis_);
    indatom_basiswise_.clear();  indatom_basiswise_.resize(nbasis_);

    indx_cellwise_.clear();indx_cellwise_.resize(ncells_);
    indy_cellwise_.clear();indy_cellwise_.resize(ncells_);

    Nbasis_.resize(lx_);
    for(int ix=0;ix<lx_;ix++){
        Nbasis_[ix].resize(ly_);
        for(int iy=0;iy<ly_;iy++){
            Nbasis_[ix][iy].resize(n_orbs_*n_atoms_);
        }
    }

    neigh_.resize(ncells_,17);


    //basis labeling
    int icount=0;
    for(int j=0;j<ly_;j++){
        for(int i=0;i<lx_;i++){
            for(int orb=0;orb<n_orbs_;orb++){
                for(int atom=0;atom<n_atoms_;atom++){
                indx_basiswise_[icount]=i;
                indy_basiswise_[icount]=j;
                indorb_basiswise_[icount]=orb;
                indatom_basiswise_[icount]=atom;
                Nbasis_[i][j][atom + n_atoms_*orb]=icount;
                icount++;
                }
            }
        }
    }



    Ncell_.resize(lx_,ly_);
    //cell labeling
    icount=0;
    for(int j=0;j<ly_;j++){
        for(int i=0;i<lx_;i++){
            indx_cellwise_[icount]=i;
            indy_cellwise_[icount]=j;
            Ncell_(i,j)=icount;
            icount++;
        }
    }


    // Neighbors for each unit cell
    for(int i=0;i<ncells_;i++){ 	// ith site
        for(int j=0;j<17;j++) {		// jth neighbor
            neigh_(i,j)=getneigh(i,j);
        }
    }


} // ----------


int Coordinates_G2dLattice::getneigh(int site,int wneigh){
    if(site>ncells_-1 || wneigh>17){perror("Coordinates_G2dLattice.h:getneigh -> ifstatement-1");}
    int nx=indx_cellwise(site);
    int ny=indy_cellwise(site);
    int mx=0;
    int my=0;
    int mx_,my_;

    // Nearest Neighbours
    if(wneigh==0){ //PX
        mx=(nx+1);
        my=ny;
    }
    if(wneigh==1){ //MX
        mx=(nx-1);
        my=ny;
    }
    if(wneigh==2){ //PY
        mx=nx;
        my=(ny+1);
    }
    if(wneigh==3){ //MY
        mx=nx;
        my=(ny-1);
    }


    // Next-Nearest!
    if(wneigh==4){ //PXPY
        mx=(nx+1);
        my=(ny+1);
    }
    if(wneigh==5){ //MXPY
        mx=(nx-1);
        my=(ny+1);
    }
    if(wneigh==6){ //MXMY
        mx=(nx-1);
        my=(ny-1);
    }
    if(wneigh==7){ //PXMY
        mx=(nx+1);
        my=(ny-1);
    }

    if(wneigh==8){ //2PX MY
        mx=(nx+2);
        my=(ny-1);
    }
    if(wneigh==9){ //MX 2PY
        mx=(nx-1);
        my=(ny+2);
    }

    if(wneigh==10){ //2PX
        mx=(nx+2);
        my=(ny);
    }
    if(wneigh==11){ //2PY
        mx=(nx);
        my=(ny+2);
    }
    if(wneigh==12){ //2PX 2MY
        mx=(nx+2);
        my=(ny-2);
    }

    if(wneigh==13){ //2MX PY
        mx=(nx-2);
        my=(ny+1);
    }

    if(wneigh==14){ //2MY PX
        mx=(nx+1);
        my=(ny-2);
    }

    if(wneigh==15){ //2PY 2MX
        mx=(nx-2);
        my=(ny+2);
    }

    if(wneigh==16){ //2MY
        mx=(nx);
        my=(ny-2);
    }


    if(mx>=lx_ || mx<0){
        HIT_X_BC=true;
    }
    else{
        HIT_X_BC=false;
    }

    if(my>=ly_ || my<0){
        HIT_Y_BC=true;
    }
    else{
        HIT_Y_BC=false;
    }


    mx = (mx+lx_)%lx_;
    my = (my+ly_)%ly_;

//    //t4
//    if(wneigh==13){ //2PX+PY
//        mx=(nx+lx_+2)%(lx_);
//        my=(ny+ly_+1)%(ly_);
//    }
//    if(wneigh==14){ //2PY-MX
//        mx=(nx+lx_-1)%(lx_);
//        my=(ny+ly_+2)%(ly_);
//    }
//    if(wneigh==15){ //2PX 2MY
//        mx=(nx+lx_+2)%(lx_);
//        my=(ny+ly_-2)%(ly_);
//    }



    return Ncell(mx,my); //Nc(mx,my);
} // ----------

#endif
