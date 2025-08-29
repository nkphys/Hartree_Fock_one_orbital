#include <iostream>
#include <math.h>
#include <time.h>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <assert.h>
//#include "tensor.h"
#include <algorithm>
using namespace std;
#include <vector>
#include <complex>

typedef vector< double >  Mat_1_doub;
typedef vector< int >  Mat_1_int;

typedef vector< complex<double> >  Mat_1_Complex_doub;
typedef vector< Mat_1_Complex_doub >  Mat_2_Complex_doub;

int main(){

int N_atoms=1;
int N_orbs=3;


double lambda_SOC=0.2;
complex<double> iota_complex = complex<double>(0.0, 1.0);
complex<double> one_complex = complex<double>(1.0, 0.0);


string t0_file_str = "t0_mat_with_SOC.txt" ;
ofstream t0_file(t0_file_str.c_str());

 /* Convention
 0==xz
 1==yz
 2==xy
 0==up
 1==dn
 index = orb + spin*3 + site*6;
yz_up(site=0),xz_up(site=0),xy_up(site=0), yz_dn(site=0),xz_dn(site=0),xy_dn(site=0)....site(n)...
*/

    Mat_2_Complex_doub A_SOC;
    A_SOC.resize(6);
    for(int i=0;i<6;i++){
        A_SOC[i].resize(6);
    }
    A_SOC[1][0]=iota_complex;A_SOC[0][1]=-1.0*iota_complex;
    A_SOC[0][5]=iota_complex;A_SOC[5][0]=-1.0*iota_complex;
    A_SOC[1][5]=-1.0*one_complex;A_SOC[5][1]=-1.0*one_complex;

    A_SOC[4][3]=-1.0*iota_complex;A_SOC[3][4]=iota_complex;
    A_SOC[3][2]=iota_complex;A_SOC[2][3]=-1.0*iota_complex;
    A_SOC[4][2]=one_complex;A_SOC[2][4]=one_complex;


complex<double> val=0.0;

//t0[2][1] c_{2}^{dag}c_{1}
//bond-1
for(int spin2=0;spin2<2;spin2++){
for(int orb2=0;orb2<N_orbs;orb2++){
for(int atom2=0;atom2<N_atoms;atom2++){

for(int spin1=0;spin1<2;spin1++){
for(int orb1=0;orb1<N_orbs;orb1++){
for(int atom1=0;atom1<N_atoms;atom1++){

val=0;
if((atom1==0) && (atom2==0) ){
val = 0.5*lambda_SOC*A_SOC[orb2 + spin2*N_orbs ][ orb1 + spin1*N_orbs ];
}

t0_file<<val<<" ";

}}}

t0_file<<endl;
}}}



return 0;
}
