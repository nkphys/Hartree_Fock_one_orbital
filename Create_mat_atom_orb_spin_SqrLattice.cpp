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


int main(){

int N_atoms=1;
int N_orbs=2;

double spin_fac=1.0;
double orb_fac=1.0;
double t1_parameter=-1.0; //square lattice connection
double t2_parameter=0.0; //diagonal connections

//e1=(sqrt(3)/2,1/2)
//e2=(-sqrt(3)/2,1/2)
//e3=(0,-1)




//atom + orb*2 + spin*4

string t0_file_str = "t0_mat.txt" ;
ofstream t0_file(t0_file_str.c_str());

string t1_plus_a1_file_str = "t1_plus_a1_mat.txt" ;
ofstream t1_plus_a1_file(t1_plus_a1_file_str.c_str());

string t1_minus_a2_file_str = "t1_minus_a2_mat.txt" ;
ofstream t1_minus_a2_file(t1_minus_a2_file_str.c_str());

string t1_plus_a1_minus_a2_file_str = "t1_plus_a1_minus_a2_mat.txt" ;
ofstream t1_plus_a1_minus_a2_file(t1_plus_a1_minus_a2_file_str.c_str());

string t1_plus_a1_plus_a2_file_str = "t1_plus_a1_plus_a2_mat.txt" ;
ofstream t1_plus_a1_plus_a2_file(t1_plus_a1_plus_a2_file_str.c_str());


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
if((atom1==0) && (atom2==0) && (spin1==spin2) ){
if(orb1==orb2){
//val=t1_parameter;
//val=val*(1.0-2.0*orb1);
}
}

t0_file<<val<<" ";

}}}

t0_file<<endl;
}}}




//t1_plus_a1[2][1] c_{2}^{dag}c_{1}
//site---->neigh
//neigh
for(int spin2=0;spin2<2;spin2++){
for(int orb2=0;orb2<N_orbs;orb2++){
for(int atom2=0;atom2<N_atoms;atom2++){

//site
for(int spin1=0;spin1<2;spin1++){
for(int orb1=0;orb1<N_orbs;orb1++){
for(int atom1=0;atom1<N_atoms;atom1++){

val=0;
if((atom1==0) && (atom2==0) && (spin1==spin2) ){
if(orb1==orb2){
val=t1_parameter;

if(orb1==1){
val=val*orb_fac;
}


if(spin1==1){
val=val*spin_fac;
}

}
}



t1_plus_a1_file<<val<<" ";

}}}

t1_plus_a1_file<<endl;
}}}



//t1_minus_a2[2][1] c_{2}^{dag}c_{1}
//site---->neigh
//neigh
for(int spin2=0;spin2<2;spin2++){
for(int orb2=0;orb2<N_orbs;orb2++){
for(int atom2=0;atom2<N_atoms;atom2++){

//site
for(int spin1=0;spin1<2;spin1++){
for(int orb1=0;orb1<N_orbs;orb1++){
for(int atom1=0;atom1<N_atoms;atom1++){

val=0;
if((atom1==0) && (atom2==0) && (spin1==spin2) ){
if(orb1==orb2){
val=t1_parameter;

if(orb1==1){
val=val*orb_fac;
}


if(spin1==1){
val=val*spin_fac;
}


}
}


t1_minus_a2_file<<val<<" ";

}}}

t1_minus_a2_file<<endl;
}}}


//t_plus_a1_pluss_a2[2][1] c_{2}^{dag}c_{1}
//site---->neigh
//neigh
for(int spin2=0;spin2<2;spin2++){
for(int orb2=0;orb2<N_orbs;orb2++){
for(int atom2=0;atom2<N_atoms;atom2++){

//site
for(int spin1=0;spin1<2;spin1++){
for(int orb1=0;orb1<N_orbs;orb1++){
for(int atom1=0;atom1<N_atoms;atom1++){

val=0;
if((atom1==0) && (atom2==0) && (spin1==spin2) ){
if(orb1==orb2){
val=t2_parameter;
//val=val*(1.0-2.0*orb1);
}
}

t1_plus_a1_plus_a2_file<<val<<" ";

}}}

t1_plus_a1_plus_a2_file<<endl;
}}}


//t_plus_a1_minus_a2[2][1] c_{2}^{dag}c_{1}
//site---->neigh
//neigh
for(int spin2=0;spin2<2;spin2++){
for(int orb2=0;orb2<N_orbs;orb2++){
for(int atom2=0;atom2<N_atoms;atom2++){

//site
for(int spin1=0;spin1<2;spin1++){
for(int orb1=0;orb1<N_orbs;orb1++){
for(int atom1=0;atom1<N_atoms;atom1++){

val=0;
if((atom1==0) && (atom2==0) && (spin1==spin2) ){
if(orb1==orb2){
val=t2_parameter;
//val=val*(1.0-2.0*orb1);
}
}




t1_plus_a1_minus_a2_file<<val<<" ";

}}}

t1_plus_a1_minus_a2_file<<endl;
}}}


return 0;
}
