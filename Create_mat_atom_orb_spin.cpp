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

string t2_file_str = "t2_mat.txt" ;
ofstream t2_file(t2_file_str.c_str());

string t3_file_str = "t3_mat.txt" ;
ofstream t3_file(t3_file_str.c_str());

complex<double> val=0.0;

//t0[2][1] c_{2}^{dag}c_{1}
//bond-1
for(int spin2=0;spin2<2;spin2++){
for(int orb2=0;orb2<2;orb2++){
for(int atom2=0;atom2<2;atom2++){

for(int spin1=0;spin1<2;spin1++){
for(int orb1=0;orb1<2;orb1++){
for(int atom1=0;atom1<2;atom1++){


val=0;
//p1= ( (sqrt(3)/2)px + (1/2)py )
//p1_0(A)^dagg p1_1(B)

// 3/4 pxA pxB  +  sqrt(3)/4 pxA pyB  + 1/4 pyA pyB  + sqrt(3)/4 pyA pxB

if((atom1==1) && (atom2==0) && (spin1==spin2) ){

if(orb1==0 && orb2==0){
val=3.0/4.0;
}
if(orb1==0 && orb2==1){
val=sqrt(3.0)/4.0;
}
if(orb1==1 && orb2==0){
val=sqrt(3.0)/4.0;
}
if(orb1==1 && orb2==1){
val=1.0/4.0;
}

}

t0_file<<val<<" ";

}}}

t0_file<<endl;
}}}




//t1_plus_a1
for(int spin2=0;spin2<2;spin2++){
for(int orb2=0;orb2<2;orb2++){
for(int atom2=0;atom2<2;atom2++){

for(int spin1=0;spin1<2;spin1++){
for(int orb1=0;orb1<2;orb1++){
for(int atom1=0;atom1<2;atom1++){


val=0;
//p2= ( (-sqrt(3)/2)px + (1/2)py )
//p2_0(A)^dagg p2_1(B)

// 3/4 pxA pxB  - sqrt(3)/4 pxA pyB  + 1/4 pyA pyB  - sqrt(3)/4 pyA pxB

if((atom1==1) && (atom2==0) && (spin1==spin2) ){

if(orb1==0 && orb2==0){
val=3.0/4.0;
}
if(orb1==0 && orb2==1){
val=-sqrt(3.0)/4.0;
}
if(orb1==1 && orb2==0){
val=-sqrt(3.0)/4.0;
}
if(orb1==1 && orb2==1){
val=1.0/4.0;
}

}


t1_plus_a1_file<<val<<" ";

}}}

t1_plus_a1_file<<endl;
}}}



//t1_minus_a2

for(int spin2=0;spin2<2;spin2++){
for(int orb2=0;orb2<2;orb2++){
for(int atom2=0;atom2<2;atom2++){

for(int spin1=0;spin1<2;spin1++){
for(int orb1=0;orb1<2;orb1++){
for(int atom1=0;atom1<2;atom1++){



val=0;
//p3= -py
//p3_0(A)^dagg p2_1(B)
//pyA pyB
//pyB pyA

if((atom1==0) && (atom2==1) && (spin1==spin2) ){

if(orb1==1 && orb2==1){
val=1.0;
}
}



t1_minus_a2_file<<val<<" ";

}}}

t1_minus_a2_file<<endl;
}}}



//t2[2][1] c_{2}^{dag}c_{1}
for(int spin2=0;spin2<2;spin2++){
for(int orb2=0;orb2<2;orb2++){
for(int atom2=0;atom2<2;atom2++){

for(int spin1=0;spin1<2;spin1++){
for(int orb1=0;orb1<2;orb1++){
for(int atom1=0;atom1<2;atom1++){


val=0.0;

t2_file<<val<<" ";

}}}

t2_file<<endl;
}}}



//t3[2][1] c_{2}^{dag}c_{1}

for(int spin2=0;spin2<2;spin2++){
for(int orb2=0;orb2<2;orb2++){
for(int atom2=0;atom2<2;atom2++){

for(int spin1=0;spin1<2;spin1++){
for(int orb1=0;orb1<2;orb1++){
for(int atom1=0;atom1<2;atom1++){


val=0.0;

t3_file<<val<<" ";

}}}

t3_file<<endl;
}}}



return 0;
}
