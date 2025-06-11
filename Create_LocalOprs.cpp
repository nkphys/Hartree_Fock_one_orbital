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



double delta(int i, int j){

double val;
if(i==j){
val=1.0;
}
else{
val=0.0;
}
return val;
}



int main(){

int N_atoms=1;
int N_orbs=3;


//atom + orb*2 + spin*4
complex<double> val;
string Sz_file_str_temp = "SzLocal" ;
string Sp_file_str_temp = "SpLocal" ;
string Sm_file_str_temp = "SmLocal" ;

string Tauz_file_str_temp = "Tauz";
string Taup_file_str_temp = "Taup";
string Taum_file_str_temp = "Taum";




//Tauz
//atom + orb*2 + spin*4

//Tz = Sigma{2,1} c_{2}^{dag}c_{1}
//bond-1

for(int atom=0;atom<N_atoms;atom++){
string Tz_file_str = Tauz_file_str_temp + "_atom" + to_string(atom) + "_total.txt";
ofstream Tz_file(Tz_file_str.c_str());

for(int spin2=0;spin2<2;spin2++){
for(int orb2=0;orb2<N_orbs;orb2++){
for(int atom2=0;atom2<N_atoms;atom2++){

for(int spin1=0;spin1<2;spin1++){
for(int orb1=0;orb1<N_orbs;orb1++){
for(int atom1=0;atom1<N_atoms;atom1++){

complex<double> val_temp=0.0;
if(orb1 !=2 && orb2!=2){ //not xy=2
val=delta(spin1,spin2);
val=val*delta(atom1,atom2)*delta(atom1,atom);
val=val*(delta(orb1,orb2))*(0.5 - orb1)*1.0;
}
else{
val=0.0;
}
val_temp =val;


Tz_file<<val_temp<<" ";
}}}
Tz_file<<endl;
}}}

}





//Taup
//atom + orb*2 + spin*4

//Tp = Sigmap{2,1} c_{2}^{dag}c_{1}
//bond-1

for(int atom=0;atom<N_atoms;atom++){
string Tp_file_str = Taup_file_str_temp + "_atom" + to_string(atom) + "_total.txt";
ofstream Tp_file(Tp_file_str.c_str());

for(int spin2=0;spin2<2;spin2++){
for(int orb2=0;orb2<N_orbs;orb2++){
for(int atom2=0;atom2<N_atoms;atom2++){

for(int spin1=0;spin1<2;spin1++){
for(int orb1=0;orb1<N_orbs;orb1++){
for(int atom1=0;atom1<N_atoms;atom1++){

complex<double> val_temp=0.0;
if(orb1 !=2 && orb2!=2){ //not xy=2
val=delta(spin1,spin2);
val=val*delta(atom1,atom2)*delta(atom1,atom);

if(orb1==0 && orb2==1){
val=val*1.0;
}
else{
val=0.0;
}

}
else{
val=0.0;
}
val_temp =val;


Tp_file<<val_temp<<" ";
}}}
Tp_file<<endl;
}}}

}






//Taum
//atom + orb*2 + spin*4

//Tm = Sigmam{2,1} c_{2}^{dag}c_{1}
//bond-1

for(int atom=0;atom<N_atoms;atom++){
string Tm_file_str = Taum_file_str_temp + "_atom" + to_string(atom) + "_total.txt";
ofstream Tm_file(Tm_file_str.c_str());

for(int spin2=0;spin2<2;spin2++){
for(int orb2=0;orb2<N_orbs;orb2++){
for(int atom2=0;atom2<N_atoms;atom2++){

for(int spin1=0;spin1<2;spin1++){
for(int orb1=0;orb1<N_orbs;orb1++){
for(int atom1=0;atom1<N_atoms;atom1++){

complex<double> val_temp=0.0;
if(orb1 !=2 && orb2!=2){ //not xy=2
val=delta(spin1,spin2);
val=val*delta(atom1,atom2)*delta(atom1,atom);

if(orb1==1 && orb2==0){
val=val*1.0;
}
else{
val=0.0;
}

}
else{
val=0.0;
}
val_temp =val;


Tm_file<<val_temp<<" ";
}}}
Tm_file<<endl;
}}}

}










for(int orb=0;orb<N_orbs;orb++){
for(int atom=0;atom<N_atoms;atom++){

string Sz_file_str = Sz_file_str_temp + "_atom" + to_string(atom) + "_orb" + to_string(orb) + ".txt";
ofstream Sz_file(Sz_file_str.c_str());

//Sz = Sigma{2,1} c_{2}^{dag}c_{1}
//bond-1
for(int spin2=0;spin2<2;spin2++){
for(int orb2=0;orb2<N_orbs;orb2++){
for(int atom2=0;atom2<N_atoms;atom2++){

for(int spin1=0;spin1<2;spin1++){
for(int orb1=0;orb1<N_orbs;orb1++){
for(int atom1=0;atom1<N_atoms;atom1++){


val=delta(orb1,orb2)*delta(orb1,orb);
val=val*delta(atom1,atom2)*delta(atom1,atom);

val=val*(delta(spin1,spin2))*(0.5 - spin1)*1.0;

Sz_file<<val<<" ";
}}}
Sz_file<<endl;
}}}

}
}




//local_total_oprs
//atom + orb*2 + spin*4

//Sz = Sigma{2,1} c_{2}^{dag}c_{1}
//bond-1

for(int atom=0;atom<N_atoms;atom++){
string Sz_file_str = Sz_file_str_temp + "_atom" + to_string(atom) + "_total.txt";
ofstream Sz_file(Sz_file_str.c_str());

for(int spin2=0;spin2<2;spin2++){
for(int orb2=0;orb2<N_orbs;orb2++){
for(int atom2=0;atom2<N_atoms;atom2++){

for(int spin1=0;spin1<2;spin1++){
for(int orb1=0;orb1<N_orbs;orb1++){
for(int atom1=0;atom1<N_atoms;atom1++){

complex<double> val_temp=0.0;
for(int orb=0;orb<N_orbs;orb++){
val=delta(orb1,orb2)*delta(orb1,orb);
val=val*delta(atom1,atom2)*delta(atom1,atom);
val=val*(delta(spin1,spin2))*(0.5 - spin1)*1.0;

val_temp +=val;
}


Sz_file<<val_temp<<" ";
}}}
Sz_file<<endl;
}}}


}




//Sp = Sigma_plus{2,1} c_{2}^{dag}c_{1}
//bond-1

for(int atom=0;atom<N_atoms;atom++){
string Sp_file_str = Sp_file_str_temp + "_atom" + to_string(atom) + "_total.txt";
ofstream Sp_file(Sp_file_str.c_str());

for(int spin2=0;spin2<2;spin2++){
for(int orb2=0;orb2<N_orbs;orb2++){
for(int atom2=0;atom2<N_atoms;atom2++){

for(int spin1=0;spin1<2;spin1++){
for(int orb1=0;orb1<N_orbs;orb1++){
for(int atom1=0;atom1<N_atoms;atom1++){

complex<double> val_temp=0.0;
for(int orb=0;orb<N_orbs;orb++){
val=delta(orb1,orb2)*delta(orb1,orb);
val=val*delta(atom1,atom2)*delta(atom1,atom);

if(spin2==1 && spin1==0){
val=val*1.0;
}
else{
val=0.0;
}

val_temp +=val;
}


Sp_file<<val_temp<<" ";
}}}
Sp_file<<endl;
}}}


}




//Sm_total

for(int atom=0;atom<N_atoms;atom++){
string Sm_file_str = Sm_file_str_temp + "_atom" + to_string(atom) + "_total.txt";
ofstream Sm_file(Sm_file_str.c_str());

for(int spin2=0;spin2<2;spin2++){
for(int orb2=0;orb2<N_orbs;orb2++){
for(int atom2=0;atom2<N_atoms;atom2++){

for(int spin1=0;spin1<2;spin1++){
for(int orb1=0;orb1<N_orbs;orb1++){
for(int atom1=0;atom1<N_atoms;atom1++){

complex<double> val_temp=0.0;
for(int orb=0;orb<N_orbs;orb++){
val=delta(orb1,orb2)*delta(orb1,orb);
val=val*delta(atom1,atom2)*delta(atom1,atom);

if(spin2==0 && spin1==1){
val=val*1.0;
}
else{
val=0.0;
}

val_temp +=val;
}


Sm_file<<val_temp<<" ";
}}}
Sm_file<<endl;
}}}


}







return 0;
}
