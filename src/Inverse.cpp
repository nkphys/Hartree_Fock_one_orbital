#include <algorithm>
#include <functional>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <random>
#include <complex>
#include <cmath>
#include <cassert>

#include "tensor_type.h"
#include "random"
#include "Matrix.h"


#define PI acos(-1.0)


using namespace std;


extern "C" void zgetri_(int *, std::complex<double> *, int *, int *, std::complex<double> *, int *, int * );
extern "C" void zgetrf_(int *, int* , std::complex<double> *, int* , int *, int *);


void Inverse(Matrix<complex<double>> & A){

int M=A.n_col();;
vector<complex<double>> work(M);
vector<int> ipiv(M);
int n,info;

n=M;

zgetrf_(&n, &n, &(A(0,0)), &n, &(ipiv[0]), &info);

if(info != 0){
	std::cerr<<"zgetrf failed"<<endl;
}

int lwork=64*M;
work.resize(lwork);
//cout<<lwork<<endl;
zgetri_(&n, &(A(0,0)), &n,  &(ipiv[0]), &(work[0]), &lwork, &info);

if(info != 0){
	std::cerr<<"zgetri failed"<<endl;
}



}




int main(){

Matrix<complex<double>> Gamma,Gamma_inv;
Gamma.resize(4,4);

Gamma(0,0)=1.5; Gamma(0,1)=2.3;Gamma(0,2)=2.3;Gamma(0,3)=6.3;
Gamma(1,0)=1.7; Gamma(1,1)=5.3;Gamma(1,2)=2.7;Gamma(1,3)=5.3;
Gamma(2,0)=4.8;Gamma(2,1)=12.3;Gamma(2,2)=6.3;Gamma(2,3)=7.3;
Gamma(3,0)=complex<double>(2.4,1.3);Gamma(3,1)=22.3;Gamma(3,2)=5.3;Gamma(3,3)=2.9;


Gamma.print();


Gamma_inv=Gamma;

Inverse(Gamma_inv);

Gamma_inv.print();


Matrix<complex<double>> Iden;
Iden.resize(Gamma.n_row(), Gamma.n_col());

for(int i=0;i<Gamma.n_row();i++){
	for(int j=0;j<Gamma.n_row();j++){
	complex<double> val=0.0;
	for(int k=0;k<Gamma.n_row();k++){
	val += Gamma(i,k)*Gamma_inv(k,j);
	}
	Iden(i,j)=val;
	}
}

Iden.print();

return 0;
}
