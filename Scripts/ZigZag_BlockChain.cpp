#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <complex>
#include <vector>
#include <math.h>
#include <iomanip>
#include <assert.h>
using namespace std;
#define PI_ 3.14159265

typedef vector< complex<double> >  Mat_1_Complex_doub;
typedef vector<Mat_1_Complex_doub> Mat_2_Complex_doub;

typedef vector<double>  Mat_1_doub;
typedef vector<Mat_1_doub> Mat_2_doub;

typedef vector<double>  Mat_1_int;
typedef vector<Mat_1_int> Mat_2_int;

int main(){


int UnitCells=6;
int Nsites;
Nsites=UnitCells*4;

int index=0;
Mat_1_doub Xpos, Ypos;
Mat_2_int Index_;
Index_.resize(UnitCells);
for(int n=0;n<UnitCells;n++){
Index_[n].resize(4);
}

for(int cell=0;cell<UnitCells;cell++){

//atom=0
Xpos.push_back( (sqrt(3)+1)*cell);
Ypos.push_back(0.5);
Index_[cell][0]=index;
index++;

//atom=1
Xpos.push_back( ((sqrt(3)+1)*cell) + sqrt(3)/2.0);
Ypos.push_back(0);
Index_[cell][1]=index;
index++;

//atom=2
Xpos.push_back( ((sqrt(3)+1)*cell) + sqrt(3)/2.0);
Ypos.push_back(1.0);
Index_[cell][2]=index;
index++;

//atom=3
Xpos.push_back( ((sqrt(3)+1)*cell) + sqrt(3));
Ypos.push_back(0.5);
Index_[cell][3]=index;
index++;

}


cout<<"index = "<<index<<endl;

Mat_1_doub sz,sx,sy;
sz.resize(Nsites);
sx.resize(Nsites);
sy.resize(Nsites);

string filein_LocalS="../LocalS_HF.txt";
ifstream LocalS_in(filein_LocalS.c_str());

string line_temp_;
int i_temp_;
            //getline(inputfile_Ansatz_file,line_temp_);

for(int line=0;line<Nsites;line++)
{
//getline(LocalS_in,line_temp_);
//stringstream line_temp_ss_(line_temp_);

LocalS_in >> i_temp_;
LocalS_in >> sz[i_temp_]>> sx[i_temp_] >> sy[i_temp_];
cout << i_temp_<<"  "<<sz[i_temp_]<<endl;
//assert(i_temp_==line);
}


string fileout_LocalS="LocalS.txt";
ofstream LocalS(fileout_LocalS.c_str());
//LocalS<<"# index   cell atom   x  y   Sz Sx Sy"<<endl;


double sz_, sx_, sy_;
for(int cell=0;cell<UnitCells;cell++){
for(int atom=0;atom<4;atom++){
sz_=sz[Index_[cell][atom]];
sx_=sx[Index_[cell][atom]];
sy_=sy[Index_[cell][atom]];
LocalS<<Index_[cell][atom]<<setw(15)<<cell<<setw(15)<<atom<<setw(15)<<Xpos[Index_[cell][atom]]<<setw(15)<<Ypos[Index_[cell][atom]]<<setw(15)<<sz_<<setw(15)<<sx_<<setw(15)<<sy_<<endl;
}
}


return 0;
}
