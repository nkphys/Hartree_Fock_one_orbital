//--------This file creates a tex file for Dice-Lattice and its Current Configurations--------//

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


string decouple(int m, int Lx){

int ix,iy;
double x,y,a;
a=0.5;

string xy_,xstr,ystr;

ix= ( m %Lx);
iy= ( (m-ix)/(Lx) );

x=1.0*ix*a+0.5*iy*a;
y=sqrt(3)*0.5*iy*a;

ostringstream convertx;
ostringstream converty;
convertx << x;
xstr = convertx.str();
converty << y;
ystr = converty.str();

xy_= "(" + xstr + "," + ystr + ")";

return xy_;
}

//----------------Main----------------------//
int main(int argc, char** argv){

int Lx=atoi(argv[1]);
int Ly=atoi(argv[2]);
string file_input;
stringstream ss;
ss<<argv[3];
ss>> file_input;

double eps=atof(argv[4]);
double Factor_=1.0;

Mat_1_int m_val,xp_val,yp_val;
Mat_1_doub S_val,Sp_val,P_val,N_val;
m_val.clear();S_val.clear();Sp_val.clear();P_val.clear();N_val.clear();
xp_val.clear();yp_val.clear();

int temp_m,temp_xp,temp_yp;
double temp_S,temp_null;

string file_in_current_updn = file_input + ".txt";
ifstream file_in(file_in_current_updn.c_str());
string line;

int site_temp; //x + y*lx
double temp_up, temp_dn;
getline(file_in,line);
while(getline(file_in,line) ){
stringstream line_ss;
line_ss << line;
line_ss>>site_temp>>temp_up>>temp_dn;
//cout<<site_temp<<"  "<<temp_up<<"   "<<temp_dn<<endl;
temp_xp=site_temp%Lx;  //site = ix + iy*lx;
temp_yp=(site_temp - temp_xp)/Lx;
//cout<<site_temp<<"  "<<temp_xp<<"  "<<temp_yp<<"  "<<temp_up<<"   "<<temp_dn<<endl;
xp_val.push_back(temp_xp);
yp_val.push_back(temp_yp);
m_val.push_back(site_temp);
//temp_S=temp_S-(0.5*0.5);
S_val.push_back(temp_up + temp_dn);
N_val.push_back(temp_null);
}

//-------------------------------------------//
P_val.resize(S_val.size());
for (int x=0;x<S_val.size();x++){
	P_val[x]=S_val[x];
}

double max_S_val;
max_S_val=0.0;

for(int x=0;x<S_val.size();x++){
	if(abs(P_val[0]) < abs(P_val[x])){
		P_val[0] = abs(P_val[x]);
	}
}
max_S_val=P_val[0];

Sp_val.resize(m_val.size());
for(int x=0;x<m_val.size();x++){
        Sp_val[x]=Factor_*abs(S_val[x])/max_S_val;
}

//------------------------------------------//
string file_out_tex = file_input + ".tex";
ofstream file_out( file_out_tex.c_str());

file_out<<"\\documentclass[tikz,border=10pt]{standalone}"<<endl;
file_out<<"\\usepackage{tikz}"<<endl;
file_out<<"\\usetikzlibrary{shapes}"<<endl;
file_out<<"\\usetikzlibrary{positioning,arrows.meta}"<<endl;
file_out<<"\\usepackage[utf8x]{inputenc}"<<endl;
file_out<<""<<endl;
file_out<<""<<endl;
file_out<<"\\begin{document}"<<endl;
file_out<<"\\begin{tikzpicture}"<<endl;
file_out<<""<<endl;

double x_val,y_val,a;
int site;
a=0.5;

//\node[text width=0.04 mm] at (1.7,0.05) {\fontsize{4}{4}\selectfont $\epsilon$=1.0};
//file_out<<"\\node[text width=0.04mm] at ("<<1.7<<","<<0.05<<") { \\fontsize{4}{4}\\selectfont $\\epsilon$="<<eps<<"};"<<endl;

for(int x=0;x<Lx;x++){
	for(int y=0;y<Ly;y++){
		x_val=1.0*x*a+0.5*y*a;
                y_val=sqrt(3)*0.5*y*a;

		if(x!=Lx-1){
			file_out<<"\\draw[line width=0.1 mm,densely dotted,black] "<<"("<<x_val<<","<<y_val<<")"<<" -- "<<"("<<x_val+1.0*a<<","<<y_val<<")"<<";"<<endl;
		}
		if(y!=Ly-1){
			file_out<<"\\draw[line width=0.1 mm,densely dotted,black] "<<"("<<x_val<<","<<y_val<<")"<<" -- "<<"("<<x_val+0.5*a<<","<<y_val+sqrt(3)*0.5*a<<")"<<";"<<endl;
		}
		if(x!=0 && y!=Ly-1){
			file_out<<"\\draw[line width=0.1 mm,densely dotted,black] "<<"("<<x_val<<","<<y_val<<")"<<" -- "<<"("<<x_val-0.5*a<<","<<y_val+sqrt(3)*0.5*a<<")"<<";"<<endl;
		}
	}
}

int site_ref;
for(int x=0;x<Lx;x++){
        for(int y=0;y<Ly;y++){
                site=x+Lx*y;
		site_ref=x+Lx*y;

                if(S_val[site_ref]>0){
                        file_out<<"\\filldraw[fill=blue,draw=blue,opacity="<<Sp_val[site_ref]*Sp_val[site_ref]*sqrt(Sp_val[site_ref])<<"]" <<decouple(m_val[site],Lx)<<" circle (0.7mm);"<<endl;
                }

                if(S_val[site_ref]<0){
                        file_out<<"\\filldraw[fill=red,draw=red,opacity="<<Sp_val[site_ref]<<"]" <<decouple(m_val[site],Lx)<<" circle (0.7mm);"<<endl;
                }

//		cout<<site<<"	"<<decouple(m_val[site],Ly)<<endl;
/*		if(S_val[site]=P_val[0]){
			file_out<<"\\filldraw[fill=green,draw=blue]" <<"("<<x_val<<","<<y_val<<") "<<" circle ("<<0.8<<"mm);"<<endl;
		}
*/
        }
}

file_out<<""<<endl;
file_out<<""<<endl;

file_out<<"\\end{tikzpicture}"<<endl;
file_out<<"\\end{document}"<<endl;

return 0;
}

