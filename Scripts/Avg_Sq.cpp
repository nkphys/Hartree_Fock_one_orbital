#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <complex>
#include <vector>
#include <math.h>
#include <assert.h>
using namespace std;
#define PI_ 3.14159265

typedef vector< complex<double> >  Mat_1_Complex_doub;
typedef vector<Mat_1_Complex_doub> Mat_2_Complex_doub;

typedef vector<double>  Mat_1_doub;
typedef vector<Mat_1_doub> Mat_2_doub;

complex<double> reading_pair(string pair_str){
	
		complex<double> value;
           double part_one, part_two;
          int pair_str_Length=pair_str.length()-1;
          int pos = pair_str.find(",");
        string one_str=pair_str.substr(1,pos-1);
        string two_str=pair_str.substr(pos+1,pair_str_Length-1-pos);

        stringstream onestream(one_str);
        stringstream twostream(two_str);
        onestream>>part_one;
        twostream>>part_two;
        
        value.real(part_one);
        value.imag(part_two);
 
	return value;
	
	}



int main(){
	

	

//Lx_Ly_16x16_PBC/NP_1024/Seed_AFM_Pi_Exc_large_U/U_15.0/JbyU_0.25/lambda_0.15/Sq.txt
//Lx_Ly_12x12_PBC_Run2/NP_576/Seed_1/U_4.0/JbyU_0.25/lambda_3.0/dis_str_6.75/dis_seed_4/Sq.txt
	
	double tmp_n,tmp_tm;
	string tmp_line;
	int tmp_site;

	int Length =12;
	char Length_char[50];
	sprintf(Length_char,"%d",Length);
	
	int N_total = 576;
	char N_total_char[50];
	sprintf(N_total_char,"%d",N_total);
	
	
	complex<double> iota;iota.real(0);iota.imag(1);
	

	double Del_pi;
	double Del_total;

	double U[1] = {4.0};//{0.049, 0.098, 0.245, 0.98, 1.96, 2.45, 4.9, 12.25, 24.5, 49, 98};
	int U_nos=1;
	char U_char[50];


	//U=15.0
	//double Lambda[20]={0.0 ,0.1 ,0.15, 0.2, 0.3 ,0.4 ,0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1 ,1.2 ,1.3 ,1.4 ,1.5, 1.6 ,2.0 ,3.0};
	//string seed[20]={"3", "5","AFM_Pi_Exc_large_U","AFM_Pi_Exc_large_U","AFM_Pi_Exc_large_U","AFM_Pi_Exc_large_U","AFM_Pi_Exc_large_U","AFM_Pi_Exc_large_U","AFM_Pi_Exc_large_U","AFM_Pi_Exc_large_U","AFM_Pi_Exc_large_U","AFM_Pi_Exc_large_U","AFM_Pi_Exc_large_U","AFM_Pi_Exc_large_U","AFM_Pi_Exc_large_U","AFM_Pi_Exc_large_U","AFM_Pi_Exc_large_U","AFM_Pi_Exc_large_U","AFM_Pi_Exc_large_U","AFM_Pi_Exc_large_U"};

	//U=2.0
	//double Lambda[17]= {0.0, 0.5, 1.0, 1.5, 1.6, 1.7, 1.8,1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 3.0, 3.5, 4.0};//,0.3,0.4,0.5};//{0.0, 0.1, 0.2 ,0.3 ,0.4 ,0.5, 0.7, 0.8, 1.0, 2.0, 3.0};
	//string seed[17]={"10", "4", "AFM_Pi_Exc_small_U", "9", "9", "9", "6", "10", "AFM_Pi_Exc_small_U", "AFM_Pi_Exc_small_U", "3", "5", "10", "9", "9", "9"};


	double Lambda[1] = {3.0};
	int Lambda_nos=1;
	char Lambda_char[50];
	

	double DisStr[29] = {0.0, 0.25,  0.5, 0.75,  1.0, 1.25,  1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0, 5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0};
        int Dis_Str_nos=29;
        char DisStr_char[50];

	ofstream file_check("file_check.txt");



	 for (int U_i=0;U_i<U_nos;U_i++){
        //cout<<"here 9"<<endl;

        if(U[U_i]<=0.245 || U[U_i]==2.165 || U[U_i]==3.897){
        sprintf(U_char,"%.3f",U[U_i]);}
        if(U[U_i]==4.33 || U[U_i]==1.96 || U[U_i]==2.45 || U[U_i]==12.25 || U[U_i]==0.98 || U[U_i]==3.25 || U[U_i]==21.65){
        sprintf(U_char,"%.2f",U[U_i]);}
        if( U[U_i]==4.9 ||U[U_i]==24.5 ) {
        sprintf(U_char,"%.1f",U[U_i]);}
        if (U[U_i]==49 || U[U_i]==98) {
        sprintf(U_char,"%.0f",U[U_i]);
        }

        sprintf(U_char,"%.1f",U[U_i]);



	for(int lambda_i=0;lambda_i<Lambda_nos;lambda_i++){	
        sprintf(Lambda_char,"%.1f",Lambda[lambda_i]);
        if(Lambda[lambda_i]==0.05 || Lambda[lambda_i]==0.15 || Lambda[lambda_i]==0.25 || Lambda[lambda_i]==0.35 || Lambda[lambda_i]==0.45 || Lambda[lambda_i]==0.55 ||
	   Lambda[lambda_i]==0.75 || Lambda[lambda_i]==1.25 || Lambda[lambda_i]==1.75 || Lambda[lambda_i]==2.25 || Lambda[lambda_i]==2.75 || Lambda[lambda_i]==3.25 ||
	   Lambda[lambda_i]==3.75 || Lambda[lambda_i]==4.25 || Lambda[lambda_i]==4.75 || Lambda[lambda_i]==5.25 || Lambda[lambda_i]==5.75 || Lambda[lambda_i]==6.25 ||
	   Lambda[lambda_i]==6.75
			){
        sprintf(Lambda_char,"%.2f",Lambda[lambda_i]);
        }
        if( Lambda[lambda_i]== 0.125 || Lambda[lambda_i]==0.175 || Lambda[lambda_i]==0.225 || Lambda[lambda_i]==0.275 || Lambda[lambda_i]==0.325 || Lambda[lambda_i]==0.425 || Lambda[lambda_i]==0.475 || Lambda[lambda_i]==0.525 || Lambda[lambda_i]==0.575){
        sprintf(Lambda_char,"%.3f",Lambda[lambda_i]);
        }       
 
	//sprintf(Lambda_char,"%.1f",Lambda[lambda_i]);

	for(int dis_str_i=0;dis_str_i<Dis_Str_nos;dis_str_i++){
	  sprintf(DisStr_char,"%.1f", DisStr[dis_str_i]);

        if(DisStr[dis_str_i]==0.05 || DisStr[dis_str_i]==0.15 || DisStr[dis_str_i]==0.25 || DisStr[dis_str_i]==0.35 || DisStr[dis_str_i]==0.45 || DisStr[dis_str_i]==0.55 ||
           DisStr[dis_str_i]==0.75 || DisStr[dis_str_i]==1.25 || DisStr[dis_str_i]==1.75 || DisStr[dis_str_i]==2.25 || DisStr[dis_str_i]==2.75 || DisStr[dis_str_i]==3.25 ||
           DisStr[dis_str_i]==3.75 || DisStr[dis_str_i]==4.25 || DisStr[dis_str_i]==4.75 || DisStr[dis_str_i]==5.25 || DisStr[dis_str_i]==5.75 || DisStr[dis_str_i]==6.25 ||
           DisStr[dis_str_i]==6.75
                        ){
        sprintf(DisStr_char,"%.2f",DisStr[dis_str_i]);
        }

	Mat_2_doub Sq_real, Sq_imag, Sq2;
	Sq_real.resize(Length);Sq_imag.resize(Length);
	Sq2.resize(Length);

	for(int qx=0;qx<Length;qx++){
		Sq_real[qx].resize(Length);
                Sq_imag[qx].resize(Length);
		Sq2[qx].resize(Length);
	}


	int NumDisConf=60;
	for (int dis_seed=1;dis_seed<=NumDisConf;dis_seed++){	

	char dis_seed_char[50];
	sprintf(dis_seed_char,"%d",dis_seed);

	//Lx_Ly_12x12_PBC_Run2/NP_576/Seed_1/U_4.0/JbyU_0.25/lambda_3.0/dis_str_6.75/dis_seed_4/Sq.txt
	
       string fl_in_Sq = "Lx_Ly_" + string(Length_char) + "x" + string(Length_char) +  "_PBC_Run2/NP_" + string(N_total_char) + "/Seed_1/U_"+ string(U_char) +"/JbyU_0.25/lambda_"+string(Lambda_char)+ "/dis_str_" + string(DisStr_char) + "/dis_seed_" + string(dis_seed_char)  + "/Lq.txt" ;

	ifstream file_Sq(fl_in_Sq.c_str());
	
	if(!file_Sq.is_open()){
	cout<<fl_in_Sq<<" not present"<<endl;
	}
	
	else{
	int line_no=1;
	string tmp_str,temp_str;
	int temp_site_2;	
	double temp_val, temp_val_real, temp_val_imag;
	int val_x, val_y;	


	getline(file_Sq,tmp_str);
	
	for (int site_i=0;site_i<Length;site_i++){
		for (int site_j=0;site_j<Length;site_j++){
	file_Sq>>temp_val>>temp_val>>val_x>>val_y>>temp_val_real>>temp_val_imag;
	Sq_real[val_x][val_y] += temp_val_real;
	Sq_imag[val_x][val_y] += temp_val_imag;
	Sq2[val_x][val_y] += (temp_val_real*temp_val_real) + (temp_val_imag*temp_val_imag);
        

	line_no++;
		}
	getline(file_Sq,tmp_str);
	}
	
	
	//-----------------------------------------------------------------------------------------------//
	
	//----------------------------------------------------------------------------------------------//

	//---------------------------------------------------------------------------------------------//
	

	//----------------------------------------------------------------------------------------------//



	
	}	

							}//dis_seed
	



	double Sq2_mean, Sq_mean_real, Sq_mean_imag, Sq_mean;
	double variance, std_dev;
	string fl_out_Sq = "Lx_Ly_" + string(Length_char) + "x" + string(Length_char) +  "_PBC_Run2/NP_" + string(N_total_char) + "/Seed_1/U_"+ string(U_char) +"/JbyU_0.25/lambda_"+string(Lambda_char)+ "/dis_str_" + string(DisStr_char) + "/Lq_Avr.txt" ;
        ofstream file_out_Sq(fl_out_Sq.c_str());
	file_out_Sq<<"#qx   qy   qx_index    qy_index   S_avg(qx,qy).real   S_avg(qx,qy).imag   std.dev.Sq(qx,qy)"<<endl;
	for (int qx_i=0;qx_i<Length;qx_i++){
                for (int qy_i=0;qy_i<Length;qy_i++){
			Sq2_mean = (Sq2[qx_i][qy_i])/(1.0*NumDisConf);

			Sq_mean_real = (Sq_real[qx_i][qy_i])/(1.0*NumDisConf);
                        Sq_mean_imag = (Sq_imag[qx_i][qy_i])/(1.0*NumDisConf);

			Sq_mean = sqrt( (Sq_mean_real*Sq_mean_real)   + (Sq_mean_imag*Sq_mean_imag)   );
			variance = Sq2_mean - (Sq_mean*Sq_mean); 

			if(variance<-10e-6){
			cout<<"variance = "<<variance<<endl;
			assert(variance >= 0);
			}

			std_dev = sqrt(variance);			

			file_out_Sq<<(2.0*PI_*qx_i)/(1.0*Length)<<"\t"<<(2.0*PI_*qy_i)/(1.0*Length)<<"\t"<<qx_i<<"\t"<<qy_i<<"\t"<<Sq_mean_real<<"\t"<<Sq_mean_imag<<
				"\t"<< std_dev<<endl;
		
		}

		file_out_Sq<<endl;	
	}



	cout << "DONE U="<<string(U_char)<<", lambda = "<<string(Lambda_char)<<", Disorder strength = "<<string(DisStr_char)<<endl;

		
			}//dis_str_i	
		}// lambda_i


		}//U_i






	
return 0;
}


