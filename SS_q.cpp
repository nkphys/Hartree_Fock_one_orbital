#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <complex>
#include <vector>
#include <math.h>
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
	

//Target_Jz/L_8/N_32/Jz_48/lambda_SO_0.1/U_2.45/m_400_sweeps_9/CdagC_up_0_dn_1_vals.txt,n_up_0.txt
//Target_Jz/L_24/N_96/Jz_144/lambda_SO_0.5/U_2.45/m_500_s15/
//Target_Jz/L_32/N_128/Jz_192/
//L_16/N_56/Jz_84/
//L_16/N_64/Jz_96/
///L_16/N_72/Jz_108/
//L_16/N_80/Jz_120
//Target_Jz_degenerate/L_16/N_56/Jz_84/lambda_SO_1.0/U_2.0/m_500_s19/SzSzcorrs_full.txt
//L_32/N_112/Jz_168/
//Target_Jz_degenerate/L_32/N_128/Jz_192/lambda_SO_0.575/U_4.0/m_500_s19/SzSzcorrs_full.txt
//Target_Jz_degenerate/L_24/N_96/Jz_144/
//Target_Jz_degenerate/L_16/N_64/Jz_96/
//Target_Jz_degenerate/L_12/N_48/Jz_72/lambda_SO_0.575/U_4.0/m_500_s19/SzSzcorrs_full.txt
//Target_Jz_degenerate/L_8/N_32/Jz_48/lambda_SO_0.575/
//Target_Jz_degenerate/L_20/N_80/Jz_120/lambda_SO_0.575/U_4.0/m_500_s19/
//Target_Jz_degenerate/L_28/N_112/Jz_168/
//Target_Jz_degenerate_onesite_PBC/L_12/N_42/Jz_63/lambda_SO_0.3/U_10.0/m_500_s19/

	//L24_Jz/lambda0.4/SzSz.txt

	double tmp_n,tmp_tm;
	string tmp_line;
	int tmp_site;

	int Length =24;
	char Length_char[50];
	sprintf(Length_char,"%d",Length);
		
	
	complex<double> iota;iota.real(0);iota.imag(1);
	

	Mat_2_Complex_doub   SzSz;
	Mat_2_Complex_doub   SpSm;
	Mat_2_Complex_doub   SmSp;
	Mat_2_Complex_doub   SS;
	Mat_1_Complex_doub   Sz;
	Mat_1_Complex_doub SS_d, SzSz_d, SS_trnvrs_d;
	Mat_1_Complex_doub SSq, SzSzq, SS_trnvrsq;
	Mat_1_doub q_vals;		

	Sz.resize(Length);
	SS_d.resize(Length);
	SzSz_d.resize(Length);
	SS_trnvrs_d.resize(Length);

	SS.resize(Length);
	SzSz.resize(Length);
	SpSm.resize(Length);
	SmSp.resize(Length);
	
	for(int i=0;i<Length;i++){
	SS[i].resize(Length);
        SzSz[i].resize(Length);
        SpSm[i].resize(Length);
        SmSp[i].resize(Length);
	}  
	
	double n_i=0;
	double L_inv;
	L_inv=(1.0/Length);

        double qbypi=(n_i)*(1.0/(Length+1.0));

        //FOR OBC
        /*
	n_i=1.0;
        while(qbypi<=1.0){
        qbypi=(n_i)*(1.0/(Length+1.0));
        q_vals.push_back(qbypi);
	  if(n_i==((Length*1.0)/2.0)){
        q_vals.push_back(0.5);
        }
        n_i=n_i+1;
        }
	*/



	n_i=0.0;
        while(qbypi<=1.0){
         qbypi = ((n_i*1.0))/((Length+1)*1.0);
		//qbypi=(2.0*n_i)*(1.0/(Length));
        q_vals.push_back(qbypi);
        n_i=n_i+1;
        }

	SSq.resize(q_vals.size());
	SzSzq.resize(q_vals.size());
	SS_trnvrsq.resize(q_vals.size());

	double U[1] ={40.0};//{0.049, 0.098, 0.245, 0.98, 1.96, 2.45, 4.9, 12.25, 24.5, 49, 98};
	//double U[1] = {2.45};
	int U_nos=1;
	char U_char[50];
	
//L24_Jz/lambda0.4/SzSz.txt

	double Lambda[12]= {0.0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
	
	//{0.21, 0.225, 0.25, 0.275, 0.3, 0.35, 0.45};//{0.1, 0.15, 0.175, 0.18, 0.185, 0.19, 0.2, 0.205}; //{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };//, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0};//,0.3,0.4,0.5};//{0.0, 0.1, 0.2 ,0.3 ,0.4 ,0.5, 0.7, 0.8, 1.0, 2.0, 3.0};
	int Lambda_nos=12;
	//double Lambda[1]= {0.8};
	char Lambda_char[50];
	
	ofstream file_check("file_check.txt");
	
	for(int lambda_i=0;lambda_i<Lambda_nos;lambda_i++){
	
	//0.05 0.1 0.14 0.18 0.2 0.21 0.22 0.23 0.24 0.25 0.27 0.28 0.3 0.32 0.32 0.38 
        if(Lambda[lambda_i]==0.02 || Lambda[lambda_i]==0.05 || Lambda[lambda_i]==0.15 ||
        Lambda[lambda_i]==0.18 || Lambda[lambda_i]==0.19 || Lambda[lambda_i]==0.21 ||
          Lambda[lambda_i]==0.25 || Lambda[lambda_i]==0.35 || Lambda[lambda_i]==0.45 ||
          Lambda[lambda_i]==0.65 || Lambda[lambda_i]==0.55  || Lambda[lambda_i]==0.75 ||
          Lambda[lambda_i]==0.36 || Lambda[lambda_i]==0.38  || Lambda[lambda_i]==0.39 ||
          Lambda[lambda_i]==0.41 || Lambda[lambda_i]==0.45 || Lambda[lambda_i]==0.43 || 
	  Lambda[lambda_i]==0.55 || Lambda[lambda_i]==0.85 ||
          Lambda[lambda_i]==0.95 || Lambda[lambda_i]==0.14 || Lambda[lambda_i]==0.22
	  || Lambda[lambda_i]==0.23 || Lambda[lambda_i]==0.24 || Lambda[lambda_i]==0.27
	  || Lambda[lambda_i]==0.32 
          ){
        sprintf(Lambda_char,"%.2f",Lambda[lambda_i]);
        }
        else if(Lambda[lambda_i]==0.125 || Lambda[lambda_i]==0.175 || Lambda[lambda_i]==0.225
                || Lambda[lambda_i]==0.185 || Lambda[lambda_i]==0.205 || Lambda[lambda_i]==0.215
                        || Lambda[lambda_i]==0.275 || Lambda[lambda_i]==0.325 || Lambda[lambda_i]==0.425 || Lambda[lambda_i]==0.475 || Lambda[lambda_i]==0.525 || Lambda[lambda_i]==0.575
        || Lambda[lambda_i]==0.625 || Lambda[lambda_i]==0.725  ||
       Lambda[lambda_i]==0.825  ||      Lambda[lambda_i]==0.875 || Lambda[lambda_i]==0.925){
        sprintf(Lambda_char,"%.3f",Lambda[lambda_i]);
        }
        else{
        sprintf(Lambda_char,"%.1f",Lambda[lambda_i]);
        }


        cout<<Lambda_char<<"   "<<Lambda[lambda_i]<<endl;
	complex<double> S_2, Sz_2, S_trnvrs_2;

	string fl_out_S_2 = "Data/S2_new_L_"+string(Length_char)+"_lambda_" + string(Lambda_char)
        + ".txt";
        ofstream file_out_S_2(fl_out_S_2.c_str());
	file_out_S_2<<"#U	S_2	Sz_2	S_transverse_2(i.e Sx2+Sy2)"<<endl;

	for (int U_i=0;U_i<U_nos;U_i++){
	//cout<<"here 9"<<endl;

	//S_2.resize(U_nos);
	//Sz_2.resize(U_nos);
	//S_trnvrs_2.resize(U_nos);
/*	
	if(U[i]==0.0 || U[i]==4.9){
		sprintf(u_val,"%1.1f",U[i]);}
	
	else if((U[i]<=1.225) && (U[i]>=0.245)){	
	sprintf(u_val,"%1.3f",U[i]);}
	
	else if(U[i]==2.45 ||U[i]==7.35){	
	sprintf(u_val,"%1.2f",U[i]);}
	 
	else if(( (U[i]==14.7) || (U[i]==24.5) ) || (U[i]==49.0)){
		sprintf(u_val,"%2.1f",U[i]);} 
		
	else if(U[i]==12.25 || U[i]==36.75){
	sprintf(u_val,"%2.2f",U[i]);}
	
	else {	
	sprintf(u_val,"%1.4f",0.0245);}
	
*/

	if(U[U_i]<=0.245 || U[U_i]==2.165){
	sprintf(U_char,"%.3f",U[U_i]);}
	if( U[U_i]==1.96 || U[U_i]==2.45 || U[U_i]==12.25 || U[U_i]==0.98 || U[U_i]==3.25 || U[U_i]==21.65){
	sprintf(U_char,"%.2f",U[U_i]);}
	if( U[U_i]==4.9 ||U[U_i]==24.5 ) {
	sprintf(U_char,"%.1f",U[U_i]);}
	if (U[U_i]==49 || U[U_i]==98) {
	sprintf(U_char,"%.0f",U[U_i]);
	}	
	sprintf(U_char,"%.1f",U[U_i]);	
	
	cout<<"Started lambda = "<<Lambda_char << ", U = "<<U_char<<endl;
	


	string fl_out_SSq = "Data/SS_q_L_"+string(Length_char)+"_lambda_" + string(Lambda_char)
        + ".txt" ;
        ofstream file_out_SSq(fl_out_SSq.c_str());

        //file_out_SS_real<<"# this file have SScorr for for given lambda,U produced by 'SS_corr.cpp' "<<endl;
        //file_out_SS_imag<<"# this file have SScorr for given lambda,U produced by 'SS_corr.cpp' "<<endl;
	file_out_SSq<<"# this file have SS(q) for given lambda,U produced by 'SS_corr.cpp' "<<endl;
	file_out_SSq<<"#q         SS(q)"<<endl;

        string str_vars ="#distance          SS";
       // file_out_SS_real<<str_vars<<endl;
        //file_out_SS_imag<<str_vars<<endl;

	//cout<<U_i<<"   "<<U[U_i]<<"  "<<U_char<<endl;
	
	//Target_Sz/L_8/N_32/Nup_16/U_12.25/n_vals.txt
	//Target_Jz/L_8/N_32/Jz_48/lambda_SO_0.1/U_2.45/m_400_sweeps_9/CdagC_up_0_dn_1_vals.txt,n_up_0.txt
	
	string fl_in_SzSz;
	string fl_in_SpSm;
	string fl_in_Sz;

	//L24_Jz/lambda0.4/SzSz.txt
        fl_in_SzSz = "L" + string(Length_char) + "_Jz/" + "lambda"+string(Lambda_char)+"/SzSz.txt";
	fl_in_SpSm = "L" + string(Length_char) + "_Jz/" + "lambda"+string(Lambda_char)+"/SpSm.txt";

//fl_in_Sz = Target + "/L_" + string(Length_char) + "/N_" + string(N_total_char) + "/Jz_" + string(symmetry_sec_char) + "/lambda_SO_"+string(Lambda_char)+"/U_"+ string(U_char) +"/m_500_s19/Sz_T_values.txt";

//	cout<<fl_in_Sz<<endl;
		
	ifstream file_SzSz(fl_in_SzSz.c_str());
	ifstream file_SpSm(fl_in_SpSm.c_str());
	//ifstream file_Sz(fl_in_Sz.c_str());
	
		
	cout<< fl_in_SzSz<<endl;
	
	//cout<< "here 2"<<endl;
	int line_no=1;
	string tmp_str,temp_str;
	int temp_site_2;	



/*
           inputfile_up>>O1_str;
          O1_str_L=O1_str.length()-1;
          pos = O1_str.find(",");
        O1_real_str=O1_str.substr(1,pos-1);
        O1_imag_str=O1_str.substr(pos+1,O1_str_L-1-pos);

        stringstream Rstream(O1_real_str);
        stringstream Istream(O1_imag_str);
        Rstream>>O1_real[i][j];
        Istream>>O1_imag[i][j];
 */

//cout<< "here 1.5"<<endl;


	file_SzSz.clear();
        file_SzSz.seekg(0, ios::beg);

        string temp_string="123ddf";

	double tmp_double;

/*
        while(temp_string != "<gs|:Sz_T.txt;:Sz_T.txt|gs>"){
        getline(file_SzSz,temp_string);
        //cout<<temp_string<<endl;
        }
        getline(file_SzSz,temp_string);	
*/
	//----------------------------------------Read SzSz----------------------------------------------------//
	//cout<<"here 2"<<endl;
	for (int site_i=0;site_i<Length;site_i++){
		for (int site_j=0;site_j<Length;site_j++){
	file_SzSz>>tmp_double;
	//cout<<tmp_str<<endl;
	line_no++;
	//cout<<"line "<<line_no<<endl;	
	SzSz[site_i][site_j]=tmp_double;
	//cout<<SzSz[site_i][site_j]<<"  ";
		}
	//cout<<endl;
	}
	//cout<<"here 3"<<endl;
	//-----------------------------------------------------------------------------------------------//
	
//cout<< "here 2"<<endl;


		

	file_SpSm.clear();
        file_SpSm.seekg(0, ios::beg);

        temp_string="123ddf";
/*
        while(temp_string != "<gs|:Splus_T.txt;:Sminus_T.txt|gs>"){
        getline(file_SpSm,temp_string);
        //cout<<temp_string<<endl;
        }
        getline(file_SpSm,temp_string);
*/
	//-------------------------------------Read SpSm-----------------------------------------------------//
        for (int site_i=0;site_i<Length;site_i++){
                for (int site_j=0;site_j<Length;site_j++){
        file_SpSm>>tmp_double;
        line_no++;

        SpSm[site_i][site_j]=tmp_double;
		}
	}
	//----------------------------------------------------------------------------------------------//


	//----------------------------Read Sz ----------------------------------------------------//


//site <gs|:Sz_T_transpose.txt|gs> time

	// ifstream file_Sz;
        // file_Sz.open(fl_in_Sz.c_str());
   
       // file_Sz.clear();
       // file_Sz.seekg(0, ios::beg);

        temp_string="123ddf";	

//	while(temp_string != "site <gs|:Sz_T_transpose.txt|gs> time"){
       // getline(file_Sz,temp_string);
        //cout<<temp_string<<endl;
  //      }
//	getline(file_Sz,temp_string);


//        line_no=1;
  //      while(line_no<91){
    //    getline(file_Sz,temp_str);
      //  line_no++;}


  //      for (int site=0;site<Length;site++){
    //    file_Sz>>tmp_site>>tmp_str>>temp_site_2;
      //  line_no++;
        //Sz[tmp_site]=reading_pair(tmp_str);
	
	//cout<<tmp_site<<"  "<<Sz[tmp_site]<<endl;
          //      }

	//-------------------------------------------------------------------------------------//


	//---------------------------------Read SmSp------------------------------------------------//
	for (int site_i=0;site_i<Length;site_i++){
        for (int site_j=0;site_j<Length;site_j++){
        //file_SmSp>>tmp_str;
        //line_no++;
	
	if(site_i==site_j){
        SmSp[site_i][site_j]= SpSm[site_i][site_j];}// + 2.0*Sz[site_i];}
	else{
	SmSp[site_i][site_j]= conj(SpSm[site_i][site_j]);
	}

                }
        }
	//---------------------------------------------------------------------------------------------//
	

	//---------------Construct S.S[i][j]-----------------------------------------------------------------//
	for (int site_i=0;site_i<Length;site_i++){
        for (int site_j=0;site_j<Length;site_j++){

        SS[site_i][site_j]=SzSz[site_i][site_j] + 0.5*(SpSm[site_i][site_j] + SmSp[site_i][site_j]);

	//cout<<SS[site_i][site_j]<<"   ";
                }
	//cout<<endl;
        }
	//----------------------------------------------------------------------------------------------//




	
	for (int site_i=0;site_i<Length;site_i++){
		SS_d[site_i]=(0,0);
                for (int site_j=0;site_j<site_i;site_j++){
       
        SS[site_i][site_j]=SS[site_j][site_i];
	SzSz[site_i][site_j]=SzSz[site_j][site_i];
	SpSm[site_i][site_j]=SpSm[site_j][site_i];
	SmSp[site_i][site_j]=SmSp[site_j][site_i];
                }
        }

	
	S_2=(0,0); Sz_2=(0,0); S_trnvrs_2=(0,0);
	for(int Dis=0;Dis<Length;Dis++){
		for(int site_i=0;site_i<Length;site_i++){
                for(int site_j=0;site_j<=site_i;site_j++){
		if(abs(site_j-site_i)==Dis){
		SS_d[Dis]=SS_d[Dis]+SS[site_i][site_j];
		}
		if(abs(site_j-site_i)==Dis && Dis==0){
                S_2 = S_2 + SS[site_i][site_j];
		Sz_2 = Sz_2 + SzSz[site_i][site_j];
		S_trnvrs_2 = S_trnvrs_2 + 0.5*(SpSm[site_i][site_j]+SmSp[site_i][site_j]);
                }

		}}
		//cout<<"here 4"<<endl;
		double x=1.0/(Length-Dis);
		//cout<<x<<endl; 
		SS_d[Dis].real((SS_d[Dis].real())*x);
		//cout<<"here 5"<<endl;
		SS_d[Dis].imag((SS_d[Dis].imag())*x);									
		//cout<<"here 6"<<endl;


        //file_out_SS_real<<Dis<<"   "<<SS_d[Dis].real()<<endl;
	//cout<<"here 7"<<endl;
        //file_out_SS_imag<<Dis<<"   "<<SS_d[Dis].imag()<<endl;
	//cout<<"here 8"<<endl;
		}

		S_2=L_inv*S_2;
                Sz_2=L_inv*Sz_2;
                S_trnvrs_2=L_inv*S_trnvrs_2;
//	file_out_S_2<<U[U_i]<<"   "<<S_2.real()<<"   "<<Sz_2.real()<<"   "<<S_trnvrs_2.real()<<endl;

	for(int q_i=0;q_i<q_vals.size();q_i++){
		SSq[q_i].real(0); SSq[q_i].imag(0);
	    for(int site_i=0;site_i<Length;site_i++){
            for(int site_j=0;site_j<Length;site_j++){
	
		//SSq[q_i].real( SSq[q_i].real() + sin(q_vals[q_i]*PI_*(site_j+1))*sin(q_vals[q_i]*PI_*(site_i+1))*SS[site_i][site_j].real()*L_inv);
		SSq[q_i].real( SSq[q_i].real() + cos(q_vals[q_i]*PI_*(site_j - site_i))*SS[site_i][site_j].real()*L_inv);

	}}
	file_out_SSq<<q_vals[q_i]<<"    "<<SSq[q_i].real()<<endl;
	}

	
	
			
	
	cout<<"Completed lambda = "<<Lambda_char << ", U = "<<U_char<<endl;	
		
							}// U_i


		}//lambda_i






	
return 0;
}


