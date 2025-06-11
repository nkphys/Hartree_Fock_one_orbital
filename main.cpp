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
using namespace std;
#include "Matrix.h"
#include "ParametersEngine.h"
#include "Coordinates.h"
#include "MFParams.h"
#include "Hamiltonian.h"
#include "Observables.h"
#include "SelfConsistencyEngine.h"

//Model files
#include "src/Models/DiceLattice/Parameters_DL.h"
#include "src/Models/DiceLattice/Coordinates_DL.h"
#include "src/Models/DiceLattice/Connections_DL.h"
#include "src/Models/DiceLattice/Observables_DL.h"
#include "src/Models/DiceLattice/Kspace_calculation_DL.h"

#include "src/Models/LiebLattice/Parameters_LL.h"
#include "src/Models/LiebLattice/Coordinates_LL.h"
#include "src/Models/LiebLattice/Connections_LL.h"
#include "src/Models/LiebLattice/Observables_LL.h"
#include "src/Models/LiebLattice/Kspace_calculation_LL.h"

#include "src/Models/SquareLattice/Parameters_SQL.h"
#include "src/Models/SquareLattice/Coordinates_SQL.h"
#include "src/Models/SquareLattice/Connections_SQL.h"
#include "src/Models/SquareLattice/Observables_SQL.h"

#include "src/Models/TriangularLattice_TBG/Parameters_TL.h"
#include "src/Models/TriangularLattice_TBG/Coordinates_TL.h"
#include "src/Models/TriangularLattice_TBG/Connections_TL.h"
#include "src/Models/TriangularLattice_TBG/Observables_TL.h"
#include "src/Models/TriangularLattice_TBG/Kspace_calculation_TL.h"


#include "src/Models/HoneycombLattice_TBG/Parameters_HC.h"
#include "src/Models/HoneycombLattice_TBG/Coordinates_HC.h"
#include "src/Models/HoneycombLattice_TBG/Connections_HC.h"
#include "src/Models/HoneycombLattice_TBG/Observables_HC.h"
#include "src/Models/HoneycombLattice_TBG/Kspace_calculation_HC.h"

#include "src/Models/HoneycombLattice_MultiOrb/Parameters_HC_MO.h"
#include "src/Models/HoneycombLattice_MultiOrb/Coordinates_HC_MO.h"
#include "src/Models/HoneycombLattice_MultiOrb/Connections_HC_MO.h"
#include "src/Models/HoneycombLattice_MultiOrb/Observables_HC_MO.h"
#include "src/Models/HoneycombLattice_MultiOrb/Kspace_calculation_HC_MO.h"


#include "src/Models/Generic2dLattice/Parameters_G2dLattice.h"
#include "src/Models/Generic2dLattice/Coordinates_G2dLattice.h"
#include "src/Models/Generic2dLattice/Connections_G2dLattice.h"
#include "src/Models/Generic2dLattice/Observables_G2dLattice.h"
#include "src/Models/Generic2dLattice/Kspace_calculation_G2dLattice.h"

#include "src/Models/Generic2dLatticeNew/Parameters_G2dLatticeNew.h"
#include "src/Models/Generic2dLatticeNew/Coordinates_G2dLatticeNew.h"
#include "src/Models/Generic2dLatticeNew/Connections_G2dLatticeNew.h"
#include "src/Models/Generic2dLatticeNew/Observables_G2dLatticeNew.h"
#include "src/Models/Generic2dLatticeNew/Kspace_calculation_G2dLatticeNew.h"


#include "random"

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char *argv[]) {

    string ex_string_original =argv[0];

    cout<<"'"<<ex_string_original<<"'"<<endl;
    string ex_string;
    ex_string = ex_string_original.substr(ex_string_original.length()-5);
    //ex_string=ex_string_original.substr (2);
    cout<<"'"<<ex_string<<"'"<<endl;




    if(ex_string=="tency"){ //k_space_SelfConsistency


        string ModelType = argv[1];
        string model_inputfile = argv[2];

        if (argc<3) { throw std::invalid_argument("USE:: executable inputfile"); }

        if(ModelType=="DiceLattice"){

            Parameters_DL Parameters_DL_;
            Parameters_DL_.Initialize(model_inputfile);

            Coordinates_DL Coordinates_DL_(Parameters_DL_.lx, Parameters_DL_.ly, Parameters_DL_.n_orbs);

            mt19937_64 Generator_(Parameters_DL_.RandomSeed);
            Kspace_calculation_DL Kspace_calculation_DL_(Parameters_DL_, Coordinates_DL_,Generator_);

            Kspace_calculation_DL_.SelfConsistency();

        }

        if(ModelType=="LiebLattice"){

            Parameters_LL Parameters_LL_;
            Parameters_LL_.Initialize(model_inputfile);

            Coordinates_LL Coordinates_LL_(Parameters_LL_.lx, Parameters_LL_.ly, Parameters_LL_.n_orbs);

            mt19937_64 Generator_(Parameters_LL_.RandomSeed);
            Kspace_calculation_LL Kspace_calculation_LL_(Parameters_LL_, Coordinates_LL_,Generator_);

            Kspace_calculation_LL_.SelfConsistency();

        }



        if(ModelType=="TriangularLattice"){

            Parameters_TL Parameters_TL_;
            Parameters_TL_.Initialize(model_inputfile);

            Coordinates_TL Coordinates_TL_(Parameters_TL_.lx, Parameters_TL_.ly, Parameters_TL_.n_orbs);
            Connections_TL Connections_TL_(Parameters_TL_, Coordinates_TL_);
            Connections_TL_.Print_Hopping();                                       //::DONE
            Connections_TL_.InteractionsCreate();
            Connections_TL_.Print_LongRangeInt();
            Connections_TL_.Print_Spin_resolved_OnsiteE();

            Coordinates_TL Coordinates_TL_UC_(Parameters_TL_.lx/Parameters_TL_.UnitCellSize_x, Parameters_TL_.ly/Parameters_TL_.UnitCellSize_y, 1);


            mt19937_64 Generator_(Parameters_TL_.RandomSeed);
            Kspace_calculation_TL Kspace_calculation_TL_(Parameters_TL_, Coordinates_TL_UC_, Connections_TL_, Generator_);

            Kspace_calculation_TL_.SelfConsistency();
            //Observables_TL_.Calculate_Akxw_ribbon();

        }


        if(ModelType=="HoneycombLattice"){

            Parameters_HC Parameters_HC_;
            Parameters_HC_.Initialize(model_inputfile);

            Coordinates_HC Coordinates_HC_(Parameters_HC_.lx, Parameters_HC_.ly, Parameters_HC_.n_orbs);
            Connections_HC Connections_HC_(Parameters_HC_, Coordinates_HC_);
            Connections_HC_.Print_Hopping();                                       //::DONE
            Connections_HC_.InteractionsCreate();
            Connections_HC_.Print_LongRangeInt();
            Connections_HC_.Print_Spin_resolved_OnsiteE();
            Connections_HC_.Interactions_Sorting();

            Coordinates_HC Coordinates_HC_UC_(Parameters_HC_.lx/Parameters_HC_.UnitCellSize_x, Parameters_HC_.ly/Parameters_HC_.UnitCellSize_y, 1);


            mt19937_64 Generator_(Parameters_HC_.RandomSeed);
            Kspace_calculation_HC Kspace_calculation_HC_(Parameters_HC_, Coordinates_HC_UC_, Connections_HC_, Generator_);
            Kspace_calculation_HC_.SelfConsistency();

        }

         if(ModelType=="HoneycombLattice_MO"){

            Parameters_HC_MO Parameters_HC_;
            Parameters_HC_.Initialize(model_inputfile);

            Coordinates_HC_MO Coordinates_HC_(Parameters_HC_.lx, Parameters_HC_.ly, Parameters_HC_.n_atoms, Parameters_HC_.n_orbs);
            Connections_HC_MO Connections_HC_(Parameters_HC_, Coordinates_HC_);
            Connections_HC_.Print_Hopping();                                       //::DONE
            Connections_HC_.InteractionsCreate();
            Connections_HC_.Print_LongRangeInt();
            //Connections_HC_.Print_Spin_resolved_OnsiteE();
            Connections_HC_.Interactions_Sorting();

            Coordinates_HC_MO Coordinates_HC_UC_(Parameters_HC_.lx/Parameters_HC_.UnitCellSize_x, Parameters_HC_.ly/Parameters_HC_.UnitCellSize_y, 1, 1);


            mt19937_64 Generator_(Parameters_HC_.RandomSeed);
            Kspace_calculation_HC_MO Kspace_calculation_HC_(Parameters_HC_, Coordinates_HC_UC_, Connections_HC_, Generator_);
            Kspace_calculation_HC_.SelfConsistency();

        }

         if(ModelType=="Generic2dLattice"){

             Parameters_G2dLattice Parameters_G2d_;
             Parameters_G2d_.Initialize(model_inputfile);

             Coordinates_G2dLattice  Coordinates_G2d_(Parameters_G2d_.lx, Parameters_G2d_.ly, Parameters_G2d_.n_atoms, Parameters_G2d_.n_orbs);
             Connections_G2dLattice  Connections_G2d_(Parameters_G2d_, Coordinates_G2d_);
             Connections_G2d_.Print_Hopping();                                       //::DONE
             // Connections_HC_.InteractionsCreate();
             // Connections_HC_.Print_LongRangeInt();
             //Connections_HC_.Print_Spin_resolved_OnsiteE();
             // Connections_HC_.Interactions_Sorting();

             Coordinates_G2dLattice Coordinates_G2d_UC_(Parameters_G2d_.lx/Parameters_G2d_.UnitCellSize_x, Parameters_G2d_.ly/Parameters_G2d_.UnitCellSize_y, 1, 1);


             mt19937_64 Generator_(Parameters_G2d_.RandomSeed);
             Kspace_calculation_G2dLattice Kspace_calculation_G2d_(Parameters_G2d_, Coordinates_G2d_UC_, Connections_G2d_, Generator_);
             Kspace_calculation_G2d_.SelfConsistency();
             // Kspace_calculation_G2d_.Get_Bare_Susceptibility();
             // Kspace_calculation_G2d_.Get_RPA_Susceptibility();

         }


         if(ModelType=="Generic2dLatticeNew"){

             Parameters_G2dLatticeNew Parameters_G2d_;
             Parameters_G2d_.Initialize(model_inputfile);


#ifdef _OPENMP
             int N_p = omp_get_max_threads();
             omp_set_num_threads(Parameters_G2d_.NProcessors);

             cout<<"Max threads which can be used parallely = "<<N_p<<endl;
             cout<<"No. of threads used parallely = "<<Parameters_G2d_.NProcessors<<endl;
#endif


             Coordinates_G2dLatticeNew  Coordinates_G2d_(Parameters_G2d_.lx, Parameters_G2d_.ly, Parameters_G2d_.n_atoms, Parameters_G2d_.n_orbs);
             Connections_G2dLatticeNew  Connections_G2d_(Parameters_G2d_, Coordinates_G2d_);
             Connections_G2d_.Print_Hopping();                                       //::DONE
             // Connections_HC_.InteractionsCreate();
             // Connections_HC_.Print_LongRangeInt();
             //Connections_HC_.Print_Spin_resolved_OnsiteE();
             // Connections_HC_.Interactions_Sorting();

             Coordinates_G2dLatticeNew Coordinates_G2d_UC_(Connections_G2d_.lx_cells, Connections_G2d_.ly_cells, 1, 1);


              mt19937_64 Generator_(Parameters_G2d_.RandomSeed);
              Kspace_calculation_G2dLatticeNew Kspace_calculation_G2d_(Parameters_G2d_, Coordinates_G2d_UC_, Connections_G2d_, Generator_);
              Kspace_calculation_G2d_.SelfConsistency();
              if(Parameters_G2d_.Calculate_Susc){
                  // Kspace_calculation_G2d_.Get_Bare_Susceptibility_New();
                  // if(Kspace_calculation_G2d_.NSites_in_MUC>0){
                  //     Kspace_calculation_G2d_.Get_Bare_Susceptibility_in_NBZ();
                  // }
                  Kspace_calculation_G2d_.Calculate_InteractionKernel();
                  Kspace_calculation_G2d_.Get_RPA_Susceptibility_Matrix();
                 // Kspace_calculation_G2d_.Get_Bare_Susceptibility_in_NBZ();
                  Kspace_calculation_G2d_.Get_RPA_Susceptibility_in_NBZ();

                  //Kspace_calculation_G2d_.Get_Bare_Susceptibility();
              }
             // Kspace_calculation_G2d_.Get_RPA_Susceptibility();

         }

    }




    else if(ex_string=="tions"){

        cout <<"Creating connections for the given model"<<endl;
        string ModelType = argv[1];
        if (argc<3) { throw std::invalid_argument("USE:: executable inputfile"); }
        string inputfile = argv[2];

        if(ModelType=="LiebLattice"){

            Parameters_LL Parameters_LL_;
            Parameters_LL_.Initialize(inputfile);

            Coordinates_LL Coordinates_LL_(Parameters_LL_.lx, Parameters_LL_.ly, Parameters_LL_.n_orbs);

            Connections_LL Connections_LL_(Parameters_LL_, Coordinates_LL_);
            Connections_LL_.Print_Hopping();                                       //::DONE
            Connections_LL_.Print_LongRangeInt();
            Connections_LL_.Print_Spin_resolved_OnsiteE();

        }
        if(ModelType=="DiceLattice"){

            Parameters_DL Parameters_DL_;
            Parameters_DL_.Initialize(inputfile);

            Coordinates_DL Coordinates_DL_(Parameters_DL_.lx, Parameters_DL_.ly, Parameters_DL_.n_orbs);

            Connections_DL Connections_DL_(Parameters_DL_, Coordinates_DL_);
            Connections_DL_.Print_Hopping();                                       //::DONE
            Connections_DL_.Print_LongRangeInt();
            Connections_DL_.Print_Spin_resolved_OnsiteE();

        }
        if(ModelType=="SquareLattice"){

            Parameters_SQL Parameters_SQL_;
            Parameters_SQL_.Initialize(inputfile);

            Coordinates_SQL Coordinates_SQL_(Parameters_SQL_.lx, Parameters_SQL_.ly, Parameters_SQL_.n_orbs);

            Connections_SQL Connections_SQL_(Parameters_SQL_, Coordinates_SQL_);
            Connections_SQL_.Print_Hopping();                                       //::DONE
            Connections_SQL_.Print_LongRangeInt();
            Connections_SQL_.Print_Spin_resolved_OnsiteE();

        }
        if(ModelType=="TriangularLattice"){

            Parameters_TL Parameters_TL_;
            Parameters_TL_.Initialize(inputfile);

            Coordinates_TL Coordinates_TL_(Parameters_TL_.lx, Parameters_TL_.ly, Parameters_TL_.n_orbs);

            Connections_TL Connections_TL_(Parameters_TL_, Coordinates_TL_);
            Connections_TL_.Print_Hopping();                                       //::DONE
            Connections_TL_.Print_Hopping2();                                       //::DONE
            Connections_TL_.Print_Hopping3();
            Connections_TL_.Print_Hopping_Snake3();
            Connections_TL_.Print_Hopping_YCn();
            Connections_TL_.Print_RingExchange();
            Connections_TL_.Print_RingExchange_YCn();

            Connections_TL_.InteractionsCreate();
            Connections_TL_.Print_LongRangeInt_SNAKE3();
            Connections_TL_.Print_LongRangeInt();
            Connections_TL_.Print_Spin_resolved_OnsiteE();

            Connections_TL_.Print_Ansatz_LocalDen_CDW();

        }

        if(ModelType=="HoneycombLattice"){

            Parameters_HC Parameters_HC_;
            Parameters_HC_.Initialize(inputfile);

            Coordinates_HC Coordinates_HC_(Parameters_HC_.lx, Parameters_HC_.ly, Parameters_HC_.n_orbs);

            Connections_HC Connections_HC_(Parameters_HC_, Coordinates_HC_);
            Connections_HC_.Print_Hopping();                                       //::DONE

            Connections_HC_.InteractionsCreate();
            Connections_HC_.Print_LongRangeInt();
            Connections_HC_.Print_Spin_resolved_OnsiteE();

        }
    }
    else if(ex_string=="serve"){

        string ModelType = argv[1];
        string inputfile = argv[2];
        string model_inputfile = argv[3];

        if (argc<4) { throw std::invalid_argument("USE:: executable inputfile"); }

        Parameters Parameters_;
        Parameters_.Initialize(inputfile);

        Coordinates Coordinates_(Parameters_.ns);

        mt19937_64 Generator_(Parameters_.RandomSeed); //for random fields
        mt19937_64 Generator2_(Parameters_.RandomDisorderSeed); //for random disorder

        MFParams MFParams_(Parameters_,Coordinates_,Generator_, Generator2_);
        Hamiltonian Hamiltonian_(Parameters_,Coordinates_,MFParams_);

        double initial_mu_guess;
        int n_states_occupied_zeroT;
        Parameters_.Dflag='V'; // flag to calculate Eigenspectrum


        Hamiltonian_.InteractionsCreate();
        Hamiltonian_.Diagonalize(Parameters_.Dflag);

        //calculating mu
        n_states_occupied_zeroT=Parameters_.Total_Particles;
        initial_mu_guess=0.5*(Hamiltonian_.eigs_[n_states_occupied_zeroT-1] + Hamiltonian_.eigs_[n_states_occupied_zeroT]);
        Parameters_.mus=Hamiltonian_.chemicalpotential(initial_mu_guess,Parameters_.Total_Particles);


        if(ModelType=="LiebLattice"){

            Parameters_LL Parameters_LL_;
            Parameters_LL_.Initialize(model_inputfile);

            Coordinates_LL Coordinates_LL_(Parameters_LL_.lx, Parameters_LL_.ly, Parameters_LL_.n_orbs);
            Connections_LL Connections_LL_(Parameters_LL_, Coordinates_LL_);

            Observables_LL Observables_LL_(Parameters_LL_, Coordinates_LL_, Connections_LL_ );

            Parameters_LL_.beta=Parameters_.beta;
            Parameters_LL_.mus=Parameters_.mus;
            Parameters_LL_.eta=Parameters_.eta_dos;
            Parameters_LL_.domega=Parameters_.dw_dos;
            Observables_LL_.Ham_ = Hamiltonian_.Ham_;
            Observables_LL_.eigs_ = Hamiltonian_.eigs_;


            Observables_LL_.Calculate_Akw();
            //        Observables_LL_.Calculate_OrbResolved_Nw();
            Observables_LL_.Calculate_Nw();

            Observables_LL_.Create_Current_Oprs();
            Observables_LL_.Hall_conductance();

            //Create and call:
            //quantum spin-spin corrs
            //<s_i^2>
            //<n_iup n_idn>
            //spin resolved local density

        }

        if(ModelType=="DiceLattice"){

            Parameters_DL Parameters_DL_;
            Parameters_DL_.Initialize(model_inputfile);

            Coordinates_DL Coordinates_DL_(Parameters_DL_.lx, Parameters_DL_.ly, Parameters_DL_.n_orbs);
            Connections_DL Connections_DL_(Parameters_DL_, Coordinates_DL_);

            Observables_DL Observables_DL_(Parameters_DL_, Coordinates_DL_, Connections_DL_ );

            Parameters_DL_.beta=Parameters_.beta;
            Parameters_DL_.mus=Parameters_.mus;
            Parameters_DL_.eta=Parameters_.eta_dos;
            Parameters_DL_.domega=Parameters_.dw_dos;
            Observables_DL_.Ham_ = Hamiltonian_.Ham_;
            Observables_DL_.eigs_ = Hamiltonian_.eigs_;


            //        Observables_DL_.Calculate_Akw();
            //        Observables_DL_.Calculate_OrbResolved_Nw();
            //        Observables_DL_.Calculate_Nw();


            Observables_DL_.Create_Current_Oprs();
            Observables_DL_.Hall_conductance();
            //Create and call:
            //quantum spin-spin corrs
            //<s_i^2>
            //<n_iup n_idn>
            //spin resolved local density

        }

        if(ModelType=="SquareLattice"){

            Parameters_SQL Parameters_SQL_;
            Parameters_SQL_.Initialize(model_inputfile);

            Coordinates_SQL Coordinates_SQL_(Parameters_SQL_.lx, Parameters_SQL_.ly, Parameters_SQL_.n_orbs);
            Connections_SQL Connections_SQL_(Parameters_SQL_, Coordinates_SQL_);

            Observables_SQL Observables_SQL_(Parameters_SQL_, Coordinates_SQL_, Connections_SQL_ );

            Parameters_SQL_.beta=Parameters_.beta;
            Parameters_SQL_.mus=Parameters_.mus;
            Parameters_SQL_.eta=Parameters_.eta_dos;
            Parameters_SQL_.domega=Parameters_.dw_dos;
            Observables_SQL_.Ham_ = Hamiltonian_.Ham_;
            Observables_SQL_.eigs_ = Hamiltonian_.eigs_;


            Observables_SQL_.Calculate_Akw();
            //Observables_SQL_.Calculate_OrbResolved_Nw();
            Observables_SQL_.Calculate_Nw();
            //Create and call:
            //quantum spin-spin corrs
            //<s_i^2>
            //<n_iup n_idn>
            //spin resolved local density

        }


        if(ModelType=="TriangularLattice"){

            Parameters_TL Parameters_TL_;
            Parameters_TL_.Initialize(model_inputfile);

            Coordinates_TL Coordinates_TL_(Parameters_TL_.lx, Parameters_TL_.ly, Parameters_TL_.n_orbs);
            Connections_TL Connections_TL_(Parameters_TL_, Coordinates_TL_);

            Observables_TL Observables_TL_(Parameters_TL_, Coordinates_TL_, Connections_TL_ );

            Parameters_TL_.beta=Parameters_.beta;
            Parameters_TL_.mus=Parameters_.mus;
            Parameters_TL_.eta=Parameters_.eta_dos;
            Parameters_TL_.domega=Parameters_.dw_dos;
            Observables_TL_.Ham_ = Hamiltonian_.Ham_;
            Observables_TL_.eigs_ = Hamiltonian_.eigs_;


            Observables_TL_.Calculate_Local_spins_resolved();
            //Observables_TL_.Calculate_Akw();
            Observables_TL_.Calculate_Akxw_ribbon();
            //Observables_TL_.Calculate_OrbResolved_Nw();
            Observables_TL_.Calculate_Nw();
            //Create and call:
            //quantum spin-spin corrs
            //<s_i^2>
            //<n_iup n_idn>
            //spin resolved local density

        }



    }
    else{


        if (argc<2) { throw std::invalid_argument("USE:: executable inputfile"); }
        string inputfile = argv[1];

        Parameters Parameters_;
        Parameters_.Initialize(inputfile);

        Coordinates Coordinates_(Parameters_.ns);

        mt19937_64 Generator_(Parameters_.RandomSeed); //for random fields
        mt19937_64 Generator2_(Parameters_.RandomDisorderSeed); //for random disorder

        MFParams MFParams_(Parameters_,Coordinates_,Generator_, Generator2_);
        Hamiltonian Hamiltonian_(Parameters_,Coordinates_,MFParams_);
        Observables Observables_(Parameters_,Coordinates_,MFParams_,Hamiltonian_);

        //SVD TEST START
//        Matrix<double> A_;
//        Matrix<double> VT_;
//        Matrix<double> U_;
//        vector<double> Sigma_;
//        A_.resize(6,5);
//        A_(0,0)=1.0;A_(0,4)=2.0;A_(1,2)=3.0;A_(3,1)=2.0;
//        A_(4,1)=4.0;A_(5,4)=-5.0;
//        cout<<"A"<<endl;
//        A_.print();
//        Observables_.Perform_SVD(A_,VT_,U_,Sigma_);
//        cout<<"A"<<endl;
//        A_.print();
//        cout<<"VT"<<endl;
//        VT_.print();
//        cout<<"U"<<endl;
//        U_.print();
//        cout<<"Printing Sigma with size ="<<Sigma_.size()<<endl;
//        for(int i=0;i<Sigma_.size();i++){
//            cout<<setprecision(10)<<Sigma_[i]<<"  ";
//        }
//        cout<<endl;
//        assert(false);
        //SVD TEST DONE
        //IT IS WORKING :)




        cout<<setprecision(9);
        SelfConsistencyEngine SelfConsistencyEngine_(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_);
        SelfConsistencyEngine_.RUN_SelfConsistencyEngine();

        cout<<"_________________________________"<<endl;
        for(int i=0;i<Hamiltonian_.eigs_.size();i++){
            cout<<Hamiltonian_.eigs_[i]<<endl;
        }
        cout<<"_________________________________"<<endl;

       // Observables_.Calculate_Local_n_orb_resolved();
        //Observables_.Calculate_SpinSpincorrelations_Smartly();
       // Observables_.Calculate_DenDencorrelations_Smartly();
        Observables_.Calculate_Nw();
       // Observables_.Calculate_Local_spins_resolved();


    }



    cout << "--------THE END--------" << endl;
} // main
