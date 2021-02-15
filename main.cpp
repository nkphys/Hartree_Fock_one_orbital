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
#include "random"


int main(int argc, char *argv[]) {

    string ex_string_original =argv[0];

    string ex_string;
    //ex_string.substr(ex_string_original.length()-5);
    ex_string=ex_string_original.substr (2);
    cout<<"'"<<ex_string<<"'"<<endl;



    if(ex_string=="k_space_SelfConsistency"){
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

    }




    else if(ex_string=="CreateConnections"){

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

            Connections_TL_.InteractionsCreate();
            Connections_TL_.Print_LongRangeInt();
            Connections_TL_.Print_Spin_resolved_OnsiteE();

        }
    }
    else if(ex_string=="observe"){

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
            Observables_SQL_.Calculate_OrbResolved_Nw();
            Observables_SQL_.Calculate_Nw();
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
        cout<<setprecision(9);
        SelfConsistencyEngine SelfConsistencyEngine_(Parameters_,Coordinates_,MFParams_,Hamiltonian_,Observables_);
        SelfConsistencyEngine_.RUN_SelfConsistencyEngine();

        Observables_.Calculate_Local_n_orb_resolved();
        Observables_.Calculate_SpinSpincorrelations_Smartly();
        Observables_.Calculate_DenDencorrelations_Smartly();
        Observables_.Calculate_Nw();
        Observables_.Calculate_Local_spins_resolved();


    }



    cout << "--------THE END--------" << endl;
} // main
