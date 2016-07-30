    // DMRG for transverse field ising model
    
    #include "itensor/all.h"

    #include <iostream>
    #include <fstream>
    #include <boost/lexical_cast.hpp>     

    using namespace itensor;
    using namespace std ;
    using boost::lexical_cast;
    using boost::bad_lexical_cast;

    int main()

    {

        double h=1.0; // h=1 is the critical point
        
        ofstream fout;
        fout.open("Energy-TFIM-APBC.dat");// Opens a file for output

        for (int N = 100; N < 1000; N +=11)
	{
     

         auto sites = SpinHalf(N);

     

         auto ampo = AutoMPO(sites);

         for(int j = 1; j < N; ++j)

         {
            ampo += -1, "Sz",j,"Sz",j+1;
            ampo += -h,"Sx",j;

         }
         
         ampo += -h,"Sx",N;

        //ampo +=  1, "Sz",1,"Sz",N;  // Anti-Periodic Boundary condition 
        // ampo +=  -1, "Sz",1,"Sz",N;  //Periodic Boundary condition 


         auto H = MPO(ampo);

         auto sweeps = Sweeps(5); //number of sweeps is 5

         sweeps.maxm() = 10,20,100,100,200;

         sweeps.cutoff() = 1E-10;

         auto psi = MPS(sites);

     

         auto energy = dmrg(psi,H,sweeps);

     

         println("Ground State Energy = ",energy);

     

         //calculate the entanglement entropy here
         //Given an MPS or IQMPS called "psi",
         //and some particular bond "b" (1 <= b < psi.N())
         //across which we want to compute the von Neumann entanglement

     
         string N_str = lexical_cast<string>(N);
	 ofstream myfile(string("APBC-tfim-" + N_str + ".dat").c_str());// Opens a file for output
        
         for (int b=1; b<psi.N(); b++){

     

            //"Gauge" the MPS to site b

            psi.position(b); 

     

            //Here assuming an MPS of ITensors, but same code works

            //for IQMPS by replacing ITensor -> IQTensor

     

            //Compute two-site wavefunction for sites (b,b+1)

            ITensor wf = psi.A(b)*psi.A(b+1);

     

            //SVD this wavefunction to get the spectrum

            //of density-matrix eigenvalues

            auto U = psi.A(b);

            ITensor S,V;

            auto spectrum = svd(wf,U,S,V);

     

            //Apply von Neumann formula

            //spectrum.eigs() is a Vector containing

            //the density matrix eigenvalues

            //(squares of the singular values)

            Real SvN = 0.;

            for(auto p : spectrum.eigs())

            {

                if(p > 1E-12) SvN += -p*log(p);

            }

            printfln("Across bond b=%d, SvN = %.10f",b,SvN);
            myfile <<b<<'\t'<<SvN<<std::endl;

        }
         myfile <<"GS energy"<<'\t'<<energy<<std::endl;
         myfile.close();

         fout <<N<<'\t'<<energy<<std::endl;
     }
       fout.close(); 
       return 0;

    }
