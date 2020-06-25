#ifndef _TestFi
#define _TestFi  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "FiOracle.h"
#include "NDOSlver.h"
#include "structures.hpp"

#include <string>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <deque>
#include <cfloat>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <list>
//#include <iomanip>



using namespace std;


/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

namespace NDO_di_unipi_it
{
	
    /*--------------------------------------------------------------------------*/
    /*----------------------------- CLASS TestFi -------------------------------*/
    /*--------------------------------------------------------------------------*/
	
	
    class TestFi : public FiOracle
    {
        
        /*--------------------------------------------------------------------------*/
        /*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
        /*--------------------------------------------------------------------------*/
        
    public:
        
        //model primal variables
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> rc; //reduced cost
        double lcost;
        
        std::vector<std::vector<double> > memx;
        std::vector<std::vector<double> > memy;

        double knapsack(int a,  const std::vector<double> & rc, std::vector<double> & x);
	    void formfinalsol();

        int reset_topology();
	    void reset( int unfixd, const std::vector<double> & hs, bool hotstart);

        //getters translaters
        void transp_h(std::vector<double> & y1, double div, int place=0);
        
        //modifiers
        int deque_concat(const std::deque<int> & out,const std::deque<int> & from);
        void insert_in_place(const std::deque<int> & in, int nplace);
        
        //randomers
        double re_rand_fix(const std::vector<int> & ya2, const std::vector<int> & ya1,  const std::vector<double> & y1, std::deque<std::string> & hist);


	
	   const Data *data;
	   int ndemands, narcs, nnodes;
	   int szopnd, szunfxd, sznon0;
	   std::deque<int> non0; 
	   std::vector<int> h1;
        
	   std::list<HeapCell> heap;  
	   std::vector<int> Iu; //index of lag multipliers

	
	  TestFi(const Data  *d);
	    
        /*--------------------------------------------------------------------------*/
        /*--------------------------- PUBLIC METHODS -------------------------------*/
        /*--------------------------------------------------------------------------*/
        /*---------------------------- CONSTRUCTOR ---------------------------------*/
        /*--------------------------------------------------------------------------*/
        
        TestFi(std::string s);
        
        /* NV is the number of variables, and L0 the "center" of the function. */
        
        /*--------------------------------------------------------------------------*/
        /*-------------------------- OTHER INITIALIZATIONS -------------------------*/
        /*--------------------------------------------------------------------------*/
        
        void SetMaxName( cIndex MxNme = 0 )
        {
            MBN = MxNme;
            memx.resize(MBN+1);
            memy.resize(MBN+1);
        }

        
        /*--------------------------------------------------------------------------*/
        /*-------------- METHODS FOR READING THE DATA OF THE PROBLEM ---------------*/
        /*--------------------------------------------------------------------------*/
         void SetNDOSolver( NDOSolver *NwSlvr = 0 )
        {
            Slvr = NwSlvr;
        }
        
        Index GetNumVar( void ) const
        {
            return( NumVar );
        }
        
        bool GetUC( cIndex i )
		{
			return( true );
		}
		
        
        /*--------------------------------------------------------------------------*/
        /*----------------------- METHODS FOR CHANGING DATA ------------------------*/
        /*--------------------------------------------------------------------------*/
        
        void SetLambda( cLMRow Lmbd = 0 );
        
        void SetLamBase( cIndex_Set LmbdB = 0 , cIndex LmbdBD = 0 );
        
        void SetGiName( cIndex Name );
        
        
        /*--------------------------------------------------------------------------*/
        /*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
        /*--------------------------------------------------------------------------*/
        
        HpNum Fi( cIndex wFi = Inf<Index>() );
        void calcul_rc();
        
        HpNum solve();
       
        
        
        /*--------------------------------------------------------------------------*/
        /*---------------------- METHODS FOR READING RESULTS -----------------------*/
        /*--------------------------------------------------------------------------*/
        
        bool NewGi( cIndex wFi = Inf<Index>() );
        
        Index GetGi( SgRow SubG , cIndex_Set &SGBse , cIndex Name = Inf<Index>() ,
                    cIndex strt = 0 , Index stp = Inf<Index>() );
        
        HpNum GetVal( cIndex Name = Inf<Index>() );
        
            
        
        /*--------------------------------------------------------------------------*/
        /*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
        /*--------------------------------------------------------------------------*/
               
		void Aggregate( cHpRow Mlt , cIndex_Set NmSt , cIndex Dm ,
                       cIndex NwNm );
        
        /*--------------------------------------------------------------------------*/
        /*------------------------------ DESTRUCTOR --------------------------------*/
        /*--------------------------------------------------------------------------*/
        
        virtual ~TestFi();
        
        /* Destructor of the class: it must be virtual. */
        
        /*--------------------------------------------------------------------------*/
        /*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
        /*--------------------------------------------------------------------------*/
        /*--                                                                      --*/
        /*--  The standard user should not care about the following part: users   --*/
        /*--  who need to extend the code by deriving a new class may use these   --*/
        /*--  methods and data structures. It is *dangerous* to *modify* the      --*/
        /*--  data structures, while it safe to read them.                        --*/
        /*--                                                                      --*/
        /*--------------------------------------------------------------------------*/
        
    protected:
        
        /*--------------------------------------------------------------------------*/
        /*-------------------------- PROTECTED METHODS -----------------------------*/
        /*--------------------------------------------------------------------------*/
        
        /*--------------------------------------------------------------------------*/
        /*----------------------- PROTECTED DATA STRUCTURES  -----------------------*/
        /*--------------------------------------------------------------------------*/
        
        Index NumVar;      // the number of variables
        LMNum Cntr;        // the minimum of the function (L0)
        
        Index MBN;         // max n. of names
       
		NDOSolver *Slvr; 

        cLMRow Lam1;       // point at which to evaluate the function
        cIndex_Set LamB;   // vector of the indices of nonzeroes in Lam1
        Index LamBd;       // lenght of LamB
        
        /*--------------------------------------------------------------------------*/
        /*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
        /*--------------------------------------------------------------------------*/
        /*--                                                                      --*/
        /*-- Nobody should ever look at this part: everything that is under this  --*/
        /*-- advice may be changed without notice in any new release of the code. --*/
        /*--                                                                      --*/
        /*--------------------------------------------------------------------------*/
        
    private:
        
        /*--------------------------------------------------------------------------*/
        /*-------------------------- PRIVATE METHODS -------------------------------*/
        /*--------------------------------------------------------------------------*/
        
        /*--------------------------------------------------------------------------*/
        /*----------------------- PRIVATE DATA STRUCTURES  -------------------------*/
        /*--------------------------------------------------------------------------*/
        
        /*--------------------------------------------------------------------------*/
        
    };  // end( class TestFi )
    
    /*--------------------------------------------------------------------------*/
    
};  // end( namespace FiOracle_NDO_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* TestFi.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File TestFi.h ------------------------------*/
/*--------------------------------------------------------------------------*/

