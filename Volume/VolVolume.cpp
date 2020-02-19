// $Id: VolVolume.cpp 281 2011-01-05 23:58:31Z lou $
// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "VolVolume.hpp"

//#############################################################################
/// Usage: v=w; where w is a VOL_dvector
VOL_dvector&
VOL_dvector::operator=(const VOL_dvector& w) {
   if (this == &w) 
      return *this;
   delete[] v;
   const int wsz = w.size();
   if (wsz == 0) {
      v = 0;
      sz = 0;
   } else {
      v = new double[sz = wsz];
      for (int i = sz ;i--;)
      v[i] = w[i];
   }
   return *this;
}
/// Usage v=w; where w is a double. It copies w in every entry of v
VOL_dvector&
VOL_dvector::operator=(const double w) {
   for (int i = sz ;i--;)
	v[i] = w;
   return *this;
}

//#############################################################################
/// Usage: v=w; where w is a VOL_ivector
VOL_ivector&
VOL_ivector::operator=(const VOL_ivector& w) {
   if (this == &w) return *this;
   delete[] v;
   const int wsz = w.size();
   if (wsz == 0) {
      v = 0;
      sz = 0;
   } else {
      v = new int[sz = wsz];
      for (int i = sz ;i--;)
      v[i] = w[i];
   }
   return *this;
}
/// Usage v=w; where w is an int. It copies w in every entry of v
VOL_ivector&
VOL_ivector::operator=(const int w) {
   for (int i = sz ;i--;)
   v[i] = w;
   return *this;
}

//############################################################################
/// Dual step. It takes a step in the direction v
// lcost is a member of VOL_dual
double 
VOL_dual::step(const double target, const double lambda, VOL_dual & trial, const VOL_dvector& v, double & t) {

   const int nc = u.size();

   int i;
   double viol = 0.0;
   for (i=nc;i--;) {
		viol += v[i] * v[i];
   }
  
   t = viol == 0.0 ? 0.0 : (target - lcost) / viol * lambda;

   for (i=nc;i--;) {
		trial.u[i] = u[i] + t * v[i];
   }
   

   return viol;
}

//############################################################################
/** Computing inner products. */
VOL_vh::VOL_vh(const VOL_dvector& v, const VOL_dvector& vstar) : 
  vh(0), vv(0)
{
	int i;
   const int nc = vstar.size();
   for (i=nc;i--;) {
      vh   += v[i] * vstar[i];
      vv +=  v[i] * v[i];
   }
}


//############################################################################

// reading parameters that control the algorithm
void
VOL_problem::read_params(const char* filename) 
{
   char s[100];
   FILE* infile = fopen(filename, "r");
   if (!infile) {
      printf("Failure to open file: %s\n", filename);
      abort();
   }
   while (fgets(s, 100, infile)) {
      const int len = strlen(s) - 1;
      if (s[len] == '\n')
	 s[len] = 0;
      std::string ss(s);

      if (ss.find("ubinit") == 0) {
	 int i = ss.find("=");  
	 parm.ubinit = atof(&s[i+1]);

      } else if (ss.find("printflag") == 0) {
	 int i = ss.find("=");  
	 parm.printflag = atoi(&s[i+1]);
	 
      } else if (ss.find("printinvl") == 0) {
	 int i = ss.find("=");  
	 parm.printinvl = atoi(&s[i+1]);

      } else if (ss.find("maxsgriters") == 0) {
	 int i = ss.find("=");  
	 parm.maxsgriters = atoi(&s[i+1]);
	 
      } else if (ss.find("greentestinvl") == 0) {
	 int i = ss.find("=");  
	 parm.greentestinvl = atoi(&s[i+1]);
	 
      } else if (ss.find("yellowtestinvl") == 0) {
	 int i = ss.find("=");  
	 parm.yellowtestinvl = atoi(&s[i+1]);
	 
      } else if (ss.find("redtestinvl") == 0) {
	 int i = ss.find("=");  
	 parm.redtestinvl = atoi(&s[i+1]);
	 
      } else if (ss.find("lambdainit") == 0) {
	 int i = ss.find("=");  
	 parm.lambdainit = atof(&s[i+1]);
	 
      } else if (ss.find("alphainit") == 0) {
	 int i = ss.find("=");  
	 parm.alphainit = atof(&s[i+1]);
	 
      } else if (ss.find("alphamin") == 0) {
	 int i = ss.find("=");  
	 parm.alphamin = atof(&s[i+1]);
	 
      } else if (ss.find("alphafactor") == 0) {
	 int i = ss.find("=");  
	 parm.alphafactor = atof(&s[i+1]);
	 
      } else if (ss.find("alphaint") == 0) {
	 int i = ss.find("=");  
	 parm.alphaint = atoi(&s[i+1]);

      }else if (ss.find("maxtime") == 0) {
      int i = ss.find("=");
      parm.maxtime = atof(&s[i+1]);
      
      }else if (ss.find("epslin") == 0) {
      int i = ss.find("=");
      parm.epslin = atof(&s[i+1]);
      
      }else if (ss.find("tstar") == 0) {
      int i = ss.find("=");
      parm.tstar = atof(&s[i+1]);
      }

   }
   fclose(infile);
}

//#############################################################################

void
VOL_problem::set_default_parm()
{
   parm.epslin = 1e-4;
   parm.tstar = 1.0;	
	
   parm.lambdainit = 0.1;
   parm.alphainit = 0.01;
   parm.alphamin = 0.001;
   parm.alphafactor = 0.5;
   parm.ubinit = DBL_MAX;
 
   parm.maxsgriters = 2000;
   parm.printflag = 3;
   parm.printinvl = 50;
   parm.greentestinvl = 1;
   parm.yellowtestinvl = 2;
   parm.redtestinvl = 10;
   parm.alphaint = 80;
    parm.maxtime = 3600.0;
}
   
//#############################################################################

VOL_problem::VOL_problem() : 
   alpha_(-1),
   lambda_(-1),
   iter_(0),
   value(-1),
   psize(-1),
   dsize(-1),
   maxlb(DBL_MAX)
{
   set_default_parm();
}

//
VOL_problem::VOL_problem(const char *filename,std::string instance) :
   alpha_(-1),
   lambda_(-1),
   iter_(0),
  value(-1),
    psize(-1),
   dsize(-1),
   maxlb(DBL_MAX)
{
   set_default_parm();
   read_params(filename);
   
}

//######################################################################

VOL_problem::~VOL_problem()
{
  // delete[] parm.temp_dualfile;
  
}

//######################################################################
/// print information about the current iteration
void
VOL_problem::print_info(const int iter,
			const VOL_primal& primal, const VOL_primal& pstar,
			const VOL_dual& dual)
{
   printf("%i. L=%f P=%f \n",
	  iter, dual.lcost, pstar.value);
}

//######################################################################
/// this is the Volume Algorithm
int
VOL_problem::solve(VOL_user_hooks& hooks, const bool use_preset_dual) 
{
   clock_t t_u;
   struct tms buff;
   times( &buff );
   t_u = buff.tms_utime;
   
   double timespent = 0;
   
   if (initialize(use_preset_dual) < 0) // initialize several parameters
      return -1;

   double best_ub = parm.ubinit;      // upper bound
   int retval = 0; 

   VOL_dvector rc(psize); // reduced costs
   VOL_dual dual(dsize); // dual vector
   dual.u = dsol;
   VOL_primal primal(psize, dsize);  // primal vector

   retval = hooks.compute_rc(dual.u, rc); // compute reduced costs
   if (retval < 0)  return -1;
   retval = hooks.solve_subproblem(dual.u, rc, dual.lcost,
				   primal.x, primal.v, primal.value);// solve relaxed problem
   if (retval < 0)  return -1;
   double target = readjust_target(-DBL_MAX/2, dual.lcost); // set target for the lagrangian value
   
   VOL_primal pstar(primal); // set pstar=primal
   VOL_dual dstar(dual); // dstar is the best dual solution so far
   VOL_swing swing;
   VOL_alpha_factor alpha_factor;
  
   
   iter_ = 0;
   if (parm.printflag)
     print_info(iter_, primal, pstar, dual);
	
	double dnorm = 0.0;
	double linerror = 0;
	double sigma = 0;
	double t = 0;
	double mult =1;
	dnorm = dstar.step(target, lambda_,dual, pstar.v, t);// take a dual step
	
   for (iter_ = 1; iter_ <= parm.maxsgriters; ++iter_) {  // main iteration
	   
      hooks.compute_rc(dual.u, rc);  // compute reduced costs
      hooks.solve_subproblem(dual.u, rc, dual.lcost,
				      primal.x, primal.v, primal.value);// solve relaxed problem
    
	 	 
      if (iter_ % parm.alphaint == 0) { // change alpha if needed
		  double best = dstar.lcost;
		  if (dual.lcost > dstar.lcost) best = dual.lcost;
		  const double fact = alpha_factor.factor(parm, best, alpha_);
		  if (fact != 1.0 && (parm.printflag & 2)) {
			printf(" ------------decreasing alpha to %f\n", alpha_*fact);
			}
		  alpha_ *= fact;
      }
       if ((iter_ % parm.printinvl == 0) && (parm.printflag<3) && (parm.printflag>0)) { // printing iteration information
	 print_info(iter_, primal, pstar, dual);
	 swing.print();
      }
      
      // compute inner product between the new subgradient and the
      // last direction. This to decide among green, yellow, red
      VOL_vh prod(primal.v, pstar.v);
      // green, yellow, red
      swing.cond(dstar, dual.lcost, prod.vh, iter_);
      // change lambda if needed
      lambda_ *= swing.lfactor(parm, lambda_, iter_);
      
     
      
      if (dual.lcost > dstar.lcost) {
		sigma += dstar.lcost - dual.lcost;
		int i;
		for(i=dsize;i--;) sigma +=(dual.u[i]-dstar.u[i])*pstar.v[i];
		dstar = dual; // update dstar
		linerror =0;
      }else{
		  linerror = dual.lcost  - dstar.lcost - prod.vh*t;
		  if(linerror<1e-30 ){
			  linerror =0;
		   }
       }
      // check if target should be updated
		target = readjust_target(target, dstar.lcost);
      // convex combination with new primal vector
      
      
      mult = power_heur(prod, linerror, sigma);
      pstar.cc(mult, primal);
	  dnorm = dstar.step(target, lambda_,dual, pstar.v, t);// take a dual step
	  sigma = mult*linerror + (1 - mult)*sigma;
      
      
      // test terminating criteria
       if(parm.tstar*dnorm + sigma <= parm.epslin*dstar.lcost){ retval=1; break;}
	  if(dstar.lcost-maxlb>=10){ retval =2; break;}
       
 	  times( &buff );
	  timespent = ( double( buff.tms_utime - t_u ) ) / double( CLK_TCK );
      
      if (timespent > parm.maxtime) {
          retval=0;
           break;
       }

   }

   if (parm.printflag)
     print_info(iter_, primal, pstar, dual);
   // set solution to return
   value = dstar.lcost;
   psol = pstar.x;
   dsol = dstar.u;

   return retval;
}


/// A function to initialize a few variables
int
VOL_problem::initialize(const bool use_preset_dual) {
  
   // setting initial values for parameters
   alpha_       = parm.alphainit;
   lambda_      = parm.lambdainit;
   // check if there is an initial dual solution
   if (use_preset_dual) {
      if (dsol.size() != dsize) {
	 printf("size inconsistent (dsol)\n");
	 return -1;
      }
   } else {
      dsol.clear();
      dsol.allocate(dsize);
      dsol = 0.0;
   }

   return 0;
}


/// Here we increase the target once we get within 5% of it
double
VOL_problem::readjust_target(const double oldtarget, const double lcost) const
{
   double target = oldtarget;
   if (lcost >= target - VolAbs(target) * 0.05) {
      if (VolAbs(lcost) < 10.0) {
	 target = 10.0;
      } else { 
	 target += 0.025 * VolAbs(target);
	 target = VolMax(target, lcost + 0.05 * VolAbs(lcost));
      }
      if (target != oldtarget && (parm.printflag & 2)) {
	 printf("     **** readjusting target!!! new target = %f *****\n",
		target);
      }
   }
   return target;
}

/** Here we decide the value of alpha_fb to be used in the convex
    combination. More details of this are in doc.ps
    IN:  alpha, primal, pstar, dual
    OUT: pstar = alpha_fb * pstar + (1 - alpha_fb) * primal
*/
double
VOL_problem::power_heur(const VOL_vh& prod, const double & linerror, const double & sigma) const 
{
   const double alpha = alpha_;
   
   double a_asc = (alpha * prod.vv - prod.vh) / (prod.vv - prod.vh);
   double alpha_fb;
   
   
   alpha_fb = alpha;
   
   if (alpha_fb < a_asc)
      alpha_fb = a_asc;
   if (alpha_fb > 1.0)
      alpha_fb = alpha;
   if (alpha_fb < 1e-8)
      alpha_fb = alpha / 10.0;

   return alpha_fb;
}

