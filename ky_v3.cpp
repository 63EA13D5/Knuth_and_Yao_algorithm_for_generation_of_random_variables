/***** BEGINNING OF DISCLAIMER *****
       
       This code is provided without warranty of any kind, either express or implied. Use at your own risk.
       
       The use of this is done at your own discretion and risk and with agreement that you will be solely responsible for any damage to your computer system or loss of data that results from the use of this code. You are solely responsible for adequate protection and backup of the data and equipment used in connection with this code, and we will not be liable for any damages that you may suffer in connection with using, modifying or distributing any of this code. No advice or information, whether oral or written, obtained by you from us or from this website shall create any warranty for this code.
       
       We make makes no warranty that:
       
       1) the code will meet your requirements
       2) the code will be uninterrupted, timely, secure or error-free
       3) the results that may be obtained from the use of the software will be effective, accurate or reliable
       4) the quality of the code will meet your expectations
       5) any errors in the code obtained from us will be corrected. 
       
       The code and its documentation made available on this website:
       
       1) could include technical or other mistakes, inaccuracies or typographical errors. We may make changes to the code or documentation made available on its web site at any time without prior-notice.
       2) may be out of date, and we make no commitment to update such materials. 
       
       We assume no responsibility for errors or omissions in the code or documentation available from its web site.
       
       In no event shall we be liable to you or any third parties for any special, punitive, incidental, indirect or consequential damages of any kind, or any damages whatsoever, including, without limitation, those resulting from loss of use, data or profits, and on any theory of liability, arising out of or in connection with the use of this code. 
       
       ***** END OF DISCLAIMER *****/

/*
  
  AUTHORS: Claude Gravel
  
  DATE: July 27th, 2020 
  
  TECHNICAL REFERENCE: Original reference by Knuth and Yao 1976

  DESCRIPTION:

  The code here implement the Knuth and Yao algorithm for sampling probability distribution.
  More precisely, given the description of a probabilitiy distribution through a vector of
     arbitrary finite length, the algorithm uses random unbiased i.i.d. bits from a source to generate
     an instance using nearly the entropy as the number of random bits which is optimal from an information theory perspective.
  This code relies on NTL library which may (or may not) depend on GMPL.

  NOTES:

  The libraries are used only to compute probabilities with arbitrary precision.
  For information about NTL (Number Theory Library) see: https://www.shoup.net/ntl/
  For information about GMPL (GNU Multiple Precision Arithmetic Library) see: https://gmplib.org/
 
  The words "level" and "depth" can sometimes be used interchangeably within comments left in this code.
*/


#include <iomanip>
#include <iostream>
#include <random>
#include <ctime>
#include <chrono>
#include <cstring>

#include <NTL/ZZ.h>
#include <NTL/RR.h>

/**** ***** ***** ***** ***** **** ***** ***** ***** *****/

std::random_device r_dev;
std::mt19937_64 mt(r_dev());
//std::knuth_b mt(r_dev());

std::bernoulli_distribution distBer(0.5);

/**** ***** ***** ***** ***** **** ***** ***** ***** *****/

/*
  General rules to set the minimal accuracy.

  1) Size of the support which is the length of the probability vector, let l1 = log2(size support)
     Note here that the entropy of the probability vector <= log2(size support) so this is 
     always a safe bet whenever one does not know the exact entropy analitically.

  2) If the probability vector represents a truncation of an infinite countable discrete 
     distribution with delta as the leftover tail probilities, let l2 = -log2(-delta). Usually 
     delta is very small here.

  3) If the probability vector represents a discretization of an absolutely continuous distribution
     with 0 < epsilon < 1 as the discretization precision, let l3 = -log2(-epsilon)
     (For singular distribution, this is more a case by case analysis. Also for an absolutely continuous, 
      it is assumed that its differential entropy is bounded above so that the discretization
      is relevant.)

  4) An extra security parameter, say l4, that depends on how much resource one can use to
     sample enough in order to distinguish from the ideal (targer) distribution.

  The general rule, take min_accuray = l1 + l2 + l3 + l4;

  min_accuracy is passed as an argument to NTL::RR::SetPrecision
*/

const long min_accuracy = 1000;//ok this is a bit overkilling perhaps

/**** ***** ***** ***** ***** **** ***** ***** ***** *****/

/*
  To count the average number of random bits needed from
  the generator per random instance in order to compare
  with the numerically computed entropy.
*/

long ct_rnd_bit=0;

/**** ***** ***** ***** ***** **** ***** ***** ***** *****/

void roll_dice(long & outcome, long nb_faces)
{
  /*
    Outcome a uniform discrete random varible over nb_faces possible outcomes.
    It uses no more than entropy + 2 bits from the source generator according Knuth and Yao's result
    The mimimum expected number of random bits from an information point of view is log_2(nb_faces).
    Based on J. Lumbroso phd thesis idea.
    This is basically Knuth and Yao 1976 algorithm for uniform discrete distribution.
  */
  
  long x = 0;
  long y = 1;

  while(true)
    {
      y = 2*y;
      x = 2*x + distBer(mt);
      ct_rnd_bit++;
      if(y>=nb_faces)
	{
	  if(x<nb_faces)
	    {
	      outcome = x;
	      break;
	    }
	  else
	    {
	      y = y - nb_faces;
	      x = x - nb_faces;
	    }
	}
    }
}

/**** ***** ***** ***** ***** **** ***** ***** ***** *****/

struct prob_dist_data{

  long nb_atoms;//size of probability distribution support or equivalently the lenght of the probability vector
  
  std::vector<NTL::RR> pmv;//pmv.size()=nb_atoms;
  /*
    pmv = probability mass vector
    pmv[i] > 0
    sum pmv[i] = 1
    pmv[i] is accurate up to min_accuracy bits
  */

  std::vector<long> nb_internal_nodes;

  std::vector<std::vector<long>> L;
  /*
    L.size = m;
    L[i].size = number of leaves at depth i, L[i] is a lookup table with outcomes
    L[i][j] = k if and only if p_{k, i} = 1 for some k
    0 <= j < L[i].size
    
    nb_internal_nodes[i]= L[i].size + (number of non-leave at depth i) which is
       computed with the formula. We don't want to recompute this quantity 
       for every outcome generated if we are generating a sample for instance.
       Thus it is faster to maintain this auxiliary info for generating i.i.d. sample.
  */
};

/**** ***** ***** ***** ***** **** ***** ***** ***** *****/

void entropy(NTL::RR & binary_entropy, struct prob_dist_data P)
{
  //Store the binary entropy of P in binary_entropy
  
  NTL::RR r=NTL::RR(0.0);

  for(long i = 0; i < P.nb_atoms ; i++)
    {

      r = r - P.pmv[i]*NTL::log(P.pmv[i]);
      
    }
  r = r/NTL::log(NTL::RR(2.0));
  binary_entropy = r;
}

void get_DDG_tree(struct prob_dist_data & T)
{
  /*
    Get the binary representations of T.pmv and then get the number of internal
    nodes for the DDG tree reprensenting T. The representation of the DDG is
    canonical in the sens that a leaf's symbol corresponds to its index from left
    to right. Hence the random outcome can be remapped. 
  */
  
  if((T.nb_atoms==0)|| T.pmv.empty() ){std::cerr << "Wrong selection number --- exit\n\n"; exit(-1);}
  
  long ct_atoms=0;
  
  NTL::RR chk_sum=NTL::RR(0.0);

  long max_ee = 0;
  
  ct_atoms=0;
  do
    {
      
      if((T.pmv[ct_atoms].mantissa()==NTL::ZZ(0)) && (T.pmv[ct_atoms].exponent()==0))
	{
	  std::cerr << "zero prob --- exit" << ct_atoms << " " << T.pmv[ct_atoms] << "\n\n";
	  exit(-1);
	}
      
      chk_sum += T.pmv[ct_atoms];
      
      //NTL::ZZ mm = T.pmv[ct_atoms].mantissa();
      long ee = T.pmv[ct_atoms].exponent();
      
      if(-ee > max_ee)
	{
	  max_ee = -ee;
	}
      ct_atoms++;
    }
  while(ct_atoms<T.nb_atoms);

  
  if( (chk_sum - NTL::RR(1.0)) >= NTL::RR(T.nb_atoms)*NTL::power(NTL::RR(2.0),-min_accuracy) )
    {
      //For the sum of two p-bit numbers, we lose one bit. The minimal accuracy is multiplied
      //  the size of the support in the if statement.
      
      std::cerr << "check sum prob vector error --- exit\n\n";
      std::cerr << chk_sum-NTL::RR(1.0) << "  " << NTL::power(NTL::RR(2.0),-min_accuracy) << "\n";
      exit(-1);
    }

  T.L.reserve(max_ee);
  for(long i = 0 ; i < max_ee+1; i++)
    {
      std::vector<long> tmp;
      T.L.push_back(tmp);
    }

  ct_atoms=0;
  do
    {
      NTL::ZZ mm = T.pmv[ct_atoms].mantissa();
      long ee = T.pmv[ct_atoms].exponent();
      
      long j = -ee-1;
       
      do
	{
	  if((( mm&(NTL::ZZ(1)<<j) ) >> j) == NTL::ZZ(1))
	    {
	      T.L[-ee-1-j].push_back(ct_atoms);
	    }
	  j--;
	}
      while(j>-1);
     
      ct_atoms++;
    }
  while(ct_atoms<T.nb_atoms);
  
  
  /****************************************/
  
  T.nb_internal_nodes.push_back(2);//(there are always 2 nodes at depth 1 indexed by (1 - 1) )
  for(long l=1;l<min_accuracy;l++) 
    {
      //T->nb_nodes[l-1]=(1L<<l);//maximum possibles of nodes at depth l indexed by (l-1) + 1, 
      long tmp = 0;
      for(long m=0;m<l;m++)//account of previous levels
	{
	  tmp = tmp + (1L<<(l-m))*(T.L[m].size());
	}
      T.nb_internal_nodes.push_back((1L<<(l+1))-tmp);
    }
}

/**** ***** ***** ***** ***** **** ***** ***** ***** *****/

void gen_rnd_var(long & index_rv, struct prob_dist_data T)
{
   /*
    index_rv is a random outcome from distribution described by T.
    0 <= index_rv < size of support of T = T.nb_atoms

    A user can transform the range of index_rv if required. For instance by an affine transformation.
    
    It uses no more than the entropy(T) + 2 bits on average from the random source.
    A simple adaptation of the roll_dice above using lists.
    This is basically Knuth and Yao 1976 algorithm for discrete distributions.
   */
  
  long x=0;
  long y=1;
  long lev=0;
  
  while(true)
    {
      x = 2*x + distBer(mt);
      ct_rnd_bit++;
      
      y = 2*y;
      if( y >= T.nb_internal_nodes[lev])//number of internal nodes = exit nodes (leaves) + non-exit nodes (continue to next level)
	{
	  
	  if(x<(long)T.L[lev].size())//number of leaves
	    {
	      
	      index_rv = T.L[lev][x];
	      return;
	    
	    }
	  else
	    {
	      y = y - T.L[lev].size();
	      x = x - T.L[lev].size();
	    }
	}
      lev++;
    }
}

/**** ***** ***** ***** ***** **** ***** ***** ***** *****/

/*
  In ky_utils.h: A module to query users for specific information
                 about certain kinds of probability distributions.

		 A simple check of the top N most probable outcomes together
		 with the empirical observed frequencies. N is ask to the user
		 together with the sample size desired.
*/

#include "ky_utils.h" 

int main(void)
{
 
  std::cout.precision(15);
  std::cout.setf( std::ios::fixed, std::ios::floatfield );

  NTL::RR::SetPrecision(min_accuracy);
  NTL::RR::SetOutputPrecision(15);

  struct prob_dist_data P;

  
  query_user(P);//Ask user for some kind of distributions.
  
  get_DDG_tree(P);//Get the Discrete Data Generator tree (list of lists representation here)

  /*
    Example of a call to the Knuth & Yao algorithm.
    Note gen_rnd_var increments the counter of the number of 
    call to the generator distBer(0.5) above.
  */
  if(false)
    {
      long rv;
      gen_rnd_var(rv,P);
      
      std::cout << "\nP{outcome = "<<rv<<"} = " << P.pmv[rv] << "\n";
    }

  /*
    Compare frequencies (empirical and theoritical)
    
    Compare average number of bits required to 
     generate an instance and the entropy.
    
    Get average time to generate a random instance.
  */
  simple_checks(P);


    
  return 0;
}
