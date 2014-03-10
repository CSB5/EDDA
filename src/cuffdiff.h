#ifndef _cuffdiffr_CUFFDIFF_H
#define _cuffdiffr_CUFFDIFF_H

#include <Rcpp.h>
#include <stdlib.h>
#include <getopt.h>
#include <string>
#include <numeric>
#include <limits>
#include <functional>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <assert.h>  
//#include <boost/shared_ptr.hpp>
//#include <boost/math/distributions/normal.hpp>

using namespace std;
//using namespace boost;
//using namespace math;

enum TestStatus {
  NOTEST,  // successful calculation, test not performed
  LOWDATA, // unsuccessful calculation due to low data, test not performed
  HIDATA,  // skipped calculation due to too many reads data, test not performed
  OK,      // successful numerical calc, test performed
  FAIL     // numerical exception, test not performed
};


enum AbundanceStatus { NUMERIC_OK, NUMERIC_FAIL, NUMERIC_LOW_DATA, NUMERIC_HI_DATA };


struct RPMContext
{
  RPMContext(double r, double v)
  : RPM(r), RPM_variance(v) {}
  double RPM;
  double RPM_variance;
  AbundanceStatus status;
};



// Stores the differential expression of an isoform or set of isoforms in two
// different samples, along with a significance test statistic for the difference.
struct SampleDifference
{
  SampleDifference() :
  sample_1(-1),
  sample_2(-1),
  value_1(0.0),
  value_2(0.0),
  test_stat(0.0),
  p_value(1.0),
  corrected_p(1.0),
  tested_group_id(-1),
  test_status(NOTEST),
  significant(false){}

  size_t sample_1;
  size_t sample_2;

  double value_1;
  double value_2;
  double differential;
  double test_stat;
  double p_value;
  double corrected_p;

  size_t tested_group_id; // which scaffolds' FPKMs contribute

  //shared_ptr<SampleDifferenceMetaData> meta_data;

  TestStatus test_status;
  bool significant;
};

// Returns the probability of x, given the distribution described by mu and sigma.
double pdf(double x, double mu, double sigma)
{
	static const double pi = 3.14159265; 
	return exp( -1 * (x - mu) * (x - mu) / (2 * sigma * sigma)) / (sigma * sqrt(2 * pi));
}

// Returns the integral from -inf to x of any function that accepts x and 2 other parameters
double cdf(double x, double arg1, double arg2, double(*pPDF)(double,double,double))
{

	double sum = 0;
	double ninf = -1e3; // Negative infinity, just use something small
	double n = 1e6; // The number of "buckets" that we'll calculate, more is more accurate but takes more time;
	//double ninf = -1e2;
	//double n = 1000.0;

	for (double k = 1.; k < n-1; k++)
	{
		sum = sum+ pPDF( x + k*(x-ninf)/n ,arg1,arg2);
		//cout << x + k*(x-ninf)/n << "\t" << sum << "\t" << pPDF( x + k*(x-ninf)/n ,arg1,arg2) << endl;
	}
/*	cout << "(x - ninf) / n = " << (x - ninf) / n << endl;
	cout << "pPDF(x,arg1,arg2) = " << pPDF(x,arg1,arg2) << endl;
	cout << "pPDF(ninf,arg1,arg2) = " << pPDF(ninf,arg1,arg2) << endl;
	cout << "sum = " << sum << endl;
	cout << "((pPDF(x) + pPDF(ninf)/2 + sum) = " << ((pPDF(x,arg1,arg2) + pPDF(ninf,arg1,arg2))/2 + sum) << endl;
	cout << "p = " << ((x - ninf) / n) * ((pPDF(x,arg1,arg2) + pPDF(ninf,arg1,arg2))/2 + sum)<<endl;*/
	double p=1-((x - ninf) / n) * ((pPDF(x,arg1,arg2) + pPDF(ninf,arg1,arg2))/2 + sum);
	if (p<0) {p=0.0;}
	if (p>1) {p=1.0;}
	return (p);
}



// This performs a between-group test on an isoform or TSS grouping, on two
// different samples.
double test_diffexp(const RPMContext curr, const RPMContext prev, SampleDifference test)
{
  bool performed_test = false;
  if (curr.RPM > 0.0 && prev.RPM > 0.0)
  {
    
    double stat = 0.0;
    double p_value = 1.0;
    
    if (curr.RPM_variance > 0.0 || prev.RPM_variance > 0.0)
    {
      double curr_log_rpm_var = (curr.RPM_variance) / (curr.RPM * curr.RPM);
      double prev_log_rpm_var = (prev.RPM_variance) / (prev.RPM * prev.RPM);
      
      double numerator = log(prev.RPM / curr.RPM);
      
      double denominator = sqrt(prev_log_rpm_var + curr_log_rpm_var);
      stat = numerator / denominator;
      
      
//     normal norm;
      double t1, t2;
      if (stat > 0.0)
      {
	t1 = stat;
	t2 = -stat;
      }
      else
      {
	t1 = -stat;
	t2 = stat;
      }
      
      if (isnan(t1) || isinf(t1) || isnan(t2) || isnan(t2))
      {
	
	//fprintf(stderr, "Warning: test statistic is NaN! %s (samples %lu and %lu)\n", test.locus_desc.c_str(), test.sample_1, test.sample_2);
	p_value = 1.0;
      }
      else
      {
	double tail_1 = cdf(t1,0,1,pdf);
	double tail_2 = cdf(t2,0,1,pdf);
	p_value = 1.0 - (tail_1 - tail_2);
      }
    }
    
    double differential = log2(curr.RPM) - log2(prev.RPM);
    
    //test = SampleDifference(sample1, sample2, prev.FPKM, curr.FPKM, stat, p_value, transcript_group_id);
    test.p_value = p_value;
    test.differential = differential;
    test.test_stat = stat;
    test.value_1 = prev.RPM;
    test.value_2 = curr.RPM;
    
    performed_test = true;
  }
  else
  {
    if (curr.RPM > 0.0)
    {
      if (curr.status != NUMERIC_LOW_DATA && curr.RPM_variance > 0.0)
      {
	test.p_value = cdf(0,curr.RPM, sqrt(curr.RPM_variance),pdf);
	performed_test = true;
	test.differential = numeric_limits<double>::max();;
	test.test_stat = numeric_limits<double>::max();
	test.value_1 = 0;
	test.value_2 = curr.RPM;
      }
      else
      {
	test.differential = -numeric_limits<double>::max();
	test.test_stat = -numeric_limits<double>::max();
	test.value_1 = prev.RPM;
	test.value_2 = 0;
	test.p_value = 1;
	performed_test = false;
      }
    }
    else if (prev.RPM > 0.0)
    {
      if (prev.status != NUMERIC_LOW_DATA && prev.RPM_variance > 0.0)
      {
	test.p_value = cdf(0,prev.RPM, sqrt(prev.RPM_variance),pdf);
	performed_test = true;
	
	test.differential = -numeric_limits<double>::max();
	test.test_stat = -numeric_limits<double>::max();
	test.value_1 = prev.RPM;
	test.value_2 = 0;
      }
      else
      {
	test.differential = -numeric_limits<double>::max();
	test.test_stat = -numeric_limits<double>::max();
	test.value_1 = prev.RPM;
	test.value_2 = 0;
	test.p_value = 1;
	performed_test = false;
      }
    }
    else
    {
      assert (prev.RPM == 0.0 && curr.RPM == 0.0);
      performed_test = false;
    }
    
  }
  
  test.test_status = performed_test ? OK : NOTEST;
  return test.p_value;
}


double cuffdiff(double RPM_A, double RPMvar_A, double RPM_B, double RPMvar_B){
  double pValue = 999.0;
  
    RPMContext conditionA( RPM_A, RPMvar_A );
    RPMContext conditionB( RPM_B, RPMvar_B );
    SampleDifference test;
    pValue = test_diffexp( conditionA, conditionB, test );
    
  return pValue;
}

RcppExport SEXP cuffdiff_wrapper(SEXP a, SEXP b, SEXP c, SEXP d);
#endif
