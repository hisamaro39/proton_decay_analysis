#include <math.h>
#include <stdio.h>
     
#include <gsl_math.h>
#include <gsl_min.h>
#include <gsl_monte_vegas.h>

/* 0=background,1=background_mc,2=exposure,3=efficiency */
#define NUM_DIM 4
#define MAX_MODE 10

#define CALL_GRAND    100
#define CALL_EACH  1000
//#define CALL_EACH  100000

double xl[NUM_DIM];
double xu[NUM_DIM];

double clevel = 0.9;
static double efficiency_estimate[MAX_MODE];
static double efficiency_uncertainty[MAX_MODE];
static double exposure_estimate[MAX_MODE];
static double exposure_uncertainty[MAX_MODE];
static double oversample_factor[MAX_MODE];
static double background_mc_uncertainty[MAX_MODE];
static double lgamma_number_of_cand[MAX_MODE];
static double lgamma_number_of_bkg[MAX_MODE];
static int number_of_candidates[MAX_MODE];
static int number_of_mc_background[MAX_MODE];
static int use_jeffreys_for_bg = 0;

double norm_fgrand;
int nmode, imode;
double xrate;

/**********************************************************************/
static void usage(char *name) {
    printf(
	"Usage: %s [16 numbers] [1 optional]\n"
	" 1) Number of candidates\n"
	" 2) Minimum decay rate\n"
	" 3) Maximum decay rate\n"
	" 4) Minimum efficiency\n"
	" 5) Maximum efficiency\n"
	" 6) Mean efficiency\n"
	" 7) Systematic error (absolute) in efficiency\n"
	" 8) Minimum exposure\n"
	" 9) Maximum exposure\n"
	"10) Mean exposure\n"
	"11) Systematic error (absolute) in exposure\n"
	"12) Minimum background (normalized)\n"
	"13) Maximum background (normalized)\n"
	"14) Number of backgrounds (un-normalized)\n"
	"15) Bkg/data Oversampling factor\n"
	"16) Systematic error (relative) in background\n"
	"17) optional: 1 if use Jefferey's prior for bg\n"
	"Use shell scripts for easier interface\n"
	, name);
    printf("\nOr,\n\n");
    printf(
	"Usage: %s [8 numbers] [1 optional]\n"
	" 1) Number of candidates\n"
	" 2) Mean efficiency\n"
	" 3) Systematic error (absolute) in efficiency\n"
	" 4) Mean exposure\n"
	" 5) Systematic error (absolute) in exposure\n"
	" 6) Number of backgrounds (un-normalized)\n"
	" 7) Bkg/data Oversampling factor\n"
	" 8) Systematic error (relative) in background\n"
	" 9) optional: 1 if use Jefferey's prior for bg\n"
	"Use shell scripts for easier interface\n"
	, name);
    exit(EXIT_FAILURE);
}

/*********************************************************************/
static void dump_ranges()
{
  fprintf(stderr,"mode #%d ====integral information====\n",imode);
  fprintf(stderr,"        bkg  = [%.3e %.3e]\n", xl[0],xu[0]);
  fprintf(stderr,"       bkgmc = [%.3e %.3e]\n", xl[1],xu[1]);
  fprintf(stderr,"        exp  = [%.3e %.3e]\n", xl[2],xu[2]);
  fprintf(stderr,"        eff  = [%.3e %.3e]\n", xl[3],xu[3]);
  fprintf(stderr,"-----------------------------------------------\n");
  fflush(stdout);
}

/**********************************************************************/
double feach (double *x, size_t dim, void *params)
{
    double bgmc, ex, ef;
    double e, m, rbkg, rcand, func;

    /* 0=background,1=background_mc,2=exposure,3=efficiency */

    if (background_mc_uncertainty[imode] == 0.)
      bgmc = 0.;
    else{
      bgmc = (x[1] - x[0]*oversample_factor[imode])/(background_mc_uncertainty[imode]*x[1]);
      bgmc *= 0.5*bgmc;
      bgmc += x[1];
    }

    if (exposure_uncertainty[imode] == 0.)
      ex = 0.;
    else{
      ex = (exposure_estimate[imode]-x[2])/exposure_uncertainty[imode];
      ex *= 0.5*ex;
    }

    if (efficiency_uncertainty[imode] == 0.)
      ef = 0.;
    else{
      ef = (efficiency_estimate[imode]-x[3])/efficiency_uncertainty[imode];
      ef *= 0.5*ef;
    }

    m = (xrate*x[2]*x[3]+x[0]);
    e = -bgmc - ex - ef - m;

    //    if (use_jeffreys_for_bg) rbkg = 1.0/x[0];

    rbkg  = number_of_mc_background[imode] * log(x[1]) - lgamma_number_of_bkg[imode];
    rcand = number_of_candidates[imode] * log(m) - lgamma_number_of_cand[imode];

    /*
    rbkg  = number_of_mc_background[imode] * log(x[1]) - lgamma((double)number_of_mc_background[imode] + 1.);
    rcand = number_of_candidates[imode] * log(m) - lgamma((double)number_of_candidates[imode] + 1.) ;
    */

    /*  log_e(x!) = lgamma(x+1)  */

    func = exp ( e + rcand + rbkg );

    return (func);
}

/*********************************************************************/
double inte_feach (void *params)
{
  double res = 0.;
  double err = 0.;
  double chisq = 0.;
  int status = 0;

  static int vegas_stage = 0;

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_monte_function pdf_func = { &feach, NUM_DIM, 0 };
  size_t dimension = NUM_DIM;
  size_t feach_calls = CALL_EACH;
  
  gsl_rng_env_setup ();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  //  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(NUM_DIM);
  gsl_monte_vegas_state *s;
  s  = gsl_monte_vegas_alloc(NUM_DIM);
    
  s->verbose = -1;        /* default, not print out */

  //  printf("%e %e %e %e %e\n",xl[0],xl[1],xl[2],xl[3],xl[4]);
  //  printf("%e %e %e %e %e\n",xu[0],xu[1],xu[2],xu[3],xu[4]);
  status = gsl_monte_vegas_integrate(&pdf_func, xl, xu, dimension, 
				     feach_calls, r, s, &res, &err);

  if (status)
    fprintf(stderr,"...error\n");
  else{
    //    fprintf(stderr,"int = %e +/- %e [%e, %e] %e\n",res,err,xl[0],xu[0],s->chisq);
  }
  /*
  printf ("converging...\n");

    do
      {
        gsl_monte_vegas_integrate (&pdf_func, xl, xu, dimension,
				   feach_calls/5., r, s, &res, &err);
        printf ("result = % .6f sigma = % .6f "
	        "chisq/dof = %.1f\n", res, err, s->chisq);
      }
    while (fabs (s->chisq - 1.0) > 0.5);
  */

  gsl_monte_vegas_free (s);
  return res;
}

/**********************************************************************/
double fgrand (double *x, size_t dim, void *params)
{
    double tot = 1.;

    /* 0=background,1=background_mc,2=exposure,3=efficiency */
    for(imode=0;imode<nmode;imode++){
      xl[0] = 0.000001;
      xu[0] = 10.*(number_of_mc_background[imode]+1)/oversample_factor[imode];
      xl[1] = xl[0];
      xu[1] = xu[0] * oversample_factor[imode] * 3.;
      xl[2] = exposure_estimate[imode] - 5.*exposure_uncertainty[imode];
      //      xl[2] = exposure_estimate[imode] - 10.*exposure_uncertainty[imode];
      if (xl[2] < 0.0) xl[2] = 0.000001;
      xu[2] = exposure_estimate[imode] + 5.*exposure_uncertainty[imode];
      //      xu[2] = exposure_estimate[imode] + 10.*exposure_uncertainty[imode];
      xl[3] = efficiency_estimate[imode] - 5.*efficiency_uncertainty[imode];
      //      xl[3] = efficiency_estimate[imode] - 10.*efficiency_uncertainty[imode];
      if (xl[3] < 0.0) xl[3] = 0.000001;
      xu[3] = efficiency_estimate[imode] + 5.*efficiency_uncertainty[imode];
      //      xu[3] = efficiency_estimate[imode] + 10.*efficiency_uncertainty[imode];
      if (xu[3] > 1.0) xu[3] = 1.0;
      xrate = x[0];

      /* dump parameter rages */
      //      dump_ranges();
      /*
      fprintf(stderr,"        mode :  %d\n", imode);
      fprintf(stderr,"Using ranges:  bg  = [%.3e %.3e]\n", xl[0],xu[0]);
      fprintf(stderr,"              bgmc = [%.3e %.3e]\n", xl[1],xu[1]);
      fprintf(stderr,"              exp  = [%.3e %.3e]\n", xl[2],xu[2]);
      fprintf(stderr,"              eff  = [%.3e %.3e]\n", xl[3],xu[3]);
      fprintf(stderr,"-----------------------------------------------\n");
      */
      tot *= inte_feach(0);
    }
    return (tot);
}

/*********************************************************************/
double inte_fgrand (double rate, void *params)
{
  double res = 0.;
  double err = 0.;
  double chisq = 0.;
  int status = 0;
  size_t rate_dim = 1;
  static int vegas_stage = 0;
  double xratel[1], xrateu[1];
  const gsl_rng_type *T;
  gsl_rng *r;

  xratel[0] = 0.000001;
  xrateu[0] = rate;

  gsl_monte_function pdf_func = { &fgrand, rate_dim, 0 };
  size_t fgrand_calls = CALL_GRAND;
  
  gsl_rng_env_setup ();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(rate_dim);
    
  s->verbose = -1;        /* default, not print out */

  status = gsl_monte_vegas_integrate(&pdf_func, xratel, xrateu, rate_dim, 
				     fgrand_calls, r, s, &res, &err);

  if (status)
    fprintf(stderr,"...error\n");
  else{
    //    fprintf(stderr,"int = %e +/- %e [%e, %e] %e\n",res,err,xratel[0],xrateu[0],s->chisq);
    //    fflush(stdout);
  }
  /*
  printf ("converging...\n");

    do
      {
        gsl_monte_vegas_integrate (&pdf_func, xratel, xrateu, rate_dim,
				   fgrand_calls/5., r, s, &res, &err);
        printf ("result = % .6f sigma = % .6f "
	        "chisq/dof = %.1f\n", res, err, s->chisq);
      }
    while (fabs (s->chisq - 1.0) > 0.5);
  */

  gsl_monte_vegas_free (s);
  return res;
}

/***************************************************************/
double minimize_inte_fgrand (double x, void *params)
{
  double integral,d;
  integral = inte_fgrand(x, 0);
  d = fabs ( clevel - integral/norm_fgrand);
  return d;
}

/*********************************************************************/
static void dump_args(int argc, char **argv)
{
    int i;
    for (i=0;i<argc;++i) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"[%d]\n",argc);
    fprintf(stderr,"-----------------------------------------------\n");
}

/*********************************************************************/
static parse_args_autoset(int argc, char **argv)
{
    double rate;

    /*
    if (argc == 10) {
	use_jeffreys_for_bg = atoi(argv[9]);
	--argc;
    }
    if (argc != 9) usage(argv[0]);
    */

    nmode = atoi(argv[1]);

    for(imode=0;imode<nmode;imode++){
      /* Number of candidates */
      number_of_candidates[imode] = atoi(argv[imode*8+2]);
      
      /* 0=rate, 1=bg, 2=bgmc, 3=ex, 4=ef */

      /* Background parameters */
      number_of_mc_background[imode] = atof(argv[imode*8+3]);
      oversample_factor[imode] = atof(argv[imode*8+4]);
      background_mc_uncertainty[imode] = atof(argv[imode*8+5]);
      /* Exposure parameters */
      exposure_estimate[imode] = atof(argv[imode*8+6]);
      exposure_uncertainty[imode] = atof(argv[imode*8+7])*exposure_estimate[imode];

      /* Efficiency parameters */
      efficiency_estimate[imode] = atof(argv[imode*8+8]); 
      efficiency_uncertainty[imode] = atof(argv[imode*8+9])*efficiency_estimate[imode];
      /* calculate lgamma in advance to save calculation time */
    lgamma_number_of_cand[imode] = lgamma((double)number_of_mc_background[imode] + 1.);
    lgamma_number_of_bkg[imode] = lgamma((double)number_of_candidates[imode] + 1.) ;
    }
}

/*********************************************************************/
static void dump_input()
{
  int jmode;
  for(jmode=0;jmode<nmode;jmode++){
    fprintf(stderr,"mode #%1d cand, eff, exp =  %3d, %f, %f\n", 
	    jmode, number_of_candidates[jmode], efficiency_estimate[jmode],
	    exposure_estimate[jmode]);
    fprintf(stderr,"        bkgmc, samfac = %d %f\n",
	    number_of_mc_background[jmode], oversample_factor[jmode]);
    fprintf(stderr,"        uncertainties(bkg, exp, eff) = %f %f %f\n",
	    background_mc_uncertainty[jmode], exposure_uncertainty[jmode],
	    efficiency_uncertainty[jmode]);
  }
  fprintf(stderr,"-----------------------------------------------\n");
  fflush(stdout);

}

/*********************************************************************/
/*********************************************************************/
int main(int argc, char **argv)
{
    double norm;
    double rate, tradrate;
    double min_rate, max_rate, max_rate_norm;
    double clevel,level, diff_level;
    double clfinal;
    double minimize_tol = 1.0E-2;
    double classic_rate_limit;

    int i;
    double min_min, min_tmp, min_a, min_b;
    double func_tmpl, func_tmpu;
    int status;
    int iter =0, max_iter = 100;
    const gsl_min_fminimizer_type *T;
    gsl_min_fminimizer *s;
    double a = 0.001, b;
    double m;
    gsl_function my_minimize;

    my_minimize.function = &minimize_inte_fgrand;
    my_minimize.params = 0;

    /* dump command line*/
    dump_args(argc,argv);

    /* input parameters */
    parse_args_autoset (argc, argv);

    dump_input();

    min_rate = 0.00001;
    max_rate = 4.;
    //    max_rate = 6.1;
    max_rate_norm = max_rate*10.;

    /* Calculate the norm, then drop down for speed */
    norm = inte_fgrand(max_rate_norm, 0);
    fprintf(stderr,"Norm factor = %e\n",norm);
    fprintf(stderr,"=============================================\n");

    T = gsl_min_fminimizer_brent;
    s = gsl_min_fminimizer_alloc (T); 

    norm_fgrand = norm;
    a = min_rate;
    b = max_rate_norm;
    m = 1.;
    b = max_rate;

    /*
    norm_fgrand = norm;
    b = max_rate_norm;
    min_min = min_a;
    m = a;
    if (min_b <= min_a){
      min_min = min_b;
      m = b;
    }
    */
    min_a = minimize_inte_fgrand(a,0);
    min_b = minimize_inte_fgrand(b,0);
    if (min_b <= min_a){
      min_min = min_b;
      m = b;
    }

    for(i=1;i<20;i++){
	min_tmp = minimize_inte_fgrand(b*(float)i/20.,0);
	printf("%f %e %e\n",b*(float)i/20.,min_tmp, inte_fgrand(b*(float)i/20.,0));
	if (min_tmp<=min_min){
	    min_min = min_tmp;
	    m = b*(float)i/20.;
	}	
	fflush(stdout);
    }

    gsl_min_fminimizer_set (s, &my_minimize, m, a, b);

    do
      {
	iter++;
	status = gsl_min_fminimizer_iterate (s);
	m = gsl_min_fminimizer_minimum (s);
	//	m = gsl_min_fminimizer_x_minimum (s);
	a = gsl_min_fminimizer_x_lower (s);
	b = gsl_min_fminimizer_x_upper (s);

	//status = gsl_min_test_interval (a, b, 0.001, 0.0);
	status = gsl_min_test_interval (a, b, 0.00001, 0.0);
	if (status == GSL_SUCCESS)
	  fprintf (stderr,"Converged:\n");

	fprintf (stderr,"%5d [%e, %e] %e\n", iter, a, b, m);
	fflush(stdout);
      }
    while (status == GSL_CONTINUE && iter < max_iter);

    rate = m;
    diff_level = minimize_inte_fgrand(m,0);
    fprintf(stderr,"=============================================\n");
    fprintf(stderr,"Result with systematic uncertainties:\n");
    fprintf(stderr,"     rate = %f\n",rate);
    fprintf(stderr,"    limit = %f *10^33 year\n",1.0/rate);
    fprintf(stderr," |0.9-CL| = %e\n",diff_level);
    fprintf(stderr,"=============================================\n");

    gsl_min_fminimizer_free (s);

    fflush(stdout);
    return 0;
}
