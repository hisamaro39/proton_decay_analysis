#include <math.h>
#include <stdio.h>
     
#include <gsl_math.h>
#include <gsl_min.h>
#include <gsl_monte_vegas.h>

#define NUM_DIM 5

#define CALL_GRAND  1000000
#define LRGMCSTAT

/* 0=rate,1=background,2=background_mc,3=exposure,4=efficiency */
double xl[NUM_DIM];
double xu[NUM_DIM];

double clevel = 0.9;
static double efficiency_estimate;
static double efficiency_uncertainty;
static double efficiency_uncertainty_low;
static double exposure_estimate;
static double exposure_uncertainty;
static double oversample_factor;
static double background_mc_uncertainty_fraction;
static int number_of_candidates;
//static int number_of_mc_background;
static double number_of_mc_background;

static double efficiency_estimate2;
static double efficiency_uncertainty2;
static double efficiency_uncertainty_low2;
static double exposure_estimate2;
static double exposure_uncertainty2;
static double oversample_factor2;
static double background_mc_uncertainty_fraction2;
static int number_of_candidates2;
static double number_of_mc_background2;

static double x1_lim[2];
static double x3_lim[2];
static double x4_lim[2];

static int use_jeffreys_for_bg = 0;
double norm_fgrand;

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

/**********************************************************************/
static double ipow(double x, int n) 
{
  /* Calculate X to the N'th with N an int. */
    int i;
    double ans = 1.0;
    for (i = 0; i < n; ++i) ans *= x;
    return ans;
}

/**********************************************************************/
double fgrand (double *x, size_t dim, void *params)
{
    double bgmc, ex, ef;
    double e, m, rbkg, rcand, func;
    double x1_2, x3_2, x4_2, rcand2, m2;
    /*
    x1_2 = number_of_mc_background2/oversample_factor2 + 
	(x[1]-number_of_mc_background/oversample_factor)
	*(background_mc_uncertainty_fraction2*number_of_mc_background2)
	/(background_mc_uncertainty_fraction*number_of_mc_background);
    */
    x1_2 = number_of_mc_background2/oversample_factor2 + 
	(x[1]-number_of_mc_background/oversample_factor)
	*(background_mc_uncertainty_fraction2*number_of_mc_background2)
	/(background_mc_uncertainty_fraction*number_of_mc_background)
	*oversample_factor/oversample_factor2;
    x3_2 = exposure_estimate2 
	+ (x[3]-exposure_estimate)*exposure_uncertainty2/exposure_uncertainty;
    if(x[4]>efficiency_estimate){
	x4_2 = efficiency_estimate2
	    + (x[4]-efficiency_estimate)*efficiency_uncertainty2/efficiency_uncertainty;
    }
    else{
	x4_2 = efficiency_estimate2
	    + (x[4]-efficiency_estimate)
	    *efficiency_uncertainty_low2/efficiency_uncertainty_low;
    }

    /*
    printf("1: %lf %lf %lf %lf\n", exposure_estimate, x[1], x[3]-exposure_estimate, x[4]-efficiency_estimate);
    printf("2: %lf %lf %lf %lf\n", exposure_estimate2, x1_2, x3_2-exposure_estimate2, x4_2-efficiency_estimate2);
    */

    if(x1_2<x1_lim[0]) x1_2 = x1_lim[0];
    if(x1_2>x1_lim[1]) x1_2 = x1_lim[1];
    if(x3_2<x3_lim[0]) x3_2 = x3_lim[0];
    if(x3_2>x3_lim[1]) x3_2 = x3_lim[1];
    if(x4_2<x4_lim[0]) x4_2 = x4_lim[0];
    if(x4_2>x4_lim[1]) x4_2 = x4_lim[1];

    //printf("3: %lf %lf %lf %lf\n", exposure_estimate2, x1_2, x3_2-exposure_estimate2, x4_2-efficiency_estimate2);

    /* 0=rate,1=background,2=background_mc,3=exposure,4=efficiency */

    if (background_mc_uncertainty_fraction == 0.)
      bgmc = 0.;
    else{
#ifdef LRGMCSTAT
      bgmc = (number_of_mc_background - x[1]*oversample_factor)/(background_mc_uncertainty_fraction*number_of_mc_background);
      bgmc *= 0.5*bgmc;
#else
      bgmc = (x[2] - x[1]*oversample_factor)/(background_mc_uncertainty_fraction*x[2]);
      bgmc *= 0.5*bgmc;
      bgmc += x[2];
#endif
    }

    if (exposure_uncertainty == 0.)
      ex = 0.;
    else{
      ex = (exposure_estimate-x[3])/exposure_uncertainty;
      ex *= 0.5*ex;
    }

    if (efficiency_uncertainty == 0.)
      ef = 0.;
    else if(x[4]-efficiency_estimate>=0.){
      ef = (efficiency_estimate-x[4])/efficiency_uncertainty;
      ef *= 0.5*ef;
    }
    else{
      ef = (efficiency_estimate-x[4])/efficiency_uncertainty_low;
      ef *= 0.5*ef;
    }

    m = (x[0]*x[3]*x[4]+x[1]);
    e = -bgmc - ex - ef - m;

    m2 = x[0]*x3_2*x4_2+x1_2;
    
    if (use_jeffreys_for_bg) rbkg = 1.0/x[1];
#ifdef LRGMCSTAT
    rbkg  = 0.;
#else
    rbkg  = number_of_mc_background * log(x[2]) - lgamma((double)number_of_mc_background + 1.);
#endif
    rcand = number_of_candidates * log(m) - lgamma((double)number_of_candidates + 1.) ;

    rcand2 = number_of_candidates2 * log(m2) - lgamma((double)number_of_candidates2 + 1.) ;

    //printf("%lf %lf %lf\n", rcand2, m, m2);

    /*  log_e(x!) = lgamma(x+1)  */

    func = exp ( e + rcand + rbkg + rcand2 - m2);
    //func = exp ( e + rcand + rbkg + rcand2 - m2 -bgmc-ex-ef+rbkg);

    if(func>0. || func < 1.){
    }
    else{
	printf("invalid :%lf\n", func);
	printf("%lf %lf %lf\n", x1_2, x3_2, x4_2);
    }

    return (func);
}

/*********************************************************************/
double inte_fgrand (double rate, void *params)
{
  double res = 0.;
  double err = 0.;
  double chisq = 0.;
  int status = 0;
  /* double s; */
  gsl_monte_vegas_state *s ;
  static int vegas_stage = 0;

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_monte_function pdf_func = { &fgrand, NUM_DIM, 0 };
  size_t dimension = NUM_DIM;
  size_t fgrand_calls = CALL_GRAND;
  
  gsl_rng_env_setup ();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  xu[0] = rate;
  /* gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(NUM_DIM); */
  s = gsl_monte_vegas_alloc(NUM_DIM); 
    
  s->verbose = -1;        /* default, not print out */

  //  printf("%e %e %e %e %e\n",xl[0],xl[1],xl[2],xl[3],xl[4]);
  //  printf("%e %e %e %e %e\n",xu[0],xu[1],xu[2],xu[3],xu[4]);
  status = gsl_monte_vegas_integrate(&pdf_func, xl, xu, dimension, 
				     fgrand_calls, r, s, &res, &err);

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
  printf("%f %e %e\n",x,d, integral);
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
static void parse_args_handset(int argc, char **argv)
{
    double rate;

    if (argc == 18) {
	use_jeffreys_for_bg = atoi(argv[17]);
	--argc;
    }
    if (argc != 17) usage(argv[0]);

    /* Number of candidates */
    number_of_candidates = atoi(argv[1]);

    /* 0=rate, 1=bg, 2=bgmc, 3=ex, 4=ef */

    /* Decay Rate parameters */
    xl[0] = atof(argv[2]);
    xu[0] = atof(argv[3]);

    /* Efficiency parameters */
    xl[4] = atof(argv[4]);
    xu[4] = atof(argv[5]);
    efficiency_estimate = atof(argv[6]); 
    efficiency_uncertainty = atof(argv[7]);

    /* Exposure parameters */
    xl[3] = atof(argv[8]);
    xu[3] = atof(argv[9]);
    exposure_estimate = atof(argv[10]);
    exposure_uncertainty = atof(argv[11]);

    /* Background parameters */
    xl[1] = atof(argv[12]);
    xu[1] = atof(argv[13]);
    number_of_mc_background = atof(argv[14]);
    oversample_factor = atof(argv[15]);
    background_mc_uncertainty_fraction = atof(argv[16]);
    xl[2] = xl[1];
    xu[2] = xu[1] * oversample_factor * 3;
}

/*********************************************************************/
static parse_args_autoset(int argc, char **argv)
{
    double rate;

    if (argc == 10) {
	use_jeffreys_for_bg = atoi(argv[9]);
	--argc;
    }
    if (argc!=17 && argc!=19) usage(argv[0]);

    /* Number of candidates */
    number_of_candidates = atoi(argv[1]);

    /* 0=rate, 1=bg, 2=bgmc, 3=ex, 4=ef */

    /* Decay Rate parameters */
    xl[0] = 0.000001;
    xu[0] = 1000.;

    /* Background parameters */
    number_of_mc_background = atof(argv[6]);
    oversample_factor = atof(argv[7]);
    background_mc_uncertainty_fraction = atof(argv[8]);
    xl[1] = 0.000001;
    xu[1] = 10.*(number_of_mc_background+1)/oversample_factor;
    xl[2] = xl[1];
    xu[2] = xu[1] * oversample_factor * 3.;

    /* Exposure parameters */
    exposure_estimate = atof(argv[4]);
    exposure_uncertainty = atof(argv[5])*exposure_estimate;
    xl[3] = exposure_estimate - 10*exposure_uncertainty;
    if (xl[3] < 0.0) xl[3] = 0.000001;
    xu[3] = exposure_estimate + 10*exposure_uncertainty;

    /* Efficiency parameters */
    efficiency_estimate = atof(argv[2]); 
    efficiency_uncertainty = atof(argv[3])*efficiency_estimate;
    xl[4] = efficiency_estimate - 10*efficiency_uncertainty;
    if (xl[4] < 0.00001) xl[4] = 0.00001;
    xu[4] = efficiency_estimate + 10*efficiency_uncertainty;
    //xu[4] = exposure_estimate + 10*exposure_uncertainty;
    if (xu[4] > 1.0) xu[4] = 1.0;

    /* SK-II parameters*/
    number_of_candidates2 = atoi(argv[9]);
    efficiency_estimate2 = atof(argv[10]); 
    efficiency_uncertainty2 = atof(argv[11])*efficiency_estimate2;
    exposure_estimate2 = atof(argv[12]);
    exposure_uncertainty2 = atof(argv[13])*exposure_estimate2;
    number_of_mc_background2 = atof(argv[14]);
    oversample_factor2 = atof(argv[15]);
    background_mc_uncertainty_fraction2 = atof(argv[16]);
    x1_lim[0] = 0.000001;
    x1_lim[1] = 10.*(number_of_mc_background2+1)/oversample_factor2;
    x3_lim[0] = exposure_estimate2 - 10.*exposure_uncertainty2;
    if(x3_lim[0]<0.000001) x3_lim[0] = 0.000001;
    x3_lim[1] = exposure_estimate2 + 10.*exposure_uncertainty2;
    x4_lim[0] = efficiency_estimate2 - 10.*efficiency_uncertainty2;
    if(x4_lim[0]<0.00001) x4_lim[0] = 0.00001;
    x4_lim[1] = efficiency_estimate2 + 10.*efficiency_uncertainty2;
    if(x4_lim[1]>1.) x4_lim[1] = 1.;

    /* asymmetric error for detection efficiency */
    if(argc==19){
	efficiency_uncertainty_low  = atof(argv[17])*efficiency_estimate;
	efficiency_uncertainty_low2 = atof(argv[18])*efficiency_estimate2;
	xl[4] = efficiency_estimate - 10.*efficiency_uncertainty_low;
	if (xl[4]   < 0.00001) xl[4]     = 0.00001;
	x4_lim[0] = efficiency_estimate2 - 10.*efficiency_uncertainty_low2;
	if(x4_lim[0]< 0.00001) x4_lim[0] = 0.00001;
    }else{
	efficiency_uncertainty_low  = efficiency_uncertainty;
	efficiency_uncertainty_low2 = efficiency_uncertainty2;
    }
    
}

/*********************************************************************/
static void dump_ranges(double max_rate)
{
    fprintf(stderr,"Using ranges: rate = [%.3e %.3e]\n", xl[0],max_rate);
    fprintf(stderr,"               bg  = [%.3e %.3e]\n", xl[1],xu[1]);
    fprintf(stderr,"              bgmc = [%.3e %.3e]\n", xl[2],xu[2]);
    fprintf(stderr,"              exp  = [%.3e %.3e]\n", xl[3],xu[3]);
    fprintf(stderr,"              eff  = [%.3e %.3e]\n", xl[4],xu[4]);
    fprintf(stderr,"-----------------------------------------------\n");
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

    parse_args_autoset (argc, argv);

    /*
    level = poisson_limit(clevel, number_of_candidates, 
                          number_of_mc_background/oversample_factor, 
                          &classic_rate_limit );

    rate = classic_rate_limit/(exposure_estimate * efficiency_estimate);
    xu[0] = 10. * rate;

    fprintf(stderr,"=============================================\n");
    fprintf(stderr,"Result without systematic uncertainties:\n");
    fprintf(stderr,"     rate = %f\n",rate);
    fprintf(stderr,"    limit = %f *10^33 year\n",1.0/rate);
    fprintf(stderr," |0.9-CL| = %e\n",fabs(level - clevel));
    fprintf(stderr,"=============================================\n");
    */
    min_rate = xl[0] + 0.00001;	/* get error if use exactly min */
    max_rate = xu[0];
    max_rate = 100.;
    max_rate = 100.;
    //max_rate = 2.;
    //    max_rate = 6.1;
    //max_rate_norm = max_rate*10.;
    max_rate_norm = 100.;
    /* dump parameter rages */
    dump_ranges(max_rate_norm);

    fprintf(stderr,"# of grand call = %d\n",CALL_GRAND);

    /* Calculate the norm, then drop down for speed */
    norm = inte_fgrand(max_rate_norm, 0);
    fprintf(stderr,"Norm factor = %e\n",norm);
    fprintf(stderr,"=============================================\n");

    T = gsl_min_fminimizer_brent;
    s = gsl_min_fminimizer_alloc (T); 

    norm_fgrand = norm;
    b = max_rate_norm;
    b = max_rate;
    m = b * 2./3.;
    min_a = minimize_inte_fgrand(a,0);
    min_min = min_a;
    m = a;
    min_b = minimize_inte_fgrand(b,0);
    if (min_b <= min_a){
      min_min = min_b;
      m = b;
    }
    for(i=1;i<20;i++){
      min_tmp = minimize_inte_fgrand(b*(float)i/20.,0);
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

	fprintf (stderr,"%5d :", iter);
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
