#include <math.h>
#include <stdio.h>
     
#include "gsl_math.h"
#include "gsl_min.h"
#include "gsl_monte_vegas.h"

#define NUM_DIM 5

//#define CALL_GRAND  1000000
#define CALL_GRAND  1000000
//#define POISSONLIMIT

/* 0=rate,1=background,2=background_mc,3=exposure,4=efficiency */
double xl[NUM_DIM];
double xu[NUM_DIM];

double clevel = 0.9;
static double efficiency_estimate;
static double efficiency_uncertainty;
static double exposure_estimate;
static double exposure_uncertainty;
static double oversample_factor;
static double background_mc_uncertainty_fraction;
static int number_of_candidates;
//static int number_of_mc_background;
static double number_of_mc_background;
static int use_baisian = 0;
double norm_fgrand;

/**********************************************************************/
static void usage(char *name) {
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
	" 9) Function: 0->poisson 1->baisian\n"
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
  //fprintf(stderr,"fgrand\n");
  //fprintf(stderr,"dim=%d\n",dim);
  //fprintf(stderr,"x %f/%f/%f/%f/%f\n",x[0],x[1],x[2],x[3],x[4]);
    double bgmc, ex, ef;
    double e, m, rbkg, rcand, func;

    /* 0=rate,1=background,2=background_mc,3=exposure,4=efficiency */

    if (background_mc_uncertainty_fraction == 0.)
      bgmc = 0.;
    else{
      bgmc = (x[2] - x[1]*oversample_factor)/(background_mc_uncertainty_fraction*x[2]);
      bgmc *= 0.5*bgmc;
      bgmc += x[2];
    }

    if (exposure_uncertainty == 0.)
      ex = 0.;
    else{
      ex = (exposure_estimate-x[3])/exposure_uncertainty;
      ex *= 0.5*ex;
    }

    if (efficiency_uncertainty == 0.)
      ef = 0.;
    else{
      ef = (efficiency_estimate-x[4])/efficiency_uncertainty;
      ef *= 0.5*ef;
    }

    m = (x[0]*x[3]*x[4]+x[1]);
    e = -bgmc - ex - ef - m;

    rbkg  = number_of_mc_background * log(x[2]) - lgamma((double)number_of_mc_background + 1.);
    rcand = number_of_candidates * log(m) - lgamma((double)number_of_candidates + 1.) ;

    /*  log_e(x!) = lgamma(x+1)  */

    func = exp ( e + rcand + rbkg );
    //fprintf(stderr,"func=%f\n",func);

    return (func);
}

/**********************************************************************/
double fgrand_pois (double *x, size_t dim, void *params)
{
  //fprintf(stderr,"fgrand_pois\n");
  //fprintf(stderr,"dim=%d\n",dim);
  //fprintf(stderr,"x %f/%f/%f/%f/%f\n",x[0],x[1],x[2],x[3],x[4]);
    double bgmc, ex, ef;
    double e, m, rbkg, rcand, func;

    /* 0=rate,1=background,2=background_mc,3=exposure,4=efficiency */

    m = (x[0]*exposure_estimate*efficiency_estimate+number_of_mc_background/oversample_factor);

    rcand = number_of_candidates * log(m) - lgamma((double)number_of_candidates + 1.) ;

    /*  log_e(x!) = lgamma(x+1)  */

    func = exp ( -m + rcand );

    return (func);
}

/*********************************************************************/
double inte_fgrand (double rate, void *params)
{
  //fprintf(stderr,"inte_fgrand\n");
  //fprintf(stderr,"rate=%f\n",rate);
  double res = 0.;
  double err = 0.;
  double chisq = 0.;
  int status = 0;
  /* double s; */
  gsl_monte_vegas_state *s ;
  static int vegas_stage = 0;

  const gsl_rng_type *T;
  gsl_rng *r;

#ifdef POISSONLIMIT
  //fprintf(stderr,"define pdf_func POISSONLIMIT\n");
  gsl_monte_function pdf_func = { &fgrand_pois, NUM_DIM, 0 };
#else
  //fprintf(stderr,"define pdf_func\n");
  gsl_monte_function pdf_func = { &fgrand, NUM_DIM, 0 };
#endif
  size_t dimension = NUM_DIM;
  size_t fgrand_calls = CALL_GRAND;
  
  gsl_rng_env_setup ();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  xu[0] = rate;
  /* gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(NUM_DIM); */
  s = gsl_monte_vegas_alloc(NUM_DIM); 
    
  s->verbose = -1;        /* default, not print out */

  //fprintf(stderr,"%e %e %e %e %e\n",xl[0],xl[1],xl[2],xl[3],xl[4]);
  //fprintf(stderr,"%e %e %e %e %e\n",xu[0],xu[1],xu[2],xu[3],xu[4]);
  //fprintf(stderr,"gsl_monte_vegas_integrate\n");
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
  //fprintf(stderr,"minimize_inte_fgrand\n");
  double integral,d;
  integral = inte_fgrand(x, 0);
  d = fabs ( clevel - integral/norm_fgrand);
  //fprintf(stderr,"rate 0.9-this integral %f %e %e\n",x,d, integral);

  return d;
}

/*********************************************************************/
static void dump_args(int argc, char **argv)
{
    fprintf(stderr,"dump_args\n");
    fprintf(stderr,"argc=%d\n",argc);
    int i;
    for (i=0;i<argc;++i) fprintf(stderr,"%s ",argv[i]);
    fprintf(stderr,"[%d]\n",argc);
    fprintf(stderr,"-----------------------------------------------\n");
}


/*********************************************************************/
static parse_args_autoset(int argc, char **argv)
{
  printf("parse_args_autoset\n");
    double rate;

    if (argc == 10) {
	use_baisian = atoi(argv[9]);
	--argc;
    }
    if (argc != 9) usage(argv[0]);

    /* Number of candidates */
    number_of_candidates = atoi(argv[1]);

    /* 0=rate, 1=bg, 2=bgmc, 3=ex, 4=ef */

    /* Decay Rate parameters */
    xl[0] = 0.000001;
    //xu[0] = 10. * rate;
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
    //if (xl[4] < 0.0) xl[4] = 0.0001;
    xu[4] = efficiency_estimate + 10*efficiency_uncertainty;
    //xu[4] = exposure_estimate + 10*exposure_uncertainty;
    if (xu[4] > 1.0) xu[4] = 1.0;
}

/*********************************************************************/
static void dump_ranges(double max_rate)
{
  printf("dump_ranges\n");
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
    double clevel,level, diff_level, tmpclevel;
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

    /* input parameters and calculate limit by traditional method */
    parse_args_autoset (argc, argv);

    min_rate = xl[0] + 0.00001;	/* get error if use exactly min */
    max_rate = xu[0];
    //max_rate = 100.;
    //max_rate = 10.;
    max_rate = 1.;
    //max_rate = 2.;
    //    max_rate = 6.1;
    //max_rate_norm = max_rate*10;
    max_rate_norm = 100.;
    /* dump parameter rages */
    dump_ranges(max_rate_norm);

    //####################### make pdf plot
    double res = 0.;
    double err = 0.;
    double chisq = 0.;
    /* double s; */
    gsl_monte_vegas_state *s_tmp ;
    static int vegas_stage = 0;

    const gsl_rng_type *T_tmp;
    gsl_rng *r_tmp;

    char table_name[256];
    #ifdef POISSONLIMIT
    gsl_monte_function pdf_func = { &fgrand_pois, NUM_DIM, 0 };
    sprintf(table_name,"output_table/pois_ncand%d_eff%g_syseff%g_nbkg%g_sysbkg%g_factor%g.txt"
        ,number_of_candidates,efficiency_estimate,atof(argv[3]),number_of_mc_background,background_mc_uncertainty_fraction,oversample_factor);
    #else
    gsl_monte_function pdf_func = { &fgrand, NUM_DIM, 0 };
    sprintf(table_name,"output_table/bais_ncand%d_eff%g_syseff%g_nbkg%g_sysbkg%g_factor%g.txt"
        ,number_of_candidates,efficiency_estimate,atof(argv[3]),number_of_mc_background,background_mc_uncertainty_fraction,oversample_factor);
    #endif

    size_t dimension = NUM_DIM;
    size_t fgrand_calls = CALL_GRAND;

    gsl_rng_env_setup ();
    T_tmp = gsl_rng_default;
    r_tmp = gsl_rng_alloc (T_tmp);

    /* gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(NUM_DIM); */
    s_tmp = gsl_monte_vegas_alloc(NUM_DIM); 

    s_tmp->verbose = -1;        /* default, not print out */

    fprintf(stderr,"table_name=%s\n",table_name);
    FILE *file;
    file = fopen(table_name,"w");
    int t;
    for(t=0;t<500;t++){
      if(t%100==0) fprintf(stderr,"%d iteration finished!!\n",t);
      xu[0] = 0.002*(t + 1);
      xl[0] = 0.002*t;
      float mean = (xl[0]+xu[0])*0.5;
      status = gsl_monte_vegas_integrate(&pdf_func, xl, xu, dimension, 
          fgrand_calls, r_tmp, s_tmp, &res, &err);
      //fprintf(stderr,"xl/xu/res=%f/%f/%f\n",xl[0],xu[0],res);
      fprintf(file,"%f %f %f\n",xl[0],xu[0],res);
    }
    fclose(file);

    gsl_monte_vegas_free (s_tmp);

    //############################ make pdf plot end

    return 0;
}
