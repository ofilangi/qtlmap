#include "nlopt.h"
#include <stdio.h>
#include <math.h>

#ifdef MANAGE_OMP
#include <omp.h>
#endif



#ifdef AIX_SYS
    void get_nb_opti_nlopt_qtlmap(int *numOptim) {
#else
	void get_nb_opti_nlopt_qtlmap_(int *numOptim) {
#endif
	*numOptim=NLOPT_NUM_ALGORITHMS;
}

#ifdef AIX_SYS
    void info_nlopt_qtlmap(char list_optim[NLOPT_NUM_ALGORITHMS][200]) {
#else
	void info_nlopt_qtlmap_(char list_optim[NLOPT_NUM_ALGORITHMS][200]) {
#endif
	int i;
	for (i = 0; i < NLOPT_NUM_ALGORITHMS; ++i) {
		strcpy(list_optim[i],nlopt_algorithm_name((nlopt_algorithm) i));
	}
}


#ifdef AIX_SYS
    void interface_minimizing_nlopt_partial(
#else
    void interface_minimizing_nlopt_partial_(
#endif
		                         int * opti_nlopt,
		                         void (*funct_value)(),
		                         void (*funct_gradient)(),
		                         void (*funct_partial)(),
		                         int *n,
		                         double *x,
		                         const double *lb,
		                         const double *ub,
		                         int    *s1,
		                         int    *iuser,
		                         int    *s2,
		                         double *user,
		                         int    *np,
		                         double *fp,
		                         double *filter_fp,
		                         int    *nm,
		                         double *fm,
                                 double *filter_fm,
		                         double *minf,
		                         double *xtol_rel,
		                         double *ftol_rel,
                                 int * maxeval,
                                 double * maxtime,
                                 int *ifail)
{

	double my_funct_nlop(int n, const double *x, double *grad, void *data)
	{
		double f = 0.0;
	    int i;


	    /*
		printf("---->my_funct_nlop\n");
		printf(" n:%d\n s1=%d\n s2=%d\n np=%d\n nm=%d\n",n,*s1,*s2,*np,*nm);
	    printf("  \n************ x ******************\n");
	    for (i=0;i<n;i++) { printf("x[%d]=%f,",i,x[i]); }
	    printf("  \n************ lb ******************\n");
	    for (i=0;i<n;i++) { printf("lb[%d]=%f,",i,lb[i]); }
	    printf("  \n************ ub ******************\n");
	    for (i=0;i<n;i++) { printf("ub[%d]=%f,",i,ub[i]); }
	    printf("  \n************ fp ******************\n");
	    for (i=0;i<*np;i++) { printf("fp[%d]=%f,",i,fp[i]); }
	    printf("  \n************ fm ******************\n");
	    for (i=0;i<*nm;i++) { printf("fm[%d]=%f,",i,fm[i]); }
*/

	    funct_value(funct_partial,&n,x,&f,s1,iuser,s2,user,np,fp,filter_fp,nm,fm,filter_fm);
#ifdef MANAGE_OMP
	   /* printf("---->my_funct_nlop:%i --  f:%f\n",omp_get_thread_num(),f);*/
#endif
/*	    printf("\n==>%f\n",f); */
	    if (grad) {
	    	/*printf("---->funct_gradient");*/
	    	funct_gradient(&f,&n,x,funct_partial,s1,iuser,s2,user,ub,lb,grad,np,fp,filter_fp,nm,fm,filter_fm);
	    }
	    return (f);
	}

	void *f_data = NULL ;
	double ftol_abs = 0 ;
	double minf_max = -HUGE_VAL;
	double *xtol_abs = NULL ;

	sub_minimizing_nlopt((*opti_nlopt)-1,*n,x,my_funct_nlop,f_data,lb,ub,minf,minf_max,(*ftol_rel),
			             ftol_abs,(*xtol_rel),xtol_abs,*maxeval,*maxtime,ifail);

/*	printf("\n==>%f\n",*minf); */
}

void sub_minimizing_nlopt(  int opti_nlopt,
		                    int n,
			                double *x,
			                nlopt_func f,
			                void *f_data,
			                const double *lb,
			                const double *ub,
			                double *minf,
			                double minf_max,
			                double ftol_rel,
			                double ftol_abs,
			                double xtol_rel,
			                const double *xtol_abs,
			                int maxeval,
			                double maxtime,
			                int *ifail)
{
    nlopt_algorithm algorithm = opti_nlopt ;
 /*   printf("NLOPT :%d\n",opti_nlopt); */
   /*
    int i =0;
       printf("N :%d\n",n);

       printf("  \n************ x ******************\n");
       for (i=0;i<n;i++) { printf("x[%d]=%f,",i,x[i]); }
       printf("  \n************ lb ******************\n");
       for (i=0;i<n;i++) { printf("lb[%d]=%f,",i,lb[i]); }
       printf("  \n************ ub ******************\n");
       for (i=0;i<n;i++) { printf("ub[%d]=%f,",i,ub[i]); }

       printf("\n****>MINF :%f",*minf);*/
    nlopt_result r = nlopt_minimize(algorithm,n,f,f_data,lb,ub,x,minf,minf_max,ftol_rel,
    		                        ftol_abs,xtol_rel,xtol_abs,maxeval,maxtime);

    switch(r) {
      case NLOPT_FAILURE  : *ifail=-4;break;  /* generic failure code */
      case NLOPT_INVALID_ARGS : *ifail=-3;break;
      case NLOPT_OUT_OF_MEMORY : *ifail=-2;break;
      case NLOPT_ROUNDOFF_LIMITED : *ifail=-1;break;

      case NLOPT_SUCCESS : *ifail=0;break;
      case NLOPT_MINF_MAX_REACHED: *ifail=1;break;
      case NLOPT_FTOL_REACHED : *ifail=2;break;
      case NLOPT_XTOL_REACHED : *ifail=3;break;
      case NLOPT_MAXEVAL_REACHED : *ifail=4;break;
      case NLOPT_MAXTIME_REACHED : *ifail=5;break;
    }

}



