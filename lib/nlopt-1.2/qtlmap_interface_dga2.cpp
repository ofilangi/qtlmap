

extern "C" void get_nb_opti_nlopt_qtlmap(int *numOptim) {
	*numOptim=0;
}

extern "C" void info_nlopt_qtlmap(char list_optim[10][200]) {
}


extern "C" void interface_minimizing_nlopt_partial(
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

   *ifail = -4;
}



