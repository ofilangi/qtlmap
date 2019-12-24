#ifndef CONST_GPU_IMPL_H
#define CONST_GPU_IMPL_H

#define safecall(call) do{\
		cudaError_t errERR = call ; \
		if (cudaSuccess != errERR){\
			fprintf(stderr, "cuda error at %s:%d, %s\n",\
					__FILE__, __LINE__, cudaGetErrorString(errERR));exit(1);\
		}\
}while(0)

#define PRECISION_FLOTTANT 4
/* ------------------------------- ACCESSORS ------------------------------------------ */
#define OPTIMIZE_MEMORY_MATSYM true

#if(!OPTIMIZE_MEMORY_MATSYM)

#define BASE_SCALE_SYMMAT data.ntniv*data.npos*data.ntniv
#define GETACC(data,isim,ipos,iniv,jniv) (iniv>=jniv?isim*BASE_SCALE_SYMMAT+ipos+iniv*data.npos+jniv*data.ntniv*data.npos:isim*BASE_SCALE_SYMMAT+ipos+jniv*data.npos+iniv*data.ntniv*data.npos)
#define SIZEMATSYM(nsim,npos,ntniv) (nsim*npos*ntniv*ntniv)
#else

#define CEIL(VARIABLE) ( (VARIABLE - (int)VARIABLE)==0 ? (int)VARIABLE : (int)VARIABLE+1 )
#define BASE_SCALE_SYMMAT ((data.npos*data.ntniv*(data.ntniv+1)/2))
#define GETACC(data,isim,ipos,iniv,jniv) ( iniv>jniv?ipos+(iniv*(iniv+1)/2+jniv)*data.npos+isim*BASE_SCALE_SYMMAT:ipos+((jniv*(jniv+1))/2+iniv)*data.npos+isim*BASE_SCALE_SYMMAT )
#define SIZEMATSYM(nsim,npos,ntniv) (nsim*npos*( ntniv*(ntniv+1) /2 ))

#endif



#ifdef CUDA_SP
#define DT float
#else
#define DT double
#endif


/* Generate print information about sbu-results */
#if ! defined(_CUDA_HOST_DEBUG_)
#define _CUDA_HOST_DEBUG_ 0
#endif

/* Generate time comsuming information about sub process */
#if ! defined(_CUDA_HOST_TIME_PROF_)
#define _CUDA_HOST_TIME_PROF_ 0
#endif
/* number of thread available by block */
#define MAX_BLOCKDIM_512 512
#define MAX_BLOCKDIM_64 64
#define MAX_BLOCKDIM_32 32
#define BLOCKDIMX_XTX_C 96
#define BLOCKDIMX 128
//#define MAXND_WORK 5000
#define MAXNP_WORK 80


#define MAXNTNIV_LOCALTHREAD_MEMORY 200

/* Definition des donnees en memoire constante sur gpu */

__constant__ int constNbLevelFix ;
__constant__ int constNbLevelVar ;
__constant__ int constNSim;
//__constant__ int constCorrIndexCol[MAXNTNIV_LOCALTHREAD_MEMORY];
__constant__ int constSizeFam[MAXNP_WORK];

#endif
