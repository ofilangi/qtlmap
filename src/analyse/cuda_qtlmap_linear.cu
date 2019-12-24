/*
 *
 * O.Filangi : linearisation de la vraissemblance -  calcul simultané des positions sur le génome
 * -deviceemu => emulation
 * nvcc -g -arch=sm_20 --cuda analyse_whole-genome.cu
 * debugging : ddd --debugger /usr/local/cuda/bin/cuda-gdb qtlmap
 * Debug   : nvcc -g -G -gencode=arch=compute_20,code=\"sm_20,compute_20\" --cuda analyse_whole-genome.cu
 * Release : nvcc -O4 --use_fast_math  -gencode=arch=compute_20,code=\"sm_20,compute_20\" --cuda analyse_whole-genome.cu
 * 
 * 
 * TODO : XX est une matrice symetrique => utiliser que la moitie superieur de la matrice
 * 
 * pour le profiling : /usr/local/cuda/computeprof/bin/computeprof
 */
#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <iostream>
using namespace std;

#include "const_gpu_impl.h"

/* Keep in memory the residual variance fitted for each Hypothesis analysis */
static DT * allVarFitted = NULL ;

void print_info_memory() {
	size_t free,total;
	cudaMemGetInfo(&free, &total) ;
	printf("===============> free: %f Mb (%lu bytes)   total:%f Mb (%lu bytes) \n",double(free)/(1024*1024),free,double(total)/(1024*1024),total);
}

class  QTLMapStructDeviceDataLinear {
public:
	//DT * contingence ; /* matrice de contingence les effets statiques sont dans les premieres colonnes. */
	DT * contingence_fix ; /*                                      [NDMAX][NFIX]                         */
	DT * contingence_fix_host ;
	DT * contingence_var ; /*                                      [ND][NVAR][NPOS]                   */  
	int * constCorrIndexCol ;
	DT * Y          ; /* Vecteur de Performance Y                  [NSIM][ND]        */
	DT * CD         ; /* Vecteur des CD                            [ND]        */
	int nqtl        ; /* hypothese courante */
	int np          ; /* nombre de famille de pere */
	int nsim        ; /* nombre de simulation */
	int nposGlobal  ; /* nombre de position teste */
	int npos        ; /* nombre de position teste par block */
	int ndmax       ; /* taille du tableau Y */
	int nd          ; /* nombre de descendants a prendre en compte */
	int nLevelFix   ; /* nombre d effet non dependant de la position */
	int nLevelVar   ; /* nombre d'effet dependant de la position */
	int ntniv       ; /* nLevelFix + nLevelVar */
	double seuil_cho    ; /* seuil pour l estimation des effets */
	DT * isigsq     ; /* variance residuelle des Hypothese < NQTL           [NSIM][NP] */
	int *corIpKd    ; /* *correspondance kd -> ip (indice du pere de la progeniture kd) */

	static const int MODEL_HOMO_POLYGENIC   = 0;
	static const int MODEL_HETERO_POLYGENIC = 1;
	static const int MODEL_HOMO_ANIMAL      = 2;
public:

	QTLMapStructDeviceDataLinear();

	void init(   
			int * mode,
			int * nqtlPtr,
			DT * sigsquare_host,
			DT * xinc_d_fix,
			DT * Y_d,
			DT * CD_d,
			int *corrLevelColPtr,
			double *seuil_choPtr,
			int * ndPtr,
			int * nkdPtr,
			int * nsimPtr,
			int * npositionGlobalPtr,
			//		int * npositionBlockPtr,
			int * nLevelFixPtr,
			int * nLevelVarPtr,
			int * npPtr,
			int * sizeFamilyNp);

	virtual ~QTLMapStructDeviceDataLinear();

	void set_contingence_var(int ndmax, int nLevelVar,int nbPositionsTest, int startPosition, DT* xinc_d_var,cudaStream_t * stream);

	void releaseDevice();

	size_t calculBlockPositionWorkSize(int mode) ;
} ;

class QTLMapStructDeviceWorkLinear {
public :
	DT * A_res      ; /* matrice d'incidence                                homoscedastic : [NFIX][NFIX]  , heteroscedastic:  not used ! */

	DT * XX         ; /* matrice d'incidence                                homoscedastic : [NPOS][NTNIV][NTNIV]  , heteroscedastic: [NPOS][NTNIV][NTNIV][NSIM] */
	DT * triang     ; /* matrice temporaire de la decomposition Cholesky    homoscedastic : [NPOS][NTNIV][NTNIV]  , heteroscedastic: [NPOS][NTNIV][NTNIV][NSIM] */
	DT * rhs        ; /* Vecteur RHS                                                         [NPOS][NSIM][NTNIV]  */
	cudaDeviceProp  prop ; /* property of the GPU card used */ 

public :
	QTLMapStructDeviceWorkLinear();
	/* init de XX, triang, rhs */
	void initResolution(int mode,const QTLMapStructDeviceDataLinear & data,DT *CD,DT * sigsquare_host,int * sizeFamily);

	virtual ~QTLMapStructDeviceWorkLinear();

	void releaseDeviceResolution();


	/* recuperation du tableau vecsol stoque sur le device */
	void getAllVecsolDevice(const QTLMapStructDeviceDataLinear & data,int nbBloc,int nbThread,int *vecsol) ;

} ;


/**
 * 
 */
class QTLMapStructDeviceSolutionLinear {
public:
	DT * osigsq     ; /* Vecteur des variances residuelles                  [NPOS][NSIM][NP]     */
	DT * bestim     ; /* Vecteur des solution des effets                    [NPOS][NSIM][NTNIV ] */
	DT * lrt        ; /* rapport de vraisemblance                           [NPOS][NSIM][NP] */
public :

	QTLMapStructDeviceSolutionLinear();

	void init(const QTLMapStructDeviceDataLinear & data);

	virtual ~QTLMapStructDeviceSolutionLinear();

	void releaseDevice();

} ;


class QTLMapGenericModelCalcul {

public: 
	virtual int getType() = 0 ;
virtual void start_analysis(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work) = 0 ;
virtual bool convergenceOk(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work,QTLMapStructDeviceSolutionLinear solution) = 0 ;
virtual void calcul_XT_X_A(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work) = 0;
virtual void calcul_XT_X_B(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work) = 0;
virtual void calcul_XT_X_C(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work) = 0;
virtual void calcul_Cholesky_Decomposition(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work) = 0;
virtual void calcul_RHS(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work) = 0;
virtual void calcul_LU(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work,QTLMapStructDeviceSolutionLinear solution) = 0;
virtual void calcul_SIGSQ(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work,QTLMapStructDeviceSolutionLinear solution) = 0;
};

class QTLMapHomoscedasticModelCalcul : public QTLMapGenericModelCalcul {

public:

	QTLMapHomoscedasticModelCalcul() {};
	virtual ~QTLMapHomoscedasticModelCalcul() {};

	virtual int getType() { return QTLMapStructDeviceDataLinear::MODEL_HOMO_POLYGENIC; };

	virtual void start_analysis(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work);
	virtual bool convergenceOk(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work,QTLMapStructDeviceSolutionLinear solution);
	virtual void calcul_XT_X_A(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work);
	virtual void calcul_XT_X_B(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work);
	virtual void calcul_XT_X_C(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work);
	virtual void calcul_Cholesky_Decomposition(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work);
	virtual void calcul_RHS(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work);
	virtual void calcul_LU(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work,QTLMapStructDeviceSolutionLinear solution);
	virtual void calcul_SIGSQ(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work,QTLMapStructDeviceSolutionLinear solution);
};

class QTLMapHomoscedasticAnimalModelCalcul : public QTLMapHomoscedasticModelCalcul {
private:
	DT * M            ; /* computat from animal matrix [ND][ND] => [ I - (I + lambda.A**-1)]*/
	DT * M_h          ;
	DT * FixM         ;
public:
	QTLMapHomoscedasticAnimalModelCalcul(int nd,int nfix,DT *M_host) : QTLMapHomoscedasticModelCalcul() {
		size_t size = nd*nd*sizeof(DT);
		safecall(cudaMalloc(&M,size));
		safecall(cudaMemcpy(M, M_host, size, cudaMemcpyHostToDevice));
		M_h = M_host;
		size = nd*nfix*sizeof(DT);
		safecall(cudaMalloc(&FixM,size));
	};

	virtual ~QTLMapHomoscedasticAnimalModelCalcul() { 
		safecall(cudaFree(M)); 
		safecall(cudaFree(FixM)); 
	};

	virtual int getType() { return QTLMapStructDeviceDataLinear::MODEL_HOMO_ANIMAL; };

	virtual void start_analysis(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work);
	virtual void calcul_XT_X_B(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work);
	virtual void calcul_XT_X_C(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work);
	virtual void calcul_RHS(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work);
	virtual void calcul_SIGSQ(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work,QTLMapStructDeviceSolutionLinear solution);
};

class QTLMapHeteroscedasticModelCalcul : public QTLMapGenericModelCalcul {
private :
	DT * lastLrtConvergenceHost ; /* save the last lrt to compute the convergence for each position/simulation */
	DT * lrtConvergenceDevice   ; /* do once an allocation to get lrt from the current analysis */
	DT * varInDevice            ; /* [NPOS][NSIM][NP]*/
public:
	/**
	 * Constructor
	 */
	QTLMapHeteroscedasticModelCalcul(int npos,int nsim,int np) { 
		size_t size=npos*nsim;
		safecall(cudaMallocHost(&lastLrtConvergenceHost,size*sizeof(DT)));
		safecall(cudaMalloc(&lrtConvergenceDevice,size*sizeof(DT)));

		for (size_t i=0;i<size;i++)
			lastLrtConvergenceHost[i]=99999.0;

		size=npos*nsim*np;
		safecall(cudaMalloc(&varInDevice,size*sizeof(DT)));
	} ;

	/**
	 * Destructor => delete lastlrt
	 */
	virtual ~QTLMapHeteroscedasticModelCalcul() { 
		safecall(cudaFreeHost(lastLrtConvergenceHost)); 
		safecall(cudaFree(lrtConvergenceDevice)); 
		safecall(cudaFree(varInDevice)); 
	} ;

	virtual int getType() { return QTLMapStructDeviceDataLinear::MODEL_HETERO_POLYGENIC; };

	virtual void start_analysis(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work);
	virtual bool convergenceOk(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work,QTLMapStructDeviceSolutionLinear solution);
	virtual void calcul_XT_X_A(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work);
	virtual void calcul_XT_X_B(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work);
	virtual void calcul_XT_X_C(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work);
	virtual void calcul_Cholesky_Decomposition(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work);
	virtual void calcul_RHS(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work);
	virtual void calcul_LU(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work,QTLMapStructDeviceSolutionLinear solution);
	virtual void calcul_SIGSQ(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work,QTLMapStructDeviceSolutionLinear solution);
};

size_t QTLMapStructDeviceDataLinear::calculBlockPositionWorkSize(int mode) {

	/* le calcul comprend 
	 * A_Res           : NFIX^2 ou 0 (heteroscedastic)
	 * XX              : NPOS*NTNIV*NTNIV ou NPOS*NTNIV*NTNIV*NSIM (heteroscedastic)
	 * Triang          : NPOS*NTNIV*NTNIV ou NPOS*NTNIV*NTNIV*NSIM (heteroscedastic)
	 * Rhs             : NPOS*NSIM*NTNIV
	 * OSIG            : NPOS*NSIM*NP
	 * BESTIM          : NPOS*NSIM*NTNIV
	 * LRT             : NQTL*NPOSGLOBAL*NSIM
	 * CONTINGENCE_VAR : NPOS*NVAR*NPOS
	 * 
	 * 
	 * Y,CD et contingence_fix ont deja ete alloue 
	 * 
	 */

	size_t free,total,num;
	cudaMemGetInfo(&free, &total) ;
	// avec la version driver 3.7 on etait obliger de diviser par 1.8 pour que ca passe.....;
	free = free / 1.05 ;


	if ( MODEL_HETERO_POLYGENIC == mode ) {
		num = ( nqtl*nposGlobal*nsim + np*nsim );
	}
	else if ( MODEL_HOMO_POLYGENIC == mode )  {
		num = ( nLevelFix*nLevelFix + nqtl*nposGlobal*nsim + np*nsim );
	} else if ( MODEL_HOMO_ANIMAL == mode )  {
		num = ( nLevelFix*(nLevelFix+ndmax) + nqtl*nposGlobal*nsim + np*nsim );
	}
	if ( free < num ) {
		cerr << "Not enough memory free:"<< free << " need:"<< num << endl ;
		exit(1);
	}
	size_t res = free - num*sizeof(DT);
	size_t denom ;
	/*  SIZE in memory for ONE POSITION :  HETEROSCEDASTIC              HOMOSCEDASTIC
	 *  contingence matrix :                                 ND x NVAR
	 *  XX                 :            NTNIV*NTNIVxNSIM                NTNIV * NTNIV
	 *  TRIANG             :                   ""                            ""
	 *  RHS                :                                 NSIM*NTNIV
	 *  OSIGSQ             :                                 NSIM * NP
	 *  
	 *  LRTCONVERGENCEDEVICE :                 NSIM                          0
	 *  VarInDevice          :               NSIM x NP                       0
	 */ 
	if ( MODEL_HETERO_POLYGENIC == mode ) {
		denom = (2*SIZEMATSYM(nsim,1,ntniv) + 2*nsim*ntniv+ 2*nsim*np + nd*nLevelVar + nsim);
	} else if ( MODEL_HOMO_POLYGENIC == mode )  {
		denom = (2*SIZEMATSYM(1,1,ntniv) + 2*nsim*ntniv+ nsim*np + nd*nLevelVar);
	} else if ( MODEL_HOMO_ANIMAL == mode )  {
		denom = (2*SIZEMATSYM(1,1,ntniv) + 2*nsim*ntniv+ nsim*np + nd*nLevelVar*2);
	}

	res = ceil(res/(denom*sizeof(DT)))+1;

	/* The CUDA methods Set_RHS_VAR have a limited block size in x : 16000 */
	if ( nLevelVar > 0) {
		int sizeCalculSet_Rhs_VAR = (16000 * MAX_BLOCKDIM_64) / nLevelVar ;
		if ( sizeCalculSet_Rhs_VAR < res ) res = sizeCalculSet_Rhs_VAR;
	}

	cout << "*******************************************************************************************************"<< endl;
	cout << "**   Maximum number of position to test in one block :"<< res << " / total number of position to test :"<< nposGlobal << " **"<< endl;
	cout << "*******************************************************************************************************"<< endl;
	return res ;

} ;

template <class T>
class Utils {
public:
	static void printFloatDeviceArray1D(int dim1,int nb1,T* array1D);
	static void printFloatDeviceArray2D(int dim1,int dim2,int nb1,int nb2,T* array2D);
	static void printFloatDeviceArray3D(int dim1,int dim2,int dim3,int nb1,int nb2,int nb3,T* array3D);
	static void printFloatDeviceArray4D(int dim1,int dim2,int dim3,int dim4,int nb1,int nb2,int nb3,int nb4,T* array4D);
	static void printFloatHostArray2D(int dim1,int dim2,int nb1,int nb2,T* array2D);
	static void printFloatHostArray3D(int dim1,int dim2,int dim3,int nb1,int nb2,int nb3,T* array3D);
	static void getArrayDeviceToHost(int nbBloc,int nbThread,int size,T * inDeviceArray, T *outHostArray,cudaStream_t * stream);
};

/*******************************************************************************************************************************************/

QTLMapStructDeviceDataLinear::QTLMapStructDeviceDataLinear() {
	contingence_fix = NULL ; 
	contingence_var = NULL ;
	Y = NULL;
	CD = NULL;
	isigsq = NULL;
}


void QTLMapStructDeviceDataLinear::init(
		int * mode,
		int * nqtlPtr,
		DT * sigsquare_host,
		DT * xinc_d_fix,
		DT * Y_d,
		DT * CD_d,
		int *corrLevelColPtr,
		double *seuil_choPtr,
		int * ndPtr,
		int * nkdPtr,
		int * nsimPtr,
		int * npositionGlobalPtr,
		//		int * npositionBlockPtr,
		int * nLevelFixPtr,
		int * nLevelVarPtr,
		int * npPtr,
		int * sizeFamilyNp) {

	seuil_cho   = *seuil_choPtr;
	seuil_cho   = sqrt(seuil_cho);
	ndmax       = *ndPtr       ; /* taille reelle de la matrice de contingence : lignes */
	nd          = *nkdPtr      ; /* taille de la population pris en compte */
	nsim        = *nsimPtr     ;
	npos        = 0;
	nposGlobal  = *npositionGlobalPtr;
	nLevelFix   = *nLevelFixPtr ;
	nLevelVar   = *nLevelVarPtr;
	np          = *npPtr;
	ntniv       = nLevelFix + nLevelVar ;
	nqtl        = *nqtlPtr;

#if(_CUDA_HOST_DEBUG_)
	printf("*********************** INIT DATA *********************************\n");
	printf("mode          : %d\n",*mode);
	printf("seuil_cho  SQ : %f\n",seuil_cho);
	printf("ndmax         : %d\n",ndmax);
	printf("nd            : %d\n",nd);
	printf("nsim          : %d\n",nsim);
	printf("np            : %d\n",np);
	printf("nposGlobal    : %d\n",nposGlobal);
	printf("npos          : %d\n",npos);
	printf("nLevelFix     : %d\n",nLevelFix);
	printf("nLevelVar     : %d\n",nLevelVar);
	for (int i=0;i<np;i++) printf("sizeFamilyNp[%d] : %d\n",i,sizeFamilyNp[i]);
	for (int i=0;i<ntniv;i++) printf("corrLevelCol[%d] : %d\n",i,corrLevelColPtr[i]);
	for (int i=0;i<nd;i++) printf("Y[%d]=%f,",i,(double)Y_d[i*nsim]); printf("\n");
	for (int i=0;i<nd;i++) printf("CD[%d]=%f,",i,(double)CD_d[i]);printf("\n");
	for (int i=0;i<np;i++) printf("sigsq[%d]         : %f\n",i,(double)sigsquare_host[i*nsim]);
	printf("LAST Y[%d]        : %f\n",(nd-1)*nsim,(double)Y_d[(nd-1)*nsim]);
	printf("DEPASSEMENT 1 - LAST Y[%d]        : %f\n",(nd)*nsim,(double)Y_d[(nd)*nsim]);
#endif

	/* Y */
	size_t size = nsim*nd*sizeof(DT);
	int err = cudaMalloc(&Y,size);
	if ( err ) {
		cerr << "Error allocation Y ERR="<< err << " SIZE=" << size<< endl;
		exit(1);
	}
	safecall(cudaMemcpy(Y, Y_d, size, cudaMemcpyHostToDevice));

	/* CD */
	size = nd*sizeof(DT);
	safecall(cudaMalloc(&CD,size));
	safecall(cudaMemcpy(CD, CD_d, size, cudaMemcpyHostToDevice));

	/* les effets fixes */
	size = ndmax*nLevelFix*sizeof(DT) ;
	safecall(cudaMalloc(&contingence_fix,size));
	safecall(cudaMemcpy(contingence_fix, xinc_d_fix, size, cudaMemcpyHostToDevice));

	safecall(cudaMallocHost(&contingence_fix_host,size));
	for (int i=0;i<ndmax*nLevelFix;i++) {
		contingence_fix_host[i] = xinc_d_fix[i];
	}

	if ( nLevelVar > 0 ) {
		/* isigsq */
		size = np*nsim*nqtl*sizeof(DT) ;
		safecall(cudaMalloc(&isigsq,size));
		safecall(cudaMemcpy(isigsq,allVarFitted, size, cudaMemcpyHostToDevice));
		//		for (int isim=0;isim<nsim;isim++) {
		//			for (int ip=0;ip<np;ip++)
		//				cout << allVarFitted[isim+ip*nsim]<<" " ;
		//			cout << endl ;
		//		}
		//exit(1);
	} 

	int ip=0;
	int subtotal=sizeFamilyNp[ip];
	int corIpKd_host[nd];
	for (int kd=0;kd<nd;kd++) {
		if ( kd < subtotal ) {
			corIpKd_host[kd]=ip;
		} else {
			ip++;
			subtotal+=sizeFamilyNp[ip];
			corIpKd_host[kd]=ip;
		}					
	}

	/* corIpKd */
	size = nd*sizeof(int) ;
	safecall(cudaMalloc(&corIpKd,size));
	safecall(cudaMemcpy(corIpKd, corIpKd_host, size, cudaMemcpyHostToDevice));
	safecall(cudaGetLastError());
}


void QTLMapStructDeviceDataLinear::set_contingence_var(int ndmax, int ,int nbPositionsTest, int startPosition, DT* xinc_d_var,cudaStream_t * stream) {
	if ( nLevelVar > 0 ) {
		size_t size = nd*nLevelVar*nbPositionsTest*sizeof(DT) ;
		/* 1er appel initialisation du block memoire */
		if ( contingence_var == NULL ) {
			cout << "Size of subset contingence matrix :"<< size/(1024*1024) << " Mo"<< endl ;
			/* les effets variables a chaque position */
			safecall(cudaMalloc(&contingence_var,size));

		}

		/* rearrangement memoire */
		/* On peut stocker le pointeur et eviter une reallocation memoire du block temporaire pour le transfert */
		DT* n_contvar;
		safecall(cudaMallocHost(& n_contvar, size));

		for (int i=0;i<(nLevelVar*nbPositionsTest);i++) {
			for (int j=0;j<nd;j++) {
				n_contvar[j*(nLevelVar*nbPositionsTest)+i] = xinc_d_var[i*ndmax+j];
#if(_CUDA_HOST_DEBUG_)
				if ( xinc_d_var[i*ndmax+j] != xinc_d_var[i*ndmax+j] ) { // Detection des NaN
					cerr << "npos:"<< nbPositionsTest << " nlevelvar:"<< nLevelVar  << " total:"<< nLevelVar*nbPositionsTest << " nd:" << nd << " ndmax:"<< ndmax << endl ;
					cerr << "Detected NaN in xinc_d_var at position " << i << ", animal :"<< (j+1) << ":"<< n_contvar[j*(nLevelVar*nbPositionsTest)+i]<< ","<< xinc_d_var[i*ndmax+j] << endl ;
					exit(1);
				}
#endif
			}
		}

		/* */
		safecall(cudaMemcpyAsync(contingence_var, n_contvar, size, cudaMemcpyHostToDevice,*stream));
		safecall(cudaFreeHost(n_contvar));
	}
}

QTLMapStructDeviceDataLinear::~QTLMapStructDeviceDataLinear() {
}

void QTLMapStructDeviceDataLinear::releaseDevice(){
	safecall(cudaFree(contingence_fix));
	if ( contingence_var != NULL ) safecall(cudaFree(contingence_var));
	safecall(cudaFree(Y));
	safecall(cudaFree(CD));
	safecall(cudaFreeHost(contingence_fix_host));
	//	safecall(cudaFree(sizeFam));
	if (isigsq != NULL ) safecall(cudaFree(isigsq));
}


QTLMapStructDeviceWorkLinear::QTLMapStructDeviceWorkLinear() {
	XX = NULL ;
	triang = NULL ;
	rhs = NULL ;
}

void QTLMapStructDeviceWorkLinear::initResolution(int mode,const QTLMapStructDeviceDataLinear & data,DT *CD,DT * sigsquare_host,int * sizeFamily) {
	//printf("*********************** INIT WORK *********************************\n");

	/* XX */
	size_t size = (mode == QTLMapStructDeviceDataLinear::MODEL_HETERO_POLYGENIC)?SIZEMATSYM(data.nsim,data.npos,data.ntniv):SIZEMATSYM(1,data.npos,data.ntniv) ;

	if ( size > 0 ) safecall(cudaMalloc(&XX,size*sizeof(DT)));

#if(_CUDA_HOST_DEBUG_)
	DT * XX_HOST ;
	safecall(cudaMallocHost(&XX_HOST,size*sizeof(DT)));
	for (size_t i=0;i<size;i++) 
		XX_HOST[i]=0.0;
	safecall(cudaMemcpy(XX,XX_HOST, size*sizeof(DT), cudaMemcpyHostToDevice));
	safecall(cudaFreeHost(XX_HOST));
#endif

	/* triang */
	if ( size > 0 ) safecall(cudaMalloc(&triang,size*sizeof(DT)));

#if(_CUDA_HOST_DEBUG_)
	DT * triang_HOST ;
	safecall(cudaMallocHost(&triang_HOST,size*sizeof(DT)));
	for (size_t i=0;i<size;i++) 
		triang_HOST[i]=0.0;
	safecall(cudaMemcpy(triang,triang_HOST, size*sizeof(DT), cudaMemcpyHostToDevice));
	safecall(cudaFreeHost(triang_HOST));
#endif

	/* A_res */
	size = (mode == QTLMapStructDeviceDataLinear::MODEL_HETERO_POLYGENIC)?0:data.nLevelFix*data.nLevelFix;
	safecall(cudaMalloc(&A_res,size*sizeof(DT)));

#if(_CUDA_HOST_DEBUG_)
	if ( size > 0 ) {
		DT * A_res_HOST ;
		safecall(cudaMallocHost(&A_res_HOST,size*sizeof(DT)));
		for (size_t i=0;i<size;i++) 
			A_res_HOST[i]=0.0;
		safecall(cudaMemcpy(A_res,A_res_HOST, size*sizeof(DT), cudaMemcpyHostToDevice));
		safecall(cudaFreeHost(A_res_HOST));
	}
#endif

	/* RHS */
	size = data.nsim*data.npos*data.ntniv;
	if ( size > 0 ) safecall(cudaMalloc(&rhs,size*sizeof(DT)));

#if(_CUDA_HOST_DEBUG_)
	DT * rhs_HOST ;
	safecall(cudaMallocHost(&rhs_HOST,size*sizeof(DT)));

	for (int i=0;i<size;i++) 
		rhs_HOST[i]=0.0;

	safecall(cudaMemcpy(rhs,rhs_HOST, size*sizeof(DT), cudaMemcpyHostToDevice));
	safecall(cudaFreeHost(rhs_HOST));
#endif
}

QTLMapStructDeviceWorkLinear::~QTLMapStructDeviceWorkLinear() {

}

void QTLMapStructDeviceWorkLinear::releaseDeviceResolution(){
	safecall(cudaFree(A_res));
	if (rhs != NULL) safecall(cudaFree(rhs));
	if (XX != NULL) safecall(cudaFree(XX));
	if (triang != NULL) safecall(cudaFree(triang));
}


QTLMapStructDeviceSolutionLinear::QTLMapStructDeviceSolutionLinear() {
	bestim = NULL;
	osigsq = NULL;
	lrt    = NULL;
}

void QTLMapStructDeviceSolutionLinear::init(const QTLMapStructDeviceDataLinear & data) {
	/* bestim */
	size_t size = data.nsim*data.npos*data.ntniv*sizeof(DT);
	safecall(cudaMalloc(&bestim,size));

	/* osigsq */
	size = data.nsim*data.npos*data.np*sizeof(DT);
	safecall(cudaMalloc(&osigsq,size));

	/* lrt */
	if ( data.nqtl > 0) {
		size = data.np*data.nqtl*data.nsim*data.nposGlobal*sizeof(DT);
		safecall(cudaMalloc(&lrt,size));
	} 
}

QTLMapStructDeviceSolutionLinear::~QTLMapStructDeviceSolutionLinear() {

}


void QTLMapStructDeviceSolutionLinear::releaseDevice(){
	if ( bestim != NULL ) safecall(cudaFree(bestim));
	if ( osigsq != NULL ) safecall(cudaFree(osigsq));
	if ( lrt != NULL ) safecall(cudaFree(lrt));
}

template <class T>
void Utils<T>::printFloatDeviceArray1D(int dim1,int nb1,T* array1D) {
	size_t size = dim1*sizeof(T);

	T* array1DHost ;
	int errMH = (int) cudaMallocHost(& array1DHost, size);
	if ( errMH ) {
		printf("Error host allocation of printFloatDeviceArray1D::array1DHost ERR=%d \n",errMH);
		exit(1);
	}

	cudaError_t err = cudaMemcpy(array1DHost, array1D, size, cudaMemcpyDeviceToHost);
	safecall(err);
	cout.precision(PRECISION_FLOTTANT);
	if ( err == 0 ) {
		cout << "****************************"<< endl;
		for (int i=0;i<nb1;i++ ) {
			cout << " " << array1DHost[i];
		}
		cout << endl <<"--"<< endl ;
		cout << "****************************"<< endl ;
	} else {
		cerr << "Error Memcpy :"<< err << endl ;
		exit(1);
	}
	safecall(cudaFreeHost(array1DHost));
}

template <class T>
void Utils<T>::printFloatDeviceArray2D(int dim1,int dim2,int nb1,int nb2,T* array2D) {
	size_t size = dim1*dim2*sizeof(T);

	T* array2DHost ;
	int errMH = (int) cudaMallocHost(& array2DHost, size);
	if ( errMH ) {
		printf("Error host allocation of printFloatDeviceArray2D::array2DHost ERR=%d \n",errMH);
		exit(1);
	}

	cudaError_t err = cudaMemcpy(array2DHost, array2D, size, cudaMemcpyDeviceToHost);
	safecall(err);
	cout.precision(PRECISION_FLOTTANT);
	if ( err == 0 ) {
		cout << "****************************"<< endl ;
		for (int i=0;i<nb1;i++ ) {
			for (int j=0;j<nb2;j++ ) {
				cout << " " << array2DHost[dim1*j+i];
			}
			cout << endl ;
		}
		cout << endl <<"--"<< endl ;
		cout << "****************************"<< endl ;
	} else {
		cerr << "Error Memcpy :"<< err << endl ;
		exit(1);
	}
	safecall(cudaFreeHost(array2DHost));
}

template <class T>
void Utils<T>::printFloatDeviceArray3D(int dim1,int dim2,int dim3,int nb1,int nb2,int nb3,T* array3D) {
	size_t size = dim1*dim2*dim3*sizeof(T);
	T* array3DHost ;
	int errMH = (int) cudaMallocHost(& array3DHost, size);
	if ( errMH ) {
		printf("Error host allocation of printFloatDeviceArray3D::array3DHost ERR=%d \n",errMH);
		exit(1);
	}

	cudaError_t err = cudaMemcpy(array3DHost, array3D, size, cudaMemcpyDeviceToHost);
	safecall(err);
	cout.precision(PRECISION_FLOTTANT);
	if ( err == 0 ) {
		cout << "****************************"<< endl ;
		for (int i=0;i<nb1;i++ ) {
			for (int j=0;j<nb2;j++ ) {
				for (int k=0;k<nb3 ;k++ ) {
					cout << " " << array3DHost[dim1*dim2*k+dim1*j+i];
				}
				cout << endl ;
			}
			cout << endl <<"--"<< endl ;
		}
		cout << "****************************"<< endl ;
	} else {
		cerr << "Error Memcpy :"<< err << endl ;
		exit(1);
	}
	safecall(cudaFreeHost(array3DHost));
}

template <class T>
void Utils<T>::printFloatDeviceArray4D(int dim1,int dim2,int dim3,int dim4,int nb1,int nb2,int nb3,int nb4,T* array4D) {
	size_t size = dim1*dim2*dim3*dim4*sizeof(T);
	T* array4DHost ;
	int errMH = (int) cudaMallocHost(& array4DHost, size);
	if ( errMH ) {
		printf("Error host allocation of printFloatDeviceArray4D::array4DHost ERR=%d \n",errMH);
		exit(1);
	}
	cout.precision(PRECISION_FLOTTANT);
	cudaError_t err = cudaMemcpy(array4DHost, array4D, size, cudaMemcpyDeviceToHost);
	safecall(err);
	if ( err == 0 ) {
		cout << "****************************"<< endl ;
		for (int i=0;i<nb1;i++ ) {
			for (int j=0;j<nb2;j++ ) {
				for (int k=0;k<nb3 ;k++ ) {
					for (int l=0;l<nb4 ;l++ ) {
						cout << " " << array4DHost[l*dim1*dim2*dim3+dim1*dim2*k+dim1*j+i];
					}
					cout << endl ;
				}
				cout << endl << endl;
			}
			cout << endl << " --*****************************************************************************-- "<<endl;
		}
		cout << "****************************"<< endl ;
	} else {
		cerr << "Error Memcpy :"<< err << endl ;
		exit(1);
	}
	safecall(cudaFreeHost(array4DHost));
}


template <class T>
void Utils<T>::printFloatHostArray2D(int dim1,int dim2,int nb1,int nb2,T* array2D) {

}

template <class T>
void Utils<T>::printFloatHostArray3D(int dim1,int dim2,int dim3,int nb1,int nb2,int nb3,T* array3D) {

}

template <class T>
void Utils<T>::getArrayDeviceToHost(int nbBloc,int nbThread,int size,T * inDeviceArray, T *outHostArray,cudaStream_t * stream) {
	cudaError_t err ;
	if ( stream != NULL ) {
		err = cudaMemcpyAsync(outHostArray,inDeviceArray, size*sizeof(T), cudaMemcpyDeviceToHost,*stream);
	} else {
		err = cudaMemcpy(outHostArray,inDeviceArray, size*sizeof(T), cudaMemcpyDeviceToHost);
	}
	safecall(err);
	if ( err ) {
		printf("Error QTLMapStructDeviceWorkLinear::getArrayDeviceToHost = %d \n",err);
		exit(1);
	}
}



void printMatSymXXTriang(int mode,QTLMapStructDeviceDataLinear data,DT * w) {

	int sSim=(mode == QTLMapStructDeviceDataLinear::MODEL_HETERO_POLYGENIC)?data.nsim:1;
	sSim=1;
#if(!OPTIMIZE_MEMORY_MATSYM)
	size_t size = (mode == QTLMapStructDeviceDataLinear::MODEL_HETERO_POLYGENIC)?data.ntniv*data.ntniv*data.npos*data.nsim:data.ntniv*data.ntniv*data.npos;
#else
	size_t size = (mode == QTLMapStructDeviceDataLinear::MODEL_HETERO_POLYGENIC)?(CEIL(data.ntniv*(data.ntniv+1))/2)*data.npos*data.nsim:((data.ntniv*(data.ntniv+1))/2)*data.npos;
#endif	
	DT* array ;
	safecall(cudaMallocHost(& array, size));
	safecall(cudaMemcpy(array, w, size, cudaMemcpyDeviceToHost));

	cout.precision(PRECISION_FLOTTANT);
	for (int isim=0;isim<sSim;isim++) {
		cout << endl << "************ ISIM:"<< isim << endl;
		for (int ipos=0;ipos<1;ipos++) {
			cout << " **> POS:"<< ipos << endl << "****"<< endl ;
			for (int iniv=0;iniv<data.ntniv;iniv++) {
				for (int jniv=0;jniv<data.ntniv;jniv++) {
					cout << array[GETACC(data,isim,ipos,iniv,jniv)]<<" ";
				}
				cout << endl ;
			}
		}
	}
	safecall(cudaFreeHost(array));
}


/**************************************************** MODEL METHODS ******************************************************************************************/

void QTLMapHomoscedasticModelCalcul::start_analysis(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work){
	/* calcul de XX' pour la partie des effets fixe a la position sur le host */
	size_t size = data.nLevelFix*data.nLevelFix;
	DT A_res_host[size];

	for (size_t i=0;i<data.nLevelFix;i++) {
		for (size_t j=i;j<data.nLevelFix;j++) {
			A_res_host[j*data.nLevelFix+i] = 0.0;
			for (size_t kd=0;kd<data.nd;kd++) {
				A_res_host[j*data.nLevelFix+i] += data.contingence_fix_host[i*data.ndmax+kd]*data.contingence_fix_host[j*data.ndmax+kd];
			}
		}
	}

	safecall(cudaMemcpy(work.A_res, A_res_host, size*sizeof(DT), cudaMemcpyHostToDevice));
}

void QTLMapHomoscedasticAnimalModelCalcul::start_analysis(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work){

	/* calcul de X.M.X' pour la partie des effets fixe a la position sur le host */
	size_t size = data.nLevelFix*data.nLevelFix;
	DT A_res_host[size];

	for (size_t i=0;i<data.nLevelFix;i++) {
		for (size_t j=i;j<data.nLevelFix;j++) {
			A_res_host[j*data.nLevelFix+i] = 0.0;
			for (size_t kd=0;kd<data.nd;kd++) {
				DT v2=0;
				// Calcul de Xt . M
				for (size_t kd2=0;kd2<data.nd;kd2++ ) {
					v2 += data.contingence_fix_host[i*data.ndmax+kd2]*M_h[kd2*data.nd+kd];
				}

				A_res_host[j*data.nLevelFix+i] += v2*data.contingence_fix_host[j*data.ndmax+kd];
			}
		}
	}
	safecall(cudaMemcpy(work.A_res, A_res_host, size*sizeof(DT), cudaMemcpyHostToDevice));	

	DT FixM_h[data.nLevelFix*data.nd];

	/* Calcul de l'entite Colonne Fix * M =>   FixM [NFIX][ND] */
	for (size_t i=0;i<data.nLevelFix;i++) {
		for (size_t kd=0;kd<data.nd;kd++) {
			FixM_h[kd*data.nLevelFix+i] =0.0;
			// Calcul de Xt . M
			for (size_t kd2=0;kd2<data.nd;kd2++ ) {
				FixM_h[kd*data.nLevelFix+i] += data.contingence_fix_host[i*data.ndmax+kd2]*M_h[kd2*data.nd+kd];
			}
			//cout << kd << ","<< i << ":"<< FixM_h[kd*data.nLevelFix+i] << endl ;
		}
	}
	safecall(cudaMemcpy(FixM, FixM_h, (data.nd*data.nLevelFix)*sizeof(DT), cudaMemcpyHostToDevice));	
}



__global__ void initialize_start_variance_H0(QTLMapStructDeviceDataLinear data, DT * varInDevice) {
	int isim = blockIdx.x * blockDim.x + threadIdx.x ;

	if ( isim < data.nsim ) {

		for (size_t ip=0;ip<data.np;ip++) {
			DT somyp = 0;
			int effp = 0;
			DT var   = 0;	


			for (size_t kd=0;kd<data.nd;kd++ ) {
				if ( ip == data.corIpKd[kd] )	{
					somyp = somyp + data.Y[isim+kd*data.nsim]*data.CD[kd];
					effp++;
				}
			}

			/* mean */
			DT mu = somyp / DT(effp);

			/* variance */
			for (size_t kd=0;kd<data.nd;kd++ ) {
				if ( ip == data.corIpKd[kd] ) {
					DT v = data.Y[isim+kd*data.nsim] - mu ;
					var = var + v*v;
				}
			}

			var = var / (DT(effp) - 1);
			varInDevice[isim+ip*data.nsim] = var;
		}
	}

}

/**
 * Initialize variance to start the model (classical variance estimator)
 * Calcul de XX pour les fixe a la position
 */
__global__ void initialize_start_variance(QTLMapStructDeviceDataLinear data,DT * varInDevice) {
	int ipos = blockIdx.x * blockDim.x + threadIdx.x ;
	int isim = blockIdx.y * blockDim.y + threadIdx.y ;	

	if ( ipos < data.npos && isim < data.nsim ) {
		for (int ip=0;ip<data.np;ip++) {
			int k = (data.nqtl-1)*data.nsim*data.np;
			varInDevice[ipos+isim*data.npos+ip*data.npos*data.nsim] = data.isigsq[isim+ip*data.nsim+k] ;
		}
	}
}

void QTLMapHeteroscedasticModelCalcul::start_analysis(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work){
	/* compute start variance for analysis under H0 */


	if ( data.nqtl == 0 ) {
		dim3 dimBlock(BLOCKDIMX);
		dim3 dimGrid(ceil(data.nsim / dimBlock.x )+1);
		initialize_start_variance_H0<<<dimGrid,dimBlock,0,stream>>>(data,varInDevice);

	} else {


		//		int nbBlockY=1;
		//		while ( work.prop.maxGridSize[1] < ceil(data.nsim / nbBlockY )+1 ) {
		//			nbBlockY++;
		//		}
		//
		//		dim3 dimBlock(BLOCKDIMX,nbBlockY);
		dim3 dimBlock(32,16);
		if ( work.prop.maxGridSize[0] < ceil(data.npos / dimBlock.x )+1 ) {
			cerr << "QTLMap can not support this number of position to tested..." << endl ;
			exit(1);
		}

		dim3 dimGrid(ceil(data.npos / dimBlock.x )+1,ceil(data.nsim / dimBlock.y )+1);
		initialize_start_variance<<<dimGrid,dimBlock,0,stream>>>(data,varInDevice);

	}

	cudaThreadSynchronize();

#if(_CUDA_HOST_DEBUG_)
	cudaThreadSynchronize();
	printf("PARTIEL IVAR0\n------------POS=1,SIM=1,..,nsim => NP----------------------------------------------\n");
	Utils<DT>::printFloatDeviceArray3D(data.npos,data.nsim,data.np,data.npos,1,data.np,varInDevice);
#endif
}

void QTLMapHomoscedasticModelCalcul::calcul_XT_X_A(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work) {
	/* nothing to do */
}




__global__ void initialize_A_res(QTLMapStructDeviceDataLinear data,DT * XX,DT * varInDevice) {


	int ipos = blockIdx.x * blockDim.x + threadIdx.x ;
	int isim = blockIdx.y * blockDim.y + threadIdx.y ;	

	if ( ipos < data.npos && isim < data.nsim ) {

		/* calcul de XX' des effets fixes a la position */
		for (int i=0;i<data.nLevelFix;i++) {
			for (int j=i;j<data.nLevelFix;j++) {
				DT v = 0;
				for (int kd=0;kd<data.nd;kd++) {
					v += data.contingence_fix[i*data.ndmax+kd]*data.contingence_fix[j*data.ndmax+kd]*(data.CD[kd]/varInDevice[ipos+isim*data.npos+data.corIpKd[kd]*data.npos*data.nsim]);
				}
				XX[GETACC(data,isim,ipos,data.constCorrIndexCol[j],data.constCorrIndexCol[i])] = v;
				//work.XX[GETACC(data,isim,ipos,data.constCorrIndexCol[i],data.constCorrIndexCol[j])] = v;
			}
		}
	}
}

void QTLMapHeteroscedasticModelCalcul::calcul_XT_X_A(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work) {

	/* sur les M2050 la limitation est de 65535, factuer limitant pour les simulations 
	 * si l utilisateur demande plus, on cree des blocks pour les simuls...*/
	//	int nbBlockY=1;
	//	while ( work.prop.maxGridSize[1] < ceil(data.nsim / nbBlockY )+1 ) {
	//		nbBlockY++;
	//	}
	//dim3 dimBlock(BLOCKDIMX,nbBlockY);
	dim3 dimBlock(32,16);

	if ( work.prop.maxGridSize[0] < ceil(data.npos / dimBlock.x )+1 ) {
		cerr << "QTLMap can not support this number of position to tested..." << endl ;
		exit(1);
	}

	dim3 dimGrid(ceil(data.npos / dimBlock.x )+1,ceil(data.nsim / dimBlock.y )+1);
	initialize_A_res<<<dimGrid,dimBlock,0,stream>>>(data,work.XX,varInDevice);
}

bool QTLMapHomoscedasticModelCalcul::convergenceOk(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work,QTLMapStructDeviceSolutionLinear solution) {
	return true;
}


__global__ void calcul_lrt_and_initialize_start_variance(QTLMapStructDeviceDataLinear data,QTLMapStructDeviceSolutionLinear solution,DT * varInDevice,DT * lrtConv) {
	int ipos = blockIdx.x * blockDim.x + threadIdx.x ;
	int isim = blockIdx.y * blockDim.y + threadIdx.y ;

	if ( ipos < data.npos && isim < data.nsim ) {
		DT v = 0.0; 
		for (int ip=0;ip<data.np;ip++) {
			v += constSizeFam[ip] * (log(solution.osigsq[ipos+isim*data.npos+ip*data.npos*data.nsim]) - log(varInDevice[ipos+isim*data.npos+ip*data.npos*data.nsim]));
			varInDevice[ipos+isim*data.npos+ip*data.npos*data.nsim] = solution.osigsq[ipos+isim*data.npos+ip*data.npos*data.nsim];
		}
		lrtConv[ipos+isim*data.npos] = v;
	}
}

bool QTLMapHeteroscedasticModelCalcul::convergenceOk(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work,QTLMapStructDeviceSolutionLinear solution) {

	//	int nbBlockY=1;
	//	while ( work.prop.maxGridSize[1] < ceil(data.nsim / nbBlockY )+1 ) {
	//		nbBlockY++;
	//	}

	//dim3 dimBlock(BLOCKDIMX,nbBlockY);
	dim3 dimBlock(32,16);

	if ( work.prop.maxGridSize[0] < ceil(data.npos / dimBlock.x )+1 ) {
		cerr << "QTLMap can not support this number of position to tested..." << endl ;
		exit(1);
	}

	dim3 dimGrid(ceil(data.npos / dimBlock.x )+1,ceil(data.nsim / dimBlock.y )+1);
	calcul_lrt_and_initialize_start_variance<<<dimGrid,dimBlock,0,stream>>>(data,solution,varInDevice,lrtConvergenceDevice);
	cudaThreadSynchronize();

	//
	size_t size=data.npos*data.nsim;
	DT  *lrt;
	//safecall(cudaMallocHost(&lrt,size*sizeof(DT)));
	lrt = new DT[size];

	safecall(cudaMemcpyAsync(lrt,lrtConvergenceDevice, size*sizeof(DT), cudaMemcpyDeviceToHost,stream));
	int countOk=0;

	for (int ipos=0;ipos<data.npos;ipos++) {
		for (int isim=0;isim<data.nsim;isim++) {
			int i = ipos + isim*data.npos ;
			if ( fabs(lastLrtConvergenceHost[i] - lrt[i]) < 0.5 )  countOk++;
			//			else
			//				cout << "ipos:"<< ipos << " isim:"<< isim << " "<< lrt[i] << " "<< lastLrtConvergenceHost[i] << " diff:" << fabs(lastLrtConvergenceHost[i] - lrt[i])<< endl ;
			lastLrtConvergenceHost[i] = lrt[i] ;
		}
	}

	delete [] lrt;
	//safecall(cudaFreeHost(lrt));

	if ( countOk ==  data.nsim*data.npos ) {
		cout << "ok" << endl ;
		return true;
	}

	cout << "convergence :" << countOk << "/"<< data.nsim*data.npos << endl ;

	return false;
}

/**********************************************************************************************************************************************
 * NFIX  : nombre de niveau fixe a la position (moyenne general, effet polygenique, effets de nuisances)
 * NVAR  : nombre de niveaux varaible a la position (effet qtl ou haplotype)
 * NTNIV : NFIX + NVAR 
 * NPOS  : nombre de position teste 
 * X     : matrice de contingence (animal par ligne / niveau exprime en colonne)
 * XT    : transpose de X
 * XT.X => Calcul independant de block A ( NFIX * NFIX operation ), B ( NFIX * NVAR * NPOS) , C (NVAR * NVAR/2 * NPOS)
 * 
 * Exemple : NFIX = 4, NVAR = 3   
 * 
 *              LEV FIX     POS1        POS2        POS3          POS4
 *           ( 1 1 0 0  0.5 0  0    0.4  0  0    0.3  0  0    0.2  0  0  )      NPOS = 4
 * TOTAL_X = ( 1 1 0 0  0.5 0  0    0.4  0  0    0.3  0  0    0.2  0  0  )
 *           ( 1 0 1 0  0  0.5 0     0  0.4  0    0  0.4  0    0  0.4  0 )    
 *           ( 1 0 1 0  0  0.5 0     0  0.4  0    0  0.4  0    0  0.4  0 )
 *           ( 1 0 0 1  0   0 0.5    0  0  0.4    0  0   0.4   0  0  0.4 )
 *           ( 1 0 0 1  0   0 0.5    0  0  0.4    0  0   0.4   0  0  0.4 )
 * 
 * pour POS = 1
 * 
 *     ( 1 1 0 0 0.5 0  0  )
 * X = ( 1 1 0 0 0.5 0  0  )
 *     ( 1 0 1 0  0 0.5 0  )     ND=6
 *     ( 1 0 1 0  0 0.5 0  )
 *     ( 1 0 0 1  0  0 0.5 )
 *     ( 1 0 0 1  0  0 0.5 )
 * 
 * 
 * XT.X = ( 1.1 1.2 1.3 1.4 1.5 1.6 1.7 )
 *        (     2.2 2.3 2.4 2.5 2.6 2.7 )
 *        (         3.3 3.4 3.5 3.6 3.7 )            ( 1.1 1.2 1.3 1.4 )         ( 1.5 1.6 1.7 )        ( 5.5 5.6 5.7 ) 
 *        (             4.4 4.5 4.6 4.7 )   avec A = (     2.2 2.3 2.4 )    B  = ( 2.5 2.6 2.7 )    C = (     6.6 6.7 )
 *        (                 5.5 5.6 5.7 )            (         3.3 3.4 )         ( 3.5 3.6 3.7 )        (         7.7 )
 *        (                     6.6 6.7 )            (             4.4 )         ( 4.5 4.6 4.7 )        
 *        (                         7.7 )
 * 
 * 
 *  - Calcul de A sur le HOST
 *  - Calcul de B sur le DEVICE
 *  - Calcul de C sur le DEVICE
 *  
 *  
 *  Calcul B :
 *  BLOCK à 2 dimensions : ivar => numero de la colonne dans la liste des niveaux variables
 *                         ifix => numero de la colonne dans la liste des niveaux fixes
 *                         
 *  
 *  Ordonnancement des resultats :
 *  [ifix] => |  1.5 pos1  | 1.6 pos1 |  1.7 pos1 |  1.5 pos2 |   | 1.6 pos2 |  1.7 pos2 ....
 *  
 */


__global__ void XT_X_B_homoscedastic(QTLMapStructDeviceDataLinear data,DT * XX) {

	extern __shared__ DT FIXVALUE[];

	int ivar = blockIdx.x * blockDim.x + threadIdx.x ;

	if ( blockIdx.y < constNbLevelFix ) {
		int i_shared = threadIdx.x;
		while(i_shared < data.nd)
		{
			FIXVALUE[i_shared] = data.contingence_fix[blockIdx.y*data.ndmax+i_shared];
			i_shared += blockDim.x;
		}
		__syncthreads();
	}


	if ( ivar < constNbLevelVar*data.npos ) {
		DT v = 0;
		for (int kd=0;kd<data.nd;kd++ ) {
			v += FIXVALUE[kd]*data.contingence_var[kd*constNbLevelVar*data.npos+ivar];
		}

		/* niveau à une position */
		int iniv = ivar % constNbLevelVar ;
		/* position teste */
		int ipos = ivar / constNbLevelVar ;

		XX[GETACC(data,0,ipos,data.constCorrIndexCol[constNbLevelFix + iniv],data.constCorrIndexCol[blockIdx.y])] = v;
	}
}
void QTLMapHomoscedasticModelCalcul::calcul_XT_X_B(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work) {

	dim3 dimBlock(BLOCKDIMX);
	dim3 dimGrid (ceil(data.npos*data.nLevelVar / dimBlock.x) + 1 , (data.nLevelFix / dimBlock.y) + 1);
	XT_X_B_homoscedastic<<<dimGrid,dimBlock,data.nd*sizeof(DT),stream>>>(data,work.XX);
	safecall(cudaGetLastError());
}

//******************************************************
__global__ void XT_X_B_homoscedastic_animal(QTLMapStructDeviceDataLinear data,DT * FixM,DT* M ,DT * XX) {

	//ivar : indice du niveau pour toute les positions (0<ivar<ntniv_var*npos)
	int ivar = blockIdx.x * blockDim.x + threadIdx.x ;
	// inivfix : iniv fix dans contingence fix
	int inivfix = blockIdx.y;
	if ( ivar < constNbLevelVar*data.npos ) {
		/* indice du niveau dans une matrice de contingence independant de la position (0<iniv<ntniv) */
		int iniv = ivar % constNbLevelVar ;
		/* position teste */
		int ipos = ivar / constNbLevelVar ;

		DT v = 0;
		for (int kd=0;kd<data.nd;kd++ ) {
			v += FixM[kd*constNbLevelFix+inivfix]*data.contingence_var[kd*constNbLevelVar*data.npos+ivar];
		}
		XX[GETACC(data,0,ipos,data.constCorrIndexCol[constNbLevelFix + iniv],data.constCorrIndexCol[inivfix])] = v;
	}
}

void QTLMapHomoscedasticAnimalModelCalcul::calcul_XT_X_B(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work) {

	dim3 dimBlock(BLOCKDIMX);
	dim3 dimGrid (ceil(data.npos*data.nLevelVar / dimBlock.x) + 1 , (data.nLevelFix / dimBlock.y) + 1);


	XT_X_B_homoscedastic_animal<<<dimGrid,dimBlock,0,stream>>>(data,FixM,M,work.XX);
	safecall(cudaGetLastError());
}



__global__ void XT_X_B_heteroscedastic(DT * varInDevice,QTLMapStructDeviceDataLinear data,DT * XX) {

	extern __shared__ DT FIXVALUE[];

	int ivar = blockIdx.x * blockDim.x + threadIdx.x ;
	int isim = blockIdx.y * blockDim.y + threadIdx.y ;
	int inivfix = blockIdx.z ;

	if ( inivfix < constNbLevelFix ) {
		int i_shared = threadIdx.x;
		while(i_shared < data.nd)
		{
			FIXVALUE[i_shared] = data.contingence_fix[inivfix*data.ndmax+i_shared];
			i_shared += blockDim.x;
		}
		__syncthreads();
	}


	if ( ivar < constNbLevelVar*data.npos && isim < data.nsim) {
		/* niveau à une position */
		int iniv = ivar % constNbLevelVar ;
		/* position teste */
		int ipos = ivar / constNbLevelVar ;

		DT v = 0.0;
		for (int kd=0;kd<data.nd;kd++ ) {
			v += FIXVALUE[kd]*data.contingence_var[kd*constNbLevelVar*data.npos+ivar]*(data.CD[kd] / varInDevice[ipos+isim*data.npos+data.corIpKd[kd]*data.nsim*data.npos]) ; 
		}

		XX[GETACC(data,isim,ipos,data.constCorrIndexCol[constNbLevelFix + iniv],data.constCorrIndexCol[inivfix])] = v;
		//work.XX[GETACC(data,isim,ipos,data.constCorrIndexCol[inivfix],data.constCorrIndexCol[constNbLevelFix + iniv])] = v;
	}
}

void QTLMapHeteroscedasticModelCalcul::calcul_XT_X_B(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work) {

	int nbBlockY=1;
	while ( work.prop.maxGridSize[1] < ceil(data.nsim / nbBlockY )+1 ) {
		nbBlockY++;
	}

	dim3 dimBlock(BLOCKDIMX,nbBlockY);

	if ( work.prop.maxGridSize[0] < ceil(data.npos*data.nLevelVar / dimBlock.x) + 1 ) {
		cerr << "QTLMap can not support this number of position ["<< data.npos << "] to tested with the number of lever ["<< data.nLevelVar <<"]" << endl ;
		exit(1);
	}

	dim3 dimGrid (ceil(data.npos*data.nLevelVar / dimBlock.x) + 1 , ceil(data.nsim / dimBlock.y )+1, data.nLevelFix);
	XT_X_B_heteroscedastic<<<dimGrid,dimBlock,data.nd*sizeof(DT),stream>>>(varInDevice,data,work.XX);
}

/* calcul d'une valeur de la diagonale et des valeurs a droite de cette diagonale de XT.X   = colonne et ligne de niveau variable =  */
/*
 * 
 *     ( 5.5 5.6 5.7 )
 * C = (     6.6 6.7 )
 *     (         7.7 )
 *     
 *     
 *     stockage du resultat dans C_res sous la forme
 *     
 *            POS 1              POS2                   POS N
 *     | 5.5 | 6.6 | 7.7 | 5.5 | 6.6 | 7.7 | ... | 5.5 | 6.6 | 7.7 |
 *     puis
 *            POS 1     POS 2            POS N
 *     | 5.6 | 6.7 | X | 5.6 | 6.7 | X | ... | 5.6 | 6.7 | X |
 *     puis
 *       POS1              POS2               POS N
 *     [ 5.7 | X | X | 5.7 | X | X |.... | 5.7 | X | X |  
 * 
 */
__global__ void XT_X_C_homoscedastic(QTLMapStructDeviceDataLinear data,DT * XX) {
	int ivar = blockIdx.x * blockDim.x + threadIdx.x ;

	/* ivar : la diagonale de XT . X */

	if ( ivar < constNbLevelVar*data.npos ) {	           
		/* iniv position entre 0 et ntlevVar */ 
		int nfix = constNbLevelFix ;
		int nvar = constNbLevelVar ;

		/* niveau à une position */
		int iniv = ivar % nvar ;
		/* position teste */
		int ipos = ivar / nvar ;

		/* initialisation */
		//DT viniv[MAXND_WORK] ;
		DT v=0.0,v2=0.0;

		for (int kd=0;kd<data.nd;kd++ ) {
			//viniv[kd] = data.contingence_var[kd*nvar*data.npos+ivar];
			//v += viniv[kd]*viniv[kd];	
			v2=data.contingence_var[kd*nvar*data.npos+ivar];
			v += v2*v2;
		}		

		/* calcul de la diagonal */

		/* les diagonales sont stockees en premieres...*/ 

		XX[GETACC(data,0,ipos,data.constCorrIndexCol[nfix + iniv],data.constCorrIndexCol[nfix + iniv])] = v;

		for (int jniv=1;jniv<(nvar-iniv);jniv++) {
			DT v = 0.0;	
			for (int kd=0;kd<data.nd;kd++ ) {
				//v +=  viniv[kd]*data.contingence_var[kd*nvar*data.npos+(ivar+jniv)];
				v += data.contingence_var[kd*nvar*data.npos+ivar]*data.contingence_var[kd*nvar*data.npos+(ivar+jniv)];
			}

			XX[GETACC(data,0,ipos,data.constCorrIndexCol[nfix + jniv + iniv],data.constCorrIndexCol[nfix + iniv])] = v;
			//work.XX[GETACC(data,0,ipos,data.constCorrIndexCol[nfix + iniv],data.constCorrIndexCol[nfix + jniv + iniv])] = v;
		}
	}
}


void QTLMapHomoscedasticModelCalcul::calcul_XT_X_C(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work) {
	dim3 dimBlock(BLOCKDIMX);
	dim3 dimGrid2 (ceil(data.npos*data.nLevelVar / dimBlock.x) + 1 , 1);
	XT_X_C_homoscedastic<<<dimGrid2,dimBlock,0,stream>>>(data,work.XX);
}


__global__ void XT_X_C_homoscedastic_animal_init(QTLMapStructDeviceDataLinear data,DT * XX) {
	int ivar = blockIdx.x * blockDim.x + threadIdx.x ;
	/* ivar : la diagonale de XT . X */
	if ( ivar < constNbLevelVar*data.npos ) {	           
		/* iniv position entre 0 et ntlevVar */ 
		int nfix = constNbLevelFix ;

		/* niveau à une position */
		int iniv = ivar % constNbLevelVar ;
		/* position teste */
		int ipos = ivar / constNbLevelVar ;

		XX[GETACC(data,0,ipos,data.constCorrIndexCol[nfix + iniv],data.constCorrIndexCol[nfix + iniv])] = 0.0;

		for (int jniv=1;jniv<(constNbLevelVar-iniv);jniv++) {
			XX[GETACC(data,0,ipos,data.constCorrIndexCol[nfix + jniv + iniv],data.constCorrIndexCol[nfix + iniv])] = 0.0;
		}
	}


}

__global__ void XT_X_C_homoscedastic_animal(QTLMapStructDeviceDataLinear data,int ikd,DT * M,DT * XX) {
	int ivar = blockIdx.x * blockDim.x + threadIdx.x ;

	/* ivar : la diagonale de XT . X */
	if ( ivar < constNbLevelVar*data.npos && ikd < data.nd ) {	           
		/* iniv position entre 0 et ntlevVar */ 
		int nfix = constNbLevelFix ;
		int nvar = constNbLevelVar ;

		/* niveau à une position */
		int iniv = ivar % nvar ;
		/* position teste */
		int ipos = ivar / nvar ;

		DT v=0;/* Calcul de XT.M */
		for (int kd=0;kd<data.nd;kd++) {
			v += data.contingence_var[kd*nvar*data.npos+ivar]*M[ikd*data.nd+kd];
		}

		for (int jniv=0;jniv<(nvar-iniv);jniv++) {
			/* fin entite constante..... */
			DT v2 =  v*data.contingence_var[ikd*nvar*data.npos+(ivar+jniv)];
			XX[GETACC(data,0,ipos,data.constCorrIndexCol[nfix + jniv + iniv],data.constCorrIndexCol[nfix + iniv])] += v2;
		}
	}
}


void QTLMapHomoscedasticAnimalModelCalcul::calcul_XT_X_C(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work) {

	dim3 dimBlock(BLOCKDIMX);
	dim3 dimGrid2 (ceil(data.npos*data.nLevelVar / dimBlock.x) + 1 , 1 );

	// initialize memory to 0
	XT_X_C_homoscedastic_animal_init<<<dimGrid2,dimBlock,0,stream>>>(data,work.XX);
	cudaStreamSynchronize(stream);

	/* for each ikd we compute the following sub value  : IVAR[ .... ]    * M[ IKD, ....]   * CONT_VAR[....,IVAR:NVAR] */
	for (int ikd=0;ikd < data.nd; ikd++) {
		XT_X_C_homoscedastic_animal<<<dimGrid2,dimBlock,0,stream>>>(data,ikd,M,work.XX);
		cudaStreamSynchronize(stream);
	}
}


__global__ void XT_X_C_heteroscedastic(DT * varInDevice,QTLMapStructDeviceDataLinear data,DT * XX) {
	int ivar = blockIdx.x * blockDim.x + threadIdx.x ;
	int isim = blockIdx.y * blockDim.y + threadIdx.y ;

	/* ivar : la diagonale de XT . X */

	if ( ivar < constNbLevelVar*data.npos && isim  < data.nsim ) {	           
		/* iniv position entre 0 et ntlevVar */ 
		size_t nfix = constNbLevelFix ;
		size_t nvar = constNbLevelVar ;

		/* niveau à une position */
		size_t iniv = ivar % nvar ;
		/* position teste */
		size_t ipos = ivar / nvar ;

		/* initialisation */
		/* On ne passse plus par un tableau local==> trop gourmand en moire du coup les grosse analyses ne passait pas (2QTL sur puce 54K)*
		   On accede directement a la globalMemory
		 */
		/* DT viniv[MAXND_WORK] ; */
		DT v=0.0,v2=0.0;

		//int scalesimXX = isim*BASE_SCALE_SYMMAT;
		v=0.0;
		for (size_t kd=0;kd<data.nd;kd++ ) {
			/* viniv[kd] = data.contingence_var[kd*nvar*data.npos+ivar];*/
			//v += viniv[kd]*viniv[kd]*(data.CD[kd] / varInDevice[ipos+isim*data.npos+data.corIpKd[kd]*data.nsim*data.npos]);
			v2 = data.contingence_var[kd*nvar*data.npos+ivar];
			v += v2*v2*(data.CD[kd] / varInDevice[ipos+isim*data.npos+data.corIpKd[kd]*data.nsim*data.npos]);	
		}		

		/* calcul de la diagonal */

		/* les diagonales sont stockees en premieres...*/ 
		XX[GETACC(data,isim,ipos,data.constCorrIndexCol[nfix + iniv],data.constCorrIndexCol[nfix + iniv])] = v;

		for (size_t jniv=1;jniv<(nvar-iniv);jniv++) {
			//scale +=  nvar*data.npos;
			DT v = 0.0;	

			for (size_t kd=0;kd<data.nd;kd++ ) {
				//v +=  viniv[kd]*data.contingence_var[kd*nvar*data.npos+(ivar+jniv)]*(data.CD[kd] / varInDevice[ipos+isim*data.npos+data.corIpKd[kd]*data.nsim*data.npos]);
				v2 = data.contingence_var[kd*nvar*data.npos+ivar];
				v +=  v2*data.contingence_var[kd*nvar*data.npos+(ivar+jniv)]*(data.CD[kd] / varInDevice[ipos+isim*data.npos+data.corIpKd[kd]*data.nsim*data.npos]);
			}
			XX[GETACC(data,isim,ipos,data.constCorrIndexCol[nfix + jniv + iniv],data.constCorrIndexCol[nfix + iniv])] = v;
		}
	}
}

__global__ void XT_X_C_heteroscedastic_opt_init(QTLMapStructDeviceDataLinear data,DT * XX) {
	int ivar = blockIdx.x * blockDim.x + threadIdx.x ;
	int isim = blockIdx.y * blockDim.y + threadIdx.y ;

	if ( ivar < constNbLevelVar*data.npos && isim  < data.nsim ) {	           
		/* iniv position entre 0 et ntlevVar */ 
		int nfix = constNbLevelFix ;
		int nvar = constNbLevelVar ;

		/* niveau à une position */
		int iniv = ivar % nvar ;
		/* position teste */
		int ipos = ivar / nvar ;

		for (int jniv=0;jniv<(nvar-iniv);jniv++) {	
			XX[GETACC(data,isim,ipos,data.constCorrIndexCol[nfix + jniv + iniv],data.constCorrIndexCol[nfix + iniv])] = 0;
		}
	}
}

__global__ void XT_X_C_heteroscedastic_opt(DT * varInDevice,int kd,QTLMapStructDeviceDataLinear data,DT * XX) {
	int ivar = blockIdx.x * blockDim.x + threadIdx.x ;
	int isim = blockIdx.y * blockDim.y + threadIdx.y ;

	/* ivar : la diagonale de XT . X */

	if ( ivar < constNbLevelVar*data.npos && isim  < data.nsim ) {	           
		/* iniv position entre 0 et ntlevVar */ 
		int nfix = constNbLevelFix ;
		int nvar = constNbLevelVar ;

		/* niveau à une position */
		int iniv = ivar % nvar ;
		/* position teste */
		int ipos = ivar / nvar ;

		DT v = data.contingence_var[kd*nvar*data.npos+ivar] ;
		for (int jniv=0;jniv<(nvar-iniv);jniv++) {	
			DT v2 =  v*data.contingence_var[kd*nvar*data.npos+(ivar+jniv)] *(data.CD[kd] / varInDevice[ipos+isim*data.npos+data.corIpKd[kd]*data.nsim*data.npos]);
			XX[GETACC(data,isim,ipos,data.constCorrIndexCol[nfix + jniv + iniv],data.constCorrIndexCol[nfix + iniv])] += v2;
		}
	}
}




void QTLMapHeteroscedasticModelCalcul::calcul_XT_X_C(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work) {
	int nbBlockY=1;
	while ( work.prop.maxGridSize[1] < ceil(data.nsim / nbBlockY )+1 ) {
		nbBlockY++;
	}

	dim3 dimBlock(BLOCKDIMX,nbBlockY);

	if ( work.prop.maxGridSize[0] < ceil(data.npos*data.nLevelVar / dimBlock.x) + 1 ) {
		cerr << "QTLMap can not support this number of position ["<< data.npos << "] to tested with the number of lever ["<< data.nLevelVar <<"]" << endl ;
		exit(1);
	}

	dim3 dimGrid (ceil(data.npos*data.nLevelVar / dimBlock.x) + 1 , ceil(data.nsim / nbBlockY )+1);
	XT_X_C_heteroscedastic<<<dimGrid,dimBlock,0,stream>>>(varInDevice,data,work.XX);
	//	// initialize memory to 0
	//	XT_X_C_heteroscedastic_opt_init<<<dimGrid,dimBlock,0,stream>>>(data,work.XX);
	//	cudaStreamSynchronize(stream);
	//
	//	/* for each ikd we compute the following sub value  : IVAR[ .... ]    * M[ IKD, ....]   * CONT_VAR[....,IVAR:NVAR] */
	//	for (int ikd=0;ikd < data.nd; ikd++) {
	//		XT_X_C_heteroscedastic_opt<<<dimGrid,dimBlock,0,stream>>>(varInDevice,ikd,data,work.XX);
	//		cudaStreamSynchronize(stream);
	//	}

}


/**
 *
 * Cholesky_Decomposition
 *
 *
 *
 */
__global__ void Cholesky_Decomposition_homoscedastic(QTLMapStructDeviceDataLinear data,DT * A_res,DT * XX, DT * triang) {

	int ipos = blockIdx.x * blockDim.x + threadIdx.x ; /* position */

	//DT * triang = work.triang ;
	int ntniv = data.ntniv;
	int nfix = constNbLevelFix;

	if ( ipos < data.npos ) {

		for (int i=0;i<ntniv;i++) {
			for (int j=i;j<ntniv;j++) {
				triang[GETACC(data,0,ipos,i,j)] = 0.0;
			}
		}
		/* a optimiser , ces valeurs sont redondantes dans les matrices XX pour chaque position  */
		for (int i=0;i<nfix;i++) {
			for (int j=i;j<nfix;j++) {
				XX[GETACC(data,0,ipos,data.constCorrIndexCol[j],data.constCorrIndexCol[i])] = A_res[j*nfix+i];
			}
		}

		for (int j=0;j<ntniv;j++) {
			DT v = XX[GETACC(data,0,ipos,j,j)]; //XX[j][j] ;

			for (int k=0;k<=j-1;k++) {
				// [j][j] = [j][j] - [k][j]*[k][j]
				//v -= triang[ipos+k*data.npos+j*ntniv*data.npos]*triang[ipos+k*data.npos+j*ntniv*data.npos];
				v -= triang[GETACC(data,0,ipos,k,j)]*triang[GETACC(data,0,ipos,k,j)];
			}
			//si estimable .. on ne met pas de structure vecseuil==> on utilise la diagonale de triangle pour savoir l effet est estimable
			//triang[ipos+j*data.npos +j*ntniv*data.npos] =  sqrt(v);
			triang[GETACC(data,0,ipos,j,j)] =  sqrt(v);
			//if ( triang[ipos+j*data.npos +j*ntniv*data.npos] != triang[ipos+j*data.npos +j*ntniv*data.npos] ) triang[ipos+j*data.npos +j*ntniv*data.npos] = 0.0;
			if ( triang[GETACC(data,0,ipos,j,j)] != triang[GETACC(data,0,ipos,j,j)] ) triang[GETACC(data,0,ipos,j,j)] = 0.0;
			if ( triang[GETACC(data,0,ipos,j,j)] > data.seuil_cho ) { 
				for (int i=j+1;i<ntniv;i++ ) {
					//triang[ipos+i*data.npos +j*ntniv*data.npos] = work.XX[ipos+i*data.npos +j*ntniv*data.npos];//XX[i][j];
					triang[GETACC(data,0,ipos,i,j)] = XX[GETACC(data,0,ipos,i,j)];//XX[i][j];
					for ( int k=0;k<=j-1;k++) {
						triang[GETACC(data,0,ipos,i,j)] -= (triang[GETACC(data,0,ipos,i,k)] * triang[GETACC(data,0,ipos,j,k)]);
					}
					triang[GETACC(data,0,ipos,i,j)] = triang[GETACC(data,0,ipos,i,j)] / triang[GETACC(data,0,ipos,j,j)];
					//triang[ipos+j*data.npos +i*ntniv*data.npos] = triang[ipos+i*data.npos +j*ntniv*data.npos] ;
				}
			}  
		}
	}

}


void QTLMapHomoscedasticModelCalcul::calcul_Cholesky_Decomposition(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work) {

	int nbblock_512 = ceil(data.npos / MAX_BLOCKDIM_512) + 1 ;
	Cholesky_Decomposition_homoscedastic<<<nbblock_512,MAX_BLOCKDIM_512,0,stream>>>(data,work.A_res,work.XX,work.triang);
}


__global__ void Cholesky_Decomposition_heterocedastic(QTLMapStructDeviceDataLinear data,DT * XX, DT * triang) {

	int ipos = blockIdx.x * blockDim.x + threadIdx.x ; /* position */
	int isim = blockIdx.y * blockDim.y + threadIdx.y ;

	//DT * triang = work.triang ;
	int ntniv = data.ntniv;

	if ( ipos < data.npos && isim < data.nsim ) {
		for (int i=0;i<ntniv;i++) {
			for (int j=i;j<ntniv;j++) {
				triang[GETACC(data,isim,ipos,i,j)] = 0.0;
			}
		}


		for (int j=0;j<ntniv;j++) {
			DT v = XX[GETACC(data,isim,ipos,j,j)]; //XX[j][j] ;

			for (int k=0;k<=j-1;k++) {
				// [j][j] = [j][j] - [k][j]*[k][j]
				v -= triang[GETACC(data,isim,ipos,j,k)]*triang[GETACC(data,isim,ipos,j,k)];
			}

			//si estimable .. on ne met pas de structure vecseuil==> on utilise la diagonale de triangle pour savoir l effet est estimable
			triang[GETACC(data,isim,ipos,j,j)] =  sqrt(v);
			if ( triang[GETACC(data,isim,ipos,j,j)] != triang[GETACC(data,isim,ipos,j,j)] ) triang[GETACC(data,isim,ipos,j,j)] = 0.0;

			if ( triang[GETACC(data,isim,ipos,j,j)] > data.seuil_cho ) { 
				for (int i=j+1;i<ntniv;i++ ) {
					triang[GETACC(data,isim,ipos,i,j)] = XX[GETACC(data,isim,ipos,i,j)];//XX[i][j];
					for ( int k=0;k<=j-1;k++) {
						triang[GETACC(data,isim,ipos,i,j)] -= (triang[GETACC(data,isim,ipos,i,k)] * triang[GETACC(data,isim,ipos,j,k)]);
					}
					triang[GETACC(data,isim,ipos,i,j)] = triang[GETACC(data,isim,ipos,i,j)] / triang[GETACC(data,isim,ipos,j,j)];
				}
			}
		}
	}

}


void QTLMapHeteroscedasticModelCalcul::calcul_Cholesky_Decomposition(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work) {
	int nbBlockY=1;
	while ( work.prop.maxGridSize[1] < ceil(data.nsim / nbBlockY )+1 ) {
		nbBlockY++;
	}

	dim3 dimBlock(BLOCKDIMX,nbBlockY);

	if ( work.prop.maxGridSize[0] < ceil(data.npos / dimBlock.x) + 1 ) {
		cerr << "QTLMap can not support this number of position ["<< data.npos << "] to tested."<< endl ;
		exit(1);
	}

	dim3 dimGrid (ceil(data.npos / dimBlock.x) + 1 , ceil(data.nsim / nbBlockY )+1);
	Cholesky_Decomposition_heterocedastic<<<dimGrid,dimBlock,0,stream>>>(data,work.XX,work.triang);
}

/*
 * A optimiser : rhs pour les effet fixe a la position ne devrait prendre que nsim*nlevelfix et non pas nsim*nlevelfix*npos
 * 
 * 
 * logique de corrIndexCol :
 * 
 *     I=1,..,NTNIV | I=1,..,NFIX,NFIX+1,..,NVAR
 *     -------------------------------------------
 *     0            |            0
 *     1            |            J
 *     2            |            J+1
 *     3            |            .
 *     .            |            1
 *     .            |            2
 *     .            |            3
 *     NTNIV        |            NTNIV
 * 
 * 
 * 
 */
__global__ void Set_RHS_FIX_homoscedatic(QTLMapStructDeviceDataLinear data,DT * rhs,DT * triang) {
	int ipos = blockIdx.x * blockDim.x + threadIdx.x ; /* position */
	int isim = blockIdx.y * blockDim.y + threadIdx.y ;// * blockDim.y + threadIdx.y ;

	int idNiv;
	if ( ipos < data.npos && isim < data.nsim ) {
		for (int iniv=0;iniv<constNbLevelFix;iniv++) {

			idNiv = data.constCorrIndexCol[iniv];

			int index = ipos + isim*data.npos + idNiv*data.nsim*data.npos ;
			DT v=0;
			if (triang[GETACC(data,0,ipos,idNiv,idNiv)]>data.seuil_cho) {	
				for (int kd = 0;kd<data.nd;kd++ ) {
					v += data.contingence_fix[iniv*data.ndmax + kd]*data.Y[isim+kd*data.nsim];
				}	
				rhs[index] = v;
			}
		}
	}
}

__global__ void Set_RHS_VAR_homoscedatic(QTLMapStructDeviceDataLinear data,DT * rhs,DT * triang) {
	int ivar = blockIdx.x * blockDim.x + threadIdx.x ; /* position */
	int isim = blockIdx.y ;// * blockDim.y + threadIdx.y ;

	if ( ivar < constNbLevelVar*data.npos ) {	
		/* iniv position entre 0 et ntlevVar */ 

		/* niveau à une position */
		int iniv = ivar % constNbLevelVar ;
		/* position teste */
		int ipos = ivar / constNbLevelVar ;

		//for (int isim=0;isim<data.nsim;isim++) {
		int index = ipos + isim*data.npos + data.constCorrIndexCol[constNbLevelFix+iniv]*data.nsim*data.npos ;
		DT v=0;
		if (triang[GETACC(data,0,ipos,data.constCorrIndexCol[constNbLevelFix+iniv],data.constCorrIndexCol[constNbLevelFix+iniv])]>data.seuil_cho) {
			for (int kd = 0;kd<data.nd;kd++ ) {
				v += data.contingence_var[kd*constNbLevelVar*data.npos+ivar]*data.Y[isim+kd*data.nsim];
			}
			rhs[index] = v;
		}
		//}
	}
}


void QTLMapHomoscedasticModelCalcul::calcul_RHS(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work) {
	dim3 dimBlock_Set_RHS_FIX(16,32);
	dim3 dimGrid_Set_RHS_FIX (ceil(data.npos*data.nLevelVar /dimBlock_Set_RHS_FIX.x ) + 1 ,ceil(data.nsim / dimBlock_Set_RHS_FIX.y )+1);
	Set_RHS_FIX_homoscedatic<<<dimGrid_Set_RHS_FIX,dimBlock_Set_RHS_FIX,0,stream>>>(data,work.rhs,work.triang);
	dim3 dimGrid_Set_RHS_VAR (ceil(data.npos*data.nLevelVar / MAX_BLOCKDIM_64) + 1 ,data.nsim);
	Set_RHS_VAR_homoscedatic<<<dimGrid_Set_RHS_VAR,MAX_BLOCKDIM_64,0,stream>>>(data,work.rhs,work.triang);
}

//***********************************************************************************************************************************************

__global__ void Set_RHS_FIX_homoscedatic_animal(QTLMapStructDeviceDataLinear data,DT *FixM,DT * M,DT * rhs,DT * triang) {
	int ipos = blockIdx.x * blockDim.x + threadIdx.x ; /* position */
	int isim = blockIdx.y * blockDim.y + threadIdx.y ;// * blockDim.y + threadIdx.y ;

	int idNiv;
	if ( ipos < data.npos && isim < data.nsim ) {
		for (int iniv=0;iniv<constNbLevelFix;iniv++) {

			idNiv = data.constCorrIndexCol[iniv];
			int index = ipos + isim*data.npos + idNiv*data.nsim*data.npos ;
			DT v2=0;
			if (triang[GETACC(data,0,ipos,idNiv,idNiv)]>data.seuil_cho) {	
				for (int kd = 0;kd<data.nd;kd++ ) {
					v2 += FixM[kd*constNbLevelFix+iniv]*data.Y[isim+kd*data.nsim];
				}	
				rhs[index] = v2 ;
			}
		}
	}
}

__global__ void Set_RHS_VAR_homoscedatic_animal_1(QTLMapStructDeviceDataLinear data,DT * M,DT * tempRHS) {
	int ivar = blockIdx.x * blockDim.x + threadIdx.x ;
	int kd   = blockIdx.y ;

	if ( ivar < constNbLevelVar*data.npos && kd < data.nd ) {
		DT v=0.0;
		for (int ikd=0;ikd<data.nd;ikd++) {
			v += data.contingence_var[ikd*constNbLevelVar*data.npos+ivar] * M[ikd*data.nd+kd];
		}
		tempRHS[kd*constNbLevelVar*data.npos + ivar] = v;
	}			
}

__global__ void Set_RHS_VAR_homoscedatic_animal_2(QTLMapStructDeviceDataLinear data,DT * tempRHS,DT * rhs,DT * triang) {
	int ivar = blockIdx.x * blockDim.x + threadIdx.x ; /* position */
	int isim = blockIdx.y ;// * blockDim.y + threadIdx.y ;

	if ( ivar < constNbLevelVar*data.npos && isim < data.nsim ) {	
		/* iniv position entre 0 et ntlevVar */ 

		/* niveau à une position */
		int iniv = ivar % constNbLevelVar ;
		/* position teste */
		int ipos = ivar / constNbLevelVar ;

		//for (int isim=0;isim<data.nsim;isim++) {
		int index = ipos + isim*data.npos + data.constCorrIndexCol[constNbLevelFix+iniv]*data.nsim*data.npos ;
		DT v=0.0;
		if (triang[GETACC(data,0,ipos,data.constCorrIndexCol[constNbLevelFix+iniv],data.constCorrIndexCol[constNbLevelFix+iniv])]>data.seuil_cho) {
			for (int kd = 0;kd<data.nd;kd++ ) {
				v += tempRHS[kd*constNbLevelVar*data.npos+ivar]*data.Y[isim+kd*data.nsim];
			}
			rhs[index] = v;
		}
		//}
	}
}


void QTLMapHomoscedasticAnimalModelCalcul::calcul_RHS(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work) {
	dim3 dimBlock_Set_RHS_FIX(16,32);
	dim3 dimGrid_Set_RHS_FIX (ceil(data.npos*data.nLevelVar /dimBlock_Set_RHS_FIX.x ) + 1 ,ceil(data.nsim / dimBlock_Set_RHS_FIX.y )+1);
	Set_RHS_FIX_homoscedatic_animal<<<dimGrid_Set_RHS_FIX,dimBlock_Set_RHS_FIX,0,stream>>>(data,FixM,M,work.rhs,work.triang);

	dim3 dimBlock_1(MAX_BLOCKDIM_64);
	dim3 dimGrid_1 (ceil(data.npos*data.nLevelVar / dimBlock_1.x) + 1 ,ceil(data.nd / dimBlock_1.y )+1);
	size_t size = data.nLevelVar*data.npos*data.nd*sizeof(DT);
	DT * tempRHS = NULL ;
	safecall(cudaMalloc(&tempRHS,size));
	Set_RHS_VAR_homoscedatic_animal_1<<<dimGrid_1,dimBlock_1,0,stream>>>(data,M,tempRHS); 

	dim3 dimBlock_2(MAX_BLOCKDIM_64);
	dim3 dimGrid_2(ceil(data.npos*data.nLevelVar / dimBlock_2.x) + 1 ,ceil(data.nsim / dimBlock_2.y )+1);
	Set_RHS_VAR_homoscedatic_animal_2<<<dimGrid_2,dimBlock_2,0,stream>>>(data,tempRHS,work.rhs,work.triang);
	safecall(cudaFree(tempRHS));

}


//***********************************************************************************************************************************************

__global__ void Set_RHS_FIX_heteroscedastic(DT * varInDevice,QTLMapStructDeviceDataLinear data,DT * rhs,DT * triang) {
	int ipos = blockIdx.x * blockDim.x + threadIdx.x ; /* position */
	int isim = blockIdx.y * blockDim.y + threadIdx.y ;// * blockDim.y + threadIdx.y ;

	int idNiv;
	if ( ipos < data.npos && isim < data.nsim ) {
		for (int iniv=0;iniv<constNbLevelFix;iniv++) {

			idNiv = data.constCorrIndexCol[iniv];

			int index = ipos + isim*data.npos + idNiv*data.nsim*data.npos ;
			rhs[index]=0;
			if (triang[GETACC(data,isim,ipos,idNiv,idNiv)]>data.seuil_cho) {
				for (int kd=0;kd<data.nd;kd++ ) {
					rhs[index] += data.contingence_fix[iniv*data.ndmax + kd]*data.Y[isim+kd*data.nsim]*(data.CD[kd] / varInDevice[ipos+isim*data.npos+data.corIpKd[kd]*data.nsim*data.npos]);
				}
			}
		}
	}
}

__global__ void Set_RHS_VAR_heteroscedastic(DT * varInDevice,QTLMapStructDeviceDataLinear data,DT * rhs,DT * triang) {
	int ivar = blockIdx.x * blockDim.x + threadIdx.x ; /* position */
	int isim = blockIdx.y * blockDim.y + threadIdx.y ;// * blockDim.y + threadIdx.y ;

	if ( ivar < constNbLevelVar*data.npos && isim < data.nsim ) {	
		/* iniv position entre 0 et ntlevVar */ 

		/* niveau à une position */
		int iniv = ivar % constNbLevelVar ;
		/* position teste */
		int ipos = ivar / constNbLevelVar ;

		int index = ipos + isim*data.npos + data.constCorrIndexCol[constNbLevelFix+iniv]*data.nsim*data.npos ;

		rhs[index]=0;

		if (triang[GETACC(data,isim,ipos,data.constCorrIndexCol[constNbLevelFix+iniv],data.constCorrIndexCol[constNbLevelFix+iniv])]>data.seuil_cho) {
			for (int kd = 0;kd<data.nd;kd++ ) {
				rhs[index] += data.contingence_var[kd*constNbLevelVar*data.npos+ivar]*data.Y[isim+kd*data.nsim]*(data.CD[kd] / varInDevice[ipos+isim*data.npos+data.corIpKd[kd]*data.nsim*data.npos]);
			}
		}
	}
}

void QTLMapHeteroscedasticModelCalcul::calcul_RHS(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work) {
	dim3 dimBlock_Set_RHS_FIX(16,32);
	dim3 dimGrid_Set_RHS_FIX (ceil(data.npos*data.nLevelVar /dimBlock_Set_RHS_FIX.x ) + 1 ,ceil(data.nsim / dimBlock_Set_RHS_FIX.y )+1);
	Set_RHS_FIX_heteroscedastic<<<dimGrid_Set_RHS_FIX,dimBlock_Set_RHS_FIX,0,stream>>>(varInDevice,data,work.rhs,work.triang);

	int nbBlockY=1;
	while ( work.prop.maxGridSize[1] < ceil(data.nsim / nbBlockY )+1 ) {
		nbBlockY++;
	}
	dim3 dimBlock(MAX_BLOCKDIM_64,nbBlockY);
	dim3 dimGrid_Set_RHS_VAR (ceil(data.npos*data.nLevelVar / MAX_BLOCKDIM_64) + 1 ,ceil(data.nsim / nbBlockY )+1);
	Set_RHS_VAR_heteroscedastic<<<dimGrid_Set_RHS_VAR,dimBlock,0,stream>>>(varInDevice,data,work.rhs,work.triang);
}

/**
 *
 * Resolution LU (Partie descendante descendante )
 *
 *
 */
__global__ void Resolve_LU_homoscedastic(QTLMapStructDeviceDataLinear data, DT * rhs, DT * triang,DT * bestim) {
	int ipos = blockIdx.x * blockDim.x + threadIdx.x ; /* position */
	int isim = blockIdx.y * blockDim.y + threadIdx.y ;

	int ntniv = data.ntniv ;
	DT temp[MAXNTNIV_LOCALTHREAD_MEMORY];

	/******* DESCENT LU *****/
	if ( ipos < data.npos && isim < data.nsim ) {
		//for (int isim=0;isim<data.nsim;isim++) {
		int scale = ipos + isim*data.npos ;
		for (int i=0;i<ntniv;i++) {
			temp[i] = 0.0;
			if (triang[GETACC(data,0,ipos,i,i)]>data.seuil_cho) {
				temp[i] = rhs[scale + i*data.nsim*data.npos];
				for (int j=i-1;j>=0;j--) {
					if (triang[GETACC(data,0,ipos,j,j)]>data.seuil_cho) {
						temp[i] = temp[i] - temp[j]*triang[GETACC(data,0,ipos,i,j)];
					}
				}
				temp[i] = temp[i] / triang[GETACC(data,0,ipos,i,i)];
			} 
		}

		for (int i=ntniv-1;i>=0;i--) {
			int indexI = scale + i*data.nsim*data.npos;
			bestim[indexI] = 0.0 ;
			if (triang[GETACC(data,0,ipos,i,i)]>data.seuil_cho) {
				bestim[indexI] = temp[i];
				for (int j=i+1;j<ntniv;j++) {
					if (triang[GETACC(data,0,ipos,j,j)]>data.seuil_cho) {
						int indexJ = scale + j*data.nsim*data.npos;
						bestim[indexI] -= bestim[indexJ]*triang[GETACC(data,0,ipos,i,j)];
					}
				}
				bestim[indexI] = bestim[indexI] / triang[GETACC(data,0,ipos,i,i)]; 
			} 
		}
		//}//fin for isim

	}
}


void QTLMapHomoscedasticModelCalcul::calcul_LU(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work,QTLMapStructDeviceSolutionLinear solution) {
	/* Solution */
	dim3 Resolve_LU_dimBlock(16,32);
	dim3 Resolve_LU_dimGrid (ceil(data.npos/Resolve_LU_dimBlock.x)+1 , ceil(data.nsim / Resolve_LU_dimBlock.y)+1);
	Resolve_LU_homoscedastic<<<Resolve_LU_dimGrid,Resolve_LU_dimBlock,0,stream>>>(data,work.rhs,work.triang,solution.bestim);
}


__global__ void Resolve_LU_heteroscedastic(QTLMapStructDeviceDataLinear data, DT * rhs, DT * triang,DT * bestim) {
	int ipos = blockIdx.x * blockDim.x + threadIdx.x ; /* position */
	int isim = blockIdx.y * blockDim.y + threadIdx.y ;

	int ntniv = data.ntniv ;

	DT temp[MAXNTNIV_LOCALTHREAD_MEMORY];

	/******* DESCENT LU *****/
	if ( ipos < data.npos && isim < data.nsim ) {
		int scale = ipos + isim*data.npos ;
		for (int i=0;i<ntniv;i++) {
			temp[i] = 0.0;
			if (triang[GETACC(data,isim,ipos,i,i)]>data.seuil_cho) {
				temp[i] = rhs[scale + i*data.nsim*data.npos];
				for (int j=i-1;j>=0;j--) {
					if (triang[GETACC(data,isim,ipos,j,j)]>data.seuil_cho) {
						temp[i] = temp[i] - temp[j]*triang[GETACC(data,isim,ipos,j,i)];
					}
				}
				temp[i] = temp[i] / triang[GETACC(data,isim,ipos,i,i)];
			} 
		}

		for (int i=ntniv-1;i>=0;i--) {
			int indexI = scale + i*data.nsim*data.npos;
			bestim[indexI] = 0.0 ;
			if (triang[GETACC(data,isim,ipos,i,i)]>data.seuil_cho) {
				bestim[indexI] = temp[i];
				for (int j=i+1;j<ntniv;j++) {
					if (triang[GETACC(data,isim,ipos,j,j)]>data.seuil_cho) {
						int indexJ = scale + j*data.nsim*data.npos;
						bestim[indexI] -= bestim[indexJ]*triang[GETACC(data,isim,ipos,i,j)];
					}
				}
				bestim[indexI] = bestim[indexI] / triang[GETACC(data,isim,ipos,i,i)]; 
			} 
		}
	}
}



void QTLMapHeteroscedasticModelCalcul::calcul_LU(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work,QTLMapStructDeviceSolutionLinear solution) {
	/* Solution */
	dim3 Resolve_LU_dimBlock(16,32);
	dim3 Resolve_LU_dimGrid (ceil(data.npos/Resolve_LU_dimBlock.x)+1 , ceil(data.nsim / Resolve_LU_dimBlock.y)+1);
	Resolve_LU_heteroscedastic<<<Resolve_LU_dimGrid,Resolve_LU_dimBlock,0,stream>>>(data,work.rhs,work.triang,solution.bestim);
}




/**
 *
 * Set_SIGSQ
 * XB = Y - ( X' . Bestim )
 * SIGSQ = SUM(XB^2)
 */
__global__ void Set_SIGSQ_homoscedastic(QTLMapStructDeviceDataLinear data,DT * bestim,DT * osigsq) {
	int ipos = blockIdx.x * blockDim.x + threadIdx.x ; /* position */
	int isim = blockIdx.y * blockDim.y + threadIdx.y ; /* simulation */

	int ntniv = data.ntniv ;

	if ( ipos < data.npos && isim < data.nsim ) {
		//	for (int isim=0;isim<data.nsim;isim++) {
		int scale = ipos + isim*data.npos ;
		DT xb = 0 ;

		for (int kd = 0;kd<data.nd;kd++ ) {
			DT v=0;

#pragma unroll 4
			for (int iniv=0;iniv<constNbLevelFix;iniv++) {
				v += data.contingence_fix[iniv*data.ndmax+kd]*bestim[scale + data.constCorrIndexCol[iniv]*data.nsim*data.npos];
			}

			for (int iniv=constNbLevelFix;iniv<ntniv;iniv++) {	
				v += data.contingence_var[kd*constNbLevelVar*data.npos+ipos*constNbLevelVar+(iniv-constNbLevelFix)]*bestim[scale + data.constCorrIndexCol[iniv]*data.nsim*data.npos];
			}

			v = data.Y[isim+kd*data.nsim] - v;
			xb += v*v ;
		}
		DT v = xb / data.nd ;

#pragma unroll 4
		for (int ip=0;ip<data.np;ip++) {
			osigsq[ipos+isim*data.npos+ip*data.npos*data.nsim] = v ; 
		}
	}
	//}
}

void QTLMapHomoscedasticModelCalcul::calcul_SIGSQ(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work,QTLMapStructDeviceSolutionLinear solution) {
	dim3 dimBlock(16,32);
	dim3 dimGrid (ceil(data.npos/dimBlock.x)+1 , ceil(data.nsim / dimBlock.y)+1);					
	Set_SIGSQ_homoscedastic<<<dimGrid,dimBlock,0,stream>>>(data,solution.bestim,solution.osigsq);
}
//****************************************************************************************************************************************************
__global__ void Set_SIGSQ_homoscedastic_animal(QTLMapStructDeviceDataLinear data,DT * M,DT * bestim,DT * osigsq) {
	int ipos = blockIdx.x * blockDim.x + threadIdx.x ; /* position */
	int isim = blockIdx.y * blockDim.y + threadIdx.y ; /* simulation */

	int ntniv = data.ntniv ;

	if ( ipos < data.npos && isim < data.nsim ) {
		//	for (int isim=0;isim<data.nsim;isim++) {
		int scale = ipos + isim*data.npos ;
		DT xb[5000]  ; /* attention ce 5000 remplace MAXND qui n est plus utilisé....ceci est pour l'experimentation du modele animale....*/

		for (int kd = 0;kd<data.nd;kd++ ) {
			DT v=0;

#pragma unroll 4
			for (int iniv=0;iniv<constNbLevelFix;iniv++) {
				v += data.contingence_fix[iniv*data.ndmax+kd]*bestim[scale + data.constCorrIndexCol[iniv]*data.nsim*data.npos];
			}

			for (int iniv=constNbLevelFix;iniv<ntniv;iniv++) {	
				v += data.contingence_var[kd*constNbLevelVar*data.npos+ipos*constNbLevelVar+(iniv-constNbLevelFix)]*bestim[scale + data.constCorrIndexCol[iniv]*data.nsim*data.npos];
			}

			v = data.Y[isim+kd*data.nsim] - v;
			xb[kd] = v;
		}

		DT tot=0;
		for ( int kd=0;kd<data.nd;kd++  ) {
			DT v=0;
			for (int kd2 = 0;kd2<data.nd;kd2++ ) {
				v += xb[kd2]*M[kd*data.nd+kd2];
			}
			tot += v*xb[kd];
		}

		DT v = tot / data.nd ;

#pragma unroll 4
		for (int ip=0;ip<data.np;ip++) {
			osigsq[ipos+isim*data.npos+ip*data.npos*data.nsim] = v ; 
			//osigsq[ipos+isim*data.npos+ip*data.npos*data.nsim] = 0.5;
		}

	}
	//}
}


//****************************************************************************************************************************************************
__global__ void Set_SIGSQ_homoscedastic_animal_XB(QTLMapStructDeviceDataLinear data,DT * XB,DT * bestim) {
	int ipos = blockIdx.x * blockDim.x + threadIdx.x ; /* position */
	int isim = blockIdx.y * blockDim.y + threadIdx.y ; /* simulation */

	int ntniv = data.ntniv ;

	if ( ipos < data.npos && isim < data.nsim ) {
		int scale = ipos + isim*data.npos ;

		for (int kd = 0;kd<data.nd;kd++ ) {
			DT v=0;

			for (int iniv=0;iniv<constNbLevelFix;iniv++) {
				v += data.contingence_fix[iniv*data.ndmax+kd]*bestim[scale + data.constCorrIndexCol[iniv]*data.nsim*data.npos];
			}

			for (int iniv=constNbLevelFix;iniv<ntniv;iniv++) {	
				v += data.contingence_var[kd*constNbLevelVar*data.npos+ipos*constNbLevelVar+(iniv-constNbLevelFix)]*bestim[scale + data.constCorrIndexCol[iniv]*data.nsim*data.npos];
			}

			XB[kd*data.nsim*data.npos+scale] = data.Y[isim+kd*data.nsim] - v;
		}		
	}
}

__global__ void Set_SIGSQ_homoscedastic_animal_2(QTLMapStructDeviceDataLinear data,DT * M,DT * XB,DT * osigsq) {
	int ipos = blockIdx.x * blockDim.x + threadIdx.x ; /* position */
	int isim = blockIdx.y * blockDim.y + threadIdx.y ; /* simulation */

	if ( ipos < data.npos && isim < data.nsim ) {
		int scale = ipos + isim*data.npos ;

		DT tot=0;
		for ( int kd=0;kd<data.nd;kd++  ) {
			DT v=0;
			for (int kd2 = 0;kd2<data.nd;kd2++ ) {
				v += XB[kd2*data.nsim*data.npos+scale]*M[kd*data.nd+kd2];
				//v += xb_loc[kd2]*M[kd*data.nd+kd2];
			}
			tot += v*XB[kd*data.nsim*data.npos+scale];
			//		tot += v*xb_loc[kd];
		}

		DT v = tot / data.nd ;

		for (int ip=0;ip<data.np;ip++) {
			osigsq[ipos+isim*data.npos+ip*data.npos*data.nsim] = v ; 
		}	
	}
}

void QTLMapHomoscedasticAnimalModelCalcul::calcul_SIGSQ(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work,QTLMapStructDeviceSolutionLinear solution) {
	dim3 dimBlock(16,32);
	dim3 dimGrid (ceil(data.npos/dimBlock.x)+1 , ceil(data.nsim / dimBlock.y)+1);	
	Set_SIGSQ_homoscedastic_animal<<<dimGrid,dimBlock,0,stream>>>(data,M,solution.bestim,solution.osigsq);

	//	int size=data.npos*data.nsim*data.nd*sizeof(DT);
	//	DT * XB = NULL ;
	//	safecall(cudaMalloc(&XB,size));
	//
	//	dim3 dimBlock(32,32);
	//	dim3 dimGrid (ceil(data.npos/dimBlock.x)+1 , ceil(data.nsim / dimBlock.y)+1);
	//
	//	Set_SIGSQ_homoscedastic_animal_XB<<<dimGrid,dimBlock,0,stream>>>(data,XB,solution.bestim);
	//
	//	Set_SIGSQ_homoscedastic_animal_2<<<dimGrid,dimBlock,0,stream>>>(data,M,XB,solution.osigsq);
	//	safecall(cudaFree(XB));

}
//****************************************************************************************************************************************************
/**
 *
 * Set_SIGSQ
 * XB = Y - ( X' . Bestim )
 * SIGSQ = SUM(XB^2)
 */
__global__ void Set_SIGSQ_heteroscedastic(QTLMapStructDeviceDataLinear data,DT * bestim,DT * osigsq) {
	int ipos = blockIdx.x * blockDim.x + threadIdx.x ; /* position */
	int isim = blockIdx.y * blockDim.y + threadIdx.y  ; 

	int ntniv = data.ntniv ;

	if ( ipos < data.npos && isim < data.nsim ) {
		int scale = ipos + isim*data.npos ;

		for (int ip=0;ip<data.np;ip++) {
			DT xb = 0 ; 
			for (int kd = 0;kd<data.nd;kd++ ) {
				if ( data.corIpKd[kd] == ip ) {
					DT v=0;

					for (int iniv=0;iniv<constNbLevelFix;iniv++) {
						v += data.contingence_fix[iniv*data.ndmax+kd]*bestim[scale + data.constCorrIndexCol[iniv]*data.nsim*data.npos];
					}

					for (int iniv=constNbLevelFix;iniv<ntniv;iniv++) {	
						v += data.contingence_var[kd*constNbLevelVar*data.npos+ipos*constNbLevelVar+(iniv-constNbLevelFix)]*bestim[scale + data.constCorrIndexCol[iniv]*data.nsim*data.npos];
					}

					v = data.Y[isim+kd*data.nsim] - v;
					xb += v*v*data.CD[kd]*data.CD[kd] ;
				}
			}
			osigsq[ipos+isim*data.npos+ip*data.npos*data.nsim] = xb / constSizeFam[ip];
		}
	}
}


void QTLMapHeteroscedasticModelCalcul::calcul_SIGSQ(cudaStream_t stream,QTLMapStructDeviceDataLinear data,QTLMapStructDeviceWorkLinear work,QTLMapStructDeviceSolutionLinear solution) {
	dim3 dimBlock(16,32);
	dim3 dimGrid (ceil(data.npos/dimBlock.x)+1 , ceil(data.nsim / dimBlock.y)+1);			
	Set_SIGSQ_heteroscedastic<<<dimGrid,dimBlock,0,stream>>>(data,solution.bestim,solution.osigsq);
	safecall(cudaGetLastError());
}

/**
 *
 *  Set_LRT
 *
 *  Modif Avril 2013 : return LRT sire contribution
 */
__global__ void Set_LRT(QTLMapStructDeviceDataLinear data,DT * osigsq, DT * lrt,int scalePosition) {
	int ipos = blockIdx.x * blockDim.x + threadIdx.x ; /* position */
	int ip   = blockIdx.y ;

	if ( data.nqtl > 0 ) {
		if ( ipos < data.npos && ip < data.np ) {
			for (int isim=0;isim<data.nsim;isim++) {

				int k = ipos + isim * data.npos ;

				for (int iqtl=0;iqtl<data.nqtl;iqtl++) {
					//					DT v = 0.0;
					//					for (int ip=0;ip<data.np;ip++) {
					//						v += constSizeFam[ip] * (log(data.isigsq[isim+ip*data.nsim+(iqtl)*data.nsim*data.np]) - log(osigsq[k+ip*data.npos*data.nsim]));
					//					}
					lrt[scalePosition+ipos+isim*data.nposGlobal+iqtl*data.nposGlobal*data.nsim+ip*data.nqtl*data.nposGlobal*data.nsim] =
							constSizeFam[ip] * (log(data.isigsq[isim+ip*data.nsim+(iqtl)*data.nsim*data.np]) - log(osigsq[k+ip*data.npos*data.nsim]));
				}
			}//fin for isim
		} 
	}
}

extern "C" void release_allocated_internal_sig_() 
{
	if (allVarFitted != NULL ) delete[] allVarFitted;
	allVarFitted = NULL;
}

/**
 * - save residual variance at the maximum position of the likelihood. this information is useful to compute LRT
 */
void add_internal_sig(DT * osigsquare,int *maxPosition,int npos,int nsim,int np,int nqtl) {
	DT * allVarFitted_buf = NULL ;
	if ( allVarFitted ) {
		allVarFitted_buf = new DT [nqtl*np*nsim];
		for (int iq=0;iq<nqtl;iq++) {
			for (int isim=0;isim<nsim;isim++) {
				for (int ip=0;ip<np;ip++)
					allVarFitted_buf[isim+ip*nsim+iq*np*nsim] = allVarFitted[isim+ip*nsim+iq*np*nsim];
			}
		}
		delete[] allVarFitted;
		allVarFitted = NULL;
	}

	allVarFitted = new DT [(nqtl+1)*np*nsim];
	if ( nqtl >= 1 ) {
		for (int iq=0;iq<nqtl;iq++) {
			for (int isim=0;isim<nsim;isim++) {
				for (int ip=0;ip<np;ip++) {
					allVarFitted[isim+ip*nsim+iq*np*nsim] = allVarFitted_buf[isim+ip*nsim+iq*np*nsim];
				}
			}
		}
	}

	if (allVarFitted_buf) {
		delete [] allVarFitted_buf ;
		allVarFitted_buf = NULL;
	}

	/* update with the new values osig computed */
	for (int isim=0;isim<nsim;isim++) {
		if ( maxPosition[isim] >= 1 ) { /* index fortran */
			int ipos = maxPosition[isim] - 1;
			for (int ip=0;ip<np;ip++) {
				//			cout << osigsquare[ipos+isim*npos+ip*nsim*npos] << endl ;
				allVarFitted[isim+ip*nsim+nqtl*np*nsim] = osigsquare[ipos+isim*npos+ip*nsim*npos];	
			}
		}
	}
}



/*
  xinc  = matrice d'incidence
  nkd   = nombre de ligne
  ntniv = nombre de colonne
  seuil_choPtr = seuil d'estimabilite pour la solution

  column major order ===> offset = row + column*NUMROWS
  XINC ACCESS => ndmax*npos*iniv+npos*kd+ipos

  XINC(construction : level fixe puis levels variables)
  Mu MuPere1 MuPere2 MuMere1 MuMere2 QTLPere1Pos1 QTLPere2Pos1 QTLMere1Pos1 QTLMere2Pos1 ... QTLPere1PosN QTLPere2PosN QTLMere1PosN QTLMere2PosN

  nLevelFixPtr : 5 (Moyenne gen + 4 effet polygenique)
  nLevelVarPtr : 4 (4 effet qtl a estimer)
  corrLevelColPtr(1) = 1    moyenne generale
  corrLevelColPtr(2) = 6    effet qtl pere 1
  corrLevelColPtr(3) = 7    effet qtl pere 2
  corrLevelColPtr(4) = 8    effet qtl mere 1
  corrLevelColPtr(5) = 9    effet qtl mere 2
  corrLevelColPtr(6) = 2    effet polygenique pere 1
  corrLevelColPtr(7) = 3    effet polygenique pere 2
  corrLevelColPtr(8) = 4    effet polygenique mere 1
  corrLevelColPtr(9) = 5    effet polygenique mere 2
 */
extern "C" void cuda_model_resolv_genome_(
		int *gpu_device_id,                   /* IN  : The id device to used */
		int *heteroscedastic_mode,            /* IN  : heteroscedastic = 1 , homoscedastic = 0,  homoscedastic_animal = 2 */
		int *nqtlPtr,                         /* IN  : Number of QTL */
		DT  *sigsquare,                       /* IN  : variance residuelle sous les hypothese 0 à NQTL-1 tableau [NSIM][NP] */
		int *nLevelFixPtr,                    /* IN  : nombre de niveaux fixe a la position : moyenne genrale, polygenic, effet fixe, covariables */
		int *nLevelVarPtr,                    /* IN  : nombre de niveau variable a la position */
		int *corrLevelColPtr,                 /* IN  : Correspondance du niveau associe a une colonne de la matrice de contingence */
		void *work_cuda,                      /* IN  : structure de travail pour la fonction get_partialXinc */
		void (*get_partialXincFix)(void*,DT*),   /* IN  : Fonction qui donne les donnees fixe a la position de la matrice d incidences... */
		void (*get_partialXincVar)(void*,int*,int*,DT*),            /* IN  : Fonction qui donne les donnees variables a la position des matrices d incidences... */
		DT *Y_d,                              /* IN  : performance [nsim,nd] */
		DT *CD_d,                             /* IN  : cd of anim  [nd]      */
		DT *M_d ,                             /* IN  : matrix quantity (I - (I -lambda. A ** -1))  [nd][nd]      */
		int *nsimPtr,                         /* IN  : nombre de simulation */
		int *ndPtr,                           /* IN  : nombre d'individu total  : nombre de ligne reelle de XINC */
		int *nkdPtr,                          /* IN  : nombre d'individu pris en compte dans les matrices de contingence */
		int *npPtr,                           /* IN  : nombre de famille de pere */
		int *sizeFamilyNp,                    /* IN  : nombre de descendant par pere tableau [NP] */
		int *ntnivmaxPtr,                     /* IN  : nombre maximum de niveau : nombre de colonne reelle de XINC */
		int *npositionPtr,                    /* IN  : number of position estimation */
		double  *seuil_choPtr,                /* IN  : seuil pour considere grace a la decomposition de cholesky si un niveau est estimable */
		//		int    * vecsol,                      /* OUT : tableau de booleen de taille [NPOSITION,NTNIV]  pour l'estimabilite de chaque niveau */
		DT * bestim,                          /* OUT : tableau solution de l'estimation de chaque niveau (taille [NPOSITION,NSIM,NTNIV]) */
		DT * osigsquare,                      /* OUT : tableau des variances residuelles pour chaque familles de pere (taille [NPOSITION,NSIM,NP])  */
		DT * lrt,                             /* OUT : LRT (taille [NPOSITION,NSIM,NQTL,NP])  */
		DT * maxLRT,                          /* OUT : Value of the Maximul LRT  [NQTL,NSIM] */
		int *maxPosition                      /* OUT : Index of the position maximum tableau [NSIM] */
) {
	/* reset the cards to remove leaks */
	safecall(cudaDeviceReset());
#if(_CUDA_HOST_DEBUG_)
	cout << "GPU_DEVICE_ID               " <<  gpu_device_id << "="<< *gpu_device_id << endl;
	cout << "HETEROSCEDACTIC_MODE        " <<  heteroscedastic_mode <<"="<< *heteroscedastic_mode << endl;
	cout << "NQTL                        " <<  nqtlPtr << "="<<*nqtlPtr << endl;
	cout << "SISQUARE                    " <<  sigsquare << "="<<*sigsquare << endl;
	cout << "nLevelFixPtr                " <<  nLevelFixPtr << "="<<*nLevelFixPtr << endl;
	cout << "nLevelnLevelVarPtrFixPtr    " <<  nLevelVarPtr << "="<<*nLevelVarPtr << endl;
	cout << "corrLevelColPtr             " <<  corrLevelColPtr << "="<<*corrLevelColPtr << endl;
	cout << "work_cuda                   " <<  work_cuda << endl;
	cout << "bestim                      " <<  bestim << "="<<*bestim << endl;
	cout << "lrt                         " <<  lrt << "="<<*lrt << endl;
	cout << "osigsquare                  " <<  osigsquare << "="<<*osigsquare << endl;
	cout << "nsimPtr                     " <<  nsimPtr << "="<<*nsimPtr << endl;
#endif
	int  mode = *heteroscedastic_mode;

	QTLMapStructDeviceDataLinear      data ;
	QTLMapStructDeviceWorkLinear      work ;
	QTLMapStructDeviceSolutionLinear  solution ;

	cout << "****************** INFO DEVICES ******************************" << endl ;
	int nbDevice ;
	cudaGetDeviceCount(&nbDevice);
	cout << "number of device    : " << nbDevice << endl ;

	cudaDeviceProp  prop ;

	for (int i=0;i<nbDevice;i++) {
		cudaGetDeviceProperties(&prop,i);
		cout << "id:"<< i << " ******* ===  "<< prop.name << "  ===  *******"<< endl ;
	}
	cout << "***************************************************************" << endl ;

	if ( *gpu_device_id >= nbDevice || *gpu_device_id < 0 ) {
		cerr << " bad gpu_divice_id : " << *gpu_device_id << endl ;
		exit(1);
	}

	/* Choose the device */
	cudaGetDeviceProperties(&work.prop,*gpu_device_id);
	safecall(cudaSetDevice(*gpu_device_id));


	cout << "Device Num:"<<  gpu_device_id << " - "<< work.prop.name << endl ;
	cout << "  totalGlobalMem : "      << prop.totalGlobalMem     << endl ;
	cout << "  sharedMemPerBlock : "   << prop.sharedMemPerBlock  << endl ;
	cout << "  regsPerBlock : "        << prop.regsPerBlock       << endl ; 
	cout << "  warpSize : "            << prop.warpSize           << endl ;
	cout << "  maxThreadsPerBlock : "  << prop.maxThreadsPerBlock << endl;
	cout << "  maxThreadsDim : "       << prop.maxThreadsDim[0]<< "," << prop.maxThreadsDim[1] << "," << prop.maxThreadsDim[2] << endl ;
	cout << "  maxGridSize : "         << prop.maxGridSize[0] << "," << prop.maxGridSize[1] << "," << prop.maxGridSize[2] << endl ;
	cout << "  clockRate : " << prop.clockRate << endl ;
	cout << "  totalConstMem : " << prop.totalConstMem << endl ;
	cout << "  multiProcessorCount : " << prop.multiProcessorCount << endl ;

	safecall(cudaMemcpyToSymbol(constNbLevelFix, nLevelFixPtr, sizeof(int)));
	safecall(cudaMemcpyToSymbol(constNbLevelVar, nLevelVarPtr, sizeof(int)));
	safecall(cudaMemcpyToSymbol(constNSim, nsimPtr, sizeof(int)));
	/*
	if ( ( (*nLevelFixPtr) + (*nLevelVarPtr)) >= MAXNTNIV_LOCALTHREAD_MEMORY ) {
		cerr << "Can not initialized cuda structure => update constant MAXNTNIV_LOCALTHREAD_MEMORY :"<<  MAXNTNIV_LOCALTHREAD_MEMORY << "<=" <<  (*nLevelFixPtr + *nLevelVarPtr) << endl ;
		exit(1);
	}	
	 */
	//	safecall(cudaMemcpyToSymbol(constCorrIndexCol, corrLevelColPtr, (*nLevelFixPtr + *nLevelVarPtr)*sizeof(DT)));
	size_t size=(*nLevelFixPtr + *nLevelVarPtr)*sizeof(DT);
	safecall(cudaMalloc(&data.constCorrIndexCol,size));
	size_t size2=(*nLevelFixPtr + *nLevelVarPtr)*sizeof(DT);
	safecall(cudaMemcpy(data.constCorrIndexCol, corrLevelColPtr, size2,cudaMemcpyHostToDevice));

	if ( (*npPtr) >=  MAXNP_WORK ) {
		cerr << "Can not initialized cuda structure => update constant MAXNP_WORK :"<<  MAXNP_WORK << "<=" <<  (*npPtr) << endl ;
		exit(1);
	}

	safecall(cudaMemcpyToSymbol(constSizeFam,sizeFamilyNp,(*npPtr)*sizeof(int)));



	int err ;

	DT *xinc_d_var,*xinc_d_fix;

	/* Initialisation des donnees fixes a la position de la matrice de contingence */
	safecall(cudaMallocHost(& xinc_d_fix, (*ndPtr)*((*nLevelFixPtr)*sizeof(DT))));
	get_partialXincFix(work_cuda,xinc_d_fix);

	data.init(heteroscedastic_mode,nqtlPtr,
			sigsquare,xinc_d_fix,Y_d,
			CD_d,corrLevelColPtr,seuil_choPtr,
			ndPtr,nkdPtr,nsimPtr,npositionPtr,
			nLevelFixPtr,nLevelVarPtr,npPtr,sizeFamilyNp);




	int nbPosBlock;
	nbPosBlock = data.calculBlockPositionWorkSize(mode);
	if ( nbPosBlock > *npositionPtr ) {
		nbPosBlock = *npositionPtr ;
	} 

	data.npos = nbPosBlock ;
	work.initResolution(mode,data,CD_d,sigsquare,sizeFamilyNp) ;
	solution.init(data);

	/*
	 * 
	 *    Boucle des traitements par block
	 *    Utilisation de stream cuda
	 * 
	 * 
	 */
	/* Variables de travail pour recolter les resultats intermediaires */
	DT *local_osigsquare;
	DT *local_bestim;//[nbPosBlock*data.nsim*data.ntniv];

	safecall(cudaMallocHost(& local_osigsquare, nbPosBlock*data.nsim*data.np*sizeof(DT)));
	safecall(cudaMallocHost(& local_bestim, nbPosBlock*data.nsim*data.ntniv*sizeof(DT)));

	/* Creation des streams */
	int nstream=2;
	cudaStream_t stream[nstream];
	for (int i=0;i<nstream;i++) {
		cudaStreamCreate(&stream[i]);
	}
	int istream = 0 ;
	/* premier appel */

	int sizeBlockXincDVar = 0; /* taille en nombre de position a teste sur le groupe de liaison */
	int nbBlockXincDVar=2; /* On gere 2 block .l'un a remplir, l autre pour le calcul */
	int currentBlockMod = 0; 

	if ((*nLevelVarPtr) > 0) {
		/* Initialisation des donnees variables a la position de la matrice de contingence */
		sizeBlockXincDVar = (*ndPtr)*((*nLevelVarPtr)*nbPosBlock) ;
		safecall(cudaMallocHost(& xinc_d_var, nbBlockXincDVar*sizeBlockXincDVar*sizeof(DT)));
	}

	int lastPosition = 1;
	int nextPosition = nbPosBlock;
	get_partialXincVar(work_cuda,&lastPosition,&nextPosition,xinc_d_var);
	print_info_memory();

	QTLMapGenericModelCalcul *model ;

	switch (mode) {
	case QTLMapStructDeviceDataLinear::MODEL_HETERO_POLYGENIC :
		model = new QTLMapHeteroscedasticModelCalcul(data.npos,data.nsim,data.np);
		break;
	case QTLMapStructDeviceDataLinear::MODEL_HOMO_POLYGENIC :
		model = new QTLMapHomoscedasticModelCalcul();
		break;
	case QTLMapStructDeviceDataLinear::MODEL_HOMO_ANIMAL :
		model = new QTLMapHomoscedasticAnimalModelCalcul(data.nd,data.nLevelFix,M_d);
		break;
	default :
		cerr << "Devel error . not implementation of mode :"<< mode << endl ;
		exit(1);
	}

	/*
	 * Initialisation pour benchmark
	 */	
#if(_CUDA_HOST_TIME_PROF_)
	cudaEvent_t start, stop;
	float elapsedTime;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
#endif

	for (int lastPosition = 1;lastPosition<=data.nposGlobal;lastPosition += nbPosBlock){   
#pragma omp sections
		{
#pragma omp section
			{
				int lastPositionNext = lastPosition + nbPosBlock ;
				if ( lastPositionNext <=data.nposGlobal ) {
					int nextPosition = lastPositionNext + nbPosBlock - 1;
					if ( nextPosition > data.nposGlobal) {
						nextPosition = data.nposGlobal;
					}

					int nextBlockMod = (currentBlockMod+1)%nbBlockXincDVar;
					get_partialXincVar(work_cuda,&lastPositionNext,&nextPosition,xinc_d_var+(nextBlockMod)*sizeBlockXincDVar);
				}
			}
#pragma omp section
			{
				nextPosition = lastPosition + nbPosBlock - 1;
				if ( nextPosition > data.nposGlobal) {
					nextPosition = data.nposGlobal;
				}
				data.npos = nextPosition - lastPosition + 1 ;

				cout << "   ** computing from position:"<< lastPosition << " at:"<< nextPosition << " (NPOS:"<< data.nposGlobal << ")   **" << endl << endl;

				/* Type de block utilise par les kernels */
				int nbblock_512 = ceil(data.npos / MAX_BLOCKDIM_512) + 1 ;
				int nbblock_64  = ceil(data.npos / MAX_BLOCKDIM_64) + 1 ;

				data.set_contingence_var(*ndPtr,*nLevelVarPtr,data.npos,0,xinc_d_var+currentBlockMod*sizeBlockXincDVar,&stream[istream]);//

				err = (int)cuInit(0);
				if ( err ) {
					printf(" cuInit error : %d \n",err);
					exit(1);
				}

				do {   /*  START IF CONVERGENCE IS OK */
					/* initializing model */
					model->start_analysis(stream[istream],data,work);
					cudaStreamSynchronize(stream[istream]);
#if(_CUDA_HOST_TIME_PROF_)
					cudaEventRecord(start, 0);
#endif		
					model->calcul_XT_X_A(stream[istream],data,work);
					/* on peut peut etre enlever cette synchronisation et la mettre plus tard */
					//cudaStreamSynchronize(stream[istream]);
#if(_CUDA_HOST_TIME_PROF_)
					cudaStreamSynchronize(stream[istream]);
					cudaEventRecord(stop, 0);
					cudaEventElapsedTime(&elapsedTime, start, stop);
					cout << "----------------------------------------------------------------------------------------------------------------"<< endl ;
					cout << "[BENCHMARK]  XT_X_A grid(x,y,z): elapsedtime:"<< elapsedTime<< endl ;
					cout << "----------------------------------------------------------------------------------------------------------------"<< endl ;
#endif
#if(_CUDA_HOST_TIME_PROF_)
					cudaEventRecord(start, 0);
#endif				
					model->calcul_XT_X_B(stream[istream],data,work);
					/* on peut peut etre enlever cette synchronisation et la mettre plus tard */
					//cudaStreamSynchronize(stream[istream]);
					safecall(cudaGetLastError());

#if(_CUDA_HOST_TIME_PROF_)
					cudaStreamSynchronize(stream[istream]);
					cudaEventRecord(stop, 0);
					cudaEventElapsedTime(&elapsedTime, start, stop);
					cout << "----------------------------------------------------------------------------------------------------------------"<< endl ;
					cout << "[BENCHMARK]  XT_X_B grid(x,y,z): elapsedtime:"<< elapsedTime<< endl ;
					cout << "----------------------------------------------------------------------------------------------------------------"<< endl ;
#endif


#if(_CUDA_HOST_TIME_PROF_)
					cudaEventRecord(start, 0);
#endif
					model->calcul_XT_X_C(stream[istream],data,work);
					cudaStreamSynchronize(stream[istream]);
					safecall(cudaGetLastError());


#if(_CUDA_HOST_TIME_PROF_)
					cudaEventRecord(stop, 0);
					cudaEventElapsedTime(&elapsedTime, start, stop);
					cout << "----------------------------------------------------------------------------------------------------------------"<< endl ;
					cout << "[BENCHMARK]  XT_X_C grid(x,y,z) elapsedtime:"<< elapsedTime<< endl ;
					cout << "----------------------------------------------------------------------------------------------------------------"<< endl ;
#endif

#if(_CUDA_HOST_DEBUG_)
					printMatSymXXTriang(model->getType(),data,work.XX);
#endif

#if(_CUDA_HOST_TIME_PROF_)
					cudaEventRecord(start, 0);
#endif
					/* Cholesky decomp */
					model->calcul_Cholesky_Decomposition(stream[istream],data,work);
					cudaStreamSynchronize(stream[istream]);
					safecall(cudaGetLastError());

					//	printMatSymXXTriang(model->getType(),data,work.triang);

#if(_CUDA_HOST_TIME_PROF_)
					cudaEventRecord(stop, 0);
					cudaEventElapsedTime(&elapsedTime, start, stop);
					cout << "----------------------------------------------------------------------------------------------------------------"<< endl ;
					cout << "[BENCHMARK]  Cholesky_Decomposition grid(x,y,z) elapsedtime:"<< elapsedTime<< endl ;
					cout << "----------------------------------------------------------------------------------------------------------------"<< endl ;
#endif

#if(_CUDA_HOST_DEBUG_)
					printf("Cholesky decomp\n----------------------------------------------------------\n");
					printMatSymXXTriang(model->getType(),data,work.triang);
#endif

#if(_CUDA_HOST_TIME_PROF_)
					cudaEventRecord(start, 0);
#endif								
					model->calcul_RHS(stream[istream],data,work);
					cudaStreamSynchronize(stream[istream]);
					safecall(cudaGetLastError());

#if(_CUDA_HOST_TIME_PROF_)
					cudaEventRecord(stop, 0);
					cudaEventElapsedTime(&elapsedTime, start, stop);
					cout << "----------------------------------------------------------------------------------------------------------------"<< endl ;
					cout << "[BENCHMARK]  Set_RHS_FIX + Set_RHS_VAR grid(x,y,z):"<< " elapsedtime:"<< elapsedTime<< endl ;
					cout << "----------------------------------------------------------------------------------------------------------------"<< endl ;
#endif


#if(_CUDA_HOST_DEBUG_)
					printf("RHS\n-------POS=1,SIM=nsim => NTNIV---------------------------------------------------\n");
					Utils<DT>::printFloatDeviceArray3D(data.npos,data.nsim,data.ntniv,10,1,data.ntniv,work.rhs);
#endif		

#if(_CUDA_HOST_TIME_PROF_)
					cudaEventRecord(start, 0);
#endif
					/* Solution */
					model->calcul_LU(stream[istream],data,work,solution);
					cudaStreamSynchronize(stream[istream]);
					safecall(cudaGetLastError());

#if(_CUDA_HOST_TIME_PROF_)
					cudaEventRecord(stop, 0);
					cudaEventElapsedTime(&elapsedTime, start, stop);
					cout << "----------------------------------------------------------------------------------------------------------------"<< endl ;
					cout << "[BENCHMARK]  Resolve_LU grid(x,y,z): elapsedtime:"<< elapsedTime<< endl ;
					cout << "----------------------------------------------------------------------------------------------------------------"<< endl ;
#endif


#if(_CUDA_HOST_DEBUG_)
					printf("Resolution LU\n------POS=1,SIM=nsim => NTNIV ----------------------------------------------------\n");
					Utils<DT>::printFloatDeviceArray3D(data.npos,data.nsim,data.ntniv,4,1,data.ntniv,solution.bestim);
#endif

					/* Recuperation des resultats */
					Utils<DT>::getArrayDeviceToHost(nbblock_512,MAX_BLOCKDIM_512,data.npos*data.nsim*data.ntniv,solution.bestim,local_bestim,&stream[istream]);

#if(_CUDA_HOST_TIME_PROF_)
					cudaEventRecord(start, 0);
#endif
					/* Variances residuelles */
					model->calcul_SIGSQ(stream[istream],data,work,solution);
					/* on synchronise pour que local_bestim soit initialise */
					cudaStreamSynchronize(stream[istream]);

#if(_CUDA_HOST_TIME_PROF_)
					cudaEventRecord(stop, 0);
					cudaEventElapsedTime(&elapsedTime, start, stop);
					cout << "----------------------------------------------------------------------------------------------------------------"<< endl ;
					cout << "[BENCHMARK]  Set_SIGSQ grid(x,y,z)  elapsedtime:"<< elapsedTime<< endl ;
					cout << "----------------------------------------------------------------------------------------------------------------"<< endl ;
#endif

					for (int ipos=0;ipos<data.npos;ipos++) {
						for (int isim=0;isim<data.nsim;isim++) {
							for (int iniv=0;iniv<data.ntniv;iniv++) {
								bestim[ipos+lastPosition-1+isim*data.nposGlobal+iniv*data.nsim*data.nposGlobal] = local_bestim[ipos+isim*data.npos+iniv*data.nsim*data.npos] ;
							}
						}
					}

					cudaStreamSynchronize(stream[istream]);
					safecall(cudaGetLastError());

#if(_CUDA_HOST_DEBUG_)
					cout << "SIGSQ------------POS=1,SIM=nsim => NP----------------------------------------------"<<endl ;
					Utils<DT>::printFloatDeviceArray3D(data.npos,data.nsim,data.np,data.npos,1,data.np,solution.osigsq);
#endif
					/* Recuperation des resultats */
					Utils<DT>::getArrayDeviceToHost(nbblock_512,MAX_BLOCKDIM_512,data.npos*data.nsim*data.np,solution.osigsq,local_osigsquare,&stream[istream]);



#if(_CUDA_HOST_TIME_PROF_)
					cudaEventRecord(start, 0);
#endif
					cudaStreamSynchronize(stream[istream]);
				} while ( ! model->convergenceOk(stream[istream],data,work,solution) ); /* FIN CONVERGENCE */

				if (data.nqtl > 0 ) {
					dim3 dimBlockLrt(32,16);
					dim3 dimGridLrt (ceil(data.npos/dimBlockLrt.x)+1 , data.np );
					/* Calcul LRT */
					Set_LRT<<<dimGridLrt,dimBlockLrt,0,stream[istream]>>>(data,solution.osigsq,solution.lrt,lastPosition-1);

#if(_CUDA_HOST_DEBUG_)
					printf("PARTIEL LRT\n------------POS=1,SIM=1,..,nsim => NP----------------------------------------------\n");
					cudaStreamSynchronize(stream[istream]);
					Utils<DT>::printFloatDeviceArray3D(data.nposGlobal,data.nsim,data.nqtl,data.nposGlobal,1,data.nqtl,solution.lrt);
#endif
					safecall(cudaGetLastError());
					cudaStreamSynchronize(stream[istream]);
				}


#if(_CUDA_HOST_TIME_PROF_)
				cudaEventRecord(stop, 0);
				cudaEventElapsedTime(&elapsedTime, start, stop);
				cout << "----------------------------------------------------------------------------------------------------------------"<< endl ;
				cout << "[BENCHMARK]  Set_LRT grid(x,y,z)  elapsedtime:"<< elapsedTime<< endl ;
				cout << "----------------------------------------------------------------------------------------------------------------"<< endl ;
#endif



				for (int ipos=0;ipos<data.npos;ipos++) {
					for (int isim=0;isim<data.nsim;isim++) {
						for (int ip=0;ip<data.np;ip++) {
							osigsquare[ipos+lastPosition-1+isim*data.nposGlobal+ip*data.nsim*data.nposGlobal] = local_osigsquare[ipos+isim*data.npos+ip*data.nsim*data.npos];//[ip+isim*data.np+ipos*data.nsim*data.np] ;
						}
					}
				}
			}
		}
		istream = (istream+1)%nstream;
		currentBlockMod = (currentBlockMod+1)%nbBlockXincDVar ;

	}

	/* Fin boucle principale  */

	for (int i=0;i<nstream;i++) {
		cudaStreamDestroy(stream[i]);
	}


	work.releaseDeviceResolution();
	safecall(cudaFreeHost(xinc_d_fix));
	if ( data.nqtl>=1 ) safecall(cudaFreeHost(xinc_d_var));
	safecall(cudaFreeHost(local_osigsquare));
	safecall(cudaFreeHost(local_bestim));

	if ( data.nqtl>=1 ) {	
		/* Recuperation des LRTs */
		int nbblock_512_global = ceil(data.nposGlobal / MAX_BLOCKDIM_512) + 1 ;
		Utils<DT>::getArrayDeviceToHost(nbblock_512_global,MAX_BLOCKDIM_512,data.nqtl*data.nposGlobal*data.nsim*data.np,solution.lrt,lrt,NULL);

//		for (int ip=0;ip<data.np;ip++) {
//			cout << "IP="<< ip << endl ;
//			for (int i=0;i<data.nposGlobal;i++) {
//				cout << lrt[i+ip*data.nqtl*data.nposGlobal*data.nsim]<< " ";
//			}
//			cout << endl ;
//		}

		for (int isim=0;isim<data.nsim;isim++) {
			for (int iq=0;iq<data.nqtl;iq++) {
				maxLRT[isim+iq*data.nsim]      = -999.9;
			}
			maxPosition[isim] = -1;


			for (int i=0;i<data.nposGlobal;i++) {
				DT lrtmax = 0 ;
				for (int ip=0;ip<data.np;ip++) {
					lrtmax += lrt[i+isim*data.nposGlobal+(data.nqtl-1)*data.nposGlobal*data.nsim+ip*data.nqtl*data.nposGlobal*data.nsim];
				}
				if ( maxLRT[isim+(data.nqtl-1)*data.nsim] < lrtmax ) {
					for (int iq=0;iq<data.nqtl;iq++) {
						DT lrtmax = 0 ;
						for (int ip=0;ip<data.np;ip++) {
							lrtmax += lrt[i+isim*data.nposGlobal+iq*data.nposGlobal*data.nsim+ip*data.nqtl*data.nposGlobal*data.nsim];
						}
						maxLRT[isim+iq*data.nsim] = lrtmax;
					}
					maxPosition[isim] = i+1 ; /* index pour fortran */
				}
			} // fin i
	} // fin isim

	add_internal_sig(osigsquare,maxPosition,data.nposGlobal,data.nsim,data.np,data.nqtl);
} else { /* residual variance under H0 */
	int *maxP;
	safecall(cudaMallocHost(&maxP,data.nsim*sizeof(DT)));
	for (int isim=0;isim<data.nsim;isim++)maxP[isim]=1;
	add_internal_sig(osigsquare,maxP,data.nposGlobal,data.nsim,data.np,data.nqtl);
	safecall(cudaFreeHost(maxP));
}
delete model;

#if(_CUDA_HOST_TIME_PROF_)
cudaEventDestroy(start);
cudaEventDestroy(stop);
#endif

#if(_CUDA_HOST_DEBUG_)
cout << "DATA RELEASED ON DEVICE" << endl ;
#endif
data.releaseDevice();

#if(_CUDA_HOST_DEBUG_)
cout << "SOLUTION RELEASED ON DEVICE" << endl ;
#endif
solution.releaseDevice();

#if(_CUDA_HOST_DEBUG_)
cout << " ** END ** " << endl ;
#endif

/* Modif OFI Mai 2012 constCorrIndexCol n est plus dans la memoire constante (ne fonctionnait pas quand le nombre de niveau atteignait 160)
 * cet appel devra etre mis dans data.releasedevice() */
safecall(cudaFree(data.constCorrIndexCol));

}
