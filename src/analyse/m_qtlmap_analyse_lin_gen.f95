!!****m* ANALYSE/m_qtlmap_analyse_lin_gen
!!  NAME
!!    m_qtlmap_analyse_lin_gen
!!  SYNOPSIS
!!
!!  DESCRIPTION
!!
!!  NOTES
!!
!!  BUGS
!!
!!  HISTORY
!!
!!  SEE ALSO
!!
!!  COPYRIGHT
!!***
!! You can use this space for remarks that should not be included
!! in the documentation.
!!/
module m_qtlmap_analyse_lin_gen
    use m_qtlmap_base
    use m_qtlmap_math
    use m_qtlmap_analyse_gen
    use m_qtlmap_types
    use m_qtlmap_output_handler
    use m_qtlmap_log

    !TODO
    !a enlever lorsque ca ne sera plus connecte a ce module....
    use m_qtlmap_haplotype_ldla

    implicit none
    private
    save

!!****v* m_qtlmap_analyse_lin_gen/ntlevp
!!  NAME
!!   ntlevp
!!  DESCRIPTION
!!   number of level for qtl sire interaction with fixed effect
!!  NOTES
!!   contingence
!!***
    integer  ,                              public       :: ntlevp
!!****v* m_qtlmap_analyse_lin_gen/ntlevm
!!  NAME
!!   ntlevm
!!  DESCRIPTION
!!   number of level for qtl dam interaction with fixed effect
!!  NOTES
!!   contingence
!!***
    integer  ,                              public       :: ntlevm
!!****v* m_qtlmap_analyse_lin_gen/nbco
!!  NAME
!!   nbco
!!  DESCRIPTION
!!   number of covariate in the model for the current trait ic
!!  NOTES
!!   contingence
!!***
    integer  ,                              public       :: nbco
!!****v* m_qtlmap_analyse_lin_gen/nbfem
!!  NAME
!!   nbfem
!!  DESCRIPTION
!!   number of female estimable
!!  NOTES
!!   contingence
!!***
    integer  ,                              public       :: nbfem
!!****v* m_qtlmap_analyse_lin_gen/nbef
!!  NAME
!!   nbef
!!  DESCRIPTION
!!   number of fixed effect in the model for the current trait ic
!!  NOTES
!!   contingence
!!***
    integer  ,                              public       :: nbef
!!****v* m_qtlmap_analyse_lin_gen/ntniv
!!  NAME
!!   ntniv
!!  DESCRIPTION
!!   number of level defined in the contingence matrix
!!  NOTES
!!   contingence
!!***
    integer  ,                              public       :: ntniv
!!****v* m_qtlmap_analyse_lin_gen/nbniv
!!  NAME
!!   nbniv
!!  DESCRIPTION
!!   number of level for all fixed effect for the current trait ic
!!  NOTES
!!   contingence
!!***
    integer  ,                              public       :: nbniv
!!****v* m_qtlmap_analyse_lin_gen/ntnivmax
!!  NAME
!!   ntnivmax
!!  DESCRIPTION
!!   maximum number of level (maximum column size) in the contingence matrix
!!  NOTES
!!  contingence
!!***
    integer  ,                              public       :: ntnivmax
!!****v* m_qtlmap_analyse_lin_gen/nteffmax
!!  NAME
!!   nteffmax
!!  DESCRIPTION
!!   maximum number of effect defined in the contingence matrix
!!  NOTES
!!  contingence
!!***
    integer  ,                              public       :: nteffmax
!!****v* m_qtlmap_analyse_lin_gen/nbnivest
!!  NAME
!!   nbnivest
!!  DESCRIPTION
!!   number of level estimable (given by the cholesky decomposition and the seuil SEUIL_CHO )
!!  NOTES
!!   nbnivest <= ntniv
!!   contingence
!!***
    integer  ,                              public       :: nbnivest
!!****v* m_qtlmap_analyse_lin_gen/ntnifix
!!  NAME
!!   ntnifix
!!  DESCRIPTION
!!   index of the latest level of fixed effect (use to get first level/effect of covariate)
!!  NOTES
!!   contingence
!!***
    integer  ,                              public       :: ntnifix
!!****v* m_qtlmap_analyse_lin_gen/meff
!!  NAME
!!   meff
!!  DESCRIPTION
!!   index of fixed effect to remove in the construction of the contingence matrix
!!  NOTES
!!    meff<=0 => no effect, see test_lin
!!***
    integer  ,                              public       :: meff   ! modlin
!!****v* m_qtlmap_analyse_lin_gen/mcov
!!  NAME
!!   mcov
!!  DESCRIPTION
!!   index of covariate  to remove in the construction of the contingence matrix
!!  NOTES
!!    mcov<=0 => no effect, see test_lin
!!   contingence
!!***
    integer  ,                              public       :: mcov   ! modlin
!!****v* m_qtlmap_analyse_lin_gen/mint
!!  NAME
!!   mint
!!  DESCRIPTION
!!   index of interaction fixed effect-qtl  to remove in the construction of the contingence matrix
!!  NOTES
!!    mint<=0 => no effect, see test_lin
!!   contingence
!!***
    integer  ,                              public       :: mint   ! modlin

    !$omp threadprivate (ntlevp,ntlevm,nbco,nbfem,nbef,ntniv,nbniv,ntnivmax,nteffmax,nbnivest,ntnifix,meff,mcov,mint)

!!****v* m_qtlmap_analyse_lin_gen/corniv
!!  NAME
!!   corniv
!!  DESCRIPTION
!!    get the corresponding estimable index of a level inside a contingence matrix
!!  DIMENSIONS
!!   ntniv
!!  NOTES
!!  contingence
!!***
    integer      ,dimension(:),   allocatable,public        :: corniv
!!****v* m_qtlmap_analyse_lin_gen/nivdir
!!  NAME
!!   nivdir
!!  DESCRIPTION
!!    get the corresponding level of a given effect for a kd progeny
!!  DIMENSIONS
!!   nd,nteff
!!  NOTES
!!  JM:
!!  le tableau nivdir donne pour chaque descendant la position
!!  dans le vecteur des effets fixes, des effets exprimes par ce descendant
!!
!!***
    integer         ,dimension(:,:),   allocatable,public   :: nivdir
!!****v* m_qtlmap_analyse_lin_gen/vecsol
!!  NAME
!!   vecsol
!!  DESCRIPTION
!!    boolean vector of estimability of level referenced in the contingence matrix
!!  DIMENSIONS
!!   ntniv
!!  NOTES
!!  contingence
!!***
    logical       ,dimension(:),   allocatable,public       :: vecsol

!!****v* m_qtlmap_analyse_lin_gen/covdir
!!  NAME
!!   covdir
!!  DESCRIPTION
!!    get the corresponding value of level of a given covariate effect for a kd progeny
!!  DIMENSIONS
!!   nd,ncov
!!  NOTES
!!  le tableau covdir donne pour chaque descendant la valeur des  covariables
!!  retenues dans le modele
!!  contingence
!!***
    real (kind=dp)  ,dimension(:,:),  allocatable, public   :: covdir
!!****v* m_qtlmap_analyse_lin_gen/precis
!!  NAME
!!   precis
!!  DESCRIPTION
!!    precision of each level
!!  DIMENSIONS
!!   ntniv
!!  NOTES
!!   contingence
!!***
    real (kind=dp)       ,dimension(:),allocatable,public   :: precis
!!****v* m_qtlmap_analyse_lin_gen/xx
!!  NAME
!!   xx
!!  DESCRIPTION
!!    the incidence matrix X'.X , X contingence matrix
!!  DIMENSIONS
!!   ntniv,ntniv
!!  NOTES
!!   contingence
!!***
    real (kind=dp)   ,dimension(:,:), allocatable, public   :: xx   ! Matrice d'incidence
!!****v* m_qtlmap_analyse_lin_gen/xxx
!!  NAME
!!   xxx
!!  DESCRIPTION
!!    the inverse of the incidence matrix  (X'.X) -1 , X contingence matrix
!!  DIMENSIONS
!!   ntniv,ntniv
!!  NOTES
!!   contingence,confusion
!!***
    real (kind=dp)   ,dimension(:,:), allocatable, private  :: xxx  ! Inverse de la matrice d incidence


    !$omp threadprivate (corniv,nivdir,xx,xxx,vecsol,covdir,precis)


!!****v* m_qtlmap_analyse_lin_gen/prbp
!!  NAME
!!   prbp
!!  DESCRIPTION
!!    probabilities cumulates (according all dam genotype probabilities) to receive the qtl from sire
!!  DIMENSIONS
!!   nd
!!  NOTES
!!    preprinc,contingence
!!***
    real (kind=dp)  ,dimension(:),    allocatable, public     :: prbp
!!****v* m_qtlmap_analyse_lin_gen/prbm
!!  NAME
!!   prbm
!!  DESCRIPTION
!!    probabilities cumulates (according all dam genotype probabilities) to receive the qtl from dam
!!  DIMENSIONS
!!   nd
!!  NOTES
!!    preprinc,contingence
!!***
    real (kind=dp)  ,dimension(:),    allocatable, public     :: prbm
!!****v* m_qtlmap_analyse_lin_gen/pp_ldla
!!  NAME
!!   pp_ldla
!!  DESCRIPTION
!!    probabilities to receive the ith haplotype sire for the kkd progenies, kkd come from ngend
!!  DIMENSIONS
!!   maxval(ngend),2
!!  NOTES
!!    preprinc,contingence_ldla
!!***
    real (kind=dp)  ,dimension(:,:),    allocatable, public   :: pp_ldla
!!****v* m_qtlmap_analyse_lin_gen/pm_ldla
!!  NAME
!!   pm_ldla
!!  DESCRIPTION
!!    probabilities to receive the ith haplotype dam for the kkd progenies according the genotype, kkd come from ngend
!!  DIMENSIONS
!!   maxval(ngenom),maxval(ngend),2
!!  NOTES
!!    preprinc,contingence_ldla
!!***
    real (kind=dp)  ,dimension(:,:,:),    allocatable, public :: pm_ldla
!!****v* m_qtlmap_analyse_lin_gen/pb_haplo
!!  NAME
!!   pb_haplo
!!  DESCRIPTION
!!
!!  DIMENSIONS
!!   nd,NB_HAPLO_PRIOR
!!  NOTES
!!    preprinc,contingence_ldla
!!***
    real (kind=dp)  ,dimension(:,:),    allocatable, public   :: pb_haplo
!!****v* m_qtlmap_analyse_lin_gen/prop_haplo_info
!!  NAME
!!   prop_haplo_info
!!  DESCRIPTION
!!
!!  DIMENSIONS
!!
!!  NOTES
!!
!!***
    real (kind=dp)  ,dimension(:),    allocatable, public     :: prop_haplo_info
!!****v* m_qtlmap_analyse_lin_gen/ppt
!!  NAME
!!   ppt
!!  DESCRIPTION
!!
!!  DIMENSIONS
!!
!!  NOTES
!!
!!***
    real (kind=dp)  ,dimension(:),    allocatable, public     :: ppt     ! modlin
!!****v* m_qtlmap_analyse_lin_gen/pmt
!!  NAME
!!   pmt
!!  DESCRIPTION
!!
!!  DIMENSIONS
!!
!!  NOTES
!!
!!***
    real (kind=dp)  ,dimension(:),    allocatable, public     :: pmt     ! modlin

    !$omp threadprivate (prbp,prbm,pp_ldla,pm_ldla,pb_haplo,prop_haplo_info,ppt,pmt)

!!****v* m_qtlmap_analyse_lin_gen/par0
!!  NAME
!!   par0
!!  DESCRIPTION
!!   Solution finded under H0
!!  DIMENSIONS
!!   nbnivest
!!  NOTES
!!
!!***
    real (kind=dp)       ,dimension(:),allocatable,public   :: par0
!!****v* m_qtlmap_analyse_lin_gen/precis0
!!  NAME
!!   precis0
!!  DESCRIPTION
!!   Precision of each level under H0
!!  DIMENSIONS
!!   ntniv
!!  NOTES
!!
!!***
    real (kind=dp)       ,dimension(:),allocatable,public   :: precis0
!!****v* m_qtlmap_analyse_lin_gen/vecsol0
!!  NAME
!!   vecsol0
!!  DESCRIPTION
!!   boolean vector of estimability level under H0
!!  DIMENSIONS
!!   ntniv
!!  NOTES
!!
!!***
    logical              ,dimension(:),allocatable,public   :: vecsol0

    !$omp threadprivate (par0,precis0,vecsol0)

!!****v* m_qtlmap_analyse_lin_gen/par1
!!  NAME
!!   par1
!!  DESCRIPTION
!!   Solution finded under H1
!!  DIMENSIONS
!!   nbnivest
!!  NOTES
!!
!!***
    real (kind=dp)       ,dimension(:),allocatable,public   :: par1
!!****v* m_qtlmap_analyse_lin_gen/precis1
!!  NAME
!!   precis1
!!  DESCRIPTION
!!   Precision of each level under H1
!!  DIMENSIONS
!!   ntniv
!!  NOTES
!!
!!***
    real (kind=dp)       ,dimension(:),allocatable,public   :: precis1
!!****v* m_qtlmap_analyse_lin_gen/vecsol1
!!  NAME
!!   vecsol1
!!  DESCRIPTION
!!   boolean vector of estimability level under H1
!!  DIMENSIONS
!!   ntniv
!!  NOTES
!!
!!***
    logical              ,dimension(:),allocatable,public   :: vecsol1


!!****v* m_qtlmap_analyse_lin_gen/xinc
!!  NAME
!!   xinc
!!  DESCRIPTION
!!   contingence matrix
!!  DIMENSIONS
!!   nd,ntniv
!!  NOTES
!!
!!***
    real (kind=dp) , dimension(:,:),allocatable   ,public   :: xinc

    !$omp threadprivate (par1,precis1,vecsol1,xinc)

    public :: init_analyse_lin_gen
    public :: init_contingence
    public :: end_contingence
    public :: contingence
    public :: contingence_cox
    public :: contingence_ldla
    public :: set_filter_optim
    public :: confusion
    public :: prepinc
    public :: precision
    public :: end_analyse_lin_gen

    contains

!!****f* m_qtlmap_analyse_lin_gen/init_analyse_lin_gen
!!  NAME
!!    init_analyse_lin_gen
!!  DESCRIPTION
!!    Initialisation/allocation of solution array
!!  INPUTS
!!     ic         : index of trait
!!    nqtl        : the hypothesis tested
!!
!!  NOTES
!!  SOURCE
      subroutine init_analyse_lin_gen(dataset,spt,ic,nqtl)
          type(QTLMAP_DATASET)       ,intent(in) :: dataset
          type(PDD_BUILD)            ,intent(in) :: spt
          integer        ,intent(in)   :: ic
          integer        ,intent(in)   :: nqtl
          integer                      :: ntniv,nteff
          integer                      :: stat,ngd,npar,chr
          type(GENEALOGY_BASE) , pointer :: dg
          type(PHENOTYPE_BASE) , pointer :: dpa
          type(DATAMODEL_BASE) , pointer :: dpm

          dpm => dataset%phenoModel
          dg => dataset%genea
          dpa => dataset%phenoAnimal

          call set_ntnivmax(dataset,ic,nqtl,ntnivmax,nteffmax,ntlevp,nbniv)
          ntlevm=ntlevp
          npar = dg%np+ntnivmax

          if ( associated (dpm%nmod) ) then
           if ( dpm%nmod(ic)>1) then
              npar = npar+dpm%nmod(ic)-1
           end if
          end if

          allocate (par0(npar),STAT=stat)
          call check_allocate(stat,'par0 [m_qtlmap_analyse_lin_gen]')
          par0=0.d0
          allocate (precis0(ntnivmax),STAT=stat)
          call check_allocate(stat,'precis0 [m_qtlmap_analyse_lin_gen]')
          precis0=0.d0
          allocate (vecsol0(ntnivmax),STAT=stat)
          call check_allocate(stat,'vecsol0 [m_qtlmap_analyse_lin_gen]')

          allocate (par1(npar),STAT=stat)
          call check_allocate(stat,'par1 [m_qtlmap_analyse_lin_gen]')
          par1=0.d0
          allocate (precis1(ntnivmax),STAT=stat)
          call check_allocate(stat,'precis1 [m_qtlmap_analyse_lin_gen]')
          precis1=0.d0
          allocate (vecsol1(ntnivmax),STAT=stat)
          call check_allocate(stat,'vecsol1 [m_qtlmap_analyse_lin_gen]')

          !Pour garder la compatibilite, on instancie pour le thread courant les
          !structure de contingence
          call init_contingence(dataset,spt)

      end subroutine init_analyse_lin_gen
!!***

!!****f* m_qtlmap_analyse_lin_gen/init_contingence
!!  NAME
!!    init_contingence
!!  DESCRIPTION
!!    Initialisation/allocation of contingence arrays
!!
!!  NOTES
!!  SOURCE
      subroutine init_contingence(dataset,spt)
          type(QTLMAP_DATASET)       ,intent(in)   :: dataset
          type(PDD_BUILD)            ,intent(in)   :: spt

          integer   :: stat,ngd,npar,chr
          type(GENEALOGY_BASE) , pointer :: dg
          type(PHENOTYPE_BASE) , pointer :: dpa
          type(DATAMODEL_BASE) , pointer :: dpm

          dpm => dataset%phenoModel
          dg => dataset%genea
          dpa => dataset%phenoAnimal

          !print *,"** init_contingence ** "
          allocate (xx(ntnivmax,ntnivmax),STAT=stat)
	      call check_allocate(stat,'xx [init_contingence]')

	      allocate (vecsol(ntnivmax),STAT=stat)
	      call check_allocate(stat,'vecsol [init_contingence]')

          allocate (corniv(ntnivmax),STAT=stat)
          call check_allocate(stat,'corniv [init_contingence]')
          allocate (nivdir(dg%nd,nteffmax),STAT=stat)
          call check_allocate(stat,'nivdir [init_contingence]')

          allocate (xxx(ntnivmax,ntnivmax),STAT=stat)
          call check_allocate(stat,'xxx [init_contingence]')

          allocate (covdir(dg%nd,dpm%ncov),STAT=stat)
          call check_allocate(stat,'covdir [init_contingence]')

          !ngd = 4*nd
!          ngd=0
!          do chr=1,nchr
!           ngd = max(ngd,ngend(chr,ngenom(chr,nm+1)+1))
!          end do
          !31/08/2010 ! correction bug OFI: on cherche l indice maximum correspondant a ngend
          ngd = maxval(spt%ngend)

          allocate (pp_ldla(ngd,2),STAT=stat)
          call check_allocate(stat,'pp_ldla [init_contingence]')
          !31/08/2010 ! correction bug OFI: size(ngend) => maxval(ngenom) => donne le nombre de genotype possible toute femelles confondues
          allocate (pm_ldla(maxval(spt%ngenom),ngd,2),STAT=stat)
          call check_allocate(stat,'pm_ldla [init_contingence]')
          allocate (pb_haplo(dg%nd,dataset%params%NB_HAPLO_PRIOR),STAT=stat)
          call check_allocate(stat,'pb_haplo [init_contingence]')
          allocate (prop_haplo_info(dg%nd),STAT=stat)
          call check_allocate(stat,'prop_haplo_info [init_contingence]')

          allocate (prbp(dg%nd),STAT=stat)
          call check_allocate(stat,'prbp [init_contingence]')
          allocate (prbm(dg%nd),STAT=stat)
          call check_allocate(stat,'prbm [init_contingence]')
          allocate (ppt(ngd),STAT=stat)
          call check_allocate(stat,'ppt [init_contingence]')
          allocate (pmt(ngd),STAT=stat)
          call check_allocate(stat,'pmt [init_contingence]')
          allocate (precis(ntnivmax),STAT=stat)
          call check_allocate(stat,'precis [init_contingence]')
          allocate (xinc(dg%nd,ntnivmax),STAT=stat)
          call check_allocate(stat,'xinc [init_contingence]')
      end subroutine init_contingence
!!***


!!****f* m_qtlmap_analyse_lin_gen/end_contingence
!!  NAME
!!    end_contingence
!!  DESCRIPTION
!!    Deallocation of contingence arrays
!!
!!  NOTES
!!  SOURCE
      subroutine end_contingence
           !print *,"** end_contingence ** "
          deallocate (corniv)
          deallocate (nivdir)
          deallocate(xx)
          deallocate(xxx)
          deallocate(vecsol)
          deallocate (covdir)
          deallocate (prbp)
          deallocate (prbm)
          deallocate (ppt)
          deallocate (pmt)
          deallocate (pp_ldla)
          deallocate (pm_ldla)
          deallocate (pb_haplo)
          deallocate (prop_haplo_info)
          deallocate (precis)
          deallocate(xinc)
      end subroutine end_contingence
!!***


!!****f* m_qtlmap_analyse_lin_gen/end_analyse_lin_gen
!!  NAME
!!    end_analyse_lin_gen
!!  DESCRIPTION
!!    Deallocation of solution arrays
!!
!!  NOTES
!!  SOURCE
      subroutine end_analyse_lin_gen
          deallocate (par0)
          deallocate (precis0)
          deallocate (vecsol0)
          deallocate (par1)
          deallocate (precis1)
          deallocate (vecsol1)
          call end_contingence
      end subroutine end_analyse_lin_gen
!!***

!!****f* m_qtlmap_analyse_lin_gen/contingence
!!  NAME
!!    contingence
!!  DESCRIPTION
!!
!!  INPUTS
!!     ic  : index of the trait
!!    ihyp : under this hypothesis
!!   itest : testing nuisances effects
!! est_moy : estimation of the general mean
!!
!!  NOTES
!!    construction de la matrice de contingence et recherche des contraintes
!!  SOURCE
      subroutine contingence (dataset,spt,ic,ihyp,itest,est_moy)
          type(QTLMAP_DATASET)       ,intent(in)   :: dataset
          type(PDD_BUILD)            ,intent(in)   :: spt
!
!  Tableaux dimensionn�s au nombre d'effets � estimer
      integer       :: indestim(dataset%genea%nm)
!
!  divers
      integer i,j,k,iniv,jniv,ki,kj,kniv,kkd,kd,ief,jef,icov,jcov
      integer km,ip,nm1,nm2,nd1,nd2,jm,nbtef,nbtco,nbint,nbtint
      integer nkd,nteff,nivintp,nbt,nivintm,nd
      logical , intent(in) :: itest
      integer , intent(in) :: ic,ihyp
      logical , intent(in) , optional :: est_moy
      real (kind=dp) , dimension(:,:),allocatable   :: triang
      integer :: stat

      logical :: estmoylocal=.true.
      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa
      type(DATAMODEL_BASE) , pointer :: dpm

      dpm => dataset%phenoModel
      dg => dataset%genea
      dpa => dataset%phenoAnimal

      ! NOTE Olivier : j ai mis est_moy optional et j initlise une variable local
      ! modlin de jm n a pas besoin d etre modifie, c est donc moin intrusif....
      if ( present(est_moy) ) then
          estmoylocal=est_moy
      end if

!cccccccccccccccccccccccc
!cccccccccccccccccccccccc
! L'analyse est faite pour le ic �me caract�re
!
!  Les premier �l�ment du vecteur des effets � estimer sont la moyenne
!  g�n�rale, les effets p�re et les effets m�re
!
!
!  initilisation des compteurs
!
!  nbfem est le nombre de m�res "estimables" d'apr�s ndmin
!  nbef est le nombre d'effets fix�s parasites
!  nbco le nombre de covariables
!  nbintp le nombre d'effet hierarchis�s dnas les effets qtl p�re
!  nbintm idem pour les effets qtl m�re
!  nbniv est le nombre total de niveaux d'effets parasites
!  nbnivp, celui des efets hierarchis�s dans les effets qtl p�re
!  nbnivm, idem pour les qtl m�re
!
!
      nbfem=0
      do jm=1,dg%nm
        if(dpa%estime(ic,jm)) then
          nbfem=nbfem+1
          indestim(jm)=nbfem
        end if
      end do
      nbef=dpm%modele(ic,1)
      nbtef=nbef
      nbco=dpm%modele(ic,2)
      nbtco=nbco
      nbint=dpm%modele(ic,3)
      nbtint=nbint
!
!  si itest est � true, il s'agit de retirer un des effet du modele
!  et de recalculer la matrice X'X et les estimabilit�
!
       if(.not.itest) then
         meff=0
         mcov=0
         mint=0
       else
         if(meff.ne.0) nbef=nbef-1
         if(mcov.ne.0) nbco=nbco-1
         if(mint.ne.0) nbint=nbint-1
       end if

!
! construction de la matrice d' incidence  xinc
!   de dimension nb de descendants x 1+nb pere + nb mere +
!   nb niveaux effets fix�s + nb c   covariables  (sous H0)
!
!
!  le tableau nivdir donne pour chaque descendant la position
!  dans le vecteur des effets fix�s, des effets exprim�s par ce descendant
!
!  le tableau covdir donne pour chaque descendant la valeur des  covariables
!  retenues dans le mod�le
!
!  ntniv est le nombre total de niveaux d'effets et covariables
!
      nkd=0
      xinc=0.d0
      nivdir = 0
      do ip=1,dg%np
      nm1=dg%nmp(ip)+1
      nm2=dg%nmp(ip+1)
      do jm=nm1,nm2
      nd1=dg%ndm(jm)+1
      nd2=dg%ndm(jm+1)
      do kd=nd1,nd2
      if(dpa%presentc(ic,kd)) then
        nkd=nkd+1
!
!  le premier effet fix� est la moyenne g�n�rale
!
       if (estmoylocal) then
	  nivdir(kd,1)=1
          xinc(nkd,1)=1
          ntniv=1
          nteff=1
	else  ! pas d esitmation de la moyenne generale
          ntniv=0
          nteff=0
	endif

!
!  puis, sous H1, les effets QTL
!
        if(ihyp.eq.1) then
!
!  si des effets fix�s sont en interaction avec l'all�le paternel au qtl
!  il faut estimer les effets qtl pour chacun des niveaux concern�s
!
!  on commence donc par cr�er un effet composites regroupant l'ensemble de
!  ces effets
!
         nivintp=1
         ntlevp=1
         ! 07/04/2010
         ! correction bug ofi : nbt=3+nbef+nbco n est pas valide
         nbt=3+nbtef+nbtco
         if(nbint.ne.0) then
           ief=0
           do jef=1,nbtint
           if(mint.ne.jef)then
             ief=ief+1
             nivintp=nivintp+ntlevp*(dpa%nivx(kd,dpm%modele(ic,nbt+jef))-1)
             ntlevp=ntlevp*dpm%nlev(dpm%modele(ic,nbt+jef))
           end if
           end do
         end if

!
!  puis on cr�e le tableau des niveaux cumul�s et la matrice d'incidence
!
         nivdir(kd,nteff+1)=ntniv+nivintp+ntlevp*(ip-1)
         xinc(nkd,nivdir(kd,nteff+1))=prbp(kd)
         ntniv=ntniv+ntlevp*dg%np
         nteff=nteff+1
!
!  m�me op�ration pour le qtl transmis par la m�re
!
        if(nbfem.ne.0)then
         nivintm=1
         ntlevm=1
          ! 07/04/2010
         ! correction bug ofi : nbt=nbt+nbint n est pas valide
         !nbt=nbt+nbint
         if(nbint.ne.0) then
           ief=0
           do jef=1,nbtint
           if(nbint.ne.0) then
             ief=ief+1
             nivintm=nivintm+ntlevm*(dpa%nivx(kd,dpm%modele(ic,nbt+jef))-1)
             ntlevm=ntlevm*dpm%nlev(dpm%modele(ic,nbt+jef))
           end if
           end do
         end if

         if(dpa%estime(ic,jm)) then
           km=dpa%iam(ic,dg%repfem(jm))
           nivdir(kd,nteff+1)=ntniv+nivintm+ntlevm*(km-1)
           xinc(nkd,nivdir(kd,nteff+1))=prbm(kd)
         end if
         ntniv=ntniv+ntlevm*dpa%namest(ic)
         nteff=nteff+1
        end if
!
!  fin de H1
!
        end if
!
!  puis les effets polyg�niques parentaux
!
        nivdir(kd,nteff+1)= ntniv+ip
        xinc(nkd,nivdir(kd,nteff+1))=1
        ntniv=ntniv+dg%np
        nteff=nteff+1

        if(nbfem.ne.0)then
        if(dpa%estime(ic,jm)) then
          nivdir(kd,nteff+1)=ntniv+indestim(jm)
          xinc(nkd,nivdir(kd,nteff+1))=1
        end if
        ntniv=ntniv+nbfem
        nteff=nteff+1
        end if

!
!  autres effets fix�s
!
        nbniv=0
        ief=0
        do jef=1,nbtef
          if(meff.ne.jef) then
          ief=ief+1
          nivdir(kd,nteff+ief)=ntniv+nbniv+dpa%nivx(kd,dpm%modele(ic,3+jef))
          xinc(nkd,nivdir(kd,nteff+ief))=1
          nbniv=nbniv+dpm%nlev(dpm%modele(ic,3+jef))
          end if
        end do

        ntnifix=ntniv+nbniv
        ! 07/04/2010
        ! correction bug ofi : nbt=3+nbef n est pas valide
        nbt=3+nbtef

!
!  covariables
!
        icov=0
        do jcov=1,nbtco
          if(mcov.ne.jcov) then
          icov=icov+1
          covdir(kd,icov)=dpa%covar(kd,dpm%modele(ic,nbt+jcov))
          xinc(nkd,ntnifix+icov)=dpa%covar(kd,dpm%modele(ic,nbt+jcov))
          end if
        end do
        ntniv=ntnifix+nbco
        nbt=nbt+nbco

        end if


      end do
      end do
      end do

      allocate (triang(ntniv,ntniv),STAT=stat)
	  call check_allocate(stat,'triang [m_qtlmap_analyse_lin_gen]')

!
!  construction de X'X
!
      xx=0.d0

      do kd=1,nkd
        do iniv =1,ntniv
          do jniv=1,ntniv
            xx(iniv,jniv)=xx(iniv,jniv)+xinc(kd,iniv)*xinc(kd,jniv)
          end do
        end do
      end do

!  application de la d�composition de choleski pour d�terminer les contraintes
!
!  Dans le vecteur vecsol les effets � estimer sont indiqu�s � TRUE et ceux
!  qui ne sont pas estimables � FALSE
!
!  La matrice triang est la matrice triangulaire sup�rieure M telle que
!  M'M = X'X
!
      do i=1,ntniv
        do j=1,ntniv
	     triang(i,j)=0.d0
	    end do
      end do

      do j=1,ntniv
       vecsol(j)=.true.
       triang(j,j)=xx(j,j)
       do k=1,j-1
         triang(j,j)=triang(j,j)- triang(j,k)*triang(j,k)
       end do

       if(triang(j,j).gt.dataset%params%SEUIL_CHO) then
         triang(j,j)=dsqrt(triang(j,j))
         do i=j+1,ntniv
           triang(i,j)=xx(i,j)
           do k=1,j-1
             triang(i,j)=triang(i,j)-triang(i,k)*triang(j,k)
           end do
           triang(i,j)=triang(i,j)/triang(j,j)
         end do

       else
         vecsol(j)=.false.
       end if
      end do

!
!  Table de correspondance entre niveaux de d�part et niveaux estimables
!  Au lieu de ntniv effets � estimer, on en a plus que nbnivest
!  corniv(i) est la position, dans la liste des effets estimables
!  du i �me effet initial
!

      nbnivest=0
      do ief=1,ntniv
        if(vecsol(ief))then
          nbnivest=nbnivest+1
          corniv(ief)=nbnivest
        end if
      end do

      return
      end subroutine contingence
!!***

!!****f* m_qtlmap_analyse_lin_gen/contingence_cox
!!  NAME
!!    contingence_cox
!!  DESCRIPTION
!!
!!  INPUTS
!!     ic  : index of the trait
!!    ihyp : under this hypothesis
!!   itest : testing nuisances effects
!! est_moy : estimation of the general mean
!!
!!  NOTES
!!    construction de la matrice de contingence et recherche des contraintes
!!  SOURCE
      subroutine contingence_cox (dataset,spt,ic,ihyp,itest,est_moy)
        type(QTLMAP_DATASET)       ,intent(in)         :: dataset
        type(PDD_BUILD)            ,intent(in)         :: spt
!
!  Tableaux dimensionn�s au nombre d'effets � estimer
      integer       :: indestim(dataset%genea%nm)
!
!  divers
      integer i,j,k,iniv,jniv,ki,kj,kniv,kkd,kd,ief,jef,icov,jcov
      integer km,ip,nm1,nm2,nd1,nd2,jm,nbtef,nbtco,nbint,nbtint
      integer nkd,nteff,nivintp,nbt,nivintm,nd
      logical , intent(in) :: itest
      integer , intent(in) :: ic,ihyp
      logical , intent(in) , optional :: est_moy
      real (kind=dp) , dimension(:,:),allocatable   :: triang
      integer :: stat

      logical :: estmoylocal=.true.
      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa
      type(DATAMODEL_BASE) , pointer :: dpm

      dpm => dataset%phenoModel
      dg => dataset%genea
      dpa => dataset%phenoAnimal


      ! NOTE Olivier : j ai mis est_moy optional et j initlise une variable local
      ! modlin de jm n a pas besoin d etre modifie, c est donc moin intrusif....
      if ( present(est_moy) ) then
          estmoylocal=est_moy
      end if

!cccccccccccccccccccccccc
!cccccccccccccccccccccccc
! L'analyse est faite pour le ic �me caract�re
!
!  Les premier �l�ment du vecteur des effets � estimer sont la moyenne
!  g�n�rale, les effets p�re et les effets m�re
!
!
!  initilisation des compteurs
!
!  nbfem est le nombre de m�res "estimables" d'apr�s ndmin
!  nbef est le nombre d'effets fix�s parasites
!  nbco le nombre de covariables
!  nbintp le nombre d'effet hierarchis�s dnas les effets qtl p�re
!  nbintm idem pour les effets qtl m�re
!  nbniv est le nombre total de niveaux d'effets parasites
!  nbnivp, celui des efets hierarchis�s dans les effets qtl p�re
!  nbnivm, idem pour les qtl m�re
!
!
      nbfem=0
      do jm=1,dg%nm
        if(dpa%estime(ic,jm)) then
          nbfem=nbfem+1
          indestim(jm)=nbfem
        end if
      end do
      nbef=dpm%modele(ic,1)
      nbtef=nbef
      nbco=dpm%modele(ic,2)
      nbtco=nbco
      nbint=dpm%modele(ic,3)
      nbtint=nbint
!
!  si itest est � true, il s'agit de retirer un des effet du modele
!  et de recalculer la matrice X'X et les estimabilit�
!
       if(.not.itest) then
         meff=0
         mcov=0
         mint=0
       else
         if(meff.ne.0) nbef=nbef-1
         if(mcov.ne.0) nbco=nbco-1
         if(mint.ne.0) nbint=nbint-1
       end if

!
! construction de la matrice d' incidence  xinc
!   de dimension nb de descendants x 1+nb pere + nb mere +
!   nb niveaux effets fix�s + nb c   covariables  (sous H0)
!
!
!  le tableau nivdir donne pour chaque descendant la position
!  dans le vecteur des effets fix�s, des effets exprim�s par ce descendant
!
!  le tableau covdir donne pour chaque descendant la valeur des  covariables
!  retenues dans le mod�le
!
!  ntniv est le nombre total de niveaux d'effets et covariables
!
      nkd=0
      xinc=0.d0
      nivdir = 0
      do ip=1,dg%np
      nm1=dg%nmp(ip)+1
      nm2=dg%nmp(ip+1)
      do jm=nm1,nm2
      nd1=dg%ndm(jm)+1
      nd2=dg%ndm(jm+1)
      do kd=nd1,nd2
      if(dpa%presentc(ic,kd)) then
        nkd=nkd+1
!
!  le premier effet fix� est la moyenne g�n�rale
!
       if (estmoylocal) then
	  nivdir(kd,1)=1
          xinc(nkd,1)=1
          ntniv=1
          nteff=1
	else  ! pas d esitmation de la moyenne generale
          ntniv=0
          nteff=0
	endif

!
!  puis, sous H1, les effets QTL
!
        if(ihyp.eq.1) then
!
!  si des effets fix�s sont en interaction avec l'all�le paternel au qtl
!  il faut estimer les effets qtl pour chacun des niveaux concern�s
!
!  on commence donc par cr�er un effet composites regroupant l'ensemble de
!  ces effets
!
         nivintp=1
         ntlevp=1
         ! 07/04/2010
         ! correction bug ofi : nbt=3+nbef+nbco n est pas valide
         nbt=3+nbtef+nbtco
         if(nbint.ne.0) then
           ief=0
           do jef=1,nbtint
           if(mint.ne.jef)then
             ief=ief+1
             nivintp=nivintp+ntlevp*(dpa%nivx(kd,dpm%modele(ic,nbt+jef))-1)
             ntlevp=ntlevp*dpm%nlev(dpm%modele(ic,nbt+jef))
           end if
           end do
         end if

!
!  puis on cr�e le tableau des niveaux cumul�s et la matrice d'incidence
!
         nivdir(kd,nteff+1)=ntniv+nivintp+ntlevp*(ip-1)
         xinc(nkd,nivdir(kd,nteff+1))=prbp(kd)
         ntniv=ntniv+ntlevp*dg%np
         nteff=nteff+1
!
!  m�me op�ration pour le qtl transmis par la m�re
!
        if(nbfem.ne.0)then
         nivintm=1
         ntlevm=1
          ! 07/04/2010
         ! correction bug ofi : nbt=nbt+nbint n est pas valide
         !nbt=nbt+nbint
         if(nbint.ne.0) then
           ief=0
           do jef=1,nbtint
           if(nbint.ne.0) then
             ief=ief+1
             nivintm=nivintm+ntlevm*(dpa%nivx(kd,dpm%modele(ic,nbt+jef))-1)
             ntlevm=ntlevm*dpm%nlev(dpm%modele(ic,nbt+jef))
           end if
           end do
         end if

         if(dpa%estime(ic,jm)) then
           km=dpa%iam(ic,dg%repfem(jm))
           nivdir(kd,nteff+1)=ntniv+nivintm+ntlevm*(km-1)
           xinc(nkd,nivdir(kd,nteff+1))=prbm(kd)
         end if
         ntniv=ntniv+ntlevm*dpa%namest(ic)
         nteff=nteff+1
        end if
!
!  fin de H1
!
        end if
!
!  puis les effets polyg�niques parentaux
!
        if(dg%np.GE.2) then
          nivdir(kd,nteff+1)= ntniv+ip
          xinc(nkd,nivdir(kd,nteff+1))=1
          ntniv=ntniv+dg%np
          nteff=nteff+1
        endif
        if(nbfem.GE.2)then
          if(dpa%estime(ic,jm)) then
            nivdir(kd,nteff+1)=ntniv+indestim(jm)
            xinc(nkd,nivdir(kd,nteff+1))=1
          end if
        ntniv=ntniv+nbfem
        nteff=nteff+1
        end if

!
!  autres effets fix�s
!
        nbniv=0
        ief=0
        do jef=1,nbtef
          if(meff.ne.jef) then
          ief=ief+1
          nivdir(kd,nteff+ief)=ntniv+nbniv+dpa%nivx(kd,dpm%modele(ic,3+jef))
          xinc(nkd,nivdir(kd,nteff+ief))=1
          nbniv=nbniv+dpm%nlev(dpm%modele(ic,3+jef))
          end if
        end do

        ntnifix=ntniv+nbniv
        ! 07/04/2010
        ! correction bug ofi : nbt=3+nbef n est pas valide
        nbt=3+nbtef

!
!  covariables
!
        icov=0
        do jcov=1,nbtco
          if(mcov.ne.jcov) then
          icov=icov+1
          covdir(kd,icov)=dpa%covar(kd,dpm%modele(ic,nbt+jcov))
          xinc(nkd,ntnifix+icov)=dpa%covar(kd,dpm%modele(ic,nbt+jcov))
          end if
        end do
        ntniv=ntnifix+nbco
        nbt=nbt+nbco

        end if


      end do
      end do
      end do

      allocate (triang(ntniv,ntniv),STAT=stat)
	  call check_allocate(stat,'triang [m_qtlmap_analyse_lin_gen]')

!
!  construction de X'X
!
      xx=0.d0

      do kd=1,nkd
        do iniv =1,ntniv
          do jniv=1,ntniv
            xx(iniv,jniv)=xx(iniv,jniv)+xinc(kd,iniv)*xinc(kd,jniv)
          end do
        end do
      end do

!  application de la d�composition de choleski pour d�terminer les contraintes
!
!  Dans le vecteur vecsol les effets � estimer sont indiqu�s � TRUE et ceux
!  qui ne sont pas estimables � FALSE
!
!  La matrice triang est la matrice triangulaire sup�rieure M telle que
!  M'M = X'X
!
      do i=1,ntniv
        do j=1,ntniv
	     triang(i,j)=0.d0
	    end do
      end do

      do j=1,ntniv
       vecsol(j)=.true.
       triang(j,j)=xx(j,j)
       do k=1,j-1
         triang(j,j)=triang(j,j)- triang(j,k)*triang(j,k)
       end do

       if(triang(j,j).gt.dataset%params%SEUIL_CHO) then
         triang(j,j)=dsqrt(triang(j,j))
         do i=j+1,ntniv
           triang(i,j)=xx(i,j)
           do k=1,j-1
             triang(i,j)=triang(i,j)-triang(i,k)*triang(j,k)
           end do
           triang(i,j)=triang(i,j)/triang(j,j)
         end do

       else
         vecsol(j)=.false.
       end if
      end do

!
!  Table de correspondance entre niveaux de d�part et niveaux estimables
!  Au lieu de ntniv effets � estimer, on en a plus que nbnivest
!  corniv(i) est la position, dans la liste des effets estimables
!  du i �me effet initial
!

      nbnivest=0
      do ief=1,ntniv
        if(vecsol(ief))then
          nbnivest=nbnivest+1
          corniv(ief)=nbnivest
        end if
      end do

      return
      end subroutine contingence_cox
!!***


!!****f* m_qtlmap_analyse_lin_gen/contingence_ldla
!!  NAME
!!    contingence_ldla
!!  DESCRIPTION
!!
!!  INPUTS
!!     ic      : index of the trait
!!    ihyp     : under this hypothesis
!!   itest     : testing nuisances effects
!! est_moy     : estimation of the general mean
!!    chr      : chromosome index
!!    n        : position tested
!! details     :
!! est_moy     :
!! option_anal : type analysis  * 'LA  ', 'LD  ', 'LDLA', 'LDJH' *
!! hsire       :
!! hdam        :
!!
!!  OUTPUTS
!! var_yQy     :
!!
!!  NOTES
!!    construction de la matrice de contingence et recherche des contraintes
!!  SOURCE
      subroutine contingence_ldla (dataset,spt,shp,ic,ihyp,itest,chr,&
          n,var_yQy,details,est_moy,option_anal,hsire,hdam)

      type(QTLMAP_DATASET)       ,intent(in)            :: dataset
      type(PDD_BUILD)            ,intent(in)            :: spt
      type(HAPLOTYPE_POSITION_BUILD)   ,intent(inout)   :: shp
!
!  Tableaux dimensionn�s au nombre d'effets � estimer
      integer       :: indestim(dataset%genea%nm)
!
!  divers
      integer i,j,k,iniv,jniv,ki,kj,kniv,kkd,kd,ief,jef,icov,jcov,chr
      integer km,ip,nm1,nm2,nd1,nd2,jm,nbtef,nbtco,nbint,nbtint
      integer nkd,nteff,nivintp,nbt,nivintm
      integer irank
      integer i_gam
      logical , intent(in) :: itest
      integer , intent(in) :: ic,ihyp,n
      logical , intent(in) :: details
      logical , intent(in) :: est_moy
      character(len=4) , intent(in)           :: option_anal
      real (kind=dp) var_yQy,ups
      real (kind=dp) :: xinc_red(dataset%genea%nd,ntnivmax)
      real (kind=dp) :: temp_q(dataset%genea%nd,ntnivmax),temp_x(dataset%genea%nd,ntnivmax)
      real (kind=dp) :: xx_red(ntnivmax,ntnivmax),temp_xx(ntnivmax,ntnivmax)
      real (kind=dp) :: mat_q(dataset%genea%nd,dataset%genea%nd),su
      real (kind=dp) , dimension(:,:),allocatable  :: triang
      integer :: stat
      integer i_haplo,geno,j_gam,iall,jp
      logical :: estmoylocal,is_la,is_ld,is_ldla,is_ldjh
      logical :: hsire,hdam
      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa
      type(DATAMODEL_BASE) , pointer :: dpm

      dpm => dataset%phenoModel
      dg => dataset%genea
      dpa => dataset%phenoAnimal

      is_la = (option_anal == 'LA  ')
      is_ld = (option_anal == 'LD  ')
      is_ldla = (option_anal == 'LDLA')
      is_ldjh = (option_anal == 'LDJH')


         estmoylocal=est_moy

!cccccccccccccccccccccccc
!cccccccccccccccccccccccc
! L'analyse est faite pour le ic �me caract�re
!
!  Les premier �l�ment du vecteur des effets � estimer sont la moyenne
!  g�n�rale, les effets p�re et les effets m�re
!
!
!  initilisation des compteurs
!
!  nbfem est le nombre de m�res "estimables" d'apr�s ndmin
!  nbef est le nombre d'effets fix�s parasites
!  nbco le nombre de covariables
!  nbintp le nombre d'effet hierarchis�s dnas les effets qtl p�re
!  nbintm idem pour les effets qtl m�re
!  nbniv est le nombre total de niveaux d'effets parasites
!  nbnivp, celui des efets hierarchis�s dans les effets qtl p�re
!  nbnivm, idem pour les qtl m�re
!
!
    ! print*,'entree dans contingence ldla'
      nbfem=0
      do jm=1,dg%nm
        if(dpa%estime(ic,jm)) then
          nbfem=nbfem+1         ! correspond a namest
          indestim(jm)=nbfem    ! correspond a iam .....
        end if
      end do
      nbef=dpm%modele(ic,1)
      nbtef=nbef
      nbco=dpm%modele(ic,2)
      nbtco=nbco
      nbint=dpm%modele(ic,3)
      nbtint=nbint
!
!  si itest est � true, il s'agit de retirer un des effet du modele
!  et de recalculer la matrice X'X et les estimabilit�
!
       if(.not.itest) then
         meff=0
         mcov=0
         mint=0
       else
         if(meff.ne.0) nbef=nbef-1
         if(mcov.ne.0) nbco=nbco-1
         if(mint.ne.0) nbint=nbint-1
       end if

    !  allocate(xinc(nd,size(xx,1)))

!
! construction de la matrice d' incidence  xinc
!   de dimension nb de descendants x 1+nb pere + nb mere +
!   nb niveaux effets fix�s + nb c   covariables  (sous H0)
!
!
!  le tableau nivdir donne pour chaque descendant la position
!  dans le vecteur des effets fix�s, des effets exprim�s par ce descendant
!
!  le tableau covdir donne pour chaque descendant la valeur des  covariables
!  retenues dans le mod�le
!
!  ntniv est le nombre total de niveaux d'effets et covariables
!
      nkd=0
      xinc=0.d0
      nivdir = 0
      do ip=1,dg%np
      nm1=dg%nmp(ip)+1
      nm2=dg%nmp(ip+1)
      do jm=nm1,nm2
      nd1=dg%ndm(jm)+1
      nd2=dg%ndm(jm+1)
      do kd=nd1,nd2
      if(dpa%presentc(ic,kd) ) then
        nkd=nkd+1
!
!  le premier effet fix� est la moyenne g�n�rale
!
       if (estmoylocal) then
	  nivdir(kd,1)=1
          xinc(nkd,1)=1
          ntniv=1
          nteff=1
	else  ! pas d estimation de la moyenne generale
          ntniv=0
          nteff=0
	endif

!
!  puis, sous H1, les effets QTL
!
        if(ihyp.eq.1) then

        if(is_ld .or. is_ldla) then

           do i_haplo=1,shp%nb_haplo_reduit
             pb_haplo(kd,i_haplo)=0.d0
           end do
!
! haplotype du pere
           if(hsire)then
           do j=1,2
             su = dble(sum(shp%pb_haplo_reduit(shp%num_haplo_pere(ip,j,1:shp%nb_gam_pere(ip,j)))))
             if ( su == 0 ) cycle
             do i_gam=1,shp%nb_gam_pere(ip,j)
               pb_haplo(kd,shp%num_haplo_pere(ip,j,i_gam))= pb_haplo(kd,shp%num_haplo_pere(ip,j,i_gam)) &
                                                       +pp_ldla(kd,j)*(shp%pb_haplo_reduit(shp%num_haplo_pere(ip,j,i_gam)) / su )
         !      print *,'igam:',i_gam,' pb_haplo_reduit:',pb_haplo_reduit(num_haplo_pere(ip,j,i_gam))&
         !     ,' pp_ldla:',pp_ldla(kd,j)
             end do !i_gam
           end do !j
           ! 07/04/2010
           ! correction bug ofi : nbt=3+nbef+nbco n est pas valide
           end if ! opt_model
! haplotype de la mere si phase connue   
           if(hdam)then
           if(dpa%estime(ic,jm)) then
             do geno=spt%ngenom(chr,jm)+1,spt%ngenom(chr,jm+1)
               do j=1,2
                 su = dble(sum(shp%pb_haplo_reduit(shp%num_haplo_mere(geno,j,1:shp%nb_gam_mere(geno,j)))))
                 if ( su == 0 ) cycle
                 do i_gam=1,shp%nb_gam_mere(geno,j)
                   pb_haplo(kd,shp%num_haplo_mere(geno,j,i_gam))  = pb_haplo(kd,shp%num_haplo_mere(geno,j,i_gam)) &
                                                 +pm_ldla(geno,kd,j)*(shp%pb_haplo_reduit(shp%num_haplo_mere(geno,j,i_gam)) / su )
               !    print *,'hdam igam:',i_gam,' pb_haplo_reduit:',pb_haplo_reduit(num_haplo_mere(geno,j,i_gam))&
             ! ,' pp_ldla:',pm_ldla(geno,kd,j)
                 end do !i_gam
               end do !j 
            end do !geno
! haplotype de la mere si phase inconnue   
           else
 !            do j_gam=1,nb_gam(kd)
 !              su = dble(sum(pb_haplo_reduit(num_haplo_desc(kd,j_gam,1:nb_gam_desc(kd,j_gam)))))
 !              if ( su == 0 ) cycle
 !              do i_gam=1,nb_gam_desc(kd,j_gam)
 !                  pb_haplo(kd,num_haplo_desc(kd,j_gam,i_gam))  = pb_haplo(kd,num_haplo_desc(kd,j_gam,i_gam)) &
  !                                               +prob_gam(kd,j_gam) * (pb_haplo_reduit(num_haplo_desc(kd,j_gam,i_gam)) /su)
   !            end do !i_gam
   !          end do !j_gam 

!            do i_gam=1,nb_gam_desc(kd)
!             pb_haplo(kd,num_haplo_desc(kd,i_gam))= pb_haplo(kd,num_haplo_desc(kd,i_gam))+&  
!               pb_haplo_desc(kd,i_gam)
!            end do  !i_gam


           do i_gam=1,shp%nb_haplo_reduit
             pb_haplo(kd,i_gam)= pb_haplo(kd,i_gam)+shp%pb_haplo_desc(kd,i_gam)
           end do  !i_gam

           end if !estime
           end if !opt_model (hdam)
!
! ligne d'incidence
!   
           do i_haplo=1,shp%nb_haplo_reduit
             nivdir(kd,nteff+i_haplo)=ntniv+i_haplo
             xinc(nkd,nivdir(kd,nteff+i_haplo))=pb_haplo(kd,i_haplo)
          !   print *,name_haplo_reduit(i_haplo)
           end do
           ntniv=ntniv+shp%nb_haplo_reduit
           nteff=nteff+shp%nb_haplo_reduit

         end if ! option_anal == LD , LDLA

         if(option_anal /= 'LD  ') then

!
!  si des effets fix�s sont en interaction avec l'all�le paternel au qtl
!  il faut estimer les effets qtl pour chacun des niveaux concern�s
!
!  on commence donc par cr�er un effet composites regroupant l'ensemble de
!  ces effets
!
         nivintp=1
         ntlevp=1
         ! 07/04/2010
         ! correction bug ofi : nbt=3+nbef+nbco n est pas valide
         nbt=3+nbtef+nbtco
         if(nbint.ne.0) then
           ief=0
           do jef=1,nbtint
           if(mint.ne.jef)then
             ief=ief+1
             nivintp=nivintp+ntlevp*(dpa%nivx(kd,dpm%modele(ic,nbt+jef))-1)
             ntlevp=ntlevp*dpm%nlev(dpm%modele(ic,nbt+jef))
           end if
           end do
         end if

!
!  puis on cr�e le tableau des niveaux cumul�s et la matrice d'incidence
!
         nivdir(kd,nteff+1)=ntniv+nivintp+ntlevp*(ip-1)
         xinc(nkd,nivdir(kd,nteff+1))=prbp(kd)
         ntniv=ntniv+ntlevp*dg%np
         nteff=nteff+1

!
!  m�me op�ration pour le qtl transmis par la m�re
!
        if(option_anal == 'LDJH') then
          if(shp%shared_haplo(kd)/=0) then
            jp=int((1+shp%shared_haplo(kd))/2)
            iall=shp%shared_haplo(kd)-2*(jp-1)
            nivdir(kd,nteff+1)=ntniv+jp
            xinc(nkd,nivdir(kd,nteff+1))=xinc(nkd,nivdir(kd,nteff+1))+(-1)**iall
          end if
          ntniv=ntniv+dg%np
          nteff=nteff+1


        else
          if(nbfem.ne.0)then
            nivintm=1
            ntlevm=1
            ! 07/04/2010
         ! correction bug ofi : nbt=nbt+nbint n est pas valide
            !nbt=nbt+nbint
            if(nbint.ne.0) then
              ief=0
              do jef=1,nbtint
                if(nbint.ne.0) then
                  ief=ief+1
                  nivintm=nivintm+ntlevm*(dpa%nivx(kd,dpm%modele(ic,nbt+jef))-1)
                  ntlevm=ntlevm*dpm%nlev(dpm%modele(ic,nbt+jef))
                end if
              end do
            end if ! nbint ne 0
            if(dpa%estime(ic,jm)) then
              km=dpa%iam(ic,dg%repfem(jm))
              nivdir(kd,nteff+1)=ntniv+nivintm+ntlevm*(km-1)
              xinc(nkd,nivdir(kd,nteff+1))=prbm(kd)
            end if
            ntniv=ntniv+ntlevm*dpa%namest(ic)
            nteff=nteff+1
          end if ! nbfem ne 0
        end if ! option LDJH
!
!  fin de l optionLA , LDLA ou LDJH
       end if
!
!  fin de H1
!
        end if
!
!  puis les effets polyg�niques parentaux
!
        nivdir(kd,nteff+1)= ntniv+ip
        xinc(nkd,nivdir(kd,nteff+1))=1
        ntniv=ntniv+dg%np
        nteff=nteff+1

        if(nbfem.ne.0)then
        if(dpa%estime(ic,jm)) then
          nivdir(kd,nteff+1)=ntniv+indestim(jm)
          xinc(nkd,nivdir(kd,nteff+1))=1
        end if
        ntniv=ntniv+nbfem
        nteff=nteff+1
        end if

!
!  autres effets fix�s
!
        nbniv=0
        ief=0
        do jef=1,nbtef
          if(meff.ne.jef) then
          ief=ief+1
          nivdir(kd,nteff+ief)=ntniv+nbniv+dpa%nivx(kd,dpm%modele(ic,3+jef))
          xinc(nkd,nivdir(kd,nteff+ief))=1
          nbniv=nbniv+dpm%nlev(dpm%modele(ic,3+jef))
          end if
        end do

        ntnifix=ntniv+nbniv
        ! 07/04/2010
        ! correction bug ofi : nbt=3+nbef n est pas valide
        nbt=3+nbtef

!
!  covariables
!
        icov=0
        do jcov=1,nbtco
          if(mcov.ne.jcov) then
          icov=icov+1
          covdir(kd,icov)=dpa%covar(kd,dpm%modele(ic,nbt+jcov))
          xinc(nkd,ntnifix+icov)=dpa%covar(kd,dpm%modele(ic,nbt+jcov))
          end if
        end do
        ntniv=ntnifix+nbco
        nbt=nbt+nbco

        end if

      end do
      end do
      end do

      allocate (triang(ntniv,ntniv),STAT=stat)
	  call check_allocate(stat,'triangle [m_qtlmap_analyse_lin_gen]')

!
!  construction de X'X
!
      xx=0.d0

      do kd=1,nkd
        do iniv =1,ntniv
          do jniv=1,ntniv
            xx(iniv,jniv)=xx(iniv,jniv)+xinc(kd,iniv)*xinc(kd,jniv)
          end do
        end do
      end do

!  application de la d�composition de choleski pour d�terminer les contraintes
!
!  Dans le vecteur vecsol les effets � estimer sont indiqu�s � TRUE et ceux
!  qui ne sont pas estimables � FALSE
!
!  La matrice triang est la matrice triangulaire sup�rieure M telle que
!  M'M = X'X
!
      do i=1,ntniv
        do j=1,ntniv
	     triang(i,j)=0.d0
	    end do
      end do

      do j=1,ntniv
       vecsol(j)=.true.
       triang(j,j)=xx(j,j)
       do k=1,j-1
         triang(j,j)=triang(j,j)- triang(j,k)*triang(j,k)
       end do

       if(triang(j,j).gt.dataset%params%SEUIL_CHO) then
         triang(j,j)=dsqrt(triang(j,j))
         do i=j+1,ntniv
           triang(i,j)=xx(i,j)
           do k=1,j-1
             triang(i,j)=triang(i,j)-triang(i,k)*triang(j,k)
           end do
           triang(i,j)=triang(i,j)/triang(j,j)
         end do

       else
         vecsol(j)=.false.
       end if
      end do
!
!  Table de correspondance entre niveaux de d�part et niveaux estimables
!  Au lieu de ntniv effets � estimer, on en a plus que nbnivest
!  corniv(i) est la position, dans la liste des effets estimables
!  du i �me effet initial
!

      nbnivest=0
      corniv=0
      do ief=1,ntniv
        if(vecsol(ief))then
          nbnivest=nbnivest+1
          corniv(ief)=nbnivest
        end if
      end do
!
      if(details) then

      deallocate (triang)

!
! compactage de X'X et xinc
!
      ki=0
      do iniv =1,ntniv
        if(vecsol(iniv)) then
          ki=ki+1
          do kd=1,nkd
            xinc_red(kd,ki)=xinc(kd,iniv)
          end do ! kd
          kj=0
          do jniv=1,ntniv
            if(vecsol(jniv)) then
              kj=kj+1
              xx_red(ki,kj)=xx(iniv,jniv)
            end if
          end do ! jniv
        end if
      end do ! iniv
!
!  calcul de la variance de la forme quadratique
!

      temp_x=xinc_red(:nkd,:nbnivest)
      temp_xx=xx_red(:nbnivest,:nbnivest)

!
!  inversion de la matrice d'incidence r�duite
!
      ups=1.d-15
      call ginv1(temp_xx,nbnivest,size(temp_xx,1),ups,irank)

      temp_q=matmul(temp_x,temp_xx)
      mat_q(:nbnivest,:nbnivest)=matmul(temp_q,transpose(temp_x))

      var_yQy=0.d0
      do iniv=1,nbnivest
        var_yQy=var_yQy+1.d0-mat_q(iniv,iniv)
      end do

      end if ! details

 !     deallocate(xinc)

      return
      end subroutine contingence_ldla
!!***

!!****f* m_qtlmap_analyse_lin_gen/set_filter_optim
!!  NAME
!!    set_filter_optim
!!  DESCRIPTION
!!    provide a boolean vector that informed which level are expressed according the half sib family
!!
!!  INPUTS
!!     ic             : index of the trait
!!    hetero          : true => heteroscedastic otherwise homoscedastic
!!   variance_haplo   : true => setting the evaluation of a residual variance with the evaluation all haplotypes effects
!!   ntnivmax         : maximum number of level
!!   ntniv            : number of level of the contingence matrix
!!   vecsol           : the boolean vector of the estimability of each effect
!!   xinc             : contingence matrix
!!
!!  OUTPUTS
!!    filter_inc      : the boolean result array
!!
!!  NOTES
!!    see m_qtlmap_optimization for further details
!!  SOURCE
    subroutine set_filter_optim(dataset,ic,hetero,variance_haplo,ntnivmax,ntniv,vecsol,xinc,filter_inc)
       type(QTLMAP_DATASET)                ,intent(in)          :: dataset
       integer                             ,intent(in)          :: ic,ntnivmax,ntniv
       logical                             ,intent(in)          :: hetero,variance_haplo
       logical  , dimension(ntnivmax)      ,intent(in)          :: vecsol
       real(kind=dp) , dimension(dataset%genea%nd,ntnivmax),intent(in)        :: xinc
       logical         , dimension(:,:,:)    ,intent(inout)     :: filter_inc

       integer :: kd1,kd2,ip,jm,ii,i,nstart
       type(GENEALOGY_BASE) , pointer :: dg
       type(PHENOTYPE_BASE) , pointer :: dpa

      dg => dataset%genea
      dpa => dataset%phenoAnimal

!Construction du filtre d optimisation
      filter_inc=.true.
      if (hetero) then
        filter_inc(:,:,1:dg%np)=.false.
        nstart = dg%np
      else
        filter_inc(:,:,1)=.true.
        nstart=1
      end if

      kd1=0
      kd2=0

      do ip=1,dg%np
         if (hetero) then
           filter_inc(ip,dg%nmp(ip)+1:dg%nmp(ip+1),ip)=.true.
         end if

         do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
            kd1=kd2+1
            kd2=kd2+count(dpa%presentc(ic,dg%ndm(jm)+1:dg%ndm(jm+1)))
            if (kd1<=kd2) then
               ii=0
               do i=1,ntniv
                if (vecsol(i)) then
                  ii=ii+1
                  filter_inc(ip,jm,nstart+ii)=any(xinc(kd1:kd2,i)/=0.d0,dim=1)
                  end if
               end do
             else
                filter_inc(ip,jm,:)=.false.
             end if
           end do
        end do




      end subroutine set_filter_optim
!!***

!!****f* m_qtlmap_analyse_lin_gen/confusion
!!  NAME
!!    confusion
!!  DESCRIPTION
!!    Print in the result file, the confusion (real values) between qtl effet and the other effect
!!
!!  INPUTS
!!     mod            : 'avant','apres'
!!     ic             : index of the trait
!!     chr            : index of the chromosome
!!    est_moy         : true => the general mean is included as an effect
!!   option_anal      : type analysis  * 'LA  ', 'LD  ', 'LDLA', 'LDJH' *
!!
!!  NOTES
!!    Mesure de la confuqion entre les effets QTL et d'autres effets
!!    THRES_CONFUSION
!!  SOURCE
     subroutine confusion(dataset,spt,mod,chr,ic,est_moy,option_anal)
      type(QTLMAP_DATASET)       ,intent(in)         :: dataset
      type(PDD_BUILD)            ,intent(in)         :: spt

      integer          , intent(in)           :: ic,chr
      character(len=*) , intent(in)           :: mod
      logical          , intent(in),optional  :: est_moy
      character(len=4) , intent(in)           :: option_anal
      real (kind=dp) :: xcor(ntnivmax,ntnivmax)
      logical , dimension(dataset%genea%nfem) :: femIsOk
     ! real (kind=dp)  ,dimension(ntniv,ntniv)             :: xxx
!
! Divers
      logical :: icas,itest
      integer :: nfirst,ief,i,j,jef,ip,ilevp,iclic,iqtl,ntot
      integer :: jniv,ilev,iniveau,ico,jm,ilevm,ifem
      real (kind=dp) ::corrmin
      logical        :: estMoyLocal = .true.
      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa
      type(DATAMODEL_BASE) , pointer :: dpm

      dpm => dataset%phenoModel
      dg => dataset%genea
      dpa => dataset%phenoAnimal

      if ( present(est_moy) ) then
           estMoyLocal = est_moy
      end if
!
      if(mod.eq.'avant') then
        itest=.false.
        call prepinc(dataset,spt,chr,1,ic,option_anal)
        call contingence(dataset,spt,ic,1,itest,estMoyLocal)
        write(nficout,601)
  601   format(//,80('*')/'Test of confusion between QTL ',  &
      'and other effects in the initial full model',/,       &
      '(test based on the correlation between columns of the ','incidence matrix)',//)
      end if

      if(mod.eq.'apres') then
        write(nficout,602)
  602   format(//,80('*')/'Test of confusion between QTL ',   &
       'and other effects in the final constained model',/,   &
       '(test based on the correlation between columns of the ','incidence matrix)',//)
      end if
!  icas est un indicateur de defaut possible d'identifiabilit�
!
      icas=.false.
!

      if (estMoyLocal) then
         nfirst= 1+ntlevp*dg%np+ntlevm*dpa%namest(ic)+1
      else
         nfirst=   ntlevp*dg%np+ntlevm*dpa%namest(ic)+1
      end if

      corrmin = 0.d0

      write(nficout,600)
  600 format(//,80('*')/'Confusion between QTL and other effects ','(final constained model)',//)

!
!  xcor(i,j), la correlation entre les estim�es i et j d'apr�s X'X-1
!
      xcor=0.d0
      ief=0
      do i=1,ntniv
         if(vecsol(i))then
            ief=ief+1
            jef=0
            do j=1,ntniv
              if(vecsol(j))then
                jef=jef+1
                if (xxx(ief,ief) == 0 .or. xxx(jef,jef) == 0) then
                  xcor(i,j)=0.d0
                else
                xcor(i,j)=xxx(ief,jef)/dsqrt(xxx(ief,ief)*xxx(jef,jef))
                end if
               end if
            end do
         end if
      end do

!
!  on caclule les correlations entre les colonnes de X'X des effets qtl p�re et de
!  tous les autres effets
!
      do ip=1,dg%np
!        print *,"ip:",ip
        do ilevp=1,ntlevp
           iclic=0
           if (estMoyLocal) then
            iqtl=1+ilevp+ntlevp*(ip-1)
          else
            iqtl=ilevp+ntlevp*(ip-1)
          end if

          if (vecsol(iqtl))then
!
!  test d'une confusion avec les effets polyg�niques p�re
!
          ntot=nfirst-1
          !print *,"iqtl:",iqtl
          do ief=1,dg%np
            jniv=ntot+ief
           ! print *,'jniv:',jniv
           ! print *,'jniv:',jniv,dabs(xcor(iqtl,jniv))
             if(dabs(xcor(iqtl,jniv)).ge.corrmin)corrmin=dabs(xcor(iqtl,jniv))
!
!  si la correlation est sup�rieure � corMax, on imprime une alerte
!
            if(dabs(xcor(iqtl,jniv)).ge.dataset%params%THRES_CONFUSION)then
              icas=.true.
              if(iclic.eq.0) then
                write(nficout,605)trim(dg%pere(ip)),ilevp
  605 format('Risk for sire ',a,' of confusion between the QTL level',i3,' and :')
                iclic=1
              end if
              write(nficout,610)trim(dg%pere(ief)),xcor(iqtl,jniv)
  610 format('the sire ',a,' polygenic effect ',' (correlation ',f5.2,')')
             end if
            end do

!
!  test d'une confusion avec les effets polyg�niques m�re
!
          ntot=ntot+dg%np

          do ief=1,nbfem
            jniv=ntot+ief
!
!  si la correlation est sup�rieur � seuil, on imprime une alerte
!
          if(dabs(xcor(iqtl,jniv)).ge.corrmin)corrmin=dabs(xcor(iqtl,jniv))
              if(dabs(xcor(iqtl,jniv)).ge.dataset%params%THRES_CONFUSION)then
              icas=.true.
               if(iclic.eq.0) then
                 write(nficout,605)trim(dg%pere(ip)),ilevp
                 iclic=1
               end if
               write(nficout,611)trim(dg%mere(ief)),xcor(iqtl,jniv)
  611 format('the dam ',a,' polygenic effect ',' (correlation ',f5.2,')')
            end if
            end do
!
!  test d'une confusion avec les effets fix�s
!

          ntot=ntot+nbfem
          if(nbniv.ne.0)then
           ilev=0
           do ief=1,nbef
             do iniveau=1,dpm%nlev(dpm%modele(ic,3+ief))
             ilev=ilev+1
             jniv=ntot+ilev
!
!  si la correlation est sup�rieur � corMax, on imprime une alerte
!
           if(dabs(xcor(iqtl,jniv)).ge.corrmin)corrmin=dabs(xcor(iqtl,jniv))
              if(dabs(xcor(iqtl,jniv)).ge.dataset%params%THRES_CONFUSION)then
              icas=.true.
               if(iclic.eq.0) then
                 write(nficout,605)trim(dg%pere(ip)),ilevp
                 iclic=1
               end if

               write(nficout,612)iniveau,trim(dpm%namefix(dpm%modele(ic,3+ief))),xcor(iqtl,jniv)
  612 format('the level ',i3,' of the effect ',a15,' (correlation ',f5.2,')')
             end if
             end do
            end do
          end if
!
!  test d'une confusion avec les covariables
!
      ntot=ntot+nbniv
      if(nbco.ne.0)then
        do ico=1,nbco
          jniv=ntot+ico
!
!  si la correlation est sup�rieur � corMax, on imprime une alerte
!
          if(dabs(xcor(iqtl,jniv)).ge.corrmin)corrmin=dabs(xcor(iqtl,jniv))
              if(dabs(xcor(iqtl,jniv)).ge.dataset%params%THRES_CONFUSION)then
               icas=.true.
               if(iclic.eq.0) then
                 write(nficout,605)trim(dg%pere(ip)),ilevp
                 iclic=1
               end if
               write(nficout,613)trim(dpm%namecov(dpm%modele(ic,3+nbef+ico))),xcor(iqtl,jniv)
  613 format('the covariable ',a15,' (correlation ',f5.2,')')
              end if
        end do
!
      end if
      end if
!
!  fin pour le niveau du qtl du pere
!
        end do
      end do

!
!  meme scenario pour les m�res
!
      femisok=.false.
      do jm=1,dg%nm
        if (.not. dpa%estime(ic,jm)) cycle
        ifem=dpa%iam(ic,dg%repfem(jm))
        if (femIsok(ifem)) cycle
        femisok(ifem)=.true.
      !  print *,"jm:",jm
        do ilevm=1,ntlevm
          iclic=0

          if (estMoyLocal) then
            iqtl=1+ntlevp*dg%np+ilevm+ntlevm*(ifem-1)
          else
            iqtl=ntlevp*dg%np+ilevm+ntlevm*(ifem-1)
          end if
       !   print *,"iqtl:",iqtl
          if (vecsol(iqtl))then
!
!  test d'une confusion avec les effets polyg�niques p�re
!
          ntot=nfirst-1
          do ief=1,dg%np
            jniv=ntot+ief
!
!  si la correlation est sup�rieure � corMax, on imprime une alerte
!
 !           print *,'jniv:',jniv
  !          print *,dabs(xcor(iqtl,jniv))
            if(dabs(xcor(iqtl,jniv)).ge.dataset%params%THRES_CONFUSION)then
              icas=.true.
              if(iclic.eq.0) then
                write(nficout,615)trim(dg%mere(jm)),ilevm
  615 format('Risk for dam ',a,' of confusion between the QTL level',i3,' and :')
                iclic=1
              end if
              write(nficout,610)trim(dg%pere(ief)),xcor(iqtl,jniv)
             end if
            end do
!
!  test d'une confusion avec les effets polyg�niques m�re
!
          ntot=ntot+dg%np
          do ief=1,nbfem
            jniv=ntot+ief
!
!  si la correlation est sup�rieur � corMax, on imprime une alerte
!
          if(dabs(xcor(iqtl,jniv)).ge.corrmin)corrmin=dabs(xcor(iqtl,jniv))
              if(dabs(xcor(iqtl,jniv)).ge.dataset%params%THRES_CONFUSION)then
               icas=.true.
               if(iclic.eq.0) then
                write(nficout,615)trim(dg%mere(jm)),ilevm
                 iclic=1
               end if
               write(nficout,611)trim(dg%mere(ief)),xcor(iqtl,jniv)
             end if
            end do
!
!  test d'une confusion avec les effets fix�s
!

          ntot=ntot+nbfem
          if(nbniv.ne.0)then
           ilev=0
           do ief=1,nbef
             do iniveau=1,dpm%nlev(dpm%modele(ic,3+ief))
             ilev=ilev+1
             jniv=ntot+ilev
!
!  si la correlation est sup�rieur � corMax, on imprime une alerte
!
           if(dabs(xcor(iqtl,jniv)).ge.corrmin)corrmin=dabs(xcor(iqtl,jniv))
              if(dabs(xcor(iqtl,jniv)).ge.dataset%params%THRES_CONFUSION)then
               icas=.true.
               if(iclic.eq.0) then
                write(nficout,615)trim(dg%mere(jm)),ilevm
                 iclic=1
               end if

             write(nficout,612)iniveau,trim(dpm%namefix(dpm%modele(ic,3+ief))),xcor(iqtl,jniv)
               end if
             end do
            end do
          end if
!
!  test d'une confusion avec les covariables
!
      ntot=ntot+nbniv
      if(nbco.ne.0)then
        do ico=1,nbco
          jniv=ntot+ico
!
!  si la correlation est sup�rieur � corMax, on imprime une alerte
!
          if(dabs(xcor(iqtl,jniv)).ge.corrmin)corrmin=dabs(xcor(iqtl,jniv))
              if(dabs(xcor(iqtl,jniv)).ge.dataset%params%THRES_CONFUSION)then
               icas=.true.
               if(iclic.eq.0) then
                write(nficout,615)trim(dg%mere(jm)),ilevm
                 iclic=1
               end if
               write(nficout,613)trim(dpm%namecov(dpm%modele(ic,3+nbef+ico))),xcor(iqtl,jniv)
              end if
        end do
!
      end if
      end if
!

!  fin pour le niveau du qtl de la mere
!
        end do
      end do
!
!  situation sans probleme
!
      if(.not.icas) then
        write(nficout,616)corrmin
  616   format(/,' No confusion detected',/, &
      ' the highest correlation is : ', F7.3,/,80('*')/)
      end if

      return
      end subroutine confusion
!!***

!!****f* m_qtlmap_analyse_lin_gen/prepinc
!!  NAME
!!    prepinc
!!  DESCRIPTION
!!    Prepare the following arrays for the build of the contingence matrix :
!!     prbp,prbm,pp_ldla,pm_ldla,pmt,ppt
!!  INPUTS
!!     chr            : index of the chromosome
!!     n              : the tested position
!!   bcar_icar        : the index of trait
!!   option_anal      : type analysis  * 'LA  ', 'LD  ', 'LDLA', 'LDJH' *
!!
!!  NOTES
!!  SOURCE
   subroutine prepinc(dataset,spt,chr,n,bcar_icar,option_anal)
      type(QTLMAP_DATASET)       ,intent(in)        :: dataset
      type(PDD_BUILD)            ,intent(in)        :: spt
      integer         , intent(in)                  :: chr,n
      integer         , intent(in)                  :: bcar_icar
      character(len=4)                  ,intent(in) :: option_anal

      integer :: ip,nm1,nm2,jm,nd1,nd2,kd,ngeno1,ngeno2,ig,kkd
      logical :: is_diff_la=.false.
      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa
      type(DATAMODEL_BASE) , pointer :: dpm

      dpm => dataset%phenoModel
      dg => dataset%genea
      dpa => dataset%phenoAnimal

      is_diff_la = (option_anal /= "LA  ")
!******************************************************************************
!  preparation de la matrice d'incidence
!
! ************* deb chgt *************
      prbp=0.d0
      prbm=0.d0
      pp_ldla=0.d0
      pm_ldla=0.d0

! ************* fin chgt *************
      do ip=1,dg%np
        nm1=dg%nmp(ip)+1
        nm2=dg%nmp(ip+1)
        do jm=nm1,nm2
          nd1=dg%ndm(jm)+1
          nd2=dg%ndm(jm+1)
           do kd=nd1,nd2
              prbp(kd)=0.d0
              prbm(kd)=0.d0
           end do
          ngeno1=spt%ngenom(chr,jm)+1
          ngeno2=spt%ngenom(chr,jm+1)
          do ig=ngeno1,ngeno2
            nd1=spt%ngend(chr,ig)+1
            nd2=spt%ngend(chr,ig+1)
            do kd=nd1,nd2
              kkd=spt%ndesc(chr,kd)
              if(dpa%presentc(bcar_icar,kkd)) then
	       if( (.not. dpa%estime(bcar_icar,jm)) .or. ( dataset%params%opt_sib.eq.OPT_SIB_HS) )then
                ppt(kd)=-spt%pdd(chr,kd,1,n)+spt%pdd(chr,kd,3,n)
! ************* deb chgt *************
                if(is_diff_la) then
                  pp_ldla(kkd,1)=spt%pdd(chr,kd,1,n)
                  pp_ldla(kkd,2)=spt%pdd(chr,kd,3,n)
                end if
! ************* fin chgt *************
                prbp(kkd)=ppt(kd)
	       else
                ppt(kd)=-spt%pdd(chr,kd,1,n)-spt%pdd(chr,kd,2,n)+spt%pdd(chr,kd,3,n)+spt%pdd(chr,kd,4,n)
                prbp(kkd)=prbp(kkd)+spt%probg(chr,ig)*ppt(kd)
                pmt(kd)=-spt%pdd(chr,kd,1,n)+spt%pdd(chr,kd,2,n)-spt%pdd(chr,kd,3,n)+spt%pdd(chr,kd,4,n)
                prbm(kkd)=prbm(kkd)+spt%probg(chr,ig)*pmt(kd)
! ************* deb chgt *************
                if(is_diff_la) then
                 pp_ldla(kkd,1)=pp_ldla(kkd,1)+spt%probg(chr,ig)*(spt%pdd(chr,kd,1,n)+spt%pdd(chr,kd,2,n))
                 pp_ldla(kkd,2)=pp_ldla(kkd,2)+spt%probg(chr,ig)*(spt%pdd(chr,kd,3,n)+spt%pdd(chr,kd,4,n))
                 pm_ldla(ig,kd,1)=spt%probg(chr,ig)*(spt%pdd(chr,kd,1,n)+spt%pdd(chr,kd,3,n))
                 pm_ldla(ig,kd,2)=spt%probg(chr,ig)*(spt%pdd(chr,kd,2,n)+spt%pdd(chr,kd,4,n))
                end if
! ************* fin chgt *************
	       end if
              end if
            end do
          end do
        end do
      end do

      return
      end subroutine prepinc
!!***

!!****f* m_qtlmap_analyse_lin_gen/precision
!!  NAME
!!    precision
!!  DESCRIPTION
!!
!!  INPUTS
!!     xx             : incidence matrix
!!  OUTPUTS
!!     precis         : vector of precision values
!!
!!  NOTES
!!  SOURCE
      subroutine precision(xx,precis)
      implicit none

      real (kind=dp)  ,dimension(ntnivmax,ntnivmax)       , intent(in)   :: xx
      real (kind=dp)  ,dimension(ntnivmax), intent(out)  :: precis

      integer                                             :: iniv,kniv,jniv,lniv,irank
      real (kind=dp)                                      :: ups
     ! real (kind=dp)  ,dimension(ntniv,ntniv)             :: xxx
      call log_mess("START precision",DEBUG_DEF)
!
!  R�duction de la matrice d'incidence
!
     call log_mess("size(xx)="//trim(str(size(xx))),DEBUG_DEF)
     call log_mess("size(precis)="//trim(str(size(precis))),DEBUG_DEF)
     call log_mess("size(vecsol)="//trim(str(size(vecsol))),DEBUG_DEF)
     call log_mess("NTNIV="//trim(str(ntniv)),DEBUG_DEF)
     call log_mess("size(xxx)="//trim(str(size(xxx))),DEBUG_DEF)
 
     xxx=0.d0
     lniv=0
     kniv=0
     
     do iniv =1,ntniv
         if(vecsol(iniv)) then
           kniv=kniv+1
           lniv=0
           do jniv=1,ntniv
             if(vecsol(jniv)) then
               lniv=lniv+1
               xxx(kniv,lniv)=xx(iniv,jniv)
             end if
           end do
          end if
       end do
!
!  inversion de la matrice d'incidence r�duite pour calcul des precisions
!

      ups=1.d-15
      call ginv1(xxx,lniv,size(xxx,1),ups,irank)

      jniv=0
      do iniv=1,ntniv
        precis(iniv)=9999.d0
        if(vecsol(iniv)) then
          jniv=jniv+1
          precis(iniv)=xxx(jniv,jniv)
        end if
      end do

      !stop
      call log_mess("END precision",DEBUG_DEF)
      return
      end subroutine precision
!!***


end module m_qtlmap_analyse_lin_gen
