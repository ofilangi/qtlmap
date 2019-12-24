!!****m* ANALYSE/m_qtlmap_analyse_multitrait_DA
!!  NAME
!!    m_qtlmap_analyse_multitrait_DA
!!  DESCRIPTION
!!    Module analyse multicaractere DA
!!    init_analyse_DA_1QTL,perform_DA_1QTL, opti_DA_0qtl, opti_DA_1qtl, end_analyse_DA_1QTL
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
module m_qtlmap_analyse_multitrait_DA
  use m_qtlmap_types
  use m_qtlmap_log
  use m_qtlmap_optimization
  use m_qtlmap_math
  use m_qtlmap_analyse_gen, only : estmum,eff,cary,somy,effdf,somydf,carydf,xmu0p,xmu0m,sig0

  implicit none
  save

  type(GENEALOGY_BASE) , pointer :: p_dg
  type(PDD_BUILD)      , pointer :: p_spt
  type(PHENOTYPE_BASE) , pointer :: p_dpa
  type(DATAMODEL_BASE) , pointer :: p_dpm


!!****v* m_qtlmap_analyse_multitrait_DA/current_chr
!!  NAME
!!   current_chr
!!  DESCRIPTION
!!   The current chromosome while the likelihood calculs
!!***
  integer                                       , private :: current_ch ! current chromosome
!!****v* m_qtlmap_analyse_multitrait_DA/std
!!  NAME
!!   std
!!  DESCRIPTION
!!   the standart deviation (residual) found under H1
!! DIMENSIONS
!!   np
!!
!!***
  real (kind=dp)       ,dimension(:),allocatable,public   :: std
!!****v* m_qtlmap_analyse_multitrait_DA/xmoyp
!!  NAME
!!   xmoyp
!!  DESCRIPTION
!!   The polygenic sire effect found under H1
!! DIMENSIONS
!!   np
!!
!!***
  real (kind=dp)       ,dimension(:),allocatable,public   :: xmoyp
!!****v* m_qtlmap_analyse_multitrait_DA/ap
!!  NAME
!!   ap
!!  DESCRIPTION
!!   The qtl sire effect found der H1
!! DIMENSIONS
!!   np
!!
!!***
  real (kind=dp)       ,dimension(:),allocatable,public   :: ap
!!****v* m_qtlmap_analyse_multitrait_DA/am
!!  NAME
!!   am
!!  DESCRIPTION
!!   The qtl dam effect found der H1
!! DIMENSIONS
!!   nm
!!
!!***
  real (kind=dp)       ,dimension(:),allocatable,public   :: am
!!****v* m_qtlmap_analyse_multitrait_DA/xmoym
!!  NAME
!!   xmoym
!!  DESCRIPTION
!!   The polygenic dam effect found under H1
!! DIMENSIONS
!!   nm
!!
!!***
  real (kind=dp)       ,dimension(:),allocatable,public   :: xmoym

!!****v* m_qtlmap_analyse_multitrait_DA/sig1
!!  NAME
!!   sig1
!!  DESCRIPTION
!!   the standart deviation (residual) found under H0
!! DIMENSIONS
!!   np
!!
!!***
  real (kind=dp)       ,dimension(:),allocatable,private   :: sig1
!!****v* m_qtlmap_analyse_multitrait_DA/xmu1p
!!  NAME
!!   xmu1p
!!  DESCRIPTION
!!   The polygenic sire effect found under H0
!! DIMENSIONS
!!   np
!!
!!***
  real (kind=dp)       ,dimension(:),allocatable,private   :: xmu1p
!!****v* m_qtlmap_analyse_multitrait_DA/xmu1m
!! NAME
!!   xmu1m
!! DESCRIPTION
!!   The polygenic dam effect found under H0
!! DIMENSIONS
!!   nm
!!
!!***
  real (kind=dp)       ,dimension(:),allocatable,private   :: xmu1m
!!****v* m_qtlmap_analyse_multitrait_DA/f0
!!  NAME
!!   f0
!!  DESCRIPTION
!!   The maximum likelihood found under H0
!!***
  real (kind=dp)                                ,private    :: f0
!!****v* m_qtlmap_analyse_multitrait_DA/fp0
!!  NAME
!!   fp0
!!  DESCRIPTION
!!   The likelihood by half sib family under H0
!! DIMENSIONS
!!   np
!!
!!***
  real (kind=dp)       ,dimension(:),allocatable,private   :: fp0
!!****v* m_qtlmap_analyse_multitrait_DA/fp1
!!  NAME
!!   fp1
!!  DESCRIPTION
!!   The likelihood by half sib family under H1
!! DIMENSIONS
!!   np
!!
!!***
  real (kind=dp)       ,dimension(:),allocatable,private   :: fp1
!!****v* m_qtlmap_analyse_multitrait_DA/fm0
!!  NAME
!!   fm0
!!  DESCRIPTION
!!   The likelihood by full sib family under H0
!! DIMENSIONS
!!   nm
!!
!!***
  real (kind=dp)       ,dimension(:),allocatable,private   :: fm0
!!****v* m_qtlmap_analyse_multitrait_DA/fm1
!!  NAME
!!   fm1
!!  DESCRIPTION
!!   The likelihood by full sib family under H1
!! DIMENSIONS
!!   nm
!!
!!***
  real (kind=dp)       ,dimension(:),allocatable,private   :: fm1
!!****v* m_qtlmap_analyse_multitrait_DA/sompp
!!  NAME
!!   sompp
!!  DESCRIPTION
!!   initialization : loop before the minimization of the likelihood : opti_DA_1qtl
!!
!! DIMENSIONS
!!   maxval(ngenom)
!!***
  real (kind=dp)       ,dimension(:),allocatable,private   :: sompp
!!****v* m_qtlmap_analyse_multitrait_DA/sompm
!!  NAME
!!   sompm
!!  DESCRIPTION
!!   initialization : loop before the minimization of the likelihood : opti_DA_1qtl
!! DIMENSIONS
!!   maxval(ngenom)
!!***
  real (kind=dp)       ,dimension(:),allocatable,private   :: sompm
!!****v* m_qtlmap_analyse_multitrait_DA/sompmy
!!  NAME
!!   sompmy
!!  DESCRIPTION
!!   initialization : loop before the minimization of the likelihood : opti_DA_1qtl
!! DIMENSIONS
!!   maxval(ngenom)
!!***
  real (kind=dp)       ,dimension(:),allocatable,private   :: sompmy
!!****v* m_qtlmap_analyse_multitrait_DA/carpp
!!  NAME
!!   carpp
!!  DESCRIPTION
!!   Sum of square of probabilities to receive the qtl dam (for dam with enough dndmin)
!!   initialization : loop before the minimization of the likelihood : opti_DA_1qtl
!! DIMENSIONS
!!   maxval(ngenom)
!!***
  real (kind=dp)       ,dimension(:),allocatable,private   :: carpp
!!****v* m_qtlmap_analyse_multitrait_DA/carpm
!!  NAME
!!   carpm
!!  DESCRIPTION
!!   Sum of square of probabilities to receive the qtl dam (for dam with enough dndmin)
!!   initialization : loop before the minimization of the likelihood : opti_DA_1qtl
!! DIMENSIONS
!!   maxval(ngenom)
!!***
  real (kind=dp)       ,dimension(:),allocatable,private   :: carpm
!!****v* m_qtlmap_analyse_multitrait_DA/somppy
!!  NAME
!!   somppy
!!  DESCRIPTION
!!   Sum of probabilities to receive the qtl dam mult with a trait
!!   initialization : loop before the minimization of the likelihood : opti_DA_1qtl
!! DIMENSIONS
!!   maxval(ngenom)
!!***
  real (kind=dp)       ,dimension(:),allocatable,private   :: somppy
!!****v* m_qtlmap_analyse_multitrait_DA/somppm
!!  NAME
!!   somppm
!!  DESCRIPTION
!!    Sum of probabilities to receive the qtl dam
!!    initialization : loop before the minimization of the likelihood : opti_DA_1qtl
!! DIMENSIONS
!!    maxval(ngenom)
!!***
  real (kind=dp)       ,dimension(:),allocatable,private   :: somppm
!!****v* m_qtlmap_analyse_multitrait_DA/somppdf
!!  NAME
!!   somppdf
!!  DESCRIPTION
!!   Sum of probabilities to receive the qtl sire
!!   initialization : loop before the minimization of the likelihood : opti_DA_1qtl
!! DIMENSIONS
!!
!!***
  real (kind=dp)       ,dimension(:),allocatable,private   :: somppdf
!!****v* m_qtlmap_analyse_multitrait_DA/carppdf
!!  NAME
!!   carppdf
!!  DESCRIPTION
!!   Sum of square of probabilities to receive the qtl sire
!!   initialization : loop before the minimization of the likelihood : opti_DA_1qtl
!! DIMENSIONS
!!  np
!!***
  real (kind=dp)       ,dimension(:),allocatable,private   :: carppdf
!!****v* m_qtlmap_analyse_multitrait_DA/somppydf
!!  NAME
!!   somppydf
!!  DESCRIPTION
!!   Sum of probabilities to receive the qtl sire mult with a trait
!!   initialization : loop before the minimization of the likelihood : opti_DA_1qtl
!! DIMENSIONS
!!  np
!!***
  real (kind=dp)       ,dimension(:),allocatable,private   :: somppydf

!!****v* m_qtlmap_analyse_multitrait_DA/fp_1
!!  NAME
!!   fp_1
!!  DESCRIPTION
!!   The maximum likelihood by half sib family
!! DIMENSIONS
!!   np
!!***
  real (kind=dp)       ,dimension(:),allocatable,private   :: fp_1
!!****v* m_qtlmap_analyse_multitrait_DA/fm_1
!!  NAME
!!   fm_1
!!  DESCRIPTION
!!   The maximum likelihood by full sib family
!! DIMENSIONS
!!   nm
!!***
  real (kind=dp)       ,dimension(:),allocatable,private   :: fm_1
!!****v* m_qtlmap_analyse_multitrait_DA/f_1
!!  NAME
!!   f_1
!!  DESCRIPTION
!!   The maximum likelihood under H1
!!***
  real (kind=dp)      ,                            private   :: f_1



  !---------------------------------------------------------------
  ! INTERFACE
  !---------------------------------------------------------------
  public :: init_analyse_DA_1QTL
  public :: perform_DA_1QTL
  public :: opti_DA_0qtl
  public :: opti_DA_1qtl
  public :: end_analyse_DA_1QTL
  public :: set_solution_hypothesis1

contains

!!****f* m_qtlmap_analyse_multitrait_DA/init_analyse_DA_1QTL
!!  NAME
!!    init_analyse_DA_1QTL
!!  DESCRIPTION
!!    allocation of buffer/solution arrays
!!  SOURCE
  subroutine init_analyse_DA_1QTL(dataset,spt)
    type(QTLMAP_DATASET)       ,intent(in) :: dataset
    type(PDD_BUILD)            ,intent(in) :: spt

    integer           :: ic,jm,stat,ng
    type(GENEALOGY_BASE) , pointer :: dg
    type(PHENOTYPE_BASE) , pointer :: dpa
    type(DATAMODEL_BASE) , pointer :: dpm

    dpm => dataset%phenoModel
    dg => dataset%genea
    dpa => dataset%phenoAnimal

    ! Get Log Debug Information about estim

    do jm=1,size(dpa%estime,2)
     if( count(dpa%estime(:,jm))== dpm%ncar )  then
         call log_mess('DAM ['//trim(dg%mere(jm))//"] estime ** OK **",DEBUG_DEF)
     else
         call log_mess('DAM ['//trim(dg%mere(jm))//"] estime -- KO --",DEBUG_DEF)
         do ic=1,dpm%ncar
          if ( .not. dpa%estime(ic,jm) ) then
            call log_mess('  --> Trait ['//trim(dpm%carac(ic))//'] is not estim ',DEBUG_DEF)
          end if
         end do
     end if
   end do

    ng = maxval(spt%ngenom)
    allocate (sig1(dg%np),STAT=stat)
    call check_allocate(stat,'sig1 [m_qtlmap_analyse_multitrait_DA]')
    allocate (xmu1p(dg%np),STAT=stat)
    call check_allocate(stat,'xmu1p [m_qtlmap_analyse_multitrait_DA]')
    allocate (xmu1m(dg%nm),STAT=stat)
    call check_allocate(stat,'xmu1m [m_qtlmap_analyse_multitrait_DA]')
    allocate (fp0(dg%np),STAT=stat)
    call check_allocate(stat,'fp0 [m_qtlmap_analyse_multitrait_DA]')
    allocate (fp1(dg%np),STAT=stat)
    call check_allocate(stat,'fp1 [m_qtlmap_analyse_multitrait_DA]')
    allocate (fm0(dg%nm),STAT=stat)
    call check_allocate(stat,'fm0 [m_qtlmap_analyse_multitrait_DA]')
    allocate (fm1(dg%nm),STAT=stat)
    call check_allocate(stat,'fm1 [m_qtlmap_analyse_multitrait_DA]')
    allocate (sompp(ng),STAT=stat)
    call check_allocate(stat,'sompp [m_qtlmap_analyse_multitrait_DA]')
    allocate (sompm(ng),STAT=stat)
    call check_allocate(stat,'sompm [m_qtlmap_analyse_multitrait_DA]')
    allocate (sompmy(ng),STAT=stat)
    call check_allocate(stat,'sompmy [m_qtlmap_analyse_multitrait_DA]')
    allocate (carpp(ng),STAT=stat)
    call check_allocate(stat,'carpp [m_qtlmap_analyse_multitrait_DA]')
    allocate (carpm(ng),STAT=stat)
    call check_allocate(stat,'carpm [m_qtlmap_analyse_multitrait_DA]')
    allocate (somppy(ng),STAT=stat)
    call check_allocate(stat,'somppy [m_qtlmap_analyse_multitrait_DA]')
    allocate (somppm(ng),STAT=stat)
    call check_allocate(stat,'somppm [m_qtlmap_analyse_multitrait_DA]')
    allocate (somppdf(dg%np),STAT=stat)
    call check_allocate(stat,'somppdf [m_qtlmap_analyse_multitrait_DA]')
    allocate (carppdf(dg%np),STAT=stat)
    call check_allocate(stat,'carppdf [m_qtlmap_analyse_multitrait_DA]')
    allocate (somppydf(dg%np),STAT=stat)
    call check_allocate(stat,'somppydf [m_qtlmap_analyse_multitrait_DA]')
    allocate (ap(dg%np),STAT=stat)
    call check_allocate(stat,'ap [m_qtlmap_analyse_multitrait_DA]')
    ap=0.d0
    allocate (xmoyp(dg%np),STAT=stat)
    call check_allocate(stat,'xmoyp [m_qtlmap_analyse_multitrait_DA]')
    xmoyp=0.d0
    allocate (std(dg%np),STAT=stat)
    call check_allocate(stat,'std [m_qtlmap_analyse_multitrait_DA]')
    std=0.d0
    allocate (am(dg%nm),STAT=stat)
    call check_allocate(stat,'am [m_qtlmap_analyse_multitrait_DA]')
    am=0.d0
    allocate (xmoym(dg%nm),STAT=stat)
    call check_allocate(stat,'xmoym [m_qtlmap_analyse_multitrait_DA]')
    xmoym=0.d0
    allocate (fp_1(dg%np),STAT=stat)
    call check_allocate(stat,'fp_1 [m_qtlmap_analyse_multitrait_DA]')
    allocate (fm_1(dg%nm),STAT=stat)
    call check_allocate(stat,'fm_1 [m_qtlmap_analyse_multitrait_DA]')

  end subroutine init_analyse_DA_1QTL
!!***


!!****f* m_qtlmap_analyse_multitrait_DA/end_analyse_DA_1QTL
!!  NAME
!!    end_analyse_DA_1QTL
!!  DESCRIPTION
!!    deallocation of buffer/solution arrays
!!  SOURCE
  subroutine end_analyse_DA_1QTL
    deallocate (sig1)
    deallocate (xmu1p)
    deallocate (xmu1m)
    deallocate (fp0)
    deallocate (fp1)
    deallocate (fm0)
    deallocate (fm1)
    deallocate (sompp)
    deallocate (sompm)
    deallocate (sompmy)
    deallocate (carpp)
    deallocate (carpm)
    deallocate (somppy)
    deallocate (somppm)
    deallocate (somppdf)
    deallocate (carppdf)
    deallocate (somppydf)
    deallocate (ap)
    deallocate (xmoyp)
    deallocate (std)
    deallocate (xmoym)
    deallocate (am)
    deallocate (fp_1)
    deallocate (fm_1)
  end subroutine end_analyse_DA_1QTL
!!***


!!****f* m_qtlmap_analyse_multitrait_DA/perform_DA_1qtl
!!  NAME
!!    perform_DA_1qtl
!!  DESCRIPTION
!!    Calcul des combinaisons lineaires des performances pour la position testee
!!  INPUTS
!!    ch : index of the chromosome tested
!!    n  : the position tested
!!  npo  : number of position tested among the chromosome
!!
!!  OUTPUTS
!!   yda  :
!!  coeff :
!!
!!  SOURCE
!!
  subroutine perform_DA_1qtl(dataset,spt,ch,n,yda,coeff,npo)
!
    type(QTLMAP_DATASET)       ,intent(in) :: dataset
    type(PDD_BUILD)            ,intent(in) :: spt
    integer         ,intent(in)                                :: ch ! chromosome
    real (kind=dp)  ,intent(out) ,dimension(dataset%genea%nd)  :: yda
    integer         ,intent(in)                          :: n,npo
    real (kind=dp)  ,intent(inout),dimension(dataset%phenoModel%ncar,npo)   :: coeff

    integer                                              :: ic,k,kkd,nm1,nm2,jm,ip,ngeno1,ngeno2,nd1,nd2,kd,ig
    integer                                              :: nit,ind,ifail,it,NCV,IRANKX,LDX1,NI,IWKI
    real(kind=dp)  , dimension(:,:),allocatable          :: Xy,pb
    real(kind=dp)  , dimension(:),allocatable            :: WT, WK
    integer        , dimension(:),allocatable            :: NIG, ING, ISX
    integer        , dimension(:,:),allocatable          :: NIGP
    integer                                              :: NGP,LDX,M,NBX,LDCVM,LDE,LDCVX,IWK

    real (kind=dp) , parameter                          :: TOL=0.01d0
    real(kind=dp)  , dimension(:,:) ,allocatable        :: CVM, E, CVX

    character (len=1)                                   :: WEIGHT='W'
    real (kind=dp)                                        ::s1,s2
    type(GENEALOGY_BASE) , pointer :: dg
    type(PHENOTYPE_BASE) , pointer :: dpa
    type(DATAMODEL_BASE) , pointer :: dpm

    dpm => dataset%phenoModel
    dg => dataset%genea
    dpa => dataset%phenoAnimal

! Initialisation des parametres
    NGP=2

! Allocation  tableau des poids (4*nd est mal ajust�)
    allocate(pb(maxval(spt%ndesc),NGP)) !

  do ic=1,dpm%ncar
! Initialisation des poids des perf dans chaque groupe genotypique
    nit=0
    do ip=1,dg%np
      nm1=dg%nmp(ip)+1
      nm2=dg%nmp(ip+1)
      do jm=nm1,nm2
         ngeno1=spt%ngenom(ch,jm)+1
         ngeno2=spt%ngenom(ch,jm+1)
         do ig=ngeno1,ngeno2
            nd1=spt%ngend(ch,ig)+1
            nd2=spt%ngend(ch,ig+1)
            do kd=nd1,nd2
               kkd=spt%ndesc(ch,kd)
               if(count(dpa%presentc(:,kkd))==dpm%ncar) then
                   if (NGP.eq.2) then
                    pb(kd,1)=spt%pdd(ch,kd,1,n)+spt%pdd(ch,kd,2,n)
                    pb(kd,2)=spt%pdd(ch,kd,3,n)+spt%pdd(ch,kd,4,n)
                   else
                     if (NGP.eq.4) then
                       pb(kd,1)=spt%pdd(ch,kd,1,n)
                       pb(kd,2)=spt%pdd(ch,kd,2,n)
                       pb(kd,3)=spt%pdd(ch,kd,3,n)
                       pb(kd,4)=spt%pdd(ch,kd,4,n)
                     else
                      if (NGP.Eq.3) then
                        pb(kd,1)=spt%pdd(ch,kd,1,n)
                        pb(kd,2)=spt%pdd(ch,kd,2,n)+spt%pdd(ch,kd,3,n)
                        pb(kd,3)=spt%pdd(ch,kd,4,n)
                      end if ! NGP 3
                    end if ! NGP 4
                   end if ! NGP 2
                   nit=nit+1
                 end if ! PRESENTC
              end do
           end do
        end do
     end do
   end do ! ic
! Initialisation des dimensions du probleme
    NI=nit*NGP

! Initialisation des parametres
    LDX=NI
    M=dpm%ncar
    NBX=dpm%ncar
    LDCVM=NGP
    if(dpm%ncar> NGP) then
      LDE=dpm%ncar
    else
      LDE=NGP
    end if
    LDCVX=dpm%ncar
    if(NI> NGP-1) then
      IWKi=NI
    else
      IWKi=NGP-1
    end if
    IWK=NI*NBX+5*(NBX-1)+(IWKi+1)*NBX+2
    IFAIL=0

! Allocation des vecteurs et tables
    allocate(Xy(LDX,dpm%ncar))
    allocate(WT(LDX))
    allocate(WK(IWK))
    allocate(NIG(NGP))
    allocate(ING(NI))
    allocate(ISX(M))
    allocate(CVM(LDCVM,NBX))
    allocate(E(LDE,6))
    allocate(CVX(LDCVX,NGP-1))

! Initialisation des vecteurs et tables
    Xy=0.d0
    WT=0.d0
    WK=0.d0
    NIG=nit
    ING=0
    ISX=1
    CVM=0.d0
    E=0.d0
    CVX=0.d0
! determination de la VC � la position
    ind=0
    LDX1=0
    ISX(:)=1 ! all trait are included in the analysis
    do while (ind < NGP)
         do ip=1,dg%np
          nm1=dg%nmp(ip)+1
          nm2=dg%nmp(ip+1)
          do jm=nm1,nm2
            nd1=dg%ndm(jm)+1
            nd2=dg%ndm(jm+1)
            ngeno1=spt%ngenom(ch,jm)+1
            ngeno2=spt%ngenom(ch,jm+1)
            do ig=ngeno1,ngeno2
              nd1=spt%ngend(ch,ig)+1
              nd2=spt%ngend(ch,ig+1)
              do kd=nd1,nd2
                kkd=spt%ndesc(ch,kd)
                if(count(dpa%presentc(:,kkd))==dpm%ncar) then
                  LDX1=LDX1+1
                  ING(LDX1)=1+ind
                  Xy(LDX1,:)=dpa%y(:,kkd)
                  WT(LDX1)=pb(kd,1+ind)
                end if
              end do
            end do
         end do
       end do
      ind=ind+1
    end do

    s1=0.d0
    s2=s1
    do kd=1,nit
       s1=s1+WT(kd)
   end do
    do kd=1+nit,nit*2
      s2=s2+WT(kd)
    end do
!    print*,'sum', s1,s2
    call MATH_QTLMAP_G03ACF(WEIGHT,NI,M,Xy,LDX,ISX,NBX,ING,NGP,WT,NIG,CVM,LDCVM,E &
                          ,LDE,NCV,CVX,LDCVX,TOL,IRANKX,IFAIL)
    if (ifail.ne.0) then
         do it=1,LDX1
!           call log_mess('LDX1 '//trim(str(LDX1))//trim(str(WT(LDX1),DEBUG_DEF)))
         end do
         call stop_application('Analyse discriminante sous H1; ifail='//trim(str(ifail)))
    end if

    do kd=1,dg%nd
         yda(kd)=0.d0
         do ic=1,dpm%ncar
          yda(kd)=CVX(ic,1)*dpa%y(ic,kd)+yda(kd)
         end do
    end do
    do ic=1,dpm%ncar
      coeff(ic,n)=CVX(ic,1)
    end do

    deallocate(pb)
    deallocate(Xy)
    deallocate(WT)
    deallocate(WK)
    deallocate(NIG)
    deallocate(ING)
    deallocate(ISX)
    deallocate(CVM)
    deallocate(E)
    deallocate(CVX)
  end subroutine perform_DA_1qtl
!!***


!!****f* m_qtlmap_analyse_multitrait_DA/opti_DA_0qtl
!!  NAME
!!    opti_DA_0qtl
!!  DESCRIPTION
!!    Calcul de la vraisemblance 0 QTL, 1 caractere
!!  SOURCE
  subroutine opti_DA_0qtl(dataset)
    type(QTLMAP_DATASET)       ,intent(in) :: dataset
    !
    ! Divers

    integer iuser(1)
    double precision  ,dimension(:),allocatable :: borni,borns,par
    real (kind=dp) :: user(1)
    integer :: npar,ibound,ip,jm,ifail
    type(GENEALOGY_BASE) , pointer :: dg
    type(PHENOTYPE_BASE) , pointer :: dpa
    type(DATAMODEL_BASE) , pointer :: dpm

    dpm => dataset%phenoModel
    dg => dataset%genea
    dpa => dataset%phenoAnimal

    p_dpm => dataset%phenoModel
    p_dg => dataset%genea
    p_dpa => dataset%phenoAnimal

    !******************************************************************************
    ! Parametres de maximisation
    npar=(2*dg%np)+dpa%nmumest(1)

    allocate (borni(npar))
    allocate (borns(npar))
    allocate (par(npar))

    ibound=0
    do ip=1,dg%np
       borni(ip)=1.d-6
       borns(ip)=1.d6
       borni(ip+dg%np)=-1d6
       borns(ip+dg%np)=1.d6
    end do
    do jm=1,dpa%nmumest(1)
       borni(2*dg%np+jm)=-1.d6
       borns(2*dg%np+jm)=1.d6
    end do
    !
    ! Point de depart
    do ip=1,dg%np
       par(ip)=sig0(ip)
       par(ip+dg%np)=xmu0p(ip)
    end do
    do jm=1,dpa%nmumest(1)
       par(2*dg%np+jm)=xmu0m(jm)
    end do
    !
    ! Optimisation de la vraisemblance
    ifail=0

    call minimizing_funct(dataset,npar,ibound,funct_0qtl,borni,borns,par,f0,iuser,user,ifail)

    do ip=1,dg%np
       sig1(ip)=par(ip)
       xmu1p(ip)=par(ip+dg%np)
    end do
    do jm=1,dpa%nmumest(1)
       xmu1m(jm)=par(2*dg%np+jm)
    end do
    !
    deallocate (borni)
    deallocate (borns)
    deallocate (par)

  end subroutine opti_DA_0qtl
!!***


!!****f* m_qtlmap_analyse_multitrait_DA/funct_0qtl
!!  NAME
!!    funct_0qtl
!!  DESCRIPTION
!!    Calcul de la vraisemblance et de ses derivees partielles sous H0
!!  SOURCE
  subroutine funct_0qtl(n,x,f,iuser,user)
    integer         , intent(in)                  :: n
    !
    ! Tableaux dimensionnes selon n, le nombre de parametres a estimer
    !
    real (kind=dp)      ,dimension(n), intent(in) :: x
    real (kind=dp)  , intent(inout)               :: f
    integer ,       dimension(1), intent(in)      :: iuser
    real (kind=dp)      ,dimension(1), intent(in) :: user

    integer :: ip,nm1,nm2,jm
    integer :: nestim,nest

    real (kind=dp) :: sig,var,vpf,xmup,xmup2
    real (kind=dp) :: vdf,somxmu,xmum,xmum2,xmupm

    !******************************************************************************
    f=0.d0
    nestim=0
    do ip=1,p_dg%np
       sig=x(ip)
       var=sig*sig
       xmup=x(ip+p_dg%np)
       xmup2=xmup*xmup
       vdf=carydf(ip)+dble(effdf(ip))*xmup2-2.d0*xmup*somydf(ip)
       vdf=0.5d0*vdf/var
       fp0(ip)=vdf+(dble(effdf(ip))*dlog(sig))
       f=f+fp0(ip)

       nm1=p_dg%nmp(ip)+1
       nm2=p_dg%nmp(ip+1)
       somxmu=0.d0
       nest=0
       do jm=nm1,nm2
          !AJOUT OFI. le tableau 'estime' depend d un caractere
          ! on ne peut donc pas faire de if estime(ic,jm)
          ! A VALIDEr PAR HELENE : si un caractere est non estimable alors tous les caracteres ne le sont pas...

          if( count(p_dpa%estime(:,jm))== p_dpm%ncar )  then
             nest=nest+1
             if(nest.le.estmum(ip)) then
                nestim=nestim+1
                xmum=x(2*p_dg%np+nestim)
                somxmu=somxmu+xmum
             else
                xmum=-somxmu
             end if
             xmum2=xmum*xmum
             xmupm=2.d0*(xmup+xmum)
             vpf=cary(jm)+dble(eff(jm))*(xmup2+xmum2+2.d0*xmup*xmum)-xmupm*somy(jm)
             vpf=0.5d0*vpf/var
             fm0(jm)=vpf+(dble(eff(jm))*dlog(sig))
             f=f+fm0(jm)
             fp0(ip)=fp0(ip)+fm0(jm)
          else
             fm0(jm)=0.d0
          end if
       end do
    end do

  end subroutine funct_0qtl
!!***

!!****f* m_qtlmap_analyse_multitrait_DA/opti_DA_1qtl
!!  NAME
!!    opti_DA_1qtl
!!  DESCRIPTION
!!    Computing statistical test across the chromosome
!!
!!  HISTORY
!!   01/04/2009 : -- imported this function into this module
!!              : -- export all printing outside from this module
!!  SOURCE
  subroutine opti_DA_1qtl(dataset,spt,ch,ix,n,yda,lrtsol)
    type(QTLMAP_DATASET)       ,intent(in) :: dataset
    type(PDD_BUILD)    ,target ,intent(in) :: spt

    integer , intent(in)                                     :: ch,n,ix
    type(TYPE_LRT_SOLUTION)  , intent(inout)                 :: lrtsol
    real (kind=dp)  ,intent(in) ,dimension(dataset%genea%nd) :: yda
    !
    ! Divers
    integer :: iuser(1)
    double precision  ,dimension(:),allocatable :: borni,borns,par
    double precision    :: user(1)

    integer :: ibound,ip,jm,ngeno1,ngeno2,ilong,nm1,nm2,ig
    integer :: npar,nd1,nd2,kd,kkd,ifail,ii
    integer :: nestim,nest,ifem,indam
    real (kind=dp) :: pp, pm, somxmu,xlrt_t,f1
    type(GENEALOGY_BASE) , pointer :: dg
    type(PHENOTYPE_BASE) , pointer :: dpa
    type(DATAMODEL_BASE) , pointer :: dpm

    dpm => dataset%phenoModel
    dg => dataset%genea
    dpa => dataset%phenoAnimal

    p_dpm => dataset%phenoModel
    p_dg => dataset%genea
    p_dpa => dataset%phenoAnimal
    p_spt => spt

    !******************************************************************************
    ! Calcul de la vraisemblance sous H1
    ! Parametres de maximisation
    ibound=0
    npar=(3*dg%np)+dpa%nmumest(1)+dpa%namest(1)

    call log_mess("**** OPTI_DA_1QTL *****",DEBUG_DEF)
    call log_mess("NP:"//str(dg%np),DEBUG_DEF)
    call log_mess("NPAR:"//str(npar),DEBUG_DEF)

    allocate (borni(npar))
    allocate (borns(npar))
    allocate (par(npar))

    current_ch = ch

    do ip=1,dg%np
       borni(ip)=1.d-6
       borns(ip)=1.d6
       borni(ip+dg%np)=-1.d6
       borns(ip+dg%np)=1.d6
       borni(2*dg%np+ip)=-1.d6
       borns(2*dg%np+ip)=1.d6
    end do
    do jm=1,dpa%nmumest(1)
       borni(3*dg%np+jm)=-1.d6
       borns(3*dg%np+jm)=1.d6
    end do
    do jm=1,dpa%namest(1)
       borni(3*dg%np+dpa%nmumest(1)+jm)=-1.d6
       borns(3*dg%np+dpa%nmumest(1)+jm)=1.d6
    end do
    par=0.d0
    !
    ! Point de depart
    do ip=1,dg%np
       par(ip)=sig1(ip)
       par(ip+dg%np)=xmu1p(ip)
       par(2*dg%np+ip)=0.d0
    end do

    do jm=1,dpa%nmumest(1)
       par(3*dg%np+jm)=xmu1m(jm)
    end do
    do jm=1,dpa%namest(1)
       par(3*dg%np+dpa%nmumest(1)+jm)=0.d0
    end do
    !
    !A la position en cours
    do ip=1,dg%np
      nm1=dg%nmp(ip)+1
      nm2=dg%nmp(ip+1)
      somppdf(ip)=0.d0
      carppdf(ip)=0.d0
      somppydf(ip)=0.d0
      do jm=nm1,nm2
         ngeno1=spt%ngenom(ch,jm)+1
         ngeno2=spt%ngenom(ch,jm+1)
         do ig=ngeno1,ngeno2
            sompp(ig)=0.d0
            carpp(ig)=0.d0
            somppy(ig)=0.d0
            sompm(ig)=0.d0
            carpm(ig)=0.d0
            sompmy(ig)=0.d0
            somppm(ig)=0.d0
            nd1=spt%ngend(ch,ig)+1
            nd2=spt%ngend(ch,ig+1)
            do kd=nd1,nd2
               kkd=spt%ndesc(ch,kd)
               if(dpa%presentc(1,kkd)) then
                  pp=-spt%pdd(ch,kd,1,n)-spt%pdd(ch,kd,2,n)+spt%pdd(ch,kd,3,n)+spt%pdd(ch,kd,4,n)
                  pm=-spt%pdd(ch,kd,1,n)+spt%pdd(ch,kd,2,n)-spt%pdd(ch,kd,3,n)+spt%pdd(ch,kd,4,n)
                  !AJOUT OFI. le tableau 'estime' depend d un caractere
                  ! on ne peut donc pas faire de if estime(ic,jm)
                  !
                  ! A VALIDEr PAR HELENE : si un caractere est non estimable alors tous les caracteres ne le sont pas...

                   if( count(dpa%estime(:,jm))== dpm%ncar )  then
                     sompp(ig)=sompp(ig)+pp
                     carpp(ig)=carpp(ig)+pp*pp
                     somppy(ig)=somppy(ig)+pp*yda(kkd)
                     sompm(ig)=sompm(ig)+pm
                     carpm(ig)=carpm(ig)+pm*pm
                     sompmy(ig)=sompmy(ig)+pm*yda(kkd)
                     somppm(ig)=somppm(ig)+pp*pm
                   else
                     somppdf(ip)=somppdf(ip)+pp
                     carppdf(ip)=carppdf(ip)+pp*pp
                     somppydf(ip)=somppydf(ip)+pp*yda(kkd)
                   end if
                 end if
              end do
           end do
        end do
     end do
     !
     ! Optimisation de la vraisemblance a la position en cours
     ifail=0

     call minimizing_funct(dataset,npar,ibound,funct_1qtl,borni,borns,par,f1,iuser,user,ifail)

     xlrt_t=-2.d0*(f1-f0)
     call log_mess('N:'//str(n)//" LRT:"//str(xlrt_t)//" F0:"//str(f0)//" F1:"//str(f1),DEBUG_DEF)
     do ii=1,dg%np
        call lrtsol%LRT_SIRES(ii)%add1p(dataset,ch,n,(-2.d0*(fp1(ii)-fp0(ii))))
        lrtsol%pater_eff(ch,ii,n)=par(2*dg%np+ii)
     end do

     do ii=1,dg%nm
        call lrtsol%LRT_DAMS(ii)%add1p(dataset,ch,n,(-2.d0*(fm1(ii)-fm0(ii))))
     end do

     do ii=1,dpa%namest(1)
        lrtsol%mater_eff(ch,ii,n)=par(3*dg%np+dpa%nmumest(1)+ii)
     end do

     if(lrtsol%lrtmax(0).le.xlrt_t) then
          lrtsol%lrtmax(0)= xlrt_t
          lrtsol%nxmax(0)= n
          lrtsol%chrmax(0)= ch
          f_1=f1       ! ic
          nestim=0
          do ip=1,dg%np
             ap(ip)=par(2*dg%np+ip)
             xmoyp(ip)=par(dg%np+ip)
             std(ip)=par(ip)
             fp_1(ip)=fp1(ip)
             nm1=dg%nmp(ip)+1
             nm2=dg%nmp(ip+1)
             nest=0
             somxmu=0.d0
             do jm=nm1,nm2
                fm_1(jm)=fm1(jm)
                !AJOUT OFI. le tableau 'estime' depend d un caractere
                  ! on ne peut donc pas faire de if estime(ic,jm)
                  !
                  ! A VALIDEr PAR HELENE : si un caractere est non estimable alors tous les caracteres ne le sont pas...
                  if( count(dpa%estime(:,jm))== dpm%ncar )  then
                !if (estime(ic,jm)) then
                   nest=nest+1
                   ifem=dg%repfem(jm)
                   indam=dpa%iam(1,ifem)
                   am(jm)=par(3*dg%np+dpa%nmumest(1)+indam)
                   if (nest.le.estmum(ip)) then
                      nestim=nestim+1
                      xmoym(jm)=par(3*dg%np+nestim)
                      somxmu=somxmu+xmoym(jm)
                   else
                      xmoym(jm)=-somxmu
                   end if
                else
                   am(jm)=0.d0
                   xmoym(jm)=0.d0
                end if
             end do
          end do
       end if

       call lrtsol%LRT%add1p(dataset,ch,n,xlrt_t)


!
    deallocate (borni)
    deallocate (borns)
    deallocate (par)

    return
  end subroutine opti_DA_1qtl
!!***

!!****f* m_qtlmap_analyse_multitrait_DA/funct_1qtl
!!  NAME
!!    funct_1qtl
!!  DESCRIPTION
!!    Calcul de la vraisemblance et de ses derivees partielles sous Hypothese 1QTL 1carac
!!
!!  SOURCE
  subroutine funct_1qtl(n,x,f,iuser,user)

    integer       , intent(in)          :: n    ! number of variables
    double precision, intent(in)        :: x(n) ! Position
    double precision, intent(out)       :: f    ! Result value at the position
    double precision      user(1)
    integer   iuser(1)


    !
    ! Divers
    integer :: ip,nestim,nm1,nm2,nest,jm,ifem,indam,i
    integer :: ngeno1,ngeno2,ig
    real (kind=dp)  :: sig,var,xmup,xmup2,vdf,somxmu,xmum,xmum2
    real (kind=dp)  :: xmupm,vmere,z,vpf,x_ap,x_ap2,x_am,x_am2


    !******************************************************************************
    f=0.d0
    nestim=0

    do ip=1,p_dg%np
       sig=x(ip)
       var=sig*sig
       xmup=x(ip+p_dg%np)
       xmup2=xmup*xmup
       x_ap=x(2*p_dg%np+ip)
       x_ap2=x_ap*x_ap
       vdf=carydf(ip)+dble(effdf(ip))*xmup2-2.d0*xmup*somydf(ip)
       vdf=vdf+x_ap2*carppdf(ip)+2.d0*xmup*x_ap*somppdf(ip)-2.d0*x_ap*somppydf(ip)
       vdf=0.5d0*vdf/var
       fp1(ip)=vdf+dble(effdf(ip))*dlog(sig)
       f=f+fp1(ip)
       nm1=p_dg%nmp(ip)+1
       nm2=p_dg%nmp(ip+1)
       somxmu=0.d0
       nest=0
       do jm=nm1,nm2
        !AJOUT OFI. le tableau 'estime' depend d un caractere
        ! on ne peut donc pas faire de if estime(ic,jm)
        !
        ! A VALIDEr PAR HELENE : si un caractere est non estimable alors tous les caracteres ne le sont pas...
          if( count(p_dpa%estime(:,jm))== p_dpm%ncar )  then
             nest=nest+1
             if(nest.le.estmum(ip)) then
                nestim=nestim+1
                xmum=x(3*p_dg%np+nestim)
                somxmu=somxmu+xmum
             else
                xmum=-somxmu
             end if
             xmum2=xmum*xmum
             xmupm=2.d0*(xmup+xmum)
             ifem=p_dg%repfem(jm)
             indam=p_dpa%iam(1,ifem)
             x_am=x(3*p_dg%np+p_dpa%nmumest(1)+indam)
             x_am2=x_am*x_am
             z=cary(jm)+dble(eff(jm))*(xmup2+xmum2+2.d0*xmup*xmum)-xmupm*somy(jm)
             vmere=0.d0
             ngeno1=p_spt%ngenom(current_ch,jm)+1
             ngeno2=p_spt%ngenom(current_ch,jm+1)
             do ig=ngeno1,ngeno2
                vpf=z+x_ap2*carpp(ig)+x_am2*carpm(ig)            &
                     +2.d0*x_ap*x_am*somppm(ig)                   &
                     +xmupm*(x_ap*sompp(ig)+x_am*sompm(ig))       &
                     -2.d0*x_ap*somppy(ig)-2.d0*x_am*sompmy(ig)
                vpf=-0.5d0*vpf/var
                vmere=vmere+p_spt%probg(current_ch,ig)*dexp(vpf)
             end do
             !               if (vmere.lt.1.e-8) then
             !                 write (nficerr,*) '** WARNING : vmere.eq.0 ',
             !   $                 'in funct_1qtl.F **'
             !               fm1(jm)=(dble(eff(jm))*dlog(sig))
             !         else
             fm1(jm)=-dlog(vmere)+(dble(eff(jm))*dlog(sig))
             !         end if
             f=f+fm1(jm)
             fp1(ip)=fp1(ip)+fm1(jm)
          else
             fm1(jm)=0.d0
          end if
       end do
    end do

  end subroutine funct_1qtl
!!***

!!****f* m_qtlmap_analyse_multitrait_DA/set_solution_hypothesis1
!!  NAME
!!    set_solution_hypothesis1
!!  DESCRIPTION
!!
!!  SOURCE
 subroutine set_solution_hypothesis1(dataset,incsol)
       type(QTLMAP_DATASET)       ,intent(in) :: dataset
       type(TYPE_INCIDENCE_SOLUTION)      ,intent(inout)    :: incsol

       integer :: nteff,maxNbPar,jm,ifem,ip,ieff
       type(GENEALOGY_BASE) , pointer :: dg
       type(PHENOTYPE_BASE) , pointer :: dpa
       type(DATAMODEL_BASE) , pointer :: dpm

       dpm => dataset%phenoModel
       dg => dataset%genea
       dpa => dataset%phenoAnimal

       allocate (incsol%sig(1,dg%np))

       incsol%hypothesis=1
       !  Mean family , Qtl effect
       nteff = 2
       if ( count(dpa%estime(1,:))  > 0 ) nteff = 4

       maxNbPar = max(dg%np,count(dpa%estime(1,:)))
       allocate (incsol%groupeName(nteff))
       allocate (incsol%nbParameterGroup(nteff))
       allocate (incsol%parameterName(nteff,maxNbPar))
       allocate (incsol%paramaterValue(nteff,maxNbPar))
       allocate (incsol%parameterVecsol(nteff,maxNbPar))
       allocate (incsol%parameterPrecis(nteff,maxNbPar))
       allocate (incsol%qtl_groupeName(1,1))

       incsol%parameterPrecis=0.d0
       incsol%parameterVecsol=.true.

       do ip=1,dg%np
            incsol%sig(1,ip) = std(ip)*dpm%sigt(1)
       end do

       ieff=1
         incsol%qtl_groupeName(1,1)=ieff
         incsol%groupeName(ieff) = 'Sire Qtl effect'
         incsol%nbParameterGroup(ieff)=dg%np

         do ip=1,dg%np
           incsol%parameterName(ieff,ip)='Sire '//trim(dg%pere(ip))
           incsol%paramaterValue(ieff,ip) = ap(ip)*dpm%sigt(1)! en DA quel sigt prendre ? *sigt(ic)
         end do

         if ( count(dpa%estime(1,:)) > 0 ) then
           ieff = ieff +1
           incsol%groupeName(ieff) = 'Dam Qtl effect'
           incsol%nbParameterGroup(ieff)=dpa%namest(1)
           ifem=0
           do jm=1,dg%nm
             if ( dpa%estime(1,jm)) then
              ifem=ifem+1
              incsol%parameterName(ieff,ifem)='Dam '//trim(dg%mere(jm))
              incsol%paramaterValue(ieff,ifem) = am(jm)*dpm%sigt(1)
              incsol%parameterVecsol(ieff,ifem) = .true.
             end if
           end do
         end if


       ieff = ieff +1
       incsol%groupeName(ieff) = 'Mean Sire'
       incsol%nbParameterGroup(ieff)=dg%np

       do ip=1,dg%np
           incsol%parameterName(ieff,ip)='Sire '//trim(dg%pere(ip))
           incsol%paramaterValue(ieff,ip) = xmoyp(ip)*dpm%sigt(1) + dpm%xmut(1)
           incsol%parameterVecsol(ieff,ip) = .true.
       end do

       if ( count(dpa%estime(1,:)) > 0 ) then
           ieff = ieff +1
           incsol%groupeName(ieff) = 'Mean dam'
           incsol%nbParameterGroup(ieff)=count(dpa%estime(1,:))
           ifem=0
         do ip=1,dg%np
           do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
             if ( dpa%estime(1,jm)) then
              ifem=ifem+1
              incsol%parameterName(ieff,ifem)='Dam '//trim(dg%mere(jm))//" [Sire "//trim(dg%pere(ip))//"]"
              incsol%paramaterValue(ieff,ifem) = xmoym(jm)*dpm%sigt(1)
              incsol%parameterVecsol(ieff,ifem) = .true.
             end if
           end do
         end do
         end if

       end subroutine set_solution_hypothesis1
!!***


end module m_qtlmap_analyse_multitrait_DA
