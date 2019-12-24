!!****m* ANALYSE/m_qtlmap_analyse_gen
!!  NAME
!!    m_qtlmap_analyse_gen
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
   module m_qtlmap_analyse_gen
      use m_qtlmap_base
      use m_qtlmap_types
      use m_qtlmap_log
      use m_qtlmap_output_handler
      use m_qtlmap_haplotype_ldla

      implicit none
      save

!!****v* m_qtlmap_analyse_gen/estmum
!!  NAME
!!   estmum
!!  DESCRIPTION
!!   number of female estimable in a half sib family ip
!!  DIMENSIONS
!!    np
!!  NOTES
!!
!!***
      integer      ,dimension(:),pointer, public          :: estmum
!!****v* m_qtlmap_analyse_gen/eff
!!  NAME
!!   eff
!!  DESCRIPTION
!!   number of progenies with a trait value in the full sib family jm
!!  DIMENSIONS
!!    nm
!!  NOTES
!!
!!***
      integer      ,dimension(:),pointer     ,public      :: eff
!!****v* m_qtlmap_analyse_gen/effp
!!  NAME
!!   effp
!!  DESCRIPTION
!!   number of progenies with a trait value in the half sib family ip
!!  DIMENSIONS
!!    np
!!  NOTES
!!
!!***
      integer       ,dimension(:),allocatable,public      :: effp

      !$omp threadprivate (estmum,eff,effp)

!!****v* m_qtlmap_analyse_gen/somcd
!!  NAME
!!   somcd
!!  DESCRIPTION
!!   sum of the censured data in the full sib family jm
!!  DIMENSIONS
!!    nm
!!  NOTES
!!
!!***
      real (kind=dp)       ,dimension(:),pointer,public   :: somcd

!!****v* m_qtlmap_analyse_gen/somy
!!  NAME
!!   somy
!!  DESCRIPTION
!!   sum of the trait value in the full sib family jm
!!  DIMENSIONS
!!    nm
!!  NOTES
!!
!!***
      real (kind=dp)       ,dimension(:),pointer,public   :: somy

!!****v* m_qtlmap_analyse_gen/cary
!!  NAME
!!   cary
!!  DESCRIPTION
!!   sum of square of the trait value in the full sib family jm
!!  DIMENSIONS
!!    nm
!!  NOTES
!!
!!***
      real (kind=dp)       ,dimension(:),pointer,public   :: cary

      !$omp threadprivate (somcd,somy,cary)

!!****v* m_qtlmap_analyse_gen/effdf
!!  NAME
!!   effdf
!!  DESCRIPTION
!!   number of progenies (only with dam estimable) with a trait value in the half sib family ip
!!  DIMENSIONS
!!    np
!!  NOTES
!!
!!***
      integer              ,dimension(:),pointer,public   :: effdf
      !$omp threadprivate (effdf)

!!****v* m_qtlmap_analyse_gen/somcddf
!!  NAME
!!   somcddf
!!  DESCRIPTION
!!   sum of the censured data in the half sib family ip (only with dam estimable)
!!  DIMENSIONS
!!    np
!!  NOTES
!!
!!***
      real (kind=dp)       ,dimension(:),pointer,public   :: somcddf

!!****v* m_qtlmap_analyse_gen/somydf
!!  NAME
!!   somydf
!!  DESCRIPTION
!!   sum of the trait value in the half sib family ip (only with dam estimable)
!!  DIMENSIONS
!!    np
!!  NOTES
!!
!!***
      real (kind=dp)       ,dimension(:),pointer,public   :: somydf

!!****v* m_qtlmap_analyse_gen/carydf
!!  NAME
!!   carydf
!!  DESCRIPTION
!!  sum of square of the trait value in the half sib family ip (only with dam estimable)
!!  DIMENSIONS
!!    np
!!  NOTES
!!
!!***
      real (kind=dp)       ,dimension(:),pointer,public   :: carydf
!!****v* m_qtlmap_analyse_gen/sig0
!!  NAME
!!   sig0
!!  DESCRIPTION
!!    variance of traits values for the half sib family ip
!!  DIMENSIONS
!!    np
!!  NOTES
!!
!!***
      real (kind=dp)       ,dimension(:),allocatable,public   :: sig0
!!****v* m_qtlmap_analyse_gen/xmu0p
!!  NAME
!!   xmu0p
!!  DESCRIPTION
!!    mean of traits values for the half sib family ip
!!  DIMENSIONS
!!    np
!!  NOTES
!!
!!***
      real (kind=dp)       ,dimension(:),allocatable,public   :: xmu0p
!!****v* m_qtlmap_analyse_gen/xmu0m
!!  NAME
!!   xmu0p
!!  DESCRIPTION
!!    mean of traits values for the full sib family jm
!!  DIMENSIONS
!!    nm
!!  NOTES
!!
!!***
      real (kind=dp)       ,dimension(:),allocatable,public   :: xmu0m
      !$omp threadprivate (somcddf,somydf,carydf,sig0,xmu0p,xmu0m)


      public :: init_analyse_gen
      public :: optinit
      public :: optinit_da
      public :: get_eff_paternal_and_total
      public :: scale_variables_dam
      public :: set_ntnivmax
      public :: set_phase_pdd
      public :: end_analyse_gen

      contains

!!****f* m_qtlmap_analyse_gen/init_analyse_gen
!!  NAME
!!   init_analyse_gen
!!  DESCRIPTION
!!   allocation of arrays provide by this module.
!!  NOTES
!!
!!  SOURCE
      subroutine init_analyse_gen(dataset)
          type(QTLMAP_DATASET)       ,intent(in)            :: dataset
          integer :: stat
          type(GENEALOGY_BASE) , pointer :: dg

          dg => dataset%genea

          allocate (estmum(dg%np),STAT=stat)
          call check_allocate(stat,'estmum [m_qtlmap_analyse_gen]')
          allocate (eff(dg%nm),STAT=stat)
          call check_allocate(stat,'eff [m_qtlmap_analyse_gen]')
          allocate (somcd(dg%nm),STAT=stat)
          call check_allocate(stat,'somcd [m_qtlmap_analyse_gen]')
          allocate (somy(dg%nm),STAT=stat)
          call check_allocate(stat,'somy [m_qtlmap_analyse_gen]')
          allocate (cary(dg%nm),STAT=stat)
          call check_allocate(stat,'cary [m_qtlmap_analyse_gen]')
          allocate (effp(dg%np))
          !call check_allocate(stat,'estime [m_qtlmap_analyse_gen]')
          allocate (effdf(dg%np),STAT=stat)
          call check_allocate(stat,'effdf [m_qtlmap_analyse_gen]')
          allocate (somcddf(dg%np),STAT=stat)
          call check_allocate(stat,'somcddf [m_qtlmap_analyse_gen]')
          allocate (somydf(dg%np),STAT=stat)
          call check_allocate(stat,'somydf [m_qtlmap_analyse_gen]')
          allocate (carydf(dg%np),STAT=stat)
          call check_allocate(stat,'carydf [m_qtlmap_analyse_gen]')
          allocate (sig0(dg%np),STAT=stat)
          call check_allocate(stat,'sig0 [m_qtlmap_analyse_gen]')
          allocate (xmu0p(dg%np),STAT=stat)
          call check_allocate(stat,'xmu0p [m_qtlmap_analyse_gen]')
          allocate (xmu0m(dg%nm),STAT=stat)
          call check_allocate(stat,'xmu0m [m_qtlmap_analyse_gen]')
      end subroutine init_analyse_gen
!!***

!!****f* m_qtlmap_analyse_gen/end_analyse_gen
!!  NAME
!!   end_analyse_gen
!!  DESCRIPTION
!!   deallocation of arrays provide by this module.
!!  NOTES
!!
!!  SOURCE
      subroutine end_analyse_gen

          deallocate (estmum)
          deallocate (eff)
          deallocate (somcd)
          deallocate (somy)
          deallocate (cary)
          deallocate (effp)
          deallocate (effdf)
          deallocate (somcddf)
          deallocate (somydf)
          deallocate (carydf)
          deallocate (sig0)
          deallocate (xmu0p)
          deallocate (xmu0m)

      end subroutine end_analyse_gen
!!***

!!****f* m_qtlmap_analyse_gen/optinit
!!  NAME
!!   optinit
!!  DESCRIPTION
!!    initialisation for half sib, full sib family :
!!      * number of progenies
!!      * sum of cd
!!      * sum of y
!!      * sum of square y
!!      * sig0
!!      * mean0
!!  INPUTS
!!     ic : index of the trait
!!  NOTES
!!   Initialise les ecart types et moyennes intra-famille
!!   Sous programme appele par analyse
!!   Version XXX _ HEG030507
!!
!!  SOURCE
      subroutine optinit(dataset,ic)
      type(QTLMAP_DATASET)       ,intent(in)        :: dataset
      integer                    , intent(in)       :: ic

! Divers
      integer         :: ifem,ip,nm1,nm2,jm,nd1,nd2,kd,imumest,namest
      real (kind=dp)  :: somyp
      logical         :: estfem(dataset%genea%nfem)
      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa

      dpa => dataset%phenoAnimal
      dg => dataset%genea
!
!******************************************************************************
      effp=0
      imumest=0
      estfem=.false.
      do ip=1,dg%np
        somyp=0.d0
        sig0(ip)=0.d0
        effdf(ip)=0
    somcddf(ip)=0.d0
        somydf(ip)=0.d0
        carydf(ip)=0.d0
        estmum(ip)=0
        nm1=dg%nmp(ip)+1
        nm2=dg%nmp(ip+1)
        do jm=nm1,nm2
          eff(jm)=0
      somcd(jm)=0.d0
          somy(jm)=0.d0
          cary(jm)=0.d0
          nd1=dg%ndm(jm)+1
          nd2=dg%ndm(jm+1)
          ifem=dg%repfem(jm)
          do kd=nd1,nd2
            if(dpa%presentc(ic,kd)) then
              eff(jm)=eff(jm)+1
              somcd(jm)=somcd(jm)+dpa%cd(ic,kd)
              somy(jm)=somy(jm)+dpa%y(ic,kd)*dpa%cd(ic,kd)
              cary(jm)=cary(jm)+(dpa%y(ic,kd)*dpa%y(ic,kd))*dpa%cd(ic,kd)
            end if
          end do

          effp(ip)=effp(ip)+eff(jm)
          somyp=somyp+somy(jm)
          if( dpa%estime(ic,jm) ) then
            imumest=imumest+1
            estmum(ip)=estmum(ip)+1
            xmu0m(imumest)=somy(jm)/dble(eff(jm))
            estfem(ifem)=.true.
          else
            effdf(ip)=effdf(ip)+eff(jm)
            somcddf(ip)=somcddf(ip)+somcd(jm)
            somydf(ip)=somydf(ip)+somy(jm)
            carydf(ip)=carydf(ip)+cary(jm)
          end if
        end do
        if (effp(ip) == 0.d0) then
           call stop_application('Father ['//trim(dg%pere(ip))//'] has got no child with trait value')
        else
           xmu0p(ip)=somyp/dble(effp(ip))
        end if
        do jm=nm1,nm2
          nd1=dg%ndm(jm)+1
          nd2=dg%ndm(jm+1)
          do kd=nd1,nd2
             if(dpa%presentc(ic,kd)) then
               sig0(ip)=sig0(ip)+(dpa%y(ic,kd)-xmu0p(ip))*(dpa%y(ic,kd)-xmu0p(ip))
             end if
          end do
        end do
        if ( effp(ip) == 1) then
           call stop_application('Father ['//trim(dg%pere(ip))//'] has got only one child with trait value')
        else
           sig0(ip)=dsqrt(sig0(ip)/dble(effp(ip)-1))
        end if

        if(estmum(ip).gt.0) then
          estmum(ip)=estmum(ip)-1
          dpa%nmumest(ic)=dpa%nmumest(ic)-1
        end if
      end do
      end subroutine optinit
!!***

!!****f* m_qtlmap_analyse_gen/optinit_da
!!  NAME
!!   optinit_da
!!  DESCRIPTION
!!  INPUTS
!!     yda :
!!  NOTES
!!
!!  SOURCE
      subroutine optinit_da(dataset,yda)
       type(QTLMAP_DATASET)       ,intent(in)        :: dataset
!      integer        , intent(in)         :: ic
       real (kind=dp), dimension (dataset%genea%nd)  :: yda

! Divers
      integer         :: effp,ifem,ip,nm1,nm2,jm,nd1,nd2,kd,ic,imumest
      real (kind=dp)  :: somyp
      logical         :: estfem(dataset%genea%nfem)
      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa
      type(DATAMODEL_BASE) , pointer :: dpm

      dpm => dataset%phenoModel
      dg => dataset%genea
      dpa => dataset%phenoAnimal
!
!******************************************************************************
      ic=1
      imumest=0

      do ifem=1,dg%nfem
        estfem(ifem)=.false.
      end do
      do ip=1,dg%np
        effp=0
        somyp=0.d0
        sig0(ip)=0.d0
        effdf(ip)=0
        somydf(ip)=0.d0
        carydf(ip)=0.d0
        estmum(ip)=0
        nm1=dg%nmp(ip)+1
        nm2=dg%nmp(ip+1)
        do jm=nm1,nm2
          eff(jm)=0
          somy(jm)=0.d0
          cary(jm)=0.d0
          nd1=dg%ndm(jm)+1
          nd2=dg%ndm(jm+1)
          ifem=dg%repfem(jm)
          do kd=nd1,nd2
            if(dpa%presentc(ic,kd)) then
              eff(jm)=eff(jm)+1
              somy(jm)=somy(jm)+yda(kd)
              cary(jm)=cary(jm)+(yda(kd)*yda(kd))
            end if
          end do

          effp=effp+eff(jm)
          somyp=somyp+somy(jm)
          if(dpa%estime(ic,jm)) then
            imumest=imumest+1
            estmum(ip)=estmum(ip)+1
            xmu0m(imumest)=somy(jm)/dble(eff(jm))
            estfem(ifem)=.true.
          else
            effdf(ip)=effdf(ip)+eff(jm)
            somydf(ip)=somydf(ip)+somy(jm)
            carydf(ip)=carydf(ip)+cary(jm)
          end if
        end do
        if (effp == 0.d0) then
           call log_mess('optinit : effp <- 0.d0',WARNING_DEF)
           xmu0p(ip)= 0.d0
        else
           xmu0p(ip)=somyp/dble(effp)
        end if
        do jm=nm1,nm2
          nd1=dg%ndm(jm)+1
          nd2=dg%ndm(jm+1)
          do kd=nd1,nd2
             if(dpa%presentc(ic,kd)) then
               sig0(ip)=sig0(ip)+(yda(kd)-xmu0p(ip))*(yda(kd)-xmu0p(ip))
             end if
          end do
        end do
        sig0(ip)=dsqrt(sig0(ip)/dble(effp-1))
         if(estmum(ip).gt.0) then
          estmum(ip)=estmum(ip)-1
          dpa%nmumest(:)=imumest-1
        end if
      end do

      end subroutine optinit_da
!!***

!!****f* m_qtlmap_analyse_gen/get_eff_paternal_and_total
!!  NAME
!!   get_eff_paternal_and_total
!!  DESCRIPTION
!!   Get the number of femal with enough progeny with perf
!!  OUTPUTS
!!     effp : number of progenies in the half sib family ip with a trait value
!!     efft : number of progenies in the dataset with a trait value
!!  NOTES
!!
!!  SOURCE
      subroutine get_eff_paternal_and_total(dataset,ic,effp,efft)
             type(QTLMAP_DATASET)       ,intent(in)    :: dataset
             integer                    , intent(in)   :: ic
             integer    , dimension(dataset%genea%np) , intent(out)  :: effp
             integer                    , intent(out)  :: efft

             integer     :: ip, nm1,nm2,jm
             type(GENEALOGY_BASE) , pointer :: dg
             type(PHENOTYPE_BASE) , pointer :: dpa

             dg => dataset%genea
             dpa => dataset%phenoAnimal

             efft=0
             do ip=1,dg%np
                effp(ip)=0
                nm1=dg%nmp(ip)+1
                nm2=dg%nmp(ip+1)
                do jm=nm1,nm2
                   efft=efft+count(dpa%presentc(ic,dg%ndm(jm)+1:dg%ndm(jm+1)))
                   effp(ip)=effp(ip)+eff(jm)
                   print *,'get-eff:',efft
                end do
              end do
      end subroutine get_eff_paternal_and_total
!!***

!!****f* m_qtlmap_analyse_gen/scale_variables_dam
!!  NAME
!!   scale_variables_dam
!!  DESCRIPTION
!!
!!  INPUTS
!!    chr       : index of the chromosome
!!    ic        : index of the trait
!!   am         : qtleffect
!!   sigt       : trait variance
!!  OUTPUTS
!!     xmoym1   :
!!     am1      :
!!    subphasm  :
!!    nmsub     :
!!    submere   :
!!    impfem    :
!!  NOTES
!!
!!  SOURCE
      subroutine scale_variables_dam(dataset,spt,chr,ic,am,xmoym,sigt,submere,nmsub,am1,xmoym1,subphasm,impfem)
             type(QTLMAP_DATASET)       ,intent(in)           :: dataset
             type(PDD_BUILD)            ,intent(in)           :: spt
             integer                             ,intent(in)  :: chr
             integer              ,dimension(:),  intent(in)  :: ic
             real (kind=dp)       ,dimension(dataset%genea%nm),intent(in)  :: xmoym
             real (kind=dp)       , dimension(dataset%genea%nm),intent(in)  :: am
             real (kind=dp)                      ,intent(in)  :: sigt
             real (kind=dp)       , dimension(dataset%genea%nm),intent(out) :: xmoym1
             real (kind=dp)       , dimension(dataset%genea%nm),intent(out) :: am1
             logical              , dimension(dataset%genea%nm),intent(out) :: subphasm
             integer              , intent(out)               :: nmsub
             character(len=LEN_DEF) , dimension(dataset%genea%nm),intent(out) :: submere
             logical              , intent(out)               :: impfem


             integer     :: ip, nm1,nm2,jm,i
             logical  ,dimension(size(ic))   :: estime_ic
             type(GENEALOGY_BASE) , pointer :: dg
             type(PHENOTYPE_BASE) , pointer :: dpa
             type(DATAMODEL_BASE) , pointer :: dpm

             dpm => dataset%phenoModel
             dg => dataset%genea
             dpa => dataset%phenoAnimal

             impfem=.false.
             nmsub = 0

             do ip=1,dg%np
                nm1=dg%nmp(ip)+1
                nm2=dg%nmp(ip+1)
                do jm=nm1,nm2
                          do i=1,size(ic)
                           estime_ic(i) = dpa%estime(ic(i),jm)
                          end do
                  if( count(estime_ic(:))== size(ic) )then
                   impfem=.true.
                   nmsub=nmsub+1
                   submere(nmsub)=dg%mere(jm)
                   subphasm(nmsub)=spt%phasm(chr,jm)
                   xmoym1(nmsub)=xmoym(jm)*sigt
                   am1(nmsub)=am(jm)*sigt
                  end if
               end do
              end do
      end subroutine
!!***

!!****f* m_qtlmap_analyse_gen/set_ntnivmax
!!  NAME
!!   set_ntnivmax
!!  DESCRIPTION
!!   get the maximum number parameter to build the incidence matrix
!!  INPUTS
!!    ic        : index of the trait
!!   nqtl       : the hypothesis
!!  OUTPUTS
!!     ntnivmax : maximum number of parameter
!!     nteffmax : maximum number of effect
!!    ntlev     : number of level interaction fixed effect-qtl
!!    nbniv     : number of level fixed effect
!!  NOTES
!!
!!  SOURCE
       subroutine set_ntnivmax(dataset,ic,nqtl,ntnivmax,nteffmax,ntlev,nbniv)
          type(QTLMAP_DATASET)       ,intent(in)            :: dataset
          integer        ,intent(in)   :: ic
          integer        ,intent(in)   :: nqtl
          integer        ,intent(out)  :: ntnivmax
          integer        ,intent(out)  :: nteffmax
          integer        ,intent(out)  :: ntlev
          integer        ,intent(out)  :: nbniv

          integer                 :: nbtef,nbtco,nbtint
          integer                 :: nbtp,nbtm,jef

          type(GENEALOGY_BASE) , pointer :: dg
          type(PHENOTYPE_BASE) , pointer :: dpa
          type(DATAMODEL_BASE) , pointer :: dpm

          dpm => dataset%phenoModel
          dg => dataset%genea
          dpa => dataset%phenoAnimal

	      nbtef=dpm%modele(ic,1)
	      nbtco=dpm%modele(ic,2)
	      nbtint=dpm%modele(ic,3)
	      ntlev=1
		  nbtp=3+nbtef+nbtco
		  nbtm=nbtp+nbtint

          !compute number of level for qtl*interaction
          if(nbtint > 0) then
		     do jef=1,nbtint
		          ntlev=ntlev*dpm%nlev(dpm%modele(ic,nbtp+jef))
		     end do
		  end if

          ! number of level for fixed effect
		  nbniv=0
		  do jef=1,nbtef
		      nbniv=nbniv+dpm%nlev(dpm%modele(ic,3+jef))
		  end do


          !General mean + polygenic sire + polygenic dam  + 2*qtl*interaction+nbniv effet fixed + number of covariate
	      ! namest : nombre de fem estim
	      ! count(estime) : nombre de dam estimable
	      ntnivmax = 1 + dg%np + count(dpa%estime(ic,:)) + nqtl*(ntlev*dg%np+ntlev*dpa%namest(ic))+nbniv+nbtco &
	       + dataset%params%NB_HAPLO_PRIOR
	      nteffmax = 3+2*nqtl+nbtef+nbtco + dataset%params%NB_HAPLO_PRIOR

          call log_mess('Value NTNIV MAX:'//trim(str(ntnivmax)),VERBOSE_DEF)
          call log_mess('Value NTEFF MAX:'//trim(str(nteffmax)),VERBOSE_DEF)

     end subroutine set_ntnivmax
!!***

!!****f* m_qtlmap_analyse_gen/courbe_lin_ldla
!! NAME
!!   courbe_lin_ldla
!! DESCRIPTION
!!   Print the information about result (curve, maximum finded, solution of the estimation)
!! INPUTS
!!    ic        : index of the trait
!!    chr       : index of the chromosome
!!   est_moy    : True if the mean is given
!!   est_var    : True if the variance is given
!!   xlrmax     : LRT maximum
!!   nbfem      : number of female estimable
!!   nbniv      : number of level fixed effect
!!   nbef       : number of fixed effect
!!   nbco       : number of covariate
!!   ntlevp     : number of level for qtl sire interaction with fixed effect
!!   ntlevm     : number of level for qtl dam interaction with fixed effect
!!   par0       : solution (bestim) under H0
!!   par1       : solution (bestim) under H1
!!   ntniv      : number of level in the incidence matrix
!!   vecsol0    : logical array : estimability of level in the incidence matrix under H0
!!   vecsol1    : logical array : estimability of level in the incidence matrix under H1
!!   precis0    : precision information about each level under H0
!!   precis1    : precision information about each level under H1
!!   ordo       : LRT curve under H1
!!   npo        : number of position tested
!!   nposx      : the maximum position finded
!! printRiskFactor : print cox case with Risk factor
!!   qtl        : qtl effect are estimated
!!  hsire       : haplotype effect are estimated
!!  hdam        : haplotype effect with consideration of female estimable are estimated
!! NOTES
!!   En fonction du profil de stat de test trace la courbe correspondante
!!   Sous programme appele par analyse
!!   Version 1
!!   Date 1er Mars 2010. JM Elsen
!! SOURCE
       subroutine courbe_lin_ldla(dataset,spt,shp,chr,est_moy,est_var,xlrmax,ic, &
                            nbfem,nbniv,nbef,nbco,ntlevp,ntlevm,par0,par1,&
                            ntniv,vecsol0,precis0,vecsol1,precis1,        &
                            ordo,npo,nposx,printRiskFactor,qtl,hsire,hdam)
      type(QTLMAP_DATASET)             ,intent(in)     :: dataset
      type(PDD_BUILD)                  ,intent(in)     :: spt
      type(HAPLOTYPE_POSITION_BUILD)   ,intent(inout)     :: shp
      logical , intent(in)                        :: est_moy,est_var
      integer                       ,intent(in)   :: ic
      real (kind=dp)                ,intent(in)   :: xlrmax
      integer                       ,intent(in)   :: chr,nbfem
      integer                       ,intent(in)   :: nbniv
      integer                       ,intent(in)   :: nbef
      integer                       ,intent(in)   :: nbco
      integer                       ,intent(in)   :: ntlevp
      integer                       ,intent(in)   :: ntlevm
      real (kind=dp),dimension(:),intent(in)   :: par0
      real (kind=dp),dimension(:),intent(in)   :: par1
       integer                       ,intent(in)   :: ntniv
      real (kind=dp),dimension(:),intent(in)   :: precis0
      logical       ,dimension(:),intent(in)   :: vecsol0
      real (kind=dp),dimension(:),intent(in)   :: precis1
      logical       ,dimension(:),intent(in)   :: vecsol1
      integer                       ,intent(in)    :: npo,nposx
      real (kind=dp)  ,intent(in) ,dimension(npo) :: ordo
      logical , intent(in) , optional             :: printRiskFactor
      real (kind=dp)                              :: par_refm,  par_refp
      real (kind=dp), dimension(nbef)             :: par_refef
      logical,intent(in)                          :: qtl,hsire,hdam

! ************* deb chgt *************
      integer :: n_est_par,i_haplo,lll
! ************* fin chgt *************

      logical :: estRiskFactor=.false.
      character*40  chrom
!
! Pour la courbe
      real    :: x_p(npo), y_p(npo)
      integer :: effp(dataset%genea%np)
      integer :: efft,ipos,ip,nm1,nm2,jm,ipar,indest,ntot,km
      integer :: ilev,ief,iniveau,ico,lp,ilp,lm,ilm,i
      logical :: impfem
      real (kind=dp),dimension(size(par0)) :: par0_t
      real (kind=dp),dimension(size(par1)) :: par1_t
      real(kind=dp) :: dxmax
      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa
      type(DATAMODEL_BASE) , pointer :: dpm

      dpm => dataset%phenoModel
      dg => dataset%genea
      dpa => dataset%phenoAnimal

       if ( present(printRiskFactor) ) then
          estRiskFactor=printRiskFactor
      end if
      dxmax=dataset%map%absi(chr,nposx)

!  comptages
!
      efft=0
      impfem=.false.
      do ip=1,dg%np
        effp(ip)=0
        nm1=dg%nmp(ip)+1
        nm2=dg%nmp(ip+1)
        do jm=nm1,nm2
          if(dpa%estime(ic,jm))impfem=.true.
          efft=efft+eff(jm)
          effp(ip)=effp(ip)+eff(jm)
        end do
      end do

!      write(nficout,3003)
! 3003 format(//1x,'Maximum likelihood ratio test :'/)
!      write(nficout,3000)dxmax*1.d2,xlrmax
! 3000 format(1x,'The maximum is reached at position ',f7.2,' cM, with value ',f8.3/)
!
!  impression des r�sultats sous H0
!
    write(nficout,600)
600 format(//80('*'),/,'Estimation of parameters under H0'//)

    par0_t=par0
! ************* deb chgt ************
    if (est_var) then
!  remise � l'�chelle des variables
!
      ! ipar-> 1 to np : std
      !            np+1: General mean (les autres analyses on un mean / pere)
      !
      par0_t=par0*dpm%sigt(ic)
    end if
! ************* fin chgt ************

    indest=0
    if (est_var) then
      indest=dg%np
      if (est_moy.and.vecsol0(1)) par0_t(dg%np+1)=par0(dg%np+1)+dpm%xmut(ic)

      write(nficout,601)
  601 format('Within sire standard deviation',/)
      do ip = 1,dg%np
        write(nficout,602)trim(dg%pere(ip)),par0_t(ip)
  602   format(' sire ',a, '  s.d. :',f10.3)
      end do
    endif

    if (.not. estRiskFactor) then
      write(nficout,603)
  603 format(//,'  parameter    ','        estimable ?    value     ','precision'/)
      if (est_moy.and.vecsol0(1)) then
        indest=indest+1
        write(nficout,604) par0_t(indest),precis0(1)
  604   format('General mean               yes ',2f10.3,1x)
      else
        write(nficout,6041)
 6041   format('General mean                no ')
      end if
    else
      write(nficout,6030)
  6030 format(//,'  parameter    ','        estimable ?    risk factor     '/)

    endif

      write(nficout,606)
  606 format(/,'Sire polygenic effects')
      par_refp=par0_t(indest+1)
      do ip =1,dg%np
         i=ip
         if (est_moy)  i=ip+1
         if (vecsol0(i)) then
           indest=indest+1
           if (.not. estRiskFactor) then
             write(nficout,607)trim(dg%pere(ip)), par0_t(indest),precis0(i)
  607        format(' sire ',a, 15x,'yes ',2f10.3,1x)
           else
             write(nficout,6070)trim(dg%pere(ip)), dexp(par0_t(indest)-par_refp)
  6070       format(' sire ',a, 15x,'yes ',f10.3)
           endif
        else
          write(nficout,608)trim(dg%pere(ip))
  608     format(' sire ',a, 15x,'no ')
        end if

      end do

      ntot=dg%np
      par_refm=par0_t(indest+1)
      if (est_moy) ntot=dg%np+1
      if(impfem) write(nficout,609)
  609 format(/,'Dam polygenic effects')
      km=0
      do jm=1,dg%nm
        if (dpa%estime(ic,jm))then
          km=km+1
          if(vecsol0(ntot+km)) then
            indest=indest+1
      if (.not. estRiskFactor) then
            write(nficout,610)trim(dg%mere(jm)), par0_t(indest),precis0(ntot+km)
  610 format(' dam  ',a, 15x,'yes ',2f10.3,1x)
      else
             write(nficout,6100)trim(dg%mere(jm)), dexp(par0_t(indest)-par_refm)
      endif
  6100 format(' dam  ',a, 15x,'yes ',f10.3)
          else
            write(nficout,611)trim(dg%mere(jm))
          end if
  611 format(' dam  ',a, 15x,'no ')
        end if
      end do

      ntot=ntot+nbfem
      if(nbniv.ne.0)then
        ilev=0
        write(nficout,612)
  612 format(/,'fixed effects')
        do ief=1,nbef
          par_refef(ief)=par0_t(indest+1)
          do iniveau=1,dpm%nlev(dpm%modele(ic,3+ief))
          ilev=ilev+1
          if (vecsol0(ntot+ilev)) then
            indest=indest+1
        if (.not. estRiskFactor) then
               write(nficout,613)trim(dpm%namefix(dpm%modele(ic,3+ief))),iniveau,par0_t(indest),precis0(ntot+ilev)
  613 format(a15,' level',i3,'  yes ',2f10.3,1x)
            else
           write(nficout,6130)trim(dpm%namefix(dpm%modele(ic,3+ief))),iniveau,dexp(par0_t(indest)-par_refef(ief))
  6130 format(a15,' level',i3,'  yes ',f10.3)
        endif
          else
            write(nficout,614)trim(dpm%namefix(dpm%modele(ic,3+ief))),iniveau
  614 format(a15,' level',i3,'  no ')
          end if
        end do
        end do
      end if

      ntot=ntot+nbniv
      if(nbco.ne.0)then
        write(nficout,615)
  615 format(/,'covariables')
        do ico=1,nbco
          if (vecsol0(ntot+ico)) then
            indest=indest+1
        if (.not. estRiskFactor) then
               write(nficout,616)trim(dpm%namecov(dpm%modele(ic,3+nbef+ico))),par0_t(indest),precis0(ntot+ico)
  616          format(a15,12x,'yes ',2f10.3,1x)
            else
               write(nficout,6160)trim(dpm%namecov(dpm%modele(ic,3+nbef+ico))),dexp(par0_t(indest))
  6160         format(a15,12x,'yes ',f10.3)
             endif
          else
            write(nficout,617)trim(dpm%namecov(dpm%modele(ic,3+nbef+ico)))
  617 format(a15,12x,'no ')
          end if
        end do
      end if
!
!  impression des resultats sous H1
      write(nficout,618)
  618 format(//80('*'),/,'Estimation of parameters under H1'//)

      par1_t=par1
!
!  remise � l'�chelle des variables
!
! ************* deb chgt *************
      n_est_par=0
      if (est_moy)  n_est_par=1
      ntot=n_est_par

      if(est_var)then
! ************* fin chgt *************
        par1_t=par1*dpm%sigt(ic)
      end if

      indest=0
      if (est_var) then
        indest=dg%np
        if(est_moy.and.vecsol1(1))par1_t(dg%np+1)=par1(dg%np+1)+dpm%xmut(ic)

        write(nficout,619)
  619   format('Within sire standard deviation ',/)
        do ip = 1,dg%np
          write(nficout,620)trim(dg%pere(ip)),par1_t(ip)
  620     format(' sire ',a, '  s.d. :',f10.3)
        end do
      endif



      if (.not. estRiskFactor) then
        write(nficout,621)
  621   format(//,'  parameter    ','        estimable ?    value     ','precision'/)
        if (est_moy.and.vecsol1(1)) then
          indest=indest+1
          write(nficout,604) par1_t(indest),precis1(1)
        else
          write(nficout,6041)
        end if
      else
        write(nficout,6210)
  6210   format(//,'  parameter    ','        estimable ?    risk factor     '/)
      end if

        if( hsire .or. hdam ) then

        write(nficout,640)
  640   format(/,'Haplotypes effects')

!        do i_haplo=1,nb_max_haplo(nposx)
 !       do i_haplo=1,nb_max_haplo-1
        do i_haplo=1,shp%nb_haplo_reduit
!
!  TRIM non disponible au CSIRO
!
          if(vecsol1(n_est_par+i_haplo)) then
            indest=indest+1
            if (.not. estRiskFactor)then
!              write(nficout,641)name_haplo_reduit(i_haplo+1)(:longhap),count_haplo(i_haplo+1),&
  !            write(nficout,641)name_haplo_reduit(i_haplo),pb_haplo_reduit(i_haplo),&
  !               par1_t(indest),precis1(n_est_par+i_haplo)
 !             print*,name_haplo_reduit(i_haplo),pb_haplo_reduit(i_haplo),&
 !                par1_t(indest),precis1(n_est_par+i_haplo)
              write(nficout,641)shp%name_haplo(i_haplo)(:dataset%params%longhap),shp%count_haplo(i_haplo),&
                 par1_t(indest),precis1(n_est_par+i_haplo)
! 641          format(' haplotype ',a, 14x,' yes ',2f10.3,9x)
 641          format(' haplotype ',a, ' freq = ',f4.2,2x,' yes ',2f10.3,9x)
            else
!              write(nficout,642)name_haplo(nposx,i_haplo)(:longhap), count_haplo(nposx,i_haplo),&
!              write(nficout,642)name_haplo_reduit(i_haplo+1)(:longhap), count_haplo(i_haplo+1),&
              write(nficout,642)shp%name_haplo_reduit(i_haplo),shp%pb_haplo_reduit_r(i_haplo),&
                        dexp(par1_t(indest))
 642          format(' haplotype ',a, ' freq = ',f4.2,2x,' yes ',f10.3,9x)
            endif
          else
 !          write(nficout,643)name_haplo(nposx,i_haplo)(:longhap),count_haplo(nposx,i_haplo)
 !           write(nficout,643)name_haplo_reduit(i_haplo+1)(:longhap),count_haplo(i_haplo+1)
              write(nficout,643)shp%name_haplo_reduit(i_haplo),shp%pb_haplo_reduit_r(i_haplo)
 643        format(' haplotype ',a, ' freq = ',f4.2,2x,' no ')
          end if
        end do

 !       n_est_par=n_est_par+nb_max_haplo(nposx)
!        n_est_par=n_est_par+nb_max_haplo-1
         n_est_par=n_est_par+shp%nb_haplo_reduit

      end if ! option

      if(qtl) then
      write(nficout,624)
  624 format(/,'Sire QTL effects')

      if(ntlevp.eq.1) then
        write(nficout,6241)
 6241   format(59x,'allelic origin'/)
        do ip =1,dg%np

          i=n_est_par+ip

  ! ****** POUR JM ********
          call shp%liste_haplo(chr,dxmax,nposx,hsire,hdam)
          call shp%sort_haplo(chr,nposx,hsire,hdam)
!          print *,'num_haplo_pere:',num_haplo_pere(ip,1),num_haplo_pere(ip,2)
          chrom=shp%name_haplo(shp%num_haplo_pere(ip,1,1))(:dataset%params%longhap)//' '//&
                shp%name_haplo(shp%num_haplo_pere(ip,2,1))(:dataset%params%longhap)

          if (vecsol1(i)) then
            indest=indest+1
            if(spt%phasp(chr,ip)) then
              if (.not. estRiskFactor)then
                 write(nficout,6251)trim(dg%pere(ip)), par1_t(indest),precis1(i),chrom
 6251            format(' sire ',a, 14x,' yes ',2f10.3,9x,'known  ',a)
              else
                 write(nficout,62510)trim(dg%pere(ip)), dexp(par1_t(indest)),chrom
 62510           format(' sire ',a, 14x,' yes ',f10.3,9x,'known  ',a)
              endif
            else
              if (.not. estRiskFactor)then
                 write(nficout,6252)trim(dg%pere(ip)), par1_t(indest),precis1(i),chrom
 6252            format(' sire ',a, 14x,' yes ',2f10.3,9x,'unknown  ',a)
              else
                 write(nficout,62520)trim(dg%pere(ip)), dexp(par1_t(indest)),chrom
 62520           format(' sire ',a, 14x,' yes ',f10.3,9x,'unknown  ',a)
              endif
            end if
          else
            write(nficout,626)trim(dg%pere(ip)),chrom
  626 format(' sire ',a, 14x,' no  ',a)
          end if
        end do
      end if

      if(ntlevp.ne.1) then
       write(nficout,6242)
 6242   format(59x,'allelic origin'/)
      do ip =1,dg%np
        do lp = 1,ntlevp
! ************* deb chgt *************
      ilp=ntlevp*(ip-1)+lp-1+ n_est_par
!       ilp=ntlevp*(ip-1)+lp-1
!      if (est_moy) ilp=ntlevp*(ip-1)+lp
! ************* fin chgt *************
          if (vecsol1(1+ilp)) then
            indest=indest+1
            if (.not. estRiskFactor)then
              if(spt%phasp(chr,ip)) then
                 write(nficout,6271)trim(dg%pere(ip)), lp,par1_t(indest),precis1(1+ilp)
 6271            format(' sire ',a, 'level ',i3,5x,' yes ',2f10.3,9x,'known')
              else
                 write(nficout,6272)trim(dg%pere(ip)), lp,par1_t(indest),precis1(1+ilp)
 6272            format(' sire ',a, 'level ',i3,5x,' yes ',2f10.3,9x,'unknown')
              end if
         else
               if(spt%phasp(chr,ip)) then
                 write(nficout,62710)trim(dg%pere(ip)), lp,dexp(par1_t(indest))
 62710           format(' sire ',a, 'level ',i3,5x,' yes ',f10.3,9x,'known')
              else
                 write(nficout,62720)trim(dg%pere(ip)), lp,dexp(par1_t(indest))
 62720           format(' sire ',a, 'level ',i3,5x,' yes ',f10.3,9x,'unknown')
              end if
         endif
         else
            write(nficout,628)trim(dg%pere(ip)),lp
  628 format(' sire ',a, 'level ',i3,5x,' no ')
          end if
        end do
      end do
      end if

! ************* deb chgt *************
      ntot=ntlevp*dg%np+ n_est_par
!      ntot=ntlevp*np
!      if (est_moy) ntot=1+ntlevp*np
! ************* fin chgt *************
      if(impfem) write(nficout,629)
  629 format(/,'Dam QTL effects')

      if(ntlevm.eq.1) then

      if(impfem) write(nficout,6241)
      do jm=1,dg%nm
        if (dpa%estime(ic,jm))then
          km=dpa%iam(ic,dg%repfem(jm))
          if(vecsol1(ntot+km)) then
            indest=indest+1

        if (.not. estRiskFactor)then
              if(spt%phasm(chr,jm)) then
                 write(nficout,6301)trim(dg%mere(jm)), par1_t(indest),precis1(ntot+km)
 6301            format(' dam  ',a, 14x,' yes ',2f10.3,9x,'known')
              else
                 write(nficout,6302)trim(dg%mere(jm)), par1_t(indest),precis1(ntot+km)
 6302            format(' dam  ',a, 14x,' yes ',2f10.3,9x,'unknown')
              end if
        else
              if(spt%phasm(chr,jm)) then
                 write(nficout,63010)trim(dg%mere(jm)), dexp(par1_t(indest))
 63010            format(' dam  ',a, 14x,' yes ',f10.3,9x,'known')
              else
                 write(nficout,63020)trim(dg%mere(jm)), dexp(par1_t(indest))
 63020            format(' dam  ',a, 14x,' yes ',f10.3,9x,'unknown')
              end if
        endif

          end if
        end if
      end do
      end if

      if(ntlevm.ne.1) then

      if(impfem) write(nficout,6242)
          do jm=1,nbfem
            if(dpa%estime(ic,jm))then
              km=dpa%iam(ic,dg%repfem(jm))
              do lm=1,ntlevm
                ilm=ntlevm*(km-1)+lm
                if(vecsol1(ntot+ilm))then
                  indest=indest+1

          if (.not. estRiskFactor)then
                    if(spt%phasm(chr,jm)) then
                       write(nficout,6321)trim(dg%mere(jm)),lm,par1_t(indest),precis1(ntot+ilm)
 6321                  format(' dam  ',a, 'level ',i3,6x,' yes ',2f10.3,9x,'known')
                    else
                       write(nficout,6322)trim(dg%mere(jm)),lm,par1_t(indest),precis1(ntot+ilm)
 6322                  format(' dam  ',a, 'level ',i3,6x,' yes ',2f10.3,9x,'unknown')
                    end if
          else
                    if(spt%phasm(chr,jm)) then
                       write(nficout,63210)trim(dg%mere(jm)),lm,dexp(par1_t(indest))
 63210                  format(' dam  ',a, 'level ',i3,6x,' yes ',f10.3,9x,'known')
                    else
                       write(nficout,63220)trim(dg%mere(jm)),lm,dexp(par1_t(indest))
 63220                  format(' dam  ',a, 'level ',i3,6x,' yes ',f10.3,9x,'unknown')
                    end if
          endif

                end if
              end do
            end if
          end do

        end if

      write(nficout,*)
      write(nficout,*)' NOTE: known allelic origin means QTL effect =  maternal - paternal allele effects'

      ntot=ntot+dpa%namest(ic)*ntlevm
      end if ! option


      write(nficout,606)
       par_refp=par1_t(indest+1)
      do ip =1,dg%np
        if (vecsol1(ntot+ip)) then
          indest=indest+1
          if (.not. estRiskFactor)then
             write(nficout,607)trim(dg%pere(ip)),par1_t(indest),precis1(ntot+ip)
          else
             write(nficout,6070)trim(dg%pere(ip)),dexp(par1_t(indest)-par_refp)
      endif
    else
          write(nficout,608)trim(dg%pere(ip))
        end if
      end do

      ntot=ntot+dg%np
      write(nficout,609)
      km=0
      par_refm=par1_t(indest+1)
      do jm=1,dg%nm
        if (dpa%estime(ic,jm)) then
          km=km+1
          if(vecsol1(ntot+km)) then
            indest=indest+1
            if (.not. estRiskFactor)then
              write(nficout,610)trim(dg%mere(jm)),par1_t(indest),precis1(ntot+km)
            else
              write(nficout,6100)trim(dg%mere(jm)),dexp(par1_t(indest)-par_refm)
        endif
         else
          write(nficout,611)trim(dg%mere(jm))
         end if
       end if
      end do

      ntot=ntot+nbfem
      if(nbniv.ne.0)then
        ilev=0
        write(nficout,612)

        do ief=1,nbef
      par_refef(ief)=par0_t(indest+1)
      do iniveau=1,dpm%nlev(dpm%modele(ic,3+ief))
          ilev=ilev+1
          if (vecsol1(ntot+ilev)) then
            indest=indest+1
            if (.not. estRiskFactor)then
              write(nficout,613)trim(dpm%namefix(dpm%modele(ic,3+ief))),iniveau,par1_t(indest),precis1(ntot+ilev)
        else
              write(nficout,6130)trim(dpm%namefix(dpm%modele(ic,3+ief))),iniveau,dexp(par1_t(indest)-par_refef(ief))
        endif

          else
            write(nficout,614)trim(dpm%namefix(dpm%modele(ic,3+ief))),iniveau
          end if
          end do
        end do
      end if

      ntot=ntot+nbniv
      if(nbco.ne.0)then
        write(nficout,615)
        do ico=1,nbco
          if (vecsol1(ntot+ico)) then
            indest=indest+1
            if (.not. estRiskFactor)then
              write(nficout,616)trim(dpm%namecov(dpm%modele(ic,3+nbef+ico))),par1_t(indest),precis1(ntot+ico)
            else
          write(nficout,6160)trim(dpm%namecov(dpm%modele(ic,3+nbef+ico))),dexp(par1_t(indest))
        endif
          else
            write(nficout,617)trim(dpm%namecov(dpm%modele(ic,3+nbef+ico)))
          end if
        end do
      end if


      return
      end subroutine courbe_lin_ldla
!!***

      !----------------------------------------------------------------------------------------
!
!   LONGUEUR  remplace la defaillance de TRIM sur l'ordi du CSIRO
!
!----------------------------------------------------------------------------------------

      integer function longueur(a)
      character*40 a
      integer i
       longueur=0
       do i=1,len(a)
         if(a(i:i) /=' ') longueur=longueur+1
       end do
       return
       end function longueur


     ! Dev 17/05/2013
     ! La numerotation "phase 1/phase 2" ne suppose pas de DL avec les alleles au QTL
     ! Cette routine est utilisée pour estimer des effets QTL communs aux differentes familles,
     ! il faut etre sur que le meme allele au QTL est porté par la phase avec la même designation (1 ou 2) dans les differentes familles
     ! => utilse pour le modele biallelique/fonction discriminante du multicaractere
     !
     ! Proposition : la phase 1 sera celle portant l'allèle diminuant la valeur du caractere (Q), la phase 2 portera l'allèle augmentant la valeur du caractere (dit q)
     !
     ! en sortie : Les PDD sont modifiées dans une famille (les proba de reception du chr 1 devient chr 2 et vice versa)
     !
     ! Algo :
     !  1) On corrige Y avec l'estimation de la solution de la position precedente
     !  2) Calcul des moyennes Mu1 et Mu2 des Ycorr porté par 1 ou 2
     !  3) si Mu1 < Mu2 => on inverse les PDD

     subroutine set_phase_pdd(dataset,spt,chr,ipos,ic,nd)!,ntnivmax,ntniv,matrix_inc,vecsol,beta)
       type(QTLMAP_DATASET)             ,intent(in)       :: dataset   ! je jeux de donnée
       type(PDD_BUILD)                  ,intent(inout)    :: spt       ! les probas de transmissions
       integer                          ,intent(in)       :: chr       ! le GL ou se trouve la position à testé
       integer                          ,intent(in)       :: ipos      ! La position
       integer                          ,intent(in)       :: ic        ! caractere etudié
       integer                          ,intent(in)       :: nd        ! dimension ligne de la matrice d incidence
     !  integer                          ,intent(in)       :: ntnivmax  ! dimension colonne de la matrice d'incidence
     !  integer                          ,intent(in)       :: ntniv     ! dimension colonne de la matrice d'incidence
     !  real(kind=dp) , dimension(nd,ntnivmax),intent(in)     :: matrix_inc ! description des echantillons
     !  logical          , dimension(ntnivmax),intent(in)     :: vecsol     ! corresponding vector for column incidence matrix / solution
     !  real(kind=dp)    , dimension(ntnivmax),intent(in)     :: beta       ! solution pour la correction de Y . devrait etre de la position precedente

       real(kind=dp) :: ycorr,buffpdd
       real(kind=dp) :: buf,sumYcorr,Mu1p,Mu2p,Mu1m,Mu2m
       integer :: ip,jm,kd,ikd,iniv,ig,kkd,iniv2
       type(GENEALOGY_BASE) , pointer :: dg
       dg => dataset%genea

       ikd=0
       do ip=1,dg%np
         ! Manage half-sib
         Mu1p = 0.d0
         Mu2p = 0.d0
         do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
          do ig=spt%ngenom(chr,jm)+1,spt%ngenom(chr,jm+1)
            do kd=spt%ngend(chr,ig)+1,spt%ngend(chr,ig+1)
             kkd=spt%ndesc(chr,kd)
             if ( dataset%phenoAnimal%presentc(ic,kd)) then
              ikd=ikd+1
              ycorr=0.d0
              iniv2=0
!              do iniv=1,ntniv
!               if ( vecsol(iniv)) then
!                iniv2=iniv2+1
!                ycorr = ycorr + matrix_inc(ikd,iniv)*beta(iniv2)
!               end if
!              end do
              ycorr = dataset%phenoAnimal%y(ic,kkd) - ycorr
              Mu1p = Mu1p + ycorr*(spt%pdd(chr,kd,1,ipos)+spt%pdd(chr,kd,2,ipos))
              Mu2p = Mu2p + ycorr*(spt%pdd(chr,kd,3,ipos)+spt%pdd(chr,kd,4,ipos))
             end if
            end do ! kd
           end do ! ig
         end do ! jm

         if ( Mu1p < Mu2p ) then
            print *,"inverse pdd....ip=",ip,' pos=',ipos," Mu1=",Mu1p," Mu2=",Mu2p
            do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
             do ig=spt%ngenom(chr,jm)+1,spt%ngenom(chr,jm+1)
              do kd=spt%ngend(chr,ig)+1,spt%ngend(chr,ig+1)

               buffpdd = spt%pdd(chr,kd,1,ipos)
               spt%pdd(chr,kd,1,ipos) = spt%pdd(chr,kd,3,ipos)
               spt%pdd(chr,kd,3,ipos) = buffpdd

               buffpdd = spt%pdd(chr,kd,2,ipos)
               spt%pdd(chr,kd,2,ipos) = spt%pdd(chr,kd,4,ipos)
               spt%pdd(chr,kd,4,ipos) = buffpdd

              end do
             end do
            end do
         end if


         !Manage full sib
         do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
          if ( dataset%phenoAnimal%estime(ic,jm) ) then
           Mu1m = 0.d0
           Mu2m = 0.d0
           do ig=spt%ngenom(chr,jm)+1,spt%ngenom(chr,jm+1)
            do kd=spt%ngend(chr,ig)+1,spt%ngend(chr,ig+1)
             kkd=spt%ndesc(chr,kd)
             if ( dataset%phenoAnimal%presentc(ic,kd)) then
              ikd=ikd+1
              ycorr=0.d0
              iniv2=0
!              do iniv=1,ntniv
!               if ( vecsol(iniv)) then
!                iniv2=iniv2+1
!                ycorr = ycorr + matrix_inc(ikd,iniv)*beta(iniv2)
!               end if
!              end do
              Mu1m = Mu1m + ycorr*(spt%pdd(chr,kd,1,ipos)+spt%pdd(chr,kd,3,ipos))
              Mu2m = Mu2m + ycorr*(spt%pdd(chr,kd,2,ipos)+spt%pdd(chr,kd,4,ipos))
             end if
            end do ! kd
           end do ! ig

           if ( Mu1m < Mu2m ) then
            print *,"inverse pdd....jm=",jm,' pos=',ipos," Mu1=",Mu1m," Mu2=",Mu2m
             do ig=spt%ngenom(chr,jm)+1,spt%ngenom(chr,jm+1)
              do kd=spt%ngend(chr,ig)+1,spt%ngend(chr,ig+1)

               buffpdd = spt%pdd(chr,kd,1,ipos)
               spt%pdd(chr,kd,1,ipos) = spt%pdd(chr,kd,2,ipos)
               spt%pdd(chr,kd,2,ipos) = buffpdd

               buffpdd = spt%pdd(chr,kd,3,ipos)
               spt%pdd(chr,kd,3,ipos) = spt%pdd(chr,kd,4,ipos)
               spt%pdd(chr,kd,4,ipos) = buffpdd
             end do
            end do
           end if
          end if ! estime
          end do ! jm
       end do ! ip

     end subroutine set_phase_pdd

 end module m_qtlmap_analyse_gen
