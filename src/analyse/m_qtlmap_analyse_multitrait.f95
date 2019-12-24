!!****m* ANALYSE/m_qtlmap_analyse_multitrait
!!  NAME
!!    m_qtlmap_analyse_multitrait
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
module m_qtlmap_analyse_multitrait
   use m_qtlmap_types
   use m_qtlmap_log
   use m_qtlmap_math
   use m_qtlmap_optimization

   implicit none
   save

   type(GENEALOGY_BASE) , pointer :: p_dg
   type(PDD_BUILD)      , pointer :: p_spt
   type(PHENOTYPE_BASE) , pointer :: p_dpa
   type(DATAMODEL_BASE) , pointer :: p_dpm

!!****v* m_qtlmap_analyse_multitrait/std
!!  NAME
!!   std
!!  DESCRIPTION
!!   the standart deviation (residual) found under H1
!! DIMENSIONS
!!   ncar,np
!!
!!***
   real (kind=dp) ,dimension(:,:),allocatable,public   :: std
!!****v* m_qtlmap_analyse_multitrait/ap
!!  NAME
!!   ap
!!  DESCRIPTION
!!   The qtl sire effect found der H1
!! DIMENSIONS
!!   ncar,np
!!
!!***
   real (kind=dp) ,dimension(:,:),allocatable,public   :: ap
!!****v* m_qtlmap_analyse_multitrait/xmoyp
!!  NAME
!!   xmoyp
!!  DESCRIPTION
!!   The polygenic sire effect found under H1
!! DIMENSIONS
!!  ncar,np
!!
!!***
   real (kind=dp) ,dimension(:,:),allocatable,public   :: xmoyp
!!****v* m_qtlmap_analyse_multitrait/am
!!  NAME
!!   am
!!  DESCRIPTION
!!   The qtl dam effect found der H1
!! DIMENSIONS
!!   ncar,nm
!!
!!***
   real (kind=dp) ,dimension(:,:),allocatable,public   :: am
!!****v* m_qtlmap_analyse_multitrait/xmoym
!!  NAME
!!   xmoym
!!  DESCRIPTION
!!   The polygenic dam effect found under H1
!! DIMENSIONS
!!   nm
!!
!!***
   real (kind=dp) ,dimension(:,:),allocatable,public   :: xmoym


!!****v* m_qtlmap_analyse_multitrait/sompp11
!!  NAME
!!   sompp11
!!  DESCRIPTION
!!   initialization : loop before the minimization of the likelihood : opti_mcar_1qtl
!! DIMENSIONS
!!
!!***
   real (kind=dp)       ,dimension(:),allocatable,private   :: sompp11
!!****v* m_qtlmap_analyse_multitrait/sompm11
!!  NAME
!!   sompm11
!!  DESCRIPTION
!!   initialization : loop before the minimization of the likelihood : opti_mcar_1qtl
!! DIMENSIONS
!!
!!***
   real (kind=dp)       ,dimension(:),allocatable,private   :: sompm11
!!****v* m_qtlmap_analyse_multitrait/carpp11
!!  NAME
!!   carpp11
!!  DESCRIPTION
!!   initialization : loop before the minimization of the likelihood : opti_mcar_1qtl
!! DIMENSIONS
!!
!!***
   real (kind=dp)       ,dimension(:),allocatable,private   :: carpp11
!!****v* m_qtlmap_analyse_multitrait/carpm11
!!  NAME
!!   carpm11
!!  DESCRIPTION
!!   initialization : loop before the minimization of the likelihood : opti_mcar_1qtl
!! DIMENSIONS
!!
!!***
   real (kind=dp)       ,dimension(:),allocatable,private   :: carpm11
!!****v* m_qtlmap_analyse_multitrait/somppm11
!!  NAME
!!   somppm11
!!  DESCRIPTION
!!   initialization : loop before the minimization of the likelihood : opti_mcar_1qtl
!! DIMENSIONS
!!
!!***
   real (kind=dp)       ,dimension(:),allocatable,private   :: somppm11
!!****v* m_qtlmap_analyse_multitrait/sompmy11
!!  NAME
!!   sompmy11
!!  DESCRIPTION
!!   initialization : loop before the minimization of the likelihood : opti_mcar_1qtl
!! DIMENSIONS
!!
!!***
   real (kind=dp)     ,dimension(:,:),allocatable,private   :: sompmy11
!!****v* m_qtlmap_analyse_multitrait/somppy11
!!  NAME
!!   somppy11
!!  DESCRIPTION
!!   initialization : loop before the minimization of the likelihood : opti_mcar_1qtl
!! DIMENSIONS
!!
!!***
   real (kind=dp)     ,dimension(:,:),allocatable,private   :: somppy11
!!****v* m_qtlmap_analyse_multitrait/somppdf11
!!  NAME
!!   somppdf11
!!  DESCRIPTION
!!   initialization : loop before the minimization of the likelihood : opti_mcar_1qtl
!! DIMENSIONS
!!
!!***
   real (kind=dp)       ,dimension(:),allocatable,private   :: somppdf11
!!****v* m_qtlmap_analyse_multitrait/carppdf11
!!  NAME
!!   carppdf11
!!  DESCRIPTION
!!   initialization : loop before the minimization of the likelihood : opti_mcar_1qtl
!! DIMENSIONS
!!
!!***
   real (kind=dp)       ,dimension(:),allocatable,private   :: carppdf11
!!****v* m_qtlmap_analyse_multitrait/somppydf11
!!  NAME
!!   somppydf11
!!  DESCRIPTION
!!   initialization : loop before the minimization of the likelihood : opti_mcar_1qtl
!! DIMENSIONS
!!
!!***
   real (kind=dp)     ,dimension(:,:),allocatable,private   :: somppydf11
!!****v* m_qtlmap_analyse_multitrait/fp2
!!  NAME
!!   fp2
!!  DESCRIPTION
!!   the likelihood by half sib family
!! DIMENSIONS
!!   np
!!***
   real (kind=dp)       ,dimension(:),allocatable,private   :: fp2
!!****v* m_qtlmap_analyse_multitrait/fm2
!!  NAME
!!   fm2
!!  DESCRIPTION
!!   the likelihood by full sib family
!! DIMENSIONS
!!   nm
!!***
   real (kind=dp)       ,dimension(:),allocatable,private   :: fm2

   !$omp threadprivate (sompp11,sompm11,carpp11,carpm11,somppm11,sompmy11,somppy11,somppdf11,carppdf11,somppydf11,fp2,fm2)


!!****v* m_qtlmap_analyse_multitrait/somydf01
!!  NAME
!!   somydf01
!!  DESCRIPTION
!!   initialization : optinit_mcar
!! DIMENSIONS
!!
!!***
   real (kind=dp)       ,dimension(:,:),allocatable,private :: somydf01
!!****v* m_qtlmap_analyse_multitrait/carydf01
!!  NAME
!!   carydf01
!!  DESCRIPTION
!!   initialization : optinit_mcar
!! DIMENSIONS
!!
!!***
   real (kind=dp)       ,dimension(:,:),allocatable,private :: carydf01
!!****v* m_qtlmap_analyse_multitrait/somyydf01
!!  NAME
!!   somyydf01
!!  DESCRIPTION
!!   initialization : optinit_mcar
!! DIMENSIONS
!!
!!***
   real (kind=dp)     ,dimension(:,:,:),allocatable,private :: somyydf01
!!****v* m_qtlmap_analyse_multitrait/xmu01p
!!  NAME
!!   xmu01p
!!  DESCRIPTION
!!
!! DIMENSIONS
!!
!!***
   real (kind=dp)      ,dimension(:,:),allocatable,private  :: xmu01p
!!****v* m_qtlmap_analyse_multitrait/sig01
!!  NAME
!!   sig01
!!  DESCRIPTION
!!
!! DIMENSIONS
!!
!!***
   real (kind=dp)      ,dimension(:,:),allocatable,private  :: sig01
!!****v* m_qtlmap_analyse_multitrait/xmu01m
!!  NAME
!!   xmu01m
!!  DESCRIPTION
!!
!! DIMENSIONS
!!
!!***
   real (kind=dp)     ,dimension(:,:),allocatable,private   :: xmu01m
!!****v* m_qtlmap_analyse_multitrait/rhoi
!!  NAME
!!   rhoi
!!  DESCRIPTION
!!
!! DIMENSIONS
!!
!!***
   real (kind=dp)     ,dimension(:,:),allocatable,private   :: rhoi
!!****v* m_qtlmap_analyse_multitrait/fm21
!!  NAME
!!   fm21
!!  DESCRIPTION
!!   maximum value of likelihood by full sib family under H0
!! DIMENSIONS
!!   nm
!!***
   real (kind=dp)       ,dimension(:),allocatable,private   :: fm21
!!****v* m_qtlmap_analyse_multitrait/fp21
!!  NAME
!!   fp21
!!  DESCRIPTION
!!   maximum value of likelihood by half sib family under H0
!! DIMENSIONS
!!   np
!!***
   real (kind=dp)       ,dimension(:),allocatable,private   :: fp21
!!****v* m_qtlmap_analyse_multitrait/somy11
!!  NAME
!!   somy11
!!  DESCRIPTION
!!
!! DIMENSIONS
!!
!!***
   real (kind=dp)     ,dimension(:,:),allocatable,private   :: somy11
!!****v* m_qtlmap_analyse_multitrait/cary11
!!  NAME
!!   cary11
!!  DESCRIPTION
!!
!! DIMENSIONS
!!
!!***
   real (kind=dp)     ,dimension(:,:),allocatable,private   :: cary11
!!****v* m_qtlmap_analyse_multitrait/somyy11
!!  NAME
!!   somyy11
!!  DESCRIPTION
!!
!! DIMENSIONS
!!
!!***
   real (kind=dp)   ,dimension(:,:,:),allocatable,private   :: somyy11
!!****v* m_qtlmap_analyse_multitrait/xmu2p
!!  NAME
!!   xmu2p
!!  DESCRIPTION
!!
!! DIMENSIONS
!!
!!***
   real (kind=dp)     ,dimension(:,:),allocatable,private   :: xmu2p
!!****v* m_qtlmap_analyse_multitrait/sig2
!!  NAME
!!   sig2
!!  DESCRIPTION
!!
!! DIMENSIONS
!!
!!***
   real (kind=dp)     ,dimension(:,:),allocatable,private   :: sig2
!!****v* m_qtlmap_analyse_multitrait/xmu2m
!!  NAME
!!   xmu2m
!!  DESCRIPTION
!!
!! DIMENSIONS
!!
!!***
   real (kind=dp)     ,dimension(:,:),allocatable,private   :: xmu2m
!!****v* m_qtlmap_analyse_multitrait/rhoi2
!!  NAME
!!   rhoi2
!!  DESCRIPTION
!!
!! DIMENSIONS
!!
!!***
   real (kind=dp)     ,dimension(:,:),allocatable,private   :: rhoi2
!!****v* m_qtlmap_analyse_multitrait/f01
!!  NAME
!!   f01
!!  DESCRIPTION
!!
!! DIMENSIONS
!!
!!***
   real (kind=dp)                                ,private   :: f01
!!****v* m_qtlmap_analyse_multitrait/rhoi2_1
!!  NAME
!!   rhoi2_1
!!  DESCRIPTION
!!
!! DIMENSIONS
!!
!!***
   real (kind=dp)     ,dimension(:,:),allocatable,private   :: rhoi2_1
!!****v* m_qtlmap_analyse_multitrait/rhoi2_2
!!  NAME
!!   rhoi2_2
!!  DESCRIPTION
!!
!! DIMENSIONS
!!
!!***
   real (kind=dp)     ,dimension(:,:),allocatable,private   :: rhoi2_2
!!****v* m_qtlmap_analyse_multitrait/corcd
!!  NAME
!!   corcd
!!  DESCRIPTION
!!
!! DIMENSIONS
!!
!!***
   real (kind=dp)   ,dimension(:,:,:),allocatable,private   :: corcd
!!****v* m_qtlmap_analyse_multitrait/effdf
!!  NAME
!!   effdf
!!  DESCRIPTION
!!
!! DIMENSIONS
!!
!!***
   integer              ,dimension(:),pointer,public        :: effdf
!!****v* m_qtlmap_analyse_multitrait/estmum
!!  NAME
!!   estmum
!!  DESCRIPTION
!!
!! DIMENSIONS
!!
!!***
   integer      ,dimension(:),pointer, public               :: estmum
!!****v* m_qtlmap_analyse_multitrait/eff
!!  NAME
!!   eff
!!  DESCRIPTION
!!
!! DIMENSIONS
!!
!!***
   integer      ,dimension(:),pointer     ,public           :: eff
!!****v* m_qtlmap_analyse_multitrait/effp
!!  NAME
!!   effp
!!  DESCRIPTION
!!
!! DIMENSIONS
!!
!!***
   integer       ,dimension(:),allocatable,public           :: effp
!!****v* m_qtlmap_analyse_multitrait/current_chr
!!  NAME
!!   current_chr
!!  DESCRIPTION
!!
!! DIMENSIONS
!!
!!***
   integer              ,private                            :: current_chr

   public :: init_analyse_multitrait
   private :: init_sub_1qtl
   private :: end_sub_1qtl

   public :: optinit_mcar
   public :: opti_mcar_0qtl
   public :: opti_mcar_1qtl
   public :: set_solution_hypothesis0
   public :: set_solution_hypothesis1
   public :: end_analyse_multitrait


   contains
!!****f* m_qtlmap_analyse_multitrait/init_analyse_multitrait
!!  NAME
!!    init_analyse_multitrait
!!  DESCRIPTION
!!    allocation of buffer/solution arrays
!!  SOURCE
      subroutine init_analyse_multitrait(dataset)
         type(QTLMAP_DATASET)       ,intent(in) :: dataset
         integer              :: stat

         type(GENEALOGY_BASE) , pointer :: dg
         type(DATAMODEL_BASE) , pointer :: dpm

         dpm => dataset%phenoModel
         dg => dataset%genea

         allocate (std(dpm%ncar,dg%np),STAT=stat)
         call check_allocate(stat,'xmoyp [m_qtlmap_analyse_multitrait]')
         allocate (ap(dpm%ncar,dg%np),STAT=stat)
         call check_allocate(stat,'xmoyp [m_qtlmap_analyse_multitrait]')
         allocate (xmoyp(dpm%ncar,dg%np),STAT=stat)
         call check_allocate(stat,'xmoyp [m_qtlmap_analyse_multitrait]')
         allocate (am(dpm%ncar,dg%nm),STAT=stat)
         call check_allocate(stat,'xmoyp [m_qtlmap_analyse_multitrait]')
         allocate (xmoym(dpm%ncar,dg%nm),STAT=stat)
         call check_allocate(stat,'xmoyp [m_qtlmap_analyse_multitrait]')
         allocate (somydf01(dpm%ncar,dg%np),STAT=stat)
         call check_allocate(stat,'somydf01 [m_qtlmap_analyse_multitrait]')
         allocate (carydf01(dpm%ncar,dg%np),STAT=stat)
         call check_allocate(stat,'carydf01 [m_qtlmap_analyse_multitrait]')
         allocate (somyydf01(dpm%ncar,dpm%ncar,dg%np),STAT=stat)
         call check_allocate(stat,'somyydf01 [m_qtlmap_analyse_multitrait]')
         allocate (sig01(dpm%ncar,dg%np),STAT=stat)
         call check_allocate(stat,'sig01 [m_qtlmap_analyse_multitrait]')
         allocate (xmu01p(dpm%ncar,dg%np),STAT=stat)
         call check_allocate(stat,'xmu01p [m_qtlmap_analyse_multitrait]')
         allocate (xmu01m(dpm%ncar,dg%nm),STAT=stat)
         call check_allocate(stat,'xmu01m [m_qtlmap_analyse_multitrait]')
         allocate (rhoi(dpm%ncar,dpm%ncar),STAT=stat)
         call check_allocate(stat,'rhoi [m_qtlmap_analyse_multitrait]')

         allocate (fp21(dg%np),STAT=stat)
         call check_allocate(stat,'fp21 [m_qtlmap_analyse_multitrait]')
         allocate (fm21(dg%nm),STAT=stat)
         call check_allocate(stat,'fm21 [m_qtlmap_analyse_multitrait]')
         allocate (somy11(dpm%ncar,dg%nm),STAT=stat)
         call check_allocate(stat,'somy11 [m_qtlmap_analyse_multitrait]')
         allocate (cary11(dpm%ncar,dg%nm),STAT=stat)
         call check_allocate(stat,'cary11 [m_qtlmap_analyse_multitrait]')
         allocate (somyy11(dpm%ncar,dpm%ncar,dg%nm),STAT=stat)
         call check_allocate(stat,'somyy11 [m_qtlmap_analyse_multitrait]')
         allocate (xmu2p(dpm%ncar,dg%np),STAT=stat)
         call check_allocate(stat,'xmu2p [m_qtlmap_analyse_multitrait]')
         allocate (sig2(dpm%ncar,dg%np),STAT=stat)
         call check_allocate(stat,'sig2 [m_qtlmap_analyse_multitrait]')
         allocate (xmu2m(dpm%ncar,dg%nm),STAT=stat)
         call check_allocate(stat,'xmu2m [m_qtlmap_analyse_multitrait]')
         xmu2m=0.d0         
         allocate (rhoi2(dpm%ncar,dpm%ncar),STAT=stat)
         call check_allocate(stat,'rhoi2 [m_qtlmap_analyse_multitrait]')
         rhoi2=0.d0
         allocate (rhoi2_2(dpm%ncar,dpm%ncar))
         rhoi2_2=0.d0
         allocate (rhoi2_1(dpm%ncar,dpm%ncar))
         rhoi2_1=0.d0

         allocate (effdf(dg%np))
         allocate (estmum(dg%np))
         allocate (eff(dg%nm))
         allocate (effp(dg%np))

         allocate (corcd(dpm%ncar,dpm%ncar,dg%nd))
         corcd=1.d0

      end subroutine init_analyse_multitrait
!!***

!!****f* m_qtlmap_analyse_multitrait/init_sub_1qtl
!!  NAME
!!    init_sub_1qtl
!!  DESCRIPTION
!!    allocation of buffer/solution arrays
!!  SOURCE
      subroutine init_sub_1qtl(dataset,spt)
         type(QTLMAP_DATASET)       ,intent(in) :: dataset
         type(PDD_BUILD)            ,intent(in) :: spt
         integer :: stat,ng
         type(GENEALOGY_BASE) , pointer :: dg
         type(DATAMODEL_BASE) , pointer :: dpm

         dg => dataset%genea
         dpm => dataset%phenoModel

         ng = spt%get_maxnbgenotypedam(dataset)

         allocate (sompp11(ng),STAT=stat)
         call check_allocate(stat,'sompp11 [m_qtlmap_analyse_multitrait]')
         allocate (sompm11(ng),STAT=stat)
         call check_allocate(stat,'sompm11 [m_qtlmap_analyse_multitrait]')
         allocate (carpp11(ng),STAT=stat)
         call check_allocate(stat,'carpp11 [m_qtlmap_analyse_multitrait]')
         allocate (carpm11(ng),STAT=stat)
         call check_allocate(stat,'carpm11 [m_qtlmap_analyse_multitrait]')
         allocate (somppm11(ng),STAT=stat)
         call check_allocate(stat,'somppm11 [m_qtlmap_analyse_multitrait]')
         allocate (sompmy11(dpm%ncar,ng),STAT=stat)
         call check_allocate(stat,'sompmy11 [m_qtlmap_analyse_multitrait]')
         allocate (somppy11(dpm%ncar,ng),STAT=stat)
         call check_allocate(stat,'somppy11 [m_qtlmap_analyse_multitrait]')
         allocate (somppdf11(ng),STAT=stat)
         call check_allocate(stat,'somppdf11 [m_qtlmap_analyse_multitrait]')
         allocate (carppdf11(ng),STAT=stat)
         call check_allocate(stat,'carppdf11 [m_qtlmap_analyse_multitrait]')
         allocate (somppydf11(dpm%ncar,ng),STAT=stat)
         call check_allocate(stat,'somppydf11 [m_qtlmap_analyse_multitrait]')

         allocate (fp2(dg%np),STAT=stat)
         call check_allocate(stat,'fp2 [m_qtlmap_analyse_multitrait]')
         allocate (fm2(dg%nm),STAT=stat)
         call check_allocate(stat,'fm2 [m_qtlmap_analyse_multitrait]')

      end subroutine init_sub_1qtl
!!***

!!****f* m_qtlmap_analyse_multitrait/end_sub_1qtl
!!  NAME
!!    end_sub_1qtl
!!  DESCRIPTION
!!    deallocation of buffer/solution arrays
!!  SOURCE
      subroutine end_sub_1qtl
         deallocate (sompp11)
         deallocate (sompm11)
         deallocate (carpp11)
         deallocate (carpm11)
         deallocate (somppm11)
         deallocate (sompmy11)
         deallocate (somppy11)
         deallocate (somppdf11)
         deallocate (carppdf11)
         deallocate (somppydf11)
         deallocate (fp2)
         deallocate (fm2)

      end subroutine end_sub_1qtl
!!***

!!****f* m_qtlmap_analyse_multitrait/end_analyse_multitrait
!!  NAME
!!    end_analyse_multitrait
!!  DESCRIPTION
!!    deallocation of buffer/solution arrays
!!  SOURCE
      subroutine end_analyse_multitrait

          deallocate (somydf01)
          deallocate (carydf01)
          deallocate (somyydf01)
          deallocate (sig01)
          deallocate (xmu01p)
          deallocate (xmu01m)
          deallocate (rhoi)
          deallocate (fp21)
          deallocate (fm21)
          deallocate (somy11)
          deallocate (cary11)
          deallocate (somyy11)
          deallocate (xmu2p)
          deallocate (sig2)
          deallocate (xmu2m)
          deallocate (rhoi2)
          deallocate (std)
          deallocate (ap)
          deallocate (xmoyp)
          deallocate (am)
          deallocate (xmoym)
          deallocate (rhoi2_2)
          deallocate (rhoi2_1)

          deallocate (effdf)
          deallocate (estmum)
          deallocate (eff)
          deallocate (effp)

          deallocate (corcd)

      end subroutine end_analyse_multitrait
!!***

!!****f* m_qtlmap_analyse_multitrait/optinit_mcar
!!  NAME
!!    optinit_mcar
!!  DESCRIPTION
!!    Initialise les ecart types et moyennes intra-famille
!!  SOURCE
      subroutine optinit_mcar(dataset)
      type(QTLMAP_DATASET)       ,intent(in) :: dataset
      integer        :: efft(dataset%phenoModel%ncar),i
      integer        :: ic,nm1,nm2,jm,nest,ifem,imumest,jc,ip,nd1,nd2,kd
      real (kind=dp) , dimension(dataset%phenoModel%ncar) :: somt,xmut,sigt
      real (kind=dp) , dimension(dataset%phenoModel%ncar,dataset%phenoModel%ncar) :: cov

      real (kind=dp) , dimension(dataset%genea%nm):: somy,cary
      real (kind=dp) :: somyp

      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa
      type(DATAMODEL_BASE) , pointer :: dpm

      dpm => dataset%phenoModel
      dg => dataset%genea
      dpa => dataset%phenoAnimal

      cov=0.d0
      somyydf01=0.d0
      somyy11=0.d0      !! que moitie inf
!
      do ic=1,dpm%ncar


      effp=0
      imumest=0

      do ip=1,dg%np
        somyp=0.d0
        sig01(ic,ip)=0.d0
        effdf(ip)=0
        somydf01(ic,ip)=0.d0
        carydf01(ic,ip)=0.d0
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
            if(count(dpa%presentc(:,kd)) == dpm%ncar ) then
              eff(jm)=eff(jm)+1
              somy(jm)=somy(jm)+dpa%y(ic,kd)*dpa%cd(ic,kd)
              cary(jm)=cary(jm)+(dpa%y(ic,kd)*dpa%y(ic,kd))*dpa%cd(ic,kd)
            end if
          end do

          effp(ip)=effp(ip)+eff(jm)
          somyp=somyp+somy(jm)
          if( dpa%estime(ic,jm) ) then
            imumest=imumest+1
            estmum(ip)=estmum(ip)+1
            xmu01m(ic,imumest)=somy(jm)/dble(eff(jm))
          else
            effdf(ip)=effdf(ip)+eff(jm)
            somydf01(ic,ip)=somydf01(ic,ip)+somy(jm)
            carydf01(ic,ip)=carydf01(ic,ip)+cary(jm)
          end if
        end do
        if (effp(ip) == 0.d0) then
           call stop_application('Father ['//trim(dg%pere(ip))//'] has got no child with trait value')
        else
           xmu01p(ic,ip)=somyp/dble(effp(ip))
        end if
        do jm=nm1,nm2
          nd1=dg%ndm(jm)+1
          nd2=dg%ndm(jm+1)
          do kd=nd1,nd2
             if(count(dpa%presentc(:,kd)) == dpm%ncar ) then
               sig01(ic,ip)=sig01(ic,ip)+(dpa%y(ic,kd)-xmu01p(ic,ip))*(dpa%y(ic,kd)-xmu01p(ic,ip))
             end if
          end do
        end do
        if ( effp(ip) == 1) then
           call stop_application('Father ['//trim(dg%pere(ip))//'] has got only one child with trait value')
        else
           sig01(ic,ip)=dsqrt(sig01(ic,ip)/dble(effp(ip)-1))
        end if

        if(estmum(ip).gt.0) then
          estmum(ip)=estmum(ip)-1
          dpa%nmumest(ic)=dpa%nmumest(ic)-1
        end if
      end do

!  Initialisation des param�tres dans des vecteurs d�pendants du nombre de caracteres
       efft(ic)=0
       somt(ic)=0.d0
       xmut(ic)=0.d0
       sigt(ic)=0.d0
       imumest=0
       do ip=1,dg%np
         nm1=dg%nmp(ip)+1
         nm2=dg%nmp(ip+1)
         nest=0
         do jm=nm1,nm2

           efft(ic)=efft(ic)+eff(jm)
           somt(ic)=somt(ic)+somy(jm)

           somy11(ic,jm)=somy(jm)
           cary11(ic,jm)=cary(jm)

           if(dpa%estime(ic,jm)) then
             nest=nest+1
           if (nest.le.estmum(ip)) then
              imumest=imumest+1
           end if
           end if
         end do             !! fin jm
        end do              !! fin ip

        xmut(ic)=somt(ic)/dble(efft(ic))

        do ip=1,dg%np
          nm1=dg%nmp(ip)+1
          nm2=dg%nmp(ip+1)
          do jm=nm1,nm2
            nd1=dg%ndm(jm)+1
            nd2=dg%ndm(jm+1)
            do kd=nd1,nd2
             if(count(dpa%presentc(:,kd)) == dpm%ncar ) then
               sigt(ic)=sigt(ic)+(dpa%y(ic,kd)-dpm%xmut(ic))*(dpa%y(ic,kd)-dpm%xmut(ic))
              if (ic.gt.1) then
               do jc=1,ic-1
                if(count(dpa%presentc(:,kd)) == dpm%ncar ) then
                 somyy11(jc,ic,jm)=somyy11(jc,ic,jm)+(dpa%y(jc,kd)*dpa%y(ic,kd))         !! que moitie inf
                 somyy11(ic,jc,jm)=somyy11(jc,ic,jm)
                 cov(jc,ic)=cov(jc,ic)+(dpa%y(jc,kd)-dpm%xmut(jc))*(dpa%y(ic,kd)-dpm%xmut(ic))
                end if
               end do
              end if
             end if
            end do
            if(.not.dpa%estime(ic,jm)) then
             if (ic.gt.1) then
              do jc=1,ic-1
               somyydf01(jc,ic,ip)=somyydf01(jc,ic,ip)+somyy11(jc,ic,jm)
               somyydf01(ic,jc,ip)=somyydf01(jc,ic,ip)
              end do
             end if
            end if
          end do
        end do
       sigt(ic)=sigt(ic)/dble(efft(ic))
       sigt(ic)=dsqrt(sigt(ic))

        if(ic.gt.1) then
          do jc=1,ic-1
           if (efft(jc).ne.efft(ic)) then
            do i=1,dg%nd
              if (count(dpa%presentc(:,i))<dpm%ncar) then
                call log_mess("*Remove animal :"//trim(dg%animal(i)))
              end if
            end do
             call log_mess('Missing trait data for traits '//trim(dpm%carac(ic))//  &
           ' and '//trim(dpm%carac(jc))//', not possible to perform multiple trait analysis',ERROR_DEF)
            call stop_application(" ** ")
           end if
           cov(jc,ic)=cov(jc,ic)/dble(efft(ic))
           cov(ic,jc)=cov(jc,ic)
           rhoi(jc,ic)=cov(jc,ic)/(sigt(jc)*sigt(ic))
           rhoi(ic,jc)=rhoi(jc,ic)
          end do
        end if
      end do
                                    !! fin ic
      return
      end subroutine optinit_mcar
!!***

!!****f* m_qtlmap_analyse_multitrait/opti_mcar_0qtl
!!  NAME
!!    opti_mcar_0qtl
!!  DESCRIPTION
!!    Calcul de la statistique de test multicaractere le long du chromosome
!!  SOURCE
      subroutine opti_mcar_0qtl(dataset)
      type(QTLMAP_DATASET)       ,intent(in) :: dataset
! Divers
      real (kind=dp) , dimension(:), allocatable :: par,borni,borns
      integer iuser(1)
      real (kind=dp)user(1)

      integer :: ibound,npar,nrho,ip,jm,j,k,nest,nm1,nm2
      integer :: r1,imumest,i,nestim,ifail,ic
      real (kind=dp) ::somxmu
      logical  , dimension(:,:),pointer        :: filter_inc
      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa
      type(DATAMODEL_BASE) , pointer :: dpm

      dpm => dataset%phenoModel
      dg => dataset%genea
      dpa => dataset%phenoAnimal

      p_dg => dataset%genea
      p_dpm => dataset%phenoModel
      p_dpa => dataset%phenoAnimal

!
!****************************************************************************
!
! Parametres de maximisation
      nrho=dpm%ncar*(dpm%ncar-1)/2
      npar=(2*dpm%ncar*dg%np)+(dpm%ncar*dpa%nmumest(1))+nrho
      allocate (par(npar))
      allocate (borni(npar))
      allocate (borns(npar))
      allocate (filter_inc(dg%np,npar))

      ibound=0
      filter_inc=.false.
      nestim=0
      do ip=1,dg%np
       do j=1,dpm%ncar
        borni((j-1)*dg%np+ip)=SIG_MIN
        borns((j-1)*dg%np+ip)=SIG_MAX
        filter_inc(ip,(j-1)*dg%np+ip)=.true.
        borni(ip+dpm%ncar*dg%np+(j-1)*dg%np)=XMU_MIN
        borns(ip+dpm%ncar*dg%np+(j-1)*dg%np)=XMU_MAX
        filter_inc(ip,ip+dpm%ncar*dg%np+(j-1)*dg%np)=.true.
       end do

       nest=0
       do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
           if(dpa%estime(1,jm)) then
            nest=nest+1
            if(nest.le.estmum(ip)) then
              nestim=nestim+1
               do j=1,dpm%ncar
                filter_inc(ip,2*dpm%ncar*dg%np+(j-1)*dpa%nmumest(1)+nestim)=.true.
               end do
            end if
           end if
        end do

      end do
      do jm=1,dpa%nmumest(1)
        do j=1,dpm%ncar
          borni(2*dpm%ncar*dg%np+(j-1)*dpa%nmumest(1)+jm)=XMU_MIN
          borns(2*dpm%ncar*dg%np+(j-1)*dpa%nmumest(1)+jm)=XMU_MAX
        end do
      end do
      do j=1,nrho
        borni(2*dpm%ncar*dg%np+dpm%ncar*dpa%nmumest(1)+j)=DEFAULT_PARAM_MIN
        borns(2*dpm%ncar*dg%np+dpm%ncar*dpa%nmumest(1)+j)=DEFAULT_PARAM_MAX
        filter_inc(:,2*dpm%ncar*dg%np+dpm%ncar*dpa%nmumest(1)+j)=.true.
      end do
!
! Point de depart
      k=0
      imumest=0
      do i=1,dpm%ncar
        nestim=0
        do ip=1,dg%np
          par((i-1)*dg%np+ip)=sig01(i,ip)
          par(dpm%ncar*dg%np+(i-1)*dg%np+ip)=xmu01p(i,ip)
          nm1=dg%nmp(ip)+1
          nm2=dg%nmp(ip+1)
          nest=0
          somxmu=0.d0
          do jm=nm1,nm2
           if (dpa%estime(i,jm)) then
            nest=nest+1
            if (nest.le.estmum(ip)) then
	     nestim=nestim+1
             par(2*dpm%ncar*dg%np+(i-1)*dpa%nmumest(1)+nestim)=xmu01m(i,nestim)
	    end if
           end if
          end do
        end do
        if (i.gt.1) then
         do j=1,i-1
          k=k+1
          r1=rhoi(j,i)
          par(2*dpm%ncar*dg%np+dpm%ncar*dpa%nmumest(1)+k)=dlog((1.d0+r1)/(1.d0-r1))
         end do
        end if
      end do
!
! Optimisation de la vraisemblance
      ifail=1
!      call minimizing_funct(npar,ibound,funct_mcar_0qtl,borni,borns,par,f01,iuser,user,ifail)
      call minimizing_funct_family_sire(dataset,npar,ibound,funct_mcar_0qtl_family,&
        filter_inc,fp21,borni,borns,par,f01,iuser,user,ifail)

      k=0
      do i=2,dpm%ncar
        do j=1,i-1
          k=k+1
          rhoi2(j,i)=par(2*dpm%ncar*dg%np+dpm%ncar*dpa%nmumest(1)+k)
          !**ajout ofi pour affichage sous H0
          rhoi2_1(j,i)=(dexp(rhoi2(j,i))-1.d0)/(dexp(rhoi2(j,i))+1.d0)
          rhoi2_1(i,j)=rhoi2_1(j,i)
          !**fin ajout
          rhoi2(i,j)=rhoi2(j,i)
        end do
      end do

      do i=1,dpm%ncar
       do ip=1,dg%np
        sig2(i,ip)=par((i-1)*dg%np+ip)
        xmu2p(i,ip)=par(dpm%ncar*dg%np+(i-1)*dg%np+ip)
       end do
       do jm=1,dpa%nmumest(1)
         xmu2m(i,jm)=par(2*dpm%ncar*dg%np+(i-1)*dpa%nmumest(1)+jm)
       end do
      end do

      deallocate (par)
      deallocate (borni)
      deallocate (borns)
      deallocate (filter_inc)
      end subroutine opti_mcar_0qtl
!!***

!!****f* m_qtlmap_analyse_multitrait/funct_mcar_0qtl
!!  NAME
!!    funct_mcar_0qtl
!!  DESCRIPTION
!!    Calcul de la statistique de test multicaractere le long du chromosome
!!  SOURCE
    subroutine funct_mcar_0qtl(n,x,f,iuser,user)
      integer         , intent(in)                  :: n
      real (kind=dp)      ,dimension(n), intent(inout) :: x
      real (kind=dp)    , intent(inout)             :: f
      integer ,       dimension(1), intent(in)      :: iuser
      real (kind=dp)      ,dimension(1), intent(in) :: user
      integer :: ip

      do ip=1,p_dg%np
        call funct_mcar_0qtl_family(ip,n,x,fp21(ip) ,iuser,user)
      end do

      f = sum(fp21)
      return

   end subroutine funct_mcar_0qtl
!!***

!!****f* m_qtlmap_analyse_multitrait/funct_mcar_0qtl_family
!!  NAME
!!    funct_mcar_0qtl_family
!!  DESCRIPTION
!!    Calcul de la statistique de test multicaractere le long du chromosome
!!  SOURCE
    subroutine funct_mcar_0qtl_family(ip,n,x,f,iuser,user)
      integer         , intent(in)                  :: ip,n
      real (kind=dp)      ,dimension(n), intent(inout) :: x
      real (kind=dp)    , intent(inout)             :: f
      integer ,       dimension(1), intent(in)      :: iuser
      real (kind=dp)      ,dimension(1), intent(in) :: user
!
! Tableaux dimensionnes selon np, le nombre de peres
      real (kind=dp) ,dimension(p_dg%np) :: det

!
! Tableaux dimensionnes selon le nombre nc de caracteres
      real (kind=dp) ,dimension(p_dpm%ncar,p_dg%np) :: sigc,xmupc,xmup2c
      real (kind=dp) ,dimension(p_dpm%ncar,p_dg%nm) :: xmumc,xmum2c,xmupmc
      real (kind=dp) ,dimension(p_dpm%ncar,p_dpm%ncar,p_dg%np) ::covc,vc
      real (kind=dp) ,dimension(p_dpm%ncar+1,p_dpm%ncar) :: A
      real (kind=dp) ,dimension(p_dpm%ncar,p_dpm%ncar)   :: rhoc,Ab
      real (kind=dp) ,dimension(p_dpm%ncar) :: wrkspce,somxmu
!
! Divers
      integer ::i,k,j,nt,nestim,ifail,ic,jc,nest
      integer :: jm
      real (kind=dp) :: rh,zdf,vdf,zpf,vpf,determ
!***********************************************************************************
      f=0.d0
! R�cup�ration des corr�lations
      k=0
      do i=2,p_dpm%ncar
        do j=1,i-1
          k=k+1
          nt=2*p_dpm%ncar*p_dg%np+p_dpm%ncar*p_dpa%nmumest(1)+k
          rh=x(nt)
         ! print *,"rh:",rh
          rhoc(j,i)=((dexp(rh)-1.d0)/(dexp(rh)+1.d0))
          rhoc(i,j)=rhoc(j,i)
         ! print *,"rhoc ",i,j,rhoc(i,j)
        end do
      end do
!
! Calcul de l'inverse et du determinant de la matrice de var-cov
      nestim=0
!
! Initialisation des parametres matrice de variance cov
       do i=1,p_dpm%ncar
        sigc(i,ip)=x((i-1)*p_dg%np+ip)
        covc(i,i,ip)=sigc(i,ip)*sigc(i,ip)
        if(i.gt.1) then
          do j=1,i-1
            covc(i,j,ip)=rhoc(i,j)*sigc(i,ip)*sigc(j,ip)
            covc(j,i,ip)=covc(i,j,ip)
          end do
        end if
       end do
!
! Si 2 caracteres, calcul a la main
       if(p_dpm%ncar.eq.2) then
         det(ip)=covc(1,1,ip)*covc(2,2,ip)-covc(1,2,ip)*covc(2,1,ip)

         vc(1,1,ip)=covc(2,2,ip)/det(ip)
         vc(2,2,ip)=covc(1,1,ip)/det(ip)
         vc(1,2,ip)=-covc(1,2,ip)/det(ip)
         vc(2,1,ip)=-covc(2,1,ip)/det(ip)
       else
! sinon appelle nag
         do i=1,p_dpm%ncar
           do j=1,p_dpm%ncar
             A(i,j)=covc(i,j,ip)
             Ab(i,j)=covc(i,j,ip)
           end do
         end do
         ifail=1
        ! call MATH_QTLMAP_F03ABF(Ab,ncar,ncar,determ,ifail)
      !   if (ifail.ne.0) then
          ! MODIFICATION COMPORTEMENT
          ! OFI:Si le determinant n'est pas calculable, on considere que les parametre de la fonction
          ! ne sont pas viable, on met une valeur tres haute au LRT

          !  fp21 = INIFINY_REAL_VALUE
         !   f  = INIFINY_REAL_VALUE
         !   fm21 = INIFINY_REAL_VALUE
          !  return
       !  end if
         CALL MATH_QTLMAP_INVDETMATSYM(p_dpm%ncar,A,p_dpm%ncar+1,determ,ifail)
         det(ip)=determ

         !ifail=1
         !call MATH_QTLMAP_F01ADF(ncar,A,ncar+1,ifail)
         if (ifail.ne.0 .or. det(ip)< 1.d-9) then
           !print*,'ifail=',ifail,' ; inversion matrice0'
           !stop
            ! MODIFICATION COMPORTEMENT
          ! OFI:Si le determinant n'est pas calculable, on considere que les parametre de la fonction
          ! ne sont pas viable, on met une valeur tres haute au LRT
            fp21 = INIFINY_REAL_VALUE
            f  = INIFINY_REAL_VALUE
            fm21 = INIFINY_REAL_VALUE
            return
         end if

         do i=1,p_dpm%ncar
           do j=1,i
            vc(i,j,ip)=A(i+1,j)
            vc(j,i,ip)=vc(i,j,ip)
           end do
         end do
       end if
!
! Initialisation polyg�nique pere
       do i=1,p_dpm%ncar
          xmupc(i,ip)=x(p_dpm%ncar*p_dg%np+(i-1)*p_dg%np+ip)
          xmup2c(i,ip)=xmupc(i,ip)*xmupc(i,ip)
       end do
!
! Calcul de la vraisemblance
!
! Vraisemblance demi-fr�res
       zdf=0.d0
       vdf=0.d0
       do ic=1,p_dpm%ncar
        vdf=carydf01(ic,ip)+dble(effdf(ip))*xmup2c(ic,ip)-2.d0*xmupc(ic,ip)*somydf01(ic,ip)
        zdf=zdf + vc(ic,ic,ip)*vdf
	do jc=1,ic-1
	 vdf=somyydf01(ic,jc,ip) - xmupc(ic,ip)*somydf01(jc,ip)  &
            - xmupc(jc,ip)*somydf01(ic,ip)+ dble(effdf(ip))*xmupc(ic,ip)*xmupc(jc,ip)
         zdf = zdf + 2.d0*vdf*vc(ic,jc,ip)
	end do
       end do
       f= 0.5d0 * zdf + 0.5*dble(effdf(ip))*dlog(det(ip))
!
! Vraisemblance plein-fr�res
       nest=0
       somxmu=0.d0
       do jm=p_dg%nmp(ip)+1,p_dg%nmp(ip+1)
        zpf=0.d0
        ! WARNING********************************
        ! On estime par rapport a quel caractere
        if(p_dpa%estime(p_dpm%ncar,jm)) then
         nest=nest+1
         if (ip>1) then
               nestim=sum(estmum(:ip-1))+nest
         else
               nestim=nest
         end if

         do ic=1,p_dpm%ncar
!          print *,"ip:",ip,nest,estmum(ip)
          if(nest.le.estmum(ip)) then
	!    if(ic.eq.1) nestim=nestim+1
	      !  print *,"polygenique mere:",x(2*ncar*np+(ic-1)*nmumest(1)+nestim),' ic:',ic
            xmumc(ic,jm)=x(2*p_dpm%ncar*p_dg%np+(ic-1)*p_dpa%nmumest(1)+nestim)
            somxmu(ic)=somxmu(ic)+xmumc(ic,jm)
          else
            xmumc(ic,jm)=-somxmu(ic)
          end if
          xmum2c(ic,jm)=xmumc(ic,jm)*xmumc(ic,jm)
          xmupmc(ic,jm)=xmupc(ic,ip)+xmumc(ic,jm)

          vpf=dble(eff(jm))*(xmup2c(ic,ip)+xmum2c(ic,jm)         &
     	  +2.d0*xmupc(ic,ip)*xmumc(ic,jm))                       &
         -2.d0*xmupmc(ic,jm)*somy11(ic,jm)+cary11(ic,jm)
          zpf=zpf + vpf*vc(ic,ic,ip)
          if (ic.gt.1) then
           do jc=1,ic-1
            vpf=dble(eff(jm))*xmupmc(ic,jm)*xmupmc(jc,jm)             &
            +somyy11(ic,jc,jm)-xmupmc(ic,jm)*somy11(jc,jm)-xmupmc(jc,jm)*somy11(ic,jm)
            zpf = zpf + 2.d0*vpf*vc(ic,jc,ip)
           end do
          end if
         end do
         fm21(jm)= 0.5d0*zpf+0.5d0*dble(eff(jm))*dlog(det(ip))
        else
	     fm21(jm)=0.d0
        end if
        f=f+fm21(jm)
       end do
      return
      end subroutine funct_mcar_0qtl_family
!!***

!!****f* m_qtlmap_analyse_multitrait/opti_mcar_1qtl
!!  NAME
!!    opti_mcar_1qtl
!!  DESCRIPTION
!!    Calcul de la statistique de test multicaractere le long du chromosome
!!  SOURCE
      subroutine opti_mcar_1qtl(dataset,spt,lrtsol)
      type(QTLMAP_DATASET)       ,intent(in) :: dataset
      type(PDD_BUILD)  ,target   ,intent(in) :: spt
      type(TYPE_LRT_SOLUTION)  , intent(inout)            :: lrtsol

      real (kind=dp) , dimension(dataset%phenoModel%ncar,dataset%phenoModel%ncar) :: r,rhof


!
! Tableaux dimensionnes selon nm, le nombre de meres

! Divers
      !real (kind=dp), dimension((3*ncar*np)+ncar*nmumest(1)+ncar*namest(1)+ncar*(ncar-1)/2) :: par,borni,borns
      real (kind=dp), dimension(:),allocatable :: par,borni,borns
      integer iuser(1)
      double precision user(1)
      integer :: nrho,npar,ifail,ibound,i,k,ip,chr,ntotal,indexchr(dataset%map%nchr),nprime
      integer :: jm,nm1,nm2,nd1,nd2,kd,kkd,ifem,j,ilong,n,ix
      integer :: ngeno1,ngeno2,ig,nestim,nest,indam

      real (kind=dp) , dimension(:,:,:),pointer  :: listepar
      real (kind=dp) :: somxmu,pp,pm,xlrtc_t,f2,f2_2
      logical  , dimension(:,:),pointer        :: filter_inc
      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa
      type(DATAMODEL_BASE) , pointer :: dpm

      real(kind=dp)       ,dimension(:,:,:)   ,pointer     :: xlrp,xlrm
      real(kind=dp)       ,dimension(:,:)   ,pointer       :: lrt1

      dpm => dataset%phenoModel
      dg => dataset%genea
      dpa => dataset%phenoAnimal

      p_dg => dataset%genea
      p_dpm => dataset%phenoModel
      p_dpa => dataset%phenoAnimal
      p_spt => spt

      call lrtsol%new(dataset,1)

      allocate (lrt1(dataset%map%nchr,dataset%map%get_maxnpo()))
      allocate (xlrp(dataset%map%nchr,dataset%map%get_maxnpo(),dg%np))
      allocate (xlrm(dataset%map%nchr,dataset%map%get_maxnpo(),dg%nm))


!****************************************************************************
! Calcul de la vraisemblance sous H2
! Parametres de maximisation
      nrho=dpm%ncar*(dpm%ncar-1)/2
      npar=(3*dpm%ncar*dg%np)+dpm%ncar*dpa%nmumest(1)+dpm%ncar*dpa%namest(1)+nrho

      allocate (filter_inc(dg%np,npar))
      allocate (par(npar),borni(npar),borns(npar))
      par=0.d0
      filter_inc=.false.

      ibound=0
      rhof = 0.d0
      k=0
      do i=1,dpm%ncar
       nestim=0
       do ip=1,dg%np
        borni((i-1)*dg%np+ip)=SIG_MIN                          ! variance de la famille ip pour le caractere i
        borns((i-1)*dg%np+ip)=SIG_MAX
        filter_inc(ip,(i-1)*dg%np+ip)=.true.
        borni(ip+dpm%ncar*dg%np+(i-1)*dg%np)=DEFAULT_PARAM_MIN        ! moyenne de la famille ip pour le caractere i
        borns(ip+dpm%ncar*dg%np+(i-1)*dg%np)=DEFAULT_PARAM_MAX
        filter_inc(ip,ip+dpm%ncar*dg%np+(i-1)*dg%np)=.true.
        borni(ip+2*dpm%ncar*dg%np+(i-1)*dg%np)=AM_MIN                 ! effet qtl de la famille ip pour le caractere i
        borns(ip+2*dpm%ncar*dg%np+(i-1)*dg%np)= AM_MAX
        filter_inc(ip,ip+2*dpm%ncar*dg%np+(i-1)*dg%np)=.true.

        nest=0
        do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
           if(dpa%estime(1,jm)) then
            nest=nest+1
            if(nest.le.estmum(ip)) then
              nestim=nestim+1
              filter_inc(ip,3*dpm%ncar*dg%np+(i-1)*dpa%nmumest(1)+nestim)=.true.
            end if
            ifem=dg%repfem(jm)
            filter_inc(ip,3*dpm%ncar*dg%np+dpm%ncar*dpa%nmumest(1)+(i-1)*dpa%namest(1)+dpa%iam(1,ifem))=.true.
           end if
        end do
       end do ! ip
       do jm=1,dpa%nmumest(1)
         borni(3*dpm%ncar*dg%np+(i-1)*dpa%nmumest(1)+jm)=DEFAULT_PARAM_MIN    ! moyenne de la famille ip/jm pour le caractere i
         borns(3*dpm%ncar*dg%np+(i-1)*dpa%nmumest(1)+jm)=DEFAULT_PARAM_MAX
       end do
       do jm=1,dpa%namest(1)
         borni(3*dpm%ncar*dg%np+dpm%ncar*dpa%nmumest(1)+(i-1)*dpa%namest(1)+jm)=AM_MIN ! effet qtl de la famille ip/ifem pour le caractere i
         borns(3*dpm%ncar*dg%np+dpm%ncar*dpa%nmumest(1)+(i-1)*dpa%namest(1)+jm)=AM_MAX
       end do
       if (i.gt.1) then
        do j=1,i-1
         k=k+1
         borni(3*dpm%ncar*dg%np+dpm%ncar*dpa%nmumest(1)+dpm%ncar*dpa%namest(1)+k)=DEFAULT_PARAM_MIN
         borns(3*dpm%ncar*dg%np+dpm%ncar*dpa%nmumest(1)+dpm%ncar*dpa%namest(1)+k)=DEFAULT_PARAM_MAX
        end do
       end if
      end do

      filter_inc(:,3*dpm%ncar*dg%np+dpm%ncar*dpa%nmumest(1)+dpm%ncar*dpa%namest(1)+1:)=.true.
!
!
! Point de depart
      k=0
      do i=1,dpm%ncar
        do ip=1,dg%np
          par((i-1)*dg%np+ip)=sig2(i,ip)
          par(ip+dpm%ncar*dg%np+(i-1)*dg%np)=xmu2p(i,ip)
          par(ip+2*dpm%ncar*dg%np+(i-1)*dg%np)=0.d0
        end do
        do jm=1,dpa%nmumest(1)
          par(3*dpm%ncar*dg%np+(i-1)*dpa%nmumest(1)+jm)=xmu2m(i,jm)
        end do
        do jm=1,dpa%namest(1)
          par(3*dpm%ncar*dg%np+dpm%ncar*dpa%nmumest(1)+(i-1)*dpa%namest(1)+jm)= 0.d0
        end do
        if (i.gt.1) then
         do j=1,i-1
          k=k+1
          par(3*dpm%ncar*dg%np+dpm%ncar*dpa%nmumest(1)+dpm%ncar*dpa%namest(1)+k)=rhoi2(j,i)
         end do
        end if
      end do


      allocate (listepar(dataset%map%nchr,dataset%map%get_maxnpo(),npar))
      ibound=0
! Marche le long du chromosome
      lrtsol%lrtmax=-1.d75

      ntotal=0
      do chr=1,dataset%map%nchr
        ntotal=ntotal+dataset%map%get_npo(chr)
        indexchr(chr)=ntotal
      end do

      !$OMP PARALLEL DEFAULT(SHARED) FIRSTPRIVATE(par)  &
      !$OMP PRIVATE(ip,jm,chr,i,n,nd1,nd2,pp,pm,kd,kkd,ig,ifail,nestim,f2)
      call init_sub_1qtl(dataset,spt)
      !$OMP DO
      do nprime=1,ntotal
       n=nprime
       do chr=dataset%map%nchr,1,-1
        if ( indexchr(chr) >= nprime) then
          current_chr = chr
        else
          n=n-dataset%map%get_npo(chr)
        end if
       end do
       chr=current_chr

        do ip=1,dg%np
          somppdf11(ip)=0.d0
          carppdf11(ip)=0.d0
          do i=1,dpm%ncar
           somppydf11(i,ip)=0.d0
          end do
          do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
            do ig=spt%ngenom(current_chr,jm)+1,spt%ngenom(current_chr,jm+1)
              sompp11(ig)=0.d0
              carpp11(ig)=0.d0
              do i=1,dpm%ncar
                somppy11(i,ig)=0.d0
                sompmy11(i,ig)=0.d0
              end do
              sompm11(ig)=0.d0
              carpm11(ig)=0.d0
              somppm11(ig)=0.d0
              nd1=spt%ngend(current_chr,ig)+1
              nd2=spt%ngend(current_chr,ig+1)
              do kd=nd1,nd2
                kkd=spt%ndesc(current_chr,kd)
               if(count(dpa%presentc(:,kkd)) == dpm%ncar )then
                  pp=-spt%pdd(current_chr,kd,1,n)-spt%pdd(current_chr,kd,2,n)+&
                   spt%pdd(current_chr,kd,3,n)+spt%pdd(current_chr,kd,4,n)
                  pm=-spt%pdd(current_chr,kd,1,n)+spt%pdd(current_chr,kd,2,n)-&
                   spt%pdd(current_chr,kd,3,n)+spt%pdd(current_chr,kd,4,n)
                  if(dpa%estime(dpm%ncar,jm)) then
                   sompp11(ig)=sompp11(ig)+pp
                   carpp11(ig)=carpp11(ig)+pp*pp
                   do i=1,dpm%ncar
                    somppy11(i,ig)=somppy11(i,ig)+pp*dpa%y(i,kkd)
                    sompmy11(i,ig)=sompmy11(i,ig)+pm*dpa%y(i,kkd)
                   end do
                   sompm11(ig)=sompm11(ig)+pm
                   carpm11(ig)=carpm11(ig)+pm*pm
                   somppm11(ig)=somppm11(ig)+pp*pm
                  else
                   somppdf11(ip)=somppdf11(ip)+pp
                   carppdf11(ip)=carppdf11(ip)+pp*pp
                   do i=1,dpm%ncar
                    somppydf11(i,ip)=somppydf11(i,ip)+pp*dpa%y(i,kkd)
                   end do
                  end if
                 end if
               end do                  !! fin kd
             end do                    !! fin ig
           end do                      !! fin jm
         end do                        !! fin ip

! Optimisation de la vraisemblance a la position dx
         ifail=1
        ! call minimizing_funct(npar,ibound,funct_mcar_1qtl,borni,borns,par,f2,iuser,user,ifail)
         call minimizing_funct_family_sire(dataset,npar,ibound,funct_mcar_1qtl_family,filter_inc,&
          fp2,borni,borns,par,f2,iuser,user,ifail)

         listepar(chr,n,:npar) = par(:npar)
         if ( f2 < INIFINY_REAL_VALUE ) then
           lrt1(chr,n)=-2.d0*(f2-f01)
         else
           lrt1(chr,n)=0
         end if

         ! seulement le dernier effet qtl (du dernier caracteres)........
         ! A corriger...
         do i=1,dpm%ncar
           do ip=1,dg%np
             xlrp(chr,n,ip)=-2.d0*(fp2(ip)-fp21(ip))
             lrtsol%pater_eff(chr,ip,n)=par(ip+2*dpm%ncar*dg%np+(i-1)*dg%np)
           end do
           do jm=1,dg%nm
             xlrm(chr,n,jm)=-2.d0*(fm2(jm)-fm21(jm))
           end do
           do jm=1,dpa%namest(1)
             lrtsol%mater_eff(chr,jm,n)=par(3*dpm%ncar*dg%np+dpm%ncar*dpa%nmumest(1)+(i-1)*dpa%namest(1)+jm)
           end do
         end do
      end do
      !$OMP END DO
      call end_sub_1qtl
      !$OMP END PARALLEL

      !recherche du maximum
      do chr=1,dataset%map%nchr
       do n=1,dataset%map%get_npo(chr)
          if(lrtsol%lrtmax(0) < lrt1(chr,n)) then
           lrtsol%lrtmax(0)=lrt1(chr,n)
           lrtsol%nxmax(0)=n
           lrtsol%chrmax(0)=chr
          end if
       end do
      end do
      !initilisation des valeurs par rapport au maximum trouve
      par(:npar)=listepar(lrtsol%chrmax(0),lrtsol%nxmax(0),:npar)
      k=0
      rhof=0.d0
      do i=2,dpm%ncar
         do j=1,i-1
            k=k+1
            r(j,i)=par(3*dpm%ncar*dg%np+dpm%ncar*dpa%nmumest(1)+dpm%ncar*dpa%namest(1)+k)
            rhof(j,i)=(dexp(r(j,i))-1.d0)/(dexp(r(j,i))+1.d0)
            rhof(i,j)=rhof(j,i)
          end do
      end do
      do i=1,dpm%ncar
         do j=1,dpm%ncar
           rhoi2_2(i,j)=rhof(i,j)
         end do
         nestim=0
         do ip=1,dg%np
          std(i,ip)=par(ip+(i-1)*dg%np)
          xmoyp(i,ip)=par(ip+dpm%ncar*dg%np+(i-1)*dg%np)
          ap(i,ip)=par(ip+2*dpm%ncar*dg%np+(i-1)*dg%np)
          nest=0
          somxmu = 0.d0
          do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
            if (count(dpa%estime(:,jm))==dpm%ncar) then
               nest=nest+1
               ifem=dg%repfem(jm)
               indam=dpa%iam(1,ifem)
               am(i,jm)=par(3*dpm%ncar*dg%np+dpm%ncar*dpa%nmumest(1)+(i-1)*dpa%namest(1)+indam)
               if (nest.le.estmum(ip)) then
                  nestim=nestim+1
                  xmoym(i,jm)=par(3*dpm%ncar*dg%np+(i-1)*dpa%nmumest(1)+nestim)
                  somxmu=somxmu+xmoym(i,jm)
                else
                  xmoym(i,jm)=-somxmu
                end if
             else
               am(i,jm)=0.d0
               xmoym(i,jm)=0.d0
             end if
           end do !jm
          end do !ip
       end do !icar

       ! 04/2013 - New structure for LRT
      call lrtsol%LRT%add1(dataset,dataset%map%nchr,dataset%map%get_maxnpo(),lrt1)
      do ip=1,dg%np
       call lrtsol%LRT_SIRES(ip)%add1(dataset,dataset%map%nchr,dataset%map%get_maxnpo(),xlrp(:,:,ip))
      end do
      do jm=1,dg%nm
       call lrtsol%LRT_DAMS(jm)%add1(dataset,dataset%map%nchr,dataset%map%get_maxnpo(),xlrm(:,:,jm))
      end do

      deallocate (lrt1)
      deallocate (xlrp)
      deallocate (xlrm)

      deallocate (listepar)
      deallocate (filter_inc)
      deallocate (par,borni,borns)

      end subroutine opti_mcar_1qtl
!!***

!!****f* m_qtlmap_analyse_multitrait/funct_mcar_1qtl
!!  NAME
!!    funct_mcar_1qtl
!!  DESCRIPTION
!!    Calcul de la vraisemblance et de ses derivees partielles sous H1
!!  SOURCE
      subroutine funct_mcar_1qtl(n,x,f,iuser,user)

      integer         , intent(in)                  :: n
      real (kind=dp)      ,dimension(n), intent(in) :: x
      real (kind=dp)    , intent(inout)             :: f
      integer ,       dimension(1) ,intent(in)               :: iuser
      real (kind=dp)      ,dimension(1),intent(in)             :: user

      integer :: ip !fp2

      do ip=1,p_dg%np
        call funct_mcar_1qtl_family(ip,n,x,fp2(ip),iuser,user)
      end do
      f = sum(fp2)

      end subroutine funct_mcar_1qtl
!!***

!!****f* m_qtlmap_analyse_multitrait/funct_mcar_1qtl_family
!!  NAME
!!    funct_mcar_1qtl_family
!!  DESCRIPTION
!!    Calcul de la vraisemblance et de ses derivees partielles sous H1
!!  SOURCE
      subroutine funct_mcar_1qtl_family(ip,n,x,f,iuser,user)

      integer         , intent(in)                  :: ip,n
      real (kind=dp)      ,dimension(n), intent(in) :: x
      real (kind=dp)    , intent(inout)             :: f
      integer ,       dimension(1) ,intent(in)               :: iuser
      real (kind=dp)      ,dimension(1),intent(in)             :: user

      real (kind=dp) ,dimension(p_dg%np) :: det
!
! Tableaux dimensionnes selon le nombre nc de caracteres
      real (kind=dp),dimension(p_dpm%ncar,p_dg%np):: sigc,xmupc,xmup2c,apc,ap2c
      real (kind=dp),dimension(p_dpm%ncar) :: wrkspce,somxmu
      real (kind=dp),dimension(p_dpm%ncar,p_dg%nm)::  xmumc,xmum2c,xmupmc,amc,am2c
      real (kind=dp),dimension(p_dpm%ncar,p_dpm%ncar,p_dg%np) ::covc,vc
      real (kind=dp),dimension(p_dpm%ncar+1,p_dpm%ncar) :: A
      real (kind=dp),dimension(p_dpm%ncar,p_dpm%ncar) :: rhoc,Ab

      real (kind=dp)::rh,zdf,vdf,vmere,zpf,vpf,determ,savezpf
      integer :: i,j,k,nt,nestim,ifail,it,nest,jm
      integer :: ifem,indam,ic,ngeno1,ngeno2,ig,ipar,itindi,jc

!***********************************************************************************
      f=0.d0
      k=0
      rhoc=0.d0
! Recuperation des correlations
      do i=2,p_dpm%ncar
        do j=1,i-1
          k=k+1
          nt=3*p_dpm%ncar*p_dg%np+p_dpm%ncar*p_dpa%nmumest(1)+p_dpm%ncar*p_dpa%namest(1)+k
          rh=x(nt)
          rhoc(j,i)=((dexp(rh)-1.d0)/(dexp(rh)+1.d0))
          rhoc(i,j)=rhoc(j,i)
        end do
      end do
!
! Calcul de l'inverse et du determinant de la matrice de var-cov
       nestim=0
!
! Initialisation des parametres
         do i=1,p_dpm%ncar
          sigc(i,ip)=dabs(x((i-1)*p_dg%np+ip))
          covc(i,i,ip)=sigc(i,ip)*sigc(i,ip)
          if(i.gt.1) then
            do j=1,i-1
              covc(i,j,ip)=rhoc(i,j)*sigc(i,ip)*sigc(j,ip)
              covc(j,i,ip)=covc(i,j,ip)
            end do
          end if
          xmupc(i,ip)=x(p_dpm%ncar*p_dg%np+(i-1)*p_dg%np+ip)
          xmup2c(i,ip)=xmupc(i,ip)*xmupc(i,ip)
          apc(i,ip)=x(ip+2*p_dpm%ncar*p_dg%np+(i-1)*p_dg%np)
          ap2c(i,ip)=apc(i,ip)*apc(i,ip)
         end do                             !!i
! Si 2 caracteres, calcul a la main
         if(p_dpm%ncar.eq.2) then
          det(ip)=covc(1,1,ip)*covc(2,2,ip)-covc(1,2,ip)*covc(2,1,ip)
          vc(1,1,ip)=covc(2,2,ip)/det(ip)
          vc(2,2,ip)=covc(1,1,ip)/det(ip)
          vc(1,2,ip)=-covc(1,2,ip)/det(ip)
          vc(2,1,ip)=-covc(2,1,ip)/det(ip)
         else
! sinon appelle nag
          do i=1,p_dpm%ncar
           do j=1,p_dpm%ncar
             A(i,j)=covc(i,j,ip)
             Ab(i,j)=covc(i,j,ip)
           end do
          end do
          ifail=1
         ! call MATH_QTLMAP_F03ABF(Ab,ncar,ncar,determ,ifail)

        !  if (ifail.ne.0) then
          ! MODIFICATION COMPORTEMENT
          ! OFI:Si le determinant n'est pas calculable, on considere que les parametre de la fonction
          ! ne sont pas viable, on met une valeur tres haute au LRT
       !     fp21 = INIFINY_REAL_VALUE
       !     f  = INIFINY_REAL_VALUE
        !    fm21 = INIFINY_REAL_VALUE
        !    return
        !  end if
          CALL MATH_QTLMAP_INVDETMATSYM(p_dpm%ncar,A,p_dpm%ncar+1,determ,ifail)
          det(ip)=determ
!          ifail=1
  !        call MATH_QTLMAP_F01ADF(ncar,A,ncar+1,ifail)
          if (ifail.ne.0 .or. det(ip)<=0 ) then
             !call stop_application('ifail='//trim(str(ifail))//' ; inversion matrice')
              ! MODIFICATION COMPORTEMENT
          ! OFI:Si le determinant n'est pas calculable, on considere que les parametre de la fonction
          ! ne sont pas viable, on met une valeur tres haute au LRT
            fp21 = INIFINY_REAL_VALUE
            f  = INIFINY_REAL_VALUE
            fm21 = INIFINY_REAL_VALUE
            return
          end if
          do i=1,p_dpm%ncar
            do j=1,i
              vc(i,j,ip)=A(i+1,j)
              vc(j,i,ip)=vc(i,j,ip)
    !          print *,vc(i,j,ip)
            end do
          end do
         end if
! Initialisation parametres des m�res
	     nest=0
	     somxmu=0.d0
         do jm=p_dg%nmp(ip)+1,p_dg%nmp(ip+1)
          !WARNING -->utilisation de estime
          if(count(p_dpa%estime(:,jm)) == p_dpm%ncar) then
           nest=nest+1
           if (ip>1) then
               nestim=sum(estmum(:ip-1))+nest
           else
               nestim=nest
           end if

           do i=1,p_dpm%ncar
            if(nest.le.estmum(ip)) then
!	      if (i.eq.1) nestim=nestim+1
              xmumc(i,jm)=x(3*p_dpm%ncar*p_dg%np+(i-1)*p_dpa%nmumest(1)+nestim)
              somxmu(i)=somxmu(i)+xmumc(i,jm)
            else
              xmumc(i,jm)=-somxmu(i)
            end if
            xmum2c(i,jm)=xmumc(i,jm)*xmumc(i,jm)
            xmupmc(i,jm)=xmupc(i,ip)+xmumc(i,jm)
            indam=p_dpa%iam(1,p_dg%repfem(jm))
            amc(i,jm)=x(3*p_dpm%ncar*p_dg%np+p_dpm%ncar*p_dpa%nmumest(1)+(i-1)*p_dpa%namest(1)+indam)
            am2c(i,jm)=amc(i,jm)*amc(i,jm)
           end do
          end if
         end do
!      end do                   !!ip
!
! Calcul de la vraisemblance
 !     do ip=1,np
! Vraisemblance demi-fr�res
        zdf=0.d0
	    vdf=0.d0
        do ic=1,p_dpm%ncar
         vdf=carydf01(ic,ip)+dble(effdf(ip))*xmup2c(ic,ip)-2.d0*xmupc(ic,ip)*somydf01(ic,ip)
         vdf = vdf + ap2c(ic,ip)*carppdf11(ip)                  &
                  +2.d0*xmupc(ic,ip)*apc(ic,ip)*somppdf11(ip)   &
                  -2.d0*apc(ic,ip)*somppydf11(ic,ip)
         zdf=zdf + vc(ic,ic,ip)*vdf
	 do jc=1,ic-1
	  vdf=somyydf01(ic,jc,ip) - xmupc(ic,ip)*somydf01(jc,ip) &
             - xmupc(jc,ip)*somydf01(ic,ip)                  &
             + dble(effdf(ip))*xmupc(ic,ip)*xmupc(jc,ip)
          vdf=vdf - apc(ic,ip)*somppydf11(jc,ip)             &
                 - apc(jc,ip)*somppydf11(ic,ip)              &
                 + apc(ic,ip)*apc(jc,ip)*carppdf11(ip)       &
                 + xmupc(ic,ip)*apc(jc,ip)*somppdf11(ip)     &
                 + xmupc(jc,ip)*apc(ic,ip)*somppdf11(ip)
	  zdf = zdf + 2.d0*vdf*vc(ic,jc,ip)
	    end do
	  end do

      f= 0.5d0 * zdf + 0.5d0*dble(effdf(ip))*dlog(det(ip))

!
! Vraisemblance plein-fr�res
       do jm=p_dg%nmp(ip)+1,p_dg%nmp(ip+1)
        itindi=0
        zpf=0.d0
        ! WARNING : utilisation de estime
        if(p_dpa%estime(p_dpm%ncar,jm)) then
         do ic=1,p_dpm%ncar
          vpf=dble(eff(jm))*(xmup2c(ic,ip)+xmum2c(ic,jm)    &
     	  +2.d0*xmupc(ic,ip)*xmumc(ic,jm))                  &
         -2.d0*xmupmc(ic,jm)*somy11(ic,jm)+cary11(ic,jm)
          zpf=zpf + vpf*vc(ic,ic,ip)
          if (ic.gt.1) then
           do jc=1,ic-1
            vpf=dble(eff(jm))*xmupmc(ic,jm)*xmupmc(jc,jm) &
            +somyy11(ic,jc,jm)-xmupmc(ic,jm)*somy11(jc,jm) &
            -xmupmc(jc,jm)*somy11(ic,jm)
            zpf = zpf + 2.d0*vpf*vc(ic,jc,ip)
           end do
          end if
         end do
         ngeno1=p_spt%ngenom(current_chr,jm)+1
         ngeno2=p_spt%ngenom(current_chr,jm+1)
! si un seul genotype possible pour la mere, ne passe pas par les exp(vpf)
        if ((ngeno2-ngeno1).eq.0) then
    !     if (.true.) then
           ig=ngeno1
	   do ic=1,p_dpm%ncar
            vpf= ap2c(ic,ip)*carpp11(ig)                  &
           +am2c(ic,jm)*carpm11(ig)                       &
           +2.d0*apc(ic,ip)*amc(ic,jm)*somppm11(ig)       &
           +2.d0*xmupmc(ic,jm)*(apc(ic,ip)*sompp11(ig)    &
           +amc(ic,jm)*sompm11(ig))                       &
           -2.d0*apc(ic,ip)*somppy11(ic,ig)               &
           -2.d0*amc(ic,jm)*sompmy11(ic,ig)
            zpf= zpf + vpf*vc(ic,ic,ip)
          !  if ((ngeno2-ngeno1)>0)  print *,'ic:',ic,zpf,vpf,vc(ic,ic,ip)
            if (ic.gt.1) then
             do jc=1,ic-1
               vpf= xmupmc(ic,jm)*(apc(jc,ip)*sompp11(ig)   &
              +amc(jc,jm)*sompm11(ig))                      &
              +xmupmc(jc,jm)*(apc(ic,ip)*sompp11(ig)        &
              +amc(ic,jm)*sompm11(ig))                      &
              -apc(jc,ip)*somppy11(ic,ig)                   &
              -amc(jc,jm)*sompmy11(ic,ig)                   &
              -apc(ic,ip)*somppy11(jc,ig)                   &
              -amc(ic,jm)*sompmy11(jc,ig)                   &
              +apc(ic,ip)*apc(jc,ip)*carpp11(ig)            &
              +amc(ic,jm)*amc(jc,jm)*carpm11(ig)            &
              +(apc(ic,ip)*amc(jc,jm)+amc(ic,jm)*apc(jc,ip)) &
              *somppm11(ig)
               zpf= zpf + 2.d0*vpf*vc(ic,jc,ip)
             end do                                            !! jc
            end if
           end do !! ic

          fm2(jm)= 0.5d0*zpf+0.5d0*dble(eff(jm))*dlog(det(ip))
         else
! si plusieurs genotypes possibles pour la mere,
          vmere=0.d0
          savezpf=zpf
          do ig=ngeno1,ngeno2 !attention NGENO2 *******************************************************************************
           zpf=savezpf
           do ic=1,p_dpm%ncar
             vpf= ap2c(ic,ip)*carpp11(ig)                  &
           +am2c(ic,jm)*carpm11(ig)                       &
           +2.d0*apc(ic,ip)*amc(ic,jm)*somppm11(ig)       &
           +2.d0*xmupmc(ic,jm)*(apc(ic,ip)*sompp11(ig)    &
           +amc(ic,jm)*sompm11(ig))                       &
           -2.d0*apc(ic,ip)*somppy11(ic,ig)               &
           -2.d0*amc(ic,jm)*sompmy11(ic,ig)
            zpf= zpf + vpf*vc(ic,ic,ip)

            if (ic.gt.1) then
             do jc=1,ic-1
                vpf= xmupmc(ic,jm)*(apc(jc,ip)*sompp11(ig) &
              +amc(jc,jm)*sompm11(ig))                     &
              +xmupmc(jc,jm)*(apc(ic,ip)*sompp11(ig)       &
              +amc(ic,jm)*sompm11(ig))                     &
              -apc(jc,ip)*somppy11(ic,ig)                  &
              -amc(jc,jm)*sompmy11(ic,ig)                  &
              -apc(ic,ip)*somppy11(jc,ig)                  &
              -amc(ic,jm)*sompmy11(jc,ig)                  &
              +apc(ic,ip)*apc(jc,ip)*carpp11(ig)           &
              +amc(ic,jm)*amc(jc,jm)*carpm11(ig)            &
              +(apc(ic,ip)*amc(jc,jm)+amc(ic,jm)*apc(jc,ip)) &
              *somppm11(ig)
               zpf= zpf + 2.d0*vpf*vc(ic,jc,ip)
             !  print *,'zpf:',zpf,ig,ic,jc
             end do                                            !! jc
            end if
          end do                       !! ic
  !        print *,'=>',ip,jm,zpf
          zpf=-0.5d0*zpf
 !         print *,'-0.5*zpf:',zpf
          zpf=dexp(zpf)
!          print *,'dexp(zpf):',zpf
          if ( zpf == 0 ) zpf = 1/huge(zpf)

          if ((zpf .ne. zpf) ) then
	         fm2=INIFINY_REAL_VALUE
             fp2=INIFINY_REAL_VALUE
             f=INIFINY_REAL_VALUE
	       !  print *,"FIN AVANT..."
	        ! stop
	         return
           end if


          vmere=vmere+p_spt%probg(current_chr,ig)*zpf
    !      print *,'vmere:',vmere

         end do!! ig


         if ( vmere /= 0) then
            fm2(jm)=-dlog(vmere)+0.5d0*dble(eff(jm))*dlog(det(ip))
      !      print *,'fm2:',fm2(jm),-dlog(vmere)
      !      stop
         else
	        fm2=INIFINY_REAL_VALUE
	        fp2=INIFINY_REAL_VALUE
	        f=INIFINY_REAL_VALUE
        	return
         end if

        end if
!
        else
	     fm2(jm)=0.d0
        end if
        f=f+fm2(jm)
       end do

      end subroutine funct_mcar_1qtl_family
!!***

!!****f* m_qtlmap_analyse_multitrait/set_solution_hypothesis0
!!  NAME
!!    set_solution_hypothesis0
!!  DESCRIPTION
!!
!!  SOURCE
      subroutine set_solution_hypothesis0(dataset,incsol)
       type(QTLMAP_DATASET)       ,intent(in) :: dataset
       type(TYPE_INCIDENCE_SOLUTION)    ,intent(inout)    :: incsol

       integer :: nteff,maxNbPar,jm,ifem,ip,ieff,ic

       type(GENEALOGY_BASE) , pointer :: dg
       type(PHENOTYPE_BASE) , pointer :: dpa
       type(DATAMODEL_BASE) , pointer :: dpm

       dpm => dataset%phenoModel
       dg => dataset%genea
       dpa => dataset%phenoAnimal

       incsol%hypothesis=0
       allocate (incsol%sig(dpm%ncar,dg%np))
       !  Mean family
       nteff = dpm%ncar
       maxNbPar = dg%np
       do ic=1,dpm%ncar
        if ( count(dpa%estime(ic,:))  > 0 ) then
          nteff = nteff + 1
          maxNbPar = max(maxNbPar,count(dpa%estime(ic,:)))
        end if
       end do

       allocate (incsol%groupeName(nteff))
       allocate (incsol%nbParameterGroup(nteff))
       allocate (incsol%parameterName(nteff,maxNbPar))
       allocate (incsol%paramaterValue(nteff,maxNbPar))
       allocate (incsol%parameterVecsol(nteff,maxNbPar))
       allocate (incsol%parameterPrecis(nteff,maxNbPar))

       incsol%parameterVecsol=.true.
       incsol%parameterPrecis=0.d0

       ieff=0
       do ic=1,dpm%ncar
         do ip=1,dg%np
            incsol%sig(ic,ip) = sig2(ic,ip)*dpm%sigt(ic)
         end do

         ieff=ieff+1
         incsol%groupeName(ieff) = 'Mean Sire ['//trim(dpm%carac(ic))//"]"
         incsol%nbParameterGroup(ieff)=dg%np

          do ip=1,dg%np
           incsol%parameterName(ieff,ip)='Sire '//trim(dg%pere(ip))
           incsol%paramaterValue(ieff,ip) = xmu2p(ic,ip)*dpm%sigt(ic) + dpm%xmut(ic)
          end do

         if ( count(dpa%estime(ic,:)) > 0 ) then
           ieff = ieff +1
           incsol%groupeName(ieff) = 'Mean dam ['//trim(dpm%carac(ic))//"]"
           incsol%nbParameterGroup(ieff)= count(dpa%estime(ic,:))
           ifem=0
           do ip=1,dg%np
           do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
             if ( dpa%estime(ic,jm)) then
              ifem=ifem+1
              incsol%parameterName(ieff,ifem)='Dam '//trim(dg%mere(jm))//" [Sire "//trim(dg%pere(ip))//"]"
              incsol%paramaterValue(ieff,ifem) = xmu2m(ic,ifem)*dpm%sigt(ic)
             end if
           end do
          end do
         end if

        end do

        allocate (incsol%rhoi(dpm%ncar,dpm%ncar))
        incsol%rhoi=rhoi2_1

       end subroutine set_solution_hypothesis0
!!***

!!****f* m_qtlmap_analyse_multitrait/set_solution_hypothesis1
!!  NAME
!!    set_solution_hypothesis1
!!  DESCRIPTION
!!
!!  HISTORY
!!   09/09/2010  * add allocation of incsol%qtl_groupeName
!!  SOURCE
      subroutine set_solution_hypothesis1(dataset,spt,incsol)
       type(QTLMAP_DATASET)       ,intent(in) :: dataset
       type(PDD_BUILD)            ,intent(in) :: spt
       type(TYPE_INCIDENCE_SOLUTION)    ,intent(inout)    :: incsol

       integer :: nteff,maxNbPar,jm,ifem,ip,ieff,ic
       type(GENEALOGY_BASE) , pointer :: dg
       type(PHENOTYPE_BASE) , pointer :: dpa
       type(DATAMODEL_BASE) , pointer :: dpm

       dpm => dataset%phenoModel
       dg => dataset%genea
       dpa => dataset%phenoAnimal

       incsol%hypothesis=1
       allocate (incsol%sig(dpm%ncar,dg%np))
       !  Mean family , Qtl effect
       nteff = 2 * dpm%ncar
       maxNbPar = dg%np
       do ic=1,dpm%ncar
        if ( count(dpa%estime(ic,:))  > 0 ) then
          nteff = nteff + 2
          maxNbPar = max(maxNbPar,count(dpa%estime(ic,:)))
        end if
       end do

       allocate (incsol%groupeName(nteff))
       allocate (incsol%nbParameterGroup(nteff))
       allocate (incsol%parameterName(nteff,maxNbPar))
       allocate (incsol%paramaterValue(nteff,maxNbPar))
       allocate (incsol%parameterVecsol(nteff,maxNbPar))
       allocate (incsol%parameterPrecis(nteff,maxNbPar))
       allocate (incsol%qtl_groupeName(dpm%ncar,1)) ! 1 qtl

       incsol%parameterVecsol=.true.
       incsol%parameterPrecis=0.d0

       ieff=0
       do ic=1,dpm%ncar

         do ip=1,dg%np
            incsol%sig(ic,ip) = std(ic,ip)*dpm%sigt(ic)
         end do

         ieff=ieff+1
         incsol%groupeName(ieff) = 'Mean Sire ['//trim(dpm%carac(ic))//"]"
         incsol%nbParameterGroup(ieff)=dg%np

          do ip=1,dg%np
           incsol%parameterName(ieff,ip)='Sire '//trim(dg%pere(ip))
           incsol%paramaterValue(ieff,ip) = xmoyp(ic,ip)*dpm%sigt(ic) + dpm%xmut(ic)
          end do

         if ( count(dpa%estime(ic,:)) > 0 ) then
           ieff = ieff +1
           incsol%groupeName(ieff) = 'Mean dam ['//trim(dpm%carac(ic))//"]"
           incsol%nbParameterGroup(ieff)= count(dpa%estime(ic,:))
           ifem=0
           do ip=1,dg%np
           do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
             if ( dpa%estime(ic,jm)) then
              ifem=ifem+1
              incsol%parameterName(ieff,ifem)='Dam '//trim(dg%mere(jm))
              incsol%paramaterValue(ieff,ifem) = xmoym(ic,jm)*dpm%sigt(ic)+ dpm%xmut(ic)
             end if
           end do
          end do
         end if

         ieff = ieff +1
         incsol%qtl_groupeName(ic,1) = ieff
         incsol%groupeName(ieff) = 'Sire Qtl effect ['//trim(dpm%carac(ic))//"]"
         incsol%nbParameterGroup(ieff)=dg%np
         do ip=1,dg%np
           incsol%parameterName(ieff,ip)='Sire '//trim(dg%pere(ip))
           incsol%paramaterValue(ieff,ip) = ap(ic,ip)*dpm%sigt(ic)
         end do

         if ( count(dpa%estime(ic,:)) > 0 ) then
           ieff = ieff +1
           incsol%groupeName(ieff) = 'Dam Qtl effect ['//trim(dpm%carac(ic))//"]"
           incsol%nbParameterGroup(ieff)=count(dpa%estime(ic,:))
           ifem=0
           do jm=1,dg%nm
             if ( dpa%estime(ic,jm)) then
              ifem=ifem+1
              incsol%parameterName(ieff,ifem)='Dam '//trim(dg%mere(jm))
              incsol%paramaterValue(ieff,ifem) = am(ic,jm)*dpm%sigt(ic)
              incsol%parameterVecsol(ieff,ifem) = .true.
             end if
           end do
         end if

        end do

        allocate (incsol%rhoi(dpm%ncar,dpm%ncar))
        incsol%rhoi=rhoi2_2

       end subroutine set_solution_hypothesis1
!!***

end module m_qtlmap_analyse_multitrait
