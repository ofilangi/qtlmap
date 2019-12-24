module m_qtlmap_analyse_biallelic_farnir
    use m_qtlmap_types
    use m_qtlmap_log
    use m_qtlmap_optimization
    implicit none

    type work_bicfarnir
      integer                                      :: ic = -1
      !! ip,jm
      integer          ,dimension(:,:)   ,pointer  :: ndsib => null()
      !! ip,jm,kd
      real(kind=dp)    ,dimension(:,:,:) , pointer :: ysib  => null()
      !! ip,jm,kd
      real(kind=dp)    ,dimension(:,:,:) , pointer :: cdsib => null()

      !! ip,jm
      integer          ,dimension(:,:)   ,pointer  :: ndfull => null()
      !! ip,jm,kd
      real(kind=dp)    ,dimension(:,:,:) , pointer :: yfull => null()
      !! ip,jm,kd
      real(kind=dp)    ,dimension(:,:,:) , pointer :: cdfull => null()
      !! ip,jm => id dam in X
      integer          ,dimension(:,:)   ,pointer  :: damfull => null()

      integer                                      :: ndamestim =0
      !correspondance au tableau yfull et ysib
      integer          ,dimension(:,:)   ,pointer  :: ngend  => null()

      integer          ,dimension(:,:)   ,pointer  :: ndesc  => null()

      integer          ,dimension(:,:)   ,pointer  :: npdd   => null()

      real(kind=dp)    ,dimension(:)   ,pointer  :: alpha1 => null()
      real(kind=dp)    ,dimension(:)   ,pointer  :: alpha2 => null()
      real(kind=dp)    ,dimension(:)   ,pointer  :: alpha3 => null()
      real(kind=dp)    ,dimension(:)   ,pointer  :: alpha4 => null()

      real(kind=dp)    ,dimension(:)   ,pointer  :: parSolH0 => null()
      real(kind=dp)                              :: f0
      real(kind=dp)    ,dimension(:)   ,pointer  :: parSolH1 => null()
      real(kind=dp)                              :: f1

    end type work_bicfarnir

    type work_pos
     integer          ,dimension(:)   ,pointer  :: n    => null()
     integer          ,dimension(:)   ,pointer  :: chr  => null()

    end type work_pos

     !probcas
     ! transformation
     ! correspondance entre
     ! cas =1,2,3,4 => P1, P2/2 , P2/2, 1-P1-P2
     !
     !
     ! cor(1) =>    (0,1,0)    => 1*0 + P1*1 + P2*1  = P1
     ! cor(2) =>    (0,0,0.5)  => 1*0 + P1*0 + P2*0.5 = P2/2
     ! cor(3) =>    (0,0,0.5)  => 1*0 + P1*0 + P2*0.5 = P2/2
     ! cor(4) =>    (1,-1,-1)  => 1*1 + P1*(-1) + P2*(-1) = 1 - P1 - P2

   !  real(kind=dp) , dimension(4) , private , parameter :: probcasOne = (/ 0.0,0.0,0.0,1.0 /)
     real(kind=dp) , dimension(4) , private , parameter :: probcasP1  = (/ 1.0,0.0,0.0,0.0 /)
     real(kind=dp) , dimension(4) , private , parameter :: probcasP2  = (/ 0.0,1.0,1.0,0.0 /)
     real(kind=dp) , dimension(4) , private , parameter :: probcasP3  = (/ 0.0,0.0,0.0,1.0 /)

     character(len=2) , dimension(4)        ,parameter  :: casGeno = (/'QQ','Qq','qQ','qq'/)

     ! alpha
     ! transformation
     ! cas 1,2,3,4
     !
     ! General
     ! ----------
     ! g1 QQ => l1
     ! g2 Qq => l2
     ! g3 qQ => l3
     ! g4 qq => l4
     !
     ! cor(1) => (1,0,0,0)
     ! cor(2) => (0,1,0,0)
     ! cor(3) => (0,0,1,0)
     ! cor(4) => (0,0,0,1)

     ! Dominance
     ! ----------
     ! g1 QQ => l1
     ! g2 Qq => l2
     ! g3 qQ => l2
     ! g4 qq => l3
     ! cor(1) => (1,0,0,0)
     ! cor(2) => (0,1,0,0)
     ! cor(3) => (0,1,0,0)
     ! cor(4) => (0,0,1,0)

     ! Additif
     ! ----------
     ! g1 QQ => l1
     ! g2 Qq => (l1+l2)/2
     ! g3 qQ => (l1+l2)/2
     ! g4 qq => l2
     ! cor(1) => (1,0,0,0)
     ! cor(2) => (0.5,0.5,0,0)
     ! cor(3) => (0.5,0.5,0,0)
     ! cor(4) => (0,1,0,0)

    real(kind=dp) , dimension(4) , private , target :: alpha1General = (/ 1.0,0.0,0.0,0.0 /)
    real(kind=dp) , dimension(4) , private , target :: alpha2General = (/ 0.0,1.0,0.0,0.0 /)
    real(kind=dp) , dimension(4) , private , target :: alpha3General = (/ 0.0,0.0,1.0,0.0 /)
    real(kind=dp) , dimension(4) , private , target :: alpha4General = (/ 0.0,0.0,0.0,1.0 /)

    real(kind=dp) , dimension(4) , private , target :: alpha1Dominance = (/ 1.0,0.0,0.0,0.0 /)
    real(kind=dp) , dimension(4) , private , target :: alpha2Dominance = (/ 0.0,1.0,1.0,0.0 /)
    real(kind=dp) , dimension(4) , private , target :: alpha3Dominance = (/ 0.0,0.0,0.0,1.0 /)
    real(kind=dp) , dimension(4) , private , target :: alpha4Dominance = (/ 0.0,0.0,0.0,0.0 /)

    real(kind=dp) , dimension(4) , private , target :: alpha1Additif = (/ 1.0,0.5,0.5,0.0 /)
    real(kind=dp) , dimension(4) , private , target :: alpha2Additif = (/ 0.0,0.5,0.5,1.0 /)
    real(kind=dp) , dimension(4) , private , target :: alpha3Additif = (/ 0.0,0.0,0.0,0.0 /)
    real(kind=dp) , dimension(4) , private , target :: alpha4Additif = (/ 0.0,0.0,0.0,0.0 /)


    real(kind=dp) , dimension(21) , private         :: valP = (/ 0.0,0.05,0.10,0.15,0.20,0.25,&
                                                                 0.30,0.35,0.40,0.45,0.50,&
                                                                 0.55,0.60,0.65,0.70,0.75,0.80,&
                                                                 0.85,0.90,0.95,1.00 /)

    ! Recuperer l'index d un alpha en fonction des CAS Pere et Mere
    integer       , dimension(4,4,4) , private         :: ia


    type (work_bicfarnir) , pointer, private :: dataw_p
    type (QTLMAP_DATASET) , pointer, private :: dataset_p
    type (PDD_BUILD)      , pointer, private :: spt_p
    type (work_pos)       , pointer, private :: wp_p


    contains

    subroutine setIA

     ! 1 -> QQ , 1-> QQ
    ia(1,1,1) =  1
    ia(1,1,2) =  1
    ia(1,1,3) =  1
    ia(1,1,4) =  1
    ! 1 -> QQ , 2-> Qq  =>
    ia(1,2,1) =  1
    ia(1,2,2) =  2
    ia(1,2,3) =  1
    ia(1,2,4) =  2
    ! 1 -> QQ , 3-> qQ  =>
    ia(1,3,1) =  2
    ia(1,3,2) =  1
    ia(1,3,3) =  2
    ia(1,3,4) =  1
    ! 1 -> QQ , 4-> qq  =>
    ia(1,4,1) =  2
    ia(1,4,2) =  2
    ia(1,4,3) =  2
    ia(1,4,4) =  2
    ! ==
    ! 2 -> Qq , 1-> QQ  =>
    ia(2,1,1) =  1
    ia(2,1,2) =  1
    ia(2,1,3) =  3
    ia(2,1,4) =  3
    ! 2 -> Qq , 2-> Qq  =>
    ia(2,2,1) =  1
    ia(2,2,2) =  2
    ia(2,2,3) =  3
    ia(2,2,4) =  4
    ! 2 -> Qq , 3-> qQ  =>
    ia(2,3,1) =  2
    ia(2,3,2) =  1
    ia(2,3,3) =  4
    ia(2,3,4) =  3
    ! 2 -> Qq , 4-> qq  =>
    ia(2,4,1) =  2
    ia(2,4,2) =  2
    ia(2,4,3) =  4
    ia(2,4,4) =  4
    ! ==
    ! 3 -> qQ , 1-> QQ  =>
    ia(3,1,1) =  3
    ia(3,1,2) =  3
    ia(3,1,3) =  1
    ia(3,1,4) =  1
    ! 3 -> qQ , 2-> Qq  =>
    ia(3,2,1) =  3
    ia(3,2,2) =  4
    ia(3,2,3) =  1
    ia(3,2,4) =  2
    ! 3 -> qQ , 3-> qQ  =>
    ia(3,3,1) =  4
    ia(3,3,2) =  3
    ia(3,3,3) =  2
    ia(3,3,4) =  1
    ! 3 -> qQ , 4-> qq  =>
    ia(3,4,1) =  4
    ia(3,4,2) =  4
    ia(3,4,3) =  2
    ia(3,4,4) =  2
    ! ==
    ! 4 -> qq , 1-> QQ  =>
    ia(4,1,1) =  3
    ia(4,1,2) =  3
    ia(4,1,3) =  3
    ia(4,1,4) =  3
    ! 4 -> qq , 2-> Qq  =>
    ia(4,2,1) =  3
    ia(4,2,2) =  4
    ia(4,2,3) =  3
    ia(4,2,4) =  4
    ! 4 -> qq , 3-> qQ  =>
    ia(4,3,1) =  4
    ia(4,3,2) =  3
    ia(4,3,3) =  4
    ia(4,3,4) =  3
    ! 4 -> qq , 4-> qq  =>
    ia(4,4,1) =  4
    ia(4,4,2) =  4
    ia(4,4,3) =  4
    ia(4,4,4) =  4

    end subroutine setIA


     subroutine init_work(dataset,spt,ic,dataw)
      type(QTLMAP_DATASET)       ,intent(in)         :: dataset
      type(PDD_BUILD)            ,intent(in)         :: spt
      integer                    ,intent(in)         :: ic
      type (work_bicfarnir)           ,intent(out)        :: dataw

      integer :: ip,jm,kd,ifem,maxnd,chr,ig,ndd
      type(GENEALOGY_BASE) , pointer :: dg
      real(kind=dp)    ,dimension(:,:,:), pointer :: y,cd
      integer          ,dimension(:,:)   ,pointer :: nd

      dg => dataset%genea
      dataw%ic = ic
      maxnd=0
      do ip=1,dg%np
       do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
        maxnd = max(maxnd,dg%ndm(jm+1)-dg%ndm(jm))
       end do
      end do

      allocate (dataw%ndsib(dg%np,dg%nm))
      allocate (dataw%ysib(dg%np,dg%nm,maxnd))
      allocate (dataw%cdsib(dg%np,dg%nm,maxnd))
      allocate (dataw%ndfull(dg%np,dg%nm))
      allocate (dataw%yfull(dg%np,dg%nm,maxnd))
      allocate (dataw%cdfull(dg%np,dg%nm,maxnd))
      allocate (dataw%damfull(dg%np,dg%nm))
      allocate (dataw%ngend(dataset%map%nchr,size(spt%ngend,2)))
      allocate (dataw%ndesc(dataset%map%nchr,maxval(spt%ngend(:,size(spt%ngend,2)))))
      allocate (dataw%npdd(dataset%map%nchr,maxval(spt%ngend(:,size(spt%ngend,2)))))

      dataw%ndfull=0
      dataw%ndsib = 0
      dataw%damfull = 0
      dataw%ysib=0
      dataw%cdsib = 0
      dataw%yfull = 0
      dataw%cdfull = 0

      dataw%ngend=0

      ifem=0
      do ip=1,dg%np
       do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
        if ( dataset%phenoAnimal%estime(ic,jm) ) then
         y => dataw%yfull
         cd => dataw%cdfull
         nd => dataw%ndfull
         ifem=ifem+1
         dataw%damfull(ip,jm) = ifem
    !     print *,"FULL"
        else
         y => dataw%ysib
         cd => dataw%cdsib
         nd => dataw%ndsib
      !   print *,"SIB"
        end if
        do kd=dg%ndm(jm)+1,dg%ndm(jm+1)
          if ( dataset%phenoAnimal%presentc(ic,kd)) then
           nd(ip,jm)=nd(ip,jm)+1
           y(ip,jm,nd(ip,jm))=dataset%phenoAnimal%y(ic,kd)
           cd(ip,jm,nd(ip,jm))=dataset%phenoAnimal%cd(ic,kd)
          end if
        end do
!
!       print *,nd(ip,jm)
!       print *,y(ip,jm,:nd(ip,jm))

       !! correspondance spt
      do chr=1,dataset%map%nchr
       do ig=spt%ngenom(chr,jm)+1,spt%ngenom(chr,jm+1)
         ndd=0
         dataw%ngend(chr,ig+1)=dataw%ngend(chr,ig)
         do kd=spt%ngend(chr,ig)+1,spt%ngend(chr,ig+1)
           if ( dataset%phenoAnimal%presentc(ic,spt%ndesc(chr,kd))) then
            !ngend correspond a l indice kd de y(ip,jm,kd)
            dataw%ngend(chr,ig+1)=dataw%ngend(chr,ig+1)+1
            ndd=ndd+1
            dataw%ndesc(chr,dataw%ngend(chr,ig+1)) = ndd
            dataw%npdd(chr,dataw%ngend(chr,ig+1)) = kd
           end if
         end do
       end do

!       do ig=spt%ngenom(chr,jm)+1,spt%ngenom(chr,jm+1)
!         print *,"ngend:",dataw%ngend(chr,ig)+1,dataw%ngend(chr,ig+1)
!         print *,"corr:",dataw%ndesc(chr,dataw%ngend(chr,ig)+1:dataw%ngend(chr,ig+1))
!         print *,"cpdd:",dataw%npdd(chr,dataw%ngend(chr,ig)+1:dataw%ngend(chr,ig+1))
!       end do

      end do

      end do!jm
      end do!ip

      dataw%ndamestim = ifem

      !test cas general
      dataw%alpha1 => alpha1General
      dataw%alpha2 => alpha2General
      dataw%alpha3 => alpha3General
      dataw%alpha4 => alpha4General


!      dataw%alpha1 => alpha1Additif
!      dataw%alpha2 => alpha2Additif
!      dataw%alpha3 => alpha3Additif
!      dataw%alpha4 => alpha4Additif


!      dataw%alpha1 => alpha1Dominance
!      dataw%alpha2 => alpha2Dominance
!      dataw%alpha3 => alpha3Dominance
!      dataw%alpha4 => alpha4Dominance


      call setIA

     end subroutine init_work

     subroutine opti_0qtl(dataset,dataw)
       type(QTLMAP_DATASET),target,intent(in)         :: dataset
       type (work_bicfarnir) ,target   ,intent(inout)    :: dataw
!
! Divers

      integer                                    :: iuser(1)
      integer        ,dimension(:),allocatable   :: iw
      real (kind=dp) ,dimension(:),allocatable   :: borni,borns,w
      real (kind=dp)                             :: user(1)
      real (kind=dp)                             :: f
      integer :: npar,ibound,ip,ix,ifail,ntot
      integer :: km,jm,kd,kd1,kd2,ii
      logical itest
      logical  , dimension(:,:,:),pointer        :: filter_inc
      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa
      type(DATAMODEL_BASE) , pointer :: dpm

      real (kind=dp) ,dimension(dataset%genea%np) :: fp0
      real (kind=dp) ,dimension(dataset%genea%nm) :: fm0

      dpm => dataset%phenoModel
      dg => dataset%genea
      dpa => dataset%phenoAnimal

      !nombre de parametre : variances heteroscedastique + mu famille
      npar=2*dg%np+dataw%ndamestim

      allocate (borni(npar))
      allocate (borns(npar))
      allocate (dataw%parSolH0(npar))

      ! OPMIZATION - QUASINEWTON
      allocate (filter_inc(dg%np,dg%nm, npar))
      filter_inc=.false.
      do ip=1,dg%np
        do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
         filter_inc(ip,jm,ip) = .true.
         filter_inc(ip,jm,dg%np+ip) = .true.
         if ( dataw%ndfull(ip,jm)>0 ) then
          filter_inc(ip,jm,2*dg%np+dataw%damfull(ip,jm)) = .true.
         end if
        end do
      end do

      ibound=0
      do ip=1,dg%np
        borni(ip)=SIG_MIN
        borns(ip)=SIG_MAX
      end do
      do ix=dg%np+1,npar
        borni(ix)=XMU_MIN
        borns(ix)=XMU_MAX
      end do

!
! Point de depart

      do ip=1,dg%np
        !Polygenic init
        dataw%parSolH0(dg%np+ip)=0.d0
        ntot=0.d0
        do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
         if ( dataw%damfull(ip,jm)>0) then
          dataw%parSolH0(2*dg%np+dataw%damfull(ip,jm)) = sum(dataw%yfull(ip,jm,:dataw%ndfull(ip,jm)))/real(dataw%ndfull(ip,jm))
         end if
         dataw%parSolH0(dg%np+ip)=dataw%parSolH0(dg%np+ip)+sum(dataw%yfull(ip,jm,:dataw%ndfull(ip,jm)))+&
                                     sum(dataw%ysib(ip,jm,:dataw%ndsib(ip,jm)))
         ntot=ntot+dataw%ndfull(ip,jm)+dataw%ndsib(ip,jm)
        end do
        dataw%parSolH0(dg%np+ip)=dataw%parSolH0(dg%np+ip)/real(ntot)
        !variance
        dataw%parSolH0(ip)=1.d0 ! a faire...
      end do

! Optimisation de la vraisemblance
      ifail=1

      dataw_p => dataw
      dataset_p => dataset

      call minimizing_funct_family(dataset,npar,ibound,funct_0qtl_fam,&
       filter_inc,fm0,fp0,borni,borns,dataw%parSolH0,f,iuser,user,ifail)

      print *,'f0:',f
      print *,"par0:",dataw%parSolH0
      dataw%f0=f


      deallocate (borni)
      deallocate (borns)
     ! deallocate (par)
      deallocate (filter_inc)


     end subroutine opti_0qtl


     subroutine funct_0qtl_fam(ip,jm,n,x,f,iuser,user)
       integer         , intent(in)                  :: ip,jm,n
       real (kind=dp)      ,dimension(n), intent(in) :: x
       real (kind=dp)  , intent(inout)               :: f
       integer ,       dimension(1), intent(in)      :: iuser
       real (kind=dp)      ,dimension(1), intent(in) :: user

       real (kind=dp) :: vpf,v,vmere
       integer        :: kkd

       vpf=0.d0
       ! Si c est une famille de plein frere(ndmin suffisemment grand), on passe dans cette boucle
       do kkd=1,dataw_p%ndfull(ip,jm)
        v=dataw_p%yfull(ip,jm,kkd)
        ! Polygenic sire
        v=v-x(dataset_p%genea%np+ip)
        ! Polygenic dam
        v=v-x(2*dataset_p%genea%np+dataw_p%damfull(ip,jm))
        !
        vpf=vpf+v*v*dataw_p%cdfull(ip,jm,kkd)
       end do

       ! Si la mere est non estimable (ndmin trop grand), on passe dans cette boucle
       do kkd=1,dataw_p%ndsib(ip,jm)
        v=dataw_p%ysib(ip,jm,kkd)
        ! Polygenic sire
        v=v-x(dataset_p%genea%np+ip)
        vpf=vpf+v*v*dataw_p%cdsib(ip,jm,kkd)
       end do

       vmere=dexp(-0.5d0*vpf/(x(ip)*x(ip)))
       f=-log(vmere)+real(max(dataw_p%ndfull(ip,jm),dataw_p%ndsib(ip,jm)))*log(x(ip))

      end subroutine funct_0qtl_fam


     subroutine opti_1qtl(dataset,spt,dataw,lrtsol)

      type(QTLMAP_DATASET)  ,target       ,intent(in)     :: dataset
      type(PDD_BUILD)       ,target       ,intent(in)     :: spt
      type (work_bicfarnir)      ,target   ,intent(inout)    :: dataw
      type(TYPE_LRT_SOLUTION)  , intent(out)              :: lrtsol

      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa
      type(DATAMODEL_BASE) , pointer :: dpm

      real (kind=dp) ,dimension(dataset%genea%np) :: fp1
      real (kind=dp) ,dimension(dataset%genea%nm) :: fm1
      integer                                    :: iuser(1)
      real (kind=dp) ,dimension(:),allocatable   :: par,borni,borns
      real (kind=dp)                             :: user(2)
      logical  , dimension(:,:),pointer        :: filter_inc

      integer :: npar,ifail,ip,jm,icas,ibound,chr,ipos,n,ip1,ip2,sc
      real (kind=dp) :: f1,sumUnp1p2,p1,p2,p3
      real (kind=dp) ,dimension(:,:,:) ,allocatable :: parall
      real (kind=dp) ,dimension(:,:) ,allocatable :: fmax
      type(work_pos) ,target :: work
      real (kind=dp) ::f1max,p1max,p2max

      dpm => dataset%phenoModel
      dg => dataset%genea
      dpa => dataset%phenoAnimal

      call lrtsol%new(dataset,1)

      !nombre de parametre :
      ! np                    : variances heteroscedastique
      ! np                    : polygenic sire
      ! maxval(dataw%damfull) : polygenic sire-dam
      ! 2                     : P1 + P2
      ! 4                     : l1,l2,l3,l4
      npar=2*dg%np+1

      sc=0

      allocate (borni(npar))
      allocate (borns(npar))
      allocate (par(npar))
      allocate (parall(npar,dataset%map%nchr,dataset%map%get_maxnpo()))
      allocate (fmax(dataset%map%nchr,dataset%map%get_maxnpo()))
      allocate (dataw%parSolH1(npar))



      ! OPMIZATION - QUASINEWTON
      allocate (filter_inc(dg%np,npar))
      filter_inc=.true.
      !filter_inc(:,dg%np+1) = .false.

      borni(:dg%np)=SIG_MIN
      borns(:dg%np)=SIG_MAX

      borni(dg%np:)=XMU_MIN
      borns(dg%np:)=XMU_MAX

!      borni(dg%np+1:dg%np+1)=0.001d0
!      borns(dg%np+1:dg%np+1)=0.999d0

      borni(2*dg%np+1+sc:2*dg%np+1+sc)=-99999.d0
      borns(2*dg%np+1+sc:2*dg%np+1+sc)=99999.d0

      !
! Point de depart
      par=0.d0
      do ip=1,dg%np
        !var
        par(ip)=dataw%parSolH0(ip)
        !Polygenic init
        par(dg%np+ip)=dataw%parSolH0(dg%np+ip)
      end do

    !  par(dg%np+1)=0.5d0
 !     par(dg%np+2)=2d0

      par(2*dg%np+1)=0d0
!      par(2*dg%np+2+dataw%ndamestim+2)=0d0
!      par(2*dg%np+2+dataw%ndamestim+3)=0d0
!      par(2*dg%np+2+dataw%ndamestim+4)=0d0

      ibound=0
      ifail=1

      allocate (work%n(1),work%chr(1))
      work%n=1
      work%chr=1

      dataw_p => dataw
      dataset_p => dataset
      spt_p => spt
      wp_p => work
      chr=1
      do ipos=1,dataset%map%get_npo(1)
       work%n=ipos

        f1max = 99999
        do ip1=1,21
        user(1)=valP(ip1)

        call minimizing_funct_family_sire(dataset,npar,ibound,funct_1qtl_fam_halfsib_add,filter_inc,&
          fp1,borni,borns,par,f1,iuser,user,ifail)

          if ( f1max>f1 ) then
              f1max=f1
              p1max=user(1)
              p2max=user(2)
              print *,work%n(1),ip1, user(1),ip1,f1
           end if

        end do

      !  print *,par(dg%np+1),par(2*dg%np+1+dataw%ndamestim+1)
      !  f1max=f1

!       f1max = dataw%f0
!       do ip1=1,21
!         do ip2=1,21
!
!           if ( (valP(ip1)+2*valP(ip2))<=1.0) then
!           user(1)=valP(ip1)
!           user(2)=valP(ip2)
!           print *,user(1),user(2)
!          call minimizing_funct_family_sire(dataset,npar,ibound,funct_1qtl_fam_halfsib,filter_inc,&
!          fp1,borni,borns,par,f1,iuser,user,ifail)
!
!           !print *,user(1),user(2),f1
!           if ( f1max<f1 ) then
!              f1max=f1
!              p1max=user(1)
!              p2max=user(2)
!              print *,work%n(1),ip1, user(1),ip2,user(2),f1
!           end if
!           end if
!         end do
!      end do

      parall(:npar,work%chr(1),work%n(1)) = par(:npar)
      fmax(work%chr(1),work%n(1)) = f1max

       call lrtsol%LRT%add(dataset,1,work%chr,work%n,(-2.d0*(f1max-dataw%f0)),1)

     !  stop
     !  if (ipos==5) stop
      end do

      lrtsol%lrtmax(0) = 0.d0

      do chr=1,dataset%map%nchr
        work%chr(1) = chr
        do n=1,dataset%map%get_npo(chr)
          work%n(1) = n
          if (lrtsol%LRT%get(dataset,1,work%chr,work%n,1)> lrtsol%lrtmax(0)) then
             dataw%parSolH1=parall(:,chr,n)
             dataw%f1 = fmax(chr,n)
             lrtsol%lrtmax(0)=lrtsol%LRT%get(dataset,1,work%chr,work%n,1)
             lrtsol%nxmax(0)=n
             lrtsol%chrmax(0)=chr
          end if
        end do
      end do

      deallocate (borni)
      deallocate (borns)
      deallocate (par)
      deallocate (parall)
      deallocate (fmax)

     end subroutine opti_1qtl

!  QQ -> 1
!  Qq -> 2
!  qQ -> 3
!  qq -> 4
!
!  gi1 => 1er allele du genotype du pere i
!  gi2 => 2eme allele du genotype du pere i
!
!  gij1 => 1er allele du genotype de la mere ij
!  gij2 => 2eme allele du genotype de la mere ij
!
!
!                      Sire,Mere
!  A           =    alpha(1,1)  =  1 => alpha1 => indice dans X du parametre alpha
!   gi1, gij1
!
     subroutine funct_1qtl_fam(ip,n,x,f,iuser,user)
       integer         , intent(in)                  :: ip,n
       real (kind=dp)      ,dimension(n), intent(in) :: x
       real (kind=dp)  , intent(inout)               :: f
       integer ,       dimension(1), intent(in)      :: iuser
       real (kind=dp)      ,dimension(1), intent(in) :: user

       real (kind=dp) :: vpere,vpf,v,vmere,sig,p1,p2,p3,var,fp,fm,fmtot
       real (kind=dp) :: sumUnp1p2,invsig,u
       integer        :: jm,kd,ig,casSire,casDam,kdid,kdpdd,i,effm

       sig = x(ip)
       invsig = (1.d0/sig)
       var = sig*sig
     !  sumUnp1p2=real(1.d0 + x(dataset_p%genea%np+1) + 2*x(dataset_p%genea%np+2))

       p1 = x(dataset_p%genea%np+1)*x(dataset_p%genea%np+1)! / sumUnp1p2
    !   p2 = x(dataset_p%genea%np+2) / sumUnp1p2
    !   p3 = 1.d0 / sumUnp1p2

       p2 = x(dataset_p%genea%np+1)*(1-x(dataset_p%genea%np+1))
       p3 = (1-x(dataset_p%genea%np+1))*(1-x(dataset_p%genea%np+1))

       fp=0
       ! QQ -> 1, Qq -> 2, qQ -> 3, qq -> 4
       do casSire=1,4
        fmtot=1.d0
        do jm=dataset_p%genea%nmp(ip)+1,dataset_p%genea%nmp(ip+1)
         fm=0.d0
         do casDam=1,4
          vmere=0.d0
          do ig=spt_p%ngenom(wp_p%chr(1),jm)+1,spt_p%ngenom(wp_p%chr(1),jm+1)
           vpf=0.d0

           !
           ! FULL
           !
          if ( dataw_p%damfull(ip,jm) > 0 ) then
           do kd=dataw_p%ngend(wp_p%chr(1),ig)+1,dataw_p%ngend(wp_p%chr(1),ig+1)
            kdid = dataw_p%ndesc(wp_p%chr(1),kd)
            kdpdd = dataw_p%npdd(wp_p%chr(1),kd)
            v=dataw_p%yfull(ip,jm,kdid)
            ! Polygenic sire
            v = v - x(dataset_p%genea%np+1+ip)
            ! Polygenic dam
            v=v-x(2*dataset_p%genea%np+1+dataw_p%damfull(ip,jm))
            !u=0
            !1|1 , 1|2, 2|1, 2|2
            do i=1,4
             v = v - spt_p%pdd(wp_p%chr(1),kdpdd,i,wp_p%n(1))*(&
              dataw_p%alpha1(ia(casSire,casDam,i))*x(2*dataset_p%genea%np+1+dataw_p%ndamestim+1)+&
              dataw_p%alpha2(ia(casSire,casDam,i))*x(2*dataset_p%genea%np+1+dataw_p%ndamestim+2)+&
              dataw_p%alpha3(ia(casSire,casDam,i))*x(2*dataset_p%genea%np+1+dataw_p%ndamestim+3)+&
              dataw_p%alpha4(ia(casSire,casDam,i))*x(2*dataset_p%genea%np+1+dataw_p%ndamestim+4))
            end do !i

            vpf=vpf+v*v*dataw_p%cdfull(ip,jm,kdid)
           end do !kd
         else
           do kd=dataw_p%ngend(wp_p%chr(1),ig)+1,dataw_p%ngend(wp_p%chr(1),ig+1)
            kdid = dataw_p%ndesc(wp_p%chr(1),kd)
            kdpdd = dataw_p%npdd(wp_p%chr(1),kd)
            v=dataw_p%ysib(ip,jm,kdid)
           ! print *,v,kdid
            ! Polygenic sire
            v = v - x(dataset_p%genea%np+1+ip)
            !1|1 , 1|2, 2|1, 2|2
            !u=0.d0
            do i=1,4
             v = v - spt_p%pdd(wp_p%chr(1),kdpdd,i,wp_p%n(1))*(&
              dataw_p%alpha1(ia(casSire,casDam,i))*x(2*dataset_p%genea%np+1+dataw_p%ndamestim+1)+&
              dataw_p%alpha2(ia(casSire,casDam,i))*x(2*dataset_p%genea%np+1+dataw_p%ndamestim+2)+&
              dataw_p%alpha3(ia(casSire,casDam,i))*x(2*dataset_p%genea%np+1+dataw_p%ndamestim+3)+&
              dataw_p%alpha4(ia(casSire,casDam,i))*x(2*dataset_p%genea%np+1+dataw_p%ndamestim+4))
            end do !i
            vpf=vpf+v*v*dataw_p%cdsib(ip,jm,kdid)
           end do !kd
         end if

           vmere=vmere+spt_p%probg(wp_p%chr(1),ig)*exp(-0.5d0*vpf/var)
          end do !ig

          fm=(p1*probcasP1(casDam)+p2*probcasP2(casDam)+p3*probcasP3(casDam))*vmere+fm
          !print *,jm,fm
         end do  !casDam
         effm=max(dataw_p%ndsib(ip,jm),dataw_p%ndfull(ip,jm))
         do i=1,effm
          fm=fm*invsig
         end do
         fmtot=fmtot*fm
       !  print *,'fmtot:',jm,fmtot
       end do !jm
          fp = (p1*probcasP1(casSire)+p2*probcasP2(casSire)+p3*probcasP3(casSire))*fmtot+fp
         !print *,'fp:',fp

       end do ! casSire
       if ( fp == 0 ) then
        f=INIFINY_REAL_VALUE
       else
        f = -dlog(fp)
       end if
   !   print *,p1,p2,f
      ! stop
      end subroutine funct_1qtl_fam

      subroutine funct_1qtl_fam2(ip,n,x,f,iuser,user)
       integer         , intent(in)                  :: ip,n
       real (kind=dp)      ,dimension(n), intent(in) :: x
       real (kind=dp)  , intent(inout)               :: f
       integer ,       dimension(1), intent(in)      :: iuser
       real (kind=dp)      ,dimension(1), intent(in) :: user

       real (kind=dp) :: vpere,vpf,v,vmere,sig,p1,p2,p3,var,fp,fm,fmtot
       real (kind=dp) :: sumUnp1p2,invsig,u
       integer        :: jm,kd,ig,casSire,casDam,i,effm,kkd,idfem

       sig = x(ip)
       invsig = (1.d0/sig)
       var = sig*sig
       sumUnp1p2=real(1.d0 + x(dataset_p%genea%np+1) + 2.d0*x(dataset_p%genea%np+2))
       p1 = x(dataset_p%genea%np+1) / sumUnp1p2
       p2 = x(dataset_p%genea%np+2) / sumUnp1p2
       p3 = 1.d0 / sumUnp1p2

       idfem=0
       fp=0
       ! QQ -> 1, Qq -> 2, qQ -> 3, qq -> 4
       do casSire=1,4
        fmtot=1.d0
        do jm=dataset_p%genea%nmp(ip)+1,dataset_p%genea%nmp(ip+1)
         if(dataset_p%phenoAnimal%estime(1,jm)) then
          idfem=idfem+1
         end if
         fm=0.d0
         do casDam=1,4
          vmere=0.d0
          do ig=spt_p%ngenom(wp_p%chr(1),jm)+1,spt_p%ngenom(wp_p%chr(1),jm+1)
           vpf=0.d0
           effm=0
           do kd=spt_p%ngend(wp_p%chr(1),ig)+1,spt_p%ngend(wp_p%chr(1),ig+1)
            kkd=spt_p%ndesc(wp_p%chr(1),kd)
            if(dataset_p%phenoAnimal%presentc(1,kkd)) then
            effm=effm+1
            v=dataset_p%phenoAnimal%y(1,kkd)
            ! Polygenic sire
            v = v - x(dataset_p%genea%np+2+ip)
            ! poly dam
            if(dataset_p%phenoAnimal%estime(1,jm)) then
             v=v-x(2*dataset_p%genea%np+2+idfem)
            end if
            !General
            !QQ QQ
            if ( casSire == 1 .and. casDam == 1 ) then
             v = v - spt_p%pdd(wp_p%chr(1),kd,1,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+1)
             v = v - spt_p%pdd(wp_p%chr(1),kd,2,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+1)
             v = v - spt_p%pdd(wp_p%chr(1),kd,3,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+1)
             v = v - spt_p%pdd(wp_p%chr(1),kd,4,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+1)
            !QQ Qq
            else if (casSire == 1 .and. casDam == 2 ) then
             v = v - spt_p%pdd(wp_p%chr(1),kd,1,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+1)
             v = v - spt_p%pdd(wp_p%chr(1),kd,2,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+2)
             v = v - spt_p%pdd(wp_p%chr(1),kd,3,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+1)
             v = v - spt_p%pdd(wp_p%chr(1),kd,4,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+2)
            !QQ qQ
            else if (casSire == 1 .and. casDam == 3 ) then
             v = v - spt_p%pdd(wp_p%chr(1),kd,1,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+2)
             v = v - spt_p%pdd(wp_p%chr(1),kd,2,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+1)
             v = v - spt_p%pdd(wp_p%chr(1),kd,3,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+2)
             v = v - spt_p%pdd(wp_p%chr(1),kd,4,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+1)
            !QQ qq
            else if (casSire == 1 .and. casDam == 4 ) then
             v = v - spt_p%pdd(wp_p%chr(1),kd,1,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+2)
             v = v - spt_p%pdd(wp_p%chr(1),kd,2,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+2)
             v = v - spt_p%pdd(wp_p%chr(1),kd,3,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+2)
             v = v - spt_p%pdd(wp_p%chr(1),kd,4,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+2)

            !Qq QQ
            else if (casSire == 2 .and. casDam == 1 ) then
             v = v - spt_p%pdd(wp_p%chr(1),kd,1,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+1)
             v = v - spt_p%pdd(wp_p%chr(1),kd,2,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+1)
             v = v - spt_p%pdd(wp_p%chr(1),kd,3,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+3)
             v = v - spt_p%pdd(wp_p%chr(1),kd,4,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+3)
            !Qq Qq
            else if (casSire == 2 .and. casDam == 2 ) then
             v = v - spt_p%pdd(wp_p%chr(1),kd,1,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+1)
             v = v - spt_p%pdd(wp_p%chr(1),kd,2,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+2)
             v = v - spt_p%pdd(wp_p%chr(1),kd,3,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+3)
             v = v - spt_p%pdd(wp_p%chr(1),kd,4,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+4)

            !Qq qQ
            else if (casSire == 2 .and. casDam == 3 ) then
             v = v - spt_p%pdd(wp_p%chr(1),kd,1,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+2)
             v = v - spt_p%pdd(wp_p%chr(1),kd,2,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+1)
             v = v - spt_p%pdd(wp_p%chr(1),kd,3,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+4)
             v = v - spt_p%pdd(wp_p%chr(1),kd,4,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+3)
            !Qq qq
            else if (casSire == 2 .and. casDam == 4 ) then
             v = v - spt_p%pdd(wp_p%chr(1),kd,1,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+2)
             v = v - spt_p%pdd(wp_p%chr(1),kd,2,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+2)
             v = v - spt_p%pdd(wp_p%chr(1),kd,3,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+4)
             v = v - spt_p%pdd(wp_p%chr(1),kd,4,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+4)
            !qQ QQ
            else if (casSire == 3 .and. casDam == 1 ) then
             v = v - spt_p%pdd(wp_p%chr(1),kd,1,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+3)
             v = v - spt_p%pdd(wp_p%chr(1),kd,2,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+3)
             v = v - spt_p%pdd(wp_p%chr(1),kd,3,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+1)
             v = v - spt_p%pdd(wp_p%chr(1),kd,4,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+1)
            !qQ Qq
            else if (casSire == 3 .and. casDam == 2 ) then
             v = v - spt_p%pdd(wp_p%chr(1),kd,1,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+3)
             v = v - spt_p%pdd(wp_p%chr(1),kd,2,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+4)
             v = v - spt_p%pdd(wp_p%chr(1),kd,3,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+1)
             v = v - spt_p%pdd(wp_p%chr(1),kd,4,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+2)
            !qQ qQ
            else if (casSire == 3 .and. casDam == 3 ) then
             v = v - spt_p%pdd(wp_p%chr(1),kd,1,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+4)
             v = v - spt_p%pdd(wp_p%chr(1),kd,2,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+3)
             v = v - spt_p%pdd(wp_p%chr(1),kd,3,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+2)
             v = v - spt_p%pdd(wp_p%chr(1),kd,4,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+1)
            !qQ qq
            else if (casSire == 3 .and. casDam == 4 ) then
             v = v - spt_p%pdd(wp_p%chr(1),kd,1,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+4)
             v = v - spt_p%pdd(wp_p%chr(1),kd,2,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+4)
             v = v - spt_p%pdd(wp_p%chr(1),kd,3,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+2)
             v = v - spt_p%pdd(wp_p%chr(1),kd,4,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+2)
            !qq QQ
            else if (casSire == 4 .and. casDam == 1 ) then
             v = v - spt_p%pdd(wp_p%chr(1),kd,1,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+3)
             v = v - spt_p%pdd(wp_p%chr(1),kd,2,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+3)
             v = v - spt_p%pdd(wp_p%chr(1),kd,3,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+3)
             v = v - spt_p%pdd(wp_p%chr(1),kd,4,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+3)
            !qq Qq
            else if (casSire == 4 .and. casDam == 2 ) then
             v = v - spt_p%pdd(wp_p%chr(1),kd,1,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+3)
             v = v - spt_p%pdd(wp_p%chr(1),kd,2,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+4)
             v = v - spt_p%pdd(wp_p%chr(1),kd,3,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+3)
             v = v - spt_p%pdd(wp_p%chr(1),kd,4,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+4)
            !qq qQ
            else if (casSire == 4 .and. casDam == 3 ) then
             v = v - spt_p%pdd(wp_p%chr(1),kd,1,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+4)
             v = v - spt_p%pdd(wp_p%chr(1),kd,2,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+3)
             v = v - spt_p%pdd(wp_p%chr(1),kd,3,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+4)
             v = v - spt_p%pdd(wp_p%chr(1),kd,4,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+3)
            !qq qq
            else if (casSire == 4 .and. casDam == 4 ) then
             v = v - spt_p%pdd(wp_p%chr(1),kd,1,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+4)
             v = v - spt_p%pdd(wp_p%chr(1),kd,2,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+4)
             v = v - spt_p%pdd(wp_p%chr(1),kd,3,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+4)
             v = v - spt_p%pdd(wp_p%chr(1),kd,4,wp_p%n(1))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+4)
            end if

            vpf=vpf+v*v*dataset_p%phenoAnimal%cd(1,kkd)
            end if
           end do !kd
           vmere=vmere+spt_p%probg(wp_p%chr(1),ig)*exp(-0.5d0*vpf/var)
          end do !ig
          fm=(p1*probcasP1(casDam)+p2*probcasP2(casDam)+p3*probcasP3(casDam))*vmere+fm
         end do  !casDam
         do i=1,effm
          fm=fm*invsig
         end do
         fmtot=fmtot*fm
       end do !jm
       fp = (p1*probcasP1(casSire)+p2*probcasP2(casSire)+p3*probcasP3(casSire))*fmtot+fp

       end do ! casSire
       if ( fp == 0 ) then
        f=INIFINY_REAL_VALUE
       else
        f = -dlog(fp)
       end if
      ! print *,p1,p2,f
      ! stop
      end subroutine funct_1qtl_fam2


     subroutine funct_1qtl_fam_halfsib(ip,n,x,f,iuser,user)
       integer         , intent(in)                  :: ip,n
       real (kind=dp)      ,dimension(n), intent(in) :: x
       real (kind=dp)  , intent(inout)               :: f
       integer ,       dimension(1), intent(in)      :: iuser
       real (kind=dp)      ,dimension(2), intent(in) :: user

       real (kind=dp) :: vpere,vpf,v,vmere,sig,p1,p2,p3,var,fp,fm,fmtot
       real (kind=dp) :: sumUnp1p2,invsig,u
       integer        :: jm,kd,ig,casSire,casDam,i,effm,kkd,idfem

       sig = x(ip)
       invsig = (1.d0/sig)
       var = sig*sig
       sumUnp1p2=real(1.d0 + x(dataset_p%genea%np+1) + 2.d0*x(dataset_p%genea%np+2))
       p1 = x(dataset_p%genea%np+1) / sumUnp1p2
       p2 = x(dataset_p%genea%np+2) / sumUnp1p2
       p3 = 1.d0 / sumUnp1p2

!        p1 = user(1)
!        p2 = user(2)
!        p3 = 1 -p1 -2.d0*p2
!
!       p1 = 0
!       p2 = 0
!       p3 = 1

       idfem=0
       fp=0
       ! QQ -> 1, Qq -> 2, qQ -> 3, qq -> 4
       do casSire=1,4
        fmtot=1.d0
        do jm=dataset_p%genea%nmp(ip)+1,dataset_p%genea%nmp(ip+1)
         if(dataset_p%phenoAnimal%estime(1,jm)) then
          idfem=idfem+1
         end if
         fm=0.d0
         do casDam=1,4
          vmere=0.d0
          do ig=spt_p%ngenom(wp_p%chr(1),jm)+1,spt_p%ngenom(wp_p%chr(1),jm+1)
           vpf=0.d0
           effm=0
           do kd=spt_p%ngend(wp_p%chr(1),ig)+1,spt_p%ngend(wp_p%chr(1),ig+1)
            kkd=spt_p%ndesc(wp_p%chr(1),kd)
            if(dataset_p%phenoAnimal%presentc(1,kkd)) then
            effm=effm+1
            v=dataset_p%phenoAnimal%y(1,kkd)
            ! Polygenic sire
            v = v - x(dataset_p%genea%np+2+ip)

            !General
            !QQ QQ
            if ( casSire == 1 ) then
             if ( casDam == 1 ) then
               v = v - x(2*dataset_p%genea%np+2+dataw_p%ndamestim+1)
             else if ( casDam == 2 .or. casDam == 3 ) then
               v = v - (x(2*dataset_p%genea%np+2+dataw_p%ndamestim+1)+x(2*dataset_p%genea%np+2+dataw_p%ndamestim+2))*0.5d0
             else if ( casDam == 4 ) then
              v = v - x(2*dataset_p%genea%np+2+dataw_p%ndamestim+2)
             end if
            !QQ Qq
            else if (casSire == 2 ) then
             if ( casDam == 1 ) then
               v = v - (spt_p%pdd(wp_p%chr(1),kd,1,wp_p%n(1))+&
                 spt_p%pdd(wp_p%chr(1),kd,2,wp_p%n(1)))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+1) &
               - (spt_p%pdd(wp_p%chr(1),kd,3,wp_p%n(1))+&
                 spt_p%pdd(wp_p%chr(1),kd,4,wp_p%n(1)))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+3)
             else if ( casDam == 2 .or. casDam == 3 ) then
               v = v - (spt_p%pdd(wp_p%chr(1),kd,1,wp_p%n(1))+spt_p%pdd(wp_p%chr(1),kd,2,wp_p%n(1)))*&
                       (x(2*dataset_p%genea%np+2+dataw_p%ndamestim+1)+x(2*dataset_p%genea%np+2+dataw_p%ndamestim+2))*0.5d0 &
                     - (spt_p%pdd(wp_p%chr(1),kd,3,wp_p%n(1))+spt_p%pdd(wp_p%chr(1),kd,4,wp_p%n(1)))*&
                        (x(2*dataset_p%genea%np+2+dataw_p%ndamestim+3)+x(2*dataset_p%genea%np+2+dataw_p%ndamestim+4))*0.5d0
             else if ( casDam == 4 ) then
              v = v - (spt_p%pdd(wp_p%chr(1),kd,1,wp_p%n(1))+&
               spt_p%pdd(wp_p%chr(1),kd,2,wp_p%n(1)))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+2) &
             - (spt_p%pdd(wp_p%chr(1),kd,3,wp_p%n(1))+&
             spt_p%pdd(wp_p%chr(1),kd,4,wp_p%n(1)))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+4)
             end if
            !QQ qQ
            else if (casSire == 3 ) then
            if ( casDam == 1 ) then
               v = v - (spt_p%pdd(wp_p%chr(1),kd,1,wp_p%n(1))+&
                spt_p%pdd(wp_p%chr(1),kd,2,wp_p%n(1)))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+3) &
             - (spt_p%pdd(wp_p%chr(1),kd,3,wp_p%n(1))+&
             spt_p%pdd(wp_p%chr(1),kd,4,wp_p%n(1)))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+1)
             else if ( casDam == 2 .or. casDam == 3 ) then
               v = v - (spt_p%pdd(wp_p%chr(1),kd,1,wp_p%n(1))+spt_p%pdd(wp_p%chr(1),kd,2,wp_p%n(1)))*&
                       (x(2*dataset_p%genea%np+2+dataw_p%ndamestim+3)+x(2*dataset_p%genea%np+2+dataw_p%ndamestim+4))*0.5d0 &
                     - (spt_p%pdd(wp_p%chr(1),kd,3,wp_p%n(1))+spt_p%pdd(wp_p%chr(1),kd,4,wp_p%n(1)))*&
                        (x(2*dataset_p%genea%np+2+dataw_p%ndamestim+1)+x(2*dataset_p%genea%np+2+dataw_p%ndamestim+2))*0.5d0
             else if ( casDam == 4 ) then
              v = v - (spt_p%pdd(wp_p%chr(1),kd,1,wp_p%n(1))+&
               spt_p%pdd(wp_p%chr(1),kd,2,wp_p%n(1)))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+4) &
            - (spt_p%pdd(wp_p%chr(1),kd,3,wp_p%n(1))+&
            spt_p%pdd(wp_p%chr(1),kd,4,wp_p%n(1)))*x(2*dataset_p%genea%np+2+dataw_p%ndamestim+2)
             end if
            !QQ qq
            else if (casSire == 4 ) then
            if ( casDam == 1 ) then
               v = v - x(2*dataset_p%genea%np+2+dataw_p%ndamestim+3)
             else if ( casDam == 2 .or. casDam == 3 ) then
               v = v - (x(2*dataset_p%genea%np+2+dataw_p%ndamestim+3)+x(2*dataset_p%genea%np+2+dataw_p%ndamestim+4))*0.5d0
             else if ( casDam == 4 ) then
              v = v - x(2*dataset_p%genea%np+2+dataw_p%ndamestim+4)
             end if
            end if

            vpf=vpf+v*v*dataset_p%phenoAnimal%cd(1,kkd)
            end if
           end do !kd
           vmere=vmere+spt_p%probg(wp_p%chr(1),ig)*exp(-0.5d0*vpf/var)
          end do !ig
          fm=(p1*probcasP1(casDam)+p2*probcasP2(casDam)+p3*probcasP3(casDam))*vmere+fm
         end do  !casDam
         do i=1,effm
          fm=fm*invsig
         end do
         fmtot=fmtot*fm
       end do !jm
       fp = (p1*probcasP1(casSire)+p2*probcasP2(casSire)+p3*probcasP3(casSire))*fmtot+fp

       end do ! casSire
       if ( fp == 0 ) then
        f=INIFINY_REAL_VALUE
       else
        f = -dlog(fp)
       end if
      ! print *,p1,p2,f
      ! stop
      end subroutine funct_1qtl_fam_halfsib


       subroutine funct_1qtl_fam_halfsib_add(ip,n,x,f,iuser,user)
       integer         , intent(in)                  :: ip,n
       real (kind=dp)      ,dimension(n), intent(in) :: x
       real (kind=dp)  , intent(inout)               :: f
       integer ,       dimension(1), intent(in)      :: iuser
       real (kind=dp)      ,dimension(2), intent(in) :: user

       real (kind=dp) :: vpere,vpf,v,vmere,sig,p1,p2,p3,var,fp,fm,fmtot
       real (kind=dp) :: sumUnp1p2,invsig,u,p
       integer        :: jm,kd,ig,casSire,casDam,i,effm,kkd,idfem

       sig = x(ip)
       invsig = (1.d0/sig)
       var = sig*sig
!       sumUnp1p2=real(1.d0 + x(dataset_p%genea%np+1) + 2.d0*x(dataset_p%genea%np+2))
!       p1 = x(dataset_p%genea%np+1) / sumUnp1p2
!       p2 = x(dataset_p%genea%np+2) / sumUnp1p2
!       p3 = 1.d0 / sumUnp1p2

      ! p=x(dataset_p%genea%np+1)
       p=user(1)
      ! print *,p
       p1 = p*p! / sumUnp1p2
       p2 = p*(1-p)
       p3 = (1-p)*(1-p)

!        p1 = user(1)
!        p2 = user(2)
!        p3 = 1 -p1 -2.d0*p2
!
!       p1 = 0
!       p2 = 0
!       p3 = 1

       idfem=0
       fp=0
       ! QQ -> 1, Qq -> 2, qQ -> 3, qq -> 4
       do casSire=1,4
        fmtot=1.d0
        do jm=dataset_p%genea%nmp(ip)+1,dataset_p%genea%nmp(ip+1)
         if(dataset_p%phenoAnimal%estime(1,jm)) then
          idfem=idfem+1
         end if
         fm=0.d0
         do casDam=1,4
          vmere=0.d0
          do ig=spt_p%ngenom(wp_p%chr(1),jm)+1,spt_p%ngenom(wp_p%chr(1),jm+1)
           vpf=0.d0
           effm=0
           do kd=spt_p%ngend(wp_p%chr(1),ig)+1,spt_p%ngend(wp_p%chr(1),ig+1)
            kkd=spt_p%ndesc(wp_p%chr(1),kd)
            if(dataset_p%phenoAnimal%presentc(1,kkd)) then
            effm=effm+1
            v=dataset_p%phenoAnimal%y(1,kkd)
            ! Polygenic sire
            v = v - x(dataset_p%genea%np+ip)

            !General
            !QQ QQ
            if ( casSire == 1 ) then
             if ( casDam == 1 ) then
               v = v - x(2*dataset_p%genea%np+1)
             else if ( casDam == 2 .or. casDam == 3 ) then
               v = v - (x(2*dataset_p%genea%np+1))*0.5d0
             end if
            !QQ Qq
            else if (casSire == 2 ) then
             if ( casDam == 1 ) then
               v = v - (spt_p%pdd(wp_p%chr(1),kd,1,wp_p%n(1))+&
                 spt_p%pdd(wp_p%chr(1),kd,2,wp_p%n(1)))*x(2*dataset_p%genea%np+1)
             else if ( casDam == 2 .or. casDam == 3 ) then
               v = v - (spt_p%pdd(wp_p%chr(1),kd,1,wp_p%n(1))+spt_p%pdd(wp_p%chr(1),kd,2,wp_p%n(1)))*&
                       (x(2*dataset_p%genea%np+1))*0.5d0 &
                     + (spt_p%pdd(wp_p%chr(1),kd,3,wp_p%n(1))+spt_p%pdd(wp_p%chr(1),kd,4,wp_p%n(1)))*&
                        (x(2*dataset_p%genea%np+1))*0.5d0
             else if ( casDam == 4 ) then
              v = v &
             + (spt_p%pdd(wp_p%chr(1),kd,3,wp_p%n(1))+&
             spt_p%pdd(wp_p%chr(1),kd,4,wp_p%n(1)))*x(2*dataset_p%genea%np+1)
             end if
            !QQ qQ
            else if (casSire == 3 ) then
            if ( casDam == 1 ) then
               v = v &
             - (spt_p%pdd(wp_p%chr(1),kd,3,wp_p%n(1))+&
             spt_p%pdd(wp_p%chr(1),kd,4,wp_p%n(1)))*x(2*dataset_p%genea%np+1)
             else if ( casDam == 2 .or. casDam == 3 ) then
               v = v + (spt_p%pdd(wp_p%chr(1),kd,1,wp_p%n(1))+spt_p%pdd(wp_p%chr(1),kd,2,wp_p%n(1)))*&
                       (x(2*dataset_p%genea%np+1))*0.5d0 &
                     - (spt_p%pdd(wp_p%chr(1),kd,3,wp_p%n(1))+spt_p%pdd(wp_p%chr(1),kd,4,wp_p%n(1)))*&
                        (x(2*dataset_p%genea%np+1))*0.5d0
             else if ( casDam == 4 ) then
              v = v + (spt_p%pdd(wp_p%chr(1),kd,1,wp_p%n(1))+&
               spt_p%pdd(wp_p%chr(1),kd,2,wp_p%n(1)))*x(2*dataset_p%genea%np+1)

             end if
            !QQ qq
            else if (casSire == 4 ) then
              if ( casDam == 2 .or. casDam == 3 ) then
               v = v + (x(2*dataset_p%genea%np+1))*0.5d0
             else if ( casDam == 4 ) then
              v = v + x(2*dataset_p%genea%np+1)
             end if
            end if

            vpf=vpf+v*v*dataset_p%phenoAnimal%cd(1,kkd)
            end if
           end do !kd
           vmere=vmere+spt_p%probg(wp_p%chr(1),ig)*exp(-0.5d0*vpf/var)
          end do !ig
          fm=(p1*probcasP1(casDam)+p2*probcasP2(casDam)+p3*probcasP3(casDam))*vmere+fm
         end do  !casDam
         do i=1,effm
          fm=fm*invsig
         end do
         fmtot=fmtot*fm
       end do !jm
       fp = (p1*probcasP1(casSire)+p2*probcasP2(casSire)+p3*probcasP3(casSire))*fmtot+fp

       end do ! casSire
       if ( fp == 0 ) then
        f=INIFINY_REAL_VALUE
       else
        f = -dlog(fp)
       end if
      ! print *,p1,p2,f
      ! stop
      end subroutine funct_1qtl_fam_halfsib_add


       subroutine set_solution_hypothesis0(dataset,dataw,incsol)
       type(QTLMAP_DATASET)       ,intent(in) :: dataset
       type(TYPE_INCIDENCE_SOLUTION)      ,intent(inout)    :: incsol
       type (work_bicfarnir)           ,intent(in)        :: dataw

       integer :: nteff,maxNbPar,jm,ifem,ip,ieff
       type(GENEALOGY_BASE) , pointer :: dg
       type(PHENOTYPE_BASE) , pointer :: dpa
       type(DATAMODEL_BASE) , pointer :: dpm

       dpm => dataset%phenoModel
       dg => dataset%genea
       dpa => dataset%phenoAnimal

       allocate (incsol%sig(1,dg%np))

       incsol%hypothesis=0
       !  Mean family
       nteff = 1
       if ( dataw%ndamestim > 0 ) nteff = 2

       maxNbPar = max(dg%np,dataw%ndamestim)
       allocate (incsol%groupeName(nteff))
       allocate (incsol%eqtl_print(nteff))
       allocate (incsol%nbParameterGroup(nteff))
       allocate (incsol%parameterName(nteff,maxNbPar))
       allocate (incsol%paramaterValue(nteff,maxNbPar))
       allocate (incsol%parameterVecsol(nteff,maxNbPar))
       allocate (incsol%parameterPrecis(nteff,maxNbPar))
       incsol%parameterPrecis=0.d0
       incsol%parameterVecsol=.true.
       incsol%eqtl_print=.true.

       do ip=1,dg%np
            incsol%sig(1,ip) = dataw%parSolH0(ip)*dpm%sigt(dataw%ic)
       end do

       ieff=1
       incsol%groupeName(ieff) = 'Mean Sire'
       incsol%nbParameterGroup(ieff)=dg%np

       do ip=1,dg%np
           incsol%parameterName(ieff,ip)='Sire '//trim(dg%pere(ip))
           incsol%paramaterValue(ieff,ip) = dataw%parSolH0(dg%np+ip)*dpm%sigt(dataw%ic) + dpm%xmut(dataw%ic)
           incsol%parameterVecsol(ieff,ip) = .true.
       end do

       if ( dataw%ndamestim > 0 ) then
           ieff = ieff +1
           incsol%groupeName(ieff) = 'Mean dam'
           incsol%nbParameterGroup(ieff)=dataw%ndamestim
           ifem=0

          do ip=1,dg%np
           do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
             if ( dpa%estime(dataw%ic,jm)) then
              ifem=ifem+1
              incsol%parameterName(ieff,ifem)='Dam '//trim(dg%mere(jm))//" [Sire "//trim(dg%pere(ip))//"]"
              incsol%paramaterValue(ieff,ifem) = dataw%parSolH0(2*dg%np+ifem)*dpm%sigt(dataw%ic)
              incsol%parameterVecsol(ieff,ifem) = .true.
             end if
           end do
         end do

         end if

       end subroutine set_solution_hypothesis0


    subroutine set_solution_hypothesis1(dataset,dataw,incsol)
       type(QTLMAP_DATASET)       ,intent(in) :: dataset
       type (work_bicfarnir)           ,intent(in)        :: dataw
       type(TYPE_INCIDENCE_SOLUTION)      ,intent(inout)    :: incsol

       integer :: nteff,maxNbPar,jm,ifem,ip,ieff,cas
       real :: p,sumUnp1p2,p1,p2,p3
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
       if ( dataw%ndamestim > 0 ) nteff = 3

       maxNbPar = max(dg%np,dataw%ndamestim)
       allocate (incsol%groupeName(nteff))
       allocate (incsol%eqtl_print(nteff))
       allocate (incsol%nbParameterGroup(nteff))
       allocate (incsol%parameterName(nteff,maxNbPar))
       allocate (incsol%paramaterValue(nteff,maxNbPar))
       allocate (incsol%parameterVecsol(nteff,maxNbPar))
       allocate (incsol%parameterPrecis(nteff,maxNbPar))
       allocate (incsol%qtl_groupeName(1,1))

       incsol%parameterPrecis=0.d0
       incsol%parameterVecsol=.true.
       incsol%eqtl_print=.true.

       do ip=1,dg%np
            incsol%sig(1,ip) = dataw%parSolH1(ip)*dpm%sigt(dataw%ic)
       end do

         ieff=1
         incsol%qtl_groupeName(1,1)=ieff
         incsol%groupeName(ieff) = 'Qtl effect'
         incsol%nbParameterGroup(ieff)=4

         sumUnp1p2=real(1.d0 + dataw%parSolH1(dg%np+1) + 2*dataw%parSolH1(dg%np+2))

         p1 = dataw%parSolH1(dg%np+1) / sumUnp1p2
         p2 = dataw%parSolH1(dg%np+2) / sumUnp1p2
         p3 = 1.d0 / sumUnp1p2

         do cas=1,4
     !    p = (p1*probcasP1(cas)+p2*probcasP2(cas)+p3*probcasP3(cas))

      !   incsol%parameterName(ieff,cas)=casGeno(cas)//" [prob="//trim(str(int(p*100.d0)))//"%]"
!         incsol%paramaterValue(ieff,cas) = (dataw_p%alpha1(cas)*dataw%parSolH1(2*dg%np+3+dataw%ndamestim+1)+&
!                                           dataw_p%alpha2(cas)*dataw%parSolH1(2*dg%np+3+dataw%ndamestim+2)+&
!                                           dataw_p%alpha3(cas)*dataw%parSolH1(2*dg%np+3+dataw%ndamestim+3)+&
!                                           dataw_p%alpha4(cas)*dataw%parSolH1(2*dg%np+3+dataw%ndamestim+4))*dpm%sigt(dataw%ic)
       !  incsol%paramaterValue(ieff,cas) = dataw%parSolH1(2*dg%np+3+dataw%ndamestim+cas)*dpm%sigt(dataw%ic)
         end do

       incsol%parameterName(ieff,1)=casGeno(1)//" [prob="//trim(str(int(p1*100.d0)))//"%]"
       incsol%paramaterValue(ieff,1) = dataw%parSolH1(2*dg%np+3+dataw%ndamestim+1)*dpm%sigt(dataw%ic)
       incsol%parameterName(ieff,2)=casGeno(2)//" [prob="//trim(str(int(p2*100.d0)))//"%]"
       incsol%paramaterValue(ieff,2) = 0
       incsol%parameterName(ieff,3)=casGeno(3)//" [prob="//trim(str(int(p2*100.d0)))//"%]"
       incsol%paramaterValue(ieff,3) = 0
       incsol%parameterName(ieff,4)=casGeno(4)//" [prob="//trim(str(int(p3*100.d0)))//"%]"
       incsol%paramaterValue(ieff,4) = - dataw%parSolH1(2*dg%np+3+dataw%ndamestim+1)*dpm%sigt(dataw%ic)

       ieff = ieff +1
       incsol%groupeName(ieff) = 'Mean Sire'
       incsol%nbParameterGroup(ieff)=dg%np

       do ip=1,dg%np
           incsol%parameterName(ieff,ip)='Sire '//trim(dg%pere(ip))
           incsol%paramaterValue(ieff,ip) = dataw%parSolH1(dg%np+2+ip)*dpm%sigt(dataw%ic) + dpm%xmut(dataw%ic)
           incsol%parameterVecsol(ieff,ip) = .true.
       end do

       if ( dataw%ndamestim > 0 ) then
           ieff = ieff +1
           incsol%groupeName(ieff) = 'Mean dam'
           incsol%nbParameterGroup(ieff)=dataw%ndamestim
           ifem=0
         do ip=1,dg%np
           do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
             if ( dpa%estime(dataw%ic,jm) ) then
              ifem=ifem+1
              incsol%parameterName(ieff,ifem)='Dam '//trim(dg%mere(jm))//" [Sire "//trim(dg%pere(ip))//"]"
              incsol%paramaterValue(ieff,ifem) = dataw%parSolH1(2*dg%np+2+ifem)*dpm%sigt(dataw%ic)
              incsol%parameterVecsol(ieff,ifem) = .true.
             end if
           end do
         end do
         end if

       end subroutine set_solution_hypothesis1


end module m_qtlmap_analyse_biallelic_farnir
