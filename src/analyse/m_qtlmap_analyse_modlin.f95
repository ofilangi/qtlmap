!!****m* ANALYSE/m_qtlmap_analyse_modlin
!!  NAME
!!    m_qtlmap_analyse_modlin
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
module m_qtlmap_analyse_modlin
    use m_qtlmap_types
    use m_qtlmap_log
    use m_qtlmap_math
    use m_qtlmap_optimization
    use m_qtlmap_analyse_gen
    use m_qtlmap_analyse_lin_gen
    use m_qtlmap_output_handler

    implicit none
    save
!!****v* m_qtlmap_analyse_modlin/sig1
!!  NAME
!!   sig1
!!  DESCRIPTION
!!   The standart deviation under H0
!!  DIMENSIONS
!!   np
!!***
    real (kind=dp)       ,dimension(:),allocatable,public   :: sig1
!!****v* m_qtlmap_analyse_modlin/xmu1p
!!  NAME
!!   xmu1p
!!  DESCRIPTION
!!   The polygenic mean for each sire family under H0
!!  DIMENSIONS
!!   np
!!***
    real (kind=dp)       ,dimension(:),allocatable,public   :: xmu1p
!!****v* m_qtlmap_analyse_modlin/xmu1m
!!  NAME
!!   xmu1m
!!  DESCRIPTION
!!   The polygenic mean for each full sib family
!!  DIMENSIONS
!!   nm
!!***
    real (kind=dp)       ,dimension(:),allocatable,public   :: xmu1m
!!****v* m_qtlmap_analyse_modlin/xmu1g
!!  NAME
!!   xmu1g
!!  DESCRIPTION
!!   The general mean
!!***
    real (kind=dp)                                ,public   :: xmu1g
!!****v* m_qtlmap_analyse_modlin/f0
!!  NAME
!!   f0
!!  DESCRIPTION
!!   value of the likelihood under H0
!!***
    real (kind=dp)                                ,public   :: f0

    !$omp threadprivate (sig1,xmu1p,xmu1m,xmu1g,f0)


!!****v* m_qtlmap_analyse_modlin/fp0
!!  NAME
!!   fp0
!!  DESCRIPTION
!!   value of the likelihood by sire family under H0
!!  DIMENSIONS
!!   np
!!***
    real (kind=dp)       ,dimension(:),pointer,public   :: fp0
!!****v* m_qtlmap_analyse_modlin/fp1
!!  NAME
!!   fp1
!!  DESCRIPTION
!!   value of the likelihood by sire family under H1
!!  DIMENSIONS
!!   np
!!***
    real (kind=dp)       ,dimension(:),pointer,private  :: fp1
!!****v* m_qtlmap_analyse_modlin/fm0
!!  NAME
!!   fm0
!!  DESCRIPTION
!!   value of the likelihood by full-sib family under H0
!!  DIMENSIONS
!!   nm
!!***
    real (kind=dp)       ,dimension(:),pointer,public   :: fm0
!!****v* m_qtlmap_analyse_modlin/fm1
!!  NAME
!!   fm1
!!  DESCRIPTION
!!   value of the likelihood by full-sib family under H1
!!  DIMENSIONS
!!   nm
!!***
    real (kind=dp)       ,dimension(:),pointer,private  :: fm1

    !$omp threadprivate (fp0,fp1,fm0,fm1)


!!****v* m_qtlmap_analyse_modlin/current_ic
!!  NAME
!!   current_ic
!!  DESCRIPTION
!!   the current trait to analyse while the likelihood calculus
!!***
    integer              ,private                           :: current_ic
!!****v* m_qtlmap_analyse_modlin/current_chr
!!  NAME
!!   current_chr
!!  DESCRIPTION
!!   the current chromosome to analyse while the likelihood calculus
!!***
    integer              ,private                           :: current_chr

    !$omp threadprivate (current_ic,current_chr)

    type(QTLMAP_DATASET) , pointer , private :: dataset_p => null()
    type(PDD_BUILD)      , pointer , private :: spt_p => null()

    public :: init_analyse_modlin
    public :: opti_0qtl_modlin
    public :: opti_1qtl_modlin
    public :: test_lin
    public :: end_analyse_modlin
    public :: set_solution_hypothesis0
    public :: set_solution_hypothesis1

    contains
!!****f* m_qtlmap_analyse_modlin/init_analyse_modlin
!!  NAME
!!    init_analyse_modlin
!!  DESCRIPTION
!!    Initialisation/allocation of solution/buffer arrays
!!
!!  NOTES
!!  SOURCE
      subroutine init_analyse_modlin(dataset,spt)
         type(QTLMAP_DATASET) ,target      ,intent(in)         :: dataset
         type(PDD_BUILD)      ,target      ,intent(in)         :: spt

         integer           :: stat

         type(GENEALOGY_BASE) , pointer :: dg


         dataset_p => dataset
         spt_p => spt
         dg => dataset%genea

         allocate (sig1(dg%np),STAT=stat)
         call check_allocate(stat,'sig1 [m_qtlmap_analyse_modlin]')
         allocate (xmu1p(dg%np),STAT=stat)
         call check_allocate(stat,'xmu1p [m_qtlmap_analyse_modlin]')
         allocate (xmu1m(dg%nm),STAT=stat)
         call check_allocate(stat,'xmu1m [m_qtlmap_analyse_modlin]')
         allocate (fp0(dg%np),STAT=stat)
         call check_allocate(stat,'fp0 [m_qtlmap_analyse_modlin]')
         allocate (fp1(dg%np),STAT=stat)
         call check_allocate(stat,'fp1 [m_qtlmap_analyse_modlin]')
         allocate (fm0(dg%nm),STAT=stat)
         call check_allocate(stat,'fm0 [m_qtlmap_analyse_modlin]')
         allocate (fm1(dg%nm),STAT=stat)
         call check_allocate(stat,'fm1 [m_qtlmap_analyse_modlin]')

     end subroutine init_analyse_modlin
!!***

!!****f* m_qtlmap_analyse_modlin/end_analyse_modlin
!!  NAME
!!    end_analyse_modlin
!!  DESCRIPTION
!!    deallocation of solution/buffer arrays
!!
!!  NOTES
!!  SOURCE
     subroutine end_analyse_modlin
         deallocate (sig1)
         deallocate (xmu1p)
         deallocate (xmu1m)
         deallocate (fp0)
         deallocate (fp1)
         deallocate (fm0)
         deallocate (fm1)
     end subroutine end_analyse_modlin
!!***

!!****f* m_qtlmap_analyse_modlin/opti_0qtl_modlin
!!  NAME
!!    opti_0qtl_modlin
!!  DESCRIPTION
!!    Calcul de la vraisemblance 0 QTL, 1 caractere , effets parasites inclus
!!
!!  NOTES
!!  SOURCE
      subroutine opti_0qtl_modlin(dataset,spt,ic)
       type(QTLMAP_DATASET)       ,intent(in)         :: dataset
       type(PDD_BUILD)            ,intent(in)         :: spt

       integer , intent(in)   ::ic
!
! Divers

      integer                                    :: iuser(1)
      integer        ,dimension(:),allocatable   :: iw
      real (kind=dp) ,dimension(:),allocatable   :: par,borni,borns,w
      real (kind=dp)                             :: user(1)
      real (kind=dp)                             :: f
      integer :: npar,ibound,ip,ix,ifail,i,indest,ntot
      integer :: km,jm,kd,kd1,kd2,ii
      logical itest
      logical  , dimension(:,:,:),pointer        :: filter_inc
      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa
      type(DATAMODEL_BASE) , pointer :: dpm

      dpm => dataset%phenoModel
      dg => dataset%genea
      dpa => dataset%phenoAnimal

!
!******************************************************************************
!
!  recherche des elements estimables � partir de la d�composition de choleski
!  comptage et recodification de ces elements (on commence par mettre � 0 les
!  compteurs destin�s aux tests des effes fix�s
!
      itest=.false.
      call contingence(dataset,spt,ic,0,itest,.true.)
      call precision(xx,precis)
      !  Parametres de maximisation
      npar=dg%np+nbnivest

      allocate (borni(npar))
      allocate (borns(npar))
      allocate (par(npar))

      !MODIF - OPMIZATION
      allocate (filter_inc(dg%np,dg%nm, npar))
      call set_filter_optim(dataset,ic,.true.,.false.,ntnivmax,ntniv,vecsol,xinc,filter_inc)
      !FIN MODIF - OPMIZATION

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
      par(dg%np+1)=0.d0
      do ip=1,dg%np
        par(ip)=sig0(ip)
        par(dg%np+1)=par(dg%np+1)+xmu0p(ip)
      end do
      par(dg%np+1)=par(dg%np+1)/dble(dg%np)
      do ix=dg%np+2,npar
        par(ix)=0.d0
      end do
!
! Optimisation de la vraisemblance
      ifail=1
      current_ic = ic ! pour le mode lineaire....
      !call minimizing_funct(npar,ibound,funct_0qtl_modlin,borni,borns,par,f,iuser,user,ifail)
      call minimizing_funct_family(dataset,npar,ibound,funct_0qtl_modlin_family,&
       filter_inc,fm0,fp0,borni,borns,par,f,iuser,user,ifail)

      f0=f
      do i=1,npar
        par0(i)=par(i)
      end do
      do i=1,ntniv
        vecsol0(i)=vecsol(i)
        precis0(i)=precis(i)
      end do

      do ip = 1,dg%np
        sig1(ip)=par(ip)
      end do

      xmu1g=par(dg%np+1)
      indest=1
      do ip = 1,dg%np
      if (vecsol(1+ip)) then
        indest=indest+1
        xmu1p(ip)=par(dg%np+indest)
      end if
      end do

      ntot=dg%np+1
      km=0
      do jm=1,dg%nm
        if (dpa%estime(ic,jm).and.dataset%params%opt_sib.eq.2)then
          km=km+1
          if(vecsol(ntot+km)) then
            indest=indest+1
            xmu1m(jm)=par(dg%np+indest)
          end if
        end if
      end do

      deallocate (borni)
      deallocate (borns)
      deallocate (par)
      deallocate (filter_inc)

      return
      end subroutine opti_0qtl_modlin
!!***

!!****f* m_qtlmap_analyse_modlin/funct_0qtl_modlin
!!  NAME
!!    funct_0qtl_modlin
!!  DESCRIPTION
!!    Calcul de la vraisemblance et de ses derivees partielles sous H0
!!
!!  NOTES
!!  SOURCE
      subroutine funct_0qtl_modlin(n,x,f,iuser,user)
      use m_qtlmap_analyse_lin_gen, only : nbfem,nbef,ntnifix,nbco,covdir,nivdir,corniv,vecsol

      implicit none
      integer         , intent(in)                  :: n
!
! Tableaux dimensionnes selon n, le nombre de parametres a estimer
!
      real (kind=dp)      ,dimension(n), intent(in) :: x
      real (kind=dp)  , intent(inout)               :: f
      integer ,       dimension(1), intent(in)      :: iuser
      real (kind=dp)      ,dimension(1), intent(in) :: user

      ! Tableaux dimensionnes selon nm, le nombre de meres
      real (kind=dp)   ,dimension(dataset_p%genea%nm)  :: effm

      integer ::jnit,ip,nm1,nm2,jm,nd1,nd2,kkd,ilev,ief,ico,icar
      real (kind=dp) :: sig,var,vmere,vpf,v
!

!******************************************************************************
!
      jnit=3
      icar = current_ic
      if (nbfem.eq.0)jnit=2
      f=0.d0
      do ip=1,dataset_p%genea%np
        sig=x(ip)
        var=sig*sig
        fp0(ip)=0.d0
        nm1=dataset_p%genea%nmp(ip)+1
        nm2=dataset_p%genea%nmp(ip+1)
        do jm=nm1,nm2
          effm(jm)=0.d0
          vmere=0.d0
          vpf=0.d0
!
! on ne consid�re que les m�res
!
          nd1=dataset_p%genea%ndm(jm)+1
          nd2=dataset_p%genea%ndm(jm+1)
!
          do kkd=nd1,nd2
!
            if(dataset_p%phenoAnimal%presentc(icar,kkd)) then
              effm(jm)=effm(jm)+1.d0
              v=dataset_p%phenoAnimal%y(icar,kkd)
!
              if(vecsol(nivdir(kkd,1))) then
                ilev=corniv(nivdir(kkd,1))
                v=v-x(dataset_p%genea%np+ilev)
              end if

!
              if(vecsol(nivdir(kkd,2))) then
                ilev=corniv(nivdir(kkd,2))
                v=v-x(dataset_p%genea%np+ilev)
              end if
!
              if(dataset_p%phenoAnimal%estime(icar,jm)) then
                if(vecsol(nivdir(kkd,3))) then
                  ilev=corniv(nivdir(kkd,3))
                  v=v-x(dataset_p%genea%np+ilev)
                end if
              end if

              do ief=jnit+1,jnit+nbef
                if(vecsol(nivdir(kkd,ief))) then
                  ilev=corniv(nivdir(kkd,ief))
                  v=v-x(dataset_p%genea%np+ilev)
                end if
              end do

              do ico=1,nbco
                if(vecsol(ntnifix+ico)) then
                  ilev=corniv(ntnifix+ico)
                  v=v-covdir(kkd,ico)*x(dataset_p%genea%np+ilev)
                end if
              end do
              
              vpf=vpf+v*v*dataset_p%phenoAnimal%cd(icar,kkd)
            end if

          end do
          vmere=dexp(-0.5d0*vpf/var)
          if (vmere == 0) then
                  fm0(jm)=INIFINY_REAL_VALUE
          else
                  fm0(jm)=-dlog(vmere)+effm(jm)*dlog(sig)
          end if
          fp0(ip)=fp0(ip)+fm0(jm)
          f=f+fm0(jm)
        end do

      end do
      return
      end subroutine funct_0qtl_modlin
!!***

!!****f* m_qtlmap_analyse_modlin/funct_0qtl_modlin_family
!!  NAME
!!    funct_0qtl_modlin_family
!!  DESCRIPTION
!!    Calcul de la vraisemblance et de ses derivees partielles sous H0
!!
!!  NOTES
!!  SOURCE
    subroutine funct_0qtl_modlin_family(ip,jm,n,x,f,iuser,user)
      use m_qtlmap_analyse_lin_gen, only : nbfem,nbef,ntnifix,nbco,covdir,nivdir,corniv,vecsol

      implicit none
      integer         , intent(in)                  :: ip,jm,n
!
! Tableaux dimensionnes selon n, le nombre de parametres a estimer
!
      real (kind=dp)      ,dimension(n), intent(in) :: x
      real (kind=dp)  , intent(inout)               :: f
      integer ,       dimension(1), intent(in)      :: iuser
      real (kind=dp)      ,dimension(1), intent(in) :: user

      ! Tableaux dimensionnes selon nm, le nombre de meres
      real (kind=dp)    :: effm

      integer ::jnit,nm1,nm2,nd1,nd2,kkd,ilev,ief,ico,icar
      real (kind=dp) :: sig,var,vmere,vpf,v
!

!******************************************************************************
!
      jnit=3
      icar = current_ic
      if (nbfem.eq.0)jnit=2
      f=0.d0
      sig=x(ip)
      var=sig*sig
      effm=0.d0
      vmere=0.d0
      vpf=0.d0
!
! on ne consid�re que les m�res
!
      nd1=dataset_p%genea%ndm(jm)+1
      nd2=dataset_p%genea%ndm(jm+1)
!
          do kkd=nd1,nd2
!
            if(dataset_p%phenoAnimal%presentc(icar,kkd)) then
              effm=effm+1.d0
              v=dataset_p%phenoAnimal%y(icar,kkd)
!
              if(vecsol(nivdir(kkd,1))) then
                ilev=corniv(nivdir(kkd,1))
                v=v-x(dataset_p%genea%np+ilev)
              end if

!
              if(vecsol(nivdir(kkd,2))) then
                ilev=corniv(nivdir(kkd,2))
                v=v-x(dataset_p%genea%np+ilev)
              end if
!
              if(dataset_p%phenoAnimal%estime(icar,jm)) then
                if(vecsol(nivdir(kkd,3))) then
                  ilev=corniv(nivdir(kkd,3))
                  v=v-x(dataset_p%genea%np+ilev)
                end if
              end if

              do ief=jnit+1,jnit+nbef
                if(vecsol(nivdir(kkd,ief))) then
                  ilev=corniv(nivdir(kkd,ief))
                  v=v-x(dataset_p%genea%np+ilev)
                end if
              end do

              do ico=1,nbco
                if(vecsol(ntnifix+ico)) then
                  ilev=corniv(ntnifix+ico)
                  v=v-covdir(kkd,ico)*x(dataset_p%genea%np+ilev)
                end if
              end do

              vpf=vpf+v*v*dataset_p%phenoAnimal%cd(icar,kkd)
            end if

          end do
          vmere=dexp(-0.5d0*vpf/var)
          if (vmere == 0) then
                  f=INIFINY_REAL_VALUE
          else
                  f=-dlog(vmere)+effm*dlog(sig)
          end if


      end subroutine funct_0qtl_modlin_family
!!***

!!****f* m_qtlmap_analyse_modlin/opti_1qtl_modlin
!!  NAME
!!    opti_1qtl_modlin
!!  DESCRIPTION
!!    Calcul de la statistique de test le long du chromosome
!!
!!  NOTES
!!  SOURCE
      subroutine opti_1qtl_modlin(dataset,spt,ic,lrtsol,fmax,supnbnivest)

      integer , intent(in)                                :: ic
      type(TYPE_LRT_SOLUTION)  , intent(out)              :: lrtsol
      real (kind=dp)  , intent(out)                       :: fmax
      integer         , intent(out)                       :: supnbnivest
      type(QTLMAP_DATASET)         ,intent(in)            :: dataset
      type(PDD_BUILD)              ,intent(in)            :: spt
!

! Divers
      real (kind=dp) ,dimension(:),allocatable :: val,par,borni,borns
      real (kind=dp) :: user(1)
      real (kind=dp) :: fm(dataset_p%genea%nm),fp(dataset_p%genea%np),f1,save_f0

      logical itest,stvecsol(size(vecsol1))
      integer :: ip,i,n,ilong,ibound,npar,ix,j,ifail,iuser(1),indexchr(dataset%map%nchr),ntotal,nprime
      integer :: ii,chr,nkd,jm,ig,ifem,iip,kd,kd1,kd2,save_nteffmax,save_ntnivmax
      logical  , dimension(:,:,:),pointer        :: filter_inc
      real(kind=dp)       ,dimension(:,:) ,pointer           :: listef
      integer             ,dimension(:,:)    ,pointer        :: listenbnivest
      real(kind=dp)       ,dimension(:,:,:)   ,pointer       :: listepar
      real(kind=dp)       ,dimension(:) ,pointer           :: save_fm0
      real(kind=dp)       ,dimension(:) ,pointer           :: save_fp0
      real(kind=dp)       ,dimension(dataset_p%genea%np)   :: save_sig1
      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa
      type(DATAMODEL_BASE) , pointer :: dpm

      real(kind=dp)       ,dimension(:,:,:)   ,pointer     :: xlrp,xlrm
      real(kind=dp)       ,dimension(:,:)   ,pointer       :: lrt1

      dpm => dataset%phenoModel
      dg => dataset%genea
      dpa => dataset%phenoAnimal


!******************************************************************************
! Calcul de la vraisemblance sous H1

      save_fp0=>fp0
      save_fm0=>fm0
      save_ntnivmax = ntnivmax
      save_nteffmax = nteffmax
      save_f0 = f0
      save_sig1 = sig1

      call lrtsol%new(dataset,1)

! initialisation
!  on utilisera  :
!  PAR (vecteur des param�tres � optimiser),
!  CORNIV (vecteur des positions, parmi NTNIV, des effets extimables),
!  SOLVEC (vecteur disant l'estimabilit� des effets)
!  VAL le vecteur complet (ntniv positions) des valeurs des niveaux des effets
!  STVAL,STVECSOL et STCORNIV les copies de VAL, VECSOL et CORNIV � l'it�ration
!  n-1 quand on analyse la position n
!
!
      lrtsol%lrtmax=-1.d75
      lrtsol%chrmax=0
      lrtsol%nxmax=0

      ntotal=0
      do chr=1,dataset%map%nchr
        ntotal=ntotal+dataset%map%get_npo(chr)
        indexchr(chr)=ntotal
      end do

      allocate (lrt1(dataset%map%nchr,dataset%map%get_maxnpo()))
      allocate (xlrp(dataset%map%nchr,dataset%map%get_maxnpo(),dg%np))
      allocate (xlrm(dataset%map%nchr,dataset%map%get_maxnpo(),dg%nm))

      allocate (listef(dataset%map%nchr,dataset%map%get_maxnpo()))
      allocate (listenbnivest(dataset%map%nchr,dataset%map%get_maxnpo()))
      allocate (listepar(dataset%map%nchr,dataset%map%get_maxnpo(),dg%np + ntnivmax))
      !on libere l espace pour le cas non openmp
      call end_contingence
!
! Marche le long du chromosome
     !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(par,val)  &
     !$OMP PRIVATE(filter_inc,i,ifem,iip,ii,stvecsol,npar,j,ifail,fp,fm,f1,n,chr,borni,borns)
     current_ic = ic
     ntnivmax = save_ntnivmax
     nteffmax = save_nteffmax
     fp0=>save_fp0
     fm0=>save_fm0
     f0 = save_f0

     allocate ( val( dg%np+ntnivmax ) )
     allocate ( par( dg%np + ntnivmax) )
     allocate ( borni( dg%np + ntnivmax) )
     allocate ( borns( dg%np + ntnivmax) )

! Point de depart
!
      do ip=1,dg%np
        par(ip)=save_sig1(ip)
      end do

      par(dg%np+1)=xmu1g

      do i=dg%np+2,size(par)
        par(i)=0.d0
      end do

      do i=dg%np+1,size(val)
        val(i-dg%np)=par(i)
      end do

     call init_contingence(dataset,spt)
     vecsol=.true.
     allocate (filter_inc(dg%np,dg%nm,ntnivmax+dg%np))
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

!
!  on stocke les conditions en n
        do i=1,ntniv
          stvecsol(i)=vecsol(i)
        end do
!
!  preparation de la matrice d'incidence
!
        call prepinc(dataset,spt,current_chr,n,ic,"LA  ")
!
!  recherche des elements estimables � partir de la d�composition de choleski
!  comptage et recodification de ces elements (on commence par mettre � 0 les
!  compteurs destin�s aux tests des effes fix�s
!
        itest=.false.
        call contingence(dataset,spt,ic,1,itest)
        !MODIF - OPMIZATION
        call set_filter_optim(dataset,ic,.true.,.false.,ntnivmax,ntniv,vecsol,xinc,filter_inc)
        !FIN MODIF - OPMIZATION

! Parametres de maximisation
      ibound=0
      npar=dg%np+nbnivest

      do ip=1,dg%np
        borni(ip)=SIG_MIN
        borns(ip)=SIG_MAX
      end do
       do i=dg%np+1,npar
        borni(i)=XMU_MIN
        borns(i)=XMU_MAX
      end do
!
! Point de depart (on reprend les point d'arriv�e pr�c�dents
!
      do i=dg%np+1,npar
        par(i)=0.d0
      end do

      j=dg%np
      do i=1,ntniv
        if(vecsol(i))then
          j=j+1
          par(j)=0.d0
          if(stvecsol(i)) par(j)=val(i)
        end if
        stvecsol(i)=vecsol(i)
      end do

! Optimisation de la vraisemblance a la position dx
        ifail=1

        call minimizing_funct_family(dataset,npar,ibound,funct_1qtl_modlin_family,&
        filter_inc,fm,fp,borni,borns,par,f1,iuser,user,ifail)

        j=dg%np
        do i=1,ntniv
          if(vecsol(i)) then
           j=j+1
           val(i)=par(j)
          else
           val(i)=9999.d0
          end if
        end do
!
! on garde les valeurs du LRT pour pouvopir dessiner la courbe de vraisemblance
!
      if ( f1 < INIFINY_REAL_VALUE ) then
         lrt1(chr,n)=-2.d0*(f1-f0)
      else
         lrt1(chr,n)=0
      end if

      listef(chr,n)=f1
      listenbnivest(chr,n)=nbnivest
      listepar(chr,n,:npar)=par(:npar)
!
!  on met les profil / prog�niteur dnas les ficheir ad hoc
!
        iip=dg%np+1
        do ii=1,dg%np
          xlrp(chr,n,ii)=-2.d0*(fp(ii)-fp0(ii))
          if ( vecsol(1+ii)) then
           iip=iip+1
           ! ATTENTION VALIDE DANS LE CAS SANS INTERACTION AVEC QTL
           lrtsol%pater_eff(chr,ii,n)=par(iip)
          end if
        end do

        ifem=0
        do ii=1,dg%nm
          xlrm(chr,n,ii)=-2.d0*(fm(ii)-fm0(ii))
          if ( dpa%estime(ic,ii)) then
            ifem=ifem+1
            if ( vecsol(dataset_p%genea%np+1+ifem) ) then
              iip=iip+1
              ! ATTENTION VALIDE DANS LE CAS SANS INTERACTION AVEC QTL
              lrtsol%mater_eff(chr,ifem,n)=par(iip)
            end if
          end if
        end do
       end do  ! nprime
       !$OMP END DO
       call end_contingence
       deallocate (filter_inc)
       deallocate ( val )
       deallocate ( par )
       deallocate ( borni )
       deallocate ( borns )
      !$OMP END PARALLEL

      !recherche du maximum
      do chr=1,dataset%map%nchr
       do n=1,dataset%map%get_npo(chr)
           if(lrtsol%lrtmax(0) < lrt1(chr,n)) then
            lrtsol%nxmax(0)=n
            lrtsol%chrmax(0)=chr
            lrtsol%lrtmax(0)=lrt1(chr,n)
            fmax=listef(chr,n)
            supnbnivest=listenbnivest(chr,n)
            par1(:supnbnivest+dataset_p%genea%np)=listepar(chr,n,:supnbnivest+dataset_p%genea%np)
           end if
       end do
      end do

      deallocate (listef)
      deallocate (listenbnivest)
      deallocate (listepar)

      ! 04/2013 - New structure for LRT
      call lrtsol%LRT%add1(dataset,dataset%map%nchr,dataset%map%get_maxnpo(),lrt1)
      do ii=1,dg%np
       call lrtsol%LRT_SIRES(ii)%add1(dataset,dataset%map%nchr,dataset%map%get_maxnpo(),xlrp(:,:,ii))
      end do
      do ii=1,dg%nm
       call lrtsol%LRT_DAMS(ii)%add1(dataset,dataset%map%nchr,dataset%map%get_maxnpo(),xlrm(:,:,ii))
      end do

      deallocate (lrt1)
      deallocate (xlrp)
      deallocate (xlrm)

      call init_contingence(dataset,spt)

!  on calcule la pr�cision des estimation au point correspondant
! au LRT maximum
!
      call prepinc(dataset,spt,lrtsol%chrmax(0),lrtsol%nxmax(0),ic,"LA  ")
      call contingence(dataset,spt,ic,1,itest)
      call precision(xx,precis)

      do i=1,ntniv
        precis1(i)=precis(i)
        vecsol1(i)=vecsol(i)
      end do

      return
      end subroutine opti_1qtl_modlin
!!***

!!****f* m_qtlmap_analyse_modlin/funct_1qtl_modlin
!!  NAME
!!    funct_1qtl_modlin
!!  DESCRIPTION
!!    Calcul de la vraisemblance et de ses derivees partielles sous Hypothese 1QTL 1carac
!!
!!  NOTES
!!  SOURCE
      subroutine funct_1qtl_modlin(n,x,f,iuser,user)

      integer         , intent(in)                  :: n
!
! Tableaux dimensionnes selon n, le nombre de parametres a estimer
!
      real (kind=dp)      ,dimension(n), intent(in) :: x
      real (kind=dp)  , intent(inout)               :: f
      integer ,       dimension(1), intent(in)      :: iuser
      real (kind=dp)      ,dimension(1), intent(in) :: user

! Tableaux dimensionnes selon nm, le nombre de meres
      real (kind=dp)   ,dimension(dataset_p%genea%nm)   :: effm

!
! Divers
      integer :: init,jnit,ip,nm1,nm2,jm,ngeno1,ngeno2,ig
      integer :: nd1,nd2,kd,kkd,ilev,ief,ico,ic
      real (kind=dp) :: sig,var,vmere,vpf,v
!******************************************************************************
!

      ic = current_ic
      init=4
      jnit=5
      if (nbfem.eq.0)then
        init=3
        jnit=3
      end if

      f=0.d0
      do ip=1,dataset_p%genea%np
        sig=x(ip)
        var=sig*sig
        fp1(ip)=0.d0
        nm1=dataset_p%genea%nmp(ip)+1
        nm2=dataset_p%genea%nmp(ip+1)
        do jm=nm1,nm2
          vmere=0.d0
          ngeno1=spt_p%ngenom(current_chr,jm)+1
          ngeno2=spt_p%ngenom(current_chr,jm+1)
          do ig=ngeno1,ngeno2
            nd1=spt_p%ngend(current_chr,ig)+1
            nd2=spt_p%ngend(current_chr,ig+1)
            vpf=0.d0
            effm(jm)=0.d0
            do kd=nd1,nd2
              kkd=spt_p%ndesc(current_chr,kd)

              if(dataset_p%phenoAnimal%presentc(ic,kkd)) then
                effm(jm)=effm(jm)+1
                v=dataset_p%phenoAnimal%y(ic,kkd)

                ilev=corniv(nivdir(kkd,1))
                v=v-x(dataset_p%genea%np+ilev)
                if(vecsol(nivdir(kkd,2))) then
                  ilev=corniv(nivdir(kkd,2))
                  v=v-x(dataset_p%genea%np+ilev)*ppt(kd)
                end if

                if(dataset_p%phenoAnimal%estime(ic,jm)) then
                  if(vecsol(nivdir(kkd,3))) then
                    ilev=corniv(nivdir(kkd,3))
                    v=v-x(dataset_p%genea%np+ilev)*pmt(kd)
                  end if
                end if

                if(vecsol(nivdir(kkd,init))) then
                  ilev=corniv(nivdir(kkd,init))
                  v=v-x(dataset_p%genea%np+ilev)
                end if

                if(dataset_p%phenoAnimal%estime(ic,jm)) then
                  if(vecsol(nivdir(kkd,5))) then
                    ilev=corniv(nivdir(kkd,5))
                    v=v-x(dataset_p%genea%np+ilev)
                  end if
                end if

                do ief=jnit+1,jnit+nbef

                  if(vecsol(nivdir(kkd,ief))) then
                    ilev=corniv(nivdir(kkd,ief))
                    v=v-x(dataset_p%genea%np+ilev)
                  end if
                end do

                do ico=1,nbco
                  if(vecsol(ntnifix+ico)) then
                    ilev=corniv(ntnifix+ico)
                    v=v-covdir(kkd,ico)*x(dataset_p%genea%np+ilev)
                  end if
                end do
                vpf=vpf+v*v*dataset_p%phenoAnimal%cd(ic,kkd)
              end if
            end do
            vmere=vmere+spt_p%probg(current_chr,ig)*dexp(-0.5d0*vpf/var)
          end do

          if (vmere == 0) then
                 fm1(jm)=INIFINY_REAL_VALUE
                 f=INIFINY_REAL_VALUE
                 return
          else
                 fm1(jm)=-dlog(vmere)+dble(effm(jm))*dlog(sig)
          end if
          fp1(ip)=fp1(ip)+fm1(jm)
          f=f+fm1(jm)

        end do
      end do

      return
      end subroutine funct_1qtl_modlin
!!***

!!****f* m_qtlmap_analyse_modlin/funct_1qtl_modlin_family
!!  NAME
!!    funct_1qtl_modlin_family
!!  DESCRIPTION
!!    Calcul de la vraisemblance et de ses derivees partielles sous Hypothese 1QTL 1carac
!!
!!  NOTES
!!  SOURCE
    subroutine funct_1qtl_modlin_family(ip,jm,n,x,f,iuser,user)

      integer         , intent(in)                  :: ip,jm,n
!
! Tableaux dimensionnes selon n, le nombre de parametres a estimer
!
      real (kind=dp)      ,dimension(n), intent(in) :: x
      real (kind=dp)  , intent(inout)               :: f
      integer ,       dimension(1), intent(in)      :: iuser
      real (kind=dp)      ,dimension(1), intent(in) :: user

! Tableaux dimensionnes selon nm, le nombre de meres
      real (kind=dp)    :: effm

!
! Divers
      integer :: init,jnit,nm1,nm2,ngeno1,ngeno2,ig
      integer :: nd1,nd2,kd,kkd,ilev,ief,ico,ic
      real (kind=dp) :: sig,var,vmere,vpf,v
!******************************************************************************
!

      ic = current_ic
      init=4
      jnit=5
      if (nbfem.eq.0)then
        init=3
        jnit=3
      end if

      f=0.d0
      sig=x(ip)
      var=sig*sig
      vmere=0.d0
      ngeno1=spt_p%ngenom(current_chr,jm)+1
      ngeno2=spt_p%ngenom(current_chr,jm+1)
      do ig=ngeno1,ngeno2
         nd1=spt_p%ngend(current_chr,ig)+1
         nd2=spt_p%ngend(current_chr,ig+1)
         vpf=0.d0
         effm=0.d0
         do kd=nd1,nd2
            kkd=spt_p%ndesc(current_chr,kd)

            if(dataset_p%phenoAnimal%presentc(ic,kkd)) then
              effm=effm+1
                v=dataset_p%phenoAnimal%y(ic,kkd)
                ilev=corniv(nivdir(kkd,1))
                v=v-x(dataset_p%genea%np+ilev)
                if(vecsol(nivdir(kkd,2))) then
                  ilev=corniv(nivdir(kkd,2))
                  v=v-x(dataset_p%genea%np+ilev)*ppt(kd)
                end if

                if(dataset_p%phenoAnimal%estime(ic,jm)) then
                  if(vecsol(nivdir(kkd,3))) then
                    ilev=corniv(nivdir(kkd,3))
                    v=v-x(dataset_p%genea%np+ilev)*pmt(kd)
                  end if
                end if

                if(vecsol(nivdir(kkd,init))) then
                  ilev=corniv(nivdir(kkd,init))
                  v=v-x(dataset_p%genea%np+ilev)
                end if

                if(dataset_p%phenoAnimal%estime(ic,jm)) then
                  if(vecsol(nivdir(kkd,5))) then
                    ilev=corniv(nivdir(kkd,5))
                    v=v-x(dataset_p%genea%np+ilev)
                  end if
                end if

                do ief=jnit+1,jnit+nbef

                  if(vecsol(nivdir(kkd,ief))) then
                    ilev=corniv(nivdir(kkd,ief))
                    v=v-x(dataset_p%genea%np+ilev)
                  end if
                end do

                do ico=1,nbco
                  if(vecsol(ntnifix+ico)) then
                    ilev=corniv(ntnifix+ico)
                    v=v-covdir(kkd,ico)*x(dataset_p%genea%np+ilev)
                  end if
                end do
                vpf=vpf+v*v*dataset_p%phenoAnimal%cd(ic,kkd)
              end if
          end do
        vmere=vmere+spt_p%probg(current_chr,ig)*dexp(-0.5d0*vpf/var)
      end do

      if (vmere == 0) then
           f=INIFINY_REAL_VALUE
      else
           f=-dlog(vmere)+dble(effm)*dlog(sig)
      end if

      end subroutine funct_1qtl_modlin_family
!!***

!!****f* m_qtlmap_analyse_modlin/test_lin
!!  NAME
!!    test_lin
!!  DESCRIPTION
!!    Test des differents effets de nuisance du modele par une LRT compare a une chi2
!!
!!  NOTES
!!  SOURCE
      subroutine test_lin(dataset,spt,chr,ic,est_moy,supnbnivest,fmax,nposx)
       type(QTLMAP_DATASET)       ,intent(in)         :: dataset
       type(PDD_BUILD)            ,intent(in)         :: spt

      logical , intent(in)                        :: est_moy
      integer                        , intent(in) :: chr,ic
      integer                        , intent(in) :: supnbnivest
      real (kind=dp)                 , intent(in) :: fmax
      integer                        , intent(in) :: nposx

!
! Divers
      logical itest
      integer iuser(1)
      real (kind=dp)  ,dimension(:) , allocatable :: par,borni,borns
      real (kind=dp) :: user(1),prob,xlrt_t,f1
      integer :: nbint,iecd,iecq,ief,ibound,npar,ip,i,ifail
      integer :: nbreduit
      logical  , dimension(:,:,:),pointer        :: filter_inc

      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa
      type(DATAMODEL_BASE) , pointer :: dpm

      dpm => dataset%phenoModel
      dg => dataset%genea
      dpa => dataset%phenoAnimal

      write(nficout,600)
  600 format(//,80('*')/'testing model effects',//)

!
!  preparation de la matrice d'incidence
!
      call prepinc(dataset,spt,chr,nposx,ic,"LA  ")
!
!  on recalcule la vraisemblance en retirant les effets du mod�les un � un
!
!
!  on commence par les effets hierarchis�s dans les effets qtl
!
      current_chr = chr
      nbef=dpm%modele(ic,1)
      nbco=dpm%modele(ic,2)
      nbint=dpm%modele(ic,3)

      iecd=0
      iecq=0

      allocate ( par( dg%np + ntnivmax  ) )
      allocate ( borni( dg%np + ntnivmax ) )
      allocate ( borns( dg%np + ntnivmax ) )
      !MODIF - OPMIZATION
      allocate (filter_inc(dg%np,dg%nm, dg%np + ntnivmax ))
      !FIN MODIF - OPMIZATION


      do ief=1,nbint+nbef+nbco
      meff=0
      mcov=0
      mint=0
      if(ief.le.nbef)meff=ief
      if(ief.gt.nbef.and.ief.le.(nbef+nbco))mcov=ief-nbef
      if(ief.gt.(nbef+nbco))mint=ief-(nbef+nbco)
!
!  recherche des elements estimables � partir de la d�composition de choleski
!  comptage et recodification de ces elements
!

        itest=.true.
        call contingence(dataset,spt,ic,1,itest,est_moy)
        call set_filter_optim(dataset,ic,.true.,.false.,ntnivmax,ntniv,vecsol,xinc,filter_inc)

!
!  la r�dustion du nombre d'effet estim�e est calcul�e
!
      nbreduit=supnbnivest-nbnivest
!
! Parametres de maximisation
      ibound=0
      npar=dg%np+nbnivest

      do ip=1,dg%np
        borni(ip)=1.d-6
        borns(ip)=1.d6
      end do
       do i=dg%np+1,npar
        borni(i)=-1d6
        borns(i)=1.d6
        par(i)=0.d0
      end do
!
! Point de depart (on reprend les point d'arriv�e pr�c�dents
!
      do i=1,dg%np
      par(i)=par1(i)
      end do

      do i=dg%np+1,npar
        par(i)=0.d0
      end do
!
! Optimisation de la vraisemblance a la position dx
        ifail=1

     !   call minimizing_funct(npar,ibound,funct_1qtl_modlin,borni,borns,par,f1,iuser,user,ifail)
        call minimizing_funct_family(dataset,npar,ibound,funct_1qtl_modlin_family,&
         filter_inc,fm1,fp1,borni,borns,par,f1,iuser,user,ifail)
!        if (ifail.ne.0) print *,'Code retour e04jyf H1 : ',ifail
!
        xlrt_t=-2.d0*(fmax-f1)
      if(xlrt_t.le.1.d-8)xlrt_t=0.d0
!
!  impression des tests des effets du mod�le
!

      nbef=dpm%modele(ic,1)
      nbco=dpm%modele(ic,2)
      nbint=dpm%modele(ic,3)

      if(ief.le.(nbef+nbco).and.iecd.eq.0) then
        iecd=1
        write(nficout,610)
  610 format('  Direct effects'/)

        write(nficout,601)
  601 format('Tested effect     df.    Likelihood     p-value'/ &
       '                         ratio                 '/)
      end if

      if(ief.le.nbef) then
        ifail=0
  !      prob=MATH_QTLMAP_G01ECF('U',xlrt_t,dble(dpm%nlev(dpm%modele(ic,3+ief))),ifail)
         if ( nbreduit == 0 ) then
            call log_mess("The effect ["//trim(dpm%namefix(dpm%modele(ic,3+ief)))//"]"//&
                          " might be confused with another effect !",WARNING_DEF)
            prob=0.d0
         else
         prob=MATH_QTLMAP_G01ECF('U',xlrt_t,dble(nbreduit),ifail)
         end if

          write(nficout,602) trim(dpm%namefix(dpm%modele(ic,3+ief))),nbreduit,xlrt_t,prob

      end if
  602 format(a15,3x,i2,6x,f8.2,7x,f5.3)

      if(ief.gt.nbef.and.ief.le.(nbef+nbco)) then
        ifail=0
        prob=MATH_QTLMAP_G01ECF('U',xlrt_t,dble(1.0),ifail)

        write(nficout,603) trim(dpm%namecov(dpm%modele(ic,3+ief))),xlrt_t,prob
 
      end if
  603 format(a15,4x,'1',6x,f8.2,7x,f5.3)

      if(ief.gt.(nbef+nbco).and.iecq.eq.0)then
        iecq=1
        write(nficout,611)
  611 format(/'  Intra qtl effects'/)
        write(nficout,601)
      end if

      if(ief.gt.(nbef+nbco)) then
        ifail=0
        prob=MATH_QTLMAP_G01ECF('U',xlrt_t,dble(nbreduit),ifail)
        write(nficout,602) trim(dpm%namefix(dpm%modele(ic,3+ief))),nbreduit,xlrt_t,prob
  !      write(*,602) trim(dpm%namefix(dpm%modele(ic,3+ief))),nbreduit,xlrt_t,prob
      end if

      end do

      deallocate ( par  )
      deallocate ( borni )
      deallocate ( borns )
      deallocate (filter_inc)

      write(nficout,*) "When this probability exceeds the standard threshold corresponding to the 5, 1 or 0.1 Pent level",&
               ", you might consider removing this effect from the model"

      return
      end subroutine test_lin
!!***

!!****f* m_qtlmap_analyse_modlin/set_solution_hypothesis0
!!  NAME
!!    set_solution_hypothesis0
!!  DESCRIPTION
!!
!!  NOTES
!!  SOURCE
        subroutine set_solution_hypothesis0(dataset,ic,incsol)
          type(QTLMAP_DATASET)       ,intent(in)            :: dataset
          integer                            ,intent(in)       :: ic
          type(TYPE_INCIDENCE_SOLUTION)      ,intent(inout)    :: incsol

       integer :: nteff,maxNbPar,jm,ifem,ip,ieff
       integer :: indest,isol,ipar,nlevel,ief,i,ife

       real(kind=dp) ,dimension(size(par0)) :: par0_t
       type(GENEALOGY_BASE) , pointer :: dg
       type(PHENOTYPE_BASE) , pointer :: dpa
       type(DATAMODEL_BASE) , pointer :: dpm

       dpm => dataset%phenoModel
       dg => dataset%genea
       dpa => dataset%phenoAnimal

       ! remise a l echelle des variables
       do ipar=1,size(par0)
          par0_t(ipar)=par0(ipar)*dpm%sigt(ic)
       end do

       allocate (incsol%sig(1,dg%np))

       incsol%hypothesis=0
       !  Mean, Polygenic family
       nteff = 2
       if ( count(dpa%estime(ic,:))  > 0 ) nteff = 3

       !Fixed effect and covariate
       if ( dpm%modele(ic,1) > 0 ) nteff = nteff+1
       if ( dpm%modele(ic,2) > 0 ) nteff = nteff+1

       maxNbPar = max(dg%np,count(dpa%estime(ic,:)))

       !max numbre de covariable ?
       maxNbPar = max(maxNbPar,dpm%modele(ic,2))

       nlevel=0
       !max nombre de niveau pour un effet fixe ?
       do i=1,dpm%modele(ic,1)
         nlevel=dpm%nlev(dpm%modele(ic,3+i))+nlevel
       end do
       maxNbPar = max(maxNbPar,nlevel)

       allocate (incsol%groupeName(nteff))
       allocate (incsol%eqtl_print(nteff))
       allocate (incsol%nbParameterGroup(nteff))
       allocate (incsol%parameterName(nteff,maxNbPar))
       allocate (incsol%paramaterValue(nteff,maxNbPar))
       allocate (incsol%parameterVecsol(nteff,maxNbPar))
       allocate (incsol%parameterPrecis(nteff,maxNbPar))
       incsol%eqtl_print=.true.
       incsol%parameterPrecis=0.d0
       incsol%parameterVecsol=.true.

       do ip=1,dg%np
            incsol%sig(1,ip) = par0_t(ip)
       end do

       ieff=1
       incsol%groupeName(ieff) = 'General Mean'
       incsol%nbParameterGroup(ieff)=1
       incsol%parameterName(ieff,1)   ='General Mean'
       incsol%paramaterValue(ieff,1)  = par0_t(dg%np+1)+dpm%xmut(ic)
       incsol%parameterVecsol(ieff,1) = vecsol0(1)
       incsol%parameterPrecis(ieff,1) = precis0(1)


       ieff=ieff+1
       incsol%groupeName(ieff) = 'Sire polygenic effects'
       incsol%nbParameterGroup(ieff)=dg%np

       indest = dg%np+1
       isol=1
       do ip=1,dg%np
           isol=isol+1
           if ( vecsol0(isol) ) then
             indest=indest+1
             incsol%paramaterValue(ieff,ip)  = par0_t(indest)
           else
             incsol%paramaterValue(ieff,ip)  = 0.d0
           end if

           incsol%parameterName(ieff,ip)   ='Sire '//trim(dg%pere(ip))
           incsol%parameterVecsol(ieff,ip) = vecsol0(isol)
           incsol%parameterPrecis(ieff,ip)  = precis0(isol)
       end do

       if ( count(dpa%estime(ic,:)) > 0 ) then
           ieff = ieff +1
           incsol%groupeName(ieff) = 'Dam polygenic effects'
           incsol%nbParameterGroup(ieff)=count(dpa%estime(ic,:))
           ifem=0
          do ip=1,dg%np
           do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
             if ( dpa%estime(ic,jm)) then
              isol=isol+1
              ifem=ifem+1
              if ( vecsol0(isol) ) then
                indest=indest+1
                incsol%paramaterValue(ieff,ifem) = par0_t(indest)
              else
                incsol%paramaterValue(ieff,ifem) = 0.d0
              end if

              incsol%parameterName(ieff,ifem)='Dam '//trim(dg%mere(jm))//" [Sire "//trim(dg%pere(ip))//"]"
              incsol%parameterVecsol(ieff,ifem) = vecsol0(isol)
              incsol%parameterPrecis(ieff,ifem)  = precis0(isol)
             end if
           end do
          end do
         end if

       !Fixed effect
       if ( dpm%modele(ic,1) > 0 ) then
         ieff=ieff+1
         incsol%eqtl_print(ieff)=.false.
         incsol%groupeName(ieff) = 'Fixed effects'
         incsol%nbParameterGroup(ieff)=nlevel
         ife=0
         do ief=1,nbef
          do i=1,dpm%nlev(dpm%modele(ic,3+ief))
           ife=ife+1
           isol=isol+1

           if ( vecsol0(isol) ) then
             indest=indest+1
             incsol%paramaterValue(ieff,ife)  = par0_t(indest)
           else
             incsol%paramaterValue(ieff,ife)  = 0.d0
           end if
           incsol%parameterName(ieff,ife)   = trim(dpm%namefix(dpm%modele(ic,3+ief)))//' Level '//trim(str(i))
           incsol%parameterVecsol(ieff,ife) = vecsol0(isol)
           incsol%parameterPrecis(ieff,ife)  = precis0(isol)
          end do
         end do
       end if

       !Covariate
       if ( dpm%modele(ic,2) > 0 ) then
         ieff=ieff+1
         incsol%eqtl_print(ieff)=.false.
         incsol%groupeName(ieff) = 'Covariates'
         incsol%nbParameterGroup(ieff)=dpm%modele(ic,2)

         do ief=1,dpm%modele(ic,2)
           isol=isol+1
           if ( vecsol0(isol) ) then
             indest=indest+1
             incsol%paramaterValue(ieff,ief)  = par0_t(indest)
           else
             incsol%paramaterValue(ieff,ief)  = 0.d0
           end if
           incsol%parameterName(ieff,ief)   = trim(dpm%namecov(dpm%modele(ic,3+dpm%modele(ic,1)+ief)))
           incsol%parameterVecsol(ieff,ief) = vecsol0(isol)
           incsol%parameterPrecis(ieff,ief)  = precis0(isol)
          end do
       end if

       end subroutine set_solution_hypothesis0
!!***

!!****f* m_qtlmap_analyse_modlin/set_solution_hypothesis1
!!  NAME
!!    set_solution_hypothesis1
!!  DESCRIPTION
!!
!!
!!  NOTES
!!  SOURCE

     subroutine set_solution_hypothesis1(dataset,ic,incsol)
       type(QTLMAP_DATASET)               ,intent(in)       :: dataset
       integer                            ,intent(in)       :: ic
       type(TYPE_INCIDENCE_SOLUTION)      ,intent(inout)    :: incsol

       integer :: nteff,maxNbPar,jm,ifem,ip,ieff
       integer :: ntlev,nbtp,jef,lp,indest,km
       integer :: isol,ipar,nlevel,ief,i,ife

       real(kind=dp) ,dimension(size(par1)) :: par1_t

       type(GENEALOGY_BASE) , pointer :: dg
       type(PHENOTYPE_BASE) , pointer :: dpa
       type(DATAMODEL_BASE) , pointer :: dpm

       dpm => dataset%phenoModel
       dg => dataset%genea
       dpa => dataset%phenoAnimal

       ! remise a l echelle des variables
       do ipar=1,size(par1)
          par1_t(ipar)=par1(ipar)*dpm%sigt(ic)
       end do

       allocate (incsol%sig(1,dg%np))

       incsol%hypothesis=1
       !  Mean, Polygenic family, QTL effect
       nteff = 3
       if ( count(dpa%estime(ic,:))  > 0 ) nteff = 5

       !Fixed effect and covariate
       if ( dpm%modele(ic,1) > 0 ) nteff = nteff+1
       if ( dpm%modele(ic,2) > 0 ) nteff = nteff+1

       maxNbPar = max(dg%np,count(dpa%estime(ic,:)))
       !max numbre de covariable ?
       maxNbPar = max(maxNbPar,dpm%modele(ic,2))

       nlevel=0
       !max nombre de niveau pour un effet fixe ?
       do i=1,dpm%modele(ic,1)
         nlevel=dpm%nlev(dpm%modele(ic,3+i))+nlevel
       end do
       maxNbPar = max(maxNbPar,nlevel)

       ntlev=1
       nbtp = 3 + dpm%modele(ic,1)+dpm%modele(ic,2)
       do jef=1,dpm%modele(ic,3)
		    ntlev=ntlev*dpm%nlev(dpm%modele(ic,nbtp+jef))
	   end do

	   !max nombre de niveau pour un effet fixe en interaction avec le qtl ?
	   maxNbPar = max(maxNbPar,ntlev*dg%np)
	   maxNbPar = max(maxNbPar,ntlev*count(dpa%estime(ic,:)))

       allocate (incsol%groupeName(nteff))
       allocate (incsol%nbParameterGroup(nteff))
       allocate (incsol%eqtl_print(nteff))
       allocate (incsol%parameterName(nteff,maxNbPar))
       allocate (incsol%paramaterValue(nteff,maxNbPar))
       allocate (incsol%parameterVecsol(nteff,maxNbPar))
       allocate (incsol%parameterPrecis(nteff,maxNbPar))
       allocate (incsol%qtl_groupeName(1,1))
       incsol%parameterPrecis=0.d0
       incsol%parameterVecsol=.true.
       incsol%eqtl_print=.true.


       do ip=1,dg%np
            incsol%sig(1,ip) = par1_t(ip)
       end do

       ieff=1
       incsol%groupeName(ieff) = 'General Mean'
       incsol%nbParameterGroup(ieff)=1
       incsol%parameterName(ieff,1)   ='General Mean'
       incsol%paramaterValue(ieff,1)  = par1_t(dg%np+1)+dpm%xmut(ic)
       incsol%parameterVecsol(ieff,1) = vecsol1(1)
       incsol%parameterPrecis(ieff,1) = precis1(1)

       ieff=ieff+1
       incsol%qtl_groupeName(1,1)=ieff
       incsol%groupeName(ieff) = 'Sire QTL effects'
       incsol%nbParameterGroup(ieff)=dg%np*ntlevp

       indest = dg%np+1
       isol=1
       ife=0
       do ip=1,dg%np
         do lp=1,ntlev
           ife=ife+1
           isol=isol+1
           if ( vecsol1(isol) ) then
             indest=indest+1
             incsol%paramaterValue(ieff,ife)  = par1_t(indest)
           else
             incsol%paramaterValue(ieff,ife)  = 0.d0
           end if

           incsol%parameterName(ieff,ife)   ='Sire '//trim(dg%pere(ip))//" "//trim(str(lp))
           incsol%parameterVecsol(ieff,ife) = vecsol1(isol)
           incsol%parameterPrecis(ieff,ife)  = precis1(isol)
         end do
       end do

       if ( count(dpa%estime(ic,:)) > 0 ) then
           ieff = ieff +1
           incsol%groupeName(ieff) = 'Dam QTL effects'
           incsol%nbParameterGroup(ieff)=dpa%namest(ic)*ntlevp
           ifem=0
           do jm=1,dg%nm
             if ( dpa%estime(ic,jm)) then
              do lp=1,ntlev
               isol=isol+1
               ifem=ifem+1
               if ( vecsol1(isol) ) then
                indest=indest+1
                incsol%paramaterValue(ieff,ifem) = par1_t(indest)
               else
                incsol%paramaterValue(ieff,ifem) = 0.d0
               end if
               incsol%parameterName(ieff,ifem)='Dam '//trim(dg%mere(jm))//" "//trim(str(lp))
               incsol%parameterVecsol(ieff,ifem) = vecsol1(isol)
               incsol%parameterPrecis(ieff,ifem)  = precis1(isol)
              end do
             end if
           end do
         end if


       ieff=ieff+1
       incsol%groupeName(ieff) = 'Sire polygenic effects'
       incsol%nbParameterGroup(ieff)=dg%np

       do ip=1,dg%np
           isol=isol+1
           if ( vecsol1(isol) ) then
             indest=indest+1
             incsol%paramaterValue(ieff,ip)  = par1_t(indest)
           else
             incsol%paramaterValue(ieff,ip)  = 0.d0
           end if

           incsol%parameterName(ieff,ip)   ='Sire '//trim(dg%pere(ip))
           incsol%parameterVecsol(ieff,ip) = vecsol1(isol)
           incsol%parameterPrecis(ieff,ip)  = precis1(isol)
       end do

       if ( count(dpa%estime(ic,:)) > 0 ) then
           ieff = ieff +1
           incsol%groupeName(ieff) = 'Dam polygenic effects'
           incsol%nbParameterGroup(ieff)=count(dpa%estime(ic,:))
           ifem=0
          do ip=1,dg%np
           do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
             if ( dpa%estime(ic,jm)) then
              isol=isol+1
              ifem=ifem+1
              if ( vecsol1(isol) ) then
                indest=indest+1
                incsol%paramaterValue(ieff,ifem) = par1_t(indest)
              else
                incsol%paramaterValue(ieff,ifem) = 0.d0
              end if

              incsol%parameterName(ieff,ifem)='Dam '//trim(dg%mere(jm))//" [Sire "//trim(dg%pere(ip))//"]"
              incsol%parameterVecsol(ieff,ifem) = vecsol1(isol)
              incsol%parameterPrecis(ieff,ifem)  = precis1(isol)
             end if
           end do
          end do
         end if

       !Fixed effect
       if ( dpm%modele(ic,1) > 0 ) then
         ieff=ieff+1
         incsol%eqtl_print(ieff)=.false.
         incsol%groupeName(ieff) = 'Fixed effects'
         incsol%nbParameterGroup(ieff)=nlevel
         ife=0
         do ief=1,nbef
          do i=1,dpm%nlev(dpm%modele(ic,3+ief))
           isol=isol+1
           ife=ife+1
           if ( vecsol1(isol) ) then
             indest=indest+1
             incsol%paramaterValue(ieff,ife)  = par1_t(indest)
           else
             incsol%paramaterValue(ieff,ife)  = 0.d0
           end if
           incsol%parameterName(ieff,ife)   = trim(dpm%namefix(dpm%modele(ic,3+ief)))//' Level '//trim(str(i))
           incsol%parameterVecsol(ieff,ife) = vecsol1(isol)
           incsol%parameterPrecis(ieff,ife)  = precis1(isol)
          end do
         end do
       end if

       !Covariate
       if ( dpm%modele(ic,2) > 0 ) then
         ieff=ieff+1
         incsol%eqtl_print(ieff)=.false.
         incsol%groupeName(ieff) = 'Covariates'
         incsol%nbParameterGroup(ieff)=dpm%modele(ic,2)

         do ief=1,dpm%modele(ic,2)
           isol=isol+1
           if ( vecsol1(isol) ) then
             indest=indest+1
             incsol%paramaterValue(ieff,ief)  = par1_t(indest)
           else
             incsol%paramaterValue(ieff,ief)  = 0.d0
           end if
           incsol%parameterName(ieff,ief)    = trim(dpm%namecov(dpm%modele(ic,3+dpm%modele(ic,1)+ief)))
           incsol%parameterVecsol(ieff,ief)  = vecsol1(isol)
           incsol%parameterPrecis(ieff,ief)  = precis1(isol)
          end do
       end if

       end subroutine set_solution_hypothesis1
!!***

end module m_qtlmap_analyse_modlin
