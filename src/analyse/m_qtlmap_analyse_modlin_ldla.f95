!!****m* ANALYSE/m_qtlmap_analyse_modlin_ldla
!!  NAME
!!    m_qtlmap_analyse_modlin_ldla
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
module m_qtlmap_analyse_modlin_ldla
    use m_qtlmap_types
    use m_qtlmap_log
    use m_qtlmap_math
    use m_qtlmap_optimization
    use m_qtlmap_analyse_gen, only : sig0,xmu0p
    use m_qtlmap_analyse_lin_gen
    use m_qtlmap_output_handler
    use m_qtlmap_haplotype_ldla

    implicit none

    save

    type(GENEALOGY_BASE) , pointer :: p_dg
    type(PDD_BUILD)      , pointer :: p_spt
    type(PHENOTYPE_BASE) , pointer :: p_dpa
    type(HAPLOTYPE_POSITION_BUILD)   ,pointer      :: p_shp

!!****v* m_qtlmap_analyse_modlin_ldla/sig1
!!  NAME
!!   sig1
!!  DESCRIPTION
!!   The standart deviation under H0
!!  DIMENSIONS
!!   np
!!***
    real (kind=dp)       ,dimension(:),allocatable,private   :: sig1
!!****v* m_qtlmap_analyse_modlin_ldla/xmu1p
!!  NAME
!!   xmu1p
!!  DESCRIPTION
!!   The polygenic mean for each sire family under H0
!!  DIMENSIONS
!!   np
!!***
    real (kind=dp)       ,dimension(:),allocatable,private   :: xmu1p
!!****v* m_qtlmap_analyse_modlin_ldla/xmu1m
!!  NAME
!!   xmu1m
!!  DESCRIPTION
!!   The polygenic mean for each full sib family
!!  DIMENSIONS
!!   nm
!!***
    real (kind=dp)       ,dimension(:),allocatable,private   :: xmu1m
!!****v* m_qtlmap_analyse_modlin_ldla/xmu1g
!!  NAME
!!   xmu1g
!!  DESCRIPTION
!!   The genral mean
!!***
    real (kind=dp)                                ,private   :: xmu1g
!!****v* m_qtlmap_analyse_modlin_ldla/f0
!!  NAME
!!   f0
!!  DESCRIPTION
!!   value of the likelihood under H0
!!***
    real (kind=dp)                                ,public    :: f0
!!****v* m_qtlmap_analyse_modlin_ldla/fp0
!!  NAME
!!   fp0
!!  DESCRIPTION
!!   value of the likelihood by sire family under H0
!!  DIMENSIONS
!!   np
!!
!!***
    real (kind=dp)       ,dimension(:),allocatable,private   :: fp0
!!****v* m_qtlmap_analyse_modlin_ldla/fp1
!!  NAME
!!   fp1
!!  DESCRIPTION
!!   value of the likelihood by sire family under H1
!!  DIMENSIONS
!!   np
!!
!!***
    real (kind=dp)       ,dimension(:),allocatable,private   :: fp1
!!****v* m_qtlmap_analyse_modlin_ldla/fm0
!!  NAME
!!   fm0
!!  DESCRIPTION
!!   value of the likelihood by full-sib family under H0
!!  DIMENSIONS
!!   nm
!!
!!***
    real (kind=dp)       ,dimension(:),allocatable,private   :: fm0
!!****v* m_qtlmap_analyse_modlin_ldla/fm1
!!  NAME
!!   fm1
!!  DESCRIPTION
!!   value of the likelihood by full-sib family under H1
!!  DIMENSIONS
!!   nm
!!
!!***
    real (kind=dp)       ,dimension(:),allocatable,private   :: fm1

    integer              ,private                            :: current_chr

    logical              , private                           :: estime_moy = .true.
    logical              , private                           :: variance_homo = .false.
    logical              , private                           :: is_la = .false.
    logical              , private                           :: is_ld = .false.
    logical              , private                           :: is_ldla = .false.
    logical              , private                           :: is_ldjh = .false.

    integer             , private                            :: mod_nb_var
    integer             , private                            :: mod_pos
    integer             , private                            :: mod_ic
    integer             , private                            :: mod_imoy
    integer             , private                            :: mod_imoyld
    integer             , private                            :: mod_init
    integer             , private                            :: mod_jnit

    public :: init_analyse_modlin_ldla
    public :: opti_0qtl_modlin_ldla
    public :: opti_1qtl_modlin_ldla
    public :: test_lin_ldla
    public :: end_analyse_modlin_ldla

!
!**********************************************************************
    contains


!!****f* m_qtlmap_analyse_modlin_ldla/init_analyse_modlin_ldla
!!  NAME
!!    init_analyse_modlin_ldla
!!  DESCRIPTION
!!
!!  SOURCE
      subroutine init_analyse_modlin_ldla(dataset)
         type(QTLMAP_DATASET)       ,intent(in) :: dataset
         integer           :: stat

         type(GENEALOGY_BASE) , pointer :: dg
         dg => dataset%genea

         allocate (sig1(dg%np),STAT=stat)
         call check_allocate(stat,'sig1 [m_qtlmap_analyse_modlin_ldla]')
         allocate (xmu1p(dg%np),STAT=stat)
         call check_allocate(stat,'xmu1p [m_qtlmap_analyse_modlin_ldla]')
         allocate (xmu1m(dg%nm),STAT=stat)
         call check_allocate(stat,'xmu1m [m_qtlmap_analyse_modlin_ldla]')
         allocate (fp0(dg%np),STAT=stat)
         call check_allocate(stat,'fp0 [m_qtlmap_analyse_modlin_ldla]')
         allocate (fp1(dg%np),STAT=stat)
         call check_allocate(stat,'fp1 [m_qtlmap_analyse_modlin_ldla]')
         allocate (fm0(dg%nm),STAT=stat)
         call check_allocate(stat,'fm0 [m_qtlmap_analyse_modlin_ldla]')
         allocate (fm1(dg%nm),STAT=stat)
         call check_allocate(stat,'fm1 [m_qtlmap_analyse_modlin_ldla]')

     end subroutine init_analyse_modlin_ldla
!!***

!!****f* m_qtlmap_analyse_modlin_ldla/end_analyse_modlin_ldla
!!  NAME
!!    end_analyse_modlin_ldla
!!  DESCRIPTION
!!
!!  SOURCE
     subroutine end_analyse_modlin_ldla
         deallocate (sig1)
         deallocate (xmu1p)
         deallocate (xmu1m)
         deallocate (fp0)
         deallocate (fp1)
         deallocate (fm0)
         deallocate (fm1)
     end subroutine end_analyse_modlin_ldla
!!***

!!****f* m_qtlmap_analyse_modlin_ldla/opti_0qtl_modlin_ldla
!!  NAME
!!    opti_0qtl_modlin_ldla
!!  DESCRIPTION
!!
!! NOTES
!!   Calcul de la vraisemblance 0 QTL, 1 caractere , effets parasites inclus
!!  SOURCE
      subroutine opti_0qtl_modlin_ldla(dataset,ic,est_moy,option_var)
      use m_qtlmap_analyse_lin_gen, only : contingence,precision,vecsol,precis, &
                                    nbnivest,ntniv,par0,precis0,vecsol0,xx

      implicit none

      integer , intent(in)   :: ic
      type(QTLMAP_DATASET)       ,intent(in) :: dataset
      logical , intent(in)   :: est_moy
      character(len=4)         , intent(in)               :: option_var
!
! Divers

      integer                                    :: iuser(2)
      integer        ,dimension(:),allocatable   :: iw
      real (kind=dp) ,dimension(:),allocatable   :: par,borni,borns,w
      real (kind=dp)                             :: user(1)
      real (kind=dp)                             :: f,var_yQy
      logical  , dimension(:,:,:),allocatable    :: filter_inc

      integer :: npar,ibound,ip,ix,ifail,i,indest,ntot
      integer :: km,jm,imoy
      integer :: nb_var,nparx
      logical itest,details

      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa
      type(DATAMODEL_BASE) , pointer :: dpm

      type(PDD_BUILD)                   ,pointer :: spt => null()
      type(HAPLOTYPE_POSITION_BUILD)   ,pointer :: shp => null()

      dpm => dataset%phenoModel
      dg => dataset%genea
      dpa => dataset%phenoAnimal

      p_dg => dataset%genea
      p_dpa => dataset%phenoAnimal
!
!******************************************************************************
!
!  recherche des elements estimables � partir de la d�composition de choleski
!  comptage et recodification de ces elements (on commence par mettre � 0 les
!  compteurs destin�s aux tests des effes fix�s
!
      nb_var=1
      variance_homo = (option_var == 'homo')

      if(option_var == 'hete') nb_var=dg%np
     ! if(option_anal == 'LD  ' .or. option_anal == 'LDLA')nb_var=2*nb_var

      estime_moy = est_moy

      itest=.false.
      details=.false.
      call contingence_ldla(dataset,spt,shp,ic,0,itest,1,0,var_yQy,details,est_moy,"LA  ",.false.,.false.)
      call precision(xx,precis)

!******************
! passage des arguments par iuser
!   iuser(1) = nb_var
!   iuser(2) = imoy
      iuser(1) = nb_var
      iuser(2) = 0
        if(est_moy) iuser(2) = 1
      imoy=iuser(2)


      ibound=0
      npar=nb_var+nbnivest

      allocate (borni(npar))
      allocate (borns(npar))
      allocate (par(npar))


      !MODIF - OPMIZATION
      allocate (filter_inc(dg%np,dg%nm, npar))
      call set_filter_optim(dataset,ic,(option_var == 'hete'),.false.,ntnivmax,ntniv,vecsol,xinc,filter_inc)
      !FIN MODIF - OPMIZATION

      do ip=1,nb_var
        borni(ip)=SIG_MIN
        borns(ip)=SIG_MAX
      end do
      do ix=nb_var+1,npar
        borni(ix)=XMU_MIN
        borns(ix)=XMU_MAX
      end do
!
! Point de depart
      par=0.d0
      par(dg%np+1)=0.d0
      do ip=1,dg%np
        if (option_var == 'hete' ) then
          par(ip)=sig0(ip)
        else
          par(1)=sig0(ip)+par(1)
        end if
        par(dg%np+1)=par(dg%np+1)+xmu0p(ip)
      end do
      par(dg%np+1)=par(dg%np+1)/dble(dg%np)
      if (option_var == 'homo' )  par(1) =  par(1) / dble(dg%np)
      do ix=nb_var+2,npar
        par(ix)=0.d0
      end do
   !   filter_inc=.true.
!
! Optimisation de la vraisemblance
      ifail=1
      mod_ic = ic ! pour le mode lineaire....
  !    call minimizing_funct(npar,ibound,funct_0qtl_modlin_ldla,borni,borns,par,f,iuser,user,ifail)
      call minimizing_funct_family(dataset,npar,ibound,funct_0qtl_modlin_ldla_family,filter_inc,fm0,fp0,&
        borni,borns,par,f,iuser,user,ifail)
     ! if (ifail.ne.0) print *,'Code retour optimizing 0 QTL : ',ifail

      f0=f
      do i=1,npar
        par0(i)=par(i)
      end do
      do i=1,ntniv
        vecsol0(i)=vecsol(i)
        precis0(i)=precis(i)
      end do

      if ( option_var == 'hete' ) then
       do ip = 1,dg%np
        sig1(ip)=par(ip)
       end do
      else
        sig1 =par(1)
      end if

      xmu1g=par(nb_var+1)
      indest=1
      do ip=1,dg%np
      if (vecsol(1+ip)) then
        indest=indest+1
        xmu1p(ip)=par(nb_var+indest)
      end if
      end do

      ntot=nb_var+1
      km=0
      do jm=1,dg%nm
        if (dpa%estime(ic,jm))then
          km=km+1
          if(vecsol(ntot+km)) then
            indest=indest+1
            xmu1m(jm)=par(nb_var+indest)
          end if
        end if
      end do

      deallocate (borni)
      deallocate (borns)
      deallocate (par)
      deallocate (filter_inc)

      end subroutine opti_0qtl_modlin_ldla
!!***

!!****if* m_qtlmap_analyse_modlin_ldla/funct_0qtl_modlin_ldla
!!  NAME
!!    funct_0qtl_modlin_ldla
!!  DESCRIPTION
!!
!! NOTES
!!   Calcul de la vraisemblance et de ses derivees partielles sous H0
!!  SOURCE
      subroutine funct_0qtl_modlin_ldla(n,x,f,iuser,user)
      use m_qtlmap_analyse_lin_gen, only : nbfem,nbef,ntnifix,nbco,covdir,nivdir,corniv,vecsol

      implicit none
      integer         , intent(in)                  :: n
!
! Tableaux dimensionnes selon n, le nombre de parametres a estimer
!
      real (kind=dp)      ,dimension(n), intent(in) :: x
      real (kind=dp)  , intent(inout)               :: f
      integer ,       dimension(2)     :: iuser
      real (kind=dp)      ,dimension(1), intent(in) :: user

! Tableaux dimensionnes selon nm, le nombre de meres
      real (kind=dp)   ,dimension(p_dg%nm)          :: effm

      integer ::jnit,ip,nm1,nm2,jm,nd1,nd2,kkd,ilev,ief,ico,imoy
      integer :: nb_var
      real (kind=dp) :: sig,var,vmere,vpf,v
!

!******************************************************************************
!

!******************
! passage des arguments par iuser
!   iuser(1) = nb_var
!   iuser(2) = imoy
      nb_var = iuser(1)
      imoy = iuser(2)

      jnit=3
      if (nbfem.eq.0)jnit=2
      if(iuser(1)== 0) jnit=jnit-1


      f=0.d0
      if(variance_homo) then
        sig=x(1)
        var=sig*sig
      end if

      do ip=1,p_dg%np
        if(.not. variance_homo) then
          sig=x(ip)
          var=sig*sig
        end if
        fp0(ip)=0.d0
        do jm=p_dg%nmp(ip)+1,p_dg%nmp(ip+1)
           fm0(jm)=0.d0
          do kkd=p_dg%ndm(jm)+1,p_dg%ndm(jm+1)
!
            if(p_dpa%presentc(mod_ic,kkd)) then
              v=p_dpa%y(mod_ic,kkd)
!
              if(imoy ==1 .and. vecsol(nivdir(kkd,1))) v=v-x(nb_var+corniv(nivdir(kkd,1)))
!
              if(vecsol(nivdir(kkd,imoy+1))) v=v-x(nb_var+corniv(nivdir(kkd,imoy+1)))
!
              if(p_dpa%estime(mod_ic,jm)) then
               if(vecsol(nivdir(kkd,imoy+2))) v=v-x(nb_var+corniv(nivdir(kkd,imoy+2)))
              end if

              do ief=jnit+1,jnit+nbef
               if(vecsol(nivdir(kkd,ief))) v=v-x(nb_var+corniv(nivdir(kkd,ief)))
             end do

              do ico=1,nbco
                if(vecsol(ntnifix+ico))  v=v-covdir(kkd,ico)*x(nb_var+corniv(ntnifix+ico))
              end do

                fm0(jm)=fm0(jm)+v*v*p_dpa%cd(mod_ic,kkd)/var + dlog(var)
            end if!presentc
          end do!kkd

          fp0(ip)=fp0(ip)+fm0(jm)
        end do !jm
        f=f+fp0(ip)
      end do ! ip
      return
      end subroutine funct_0qtl_modlin_ldla
!!***

!!****f* m_qtlmap_analyse_modlin_ldla/funct_0qtl_modlin_ldla_family
!!  NAME
!!    funct_0qtl_modlin_ldla_family
!!  DESCRIPTION
!!
!! NOTES
!!   Calcul de la vraisemblance et de ses derivees partielles sous H0
!!  SOURCE
     subroutine funct_0qtl_modlin_ldla_family(ip,jm,n,x,f,iuser,user)
      use m_qtlmap_analyse_lin_gen, only : nbfem,nbef,ntnifix,nbco,covdir,nivdir,corniv,vecsol

      implicit none
      integer         , intent(in)                  :: ip,jm,n
!
! Tableaux dimensionnes selon n, le nombre de parametres a estimer
!
      real (kind=dp)      ,dimension(n), intent(in) :: x
      real (kind=dp)  , intent(inout)               :: f
      integer ,       dimension(2)                  :: iuser
      real (kind=dp)      ,dimension(1), intent(in) :: user

! Tableaux dimensionnes selon nm, le nombre de meres
      real (kind=dp)   ,dimension(p_dg%nm)               :: effm

      integer ::jnit,nm1,nm2,nd1,nd2,kkd,ilev,ief,ico,imoy
      integer :: nb_var
      real (kind=dp) :: sig,var,vmere,vpf,v
!

!******************************************************************************

!******************
! passage des arguments par iuser
!   iuser(1) = nb_var
!   iuser(2) = imoy
      nb_var = iuser(1)
      imoy = iuser(2)

      jnit=3
      if (nbfem.eq.0)jnit=2
      if(iuser(1)== 0) jnit=jnit-1


      f=0.d0
      if(variance_homo) then
        sig=x(1)
        var=sig*sig
      else
        sig=x(ip)
        var=sig*sig
      end if

      do kkd=p_dg%ndm(jm)+1,p_dg%ndm(jm+1)
!
            if(p_dpa%presentc(mod_ic,kkd) ) then
              v=p_dpa%y(mod_ic,kkd)
!
              if(imoy ==1 .and. vecsol(nivdir(kkd,1))) v=v-x(nb_var+corniv(nivdir(kkd,1)))
!
              if(vecsol(nivdir(kkd,imoy+1))) v=v-x(nb_var+corniv(nivdir(kkd,imoy+1)))
!
              if(p_dpa%estime(mod_ic,jm)) then
                if(vecsol(nivdir(kkd,imoy+2))) v=v-x(nb_var+corniv(nivdir(kkd,imoy+2)))
              end if

              do ief=jnit+1,jnit+nbef
                if(vecsol(nivdir(kkd,ief))) v=v-x(nb_var+corniv(nivdir(kkd,ief)))
              end do

              do ico=1,nbco
                if(vecsol(ntnifix+ico))  v=v-covdir(kkd,ico)*x(nb_var+corniv(ntnifix+ico))
              end do

                f=f+v*v*p_dpa%cd(mod_ic,kkd)/var + dlog(var)
            end if!presentc
          end do!kkd
      return
      end subroutine funct_0qtl_modlin_ldla_family
!!***


!!****f* m_qtlmap_analyse_modlin_ldla/opti_1qtl_modlin_ldla
!!  NAME
!!    opti_1qtl_modlin_ldla
!!  DESCRIPTION
!!    Calcul de la statistique de test le long du chromosome
!! NOTES
!!
!!  SOURCE
   subroutine opti_1qtl_modlin_ldla(dataset,spt,shp,ic,lrtsol,fmax,supnbnivest,est_moy,option_var,option_anal,hdam)
        integer , intent(in)                                :: ic
        type(QTLMAP_DATASET)              ,intent(in)       :: dataset
        type(HAPLOTYPE_POSITION_BUILD) ,target,intent(inout):: shp
        type(PDD_BUILD)        ,target    ,intent(in)       :: spt

        type(TYPE_LRT_SOLUTION)  , intent(out)              :: lrtsol
        real (kind=dp)  , intent(out)                       :: fmax
        integer         , intent(out)                       :: supnbnivest
        logical         , intent(in)                        :: est_moy
        character(len=4)         , intent(in)               :: option_var,option_anal
        logical                  ,intent(in)                :: hdam
!

! Divers
      real (kind=dp) ,dimension(:),allocatable :: val,par,borni,borns,par2
      real (kind=dp) :: user(1)
      real (kind=dp) :: xlrt_t,dx,f1,var_yQy,weight
      real (kind=dp) :: xlrt_t21,xlrt_t20,f2,f21

      logical itest,stvecsol(size(vecsol1))
      integer :: ip,i,n,ilong,ibound,npar,nparx,ix,j,ifail,iuser(0),chr
      integer :: ii,imoy,kd,nstart,nend
      logical :: details,hsire
      integer :: nb_var,jm

      integer :: i_haplo,ig,kkd

      logical  , dimension(:,:,:),allocatable        :: filter_inc
      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa
      type(DATAMODEL_BASE) , pointer :: dpm

      real(kind=dp)       ,dimension(:,:,:)   ,pointer     :: xlrp,xlrm
      real(kind=dp)       ,dimension(:,:)   ,pointer       :: lrt1

      dpm => dataset%phenoModel
      dg => dataset%genea
      dpa => dataset%phenoAnimal

      p_dpa => dataset%phenoAnimal
      p_dg => dataset%genea
      p_spt => spt
      p_shp => shp

      call p_shp%set(dataset,spt)

      allocate (lrt1(dataset%map%nchr,dataset%map%get_maxnpo()))
      allocate (xlrp(dataset%map%nchr,dataset%map%get_maxnpo(),dg%np))
      allocate (xlrm(dataset%map%nchr,dataset%map%get_maxnpo(),dg%nm))

      variance_homo = (option_var == 'homo')
      is_la = (option_anal == 'LA  ')
      is_ld = (option_anal == 'LD  ')
      is_ldla = (option_anal == 'LDLA')
      is_ldjh = (option_anal == 'LDJH')
      hsire = ( is_ld .or. is_ldla .or. is_ldjh )

     ! open(1212,file='lrt_LDJH')
      nb_var=1
      if(option_var == 'hete') nb_var=dg%np
  !    if(option_anal == 'LD  ' .or. option_anal == 'LDLA')nb_var=2*nb_var

      nparx = 12*(nb_var+ntnivmax)+((nb_var+ntnivmax)*((nb_var+ntnivmax)-1)/2)
      allocate ( val(nb_var+ntnivmax ) )
      allocate ( par(nb_var + ntnivmax) )
      allocate ( par2(nb_var + ntnivmax) )
      allocate ( borni(nb_var + ntnivmax) )
      allocate ( borns(nb_var + ntnivmax) )
      iuser=0

      estime_moy = est_moy

   !  open(66,file='compALvsJM')
      !MODIF - OPMIZATION
      allocate (filter_inc(dg%np,dg%nm, nb_var + ntnivmax))
      !FIN MODIF - OPMIZATION


!
!******************************************************************************
! Calcul de la vraisemblance sous H1

!
! initialisation
!  on utilisera  :
!  PAR (vecteur des param�tres � optimiser),
!  CORNIV (vecteur des positions, parmi NTNIV, des effets extimables),
!  SOLVEC (vecteur disant l'estimabilit� des effets)
!  VAL le vecteur complet (ntniv positions) des valeurs des niveaux des effets
!  STVAL,STVECSOL et STCORNIV les copies de VAL, VECSOL et CORNIV � l'it�ration
!  n-1 quand on analyse la position n
!
! Point de depart
!

      do ip=1,dg%np
        par(ip)=sig1(ip)
      end do

      par(dg%np+1)=xmu1g

      do i=dg%np+2,size(par)
        par(i)=0.d0
      end do

      do i=dg%np+1,size(val)
        val(i-dg%np)=par(i)
      end do
!
      do i=1,size(vecsol)
        vecsol(i)=.true.
      end do

      call lrtsol%new(dataset,1)
      lrtsol%lrtmax(0)=-1.d75
      lrtsol%nxmax(0)=1
      lrtsol%chrmax(0)=1
      lrt1=0.d0

      mod_ic = ic
      do chr=1,dataset%map%nchr
       current_chr = chr
      ! call  set_tab_ibs(chr) !peut etre a faire que pour LDJH
!
! Marche le long du chromosome
      !n=0

      nstart=1
      do while (dataset%map%absi(chr,nstart)<=dataset%map%posi(chr,ceiling(dataset%params%longhap/2.d0)))
        nstart=nstart+1
      end do



      nend=dataset%map%get_npo(chr)
      do while (dataset%map%absi(chr,nend)>= &
       dataset%map%posi(chr,dataset%map%nmk(chr)-ceiling(dataset%params%longhap/2.d0)+1))
          nend = nend - 1
      end do

 !     nend=nstart
      do n=nstart,nend
        
        dx=dataset%map%absi(chr,n)
!
!  recherche des haplotypes possibles
!
      if(option_anal == 'LD  ' .or. option_anal == 'LDLA' ) &
       call shp%set_haplo_for_ldla(chr,dx,n,hsire,hdam)

      if(option_anal == 'LDJH' )then
            call shp%liste_shared_haplo(chr,dx,n)
            do ip = 1,dg%np
       !       call bilan_shared_haplo(current_chr,dx,n,ip)
            end do
      end if
!
!  on stocke les conditions en n
        do i=1,ntniv
          stvecsol(i)=vecsol(i)
        end do
!
!  preparation de la matrice d'incidence
!
        call prepinc(dataset,spt,current_chr,n,ic,option_anal)
!  recherche des elements estimables � partir de la d�composition de choleski
!  comptage et recodification de ces elements (on commence par mettre � 0 les
!  compteurs destin�s aux tests des effes fix�s
!
        itest=.false.
        details=.false.
        call contingence_ldla(dataset,spt,shp,ic,1,itest,chr,n,var_yQy,details,est_moy,option_anal,hsire,hdam)
       !MODIF - OPMIZATION
       call set_filter_optim(dataset,ic,(option_var == 'hete'),(option_anal == 'LD  ' .or. option_anal == 'LDLA'),&
                ntnivmax,ntniv,vecsol,xinc,filter_inc)
       !FIN MODIF - OPMIZATION


 !       weight=dsqrt(dlog(var_yQy))
!        print*,dx,weight
!        go to 111
!******************
! passage des arguments par iuser
!   iuser(1) = nb_var
!   iuser(2) = n
!   iuser(3) = imoy
!   iuser(4) = imoyld
!   iuser(5) = init
!   iuser(6) = jnit
       mod_nb_var =  nb_var
       mod_pos = n
       mod_imoy = 0

       if(est_moy) mod_imoy = 1

       mod_imoyld=mod_imoy

       if(option_anal == 'LD  ' .or. option_anal == 'LDLA')mod_imoyld=mod_imoy+shp%nb_haplo_reduit

       mod_init=1  ! +1 par rapport au dernier effet estime

       if(mod_imoy==1)mod_init=mod_init+1 ! moyenne

       if(option_anal /= 'LD  ')mod_init=mod_init+1 ! qtl pere

       if(option_anal == 'LDJH')mod_init=mod_init+1 ! qtl pere transmis par la mere

       if(nbfem /= 0  .and. (option_anal == 'LA  ' .or. option_anal == 'LDLA'))mod_init=mod_init+1 !qtl mere
!        if(option_anal == 'LD  ' .or. option_anal == 'LDLA')iuser(5)=iuser(5)+nb_max_haplo-1 ! haplotype
       if(option_anal == 'LD  ' .or. option_anal == 'LDLA')mod_init=mod_init+shp%nb_haplo_reduit ! haplotype
       mod_jnit=mod_init
       if(nbfem /= 0 )mod_jnit=mod_jnit+1




! Parametres de maximisation
      ibound=0
      npar=nb_var+nbnivest
      do ip=1,nb_var
        borni(ip)=SIG_MIN
        borns(ip)=SIG_MAX
      end do
      do i=nb_var+1,npar
        borni(i)=XMU_MIN
        borns(i)=XMU_MAX
      end do
!
! Point de depart (on reprend les point d'arriv�e pr�c�dents
!
      do i=nb_var+1,npar
        par(i)=0.d0
      end do

  !    j=nb_var
  !    do i=1,ntniv
  !      if(vecsol(i))then
  !        j=j+1
  !        par(j)=0.d0
  !        if(stvecsol(i)) par(j)=val(i)
  !      end if
  !    end do


!      write(66,*)' a la position dx = ', dx
!      do ip=1,np
!          do jm=nmp(ip)+1,nmp(ip+1)
!            do ig=ngenom(current_chr,jm)+1,ngenom(current_chr,jm+1)
!              do kd=ngend(current_chr,ig)+1,ngend(current_chr,ig+1)
!              kkd=ndesc(current_chr,kd)
!              if(presentc(mod_ic,kkd)) then
!       write(66,*) n,ip,jm,ig,kd,animal(kkd),(trim(name_haplo_reduit(i_haplo)),&
!                     pb_haplo(kkd,i_haplo),i_haplo=1,nb_haplo_reduit)
!              end if
!              end do !kd
!             end do !ig
!           end do !jm
!       end do !ip

! Optimisation de la vraisemblance a la position dx pour le modele sans transmission par la mere
        if(option_anal == 'LDJH') then
        ifail=1
!        iuser(7)=0
!        iuser(8)=0
        call minimizing_funct(dataset,npar,ibound,funct_1qtl_modlin_ldla,borni,borns,par,f2,iuser,user,ifail)
        end if
! Optimisation de la vraisemblance a la position dx
        ifail=1
     !   iuser(7)=1
     !   call minimizing_funct(npar,ibound,funct_1qtl_modlin_ldla,borni,borns,par,f1,iuser,user,ifail)
        call minimizing_funct_family(dataset,npar,ibound,funct_1qtl_modlin_ldla_family,filter_inc(:,:,:npar),&
                                    fm1,fp1,borni,borns,par,f1,iuser,user,ifail)

        xlrt_t=f0-f1
        print*,n,dx,xlrt_t

 !       if(xlrt_t  < 0.d0) go to 111

       if(xlrt_t  < 0.d0) then
         print*,'n, nb_haplo_reduit',n, shp%nb_haplo_reduit
         print *,'par',(par(i),i=nb_var+2,nb_var+1+shp%nb_haplo_reduit)
      end if

      j=nb_var
      do i=1,ntniv
        if(vecsol(i)) then
          j=j+1
          val(i)=par(j)
        else
          val(i)=9999.d0
        end if
      end do
!
!
!  on met les profil / progeniteur dans les fichier ad hoc
!
        do ii=1,dg%np
          xlrp(chr,n,ii)=-1.d0*(fp1(ii)-fp0(ii))
        end do
        do ii=1,dg%nm
          xlrm(chr,n,ii)=-1.d0*(fm1(ii)-fm0(ii))
        end do

!  on stocke la position si elle est meilleure que les pr�c�dents
!
        if(lrtsol%lrtmax(0) < xlrt_t) then
!          dxmax=dx
          lrtsol%nxmax(0)=n
          lrtsol%chrmax(0)=chr
          lrtsol%lrtmax(0)=xlrt_t
          fmax=f1
          supnbnivest=nbnivest
          do i=1,npar
            par1(i)=par(i)
          end do

        end if
!
! on garde les valeurs du LRT pour pouvopir dessiner la courbe de vraisemblance
!
        lrt1(chr,n)=xlrt_t
  111 continue
      end do

      end do

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

!
!  on calcule la pr�cision des estimation au point correspondant
! au LRT maximum
!
      ! hsire=  true pour recuperer les haplotypes des peres
      call shp%set_haplo_for_ldla(lrtsol%chrmax(0),dataset%map%absi(lrtsol%chrmax(0),&
      lrtsol%nxmax(0)),lrtsol%nxmax(0),.true.,hdam)

      call prepinc(dataset,spt,lrtsol%chrmax(0),lrtsol%nxmax(0),ic,option_anal)
      details=.false.
      call contingence_ldla(dataset,spt,shp,ic,1,itest,lrtsol%chrmax(0),lrtsol%nxmax(0),&
      var_yQy,details,est_moy,option_anal,hsire,hdam)
      call precision(xx,precis)

      do i=1,ntniv
        precis1(i)=precis(i)
        vecsol1(i)=vecsol(i)
      end do

      deallocate ( val )
      deallocate ( par )
      deallocate ( borni )
      deallocate ( borns )
      deallocate (filter_inc)

      return
      end subroutine opti_1qtl_modlin_ldla
!!***


!!****f* m_qtlmap_analyse_modlin_ldla/funct_1qtl_modlin_ldla
!!  NAME
!!    funct_1qtl_modlin_ldla
!!  DESCRIPTION
!!
!! NOTES
!!   Calcul de la vraisemblance et de ses derivees partielles sous Hypoth�se 1QTL 1carac
!!  SOURCE
      subroutine  funct_1qtl_modlin_ldla(n,x,f,iuser,user)
      implicit none

      integer         , intent(in)                  :: n
!
! Tableaux dimensionnes selon n, le nombre de parametres a estimer
!
      real (kind=dp)      ,dimension(n), intent(in) :: x
      real (kind=dp)  , intent(inout)               :: f
      integer ,       dimension(8)     :: iuser
      real (kind=dp)      ,dimension(1), intent(in) :: user

! Tableaux dimensionnes selon nm, le nombre de meres
      real (kind=dp)   ,dimension(p_dg%nm)   :: effm


!
! Divers
      real (kind=dp) :: sig,var,vmere,vpf,v
      real (kind=dp) :: sigc,varc,vart,det
      integer :: ip,nm1,nm2,jm,ngeno1,ngeno2,ig
      integer :: nd1,nd2,kd,kkd,ilev,ief,ico
      integer :: i_haplo
      integer :: imp_funct,i

!******************************************************************************
!
!******************
! passage des arguments par iuser
!   iuser(1) = nb_var
!   iuser(2) = n
!   iuser(3) = imoy
!   iuser(4) = imoyld
!   iuser(5) = init
!   iuser(6) = jnit


      f=0.d0
      if(variance_homo) then
        sig=x(1)      
        var=sig*sig
      end if

      do ip=1,p_dg%np
        if(.not. variance_homo) then
          sig=x(ip)
          var=sig*sig
        end if

        fp1(ip)=0.d0
         do jm=p_dg%nmp(ip)+1,p_dg%nmp(ip+1)
          vmere=0.d0
          det=0.d0
          do ig=p_spt%ngenom(current_chr,jm)+1,p_spt%ngenom(current_chr,jm+1)
            vpf=0.d0
            do kd=p_spt%ngend(current_chr,ig)+1,p_spt%ngend(current_chr,ig+1)
              kkd=p_spt%ndesc(current_chr,kd)

              if(p_dpa%presentc(mod_ic,kkd)) then

                v=p_dpa%y(mod_ic,kkd)
                if(estime_moy) v=v-x(mod_nb_var+corniv(nivdir(kkd,1)))
!
! effets haplotypes
                if(is_ld .or. is_ldla ) then
                 do i_haplo=1,p_shp%nb_haplo_reduit
                   if(vecsol(nivdir(kkd,mod_imoy+i_haplo))) &
                     v=v-pb_haplo(kkd,i_haplo)*x(mod_nb_var+corniv(nivdir(kkd,mod_imoy+i_haplo)))
                 end do
                end if ! option

                if(.not. is_ld)then
! effet qtl pere
                  if(vecsol(nivdir(kkd,mod_imoyld+1))) v=v-x(mod_nb_var+corniv(nivdir(kkd,mod_imoyld+1)))*ppt(kd)
!
! effet qtl mere
                  if(.not. is_ldjh) then
                     if (p_dpa%estime(mod_ic,jm)) then
                       if (vecsol(nivdir(kkd,mod_imoyld+2))) v=v-x(mod_nb_var+corniv(nivdir(kkd,mod_imoyld+2)))*pmt(kd)

                     end if
                  end if
! effet qtl mere ou  qtl pere transmis par la mere (modele JH)

                  if(is_ldjh .and. iuser(7) ==1) then
                    if(vecsol(nivdir(kkd,mod_imoy+2)).and. nivdir(kkd,mod_imoy+2) /= 0) &
                      v=v-x(mod_nb_var+corniv(nivdir(kkd,mod_imoy+2)))*(1-2*modulo(p_shp%shared_haplo(kkd),2))
                   end if
                end if ! option
!
! effet polygenique pere
                if(vecsol(nivdir(kkd,mod_init))) v=v-x(mod_nb_var+corniv(nivdir(kkd,mod_init)))
!
! effet polygenique mere
               if(p_dpa%estime(mod_ic,jm) ) then
                 if ( vecsol(nivdir(kkd,mod_init+1))) v=v-x(mod_nb_var+corniv(nivdir(kkd,mod_init+1)))
               end if
!
! effets fixes
                do ief=mod_jnit+1,mod_jnit+nbef
                 if(vecsol(nivdir(kkd,ief))) v=v-x(mod_nb_var+corniv(nivdir(kkd,ief)))
                end do
!
! covariables
                do ico=1,nbco
                  if(vecsol(ntnifix+ico))  v=v-covdir(kkd,ico)*x(mod_nb_var+corniv(ntnifix+ico))
                end do

                vpf=vpf+v*v*p_dpa%cd(mod_ic,kkd)/var
                det=det+dlog(var)
              end if !presentc
            end do ! kd

            vmere=vmere+p_spt%probg(current_chr,ig)*dexp(-0.5d0*vpf)
          end do !ig


          if (vmere == 0) then
                 fm1(jm)=INIFINY_REAL_VALUE
          else
                 fm1(jm)=-2.d0*dlog(vmere)+det
          end if

          fp1(ip)=fp1(ip)+fm1(jm)
        end do !jm
        f=f+fp1(ip)
      end do !ip

      return
      end subroutine funct_1qtl_modlin_ldla
!!***


!!****f* m_qtlmap_analyse_modlin_ldla/funct_1qtl_modlin_ldla_family
!!  NAME
!!    funct_1qtl_modlin_ldla_family
!!  DESCRIPTION
!!
!! NOTES
!!
!!  SOURCE
      subroutine  funct_1qtl_modlin_ldla_family(ip,jm,n,x,f,iuser,user)
      implicit none

      integer         , intent(in)                  :: n,ip,jm
!
! Tableaux dimensionnes selon n, le nombre de parametres a estimer
!
      real (kind=dp)      ,dimension(n), intent(in) :: x
      real (kind=dp)  , intent(inout)               :: f
      integer ,       dimension(8)     :: iuser
      real (kind=dp)      ,dimension(1), intent(in) :: user


!
! Divers
      real (kind=dp) :: sig,var,vmere,vpf,v
      real (kind=dp) :: sigc,varc,vart,det
      integer :: nm1,nm2,ngeno1,ngeno2,ig
      integer :: nd1,nd2,kd,kkd,ilev,ief,ico
      integer :: i_haplo
      integer :: imp_funct,i
!******************************************************************************
!
!******************
! passage des arguments par iuser
!   iuser(1) = nb_var
!   iuser(2) = n
!   iuser(3) = imoy
!   iuser(4) = imoyld
!   iuser(5) = init
!   iuser(6) = jnit

      f=0.d0
      if(variance_homo) then
        sig=x(1)      
        var=sig*sig
      end if

      
      if(.not. variance_homo) then
          sig=x(ip)
          var=sig*sig
        end if
         
          vmere=0.d0
          det=0.d0
          do ig=p_spt%ngenom(current_chr,jm)+1,p_spt%ngenom(current_chr,jm+1)
            vpf=0.d0
            do kd=p_spt%ngend(current_chr,ig)+1,p_spt%ngend(current_chr,ig+1)
              kkd=p_spt%ndesc(current_chr,kd)

              if(p_dpa%presentc(mod_ic,kkd)) then

                v=p_dpa%y(mod_ic,kkd)
                if(estime_moy) v=v-x(mod_nb_var+corniv(nivdir(kkd,1)))
!
! effets haplotypes
                if(is_ld .or. is_ldla ) then
                 do i_haplo=1,p_shp%nb_haplo_reduit
                   if(vecsol(nivdir(kkd,mod_imoy+i_haplo))) &
                     v=v-pb_haplo(kkd,i_haplo)*x(mod_nb_var+corniv(nivdir(kkd,mod_imoy+i_haplo)))
                 end do
                end if ! option

                if(.not. is_ld)then
! effet qtl pere
                  if(vecsol(nivdir(kkd,mod_imoyld+1))) v=v-x(mod_nb_var+corniv(nivdir(kkd,mod_imoyld+1)))*ppt(kd)
!
! effet qtl mere
                  if(.not. is_ldjh) then
                     if (p_dpa%estime(mod_ic,jm)) then
                       if (vecsol(nivdir(kkd,mod_imoyld+2))) v=v-x(mod_nb_var+corniv(nivdir(kkd,mod_imoyld+2)))*pmt(kd)

                     end if
                  end if
! effet qtl mere ou  qtl pere transmis par la mere (modele JH)

                  if(is_ldjh .and. iuser(7) ==1) then
                    if(vecsol(nivdir(kkd,mod_imoy+2)).and. nivdir(kkd,mod_imoy+2) /= 0) &
                      v=v-x(mod_nb_var+corniv(nivdir(kkd,mod_imoy+2)))*(1-2*modulo(p_shp%shared_haplo(kkd),2))
                   end if
                end if ! option
!
! effet polygenique pere
                if(vecsol(nivdir(kkd,mod_init))) v=v-x(mod_nb_var+corniv(nivdir(kkd,mod_init)))
!
! effet polygenique mere
               if(p_dpa%estime(mod_ic,jm) ) then
                 if ( vecsol(nivdir(kkd,mod_init+1))) v=v-x(mod_nb_var+corniv(nivdir(kkd,mod_init+1)))
               end if
!
! effets fixes
                do ief=mod_jnit+1,mod_jnit+nbef
                 if(vecsol(nivdir(kkd,ief))) v=v-x(mod_nb_var+corniv(nivdir(kkd,ief)))
                end do
!
! covariables
                do ico=1,nbco
                  if(vecsol(ntnifix+ico))  v=v-covdir(kkd,ico)*x(mod_nb_var+corniv(ntnifix+ico))
                end do

                vpf=vpf+v*v*p_dpa%cd(mod_ic,kkd)/var
                det=det+dlog(var)
              end if !presentc
            end do ! kd

            vmere=vmere+p_spt%probg(current_chr,ig)*dexp(-0.5d0*vpf)
          end do !ig


          if (vmere == 0) then
                 f=INIFINY_REAL_VALUE
          else
                 f=-2.d0*dlog(vmere)+det
          end if

      return
      end subroutine funct_1qtl_modlin_ldla_family
!!***


!!****f* m_qtlmap_analyse_modlin_ldla/test_lin_ldla
!!  NAME
!!    test_lin_ldla
!!  DESCRIPTION
!!    Test des differents effets de nuisance du modele par une LRT compere a une chi2
!! NOTES
!!
!!  SOURCE
      subroutine test_lin_ldla (dataset,spt,shp,chr,ic,est_moy,supnbnivest,fmax,nposx,par1,option_var,option_anal)
      use m_qtlmap_analyse_lin_gen, only : prepinc,contingence,nbef,nbco,meff,mcov,mint,nbnivest
      type(QTLMAP_DATASET)              ,intent(in)         :: dataset
      type(PDD_BUILD)                   ,intent(in)         :: spt
      type(HAPLOTYPE_POSITION_BUILD)   ,intent(in)         :: shp
      logical , intent(in)                        :: est_moy
      integer                        , intent(in) :: chr,ic
      integer                        , intent(in) :: supnbnivest
      real (kind=dp)                 , intent(in) :: fmax
      integer                        , intent(in) :: nposx
      real (kind=dp),dimension(3*dataset%genea%np+2*dataset%genea%nm), intent(in) :: par1
      character(len=4)         , intent(in)               :: option_var,option_anal

!
! Divers
      logical itest
      integer , dimension((3*dataset%genea%np+2*dataset%genea%nm)+2) :: iw
      real (kind=dp), dimension(3*dataset%genea%np+2*dataset%genea%nm) :: par,borni,borns
      real (kind=dp) , allocatable , dimension(:)::w
      real (kind=dp) :: user(1),prob,xlrt_t,f1
      integer :: nbint,iecd,iecq,ief,ibound,npar,ip,i,ifail,liw,liwx,iuser(1)
      integer :: lw,lwx,nbreduit

      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa
      type(DATAMODEL_BASE) , pointer :: dpm


      dpm => dataset%phenoModel
      dg => dataset%genea
      dpa => dataset%phenoAnimal

      lwx = (12*(3*dg%np+2*dg%nm))+((3*dg%np+2*dg%nm)*((3*dg%np+2*dg%nm)-1)/2)
      liwx =  (3*dg%np+2*dg%nm) +2

      allocate (w(lwx))

      write(nficout,600)
  600 format(//,80('*')/'testing model effects',//)

      variance_homo = (option_var == 'homo')
      is_la = (option_anal == 'LA  ')
      is_ld = (option_anal == 'LD  ')
      is_ldla = (option_anal == 'LDLA')
      is_ldjh = (option_anal == 'LDJH')

!
!  preparation de la matrice d'incidence
!
      call prepinc(dataset,spt,chr,nposx,ic,option_anal)
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
        liw=liwx
        lw=lwx

        call minimizing_funct(dataset,npar,ibound,funct_1qtl_modlin_ldla,borni,borns,par,f1,iuser,user,ifail)
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
  !      prob=MATH_QTLMAP_G01ECF('U',xlrt_t,dble(nlev(modele(ic,3+ief))),ifail)
         prob=MATH_QTLMAP_G01ECF('U',xlrt_t,dble(nbreduit),ifail)

          write(nficout,602) trim(dpm%namefix(dpm%modele(ic,3+ief))),nbreduit,xlrt_t,prob
!     &       nlev(modele(ic,3+ief)),xlrt,prob
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

      end if

      end do

      write(nficout,*) "When this probability exceeds the standard threshold corresponding to the 5, 1 or 0.1 Pent level",&
               ", you might consider removing this effect from the model"

      deallocate (w)
      end subroutine test_lin_ldla
!!***


!!****f* m_qtlmap_analyse_modlin_ldla/set_solution_hypothesis0
!!  NAME
!!    set_solution_hypothesis0
!!  DESCRIPTION
!!
!! NOTES
!!
!!  SOURCE

     subroutine set_solution_hypothesis0(dataset,ic,is_hetero,incsol)
       type(QTLMAP_DATASET)       ,intent(in) :: dataset
       integer                            ,intent(in)       :: ic
       type(TYPE_INCIDENCE_SOLUTION)      ,intent(inout)    :: incsol
       logical                            ,intent(in)       :: is_hetero

       integer :: nteff,maxNbPar,jm,ifem,ip,ieff
       integer :: indest,isol,ipar,nlevel,ief,i,ife,start

       real(kind=dp) ,dimension(size(par0)) :: par0_t
       type(GENEALOGY_BASE) , pointer :: dg
       type(PHENOTYPE_BASE) , pointer :: dpa
       type(DATAMODEL_BASE) , pointer :: dpm

       dg => dataset%genea
       dpa => dataset%phenoAnimal
       dpm => dataset%phenoModel

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

       if ( is_hetero ) then
        do ip=1,dg%np
            incsol%sig(1,ip) = par0_t(ip)
        end do
        start = dg%np
       else
        incsol%sig(1,:dg%np) = par0_t(1)
        start = 1
       end if

       ieff=1
       incsol%groupeName(ieff) = 'General Mean'
       incsol%nbParameterGroup(ieff)=1
       incsol%parameterName(ieff,1)   ='General Mean'
       incsol%paramaterValue(ieff,1)  = par0_t(start+1)+dpm%xmut(ic)
       incsol%parameterVecsol(ieff,1) = vecsol0(1)
       incsol%parameterPrecis(ieff,1) = precis0(1)


       ieff=ieff+1
       incsol%groupeName(ieff) = 'Sire polygenic effects'
       incsol%nbParameterGroup(ieff)=dg%np

       indest = start+1
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


!!****f* m_qtlmap_analyse_modlin_ldla/set_solution_hypothesis1
!!  NAME
!!    set_solution_hypothesis1
!!  DESCRIPTION
!!
!! NOTES
!!
!!  SOURCE
     subroutine set_solution_hypothesis1(dataset,shp,ic,is_hetero,qtl,hsire,hdam,incsol)
       type(QTLMAP_DATASET)               ,intent(in)       :: dataset
       type(HAPLOTYPE_POSITION_BUILD)     ,intent(inout)    :: shp
       integer                            ,intent(in)       :: ic
       type(TYPE_INCIDENCE_SOLUTION)      ,intent(inout)    :: incsol
       logical                            ,intent(in)       :: is_hetero,qtl,hsire,hdam

       integer :: nteff,maxNbPar,jm,ifem,ip,ieff
       integer :: ntlev,nbtp,jef,lp,indest,km
       integer :: isol,ipar,nlevel,ief,i,ife,start,i_haplo
       integer :: jhr,j

       real(kind=dp) ,dimension(size(par1)) :: par1_t
       character(LEN=400)  nhr(2)

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
       if(hsire) allocate (incsol%unknown_dam_sig(1,dg%np))

       incsol%hypothesis=1
       !  Mean, Polygenic family
       nteff = 2
       if ( count(dpa%estime(ic,:))  > 0 ) nteff = nteff +1

       ! Haplotype sire
       if ( hsire ) then
        nteff = nteff + 1
       end if

       ! QTLEffect
       if ( qtl ) then
        nteff = nteff + 1
        if ( count(dpa%estime(ic,:))  > 0 ) nteff = nteff + 1
       end if

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
       nbtp = 3 + dpm%modele(ic,1)+dpm%modele(ic,2)+dpm%modele(ic,3)
       do jef=1,dpm%modele(ic,3)
		      ntlev=ntlev*dpm%nlev(dpm%modele(ic,nbtp+jef))
	   end do

	   !max nombre de niveau pour un effet fixe en interaction avec le qtl ?
	   maxNbPar = max(maxNbPar,ntlev*dg%np)
	   maxNbPar = max(maxNbPar,ntlev*count(dpa%estime(ic,:)))
	   maxNbPar = max(maxNbPar,shp%nb_haplo_reduit)

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
      
       if ( is_hetero ) then
        do ip=1,dg%np
            incsol%sig(1,ip) = par1_t(ip)
 !           if(hsire) incsol%unknown_dam_sig(1,ip) = par1_t(np+ip)
        end do
 !       start = 2*np
         start=dg%np
       else
          incsol%sig(1,:dg%np) = par1_t(1)
  !        if(hsire) incsol%unknown_dam_sig(1,:np) = par1_t(2)
  !        start = 2
           start=1
       end if

       ieff=1
       incsol%groupeName(ieff) = 'General Mean'
       incsol%nbParameterGroup(ieff)=1
       incsol%parameterName(ieff,1)   ='General Mean'
       incsol%paramaterValue(ieff,1)  = par1_t(start+1)+dpm%xmut(ic)
       incsol%parameterVecsol(ieff,1) = vecsol1(1)
       incsol%parameterPrecis(ieff,1) = precis1(1)

       indest = start+1
       isol=1

       if ( hsire ) then
         ieff=ieff+1
         incsol%qtl_groupeName(1,1)=ieff ! a modifier..
         incsol%groupeName(ieff) = 'Haplotypes effects'
         incsol%nbParameterGroup(ieff)=shp%nb_haplo_reduit

         do i_haplo=1,shp%nb_haplo_reduit
          isol=isol+1
           if(vecsol1(isol)) then
            indest=indest+1
            incsol%paramaterValue(ieff,i_haplo)    = par1_t(indest)
           else
              incsol%paramaterValue(ieff,i_haplo)  = 0.d0
           end if
          incsol%parameterName(ieff,i_haplo)=trim(shp%name_haplo_reduit(i_haplo))//"(race="//&
           trim(shp%race_haplo_reduit(i_haplo)) &
          //", freq="//trim(str(shp%pb_haplo_reduit(i_haplo)))//")"
          incsol%parameterVecsol(ieff,i_haplo) = vecsol1(isol)
          incsol%parameterPrecis(ieff,i_haplo) = precis1(isol)

          allocate (incsol%haplotypes(dg%np,incsol%hypothesis,2))
          allocate (incsol%races_haplotypes(dg%np,incsol%hypothesis,2))
          do ip=1,dg%np
            do jhr=1,2
              if(shp%num_haplo_pere(ip,jhr,1) /= 0) then
                incsol%haplotypes(ip,1,jhr)=shp%name_haplo_reduit(shp%num_haplo_pere(ip,jhr,1))
                incsol%races_haplotypes(ip,1,jhr)=shp%race_haplo_reduit(shp%num_haplo_pere(ip,jhr,1))
              end if
            end do
          end do
        end do
       end if

       if ( qtl ) then
       ieff=ieff+1
       incsol%qtl_groupeName(1,1)=ieff
       incsol%groupeName(ieff) = 'Sire QTL effects'
       incsol%nbParameterGroup(ieff)=dg%np*ntlevp

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

           nhr= ' '
           do jhr=1,2
            if(shp%num_haplo_pere(ip,jhr,1) /= 0) nhr(jhr) = shp%name_haplo_reduit(shp%num_haplo_pere(ip,jhr,1))
          end do


 !          incsol%parameterName(ieff,ife)   ='Sire '//trim(pere(ip))//" "//trim(str(lp))//" - "&
 !                     //'['//trim(name_haplo_reduit(num_haplo_pere(ip,1,1)))//'/'&
 !                     //trim(name_haplo_reduit(num_haplo_pere(ip,2,1)))//']'

         incsol%parameterName(ieff,ife)   ='Sire '//trim(dg%pere(ip))//" "//trim(str(lp))//" - "&
                      //'['//trim(nhr(1))//'/'//trim(nhr(2))//']'


 !           print*,'ip,num_haplo_pere(ip,2,1),name_haplo_reduit(num_haplo_pere(ip,2,1))',&
 !                   ip,num_haplo_pere(ip,2,1),name_haplo_reduit(num_haplo_pere(ip,2,1))

!
!  a voir avec ofi en fait il y a maintenant > 1 phases possibles
!
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

end module m_qtlmap_analyse_modlin_ldla
