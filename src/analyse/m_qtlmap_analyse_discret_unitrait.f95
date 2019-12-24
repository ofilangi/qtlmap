!!****m* ANALYSE/m_qtlmap_analyse_discret_unitrait
!!  NAME
!!    m_qtlmap_analyse_discret_unitrait
!!  SYNOPSIS
!!
!!  DESCRIPTION
!!   Analysis module for LA analysis with discrete data and model description
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
module m_qtlmap_analyse_discret_unitrait
    use m_qtlmap_types
    use m_qtlmap_log
    use m_qtlmap_base
    use m_qtlmap_math
    use m_qtlmap_optimization
    use m_qtlmap_analyse_gen
    use m_qtlmap_analyse_lin_gen
    use m_qtlmap_output_handler

    implicit none
    save
    ! dim : np
      real (kind=dp)       ,dimension(:),allocatable,private   :: sig1
      real (kind=dp)       ,dimension(:),allocatable,private   :: xmu1p
      real (kind=dp)       ,dimension(:),allocatable,private   :: xmu1m
      real (kind=dp)                                ,private   :: xmu1g
      real (kind=dp)                                ,public    :: f0

      real (kind=dp)       ,dimension(:),allocatable,private   :: fp0
      real (kind=dp)       ,dimension(:),allocatable,private   :: fp1
      real (kind=dp)       ,dimension(:),allocatable,private   :: fm0
      real (kind=dp)       ,dimension(:),allocatable,private   :: fm1

    integer              ,private                            :: current_chr
    integer              ,private                            :: current_ic
    real (kind=dp)       ,private, parameter                 :: GRAND=1.d6

    type(QTLMAP_DATASET) , pointer , private :: dataset_p => null()
    type(PDD_BUILD)      , pointer , private :: spt_p => null()


    public :: init_analyse_discret_unitrait
    public :: opti_0qtl_discret_unitrait
    public :: opti_1qtl_discret_unitrait
    public :: test_lin
    public :: end_analyse_discret_unitrait
    public :: set_solution_hypothesis0
    public :: set_solution_hypothesis1


    contains

    subroutine init_analyse_discret_unitrait(dataset,spt)
         type(QTLMAP_DATASET)   ,target    ,intent(in)         :: dataset
         type(PDD_BUILD)        ,target    ,intent(in)         :: spt

         integer           :: stat
         type(GENEALOGY_BASE) , pointer :: dg

         dg => dataset%genea

         allocate (sig1(dg%np),STAT=stat)
         call check_allocate(stat,'sig1 [m_qtlmap_analyse_discret_unitrait]')
         allocate (xmu1p(dg%np),STAT=stat)
         call check_allocate(stat,'xmu1p [m_qtlmap_analyse_discret_unitrait]')
         allocate (xmu1m(dg%nm),STAT=stat)
         call check_allocate(stat,'xmu1m [m_qtlmap_analyse_discret_unitrait]')
         allocate (fp0(dg%np),STAT=stat)
         call check_allocate(stat,'fp0 [m_qtlmap_analyse_discret_unitrait]')
         allocate (fp1(dg%np),STAT=stat)
         call check_allocate(stat,'fp1 [m_qtlmap_analyse_discret_unitrait]')
         allocate (fm0(dg%nm),STAT=stat)
         call check_allocate(stat,'fm0 [m_qtlmap_analyse_discret_unitrait]')
         allocate (fm1(dg%nm),STAT=stat)
         call check_allocate(stat,'fm1 [m_qtlmap_analyse_discret_unitrait]')

         dataset_p => dataset
         spt_p => spt

     end subroutine init_analyse_discret_unitrait

     subroutine end_analyse_discret_unitrait
         deallocate (sig1)
         deallocate (xmu1p)
         deallocate (xmu1m)
         deallocate (fp0)
         deallocate (fp1)
         deallocate (fm0)
         deallocate (fm1)

         dataset_p => null()
         spt_p => null()

     end subroutine end_analyse_discret_unitrait


!!****f* m_qtlmap_analyse_discret_unitrait/opti_0qtl_discret_unitrait
!! NAME
!!    opti_0qtl_discret_unitrait
!! DESCRIPTION
!!   Compute the likelihood 0 QTL, 1 trait , fixed effect and covaruiate included
!! INPUTS
!!    ic    : index trait
!! NOTES
!!    contingence,precision,set_filter_optim,minimizing_funct_family,funct_0qtl_discret_family
!! SOURCE
      subroutine opti_0qtl_discret_unitrait(dataset,spt,ic)
      use m_qtlmap_analyse_lin_gen, only : contingence,precision,vecsol,precis, &
                                    nbnivest,ntniv,par0,precis0,vecsol0,xx

       type(QTLMAP_DATASET)       ,intent(in)         :: dataset
       type(PDD_BUILD)            ,intent(in)         :: spt

      integer , intent(in)   :: ic
!
! Divers

      integer                                    :: iuser(1)
      real (kind=dp) ,dimension(:),allocatable   :: par,borni,borns
      real (kind=dp)                             :: user(1)
      real (kind=dp)                             :: f,eff,cumul

      integer :: npar,ibound,ideb,ip,ix,ifail,liw,lw,i,indest,ntot
      integer :: km,jm,liwx,lwx, j, m, m1, temp, k, kkd, ii
      logical itest
      logical found
      logical  , dimension(:,:,:),pointer        :: filter_inc
      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa
      type(DATAMODEL_BASE) , pointer :: dpm

      dpm => dataset%phenoModel
      dg => dataset%genea
      dpa => dataset%phenoAnimal

!******************************************************************************
!
!  recherche des elements estimables a partir de la decomposition de choleski
!  comptage et recodification de ces elements (on commence par mettre a 0 les
!  compteurs destines aux tests des effes fixes
!

      itest=.false.
      call contingence(dataset,spt,ic,0,itest,.false.)
      call precision(xx,precis)

      current_ic=ic
      npar=dg%np+nbnivest+dpm%nmod(ic)-1
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
      do ix=dg%np+1,npar-dpm%nmod(ic)+2
        borni(ix)=XMU_MIN
        borns(ix)=XMU_MAX
      end do
      do ix=npar-dpm%nmod(ic)+3,npar
        borni(ix)=SIG_MIN
        borns(ix)=XMU_MAX
      end do
!
! Point de depart
      do ip=1,dg%np
        par(ip)=sig0(ip)
       end do
      do ix=dg%np+1,dg%np+nbnivest
        par(ix)=0.d0
      end do
      par(dg%np+nbnivest+1)=dpm%seuil(ic,1)
      do ix=2,dpm%nmod(ic)-1
        par(dg%np+nbnivest+ix)=dpm%seuil(ic,ix)-dpm%seuil(ic,ix-1)
      end do
!
!
! Optimisation de la vraisemblance
      ifail=1
   !   call minimizing_funct(npar,ibound,funct_0qtl_discret_unitrait,borni,borns,par,f,iuser,user,ifail)
    call minimizing_funct_family(dataset,npar,ibound,funct_0qtl_discret_family,&
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

      end subroutine opti_0qtl_discret_unitrait
!!***


!!****f* m_qtlmap_analyse_discret_unitrait/funct_0qtl_discret_unitrait
!! NAME
!!    funct_0qtl_discret_unitrait
!! DESCRIPTION
!!   likelihood function 0 QTL, 1 trait
!! INPUTS
!!    n     : number of parameter of the vector x
!!    x     : vector inputs
!!   iuser  : user parameters of integer type
!!   user   : user parameters of real type
!! OUTPUTS
!!    f     : result of the function
!! NOTES
!!    MATH_QTLMAP_G01EAF
!! SOURCE
      subroutine funct_0qtl_discret_unitrait(n,x,f,iuser,user)
      use m_qtlmap_analyse_lin_gen, only : nbnivest, nbfem,nbef,ntnifix,nbco,covdir,nivdir,corniv,vecsol

      implicit none
      integer         , intent(in)                       :: n
!
! Tableaux dimensionnes selon n, le nombre de parametres a estimer
!
      real (kind=dp)   ,dimension(n), intent(in)         :: x
      real (kind=dp)   ,intent(inout)                    :: f
      integer          ,dimension(1), intent(in)      :: iuser
      real (kind=dp)   ,dimension(1), intent(in)         :: user

! Tableaux dimensionnes selon nm, le nombre de meres
      real (kind=dp)   ,dimension(dataset_p%phenoModel%nmod(current_ic)) :: threshold

!modif mkw
! declaration de l comme entier et on a plus besoin de la variable vpf

      integer :: jm, i, kkd, ifail
      integer ::jnit,ip,neffet,ief,ico,ic
      real (kind=dp) :: v, vpf,tig
!mkw
!
   neffet=dataset_p%genea%np+nbnivest

!******************************************************************************
!
    ifail=1
    ic = current_ic

    threshold(1)=x(neffet+1)
    do i=2,dataset_p%phenoModel%nmod(ic)-1
      threshold(i)=threshold(i-1)+x(neffet+i)
    end do
!

      jnit=2
      if (nbfem.eq.0)jnit=1

      f=0.d0
      do ip=1,dataset_p%genea%np
        tig=grand
        if (x(ip)> 0)tig=1.d0/x(ip)
        fp0(ip)=0.d0

        do jm=dataset_p%genea%nmp(ip)+1,dataset_p%genea%nmp(ip+1)
        fm0(jm)=0.d0
!
! on ne considere que les meres

          do kkd=dataset_p%genea%ndm(jm)+1,dataset_p%genea%ndm(jm+1)
!
            if(dataset_p%phenoAnimal%presentc(ic,kkd)) then

             v=0.d0
             if(vecsol(nivdir(kkd,1))) v=v+x(dataset_p%genea%np+corniv(nivdir(kkd,1)))

             if(dataset_p%phenoAnimal%estime(ic,jm)) then
               if(vecsol(nivdir(kkd,2)))  v=v+x(dataset_p%genea%np+corniv(nivdir(kkd,2)))
             end if

             do ief=jnit+1,jnit+nbef
               if(vecsol(nivdir(kkd,ief))) v=v+x(dataset_p%genea%np+corniv(nivdir(kkd,ief)))
             end do

             do ico=1,nbco
               if(vecsol(ntnifix+ico)) v=v+covdir(kkd,ico)*x(dataset_p%genea%np+corniv(ntnifix+ico))
             end do

             if (dataset_p%phenoAnimal%ydiscretord(ic,kkd)==1)  then
	      vpf=MATH_QTLMAP_G01EAF('L',((threshold(1)-v)*tig),ifail)
             else
              if (dataset_p%phenoAnimal%ydiscretord(ic,kkd)==dataset_p%phenoModel%nmod(ic)) then
        	vpf=1-MATH_QTLMAP_G01EAF('L',(threshold(dataset_p%phenoModel%nmod(ic)-1)-v)*tig,ifail)
              else
        	 vpf=MATH_QTLMAP_G01EAF('L',(threshold(dataset_p%phenoAnimal%ydiscretord(ic,kkd))-v)*tig,ifail)-&
        	     MATH_QTLMAP_G01EAF('L',(threshold(dataset_p%phenoAnimal%ydiscretord(ic,kkd)-1)-v)*tig,ifail)
              end if
	     end if

	      if (vpf <= 0) then
                fm0(jm)=INIFINY_REAL_VALUE
              else
                fm0(jm)=fm0(jm)-dlog(vpf)
	      end if

	    end if	! si presentc
	  end do	! sur kkd

	  fp0(ip)=fp0(ip)+fm0(jm)
	  f=f+fm0(jm)

        end do ! sur jm
      end do ! sur ip

      end subroutine funct_0qtl_discret_unitrait
!!***

!!****f* m_qtlmap_analyse_discret_unitrait/funct_0qtl_discret_family
!! NAME
!!    funct_0qtl_discret_family
!! DESCRIPTION
!!   likelihood function 0 QTL, 1 trait for one family
!! INPUTS
!!    ip    : index sire
!!    jm    : idex dam
!!    n     : number of parameter of the vector x
!!    x     : vector inputs
!!   iuser  : user parameters of integer type
!!   user   : user parameters of real type
!! OUTPUTS
!!    f     : result of the function
!! NOTES
!!    MATH_QTLMAP_G01EAF
!! SOURCE
      subroutine funct_0qtl_discret_family(ip,jm,n,x,f,iuser,user)
      integer         , intent(in)                       :: ip,jm,n
!
! Tableaux dimensionnes selon n, le nombre de parametres a estimer
!
      real (kind=dp)   ,dimension(n), intent(in)         :: x
      real (kind=dp)   ,intent(inout)                    :: f
      integer          ,dimension(1), intent(in)         :: iuser
      real (kind=dp)   ,dimension(1), intent(in)         :: user

! Tableaux dimensionnes selon nm, le nombre de meres
      real (kind=dp)   ,dimension(dataset_p%phenoModel%nmod(current_ic))      :: threshold

!modif mkw
! declaration de l comme entier et on a plus besoin de la variable vpf

      integer :: i, kkd, ifail
      integer ::jnit,neffet,ief,ico,ic
      real (kind=dp) :: v, vpf,tig
!mkw
!
      neffet=dataset_p%genea%np+nbnivest

!******************************************************************************
!
    ifail=1
    ic = current_ic

    threshold(1)=x(neffet+1)
    do i=2,dataset_p%phenoModel%nmod(ic)-1
      threshold(i)=threshold(i-1)+x(neffet+i)
    end do
!

    jnit=2
    if (nbfem.eq.0)jnit=1

    f=0.d0
    tig=grand
    if (x(ip)> 0)tig=1.d0/x(ip)
!
! on ne consid�e que les m�res

    do kkd=dataset_p%genea%ndm(jm)+1,dataset_p%genea%ndm(jm+1)
       if(dataset_p%phenoAnimal%presentc(ic,kkd)) then
             v=0.d0
             if(vecsol(nivdir(kkd,1))) v=v+x(dataset_p%genea%np+corniv(nivdir(kkd,1)))

             if(dataset_p%phenoAnimal%estime(ic,jm)) then
               if(vecsol(nivdir(kkd,2)))  v=v+x(dataset_p%genea%np+corniv(nivdir(kkd,2)))
             end if

             do ief=jnit+1,jnit+nbef
               if(vecsol(nivdir(kkd,ief))) v=v+x(dataset_p%genea%np+corniv(nivdir(kkd,ief)))
             end do

             do ico=1,nbco
               if(vecsol(ntnifix+ico)) v=v+covdir(kkd,ico)*x(dataset_p%genea%np+corniv(ntnifix+ico))
             end do

             if (dataset_p%phenoAnimal%ydiscretord(ic,kkd)==1)  then
	          vpf=MATH_QTLMAP_G01EAF('L',((threshold(1)-v)*tig),ifail)
             else
              if (dataset_p%phenoAnimal%ydiscretord(ic,kkd)==dataset_p%phenoModel%nmod(ic)) then
         	  vpf=1-MATH_QTLMAP_G01EAF('L',(threshold(dataset_p%phenoModel%nmod(ic)-1)-v)*tig,ifail)
              else
        	  vpf=MATH_QTLMAP_G01EAF('L',(threshold(dataset_p%phenoAnimal%ydiscretord(ic,kkd))-v)*tig,ifail)-&
        	     MATH_QTLMAP_G01EAF('L',(threshold(dataset_p%phenoAnimal%ydiscretord(ic,kkd)-1)-v)*tig,ifail)
              end if
	     end if

	      if (vpf <= 0) then
                f=INIFINY_REAL_VALUE
                return
              else
                f=f-dlog(vpf)
	      end if

	    end if	! si presentc
	  end do	! sur kkd

   end subroutine funct_0qtl_discret_family
!!***

!mkw

!!****f* m_qtlmap_analyse_discret_unitrait/opti_1qtl_discret_unitrait
!! NAME
!!    opti_1qtl_discret_unitrait
!! DESCRIPTION
!!   Compute the test statistic among the chromosome : 1 QTL, 1 trait
!! INPUTS
!!    ic        : index trait
!! OUTPUTS
!!    lrtsol          : likelihood ratio test information (see TYPE_LRT_SOLUTION)
!!    fmax            : value of the likelihood function at the maximum
!!    supnbnivest     : number of paramter estimated in the solution (see contingence)
!! NOTES
!!    prepinc,contingence,precision,set_filter_optim,minimizing_funct_family,funct_1qtl_discret_family
!! SOURCE
      subroutine opti_1qtl_discret_unitrait(dataset,spt,ic,lrtsol,fmax,supnbnivest)
      type(QTLMAP_DATASET)       ,intent(in)            :: dataset
      type(PDD_BUILD)            ,intent(in)            :: spt

      integer , intent(in)                                :: ic
      type(TYPE_LRT_SOLUTION)  , intent(out)              :: lrtsol
      real (kind=dp)  , intent(out)                       :: fmax
      integer         , intent(out)                       :: supnbnivest
!

! Divers
      integer                                    :: iuser(1)
      real (kind=dp) ,dimension(:),allocatable :: par,borni,borns
      real (kind=dp),dimension(:),allocatable       :: val
      real (kind=dp)                             :: user(1)
      real (kind=dp)                             :: f,eff,cumul
      real (kind=dp)                             :: xlrt_t,f1
      logical                                    :: found

      logical itest,stvecsol(size(vecsol1))
      integer :: ip,i,ideb,n,ilong,ibound,npar,ix,j,ifail,liw
      integer :: ii,m, m1, k,temp,chr
      logical  , dimension(:,:,:),pointer        :: filter_inc

      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa
      type(DATAMODEL_BASE) , pointer :: dpm

      dpm => dataset%phenoModel
      dg => dataset%genea
      dpa => dataset%phenoAnimal


      call lrtsol%new(dataset,1)

      current_ic  = ic

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
!
!******************************************************************************
!******************************************************************************
!
! transformation des donn�es discret � des donn�e utilisable par QTLMAP
!
!
!*****************************************************************************
!******************************************************************************
!
!
!
!     integer        ,dimension(:),allocatable   :: ydiscretord, indicemod
!     integer ::i, j, m, m1, temp, k, nmod, ii
! Point de depart

      allocate ( val( dg%np+ntnivmax +dpm%nmod(ic)-1 ) )
      allocate ( par( dg%np+ntnivmax +dpm%nmod(ic)-1) )
      allocate ( borni( dg%np+ntnivmax +dpm%nmod(ic)-1) )
      allocate ( borns( dg%np+ntnivmax +dpm%nmod(ic)-1) )

      !MODIF - OPMIZATION OFI
      allocate (filter_inc(dg%np,dg%nm, dg%np+ntnivmax +dpm%nmod(ic)-1))
      !FIN MODIF - OPMIZATION OFI

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


      lrtsol%lrtmax(0)=-1.d75

      do chr=1,dataset%map%nchr
      n=0
      current_chr = chr
      ilong=dataset%map%get_ilong(chr)
      do ix=0,ilong,dataset%map%pas
        n=n+1
!
!  on stocke les conditions en n
        do i=1,ntniv
          stvecsol(i)=vecsol(i)
        end do

        call prepinc(dataset,spt,1,n,ic,"LA  ")
	    itest=.false.
        call contingence(dataset,spt,ic,1,itest,.false.)
	    npar=dg%np+nbnivest+dpm%nmod(ic)-1

	    !MODIF - OPMIZATION
        call set_filter_optim(dataset,ic,.true.,.false.,ntnivmax,ntniv,vecsol,xinc,filter_inc)
        !FIN MODIF - OPMIZATION


        ibound=0
        do ip=1,dg%np
         borni(ip)=SIG_MIN
         borns(ip)=SIG_MAX
        end do
        do i=dg%np+1,npar-dpm%nmod(ic)+2
         borni(i)=XMU_MIN
         borns(i)=XMU_MAX
        end do
        do i=npar-dpm%nmod(ic)+3,npar
         borni(i)=SIG_MIN
         borns(i)=XMU_MAX
        end do
!
! Point de depart
        do ip=1,dg%np
         par(ip)=sig1(ip)
        end do
        do i=dg%np+1,dg%np+nbnivest
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

      par(dg%np+nbnivest+1)=dpm%seuil(ic,1)
      do i=2,dpm%nmod(ic)-1
        par(dg%np+nbnivest+i)=dpm%seuil(ic,i)-dpm%seuil(ic,i-1)
      end do


!  on stocke les conditions en n
	!call prepinc(n,ic)
	!itest=.false.
        !call contingence(ic,1,itest)

! Optimisation de la vraisemblance a la position dx

        ifail=1

        !call minimizing_funct(npar,ibound,funct_1qtl_discret_unitrait,borni,borns,par,f1,iuser,user,ifail)
        call minimizing_funct_family(dataset,npar,ibound,funct_1qtl_discret_family,&
         filter_inc,fm1,fp1,borni,borns,par,f1,iuser,user,ifail)

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
!  preparation de la matrice d'incidence
!
!
!  recherche des elements estimables � partir de la d�composition de choleski
!  comptage et recodification de ces elements (on commence par mettre � 0 les
!  compteurs destin�s aux tests des effes fix�s

     if ( f1 < INIFINY_REAL_VALUE ) then
         xlrt_t=-2.d0*(f1-f0)
      else
         xlrt_t=0
      end if


!  on met les profil / prog�niteur dnas les ficheir ad hoc

        do ii=1,dg%np
          call lrtsol%LRT_SIRES(ii)%add1p(dataset,chr,n,(-2.d0*(fp1(ii)-fp0(ii))))
        end do
        do ii=1,dg%nm
          call lrtsol%LRT_DAMS(ii)%add1p(dataset,chr,n,(-2.d0*(fm1(ii)-fm0(ii))))
        end do

!  on stocke la position si elle est meilleure que les pr�c�dents
!
        if(lrtsol%lrtmax(0) < xlrt_t) then

          lrtsol%chrmax(0)=chr
          lrtsol%nxmax(0)=n
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
        call lrtsol%LRT%add1p(dataset,chr,n,xlrt_t)
      end do ! ix
     end do ! chr
!
!  on calcule la pr�cision des estimation au point correspondant
! au LRT maximum
!
      call prepinc(dataset,spt,lrtsol%chrmax(0),lrtsol%nxmax(0),ic,"LA  ")
      call contingence(dataset,spt,ic,1,itest,.false.)
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

      end subroutine opti_1qtl_discret_unitrait
!!***

!!****f* m_qtlmap_analyse_discret_unitrait/funct_1qtl_discret_unitrait
!! NAME
!!    funct_1qtl_discret_unitrait
!! DESCRIPTION
!!   likelihood function 1 QTL, 1 trait
!! INPUTS
!!    n     : number of parameter of the vector x
!!    x     : vector inputs
!!   iuser  : user parameters of integer type
!!   user   : user parameters of real type
!! OUTPUTS
!!    f     : result of the function
!! NOTES
!!    structure used : vecsol,nivdir,corniv,covdir,presentc,ydiscretord,nmod,ngenom,ndesc,ngend
!!    function       : estime,MATH_QTLMAP_G01EAF
!! SOURCE
      subroutine funct_1qtl_discret_unitrait(n,x,f,iuser,user)
      use m_qtlmap_analyse_lin_gen, only : nbnivest

      implicit none

      integer         , intent(in)                  :: n
!
! Tableaux dimensionnes selon n, le nombre de parametres a estimer
!
      real (kind=dp)      ,dimension(n), intent(in) :: x
      real (kind=dp)  , intent(inout)               :: f
      integer ,       dimension(1), intent(in)      :: iuser
      real (kind=dp)      ,dimension(1), intent(in) :: user
      real (kind=dp)   ,dimension(dataset_p%phenoModel%nmod(current_ic))      :: threshold

! Tableaux dimensionnes selon nm, le nombre de meres
      real (kind=dp)   ,dimension(dataset_p%genea%nm)   :: effm


!
! Divers
      integer :: init,jnit,ip,nm1,nm2,jm,ngeno1,ngeno2,ig,chr,ic
      integer :: nd1,nd2,neffet,kd,kkd,ilev,ief,ico,lambda
      integer :: ntnivmax,neffmax, m, i, j, m1, temp, k, ii, ifail


      real (kind=dp) :: sig,var,vmere,vpf,v,wpf,tig

!******************************************************************************
!
      chr =current_chr
      ic = current_ic
      neffet=dataset_p%genea%np+nbnivest
      ifail=1

      threshold(1)=x(neffet+1)
      do i=2,dataset_p%phenoModel%nmod(ic)-1
        threshold(i)=threshold(i-1)+x(neffet+i)
      end do

     init=3
     jnit=4
      if (nbfem.eq.0)then
       init=2
       jnit=2
      end if

      f=0.d0
      do ip=1,dataset_p%genea%np
      tig=grand
      if (x(ip)> 0)tig=1.d0/x(ip)
        fp1(ip)=0.d0
        do jm=dataset_p%genea%nmp(ip)+1,dataset_p%genea%nmp(ip+1)
	vmere=0.d0
	fm1(jm)=0.d0
          do ig=spt_p%ngenom(chr,jm)+1,spt_p%ngenom(chr,jm+1)
	    vpf=grand
            do kd=spt_p%ngend(chr,ig)+1,spt_p%ngend(chr,ig+1)
              kkd=spt_p%ndesc(chr,kd)

              if(dataset_p%phenoAnimal%presentc(ic,kkd)) then
             	v=0.d0

             	if(vecsol(nivdir(kkd,1))) v=v+x(dataset_p%genea%np+corniv(nivdir(kkd,1)))*ppt(kd)

             	if(vecsol(nivdir(kkd,init))) v=v+x(dataset_p%genea%np+corniv(nivdir(kkd,init)))

             	if(dataset_p%phenoAnimal%estime(ic,jm)) then
             	  if(vecsol(nivdir(kkd,2))) v=v+x(dataset_p%genea%np+corniv(nivdir(kkd,2)))*pmt(kd)
             	  if(vecsol(nivdir(kkd,4))) v=v+x(dataset_p%genea%np+corniv(nivdir(kkd,4)))
             	end if

             	do ief=jnit+1,jnit+nbef
             	  if(vecsol(nivdir(kkd,ief))) v=v+x(dataset_p%genea%np+corniv(nivdir(kkd,ief)))
                end do

             	do ico=1,nbco
             	  if(vecsol(ntnifix+ico)) v=v+covdir(kkd,ico)*x(dataset_p%genea%np+corniv(ntnifix+ico))
             	end do

    	  	if (dataset_p%phenoAnimal%ydiscretord(ic,kkd)==1) then
    	  	   wpf=MATH_QTLMAP_G01EAF('L',((threshold(1)-v)*tig),ifail)
    	  	else
    	  	   if (dataset_p%phenoAnimal%ydiscretord(ic,kkd)==dataset_p%phenoModel%nmod(ic)) then
    	  	     wpf=1.d0-MATH_QTLMAP_G01EAF('L',(threshold(dataset_p%phenoModel%nmod(ic)-1)-v)*tig,ifail)
    	  	   else
    	  	     wpf=MATH_QTLMAP_G01EAF('L',(threshold(dataset_p%phenoAnimal%ydiscretord(ic,kkd))-v)*tig,ifail)-&
    	  		 MATH_QTLMAP_G01EAF('L',(threshold(dataset_p%phenoAnimal%ydiscretord(ic,kkd)-1)-v)*tig,ifail)
    	  	   end if
    	    	end if

		vpf=vpf*wpf

	      end if ! presentc
            end do ! kd
	     vmere=vmere+spt_p%probg(chr,ig)*vpf
           end do !ig
     	    if (vmere <= 0) then
     		     fm1(jm)=INIFINY_REAL_VALUE
     	      else
     		    fm1(jm)=dlog(grand)-dlog(vmere)
     	     end if
	   fp1(ip)=fp1(ip)+fm1(jm)
	   f=f+fm1(jm)

         end do ! jm
      end do ! ip

      return
      end subroutine funct_1qtl_discret_unitrait
!!***



!!****f* m_qtlmap_analyse_discret_unitrait/funct_1qtl_discret_family
!! NAME
!!    funct_1qtl_discret_family
!! DESCRIPTION
!!   likelihood function 1 QTL, 1 trait for a family
!! INPUTS
!!    ip    : index sire
!!    jm    : idex dam
!!    n     : number of parameter of the vector x
!!    x     : vector inputs
!!   iuser  : user parameters of integer type
!!   user   : user parameters of real type
!! OUTPUTS
!!    f     : result of the function
!! NOTES
!!    structure used : vecsol,nivdir,corniv,covdir,presentc,ydiscretord,nmod,ngenom,ndesc,ngend
!!    function       : estime,MATH_QTLMAP_G01EAF
!! SOURCE
       subroutine funct_1qtl_discret_family(ip,jm,n,x,f,iuser,user)

      integer         , intent(in)                  :: ip,jm,n
!
! Tableaux dimensionnes selon n, le nombre de parametres a estimer
!
      real (kind=dp)      ,dimension(n), intent(in) :: x
      real (kind=dp)  , intent(inout)               :: f
      integer ,       dimension(1), intent(in)      :: iuser
      real (kind=dp)      ,dimension(1), intent(in) :: user
      real (kind=dp)   ,dimension(dataset_p%phenoModel%nmod(current_ic))      :: threshold

! Tableaux dimensionnes selon nm, le nombre de meres
      real (kind=dp)    :: effm


!
! Divers
      integer :: init,jnit,nm1,nm2,ngeno1,ngeno2,ig,chr,ic
      integer :: nd1,nd2,neffet,kd,kkd,ilev,ief,ico,lambda
      integer :: ntnivmax,neffmax, m, i, j, m1, temp, k, ii, ifail


      real (kind=dp) :: sig,var,vmere,vpf,v,wpf,tig

!******************************************************************************
!
      chr =current_chr
      ic = current_ic
      neffet=dataset_p%genea%np+nbnivest
      ifail=1
      threshold(1)=x(neffet+1)
      do i=2,dataset_p%phenoModel%nmod(ic)-1
        threshold(i)=threshold(i-1)+x(neffet+i)
      end do

     init=3
     jnit=4
      if (nbfem.eq.0)then
       init=2
       jnit=2
      end if

      f=0.d0
      tig=grand
      if (x(ip)> 0)tig=1.d0/x(ip)
   	  vmere=0.d0
      do ig=spt_p%ngenom(chr,jm)+1,spt_p%ngenom(chr,jm+1)
	    vpf=grand
            do kd=spt_p%ngend(chr,ig)+1,spt_p%ngend(chr,ig+1)
              kkd=spt_p%ndesc(chr,kd)

              if(dataset_p%phenoAnimal%presentc(ic,kkd)) then
             	v=0.d0

             	if(vecsol(nivdir(kkd,1))) v=v+x(dataset_p%genea%np+corniv(nivdir(kkd,1)))*ppt(kd)

             	if(vecsol(nivdir(kkd,init))) v=v+x(dataset_p%genea%np+corniv(nivdir(kkd,init)))

             	if(dataset_p%phenoAnimal%estime(ic,jm)) then
             	  if(vecsol(nivdir(kkd,2))) v=v+x(dataset_p%genea%np+corniv(nivdir(kkd,2)))*pmt(kd)
             	  if(vecsol(nivdir(kkd,4))) v=v+x(dataset_p%genea%np+corniv(nivdir(kkd,4)))
             	end if

             	do ief=jnit+1,jnit+nbef
             	  if(vecsol(nivdir(kkd,ief))) v=v+x(dataset_p%genea%np+corniv(nivdir(kkd,ief)))
                end do

             	do ico=1,nbco
             	  if(vecsol(ntnifix+ico)) v=v+covdir(kkd,ico)*x(dataset_p%genea%np+corniv(ntnifix+ico))
             	end do

    	  	if (dataset_p%phenoAnimal%ydiscretord(ic,kkd)==1) then
    	  	   wpf=MATH_QTLMAP_G01EAF('L',((threshold(1)-v)*tig),ifail)
    	  	else
    	  	   if (dataset_p%phenoAnimal%ydiscretord(ic,kkd)==dataset_p%phenoModel%nmod(ic)) then
    	  	     wpf=1.d0-MATH_QTLMAP_G01EAF('L',(threshold(dataset_p%phenoModel%nmod(ic)-1)-v)*tig,ifail)
    	  	   else
    	  	     wpf=MATH_QTLMAP_G01EAF('L',(threshold(dataset_p%phenoAnimal%ydiscretord(ic,kkd))-v)*tig,ifail)-&
    	  		 MATH_QTLMAP_G01EAF('L',(threshold(dataset_p%phenoAnimal%ydiscretord(ic,kkd)-1)-v)*tig,ifail)
    	  	   end if
    	    	end if

	   	  vpf=vpf*wpf

	      end if ! presentc
            end do ! kd
	     vmere=vmere+spt_p%probg(chr,ig)*vpf
           end do !ig
     	 if (vmere <= 0) then
             f=INIFINY_REAL_VALUE
         else
     		 f=dlog(grand)-dlog(vmere)
     	 end if

      end subroutine funct_1qtl_discret_family
!!***


!!****f* m_qtlmap_analyse_discret_unitrait/test_lin
!! NAME
!!    test_lin
!! DESCRIPTION
!!   Test des differents effets de nuisance du modele par une LRT compare a une chi2 and write in the output result file.
!! INPUTS
!!    chr          : chromosome index
!!    ic           : index trait
!!   est_moy       : if the general mean is estimated
!!   supnbnivest   : number of paramter of the solution
!!   fmax          : value of the likelihood at the maximum position
!!   nposx         : index of the position (see absi)
!! NOTES
!!***
      subroutine test_lin(dataset,spt,chr,ic,est_moy,supnbnivest,fmax,nposx)
       type(QTLMAP_DATASET)       ,intent(in)            :: dataset
       type(PDD_BUILD)            ,intent(in)         :: spt

      logical , intent(in)                        :: est_moy
      integer                        , intent(in) :: chr,ic
      integer                        , intent(in) :: supnbnivest
      real (kind=dp)                 , intent(in) :: fmax
      integer                        , intent(in) :: nposx

!
! Divers
      logical :: itest
      integer :: iuser(1)
       real (kind=dp) ,dimension(:),allocatable :: par,borni,borns
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

      allocate ( par( dg%np + ntnivmax +dpm%nmod(ic)-1 ) )
      allocate ( borni( dg%np + ntnivmax +dpm%nmod(ic)-1) )
      allocate ( borns( dg%np + ntnivmax +dpm%nmod(ic)-1) )
      !MODIF - OPMIZATION
      allocate (filter_inc(dg%np,dg%nm, dg%np + ntnivmax+dpm%nmod(ic)-1 ))
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
        call contingence(dataset,spt,ic,1,itest,.false.)

!
!  la r�dustion du nombre d'effet estim�e est calcul�e
!
      nbreduit=supnbnivest-nbnivest
!
! Parametres de maximisation
      ibound=0
	npar=dg%np+nbnivest+dpm%nmod(ic)-1

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

        !call minimizing_funct(npar,ibound,funct_1qtl_discret_unitrait,borni,borns,par,f1,iuser,user,ifail)
        call minimizing_funct_family(dataset,npar,ibound,funct_1qtl_discret_family,filter_inc,&
          fm1,fp1,borni,borns,par,f1,iuser,user,ifail)
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
!     &       dpm%nlev(dpm%modele(ic,3+ief)),xlrt,prob
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
      deallocate ( par  )
      deallocate ( borni )
      deallocate ( borns )
      deallocate (filter_inc)

       write(nficout,*) "When this probability exceeds the standard threshold corresponding to the 5, 1 or 0.1 Pent level",&
               ", you might consider removing this effect from the model"

      return
      end subroutine test_lin

!!****f* m_qtlmap_analyse_discret_unitrait/set_solution_hypothesis0
!! NAME
!!    set_solution_hypothesis0
!! DESCRIPTION
!!   fill the solution under H0 in the structure TYPE_INCIDENCE_SOLUTION
!! INPUTS
!!    ic           : index trait
!! OUTPUTS
!!   incsol        : the variable to fill
!! NOTES
!!***
     subroutine set_solution_hypothesis0(dataset,ic,incsol)
       type(QTLMAP_DATASET)       ,intent(in)            :: dataset
       integer                            ,intent(in)    :: ic
       type(TYPE_INCIDENCE_SOLUTION)      ,intent(inout) :: incsol

       integer :: nteff,maxNbPar,jm,ifem,ip,ieff,nparx
       integer :: indest,isol,ipar,nlevel,ief,i,ife
       type(GENEALOGY_BASE) , pointer :: dg
       type(PHENOTYPE_BASE) , pointer :: dpa
       type(DATAMODEL_BASE) , pointer :: dpm

       dpm => dataset%phenoModel
       dpa => dataset%phenoAnimal
       dg => dataset%genea

       allocate (incsol%sig(1,dg%np))

       incsol%hypothesis=0
       !  Polygenic family
       nteff = 1
       if ( count(dpa%estime(ic,:))  > 0 ) nteff = 2

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
       allocate (incsol%nbParameterGroup(nteff))
       allocate (incsol%parameterName(nteff,maxNbPar))
       allocate (incsol%paramaterValue(nteff,maxNbPar))
       allocate (incsol%parameterVecsol(nteff,maxNbPar))
       allocate (incsol%parameterPrecis(nteff,maxNbPar))
       incsol%parameterPrecis=0.d0
       incsol%parameterVecsol=.true.

       do ip=1,dg%np
            incsol%sig(1,ip) = par0(ip)
       end do

       ieff=1
       incsol%groupeName(ieff) = 'Sire polygenic effects'
       incsol%nbParameterGroup(ieff)=dg%np

       indest = dg%np
       isol=1
       do ip=1,dg%np
           isol=isol+1
           if ( vecsol0(isol) ) then
             indest=indest+1
             incsol%paramaterValue(ieff,ip)  = par0(indest)
           else
             incsol%paramaterValue(ieff,ip)  = 0.d0
           end if

           incsol%parameterName(ieff,ip)   = 'Sire '//trim(dg%pere(ip))
           incsol%parameterVecsol(ieff,ip) = vecsol0(isol)
           incsol%parameterPrecis(ieff,ip) = precis0(isol)
       end do

       if ( count(dpa%estime(ic,:)) > 0 ) then
           ieff = ieff +1
           incsol%groupeName(ieff) = 'Dam polygenic effects'
           incsol%nbParameterGroup(ieff)=count(dpa%estime(ic,:))
           ifem=0
           do jm=1,dg%nm
             if ( dpa%estime(ic,jm)) then
              isol=isol+1
              ifem=ifem+1
              if ( vecsol0(isol) ) then
                indest=indest+1
                incsol%paramaterValue(ieff,ifem) = par0(indest)
              else
                incsol%paramaterValue(ieff,ifem) = 0.d0
              end if

              incsol%parameterName(ieff,ifem)='Dam '//trim(dg%mere(jm))
              incsol%parameterVecsol(ieff,ifem) = vecsol0(isol)
              incsol%parameterPrecis(ieff,ifem)  = precis0(isol)
             end if
           end do
         end if

       !Fixed effect
       if ( dpm%modele(ic,1) > 0 ) then
         ieff=ieff+1
         incsol%groupeName(ieff) = 'fixed effects'
         incsol%nbParameterGroup(ieff)=nlevel
         ife=0
         do ief=1,nbef
          do i=1,dpm%nlev(dpm%modele(ic,3+ief))
           ife=ife+1
           isol=isol+1

           if ( vecsol0(isol) ) then
             indest=indest+1
             incsol%paramaterValue(ieff,ife)  = par0(indest)
           else
             incsol%paramaterValue(ieff,ife)  = 0.d0
           end if
           incsol%parameterName(ieff,ife)   = trim(dpm%namefix(dpm%modele(ic,3+ief)))//' level '//trim(str(i))
           incsol%parameterVecsol(ieff,ife) = vecsol0(isol)
           incsol%parameterPrecis(ieff,ife)  = precis0(isol)
          end do
         end do
       end if

       !Covariate
       if ( dpm%modele(ic,2) > 0 ) then
         ieff=ieff+1
         incsol%groupeName(ieff) = 'Covariates'
         incsol%nbParameterGroup(ieff)=dpm%modele(ic,2)

         do ief=1,dpm%modele(ic,2)
           isol=isol+1
           if ( vecsol0(isol) ) then
             indest=indest+1
             incsol%paramaterValue(ieff,ief)  = par0(indest)
           else
             incsol%paramaterValue(ieff,ief)  = 0.d0
           end if
           incsol%parameterName(ieff,ief)   = trim(dpm%namecov(dpm%modele(ic,3+dpm%modele(ic,1)+ief)))
           incsol%parameterVecsol(ieff,ief) = vecsol0(isol)
           incsol%parameterPrecis(ieff,ief)  = precis0(isol)
          end do
       end if

       end subroutine set_solution_hypothesis0

!!****f* m_qtlmap_analyse_discret_unitrait/set_solution_hypothesis1
!! NAME
!!    set_solution_hypothesis1
!! DESCRIPTION
!!   fill the solution under H1 in the structure TYPE_INCIDENCE_SOLUTION
!! INPUTS
!!    ic           : index trait
!! OUTPUTS
!!   incsol        : the variable to fill
!! NOTES
!!***
     subroutine set_solution_hypothesis1(dataset,ic,incsol)
       type(QTLMAP_DATASET)               ,intent(in)       :: dataset
       integer                            ,intent(in)       :: ic
       type(TYPE_INCIDENCE_SOLUTION)      ,intent(inout)    :: incsol

       integer :: nteff,maxNbPar,jm,ifem,ip,ieff
       integer :: ntlev,nbtp,jef,lp,nparx,indest
       integer :: isol,ipar,nlevel,ief,i,ife

       type(GENEALOGY_BASE) , pointer :: dg
       type(PHENOTYPE_BASE) , pointer :: dpa
       type(DATAMODEL_BASE) , pointer :: dpm

       dpm => dataset%phenoModel
       dg => dataset%genea
       dpa => dataset%phenoAnimal

       allocate (incsol%sig(1,dg%np))

       incsol%hypothesis=1
       !  Polygenic family, QTL effect
       nteff = 2
       if ( count(dpa%estime(ic,:))  > 0 ) nteff = 4

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
       nbtp = 3 + dpm%modele(ic,1)+dpm%modele(ic,2)!+dpm%modele(ic,3)
       do jef=1,dpm%modele(ic,3)
		      ntlev=ntlev*dpm%nlev(dpm%modele(ic,nbtp+jef))
	   end do

	   !max nombre de niveau pour un effet fixe en interaction avec le qtl ?
	   maxNbPar = max(maxNbPar,ntlev*dg%np)
	   maxNbPar = max(maxNbPar,ntlev*count(dpa%estime(ic,:)))

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
            incsol%sig(1,ip) = par1(ip)
       end do

       ieff=1
       incsol%qtl_groupeName(1,1)=ieff
       incsol%groupeName(ieff) = 'Sire QTL effects'
       incsol%nbParameterGroup(ieff)=dg%np*ntlevp

       indest = dg%np
       isol=1
       ife=0
       do ip=1,dg%np
         do lp=1,ntlev
           ife=ife+1
           isol=isol+1
           if ( vecsol1(isol) ) then
             indest=indest+1
             incsol%paramaterValue(ieff,ife)  = par1(indest)
           else
             incsol%paramaterValue(ieff,ife)  = 0.d0
           end if

           incsol%parameterName(ieff,ife)   ='Sire '//trim(dg%pere(ip))
           incsol%parameterVecsol(ieff,ife) = vecsol1(isol)
           incsol%parameterPrecis(ieff,ife)  = precis1(isol)
         end do
       end do

       if ( count(dpa%estime(ic,:)) > 0 ) then
           ieff = ieff +1
           incsol%groupeName(ieff) = 'Dam QTL effects'
           incsol%nbParameterGroup(ieff)=count(dpa%estime(ic,:))
           ifem=0
           do jm=1,dg%nm
             if ( dpa%estime(ic,jm)) then
              do lp=1,ntlev
               isol=isol+1
               ifem=ifem+1
               if ( vecsol1(isol) ) then
                indest=indest+1
                incsol%paramaterValue(ieff,ifem) = par1(indest)
               else
                incsol%paramaterValue(ieff,ifem) = 0.d0
               end if
               incsol%parameterName(ieff,ifem)='Dam '//trim(dg%mere(jm))
               incsol%parameterVecsol(ieff,ifem) = vecsol1(isol)
               incsol%parameterPrecis(ieff,ifem)  = precis1(isol)
              end do
             end if
           end do
         end if


       ieff=ieff+1
       incsol%groupeName(ieff) = 'Sire polygenic effects'
       incsol%nbParameterGroup(ieff)=dg%np

       isol=1
       do ip=1,dg%np
           isol=isol+1
           if ( vecsol1(isol) ) then
             indest=indest+1
             incsol%paramaterValue(ieff,ip)  = par1(indest)
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
           do jm=1,dg%nm
             if ( dpa%estime(ic,jm)) then
              isol=isol+1
              ifem=ifem+1
              if ( vecsol1(isol) ) then
                indest=indest+1
                incsol%paramaterValue(ieff,ifem) = par1(indest)
              else
                incsol%paramaterValue(ieff,ifem) = 0.d0
              end if

              incsol%parameterName(ieff,ifem)='Dam '//trim(dg%mere(jm))
              incsol%parameterVecsol(ieff,ifem) = vecsol1(isol)
              incsol%parameterPrecis(ieff,ifem)  = precis1(isol)
             end if
           end do
         end if

       !Fixed effect
       if ( dpm%modele(ic,1) > 0 ) then
         ieff=ieff+1
         incsol%groupeName(ieff) = 'fixed effects'
         incsol%nbParameterGroup(ieff)=nlevel
         ife=0
         do ief=1,nbef
          do i=1,dpm%nlev(dpm%modele(ic,3+ief))
           isol=isol+1
           ife=ife+1
           if ( vecsol1(isol) ) then
             indest=indest+1
             incsol%paramaterValue(ieff,ife)  = par1(indest)
           else
             incsol%paramaterValue(ieff,ife)  = 0.d0
           end if
           incsol%parameterName(ieff,ife)   = trim(dpm%namefix(dpm%modele(ic,3+ief)))//' level '//trim(str(i))
           incsol%parameterVecsol(ieff,ife) = vecsol1(isol)
           incsol%parameterPrecis(ieff,ife)  = precis1(isol)
          end do
         end do
       end if

       !Covariate
       if ( dpm%modele(ic,2) > 0 ) then
         ieff=ieff+1
         incsol%groupeName(ieff) = 'Covariates'
         incsol%nbParameterGroup(ieff)=dpm%modele(ic,2)

         do ief=1,dpm%modele(ic,2)
           isol=isol+1
           if ( vecsol1(isol) ) then
             indest=indest+1
             incsol%paramaterValue(ieff,ief)  = par1(indest)
           else
             incsol%paramaterValue(ieff,ief)  = 0.d0
           end if
           incsol%parameterName(ieff,ief)    = trim(dpm%namecov(dpm%modele(ic,3+dpm%modele(ic,1)+ief)))
           incsol%parameterVecsol(ieff,ief)  = vecsol1(isol)
           incsol%parameterPrecis(ieff,ief)  = precis1(isol)
          end do
       end if

       end subroutine set_solution_hypothesis1

end module m_qtlmap_analyse_discret_unitrait
