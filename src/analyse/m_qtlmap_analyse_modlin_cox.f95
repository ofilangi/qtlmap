!!****m* ANALYSE/m_qtlmap_analyse_modlin_cox
!!  NAME
!!    m_qtlmap_analyse_modlin_cox
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
module m_qtlmap_analyse_modlin_cox
    use m_qtlmap_types
    use m_qtlmap_log
    use m_qtlmap_math
    use m_qtlmap_optimization
    use m_qtlmap_analyse_gen
    use m_qtlmap_analyse_lin_gen
    use m_qtlmap_output_handler

    implicit none

    type(GENEALOGY_BASE) , pointer :: p_dg
    type(PDD_BUILD)      , pointer :: p_spt
    type(PHENOTYPE_BASE) , pointer :: p_dpa

!!****v* m_qtlmap_analyse_modlin_cox/sig1
!!  NAME
!!   sig1
!!  DESCRIPTION
!!   The standart deviation under H0
!!  DIMENSIONS
!!   np
!!***
     real (kind=dp),save       ,dimension(:),allocatable,private   :: sig1
!!****v* m_qtlmap_analyse_modlin_cox/xmu1p
!!  NAME
!!   xmu1p
!!  DESCRIPTION
!!   The polygenic mean for each sire family under H0
!!  DIMENSIONS
!!   np
!!***
     real (kind=dp),save       ,dimension(:),allocatable,private   :: xmu1p
!!****v* m_qtlmap_analyse_modlin_cox/xmu1m
!!  NAME
!!   xmu1m
!!  DESCRIPTION
!!   The polygenic mean for each full sib family
!!  DIMENSIONS
!!   np
!!***
     real (kind=dp),save       ,dimension(:),allocatable,private   :: xmu1m

!!****v* m_qtlmap_analyse_modlin_cox/f0
!!  NAME
!!   f0
!!  DESCRIPTION
!!   value of the likelihood under H0
!!***
     real (kind=dp) ,save                               ,public    :: f0
!!****v* m_qtlmap_analyse_modlin_cox/fp0
!!  NAME
!!   fp0
!!  DESCRIPTION
!!   value of the likelihood by sire family under H0
!!  DIMENSIONS
!!   np
!!
!!***
     real (kind=dp) ,save      ,dimension(:),allocatable,public    :: fp0
!!****v* m_qtlmap_analyse_modlin_cox/fp1
!!  NAME
!!   fp1
!!  DESCRIPTION
!!   value of the likelihood by sire family under H1
!!  DIMENSIONS
!!   np
!!***
     real (kind=dp),save       ,dimension(:),allocatable,private   :: fp1
!!****v* m_qtlmap_analyse_modlin_cox/fm0
!!  NAME
!!   fm0
!!  DESCRIPTION
!!   value of the likelihood by full-sib family under H0
!!  DIMENSIONS
!!   np
!!***
     real (kind=dp) ,save      ,dimension(:),allocatable,public    :: fm0
!!****v* m_qtlmap_analyse_modlin_cox/fm1
!!  NAME
!!   fm1
!!  DESCRIPTION
!!   value of the likelihood by full-sib family under H1
!!  DIMENSIONS
!!   np
!!***
     real (kind=dp) ,save      ,dimension(:),allocatable,private   :: fm1

!!****v* m_qtlmap_analyse_modlin_cox/yr
!!  NAME
!!   fm1
!!  DESCRIPTION
!!
!!  DIMENSIONS
!!
!!***
      real (kind=dp) ,save,dimension(:),allocatable,private        :: yr
!!****v* m_qtlmap_analyse_modlin_cox/yt
!!  NAME
!!   yt
!!  DESCRIPTION
!!
!!  DIMENSIONS
!!
!!***
      real (kind=dp) ,save,dimension(:),allocatable,private        :: yt
!!****v* m_qtlmap_analyse_modlin_cox/rangy
!!  NAME
!!   rangy
!!  DESCRIPTION
!!
!!  DIMENSIONS
!!
!!***
      integer ,save,       dimension(:),allocatable,private      :: rangy
!!****v* m_qtlmap_analyse_modlin_cox/Hy
!!  NAME
!!   Hy
!!  DESCRIPTION
!!
!!  DIMENSIONS
!!
!!***
      real (kind=dp)   ,dimension(:),allocatable, private         ::Hy
!!****v* m_qtlmap_analyse_modlin_cox/Sy
!!  NAME
!!   Sy
!!  DESCRIPTION
!!
!!  DIMENSIONS
!!
!!***
      real (kind=dp)   ,dimension(:),allocatable, private         ::Sy
!!****v* m_qtlmap_analyse_modlin_cox/Ssy
!!  NAME
!!   Ssy
!!  DESCRIPTION
!!
!!  DIMENSIONS
!!
!!***
      real (kind=dp)   ,dimension(:),allocatable, private         ::Ssy

!!****v* m_qtlmap_analyse_modlin_cox/ntot
!!  NAME
!!   ntot
!!  DESCRIPTION
!!
!!***
      integer ,   save,                                 private    :: ntot
!!****v* m_qtlmap_analyse_modlin_cox/posi_h1
!!  NAME
!!   posi_h1
!!  DESCRIPTION
!!
!!***
      integer,    save,               private                      :: posi_h1
!!****v* m_qtlmap_analyse_modlin_cox/current_ic
!!  NAME
!!   current_ic
!!  DESCRIPTION
!!   The current trait to analyse
!!***
    integer ,   save           ,private                            :: current_ic
!!****v* m_qtlmap_analyse_modlin_cox/current_chr
!!  NAME
!!   current_chr
!!  DESCRIPTION
!!   The current chromosome while the likelihood calculs
!!***
    integer ,  save           ,private                            :: current_chr



    !$omp threadprivate (Ssy,Sy,Hy,posi_h1,current_chr,fp1,fm1)

    public :: init_analyse_modlin_cox
    public :: opti_0qtl_modlin_cox
    public :: opti_1qtl_modlin_cox
    public :: test_lin_cox
    public :: end_analyse_modlin_cox
    public :: set_solution_hypothesis0
    public :: set_solution_hypothesis1

    contains
!!****f* m_qtlmap_analyse_modlin_cox/init_analyse_modlin_cox
!!  NAME
!!    init_analyse_modlin_cox
!!  DESCRIPTION
!!    Initialisation/allocation of solution/buffer arrays
!!
!!  NOTES
!!  SOURCE
     subroutine init_analyse_modlin_cox(dataset)
         type(QTLMAP_DATASET)       ,intent(in) :: dataset

         integer           :: stat
         type(GENEALOGY_BASE) , pointer :: dg

         dg => dataset%genea

         allocate (sig1(dg%np),STAT=stat)
         call check_allocate(stat,'sig1 [m_qtlmap_analyse_modlin]')
         allocate (xmu1p(dg%np),STAT=stat)
         call check_allocate(stat,'xmu1p [m_qtlmap_analyse_modlin]')
         allocate (xmu1m(dg%nm),STAT=stat)
         call check_allocate(stat,'xmu1m [m_qtlmap_analyse_modlin]')
         allocate (fp0(dg%np),STAT=stat)
         call check_allocate(stat,'fp0 [m_qtlmap_analyse_modlin_cox]')
        ! allocate (fp1(np),STAT=stat)
        ! call check_allocate(stat,'fp1 [m_qtlmap_analyse_modlin_cox]')
         allocate (fm0(dg%nm),STAT=stat)
         call check_allocate(stat,'fm0 [m_qtlmap_analyse_modlin_cox]')
        ! allocate (fm1(nm),STAT=stat)
        ! call check_allocate(stat,'fm1 [m_qtlmap_analyse_modlin_cox]')
! allocations used in calcul_rang
         allocate (yr(dg%nd),STAT=stat)
         call check_allocate(stat,'yr    [m_qtlmap_analyse_modlin_cox]')
         allocate (yt(dg%nd),STAT=stat)
         call check_allocate(stat,'yt    [m_qtlmap_analyse_modlin_cox]')
         allocate (rangy(dg%nd),STAT=stat)
         call check_allocate(stat,'rangy [m_qtlmap_analyse_modlin_cox]')

! allocation used in opti_0qtl_modlin_cox
       !  allocate (Hy(nd),STAT=stat)
       !  call check_allocate(stat,'Hy    [m_qtlmap_analyse_modlin_cox]')
       !  allocate (Sy(nd),STAT=stat)
       ! call check_allocate(stat,'Sy    [m_qtlmap_analyse_modlin_cox]')
       !  allocate (Ssy(nd),STAT=stat)
       !  call check_allocate(stat,'Ssy   [m_qtlmap_analyse_modlin_cox]')

     end subroutine init_analyse_modlin_cox
!!***





!!****f* m_qtlmap_analyse_modlin_cox/end_analyse_modlin_cox
!!  NAME
!!    end_analyse_modlin_cox
!!  DESCRIPTION
!!    deallocation of solution/buffer arrays
!!
!!  NOTES
!!  SOURCE
     subroutine end_analyse_modlin_cox
         deallocate (sig1)
         deallocate (xmu1p)
         deallocate (xmu1m)
         deallocate (fp0)
         !deallocate (fp1)
         deallocate (fm0)
         !deallocate (fm1)
         deallocate (yt)
         deallocate (yr)
         deallocate (rangy)
        ! deallocate (Hy)
        ! deallocate (Sy)
        ! deallocate (Ssy)
     end subroutine end_analyse_modlin_cox
!!***

!!****f* m_qtlmap_analyse_modlin_cox/calcul_rang
!!  NAME
!!    calcul_rang
!!  DESCRIPTION
!!     Calcul du rang de y
!!
!!  INPUTS
!!    ic : index of the trait
!!  SOURCE
   subroutine calcul_rang(dataset,ic)
      integer                    ,intent(in) :: ic
      type(QTLMAP_DATASET)       ,intent(in) :: dataset

      integer ,	dimension(:),allocatable      :: nryt
      integer :: nm1,nm2,nd1,nd2, ifail,kd,jm,ip, nt
      integer :: stat
      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa

      dg => dataset%genea
      dpa => dataset%phenoAnimal
!
      allocate (nryt(dg%nd),STAT=stat)
      call check_allocate(stat,'nryt  [m_qtlmap_analyse_modlin_cox]')

      ntot=0
      do ip=1,dg%np
          nm1=dg%nmp(ip)+1; nm2=dg%nmp(ip+1)
          do jm=nm1,nm2
             nd1=dg%ndm(jm)+1; nd2=dg%ndm(jm+1)
             do kd=nd1,nd2
                if(dpa%presentc(ic,kd))  then
		    ntot=ntot+1
                    yt(ntot)=dpa%y(ic,kd)
                    yr(ntot)=0.d0
		    nryt(ntot)=1
                endif
             enddo
          enddo
       enddo

      ifail=0
      CALL MATH_QTLMAP_M01DAF(yt,1,Ntot,'Descending',nryt,ifail)
!      write(nficout,*) (i,yt(i),nryt(i),i=1,nt)
        nt=0
         do ip=1,dg%np
           nm1=dg%nmp(ip)+1; nm2=dg%nmp(ip+1)
           do jm=nm1,nm2
              nd1=dg%ndm(jm)+1; nd2=dg%ndm(jm+1)
              do kd=nd1,nd2
        	if(dpa%presentc(ic,kd))  then
        	   nt=nt+1
        	   rangy(kd)=nryt(nt)
        	   yr(rangy(kd))=dpa%y(ic,kd)
        	!   write(nficout,*) ic,kd, nt, y(ic,kd), rangy(kd)
        	endif
              enddo
           enddo
        enddo

	deallocate (nryt)

  end subroutine calcul_rang
!!***

!!****f* m_qtlmap_analyse_modlin_cox/opti_0qtl_modlin_cox
!!  NAME
!!    opti_0qtl_modlin_cox
!!  DESCRIPTION
!!     Calcul de la vraisemblance 0 QTL, 1 caractere , effets parasites inclus
!!
!!  INPUTS
!!    ic       : index of the trait
!!   est_moy   : general mean
!!  SOURCE
      subroutine opti_0qtl_modlin_cox(dataset,spt,ic,est_moy)
      type(QTLMAP_DATASET)       ,intent(in)         :: dataset
      type(PDD_BUILD)            ,intent(in)         :: spt

      integer , intent(in)   :: ic
      logical , intent(in)   :: est_moy
!
! Divers

      integer                                    :: iuser(1)
      integer        ,dimension(:),allocatable   :: iw
      real (kind=dp) ,dimension(:),allocatable   :: par,borni,borns,w
      real (kind=dp)                             :: user(1)
      real (kind=dp)                             :: f
      integer :: npar,ibound,ip,ix,ifail,i,indest,ntotp
      integer :: km,jm
      logical itest
      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa
      type(DATAMODEL_BASE) , pointer :: dpm

      dpm => dataset%phenoModel
      dg => dataset%genea
      dpa => dataset%phenoAnimal

      p_dg=> dataset%genea
      p_dpa => dataset%phenoAnimal

      call log_mess("START opti_0qtl_modlin_cox",DEBUG_DEF)
!
!******************************************************************************
!
!  recherche des elements estimables � partir de la d�composition de choleski
!  comptage et recodification de ces elements (on commence par mettre � 0 les
!  compteurs destin�s aux tests des effes fix�s
!
      !  allocate (Hy(nd),STAT=stat)
       !  call check_allocate(stat,'Hy    [m_qtlmap_analyse_modlin_cox]')
       allocate (Sy(dg%nd))
       allocate (Ssy(dg%nd))
       !  call check_allocate(stat,'Ssy   [m_qtlmap_analyse_modlin_cox]')

      itest=.false.
      !print *, 'est_moy',est_moy
      call log_mess("contingence cox",DEBUG_DEF)
      call contingence_cox(dataset,spt,ic,0,itest,est_moy)
      call log_mess("precision",DEBUG_DEF)
      call precision(xx,precis)
      !  Parametres de maximisation
      !npar=np+nbnivest
       !on n'a pas � prendre en compte un ecart type par pere avec cox
       npar=nbnivest
      call log_mess("NBNIVEST:"//trim(str(npar)),DEBUG_DEF)
      allocate (borni(npar))
      allocate (borns(npar))
      allocate (par(npar))

      ibound=0
      ! pas besoin d'initialiser les ecart types par pere comme il n' y en a pas
      !do ip=1,np
      !  borni(ip)=1.d-6
       ! borns(ip)=1.d6
      !end do
      !do ix=np+1,npar
! on a passé les bornes de 1.d4 à 50 car sinon on risque un plantage avec exp()
! un risque de exp(100)=5.d21 s'est largement suffisant /ref
      do ix=1,npar
        borni(ix)=1.d-6
        borns(ix)=1.d4
      end do
!
! Point de depart
! sous cox il y a ni moyenne ni ecart type intra pere
      !par(np+1)=1.d0
      !do ip=1,np
        !par(ip)=sig0(ip)
        !par(np+1)=par(np+1)+xmu0p(ip)
      !end do
      !par(np+1)=par(np+1)/dble(np)
      !do ix=np+2,npar
      do ix= 1, npar
        par(ix)=1.d0
      end do
!
! Optimisation de la vraisemblance
      ifail=1
      current_ic = ic ! pour le mode lineaire....
!print *, 'npar sous H0',npar
      call log_mess("minimizing_funct",DEBUG_DEF)
      if (npar==0) then
         call funct_0qtl_modlin_cox(npar,par,f,iuser,user)
      else
         call minimizing_funct(dataset,npar,ibound,funct_0qtl_modlin_cox,borni,borns,par,f,iuser,user,ifail)
      endif
      f0=f
     ! print *, 'Likelihood sous HO=', f0
      do i=1,npar
        par0(i)=par(i)
     ! print *, 'sous H0: par(i) =', par(i), 'npar=',npar
      end do
      do i=1,ntniv
        vecsol0(i)=vecsol(i)
        precis0(i)=precis(i)
      end do

      do ip = 1,dg%np
       ! sig1(ip)=par(ip)
        sig1(ip)=0.d0
      end do

    !  indest=1
      indest=0
      if (dg%np.ge.2) then
        do ip = 1,dg%np
    !     if (vecsol(1+ip)) then
          if (vecsol(ip)) then
            indest=indest+1
    !       xmu1p(ip)=par(np+indest)
            xmu1p(ip)=par(indest)
          end if
        end do
      endif

   !   ntot=np+1
      ntotp=0
      if (dg%np.ge.2)  ntotp=dg%np
      km=0
      if (dg%nm.ge.2) then
        do jm=1,dg%nm
          if (dpa%estime(ic,jm).and.dataset%params%opt_sib.eq.2)then
            km=km+1
            if(vecsol(ntotp+km)) then
              indest=indest+1
             ! xmu1m(jm)=par(np+indest)
              xmu1m(jm)=par(indest)
            end if
          end if
        end do
      endif
      deallocate (borni)
      deallocate (borns)
      deallocate (par)
      deallocate (Sy)
      deallocate (Ssy)

      call log_mess("END opti_0qtl_modlin_cox",DEBUG_DEF)
      return
      end subroutine opti_0qtl_modlin_cox
!!***

!!****f* m_qtlmap_analyse_modlin_cox/funct_0qtl_modlin_cox
!!  NAME
!!    funct_0qtl_modlin_cox
!!  DESCRIPTION
!!     Calcul de la vraisemblance et de ses derivees partielles sous H0
!!
!!  NOTES
!! ------------------- MODELE DE COX ET METHODE SANS DERIVES-----------------
!!   Calcul du max de vraisemblance sous H1 et H0
!!   Sous programme appele par e04jyf, optimiseur utilise par cherche
!!   LE DENOMINATEUR DE l(kd,i) est approxime:
!!        LE CALCUL DE LA VRAISEMBLANCE EST REALISE EN 3 phases:
!!   Dans la 1ere boucle, on calcul pour chaque animal de Hy(kd) utilise dans la
!!   3eme boucle et de Sy(rang(ic,kkd)) utilise dans la 2eme boucle, pour calculer
!!   les contributions de chaque animal.
!!   la 3eme boucle calcule la vraisemblance proprement dit.
!!
!!  INPUTS
!!    ic       : index of the trait
!!   est_moy   : general mean
!!  SOURCE
      subroutine funct_0qtl_modlin_cox(n,x,f,iuser,user)
      use m_qtlmap_analyse_lin_gen, only : nbfem,nbef,ntnifix,nbco,covdir,nivdir,corniv,vecsol

      integer         , intent(in)                  :: n
!
! Tableaux dimensionnes selon n, le nombre de parametres a estimer
!
      real (kind=dp)      ,dimension(n), intent(in) :: x
      real (kind=dp)  , intent(inout)               :: f
      integer ,       dimension(1), intent(in)      :: iuser
      real (kind=dp)      ,dimension(1), intent(in) :: user

! Tableaux dimensionnes selon nm, le nombre de meres
      real (kind=dp)   ,dimension(p_dg%nm)    :: effm


      integer ::jnit,ip,nm1,nm2,jm,nd1,nd2,kkd,ilev,ief,ico,icar,iid,i
      real (kind=dp) :: vmere,vpf,v,l,dv
!
!
!   INITIALISATION
      icar = current_ic
      if (p_dg%np.lt.2.and.(count(p_dpa%estime(icar,:))).lt.2) jnit=0
      if (p_dg%np.ge.2.and.(count(p_dpa%estime(icar,:))).lt.2) jnit=1
      if (p_dg%np.ge.2.and.(count(p_dpa%estime(icar,:))).ge.2) jnit=2

      f=0.d0

!   CALCULS pour chaque descendant du terme de participation � la contribution
      Sy=0.d0
      do ip=1,p_dg%np
         nm1=p_dg%nmp(ip)+1
         nm2=p_dg%nmp(ip+1)
         do jm=nm1,nm2
! on ne consid�re que les m�res
          nd1=p_dg%ndm(jm)+1
          nd2=p_dg%ndm(jm+1)
             do kkd=nd1,nd2
		IF(p_dpa%presentc(icar,kkd)) then
		  v=1.d0
! effet du pere
                  if (p_dg%np.ge.2) then
                    if(vecsol(nivdir(kkd,1))) then
                      ilev=corniv(nivdir(kkd,1))
                      v=v*x(ilev)
                    endif
                  end if
!effet de la m�re
                  if (nbfem.ge.2) then
                    if(p_dpa%estime(icar,jm)) then
                      if(vecsol(nivdir(kkd,2))) then
                        ilev=corniv(nivdir(kkd,2))
                        v=v*x(ilev)
                      end if
                    end if
                  endif
!effet fixes
                  do ief=jnit+1,jnit+nbef
                    if(vecsol(nivdir(kkd,ief))) then
                      ilev=corniv(nivdir(kkd,ief))
                      v=v*x(ilev)
                    end if
                  end do
! effet des covariables
                  do ico=1,nbco
                    if(vecsol(ntnifix+ico)) then
                      ilev=corniv(ntnifix+ico)
                      v=v*(x(ilev)**covdir(kkd,ico))
                    end if
                  end do

		  dv = (v)
                   ! if ( dv /= dv )
		   if ( dv > huge(dv)) then
		       call log_mess('algorithm goes out of bound, iteration was left!!!!!!',DEBUG_DEF)
		       fm1=INIFINY_REAL_VALUE
	               fp1=INIFINY_REAL_VALUE
	               f=INIFINY_REAL_VALUE
	               return
		   end if
                  Sy(rangy(kkd))=dv
                ENDIF                   !fin de la condition sur les perf
	     enddo                      !fin boucle sur kd
	enddo                          !fin de la boucle sur jm
      enddo                           !fin de la boucle sur ip

!   CALCUL DES CONTRIBUTIONS
      do  ip=1,p_dg%np
         fp0(ip)=0.d0
         nm1=p_dg%nmp(ip)+1
         nm2=p_dg%nmp(ip+1)
         do  jm=nm1,nm2
	     fm0(jm)=0.d0
             nd1=p_dg%ndm(jm)+1
	     nd2=p_dg%ndm(jm+1)
             do kkd=nd1,nd2
		IF(p_dpa%presentc(icar,kkd)) then
                  Ssy(kkd)=0.d0
      ind_sort :  do iid=1,ntot
		     if (yr(rangy(kkd)).GT.yr(iid)) exit ind_sort
		     Ssy(kkd)=Ssy(kkd)+Sy(iid)
		  enddo ind_sort
		  IF ((p_dpa%ndelta(icar,kkd).EQ.1))  then
		     if (Sy(rangy(kkd)) == 0.d0.or.SSy(kkd)==0.d0) then
                         fm0(jm)=INIFINY_REAL_VALUE
                         fp0=INIFINY_REAL_VALUE
                        f=INIFINY_REAL_VALUE
                         call log_mess('algorithm goes out of bound, iteration was left!!!!!!',DEBUG_DEF)
			 return
                     else
			! print *,kkd, 'ssy:',SSy(kkd), 'sy:', sy(rangy(kkd))
                         fm0(jm)= fm0(jm)-(dlog(Sy(rangy(kkd)))-dlog(SSy(kkd)))
                     !print *,'SOUS H0', ndelta(icar,kkd),kkd,Sy(rangy(kkd)),SSy(kkd), 'fm0', fm0(jm)
		     endif
                 ENDIF
		ENDIF                  !fin de la condition sur les perf
               ! write(nficout,*) icar,kkd, rangy(kkd),ndelta(icar,kkd), l
	     enddo                     !fin boucle sur kkd
             fp0(ip)=fp0(ip)+fm0(jm)
             f=f+fm0(jm)
!print *, 'FM0', fm0(jm)
          enddo                         !fin boucle sur jm
!print *, 'FP0', fp0(ip)
       enddo                         !fin boucle sur ip

      ! print *,f
      return
      end subroutine funct_0qtl_modlin_cox
!!***

!!****f* m_qtlmap_analyse_modlin_cox/opti_1qtl_modlin_cox
!!  NAME
!!    opti_1qtl_modlin_cox
!!  DESCRIPTION
!!     Calcul de la vraisemblance 1 QTL, 1 caractere , effets parasites inclus
!!
!!
!!  INPUTS
!!    ic       : index of the trait
!!  OUTPUTS
!!   lrtsol     : LRT : curves and maximum under H(nqtl)
!!   fmax       : maximum value of the likelihood
!!  supnbnivest : number of parameter estimated at the maximum value
!!
!!  SOURCE
      subroutine opti_1qtl_modlin_cox(dataset,spt,ic,lrtsol,fmax,supnbnivest)

      integer                              , intent(in)   :: ic
      type(TYPE_LRT_SOLUTION)  , intent(out)              :: lrtsol
      real (kind=dp)  , intent(out)                       :: fmax
      integer         , intent(out)                       :: supnbnivest
      type(QTLMAP_DATASET)         ,intent(in)            :: dataset
      type(PDD_BUILD)          ,target   ,intent(in)      :: spt
!
!
! Divers

      integer                                    :: iuser(1)
      integer        ,dimension(:),allocatable   :: iw
      real (kind=dp) ,dimension(:),allocatable   :: val,par,borni,borns,w
      real (kind=dp)                             :: user(1)
      real (kind=dp)                             :: f
      integer :: npar,indest
      integer :: km,jm,chr,ifem
      real(kind=dp)       ,dimension(:,:) ,pointer           :: listef
      integer             ,dimension(:,:)    ,pointer        :: listenbnivest
      real(kind=dp)       ,dimension(:,:,:)   ,pointer       :: listepar

      real (kind=dp) :: xlrt_t,f1

      logical itest,stvecsol(size(vecsol1))
      integer :: ip,i,n,ilong,ibound,j,ifail,ii,iip,indexchr(dataset%map%nchr),ntotal,nprime
      integer :: save_nteffmax,save_ntnivmax
      type(GENEALOGY_BASE) , pointer :: dg
      type(PHENOTYPE_BASE) , pointer :: dpa
      real(kind=dp)       ,dimension(:,:,:)   ,pointer     :: xlrp,xlrm
      real(kind=dp)       ,dimension(:,:)   ,pointer       :: lrt1

      dg => dataset%genea
      dpa => dataset%phenoAnimal

      p_dg => dataset%genea
      p_spt => spt
      p_dpa => dataset%phenoAnimal

      call log_mess("START opti_1qtl_modlin_cox",DEBUG_DEF)
      call lrtsol%new(dataset,1)
!******************************************************************************
!
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

      ibound=0
      ! pas besoin d'initialiser les ecart types par pere comme il n' y en a pas
      !do ip=1,np
      !  borni(ip)=1.d-6
       ! borns(ip)=1.d6
      !end do
      !do ix=np+1,npar

      do i=1,size(vecsol)
        vecsol(i)=.true.
      end do

      lrtsol%lrtmax=-1.d75
      lrtsol%chrmax=0
      lrtsol%nxmax=0

      current_ic = ic

      ntotal=0
      do chr=1,dataset%map%nchr
        ntotal=ntotal+dataset%map%get_npo(chr)
        indexchr(chr)=ntotal
      end do

      allocate (listef(dataset%map%nchr,dataset%map%get_maxnpo()))
      allocate (listenbnivest(dataset%map%nchr,dataset%map%get_maxnpo()))
      allocate (listepar(dataset%map%nchr,dataset%map%get_maxnpo(),dg%np + ntnivmax))

     call end_contingence
     save_ntnivmax = ntnivmax
     save_nteffmax = nteffmax

     allocate (lrt1(dataset%map%nchr,dataset%map%get_maxnpo()))
     allocate (xlrp(dataset%map%nchr,dataset%map%get_maxnpo(),dg%np))
     allocate (xlrm(dataset%map%nchr,dataset%map%get_maxnpo(),dg%nm))

     !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(par,val)  &
     !$OMP PRIVATE(i,ifem,iip,ii,stvecsol,npar,j,ifail,f1,n,chr,borni,borns)
     ntnivmax = save_ntnivmax
     nteffmax = save_nteffmax
     allocate ( val( ntnivmax ) )
     allocate ( par( ntnivmax) )
     allocate ( borni( ntnivmax) )
     allocate ( borns( ntnivmax) )
     allocate (fp1(dg%np))
     allocate (fm1(dg%nm))
     allocate (Hy(dg%nd))
     allocate (Sy(dg%nd))
     allocate (Ssy(dg%nd))

     !
     ! Point de depart
     !
      do i=1,size(par)
        par(i)=1.d0
      end do

      do i=1,size(val)
        val(i)=par(i)
      end do


     call init_contingence(dataset,spt)
!
! Marche le long du chromosome
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
      call contingence_cox(dataset,spt,ic,1,itest,.false.)
      !  Parametres de maximisation
       npar=nbnivest
! on a passé les bornes de 1.d4 à 50 car sinon on risque un plantage avec exp()
! un risque de exp(100)=5.d21 s'est largement suffisant /ref
     do i=1,npar
        borni(i)=1.d-5
        borns(i)=1.d4
      end do

!
! Point de depart
! sous cox il y a ni moyenne ni ecart type intra pere

      do i= 1, npar
        par(i)=1.d0
      end do

      do i=1,ntniv
        if(vecsol(i))then
          par(i)=1.d0
          if(stvecsol(i)) par(i)=val(i)
        end if
        stvecsol(i)=vecsol(i)
      end do
!
! Optimisation de la vraisemblance
      ifail=1
      posi_h1=n
      call minimizing_funct(dataset,npar,ibound,funct_1qtl_modlin_cox,borni,borns,par,f1,iuser,user,ifail)
      print *,n,f1
    ! print *, -2.d0*(fp1(ii)-fp0(ii))
      do i=1,ntniv
        if(vecsol(i)) then
          val(i)=par(i)
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
           lrt1(chr,n)=0.d0
        end if

        listef(chr,n)=f1
        listenbnivest(chr,n)=nbnivest
        listepar(chr,n,:npar)=par(:npar)
!
!  on met les profil et effets QTL / prog�niteur dnas les ficheir ad hoc
!
        iip=0
        do ii=1,dg%np
          xlrp(chr,n,ii)=-2.d0*(fp1(ii)-fp0(ii))
         if ( vecsol(ii)) then
           iip=iip+1
           ! ATTENTION VALIDE DANS LE CAS SANS INTERACTION AVEC QTL
           lrtsol%pater_eff(chr,ii,n)=par(iip)

          end if
        end do

        ifem=0
        do ii=1,dg%nm
          xlrm(chr,n,ii)=-2.d0*(fm1(ii)-fm0(ii))
          if ( dpa%estime(ic,ii)) then
            ifem=ifem+1
            if ( vecsol(dg%np+ifem) ) then
              iip=iip+1
              ! ATTENTION VALIDE DANS LE CAS SANS INTERACTION AVEC QTL
              lrtsol%mater_eff(chr,ifem,n)=par(iip)
            end if
          end if
        end do

      end do
      !$OMP END DO
      call end_contingence
      deallocate ( val )
      deallocate ( par )
      deallocate (fp1)
      deallocate (fm1)
      deallocate (Hy)
      deallocate (Sy)
      deallocate (Ssy)
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
            par1(:supnbnivest+dg%np)=listepar(chr,n,:supnbnivest+dg%np)
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


!
!  on calcule la pr�cision des estimation au point correspondant
! au LRT maximum
!
      call init_contingence(dataset,spt)
      call prepinc(dataset,spt,lrtsol%chrmax(0),lrtsol%nxmax(0),ic,"LA  ")
      call contingence_cox(dataset,spt,ic,1,itest,.false.)
      call precision(xx,precis)

      do i=1,ntniv
        precis1(i)=precis(i)
        vecsol1(i)=vecsol(i)
      end do

      call log_mess("END opti_1qtl_modlin_cox",DEBUG_DEF)
      return
      end subroutine opti_1qtl_modlin_cox
!!***

!!****f* m_qtlmap_analyse_modlin_cox/funct_1qtl_modlin_cox
!!  NAME
!!    funct_1qtl_modlin_cox
!!  DESCRIPTION
!!    Calcul de la vraisemblance et de ses derivees partielles sous H1
!!
!!  NOTES
!!    ------------------- MODELE DE COX ET METHODE SANS DERIVES-----------------
!!             Calcul du max de vraisemblance sous H1
!!     Sous programme appele par e04jyf, o
!!         LE DENOMINATEUR DE l(kd,i) est approxim�:
!!        LE CALCUL DE LA VRAISEMBLANCE EST REALISE EN 3 phases:
!!   Dans la 1�boucle, on calcul pour chaque animal de Hy(kd) utilis� dans la
!!   3� boucle et de Sy(rang(ic,kkd)) utilis� dans la 2�boucle, pour calculer
!!   les contributions de chaque animal.
!!   la 3�boucle calcule la vraisemblance proprement dit.
!!
!!  INPUTS
!!    n         : number of parameter
!!    x         : the parameter's values
!!   iuser      :
!!   user       :
!!  OUTPUTS
!!   f          : the value of the likelihood
!!
!!  SOURCE
      subroutine funct_1qtl_modlin_cox(n,x,f,iuser,user)
      use m_qtlmap_analyse_lin_gen, only : nbfem,nbef,ntnifix,nbco,covdir,nivdir,corniv,vecsol

      integer         , intent(in)                  :: n
!
! Tableaux dimensionnes selon n, le nombre de parametres a estimer
!
      real (kind=dp)      ,dimension(n), intent(in) :: x
      real (kind=dp)  , intent(inout)               :: f
      integer ,       dimension(1), intent(in)      :: iuser
      real (kind=dp)      ,dimension(1), intent(in) :: user

      real (kind=dp)                                :: my, gy,dvq


      integer ::init,jnit,ip,nm1,nm2,jm,nd1,nd2,kd,kkd,ilev,ief,ico,icar,nt,iid,i,ig, ngeno1, ngeno2, nbdv
      real (kind=dp) :: vmere,vpf,v, vq,l
!
!     
      fp1=0.d0
      fm1=0.d0
      f=0.d0
      icar = current_ic
! jnit: numero du 1° effet fixe, init: numero du 1° effet polygénique père
      if (p_dg%np==1.and.(count(p_dpa%estime(icar,:))).lt.2)   jnit=2
      if (p_dg%np==1.and.(count(p_dpa%estime(icar,:)))==0)     jnit=1
      if (p_dg%np==1.and.(count(p_dpa%estime(icar,:))).ge.2)   jnit=3
      init=3
      if (p_dg%np.ge.2.and.(count(p_dpa%estime(icar,:))).lt.2) jnit=3
      if (p_dg%np.ge.2.and.(count(p_dpa%estime(icar,:)))==0  ) then;jnit=2; init=2; endif
      if (p_dg%np.ge.2.and.(count(p_dpa%estime(icar,:))).ge.2) jnit=4

      f=0.d0
      Sy = 0.d0
!   CALCULS pour chaque descendant du terme de participation � la contribution
b_pe: do ip=1,p_dg%np
        nm1=p_dg%nmp(ip)+1
        nm2=p_dg%nmp(ip+1)
b_me:   do jm=nm1,nm2
          ngeno1=p_spt%ngenom(current_chr,jm)+1
          ngeno2=p_spt%ngenom(current_chr,jm+1)
b_gme:    do ig=ngeno1,ngeno2
            nd1=p_spt%ngend(current_chr,ig)+1
            nd2=p_spt%ngend(current_chr,ig+1)
            vpf=0.d0
b_indgm:    do kd=nd1,nd2
              Hy(kd)=0.d0
              kkd=p_spt%ndesc(current_chr,kd)
              if(p_dpa%presentc(icar,kkd)) then
		  v=1.d0
! effet du pere
                  if (p_dg%np.ge.2) then
                    if(vecsol(nivdir(kkd,init))) then
                       ilev=corniv(nivdir(kkd,init))
                       v=v*x(ilev)
                    end if
                  endif
!effet de la m�re

                  !print *,'nivdir(kkd,4):',nivdir(kkd,4),' s1:',size(vecsol)
                
                  if (p_dg%nfem.ge.2) then
                    if(p_dpa%estime(icar,jm)) then
                      if(vecsol(nivdir(kkd,4))) then
                        ilev=corniv(nivdir(kkd,4))
                        v=v*x(ilev)
                      end if
                    end if
                  endif
!effet fixes
                  do ief=jnit+1,jnit+nbef
                    if(vecsol(nivdir(kkd,ief))) then
                      ilev=corniv(nivdir(kkd,ief))
                      v=v*x(ilev)
                    end if
                  end do
! effet des covariables
                  do ico=1,nbco
                    if(vecsol(ntnifix+ico)) then
                      ilev=corniv(ntnifix+ico)
                      v=v*(x(ilev))**(covdir(kkd,ico))
                    end if
                  end do
! calcul des contribution en fonction des effets QTL pere et mere
b_pdd :        do i=1,4
 ! effet QTL du pere
                  vq=v
		  if(vecsol(nivdir(kkd,1))) then
                     ilev=corniv(nivdir(kkd,1))
                     if (i==1.or.i==2) then
		         vq=vq
		     else
		         vq=vq*x(ilev)
	             endif
                  end if
!effet QTL de la m�re
                   if(p_dpa%estime(icar,jm)) then
		     if (vecsol(nivdir(kkd,2))) then
                        ilev=corniv(nivdir(kkd,2))
                        if (i==1.or.i==3) then
		           vq=vq
		        else
		           vq=vq*x(ilev)
	                endif
                     end if
                   end if
		 ! if (vq.GT.1000.d0 .or. vq.LT.-1000.d0) then
                  !  print *, 'vq=',vq, 'Vq est trop GRAND!'
                   ! stop
                 ! endif
		   dvq = (vq)
                   !print *, 'kd=', kd, 'dvq=', dvq
		   if ( dvq > huge(dvq)) then
                      ! print *, 'OUT kd=', kd, 'dvq=', dvq
                       call log_mess('algorithm goes out of bound, iteration was left!!!!!!',DEBUG_DEF)
                       print *, 'gros dvq=', dvq
		       fm1=INIFINY_REAL_VALUE
	               fp1=INIFINY_REAL_VALUE
	               f=INIFINY_REAL_VALUE
	               return
		   endif

		   Hy(kd)=(p_spt%pdd(current_chr,kd,i,posi_h1)*(dvq))+Hy(kd)
	         enddo b_pdd
                 Sy(rangy(kkd))=Sy(rangy(kkd))+(p_spt%probg(current_chr,ig)*Hy(kd))

               ENDIF
              enddo b_indgm
            enddo   b_gme
           enddo    b_me
         enddo      b_pe
 !   CALCUL DES DENOMINATEUR POUR CHAQUE KKD
      do  ip=1,p_dg%np
         nm1=p_dg%nmp(ip)+1
         nm2=p_dg%nmp(ip+1)
         do jm=nm1,nm2
             nd1=p_dg%ndm(jm)+1
	     nd2=p_dg%ndm(jm+1)
	     fm0(jm)=0.d0
             do kkd=nd1,nd2
                  IF(p_dpa%presentc(icar,kkd)) then
                    Ssy(kkd)=0.d0
 ind_sort :	    do iid=1,ntot
		       if (yr(rangy(kkd)).GT.yr(iid)) exit ind_sort
		       Ssy(kkd)=Ssy(kkd)+Sy(iid)
		    enddo ind_sort
                  ENDIF                      !fin de la condition sur les perf
	     enddo                         !fin boucle sur kkd
! calcul de la vraisemblance
             vmere=0.d0
             nbdv=0
	     ngeno1=p_spt%ngenom(current_chr,jm)+1
             ngeno2=p_spt%ngenom(current_chr,jm+1)
             do ig=ngeno1,ngeno2
               gy=1.d0
	       nd1=p_spt%ngend(current_chr,ig)+1
               nd2=p_spt%ngend(current_chr,ig+1)
               vpf=0.d0
               do kd=nd1,nd2
                 kkd=p_spt%ndesc(current_chr,kd)
	   	 if(p_dpa%presentc(icar,kkd).and.p_dpa%ndelta(icar,kkd).EQ.1) then
                     nbdv=nbdv+1
! condition pour éviter des erreurs arithmetiques du genre division par 0 huge= plus grande valeur possible et tiny=plus petite valeur possible
                   !  if (Ssy(kkd).GE.huge(Ssy(kkd))) Ssy(kkd)=huge(Ssy(kkd))
                  !   if (Ssy(kkd).LE.tiny(Ssy(kkd))) Ssy(kkd)=tiny(Ssy(kkd))
		     gy=gy*(Hy(kd)*ntot)/(Ssy(kkd))
                    ! print *,ntot,'sous H1', kd,Hy(kd),Ssy(kkd), 'gy',gy, 'probm', probg(current_chr,ig)
		 endif
	       enddo !boucle kkd
               vmere= vmere+(p_spt%probg(current_chr,ig)*gy)
               !print *, 'vmere=', vmere
              
	     enddo  !boucle ig
             if (vmere == 0) then
               call log_mess('algorithm goes out of bound, iteration was left!!!!!!',DEBUG_DEF)
               fm1=INIFINY_REAL_VALUE
	       fp1=INIFINY_REAL_VALUE
	       f=INIFINY_REAL_VALUE
	       return
             else
	       fm1(jm)=-dlog(vmere)+(nbdv*(dlog(ntot*1.d0)))
              ! print *, 'fm', fm1(jm), vmere
	     endif
             fp1(ip)=fp1(ip)+fm1(jm)
             f=f+fm1(jm)
	   enddo                         !fin boucle sur jm
        enddo                         !fin boucle sur np
      !print *,'f:',f,'effet QTL', x(1)
      return
      end subroutine funct_1qtl_modlin_cox
!!***

!!****f* m_qtlmap_analyse_modlin_cox/test_lin_cox
!!  NAME
!!    test_lin_cox
!!  DESCRIPTION
!!    Test des differents effets de nuisance du mod�le par une LRT comp�r� � une chi2
!!
!!  INPUTS
!!
!!  OUTPUTS
!!
!!
!!  SOURCE
    subroutine test_lin_cox(dataset,spt,chr,ic,est_moy,supnbnivest,fmax,nposx)
      type(QTLMAP_DATASET)           ,intent(in)  :: dataset
      type(PDD_BUILD)                ,intent(in)  :: spt
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
      type(GENEALOGY_BASE) , pointer :: dg
      type(DATAMODEL_BASE) , pointer :: dpm

      dg => dataset%genea
      dpm => dataset%phenoModel

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

      nbef=dpm%modele(ic,1)
      nbco=dpm%modele(ic,2)
      nbint=dpm%modele(ic,3)

      iecd=0
      iecq=0

      allocate ( par( dg%np + ntnivmax  ) )
      allocate ( borni( dg%np + ntnivmax ) )
      allocate ( borns( dg%np + ntnivmax ) )

     allocate (fp1(dg%np))
     allocate (fm1(dg%nm))
     allocate (Hy(dg%nd))
     allocate (Sy(dg%nd))
     allocate (Ssy(dg%nd))

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
        call contingence_cox(dataset,spt,ic,1,itest,est_moy)

!
!  la r�dustion du nombre d'effet estim�e est calcul�e
!
      nbreduit=supnbnivest-nbnivest

!
! Parametres de maximisation
      ibound=0
      npar=nbnivest
   print *,'npar max', npar
       do i=1,npar
        borni(i)=1.d-6
        borns(i)=1.d20
        par(i)=1.d0
      end do
!
! Point de depart (on reprend les point d'arriv�e pr�c�dents
!
      do i=1,npar
        par(i)=1.d0
      end do
!
! Optimisation de la vraisemblance a la position dx
        ifail=1

        call minimizing_funct(dataset,npar,ibound,funct_1qtl_modlin_cox,borni,borns,par,f1,iuser,user,ifail)
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
         if ( nbreduit == 0 ) then
            call log_mess("The effect ["//trim(dpm%namefix(dpm%modele(ic,3+ief)))//"]"//&
                          " might be confused with another effect !",WARNING_DEF)
            prob=0.d0
         else
           prob=MATH_QTLMAP_G01ECF('U',xlrt_t,dble(nbreduit),ifail)
         end if

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

     deallocate (fp1)
     deallocate (fm1)
     deallocate (Hy)
     deallocate (Sy)
     deallocate (Ssy)
     deallocate ( par  )
     deallocate ( borni )
     deallocate ( borns )

      write(nficout,*) "When this probability exceeds the standard threshold corresponding to the 5, 1 or 0.1 Pent level",&
               ", you might consider removing this effect from the model"

      return
      end subroutine test_lin_cox
!!***

!!****f* m_qtlmap_analyse_modlin_cox/set_solution_hypothesis0
!!  NAME
!!    set_solution_hypothesis0
!!  DESCRIPTION
!!
!!
!!  INPUTS
!!
!!  OUTPUTS
!!
!!
!!  SOURCE
       subroutine set_solution_hypothesis0(dataset,ic,incsol)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        integer                            ,intent(in)       :: ic
        type(TYPE_INCIDENCE_SOLUTION)      ,intent(inout)    :: incsol

       integer :: nteff,maxNbPar,jm,ifem,ip,ieff
       integer :: indest,isol,ipar,nlevel,ief,i,ife

       real(kind=dp) :: par_ref
       type(GENEALOGY_BASE) , pointer :: dg
       type(PHENOTYPE_BASE) , pointer :: dpa
       type(DATAMODEL_BASE) , pointer :: dpm

       dpm => dataset%phenoModel
       dg => dataset%genea
       dpa => dataset%phenoAnimal

       !OFI 09/2013 => pour le summary il faut le sig
       allocate (incsol%sig(1,dg%np))
       incsol%sig=1.0 !

       incsol%hypothesis=0
       !  Polygenic family
       nteff=0
       if (dg%np.ge.2) nteff = 1
       if ( count(dpa%estime(ic,:))  > 1 ) nteff = 2

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

       incsol%parameterVecsol=.true.
       incsol%eqtl_print=.true.

       ieff=0
       indest = 0
       isol=0

       if (dg%np.ge.2) then
       ieff=ieff+1
       par_ref=par0(indest+1)
       incsol%groupeName(ieff) = 'Sire polygenic effects'
       incsol%nbParameterGroup(ieff)=dg%np
          do ip=1,dg%np
             isol=isol+1
             if ( vecsol0(isol) ) then
               indest=indest+1
               incsol%paramaterValue(ieff,ip)  = (par0(indest)/par_ref)
             else
               incsol%paramaterValue(ieff,ip)  = 1.d0
             end if

             incsol%parameterName(ieff,ip)   ='Sire '//trim(dg%pere(ip))
             incsol%parameterVecsol(ieff,ip) = vecsol0(isol)
         end do
       endif
       if ( count(dpa%estime(ic,:)) > 1 ) then
           par_ref=par0(indest+1)
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
                incsol%paramaterValue(ieff,ifem) = (par0(indest)/par_ref)
              else
                incsol%paramaterValue(ieff,ifem) = 1.d0
              end if

              incsol%parameterName(ieff,ifem)='Dam '//trim(dg%mere(jm))//&
               " [Sire "//trim(dg%pere(ip))//"]"
              incsol%parameterVecsol(ieff,ifem) = vecsol0(isol)

             end if
           end do
          end do
         end if

       !Fixed effect
       if ( dpm%modele(ic,1) > 0 ) then
         ieff=ieff+1
         incsol%eqtl_print(ieff)=.false.
         incsol%groupeName(ieff) = 'fixed effects'
         incsol%nbParameterGroup(ieff)=nlevel
         ife=0
         do ief=1,nbef
          par_ref=par0(indest+1)
          do i=1,dpm%nlev(dpm%modele(ic,3+ief))
           ife=ife+1
           isol=isol+1

           if ( vecsol0(isol) ) then
             indest=indest+1
             incsol%paramaterValue(ieff,ife)  = (par0(indest)/par_ref)
           else
             incsol%paramaterValue(ieff,ife)  = 1.d0
           end if
           incsol%parameterName(ieff,ife)   = trim(dpm%namefix(dpm%modele(ic,3+ief)))&
             //' level '//trim(str(i))
           incsol%parameterVecsol(ieff,ife) = vecsol0(isol)

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
             incsol%paramaterValue(ieff,ief)  = (par0(indest))
           else
             incsol%paramaterValue(ieff,ief)  = 1.d0
           end if
           incsol%parameterName(ieff,ief)   = trim(dpm%namecov(dpm%modele(ic,3+dpm%modele(ic,1)+ief)))
           incsol%parameterVecsol(ieff,ief) = vecsol0(isol)
          end do
       end if

       end subroutine set_solution_hypothesis0
!!***

!!****f* m_qtlmap_analyse_modlin_cox/set_solution_hypothesis1
!!  NAME
!!    set_solution_hypothesis1
!!  DESCRIPTION
!!
!!
!!  INPUTS
!!
!!  OUTPUTS
!!
!!
!!  SOURCE
     subroutine set_solution_hypothesis1(dataset,ic,incsol)
       integer                            ,intent(in)       :: ic
       type(QTLMAP_DATASET)               ,intent(in)       :: dataset
       type(TYPE_INCIDENCE_SOLUTION)      ,intent(inout)    :: incsol

       integer :: nteff,maxNbPar,jm,ifem,ip,ieff
       integer :: ntlev,nbtp,jef,lp,indest
       integer :: isol,ipar,nlevel,ief,i,ife
       real(kind=dp) :: par_ref

       type(GENEALOGY_BASE) , pointer :: dg
       type(PHENOTYPE_BASE) , pointer :: dpa
       type(DATAMODEL_BASE) , pointer :: dpm

       dg => dataset%genea
       dpa => dataset%phenoAnimal
       dpm => dataset%phenoModel

       !OFI 09/2013 => pour le summary il faut le sig
       allocate (incsol%sig(1,dg%np))
       incsol%sig=1.0 !

       incsol%hypothesis=1
       !  Polygenic family, QTL effect
       nteff=1
       if (dg%np.ge.2) nteff = 2
       if ((dg%np==1).and. count(dpa%estime(ic,:))  ==1 ) nteff = 2
       if ((dg%np==1).and. count(dpa%estime(ic,:))  > 1 ) nteff = 3
       if ((dg%np.ge.2).and. count(dpa%estime(ic,:))==1 ) nteff = 3
       if ((dg%np.ge.2).and. count(dpa%estime(ic,:))> 1 ) nteff = 4
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
       nbtp = 3 + dpm%modele(ic,1)+dpm%modele(ic,2)!+modele(ic,3)
       do jef=1,dpm%modele(ic,3)
		      ntlev=ntlev*dpm%nlev(dpm%modele(ic,nbtp+jef))
	   end do

	   !max nombre de niveau pour un effet fixe en interaction avec le qtl ?
	   maxNbPar = max(maxNbPar,ntlev*dg%np)
	   maxNbPar = max(maxNbPar,ntlev*count(dpa%estime(ic,:)))

       allocate (incsol%groupeName(nteff))
       allocate (incsol%eqtl_print(nteff))
       allocate (incsol%nbParameterGroup(nteff))
       allocate (incsol%parameterName(nteff,maxNbPar))
       allocate (incsol%paramaterValue(nteff,maxNbPar))
       allocate (incsol%parameterVecsol(nteff,maxNbPar))
       allocate (incsol%qtl_groupeName(1,1))

       incsol%parameterVecsol=.true.
       incsol%eqtl_print=.true.

       ieff=1
       incsol%qtl_groupeName(1,1)=ieff
       incsol%groupeName(ieff) = 'Sire QTL effects'
       incsol%nbParameterGroup(ieff)=dg%np*ntlevp

       indest = 0
       isol=0
       ife=0
       do ip=1,dg%np
         do lp=1,ntlev
           ife=ife+1
           isol=isol+1
          ! print *, 'vecsol',vecsol1(1), 'par', par1(1)
           if ( vecsol1(isol) ) then
             indest=indest+1
             incsol%paramaterValue(ieff,ife)  = (par1(indest))
           else
             incsol%paramaterValue(ieff,ife)  = 1.d0
           end if

           incsol%parameterName(ieff,ife)   ='Sire '//trim(dg%pere(ip))
           incsol%parameterVecsol(ieff,ife) = vecsol1(isol)
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
                incsol%paramaterValue(ieff,ifem) = (par1(indest))
               else
                incsol%paramaterValue(ieff,ifem) = 1.d0
               end if
               incsol%parameterName(ieff,ifem)='Dam '//trim(dg%mere(jm))
               incsol%parameterVecsol(ieff,ifem) = vecsol1(isol)
              end do
             end if
           end do
         end if



      if (dg%np.ge.2) then
         ieff=ieff+1
         incsol%groupeName(ieff) = 'Sire polygenic effects'
         incsol%nbParameterGroup(ieff)=dg%np
         par_ref=par1(indest+1)
         do ip=1,dg%np
           isol=isol+1
           if ( vecsol1(isol) ) then
             indest=indest+1
             incsol%paramaterValue(ieff,ip)  = (par1(indest)/par_ref)
           else
             incsol%paramaterValue(ieff,ip)  = 1.d0
           end if

           incsol%parameterName(ieff,ip)   ='Sire '//trim(dg%pere(ip))
           incsol%parameterVecsol(ieff,ip) = vecsol1(isol)
          end do
      endif

       if ( count(dpa%estime(ic,:)) > 1 ) then
           ieff = ieff +1
           incsol%groupeName(ieff) = 'Dam polygenic effects'
           incsol%nbParameterGroup(ieff)=count(dpa%estime(ic,:))
           ifem=0
           par_ref=par1(indest+1)
          do ip=1,dg%np
           do jm=dg%nmp(ip)+1,dg%nmp(ip+1)
             if ( dpa%estime(ic,jm)) then
              isol=isol+1
              ifem=ifem+1
              if ( vecsol1(isol) ) then
                indest=indest+1
                incsol%paramaterValue(ieff,ifem) = (par1(indest)/par_ref)
              else
                incsol%paramaterValue(ieff,ifem) = 1.d0
              end if

              incsol%parameterName(ieff,ifem)='Dam '//trim(dg%mere(jm))//" [Sire "//trim(dg%pere(ip))//"]"
              incsol%parameterVecsol(ieff,ifem) = vecsol1(isol)
             end if
           end do
          end do
         end if

       !Fixed effect
       if ( dpm%modele(ic,1) > 0 ) then
         ieff=ieff+1
         incsol%eqtl_print(ieff)=.false.
         incsol%groupeName(ieff) = 'fixed effects'
         incsol%nbParameterGroup(ieff)=nlevel
         ife=0
         do ief=1,nbef
          par_ref=par1(indest+1)
          do i=1,dpm%nlev(dpm%modele(ic,3+ief))
           isol=isol+1
           ife=ife+1
           if ( vecsol1(isol) ) then
             indest=indest+1
             incsol%paramaterValue(ieff,ife)  = (par1(indest)/par_ref)
           else
             incsol%paramaterValue(ieff,ife)  = 1.d0
           end if
           incsol%parameterName(ieff,ife)   = trim(dpm%namefix(dpm%modele(ic,3+ief)))//' level '//trim(str(i))
           incsol%parameterVecsol(ieff,ife) = vecsol1(isol)
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
             incsol%paramaterValue(ieff,ief)  =(par1(indest))
           else
             incsol%paramaterValue(ieff,ief)  = 1.d0
           end if
           incsol%parameterName(ieff,ief)    = trim(dpm%namecov(dpm%modele(ic,3+dpm%modele(ic,1)+ief)))
           incsol%parameterVecsol(ieff,ief)  = vecsol1(isol)
          end do
       end if
       end subroutine set_solution_hypothesis1
!!***

end module m_qtlmap_analyse_modlin_cox
