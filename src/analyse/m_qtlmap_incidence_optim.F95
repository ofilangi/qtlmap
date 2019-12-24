!!****m* ANALYSE/m_qtlmap_incidence_optim
!!  NAME
!!    m_qtlmap_incidence_optim
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
module m_qtlmap_incidence_optim
    !use ieee_arithmetic
    use m_qtlmap_types
    use m_qtlmap_log
    use m_qtlmap_incidence
    implicit none

    save

    !Definition des constantes

!!****d* m_qtlmap_incidence_optim/grand
!! NAME
!!   grand
!! DESCRIPTION
!!   Constant used by the discrete data analysis optimization
!! NOTES
!!
!! SOURCE
    real (kind=dp)                         ,parameter   :: grand     = 1.d6
!!***

!!****v* m_qtlmap_incidence_optim/dloggrand
!! NAME
!!   dloggrand
!! DESCRIPTION
!!   Constant used by the discrete data analysis optimization (log of grand)
!! NOTES
!!
!!***
    real (kind=dp)                                      :: dloggrand
    !$omp threadprivate (dloggrand)

!!****v* m_qtlmap_incidence_optim/my_desc
!! NAME
!!   my_desc
!! DESCRIPTION
!!   pointer of the incidence description to acces into the likelihood function (objective function to optimized)
!! NOTES
!!  pointer of INCIDENCE_TYPE
!!***
    type(INCIDENCE_TYPE)     ,pointer    ,private       :: my_desc
!!****v* m_qtlmap_incidence_optim/my_xincreduit
!! NAME
!!   my_xincreduit
!! DESCRIPTION
!!   the reduction of the originale incidence matrix (without level which are not estimable according the cholesky decomposition)
!! DIMENSIONS
!!   nd, nbnivest
!! NOTES
!!  used by the likelihood function
!!***
    real (kind=dp),dimension(:,:),allocatable,private   :: my_xincreduit

!!****v* m_qtlmap_incidence_optim/current_chr
!! NAME
!!   current_chr
!! DESCRIPTION
!!   The currents chromosomes analysed (one for each qtl added in the contigence)
!! DIMENSIONS
!!   nd, nbnivest
!! NOTES
!!  used by the likelihood function
!!***
    integer , dimension(:)   ,pointer    ,private       :: current_chr

    !$omp threadprivate (My_DESC,my_xincreduit,current_chr)


  contains

!!****f* m_qtlmap_incidence_optim/likelihood_h0_family
!! NAME
!!   likelihood_h0_family
!! DESCRIPTION
!!   Calculus of the likelihood under H0 (without qtl effect) -- continue data
!! INPUTS
!!   ip  : the index of sire
!!   jm  : the index of dam
!!   n   : number of parameter
!!   x   : the values of parameters
!! iuser : integer array (defined by the user)
!! user  : real array (defined by the user)
!! OUTPUTS
!!   f   : the value of the likelihood
!! NOTES
!!
!! SOURCE
  subroutine likelihood_h0_family(ip,jm,n,x,f,iuser,user)
      integer         , intent(in)                  :: ip,jm,n
      real (kind=dp)      ,dimension(n), intent(in) :: x
      real (kind=dp)  , intent(inout)               :: f
      integer ,       dimension(1), intent(in)      :: iuser
      real (kind=dp)      ,dimension(1), intent(in) :: user

      real (kind=dp) :: var,sig,vmere,vpf
      integer        :: ifem,kd1,kd2
      real(kind=dp),dimension(my_desc%dataset%nkd) :: V

      f = 0.d0
      V = 0.d0
      ifem=jm-my_desc%dataset%data%genea%nmp(ip)

      kd2=my_desc%dataset%lSires(ip)%full_sib(ifem)%lastkd

      if (kd2<0) then
        return
      end if
      kd1=my_desc%dataset%lSires(ip)%full_sib(ifem)%firstkd

      !   V = matmul(my_xincreduit(:my_desc%dataset%nkd,:My_DESC%nbnivest),X(np+1:))
      call matmul_incidence(kd1,kd2,my_desc,my_desc%dataset%data%genea%nd,&
       my_desc%ntnivmax,my_xincreduit,X(my_desc%dataset%data%genea%np+1:),my_desc%dataset%nkd,kd1,V)
      V(kd1:kd2) = my_desc%dataset%Y(1,kd1:kd2)-V(kd1:kd2)
      V(kd1:kd2) = V(kd1:kd2)*V(kd1:kd2)*my_desc%dataset%CD(1,kd1:kd2)
      sig=x(ip)
      var=sig*sig
      vpf = sum(V(kd1:kd2))
      vmere=dexp(-0.5d0*vpf/var)
      if (vmere == 0) then
          f=INIFINY_REAL_VALUE
          return
      end if

      f=-dlog(vmere)+(kd2-kd1+1)*dlog(sig)

  end subroutine likelihood_h0_family
!!***

!!****f* m_qtlmap_incidence_optim/likelihood_hn_family
!! NAME
!!   likelihood_hn_family
!! DESCRIPTION
!!    Calculus of the likelihood under H(nqtl) (with nqtl qtl effect)  -- continue data
!! INPUTS
!!   ip  : the index of sire
!!   jm  : the index of dam
!!   n   : number of parameter
!!   x   : the values of parameters
!! iuser : integer array (defined by the user)
!! user  : real array (defined by the user)
!! OUTPUTS
!!   f   : the value of the likelihood
!! NOTES
!!
!! SOURCE
   subroutine likelihood_hn_family(ip,jm,n,x,f,iuser,user)
!      use OMP_LIB

      integer         , intent(in)                  :: ip,jm,n
      real (kind=dp)      ,dimension(n), intent(in) :: x
      real (kind=dp)  , intent(inout)               :: f
      integer ,       dimension(1), intent(inout)      :: iuser
      real (kind=dp)      ,dimension(1), intent(inout) :: user

      real (kind=dp) :: var,sig,vmere,vpf,pbr,effm
      integer        :: jjm,ig(my_desc%nqtl),ifem,indf,indm,iq,ngg,iig,z
      real(kind=dp),dimension(my_desc%dataset%nkd) :: V
      integer        :: kd1,kd2,jj,adp,adm
      logical        :: ok

      f=0.d0
      sig=x(ip)
      var=sig*sig

      !scale where the qtl is biallelic or not...
      if ( my_desc%isBiAllellic ) then
        adp=0
      else
        adp=ip-1
      end if

    !  print *,'ip:',ip,' jm:',jm,' nmp(ip):',nmp(ip),' nmp(ip+1)',nmp(ip+1)
      jjm=jm-my_desc%dataset%data%genea%nmp(ip)
!      print *,'jjm:',jjm
      if (my_desc%dataset%lSires(ip)%full_sib(jjm)%lastKd>0) then
          kd1=my_desc%dataset%lSires(ip)%full_sib(jjm)%firstKd
          kd2=my_desc%dataset%lSires(ip)%full_sib(jjm)%lastKd
          effm=dble(kd2-kd1+1)
          vmere=0.d0

          if (my_desc%dataset%data%phenoAnimal%estime(my_desc%ic,jm)) then
            ifem=my_desc%dataset%data%phenoAnimal%iam(my_desc%ic,my_desc%dataset%data%genea%repfem(jm))
            if ( my_desc%isBiAllellic ) then
             adm=0
            else
             adm=ifem-1
            end if
          end if
          !ngg : le nombre de genotype possible sur tous les qtls...
          ngg=1
          do iq=1,my_desc%nqtl
             ig(iq)=my_desc%dataset%spt%ngenom(current_chr(iq),jm)+1
             ngg=ngg*(my_desc%dataset%spt%ngenom(current_chr(iq),jm+1)-my_desc%dataset%spt%ngenom(current_chr(iq),jm))
          end do
          !pour toutes les combinaisons possibles des genotypes
          do iig=1,ngg
             pbr=1
             !on modifie la matrice d incidence pour les n qtl
             do iq=1,my_desc%nqtl
              indm=my_desc%ntniv_qtlsires(iq)+adp
              if ( my_desc%dataset%data%phenoAnimal%estime(my_desc%ic,jm) ) then
               !Si la mere est estimable, on place dans la matrice les
               !pdds dam et pdds sires (si ces effets sont estimables)
               indf=my_desc%ntniv_qtldams(iq)+adm

               if ( My_desc%vecsol(indf)) then
                 indf=my_desc%corr_niv_nivb(indf) ! vrai correspondance dans le tableau si estimable
                 my_xincreduit(kd1:kd2,indf)=my_desc%dataset%lSires(ip)%full_sib(jjm)%pmt(iq,(ig(iq)-&
                  my_desc%dataset%spt%ngenom(current_chr(iq),jm)),:)
               end if
               if ( My_desc%vecsol(indm) ) then
                indm=my_desc%corr_niv_nivb(indm) ! vrai correspondance dans le tableau si estimable
                if ( indf == indm ) then ! biallelic case
                  my_xincreduit(kd1:kd2,indm)=my_xincreduit(kd1:kd2,indm) + &
                     my_desc%dataset%lSires(ip)%full_sib(jjm)%ppt(iq,(ig(iq)-&
                      my_desc%dataset%spt%ngenom(current_chr(iq),jm)),:)
                else
                  my_xincreduit(kd1:kd2,indm)=my_desc%dataset%lSires(ip)%full_sib(jjm)%ppt(iq,(ig(iq)-&
                   my_desc%dataset%spt%ngenom(current_chr(iq),jm)),:)
                end if
               end if
              else
               !la mere n est pas estimable, on place seulement les pdds males
               if ( My_desc%vecsol(indm) ) then
                indm=my_desc%corr_niv_nivb(indm) ! vrai correspondance dans le tableau si estimable
                my_xincreduit(kd1:kd2,indm)=my_desc%dataset%lSires(ip)%ppt(iq,kd1:kd2)
               end if
              end if
              !print *,ig(iq)-ngenom(current_chr(iq),jm),ig(iq)
              pbr=pbr*my_desc%dataset%spt%probg(current_chr(iq),ig(iq))
             end do ! iq

            !V(kd1:kd2) = matmul(my_xincreduit(kd1:kd2,:My_DESC%nbnivest),X(np+1:))
            call matmul_incidence(kd1,kd2,my_desc,my_desc%dataset%data%genea%nd,&
             my_desc%ntnivmax,my_xincreduit,X(my_desc%dataset%data%genea%np+1:),my_desc%dataset%nkd,kd1,V)
            V(kd1:kd2) = my_desc%dataset%Y(1,kd1:kd2) - V(kd1:kd2)
            V(kd1:kd2) = V(kd1:kd2)*V(kd1:kd2)*my_desc%dataset%CD(1,kd1:kd2)
            vpf = sum(V(kd1:kd2))

            vmere=vmere+pbr*dexp(-0.5d0*vpf/var)
            ! on increment
            ok=.true.
            do iq=1,my_desc%nqtl
             if (ok) then
               if ((ig(iq) < my_desc%dataset%spt%ngenom(current_chr(iq),jm+1))) then
                 ig(iq)=ig(iq)+1
                 ok=.false.
               end if
             end if
            end do
         end do ! iig

          if (vmere == 0) then
             f=INIFINY_REAL_VALUE
          else
             f=-dlog(vmere)+effm*dlog(sig)
          end if
         end if

  end subroutine likelihood_hn_family
!!***

!!****f* m_qtlmap_incidence_optim/likelihood_discret_h0_family
!! NAME
!!   likelihood_discret_h0_family
!! DESCRIPTION
!!   Calculus of the likelihood under H0 (without qtl effect) -- discrete data
!! INPUTS
!!   ip  : the index of sire
!!   jm  : the index of dam
!!   n   : number of parameter
!!   x   : the values of parameters
!! iuser : integer array (defined by the user)
!! user  : real array (defined by the user)
!! OUTPUTS
!!   f   : the value of the likelihood
!! NOTES
!!
!! SOURCE
     subroutine likelihood_discret_h0_family(ip,jm,n,x,f,iuser,user)

      implicit none
      integer         , intent(in)                       :: ip,jm,n
!
! Tableaux dimensionnes selon n, le nombre de parametres a estimer
!
      real (kind=dp)   ,dimension(n), intent(in)         :: x
      real (kind=dp)   ,intent(inout)                    :: f
      integer          ,dimension(1), intent(in)         :: iuser
      real (kind=dp)   ,dimension(1), intent(in)         :: user

! Tableaux dimensionnes selon nm, le nombre de meres
      real (kind=dp)   ,dimension(my_desc%dataset%data%phenoModel%nmod(my_desc%ic))      :: threshold

!modif mkw
! declaration de l comme entier et on a plus besoin de la variable vpf

      integer :: i, ifail
      integer :: neffet,ief,ifem,kd,kd1,kd2
      real (kind=dp) :: vpf,tig
      real(kind=dp),dimension(my_desc%dataset%nkd) :: V
!mkw
!
   neffet=my_desc%dataset%data%genea%np+my_desc%nbnivest

!******************************************************************************
!
    ifail=1

    threshold(1)=x(neffet+1)
    do i=2,my_desc%dataset%data%phenoModel%nmod(my_desc%ic)-1
      threshold(i)=threshold(i-1)+x(neffet+i)
    end do
    f=0
    ifem=jm-my_desc%dataset%data%genea%nmp(ip)
    !V = matmul(my_xincreduit(:my_desc%dataset%nkd,:My_DESC%nbnivest),X(np+1:neffet))
    kd1=my_desc%dataset%lSires(ip)%full_sib(ifem)%firstkd
    kd2=my_desc%dataset%lSires(ip)%full_sib(ifem)%lastkd

    if (kd2<0) then
      return
    end if

    call matmul_incidence(kd1,kd2,my_desc,my_desc%dataset%data%genea%nd,&
     my_desc%ntnivmax,my_xincreduit,X(my_desc%dataset%data%genea%np+1:neffet),my_desc%dataset%nkd,kd1,V)

    tig=grand
    if (x(ip)> 0)tig=1.d0/x(ip)

    if (my_desc%dataset%lSires(ip)%full_sib(ifem)%lastKd<0) then
      return
    end if

    do kd=kd1,kd2
         if (my_desc%dataset%ydiscretord(1,kd)==1)  then
	             vpf=MATH_QTLMAP_G01EAF('L',((threshold(1)-v(kd))*tig),ifail)
          else
          if (my_desc%dataset%ydiscretord(1,kd)==my_desc%dataset%data%phenoModel%nmod(my_desc%ic)) then
              vpf=1-MATH_QTLMAP_G01EAF('L',(threshold(my_desc%dataset%data%phenoModel%nmod(my_desc%ic)-1)-v(kd))*tig,ifail)
           else
        	  vpf=MATH_QTLMAP_G01EAF('L',(threshold(my_desc%dataset%ydiscretord(1,kd))-v(kd))*tig,ifail)-&
        	    MATH_QTLMAP_G01EAF('L',(threshold(my_desc%dataset%ydiscretord(1,kd)-1)-v(kd))*tig,ifail)
          end if
	     end if

  	     if (vpf <= 0) then
             f=INIFINY_REAL_VALUE
         else
             f=f-dlog(vpf)
	     end if
    end do ! kd

    end subroutine likelihood_discret_h0_family
!!***

!!****f* m_qtlmap_incidence_optim/likelihood_discret_hn_family
!! NAME
!!   likelihood_discret_hn_family
!! DESCRIPTION
!!    Calculus of the likelihood under H(nqtl) (with nqtl qtl effect) -- discrete data
!! INPUTS
!!   ip  : the index of sire
!!   jm  : the index of dam
!!   n   : number of parameter
!!   x   : the values of parameters
!! iuser : integer array (defined by the user)
!! user  : real array (defined by the user)
!! OUTPUTS
!!   f   : the value of the likelihood
!! NOTES
!!
!! SOURCE
     subroutine likelihood_discret_hn_family(ip,jm,n,x,f,iuser,user)

      implicit none

      integer         , intent(in)                  :: n,ip,jm
!
! Tableaux dimensionnes selon n, le nombre de parametres a estimer
!
      real (kind=dp)      ,dimension(n), intent(in) :: x
      real (kind=dp)  , intent(inout)               :: f
      integer ,       dimension(1), intent(in)      :: iuser
      real (kind=dp)      ,dimension(1), intent(in) :: user
      real (kind=dp)   ,dimension(my_desc%dataset%data%phenoModel%nmod(my_desc%ic))      :: threshold

! Tableaux dimensionnes selon nm, le nombre de meres
      real (kind=dp)   ,dimension(my_desc%dataset%data%genea%nm)   :: effm


!
! Divers
      integer :: neffet,kd,kkd,ilev,lambda
      integer :: ntnivmax,neffmax, m, i, j, m1, temp, k, ii, ifail


      real (kind=dp) :: sig,var,vmere,vpf,wpf,tig
      real (kind=dp) :: pbr

      integer        :: jjm,ig(my_desc%nqtl),ifem,indf,indm,iq,ngg,iig
      real(kind=dp),dimension(my_desc%dataset%nkd) :: V
      integer        :: kd1,kd2,jj,adp,adm
      logical        :: ok

!******************************************************************************
!
      neffet=my_desc%dataset%data%genea%np+my_desc%nbnivest
      ifail=1

      !scale where the qtl is biallelic or not...
      if ( my_desc%isBiAllellic ) then
        adp=0
      else
        adp=ip-1
      end if

      threshold(1)=x(neffet+1)
      do i=2,my_desc%dataset%data%phenoModel%nmod(my_desc%ic)-1
        threshold(i)=threshold(i-1)+x(neffet+i)
      end do
      f=0
      ok=.false.

      tig=grand
      if (x(ip)> 0)tig=1.d0/x(ip)
      jjm=jm-my_desc%dataset%data%genea%nmp(ip)

      if (my_desc%dataset%lSires(ip)%full_sib(jjm)%lastKd>0) then
          vmere=0.d0
           if (my_desc%dataset%data%phenoAnimal%estime(my_desc%ic,jm)) &
            ifem=my_desc%dataset%data%phenoAnimal%iam(my_desc%ic,my_desc%dataset%data%genea%repfem(jm))
           kd1=my_desc%dataset%lSires(ip)%full_sib(jjm)%firstKd
           kd2=my_desc%dataset%lSires(ip)%full_sib(jjm)%lastKd
         !  igg=0
           !ngg : le nombre de genotype possible sur tous les qtls...
           ngg=1
           do iq=1,my_desc%nqtl
             ig(iq)=my_desc%dataset%spt%ngenom(current_chr(iq),jm)+1
             ngg=ngg*(my_desc%dataset%spt%ngenom(current_chr(iq),jm+1)-my_desc%dataset%spt%ngenom(current_chr(iq),jm))
           end do
           do iig=1,ngg
         !  igg=0
           pbr=1
           !on modifie la matrice d incidence pour les n qtl
           do iq=1,my_desc%nqtl
            indm=my_desc%ntniv_qtlsires(iq)+adp
            if ( my_desc%dataset%data%phenoAnimal%estime(my_desc%ic,jm) ) then
               if ( my_desc%isBiAllellic ) then
                 adm=0
               else
                 adm=ifem-1
               end if
              !Si la mere est estimable, on place dans la matrice les
              !pdds dam et pdds sires (si ces effets sont estimables)
              indf=my_desc%ntniv_qtldams(iq)+adm
              if ( My_desc%vecsol(indf)) then
               indf=my_desc%corr_niv_nivb(indf) ! vrai correspondance dans le tableau si estimable
               my_xincreduit(kd1:kd2,indf)=my_desc%dataset%lSires(ip)%full_sib(jjm)%pmt(iq,ig(iq)-&
                my_desc%dataset%spt%ngenom(current_chr(iq),jm),:)
              end if
              if ( My_desc%vecsol(indm) ) then
               indm=my_desc%corr_niv_nivb(indm) ! vrai correspondance dans le tableau si estimable
               if ( indm == indf ) then !biallelic case
                my_xincreduit(kd1:kd2,indm)=my_xincreduit(kd1:kd2,indm)+&
                  my_desc%dataset%lSires(ip)%full_sib(jjm)%ppt(iq,ig(iq)-&
                   my_desc%dataset%spt%ngenom(current_chr(iq),jm),:)
               else
                my_xincreduit(kd1:kd2,indm)=my_desc%dataset%lSires(ip)%full_sib(jjm)%ppt(iq,ig(iq)-&
                 my_desc%dataset%spt%ngenom(current_chr(iq),jm),:)
               end if
              end if
            else
               !la mere n est pas estimable, on place seulement les pdds males
              if ( My_desc%vecsol(indm) ) then
               indm=my_desc%corr_niv_nivb(indm) ! vrai correspondance dans le tableau si estimable
               my_xincreduit(kd1:kd2,indm)=my_desc%dataset%lSires(ip)%ppt(iq,kd1:kd2)
              end if
            end if
           pbr=pbr*my_desc%dataset%spt%probg(current_chr(iq),ig(iq))
          end do ! iq

           !V(kd1:kd2) = matmul(my_xincreduit(kd1:kd2,:My_DESC%nbnivest),X(np+1:neffet))
           call matmul_incidence(kd1,kd2,my_desc,my_desc%dataset%data%genea%nd,my_desc%ntnivmax,&
            my_xincreduit,X(my_desc%dataset%data%genea%np+1:neffet),my_desc%dataset%nkd,kd1,V)

           vpf=grand
           do kkd=kd1,kd2
             if (my_desc%dataset%ydiscretord(1,kkd)==1) then
    	  	   wpf=MATH_QTLMAP_G01EAF('L',((threshold(1)-v(kkd))*tig),ifail)
    	     else
    	  	   if (my_desc%dataset%ydiscretord(1,kkd)==my_desc%dataset%data%phenoModel%nmod(my_desc%ic)) then
    	  	      wpf=1.d0-MATH_QTLMAP_G01EAF('L',(threshold(my_desc%dataset%data%phenoModel%nmod(my_desc%ic)-1)-v(kkd))*tig,ifail)
    	  	   else
    	  	     wpf=MATH_QTLMAP_G01EAF('L',(threshold(my_desc%dataset%ydiscretord(1,kkd))-v(kkd))*tig,ifail)-&
    	  		 MATH_QTLMAP_G01EAF('L',(threshold(my_desc%dataset%ydiscretord(1,kkd)-1)-v(kkd))*tig,ifail)
    	  	   end if
    	     end if
             vpf=vpf*wpf
    	   end do

            vmere=vmere+pbr*vpf
            ! on increment
            do iq=1,my_desc%nqtl
             ! print *,"test:",ig(iq),ngenom(current_chr(iq),jm),iig,ngg
             if ((ig(iq) < my_desc%dataset%spt%ngenom(current_chr(iq),jm+1))) then
                 ig(iq)=ig(iq)+1
                 exit
             end if
            end do
            end do !ig

     	    if (vmere <= 0) then
     		     f=INIFINY_REAL_VALUE
     		     return
            endif

            f=dloggrand-dlog(vmere)
       end if

      return
      end subroutine likelihood_discret_hn_family
!!***

!!****f* m_qtlmap_incidence_optim/model_optim_family
!!  NAME
!!    model_optim_family
!!  DESCRIPTION
!!    Generic model for an optimization of a likelihood (discrete and continue data)
!!    process :
!!    * Cholesky decomposition (to have the estimability of each level of the contingence)
!!    * build a matrix without no-null element to improve the multiplication matrix during the likelihood calculus
!!    * compute threhold in discrete data case
!!    * build a filter to improve the optimization (see m_qtlmap_optimization)
!!
!!  INPUTS
!!    xinc            : the contigence matrix
!!   incidenceDesc    : description of the contingence and incidence matrix
!!   performPrecision : indicates if the precision have to be computed
!!   FUNCT_PART       : the likelihood function to optimize
!!   tConf            : indicates if the caller want to keep a buffer array in order to compute the confusion between parameter
!!
!!  OUTPUTS
!!     f              : the maximum value of the likelihood
!!   sigsquareEstime  : residual variance
!!   Bestim           : solution values of each level of the contingence matrix
!!   tempForConfusion : a buffer array for the confusion calculus
!!
!!  HISTORY
!!    10/09/2010      * correction , my_desc utilise pour la dimension d une variable locale (my_desc pas initialise a ce stade)
!!  NOTES
!!  SOURCE
    subroutine model_optim_family(xinc,incidenceDesc,osig,sigsquareEstime,Bestim,&
    performPrecision,FUNCT_PART,tConf,tempForConfusion)
         type(INCIDENCE_TYPE),target           , intent(inout)     :: incidenceDesc
  real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax), intent(in) :: xinc
!         type(INCIDENCE_GEN_STRUCT)              ,intent(in)       :: workstruct
         real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%np)        ,intent(out)        :: osig
         real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%np)        ,intent(out) :: sigsquareEstime
         real (kind=dp) , dimension(incidenceDesc%ntnivmax)  ,intent(out)        :: Bestim
         logical                               ,intent(in)         :: performPrecision,tConf
         real (kind=dp),dimension(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax) ,intent(out)  :: tempForConfusion
         external                                                  :: FUNCT_PART

         integer :: j,i,ip,ifail,ibound,npar
!         real(kind=dp) ,dimension(np+ntnivmax) :: par
         real (kind=dp) , dimension(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax)   :: XX,xxx
         real (kind=dp) , dimension(incidenceDesc%ntniv,incidenceDesc%ntniv)   :: triang
         integer                                    :: iuser(1),ix,kd1,kd2,jm,jjm
         real (kind=dp)                             :: user(1),f
         logical       ,dimension(incidenceDesc%ntnivmax)         :: lastvecsol
         logical  , dimension(:,:,:),pointer        :: filter_inc
   !      real (kind=dp)       ,dimension(incidenceDesc%dataset%data%genea%np)        :: fp
   !      real (kind=dp)       ,dimension(incidenceDesc%dataset%data%genea%nm)        :: fm

         my_desc=>incidenceDesc
         lastvecsol=my_desc%vecsol
         ! create X'.X matrix from incidence matrix
         call model_XT_X(xinc,my_desc,XX)
         ! Check all parameters to remove from the estimation
         call estim_cholesky(XX,my_desc,my_desc%ntniv,triang)
         ! compute the precision of each parameter
         if (performPrecision) call get_precision(XX,xxx,my_desc)
         if (tConf) tempForConfusion=xxx
         allocate (my_xincreduit(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax))
         my_xincreduit=0.d0
         call set_corrxinc(xinc,my_desc,my_xincreduit)
         !optimisation : on sauvegarde les index des elements non nul de la matrice reduite
         call fill_nonull_elements(my_desc,size(my_xincreduit,1),size(my_xincreduit,2),my_xincreduit)
         ! Optimisation de la vraisemblance a la position dx
         ifail=1
         ibound=0
         npar=incidenceDesc%dataset%data%genea%np+my_desc%nbnivest

         if (count(lastvecsol)>0) then ! autre que initialisation
           j=incidenceDesc%dataset%data%genea%np
           do i=1,my_desc%ntniv
            if(my_desc%vecsol(i))then
             j=j+1
             if(.not. lastvecsol(i)) incidenceDesc%par(j)=0.d0
            end if
           end do
         else
         ! print *,'initialisation......'
          incidenceDesc%par=0.d0
         do ip=1,incidenceDesc%dataset%data%genea%np
          incidenceDesc%par(ip)=incidenceDesc%dataset%lSires(ip)%sig0(1)
         end do
        end if

          if ( incidenceDesc%dataset%data%phenoModel%nmod(incidenceDesc%ic) > 1) then
          npar=npar + incidenceDesc%dataset%data%phenoModel%nmod(incidenceDesc%ic)-1
          incidenceDesc%par(incidenceDesc%dataset%data%genea%np+my_desc%nbnivest+1)=&
           incidenceDesc%dataset%data%phenoModel%seuil(incidenceDesc%ic,1)
          do ix=2,incidenceDesc%dataset%data%phenoModel%nmod(incidenceDesc%ic)-1
           incidenceDesc%par(incidenceDesc%dataset%data%genea%np+my_desc%nbnivest+ix)=&
            incidenceDesc%dataset%data%phenoModel%seuil(incidenceDesc%ic,ix)-&
            incidenceDesc%dataset%data%phenoModel%seuil(incidenceDesc%ic,ix-1)
          end do
         end if
         !call debug_write_incidence(xinc,incidenceDesc)
              !Preparation de l optimisation....
         !call init_optim(incidenceDesc%ic)

         allocate (filter_inc(incidenceDesc%dataset%data%genea%np,&
                              incidenceDesc%dataset%data%genea%nm,&
                              npar))
         filter_inc=.true.
         filter_inc(:,:,1:incidenceDesc%dataset%data%genea%np+1)=.false.
         do ip=1,incidenceDesc%dataset%data%genea%np
           filter_inc(ip,incidenceDesc%dataset%data%genea%nmp(ip)+1:incidenceDesc%dataset%data%genea%nmp(ip+1),ip)=.true.
           jjm=0
           do jm=incidenceDesc%dataset%data%genea%nmp(ip)+1,incidenceDesc%dataset%data%genea%nmp(ip+1)
             jjm=jjm+1
             if (my_desc%dataset%lSires(ip)%full_sib(jjm)%lastKd>0) then
             kd1=my_desc%dataset%lSires(ip)%full_sib(jjm)%firstkd
             kd2=my_desc%dataset%lSires(ip)%full_sib(jjm)%lastkd
             filter_inc(ip,jm,incidenceDesc%dataset%data%genea%np+1:incidenceDesc%dataset%data%genea%np+my_desc%nbnivest)=&
              any(my_xincreduit(kd1:kd2,:my_desc%nbnivest)/=0.d0,dim=1)
             else
          !     print *,'pas de desc....'
               filter_inc(ip,jm,incidenceDesc%dataset%data%genea%np+1:incidenceDesc%dataset%data%genea%np+my_desc%nbnivest)=.false.
             end if
           end do
         end do
     !    filter_inc=.true.
         call minimizing_funct_family(incidenceDesc%dataset%data,npar,ibound,&
          FUNCT_PART,filter_inc,incidenceDesc%fmeremax,incidenceDesc%fperemax,&
          incidenceDesc%borni,incidenceDesc%borns,&
         incidenceDesc%par,f,iuser,user,ifail)

         deallocate (filter_inc)

         !getting standart deviation
         do ip=1,incidenceDesc%dataset%data%genea%np
           osig(ip)=incidenceDesc%par(ip)
           sigsquareEstime(ip)=osig(ip)*osig(ip)
         end do
         !The solution
         Bestim(:my_desc%nbnivest)=&
          incidenceDesc%par(incidenceDesc%dataset%data%genea%np+1:incidenceDesc%dataset%data%genea%np+my_desc%nbnivest)

         deallocate (my_xincreduit)
!
      end subroutine model_optim_family
!!***

!!****f* m_qtlmap_incidence_optim/model_optim_h0
!! NAME
!!   model_optim_h0
!! DESCRIPTION
!!
!!  INPUTS
!!    xinc            : the contigence matrix
!!   incidenceDesc    : description of the contingence and incidence matrix
!!   performPrecision : indicates if the precision of each parameter have to be computed
!!
!!  INPUTS/OUTPUTS
!!   workstruct       : result information about current point (not available in this function) and the specific method to used (linear or not, homo,heteroscedastic)
!!
!!  OUTPUTS
!!   sigsquareEstime  : residual variance
!!   Bestim           : solution values of each level of the contingence matrix
!!
!! NOTES
!!
!! SOURCE
      subroutine model_optim_h0(xinc,incidenceDesc,workstruct,sigsquareEstime,Bestim,performPrecision)
         type(INCIDENCE_TYPE)                                 , intent(inout)    :: incidenceDesc
         real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax), intent(in) :: xinc
         type(INCIDENCE_GEN_STRUCT)                           , intent(in)       :: workstruct
         real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%np)                       , intent(out):: sigsquareEstime
         real (kind=dp) , dimension(incidenceDesc%ntnivmax)   , intent(out)      :: Bestim
         logical                                              , intent(in)        :: performPrecision

         real (kind=dp) ,dimension(incidenceDesc%dataset%data%genea%np) :: osig
         real(kind=dp) , dimension(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax) :: xxx

         integer :: ip

         if ( incidenceDesc%dataset%data%phenoModel%natureY(incidenceDesc%ic) == 'r') then
            call model_optim_family(xinc,incidenceDesc,&
            osig,sigsquareEstime,Bestim,performPrecision,likelihood_h0_family,.false.,xxx)

         else if ( incidenceDesc%dataset%data%phenoModel%natureY(incidenceDesc%ic) == 'i' ) then
            dloggrand= dlog(grand)
            call model_optim_family(xinc,incidenceDesc,osig,&
             sigsquareEstime,Bestim,performPrecision,likelihood_discret_h0_family,.false.,xxx)
         end if

      end subroutine model_optim_h0
!!***

!!****f* m_qtlmap_incidence_optim/model_optim_hn
!! NAME
!!   model_optim_hn
!! DESCRIPTION
!!
!!  INPUTS
!!    xinc            : the contigence matrix
!!   incidenceDesc    : description of the contingence and incidence matrix
!!   curPos           : information about the current position
!!   performPrecision : indicates if the precision of each parameter have to be computed
!!   tConf            : indicates if the caller want to keep a buffer array in order to compute the confusion between parameter
!!   invlrt           : true => 2* (F0max - Fnqtlmax), false =>   2* ( Fnqtlmax - F0max )
!!
!!  INPUTS/OUTPUTS
!!   workstruct       : result information about current point and the specific method to used (linear or not, homo,heteroscedastic)
!!  OUTPUTS
!!   sigsquareEstime  : residual variance
!!   Bestim           : solution values of each level of the contingence matrix
!!   tempForConfusion : a buffer array for the confusion calculus
!!
!! NOTES
!!
!! SOURCE
      subroutine model_optim_hn(xinc,incidenceDesc,curPos,workstruct,sigsquareEstime,Bestim,&
                                performPrecision,tConf,tempForConfusion,invlrt)
         type(INCIDENCE_TYPE)                  , intent(inout)     :: incidenceDesc
         real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax), intent(in)       :: xinc
         type(POSITION_LRT_INCIDENCE)            ,intent(inout)    :: curPos
         type(INCIDENCE_GEN_STRUCT)              ,intent(inout)    :: workstruct
         real (kind=dp) , dimension(incidenceDesc%ntnivmax)  ,intent(out)        :: Bestim
         real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%np)        ,intent(out)        :: sigsquareEstime
         logical                               ,intent(in)         :: performPrecision,tConf
         real (kind=dp),dimension(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax) ,intent(out)  :: tempForConfusion
         logical                                     ,intent(in)   :: invlrt

         real (kind=dp) :: osig(incidenceDesc%dataset%data%genea%np),t1
         integer :: hypothesis,i
         character(len=LEN_L) :: dx

         allocate (current_chr(workstruct%nqtl))
         current_chr=curPos%listChr

         incidenceDesc%par(:incidenceDesc%dataset%data%genea%np)=sqrt(workstruct%sigsquare(workstruct%hypothesis-1,:,1))

         if ( incidenceDesc%dataset%data%phenoModel%natureY(incidenceDesc%ic) == 'r') then
            call model_optim_family(xinc,incidenceDesc,osig,sigsquareEstime,&
             Bestim,performPrecision,likelihood_hn_family,tConf,tempForConfusion)
         else if ( incidenceDesc%dataset%data%phenoModel%natureY(incidenceDesc%ic) == 'i') then
            dloggrand= dlog(grand)
            call model_optim_family(xinc,incidenceDesc,osig,&
            sigsquareEstime,Bestim,performPrecision,likelihood_discret_hn_family,tConf,tempForConfusion)

         end if

         if ((any(incidenceDesc%fperemax == INIFINY_REAL_VALUE)) ) then
            dx=""
            do i=1,workstruct%nqtl-1
              dx=trim(dx)//trim(str(incidenceDesc%dataset%data%map%absi(curPos%listChr(i),curPos%listN(i))))//","
            end do

            dx=trim(dx)//str(incidenceDesc%dataset%data%map%absi(curPos%listChr(workstruct%nqtl),curPos%listN(workstruct%nqtl)))
            call log_mess("dx ["//trim(dx)// "]. Can not optimize likelihood....The start point is reinitializing",WARNING_DEF)
            incidenceDesc%par=0.d0
            curPos%lrtSires=0.d0
            curPos%lrtDams=0.d0
            incidenceDesc%fperemax=0.d0
            incidenceDesc%fmeremax=0.d0
            return
         end if

         !compute LRT
         if (invLrt) then
            do hypothesis=1,workstruct%hypothesis
              curPos%lrtSires(hypothesis,:)=-2.d0*(workstruct%fnqtlsires(hypothesis-1,:)-incidenceDesc%fperemax)
              curPos%lrtDams(hypothesis,:)=-2.d0*(workstruct%fnqtldams(hypothesis-1,:)-incidenceDesc%fmeremax)
            end do
         else
           do hypothesis=1,workstruct%hypothesis
            curPos%lrtSires(hypothesis,:)=-2.d0*(incidenceDesc%fperemax-workstruct%fnqtlsires(hypothesis-1,:))
            curPos%lrtDams(hypothesis,:)=-2.d0*(incidenceDesc%fmeremax-workstruct%fnqtldams(hypothesis-1,:))
           end do
         end if


         deallocate (current_chr)

      end subroutine model_optim_hn
!!***

end module m_qtlmap_incidence_optim
