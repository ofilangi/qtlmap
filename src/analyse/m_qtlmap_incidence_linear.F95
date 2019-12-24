!!****m* ANALYSE/m_qtlmap_incidence_linear
!!  NAME
!!    m_qtlmap_incidence_linear
!!  SYNOPSIS
!!    Model Homoscedastic and Heteroscedastic Linear analysis
!!    Implementation of a likelihood ratio test linearised for Homoscedastic and Heteroscedastic Linear
!!
!!  DESCRIPTION
!!
!!    i   : sire index
!!    j   : Dam index
!!    k   : Progeny index
!!    N i : Number of progeny of Sire i
!!
!!     ~
!!    SIG i : Sqrt Variance of Sire family i under H0
!!     ^
!!    SIG i : Sqrt Variance of Sire family i under H1 at a position dx
!!
!!    General Case (Heteroscedastic model) :
!!    -----------------------------------
!!
!!          |   B = (X'.V-1.X)-1 . X' . V-1.Y
!!   System |    ^
!!          |   SIG² i = SUM [ Y ijk  - X ijk . B]² CD² ijk  / Ni
!!
!!    becomes in homoscedastic model :
!!    ------------------------------
!!
!!          |  B = (X'.X)-1 . X'.Y
!!   System |   ^
!!          |  SIG² = SUM [ Y ijk  - X ijk . B]² CD² ijk  / Ni
!!                  ~           ^                       ~           ^
!!     LRT = -2*Log(VL) + 2*Log(VL) = SUM i[ Ni * (Log(SIG²) - Log(SIG²))  ]
!!
!!  NOTES
!!
!!  BUGS
!!
!!  SEE ALSO
!!
!!  COPYRIGHT
!!***
!! You can use this space for remarks that should not be included
!! in the documentation.
!!/


module m_qtlmap_incidence_linear
    use m_qtlmap_types
    use m_qtlmap_log
    use m_qtlmap_solve_linear
    use m_qtlmap_optimization
    use m_qtlmap_output_handler
    use m_qtlmap_incidence

    implicit none

    contains

!!****f* m_qtlmap_incidence_linear/modele_homoscedastic
!! NAME
!!    modele_homoscedastic
!! DESCRIPTION
!!  From a incidence matrix, get the B estimation and the variance of each sires in a homoscedastic model
!!
!!   Bs   = (X'.X)-1 . X'.Y
!!   SIGs = SUM   [ Y    - X   . Bs]² / N
!!      i      jk    ijk    ijk          i
!! NOTES
!! SOURCE
     subroutine modele_homoscedastic(xinc,incidenceDesc,osigsquare,Bestim,performPrecision,tempForConfusion)
         type(INCIDENCE_TYPE)                  , intent(inout)     :: incidenceDesc
         real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax), intent(in) :: xinc
         real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%np)        ,intent(out)   :: osigsquare
         real (kind=dp) , dimension(incidenceDesc%ntnivmax)  ,intent(out)        :: Bestim
         logical                               ,intent(in)         :: performPrecision
         real (kind=dp),dimension(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax) ,intent(out)  :: tempForConfusion

         !local
         real (kind=dp) , dimension(incidenceDesc%ntnivmax)           :: RHS,YTemp
         real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax) :: xincreduit
         real (kind=dp) , dimension(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax)  :: XX!,Mr
         real (kind=dp) , dimension(incidenceDesc%ntniv,incidenceDesc%ntniv)  :: triang
         integer                                                              :: i,j
         real (kind=dp)                                                       :: f,XB(incidenceDesc%dataset%nkd),sigsq

         ! create X'.X matrix from incidence matrix
         call model_XT_X(xinc,incidenceDesc,XX)

         ! Check all parameters to remove from the estimation
         call estim_cholesky(XX,incidenceDesc,incidenceDesc%ntniv,triang)
         ! compute the precision of each parameter
         if (performPrecision) call get_precision(XX,tempForConfusion,incidenceDesc)
        ! call debug_write_incidence(xinc,incidenceDesc)

         call set_corrxinc(xinc,incidenceDesc,xincreduit)
         call set_RHS(xincreduit,incidenceDesc,RHS)

         do i=1,incidenceDesc%nbnivest
           YTemp(i) = RHS(i)
           do j=i-1,1,-1
            YTemp(i) = YTemp(i) - YTemp(j)*triang(i,j)
           end do
           YTemp(i) = YTemp(i) / triang(i,i)
         end do

         do i=incidenceDesc%nbnivest,1,-1
          Bestim(i) = YTemp(i)
          do j=i+1,incidenceDesc%nbnivest
             Bestim(i)  = Bestim(i) - Bestim(j)*triang(j,i)
          end do
          Bestim(i) = Bestim(i) / triang(i,i)
         end do

         XB = matmul(xincreduit(:incidenceDesc%dataset%nkd,:incidenceDesc%nbnivest),Bestim(:incidenceDesc%nbnivest))
         XB = incidenceDesc%dataset%Y(1,:)-XB
         XB = XB*XB
         sigsq=sum(XB)

         sigsq = sigsq / incidenceDesc%dataset%nkd
         osigsquare = sigsq

     end subroutine modele_homoscedastic
!!***

!!****f* m_qtlmap_incidence_linear/modele_MA_homoscedastic
!! NAME
!!    modele_MA_homoscedastic
!! DESCRIPTION
!! NOTES
!! SOURCE
     subroutine modele_MA_homoscedastic(xinc,incidenceDesc,workstruct,osigsquare,Bestim,performPrecision,tempForConfusion)
         type(INCIDENCE_TYPE)                  , intent(inout)     :: incidenceDesc
         real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax), intent(in) :: xinc
         type(INCIDENCE_GEN_STRUCT)              ,intent(in)       :: workstruct
         real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%np)        ,intent(out)        :: osigsquare
         real (kind=dp) , dimension(incidenceDesc%ntnivmax)  ,intent(out)        :: Bestim
         logical                               ,intent(in)         :: performPrecision
         real (kind=dp),dimension(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax) ,intent(out)  :: tempForConfusion

         !local
         real (kind=dp) ,dimension(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax)   :: XMX
         real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax) :: xincreduit
         real (kind=dp) , dimension(incidenceDesc%ntnivmax)                         :: RHS,YTemp
         real (kind=dp) , dimension(incidenceDesc%ntniv,incidenceDesc%ntniv)        :: triang
         integer                                                                    :: i,j
         real (kind=dp)                                                             :: f,XB(incidenceDesc%dataset%nkd)

         call model_XT_V_X(xinc,incidenceDesc,workstruct%M,XMX)
         ! Check all parameters to remove from the estimation
         call estim_cholesky(XMX,incidenceDesc,incidenceDesc%ntniv,triang)
         ! compute the precision of each parameter
         if (performPrecision) call get_precision(XMX,tempForConfusion,incidenceDesc)

         !call debug_write_incidence(xinc,incidenceDesc)

         call set_corrxinc(xinc,incidenceDesc,xincreduit)
         call set_RHS_V(xincreduit,workstruct%M,incidenceDesc,RHS)

         do i=1,incidenceDesc%nbnivest
           YTemp(i) = RHS(i)
           do j=i-1,1,-1
            YTemp(i) = YTemp(i) - YTemp(j)*triang(i,j)
           end do
           YTemp(i) = YTemp(i) / triang(i,i)
         end do

         do i=incidenceDesc%nbnivest,1,-1
          Bestim(i) = YTemp(i)
          do j=i+1,incidenceDesc%nbnivest
             Bestim(i)  = Bestim(i) - Bestim(j)*triang(j,i)
          end do
          Bestim(i) = Bestim(i) / triang(i,i)
         end do

         XB = matmul(xincreduit(:incidenceDesc%dataset%nkd,:incidenceDesc%nbnivest),Bestim(:incidenceDesc%nbnivest))
         XB = incidenceDesc%dataset%Y(1,:)-XB
         osigsquare = dot_product(matmul(XB,workstruct%M),XB)
         osigsquare = osigsquare / incidenceDesc%dataset%nkd

         !stop

     end subroutine modele_MA_homoscedastic
!!***


!!****f* m_qtlmap_incidence_linear/modele_heteroscedastic
!!  NAME
!!    modele_heteroscedastic
!!  DESCRIPTION
!!  From a incidence matrix, get the B estimation and the variance of each sires in a heteroscedastic model
!!
!!   Bs   = (X'.V-1.X)-1 . X'.Vs -1 .Y
!!   SIGs = SUM   CD * [ Y    - X   . Bs]² / N
!!      i      jk   ijk   ijk    ijk          i
!!  NOTES
!!  SOURCE
     subroutine modele_heteroscedastic(xinc,incidenceDesc,startsig,osigsquare,Bestim,performPrecision,tempForConfusion)
         type(INCIDENCE_TYPE)                  , intent(inout)     :: incidenceDesc
         real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax), intent(in) :: xinc
         real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%np)        ,intent(in)         :: startsig
         real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%np)        ,intent(out)        :: osigsquare
         real (kind=dp) , dimension(incidenceDesc%ntnivmax)  ,intent(out)        :: Bestim
         logical                               ,intent(in)         :: performPrecision
         real (kind=dp),dimension(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax) ,intent(out)  :: tempForConfusion

         !local
         real (kind=dp) , dimension(incidenceDesc%ntnivmax)           :: RHS,YTemp
         real (kind=dp) ,dimension(incidenceDesc%dataset%nkd,incidenceDesc%dataset%nkd) :: Vinv
         real (kind=dp) ,dimension(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax)   :: XVX
         real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax)        :: xincreduit
         real (kind=dp) , dimension(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax)  :: XX,Mr
         real (kind=dp) , dimension(incidenceDesc%ntniv,incidenceDesc%ntniv)  :: triang
         integer                  :: ip,ERROR,it,i,j,kd1,kd2
         real (kind=dp)           :: f,XB(incidenceDesc%dataset%nkd),lrtv,diff,lrt,oldlrt,sig(incidenceDesc%dataset%data%genea%np)

         RHS=0.d0
         sig=startsig
         diff = 2*incidenceDesc%dataset%data%params%EPS_LINEAR_HETEROSCEDASTIC
         lrt=INIFINY_REAL_VALUE
         it=0

         do while ( diff > incidenceDesc%dataset%data%params%EPS_LINEAR_HETEROSCEDASTIC )
           if ( it < incidenceDesc%dataset%data%params%MAX_LINEAR_ITERATION) then
            it=it+1
            call set_VInv( incidenceDesc, sig,Vinv )
            call model_XT_V_X(xinc,incidenceDesc,Vinv,XVX)

            call estim_cholesky(XVX,incidenceDesc,incidenceDesc%ntniv,triang)
            call set_corrxinc(xinc,incidenceDesc,xincreduit)

            call set_RHS_V(xincreduit,Vinv,incidenceDesc,RHS)

            do i=1,incidenceDesc%nbnivest
              YTemp(i) = RHS(i)
              do j=i-1,1,-1
               YTemp(i) = YTemp(i) - YTemp(j)*triang(i,j)
              end do
              YTemp(i) = YTemp(i) / triang(i,i)
            end do

            do i=incidenceDesc%nbnivest,1,-1
             Bestim(i) = YTemp(i)
             do j=i+1,incidenceDesc%nbnivest
               Bestim(i)  = Bestim(i) - Bestim(j)*triang(j,i)
             end do
               Bestim(i) = Bestim(i) / triang(i,i)
            end do

            XB = matmul(xincreduit(:incidenceDesc%dataset%nkd,:incidenceDesc%nbnivest),Bestim(:incidenceDesc%nbnivest))
            XB = incidenceDesc%dataset%Y(1,:)-XB
            XB = XB*XB*incidenceDesc%dataset%CD(1,:)*incidenceDesc%dataset%CD(1,:)

            do ip=1,incidenceDesc%dataset%data%genea%np
             kd1=incidenceDesc%dataset%lSires(ip)%half_sib%firstKd
             kd2=incidenceDesc%dataset%lSires(ip)%half_sib%lastkd
             osigsquare(ip)=sum(XB(kd1:kd2)) / (kd2-kd1+1)
            end do

            oldlrt=lrt
            lrt = 0
            do ip=1,incidenceDesc%dataset%data%genea%np
             kd1=incidenceDesc%dataset%lSires(ip)%half_sib%firstKd
             kd2=incidenceDesc%dataset%lSires(ip)%half_sib%lastkd
             lrt = lrt +(kd2-kd1+1)*(log(osigsquare(ip)) - log(sig(ip)*sig(ip)))
             sig(ip)=sqrt(osigsquare(ip))
            end do
            diff = abs(lrt-oldlrt)
          else
           diff = incidenceDesc%dataset%data%params%EPS_LINEAR_HETEROSCEDASTIC
          end if
         end do


         if (performPrecision) then
          call get_precision(XVX,tempForConfusion,incidenceDesc)
         end if

      end subroutine modele_heteroscedastic
!!***


!!****f* m_qtlmap_incidence_linear/model_lin_h0
!!  NAME
!!    model_lin_h0
!!  DESCRIPTION
!!    compute under H0 the LRT and the solution
!!  INPUTS
!!    xinc             : incidence matrix
!!    incidenceDesc    : description of incidence matrix (see INCIDENCE_TYPE)
!!    workstruct       : configuration and information about analysis (see INCIDENCE_GEN_STRUCT)
!!    performPrecision : boolean to compute for each parameter the precision (fill the variable embeded INCIDENCE_TYPE%precis)
!!  OUTPUTS
!!    sigsquareEstime  : the variance of each sire family
!!    Bestim           : the solution of the estimation
!!
!!  NOTES
!!  SOURCE
      subroutine model_lin_h0(xinc,incidenceDesc,workstruct,sigsquareEstime,Bestim,performPrecision)
         type(INCIDENCE_TYPE)                  , intent(inout)     :: incidenceDesc
         real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax), intent(in) :: xinc
         type(INCIDENCE_GEN_STRUCT)              ,intent(in)       :: workstruct
         real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%np)        ,intent(out) :: sigsquareEstime
         real (kind=dp) , dimension(incidenceDesc%ntnivmax)  ,intent(out)        :: Bestim
         logical                               ,intent(in)         :: performPrecision!,tConf
      !   real (kind=dp),dimension(ntnivmax,ntnivmax) ,intent(out)  :: tempForConfusion

         real (kind=dp),dimension(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax) :: temp
         real (kind=dp)                              :: startsig(incidenceDesc%dataset%data%genea%np)
         integer :: hypothesis,ip

         if(is_parameter(KIND_HOMOSCEDASTIC,workstruct%type_model).and.&
            is_parameter(MODEL_UNITRAIT_POLYGENIC,workstruct%type_model)) then
             call modele_homoscedastic(xinc,incidenceDesc,sigsquareEstime,Bestim,performPrecision,temp)
         else if (is_parameter(KIND_HOMOSCEDASTIC,workstruct%type_model).and.&
                  is_parameter(MODEL_UNITRAIT_RELATIONSHIP,workstruct%type_model)) then
             call modele_MA_homoscedastic(xinc,incidenceDesc,workstruct,sigsquareEstime,Bestim,performPrecision,temp)
         else if (is_parameter(KIND_HETEROSCEDASTIC,workstruct%type_model).and.&
                  is_parameter(MODEL_UNITRAIT_POLYGENIC,workstruct%type_model)) then
             do ip=1,incidenceDesc%dataset%data%genea%np
               startsig(ip) = incidenceDesc%dataset%lSires(ip)%sig0(1)
             end do
             call modele_heteroscedastic(xinc,incidenceDesc,startsig,sigsquareEstime,Bestim,performPrecision,temp)
         else
             call stop_application("Devel Error : Unknown linear type : "//trim(str(workstruct%type_model)))
         end if

      end subroutine model_lin_h0
!!***


!!****f* m_qtlmap_incidence_linear/model_lin_hn
!!  NAME
!!    model_lin_hn
!!  DESCRIPTION
!!    compute under Hn the LRT and the solution
!!  INPUTS
!!    xinc             : incidence matrix
!!    incidenceDesc    : description of incidence matrix (see INCIDENCE_TYPE)
!!    curPos           : current point to test (see POSITION_LRT_INCIDENCE)
!!    workstruct       : configuration and information about analysis (see INCIDENCE_GEN_STRUCT)
!!    performPrecision : boolean to compute for each parameter the precision (fill the variable embeded INCIDENCE_TYPE%precis)
!!    tConf            : boolean to get an array used by the confusion method
!!    invlrt           : boolean to compute the LRT
!!  OUTPUTS
!!    sigsquareEstime  : the variance of each sire family
!!    Bestim           : the solution of the estimation
!!    tempForConfusion : the buffer array for the confusion function
!!  NOTES
!!  SOURCE
      subroutine model_lin_hn(xinc,incidenceDesc,curPos,workstruct,sigsquareEstime,Bestim,&
                              performPrecision,tConf,tempForConfusion,invlrt)
         type(INCIDENCE_TYPE)                  , intent(inout)     :: incidenceDesc
         real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%nd,incidenceDesc%ntnivmax), intent(in) :: xinc
         type(POSITION_LRT_INCIDENCE)            ,intent(inout)    :: curPos
         type(INCIDENCE_GEN_STRUCT)              ,intent(in)       :: workstruct
         real (kind=dp) , dimension(incidenceDesc%ntnivmax)  ,intent(out)        :: Bestim
         real (kind=dp) , dimension(incidenceDesc%dataset%data%genea%np)        ,intent(out) :: sigsquareEstime
         logical                               ,intent(in)         :: performPrecision,tConf
         real (kind=dp),dimension(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax) ,intent(out)  :: tempForConfusion
         logical                                     ,intent(in)   :: invlrt

         real (kind=dp),dimension(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax) :: temp
         real (kind=dp)                              :: startsig(incidenceDesc%dataset%data%genea%np)
         integer :: hypothesis,ip,kd1,kd2,ifem,jm

         !print *,ic
         !call debug_write_incidence(xinc,incidenceDesc)

         if(is_parameter(KIND_HOMOSCEDASTIC,workstruct%type_model).and.&
            is_parameter(MODEL_UNITRAIT_POLYGENIC,workstruct%type_model)) then
             call modele_homoscedastic(xinc,incidenceDesc,sigsquareEstime,Bestim,performPrecision,temp)
         else if (is_parameter(KIND_HOMOSCEDASTIC,workstruct%type_model).and.&
                  is_parameter(MODEL_UNITRAIT_RELATIONSHIP,workstruct%type_model)) then
             call modele_MA_homoscedastic(xinc,incidenceDesc,workstruct,sigsquareEstime,Bestim,performPrecision,temp)
         else if (is_parameter(KIND_HETEROSCEDASTIC,workstruct%type_model).and.&
                  is_parameter(MODEL_UNITRAIT_POLYGENIC,workstruct%type_model)) then
             startsig=sqrt(workstruct%sigsquare(workstruct%hypothesis-1,:,1))
             call modele_heteroscedastic(xinc,incidenceDesc,startsig,sigsquareEstime,Bestim,performPrecision,temp)
         else
             call stop_application("Devel Error : Unknown linear type : "//trim(str(workstruct%type_model)))
         end if

         curPos%lrtSires=0.d0
         curPos%lrtDams=0

         !compute LRT
         if (invLrt) then
             do ip=1,incidenceDesc%dataset%data%genea%np
              kd1=incidenceDesc%dataset%lSires(ip)%half_sib%firstKd
              kd2=incidenceDesc%dataset%lSires(ip)%half_sib%lastkd
              if (kd2>kd1) then
               do hypothesis=1,workstruct%hypothesis
                curPos%lrtSires(hypothesis,ip)=&
                (kd2-kd1+1)*(log(sigsquareEstime(ip))-log(workstruct%sigsquare(hypothesis-1,ip,1))) !  (IQ-1) / NQTL QTL
               end do
              end if

             !lrt pour les femmes (estimables)
             ifem=0
             do jm=incidenceDesc%dataset%data%genea%nmp(ip)+1,incidenceDesc%dataset%data%genea%nmp(ip+1)
               if ( .not. incidenceDesc%dataset%data%phenoAnimal%estime(incidenceDesc%ic,jm) ) cycle
                ifem=ifem+1
                kd1=incidenceDesc%dataset%lSires(ip)%full_sib(ifem)%firstKd
                kd2=incidenceDesc%dataset%lSires(ip)%full_sib(ifem)%lastkd
                if (kd2>kd1) then
                 do hypothesis=1,workstruct%hypothesis
                   curPos%lrtDams(hypothesis,jm)=&
                    (kd2-kd1+1)*(log(sigsquareEstime(ip))-log(workstruct%sigsquare(hypothesis-1,ip,1))) !  (IQ-1) / NQTL QTL
                 end do
                end if
             end do ! jm
            end do!ip


         else
          do ip=1,incidenceDesc%dataset%data%genea%np
           kd1=incidenceDesc%dataset%lSires(ip)%half_sib%firstKd
           kd2=incidenceDesc%dataset%lSires(ip)%half_sib%lastkd
           if (kd2>kd1) then
            do hypothesis=1,workstruct%hypothesis
             curPos%lrtSires(hypothesis,ip)= &
             (kd2-kd1+1)*(log(workstruct%sigsquare(hypothesis-1,ip,1)) - log(sigsquareEstime(ip))) !  (IQ-1) / NQTL QTL
            end do
           end if
             !lrt pour les femmes (estimables)

             ifem=0
             do jm=incidenceDesc%dataset%data%genea%nmp(ip)+1,incidenceDesc%dataset%data%genea%nmp(ip+1)
               if ( .not. incidenceDesc%dataset%data%phenoAnimal%estime(incidenceDesc%ic,jm) ) cycle
                ifem=ifem+1
                kd1=incidenceDesc%dataset%lSires(ip)%full_sib(ifem)%firstKd
                kd2=incidenceDesc%dataset%lSires(ip)%full_sib(ifem)%lastkd
                if (kd2>kd1) then
                 do hypothesis=1,workstruct%hypothesis
                   curPos%lrtDams(hypothesis,jm)=&
                    (kd2-kd1+1)*(log(workstruct%sigsquare(hypothesis-1,ip,1)) - log(sigsquareEstime(ip))) !  (IQ-1) / NQTL QTL
                 end do
                end if
             end do ! jm

          end do
         end if

         if (tConf) then
           tempForConfusion = temp
         end if
      end subroutine model_lin_hn
!!***


!****f* m_qtlmap_incidence_linear/opti_2qtl_linear_interaction
!  NAME
!    opti_2qtl_linear_interaction
!  DESCRIPTION
!  <DEV>
!
!  NOTES
!  SOURCE
!      subroutine opti_2qtl_linear_interaction(ic,workstruct,incsol,lrtsol,type)
!
!      integer                             , intent(in)      :: ic,type
!      type(INCIDENCE_GEN_STRUCT)           , intent(in)     :: workstruct
!      type(TYPE_INCIDENCE_SOLUTION)           ,intent(out)      :: incsol
!      type(TYPE_LRT_SOLUTION)                 ,intent(out)      :: lrtsol
!
!      integer :: n,init_n1
!      integer :: n1,ii,iq,npar2,ip,chr,chr2,nteff1,nteff2,nteffinter,listnteff(2),kd1,kd2
!      real (kind=dp) :: dx,dx1,f,lrt,sigsquare_t(np)
!
!      real (kind=dp) ,dimension(:,:),pointer      :: xinc,tempConfusion
!      real (kind=dp) , dimension(:) ,pointer      :: Bestim
!      real (kind=dp) , dimension(np)              :: sigsquareEstime
!      type(INCIDENCE_TYPE)                        :: incidenceDesc
!
!      type(POSITION_LRT_INCIDENCE)                :: curPos
!
!      call new(2,lrtsol)
!
!      !initialisation of incidence matrix
!      call init_incidence(ic,2,incidenceDesc,workstruct)
!
!      allocate (xinc(nd,incidenceDesc%ntnivmax))
!      xinc=0.d0
!      !Position Allocation
!      allocate (curPos%listN(2))
!      allocate (curPos%listChr(2))
!      allocate (curPos%lrt(2))
!
!      lrtsol%lrtmax=-INIFINY_REAL_VALUE
!
!      allocate (Bestim(incidenceDesc%ntnivmax))
!      allocate (tempConfusion(incidenceDesc%ntnivmax,incidenceDesc%ntnivmax))
!      Bestim=0.d0
!      chr=1;n=1;n1=n+1
!      !add general mean to estim
!      call add_general_mean(xinc,incidenceDesc)
!      nteff1=2
!      !add qtl effect/interaction at position n to estim
!      call add_qtleffect(n,chr,n,xinc,incidenceDesc,0,.true.)
!      nteff2=incidenceDesc%nteff+1
!      !add qtl effect/interaction at position n+1 to estim
!      call add_qtleffect(n1,chr,n+1,xinc,incidenceDesc,0,.true.)
!
!      !add interaction
!      nteffinter=incidenceDesc%nteff+1
!      call add_qtlinteraction(nteffinter,chr,chr,n,n1,xinc,incidenceDesc)
!      !call debug_write_incidence(xinc,incidenceDesc)
!      !add polygenic effect
!      call add_polygenic_effect(xinc,incidenceDesc)
!      !add fixed effect and covariate
!      call add_effcov(xinc,incidenceDesc,0,0)
!
!
!      !Initalisation of the maximum finded
!
!      do chr=1,nchr
!        curPos%listChr(1)=chr
!        do n=1,get_npo(chr)
!        curPos%listN(1)=n
!        !change qtl effect/interaction at position on the 2nd effect (the third is qtl dam) to estim
!        call change_qtleffect(nteff1,1,chr,n,xinc,incidenceDesc,0,.true.)
!        do chr2=chr,nchr
!         curPos%listChr(2)=chr2
!         init_n1=n+1
!         if ( chr2 /= chr) then
!           init_n1=1
!         end if
!         do n1=init_n1,get_npo(chr2)
!          curPos%listN(2)=n1
!          !change qtl effect/interaction at the nteff2 th position to estim
!          call change_qtleffect(nteff2,2,chr2,n1,xinc,incidenceDesc,0,.true.)
!          call change_interaction_effect(nteffinter,chr,chr2,n,n1,xinc,incidenceDesc)
!
!         ! call debug_write_incidence(xinc,incidenceDesc)
!        !  print *,"pos:",n
!          !call model
!          call model_lin_hn(xinc,incidenceDesc,curPos,workstruct,sigsquareEstime,Bestim,.false.,.false.,tempConfusion,.false.)
          !compute LRT
!          lrt = 0.d0
!          !print *,sigsquareEstime
!          do ip=1,np
!           kd1=incidenceDesc%dataset%lSires(ip)%half_sib%firstKd
!           kd2=incidenceDesc%dataset%lSires(ip)%half_sib%lastkd
!           if (kd2>kd1) then
!            lrt=lrt + (kd2-kd1+1)*(log(workstruct%sigsquare(2,ip,1)) - log(sigsquareEstime(ip)))
!           end if
!          end do
!        !  print *,n,n1,"LRT :",lrt
!        !  stop
!
!          lrtsol%lrt1_2(chr,chr2,n,n1)=lrt
!          lrtsol%lrt0_2(chr,chr2,n,n1)=lrt
!
!         ! print *,absi(chr,n),absi(chr2,n1),lrt
!          !print *,incidenceDesc%vecsol
!
!          !print *,incidenceDesc%vecsol(:incidenceDesc%ntniv)
!         ! stop
!
!
!          if ( lrtsol%lrtmax(0) < lrt ) then
!            lrtsol%nxmax(0)=n
!            lrtsol%nxmax(1)=n1
!            lrtsol%chrmax(0)=chr
!            lrtsol%chrmax(1)=chr2
!            lrtsol%lrtmax=lrt
!            print *,"LRT MAX:",lrt," pos1:",n," pos2:",n1,absi(chr,n),absi(chr2,n1)
!          end if
!       end do
!        end do
!      end do
!     end do
!
!      call change_qtleffect(nteff1,1,lrtsol%chrmax(0),lrtsol%nxmax(0),xinc,incidenceDesc,0,.true.)
!      call change_qtleffect(nteff2,2,lrtsol%chrmax(1),lrtsol%nxmax(1),xinc,incidenceDesc,0,.true.)
!      call change_interaction_effect(nteffinter,lrtsol%chrmax(0),lrtsol%chrmax(1),lrtsol%nxmax(0),lrtsol%nxmax(1),&
!                                     xinc,incidenceDesc)
!
!      call model_lin_hn(xinc,incidenceDesc,curPos,workstruct,sigsquareEstime,Bestim,.false.,.false.,tempConfusion,.false.)
!   !   call debug_write_incidence(xinc,incidenceDesc)
!
!   !   print *,"Bestime:"
!
!
!      listnteff(1)=nteff1
!      listnteff(2)=nteff2
!      call set_solution(2,workstruct,sigsquareEstime,Bestim,incidenceDesc,incsol,2,listnteff)
!
!      call end_incidence(incidenceDesc)
!      deallocate (xinc)
!      deallocate(Bestim)
!      deallocate (tempConfusion)
!      deallocate (curPos%listN)
!      deallocate (curPos%listChr)
!      deallocate (curPos%lrt)
!
!      end subroutine opti_2qtl_linear_interaction
!***


end module m_qtlmap_incidence_linear
