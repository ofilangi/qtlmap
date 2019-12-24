!!****m* HAPLOTYPE/m_qtlmap_haplotype_V1
!!  NAME
!!    m_qtlmap_haplotype_V1
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
module m_qtlmap_haplotype_V1
   use m_qtlmap_types
   use m_qtlmap_base
   use m_qtlmap_types
   use m_qtlmap_log
   use m_qtlmap_haplotype_util, only : transmi,gammapf

   implicit none
   save

   integer, private, parameter                              :: MAX_MARKER     = 39

!   real (kind=dp)   ,dimension(:,:,:,:),allocatable,private   :: prot
   !Haplotype possible dimmension lx,llx
   integer , dimension (:,:,:), allocatable,      private     :: h
   !Prob Father haplotype : dim llx,npx
   real (kind=dp), dimension (:,:,:), allocatable,private     :: prohp
   !Prob Mother haplotype : dim llx,nmx
   real (kind=dp), dimension (:,:,:), allocatable,private     :: prohm
   ! dim : lx,ngdx,4
   real (kind=dp)   ,dimension(:,:,:,:),allocatable           :: ptfin


   public :: haplotype_V1

   contains

!!****f* m_qtlmap_haplotype_V1/haplotype_V1
!! NAME
!!    haplotype_V1
!! DESCRIPTION
!!
!! NOTES
!!
!! SOURCE
   subroutine haplotype_V1(dataset,spt)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        type(PDD_BUILD)            ,intent(inout)         :: spt

        integer           :: stat,ll,c
        integer ,dimension(dataset%map%nchr)       :: valnpo
        type(GENEALOGY_BASE) , pointer :: dgenea
        type(GENOTYPE_BASE) , pointer :: dga

        dga => dataset%genoAnimal
        dgenea => dataset%genea

        if (maxval(dataset%map%nmk) > MAX_MARKER) then
          call stop_application("You can not used this version of haplotype calculation.")
        end if

        ! BUFFER STRUCTURES...
        !-------------------------
        ll = 2**maxval(dataset%map%nmk)

        allocate (h(dataset%map%nchr,maxval(dataset%map%nmk),ll),STAT=stat)
        call check_allocate(stat,'h [m_qtlmap_haplotype_V1]')
        allocate (prohp(dataset%map%nchr,ll,dgenea%np),STAT=stat)
        call check_allocate(stat,'prohp [m_qtlmap_haplotype_V1]')
        allocate (prohm(dataset%map%nchr,ll,dgenea%nm),STAT=stat)
        allocate (ptfin(dataset%map%nchr,maxval(dataset%map%nmk),(dgenea%nd*4),4),STAT=stat)
        allocate (spt%prot(dataset%map%nchr,maxval(dataset%map%nmk),dgenea%nd,4),STAT=stat)
        call check_allocate(stat,'prot [m_qtlmap_haplotype_V1]')

        call check_allocate(stat,'ptfin [m_qtlmap_haplotype_V1]')
        !HAPLOTYPE_DATA --> DONNEE PERSITENTE A L APPLI
        allocate (spt%genotyp(dataset%map%nchr,maxval(dataset%map%nmk),(size(dga%numero)),2),STAT=stat)
        call check_allocate(stat,'genotyp [m_qtlmap_haplotype_V2]')
        spt%genotyp=dga%nmanque
        allocate (spt%ngenom(dataset%map%nchr,dgenea%nm+1),STAT=stat)
        call check_allocate(stat,'ngenom [m_qtlmap_haplotype_V2]')
        allocate (spt%phasp(dataset%map%nchr,dgenea%np),STAT=stat)
        call check_allocate(stat,'phasp [m_qtlmap_haplotype_V2]')
        allocate (spt%phasm(dataset%map%nchr,dgenea%nm),STAT=stat)
        call check_allocate(stat,'phasm [m_qtlmap_haplotype_V2]')
        spt%phasp = .false.
        spt%phasm = .false.

        call combine(dataset)
        call ancetre(dataset)
        call gammapf(dataset,spt)
        call pdegp(dataset,spt)
        call pdegm(dataset,spt)

        call log_mess('Second dim of pdd:'//str(spt%ngend(1,size(spt%ngend,2))),DEBUG_DEF)
        do c=1,dataset%map%nchr
          valnpo(c)=dataset%map%get_npo(c)
        end do

        allocate( spt%pdd(dataset%map%nchr,maxval(spt%ngend),4,maxval(valnpo)) )
        spt%pdd=0.d0

        call pded(dataset,spt)

        deallocate(prohp)
        deallocate(prohm)
        deallocate(h)
        deallocate (ptfin)

   end subroutine
!!***

!******************************************************************************
!******************************************************************************
!
! Etablissement du tableau des haplotypes possibles avec l marqueurs
!
! Sous programme appele par qtl
!
! Version 1.0.
!
!******************************************************************************
!******************************************************************************
!
   subroutine combine(dataset)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        integer           :: ll,i,j,lh,nj,c

     do c=1,dataset%map%nchr
        ll=2**dataset%map%nmk(c)

        do i=1,dataset%map%nmk(c)
        do j=1,ll
            h(c,i,j)=1
          end do
        end do
        h(c,dataset%map%nmk(c),2)=2
        lh=1
        nj=1
    1   lh=lh+1
        nj=nj*2
        do j=1,nj
          do i=1,dataset%map%nmk(c)
            h(c,i,j+nj)=h(c,i,j)
          end do
          h(c,dataset%map%nmk(c)+1-lh,j+nj)=2
        end do
        if(lh.le.(dataset%map%nmk(c)-1))go to 1
     end do

      end subroutine combine



!******************************************************************************
!******************************************************************************
!
! Mise a zero des probabilites des phases impossibles pour les peres et
! pour les meres d apres les phenotypes des grands parents
!
! Sous programme appele par qtl
!
! Version 1.3.
!
!******************************************************************************
!******************************************************************************
!
      subroutine ancetre(dataset)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
! Divers
        logical                                     :: typgp,typgm
        integer(kind=KIND_PHENO)                    :: m1,m2,mc1,mc2,md1,md2
        real (kind=dp) , dimension(:,:),allocatable :: probh
        integer                                     :: ordre_t,l2,il,kr,ip,jm,ll
        integer                                     :: igp,ngm1,ngm2,jgm,nr1,nr2,stat,c
        type(GENEALOGY_BASE) , pointer :: dgenea
        type(GENOTYPE_BASE) , pointer :: dga

        dga => dataset%genoAnimal
        dgenea => dataset%genea

!
     do c=1,dataset%map%nchr
        l2=2**dataset%map%nmk(c)

        allocate (probh(l2,dgenea%nr),STAT=stat)
        call check_allocate(stat,'h [m_qtlmap_haplotype_V1]')
!
! Initialisation
        probh=1.d0
        prohp=1.d0
        prohm=1.d0

!
! Recherche des haplotypes impossibles pour les parents
! Recherche des phenotypes des grands parents
      do ll=1,dataset%map%nmk(c)
        do 10 igp=1,dgenea%ngp
          typgp=.false.
          if(dga%corregp(igp).ne.9999)then
           if(dga%pheno(c,ll,dga%corregp(igp),1).ne.dga%nmanque) then
            typgp=.true.
            m1=dga%pheno(c,ll,dga%corregp(igp),1)
            m2=dga%pheno(c,ll,dga%corregp(igp),2)
            end if
       end if
          ngm1=dgenea%ngmgp(igp)+1
          ngm2=dgenea%ngmgp(igp+1)
          do 10 jgm=ngm1,ngm2
          typgm=.false.
          if(dga%corregm(jgm).ne.9999)then
           if(dga%pheno(c,ll,dga%corregm(jgm),1).ne.dga%nmanque) then
            typgm=.true.
            mc1=dga%pheno(c,ll,dga%corregm(jgm),1)
            mc2=dga%pheno(c,ll,dga%corregm(jgm),2)
            end if
       end if
!
! Passage en revue des parents
          nr1=dgenea%nrgm(jgm)+1
          nr2=dgenea%nrgm(jgm+1)
          do 10 kr=nr1,nr2
          ordre_t=0
          if(dga%correr(kr).eq.9999) go to 10
          if((dga%pheno(c,ll,dga%correr(kr),1).eq.dga%nmanque)) go to 10
          md1=dga%pheno(c,ll,dga%correr(kr),1)
          md2=dga%pheno(c,ll,dga%correr(kr),2)
!
! Parent heterozygote
! Un allele different chez un des grands parents
          if((typgp).and.(md1.ne.m1.and.md1.ne.m2)) ordre_t=21
          if((typgp).and.(md2.ne.m1.and.md2.ne.m2)) ordre_t=12
          if((typgm).and.(md1.ne.mc1.and.md1.ne.mc2)) ordre_t=12
          if((typgm).and.(md2.ne.mc1.and.md2.ne.mc2)) ordre_t=21
!
! Mise a zero des probabilites des haplotypes impossibles
          do il=1,l2
            if((ordre_t.eq.12).and.(h(c,ll,il).eq.2)) probh(il,kr)=0.d0
            if((ordre_t.eq.21).and.(h(c,ll,il).eq.1)) probh(il,kr)=0.d0
          end do
   10   continue
      end do
!
! Stockage des probabilites pour les peres et les meres
      do il=1,l2
        do ip=1,dgenea%np
          if(dgenea%reppere(ip).ne.9999) prohp(c,il,ip)=probh(il,dgenea%reppere(ip))
        end do
        do jm=1,dgenea%nm
          if(dgenea%repmere(jm).ne.9999) prohm(c,il,jm)=probh(il,dgenea%repmere(jm))
        end do
      end do
        deallocate(probh)
    end do

   end subroutine ancetre


!******************************************************************************
!******************************************************************************
!
! Determination de la phase la plus probable pour les peres
!
! Sous programme appele par qtl
!
! Version 1.3.
!
!******************************************************************************
!******************************************************************************
!
      subroutine pdegp(dataset,spt)
         type(QTLMAP_DATASET)       ,intent(in)            :: dataset
         type(PDD_BUILD)            ,intent(inout)         :: spt

!
! Divers
        integer :: geno,t,tder,l2,ip,nm1,nm2,jm,nd1,nd2,kd,ll,ll1,ll2
        integer :: i,jmax1,jmax2,lll,gamma(maxval(dataset%map%nmk)),genp(dataset%genea%np),c
        real (kind=dp) :: sprob,xlvi,pmax,p1,p2,pr,pdeg
        real (kind=dp) :: prob(2**maxval(dataset%map%nmk)),probmax_t

        type(GENEALOGY_BASE) , pointer :: dgenea
        type(GENOTYPE_BASE) , pointer :: dga

        dga => dataset%genoAnimal
        dgenea => dataset%genea

      do c=1,dataset%map%nchr
!
! Tous les marqueurs doivent etre consideres
        l2=2**dataset%map%nmk(c)
!
! Etablissement de la phase la plus probable
      do 100 ip=1,dgenea%np
        nm1=dgenea%nmp(ip)+1
        nm2=dgenea%nmp(ip+1)
!
! Calcul de la probabilite a posteriori de la phase
        sprob=0.d0
        probmax_t=-1.d6
        do 200 geno=1,l2
          prob(geno)=0.d0
!	  write(6,*)'ip,geno,prohp(c,geno,ip)',ip,geno,prohp(c,geno,ip)
            if(prohp(c,geno,ip).eq.0.d0)go to 200
            xlvi=0.d0
          do jm=nm1,nm2
            nd1=dgenea%ndm(jm)+1
            nd2=dgenea%ndm(jm+1)
            do kd=nd1,nd2
!
! Reconstitution de l evenement de transmission
! a partir du tableau des probabilites de transmission
                do ll=1,dataset%map%nmk(c)
                  gamma(ll)=0
                  p1=spt%prot(c,ll,kd,1)+spt%prot(c,ll,kd,2)
                  if(p1.eq.1.d0)gamma(ll)=1
                  p2=spt%prot(c,ll,kd,3)+spt%prot(c,ll,kd,4)
                  if(p2.eq.1.d0)gamma(ll)=2
                end do
!
! Calcul de la probabilite du phenotype du descendant
              ll1=1
              do while (gamma(ll1).eq.0.and.ll1.lt.dataset%map%nmk(c))
                  ll1=ll1+1
              end do
              tder=gamma(ll1)
              if (tder.ne.0) then
                if (h(c,ll1,geno).eq.2) tder=3-tder
                  lll=ll1+1
                do ll2=lll,dataset%map%nmk(c)
                  t=gamma(ll2)
                    if (t.ne.0) then
                    if (h(c,ll2,geno).eq.2) t=3-t
                    if (t.eq.tder) then
                      pr=1.d0-dataset%map%rm(c,ll1,ll2)
                        else
                      pr=dataset%map%rm(c,ll1,ll2)
                    end if
                    tder=t
                    ll1=ll2
                    xlvi=xlvi+dlog(pr)
                  end if
                end do !ll2
              end if
	 !     write(6,*)'ip,jm,kd,(gamma(ll),ll=1,dataset%map%nmk(c)),xlvi',ip,jm,kd,(gamma(ll),ll=1,dataset%map%nmk(c)),xlvi
            end do !kd
          end do
          prob(geno)=xlvi
          if(probmax_t.le.xlvi) probmax_t=xlvi
  200   continue
!
! Standardisation
      do geno=1,l2
   !   write(6,*)'standardisation , prob(geno)',prob(geno)
        if (prob(geno).ne.0.d0) then
          prob(geno)=dexp(prob(geno)-probmax_t)
          sprob=sprob+prob(geno)
        end if
      end do
!
! Recherche de la phase la plus probable
        pmax=0.d0
        if ( sprob /= 0.d0 ) then
          do geno=1,l2
            pdeg=prob(geno)/sprob
	!    if(geno.eq.14)then
!	    write(6,*)'dans pdegp ip,geno,pdeg :',ip,geno,pdeg
!	write(6,*)(trim(get_pheno(dga,c,dga%pheno(c,i,dga%correp(ip),h(c,i,geno)))),i=1,dataset%map%nmk(c))
!	write(6,*)(trim(get_pheno(dga,c,dga%pheno(c,i,dga%correp(ip),h(c,i,l2+1-geno)))),i=1,dataset%map%nmk(c))
	!       end if
            if(pmax.lt.pdeg) then
              pmax=pdeg
              genp(ip)=geno
            end if
          end do
        end if
        if (pmax.eq.0.d0) then
          call stop_application('Genotype incompatibilities in sire family '//trim(dgenea%pere(ip)))
        end if
        spt%phasp(c,ip)=.true.
        jmax1=genp(ip)
        jmax2=l2+1-jmax1
        if(prob(jmax1).eq.prob(jmax2))spt%phasp(c,ip)=.false.
        do i=1,dataset%map%nmk(c)
            spt%genotyp(c,i,dga%correp(ip),1)=dga%pheno(c,i,dga%correp(ip),h(c,i,jmax1))
            spt%genotyp(c,i,dga%correp(ip),2)=dga%pheno(c,i,dga%correp(ip),h(c,i,jmax2))
        end do
!	write(6,*)'ip,pere(ip),correp(ip)',ip,trim(dgenea%pere(ip)),dga%correp(ip)
!	write(6,*)(trim(get_pheno(dga,c,spt%genotyp(c,i,dga%correp(ip),1))),i=1,dataset%map%nmk(c))
!	write(6,*)(trim(get_pheno(dga,c,spt%genotyp(c,i,dga%correp(ip),2))),i=1,dataset%map%nmk(c))
!
! Reorganisation du tableau des probabilites de transmission
          do jm=nm1,nm2
          nd1=dgenea%ndm(jm)+1
          nd2=dgenea%ndm(jm+1)
          do kd=nd1,nd2
              do ll=1,dataset%map%nmk(c)
                if(h(c,ll,genp(ip)).eq.2) then
                    p1=spt%prot(c,ll,kd,1)
                    p2=spt%prot(c,ll,kd,2)
                  spt%prot(c,ll,kd,1)=spt%prot(c,ll,kd,3)
                  spt%prot(c,ll,kd,2)=spt%prot(c,ll,kd,4)
                  spt%prot(c,ll,kd,3)=p1
                  spt%prot(c,ll,kd,4)=p2
                end if
              end do
            end do
          end do
  100 continue
    end do
      end subroutine pdegp

!******************************************************************************
!******************************************************************************
!
! Determination des probabilites des phases pour les meres
!
! Sous programme appele par qtl
! Appelle transmi
!
! Version 3.0.
!
!******************************************************************************
!******************************************************************************
!
      subroutine pdegm(dataset,spt)
      type(QTLMAP_DATASET)       ,intent(in)            :: dataset
      type(PDD_BUILD)            ,intent(inout)         :: spt

      integer   ,dimension(:,:) ,allocatable               :: ngend_t,ndesc_t
      real (kind=dp),dimension(:,:),allocatable            :: probg_t
      real (kind=dp), dimension(2**maxval(dataset%map%nmk))    :: prob
      real (kind=dp), dimension(2**maxval(dataset%map%nmk),dataset%genea%nm) :: probfem
      real (kind=dp) , dimension(dataset%genea%nm)       :: probmax,sprobfem

      real (kind=dp)  :: r1,r2,pr,spdeg,pmax
      integer ::i,j,l0,l1,l2,l3,ifem,jm,nd1,nd2,kd,ll,ll1,ll2,nombfem
      integer ::kont,lhomo,ngeno,ngeno1,ngeno2,imax,j1,j2,ii,il,c,s3

      logical ,  dimension(2**maxval(dataset%map%nmk))        :: garde
      integer(kind=KIND_PHENO)    ,dimension(:,:,:,:),allocatable    :: genotypm_t
!
! Divers
      integer geno,sym,t1,t2,homo,hetero,eff_t,t,genm(4*dataset%genea%nm)
      real (kind=dp) trans(4,4),s(4),q(4), p(maxval(dataset%map%nmk),dataset%genea%nd,4)
      integer(kind=KIND_PHENO) :: mc1,mc2

      type(GENEALOGY_BASE) , pointer :: dgenea
      type(GENOTYPE_BASE) , pointer :: dga

      dga => dataset%genoAnimal
      dgenea => dataset%genea
!
!

        allocate ( spt%ngend(dataset%map%nchr,1) )
        allocate ( spt%probg(dataset%map%nchr,1) )
        allocate ( spt%ndesc(dataset%map%nchr,1) )
        allocate ( spt%genotypm(dataset%map%nchr,maxval(dataset%map%nmk),dgenea%nm,2) )

        spt%ngenom(:,1)=0
        spt%ngend(:,1)=0
        spt%probg=0.d0

       do c=1,dataset%map%nchr
!
! Tous les marqueurs doivent etre consideres
        l2=2**dataset%map%nmk(c)
        l1=l2/2

!
! Initialisation des probabilites de genotype par femelle
        do ifem=1,dgenea%nfem
          sprobfem(ifem)=0.d0
          probmax(ifem)=-1.d6
          do geno=1,l2
            probfem(geno,ifem)=0.d0
          end do
        end do

!
! Calcul de la probabilite a posteriori de la phase par femelle
        do 101 jm=1,dgenea%nm
          nd1=dgenea%ndm(jm)+1
          nd2=dgenea%ndm(jm+1)
          ifem=dgenea%repfem(jm)
          do 100 geno=1,l2
            if(prohm(c,geno,jm).eq.0.d0)go to 100
            do kd=nd1,nd2
!
! Reorganisation du tableau des probabilites de transmission
              do ll=1,dataset%map%nmk(c)
                do i=1,4
                  p(ll,kd,i)=spt%prot(c,ll,kd,i)
                end do
                if(h(c,ll,geno).eq.2) then
                p(ll,kd,1)=spt%prot(c,ll,kd,2)
                p(ll,kd,2)=spt%prot(c,ll,kd,1)
                p(ll,kd,3)=spt%prot(c,ll,kd,4)
                p(ll,kd,4)=spt%prot(c,ll,kd,3)
                end if
              end do
!
! Initialisation
              do i=1,4
                q(i)=p(1,kd,i)
              end do
!
! Calcul de la probabilite du phenotype du descendant
            do ll2=2,dataset%map%nmk(c)
                ll1=ll2-1
              r1=dataset%map%rm(c,ll1,ll2)
                r2=dataset%map%rf(c,ll1,ll2)
              call transmi(r1,r2,trans)
                do t2=1,4
                s(t2)=0.d0
                do t1=1,4
                    s(t2)=s(t2)+trans(t1,t2)*q(t1)
                  end do
                end do
                do i=1,4
                q(i)=p(ll2,kd,i)*s(i)
                end do
              end do
            pr=0.d0
              do t=1,4
                pr=pr+q(t)
              end do
              probfem(geno,ifem)=probfem(geno,ifem)+dlog(pr)
          end do
          if(probmax(ifem).le.probfem(geno,ifem))probmax(ifem)=probfem(geno,ifem)
  100   continue
  101 continue
!
! Cumul des probabilites d'une meme femelle
      do ifem=1,dgenea%nfem
        do geno=1,l2
          if (probfem(geno,ifem).ne.0.d0) then
            probfem(geno,ifem)=dexp(probfem(geno,ifem)-probmax(ifem))
            sprobfem(ifem)=sprobfem(ifem)+probfem(geno,ifem)
          end if
        end do
        do geno=1,l2
          probfem(geno,ifem)=probfem(geno,ifem)/sprobfem(ifem)
        end do
      end do
!
! Stockage des probabilites par mere
        nombfem=0
        do 200 jm=1,dgenea%nm
          ifem=dgenea%repfem(jm)
          do geno=1,l2
            prob(geno)=probfem(geno,ifem)
          end do
!
! Si la phase est inconnue, les deux haplotypes symetriques ont
! la meme probabilite
        spt%phasm(c,jm)=.true.
          do geno=1,l1
            sym=l2+1-geno
            if((prob(geno).ne.0.d0).and.(prob(sym).ne.0.d0)) then
            prob(geno)=prob(geno)+prob(sym)
              prob(sym)=0.d0
              spt%phasm(c,jm)=.false.
            end if
          end do
!
! Cas des marqueurs homozygotes
! Regle de decision : le premier allele lu vient du grand pere
          homo=0
          do i=1,l2
            garde(i)=.true.
          end do
          lhomo=l2
          do ll=1,dataset%map%nmk(c)
            if(dga%correm(jm).eq.9999) then
              mc1=dga%nmanque
            else
              mc1=dga%pheno(c,ll,dga%correm(jm),1)
              mc2=dga%pheno(c,ll,dga%correm(jm),2)
            end if
            if(mc1.eq.dga%nmanque) mc2=mc1
            if(mc1.eq.mc2) then
              homo=homo+1
              l3=2**(dataset%map%nmk(c)-ll)
              lhomo=lhomo-l3
              do i=1,l2
                if(h(c,ll,i).eq.2) then
                  j=i-l3
                  prob(j)=prob(i)+prob(j)
                  prob(i)=0.d0
                  garde(i)=.false.
                end if
              end do
            end if
          end do
          if(.not.spt%phasm(c,jm))then
            hetero=dataset%map%nmk(c)-homo
            if(hetero.eq.0)then
              prob(1)=1.d0
              do i=2,l2
                prob(i)=0.d0
              end do
              go to 300
            end if
            l0=2**(hetero-1)
            kont=0
            do i=1,l2
              if(kont.lt.l0)then
                if(garde(i)) then
                  j=lhomo-i+1
                  prob(i)=prob(i)+prob(j)
                  prob(j)=0.d0
                  kont=kont+1
                end if
              end if
            end do
          end if
!
! Recherche des haplotypes les plus probables
  300   continue
        spdeg=0.d0
        nd1=dgenea%ndm(jm)+1
        nd2=dgenea%ndm(jm+1)
      !  eff_t=nd2-nd1+1
      !  if(eff_t.ge.ndmin) then
        if(dga%estfem(ifem))then
          do geno=1,l2
            if(prob(geno).gt.dataset%params%prseuil) then
              spdeg=spdeg+prob(geno)
            else
              prob(geno)=0.d0
            end if
          end do
!
! Si le nombre de pleins freres est trop faible, le genotype de la mere ne sera pas considere
! Dans ce cas le premier genotype rencontre possible est considere comme le bon
        else
          do geno=1,l2
            if(spdeg.eq.0.d0) then
              if(prohm(c,geno,jm).gt.0.d0) then
                prob(geno)=1.d0
                spdeg=prob(geno)
          !      force(jm)=.true.
              end if
            else
              prob(geno)=0.d0
            end if
          end do
        end if
!
!
        if(spdeg.eq.0.d0)then
           call stop_application('Dam '//trim(dgenea%mere(jm))//' has no possible haplotype  '//   &
          '(none has probability greater than '// trim(str(dataset%params%prseuil))//') : '//        &
          'check the genotype compatibility in the all sire family')
        end if
        do geno=1,l2
          prob(geno)=prob(geno)/spdeg
        end do
        spt%ngenom(c,jm+1)=spt%ngenom(c,jm)
        do geno=1,l2
          if(prob(geno).ne.0.d0) then
             allocate ( probg_t(dataset%map%nchr,size(spt%probg,2)) )
             probg_t(:,1:size(spt%probg,2)) = spt%probg
             deallocate ( spt%probg )
             allocate ( spt%probg(dataset%map%nchr,size(probg_t,2)+1) )
             spt%probg=0.d0
             spt%probg(:,1:size(probg_t,2))  = probg_t
             deallocate(probg_t)

            spt%ngenom(c,jm+1)=spt%ngenom(c,jm+1)+1
            genm(spt%ngenom(c,jm+1))=geno
            spt%probg(c,spt%ngenom(c,jm+1))=prob(geno)
          end if
        end do
        ngeno1=spt%ngenom(c,jm)+1
        ngeno2=spt%ngenom(c,jm+1)
        ngeno=ngeno2-ngeno1+1

        if ( ngeno2 > size(spt%genotypm,3)) then
         s3=size(spt%genotypm,3)
         allocate ( genotypm_t(dataset%map%nchr,maxval(dataset%map%nmk),s3,2) )
         genotypm_t = spt%genotypm
         deallocate (spt%genotypm)
         allocate ( spt%genotypm(dataset%map%nchr,maxval(dataset%map%nmk),ngeno2,2) )
         spt%genotypm = genotypm_t(:,:,:s3,:)
         deallocate (genotypm_t)
        end if

        if(dga%correm(jm).ne.9999)then
         !if (eff_t.ge.ndmin) then
         if(dga%estfem(ifem))then
          nombfem=nombfem+1
          imax=0;pmax=0
          do i=ngeno1,ngeno2
            j1=genm(i)
            j2=l2+1-j1
            if(spt%probg(c,i).gt.pmax) then
                 imax=i
                 pmax = spt%probg(c,i)
            end if
            do il=1,dataset%map%nmk(c)
             spt%genotypm(c,il,i,1)=dga%pheno(c,il,dga%correm(jm),h(c,il,j1))
             spt%genotypm(c,il,i,2)=dga%pheno(c,il,dga%correm(jm),h(c,il,j2))
          end do
          end do
!
!          j1=genm(imax)
!          j2=l2+1-j1
!          do il=1,dataset%map%nmk(c)
!           genotyp(c,il,correm(jm),1)=dga%pheno(c,il,correm(jm),h(c,il,j1))
!           genotyp(c,il,correm(jm),2)=dga%pheno(c,il,correm(jm),h(c,il,j2))
!          end do
	 end if
    end if
!
!  rep�rage de la phase la plus probable
!
!          pmax=0.d0
!          do geno=ngeno1,ngeno2
!            if(spt%probg(c,geno).gt.pmax) phasemax(jm)=geno
!          end do

!
! Reorganisation du tableau des probabilites de transmission
        allocate ( ngend_t(dataset%map%nchr,size(spt%ngend,2)) )
        ngend_t = spt%ngend
        deallocate( spt%ngend )
        allocate ( spt%ngend(dataset%map%nchr,size(ngend_t,2)+(ngeno2-ngeno1)+1) )
        spt%ngend=0.d0
        spt%ngend(:,1:size(ngend_t,2))=ngend_t
        deallocate(ngend_t)

        do geno=ngeno1,ngeno2
          spt%ngend(c,geno+1)=spt%ngend(c,geno)
          allocate ( ndesc_t(dataset%map%nchr,size(spt%ndesc,2)) )
          ndesc_t = spt%ndesc
          deallocate( spt%ndesc )
          allocate ( spt%ndesc(dataset%map%nchr,size(ndesc_t,2)+(nd2-nd1)+1) )
          spt%ndesc(:,1:size(ndesc_t,2))=ndesc_t
          deallocate( ndesc_t )
          do kd=nd1,nd2
            spt%ngend(c,geno+1)=spt%ngend(c,geno+1)+1
            spt%ndesc(c,spt%ngend(c,geno+1))=kd
            do ll=1,dataset%map%nmk(c)
              do i=1,4
                ptfin(c,ll,spt%ngend(c,geno+1),i)=spt%prot(c,ll,kd,i)
              end do
              if(h(c,ll,genm(geno)).eq.2) then
                ptfin(c,ll,spt%ngend(c,geno+1),1)=spt%prot(c,ll,kd,2)
                ptfin(c,ll,spt%ngend(c,geno+1),2)=spt%prot(c,ll,kd,1)
                ptfin(c,ll,spt%ngend(c,geno+1),3)=spt%prot(c,ll,kd,4)
                ptfin(c,ll,spt%ngend(c,geno+1),4)=spt%prot(c,ll,kd,3)
              end if
            end do
          end do
        end do
  200 continue
    end do

      return
      end subroutine pdegm

!******************************************************************************
!******************************************************************************
!
! Calcul des probabilites de transmission le long du chromosome
!
! Sous programme appele par qtl
! Appelle transmi, pdet, haldane et kosambi
!
! Version 2.0.
!
!******************************************************************************
!******************************************************************************
!
      subroutine pded(dataset,spt)
      type(QTLMAP_DATASET)       ,intent(in)            :: dataset
      type(PDD_BUILD)            ,intent(inout)         :: spt

!
! Tableaux dimensionnes selon nt, le nombre d'evenements de transmission
      real (kind=dp) , dimension(:),allocatable :: prob
      ! dim : ntx,l
      integer   , dimension(:,:),allocatable :: t
!
! Divers
      integer tm1(4,4),tm2(4,4)
      integer tf1(4,4),tf2(4,4)
      real (kind=dp) rpere(4),rmere(4)
      real (kind=dp) tab(4,4)

      integer :: ilong,ll,nm1,nm2,ngeno1,ngeno2,nd1,nd2,ip,n,i,j,geno
      integer :: kd,ix,it,jm,nlt,itg,itd,alloc_stat,c
      real (kind=dp) :: dx,r1,r2,xint,xintm,xintf,dg,dd,dgm,ddm,dgf,ddf
      real (kind=dp) :: rgm,rdm,rgf,rdf,rtm,rtf
      real (kind=dp) ,dimension(maxval(dataset%map%nmk),4,4)   :: tabr
      type(GENEALOGY_BASE) , pointer :: dgenea
      type(GENOTYPE_BASE) , pointer :: dga

      dga => dataset%genoAnimal
      dgenea => dataset%genea

!
      data tm1 /1,1,3,3,1,1,3,3,2,2,4,4,2,2,4,4/
      data tm2 /4,4,2,2,4,4,2,2,3,3,1,1,3,3,1,1/
      data tf1 /1,3,1,3,2,4,2,4,1,3,1,3,2,4,2,4/
      data tf2 /4,2,4,2,3,1,3,1,4,2,4,2,3,1,3,1/


   do c=1,dataset%map%nchr
!******************************************************************************
!
! Taille du segment explore
      ilong=dataset%map%get_ilong(c)
!
! Tableau des probabilites des evenements de recombinaison
        do ll=2,dataset%map%nmk(c)
          r1=dataset%map%rm(c,ll-1,ll)
          r2=dataset%map%rf(c,ll-1,ll)
          call transmi(r1,r2,tab)
          do i=1,4
            do j=1,4
              tabr(ll,i,j)=tab(i,j)
            end do
          end do
        end do
!
! Calcul pour chaque descendant de la probabilite de transmission
        do 100 ip=1,dgenea%np
        nm1=dgenea%nmp(ip)+1
        nm2=dgenea%nmp(ip+1)
        do 100 jm=nm1,nm2
        ngeno1=spt%ngenom(c,jm)+1
        ngeno2=spt%ngenom(c,jm+1)
        do 100 geno=ngeno1,ngeno2
        nd1=spt%ngend(c,geno)+1
        nd2=spt%ngend(c,geno+1)
        do 100 kd=nd1,nd2
          nlt = get_nlt(dataset,spt,c,kd)

          if ( nlt >= 2**20 ) then
            call log_mess(" ** progeny :"//trim(dgenea%animal(kd)),WARNING_DEF)
            call log_mess(" ** transmission possibility :"//trim(str(nlt)),WARNING_DEF)
          end if

          if (allocated (t)) then
            deallocate (t)
          end if
          if (allocated(prob)) then
            deallocate (prob)
          end if
          !allocation
          allocate (t(nlt,dataset%map%nmk(c)),stat = alloc_stat)
          call check_allocate(alloc_stat,'Not enough memory to use this computation of transmission probabilities ')

          allocate (prob(nlt),stat = alloc_stat)
          call check_allocate(alloc_stat,'Not enough memory to use this computation of transmission probabilities ')
          prob = 0.d0
!
! Calcul des probabilites des evenements de transmission du segment
          call pdet(dataset,spt,c,kd,nlt,tabr,t,prob)
!
! Definition de la position dx
          n=0
          do 200 ix=0,ilong,dataset%map%pas
            n=n+1
            dx=dataset%map%absi(c,n)
            ll=1
   2       if(dataset%map%posi(c,ll+1).ge.dx) then
              xint=dataset%map%posi(c,ll+1)-dataset%map%posi(c,ll)
              dg=dx-dataset%map%posi(c,ll)
              dd=dataset%map%posi(c,ll+1)-dx
              xintm=dataset%map%posim(c,ll+1)-dataset%map%posim(c,ll)
              dgm=dg*xintm/xint
              ddm=dd*xintm/xint
              xintf=dataset%map%posif(c,ll+1)-dataset%map%posif(c,ll)
              dgf=dg*xintf/xint
              ddf=dd*xintf/xint
              rgm=xaldane(dgm)
              rdm=xaldane(ddm)
              rgf=xaldane(dgf)
              rdf=xaldane(ddf)
            else
              ll=ll+1
              if(ll.lt.( dataset%map%nmk(c)-1) )go to 2
            end if
!
! On se place sous l'hypothese d'absence d'interference pour la detection de QTL
            rtm=rgm+rdm-(2.d0*rgm*rdm)
            rtf=rgf+rdf-(2.d0*rgf*rdf)
            rpere(1)=(1.d0-rgm)*(1.d0-rdm)/(1.d0-rtm)
            rpere(2)=(1.d0-rgm)*rdm/rtm
            rpere(3)=rgm*(1.d0-rdm)/rtm
            rpere(4)=rgm*rdm/(1.d0-rtm)
            rmere(1)=(1.d0-rgf)*(1.d0-rdf)/(1.d0-rtf)
            rmere(2)=(1.d0-rgf)*rdf/rtf
            rmere(3)=rgf*(1.d0-rdf)/rtf
            rmere(4)=rgf*rdf /(1.d0-rtf)
!
!
! Initialisation
          do it=1,4
            spt%pdd(c,kd,it,n)=0.d0
          end do
!
! Passage en revue des transmissions possibles
          do i=1,nlt
            itg=t(i,ll)
            itd=t(i,ll+1)
            spt%pdd(c,kd,1,n)= spt%pdd(c,kd,1,n)+rpere(tm1(itg,itd))*rmere(tf1(itg,itd))*prob(i)
            spt%pdd(c,kd,2,n)= spt%pdd(c,kd,2,n)+rpere(tm1(itg,itd))*rmere(tf2(itg,itd))*prob(i)
            spt%pdd(c,kd,3,n)= spt%pdd(c,kd,3,n)+rpere(tm2(itg,itd))*rmere(tf1(itg,itd))*prob(i)
            spt%pdd(c,kd,4,n)= spt%pdd(c,kd,4,n)+rpere(tm2(itg,itd))*rmere(tf2(itg,itd))*prob(i)
          end do!

  200   continue
  100   continue

      if (allocated (t)) then
         deallocate (t)
      end if
      if (allocated(prob)) then
         deallocate (prob)
      end if
   end do
      end subroutine pded

!************************************************
! SUBROUTINE : get_nlt : transmition number
!************************************************
      function get_nlt(dataset,spt,c,kd) result(nlt)
            type(QTLMAP_DATASET)       ,intent(in)         :: dataset
            type(PDD_BUILD)            ,intent(in)         :: spt
            integer , intent(in)          :: c,kd
            integer                       :: nlt,i,ll,kont

             nlt = 0
             do i=1,4
              if(ptfin(c,1,kd,i).ne.0.d0) then
                nlt=nlt+1
              end if
             end do

              do ll=2,dataset%map%nmk(c)
                kont=0
                do i=1,4
                  if(ptfin(c,ll,kd,i).ne.0.d0) then
                    kont=kont+1
                  end if
                end do
                nlt=nlt*kont
              end do
            return
          end function get_nlt

!******************************************************************************
!******************************************************************************
!
! Calcul des probabilites de transmission des segments chromosomiques � un
! descendant sachant les haplotypes de ses parents
!
! Sous programme appele par pded
!
! Version 2.0.
!
!******************************************************************************
!******************************************************************************
!

      subroutine pdet(dataset,spt,c,kd,nlt,tabr,t,prob)
      type(QTLMAP_DATASET)       ,intent(in)            :: dataset
      type(PDD_BUILD)            ,intent(inout)         :: spt

      integer          , intent(in)                 :: c,kd
      integer          , intent(in)                 :: nlt
      real (kind=dp) ,intent(in),dimension(:,:,:)   :: tabr

      integer        ,intent(inout), dimension(:,:) :: t
      real (kind=dp) ,intent(inout),dimension(nlt)  :: prob

      ! dim : ntx,l
      integer  , dimension (:,:),allocatable       :: gardt
      real (kind=dp) , dimension (:)  ,allocatable :: gardp
!
! Divers
      integer :: gamma(4),i,ll,kont,ii,iii,j,val_nlt
      real (kind=dp) :: sprob

      type(GENEALOGY_BASE) , pointer :: dgenea
      type(GENOTYPE_BASE) , pointer :: dga

      dga => dataset%genoAnimal
      dgenea => dataset%genea

      allocate (gardt(nlt,dataset%map%nmk(c)))
      allocate (gardp(nlt))
      gardt = 0.d0
      gardp = 0.d0

!
! Creation pour le descendant du tableau des transmissions possibles
          val_nlt = 0

! Premier marqueur
          do i=1,4
            if(ptfin(c,1,kd,i).ne.0.d0) then
              val_nlt = val_nlt + 1
              t(val_nlt,1)=i
              prob(val_nlt)=ptfin(c,1,kd,i)
          end if
          end do
!
! Marqueurs suivants
          do ll=2,dataset%map%nmk(c)
!
! Transmissions possibles
          kont=0
          do i=1,4
              if(ptfin(c,ll,kd,i).ne.0.d0) then
                kont=kont+1
                gamma(kont)=i
              end if
          end do
!
! Repetition des lignes pour les marqueurs precedents
! et remplissage de la colonne correspondant au marqueur en cours
          do i=1,val_nlt
              do ii=1,kont
                iii=(kont*(i-1))+ii
                do j=1,ll-1
                  gardt(iii,j)=t(i,j)
              end do
              gardt(iii,ll)=gamma(ii)
              gardp(iii)=prob(i)*tabr(ll,gardt(iii,ll-1),gardt(iii,ll))
              end do
            end do
            val_nlt=val_nlt*kont
            do i=1,val_nlt
              do j=1,ll
                t(i,j)=gardt(i,j)
              end do
              prob(i)=gardp(i)
            end do
        end do
!
          sprob=0.d0
          do i=1,val_nlt
            sprob=sprob+prob(i)
          end do
          do i=1,val_nlt
            prob(i)=prob(i)/sprob
          end do

          deallocate (gardt)
          deallocate (gardp)
      return
      end subroutine pdet


end module m_qtlmap_haplotype_V1
