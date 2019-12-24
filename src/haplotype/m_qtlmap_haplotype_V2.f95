!!****m* HAPLOTYPE/m_qtlmap_haplotype_V2
!!  NAME
!!    m_qtlmap_haplotype_V2
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
module m_qtlmap_haplotype_V2
    use m_qtlmap_base
    use m_qtlmap_types
    use m_qtlmap_log
    use m_qtlmap_haplotype_util, only : transmi,gammapf
    use m_qtlmap_isymmax2sat

    implicit none
    save
    private

!! DESCRIPTION
!!   The maximum number of marker to use the pded_V5 analysis
!! NOTES
!!   there are a lot of allocation 2**nmk ==> overflow if nmk is too big
    integer, private, parameter                      :: MAX_MARKER     = 39

!! DESCRIPTION
!!   The maximum number of marker to use the pded_V5 analysis
!! NOTES
!!  Constante pour eviter les segmentation fault du au overflow 2 < MAX_CAS_ALLOC < 50
    integer , private,    parameter                  :: MAX_CAS_ALLOC = 22
    integer , dimension (:,:,:), allocatable      ,private   :: ordrep
    integer , dimension (:,:,:), allocatable      ,private   :: ordrem
    integer , dimension (:,:,:), allocatable       ,private  :: ordre
    integer , dimension (:,:,:), allocatable       ,private  :: ordref
    real (kind=dp),dimension(:,:,:,:),allocatable  ,private  :: ptfin

    public :: set_allocation
    public :: free_internal_struct
    public :: setting_alloc_pdd
    public :: haplotype_V2
    public :: haplotype_V3
    public :: haplotype_SNP
    public :: haplotype_SYMMAX2SAT_SNP

    public :: pded_v5_optim
    public :: check_genotype_mod
    public :: check_recombination_sire

    public :: internal_get_information_informative_marker

   contains

!!****f* m_qtlmap_haplotype_V2/set_allocation
!! NAME
!!    set_allocation
!! DESCRIPTION
!!
!! NOTES
!!
!! SOURCE
    subroutine set_allocation(dataset,spt)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        type(PDD_BUILD)            ,intent(inout)         :: spt

        integer  :: stat
        type(GENEALOGY_BASE) , pointer :: dgenea
        type(GENOTYPE_BASE) , pointer :: dga

        dga => dataset%genoAnimal
        dgenea => dataset%genea

        allocate (ordrep(dataset%map%nchr,dgenea%np,maxval(dataset%map%nmk)),STAT=stat)
        call check_allocate(stat,'ordrep [m_qtlmap_haplotype_V2]')
        allocate (ordrem(dataset%map%nchr,dgenea%nm,maxval(dataset%map%nmk)),STAT=stat)
        call check_allocate(stat,'ordrem [m_qtlmap_haplotype_V2]')
        allocate (ordre(dataset%map%nchr,dgenea%nr,maxval(dataset%map%nmk)),STAT=stat)
        call check_allocate(stat,'ordre  [m_qtlmap_haplotype_V2]')
        allocate (ordref(dataset%map%nchr,dgenea%nm,maxval(dataset%map%nmk)),STAT=stat)
        call check_allocate(stat,'ordref [m_qtlmap_haplotype_V2]')
        allocate (spt%prot(dataset%map%nchr,maxval(dataset%map%nmk),dgenea%nd,4),STAT=stat)
        call check_allocate(stat,'prot [m_qtlmap_haplotype_V2]')
        allocate (ptfin(dataset%map%nchr,maxval(dataset%map%nmk),dgenea%nd,4),STAT=stat)
        call check_allocate(stat,'ptfin [m_qtlmap_haplotype_V2]')

       ! a voir allocation genotyp et genotypm...... : la 3eme dim est approximatif...

        allocate (spt%genotyp(dataset%map%nchr,maxval(dataset%map%nmk),(size(dga%numero)),2),STAT=stat)
        call check_allocate(stat,'genotyp [m_qtlmap_haplotype_V2]')
        allocate (spt%genotypm(dataset%map%nchr,maxval(dataset%map%nmk),(size(dga%numero)),2),STAT=stat)
        call check_allocate(stat,'genotypm [m_qtlmap_haplotype_V2]')

        spt%genotyp=dga%nmanque
        allocate (spt%ngenom(dataset%map%nchr,dgenea%nm+1),STAT=stat)
        call check_allocate(stat,'spt%ngenom( [m_qtlmap_haplotype_V2]')
        allocate (spt%phasp(dataset%map%nchr,dgenea%np),STAT=stat)
        call check_allocate(stat,'phasp [m_qtlmap_haplotype_V2]')
        allocate (spt%phasm(dataset%map%nchr,dgenea%nm),STAT=stat)
        call check_allocate(stat,'phasm [m_qtlmap_haplotype_V2]')

    end subroutine set_allocation
!!***

!!****f* m_qtlmap_haplotype_V2/free_internal_struct
!! NAME
!!    free_internal_struct
!! DESCRIPTION
!!
!! NOTES
!!
!! SOURCE
    subroutine free_internal_struct

        deallocate (ordrep)
        deallocate (ordrem)
        deallocate (ordre)
        deallocate (ordref)
        deallocate (ptfin)

    end subroutine free_internal_struct
!!***

!!****f* m_qtlmap_haplotype_V2/check_genotype_mod
!! NAME
!!    check_genotype_mod
!! DESCRIPTION
!!
!! NOTES
!!
!! SOURCE
    subroutine check_genotype_mod(dataset,spt)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        type(PDD_BUILD)            ,intent(inout)         :: spt

        call set_allocation(dataset,spt)
        call ancetre(dataset)
        call set_phasp(dataset,spt)
        call set_phasm(dataset,spt)
        call gammapf(dataset,spt)
        call check_recombination_sire(dataset,spt)

        call free_internal_struct

   end subroutine check_genotype_mod
!!***

!!****f* m_qtlmap_haplotype_V2/setting_alloc_pdd
!! NAME
!!    setting_alloc_pdd
!! DESCRIPTION
!!
!! NOTES
!!
!! SOURCE
    subroutine setting_alloc_pdd(dataset,spt)
       type(QTLMAP_DATASET)        ,intent(in)       :: dataset
       type(PDD_BUILD)            ,intent(inout)     :: spt
       integer                        :: stat,c
       integer ,dimension(dataset%map%nchr)       :: valnpo

       call log_mess('Second dim of pdd:'//trim(str(maxval(spt%ngend))),DEBUG_DEF)
        do c=1,dataset%map%nchr
          valnpo(c)=dataset%map%get_npo(c)
        end do
        allocate( spt%pdd(dataset%map%nchr,maxval(spt%ngend),4,maxval(valnpo)) )
        spt%pdd=0.d0
    end subroutine setting_alloc_pdd
!!***



!!****f* m_qtlmap_haplotype_V2/get_information_informative_marker
!! NAME
!!    get_information_informative_marker
!! DESCRIPTION
!!    get information about informative marker.
!!    countH  => number of heterozygote marker by sire and linkage group
!!    markerH => array of boolean, true if for one marker and a given sire and a given linkage group the marker are heterozyte
!! NOTES
!!
!! SOURCE
    subroutine internal_get_information_informative_marker(dataset,spt,countH,markH)
      type(QTLMAP_DATASET)       ,intent(in)            :: dataset
      type(PDD_BUILD)            ,intent(inout)         :: spt
      real (kind=dp) , dimension(dataset%genea%np,dataset%map%nchr), intent(out)                :: countH
      real (kind=dp) , dimension(dataset%genea%np,dataset%map%nchr,maxval(dataset%map%nmk)) ,intent(out)    :: markH

        call set_allocation(dataset,spt)
        call ancetre(dataset)
        call gammapf(dataset,spt)

        call count_nbheterozygote_mark_byChr(dataset,spt,countH,markH)

        call free_internal_struct

   end subroutine internal_get_information_informative_marker
!!***

!!****f* m_qtlmap_haplotype_V2/haplotype_V2
!! NAME
!!    haplotype_V2
!! DESCRIPTION
!!
!! NOTES
!!
!! SOURCE
    subroutine haplotype_V2(dataset,spt)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        type(PDD_BUILD)            ,intent(inout)         :: spt

        call set_allocation(dataset,spt)
        call ancetre(dataset)
        call set_phasp(dataset,spt)
        call set_phasm(dataset,spt)
        call gammapf(dataset,spt)
        call pdegp(dataset,spt)
        call pdegm(dataset,spt)
        call setting_alloc_pdd(dataset,spt)
        call pded(dataset,spt)
        
        call free_internal_struct

   end subroutine haplotype_V2
!!***

!!****f* m_qtlmap_haplotype_V2/haplotype_V3
!! NAME
!!    haplotype_V3
!! DESCRIPTION
!!
!! NOTES
!!
!! SOURCE
    subroutine haplotype_V3(dataset,spt)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        type(PDD_BUILD)            ,intent(inout)         :: spt
        integer :: kd,c,ll

        call set_allocation(dataset,spt)
        call ancetre(dataset)
        call set_phasp(dataset,spt)
        call set_phasm(dataset,spt)
        call gammapf(dataset,spt)

        call pdegp(dataset,spt)
        call pdegm(dataset,spt)
        call setting_alloc_pdd(dataset,spt)
        !call pded_v5()
        call pded_v5_optim(dataset,spt)
        
        call free_internal_struct

   end subroutine haplotype_V3
!!***

!!****f* m_qtlmap_haplotype_V2/haplotype_SNP
!! NAME
!!    haplotype_SNP
!! DESCRIPTION
!!
!! NOTES
!!
!! SOURCE
    subroutine haplotype_SNP(dataset,spt)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        type(PDD_BUILD)            ,intent(inout)         :: spt

        call set_allocation(dataset,spt)
        call ancetre(dataset)
        call set_phasp(dataset,spt)
        call set_phasm(dataset,spt)
        call gammapf(dataset,spt)

        call pdegp_snp(dataset,spt)
        call pdegm_snp(dataset,spt)
        call setting_alloc_pdd(dataset,spt)
        !call pded_v5()
        call pded_v5_optim(dataset,spt)
        
        call free_internal_struct

   end subroutine haplotype_SNP
!!***


!!****f* m_qtlmap_haplotype_V2/haplotype_SYMMAX2SAT_SNP
!! NAME
!!    haplotype_SNP
!! DESCRIPTION
!!
!! NOTES
!!
!! SOURCE
    subroutine haplotype_SYMMAX2SAT_SNP(dataset,spt)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        type(PDD_BUILD)            ,intent(inout)         :: spt

        call set_allocation(dataset,spt)
        call ancetre(dataset)
        call set_phasp(dataset,spt)
        call set_phasm(dataset,spt)
        call gammapf(dataset,spt)

        call calcul_phases_symmax2sat_sire(dataset,spt)
        call calcul_phases_symmax2sat_dam(dataset,spt)
        call setting_alloc_pdd(dataset,spt)
        call pded_v5_optim(dataset,spt)
       
        call free_internal_struct

   end subroutine haplotype_SYMMAX2SAT_SNP
!!***
!!****f* m_qtlmap_haplotype_V2/ancetre
!! NAME
!!    ancetre
!! DESCRIPTION
!!
!! NOTES
!!   Mise a zero des probabilites des phases impossibles pour les peres et
!!   pour les meres d apres les phenotypes des grands parents
!! SOURCE
      subroutine ancetre (dataset)
       type(QTLMAP_DATASET)       ,intent(in)            :: dataset
!
! Divers
      logical typgp,typgm
      integer(kind=KIND_PHENO) :: m1,m2,mc1,mc2,md1,md2
!
      integer ll,igp,jgm,ip,jm,kr,c
      integer ngm1,ngm2,nr1,nr2
      type(GENEALOGY_BASE) , pointer :: dgenea
      type(GENOTYPE_BASE) , pointer :: dga

      dga => dataset%genoAnimal
      dgenea => dataset%genea
!
! Recherche des haplotypes impossibles pour les parents
! Recherche des phenotypes des grands parents
    do c=1,dataset%map%nchr
      do ll=1,dataset%map%nmk(c)
        do 10 igp=1,dgenea%ngp
          typgp=.false.
          if(dga%corregp(igp).ne.INT_NOT_DEFINED) then
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
            if(dga%corregm(jgm).ne.INT_NOT_DEFINED)then
             if (dga%pheno(c,ll,dga%corregm(jgm),1).ne.dga%nmanque) then
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
            ordre(c,kr,ll)=0
            if(dga%correr(kr).eq.INT_NOT_DEFINED)  go to 10
            if(dga%pheno(c,ll,dga%correr(kr),1).eq.dga%nmanque)  go to 10
            md1=dga%pheno(c,ll,dga%correr(kr),1)
            md2=dga%pheno(c,ll,dga%correr(kr),2)
!
! Parent heterozygote
! Un allele different chez un des grands parents
            if((typgp).and.(md1.ne.m1.and.md1.ne.m2))   ordre(c,kr,ll)=21
            if((typgp).and.(md2.ne.m1.and.md2.ne.m2))   ordre(c,kr,ll)=12
            if((typgm).and.(md1.ne.mc1.and.md1.ne.mc2)) ordre(c,kr,ll)=12
            if((typgm).and.(md2.ne.mc1.and.md2.ne.mc2)) ordre(c,kr,ll)=21
!
      10   continue
      end do
!
! Stockage des probabilites pour les peres et les meres
!  La gestion des valeurs de phasp et phasm (d�crite ici en commentaires) est
!  d�port�e dans get_phasp et get_phasm
      do ll=1,dataset%map%nmk(c)
        do ip=1,dgenea%np
         ! phasp(c,ip)=.false.
          ordrep(c,ip,ll)=0
          if(dgenea%reppere(ip).ne.INT_NOT_DEFINED) &
            ordrep(c,ip,ll)=ordre(c,dgenea%reppere(ip),ll)
          ! if(ordrep(ip,ll).ne.0)phasp(c,ip)=.true.
        end do
        do jm=1,dgenea%nm
          !phasm(c,jm)=.false.
          ordrem(c,jm,ll)=0
          ordref(c,dgenea%repfem(jm),ll)=0
          if(dgenea%repmere(jm).ne.INT_NOT_DEFINED) then
            ordrem(c,jm,ll)=ordre(c,dgenea%repmere(jm),ll)
            ordref(c,dgenea%repfem(jm),ll)=ordre(c,dgenea%repmere(jm),ll)
          end if
         ! if(ordrem(jm,ll).ne.0)phasm(c,jm)=.true.
        end do
      end do
     end do
   end subroutine ancetre
!!***

!!****f* m_qtlmap_haplotype_V2/count_nbheterozygote_mark_byChr
!! NAME
!!    count_nbheterozygote_mark_byChr
!! DESCRIPTION
!!    countH : informativité sur le nombre de marqueur informatif pour un pere et un groupe de liaison donnés
!!    markH  : informativité pour chaque marqueur pour un pere et un groupe de liaison donnés
!!
!! NOTES
!!   calcul de l'informativité :
!!   - +1 pour chaque descendant et chaque marqueur informatif.
!!   - on divise ce nombre par /nbdesc+nmk
!!   le resultat est mis dans le tableau counth
!!
!!      Utilisation des structures de donnees
!!      spt%prot(c,desc,ll,cas) : probabilite de transmission du marqueur ll pour le descenant desc
!! SOURCE
   subroutine count_nbheterozygote_mark_byChr(dataset,spt,countH,markH)
    type(QTLMAP_DATASET)       ,intent(in)            :: dataset
    type(PDD_BUILD)            ,intent(inout)         :: spt
    real (kind=dp) , dimension(dataset%genea%np,dataset%map%nchr), intent(out)                :: countH
    real (kind=dp) , dimension(dataset%genea%np,dataset%map%nchr,maxval(dataset%map%nmk)),intent(out)     :: markH

    !local
    integer :: c,ip,jm,kd,ll,kkd
    real (kind=dp) :: p1,p2,mHnb,cHnb
    type(GENEALOGY_BASE) , pointer :: dgenea

    dgenea => dataset%genea

!   nombre de marqueur heterozygote par pere et par groupe de liaison
    countH=0
!   marker heterozygote par pere/groupe de liaison
    markH=0

!
!  creation du tableau des transmissions observees
!
   do c=1,dataset%map%nchr
     do ip=1,dgenea%np
       cHnb=0
       mHnb=0
       do jm=dgenea%nmp(ip)+1,dgenea%nmp(ip+1)
         do kd=dgenea%ndm(jm)+1,dgenea%ndm(jm+1)
           mHnb=mHnb+1
           do ll=1,dataset%map%nmk(c)
              cHnb=cHnb+1
              p2=spt%prot(c,ll,kd,3)+spt%prot(c,ll,kd,4)
              if ( p2 == 1.d0 ) then
                countH(ip,c) =countH(ip,c)+1.d0 ! marqueur informatif
                markH(ip,c,ll)= markH(ip,c,ll)+1.d0
                cycle
              end if
              p1=spt%prot(c,ll,kd,1)+spt%prot(c,ll,kd,2)
              if ( p1 == 1.d0 ) then
                countH(ip,c) =countH(ip,c)+1.d0 ! marqueur informatif
                markH(ip,c,ll)= markH(ip,c,ll)+1.d0
                cycle
              end if
           end do !ll
          end do !!kd
        end do !!jm

        markH(ip,c,:) = markH(ip,c,:) / real(mHnb)
        countH(ip,c)=countH(ip,c)/real(cHnb)
      end do !!np
    end do !!c

   end subroutine count_nbheterozygote_mark_byChr


!!****f* m_qtlmap_haplotype_V2/pdegp
!! NAME
!!    pdegp
!! DESCRIPTION
!!
!! NOTES
!!   Determination de la phase la plus probable pour les peres
!!   Tableaux Arbitraire des 4 cas de transmission
!!
!!      Male \ Femelle   Chr1   Chr2
!!      Chr1              1       2
!!      Chr2              3       4
!!
!!      Utilisation des structures de donnees
!!      spt%prot(c,desc,ll,cas) : probabilite de transmission du marqueur ll pour le descenant desc
!!      du Cas 'cas' (1,2,3 ou 4)
!! SOURCE
      subroutine pdegp(dataset,spt)
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        type(PDD_BUILD)            ,intent(inout)         :: spt
! Divers
      integer geno,t,tder
      integer(kind=KIND_PHENO) :: m1,m2

      integer :: i,ii,imark,ip,j,jm,jmax1,jnd,kd,kkd
      integer :: linf,ll,ll1,ll2,lmark,lnd,alloc_stat,c
      integer :: trans(dataset%genea%nd,maxval(dataset%map%nmk)),hcon(maxval(dataset%map%nmk))
      integer :: indcon(maxval(dataset%map%nmk)),indinc(maxval(dataset%map%nmk))
      integer :: comptrans(maxval(dataset%map%nmk))
      integer :: nbcon,nbinc,nd1,nd2,ngeno,nhomo,nlvi,nm1,nm2
      integer :: nmark,nmkinf,nmknmk,nphp,value_alloc
      real (kind=dp) :: p1,p2,pmax,pr,probmax_t,sprob,xlvi
      logical :: hetero(maxval(dataset%map%nmk))
      !Les valeurs sont [1,2] un octet suffit
      integer(kind=1), dimension (:,:),allocatable :: h,tphp
      real (kind=dp) , dimension (:),allocatable ::prob,probphp
      type(GENEALOGY_BASE) , pointer :: dgenea
      type(GENOTYPE_BASE) , pointer :: dga

      dga => dataset%genoAnimal

      dgenea => dataset%genea

      call log_mess('calling...PDEGP_V4',DEBUG_DEF)
!
! Etablissement de la phase la plus probable
!
    do c=1,dataset%map%nchr
      do 300 ip=1,dgenea%np
        nm1=dgenea%nmp(ip)+1
        nm2=dgenea%nmp(ip+1)
!
!  on cree le tableau hetero qui liste l'heterozygotie des marqueurs
!
	    nmark=dataset%map%nmk(c)
        nhomo=0

        do imark=1,dataset%map%nmk(c)
          hetero(imark)=.true.
          if(dga%correp(ip).eq.9999) then
            m1=dga%nmanque
          else
            m1=dga%pheno(c,imark,dga%correp(ip),1)
            m2=dga%pheno(c,imark,dga%correp(ip),2)
          end if
          if(m1.eq.dga%nmanque) m2=m1
          if(m1.eq.m2) then
            nhomo=nhomo+1    ! le pere est homozygote au marqueur imark
            hetero(imark)=.false.
          end if
        end do
!
!  tableau de comptage des transmissions connues
!
        do ll=1,dataset%map%nmk(c)
	  comptrans(ll)=0
	end do
!
!  creation du tableau des transmissions observees
!
        kkd=0
        do jm=nm1,nm2
          nd1=dgenea%ndm(jm)+1
          nd2=dgenea%ndm(jm+1)
          do kd=nd1,nd2
            kkd=kkd+1
            do ll=1,dataset%map%nmk(c)
              trans(kkd,ll)=0
              p1=spt%prot(c,ll,kd,1)+spt%prot(c,ll,kd,2)
              if(p1.eq.1.d0)trans(kkd,ll)=1 ! recu 1er Chr du pere
              p2=spt%prot(c,ll,kd,3)+spt%prot(c,ll,kd,4)
              if(p2.eq.1.d0)trans(kkd,ll)=2 ! recu 2eme Chr du pere
	      if(trans(kkd,ll).ne.0)comptrans(ll)=comptrans(ll)+1
            end do
          end do
        end do
!
!  comptage des marqueurs informatifs
!
        nmkinf=0
        do ll=1,dataset%map%nmk(c)
	  if(comptrans(ll).ne.0)nmkinf=nmkinf+1
	end do
!
!  creation du tableau h (anciennement dans combine)
!  en se restreignant aux phases utiles
!
!  on commence par lister les marqueurs dont l'origine est connue
!  soit par l'ascendance, soit parcequ'ils sont homozygotes
!  (dans ce cas on donne l'ordrep=12)
!

        nbcon=0
        nbinc=0
        do lmark=1,dataset%map%nmk(c)

          !1er allele paternel / 2eme allele maternel
          if(ordrep(c,ip,lmark).eq.12)then
            nbcon=nbcon+1
            hcon(nbcon)=1
            indcon(nbcon)=lmark
          !1er allele maternel / 2eme allele paternel
          else if(ordrep(c,ip,lmark).eq.21)then
            nbcon=nbcon+1
            hcon(nbcon)=2
            indcon(nbcon)=lmark
          !Le pere est homozygote sur ce marqueur
          else if(.not.hetero(lmark))then
            nbcon=nbcon+1
            hcon(nbcon)=1
            indcon(nbcon)=lmark
          else
            nbinc=nbinc+1
            indinc(nbinc)=lmark
          end if

        end do

        !--------------------- ALLOC DYNAMIQUE ---------------------
        !   ---- A FAIRE VALIDER PAR JM ----------------------------
        if (allocated(h)) deallocate(h)
        if (allocated(tphp)) deallocate(tphp)
        if (allocated(prob)) deallocate(prob)
        if (allocated(probphp)) deallocate(probphp)

        if (nbinc>0) then
           if ( (nbinc*2) > MAX_MARKER) then
             call stop_application("You can not used this version of haplotype calculation. Please use --snp option")
           end if

          value_alloc= (2**(nbinc-1))*2
          allocate (h(dataset%map%nmk(c),value_alloc),STAT=alloc_stat)
          allocate (tphp(dataset%map%nmk(c),value_alloc),STAT=alloc_stat)
          allocate (prob(value_alloc),STAT=alloc_stat)
          allocate (probphp(value_alloc),STAT=alloc_stat)
        else
          allocate (h(dataset%map%nmk(c),1),STAT=alloc_stat)
          allocate (tphp(dataset%map%nmk(c),1),STAT=alloc_stat)
          allocate (prob(1),STAT=alloc_stat)
          allocate (probphp(1),STAT=alloc_stat)
        end if

         if (alloc_stat/=0) then
          call stop_application('Not enough memory to compute '// &
      'transmission probabilities. Deacrease number of marker.')
         end if

        !--------------------------- END OFI ---------------------

        !***************** WARNING
        !***************** WARNING
        !***************** WARNING
        !***************** WARNING
        h = 1 ! AJOUT OFI, sinon il y a des cas ou h n est pas initialise et ca fait planter l application....
        !***************** WARNING
        !***************** WARNING
        !***************** WARNING
        !***************** WARNING
!
!  traitement des cas limites (nmk(c) = 0 ou 1 marqueur)
!
        if(dataset%map%nmk(c).eq.0.or.nmkinf.eq.0)then
!          call stop_application('Chromosome '//trim(chromo(c))//&
!             ' Sire '//trim(pere(ip))// &
!             ' had no informative progeny '//&
!             ' for any of the markers'//' It should not be'//&
!             ' included in the analysis nmk='//trim(str(nmk(c)))//&
!             ' nmkinf='//trim(str(nmkinf)))

! ***************** OFI septembre 2012
!        le pere est completement homozygote

          do i=1,dataset%map%nmk(c)
           spt%genotyp(c,i,dga%correp(ip),1)=dga%pheno(c,i,dga%correp(ip),1)
           spt%genotyp(c,i,dga%correp(ip),2)=dga%pheno(c,i,dga%correp(ip),2)
          end do
          !on passe au pere suivant
          cycle
        end if


    ! Il y a au moin 1 marqueur du pere informatif !

	if(dataset%map%nmk(c).eq.1.or.nmkinf.eq.1) then
	  ngeno=1
	  prob(1)=1.d0
	  if(nbcon.eq.1)then
	    h(1,1)=hcon(1)
	  else
	    h(1,1)=1
	  end if
	  go to 200
	end if


!
!  s'il reste des inconnues
!  on remplit le tableau h, en completant les elements connus issus
!  de indcon et hcon par les phases possibles des elements inconnus
!
        if(nbinc.ne.0) then
        nmknmk=2**(nbinc-1)
        ngeno=nmknmk

        do geno=1,nmknmk
          do lmark=1,nbcon
            h(indcon(lmark),geno)=hcon(lmark)
          end do

          do lmark=1,nbinc
            lnd=2**(lmark-1)
            jnd=floor(float(lnd-1+geno)/float(lnd))
            h(indinc(nbinc+1-lmark),geno)=1+mod(jnd+1,2)
          end do

        end do


!
!  quand le genotype est oriente (phasp(ip) e true), il faut aussi calculer les
!  probabilites des phases miroirs
!

      if(spt%phasp(c,ip)) then
        ngeno=ngeno+nmknmk
        do geno=1,nmknmk

          do lmark=1,nbcon
            h(indcon(lmark),nmknmk+geno)=hcon(lmark)
          end do

          do lmark=1,nbinc
            lnd=2**(lmark-1)
            jnd=floor(float(lnd-1+geno)/float(lnd))
            h(indinc(nbinc+1-lmark),nmknmk+geno)=2-mod(jnd+1,2)
          end do
        end do
      end if
!
!  cas ou tout est connu
!
      else
        ngeno=1
        do lmark=1,nbcon
          h(indcon(lmark),ngeno)=hcon(lmark)
        end do
      end if
!
!  traitement des cas limites (0 ou 1 marqueur heterozygote)
!

      if(dataset%map%nmk(c)-nhomo.le.1)then
        ngeno=1
        prob(1)=1.d0
        go to 200
      end if
      probmax_t=-1.d6
!
!
      do 100 geno=1,ngeno
!
!  calcul de la probabilite de la phase d'apres la descendance
!
            xlvi=0.d0
            kkd=0
            do jm=nm1,nm2
              nd1=dgenea%ndm(jm)+1
              nd2=dgenea%ndm(jm+1)
              nlvi=0
              do kd=nd1,nd2
                kkd=kkd+1
!
! Calcul de la probabilite du phenotype du descendant
                ll1=1
                do while (trans(kkd,ll1).eq.0.and.ll1.lt.dataset%map%nmk(c))
                  ll1=ll1+1
                end do
                tder=trans(kkd,ll1)
                if (tder.ne.0) then
                  if (h(ll1,geno).eq.2) tder=3-tder
                    linf=ll1+1
                    do ll2=linf,dataset%map%nmk(c)
                      t=trans(kkd,ll2)
                      pr=0.d0
                      if (t.ne.0) then
                        if (h(ll2,geno).eq.2) t=3-t
                        if (t.eq.tder) then
                          pr=1.d0-dataset%map%rm(c,ll1,ll2)
                        else
                          pr=dataset%map%rm(c,ll1,ll2)
                        end if
                        tder=t
                        ll1=ll2
                        xlvi=xlvi+dlog(pr)
                      end if
                    end do
                  end if

                end do
              end do

          prob(geno)=xlvi
          if(probmax_t.le.xlvi) probmax_t=xlvi

  100   continue

!
! Standardisation
!
      sprob=0.d0
      do geno=1,ngeno
        if (prob(geno).ne.0.d0) then
          prob(geno)=dexp(prob(geno)-probmax_t)
          sprob=sprob+prob(geno)
        end if
      end do

      !AJOUT OFI : sprob peut valoir 0
      if ( sprob /= 0 ) then
       do geno=1,ngeno
         prob(geno)=prob(geno)/sprob
       end do
      end if
!
! elimination des phases improbables
!
  200 continue
            nphp=0
            do geno=1,ngeno
              if(prob(geno).ge.dataset%params%PHPSEUIL)then
                nphp=nphp+1
                 do imark=1,nmark
                  tphp(imark,nphp)=h(imark,geno)
                end do
                probphp(nphp)=prob(geno)
              end if
	    end do


            if(nphp.eq.0) then
              call stop_application(' Sire '//trim(dgenea%pere(ip))// &
              ' had no phase with a probability  above the threshold '//trim(str(dataset%params%PHPSEUIL))// &
              ' The model is not appropriate to your data  We suggest the use of the linealised model')
            end if

!
!  cas des inidvidus totalement homozygotes
!

          if(.not. spt%phasp(c,ip).and.nhomo.eq.dataset%map%nmk(c))then
            prob(1)=1.d0
            do i=2,ngeno
              prob(i)=0.d0
            end do
          end if

!
! Recherche de la phase la plus probable
        pmax=0.d0
        do geno=1,nphp
          if(pmax.lt.probphp(geno)) then
            pmax=probphp(geno)
            jmax1=geno
          end if
        end do
!
!  stockage du genotype
!
        do i=1,dataset%map%nmk(c)
          spt%genotyp(c,i,dga%correp(ip),1)=dga%pheno(c,i,dga%correp(ip),tphp(i,jmax1))
          spt%genotyp(c,i,dga%correp(ip),2)=dga%pheno(c,i,dga%correp(ip),3-tphp(i,jmax1))
        end do
!
! Reorganisation du tableau des probabilites de transmission
          do jm=nm1,nm2
          nd1=dgenea%ndm(jm)+1
          nd2=dgenea%ndm(jm+1)
          do kd=nd1,nd2
              do ll=1,dataset%map%nmk(c)
                if(tphp(ll,jmax1).eq.2) then
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
  300 continue
    end do

      if (allocated(h)) deallocate (h)
      if (allocated(tphp)) deallocate(tphp)
      if (allocated(prob)) deallocate (prob)
      if (allocated(probphp)) deallocate (probphp)
      call log_mess('end PDEGP_V4',DEBUG_DEF)

      end subroutine pdegp
!!***

!!****f* m_qtlmap_haplotype_V2/pdegm
!! NAME
!!    pdegm
!! DESCRIPTION
!!
!! NOTES
!!   Determination des probabilites des phases pour les meres
!! SOURCE
      subroutine pdegm(dataset,spt)
      type(QTLMAP_DATASET)       ,intent(in)            :: dataset
      type(PDD_BUILD)            ,intent(inout)         :: spt
!
! Divers
      integer   ,dimension(:,:) ,allocatable               :: ngend_t,ndesc_t
      real (kind=dp),dimension(:,:),allocatable            :: probg_t
      integer geno,t1,t2,t
      real (kind=dp) ::  trans(4,4),s(4),q(4)
      integer(kind=KIND_PHENO) ::  m1,m2
      integer i,ii,icon,ifem,igeno,imark,j,jj,jm,jnd,kd,c
      integer ld,ll,ll1,ll2,lmark,lnd,nbcon,nbinc,nd1,nd2
      integer ngeno,ngeno1,ngeno2,nhomo,nmark,alloc_stat,value_alloc
      integer nmknmk,nombfem,nphm,ndf(dataset%genea%nm+1),repdes(dataset%genea%nm,dataset%genea%nd),corref(dataset%genea%nm)
      integer , dimension(maxval(dataset%map%nmk)+1) :: hcon,indcon,indinc
      real (kind=dp) :: pmax,pr,probmax_t,r1,r2,sprob,sprobphm,xlvi
      logical  :: hetero(maxval(dataset%map%nmk)+1)

      integer , dimension(:,:),allocatable ::tphm,h
      real (kind=dp) , dimension(:),allocatable :: probphm,prob
      real (kind=dp) , dimension(:,:,:),allocatable ::p
      real (kind=dp) , dimension(:,:,:,:),allocatable ::ptfin_t
      type(GENOTYPE_BASE) , pointer :: dga
      type(GENEALOGY_BASE) , pointer :: dgenea

      call log_mess('calling...PDEGM_V4',DEBUG_DEF)
      allocate (p(maxval(dataset%map%nmk),dataset%genea%nd,4),stat=alloc_stat)

      if (alloc_stat/=0) then
          call stop_application('Not enough memory to compute transmission probabilities. Deacrease number of marker.')
      end if

      dga => dataset%genoAnimal
      dgenea => dataset%genea
!
!******************************************************************************
!
!
!  Creation des tableaux
!    ndf (nombre total de descendants/femelle)
!    repdes (identification des descendants de chaque femelles selon son nuero d'ordre
!  en lecture)
!    corref (l'equivalent de correm pour les femelles)
!
!  cette partie sera plutot e positionner dans lect_genea et lect_typ
!
      do ifem=1,dgenea%nfem
        ndf(ifem)=0
      end do

      do jm=1,dgenea%nm
        nd1=dgenea%ndm(jm)+1
        nd2=dgenea%ndm(jm+1)
        ifem=dgenea%repfem(jm)
        do kd=nd1,nd2
          ndf(ifem)=ndf(ifem)+1
          repdes(ifem,ndf(ifem))=kd
        end do

        corref(ifem)=dga%correm(jm)
      end do

!
!  le tableau force indique si une femelle n'a pas assez de descendants
!  pour etre estimee
!
  !    do jm=1,nm
   !     force(jm)=.false.
   !   end do
!
      allocate ( spt%ngend(dataset%map%nchr,1) )
      allocate ( spt%probg(dataset%map%nchr,1) )
      allocate ( spt%ndesc(dataset%map%nchr,1) )

      spt%ngenom(:,1)=0
      spt%ngend(:,1)=0
      spt%probg=0.d0
!
!******************************************************************************
!
!  on traite les femelles (et non les meres) une e une
!
!  nombfem est le nombre de femelles dont on peut chercher les phase
!
    do c=1,dataset%map%nchr
      nombfem=0
      do 500 jm=1,dgenea%nm
        nd1=dgenea%ndm(jm)+1
        nd2=dgenea%ndm(jm+1)
        ifem=dgenea%repfem(jm)
!
!******************************************************************************
!
! Si le nombre de pleins freres est trop faible, le genotype de la mere ne sera pas considere
! Dans ce cas le premier genotype rencontre possible est considere comme le bon

        !if(ndf(ifem).lt.ndmin) then
         ! force(jm)=.true.
       if(.not.dga%estfem(ifem))then

          if (allocated(h)) deallocate(h)
          if (allocated(tphm)) deallocate(tphm)
          if (allocated(prob)) deallocate(prob)
          if (allocated(probphm)) deallocate(probphm)
          allocate (h(dataset%map%nmk(c),1))
          allocate (tphm(dataset%map%nmk(c),1))
          allocate (prob(1))
          allocate (probphm(1))

!  on cherche le premier genotype possible
!
          do nmark=1,dataset%map%nmk(c)
            tphm(nmark,1)=1
            if(ordref(c,ifem,nmark).eq.21)tphm(nmark,1)=2
          end do
          nphm=1

          probphm(1)=1.d0
          sprobphm=probphm(1)

          go to 400

        end if
!
!******************************************************************************
!
!
!  on cree le tableau hetero qui liste l'heterozygotie des marqueurs
!

        nhomo=0
        do nmark=1,dataset%map%nmk(c)
          hetero(nmark)=.true.
          if(corref(ifem).eq.9999) then
            m1=dga%nmanque
          else
            m1=dga%pheno(c,nmark,corref(ifem),1)
            m2=dga%pheno(c,nmark,corref(ifem),2)
          end if

          if(m1.eq.dga%nmanque) m2=m1
          if(m1.eq.m2) then
            nhomo=nhomo+1
            hetero(nmark)=.false.
          end if
        end do
!
!  creation du tableau h (anciennement dans combine)
!  en se restreignant aux phases utiles
!
!  on commence par lister les marqueurs dont l'origine est connue
!  soit par l'ascendance, soit parcequ'ils sont homozygotes
!  (dans ce cas on donne l'ordrep=12)
!

        nbcon=0
        nbinc=0
        do lmark=1,dataset%map%nmk(c)

          icon=0
          if(ordref(c,ifem,lmark).eq.12)then
            nbcon=nbcon+1
            hcon(nbcon)=1
            indcon(nbcon)=lmark
            icon=1
          end if

          if(ordref(c,ifem,lmark).eq.21)then
            nbcon=nbcon+1
            hcon(nbcon)=2
            indcon(nbcon)=lmark
            icon=1
          end if

          if(.not.hetero(lmark))then
            nbcon=nbcon+1
            hcon(nbcon)=1
            indcon(nbcon)=lmark
            icon=1
          end if

          if(icon.eq.0) then
            nbinc=nbinc+1
            indinc(nbinc)=lmark
          end if

        end do

         !--------------------- ALLOC DYNAMIQUE ---------------------
        !   ---- A FAIRE VALIDER PAR JM ----------------------------
        if (allocated(h)) deallocate(h)
        if (allocated(tphm)) deallocate(tphm)
        if (allocated(prob)) deallocate(prob)
        if (allocated(probphm)) deallocate(probphm)

        if (nbinc>0) then
          value_alloc= (2**(nbinc-1))*2
          allocate (h(dataset%map%nmk(c),value_alloc),STAT=alloc_stat)
          allocate (tphm(dataset%map%nmk(c),value_alloc),STAT=alloc_stat)
          allocate (prob(value_alloc),STAT=alloc_stat)
          allocate (probphm(value_alloc),STAT=alloc_stat)
        else
          allocate (h(dataset%map%nmk(c),1),STAT=alloc_stat)
          allocate (tphm(dataset%map%nmk(c),1),STAT=alloc_stat)
          allocate (prob(1),STAT=alloc_stat)
          allocate (probphm(1),STAT=alloc_stat)
        end if

         if (alloc_stat/=0) then
          call stop_application('Not enough memory to compute transmission probabilities. Deacrease number of marker.')
         end if

        !--------------------------- END OFI ---------------------
!
!  s'il reste des inconnues
!  on remplit le tableau h, en completant les elements connus issus
!  de indcon et hcon par les phases possibles des elements inconnus
!
        if(nbinc.ne.0) then
          nmknmk=2**(nbinc-1)
          ngeno=nmknmk
          do geno=1,nmknmk

            do lmark=1,nbcon
              h(indcon(lmark),geno)=hcon(lmark)
            end do

            do lmark=1,nbinc
              lnd=2**(lmark-1)
              jnd=floor(float(lnd-1+geno)/float(lnd))
              h(indinc(nbinc+1-lmark),geno)=1+mod(jnd+1,2)
            end do

          end do
!
!  quand le genotype est oriente (phasm(c,jm) e true), il faut aussi calculer les
!  probabilites des phases miroirs
!

          if(spt%phasm(c,jm)) then
            ngeno=ngeno+nmknmk
            do geno=1,nmknmk

              do lmark=1,nbcon
                h(indcon(lmark),nmknmk+geno)=hcon(lmark)
              end do

              do lmark=1,nbinc
                lnd=2**(lmark-1)
                jnd=floor(float(lnd-1+geno)/float(lnd))
                h(indinc(nbinc+1-lmark),nmknmk+geno)=2-mod(jnd+1,2)
              end do

            end do
          end if
!
!  cas ou tout est connu
!
        else
          ngeno=1
          do lmark=1,nbcon
            h(indcon(lmark),ngeno)=hcon(lmark)
          end do
        end if

!
!  on va calculer, pour les ngeno configurations possibles
!  la probabilite de la phase d'apres la descendance
!
        probmax_t=-1.d6
!
        do 200 geno=1,ngeno
!
          xlvi=0.d0
          do 100 ld=1,ndf(ifem)
            kd=repdes(ifem,ld)
!
! Reorganisation du tableau des probabilites de transmission

            do ll=1,dataset%map%nmk(c)
              do i=1,4
                p(ll,kd,i)=spt%prot(c,ll,kd,i)
              end do
              if(h(ll,geno).eq.2) then
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
            xlvi=xlvi+log(pr)

  100     continue
          prob(geno)=xlvi
          if(probmax_t.le.xlvi) probmax_t=xlvi

  200   continue
!
!******************************************************************************

!
!  on rescale les prob(geno) en retranchant pmax e toutes les valeurs
!  ce qui fait que la plus grande vaut 0 et les autres sont augmentees
!  de -pmax.  On prend l'exponentielle de ce resultat
!
      sprob=0.d0
      do geno=1,ngeno
        if (prob(geno).ne.0.d0) then
          prob(geno)=dexp(prob(geno)-probmax_t)
          sprob=sprob+prob(geno)
        end if
      end do

       do geno=1,ngeno
         prob(geno)=prob(geno)/sprob
       end do
!
! elimination des phases improbables
!
       nphm=0
       sprobphm=0.d0
       do geno=1,ngeno
         if(prob(geno).ge.dataset%params%PRSEUIL)then
           nphm=nphm+1
           do imark=1,dataset%map%nmk(c)
             tphm(imark,nphm)=h(imark,geno)
           end do
           probphm(nphm)=prob(geno)
           sprobphm=sprobphm+prob(geno)
         end if
       end do

       do geno=1,nphm
           probphm(geno)=probphm(geno)/sprobphm
       end do
  300  continue

!******************************************************************************
!
! Impression des resultats
!
  400 continue

   if(nphm.eq.0)then
         call stop_application('Dam '//trim(dgenea%femelle(ifem))//' has no possible haplotype  '//   &
          '(none has probability greater than '// trim(str(dataset%params%PRSEUIL))//') : '//  &
          'check the genotype compatibility in the all sire family')
   end if

 !  if (corref(ifem).ne.9999) then
 !        if (ndf(ifem).ge.ndmin) then
 !           nombfem=nombfem+1
 !        end if
 !  end if
!
!   stockage des resultats par mere (et non par femelle)
!    et reperage de la phase la plus probable
!

        pmax=0.d0
        allocate ( probg_t(dataset%map%nchr,size(spt%probg,2)) )
        probg_t(:,1:size(spt%probg,2)) = spt%probg

        deallocate ( spt%probg )
        allocate ( spt%probg(dataset%map%nchr,size(probg_t,2)+nphm) )
        spt%probg=0.d0
        spt%probg(:,1:size(probg_t,2))  = probg_t

        deallocate(probg_t)

        spt%ngenom(c,jm+1)=spt%ngenom(c,jm)
        do geno=1,nphm
          spt%ngenom(c,jm+1)=spt%ngenom(c,jm+1)+1
!
!  ajout de genotypm
          do ii=1,dataset%map%nmk(c)
            if (dga%correm(jm) <= size(spt%genotyp,3).and. corref(ifem) <= size(dga%pheno,3)) then
              spt%genotypm(c,ii,spt%ngenom(c,jm+1),1)=dga%pheno(c,ii,corref(ifem),tphm(ii,geno))
              spt%genotypm(c,ii,spt%ngenom(c,jm+1),2)=dga%pheno(c,ii,corref(ifem),3-tphm(ii,geno))
            end if
          end do
!
!
          spt%probg(c,spt%ngenom(c,jm+1))=probphm(geno)
          end do
          
          ngeno1=spt%ngenom(c,jm)+1
          ngeno2=spt%ngenom(c,jm+1)
          igeno=0

          allocate ( ngend_t(dataset%map%nchr,size(spt%ngend,2)) )
          ngend_t(:,1:size( spt%ngend,2)) =  spt%ngend
          deallocate ( spt%ngend )
          allocate ( spt%ngend(dataset%map%nchr,size(ngend_t,2)+(ngeno2-ngeno1)+1) )
           spt%ngend=0.d0
          spt%ngend(:,1:size(ngend_t,2))=ngend_t
          deallocate(ngend_t)

          if ( size( ptfin,3) < maxval( spt%ngend)+(ngeno2-ngeno1+1)*(nd2-nd1+1) ) then
            allocate ( ptfin_t(dataset%map%nchr,maxval(dataset%map%nmk),size(ptfin,3),4))
            ptfin_t= ptfin
            deallocate ( ptfin)
            !pour eviter de reallouer sans cesse et faire des copies de tableau, on alloue nd valeurs supplementaires
            !normalement ca devrait etre : maxval(ngend)+(ngeno2-ngeno1+1)*(nd2-nd1+1)
            allocate (ptfin(dataset%map%nchr,maxval(dataset%map%nmk),maxval(spt%ngend)+dgenea%nd,4))
            ptfin(:,:,1:size(ptfin_t,3),:)=ptfin_t
            deallocate (ptfin_t)
          end if

          do geno=ngeno1,ngeno2
            igeno=igeno+1
            spt%ngend(c,geno+1)=spt%ngend(c,geno)
            allocate ( ndesc_t(dataset%map%nchr,size(spt%ndesc,2)) )
            ndesc_t(:,1:size(spt%ndesc,2)) = spt%ndesc
            deallocate ( spt%ndesc )
            allocate ( spt%ndesc(dataset%map%nchr,size(ndesc_t,2)+(nd2-nd1)+1) )
            spt%ndesc=0
            spt%ndesc(:,1:size(ndesc_t,2))=ndesc_t
            deallocate( ndesc_t )
            do kd=nd1,nd2
              spt%ngend(c,geno+1)=spt%ngend(c,geno+1)+1
              spt%ndesc(c,spt%ngend(c,geno+1))=kd
              do ll=1,dataset%map%nmk(c)
                do i=1,4
                  ptfin(c,ll,spt%ngend(c,geno+1),i)=spt%prot(c,ll,kd,i)
                end do
                if(tphm(ll,igeno).eq.2) then
                  ptfin(c,ll,spt%ngend(c,geno+1),1)=spt%prot(c,ll,kd,2)
                  ptfin(c,ll,spt%ngend(c,geno+1),2)=spt%prot(c,ll,kd,1)
                  ptfin(c,ll,spt%ngend(c,geno+1),3)=spt%prot(c,ll,kd,4)
                  ptfin(c,ll,spt%ngend(c,geno+1),4)=spt%prot(c,ll,kd,3)
                end if

              end do
            end do
          end do
  500 continue
    end do

      deallocate (p)
      if (allocated(h)) deallocate(h)
      if (allocated(tphm)) deallocate(tphm)
      if (allocated(prob)) deallocate(prob)
      if (allocated(probphm)) deallocate(probphm)
      call log_mess('end PDEGM_V4',DEBUG_DEF)
      end subroutine pdegm
!!***

!!****f* m_qtlmap_haplotype_V2/pded
!! NAME
!!    pded
!! DESCRIPTION
!!
!! NOTES
!!   Calcul des probabilites de transmission le long du chromosome
!! SOURCE
      subroutine pded(dataset,spt)
       type(QTLMAP_DATASET)       ,intent(in)            :: dataset
       type(PDD_BUILD)            ,intent(inout)         :: spt

!
! Tableaux dimensionnes selon l, le nombre de marqueurs
      integer :: cas(maxval(dataset%map%nmk)),trans(maxval(dataset%map%nmk)),infor(maxval(dataset%map%nmk))
      ! dim : nbt*4+nbt,nbt
      integer  , dimension(:,:),allocatable :: last
      real (kind=dp), dimension(:,:),allocatable :: pbt

      integer :: index_mark(maxval(dataset%map%nmk))
      logical :: index_mark_valide(maxval(dataset%map%nmk))

!
! Divers
      integer(kind=1):: new(11,8)
      integer geno,npas,removeCas4,maxcas4,removeCas2,maxCas2
      integer ilong,ip,it,jm,jt,rc4,rc2,lit,lkkk,c
      integer n,kd,kt,ii
      integer mark_selec,im1,im2
      integer linf,lk,llk,lmax,lq,lt,maxCas
      integer nbt,nd1,nd2,ngeno1,ngeno2,nm1,nm2,nmkinf
      real (kind=dp) pos_m,posi_lk,posi_lkinf
      real (kind=dp) dg,dd,ddf,ddm,dgf,dgm
      real (kind=dp) dx,pbtt,val_min
      real (kind=dp) recm,recf,recmlq,recflq,recmqr,recfqr
      real (kind=dp) xint,xintm,xintf

      type(GENEALOGY_BASE) , pointer :: dgenea

      dgenea => dataset%genea

      !external  xaldane
      !double precision xaldane
      call log_mess('calling...PDED_V4',DEBUG_DEF)
!
!  tableau des transmission
!
      data new / 1,1,2,2,1,1,1,2,1,1,1,   &
                1,2,1,2,1,2,1,1,1,2,1,    &
                0,0,0,0,2,2,1,2,2,2,1,    &
                0,0,0,0,2,1,2,2,1,2,2,    &
                0,0,0,0,0,0,0,0,0,0,2,    &
                0,0,0,0,0,0,0,0,0,0,1,    &
                0,0,0,0,0,0,0,0,0,0,2,    &
                0,0,0,0,0,0,0,0,0,0,2/


      ! init pour la suite
      allocate (last(2**(((MAX_CAS_ALLOC)/2)+1),2))
      allocate (pbt(4,2**(((MAX_CAS_ALLOC)/2)+1)))
   do c=1,dataset%map%nchr
!
! Taille du segment explore
      ilong=dataset%map%get_ilong(c)
      infor = 0
!
! Calcul pour chaque descendant de la probabilite de transmission
       do 997 ip=1,dgenea%np
       !print *,'***************************PERE:',trim(pere(ip))
       nm1=dgenea%nmp(ip)+1
       nm2=dgenea%nmp(ip+1)
       do 998 jm=nm1,nm2
       !print *,'**********************MERE:',trim(mere(jm))
       ngeno1=spt%ngenom(c,jm)+1
       ngeno2=spt%ngenom(c,jm+1)
       do 999 geno=ngeno1,ngeno2
       !  print *,'**********************GENO:',geno
       nd1=spt%ngend(c,geno)+1
       nd2=spt%ngend(c,geno+1)
        !print *,'ND1:',nd1,' ND2:',nd2
       do 1000 kd=nd1,nd2
!
!  le vecteur trans(lk) indique l'information jointe sur la transmission au locus lk
!  le vecteur infor contient la liste des marqeurs doublement informatifs
!
!  l'origine grand parentale des alleles recu par le descendant k au marqueur lk
!  depend de trans (tableau new)
!  soit 1 pour grand pere et 2 pour grand mere
!
!	           cas 1	   cas 2	   cas 3	   cas 4
! trans(lk)	pere  mere	pere  mere	pere  mere	pere  mere
!  	1	 1	1
!	2	 1	2
!	3	 2	1
!	4	 2	2
!	5	 1	1	 2	2
!	6	 1	2	 2	1
!	7	 1	1	 1	2
!	8	 2	1	 2	2
!	9	 1	1	 2	1
!	10	 1	2	 2	2
!	11	 1	1	 1	2	 2	1	 2	2
!
!
       llk=0
       do lk=1,dataset%map%nmk(c)

!
!  cas(lk) est le nombre d'evenements de transmission possible pour le marqueur lk
!  cas(lk) vaut 1 pour trans(lk) = 1 e 4,
!		2 pour trans(lk) = 5 e 10,
!		4 pour trans(lk) = 11
!
         cas(lk)=2
!
! incertitude totale (donnees manquantes)
!
         trans(lk)=11
!
!  triplet heterozygotes
!
         if(ptfin(c,lk,kd,2).eq.0.and.ptfin(c,lk,kd,3).eq.0)trans(lk)=5
         if(ptfin(c,lk,kd,1).eq.0.and.ptfin(c,lk,kd,4).eq.0)trans(lk)=6
!
!  incertitude sur un parent
!
         if(ptfin(c,lk,kd,3).eq.0.and.ptfin(c,lk,kd,4).eq.0)trans(lk)=7
         if(ptfin(c,lk,kd,1).eq.0.and.ptfin(c,lk,kd,2).eq.0)trans(lk)=8
         if(ptfin(c,lk,kd,2).eq.0.and.ptfin(c,lk,kd,4).eq.0)trans(lk)=9
         if(ptfin(c,lk,kd,1).eq.0.and.ptfin(c,lk,kd,3).eq.0)trans(lk)=10
!
!  cas certains
!
         if(ptfin(c,lk,kd,1).eq.1.d0) trans(lk)=1
         if(ptfin(c,lk,kd,2).eq.1.d0) trans(lk)=2
         if(ptfin(c,lk,kd,3).eq.1.d0) trans(lk)=3
         if(ptfin(c,lk,kd,4).eq.1.d0) trans(lk)=4

         if(trans(lk).le.4) then
           llk=llk+1
           infor(llk)=lk
           cas(lk)=1
         end if

!
!  reajustement de cas
!
         if(trans(lk).eq.11) cas(lk)=4
       end do

       nmkinf=llk

!
! exploration du groupe de liaison pour la recherche des transmissions possibles
! le groupe de liaison est separe en segments flanques de MIF
! avec deux segments externes pouvant ne comporter qu'un seul Marqueur Informatif
! e droite pour celui de gauche et e gauche pour celui de droite
!
       n=1
       linf=1
  100  dx=dataset%map%absi(c,n)

!
!  cas des descendants sans information
!
       if(nmkinf.eq.0) then
         linf=0
         do ii=1,4
           spt%pdd(c,kd,ii,n)=0.25d0
         end do
         go to 125
       end if


       if(dx.le.dataset%map%posi(c,infor(linf)).or.linf.eq.nmkinf) go to 105
       linf=linf+1
       go to 100

  105  continue
!
!  pour des problemes d'arrondis, la position testee pour la qtl
!  peut depasser le bout du groupe de liaison... on l'y ramene
!
       if(dx.ge.dataset%map%posi(c,dataset%map%nmk(c)))dx=dataset%map%posi(c,dataset%map%nmk(c))
!
!  cas des positions sur des marqueurs informatifs
!
       if(dx.eq.dataset%map%posi(c,infor(linf))) then

         nbt=1

         do lq=1,4
	      pbt(lq,1)=0.d0
         end do

         pbt(trans(infor(linf)),1)=1.d0
         go to 120

       end if
!
!  cas des positions entre des marqueurs informatifs
!  trois situations
!	si dx est avant le premier MIF, on explore les transmissions entre le premier
!         marqueur et le premier MIF
!	si dx est entre deux MIF (cas standard)
!	si dx est apres le dernier MIF, on explore les transmissions entre le dernier MIF et le
!         dernier marqueur
!
       if(dx.lt.dataset%map%posi(c,infor(1))) then
         lk=1
	     lmax=infor(1)
       end if
       if(dx.gt.dataset%map%posi(c,infor(1)).and.dx.lt.dataset%map%posi(c,infor(nmkinf))) then
         lk=infor(linf-1)
	     lmax=infor(linf)
       end if
       if(dx.ge.dataset%map%posi(c,infor(nmkinf))) then
         lk=infor(nmkinf)
	     lmax=dataset%map%nmk(c)
       end if

       ! ---------- AJOUT OFI ALLOCATION STRUCTURES TEMPORAIRES ----------
       !taille des tableaux Max : cas(linf)*cas(linf+1)*..*cas(lmax)=nbt
       maxCas = 0
       maxCas = MAX_CAS_ALLOC
       ! On recherche un echantillonage
       removeCas2=0
       removeCas4=0
       maxcas2=0
       maxcas4=0

       ! print *,'LK:',lk,'lMAX:',lmax
       do lq=lk,lmax
         !  print *,lq,':',cas(lq)
           if (cas(lq) == 2) maxcas2 = maxcas2 + 1
           if (cas(lq) == 4) maxcas4 = maxcas4 + 1
       end do
       !print *,'MAX CAS:',maxCas, ' nb4:',maxcas4,' nb2',maxcas2
       do while ( maxCas >= MAX_CAS_ALLOC)
         rc2=removeCas2
         rc4=removeCas4
         maxCas=0
         do lq=lk,lmax
           if (cas(lq) > 1) then
             if (cas(lq)==4 .and. rc4>0) then
                rc4 = rc4-1
                cycle
             end if
              if (cas(lq)==2 .and. rc2>0) then
                rc2 = rc2-1
                cycle
             end if
             maxCas = maxCas + cas(lq)
           end if
         end do

         if (maxCas>=MAX_CAS_ALLOC .and. removeCas4<maxcas4) then
            removeCas4 = removeCas4 +1
         else if (maxCas>=MAX_CAS_ALLOC .and. removeCas2<maxcas2) then
            removeCas2 = removeCas2 +1
         else if (maxCas>=MAX_CAS_ALLOC) then
            call stop_application("-- error dev --")
         end if
       end do
      ! print *,'MAX CAS:',maxCas, ' r4:',removeCas4,' r2',removeCas2

       index_mark_valide = .false.
       do lq=lk,lmax
         index_mark(lq-lk+1)=lq
         index_mark_valide(lq-lk+1) = .true.
       end do
       lq=lk

       do while ( removeCas4>0 )
          val_min = 0.d0
          mark_selec = -1
         ! print *,'lk:',lk
          do  lit=lk,lmax
             if (index_mark_valide(lit-lk+1) .and. cas(lit)==4) then
               pos_m=abs(dx-dataset%map%posi(c,lit))
               !print *,'lit ',lit,' pos:',pos_m,' val_min:',val_min
               if (pos_m>=val_min) then
               !  print *,'min!!! lit:',lit
                 val_min = pos_m
                 mark_selec = lit
               end if
             end if
           end do

           if (mark_selec == -1) call stop_application("ERROR");
          ! print *,'remov4 :',mark_selec
           index_mark_valide(mark_selec-lk+1) = .false.
           removeCas4 = removeCas4 -1
       end do

       do while ( removeCas2>0 )
          val_min = 0.d0
          mark_selec = -1
          do  lit=lk,lmax
             if (index_mark_valide(lit-lk+1) .and. cas(lit)==2) then
               pos_m=abs(dx-dataset%map%posi(c,lit))
               if (pos_m>=val_min) then
                 val_min = pos_m
                 mark_selec = lit
               end if
             end if
           end do

           if (mark_selec == -1) call stop_application("ERROR");
          ! print *,'remov2 :',mark_selec
           index_mark_valide(mark_selec-lk+1) = .false.
           removeCas2 = removeCas2 -1
       end do


       lkkk=1
       do lq=lk,lmax
          if ( index_mark_valide(lq-lk+1)) then
            index_mark(lkkk)=lq
            !print *,lkkk,':',lq,':valide  pos:',posi(c,lq)
            lkkk=lkkk+1
          !else
            !print *,lq,':non valide pos:',posi(c,lq)
          end if
       end do
       !stop

       lk=1
       lmax=lkkk-1

!
!  nbt est le nombre d'evenements de transmissions possibles entre le premier marqueur
!  et le premier MIF
!  pbt(lq,it) est la probabilite du it emme de ces evenements (PBT), quand la transmission
!  au qtl est lq (1 e 4, correspondant aux 4 transmission grand parentales)
!
!  on etablit ces donnees iterativement en partant de l'extremite du chromosome
!
!  Premier marqueur
!

      nbt=cas(index_mark(lk))

	  do it=1,nbt
	   last(it,1)=new(trans(index_mark(lk)),2*(it-1)+1)
	   last(it,2)=new(trans(index_mark(lk)),2*(it-1)+2)
	   do lq=1,4
	     pbt(lq,it)=1.d0/dble(cas(index_mark(lk)))
	   end do
	  end do

!
!  marqueurs suivants
!
  110  	 lk=lk+1

  	  if(lk.gt.lmax) go to 120
!
!  le qtl est entre lk-1 et lk, il faut considerer les 4 evenements de tranmissions
!  au qtl
!
           im1 = index_mark(lk-1)
           im2 = index_mark(lk)
	       if(dx.ge.dataset%map%posi(c,im1).and. dx.le.dataset%map%posi(c,im2)) then

              xint=dataset%map%posi(c,index_mark(lk))-dataset%map%posi(c,index_mark(lk-1))
              dg=dx-dataset%map%posi(c,index_mark(lk-1))
              dd=dataset%map%posi(c,index_mark(lk))-dx
              xintm=dataset%map%posim(c,index_mark(lk))-dataset%map%posim(c,index_mark(lk-1))
              dgm=dg*xintm/xint
              ddm=dd*xintm/xint
              xintf=dataset%map%posif(c,index_mark(lk))-dataset%map%posif(c,index_mark(lk-1))
              dgf=dg*xintf/xint
              ddf=dd*xintf/xint
           recmlq=xaldane(dgm)
           recflq=xaldane(dgf)
           recmqr=xaldane(ddm)
           recfqr=xaldane(ddf)

           do lq=1,4
	     do it=1,nbt
	       if(new(11,1+2*(lq-1)).eq.last(it,1)) then
	         pbt(lq,it)=pbt(lq,it)*(1.d0-recmlq)
	       else
	         pbt(lq,it)=pbt(lq,it)*recmlq
	       end if

	       if(new(11,2+2*(lq-1)).eq.last(it,2)) then
	         pbt(lq,it)=pbt(lq,it)*(1.d0-recflq)
	       else
	         pbt(lq,it)=pbt(lq,it)*recflq
	       end if
	     end do
	   end do

           do lq=1,4
	     do it=1,nbt
               last(it+nbt*(lq-1),1)=new(11,1+2*(lq-1))
               last(it+nbt*(lq-1),2)=new(11,2+2*(lq-1))
	     end do
	   end do

!
!  PBTQ est incremente pour les evenements de transmision entre le QTL et le marqueur lk
!

           do lq=1,4

             do jt=1,cas(index_mark(lk))-1
	       do it=1,nbt
	         pbt(lq,it+nbt*jt)=pbt(lq,it)
	       end do
             end do

             do jt=1,cas(index_mark(lk))
	       do it=1,nbt
	         kt=it+nbt*(jt-1)
	         lt=it+nbt*(lq-1)
	     if(new(trans(index_mark(lk)),1+2*(jt-1)).eq.last(lt,1)) then
	           pbt(lq,kt)=pbt(lq,kt)*(1.d0-recmqr)
	         else
	           pbt(lq,kt)=pbt(lq,kt)*recmqr
	         end if

	     if(new(trans(index_mark(lk)),2+2*(jt-1)).eq.last(lt,2)) then
	           pbt(lq,kt)=pbt(lq,kt)*(1.d0-recfqr)
	         else
	           pbt(lq,kt)=pbt(lq,kt)*recfqr
	         end if
	       end do
	     end do

	   end do

           do jt=1,cas(index_mark(lk))
	     do it=1,nbt
	       kt=it+nbt*(jt-1)
	       lt=it+nbt*(lq-1)
	       last(kt,1)=new(trans(index_mark(lk)),1+2*(jt-1))
	       last(kt,2)=new(trans(index_mark(lk)),2+2*(jt-1))
	     end do
	   end do

	   nbt=nbt*cas(index_mark(lk))
       end if

!
!  le qtl n'est pas entre lk-1 et lk
!
      posi_lkinf = dataset%map%posi(c,index_mark(lk-1))
      posi_lk =  dataset%map%posi(c,index_mark(lk))
	 if(.not.(dx.ge.posi_lkinf.and.dx.le.posi_lk)) then
!
!  incrementation de la probabilites des transmission de marqueurs PBT
!

           recm=xaldane(dataset%map%posim(c,index_mark(lk))-dataset%map%posim(c,index_mark(lk-1)))
           recf=xaldane(dataset%map%posif(c,index_mark(lk))-dataset%map%posif(c,index_mark(lk-1)))

           do lq=1,4
             do jt=1,cas(index_mark(lk))-1
	       do it=1,nbt
	         pbt(lq,it+nbt*jt)=pbt(lq,it)
               end do
	     end do

             do jt=1,cas(index_mark(lk))

	       do it=1,nbt
	         kt=it+nbt*(jt-1)
	     if(new(trans(index_mark(lk)),1+2*(jt-1)).eq.last(it,1)) then
	           pbt(lq,kt)=pbt(lq,kt)*(1.d0-recm)
	         else
	           pbt(lq,kt)=pbt(lq,kt)*recm
	         end if

	     if(new(trans(index_mark(lk)),2+2*(jt-1)).eq.last(it,2)) then
	           pbt(lq,kt)=pbt(lq,kt)*(1.d0-recf)
	         else
	           pbt(lq,kt)=pbt(lq,kt)*recf
	         end if
	       end do

	     end do
	   end do

           do jt=1,cas(index_mark(lk))
	     do it=1,nbt
	       kt=it+nbt*(jt-1)
	       !lt=it+nbt*(lq-1)
	       last(kt,1)=new(trans(index_mark(lk)),1+2*(jt-1))
	       last(kt,2)=new(trans(index_mark(lk)),2+2*(jt-1))
	     end do
	   end do

	   nbt=nbt*cas(index_mark(lk))
	  end if

	  go to 110

  120    continue
!
!  on stocke la probabilite de transmission e la position dx
!
         pbtt=0.d0
         do lq=1,4
	   spt%pdd(c,kd,lq,n)=0.d0
           do kt=1,nbt
	     pbtt=pbtt+pbt(lq,kt)
	     spt%pdd(c,kd,lq,n)=spt%pdd(c,kd,lq,n)+pbt(lq,kt)
	       end do
	     end do
         if (pbtt /= 0) then
           do lq=1,4
	         spt%pdd(c,kd,lq,n)=spt%pdd(c,kd,lq,n)/pbtt
	       end do
	     end if

  125    continue
!
!  on avance d'un pas (en verifiant qu'on ne depasse pas le bout droit) et on recommence
!
         n=n+1
         npas=n*dataset%map%pas
	 if(npas.gt.ilong+1) go to 1000

	 if(linf.le.nmkinf) go to 100

!
!  bout du chromosome
!
 1000 continue
  999 continue
  998 continue
  997 continue
  end do

      if (allocated(last)) deallocate(last)
      if (allocated(pbt)) deallocate(pbt)
      call log_mess('end PDED_V4',DEBUG_DEF)

      end subroutine pded
!!***

!!****f* m_qtlmap_haplotype_V2/pdegp_snp
!! NAME
!!    pdegp_snp
!! DESCRIPTION
!!
!! NOTES
!!   Determination de la phase la plus probable pour les peres
!! SOURCE
      subroutine pdegp_snp(dataset,spt)
       type(QTLMAP_DATASET)       ,intent(in)            :: dataset
       type(PDD_BUILD)            ,intent(inout)         :: spt
! Divers
      integer(kind=KIND_PHENO) ::  m1,m2

      integer ii,imark,ip,indmin,jm,kd,kkd,c
      integer ll,ll1,ll2,lmark
      integer trans(dataset%genea%nd,maxval(dataset%map%nmk)),hcon(maxval(dataset%map%nmk)),indcon(maxval(dataset%map%nmk))
      integer oppos,phase,php(maxval(dataset%map%nmk))
      integer nbcon,nd1,nd2,nm1,nm2
      real (kind=dp) :: p1,p2,ratio,ratmin
      logical :: hetero(maxval(dataset%map%nmk)), mafor(maxval(dataset%map%nmk))
      type(GENEALOGY_BASE) , pointer :: dgenea
      type(GENOTYPE_BASE) , pointer :: dga

      dga => dataset%genoAnimal
      dgenea => dataset%genea

      call log_mess('calling...PDEGP_FAST',DEBUG_DEF)

!
! Etablissement de la phase la plus probable
!
    do c=1,dataset%map%nchr
      do 300 ip=1,dgenea%np
        nm1=dgenea%nmp(ip)+1
        nm2=dgenea%nmp(ip+1)
!
!  on cree le tableau hetero qui liste l'heterozygotie des marqueurs
!
        do imark=1,dataset%map%nmk(c)
          hetero(imark)=.true.
          if(dga%correp(ip).eq.9999) then
            m1=dga%nmanque
          else
            m1=dga%pheno(c,imark,dga%correp(ip),1)
            m2=dga%pheno(c,imark,dga%correp(ip),2)
          end if
          if(m1.eq.dga%nmanque) m2=m1
          if(m1.eq.m2) then
            hetero(imark)=.false.
          end if
        end do
!
!  on reduit le jeu de marqueurs e ceux pour lesquels le pere est heterozygote
!
!  creation du tableau des transmissions observees
!
        kkd=0
        do jm=nm1,nm2
          nd1=dgenea%ndm(jm)+1
          nd2=dgenea%ndm(jm+1)
          do kd=nd1,nd2
            kkd=kkd+1
            do ll=1,dataset%map%nmk(c)
              trans(kkd,ll)=0
              p1=spt%prot(c,ll,kd,1)+spt%prot(c,ll,kd,2)
              if(p1.eq.1.d0)trans(kkd,ll)=1
              p2=spt%prot(c,ll,kd,3)+spt%prot(c,ll,kd,4)
              if(p2.eq.1.d0)trans(kkd,ll)=2
            end do
          end do
        end do
!
!  Liste des marqueurs dont l'origine est connue par l'ascendance
!
        nbcon=0
        do lmark=1,dataset%map%nmk(c)

	  mafor(lmark)=.false.
	  hcon(lmark)=0
	  php(lmark)=0

          if(ordrep(c,ip,lmark).eq.12)then
            nbcon=nbcon+1
            hcon(lmark)=1
            indcon(nbcon)=lmark
	    mafor(lmark)=.true.
          end if

          if(ordrep(c,ip,lmark).eq.21)then
            nbcon=nbcon+1
            hcon(lmark)=2
            indcon(nbcon)=lmark
	    mafor(lmark)=.true.
          end if

        end do
!
!  reconstruction de la phase la plus probable par approximations successives
!  en partant du premier marqueur heterozygote chez le pere
!
      do ll1=1,dataset%map%nmk(c)
       if(hetero(ll1)) go to 10
      end do
!
!  le pere est entierement homozygote
!
      do lmark=1,dataset%map%nmk(c)
        php(lmark)=1
      end do

      go to 100
!
  10  php(ll1)=1
      if(nbcon.ne.0)php(ll1)=hcon(indcon(1))
!
!  si le premier marqueur informatif n'est pas le premier, on pose que
!  les  phases locales e sa gauche sont identiques
!
	 if(ll1.gt.1) then
	   do ll=ll1-1,1,-1
	       php(ll)=php(ll1)
	   end do
         end if
!
!  on part du premier marqueur informatif (s'il existe) et on construit
!  les phases locales en glissant vers la droite
!  ll1 et ll2 etant les indices de 2 marqueurs heterozygotes successifs
!  (tous les marqueurs intermediaires sont homozygotes), on compare
!  les transmission en ll1 et ll2, et on en deduit la phase locale
!
!  en cas d'incoherence entre cette reconstruction et la phase donnee par
!  l'ascendance au prochain marqueur informatif, il faudra inverser le resultat
!  e partir du marqueur (indmin) le moins clair (rapport d'effectifs en phase
!  et en opposition ratmin) le plus proche de 0.5
!
   40    continue

         ratmin=0.d0
	 indmin=ll1
	 ll2=ll1
   50	 ll2=ll2+1
         if ( ll2 <= dataset%map%nmk(c) ) then
	     if(hetero(ll2)) go to 60
         end if
!
!  cas du dernier marqueur heterozygote
!
	 if(ll2.ge.dataset%map%nmk(c)) then

	   do ll=ll1+1,dataset%map%nmk(c)
	     php(ll)=php(ll1)
	   end do
	   go to 100
	 end if

	 go to 50
   60   continue
!
!  comptage des inversions apparentes des phases locales
!
           phase=0
	   oppos=0
	   do kkd=1,dgenea%ndm(nm2+1)-dgenea%ndm(nm1)
             if(trans(kkd,ll1).ne.0.and.trans(kkd,ll2).ne.0) then
	       if(trans(kkd,ll1).eq.trans(kkd,ll2))then
		 phase=phase+1
	       else
		 oppos=oppos+1
	       end if
	     end if
           end do
	   if((phase+oppos).eq.0)go to 50
!
!  etablissement des phases locales
!
           if(phase.gt.oppos)then
	     do ll=ll1+1,ll2
	       php(ll)=php(ll1)
	     end do
	     ratio=dble(oppos)/dble(phase)
	   else
	     do ll=ll1+1,ll2
	       php(ll)=3-php(ll1)
	     end do
	     ratio=dble(phase)/dble(oppos)
	   end if

           if(ratio.ge.ratmin) then
	      ratmin=ratio
	      indmin=ll2
	   end if

           ll1=ll2
	   if(ll2.eq.dataset%map%nmk(c)) go to 100
	   if(.not.mafor(ll2)) go to 50
!
!  test de la coherence entre phase reconstruite d'apres la descendance
!  et information de l'ascendance
!  en cas d'incoherence on inverse e partir de indmin
!
	 if(php(ll2).ne.hcon(ll2)) then
	   do ll = indmin,ll2
	     php(ll)=3-php(ll)
	   end do
	   go to 40
	 end if
	 go to 50
!
  100  continue

!
!  stockage du genotype
!
        do ll=1,dataset%map%nmk(c)
          spt%genotyp(c,ll,dga%correp(ip),1)=dga%pheno(c,ll,dga%correp(ip),php(ll))
          spt%genotyp(c,ll,dga%correp(ip),2)=dga%pheno(c,ll,dga%correp(ip),3-php(ll))
        end do

!
! Reorganisation du tableau des probabilites de transmission
          do jm=nm1,nm2
          nd1=dgenea%ndm(jm)+1
          nd2=dgenea%ndm(jm+1)
          do kd=nd1,nd2
              do ll=1,dataset%map%nmk(c)
                if(php(ll).eq.2) then
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
  300 continue
     end do
      call log_mess('end PDEGP_FAST',DEBUG_DEF)
      return
      end subroutine pdegp_snp
!!***

!!****f* m_qtlmap_haplotype_V2/pdegm_snp
!! NAME
!!    pdegm_snp
!! DESCRIPTION
!!
!! NOTES
!!   Determination des probabilites des phases pour les meres
!! SOURCE
      subroutine pdegm_snp(dataset,spt)
      type(QTLMAP_DATASET)       ,intent(in)            :: dataset
       type(PDD_BUILD)            ,intent(inout)         :: spt
!
! Divers
      integer          :: geno,igeno
      integer(kind=KIND_PHENO) :: m1,m2
      integer i,ii,icon,ifem,indmin,jm,kd,kkd,c
      integer ll,ll1,ll2,lmark,nbcon,nd1,nd2
      integer ngeno1,ngeno2,nmark
      integer nombfem,ndf(dataset%genea%nm+1),repdes(dataset%genea%nm,dataset%genea%nd),corref(dataset%genea%nm)
      integer oppos,phase,phm(maxval(dataset%map%nmk))
      integer trans(dataset%genea%nd,maxval(dataset%map%nmk)),hcon(maxval(dataset%map%nmk)),indcon(maxval(dataset%map%nmk))
      logical  :: hetero(maxval(dataset%map%nmk)), mafor(maxval(dataset%map%nmk))
      real (kind=dp) :: p1,p2,ratio,ratmin
      type(GENEALOGY_BASE) , pointer :: dgenea
      type(GENOTYPE_BASE) , pointer :: dga

      dga => dataset%genoAnimal
      dgenea => dataset%genea
      call log_mess('calling...PDEGM_FAST',DEBUG_DEF)
!
!******************************************************************************
!
!
!  Creation des tableaux
!    ndf (nombre total de descendants/femelle)
!    repdes (identification des descendants de chaque femelles selon son nuero d'ordre
!  en lecture)
!    corref (l'equivalent de correm pour les femelles)
!
!  cette partie sera plutot e positionner dans lect_genea et lect_typ
!

      !init OFI :
      trans = 0

      do ifem=1,dgenea%nfem
        ndf(ifem)=0
      end do

      do jm=1,dgenea%nm
        nd1=dgenea%ndm(jm)+1
        nd2=dgenea%ndm(jm+1)
        ifem=dgenea%repfem(jm)
        do kd=nd1,nd2
          ndf(ifem)=ndf(ifem)+1
          repdes(ifem,ndf(ifem))=kd
        end do

        corref(ifem)=dga%correm(jm)
      end do

!
!  le tableau force indique si une femelle n'a pas assez de descendants
!  pour etre estimee
!
    !  do jm=1,nm
    !    force(jm)=.false.
    !  end do
!
      allocate ( spt%ngend(dataset%map%nchr,dgenea%nm+1) )
      allocate ( spt%probg(dataset%map%nchr,dgenea%nm) )
      allocate ( spt%ndesc(dataset%map%nchr,dgenea%nd) )

      spt%ngenom(:,1)=0
      spt%ngend=0
      spt%probg=0.d0
!
!******************************************************************************
!
!  on traite les femelles (et non les meres) une a une
!
!  nombfem est le nombre de femelles dont on peut chercher les phases
!
    do c=1,dataset%map%nchr
      nombfem=0
      do 500 jm=1,dgenea%nm
        nd1=dgenea%ndm(jm)+1
        nd2=dgenea%ndm(jm+1)
        ifem=dgenea%repfem(jm)
!
!******************************************************************************
!
! Si le nombre de pleins freres est trop faible, le genotype de la mere ne sera pas considere
! Dans ce cas le premier genotype rencontre possible est considere comme le bon

        !if(ndf(ifem).lt.ndmin) then
         ! force(jm)=.true.
         if(dga%estfem(ifem))then
!  on cherche le premier genotype possible
!
          do nmark=1,dataset%map%nmk(c)
            phm(nmark)=1
            if(ordref(c,ifem,nmark).eq.21)phm(nmark)=2
          end do

          go to 400

        end if
!
!******************************************************************************
!
!
!  on cree le tableau hetero qui liste l'heterozygotie des marqueurs
!
        do nmark=1,dataset%map%nmk(c)
          hetero(nmark)=.true.
          if(corref(ifem).eq.9999) then
            m1=dga%nmanque
          else
            m1=dga%pheno(c,nmark,corref(ifem),1)
            m2=dga%pheno(c,nmark,corref(ifem),2)
          end if
          if(m1.eq.dga%nmanque) m2=m1
          if(m1.eq.m2) then
            hetero(nmark)=.false.
          end if
        end do
!
!  creation du tableau des transmissions observees
!
        kkd=0
        do kd=nd1,nd2
	  kkd=kkd+1
          do ll=1,dataset%map%nmk(c)
            trans(kkd,ll)=0
            p1=spt%prot(c,ll,kd,1)+spt%prot(c,ll,kd,3)
            if(p1.eq.1.d0)trans(kkd,ll)=1
            p2=spt%prot(c,ll,kd,2)+spt%prot(c,ll,kd,4)
            if(p2.eq.1.d0)trans(kkd,ll)=2
          end do
        end do
!
!  on commence par lister les marqueurs dont l'origine est connue
!  soit par l'ascendance, soit parcequ'ils sont homozygotes
!  (dans ce cas on donne l'ordrem=12)
!

        nbcon=0
        do lmark=1,dataset%map%nmk(c)

          mafor(lmark)=.false.
	  hcon(lmark)=0
	  phm(lmark)=0

          icon=0
          if(ordref(c,ifem,lmark).eq.12)then
            nbcon=nbcon+1
            hcon(lmark)=1
            indcon(nbcon)=lmark
            mafor(lmark)=.true.
          end if

          if(ordref(c,ifem,lmark).eq.21)then
            nbcon=nbcon+1
            hcon(lmark)=2
            indcon(nbcon)=lmark
            mafor(lmark)=.true.
          end if

        end do
!
!  reconstruction de la phase la plus probable par approximations successives
!  en partant du premier marqueur heterozygote chez le pere
!
      do ll1=1,dataset%map%nmk(c)
       if(hetero(ll1)) go to 10
      end do
!
!  le pere est entierement homozygote
!
      do lmark=1,dataset%map%nmk(c)
        phm(lmark)=1
      end do

      go to 100

  10  phm(ll1)=1
      if(nbcon.ne.0)phm(ll1)=hcon(indcon(1))
!
!  si le premier marqueur informatif n'est pas le premier, on pose que
!  les  phases locales e sa gauche sont identiques
!
	 if(ll1.gt.1) then
	   do ll=ll1-1,1,-1
	       phm(ll)=phm(ll1)
	   end do
         end if
!
!  on part du premier marqueur informatif (s'il existe) et on construit
!  les phases locales en glissant vers la droite
!  ll1 et ll2 etant les indices de 2 marqueurs heterozygotes successifs
!  (tous les marqueurs intermediaires sont homozygotes), on compare
!  les transmission en ll1 et ll2, et on en deduit la phase locale
!
!  en cas d'incoherence entre cette reconstruction et la phase donnee par
!  l'ascendance au prochain marqueur informatif, il faudra inverser le resultat
!  e partir du marqueur (indmin) le moins clair (rapport d'effectifs en phase
!  et en opposition ratmin) le plus proche de 0.5
!
   40    continue

         ratmin=0.d0
	 indmin=ll1
	 ll2=ll1
   50	 ll2=ll2+1

         if ( ll2 <= dataset%map%nmk(c) ) then
	     if(hetero(ll2)) go to 60
         end if
!
!  cas du dernier marqueur heterozygote
!
	 if(ll2.ge.dataset%map%nmk(c)) then

	   do ll=ll1+1,dataset%map%nmk(c)
	     phm(ll)=phm(ll1)
	   end do
	   go to 100
	 end if

	 go to 50
   60   continue
!
!  comptage des inversions apparentes des phases locales
!
           phase=0
	   oppos=0
	   do kkd=1,ndf(ifem)
             if(trans(kkd,ll1).ne.0.and.trans(kkd,ll2).ne.0) then
	       if(trans(kkd,ll1).eq.trans(kkd,ll2))then
		 phase=phase+1
	       else
		 oppos=oppos+1
	       end if
	     end if
           end do
	   if((phase+oppos).eq.0)go to 50
!
!  etablissement des phases locales
!
           if(phase.gt.oppos)then
	     do ll=ll1+1,ll2
	       phm(ll)=phm(ll1)
	     end do
	     ratio=dble(oppos)/dble(phase)
	   else
	     do ll=ll1+1,ll2
	       phm(ll)=3-phm(ll1)
	     end do
	     ratio=dble(phase)/dble(oppos)
	   end if

           if(ratio.ge.ratmin) then
	      ratmin=ratio
	      indmin=ll2
	   end if


           ll1=ll2
	   if(ll2.eq.dataset%map%nmk(c)) go to 100
	   if(.not.mafor(ll2)) go to 50
!
!  test de la coherence entre phase reconstruite d'apres la descendance
!  et information de l'ascendance
!  en cas d'incoherence on inverse e partir de indmin
!
	 if(phm(ll2).ne.hcon(ll2)) then
	   do ll = indmin,ll2
	     phm(ll)=3-phm(ll)
	   end do
	   go to 40
	 end if
	 go to 50
!
  100  continue

  400 continue

!
!   stockage des resultats par mere (et non par femelle)
!    et reperage de la phase la plus probable
!
        spt%ngenom(c,jm+1)=spt%ngenom(c,jm)+1
!
!
        spt%probg(c,spt%ngenom(c,jm+1))=1.d0

        do ii=1,dataset%map%nmk(c)
          if (dga%correm(jm) <= size(spt%genotypm,3).and. corref(ifem) <= size(dga%pheno,3)) then
            spt%genotypm(c,ii,spt%ngenom(c,jm+1),1)=dga%pheno(c,ii,corref(ifem),phm(ii))
            spt%genotypm(c,ii,spt%ngenom(c,jm+1),2)=dga%pheno(c,ii,corref(ifem),3-phm(ii))

            if ( dga%correm(jm) > size(spt%genotypm,3) ) then
                call log_mess('dam:'//trim(dgenea%mere(jm))//' have no genotype.',WARNING_DEF)
             end if
             if ( corref(ifem) > size(dga%pheno,3) ) then
                call log_mess('dam:'//trim(dgenea%mere(jm))//' have no genotype.',WARNING_DEF)

            exit
             end if

           end if
          end do


        ngeno1=spt%ngenom(c,jm)+1
        ngeno2=spt%ngenom(c,jm+1)
        igeno=0
          do geno=ngeno1,ngeno2
            igeno=igeno+1
            spt%ngend(c,geno+1)=spt%ngend(c,geno)
            do kd=nd1,nd2
              spt%ngend(c,geno+1)=spt%ngend(c,geno+1)+1
              spt%ndesc(c,spt%ngend(c,geno+1))=kd
              do ll=1,dataset%map%nmk(c)
                do i=1,4
                  ptfin(c,ll,spt%ngend(c,geno+1),i)=spt%prot(c,ll,kd,i)
                end do
                if(phm(ll).eq.2) then
                  ptfin(c,ll,spt%ngend(c,geno+1),1)=spt%prot(c,ll,kd,2)
                  ptfin(c,ll,spt%ngend(c,geno+1),2)=spt%prot(c,ll,kd,1)
                  ptfin(c,ll,spt%ngend(c,geno+1),3)=spt%prot(c,ll,kd,4)
                  ptfin(c,ll,spt%ngend(c,geno+1),4)=spt%prot(c,ll,kd,3)
                end if

              end do
            end do
          end do

  500 continue
    end do
      call log_mess('end PDEGM_FAST',DEBUG_DEF)
      return
      end subroutine pdegm_snp
!!***

!!****f* m_qtlmap_haplotype_V2/pded_v5
!! NAME
!!    pded_v5
!! DESCRIPTION
!!
!! NOTES
!!   Calcul des probabilites de transmission le long du chromosome
!! SOURCE
      subroutine pded_v5(dataset,spt)
       type(QTLMAP_DATASET)       ,intent(in)            :: dataset
       type(PDD_BUILD)            ,intent(inout)         :: spt

!
! Tableaux dimensionnes selon l, le nombre de marqueurs
      integer trans(maxval(dataset%map%nmk))
      integer indp(maxval(dataset%map%nmk)),indm(maxval(dataset%map%nmk))
      integer ingp(maxval(dataset%map%nmk)),ingm(maxval(dataset%map%nmk))
      integer indi(maxval(dataset%map%nmk)),ingi(maxval(dataset%map%nmk))
      integer hdp(maxval(dataset%map%nmk)),hdm(maxval(dataset%map%nmk))
      integer hgp(maxval(dataset%map%nmk)),hgm(maxval(dataset%map%nmk))
      integer nbdp(maxval(dataset%map%nmk)),nbdm(maxval(dataset%map%nmk))
      integer nbgp(maxval(dataset%map%nmk)),nbgm(maxval(dataset%map%nmk))
      integer transp(maxval(dataset%map%nmk)),inforp(maxval(dataset%map%nmk))

      double precision trecm(maxval(dataset%map%nmk),maxval(dataset%map%nmk))
      double precision trecf(maxval(dataset%map%nmk),maxval(dataset%map%nmk))
      double precision srecm(maxval(dataset%map%nmk),maxval(dataset%map%nmk))
      double precision srecf(maxval(dataset%map%nmk),maxval(dataset%map%nmk))
      double precision pbt(4),pbtg(4),pbtd(4),pcum(2),pbante(4,2)

      logical infor(maxval(dataset%map%nmk))

! Divers
      integer new(11,8)
      integer geno
      integer hdi,hdifin,hdideb,hgi,hgideb,hgifin,hlocp,hlocm
      integer idm,igm,idp,igp
      integer icoloc,ii,icas,idi,igi,igeno,ilong,ip,iqtl,iqtlp,iqtlm
      integer j,jdi,jgi,jj,jk,jm,jn,jnd,jma
      integer kcas,kd,kkd
      integer lecas,linf,linfp,lk,llk,llkp,lma,lma1,lma2,ln,lq,lqtl
      integer maxcas
      integer n,nbdi,nbgi
      integer nd1,nd2,ndeb,nfin,ngeno1,ngeno2,nm1,nm2,c
      integer nmkinf,npas,nposi,nmkinfp

      integer kcutopt,temps,tempsmax
      integer kcutd,lcutd,ncutd,rcutd,kcutg,lcutg,ncutg,rcutg

      double precision apm,app
      double precision ddp,ddm,dgp,dgm,dlq,dqr
      double precision dx,pbtt,pm,pp,prob
      double precision recm,recf,reclq,recqr
      double precision recdp,recdm,trecdp,trecdm,srecdp,srecdm
      double precision recgp,recgm,trecgp,trecgm,srecgp,srecgm
      double precision xint,xintp,xintm
      logical coloc,colocdp,colocdm,colocgp,colocgm,finp,finm,possible

      type(GENEALOGY_BASE) , pointer :: dgenea
      type(GENOTYPE_BASE) , pointer :: dga

      dga => dataset%genoAnimal
      dgenea => dataset%genea

!
!  tableau des transmission
!
      data new / 1,1,2,2,1,1,1,2,1,1,1,&
        	 1,2,1,2,1,2,1,1,1,2,1,&
        	 0,0,0,0,2,2,1,2,2,2,1,&
        	 0,0,0,0,2,1,2,2,1,2,2,&
        	 0,0,0,0,0,0,0,0,0,0,2,&
        	 0,0,0,0,0,0,0,0,0,0,1,&
        	 0,0,0,0,0,0,0,0,0,0,2,&
        	 0,0,0,0,0,0,0,0,0,0,2/

  do c=1,dataset%map%nchr
!
! Taille du segment explore
      ilong=dataset%map%get_ilong(c)
!
!  tableau des recombinaison
!
      do lk=1,dataset%map%nmk(c)-1
        do jk=lk+1,dataset%map%nmk(c)
         recm=xaldane(dataset%map%posim(c,jk)-dataset%map%posim(c,lk))
         recf=xaldane(dataset%map%posif(c,jk)-dataset%map%posif(c,lk))
         trecm(lk,jk)=dlog(recm)
         trecf(lk,jk)=dlog(recf)
	     srecm(lk,jk)=dlog(1.d0-recm)
	     srecf(lk,jk)=dlog(1.d0-recf)
         trecm(jk,lk)=trecm(lk,jk)
         trecf(jk,lk)=trecf(lk,jk)
	     srecm(jk,lk)=srecm(lk,jk)
	     srecf(jk,lk)=srecf(lk,jk)
	    end do
      end do
!
! Calcul pour chaque descendant de la probabilite de transmission
       do 1000 ip=1,dgenea%np
       nm1=dgenea%nmp(ip)+1
       nm2=dgenea%nmp(ip+1)
       do 1000 jm=nm1,nm2
       ngeno1=spt%ngenom(c,jm)+1
       ngeno2=spt%ngenom(c,jm+1)
       do 1000 geno=ngeno1,ngeno2
       nd1=spt%ngend(c,geno)+1
       nd2=spt%ngend(c,geno+1)
       do 1000 kd=nd1,nd2
        kkd=spt%ndesc(c,kd)
!
!  le vecteur trans(lk) indique l'information jointe sur la transmission au locus lk
!  le vecteur infor contient la liste des marqeurs doublement informatifs
!
!  l'origine grand parentale des all�les recu par le descendant k au marqueur lk
!  d�pend de trans (tableau new)
!  soit 1 pour grand p�re et 2 pour grand m�re
!
!	           cas 1	   cas 2	   cas 3	   cas 4
!     trans(lk)	pere  mere	pere  mere	pere  mere	pere  mere
!	1	 1	1
!	2	 1	2
!	3	 2	1
!	4	 2	2
!	5	 1	1	 2	2
!	6	 1	2	 2	1
!	7	 1	1	 1	2
!	8	 2	1	 2	2
!	9	 1	1	 2	1
!	10	 1	2	 2	2
!	11	 1	1	 1	2	 2	1	 2	2
!
!

       llkp=0
       do lk=1,dataset%map%nmk(c)

!
!  infor(lk) est � true si le marqueur lk est informatif pour les deux parents
!
         infor(lk)=.false.
!
! incertitude totale (donnees manquantes)
!
         trans(lk)=11
!
!  triplet h�t�rozygotes
!
         if(ptfin(c,lk,kd,2).eq.0.and.ptfin(c,lk,kd,3).eq.0)trans(lk)=5
         if(ptfin(c,lk,kd,1).eq.0.and.ptfin(c,lk,kd,4).eq.0)trans(lk)=6
!
!  incertitude sur un parent
!
         if(ptfin(c,lk,kd,3).eq.0.and.ptfin(c,lk,kd,4).eq.0)trans(lk)=7
         if(ptfin(c,lk,kd,1).eq.0.and.ptfin(c,lk,kd,2).eq.0)trans(lk)=8
         if(ptfin(c,lk,kd,2).eq.0.and.ptfin(c,lk,kd,4).eq.0)trans(lk)=9
         if(ptfin(c,lk,kd,1).eq.0.and.ptfin(c,lk,kd,3).eq.0)trans(lk)=10
!
!  cas certains
!
         if(ptfin(c,lk,kd,1).eq.1.d0) trans(lk)=1
         if(ptfin(c,lk,kd,2).eq.1.d0) trans(lk)=2
         if(ptfin(c,lk,kd,3).eq.1.d0) trans(lk)=3
         if(ptfin(c,lk,kd,4).eq.1.d0) trans(lk)=4

	 if(trans(lk).le.4) infor(lk)=.true.
!
!
!  cas des familles de demi fr�res
!
         !if(force(jm).or.opt_sib.eq.OPT_SIB_HS) then
        if(.not.dga%estfem(dgenea%repfem(jm))&
         .or. (dataset%params%opt_sib.eq.OPT_SIB_HS) )then
	   transp(lk)=0
	   if(trans(lk).le.2.or.trans(lk).eq.7)transp(lk)=1
	   if(trans(lk).eq.3.or.trans(lk).eq.4.or.trans(lk).eq.8)transp(lk)=2
           if(transp(lk).ne.0) then
             llkp=llkp+1
	     inforp(llkp)=lk
	   end if
	 end if
       end do
       nmkinfp=llkp
!
! exploration du groupe de liaison pour la recherche des transmissions possibles
!
       n=1
  100  dx=dataset%map%absi(c,n)
!
!  pour des probl�mes d'arrondis, la position test�e pour la qtl
!  peut d�passer le bout du groupe de liaison... on l'y ramene
!
       if(dx.ge.dataset%map%posi(c,dataset%map%nmk(c)))dx=dataset%map%posi(c,dataset%map%nmk(c))
!
!  cas des demi fr�res
!
      ! if(force(jm).or.opt_sib.eq.OPT_SIB_HS) then
        if(.not.dga%estfem(dgenea%repfem(jm)).or.dataset%params%opt_sib.eq.OPT_SIB_HS) then
	do ii=1,4
	  pbt(ii)=0.d0
	end do
!
!  cas des descendants sans information
!
        if(nmkinfp.eq.0) then
          do ii=1,3,2
            pbt(ii)=0.5d0
          end do
          go to 130
        end if
!
!  descendants avec information
!
        linfp=1
!
!  position du qtl avant le premier marqueur informatif
!
	if(dx.lt.dataset%map%posi(c,inforp(1)))then
	  dqr=dataset%map%posi(c,inforp(1))-dx
	  if(inforp(1).eq.1) then
	    if(dataset%map%posi(c,1).ne.0.)dqr=dqr*dataset%map%posim(c,1)/dataset%map%posi(c,1)
	  else
	    xint=dataset%map%posi(c,inforp(1))-dataset%map%posi(c,inforp(1)-1)
	    xintm=dataset%map%posim(c,inforp(1))-dataset%map%posim(c,inforp(1)-1)
	    if(xint.ne.0.d0)dqr=dqr*xintm/xint
	  end if
	  recqr=xaldane(dqr)

	  if(transp(inforp(1)).eq.1) then
	    pbt(1)=1.d0-recqr
	    pbt(3)=recqr
	  else
	    pbt(1)=recqr
	    pbt(3)=1.d0-recqr
	  end if

	   go to 130
	 end if

    5	 linfp=linfp+1
         if(linfp.gt.nmkinfp) then
!
!  qtl apr�s le dernier marqueur informatif
!
	  dlq=dx-dataset%map%posi(c,inforp(nmkinfp))
	  if(inforp(nmkinfp).eq.dataset%map%nmk(c)) then
	    xint=dataset%map%posi(c,inforp(nmkinfp))-dataset%map%posi(c,inforp(nmkinfp)-1)
	    xintm=dataset%map%posim(c,inforp(nmkinfp))-dataset%map%posim(c,inforp(nmkinfp)-1)
	    if(xint.ne.0.d0) dlq=dlq*xintm/xint
	  else
	    xint=dataset%map%posi(c,inforp(nmkinfp)+1)-dataset%map%posi(c,inforp(nmkinfp))
	    xintm=dataset%map%posim(c,inforp(nmkinfp)+1)-dataset%map%posim(c,inforp(nmkinfp))
	    if(xint.ne.0.d0) dlq=dlq*xintm/xint
	  end if
	  reclq=xaldane(dlq)
	   if(transp(inforp(nmkinfp)).eq.1) then
	     pbt(1)=1.d0-reclq
	     pbt(3)=reclq
	   else
	     pbt(1)=reclq
	     pbt(3)=1.d0-reclq
	   end if

	   go to 130

	  else

	  if(.not.(dx.ge.dataset%map%posi(c,inforp(linfp-1)) .and. dx.le.dataset%map%posi(c,inforp(linfp)))) go to 5
!
!  qtl entre deux marqueurs informatifs
!

      	  dqr=dataset%map%posi(c,inforp(linfp))-dx
     	  dlq=dx-dataset%map%posi(c,inforp(linfp-1))
	  xint=dataset%map%posi(c,inforp(linfp))-dataset%map%posi(c,inforp(linfp-1))
	  xintm=dataset%map%posim(c,inforp(linfp))-dataset%map%posim(c,inforp(linfp-1))
	  if(xint.ne.0.d0)then
	    dqr=dqr*xintm/xint
     	    recqr=xaldane(dqr)
	    dlq=dlq*xintm/xint
     	    reclq=xaldane(dlq)
	  end if

     	  if(transp(inforp(linfp-1)).eq.transp(inforp(linfp)))then
     	    if(transp(inforp(linfp)).eq.1)then
     	      pbt(1)=(1.d0-reclq)*(1.d0-recqr)
     	      pbt(3)=reclq*recqr
     	    else
     	      pbt(1)=reclq*recqr
     	      pbt(3)=(1.d0-reclq)*(1.d0-recqr)
     	    end if
     	  else
     	    if(transp(inforp(linfp)).eq.1)then
     	      pbt(1)=reclq*(1.d0-recqr)
     	      pbt(3)=(1.d0-reclq)*recqr
     	    else
     	      pbt(1)=(1.d0-reclq)*recqr
     	      pbt(3)=reclq*(1.d0-recqr)
     	    end if
     	  end if

	  go to 130

         end if
       end if
!
!  cas des m�lange plein / demi fr�res
!
!
!  cas des positions sur des marqueurs informatifs
!
!
!  si (coloc=true) le QTL est sur exactement sur le marqueur lma (icoloc =lma)
!   le cas sera trait� sp�cifiquement pour �viter de prendre Log(z�ro)
!  si de surcroit le marqueur est informatif pour les deux parents, on connait
!   directement la probabilit� de transmission
!
      coloc=.false.
      do lk=1,dataset%map%nmk(c)
        if(int(1000.d0*dx).eq.int(1000.d0*dataset%map%posi(c,lk)))then
	  coloc=.true.
	  icoloc=lk
	  if(infor(lk))then
	    do lq=1,4
	      pbt(lq)=0.d0
	    end do
	    pbt(trans(lk))=1.d0
	    go to 130
          end if
	end if
      end do
!
!  localisation de la position test�e du qtl
!
      do lk=2,dataset%map%nmk(c)
        if(dx.ge.dataset%map%posi(c,lk-1).and.dx.lt.dataset%map%posi(c,lk)) lqtl=lk-1
      end do
        if(dx.eq.dataset%map%posi(c,dataset%map%nmk(c))) lqtl=dataset%map%nmk(c)-1
!
!  Recherche des zones informatives pour la transmission
!  La zone du p�re (resp. de la m�re) doit se terminer � gauche et � droite
!  par des marqueurs informatifs (dont l'all�le transmis au descendant est
!  connu), sans la pr�sence de marqueurs incertain (trans = 5 ou 6) entre
! les bornes � gauche (resp � droite ) des deux parents
!
!  nbdp est le nombre de marqueurs connu pour le pere � droite du qtl
!  nbdm, pour la mere
!  nbdi, nombre de marqueurs incertains
!
!  hdp est le vecteur des transmision certaines pour le pere � droite
!  hdm, pour la mere
!
!  indp est le vecteur  des indices (dans la lise initiale) des marqueurs
!  connus pour le pere � droite
!  indm, pour la mere
!  indi, le tableau des indices des marqueurs incertains (colonne 1 : le pere, 2
!  : la mere)
!
       nbdi=1
       do lk=1,dataset%map%nmk(c)
         nbdp(lk)=0
	 nbdm(lk)=0
       end do
       jma=lqtl
       finp=.false.
       finm=.false.

   10  jma=jma+1
       if(jma.gt.dataset%map%nmk(c)) go to 15
       if(trans(jma).le.4) then
         nbdp(nbdi)=nbdp(nbdi)+1
         indp(nbdp(nbdi))=jma
	 hdp(nbdp(nbdi))=floor(real(trans(jma)+1)/2.)
         nbdm(nbdi)=nbdm(nbdi)+1
         indm(nbdm(nbdi))=jma
	 hdm(nbdm(nbdi))=trans(jma)-2*(hdp(nbdp(nbdi))-1)
	 finp=.true.
	 finm=.true.
       end if

       if(trans(jma).eq.5.or.trans(jma).eq.6) then
         nbdp(nbdi)=nbdp(nbdi)+1
         indp(nbdp(nbdi))=jma
         hdp(nbdp(nbdi))=0
         nbdm(nbdi)=nbdm(nbdi)+1
         indm(nbdm(nbdi))=jma
	 hdm(nbdm(nbdi))=0
	 indi(nbdi)=jma
	 nbdi=nbdi+1
	 nbdp(nbdi)=nbdp(nbdi-1)
	 nbdm(nbdi)=nbdm(nbdi-1)
	 finp=.false.
	 finm=.false.
       end if

       if(trans(jma).eq.7.or.trans(jma).eq.8) then
         nbdp(nbdi)=nbdp(nbdi)+1
         indp(nbdp(nbdi))=jma
	 hdp(nbdp(nbdi))=trans(jma)-6
	 finp=.true.
       end if


       if(trans(jma).eq.9.or.trans(jma).eq.10) then
         nbdm(nbdi)=nbdm(nbdi)+1
         indm(nbdm(nbdi))=jma
	 hdm(nbdm(nbdi))=trans(jma)-8
	 finm=.true.
       end if

       if(((.not.finp).or.(.not.finm)).and.jma.lt.dataset%map%nmk(c)) go to 10

!
!  on cherche les zones � gauche du qtl
!  les elements sont comme � droite , en remplacant d par g
!
   15  nbgi=1
       do lk=1,dataset%map%nmk(c)
         nbgp(lk)=0
	 nbgm(lk)=0
       end do
       jma=lqtl+1
       finp=.false.
       finm=.false.

   20  jma=jma-1
       if(jma.eq.0) go to 25
       if(trans(jma).le.4) then
         nbgp(nbgi)=nbgp(nbgi)+1
         ingp(nbgp(nbgi))=jma
	 hgp(nbgp(nbgi))=floor(real(trans(jma)+1)/2.)
         nbgm(nbgi)=nbgm(nbgi)+1
         ingm(nbgm(nbgi))=jma
	 hgm(nbgm(nbgi))=trans(jma)-2*(hgp(nbgp(nbgi))-1)
	 finp=.true.
	 finm=.true.
       end if

       if(trans(jma).eq.5.or.trans(jma).eq.6) then
         nbgp(nbgi)=nbgp(nbgi)+1
         ingp(nbgp(nbgi))=jma
	 hgp(nbgp(nbgi))=0
         nbgm(nbgi)=nbgm(nbgi)+1
         ingm(nbgm(nbgi))=jma
	 hgm(nbgm(nbgi))=0
	 ingi(nbgi)=jma
	 nbgi=nbgi+1
         nbgp(nbgi)=nbgp(nbgi-1)
         nbgm(nbgi)=nbgm(nbgi-1)
	 finp=.false.
	 finm=.false.
       end if

       if(trans(jma).eq.7.or.trans(jma).eq.8) then
         nbgp(nbgi)=nbgp(nbgi)+1
         ingp(nbgp(nbgi))=jma
	 hgp(nbgp(nbgi))=trans(jma)-6
	 finp=.true.
       end if


       if(trans(jma).eq.9.or.trans(jma).eq.10) then
         nbgm(nbgi)=nbgm(nbgi)+1
         ingm(nbgm(nbgi))=jma
	 hgm(nbgm(nbgi))=trans(jma)-8
	 finm=.true.
       end if

       if(((.not.finp).or.(.not.finm)).and.jma.gt.1) go to 20

!
!   colocgp et colocgm indique (true) si la position du qtl est colocalis�e
!   avec un marqueur (coloc = true) et que ce marqueur est le premier �
!   sa gauche. Dans ce cas, il n'y a pas de recombinaison possible entre
! le qtl et ce marqueur
!
!  hlocp et hlocm sont les informations sur la transmission en ce locus (� la
!  fois marqueur et qtl)
!
   25   colocgp=.false.
        colocgm=.false.
        colocdp=.false.
        colocdm=.false.
	if(coloc.and.nbgp(1).ne.0)then
	  if(int(1000.d0*dx).eq.int(1000.d0*dataset%map%posi(c,ingp(1))))colocgp=.true.
	end if
	if(coloc.and.nbgm(1).ne.0)then
	  if(int(1000.d0*dx).eq.int(1000.d0*dataset%map%posi(c,ingm(1))))colocgm=.true.
	end if
	if(coloc.and.nbdp(1).ne.0)then
	  if(int(1000.d0*dx).eq.int(1000.d0*dataset%map%posi(c,indp(1))))colocdp=.true.
	end if
	if(coloc.and.nbdm(1).ne.0)then
	  if(int(1000.d0*dx).eq.int(1000.d0*dataset%map%posi(c,indm(1))))colocdm=.true.
	end if

!
!  les taux de recombinaison avant et apres le QTL sont stock�s dans
!   (t ou s)rec(d ou g)(p ou m)
!
        idp=indp(1)
        if(nbdp(1).eq.0)idp=lqtl+1
        igp=ingp(1)
        if(nbgp(1).eq.0)igp=lqtl

	xint=dataset%map%posi(c,lqtl+1)-dataset%map%posi(c,lqtl)
        xintp=dataset%map%posim(c,lqtl+1)-dataset%map%posim(c,lqtl)
!
!  dans le cas de colocalisation, il ne faut pas calculer le taux de
!  recombinaison
!
        dgp=(dx-dataset%map%posi(c,igp))*xintp/xint
        recgp=xaldane(dgp)
	ddp=(dataset%map%posi(c,idp)-dx)*xintp/xint
	recdp=xaldane(ddp)
!
!
!  cas des m�res
!
        idm=indm(1)
        if(nbdm(1).eq.0)idm=lqtl+1
        igm=ingm(1)
        if(nbgm(1).eq.0)igm=lqtl

	xint=dataset%map%posi(c,lqtl+1)-dataset%map%posi(c,lqtl)
	xintm=dataset%map%posif(c,lqtl+1)-dataset%map%posif(c,lqtl)
!
!  dans le cas de colocalisation, il ne faut pas calculer le taux de
!  recombinaison
!
         dgm=(dx-dataset%map%posi(c,igm))*xintm/xint
	 recgm=xaldane(dgm)
	 ddm=(dataset%map%posi(c,idm)-dx)*xintm/xint
         recdm=xaldane(ddm)

!
!  calcul des probabilit�s de transmission pour le descendant kd, position dx
!
!
!  on va envisager toutes les possibilit�s pour les marqueurs incertains
!  il y a 2**(nbdi+nbgi) possibilit�s
!  ainsi que les 4 possibilit�s pour la transmission du QTL
!
!  les hdp, hdm, hgp et hgm sont compl�t�s par des 1 ou 2 pour les marqueurs
!  incertains
!
!  les boucles en jdi et jgi ont pour fonction de g�n�rer les 2**(nbdi+nbgi)
!  vecteurs d'origines correspondants aux marqueurs incertains
!
!  initialisation
      pbtt=0.d0
      do lq=1,4
	pbtg(lq)=0.d0
	pbtd(lq)=0.d0
      end do

!
!  on va consid�rer les 4 �v�nements de transmission possible au QTL, d�finis
!  pas iqtlp et iqtlm
!
!  pour chacun des �venements de transmision au qtl on v�rifie la coh�rence
!  avec un �ventuel marqueur colocalis�
!
      iqtlp=0
   30 iqtlp=iqtlp+1
      if(iqtlp.eq.3)go to 120

      iqtlm=0
   40 iqtlm=iqtlm+1
      if(iqtlm.eq.3)go to 30
!
      iqtl= iqtlm+2*(iqtlp-1)

!
!  le cas o� il n'y a pas d'incertitude � droite est trait� � part
!  dans ce cas, il faut distinguer les situations o� il n'y a que des
!  transmissions inconnues
!
      if(nbdi.eq.1) then
!
!  on a nbdp=0 si aucun marqueur � droite n'a trans diff�rent de 9, 10 ou 11
!  on a nbdm=0 si aucun marqueur � droite n'a trans diff�rent de 7,  8 ou 11
!  dans ce cas,les 2 �v�nements iqtlp = 1 ou 2 sont �quiprobables
!
         if(nbdp(1).eq.0.and.nbdm(1).eq.0) then
	   pbtd(iqtl)=0.25d0
         end if

                               prob=1.d0
         if(nbdm(1).ne.0) then
	   if(colocdm.and.iqtlm.ne.hdm(1))go to 50
	   if(.not.colocdm)then
	     if(nbdp(1).eq.0)   prob=0.5d0
             if(iqtlm.eq.hdm(1))prob=prob*(1.d0-recdm)
             if(iqtlm.ne.hdm(1))prob=prob*recdm
	   end if
	 end if

          if(nbdp(1).ne.0) then
	   if(colocdp.and.iqtlp.ne.hdp(1))go to 50
	   if(.not.colocdp)then
	     if(nbdm(1).eq.0)   prob=0.5d0
             if(iqtlp.eq.hdp(1))prob=prob*(1.d0-recdp)
             if(iqtlp.ne.hdp(1))prob=prob*recdp
	   end if
	 end if

         pbtd(iqtl)=prob
	 go to 50

       end if
!
!  cas ou il y a des incertitudes � droite
!
!
!
!  on traite diff�remment les intervalles entre marqueurs incertains : le 1er, les
!  suivants, le dernier
!
!  cas du premier intervalle
!  nfin est l'indice du marqueur incertain � droite
!
         idi=1
         possible=.false.
         do hdifin=1,2
!
!  on a deux possibilit�s couvrant l'incertitude, donnant les derniers
!  hdp et hdm du premier intervalle
!
	   nfin=indi(idi)
	   hdp(nbdp(1))=new(trans(nfin),2*hdifin-1)
	   hdm(nbdm(1))=new(trans(nfin),2*hdifin)

	   if((colocdm.and.iqtlm.ne.hdm(1)).or.(colocdp.and.iqtlp.ne.hdp(1)))then

              pbante(iqtl,hdifin)=0.d0

	   else

	     possible=.true.
             prob=0.d0

             if(.not.colocdm)then
               if(iqtlm.eq.hdm(1))prob=dlog(1.d0-recdm)
               if(iqtlm.ne.hdm(1))prob=dlog(recdm )

               do jma=2,nbdm(1)
                 if(hdm(jma-1).eq.hdm(jma))then
	           prob=prob+srecf(indm(jma-1),indm(jma))
	         else
                   prob=prob+trecf(indm(jma-1),indm(jma))
	         end if
               end do
	     end if

             if(.not.colocdp)then
               if(iqtlp.eq.hdp(1))prob=prob+dlog(1.d0-recdp)
               if(iqtlp.ne.hdp(1))prob=prob+dlog(recdp)
               do jma=2,nbdp(1)
                 if(hdp(jma-1).eq.hdp(jma))then
	           prob=prob+srecm(indp(jma-1),indp(jma))
	         else
                   prob=prob+trecm(indp(jma-1),indp(jma))
	         end if
               end do
	     end if

             pbante(iqtl,hdifin)=dexp(prob)

           end if

         end do

	 if(.not.possible) go to 50
!
!  cas des intervalles suivants
!
         do idi=2,nbdi-1

           do hdifin=1,2
!
!  on a deux possibilit�s couvrant l'incertitude, donnant les derniers
!  hdp et hdm du idi �me intervalle
!
	    nfin=indi(idi)
	    hdp(nbdp(idi))=new(trans(nfin),2*hdifin-1)
	    hdm(nbdm(idi))=new(trans(nfin),2*hdifin)
            pcum(hdifin)=0.d0

             do hdideb=1,2
!
!  il faut consid�rer les  deux possibilit�s couvrant l'incertitude,
!  de l'intervalle pr�c�dent
!
	      ndeb=indi(idi-1)
	      hdp(nbdp(idi-1))=new(trans(ndeb),2*hdideb-1)
	      hdm(nbdm(idi-1))=new(trans(ndeb),2*hdideb)

              prob=0.d0
	      do jma=nbdp(idi-1)+1,nbdp(idi)
                if(hdp(jma-1).eq.hdp(jma))then
	          prob=prob+srecm(indp(jma-1),indp(jma))
	        else
                  prob=prob+trecm(indp(jma-1),indp(jma))
	        end if
              end do

	      do jma=nbdm(idi-1)+1,nbdm(idi)
                if(hdm(jma-1).eq.hdm(jma))then
	          prob=prob+srecf(indm(jma-1),indm(jma))
	        else
                  prob=prob+trecf(indm(jma-1),indm(jma))
	        end if
              end do
!
	      pcum(hdifin)=pcum(hdifin)+pbante(iqtl,hdideb)*dexp(prob)

             end do
	   end do

	   do hdifin=1,2
	     pbante(iqtl,hdifin)=pcum(hdifin)
	   end do

	 end do

!
!  cas du dernier intervalle
!
         idi=nbdi
	 pcum(1)=0.d0
         do hdideb=1,2
!
!  il faut consid�rer les  deux possibilit�s couvrant l'incertitude,
!  de l'intervalle pr�c�dent
!
           ndeb=indi(idi-1)
           hdp(nbdp(idi-1))=new(trans(ndeb),2*hdideb-1)
           hdm(nbdm(idi-1))=new(trans(ndeb),2*hdideb)

           prob=0.d0
           if(nbdp(idi).ne.nbdp(idi-1))then
             do jma=nbdp(idi-1)+1,nbdp(idi)
               if(hdp(jma-1).eq.hdp(jma))then
        	 prob=prob+srecm(indp(jma-1),indp(jma))
               else
        	 prob=prob+trecm(indp(jma-1),indp(jma))
               end if
             end do
           end if

           if(nbdm(idi).ne.nbdm(idi-1))then
             do jma=nbdm(idi-1)+1,nbdm(idi)
               if(hdm(jma-1).eq.hdm(jma))then
        	 prob=prob+srecf(indm(jma-1),indm(jma))
               else
        	 prob=prob+trecf(indm(jma-1),indm(jma))
               end if
             end do
           end if
!
	   pcum(1)=pcum(1)+pbante(iqtl,hdideb)*dexp(prob)

         end do

         pbtd(iqtl)=pcum(1)

   50  continue
!
!  on �num�re les possibilit�s pour les transmissions incertaines � gauche
!  du qtl.
!
!
!  le cas o� il n'y a pas d'incertitude � gauche est trait� � part
!  dans ce cas, il faut distinguer les situations o� il n'y a que des
!  transmissions inconnues
!

      if(nbgi.eq.1) then
!
!  on a nbgp=0 si aucun marqueur � gauche n'a trans diff�rent de 9, 10 ou 11
!  on a nbgm=0 si aucun marqueur � gauche n'a trans diff�rent de 7,  8 ou 11
!  dans ce cas,les 2 �v�nements iqtlp = 1 ou 2 sont �quiprobables
!
         if(nbgp(1).eq.0.and.nbgm(1).eq.0) then
	   pbtg(iqtl)=0.25d0
         end if

                                prob=1.d0
         if(nbgm(1).ne.0) then
	   if(colocgm.and.iqtlm.ne.hgm(1))go to 60
	   if(.not.colocgm)then
	     if(nbgp(1).eq.0)   prob=0.5d0
             if(iqtlm.eq.hgm(1))prob=prob*(1.d0-recgm)
             if(iqtlm.ne.hgm(1))prob=prob*recgm
	   end if
	 end if

         if(nbgp(1).ne.0) then
	   if(colocgp.and.iqtlp.ne.hgp(1))go to 60
	   if(.not.colocgp)then
	     if(nbgm(1).eq.0)   prob=0.5d0
             if(iqtlp.eq.hgp(1))prob=prob*(1.d0-recgp)
             if(iqtlp.ne.hgp(1))prob=prob*recgp
	   end if
	 end if

         pbtg(iqtl)=prob
 	 go to 60

       end if
!
!  cas ou il y a des incertitudes � gauche
!
!
!
!  on traite diff�remment les intervalles entre marqueurs incertains : le 1er, les
!  suivants, le dernier
!
!  cas du premier intervalle
!  nfin est l'indice du marqueur incertain � gauche
!
         igi=1
         possible=.false.
         do hgifin=1,2

!
!  on a deux possibilit�s couvrant l'incertitude, donnant les derniers
!  hgp et hgm du premier intervalle
!
	   nfin=ingi(igi)
	   hgp(nbgp(1))=new(trans(nfin),2*hgifin-1)
	   hgm(nbgm(1))=new(trans(nfin),2*hgifin)

	   if(  (colocgm.and.iqtlm.ne.hgm(1)).or.(colocgp.and.iqtlp.ne.hgp(1)))then

              pbante(iqtl,hgifin)=0.d0

	   else

	    possible=.true.
             prob=0.d0

             if(.not.colocgm)then
               if(iqtlm.eq.hgm(1))prob=dlog(1.d0-recgm)
               if(iqtlm.ne.hgm(1))prob=dlog(recgm)

               do jma=2,nbgm(1)
         	 if(hgm(jma-1).eq.hgm(jma))then
         	   prob=prob+srecf(ingm(jma-1),ingm(jma))
         	 else
         	   prob=prob+trecf(ingm(jma-1),ingm(jma))
         	 end if
               end do
             end if

             if(.not.colocgp)then
               if(iqtlp.eq.hgp(1))prob=prob+dlog(1.d0-recgp)
               if(iqtlp.ne.hgp(1))prob=prob+dlog(recgp)
               do jma=2,nbgp(1)
         	 if(hgp(jma-1).eq.hgp(jma))then
         	   prob=prob+srecm(ingp(jma-1),ingp(jma))
         	 else
         	   prob=prob+trecm(ingp(jma-1),ingp(jma))
         	 end if
               end do
             end if

             pbante(iqtl,hgifin)=dexp(prob)

	   end if

         end do

	 if(.not.possible)go to 60
!
!  cas des intervalles suivants
!
         do igi=2,nbgi-1

           do hgifin=1,2
!
!  on a deux possibilit�s couvrant l'incertitude, donnant les derniers
!  hgp et hgm du igi �me intervalle
!
	    nfin=ingi(igi)
	    hgp(nbgp(igi))=new(trans(nfin),2*hgifin-1)
	    hgm(nbgm(igi))=new(trans(nfin),2*hgifin)
            pcum(hgifin)=0.d0

             do hgideb=1,2
!
!  il faut consid�rer les  deux possibilit�s couvrant l'incertitude,
!  de l'intervalle pr�c�dent
!
	      ndeb=ingi(igi-1)
	      hgp(nbgp(igi-1))=new(trans(ndeb),2*hgideb-1)
	      hgm(nbgm(igi-1))=new(trans(ndeb),2*hgideb)

              prob=0.d0
	      do jma=nbgp(igi-1)+1,nbgp(igi)
                if(hgp(jma-1).eq.hgp(jma))then
	          prob=prob+srecm(ingp(jma-1),ingp(jma))
	        else
                  prob=prob+trecm(ingp(jma-1),ingp(jma))
	        end if
              end do

	      do jma=nbgm(igi-1)+1,nbgm(igi)
                if(hgm(jma-1).eq.hgm(jma))then
	          prob=prob+srecf(ingm(jma-1),ingm(jma))
	        else
                  prob=prob+trecf(ingm(jma-1),ingm(jma))
	        end if
              end do
!
	      pcum(hgifin)=pcum(hgifin)+pbante(iqtl,hgideb)*dexp(prob)

             end do
	   end do

	   do hgifin=1,2
	     pbante(iqtl,hgifin)=pcum(hgifin)
	   end do

	 end do

!
!  cas du dernier intervalle
!
         igi=nbgi
	 pcum(1)=0.d0
         do hgideb=1,2
!
!  il faut consid�rer les  deux possibilit�s couvrant l'incertitude,
!  de l'intervalle pr�c�dent
!
           ndeb=ingi(igi-1)
           hgp(nbgp(igi-1))=new(trans(ndeb),2*hgideb-1)
           hgm(nbgm(igi-1))=new(trans(ndeb),2*hgideb)

           prob=0.d0
           if(nbgp(igi).ne.nbgp(igi-1))then
             do jma=nbgp(igi-1)+1,nbgp(igi)
               if(hgp(jma-1).eq.hgp(jma))then
        	 prob=prob+srecm(ingp(jma-1),ingp(jma))
               else
        	 prob=prob+trecm(ingp(jma-1),ingp(jma))
               end if
             end do
           end if

           if(nbgm(igi).ne.nbgm(igi-1))then
             do jma=nbgm(igi-1)+1,nbgm(igi)
               if(hgm(jma-1).eq.hgm(jma))then
        	 prob=prob+srecf(ingm(jma-1),ingm(jma))
               else
        	 prob=prob+trecf(ingm(jma-1),ingm(jma))
               end if
             end do
           end if
!
	   pcum(1)=pcum(1)+pbante(iqtl,hgideb)*dexp(prob)

         end do

         pbtg(iqtl)= pcum(1)

   60 continue
!
!  on teste si d'autres combinaisons au qtl sont � �valuer
!
       if(iqtlp.eq.2.and.iqtlm.eq.2) go to 120
       if(iqtlm.eq.1)go to 40
       if(iqtlm.eq.2)go to 30


  120  continue

!
!  on stocke la probabilit� de transmission � la position dx
!
         do lq=1,4
	   pbt(lq)=pbtg(lq)*pbtd(lq)
	 end do

  130   continue

	 pbtt=0.d0
         do lq=1,4
	   pbtt=pbtt+pbt(lq)
	 end do

         do lq=1,4
	   spt%pdd(c,kd,lq,n)=pbt(lq)/pbtt
	 end do

!
!
!  on avance d'un pas (en verifiant qu'on ne depasse pas le bout droit) et on recommence
!
         n=n+1
     if ( n > dataset%map%get_npo(c) ) go to 1000
         dx=dataset%map%absi(c,n)
	 if( dx.gt.dataset%map%posi(c,dataset%map%nmk(c)) ) go to 1000
	 go to 100
 1000 continue

   end do
       call log_mess('END PDED V5....',DEBUG_DEF)
      return
   end subroutine pded_v5
!!***

!!****f* m_qtlmap_haplotype_V2/pded_v5_optim
!! NAME
!!    pded_v5_optim
!! DESCRIPTION
!!
!! NOTES
!!   Calcul des probabilites de transmission le long du chromosome
!! SOURCE
    subroutine pded_v5_optim(dataset,spt)
      type(QTLMAP_DATASET)       ,intent(in)            :: dataset
      type(PDD_BUILD)            ,intent(inout)         :: spt
        
! Divers
      integer          :: c,ilong,lk,jk,ip,jm,geno,kd,iITer,nIter,sizeT,io
      real(kind=dp)    :: recm,recf
      integer ,dimension(:),allocatable  :: ipT,jmT,genoT,kdT
      real(kind=dp) , dimension(:,:), allocatable::  trecm,srecm,trecf,srecf
      type(GENEALOGY_BASE) , pointer :: dgenea

      dgenea => dataset%genea

      allocate (trecm(maxval(dataset%map%nmk),maxval(dataset%map%nmk)),&
      trecf(maxval(dataset%map%nmk),maxval(dataset%map%nmk)),stat=io)
      if ( io /= 0 ) call stop_application("calcul prob of transmission : not enough memory...")      
      allocate (srecm(maxval(dataset%map%nmk),maxval(dataset%map%nmk)),&
      srecf(maxval(dataset%map%nmk),maxval(dataset%map%nmk)),stat=io)
      if ( io /= 0 ) call stop_application("calcul prob of transmission : not enough memory...")      

      sizeT = maxval(spt%ngend(:,(maxval(spt%ngenom(:,dgenea%nm+1)))+1))
      allocate ( ipT( sizeT ),stat=io )
      if ( io /= 0 ) call stop_application("calcul prob of transmission : not enough memory...")      
      allocate ( jmT( sizeT ),stat=io )
      if ( io /= 0 ) call stop_application("calcul prob of transmission : not enough memory...")      
      allocate ( genoT( sizeT ),stat=io )
      if ( io /= 0 ) call stop_application("calcul prob of transmission : not enough memory...")      
      allocate ( kdT( sizeT ) ,stat=io)
      if ( io /= 0 ) call stop_application("calcul prob of transmission : not enough memory...")

     do c=1,dataset%map%nchr
!
! Taille du segment explore
      ilong=dataset%map%get_ilong(c)
!
!  tableau des recombinaison
!
      do lk=1,dataset%map%nmk(c)-1
        do jk=lk+1,dataset%map%nmk(c)
         recm=xaldane(dataset%map%posim(c,jk)-dataset%map%posim(c,lk))
         recf=xaldane(dataset%map%posif(c,jk)-dataset%map%posif(c,lk))
         trecm(lk,jk)=dlog(recm)
         trecf(lk,jk)=dlog(recf)
	     srecm(lk,jk)=dlog(1.d0-recm)
	     srecf(lk,jk)=dlog(1.d0-recf)
         trecm(jk,lk)=trecm(lk,jk)
         trecf(jk,lk)=trecf(lk,jk)
	     srecm(jk,lk)=srecm(lk,jk)
	     srecf(jk,lk)=srecf(lk,jk)
	    end do
      end do
!
! Calcul pour chaque descendant de la probabilite de transmission

     !calcul du nopmbre nombre d iteration
      nIter = 0
       do ip=1,dgenea%np
        do jm=dgenea%nmp(ip)+1,dgenea%nmp(ip+1)
          do geno=spt%ngenom(c,jm)+1,spt%ngenom(c,jm+1)
           do kd=spt%ngend(c,geno)+1,spt%ngend(c,geno+1)
            nIter = nIter + 1
            ipT(nIter) = ip
            jmT(nIter) = jm
            genoT(nIter) = geno
            kdT(nIter) = kd
           end do
          end do
        end do
       end do

       !Calcul + parallelisme
       !$OMP PARALLEL DO DEFAULT(SHARED)
       do iITer = 1, nIter
         call pded_v5_kd(dataset,spt,c,ipT(iIter),jmT(iIter),genoT(iIter),kdT(iIter),&
                         ilong,recm,recf,trecm,trecf,srecm,srecf)
       end do
       !$OMP END PARALLEL DO

      end do !c

      deallocate ( ipT )
      deallocate ( jmT )
      deallocate ( genoT )
      deallocate ( kdT )
      deallocate (trecm,srecm,trecf,srecf)

     end  subroutine pded_v5_optim
!!***

!!****f* m_qtlmap_haplotype_V2/pded_v5_kd
!! NAME
!!    pded_v5_kd
!! DESCRIPTION
!!
!! NOTES
!!   Calcul des probabilites de transmission le long du chromosome
!! SOURCE
      subroutine pded_v5_kd(dataset,spt,c,ip,jm,geno,kd,ilong,recm,recf,trecm,trecf,srecm,srecf)
      type(QTLMAP_DATASET)       ,intent(in)            :: dataset
      type(PDD_BUILD)            ,intent(inout)         :: spt
      integer       , intent(in) :: c,ip,jm,geno,kd,ilong
      real(kind=dp) , intent(in) :: recm,recf
      real(kind=dp) , intent(in) :: trecm(maxval(dataset%map%nmk),maxval(dataset%map%nmk))
      real(kind=dp) , intent(in) :: trecf(maxval(dataset%map%nmk),maxval(dataset%map%nmk))
      real(kind=dp) , intent(in) :: srecm(maxval(dataset%map%nmk),maxval(dataset%map%nmk))
      real(kind=dp) , intent(in) :: srecf(maxval(dataset%map%nmk),maxval(dataset%map%nmk))
!
! Tableaux dimensionnes selon l, le nombre de marqueurs
      integer trans(maxval(dataset%map%nmk))
      integer indp(maxval(dataset%map%nmk)),indm(maxval(dataset%map%nmk))
      integer ingp(maxval(dataset%map%nmk)),ingm(maxval(dataset%map%nmk))
      integer indi(maxval(dataset%map%nmk)),ingi(maxval(dataset%map%nmk))
      integer hdp(maxval(dataset%map%nmk)),hdm(maxval(dataset%map%nmk))
      integer hgp(maxval(dataset%map%nmk)),hgm(maxval(dataset%map%nmk))
      integer nbdp(maxval(dataset%map%nmk)),nbdm(maxval(dataset%map%nmk))
      integer nbgp(maxval(dataset%map%nmk)),nbgm(maxval(dataset%map%nmk))
      integer transp(maxval(dataset%map%nmk)),inforp(maxval(dataset%map%nmk))
      real(kind=dp) :: pbt(4),pbtg(4),pbtd(4),pcum(2),pbante(4,2)

      logical infor(maxval(dataset%map%nmk))

! Divers
      integer new(11,8)
      integer hdi,hdifin,hdideb,hgi,hgideb,hgifin,hlocp,hlocm
      integer idm,igm,idp,igp
      integer icoloc,ii,icas,idi,igi,igeno,iqtl,iqtlp,iqtlm
      integer j,jdi,jgi,jj,jk,jn,jnd,jma
      integer kcas,kkd
      integer lecas,linf,linfp,lk,llk,llkp,lma,lma1,lma2,ln,lq,lqtl
      integer maxcas
      integer n,nbdi,nbgi
      integer ndeb,nfin
      integer nmkinf,npas,nposi,nmkinfp

      integer kcutopt,temps,tempsmax
      integer kcutd,lcutd,ncutd,rcutd,kcutg,lcutg,ncutg,rcutg

      real(kind=dp) :: apm,app
      real(kind=dp) :: ddp,ddm,dgp,dgm,dlq,dqr
      real(kind=dp) :: dx,pbtt,pm,pp,prob
      real(kind=dp) :: reclq,recqr
      real(kind=dp) :: recdp,recdm,trecdp,trecdm,srecdp,srecdm
      real(kind=dp) :: recgp,recgm,trecgp,trecgm,srecgp,srecgm
      real(kind=dp) :: xint,xintp,xintm
      logical coloc,colocdp,colocdm,colocgp,colocgm,finp,finm,possible
      type(GENEALOGY_BASE) , pointer :: dgenea
      type(GENOTYPE_BASE) , pointer :: dga

      dga => dataset%genoAnimal
      dgenea => dataset%genea

!
!  tableau des transmission
!
      data new / 1,1,2,2,1,1,1,2,1,1,1,&
        	 1,2,1,2,1,2,1,1,1,2,1,&
        	 0,0,0,0,2,2,1,2,2,2,1,&
        	 0,0,0,0,2,1,2,2,1,2,2,&
        	 0,0,0,0,0,0,0,0,0,0,2,&
        	 0,0,0,0,0,0,0,0,0,0,1,&
        	 0,0,0,0,0,0,0,0,0,0,2,&
        	 0,0,0,0,0,0,0,0,0,0,2/


     kkd=spt%ndesc(c,kd)
    ! print *,"new kkd:",kkd
!
!  le vecteur trans(lk) indique l'information jointe sur la transmission au locus lk
!  le vecteur infor contient la liste des marqeurs doublement informatifs
!
!  l'origine grand parentale des all�les recu par le descendant k au marqueur lk
!  d�pend de trans (tableau new)
!  soit 1 pour grand p�re et 2 pour grand m�re
!
!	           cas 1	   cas 2	   cas 3	   cas 4
!     trans(lk)	pere  mere	pere  mere	pere  mere	pere  mere
!	1	 1	1
!	2	 1	2
!	3	 2	1
!	4	 2	2
!	5	 1	1	 2	2
!	6	 1	2	 2	1
!	7	 1	1	 1	2
!	8	 2	1	 2	2
!	9	 1	1	 2	1
!	10	 1	2	 2	2
!	11	 1	1	 1	2	 2	1	 2	2
!
!

       llkp=0
       do lk=1,dataset%map%nmk(c)

!
!  infor(lk) est � true si le marqueur lk est informatif pour les deux parents
!
         infor(lk)=.false.
!
! incertitude totale (donnees manquantes)
!
         trans(lk)=11
!
!  triplet h�t�rozygotes
!
         if(ptfin(c,lk,kd,2).eq.0.and.ptfin(c,lk,kd,3).eq.0)trans(lk)=5
         if(ptfin(c,lk,kd,1).eq.0.and.ptfin(c,lk,kd,4).eq.0)trans(lk)=6
!
!  incertitude sur un parent
!
         if(ptfin(c,lk,kd,3).eq.0.and.ptfin(c,lk,kd,4).eq.0)trans(lk)=7
         if(ptfin(c,lk,kd,1).eq.0.and.ptfin(c,lk,kd,2).eq.0)trans(lk)=8
         if(ptfin(c,lk,kd,2).eq.0.and.ptfin(c,lk,kd,4).eq.0)trans(lk)=9
         if(ptfin(c,lk,kd,1).eq.0.and.ptfin(c,lk,kd,3).eq.0)trans(lk)=10
!
!  cas certains
!
         if(ptfin(c,lk,kd,1).eq.1.d0) trans(lk)=1
         if(ptfin(c,lk,kd,2).eq.1.d0) trans(lk)=2
         if(ptfin(c,lk,kd,3).eq.1.d0) trans(lk)=3
         if(ptfin(c,lk,kd,4).eq.1.d0) trans(lk)=4

	 if(trans(lk).le.4) infor(lk)=.true.
!
!
!  cas des familles de demi fr�res
!
         !if(force(jm).or.opt_sib.eq.OPT_SIB_HS) then
        if(.not.(dga%estfem(dgenea%repfem(jm))).or.&
         (dataset%params%opt_sib.eq.OPT_SIB_HS))then
	   transp(lk)=0
	   if(trans(lk).le.2.or.trans(lk).eq.7)transp(lk)=1
	   if(trans(lk).eq.3.or.trans(lk).eq.4.or.trans(lk).eq.8)transp(lk)=2
           if(transp(lk).ne.0) then
             llkp=llkp+1
	     inforp(llkp)=lk
	   end if
	 end if
       end do
       nmkinfp=llkp
!
! exploration du groupe de liaison pour la recherche des transmissions possibles
!
       n=1
  100  dx=dataset%map%absi(c,n)
!
!  pour des probl�mes d'arrondis, la position test�e pour la qtl
!  peut d�passer le bout du groupe de liaison... on l'y ramene
!
       if(dx.ge.dataset%map%posi(c,dataset%map%nmk(c)))dx=dataset%map%posi(c,dataset%map%nmk(c))
!
!  cas des demi fr�res
!
      ! if(force(jm).or.opt_sib.eq.OPT_SIB_HS) then
        if(.not.dga%estfem(dgenea%repfem(jm)).or.&
          dataset%params%opt_sib.eq.OPT_SIB_HS) then
	do ii=1,4
	  pbt(ii)=0.d0
	end do
!
!  cas des descendants sans information
!
        if(nmkinfp.eq.0) then
          do ii=1,3,2
            pbt(ii)=0.5d0
          end do
          go to 130 ! condition de sortie
        end if
!
!  descendants avec information
!
        linfp=1
!
!  position du qtl avant le premier marqueur informatif
!
	if(dx.lt.dataset%map%posi(c,inforp(1)))then
	  dqr=dataset%map%posi(c,inforp(1))-dx
	  if(inforp(1).eq.1) then
	    if(dataset%map%posi(c,1).ne.0.)dqr=dqr*dataset%map%posim(c,1)/dataset%map%posi(c,1)
	  else
	    xint=dataset%map%posi(c,inforp(1))-dataset%map%posi(c,inforp(1)-1)
	    xintm=dataset%map%posim(c,inforp(1))-dataset%map%posim(c,inforp(1)-1)
	    if(xint.ne.0.d0)dqr=dqr*xintm/xint
	  end if
	  recqr=xaldane(dqr)

	  if(transp(inforp(1)).eq.1) then
	    pbt(1)=1.d0-recqr
	    pbt(3)=recqr
	  else
	    pbt(1)=recqr
	    pbt(3)=1.d0-recqr
	  end if

	   go to 130
	 end if

    5	 linfp=linfp+1
         if(linfp.gt.nmkinfp) then
!
!  qtl apr�s le dernier marqueur informatif
!
	  dlq=dx-dataset%map%posi(c,inforp(nmkinfp))
	  if(inforp(nmkinfp).eq.dataset%map%nmk(c)) then
	    xint=dataset%map%posi(c,inforp(nmkinfp))-dataset%map%posi(c,inforp(nmkinfp)-1)
	    xintm=dataset%map%posim(c,inforp(nmkinfp))-dataset%map%posim(c,inforp(nmkinfp)-1)
	    if(xint.ne.0.d0) dlq=dlq*xintm/xint
	  else
	    xint=dataset%map%posi(c,inforp(nmkinfp)+1)-dataset%map%posi(c,inforp(nmkinfp))
	    xintm=dataset%map%posim(c,inforp(nmkinfp)+1)-dataset%map%posim(c,inforp(nmkinfp))
	    if(xint.ne.0.d0) dlq=dlq*xintm/xint
	  end if
	  reclq=xaldane(dlq)
	   if(transp(inforp(nmkinfp)).eq.1) then
	     pbt(1)=1.d0-reclq
	     pbt(3)=reclq
	   else
	     pbt(1)=reclq
	     pbt(3)=1.d0-reclq
	   end if

	   go to 130

	  else

	  if(.not.(dx.ge.dataset%map%posi(c,inforp(linfp-1)) .and. dx.le.dataset%map%posi(c,inforp(linfp)))) go to 5
!
!  qtl entre deux marqueurs informatifs
!

      	  dqr=dataset%map%posi(c,inforp(linfp))-dx
     	  dlq=dx-dataset%map%posi(c,inforp(linfp-1))
	  xint=dataset%map%posi(c,inforp(linfp))-dataset%map%posi(c,inforp(linfp-1))
	  xintm=dataset%map%posim(c,inforp(linfp))-dataset%map%posim(c,inforp(linfp-1))
	  if(xint.ne.0.d0)then
	    dqr=dqr*xintm/xint
     	    recqr=xaldane(dqr)
	    dlq=dlq*xintm/xint
     	    reclq=xaldane(dlq)
	  end if

     	  if(transp(inforp(linfp-1)).eq.transp(inforp(linfp)))then
     	    if(transp(inforp(linfp)).eq.1)then
     	      pbt(1)=(1.d0-reclq)*(1.d0-recqr)
     	      pbt(3)=reclq*recqr
     	    else
     	      pbt(1)=reclq*recqr
     	      pbt(3)=(1.d0-reclq)*(1.d0-recqr)
     	    end if
     	  else
     	    if(transp(inforp(linfp)).eq.1)then
     	      pbt(1)=reclq*(1.d0-recqr)
     	      pbt(3)=(1.d0-reclq)*recqr
     	    else
     	      pbt(1)=(1.d0-reclq)*recqr
     	      pbt(3)=reclq*(1.d0-recqr)
     	    end if
     	  end if

	  go to 130

         end if
       end if
!
!  cas des m�lange plein / demi fr�res
!
!
!  cas des positions sur des marqueurs informatifs
!
!
!  si (coloc=true) le QTL est sur exactement sur le marqueur lma (icoloc =lma)
!   le cas sera trait� sp�cifiquement pour �viter de prendre Log(z�ro)
!  si de surcroit le marqueur est informatif pour les deux parents, on connait
!   directement la probabilit� de transmission
!
      coloc=.false.
      do lk=1,dataset%map%nmk(c)
        if(int(1000.d0*dx).eq.int(1000.d0*dataset%map%posi(c,lk)))then
	  coloc=.true.
	  icoloc=lk
	  if(infor(lk))then
	    do lq=1,4
	      pbt(lq)=0.d0
	    end do
	    pbt(trans(lk))=1.d0
	    go to 130
          end if
	end if
      end do
!
!  localisation de la position test�e du qtl
!
      do lk=2,dataset%map%nmk(c)
        if(dx.ge.dataset%map%posi(c,lk-1).and.dx.lt.dataset%map%posi(c,lk)) lqtl=lk-1
      end do
        if(dx.eq.dataset%map%posi(c,dataset%map%nmk(c))) lqtl=dataset%map%nmk(c)-1
!
!  Recherche des zones informatives pour la transmission
!  La zone du p�re (resp. de la m�re) doit se terminer � gauche et � droite
!  par des marqueurs informatifs (dont l'all�le transmis au descendant est
!  connu), sans la pr�sence de marqueurs incertain (trans = 5 ou 6) entre
! les bornes � gauche (resp � droite ) des deux parents
!
!  nbdp est le nombre de marqueurs connu pour le pere � droite du qtl
!  nbdm, pour la mere
!  nbdi, nombre de marqueurs incertains
!
!  hdp est le vecteur des transmision certaines pour le pere � droite
!  hdm, pour la mere
!
!  indp est le vecteur  des indices (dans la lise initiale) des marqueurs
!  connus pour le pere � droite
!  indm, pour la mere
!  indi, le tableau des indices des marqueurs incertains (colonne 1 : le pere, 2
!  : la mere)
!
       nbdi=1
       do lk=1,dataset%map%nmk(c)
         nbdp(lk)=0
	 nbdm(lk)=0
       end do
       jma=lqtl
       finp=.false.
       finm=.false.

   10  jma=jma+1
       if(jma.gt.dataset%map%nmk(c)) go to 15
       if(trans(jma).le.4) then
         nbdp(nbdi)=nbdp(nbdi)+1
         indp(nbdp(nbdi))=jma
	 hdp(nbdp(nbdi))=floor(real(trans(jma)+1)/2.)
         nbdm(nbdi)=nbdm(nbdi)+1
         indm(nbdm(nbdi))=jma
	 hdm(nbdm(nbdi))=trans(jma)-2*(hdp(nbdp(nbdi))-1)
	 finp=.true.
	 finm=.true.
       end if

       if(trans(jma).eq.5.or.trans(jma).eq.6) then
         nbdp(nbdi)=nbdp(nbdi)+1
         indp(nbdp(nbdi))=jma
         hdp(nbdp(nbdi))=0
         nbdm(nbdi)=nbdm(nbdi)+1
         indm(nbdm(nbdi))=jma
	 hdm(nbdm(nbdi))=0
	 indi(nbdi)=jma
	 nbdi=nbdi+1
	 nbdp(nbdi)=nbdp(nbdi-1)
	 nbdm(nbdi)=nbdm(nbdi-1)
	 finp=.false.
	 finm=.false.
       end if

       if(trans(jma).eq.7.or.trans(jma).eq.8) then
         nbdp(nbdi)=nbdp(nbdi)+1
         indp(nbdp(nbdi))=jma
	 hdp(nbdp(nbdi))=trans(jma)-6
	 finp=.true.
       end if


       if(trans(jma).eq.9.or.trans(jma).eq.10) then
         nbdm(nbdi)=nbdm(nbdi)+1
         indm(nbdm(nbdi))=jma
	 hdm(nbdm(nbdi))=trans(jma)-8
	 finm=.true.
       end if

       if(((.not.finp).or.(.not.finm)).and.jma.lt.dataset%map%nmk(c)) go to 10

!
!  on cherche les zones � gauche du qtl
!  les elements sont comme � droite , en remplacant d par g
!
   15  nbgi=1
       do lk=1,dataset%map%nmk(c)
         nbgp(lk)=0
	 nbgm(lk)=0
       end do
       jma=lqtl+1
       finp=.false.
       finm=.false.

   20  jma=jma-1
       if(jma.eq.0) go to 25
       if(trans(jma).le.4) then
         nbgp(nbgi)=nbgp(nbgi)+1
         ingp(nbgp(nbgi))=jma
	 hgp(nbgp(nbgi))=floor(real(trans(jma)+1)/2.)
         nbgm(nbgi)=nbgm(nbgi)+1
         ingm(nbgm(nbgi))=jma
	 hgm(nbgm(nbgi))=trans(jma)-2*(hgp(nbgp(nbgi))-1)
	 finp=.true.
	 finm=.true.
       end if

       if(trans(jma).eq.5.or.trans(jma).eq.6) then
         nbgp(nbgi)=nbgp(nbgi)+1
         ingp(nbgp(nbgi))=jma
	 hgp(nbgp(nbgi))=0
         nbgm(nbgi)=nbgm(nbgi)+1
         ingm(nbgm(nbgi))=jma
	 hgm(nbgm(nbgi))=0
	 ingi(nbgi)=jma
	 nbgi=nbgi+1
         nbgp(nbgi)=nbgp(nbgi-1)
         nbgm(nbgi)=nbgm(nbgi-1)
	 finp=.false.
	 finm=.false.
       end if

       if(trans(jma).eq.7.or.trans(jma).eq.8) then
         nbgp(nbgi)=nbgp(nbgi)+1
         ingp(nbgp(nbgi))=jma
	 hgp(nbgp(nbgi))=trans(jma)-6
	 finp=.true.
       end if


       if(trans(jma).eq.9.or.trans(jma).eq.10) then
         nbgm(nbgi)=nbgm(nbgi)+1
         ingm(nbgm(nbgi))=jma
	 hgm(nbgm(nbgi))=trans(jma)-8
	 finm=.true.
       end if

       if(((.not.finp).or.(.not.finm)).and.jma.gt.1) go to 20

!
!   colocgp et colocgm indique (true) si la position du qtl est colocalis�e
!   avec un marqueur (coloc = true) et que ce marqueur est le premier �
!   sa gauche. Dans ce cas, il n'y a pas de recombinaison possible entre
! le qtl et ce marqueur
!
!  hlocp et hlocm sont les informations sur la transmission en ce locus (� la
!  fois marqueur et qtl)
!
   25   colocgp=.false.
        colocgm=.false.
        colocdp=.false.
        colocdm=.false.
	if(coloc.and.nbgp(1).ne.0)then
	  if(int(1000.d0*dx).eq.int(1000.d0*dataset%map%posi(c,ingp(1))))colocgp=.true.
	end if
	if(coloc.and.nbgm(1).ne.0)then
	  if(int(1000.d0*dx).eq.int(1000.d0*dataset%map%posi(c,ingm(1))))colocgm=.true.
	end if
	if(coloc.and.nbdp(1).ne.0)then
	  if(int(1000.d0*dx).eq.int(1000.d0*dataset%map%posi(c,indp(1))))colocdp=.true.
	end if
	if(coloc.and.nbdm(1).ne.0)then
	  if(int(1000.d0*dx).eq.int(1000.d0*dataset%map%posi(c,indm(1))))colocdm=.true.
	end if

!
!  les taux de recombinaison avant et apres le QTL sont stock�s dans
!   (t ou s)rec(d ou g)(p ou m)
!
        idp=indp(1)
        if(nbdp(1).eq.0)idp=lqtl+1
        igp=ingp(1)
        if(nbgp(1).eq.0)igp=lqtl

	xint=dataset%map%posi(c,lqtl+1)-dataset%map%posi(c,lqtl)
        xintp=dataset%map%posim(c,lqtl+1)-dataset%map%posim(c,lqtl)
!
!  dans le cas de colocalisation, il ne faut pas calculer le taux de
!  recombinaison
!
        dgp=(dx-dataset%map%posi(c,igp))*xintp/xint
        recgp=xaldane(dgp)
	ddp=(dataset%map%posi(c,idp)-dx)*xintp/xint
	recdp=xaldane(ddp)
!
!
!  cas des m�res
!
        idm=indm(1)
        if(nbdm(1).eq.0)idm=lqtl+1
        igm=ingm(1)
        if(nbgm(1).eq.0)igm=lqtl

	xint=dataset%map%posi(c,lqtl+1)-dataset%map%posi(c,lqtl)
	xintm=dataset%map%posif(c,lqtl+1)-dataset%map%posif(c,lqtl)
!
!  dans le cas de colocalisation, il ne faut pas calculer le taux de
!  recombinaison
!
         dgm=(dx-dataset%map%posi(c,igm))*xintm/xint
	 recgm=xaldane(dgm)
	 ddm=(dataset%map%posi(c,idm)-dx)*xintm/xint
         recdm=xaldane(ddm)

!
!  calcul des probabilit�s de transmission pour le descendant kd, position dx
!
!
!  on va envisager toutes les possibilit�s pour les marqueurs incertains
!  il y a 2**(nbdi+nbgi) possibilit�s
!  ainsi que les 4 possibilit�s pour la transmission du QTL
!
!  les hdp, hdm, hgp et hgm sont compl�t�s par des 1 ou 2 pour les marqueurs
!  incertains
!
!  les boucles en jdi et jgi ont pour fonction de g�n�rer les 2**(nbdi+nbgi)
!  vecteurs d'origines correspondants aux marqueurs incertains
!
!  initialisation
      pbtt=0.d0
      do lq=1,4
	pbtg(lq)=0.d0
	pbtd(lq)=0.d0
      end do

!
!  on va consid�rer les 4 �v�nements de transmission possible au QTL, d�finis
!  pas iqtlp et iqtlm
!
!  pour chacun des �venements de transmision au qtl on v�rifie la coh�rence
!  avec un �ventuel marqueur colocalis�
!
      iqtlp=0
   30 iqtlp=iqtlp+1
      if(iqtlp.eq.3)go to 120

      iqtlm=0
   40 iqtlm=iqtlm+1
      if(iqtlm.eq.3)go to 30
!
      iqtl= iqtlm+2*(iqtlp-1)

!
!  le cas o� il n'y a pas d'incertitude � droite est trait� � part
!  dans ce cas, il faut distinguer les situations o� il n'y a que des
!  transmissions inconnues
!
      if(nbdi.eq.1) then
!
!  on a nbdp=0 si aucun marqueur � droite n'a trans diff�rent de 9, 10 ou 11
!  on a nbdm=0 si aucun marqueur � droite n'a trans diff�rent de 7,  8 ou 11
!  dans ce cas,les 2 �v�nements iqtlp = 1 ou 2 sont �quiprobables
!
         if(nbdp(1).eq.0.and.nbdm(1).eq.0) then
	   pbtd(iqtl)=0.25d0
         end if

                               prob=1.d0
         if(nbdm(1).ne.0) then
	   if(colocdm.and.iqtlm.ne.hdm(1))go to 50
	   if(.not.colocdm)then
	     if(nbdp(1).eq.0)   prob=0.5d0
             if(iqtlm.eq.hdm(1))prob=prob*(1.d0-recdm)
             if(iqtlm.ne.hdm(1))prob=prob*recdm
	   end if
	 end if

          if(nbdp(1).ne.0) then
	   if(colocdp.and.iqtlp.ne.hdp(1))go to 50
	   if(.not.colocdp)then
	     if(nbdm(1).eq.0)   prob=0.5d0
             if(iqtlp.eq.hdp(1))prob=prob*(1.d0-recdp)
             if(iqtlp.ne.hdp(1))prob=prob*recdp
	   end if
	 end if

         pbtd(iqtl)=prob
	 go to 50

       end if
!
!  cas ou il y a des incertitudes � droite
!
!
!
!  on traite diff�remment les intervalles entre marqueurs incertains : le 1er, les
!  suivants, le dernier
!
!  cas du premier intervalle
!  nfin est l'indice du marqueur incertain � droite
!
         idi=1
         possible=.false.
         do hdifin=1,2
!
!  on a deux possibilit�s couvrant l'incertitude, donnant les derniers
!  hdp et hdm du premier intervalle
!
	   nfin=indi(idi)
	   hdp(nbdp(1))=new(trans(nfin),2*hdifin-1)
	   hdm(nbdm(1))=new(trans(nfin),2*hdifin)

	   if((colocdm.and.iqtlm.ne.hdm(1)).or.(colocdp.and.iqtlp.ne.hdp(1)))then

              pbante(iqtl,hdifin)=0.d0

	   else

	     possible=.true.
             prob=0.d0

             if(.not.colocdm)then
               if(iqtlm.eq.hdm(1))prob=dlog(1.d0-recdm)
               if(iqtlm.ne.hdm(1))prob=dlog(recdm )

               do jma=2,nbdm(1)
                 if(hdm(jma-1).eq.hdm(jma))then
	           prob=prob+srecf(indm(jma-1),indm(jma))
	         else
                   prob=prob+trecf(indm(jma-1),indm(jma))
	         end if
               end do
	     end if

             if(.not.colocdp)then
               if(iqtlp.eq.hdp(1))prob=prob+dlog(1.d0-recdp)
               if(iqtlp.ne.hdp(1))prob=prob+dlog(recdp)
               do jma=2,nbdp(1)
                 if(hdp(jma-1).eq.hdp(jma))then
	           prob=prob+srecm(indp(jma-1),indp(jma))
	         else
                   prob=prob+trecm(indp(jma-1),indp(jma))
	         end if
               end do
	     end if

             pbante(iqtl,hdifin)=dexp(prob)

           end if

         end do

	 if(.not.possible) go to 50
!
!  cas des intervalles suivants
!
         do idi=2,nbdi-1

           do hdifin=1,2
!
!  on a deux possibilit�s couvrant l'incertitude, donnant les derniers
!  hdp et hdm du idi �me intervalle
!
	    nfin=indi(idi)
	    hdp(nbdp(idi))=new(trans(nfin),2*hdifin-1)
	    hdm(nbdm(idi))=new(trans(nfin),2*hdifin)
            pcum(hdifin)=0.d0

             do hdideb=1,2
!
!  il faut consid�rer les  deux possibilit�s couvrant l'incertitude,
!  de l'intervalle pr�c�dent
!
	      ndeb=indi(idi-1)
	      hdp(nbdp(idi-1))=new(trans(ndeb),2*hdideb-1)
	      hdm(nbdm(idi-1))=new(trans(ndeb),2*hdideb)

              prob=0.d0
	      do jma=nbdp(idi-1)+1,nbdp(idi)
                if(hdp(jma-1).eq.hdp(jma))then
	          prob=prob+srecm(indp(jma-1),indp(jma))
	        else
                  prob=prob+trecm(indp(jma-1),indp(jma))
	        end if
              end do

	      do jma=nbdm(idi-1)+1,nbdm(idi)
                if(hdm(jma-1).eq.hdm(jma))then
	          prob=prob+srecf(indm(jma-1),indm(jma))
	        else
                  prob=prob+trecf(indm(jma-1),indm(jma))
	        end if
              end do
!
	      pcum(hdifin)=pcum(hdifin)+pbante(iqtl,hdideb)*dexp(prob)

             end do
	   end do

	   do hdifin=1,2
	     pbante(iqtl,hdifin)=pcum(hdifin)
	   end do

	 end do

!
!  cas du dernier intervalle
!
         idi=nbdi
	 pcum(1)=0.d0
         do hdideb=1,2
!
!  il faut consid�rer les  deux possibilit�s couvrant l'incertitude,
!  de l'intervalle pr�c�dent
!
           ndeb=indi(idi-1)
           hdp(nbdp(idi-1))=new(trans(ndeb),2*hdideb-1)
           hdm(nbdm(idi-1))=new(trans(ndeb),2*hdideb)

           prob=0.d0
           if(nbdp(idi).ne.nbdp(idi-1))then
             do jma=nbdp(idi-1)+1,nbdp(idi)
               if(hdp(jma-1).eq.hdp(jma))then
        	 prob=prob+srecm(indp(jma-1),indp(jma))
               else
        	 prob=prob+trecm(indp(jma-1),indp(jma))
               end if
             end do
           end if

           if(nbdm(idi).ne.nbdm(idi-1))then
             do jma=nbdm(idi-1)+1,nbdm(idi)
               if(hdm(jma-1).eq.hdm(jma))then
        	 prob=prob+srecf(indm(jma-1),indm(jma))
               else
        	 prob=prob+trecf(indm(jma-1),indm(jma))
               end if
             end do
           end if
!
	   pcum(1)=pcum(1)+pbante(iqtl,hdideb)*dexp(prob)

         end do

         pbtd(iqtl)=pcum(1)

   50  continue
!
!  on �num�re les possibilit�s pour les transmissions incertaines � gauche
!  du qtl.
!
!
!  le cas o� il n'y a pas d'incertitude � gauche est trait� � part
!  dans ce cas, il faut distinguer les situations o� il n'y a que des
!  transmissions inconnues
!

      if(nbgi.eq.1) then
!
!  on a nbgp=0 si aucun marqueur � gauche n'a trans diff�rent de 9, 10 ou 11
!  on a nbgm=0 si aucun marqueur � gauche n'a trans diff�rent de 7,  8 ou 11
!  dans ce cas,les 2 �v�nements iqtlp = 1 ou 2 sont �quiprobables
!
         if(nbgp(1).eq.0.and.nbgm(1).eq.0) then
	   pbtg(iqtl)=0.25d0
         end if

                                prob=1.d0
         if(nbgm(1).ne.0) then
	   if(colocgm.and.iqtlm.ne.hgm(1))go to 60
	   if(.not.colocgm)then
	     if(nbgp(1).eq.0)   prob=0.5d0
             if(iqtlm.eq.hgm(1))prob=prob*(1.d0-recgm)
             if(iqtlm.ne.hgm(1))prob=prob*recgm
	   end if
	 end if

         if(nbgp(1).ne.0) then
	   if(colocgp.and.iqtlp.ne.hgp(1))go to 60
	   if(.not.colocgp)then
	     if(nbgm(1).eq.0)   prob=0.5d0
             if(iqtlp.eq.hgp(1))prob=prob*(1.d0-recgp)
             if(iqtlp.ne.hgp(1))prob=prob*recgp
	   end if
	 end if

         pbtg(iqtl)=prob
 	 go to 60

       end if
!
!  cas ou il y a des incertitudes � gauche
!
!
!
!  on traite diff�remment les intervalles entre marqueurs incertains : le 1er, les
!  suivants, le dernier
!
!  cas du premier intervalle
!  nfin est l'indice du marqueur incertain � gauche
!
         igi=1
         possible=.false.
         do hgifin=1,2

!
!  on a deux possibilit�s couvrant l'incertitude, donnant les derniers
!  hgp et hgm du premier intervalle
!
	   nfin=ingi(igi)
	   hgp(nbgp(1))=new(trans(nfin),2*hgifin-1)
	   hgm(nbgm(1))=new(trans(nfin),2*hgifin)

	   if(  (colocgm.and.iqtlm.ne.hgm(1)).or.(colocgp.and.iqtlp.ne.hgp(1)))then

              pbante(iqtl,hgifin)=0.d0

	   else

	    possible=.true.
             prob=0.d0

             if(.not.colocgm)then
               if(iqtlm.eq.hgm(1))prob=dlog(1.d0-recgm)
               if(iqtlm.ne.hgm(1))prob=dlog(recgm)

               do jma=2,nbgm(1)
         	 if(hgm(jma-1).eq.hgm(jma))then
         	   prob=prob+srecf(ingm(jma-1),ingm(jma))
         	 else
         	   prob=prob+trecf(ingm(jma-1),ingm(jma))
         	 end if
               end do
             end if

             if(.not.colocgp)then
               if(iqtlp.eq.hgp(1))prob=prob+dlog(1.d0-recgp)
               if(iqtlp.ne.hgp(1))prob=prob+dlog(recgp)
               do jma=2,nbgp(1)
         	 if(hgp(jma-1).eq.hgp(jma))then
         	   prob=prob+srecm(ingp(jma-1),ingp(jma))
         	 else
         	   prob=prob+trecm(ingp(jma-1),ingp(jma))
         	 end if
               end do
             end if

             pbante(iqtl,hgifin)=dexp(prob)

	   end if

         end do

	 if(.not.possible)go to 60
!
!  cas des intervalles suivants
!
         do igi=2,nbgi-1

           do hgifin=1,2
!
!  on a deux possibilit�s couvrant l'incertitude, donnant les derniers
!  hgp et hgm du igi �me intervalle
!
	    nfin=ingi(igi)
	    hgp(nbgp(igi))=new(trans(nfin),2*hgifin-1)
	    hgm(nbgm(igi))=new(trans(nfin),2*hgifin)
            pcum(hgifin)=0.d0

             do hgideb=1,2
!
!  il faut consid�rer les  deux possibilit�s couvrant l'incertitude,
!  de l'intervalle pr�c�dent
!
	      ndeb=ingi(igi-1)
	      hgp(nbgp(igi-1))=new(trans(ndeb),2*hgideb-1)
	      hgm(nbgm(igi-1))=new(trans(ndeb),2*hgideb)

              prob=0.d0
	      do jma=nbgp(igi-1)+1,nbgp(igi)
                if(hgp(jma-1).eq.hgp(jma))then
	          prob=prob+srecm(ingp(jma-1),ingp(jma))
	        else
                  prob=prob+trecm(ingp(jma-1),ingp(jma))
	        end if
              end do

	      do jma=nbgm(igi-1)+1,nbgm(igi)
                if(hgm(jma-1).eq.hgm(jma))then
	          prob=prob+srecf(ingm(jma-1),ingm(jma))
	        else
                  prob=prob+trecf(ingm(jma-1),ingm(jma))
	        end if
              end do
!
	      pcum(hgifin)=pcum(hgifin)+pbante(iqtl,hgideb)*dexp(prob)

             end do
	   end do

	   do hgifin=1,2
	     pbante(iqtl,hgifin)=pcum(hgifin)
	   end do

	 end do

!
!  cas du dernier intervalle
!
         igi=nbgi
	 pcum(1)=0.d0
         do hgideb=1,2
!
!  il faut consid�rer les  deux possibilit�s couvrant l'incertitude,
!  de l'intervalle pr�c�dent
!
           ndeb=ingi(igi-1)
           hgp(nbgp(igi-1))=new(trans(ndeb),2*hgideb-1)
           hgm(nbgm(igi-1))=new(trans(ndeb),2*hgideb)

           prob=0.d0
           if(nbgp(igi).ne.nbgp(igi-1))then
             do jma=nbgp(igi-1)+1,nbgp(igi)
               if(hgp(jma-1).eq.hgp(jma))then
        	 prob=prob+srecm(ingp(jma-1),ingp(jma))
               else
        	 prob=prob+trecm(ingp(jma-1),ingp(jma))
               end if
             end do
           end if

           if(nbgm(igi).ne.nbgm(igi-1))then
             do jma=nbgm(igi-1)+1,nbgm(igi)
               if(hgm(jma-1).eq.hgm(jma))then
        	 prob=prob+srecf(ingm(jma-1),ingm(jma))
               else
        	 prob=prob+trecf(ingm(jma-1),ingm(jma))
               end if
             end do
           end if
!
	   pcum(1)=pcum(1)+pbante(iqtl,hgideb)*dexp(prob)

         end do

         pbtg(iqtl)= pcum(1)

   60 continue
!
!  on teste si d'autres combinaisons au qtl sont � �valuer
!
       if(iqtlp.eq.2.and.iqtlm.eq.2) go to 120
       if(iqtlm.eq.1)go to 40
       if(iqtlm.eq.2)go to 30


  120  continue

!
!  on stocke la probabilit� de transmission � la position dx
!
         do lq=1,4
	   pbt(lq)=pbtg(lq)*pbtd(lq)
	 end do

  130   continue

	 pbtt=0.d0
         do lq=1,4
	   pbtt=pbtt+pbt(lq)
	 end do

       if ( pbtt /= 0 ) then
       do lq=1,4
          spt%pdd(c,kd,lq,n)=pbt(lq)/pbtt
	   end do
	   end if

!
!
!  on avance d'un pas (en verifiant qu'on ne depasse pas le bout droit) et on recommence
!
         n=n+1
     if ( n > dataset%map%get_npo(c) ) return
     dx=dataset%map%absi(c,n)
	 if( dx.gt.dataset%map%posi(c,dataset%map%nmk(c)) ) return
	 go to 100

      return
   end subroutine pded_v5_kd
!!***

!!****f* m_qtlmap_haplotype_V2/set_phasp
!! NAME
!!    set_phasp
!! DESCRIPTION
!!
!! NOTES
!!
!! SOURCE
    subroutine set_phasp(dataset,spt)
       type(QTLMAP_DATASET)       ,intent(in)            :: dataset
       type(PDD_BUILD)            ,intent(inout)         :: spt
       integer                                :: ll,ip,c
       type(GENEALOGY_BASE) , pointer :: dgenea

       dgenea => dataset%genea

       if (size(spt%phasp,2) /= size(ordrep,2)) then
         call stop_application("Devel error: using set_phasp with dim(phasp) :dim(ordrep,2)")
       end if
      spt%phasp=.false.
      do c=1,dataset%map%nchr
       do ll=1,dataset%map%nmk(c)
        do ip=1,dgenea%np
           if(ordrep(c,ip,ll) /= 0 ) then
               spt%phasp(c,ip)=.true.
           end if
        end do
       end do
      end do
    end subroutine set_phasp
!!***

!!****f* m_qtlmap_haplotype_V2/set_phasm
!! NAME
!!    set_phasm
!! DESCRIPTION
!!
!! NOTES
!!
!! SOURCE
    subroutine set_phasm(dataset,spt)
       type(QTLMAP_DATASET)       ,intent(in)            :: dataset
       type(PDD_BUILD)            ,intent(inout)         :: spt
       integer                                :: ll,im,c
       type(GENEALOGY_BASE) , pointer :: dgenea

       dgenea => dataset%genea

       if (size(spt%phasm,2) /= size(ordrem,2)) then
         call stop_application("Devel error: using set_phasm with dim(phasm) :dim(ordrem,2)")
       end if
      spt%phasm=.false.
      do c=1,dataset%map%nchr
       do ll=1,dataset%map%nmk(c)
        do im=1,dgenea%nm
           if(ordrem(c,im,ll) /= 0) spt%phasm(c,im)=.true.
        end do
       end do
      end do
    end subroutine set_phasm
!!***

!!****f* m_qtlmap_haplotype_V2/check_recombination_sire
!! NAME
!!    check_recombination_sire
!! DESCRIPTION
!!
!! NOTES
!!   verification de la coherence entre la carte et les recombinaisons observees intra pere
!!
!!check_recombination_sire  est  construite sur le mode de haplotype_snp.
!!Entre chaque marqueur pour lequel le pere est informative (heterozygote), on compte le nombre de transmission “en phase”
!! (meme origine grand parentale) et en opposition. Le taux de recombinaison observe est  le rapport nbr en opposition / nbr total>.
!! Connaissant ( d’apres la carte) la distance entre ces deux marqueurs ( donc le taux de recombinaison theorique) , !
!! on peut calculer la probabilite que ce taux de recombinaison observe soit obtenu
!! ( la loi est binomial de parameters  (nombre total de descendants dont on connait la transmission aux deux marqueurs
!! , taux de recombinaison theorique) . Si cette proba est trop faible, on lance un warning
!! SOURCE
      subroutine check_recombination_sire(dataset,spt)
       type(QTLMAP_DATASET)       ,intent(in)            :: dataset
       type(PDD_BUILD)            ,intent(inout)         :: spt

! Divers
      integer(kind=KIND_PHENO) ::  m1,m2
      integer imark,ip,jm,kd,c
      integer ll
      integer ,allocatable,dimension(:,:) :: trans
      integer oppos,phase,total,ig,ig2,ll1,ll2,sb,nb,nk2
      integer nd1,nd2,nm1,nm2,io,iip,opposAll,totalAll,stat,unknown,llkk
      real (kind=dp) :: p1,p2,recomb,pp,pr
      real(kind=8) :: bound,p,q,s,xn,ompr
      logical :: ischanged
      logical,dimension(:,:),allocatable :: hetero,countFalseMark
      logical, dimension(:),allocatable :: checkblock
      integer,dimension(dataset%map%nchr) :: nk

      integer, dimension(:,:),allocatable :: tab
      integer, dimension(:),allocatable :: tab2
      real (kind=dp), dimension(:,:),allocatable :: recombPP 
      character(len=20) :: ext
      type(GENEALOGY_BASE) , pointer :: dgenea
      type(GENOTYPE_BASE) , pointer :: dga

      dga => dataset%genoAnimal
      dgenea => dataset%genea

      call log_mess('checking recombination in sires...',INFO_DEF)
!
      allocate (trans(dgenea%nd,maxval(dataset%map%nmk)))
      allocate (hetero(dgenea%np,maxval(dataset%map%nmk)))
      allocate (countFalseMark(dataset%map%nchr,maxval(dataset%map%nmk)))
      allocate (checkblock(maxval(dataset%map%nmk)))
      allocate (tab(dataset%map%nchr,maxval(dataset%map%nmk)),tab2(maxval(dataset%map%nmk)))

      nk=0
      countFalseMark = .false.
      do c=1,dataset%map%nchr
       trans=0
       do ip=1,dgenea%np
        nm1=dgenea%nmp(ip)+1
        nm2=dgenea%nmp(ip+1)
!
!  on cree le tableau hetero qui liste l'heterozygotie des marqueurs
!
        do imark=1,dataset%map%nmk(c)
          hetero(ip,imark)=.true.
          if(dga%correp(ip).eq.9999) then
            m1=dga%nmanque
          else
            m1=dga%pheno(c,imark,dga%correp(ip),1)
            m2=dga%pheno(c,imark,dga%correp(ip),2)
          end if
          if(m1.eq.dga%nmanque) m2=m1
          if(m1.eq.m2) then
            hetero(ip,imark)=.false.
          end if
        end do ! imark
!
!  on reduit le jeu de marqueurs a ceux pour lesquels le pere est heterozygote
!
!  creation du tableau des transmissions observees
!
        do jm=nm1,nm2
          nd1=dgenea%ndm(jm)+1
          nd2=dgenea%ndm(jm+1)
          do kd=nd1,nd2
            do ll=1,dataset%map%nmk(c)
              p1=spt%prot(c,ll,kd,1)+spt%prot(c,ll,kd,2)
              if(p1.eq.1.d0)trans(kd,ll)=1
              p2=spt%prot(c,ll,kd,3)+spt%prot(c,ll,kd,4)
              if(p2.eq.1.d0)trans(kd,ll)=2
            end do !! ll
          end do!kd
        end do!jm
     end do ! ip
!
!  comptage des inversions apparentes des phases locales
!

     ischanged=.true.
     nk(c)=dataset%map%nmk(c)
     countFalseMark(c,:dataset%map%nmk(c))=.true. !all marker is valid !

   do ll=1,nk(c)
    tab(c,ll)=ll
   end do
   do while ( ischanged )
     ischanged=.false.
     write (*,*) "chromosome ",trim(dataset%map%chromo(c))," => nb marker:",nk(c)
    do ll=1,nk(c)-1
      totalAll=0
      opposAll=0
      iip=0
      
      ll1 = tab(c,ll)
      ll2 = tab(c,ll+1)
    
      !print *,"************* LL1,LL2=",ll1,ll2,"*****************************",trim(mark(c,ll))
      do ip=1,dgenea%np
       if ( hetero(ip,ll1) .and. hetero(ip,ll2)) then
        phase=0
	oppos=0
        do jm=dgenea%nmp(ip)+1,dgenea%nmp(ip+1)
	 do kd=dgenea%ndm(jm)+1,dgenea%ndm(jm+1)
           if(trans(kd,ll1).ne.0.and.trans(kd,ll2).ne.0) then
	     if(trans(kd,ll1).eq.trans(kd,ll2))then
		 phase=phase+1
	       else
		 oppos=oppos+1
	       end if
	     end if
          end do !kd
         end do !jm 

         total=phase+oppos
         if(phase < oppos)oppos=phase
         !comptage globale
         totalAll=totalAll+total
         opposAll=opposAll+oppos
         iip=iip+1
         !print *,ip,total,totalAll
        end if
       end do ! ip

       if ( opposAll == 0 ) cycle
         
       recomb=xaldane(dataset%map%posim(c,ll2)-dataset%map%posim(c,ll1))
      ! call probin(totalAll,recomb,opposAll,pp)
       s    = dble(opposAll - 1)
       xn   = dble(totalAll)
       pr   = recomb
       ompr = 1.d0 - recomb
       
       call cdfbin ( 1, p, q, s, xn, pr, ompr, stat, bound )

       if ( stat /= 0 ) then
         call stop_application("error with cdfbin:"//trim(str(stat))//" bound:"//trim(str(bound)))
       end if

       pp = q
       
       if(pp < dataset%params%prob_seuil_recomb) then
        countFalseMark(c,ll)=.false. ! inconsistent   
       end if
      end do ! ll
      !!-------------------------------------------------
      !! VERIFICATION SI ON ENLEVE DES MARQUEURS QUI ON TROP DE RECOMB A DROITE ET A GAUCHE
      !! SI VRAI => ON RECOMENCE LE PROCESS JUSQU4A STABILITE

      nk2=1
      tab2(nk2)=1

      do ll=2,nk(c)
       !print *,countFalseMark(c,ll-1),countFalseMark(c,ll)
       if ( (countFalseMark(c,ll-1) .or.  countFalseMark(c,ll)) ) then
         nk2=nk2+1
         tab2(nk2)=tab(c,ll)
       else
         ischanged=.true.
         write (*,*) "remove marker:",tab(c,ll),dataset%map%mark(c,tab(c,ll))
       end if
      end do

      tab(c,:nk2)=tab2(:nk2)
      countFalseMark(c,:nk(c))=.true.
      nk(c)=nk2

    end do !while
   end do !nchr

  if (sum(nk) /= sum(dataset%map%nmk)) then
      ext=".new"

      call log_mess("A new map file has been created :"//&
       trim(dataset%params%get_file_val(K_MAP))//trim(ext)//" number of snp to remove :"//&
       trim(str(sum(dataset%map%nmk)-sum(nk))),INFO_DEF)
      open(unit=987,file=trim(dataset%params%get_file_val(K_MAP))//trim(ext), form="formatted",iostat=io)
      if ( io /= 0 ) then
       call stop_application("can not create file"//trim(dataset%params%get_file_val(K_MAP))//trim(ext))
      end if
      do c=1,dataset%map%nchr
       ll=1
       do llkk=1,dataset%map%nmk(c)
        
        if ( tab(c,ll) /= llkk ) then
           write(987,fmt="(1a50,1x,a20,1x,f10.8,1x,f10.8,1x,f10.8,1x,i5)")&
            trim(dataset%map%mark(c,llkk)),trim(dataset%map%chromo(c)),&
                                      dataset%map%posi(c,llkk),dataset%map%posim(c,llkk),dataset%map%posif(c,llkk),0
        else
!            write(987,fmt="(1a50,1x,a20,1x,f10.8,1x,f10.8,1x,f10.8,1x,i5)") trim(mark(c,tab(c,ll))),trim(chromo(c)),&
!                                            posi(c,tab(c,ll)),posim(c,tab(c,ll)),posif(c,tab(c,ll)),1
            write(987,fmt="(1a50,1x,a20,1x,f10.8,1x,f10.8,1x,f10.8,1x,i5)") &
            trim(dataset%map%mark(c,llkk)),trim(dataset%map%chromo(c)),&
                                      dataset%map%posi(c,llkk),dataset%map%posim(c,llkk),dataset%map%posif(c,llkk),1
            ll = ll + 1
        end if
       !  else
          !  write(nficout,*)'Too many recombination marker (inconsistent with the next marker): '//trim(mark(c,ll))
          !  write(nficout,*) recombo(c,ll),'/',recombt(c,ll),'unk:',recombu(c,ll),' - ',recombPP(c,ll),' np=',recombiip(c,ll)
       !  end if
       end do ! ll
      end do! nchr

      close (987)
    else
      call log_mess("map :ok ")
    end if
     
     deallocate (tab,tab2)
     deallocate (trans)
     deallocate (hetero)
     deallocate (countFalseMark)
     call log_mess('end check_recombination_sire_V2',INFO_DEF)
     return
     end subroutine check_recombination_sire

!!****f* m_qtlmap_haplotype_V2/probin
!! NAME
!!    probin
!! DESCRIPTION
!!
!! NOTES
!!  probin calcule la fonction de repartition d une binomial de parametres n et r
!!  pp est la probabilite de trouver au moins m individus recombinants parmi n
!!  pour un tqux de recombinaison r
!! SOURCE
      subroutine probin(n,r,m,pp)
       double precision r,pp,c,d
       integer n,m,l

       if(r**n == 0.d0) then
         if(m == 0) then
           pp = 1.d0
         else
           pp = 0.d0
         end if
         return
       end if

       if(r == 1.d0) then
         if(m == n) then
           pp = 1.d0
         else
           pp = 0.d0
         end if
         return
       end if

       c=1.d0
       d=r**n
       pp=c*d
       do l=n,m,-1
         c=c *  dble(l) / dble(n-l+1)
         d=d*(1.d0-r)/ r
         pp=pp+c*d
       end do
     end subroutine probin
!!***

!!****f* m_qtlmap_haplotype_V2/calcul_phases_symmax2sat_sire
!! NAME
!!    probin
!! DESCRIPTION
!!
!! NOTES
!!
!! SOURCE
     subroutine calcul_phases_symmax2sat_sire(dataset,spt)
       type(QTLMAP_DATASET)       ,intent(in)            :: dataset
       type(PDD_BUILD)            ,intent(inout)         :: spt

       integer :: c,nm1,nm2,m1,m2,ip,imark,kkd,nd1,nd2,ik,ik2,ll,jm,kd,m
       real(kind=dp)  :: p1,p2,wmax
       real(kind=dp) , dimension(:,:),allocatable             :: trans,W,r
       integer       , dimension(:),allocatable                :: H
       integer       , dimension(:,:),allocatable              :: nbp,nbm

       logical :: ok
       integer :: c_11,c_12,c_21,c_22,ll1,ll2,io
       type(GENEALOGY_BASE) , pointer :: dgenea
       type(GENOTYPE_BASE) , pointer :: dga

       dga => dataset%genoAnimal
       dgenea => dataset%genea

       allocate (trans(dgenea%nd,maxval(dataset%map%nmk)+1),stat=io)
       if ( io /= 0 ) call stop_application("call phases symmax2sat : not enough memory...")
       allocate (H(maxval(dataset%map%nmk)),stat=io)
       if ( io /= 0 ) call stop_application("call phases symmax2sat : not enough memory...")
       allocate (nbp(maxval(dataset%map%nmk),maxval(dataset%map%nmk)),&
                 nbm(maxval(dataset%map%nmk),maxval(dataset%map%nmk)),&
                 W(maxval(dataset%map%nmk),maxval(dataset%map%nmk)),&
                 r(maxval(dataset%map%nmk),maxval(dataset%map%nmk)),&
                 stat=io)
       if ( io /= 0 ) call stop_application("call phases symmax2sat : not enough memory...")

       trans=0
!
! Etablissement de la phase la plus probable
!
    do c=1,dataset%map%nchr
      do ik=1,dataset%map%nmk(c)-1
         do ik2=ik+1,dataset%map%nmk(c)
            r(ik,ik2) = log( (1.d0-dataset%map%rm(c,ik,ik2))/dataset%map%rm(c,ik,ik2) )
         end do !ik2
       end do !ik

      do ip=1,dgenea%np

        W=0
!
!  on reduit le jeu de marqueurs e ceux pour lesquels le pere est heterozygote
!
!  creation du tableau des transmissions observees
!
        kkd=0
        do jm=dgenea%nmp(ip)+1,dgenea%nmp(ip+1)
          do kd=dgenea%ndm(jm)+1,dgenea%ndm(jm+1)
            kkd=kkd+1
            do ll=1,dataset%map%nmk(c)
              trans(kkd,ll)=0
              p1=spt%prot(c,ll,kd,1)+spt%prot(c,ll,kd,2)
              if(p1.eq.1.d0)trans(kkd,ll)=1
              p2=spt%prot(c,ll,kd,3)+spt%prot(c,ll,kd,4)
              if(p2.eq.1.d0)trans(kkd,ll)=2
            end do !ll
          end do !kd
        end do !jm

        nbp=0
        nbm=0
        kkd=0
        do jm=dgenea%nmp(ip)+1,dgenea%nmp(ip+1)
          do kd=dgenea%ndm(jm)+1,dgenea%ndm(jm+1)
              kkd=kkd+1
              ik=1
             do while (ik<dataset%map%nmk(c))
              do while (ik<=(dataset%map%nmk(c)-1) .and.  trans(kkd,ik) == 0 )
               ik=ik+1
              end do
              if ( ik >= dataset%map%nmk(c)) exit
              ik2=ik+1
              do while (ik2<=dataset%map%nmk(c) .and.  trans(kkd,ik2) == 0 )
                ik2=ik2+1
              end do
              if ( ik2 > dataset%map%nmk(c)) exit
              if ( trans(kkd,ik) == trans(kkd,ik2) ) then
                    nbp(ik,ik2) = nbp(ik,ik2) + 1
              else
                    nbm(ik,ik2) = nbm(ik,ik2) + 1
              end if

             ik=ik2
            end do !while(ik<nmj(c)
          end do !kd
        end do !jm

        !write(6,*)'ordrep',(ordrep(c,ip,ik),ik=1,dataset%map%nmk(c))

        wmax=0.d0
        do ik=1,dataset%map%nmk(c)-1
         do ik2=ik+1,dataset%map%nmk(c)
            if ( nbp(ik,ik2) == nbm(ik,ik2) ) cycle
            W(ik,ik2) = 0.25d0*dble(nbp(ik,ik2)-nbm(ik,ik2))*r(ik,ik2)
	    if(wmax < dabs(W(ik,ik2)))wmax=dabs(W(ik,ik2))
          end do !ik2
        end do !ik



	do ik=1,dataset%map%nmk(c)
	  if(ordrep(c,ip,ik)==12)W(ik,ik)=wmax+1.d0
	  if(ordrep(c,ip,ik)==21)W(ik,ik)=-wmax-1.d0
	 ! write(6,*)(W(ik,ik2),ik2=1,dataset%map%nmk(c))
	end do

        ok = get_h_from_w(maxval(dataset%map%nmk),W,H,m)
               
! 
!  stockage du genotype
!
        do ll=1,dataset%map%nmk(c)
          if(h(ll) == 1) then
            spt%genotyp(c,ll,dga%correp(ip),1)=dga%pheno(c,ll,dga%correp(ip),1)
            spt%genotyp(c,ll,dga%correp(ip),2)=dga%pheno(c,ll,dga%correp(ip),2)
          else
            spt%genotyp(c,ll,dga%correp(ip),1)=dga%pheno(c,ll,dga%correp(ip),2)
            spt%genotyp(c,ll,dga%correp(ip),2)=dga%pheno(c,ll,dga%correp(ip),1)
          end if 
        end do !ll
!
! Reorganisation du tableau des probabilites de transmission
          do jm=dgenea%nmp(ip)+1,dgenea%nmp(ip+1)
          do kd=dgenea%ndm(jm)+1,dgenea%ndm(jm+1)
              do ll=1,dataset%map%nmk(c)
               if(h(ll).ne.1) then
                  p1=spt%prot(c,ll,kd,1)
                  p2=spt%prot(c,ll,kd,2)
                  spt%prot(c,ll,kd,1)=spt%prot(c,ll,kd,3)
                  spt%prot(c,ll,kd,2)=spt%prot(c,ll,kd,4)
                  spt%prot(c,ll,kd,3)=p1
                  spt%prot(c,ll,kd,4)=p2
                end if
              end do !ll
            end do !kd
          end do !jm
        end do ! ip
       end do ! chromosome

       deallocate(trans,H,nbp,nbm,W,r)

     end subroutine calcul_phases_symmax2sat_sire
!!***

!!****f* m_qtlmap_haplotype_V2/calcul_phases_symmax2sat_dam
!! NAME
!!    probin
!! DESCRIPTION
!!
!! NOTES
!!
!! SOURCE
     subroutine calcul_phases_symmax2sat_dam(dataset,spt)
       type(QTLMAP_DATASET)       ,intent(in)            :: dataset
       type(PDD_BUILD)            ,intent(inout)         :: spt

       integer :: io,c,nm1,nm2,m1,m2,ip,imark,kkd,nd1,nd2,ik,ik2,ll,jm,kd,m,kdIfem,ifem,ll1,ll2
       real(kind=dp)  :: p1,p2,wmax
       real(kind=dp) , dimension(:,:),allocatable     :: trans
       integer       , dimension(:)    ,allocatable   :: H
       integer       , dimension(:,:)  ,allocatable   :: nbp,nbm
       real(kind=dp) , dimension(:,:)  ,allocatable   :: W
       integer       , dimension(:)    ,allocatable   :: ndf,corref
       integer       , dimension(:,:)   ,allocatable  :: repdes
       integer       , dimension(:,:)   ,allocatable  :: phm
       real(kind=dp) , dimension(:,:)   ,allocatable  :: r
       logical :: ok

       type(GENEALOGY_BASE) , pointer :: dgenea
       type(GENOTYPE_BASE) , pointer :: dga

       dga => dataset%genoAnimal
       dgenea => dataset%genea

       allocate (trans(dgenea%nd,maxval(dataset%map%nmk)+1),&
                 H(maxval(dataset%map%nmk)),&
                 nbp(maxval(dataset%map%nmk),maxval(dataset%map%nmk)),&
                 nbm(maxval(dataset%map%nmk),maxval(dataset%map%nmk)),&
                 W(maxval(dataset%map%nmk),maxval(dataset%map%nmk)),&
                 ndf(dgenea%nfem),&
                 corref(dgenea%nfem),&
                 repdes(dgenea%nfem,dgenea%ndm(dgenea%nm+1)),&
                 phm(dgenea%nm,maxval(dataset%map%nmk)),&
                 r(maxval(dataset%map%nmk),maxval(dataset%map%nmk)),&
                 stat=io)
       if ( io/= 0 ) call stop_application("call phases symmax2sat : not enough memory...")

       trans = 0

      repdes=0
      ndf=0
      W=0.d0
      do jm=1,dgenea%nm
        ifem=dgenea%repfem(jm)
        do kd=dgenea%ndm(jm)+1,dgenea%ndm(jm+1)
          ndf(ifem)=ndf(ifem)+1
          repdes(ifem,ndf(ifem))=kd
        end do
        corref(ifem)=dga%correm(jm)
      end do

!******************************************************************************
!
!  on traite les femelles (et non les meres) une a une
!
!  nombfem est le nombre de femelles dont on peut chercher les phases
!
    do c=1,dataset%map%nchr
      do ik=1,dataset%map%nmk(c)-1
         do ik2=ik+1,dataset%map%nmk(c)
            r(ik,ik2) = log( (1.d0-dataset%map%rm(c,ik,ik2))/dataset%map%rm(c,ik,ik2) )
         end do
      end do

      do ifem=1,dgenea%nfem
      if(dga%estfem(ifem))then
         kkd=0
         do kdIfem=1,ndf(ifem)
!
!******************************************************************************
!
! Si le nombre de pleins freres est trop faible, le genotype de la mere ne sera pas considere
! Dans ce cas le premier genotype rencontre possible est considere comme le bon

        !if(ndf(ifem).lt.ndmin) then
         ! force(jm)=.true.

!
!  creation du tableau des transmissions observees
!
           kd = repdes(ifem,kdIfem)
   	       kkd=kkd+1
           do ll=1,dataset%map%nmk(c)
            trans(kkd,ll)=0
            p1=spt%prot(c,ll,kd,1)+spt%prot(c,ll,kd,3)
            if(p1.eq.1.d0)trans(kkd,ll)=1
            p2=spt%prot(c,ll,kd,2)+spt%prot(c,ll,kd,4)
            if(p2.eq.1.d0)trans(kkd,ll)=2
!             if ( ll == 1 .and. ifem == 4 ) then
!
!               print *,'kkd:',kkd,' kd:',kd,'p1:',p1,' p2:',p2,'=>',trans(kkd,ll),' prot:',spt%prot(c,ll,kd,:),animal(kd),&
!            get_dga%pheno(c,dga%pheno(c,ll,corred(kd),1)),get_dga%pheno(c,dga%pheno(c,ll,corred(kd),2))
!              end if
          end do ! fin ll
         end do ! fin kdIfem
         !*******************************

        nbp=0
        nbm=0
        kkd=0

        do kdIfem=1,ndf(ifem)
            kd = repdes(ifem,kdIfem)
            kkd=kkd+1
            ik=1
            do while (ik<dataset%map%nmk(c))
              do while (ik<=(dataset%map%nmk(c)-1) .and.  trans(kkd,ik) == 0 )
               ik=ik+1
              end do
              if ( ik >= dataset%map%nmk(c)) exit
              ik2=ik+1
              do while (ik2<=dataset%map%nmk(c) .and.  trans(kkd,ik2) == 0 )
                ik2=ik2+1
              end do
              if ( ik2 > dataset%map%nmk(c)) exit
              if ( trans(kkd,ik) == trans(kkd,ik2) ) then
                    nbp(ik,ik2) = nbp(ik,ik2) + 1
              else
                    nbm(ik,ik2) = nbm(ik,ik2) + 1
              end if

             ik=ik2
            end do
         end do

          wmax=0.d0
          W=0.d0
          do ik=1,dataset%map%nmk(c)-1
           do ik2=ik+1,dataset%map%nmk(c)
            if ( nbp(ik,ik2) == nbm(ik,ik2) ) cycle
              W(ik,ik2) = 0.25d0*dble(nbp(ik,ik2)-nbm(ik,ik2))*r(ik,ik2)
	    if(wmax < dabs(W(ik,ik2)))wmax=dabs(W(ik,ik2))
          end do !ik2
        end do !ik
	
	do ik=1,dataset%map%nmk(c)
	  if(ordref(c,ifem,ik)==12)W(ik,ik)=wmax+1.d0
	  if(ordref(c,ifem,ik)==21)W(ik,ik)=-wmax-1.d0
	 ! write(6,*)(W(ik,ik2),ik2=1,dataset%map%nmk(c))
	end do

          ok = get_h_from_w(maxval(dataset%map%nmk),W,H,m)
        !
        !
        !   stockage des resultats par mere (et non par femelle)
        !    et reperage de la phase la plus probable
        !


!
!      Update the genotypm structure
        do ll=1,dataset%map%nmk(c)
          if(h(ll) == 1) then
           do jm=1,dgenea%nm
            if ( dgenea%repfem(jm) == ifem ) then
              spt%genotypm(c,ll,jm,1)=dga%pheno(c,ll,corref(ifem),1)
              spt%genotypm(c,ll,jm,2)=dga%pheno(c,ll,corref(ifem),2)
              phm(jm,ll)=1
            end if
           end do
          else
            do jm=1,dgenea%nm
            if ( dgenea%repfem(jm) == ifem ) then
              spt%genotypm(c,ll,jm,1)=dga%pheno(c,ll,corref(ifem),2)
              spt%genotypm(c,ll,jm,2)=dga%pheno(c,ll,corref(ifem),1)
              phm(jm,ll)=2
            end if
           end do
          end if
        end do ! ll
        end if
       end do ! ifem
      end do ! chromosome

      allocate ( spt%ngend(dataset%map%nchr,dgenea%nm+1) )
      allocate ( spt%probg(dataset%map%nchr,dgenea%nm) )
      allocate ( spt%ndesc(dataset%map%nchr,dgenea%nd) )

      do c=1,dataset%map%nchr
       spt%ngenom(c,1)=0
       do jm=1,dgenea%nm
        spt%probg(c,jm)=1.d0
        spt%ngenom(c,jm+1)=spt%ngenom(c,jm)+1
        spt%ngend(c,jm)=dgenea%ndm(jm)
        spt%ngend(c,jm+1)=dgenea%ndm(jm+1)
        do kd=dgenea%ndm(jm)+1,dgenea%ndm(jm+1)
          spt%ndesc(c,kd)=kd
        end do
        do ll=1,dataset%map%nmk(c)
          if ( phm(jm,ll) == 1 ) then
            ptfin(c,ll,spt%ngend(c,jm)+1:spt%ngend(c,jm+1),:)=spt%prot(c,ll,spt%ngend(c,jm)+1:spt%ngend(c,jm+1),:)
          else
            ptfin(c,ll,spt%ngend(c,jm)+1:spt%ngend(c,jm+1),1)=spt%prot(c,ll,spt%ngend(c,jm)+1:spt%ngend(c,jm+1),2)
            ptfin(c,ll,spt%ngend(c,jm)+1:spt%ngend(c,jm+1),2)=spt%prot(c,ll,spt%ngend(c,jm)+1:spt%ngend(c,jm+1),1)
            ptfin(c,ll,spt%ngend(c,jm)+1:spt%ngend(c,jm+1),3)=spt%prot(c,ll,spt%ngend(c,jm)+1:spt%ngend(c,jm+1),4)
            ptfin(c,ll,spt%ngend(c,jm)+1:spt%ngend(c,jm+1),4)=spt%prot(c,ll,spt%ngend(c,jm)+1:spt%ngend(c,jm+1),3)
          end if
        end do !!
       end do !jm
      end do !c

      deallocate (trans,H,nbp,nbm,W,ndf,corref,repdes,phm,r)

     end subroutine calcul_phases_symmax2sat_dam
!!***

end module m_qtlmap_haplotype_V2
