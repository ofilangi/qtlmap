!!****m* HAPLOTYPE/m_qtlmap_phase_offspring
!!  NAME
!!    m_qtlmap_phase_offspring
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

module m_qtlmap_phase_offspring
    use m_qtlmap_types
    use m_qtlmap_base
    use m_qtlmap_types
    use m_qtlmap_log
    implicit none

    public :: haplotype_offspring_v1

    contains

!!****f* m_qtlmap_phase_offspring/haplotype_offspring_v1
!! NAME
!!    haplotype_offspring_v1
!! DESCRIPTION
!!  Prediction of haplotypes in offspring:
!!  assuming a known sire phase 
!!  predicting the phase of offspring using flanking markers for marker with a unknown phase
!!  (inspired from the linkphase software)
!! NOTES
!!   creation : 14/12/2010 : carole moreno
!! 
!! SOURCE
subroutine haplotype_offspring_v1(dataset,spt)
 type(QTLMAP_DATASET)       ,intent(in)            :: dataset
 type(PDD_BUILD)            ,intent(inout)         :: spt

!A paramétricer par OLIVIER --> message de CMO

!A paramétricer par OLIVIER --> message de CMO
integer             :: mktot1, mktot2
real   (kind=dp)     :: ph1,ph2, ph
 integer :: i, j, kkd, lk, mk1, mk2, ori1, ori2, k1, k2,c, ip,jm, k,kd,nd1, nd2, nm1, nm2
  logical  , dimension (:,:),   allocatable        :: haplotyped, genotyped
 integer  :: stat
 type(GENEALOGY_BASE) , pointer :: dg
 type(GENOTYPE_BASE) , pointer :: dga

 dga => dataset%genoAnimal
 dg => dataset%genea


 allocate (haplotyped(size(dga%numero),maxval(dataset%map%nmk)),STAT=stat)
 call check_allocate(stat,'haplotyped')
 allocate (genotyped(size(dga%numero),maxval(dataset%map%nmk)),STAT=stat)
 call check_allocate(stat,'genotyped')
 
 allocate (spt%reconstructed(dataset%map%nchr,size(dga%numero),maxval(dataset%map%nmk)),STAT=stat)
 call check_allocate(stat,'reconstructed')

 DO c=1,dataset%map%nchr
!A paramétricer par OLIVIER --> message de CMO
  mktot1=1
  mktot2=dataset%map%nmk(c)
!A paramétricer par OLIVIER --> message de CMO
  haplotyped=.false.
  genotyped=.false.
  spt%reconstructed=.false.
  do ip=1,dg%np
         !print *, ip, ' pere', ' ', dga%correp(ip)
    DO k=mktot1,mktot2 !boucle sur les marqueurs
! HYPOTHESE: on suppose que genotyp est phasé pour les pères
! Definition des variables:	* genotyped pour les pères 
!				* haplotyped pour les pères
      if (dga%pheno(c,k,dga%correp(ip),1)/=dga%nmanque.and.dga%pheno(c,k,dga%correp(ip),2)/=dga%nmanque) &
         genotyped(dga%correp(ip),k) =.true.
      if (spt%genotyp(c,k,dga%correp(ip),1)/=dga%nmanque.and.spt%genotyp(c,k,dga%correp(ip),2)/=dga%nmanque) &
         haplotyped(dga%correp(ip),k)=.true.
    ENDDO
    nm1=dg%nmp(ip)+1
    nm2=dg%nmp(ip+1)
    do jm=nm1,nm2
      nd1=dg%ndm(jm)+1
      nd2=dg%ndm(jm+1)
      do kd=nd1,nd2
        ph=0.d0;ph1=0.d0;ph2=0.d0
        if (dga%corred(kd) == 9999 ) cycle
        spt%genotyp(c,:,dga%corred(kd),:)=dga%nmanque
        IF(count(dga%presentg(:,kd))>0) then !condition pour sélectionner les individus génotypés

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!1) PHASAGE DES PERES ET DES DESCENDANTS LORSQUE LE DESCENDANT ET LE PERE SONT GENOTYPES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO k=mktot1,mktot2 !boucle sur les marqueurs
! Definition de la variable genotyped pour les descendants
   if (dga%pheno(c,k,dga%corred(kd),1)/=dga%nmanque.and.dga%pheno(c,k,dga%corred(kd),2)/=dga%nmanque)&
      genotyped(dga%corred(kd),k)=.true.
   !print *, k, animal(kd), dga%pheno(c,k,i,1),dga%pheno(c,k,i,2), pere(ip), dga%pheno(c,k,dga%correp(ip),1), dga%pheno(c,k,dga%correp(ip),2)

   IF (genotyped(dga%corred(kd),k).and.genotyped(dga%correp(ip),k)) then ! condition sur les couple père-des génotypés
 ! Le père est homozygote
      if ((dga%pheno(c,k,dga%correp(ip),1)==dga%pheno(c,k,dga%correp(ip),2)) &
     .and.(dga%pheno(c,k,dga%correp(ip),1)==dga%pheno(c,k,dga%corred(kd),1)))  then
          spt%genotyp(c,k,dga%corred(kd),1)=dga%pheno(c,k,dga%corred(kd),1)
          spt%genotyp(c,k,dga%corred(kd),2)=dga%pheno(c,k,dga%corred(kd),2)
          haplotyped(dga%corred(kd),k)=.true.
          spt%genotyp(c,k,dga%correp(ip),1)=dga%pheno(c,k,dga%correp(ip),1)
          spt%genotyp(c,k,dga%correp(ip),2)=dga%pheno(c,k,dga%correp(ip),2)
          haplotyped(dga%correp(ip),k)=.true.
! print *, animal(kd), k,'OK'
      endif
      if ((dga%pheno(c,k,dga%correp(ip),1)==dga%pheno(c,k,dga%correp(ip),2)) &
     .and.(dga%pheno(c,k,dga%correp(ip),1)==dga%pheno(c,k,dga%corred(kd),2)))  then
          spt%genotyp(c,k,dga%corred(kd),1)=dga%pheno(c,k,dga%corred(kd),2)
          spt%genotyp(c,k,dga%corred(kd),2)=dga%pheno(c,k,dga%corred(kd),1)
          haplotyped(dga%corred(kd),k)=.true.
          spt%genotyp(c,k,dga%correp(ip),1)=dga%pheno(c,k,dga%correp(ip),2)
          spt%genotyp(c,k,dga%correp(ip),2)=dga%pheno(c,k,dga%correp(ip),1)
          haplotyped(dga%correp(ip),k)=.true.
      endif
!LE père est héterozygote l'allèle de l'haplotype 1 ou 2 est transmis
      if     ( (spt%genotyp(c,k,dga%correp(ip),1)/=spt%genotyp(c,k,dga%correp(ip),2)        )  &
        .and.  (dga%pheno(c,k,dga%corred(kd),1)==dga%pheno(c,k,dga%corred(kd),2)) &
        .and.  (spt%genotyp(c,k,dga%correp(ip),1)==dga%pheno(c,k,dga%corred(kd),1    ))) then
          spt%genotyp(c,k,dga%corred(kd),1)=dga%pheno(c,k,dga%corred(kd),1)
          spt%genotyp(c,k,dga%corred(kd),2)=dga%pheno(c,k,dga%corred(kd),2)
          haplotyped(dga%corred(kd),k)=.true.
      endif
      if    (  (spt%genotyp(c,k,dga%correp(ip),1)/=spt%genotyp(c,k,dga%correp(ip),2) )   &
        .and.  (dga%pheno(c,k,dga%corred(kd),1)==dga%pheno(c,k,dga%corred(kd),2)) &
        .and.  (spt%genotyp(c,k,dga%correp(ip),2)==dga%pheno(c,k,dga%corred(kd),2  ) )) then
          spt%genotyp(c,k,dga%corred(kd),1)=dga%pheno(c,k,dga%corred(kd),2)
          spt%genotyp(c,k,dga%corred(kd),2)=dga%pheno(c,k,dga%corred(kd),1)
          haplotyped(dga%corred(kd),k)=.true.
      endif

   ENDIF ! FIN de la condition sur les couple père-des génotypés
  ENDDO ! FIn de la boucle sur les marqueurs

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 2) PHASAGE POUR LES MARQUEURS sans phase en fonction des phases des marqueurs voisins
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 DO k=mktot1,mktot2 ! boucle marker mktot1 mktot2
    IF(.not.haplotyped(dga%corred(kd),k))then
     mk1=0;mk2=0;ori1=0;ori2=0
     
! search for informative marker on the left side   
     if (k==mktot1) then
         mk1=0
     else
     do k1=k-1,mktot1,-1
!       if(dga%corred(kd)==5)  print *, dga%corred(kd), k, k1,haplotyped(dga%corred(kd),k1), haplotyped(dga%correp(ip),k1)
        if(haplotyped(dga%corred(kd),k1) .and. haplotyped(dga%correp(ip),k1)) then
         if(spt%genotyp(c,k1,dga%correp(ip),1)/=spt%genotyp(c,k1,dga%correp(ip),2))then ! search for informative marker
           if(spt%genotyp(c,k1,dga%corred(kd),1)==spt%genotyp(c,k1,dga%correp(ip),1))ori1=1
           if(spt%genotyp(c,k1,dga%corred(kd),1)==spt%genotyp(c,k1,dga%correp(ip),2))ori1=2
           mk1=k1;exit
         endif
        endif
       enddo
     endif

! search for informative marker on the right side   
     if (k==mktot2) then
         mk2=0
     else
       do k2=k+1,mktot2
        if(haplotyped(dga%corred(kd),k2) .and. haplotyped(dga%correp(ip),k2)) then
         if(spt%genotyp(c,k2,dga%correp(ip),1)/=spt%genotyp(c,k2,dga%correp(ip),2))then ! search for informative marker
           if(spt%genotyp(c,k2,dga%corred(kd),1)==spt%genotyp(c,k2,dga%correp(ip),1))ori2=1
           if(spt%genotyp(c,k2,dga%corred(kd),1)==spt%genotyp(c,k2,dga%correp(ip),2))ori2=2
           mk2=k2;exit
         endif
        endif
       enddo
     endif

!! reconstruction of phases for markers with unknown phase using flanking phased markers
     if(mk1/=0 .and. mk2/=0)then !mk1, mk2
        if(ori1==ori2)then ! no recombination
         ph=(1.00-(xaldane(dataset%map%posim(c,k)-&
         dataset%map%posim(c,k1))))*(1.00-xaldane(dataset%map%posim(c,k2)-dataset%map%posim(c,k)))/&
         (1.00-xaldane(dataset%map%posim(c,k2)-dataset%map%posim(c,k1)))
         if(ph>dataset%params%PROB_PHASE_DESC)then
           spt%reconstructed(c,dga%corred(kd),k)=.true.
           if(spt%genotyp(c,k,dga%correp(ip),ori1)==dga%pheno(c,k,dga%corred(kd),1)) then
             spt%genotyp(c,k,dga%corred(kd),1)=dga%pheno(c,k,dga%corred(kd),1)
             spt%genotyp(c,k,dga%corred(kd),2)=dga%pheno(c,k,dga%corred(kd),2)
           else
             spt%genotyp(c,k,dga%corred(kd),2)=dga%pheno(c,k,dga%corred(kd),1)
             spt%genotyp(c,k,dga%corred(kd),1)=dga%pheno(c,k,dga%corred(kd),2)
         endif
        endif
       endif

      else if(mk1/=0 .and. mk2==0)then
        ph1=1.00-xaldane(dataset%map%posim(c,k)-dataset%map%posim(c,k1))
        if(ph1>dataset%params%PROB_PHASE_DESC)then
          spt%reconstructed(c,dga%corred(kd),k)=.true.
         if(spt%genotyp(c,k,dga%correp(ip),ori1)==dga%pheno(c,k,dga%corred(kd),1))then
           spt%genotyp(c,k,dga%corred(kd),1)=dga%pheno(c,k,dga%corred(kd),1)
           spt%genotyp(c,k,dga%corred(kd),2)=dga%pheno(c,k,dga%corred(kd),2)
         else
           spt%genotyp(c,k,dga%corred(kd),2)=dga%pheno(c,k,dga%corred(kd),1)
           spt%genotyp(c,k,dga%corred(kd),1)=dga%pheno(c,k,dga%corred(kd),2)
         endif
        endif

       else if(mk1==0 .and. mk2/=0)then
         ph2=1.00-xaldane(dataset%map%posim(c,k2)-dataset%map%posim(c,k))
         if(ph2>dataset%params%PROB_PHASE_DESC)then
          spt%reconstructed(c,dga%corred(kd),k)=.true.
          if(spt%genotyp(c,k,dga%correp(ip),ori2)==dga%pheno(c,k,dga%corred(kd),1))then
            spt%genotyp(c,k,dga%corred(kd),1)=dga%pheno(c,k,dga%corred(kd),1)
            spt%genotyp(c,k,dga%corred(kd),2)=dga%pheno(c,k,dga%corred(kd),2)
          else
            spt%genotyp(c,k,dga%corred(kd),2)=dga%pheno(c,k,dga%corred(kd),1)
            spt%genotyp(c,k,dga%corred(kd),1)=dga%pheno(c,k,dga%corred(kd),2)
          endif
        endif

      endif !mk1, mk2

     ENDIF ! not.haplotyped

  END DO

ENDIF! fin de la condition pour sélectionner les individus génotypés
ENDDO ! end animal
ENDDO ! end mere
ENDDO ! end pere
ENDDO ! end chromosome

 deallocate(haplotyped)
 deallocate(genotyped)

end subroutine haplotype_offspring_v1
!!***


end module m_qtlmap_phase_offspring
