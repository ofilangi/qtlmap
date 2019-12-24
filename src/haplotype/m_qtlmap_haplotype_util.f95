module m_qtlmap_haplotype_util
  use m_qtlmap_types
  implicit none

  public :: transmi
  public :: gammapf

  contains

!******************************************************************************
!******************************************************************************
!
! Etablissement du tableau des probabilites d'un couple d evenements de
! transmission
!
! Sous programme applele par pdegm et pdet
!
! Version 1.0.
!
!******************************************************************************
!******************************************************************************
!
        subroutine transmi(r1,r2,tab)
!
        real (kind=dp)   ,intent(in)     :: r1
        real (kind=dp)   ,intent(in)     :: r2
        real (kind=dp)   ,intent(out),dimension(4,4)  :: tab

        real (kind=dp) :: r1r2,r1mr2,r1r2m,r1mr2m
        integer        :: i,j
!
        r1r2=r1*r2
        r1mr2=(1.d0-r1)*r2
        r1r2m=r1*(1.d0-r2)
        r1mr2m=(1.d0-r1)*(1.d0-r2)
!
        do i=1,4
          tab(i,i)=r1mr2m
        do j=i+1,4
            tab(j,j-1)=r1mr2
            tab(j-1,j)=r1mr2
          end do
          do j=i+2,4
            tab(j,j-2)=r1r2m
            tab(j-2,j)=r1r2m
          end do
          tab(i,5-i)=r1r2
        end do
      return
      end subroutine transmi


!******************************************************************************
!******************************************************************************
!
! Etablissement de l'origine grand parentale des alleles portes
! par les descendants
!
! Sous programme appele par qtl
!
! Version 2.0.
!
!******************************************************************************
!******************************************************************************
!
      subroutine gammapf(dataset,spt)
!
        type(QTLMAP_DATASET)       ,intent(in)            :: dataset
        type(PDD_BUILD)            ,intent(inout)         :: spt
!
! Divers
        integer                                  :: c ! chromosome
        logical typp,typm
        integer(kind=KIND_PHENO) :: m1,m2,mc1,mc2,md1,md2
        integer :: ll,ip,nm1,nm2,jm,nd1,nd2,kd,i
        type(GENEALOGY_BASE) , pointer :: dgenea
        type(GENOTYPE_BASE) , pointer :: dga

        dga => dataset%genoAnimal
        dgenea => dataset%genea

        spt%prot=0.d0
!
!******************************************************************************
! Etablissement de la liste des transmissions possibles marqueur par marqueur
!******************************************************************************
!
     do c=1,dataset%map%nchr
      do 100 ll=1,dataset%map%nmk(c)
        do 100 ip=1,dgenea%np
        typp=.false.
        if(dga%correp(ip).ne.9999)then
         if(dga%pheno(c,ll,dga%correp(ip),1).ne.dga%nmanque) then
            typp=.true.
            m1=dga%pheno(c,ll,dga%correp(ip),1)
          m2=dga%pheno(c,ll,dga%correp(ip),2)
          endif
	 end if
        nm1=dgenea%nmp(ip)+1
        nm2=dgenea%nmp(ip+1)
        do 100 jm=nm1,nm2
        typm=.false.
        if(dga%correm(jm).ne.9999)then
         if(dga%pheno(c,ll,dga%correm(jm),1).ne.dga%nmanque) then
            typm=.true.
            mc1=dga%pheno(c,ll,dga%correm(jm),1)
          mc2=dga%pheno(c,ll,dga%correm(jm),2)
          endif
	 end if
        nd1=dgenea%ndm(jm)+1
        nd2=dgenea%ndm(jm+1)
!
! Initialisation : toutes les transmissions sont possibles
        do kd=nd1,nd2
            do i=1,4
              spt%prot(c,ll,kd,i)=0.25d0
            end do
          end do

!
! Pere inconnu, mere heterozygote
        if((.not.typp).and.(typm.and.mc1.ne.mc2)) then
          do 200 kd=nd1,nd2
            if(dga%corred(kd).eq.9999) go to 200
            if((dga%pheno(c,ll,dga%corred(kd),1).eq.dga%nmanque))go to 200
            md1=dga%pheno(c,ll,dga%corred(kd),1)
            md2=dga%pheno(c,ll,dga%corred(kd),2)
            if(((md1.eq.mc1).and.(md2.eq.mc2)).or.                   &
              ((md1.eq.mc2).and.(md2.eq.mc1))) go to 200
              if((md1.eq.mc1).or.(md2.eq.mc1)) then
              spt%prot(c,ll,kd,1)=0.5d0
                spt%prot(c,ll,kd,3)=0.5d0
                spt%prot(c,ll,kd,2)=0.d0
              spt%prot(c,ll,kd,4)=0.d0
            else
                spt%prot(c,ll,kd,1)=0.d0
              spt%prot(c,ll,kd,3)=0.d0
              spt%prot(c,ll,kd,2)=0.5d0
              spt%prot(c,ll,kd,4)=0.5d0
              end if
  200       continue
          go to 100
          end if
!
! Pere homozygote, mere heterozygote
        if((typp.and.m1.eq.m2).and.(typm.and.mc1.ne.mc2)) then
! Pere AA et mere AB
            if(mc1.eq.m1) then
            do 101 kd=nd1,nd2
              if(dga%corred(kd).eq.9999) go to 101
              if((dga%pheno(c,ll,dga%corred(kd),1).eq.dga%nmanque))go to 101
              md1=dga%pheno(c,ll,dga%corred(kd),1)
              md2=dga%pheno(c,ll,dga%corred(kd),2)
              if(md1.eq.md2) then
                spt%prot(c,ll,kd,1)=0.5d0
                  spt%prot(c,ll,kd,3)=0.5d0
                  spt%prot(c,ll,kd,2)=0.d0
                spt%prot(c,ll,kd,4)=0.d0
              else
                  spt%prot(c,ll,kd,1)=0.d0
                spt%prot(c,ll,kd,3)=0.d0
                spt%prot(c,ll,kd,2)=0.5d0
                spt%prot(c,ll,kd,4)=0.5d0
                end if
  101         continue
              go to 100
            end if
! Pere AA et mere BA
            if(mc2.eq.m1) then
            do 102 kd=nd1,nd2
              if(dga%corred(kd).eq.9999) go to 102
              if((dga%pheno(c,ll,dga%corred(kd),1).eq.dga%nmanque))go to 102
              md1=dga%pheno(c,ll,dga%corred(kd),1)
              md2=dga%pheno(c,ll,dga%corred(kd),2)
              if(md1.eq.md2) then
                  spt%prot(c,ll,kd,1)=0.d0
                spt%prot(c,ll,kd,3)=0.d0
                spt%prot(c,ll,kd,2)=0.5d0
                spt%prot(c,ll,kd,4)=0.5d0
              else
                spt%prot(c,ll,kd,1)=0.5d0
                  spt%prot(c,ll,kd,3)=0.5d0
                spt%prot(c,ll,kd,2)=0.d0
                spt%prot(c,ll,kd,4)=0.d0
                end if
  102         continue
              go to 100
          end if
! Pere AA et mere BC
            if((mc1.ne.m1).and.(mc2.ne.m1)) then
            do 103 kd=nd1,nd2
              if(dga%corred(kd).eq.9999) go to 103
              if((dga%pheno(c,ll,dga%corred(kd),1).eq.dga%nmanque))go to 103
              md1=dga%pheno(c,ll,dga%corred(kd),1)
              md2=dga%pheno(c,ll,dga%corred(kd),2)
                if((md1.eq.mc1).or.(md2.eq.mc1)) then
                spt%prot(c,ll,kd,1)=0.5d0
                spt%prot(c,ll,kd,3)=0.5d0
                  spt%prot(c,ll,kd,2)=0.d0
                spt%prot(c,ll,kd,4)=0.d0
              else
                spt%prot(c,ll,kd,1)=0.d0
                  spt%prot(c,ll,kd,3)=0.d0
                spt%prot(c,ll,kd,2)=0.5d0
                spt%prot(c,ll,kd,4)=0.5d0
                end if
  103         continue
            go to 100
          end if
          end if
!
! Pere heterozygote, mere inconnue
        if((typp.and.m1.ne.m2).and.(.not.typm)) then
          do 201 kd=nd1,nd2
            if(dga%corred(kd).eq.9999)go to 201
            if((dga%pheno(c,ll,dga%corred(kd),1).eq.dga%nmanque))go to 201
            md1=dga%pheno(c,ll,dga%corred(kd),1)
            md2=dga%pheno(c,ll,dga%corred(kd),2)
            if(((md1.eq.m1).and.(md2.eq.m2)).or.           &
              ((md1.eq.m2).and.(md2.eq.m1))) go to 201
              if((md1.eq.m1).or.(md2.eq.m1)) then
              spt%prot(c,ll,kd,1)=0.5d0
                spt%prot(c,ll,kd,2)=0.5d0
                spt%prot(c,ll,kd,3)=0.d0
              spt%prot(c,ll,kd,4)=0.d0
            else
                spt%prot(c,ll,kd,1)=0.d0
              spt%prot(c,ll,kd,2)=0.d0
              spt%prot(c,ll,kd,3)=0.5d0
              spt%prot(c,ll,kd,4)=0.5d0
              end if
  201       continue
          go to 100
          end if
!
! Pere heterozygote, mere homozygote
        if((typp.and.m1.ne.m2).and.(typm.and.mc1.eq.mc2)) then
! Pere AB et mere AA
            if(mc1.eq.m1) then
            do 104 kd=nd1,nd2
              if(dga%corred(kd).eq.9999) go to 104
              if((dga%pheno(c,ll,dga%corred(kd),1).eq.dga%nmanque))go to 104
              md1=dga%pheno(c,ll,dga%corred(kd),1)
              md2=dga%pheno(c,ll,dga%corred(kd),2)
                if(md1.eq.md2) then
                spt%prot(c,ll,kd,1)=0.5d0
                spt%prot(c,ll,kd,2)=0.5d0
                  spt%prot(c,ll,kd,3)=0.d0
                spt%prot(c,ll,kd,4)=0.d0
              else
                spt%prot(c,ll,kd,1)=0.d0
                  spt%prot(c,ll,kd,2)=0.d0
                spt%prot(c,ll,kd,3)=0.5d0
                spt%prot(c,ll,kd,4)=0.5d0
                end if
  104         continue
              go to 100
          end if
! Pere AB et mere BB
            if(mc1.eq.m2) then
            do 105 kd=nd1,nd2
              if(dga%corred(kd).eq.9999) go to 105
              if((dga%pheno(c,ll,dga%corred(kd),1).eq.dga%nmanque))go to 105
              md1=dga%pheno(c,ll,dga%corred(kd),1)
              md2=dga%pheno(c,ll,dga%corred(kd),2)
                if(md1.eq.md2) then
                spt%prot(c,ll,kd,1)=0.d0
                spt%prot(c,ll,kd,2)=0.d0
                  spt%prot(c,ll,kd,3)=0.5d0
                spt%prot(c,ll,kd,4)=0.5d0
              else
                spt%prot(c,ll,kd,1)=0.5d0
                  spt%prot(c,ll,kd,2)=0.5d0
                spt%prot(c,ll,kd,3)=0.d0
                spt%prot(c,ll,kd,4)=0.d0
                end if
  105         continue
              go to 100
          end if
! Pere AB et mere CC
            if((mc1.ne.m1).and.(mc1.ne.m2)) then
            do 106 kd=nd1,nd2
              if(dga%corred(kd).eq.9999) go to 106
              if((dga%pheno(c,ll,dga%corred(kd),1).eq.dga%nmanque))go to 106
                 md1=dga%pheno(c,ll,dga%corred(kd),1)
                 md2=dga%pheno(c,ll,dga%corred(kd),2)
                if((md1.eq.m1).or.(md2.eq.m1)) then
                spt%prot(c,ll,kd,1)=0.5d0
                spt%prot(c,ll,kd,2)=0.5d0
                  spt%prot(c,ll,kd,3)=0.d0
                spt%prot(c,ll,kd,4)=0.d0
              else
                spt%prot(c,ll,kd,1)=0.d0
                  spt%prot(c,ll,kd,2)=0.d0
                spt%prot(c,ll,kd,3)=0.5d0
                spt%prot(c,ll,kd,4)=0.5d0
                end if
  106         continue
              go to 100
          end if
          end if
!
! Pere heterozygote, mere heterozygote
        if((typp.and.m1.ne.m2).and.(typm.and.mc1.ne.mc2)) then
! Pere AB et mere AB
            if((mc1.eq.m1).and.(mc2.eq.m2)) then
            do 107 kd=nd1,nd2
              if(dga%corred(kd).eq.9999) go to 107
              if((dga%pheno(c,ll,dga%corred(kd),1).eq.dga%nmanque))go to 107
              md1=dga%pheno(c,ll,dga%corred(kd),1)
              md2=dga%pheno(c,ll,dga%corred(kd),2)
                if(md1.eq.md2) then
                do i=1,4
                  spt%prot(c,ll,kd,i)=0.d0
                end do
                  if(md1.eq.m1) spt%prot(c,ll,kd,1)=1.d0
                if(md1.eq.m2) spt%prot(c,ll,kd,4)=1.d0
              else
                  spt%prot(c,ll,kd,1)=0.d0
                spt%prot(c,ll,kd,3)=0.5d0
                spt%prot(c,ll,kd,2)=0.5d0
                spt%prot(c,ll,kd,4)=0.d0
                end if
  107         continue
              go to 100
          end if
! Pere AB et mere BA
            if((mc1.eq.m2).and.(mc2.eq.m1)) then
            do 108 kd=nd1,nd2
              if(dga%corred(kd).eq.9999) go to 108
              if((dga%pheno(c,ll,dga%corred(kd),1).eq.dga%nmanque))go to 108
              md1=dga%pheno(c,ll,dga%corred(kd),1)
              md2=dga%pheno(c,ll,dga%corred(kd),2)
                if(md1.eq.md2) then
                do i=1,4
                  spt%prot(c,ll,kd,i)=0.d0
                end do
                  if(md1.eq.m1) spt%prot(c,ll,kd,2)=1.d0
                if(md1.eq.m2) spt%prot(c,ll,kd,3)=1.d0
              else
                  spt%prot(c,ll,kd,1)=0.5d0
                spt%prot(c,ll,kd,3)=0.d0
                spt%prot(c,ll,kd,2)=0.d0
                spt%prot(c,ll,kd,4)=0.5d0
                end if
  108         continue
              go to 100
          end if
! Pere AB et mere AC
            if(mc1.eq.m1) then
            do 109 kd=nd1,nd2
              if(dga%corred(kd).eq.9999) go to 109
              if((dga%pheno(c,ll,dga%corred(kd),1).eq.dga%nmanque))go to 109
              md1=dga%pheno(c,ll,dga%corred(kd),1)
              md2=dga%pheno(c,ll,dga%corred(kd),2)
              do i=1,4
                spt%prot(c,ll,kd,i)=0.d0
              end do
                if(md1.eq.md2) then
                  spt%prot(c,ll,kd,1)=1.d0
                else
                  if((md1.eq.m1).or.(md2.eq.m1)) then
                    if((md1.eq.m2).or.(md2.eq.m2)) spt%prot(c,ll,kd,3)=1.d0
                    if((md1.ne.m2).and.(md2.ne.m2)) spt%prot(c,ll,kd,2)=1.d0
                  else
                    spt%prot(c,ll,kd,4)=1.d0
                  end if
                end if
  109         continue
              go to 100
          end if
! Pere AB et mere CA
            if(mc2.eq.m1) then
            do 110 kd=nd1,nd2
              if(dga%corred(kd).eq.9999) go to 110
              if((dga%pheno(c,ll,dga%corred(kd),1).eq.dga%nmanque))go to 110
              md1=dga%pheno(c,ll,dga%corred(kd),1)
              md2=dga%pheno(c,ll,dga%corred(kd),2)
              do i=1,4
                spt%prot(c,ll,kd,i)=0.d0
              end do
                if(md1.eq.md2) then
                  spt%prot(c,ll,kd,2)=1.d0
                else
                  if((md1.eq.m1).or.(md2.eq.m1)) then
                    if((md1.eq.m2).or.(md2.eq.m2)) spt%prot(c,ll,kd,4)=1.d0
                    if((md1.ne.m2).and.(md2.ne.m2)) spt%prot(c,ll,kd,1)=1.d0
                  else
                    spt%prot(c,ll,kd,3)=1.d0
                  end if
                end if
  110         continue
              go to 100
          end if
! Pere AB et mere BC
            if(mc1.eq.m2) then
            do 111 kd=nd1,nd2
              if(dga%corred(kd).eq.9999) go to 111
              if((dga%pheno(c,ll,dga%corred(kd),1).eq.dga%nmanque))go to 111
              md1=dga%pheno(c,ll,dga%corred(kd),1)
              md2=dga%pheno(c,ll,dga%corred(kd),2)
              do i=1,4
                  spt%prot(c,ll,kd,i)=0.d0
              end do
                if(md1.eq.md2) then
                  spt%prot(c,ll,kd,3)=1.d0
                else
                  if((md1.eq.m1).or.(md2.eq.m1)) then
                    if((md1.eq.m2).or.(md2.eq.m2)) spt%prot(c,ll,kd,1)=1.d0
                    if((md1.ne.m2).and.(md2.ne.m2)) spt%prot(c,ll,kd,2)=1.d0
                  else
                    spt%prot(c,ll,kd,4)=1.d0
                  end if
                end if
  111         continue
              go to 100
          end if
! Pere AB et mere CB
            if(mc2.eq.m2) then
            do 112 kd=nd1,nd2
              if(dga%corred(kd).eq.9999) go to 112
              if((dga%pheno(c,ll,dga%corred(kd),1).eq.dga%nmanque))go to 112
              md1=dga%pheno(c,ll,dga%corred(kd),1)
              md2=dga%pheno(c,ll,dga%corred(kd),2)
              do i=1,4
                  spt%prot(c,ll,kd,i)=0.d0
              end do
                if(md1.eq.md2) then
                  spt%prot(c,ll,kd,4)=1.d0
                else
                  if((md1.eq.m1).or.(md2.eq.m1)) then
                    if((md1.eq.m2).or.(md2.eq.m2)) spt%prot(c,ll,kd,2)=1.d0
                    if((md1.ne.m2).and.(md2.ne.m2)) spt%prot(c,ll,kd,1)=1.d0
                  else
                    spt%prot(c,ll,kd,3)=1.d0
                  end if
                end if
  112         continue
              go to 100
          end if
! Pere AB et mere CD
          do 113 kd=nd1,nd2
            if(dga%corred(kd).eq.9999) go to 113
            if((dga%pheno(c,ll,dga%corred(kd),1).eq.dga%nmanque))go to 113
            md1=dga%pheno(c,ll,dga%corred(kd),1)
            md2=dga%pheno(c,ll,dga%corred(kd),2)
            do i=1,4
                spt%prot(c,ll,kd,i)=0.d0
            end do
              if(((md1.eq.m1).and.(md2.eq.mc1)).or.                   &
              ((md2.eq.m1).and.(md1.eq.mc1))) spt%prot(c,ll,kd,1)=1.d0
            if(((md1.eq.m2).and.(md2.eq.mc1)).or.                     &
              ((md2.eq.m2).and.(md1.eq.mc1))) spt%prot(c,ll,kd,3)=1.d0
            if(((md1.eq.m1).and.(md2.eq.mc2)).or.                     &
              ((md2.eq.m1).and.(md1.eq.mc2))) spt%prot(c,ll,kd,2)=1.d0
            if(((md1.eq.m2).and.(md2.eq.mc2)).or.                     &
              ((md2.eq.m2).and.(md1.eq.mc2))) spt%prot(c,ll,kd,4)=1.d0
  113       continue
            go to 100
        end if
  100 continue
     end do
      end subroutine gammapf


end module m_qtlmap_haplotype_util
