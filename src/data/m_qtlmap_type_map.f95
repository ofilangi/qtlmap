module m_qtlmap_type_map
    use m_qtlmap_constant
    use m_qtlmap_log
    use m_qtlmap_base
    implicit none


!!****t* m_qtlmap_type_map/MAP_BASE
!!  NAME
!!     MAP_BASE
!!  DESCRIPTION
!!
!!  NOTES
!!
!! SOURCE
 type  :: MAP_BASE
   real                                 ,      private     :: BASE_STEP  = 1.d2
   ! Step analysis
   integer                              ,      public      :: PAS   = 0

   ! Number of chromosome to analyse
   integer                              ,      public      :: nchr  = 0
   ! List of chromosome name
   character(len=LEN_S) ,dimension(:) , pointer, public    :: chromo => null()
   ! Number of the marker by chromosome
   integer ,  dimension (:),   pointer, public             :: nmk   => null()
   ! Name of the markers by chromosome
   character(len=LEN_DEF), dimension (:,:), pointer,public :: mark  => null()
   ! Average position in Morgan by chromosome and marker
   real (kind=dp), dimension (:,:), pointer,public         :: posi  => null()
   ! Mal position in Morgan by chromosome and marker
   real (kind=dp), dimension (:,:), pointer,public         :: posim => null()
   ! Female position in Morgan by chromosome and marker
   real (kind=dp), dimension (:,:), pointer,public         :: posif => null()
   ! Xaldane distance computed with posim position
   real (kind=dp), dimension (:,:,:), pointer,public       :: rm    => null()
   ! Xaldane distance computed with posif position
   real (kind=dp), dimension (:,:,:), pointer,public       :: rf    => null()
   ! get abscis in Morgan corresponding to a a position N on the chromosome CHR
   real (kind=dp)   ,dimension(:,:),pointer   ,public      :: absi  => null()

   contains
     procedure, public :: copy      => copy_map_base
     procedure, public :: link      => link_map_base
     procedure, public :: release   => release_map_base
     procedure, public :: add_marker
     procedure, public :: get_maxnpo
     procedure, public :: set_base_and_step
     procedure, public :: get_long_step_morgan
     procedure, public :: set_absi
     procedure, public :: get_ilong
     procedure, public :: get_dx
     procedure, public :: get_pos
     procedure, public :: get_npo
     procedure, public :: get_flanking_marker
     procedure, public :: get_near_marker ! get the marker near the position
     procedure, public :: print       => print_map_base
     procedure, public :: getChridPosId

 end type MAP_BASE
!!***

    contains
      ! add a marker in the structure
      subroutine add_marker(map,chromosome,marker,posi,posim,posif)
        class(MAP_BASE)        , intent(inout)       :: map
        character(len=*)       , intent(in)          :: chromosome
        character(len=*)       , intent(in)          :: marker
        real(kind=dp)          , intent(in)          :: posi,posim,posif

        integer, parameter :: BUFFER_SIZE=10
        type(MAP_BASE)  :: map_save
        integer :: ic,l1,l2,buffer_size_tot,i
        logical :: found
        real(kind=dp) :: dm,df

        if ( posi<=0.0 .or. posim<=0.0 .or. posif<=0.0 ) then
          call log_mess("can not add marker in the map with the following position :"//&
          trim(str(posi))//" "//trim(str(posif))//" "//trim(str(posif)))
          return
        end if

        !first call
        if ( .not. associated(map%chromo) ) then
         allocate (map%chromo(1))
         map%chromo = chromosome
         map%nchr = 1
         allocate (map%mark(1,BUFFER_SIZE))
         map%mark(1,1)=marker
         allocate (map%nmk(1))
         map%nmk(1) = 1
         allocate (map%posi(1,BUFFER_SIZE),map%posim(1,BUFFER_SIZE),map%posif(1,BUFFER_SIZE))
         map%posi(1,1) = posi
         map%posim(1,1) = posim
         map%posif(1,1) = posif
         allocate(map%rm(1,BUFFER_SIZE,BUFFER_SIZE),map%rf(1,BUFFER_SIZE,BUFFER_SIZE))
         map%rm = 0.d0
         map%rf = 0.d0
         return
        end if

        found=.false.
        do ic=1,map%nchr
         if ( map%chromo(ic) == chromosome) then
          found=.true.
          exit
         end if
        end do

        ! new chromosome
        if ( .not. found ) then
          call map%copy(map_save)
          call map%release()
          buffer_size_tot = size(map_save%mark,2)
          map%nchr = map_save%nchr+1
          allocate(map%chromo(map%nchr))
          map%chromo(:map_save%nchr) = map_save%chromo
          map%chromo(map%nchr) = chromosome
          allocate (map%mark(map%nchr,buffer_size_tot))
          map%mark(:map_save%nchr,:) = map_save%mark
          map%mark(map%nchr,1)=marker
          allocate (map%nmk(map%nchr))
          map%nmk(:map_save%nchr)=map_save%nmk
          map%nmk(map%nchr) = 1
          allocate (map%posi(map%nchr,buffer_size_tot),&
                    map%posim(map%nchr,buffer_size_tot),&
                    map%posif(map%nchr,buffer_size_tot))
          map%posi(:map_save%nchr,:) = map_save%posi
          map%posim(:map_save%nchr,:) = map_save%posim
          map%posif(:map_save%nchr,:) = map_save%posif
          map%posi(map%nchr,1) = posi
          map%posim(map%nchr,1) = posim
          map%posif(map%nchr,1) = posif
          allocate(map%rm(map%nchr,buffer_size_tot,buffer_size_tot),&
                   map%rf(map%nchr,buffer_size_tot,buffer_size_tot))
          map%rm(:map_save%nchr,:,:) = map_save%rm
          map%rf(:map_save%nchr,:,:) = map_save%rf
          map%rm(map%nchr,:,:) = 0.d0
          map%rf(map%nchr,:,:) = 0.d0
          call map_save%release()
          return
        end if

        !check name marker
        do i=1,map%nmk(ic)
         if ( map%mark(ic,i) == marker ) then
           call log_mess("cannot add twice the marker ["//trim(marker)//"] in the map!")
           return
         end if
        end do


        !if new marker with not enough allocation
        if ( map%nmk(ic)>=BUFFER_SIZE) then
          call map%copy(map_save)
          call map%release();
          buffer_size_tot = size(map_save%mark,2)+BUFFER_SIZE

          map%nchr = map_save%nchr
          allocate (map%chromo(map%nchr))
          map%chromo = map_save%chromo
          allocate (map%nmk(map%nchr))
          map%nmk = map_save%nmk

          allocate (map%mark(map%nchr,buffer_size_tot))
          map%mark(:,:size(map_save%mark,2)) = map_save%mark

          allocate (map%posi(map%nchr,buffer_size_tot),&
                    map%posim(map%nchr,buffer_size_tot),&
                    map%posif(map%nchr,buffer_size_tot))
          map%posi(:,:size(map_save%mark,2)) = map_save%posi
          map%posim(:,:size(map_save%mark,2)) = map_save%posim
          map%posif(:,:size(map_save%mark,2)) = map_save%posif
          allocate(map%rm(map%nchr,buffer_size_tot,buffer_size_tot),&
                   map%rf(map%nchr,buffer_size_tot,buffer_size_tot))
          map%rm(:,:size(map_save%mark,2),:size(map_save%mark,2)) = map_save%rm
          map%rf(:,:size(map_save%mark,2),:size(map_save%mark,2)) = map_save%rf
          call map_save%release()
        end if

        map%nmk(ic) = map%nmk(ic) + 1! add the element
        map%mark(ic,map%nmk(ic))  = marker
        map%posi(ic,map%nmk(ic))  = posi
        map%posim(ic,map%nmk(ic)) = posim
        map%posif(ic,map%nmk(ic)) = posif

        l1 = map%nmk(ic)-1
        l2 = map%nmk(ic)
        dm=map%posim(ic,l2)-map%posim(ic,l1)
        df=map%posim(ic,l2)-map%posim(ic,l1)
        map%rm(ic,l1,l2)=xaldane(dm)
        map%rf(ic,l1,l2)=xaldane(df)

      end subroutine add_marker

      ! Transforme N points (Chr(N),Pos(N)) en un identidiant ChrId, PosId.
    ! Pour eviter de stoquer toutes les combinaisons de points (ex N,N1 et N1,N en BI-qtl)
    ! On reordonne les indice chr et pos du plus petit au plus grand

    subroutine getChridPosId(map,n,chr,pos,chrid,posid)
        class(MAP_BASE)       ,intent(in)    :: map
        integer               ,intent(in)    :: n
        integer ,dimension(n) ,intent(in)    :: chr
        integer ,dimension(n) ,intent(in)    :: pos
        integer               ,intent(inout) :: chrid
        integer               ,intent(inout) :: posid

        integer :: i,dimid,j,buf
        logical :: ok
        integer ,dimension(n)    :: mychr
        integer ,dimension(n)    :: mypos

        mychr=chr
        mypos = pos
        ok=.false.
        !On reordonne

        !Par chromosome
        do while (.not. ok)
            ok=.true.
            do i=1,n
                do j=i+1,n
                    if ( mychr(i)>mychr(j) ) then
                        ok=.false.
                        exit
                    end if
                end do
                if ( .not. ok ) exit
            end do
            if ( .not. ok ) then
                buf = mychr(i)
                mychr(i)=mychr(j)
                mychr(j)=buf
                buf = mypos(i)
                mypos(i)=mypos(j)
                mypos(j)=buf
            end if
        end do

        ok=.false.
        !Par Position
        do while (.not. ok)
            ok=.true.
            do i=1,n
                do j=i+1,n
                    if ( (mychr(i)==mychr(j)) .and. mypos(i)>mypos(j) ) then
                        ok=.false.
                        exit
                    end if
                end do
                if ( .not. ok ) exit
            end do
            if ( .not. ok ) then
                buf = mypos(i)
                mypos(i)=mypos(j)
                mypos(j)=buf
            end if
        end do

        ! Calcul de l'index Chr
        chrid=0
        dimid=1
        do i=1,n
            chrid = chrid + (mychr(i)-1)*dimid
            dimid = dimid*map%nchr
        end do
        chrid = chrid + 1

        ! Calcul de l'index Pos
        posid=0
        dimid=1
        do i=1,n
            posid = posid + (mypos(i)-1)*dimid
            dimid = dimid*map%get_npo(mychr(i))
        end do
        posid = posid + 1

    end subroutine getChridPosId



      subroutine release_map_base(map)
        class(MAP_BASE), intent(inout)  :: map

        if (associated(map%absi)) deallocate (map%chromo)
        if (associated(map%nmk)) deallocate (map%nmk)
        if (associated(map%mark)) deallocate (map%mark)
        if (associated(map%posi)) deallocate (map%posi)
        if (associated(map%posif)) deallocate (map%posif)
        if (associated(map%posim)) deallocate (map%posim)
        if (associated(map%rm)) deallocate (map%rm)
        if (associated(map%rf)) deallocate (map%rf)
        if (associated(map%absi)) deallocate (map%absi)

      end subroutine release_map_base




     subroutine link_map_base(map,copyMap)
       class(MAP_BASE), intent(in)  :: map
       type(MAP_BASE), intent(inout)  :: copyMap


        copyMap%nmk => map%nmk
        copyMap%chromo => map%chromo
        copyMap%mark => map%mark
        copyMap%posi => map%posi
        copyMap%posim => map%posim
        copyMap%posif => map%posif
        copyMap%rm => map%rm
        copyMap%rf => map%rf
        copyMap%absi => map%absi

       copyMap%BASE_STEP  = map%BASE_STEP
       copyMap%PAS        = map%PAS
       copyMap%nchr       = map%nchr

    end subroutine link_map_base


    subroutine copy_map_base(map,copyMap)
       class(MAP_BASE), intent(in)  :: map
       type(MAP_BASE), intent(inout)  :: copyMap

       if ( associated(map%nmk)) then
        allocate(copyMap%nmk(size(map%nmk)))
        copyMap%nmk = map%nmk
       end if

       if ( associated(map%chromo)) then
        allocate(copyMap%chromo(size(map%chromo)))
        copyMap%chromo = map%chromo
       end if


       if ( associated(map%mark)) then
        allocate(copyMap%mark(size(map%mark,1),size(map%mark,2)))
        copyMap%mark = map%mark
       end if

       if ( associated(map%posi)) then
        allocate(copyMap%posi(size(map%posi,1),size(map%posi,2)))
        copyMap%posi = map%posi
       end if

       if ( associated(map%posim)) then
        allocate(copyMap%posim(size(map%posim,1),size(map%posim,2)))
        copyMap%posim = map%posim
       end if

       if ( associated(map%posif)) then
        allocate(copyMap%posif(size(map%posif,1),size(map%posif,2)))
        copyMap%posif = map%posif
       end if

       if ( associated(map%rm)) then
        allocate(copyMap%rm(size(map%rm,1),size(map%rm,2),size(map%rm,3)))
        copyMap%rm = map%rm
       end if

       if ( associated(map%rf)) then
        allocate(copyMap%rf(size(map%rf,1),size(map%rf,2),size(map%rf,3)))
        copyMap%rf = map%rf
       end if

       if ( associated(map%absi)) then
        allocate(copyMap%absi(size(map%absi,1),size(map%absi,2)))
        copyMap%absi = map%absi
       end if

       copyMap%BASE_STEP  = map%BASE_STEP
       copyMap%PAS        = map%PAS
       copyMap%nchr       = map%nchr

    end subroutine copy_map_base


!!****f* m_qtlmap_tools/set_base_and_step
!! NAME
!!    set_base_and_step
!! DESCRIPTION
!!    Set the step (variable PAS) and the base dimension (variable BASE_STEP) according a Morgan step.
!!
!! INPUTS
!!     stepMorgan : the step in Morgan
!! SOURCE
!!
        subroutine set_base_and_step(map,stepMorgan)
         class(MAP_BASE)   ,intent(inout) :: map
         character(len=*)  ,intent(in)    :: stepMorgan
         integer :: i,s
         character(len=LEN_W) :: buf

         buf=""
         buf=adjustl(stepMorgan)
         s=len(trim(buf))
         !on cherche le point.....
         do i=1,s
           if (buf(i:i) == '.' ) then
              exit
           end if
         end do
         ! pas de point...
         if ( i > s ) then
            map%BASE_STEP = 1
            map%PAS = get_int(buf)
         else
            if ( (s-i) >= 10 ) then
              call stop_application("** Bad definition of step **")
            end if
            map%BASE_STEP = 10**(s-i)
            buf=trim(buf(:i-1))//trim(buf(i+1:))
            do while ( buf(1:1)=='0' )
               buf=trim(buf(2:))
            end do
            map%PAS = get_int(buf)
         end if

         call log_mess("BASE STEP : "//trim(str(map%BASE_STEP)),DEBUG_DEF)
         call log_mess("     STEP : "//trim(str(map%PAS)),DEBUG_DEF)

       end subroutine set_base_and_step
!!***

!!****f* m_qtlmap_tools/get_long_step_morgan
!!  NAME
!!    get_long_step_morgan
!!  DESCRIPTION
!!    Get the the step in Morgan
!!
!! SOURCE
       function get_long_step_morgan(map) result(l)
         class(MAP_BASE), intent(in):: map
         real(kind=dp) :: l

         l = real(map%pas) / map%BASE_STEP

       end function get_long_step_morgan
!!***

!!****f* m_qtlmap_tools/set_absi
!! NAME
!!    set_absi
!! DESCRIPTION
!!    build the global array absi. absi is a corresponding vector with a sampled positon.
!!    absi(ch,iposi)=<position in morgan of iposi>, iposi is an integer  1<=iposi<=get_ilong(ch)
!! HISTORY
!!  07/09/2010 - add correspond vector with marker position
!!  04/10/2011 - add support position on each marker
!! SOURCE
       subroutine set_absi(map)
         class(MAP_BASE), intent(inout):: map
         integer :: i,chr,l,n,ix,ilk,nlong,k
         real :: v,lastv,lastvm,lastvf
         logical :: vOk,vmOk,vfOk

         integer ,allocatable , dimension(:,:) :: corr_average_absi,corr_mal_absi,corr_femal_absi

         if ( map%nchr == 0 ) then
           call stop_application("Devel error: Can not set absci before the map structure!")
         end if

         allocate (map%absi(map%nchr,map%get_maxnpo()))

         if ( map%pas == 0 ) then
           do chr=1,map%nchr
             do k=1,map%nmk(chr)
               map%absi(chr,k)=map%posi(chr,k)
             end do
           end do

         else
           allocate (corr_average_absi(map%nchr,maxval(map%nmk)))
           allocate (corr_mal_absi(map%nchr,maxval(map%nmk)))
           allocate (corr_femal_absi(map%nchr,maxval(map%nmk)))
           corr_average_absi=0
           corr_mal_absi=0
           corr_femal_absi=0

          do chr=1,map%nchr
           l = map%get_ilong(chr)
           n=0
           do ix=0,l,map%pas
              n=n+1
              map%absi(chr,n)=map%get_dx(chr,ix)
           end do

          nlong=n
          do ilk=1,map%nmk(chr)
            lastv = 10000.d0;lastvm = 10000.d0;lastvf = 10000.d0
            vOk = .false.
            vmOk = .false.
            vfOk = .false.
            do n=1,nlong
             !! average map

               v = ( map%posi(chr,ilk) - map%absi(chr,n) ) / map%posi(chr,ilk)

               if (  v <= 0.0 .or. n == nlong ) then

                 if ( v < 0.d0 ) then ! maybe the last value was nearest....
                   if ( abs(v) < abs(lastv) )  then
                      corr_average_absi(chr,ilk)=n
                   else
                      corr_average_absi(chr,ilk)=n-1
                   end if
                 else
                     corr_average_absi(chr,ilk)=n
                 end if
               vOk = .true.
              end if
              lastv = v
              !! mal map

               v = ( map%posim(chr,ilk) - map%absi(chr,n) ) / map%posim(chr,ilk)

               if (  v <= 0.0 .or. n == nlong ) then

                 if ( v < 0.d0 ) then ! maybe the last value was nearest....
                   if ( abs(v) < abs(lastvm) )  then
                      corr_mal_absi(chr,ilk)=n
                   else
                      corr_mal_absi(chr,ilk)=n-1
                   end if
                 else
                     corr_mal_absi(chr,ilk)=n
                 end if
               vmOk = .true.
              end if
              lastvm = v

              !! femal map

               v = ( map%posif(chr,ilk) - map%absi(chr,n) ) / map%posif(chr,ilk)

               if (  v <= 0.0 .or. n == nlong ) then

                 if ( v < 0.d0 ) then ! maybe the last value was nearest....
                   if ( abs(v) < abs(lastvf) )  then
                      corr_femal_absi(chr,ilk)=n
                   else
                      corr_femal_absi(chr,ilk)=n-1
                   end if
                 else
                     corr_femal_absi(chr,ilk)=n
                 end if
               vfOk = .true.
              end if
              lastvf = v


            if ( vOk .and. vmOk .and. vfOk ) exit

            end do
!            call log_mess("CH "//str(chr)//" m"//str(ilk)//":"//&
!             str(corr_average_absi(chr,ilk))//" "//&
!             str(corr_mal_absi(chr,ilk))//" "//&
!             str(corr_femal_absi(chr,ilk)),DEBUG_DEF);
           end do
          end do

         deallocate (corr_average_absi)
         deallocate (corr_mal_absi)
         deallocate (corr_femal_absi)
   !      stop
        end if

       end subroutine set_absi
!!***

!!****f* m_qtlmap_tools/get_ilong
!!  NAME
!!    get_ilong
!!  DESCRIPTION
!!    Get the number of sampled position among the selected chromosome
!!  INPUTS
!!    index_chr  : index of chromosome
!!  RESULT
!!    The chromosome size in number of position sampled
!!  NOTES
!!  SEE ALSO
!!  m_qtlmap_data/nmk
!!  m_qtlmap_data/posi
!!
!! SOURCE
       function get_ilong(map,index_chr) result(ilong)
           class(MAP_BASE), intent(in) :: map
           integer, intent(in)         :: index_chr
           integer                     :: ilong

           if ( index_chr > map%nchr ) then
              call stop_application("dev error: get_ilong can not access to index chr::"&
             //trim(str(index_chr)))
           end if
           ilong=nint(map%BASE_STEP*(map%posi(index_chr,map%nmk(index_chr))-map%posi(index_chr,1)))
           return

       end function get_ilong
!!***

!!****f* m_qtlmap_tools/get_dx
!!  NAME
!!    get_dx
!!  DESCRIPTION
!!    Get the position in Morgan
!!
!!  INPUTS
!!    index_chr   : index of chromosome
!!    ix          : the position index
!!
!!  RETURN
!!   the position in Morgan
!!  NOTES
!!
!!  BUGS
!!
!! SOURCE

        function get_dx(map,index_chr,ix) result(dx)
           class(MAP_BASE), intent(in) :: map
           integer, intent(in)         :: index_chr,ix
           real(kind=dp)               :: dx

           if ( index_chr > map%nchr ) then
             call stop_application("dev error: get_ilong can not access to index chr::"&
             //trim(str(index_chr)))
           end if

           dx=map%posi(index_chr,1)+(dble(ix)/map%BASE_STEP)
          ! print *,"ix:",ix,"dx:",dx
           return

        end function get_dx
!!***

!!****f* m_qtlmap_tools/get_pos
!!  NAME
!!    get_pos
!!  DESCRIPTION
!!    Get the index of the sampled position
!!
!!  INPUTS
!!    posInMorgan : position in Morgan
!!
!!  RETURN
!!    the index of sampled position
!!  NOTES
!!   04/10/2011
!!     add support position on each marker
!!
!!  BUGS
!!
!! SOURCE
        function get_pos(map,chr,posInMorgan) result(npo)
           class(MAP_BASE), intent(inout)    :: map
           integer ,intent(in)               :: chr
           real(kind=dp), intent(in)         :: posInMorgan
           integer                           :: npo,i

           real(kind=dp) :: minv=9999.9

           if ( map%pas == 0 ) then
            do i=1,map%nmk(chr)
              if (minv>abs(map%posi(chr,i) - posInMorgan)) then
                minv = abs(map%posi(chr,i) - posInMorgan)
                npo = i
              end if
            end do
           else
            npo=int(map%BASE_STEP*posInMorgan)
            npo=nint(real(npo)/real(map%PAS))
           end if

           return

        end function get_pos
!!***

!!****f* m_qtlmap_tools/get_npo
!!  NAME
!!    get_npo
!!  DESCRIPTION
!!    Getting the number of postion depending the chromosome , step analysis and ilong (size of the chromosome)
!!  INPUTS
!!    * index_chr     -- Chromosome number
!!  RESULT
!!    The number of position to computing Probabilities of transmission and likelihood
!!  EXAMPLE
!!
!!  NOTES
!!   04/10/2011
!!     add support position on each marker
!!
!!  BUGS
!!
!!  SEE ALSO
!!  DATA/m_qtlmap_data/get_ilong
!! SOURCE
       function get_npo(map,index_chr) result(npo)
           class(MAP_BASE), intent(in) :: map
           integer, intent(in)         :: index_chr
           integer                     :: npo

           if ( map%nchr == 0 ) then
              call stop_application("Devel error: none map defined !");
           end if

           if ( index_chr > size(map%nmk) ) then
              call stop_application("Devel error: index_chr too big !");
           end if

           if ( map%pas == 0 ) then
            npo = map%nmk(index_chr)
           else
            npo =  (map%get_ilong(index_chr) / map%pas)+1
           end if
           return

       end function get_npo

!!    Getting the maximum npo
!!  RESULT
!!    maximum number of position on a chromosome
!!  NOTES
!!   04/10/2011
!!     add support position on each marker
       function get_maxnpo(map) result(maxnpo)
           class(MAP_BASE), intent(in)   :: map
           integer                       :: chr,npo,maxnpo

           if ( map%pas == 0 ) then
            maxnpo=sum(map%nmk)
           else
            maxnpo=0

            do chr=1,map%nchr
             npo = map%get_npo(chr)
             maxnpo = max(maxnpo,npo)
            end do
           end if

           return

       end function get_maxnpo

!!   get the flanking marker from a position (index sampling) according the average map
!!  INPUTS/OUTPUTS
!!    the left and right flanking marker
   subroutine  get_flanking_marker(map,chr,nx,flleft,flright)
      class(MAP_BASE), intent(in):: map
      integer , intent(in)  :: chr,nx
      integer , intent(out) :: flleft,flright

      if ( nx < 0 .or. map%get_npo(chr)<nx ) then
        call log_mess("m_qtlmap_tools/get_flanking_marker :: can not give the flanking marker, bad position :"&
        //trim(str(nx))//" chr:"//trim(str(chr)),WARNING_DEF);
        return
      end if

      flleft=1

      do while (map%absi(chr,nx) > map%posi(chr,flleft))
         flleft=flleft+1
         if ( size(map%posi,2)<flleft) exit
      end do

      if (flleft>1) flleft = flleft - 1

      if ( map%absi(chr,nx) /= map%posi(chr,flleft) .and. flleft < map%nmk(chr)) then
       flright = flleft+1
      else
       flright=flleft
      end if

   end subroutine  get_flanking_marker


   function get_near_marker(map,chr,nx) result(flk)
    class(MAP_BASE), intent(in):: map
      integer , intent(in)  :: chr,nx

      integer  :: flleft,flright,flk

      call map%get_flanking_marker(chr,nx,flleft,flright)
      if (abs(map%posi(chr,flleft)-map%absi(chr,nx)) > &
       abs(map%posi(chr,flright)-map%absi(chr,nx))) then
       flk = flright
      else
       flk = flleft
      end if

   end function  get_near_marker


   subroutine print_map_base(map,unit)
        class(MAP_BASE) ,intent(in)  :: map
        integer         ,intent(in)  :: unit
        integer :: c,ik

        do c=1,map%nchr
         do ik=1,map%nmk(c)
           write (unit,*) trim(map%chromo(c))//' '// &
                          trim(map%mark(c,ik))//' '// &
                          trim(str(map%posi(c,ik)))//' '// &
                          trim(str(map%posim(c,ik)))//' '// &
                          trim(str(map%posif(c,ik)))
         end do
        end do

   end subroutine print_map_base

end module m_qtlmap_type_map
