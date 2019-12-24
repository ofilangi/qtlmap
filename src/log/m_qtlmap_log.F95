
!!  SYNOPSIS
!!  This module manage the log console output.
!!  CREATION DATE
!!  01/01/2009

!! achar(27) est l'equivalent en C à "\x1B"

!! achar(27)//"[0m" // Normal
!achar(27)//"[1m" // Gras
!achar(27)//"[4m" // Souligné
!achar(27)//"[5m" // Clignotant
!achar(27)//"[7m" // Inverse vidéo
!
!/* VT100 : Couleurs */
!"\033[22;30m" // Noir
!"\033[22;31m" // Rouge
!"\033[22;32m" // Vert
!"\033[22;33m" // Brun
!"\033[22;34m" // Bleu
!"\033[22;35m" // Magenta
!"\033[22;36m" // Cyan
!"\033[22;37m" // Gris
!"\033[01;30m" // Gris foncé
!"\033[01;31m" // Rouge clair
!"\033[01;32m" // Vert clair
!"\033[01;33m" // Jaune
!"\033[01;34m" // Bleu clair
!"\033[01;35m" // Magenta clair
!"\033[01;36m" // Cyan clair
!"\033[01;37m" // Blanc

MODULE m_qtlmap_log

#ifdef MPI_ACTIVE
       use mpi
#endif

#ifdef MANAGE_OMP
       use omp_lib
#endif
    IMPLICIT NONE
    SAVE

    integer, public, parameter             :: DEBUG_DEF     = 0
    integer, public, parameter             :: VERBOSE_DEF   = 1
    integer, public, parameter             :: INFO_DEF      = 2
    integer, public, parameter             :: WARNING_DEF   = 3
    integer, public, parameter             :: ERROR_DEF     = 4

    character(len=4)  , private, parameter  :: default_col = achar(27)//"[0m"
    character(len=8)  , private, parameter  :: debug_col   = achar(27)//"[22;36m"
    character(len=8)  , private, parameter  :: verbose_col = achar(27)//"[22;32m"
    character(len=8)  , private, parameter  :: info_col    = achar(27)//"[22;32m"
    character(len=8)  , private, parameter  :: warning_col = achar(27)//"[22;35m"
    character(len=8)  , private, parameter  :: error_col   = achar(27)//"[31;01m"

    ! internal level to handle messages
    integer             ,private                   :: internal_level = INFO_DEF

    !identifiant du processus mpi courant
    integer                   ,private     :: thread_rang = 0

    public :: init_log
    public :: get_log_level
    public :: log_mess
    public :: log_vector
    public :: log_matrix
    public :: log_matrix_sparse
    public :: end_log


    interface log_matrix
        module procedure log_matrix_double
        module procedure log_matrix_int
        module procedure log_matrix_string
    end interface log_matrix

    interface log_vector
        module procedure log_vector_double
        module procedure log_vector_int
        module procedure log_vector_string
    end interface log_vector

    interface log_matrix_sparse
        module procedure log_matrix_sparse_double
    end interface log_matrix_sparse

CONTAINS

    function strint(value) result(res)
        integer ,intent(in)  :: value
        character(len=20)                    :: res
        integer                              :: io

        res=''
        write (res,fmt='(i11)',iostat=io) value
        if ( io /= 0 ) then
            res='<NO-DEF>'
        else
            res=adjustl(res)
        end if
        return


    end function strint

    function strfloat(value) result(res)
        double precision ,intent(in)           :: value
        character(len=20)                    :: res
        integer                              :: io

        res=''
        write (res,fmt='(f11.4)',iostat=io) value
        if ( io /= 0 ) then
            res='<NO-DEF>'
        else
            res=adjustl(res)
        end if
        return

    end function strfloat

    function strbool(value) result(res)
        logical ,intent(in)           :: value
        character(len=20)                    :: res

        if ( value ) then
            res = 'True'
        else
            res = 'False'
        end if

    end function strbool

    !-----------------------------
    ! subroutine : init_log
    !
    SUBROUTINE init_log(level_info)
        integer,optional,intent(in)                :: level_info
        integer  :: code

        thread_rang = 0
        IF ( .NOT. present(level_info) ) THEN
            internal_level = WARNING_DEF
        ELSE
            internal_level = level_info
        ENDIF


#ifdef MPI_ACTIVE
     call MPI_COMM_RANK (MPI_COMM_WORLD, thread_rang , code)
     IF ( thread_rang /= 0 ) return
#endif

    END SUBROUTINE init_log

    !-----------------------------
    ! subroutine : init_log
    !
    FUNCTION get_log_level() result(v)
        integer :: v

        v = internal_level

    END FUNCTION get_log_level

    !----------------------------------
    !subroutine end_log
    SUBROUTINE end_log()

    END SUBROUTINE end_log

    !
    SUBROUTINE log_mess(message,level)
        character (len = *) , intent(in)   :: message
        integer,optional,intent(in)        :: level

        IF ( present(level) ) THEN
            IF ( level >= internal_level ) THEN
                CALL WRITE_LOG(message,level)
                RETURN
            ENDIF
        ELSE
            IF ( INFO_DEF >= internal_level ) THEN
                CALL WRITE_LOG(message,INFO_DEF)
                RETURN
            ENDIF
        ENDIF
    END SUBROUTINE log_mess

    SUBROUTINE log_vector_double(message,vector,sn,en,level)
        character (len = *)             ,intent(in)   :: message
        double precision , dimension(:) ,intent(in)   :: vector
        integer                         , intent(in)  :: sn,en ! start index, end index

        integer,optional,intent(in)        :: level
        integer :: lev

        lev=INFO_DEF
        IF ( present(level) ) THEN
            lev = level
        END IF

        IF ( lev >= internal_level ) THEN
            CALL WRITE_LOG(message,lev)
            write (*,fmt="(1000(F7.3,1x))") vector(sn:en)
        ENDIF
    END SUBROUTINE log_vector_double

    SUBROUTINE log_vector_int(message,vector,sn,en,level)
        character (len = *)             ,intent(in)   :: message
        integer , dimension(:) ,intent(in)   :: vector
        integer                         , intent(in)  :: sn,en ! start index, end index

        integer,optional,intent(in)        :: level
        integer :: lev

        lev=INFO_DEF
        IF ( present(level) ) THEN
            lev = level
        END IF

        IF ( lev >= internal_level ) THEN
            CALL WRITE_LOG(message,lev)
            write (*,fmt="(1000(i7,1x))") vector(sn:en)
        ENDIF
    END SUBROUTINE log_vector_int

    SUBROUTINE log_vector_string(message,vector,sn,en,level)
        character (len = *)             ,intent(in)   :: message
        character (len = *) , dimension(:) ,intent(in)   :: vector
        integer                         , intent(in)  :: sn,en ! start index, end index

        integer,optional,intent(in)        :: level
        integer :: lev

        lev=INFO_DEF
        IF ( present(level) ) THEN
            lev = level
        END IF

        IF ( lev >= internal_level ) THEN
            CALL WRITE_LOG(message,lev)
            write (*,fmt="(1000(a7,1x))") vector(sn:en)
        ENDIF
    END SUBROUTINE log_vector_string


    SUBROUTINE log_matrix_double(message,matrix,sn,en,sm,em,level)
        character (len = *)             ,intent(in)   :: message
        double precision , dimension(:,:) ,intent(in) :: matrix
        integer                         , intent(in)  :: sn,en ! start line, end line
        integer                         , intent(in)  :: sm,em ! start column, end column

        integer,optional,intent(in)        :: level
        integer :: lev,i

        lev=INFO_DEF
        IF ( present(level) ) THEN
            lev = level
        END IF

        IF ( lev >= internal_level ) THEN
            CALL WRITE_LOG(message,lev)
            do i=sn,en
                write (*,fmt="(i7,':',1x,1000(F7.3,1x))") i,matrix(i,sm:em)
            end do
        ENDIF
    END SUBROUTINE log_matrix_double


    SUBROUTINE log_matrix_int(message,matrix,sn,en,sm,em,level)
        character (len = *)             ,intent(in)   :: message
        integer , dimension(:,:) ,intent(in) :: matrix
        integer                         , intent(in)  :: sn,en ! start line, end line
        integer                         , intent(in)  :: sm,em ! start column, end column

        integer,optional,intent(in)        :: level
        integer :: lev,i

        lev=INFO_DEF
        IF ( present(level) ) THEN
            lev = level
        END IF

        IF ( lev >= internal_level ) THEN
            CALL WRITE_LOG(message,lev)
            do i=sn,en
                write (*,fmt="(i7,':',1x,1000(i7,1x))") i,matrix(i,sm:em)
            end do
        ENDIF
    END SUBROUTINE log_matrix_int



    SUBROUTINE log_matrix_string(message,matrix,sn,en,sm,em,level)
        character (len = *)             ,intent(in)   :: message
        character (len = *), dimension(:,:) ,intent(in) :: matrix
        integer                         , intent(in)  :: sn,en ! start line, end line
        integer                         , intent(in)  :: sm,em ! start column, end column

        integer,optional,intent(in)        :: level
        integer :: lev,i

        lev=INFO_DEF
        IF ( present(level) ) THEN
            lev = level
        END IF

        IF ( lev >= internal_level ) THEN
            CALL WRITE_LOG(message,lev)
            do i=sn,en
                write (*,fmt="(i7,':',1x,1000(a7,1x))") i,matrix(i,sm:em)
            end do
        ENDIF
    END SUBROUTINE log_matrix_string


    SUBROUTINE log_matrix_sparse_double(message,matrix,sn,en,sm,em,level)
        character (len = *)             ,intent(in)   :: message
        double precision , dimension(:,:) ,intent(in) :: matrix
        integer                         , intent(in)  :: sn,en ! start line, end line
        integer                         , intent(in)  :: sm,em ! start column, end column

        integer,optional,intent(in)        :: level
        integer :: lev,i,j

        lev=INFO_DEF
        IF ( present(level) ) THEN
            lev = level
        END IF

        IF ( lev >= internal_level ) THEN
            CALL WRITE_LOG(message,lev)
            write (*,*) "** number of values :",count(matrix /= 0.d0)
            do i=sn,en
                write (*,*) "line ",i," number of values :",count(matrix(i,:) /= 0.d0)
                do j=sm,em
                    if ( matrix(i,j) /= 0.d0 ) then
                        write (*,fmt="(i7,',',i7,':',1x,(F7.3,1x))") i,j,matrix(i,j)
                    end if
                end do
            end do
        ENDIF
    END SUBROUTINE log_matrix_sparse_double



    !----------------------------------

    SUBROUTINE WRITE_LOG(message,level)
        character (len = *) , intent(in)   :: message
        integer,         intent(in)        :: level
        character(len=8)               :: date
        character(len=10)              :: time
        character(len=50)              :: addinfo
        character(len=800)             :: line
        character(len=5)               :: sline

        addinfo=''

#ifdef MPI_ACTIVE
        addinfo='[mpi_proc '//trim(strint(thread_rang))//']'
#endif

#ifdef MANAGE_OMP
        if ( omp_in_parallel() ) then
         addinfo=trim(addinfo)//'[omp_thread '//trim(strint(omp_get_thread_num()+1))//']'
        end if
#endif
        if ( level == WARNING_DEF ) then
            write (*,fmt="(a,a,a,a,a)") achar(27),"[22;35m",'WARNING:',message,default_col
        else if ( level == ERROR_DEF ) then
            write (0,fmt="(a,a,a,a,a)") achar(27),"[31;01m",'ERROR:',message,default_col
        else
            call DATE_AND_TIME(date,time)

            line='qtlmap:['// &
                date(7:8)//'/'//date(5:6)//'/'//date(1:4)//' '//&
                time(1:2)//':'//time(3:4)//':'//time(5:6)// &
                ']['//trim(TOSTRING(level))//']'//&
                trim(addinfo)//':'//trim(message)

           if ( level == DEBUG_DEF ) then
            write (*,fmt="(a)") debug_col//trim(line)//default_col
           else if ( level == VERBOSE_DEF ) then
            write (*,fmt="(a)") verbose_col//trim(line)//default_col
           else
            write (*,fmt="(a)") trim(line)
           end if
        end if

    END SUBROUTINE WRITE_LOG


    !----------------------------------

    FUNCTION TOSTRING(level) result(slevel)
        integer,         intent(in)        :: level
        character(len=30)                  :: slevel
        SELECT CASE(level)
            CASE(DEBUG_DEF)
                slevel = "debug  "
            CASE(VERBOSE_DEF)
                slevel = "verbose"
            CASE(INFO_DEF)
                slevel = "info   "
            CASE(WARNING_DEF)
                slevel = "warning"
            CASE(ERROR_DEF)
                slevel = "error  "
            CASE DEFAULT
                slevel = "*******"
        END SELECT

        RETURN
    END FUNCTION TOSTRING


END MODULE m_qtlmap_log
