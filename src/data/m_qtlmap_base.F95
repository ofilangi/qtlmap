!!****m* DATA/m_qtlmap_base
!!  NAME
!!    m_qtlmap_tools
!!  DESCRIPTION
!!   misc. routines.
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
module m_qtlmap_base
    use m_qtlmap_constant
    use m_qtlmap_log

    implicit none
    private

    public :: file_exist
    public :: get_real
    public :: get_int
    public :: check_allocate
    public :: next_word
    public :: parse_real_array
    public :: xaldane
    public :: xosambi
    public :: stop_on_error
    public :: stop_application
    public :: str

#ifndef VARYING_STRING_MOD
    public :: REMOVE
    public :: EXTRACT
    public :: GET
#endif


    !!****f* m_qtlmap_tools/str
    !! NAME
    !!   str
    !! DESCRIPTION
    !!   generic interface to convert a float/double or an integer to a string
    !! SOURCE
    interface str
        module procedure iget_str, &
            fget_str, &
            dget_str
    end interface str
    !!***

    !!****f* m_qtlmap_tools/get_real
    !! NAME
    !!   get_real
    !! DESCRIPTION
    !!   generic interface to convert a word into a float
    !! SOURCE
    interface get_real
        module procedure get_real_char
    end interface get_real
    !!***

    !!****f* m_qtlmap_tools/get_int
    !! NAME
    !!   get_int
    !! DESCRIPTION
    !!   generic interface to convert a word into an integer
    !! SOURCE
    interface get_int
        module procedure get_int_char
    end interface get_int
    !!***

    !!****f* m_qtlmap_tools/next_word
    !! NAME
    !!   next_word
    !! DESCRIPTION
    !!   generic interface to read a word from line (parsing context)
    !! SOURCE
    interface next_word
        module procedure next_word_char
    end interface next_word
!!***

contains

    !!****f* m_qtlmap_tools/file_exist
    !! NAME
    !!   file_exist
    !! DESCRIPTION
    !!   Test if a local file exist (the programe stop if not)
    !! INPUTS
    !!  file_in              : the pathname file
    !!
    !! SOURCE
    !!
    subroutine file_exist(file_in)
        character(len=*),intent(in) ::file_in
        logical::existe
        inquire( file=file_in, exist=existe)

        if(.not. existe) then
            call stop_application('The file ['//trim(file_in)//'] does not exist.')
        endif
    end subroutine file_exist

    !! DESCRIPTION
    !!   stop the application with an error message on a read file context
    !! INPUTS
    !!  iostat              : iostat error code
    !!  file                : the pathname file
    !!  lineNumber          : the number of the current line
    !!  message             : a message error
    subroutine stop_on_error(iostat,file,lineNumber,message)
        integer , intent(in)        ::iostat
        character(len=*),intent(in) ::file
        integer , intent(in)        ::lineNumber
        character(len=*),intent(in) ::message

        if (iostat /= 0 .and. iostat/=-2) then
            call log_mess('*** Read error ['//trim(file)//'] ' // &
                'line '//trim(iget_str(lineNumber))//':bad specification [message:'//adjustr(message)//'] ***',ERROR_DEF)
            stop 'read error'
        endif

    end subroutine stop_on_error
    !! DESCRIPTION
    !!   stop the application with an error message
    !! INPUTS
    !!  message              : an error message
    subroutine stop_application(message)
        character(len=*),intent(in) ::message
        call log_mess('*** message error :'//message//' ***',ERROR_DEF)
        stop ' - Program terminated'
    end subroutine stop_application

    !! DESCRIPTION
    !!   get a set of real value from a line
    !! INPUTS
    !!  nvalue               : number of value
    !! INPUTS/OUTPUT
    !!  line_to_parse        : the line
    !! OUTPUTS
    !!  tabV                 : a vector of double precision
    !!  is_ok                : a flag error
    subroutine parse_real_array(line_to_parse,nvalue,tabV,is_ok)
        character(len=*),intent(inout)                 :: line_to_parse
        integer, intent(in)                            :: nvalue
        real(kind=dp) ,dimension(nvalue),intent(out)   :: tabV
        logical, intent(inout),optional                :: is_ok
        character(len=LEN_BUFFER_WORD)                 :: w_token
        integer                                        :: i,ns,ne
        character(len=1000)                            :: endline

        if ( trim(line_to_parse)=='' ) then
            if ( present(is_ok) ) is_ok = .false.
            return
        end if

        ns=index(line_to_parse,"[")
        if ( ns == 0 ) then
            if ( present(is_ok) ) is_ok = .false.
            return
        end if
        ne=index(line_to_parse,"]")

        if ( ne == 0 ) then
            if ( present(is_ok) ) is_ok = .false.
            return
        end if

        endline = line_to_parse(ne+1:)
        line_to_parse = line_to_parse(ns+1:ne-1)

        do i=1,nvalue
            w_token=trim(next_word_char(line_to_parse,is_ok))
            if (is_ok) tabV(i)=get_real(w_token,is_ok)
            if ( .not. is_ok ) then
                if ( present(is_ok) ) is_ok = .false.
                return
            end if
        end do

        if ( trim(line_to_parse) /= '' ) then
            if ( present(is_ok) ) is_ok = .false.
            return
        end if

        line_to_parse = endline

    end subroutine parse_real_array
    !!***

    !!****f* m_qtlmap_tools/next_word_char
    !! NAME
    !!    next_word_char
    !! DESCRIPTION
    !!    get the last word at the right side of a line and remove it.
    !!
    !! INPUTS/OUTPUT
    !!  line_to_parse        : the line.
    !! OUTPUTS
    !!  is_ok                : a flag error
    !! RETURN
    !!  the word
    !! SOURCE
    !!
    FUNCTION next_word_char(line_to_parse,is_ok) result(w_token)

        character(len=*) ,intent(inout)        :: line_to_parse
        logical, intent(inout),optional        :: is_ok
        character(len=LEN_BUFFER_WORD)         :: w_token
        integer                                :: nd

        nd = 1
        w_token = ''
        line_to_parse=adjustl(line_to_parse)

        if ( trim(line_to_parse)=='' ) then
            if ( present(is_ok) ) is_ok = .false.
            return
        end if


        do while (nd == 1 .and. w_token == '')
            !! char(9) pour les tabulation.....
            nd=SCAN(line_to_parse," ;"//char(9))

            if ( nd == 0) then
                w_token = trim(line_to_parse)
            else
                w_token = trim(EXTRACT(line_to_parse,1,nd-1))
            end if


            if (w_token == '') then
                line_to_parse = REMOVE(line_to_parse,1,nd)
            end if
        end do

        if  (nd == 0) then ! END OF LINE
            w_token = trim(EXTRACT(line_to_parse,1))
            line_to_parse = ''
            if ( present(is_ok) ) is_ok = .true.
            return
        end if

        if (nd > 0) then  ! at least one non-terminator character in the word
            w_token = trim(EXTRACT(line_to_parse,1,nd-1))
            if ( present(is_ok) ) is_ok = .true.
        else
            w_token = ''
            if ( present(is_ok) ) is_ok = .false.
        end if


        line_to_parse = adjustl(REMOVE(line_to_parse,1,nd))
        return
    end function next_word_char
    !!***

    !!****f* m_qtlmap_tools/check_allocate
    !! NAME
    !!    check_allocate
    !! DESCRIPTION
    !!    check a iostat error code. print an error message and stop the programe if the iostat error /= 0
    !! INPUTS
    !!  iostat        : the iostat error code
    !!  error_message : an error message
    !!
    !! OUTPUTS
    !!  is_ok             : boolean
    !! RETURN
    !!   a double precision
    !! SOURCE
    !!
    subroutine check_allocate(iostat,error_message)
        !ios_id where read the parameter
        integer,intent(in)::iostat
        !message to print if iostat /= 0
        character(len=*),optional,intent(in)::error_message

        if (iostat /= 0) then
            print *, "Allocation error:[iostat:",iostat,"]"
            if ( present (error_message) ) then
                print *,'**'//error_message//'**'
            endif
            stop 'allocation error'
        end if
    end subroutine check_allocate
    !!***

    !!****f* m_qtlmap_tools/get_real_char
    !! NAME
    !!    get_real_char
    !! DESCRIPTION
    !!    conversion of a character word in a double precision using the fortran read routine
    !! INPUTS
    !!  string_to_convert : the word
    !!
    !! OUTPUTS
    !!  is_ok             : boolean
    !! RETURN
    !!   a double precision
    !! SOURCE
    !!
    function get_real_char(string_to_convert,is_ok) result(real_value)
        character(len=*),intent(in)      :: string_to_convert
        logical,optional,intent(inout)       :: is_ok
        real (kind=dp)                       :: real_value
        character(len=LEN_BUFFER_WORD)       :: word_token
        integer                              :: ios
        character(len=LEN_DEF)               :: c,string
        logical                              :: is_int

        string = trim(string_to_convert)
        is_int=.true.

        do while ( string /='' )
            c = EXTRACT (string , 1, 1)
            if (c == '.') then
                is_int=.false.
                exit
            end if
            string = REMOVE(string,1,1)
        end do

        if ( .not. is_int ) then
            word_token = string_to_convert!CHAR(string_to_convert)
            read (word_token,FMT_REAL,IOSTAT=ios) real_value
            if ( present(is_ok) ) then
                is_ok = ( ios == 0 )
            end if
            if (  ios /= 0 ) real_value = REAL_NOT_DEFINED
        else
            if (present(is_ok)) then
                real_value = get_int(string_to_convert,is_ok)
            else
                real_value = get_int(string_to_convert)
            end if
        end if
        string=''
        return

    end function get_real_char
    !!***

    !!****f* m_qtlmap_tools/get_int_char
    !! NAME
    !!    get_int_char
    !! DESCRIPTION
    !!    conversion of a character word in an integer using the fortran read routine
    !! INPUTS
    !!  string_to_convert : the word
    !!
    !! OUTPUTS
    !!  is_ok             : boolean
    !! RETURN
    !!   an integer
    !! SOURCE
    !!
    function get_int_char(string_to_convert,is_ok) result(int_value)
        character(len=*),intent(in)      :: string_to_convert
        logical,optional,intent(inout)       :: is_ok
        integer                              :: int_value
        character(len=LEN_BUFFER_WORD)       :: word_token
        integer                              :: ios

        word_token = trim(string_to_convert)!CHAR(string_to_convert)
        read (word_token,FMT_INT,IOSTAT=ios) int_value

        if ( present(is_ok) ) then
            is_ok = ( ios == 0 )
        endif

        if (  ios /= 0 ) int_value = INT_NOT_DEFINED
        word_token =''
        return

    end function get_int_char
    !!***

    !!****f* m_qtlmap_tools/iget_str
    !! NAME
    !!    iget_str
    !! DESCRIPTION
    !!    conversion of an integer into a character word
    !! INPUTS
    !!  value_int : the word
    !!
    !! RETURN
    !!   an integer
    !! SOURCE
    !!
    function iget_str(value_int) result(res)
        integer,intent(in)                   :: value_int
        character(len=20)                    :: res
        character(len=LEN_DEF)                 :: str_value,sign_value,tempvs,tempvs2
        integer                              :: v,reste,ascii_zero

        if (value_int<0) then
            sign_value='-'
            v = -value_int
        else
            sign_value=''
            v = value_int
        end if

        ascii_zero = iachar('0')
        str_value=''
        do while (v /=0 )
            reste = mod(v,10)
            tempvs=achar(ascii_zero+reste)
            tempvs2=trim(str_value)
            str_value = trim(tempvs)//trim(tempvs2)
            v = v / 10
        end do

        str_value = trim(sign_value)//trim(str_value)
        if (str_value=='' .or. str_value=='-') str_value='0'
        res = trim(str_value)

        return

    end function iget_str
    !!***

    !!****f* m_qtlmap_tools/fget_str
    !! NAME
    !!    fget_str
    !! DESCRIPTION
    !!    conversion of a a float into a character word
    !! INPUTS
    !!  value : the word
    !!
    !! RETURN
    !!   a float
    !! SOURCE
    !!
    function fget_str(value) result(res)
        real           ,intent(in)           :: value
        character(len=20)                    :: res
        real (kind=dp)                       :: v

        v = value
        res = dget_str(v)

    end function fget_str
    !!***

    !!****f* m_qtlmap_tools/dget_str
    !! NAME
    !!    dget_str
    !! DESCRIPTION
    !!    conversion of a double precision into a character word
    !! INPUTS
    !!  value : the word
    !!
    !! RETURN
    !!   a double precision
    !! SOURCE
    !!
    function dget_str(value) result(res)
        real (kind=dp) ,intent(in)           :: value
        character(len=20)                    :: res
        character(len=LEN_DEF)               :: str_value,c,sign_value,tempvs,tempvs2,temppoint
        integer                              :: ascii_zero,reste,decimal
        integer                              :: val_int,io
        real(kind=dp)                        :: v
        !real                                 :: v
        logical                              :: have_dec

        res=''
        write (res,fmt='(f13.2)',iostat=io) value
        if ( io /= 0 ) then
            res='<NO-DEF>'
        else
            res=adjustl(res)
        end if
        return

    end function dget_str
    !!***

    !!****f* m_qtlmap_tools/xaldane
    !! NAME
    !!    xaldane
    !! DESCRIPTION
    !!
    !! INPUTS
    !!  d :
    !!
    !! RETURN
    !!
    !! NOTES
    !!     Calcul du taux de recombinaison en fonction de la distance de xaldane
    !! SOURCE
    !!
    function xaldane(d) result (val)
        real (kind=dp)   , intent(in)  :: d
        real  (kind=dp)                :: val
        val=0.5d0*(1.d0-dexp(-2.d0*d))
        return
    end function xaldane
    !!***

    !!****f* m_qtlmap_tools/xosambi
    !! NAME
    !!    xosambi
    !! DESCRIPTION
    !!
    !! INPUTS
    !!  d :
    !!
    !! RETURN
    !!
    !! NOTES
    !!     Calcul du taux de recombinaison en fonction de la distance de Kosambi
    !! SOURCE
    !!
    function xosambi(d) result (val)
        real (kind=dp)   , intent(in)  :: d
        real  (kind=dp)                :: val,z

        z=dexp(4.d0*d)
        val=0.5d0*((z-1.d0)/(z+1.d0))
        return
    end function xosambi
    !!***

    !!****f* m_qtlmap_tools/REMOVE
    !! NAME
    !!    REMOVE
    !! DESCRIPTION
    !!    Remove a subset of a line
    !! INPUTS
    !!  string : the line
    !!  start  : start index of the subset
    !!  finish : end index of the subset
    !!
    !! RETURN
    !!  the new line without the subset
    !!
    !! SOURCE
    !!
    function REMOVE (string,start,finish) result(out)
        character(len=*) , intent(in)         :: string
        integer , intent(in) ,optional        :: start
        integer , intent(in) ,optional        :: finish

        character(len=len(string))                :: out
        integer                               :: s,e

        s = 1
        e = len(string)
        if (present(start)) then
            if (start > 1 .and. start <= len(string)) then
                s = start
            end if
        end if

        if (present(finish)) then
            if (finish >= start  .and. finish <= len(string)) then
                e = finish
            end if
        end if

        out=''

        if ( s > 1) then
            out=trim(string(1:s))
        end if

        if ( e < len(string) ) then
            out = trim(out) // string(e+1:len(string))
        end if

    end function REMOVE
    !!***

    !!****f* m_qtlmap_tools/EXTRACT
    !! NAME
    !!    EXTRACT
    !! DESCRIPTION
    !!    Get a subset of a character line
    !! INPUTS
    !!  string : the line
    !!  start  : start index of the subset
    !!  finish : end index of the subset
    !!
    !! RETURN
    !!  the subset of the line
    !!
    !! SOURCE
    !!
    function EXTRACT (string,start,finish) result(out)
        character(len=*) , intent(in)         :: string
        integer , intent(in) ,optional        :: start
        integer , intent(in) ,optional        :: finish

        character(len=len(string))                :: out
        integer                               :: s,e

        s=1
        e=len(string)

        if (present(start)) then
            if (start > 1 .and. start <= len(string)) then
                s = start
            end if
        end if

        if (present(finish)) then
            if (finish >= start  .and. finish <= len(string)) then
                e = finish
            end if
        end if

        out=trim(string(s:e))

    end function
    !!***

    !!****f* m_qtlmap_tools/GET
    !! NAME
    !!    GET
    !! DESCRIPTION
    !!    Read a line from an unit file
    !! INPUTS/OUTPUTS
    !!   string : the line read
    !! INPUTS
    !!   unit   : the unit record
    !!  maxlen  : maximum length of the line
    !!
    !! OUTPUTS
    !!   iostat : iostat error
    !! SOURCE
    !!
    subroutine GET (unit,string,maxlen,iostat)
        integer                , intent(in)      :: unit
        integer , intent(in)                     :: maxlen
        character(len=maxlen)  , intent(inout)   :: string
        integer , intent(out),optional           :: iostat

        character(len=LEN_LINE)                   :: word
        integer :: eof

        eof = 0
        string=''
        read (unit,iostat=eof,fmt="(a)/") string
        if ( present(iostat) ) iostat = eof

    end subroutine GET
!!***

end module m_qtlmap_base
