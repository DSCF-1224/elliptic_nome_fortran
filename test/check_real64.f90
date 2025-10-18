program check_real64

    use, intrinsic :: iso_fortran_env , only: compiler_options
    use, intrinsic :: iso_fortran_env , only: compiler_version
    use, intrinsic :: iso_fortran_env , only: error_unit
    use, intrinsic :: iso_fortran_env , only: real64

    use, intrinsic :: ieee_arithmetic , only: operator(.ne.)
    use, intrinsic :: ieee_arithmetic , only: ieee_class
    use, intrinsic :: ieee_arithmetic , only: ieee_positive_zero

    use, non_intrinsic :: elliptic_nome_fortran



    implicit none



    integer, parameter :: len_str = 44

    character(*), parameter :: file_k = 'test/k.dat'

    character(*), parameter :: file_q = 'test/q.dat'



    call check()



    contains



    subroutine check

        integer :: cnt_k, cnt_q

        integer :: file_unit

        integer :: itr

        integer :: stat

        real(real64), allocatable, dimension(:) :: num_k, num_q_ref

        real(real64), allocatable, dimension(:,:) :: num_q_err, num_q_res

        character(256) :: msg

        character(len_str), allocatable, dimension(:) :: str_k, str_q

        character(:), allocatable :: str_compiler_options



        call count_values( &!
        file    = file_k , &!
        iostat  = stat   , &!
        iomsg   = msg    , &!
        counter = cnt_k    &!
        )

        call count_values( &!!
        file    = file_q , &!
        iostat  = stat   , &!
        iomsg   = msg    , &!
        counter = cnt_q    &!
        )



        if (cnt_k .ne. cnt_q) then

            write( unit = error_unit, fmt = * ) 'count k: ', cnt_k
            write( unit = error_unit, fmt = * ) 'count q: ', cnt_q
            write( unit = error_unit, fmt = * ) 'The size of k & q must be equal.'
            error stop

        else if ( cnt_k .lt. 1 ) then

            write( unit = error_unit, fmt = * ) 'counter: ', cnt_k
            write( unit = error_unit, fmt = * ) 'The size of vals must be greater than zero.'
            error stop

        end if



        allocate( num_q_err(10, cnt_q), stat = stat, errmsg = msg )

        call export_msg(stat, msg)

        allocate( num_q_res(10, cnt_q), stat = stat, errmsg = msg )

        call export_msg(stat, msg)



        call read_values( &!
        file    = file_k , &!
        counter = cnt_k  , &!
        stat    = stat   , &!
        msg     = msg    , &!
        num     = num_k  , &!
        str     = str_k    &!
        )

        call read_values( &!
        file    = file_q    , &!
        counter = cnt_q     , &!
        stat    = stat      , &!
        msg     = msg       , &!
        num     = num_q_ref , &!
        str     = str_q       &!
        )



        do itr = 1, cnt_k

            call check_kernel(&!
            k     = num_k     (     itr ) , &!
            q_ref = num_q_ref (     itr ) , &!
            q_res = num_q_res ( : , itr ) , &!
            q_err = num_q_err ( : , itr )   &!
            )

        end do



        open( &!
        file    = 'check_real64.md' , &!
        newunit = file_unit , &!
        action  = 'write'   , &!
        status  = 'unknown' , &!
        iostat  = stat      , &!
        iomsg   = msg         &!
        )

        call export_msg(stat, msg)



        str_compiler_options = compiler_options()

        block

            integer :: i

            do

                i = index(str_compiler_options, ' -')

                print *, i

                if (i .eq. 0) exit

                str_compiler_options(i:i) = new_line('')

            end do

        end block



        write( &!
        unit   = file_unit , &!
        fmt    = '(A)'     , &!
        iostat = stat      , &!
        iomsg  = msg         &!
        ) &!
        '# check'            // new_line('') , &!
        '[compiler version]' // new_line('') , &!
        '```'                                , &!
        compiler_version()                   , &!
        '```'                // new_line('') , &!
        '[compiler options]' // new_line('') , &!
        '```'                                , &!
        str_compiler_options                 , &!
        '```'

        call export_msg(stat, msg)



        do itr = 1, cnt_k

            write( &!
            unit   = file_unit , &!
            fmt    = '(A)'     , &!
            iostat = stat      , &!
            iomsg  = msg         &!
            ) &!
            new_line('')         , &!
            '## `k = '        //   &!
            trim(str_k(itr))  //   &!
            '_real64`'

            call export_msg(stat, msg)



            if (ieee_class(num_q_ref(itr)) .ne. ieee_positive_zero) then

                write( &!
                unit   = file_unit , &!
                fmt    = '(A)'     , &!
                iostat = stat      , &!
                iomsg  = msg         &!
                ) &!
                new_line('') // '|degree|res|(relative err)/`machine epsilon`|' , &!
                                '|:----:|:--|:-------------------------------|'

            else

                write( &!
                unit   = file_unit , &!
                fmt    = '(A)'     , &!
                iostat = stat      , &!
                iomsg  = msg         &!
                ) &!
                new_line('') // '|degree|res|err|' , &!
                                '|:----:|:--|:--|'
            end if

            call export_msg(stat, msg)



            write( &!
            unit   = file_unit , &!
            fmt    = '(SP,A,A,A,ES24.16E3,A,ES24.16E3,A)' , &!
            iostat = stat      , &!
            iomsg  = msg         &!
            ) &!
            '|' , '  01' , '|' , num_q_res( 1,itr) , '|' , num_q_err( 1,itr) , '|' , &!
            '|' , '  05' , '|' , num_q_res( 2,itr) , '|' , num_q_err( 2,itr) , '|' , &!
            '|' , '  09' , '|' , num_q_res( 3,itr) , '|' , num_q_err( 3,itr) , '|' , &!
            '|' , '  13' , '|' , num_q_res( 4,itr) , '|' , num_q_err( 4,itr) , '|' , &!
            '|' , '  17' , '|' , num_q_res( 5,itr) , '|' , num_q_err( 5,itr) , '|' , &!
            '|' , '  21' , '|' , num_q_res( 6,itr) , '|' , num_q_err( 6,itr) , '|' , &!
            '|' , '  25' , '|' , num_q_res( 7,itr) , '|' , num_q_err( 7,itr) , '|' , &!
            '|' , '  29' , '|' , num_q_res( 8,itr) , '|' , num_q_err( 8,itr) , '|' , &!
            '|' , '  33' , '|' , num_q_res( 9,itr) , '|' , num_q_err( 9,itr) , '|' , &!
            '|' , 'auto' , '|' , num_q_res(10,itr) , '|' , num_q_err(10,itr) , '|' , &!
            '|' , 'ref ' , '|' , num_q_ref(   itr) , '|' , 0.0_real64  , '|'

            call export_msg(stat, msg)

        end do

    end subroutine check



    subroutine check_kernel(k, q_ref, q_res, q_err)

        real(real64), intent(in) :: k

        real(real64), intent(in) :: q_ref

        real(real64), intent(out), dimension(10) :: q_res, q_err



        q_res( 1) = elliptic_nome_01   (k)
        q_res( 2) = elliptic_nome_05   (k)
        q_res( 3) = elliptic_nome_09   (k)
        q_res( 4) = elliptic_nome_13   (k)
        q_res( 5) = elliptic_nome_17   (k)
        q_res( 6) = elliptic_nome_21   (k)
        q_res( 7) = elliptic_nome_25   (k)
        q_res( 8) = elliptic_nome_29   (k)
        q_res( 9) = elliptic_nome_33   (k)
        q_res(10) = elliptic_nome_auto (k)

        q_err(:) = q_res(:) - q_ref

        if (ieee_class(q_ref) .ne. ieee_positive_zero) then
            q_err(:) = q_err(:) / q_ref / epsilon(k)
        end if

    end subroutine check_kernel



    subroutine count_values(file, iostat, iomsg, counter)

        character(*), intent(in) :: file

        integer, intent(out) :: iostat

        character(*), intent(inout) :: iomsg

        integer, intent(out) :: counter



        integer :: file_unit

        real(real64) :: dummy



        open( &!
        file    = file        , &!
        newunit = file_unit   , &!
        action  = 'read'      , &!
        form    = 'formatted' , &!
        status  = 'old'       , &!
        iostat  = iostat      , &!
        iomsg   = iomsg         &!
        )

        call export_msg(iostat, iomsg)



        counter = 0

        do

            read( &!
            unit   = file_unit , &!
            fmt    = *         , &!
            iostat = iostat    , &!
            iomsg  = iomsg       &!
            ) &!
            dummy

            if ( iostat .eq. 0 ) then
                counter = &!
                counter + 1
            else if ( is_iostat_end(iostat) ) then
                return
            else
                call export_msg(iostat, iomsg)
            end if

        end do



        close(file_unit)

    end subroutine count_values



    subroutine export_msg(stat, msg)

        integer, intent(in) :: stat

        character(*), intent(in) :: msg



        if ( stat .ne. 0 ) then
            write( unit = error_unit, fmt = '(A,I0)' ) 'stat : ', stat
            write( unit = error_unit, fmt = '(A,A )' ) 'msg  : ', trim(msg)
            error stop
        end if

    end subroutine export_msg



    subroutine read_values(file, counter, stat, msg, num, str)

        character(*), intent(in) :: file

        integer, intent(in) :: counter

        integer, intent(out) :: stat

        character(*), intent(inout) :: msg

        real(real64), intent(inout), allocatable :: num(:)

        character(len_str), intent(inout), allocatable :: str(:)



        integer :: file_unit

        integer :: itr



        allocate( num(counter), stat = stat, errmsg = msg )

        call export_msg(stat, msg)

        allocate( str(counter), stat = stat, errmsg = msg )

        call export_msg(stat, msg)



        open( &!
        file    = file        , &!
        newunit = file_unit   , &!
        action  = 'read'      , &!
        form    = 'formatted' , &!
        status  = 'old'       , &!
        iostat  = stat        , &!
        iomsg   = msg           &!
        )

        if ( stat .ne. 0 ) then
            write( unit = error_unit, fmt = '(A)' ) trim(msg)
            error stop
        end if



        do itr = 1, counter

            read( &!
            unit   = file_unit , &!
            fmt    = '(A)'     , &!
            iostat  = stat     , &!
            iomsg   = msg        &!
            ) &!
            str(itr)

            if ( stat .ne. 0 ) then
                write( unit = error_unit, fmt = '(A)' ) trim(msg)
                error stop
            end if

            read( &!
            unit   = str(itr) , &!
            fmt    = *        , &!
            iostat  = stat    , &!
            iomsg   = msg       &!
            ) &!
            num(itr)

            if ( stat .ne. 0 ) then
                write( unit = error_unit, fmt = '(A)' ) trim(msg)
                error stop
            end if

        end do



        close(file_unit)

    end subroutine read_values

end program check_real64
