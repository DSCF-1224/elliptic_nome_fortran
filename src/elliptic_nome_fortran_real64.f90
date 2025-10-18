module elliptic_nome_fortran_real64
    !! Compute nome \( q \) for elliptic integrals corresponding to modulus \( k \)
    !! $$
    !! \begin{aligned}
    !! K (k) &:= \int_{ 0 }^{ \pi / 2 } \frac{ d \theta }{ \sqrt{ 1 - { k }^{ 2 } \sin^{ 2 } \theta } } \\
    !! { k }^{ \prime } &:= \sqrt{ 1 - { k }^{ 2 } } \\
    !! { K }^{ \prime } (k) &:= K ( { k }^{ \prime } ) \\
    !! q &:= e^{ - \pi { K }^{ \prime } (k) / K (k) } \\
    !! \varepsilon &:= \frac{ 1 }{ 2 } \frac{ 1 - \sqrt{ { k }^{ \prime } } }{ 1 + \sqrt{ { k }^{ \prime } } }
    !! \end{aligned}
    !! $$
    !! @note
    !! - 山内二郎, 宇野利雄, 一松信 共編  
    !!   電子計算機のための数値計算法 3  
    !!   培風館, 1972.  
    !!   数理科学シリーズ ; 5  
    !!   [NDLサーチ](https://ndlsearch.ndl.go.jp/books/R100000039-I2422322)
    !! - [A002103 - OEIS](https://oeis.org/A002103)
    !! @endnote

    use, intrinsic :: iso_fortran_env , only: real64

    use, intrinsic :: ieee_arithmetic , only: ieee_quiet_nan
    use, intrinsic :: ieee_arithmetic , only: ieee_value



    implicit none



    private



    public :: elliptic_nome_01_real64
    public :: elliptic_nome_05_real64
    public :: elliptic_nome_09_real64
    public :: elliptic_nome_13_real64
    public :: elliptic_nome_17_real64
    public :: elliptic_nome_21_real64
    public :: elliptic_nome_25_real64
    public :: elliptic_nome_29_real64
    public :: elliptic_nome_33_real64
    public :: elliptic_nome_auto_real64



    real(real64), parameter :: c_01 =        1.0_real64
    real(real64), parameter :: c_05 =        2.0_real64
    real(real64), parameter :: c_09 =       15.0_real64
    real(real64), parameter :: c_13 =      150.0_real64
    real(real64), parameter :: c_17 =     1707.0_real64
    real(real64), parameter :: c_21 =    20910.0_real64
    real(real64), parameter :: c_25 =   268616.0_real64
    real(real64), parameter :: c_29 =  3567400.0_real64
    real(real64), parameter :: c_33 = 48555069.0_real64

    real(real64), parameter :: q_err = tiny(0.0_real64)



    contains



    elemental function elliptic_nome_01_real64(k) result(q)
        !! Calculate the elliptic nome \( q \) using
        !! the Horner's method and
        !! the following polynomial:
        !! $$
        !! \begin{align*}
        !! \varepsilon &:= \frac{ 1 }{ 2 } \frac{ 1 - \sqrt{ { k }^{ \prime } } }{ 1 + \sqrt{ { k }^{ \prime } } } \\
        !! q(\varepsilon) &:=            \varepsilon
        !! \end{align*}
        !! $$
        !! @note
        !! - The elliptic modulus \( k \) should satisfy \( { k }^{ 2 } \le 1/2 \).
        !! @endnote

        real(real64), intent(in) :: k !! elliptic modulus \( k \)



        real(real64) :: q !! elliptic nome \( q \)



        real(real64) :: comp_k !! \( { k }^{ \prime } \)

        real(real64) :: pw02_k !! \( { k }^{ 2 } \)

        real(real64) :: sqrt_comp_k !! \( \sqrt{ { k }^{ \prime } } \)



        call evaluate_modulus(k = k, pw02_k = pw02_k, comp_k = comp_k)

        call calculate_pw01_epsilon(    &!
        &        pw02_k   =      pw02_k , &!
        &        comp_k   =      comp_k , &!
        &   sqrt_comp_k   = sqrt_comp_k , &!
        &        pw01_eps =           q   &!
        )

    end function elliptic_nome_01_real64



    elemental function elliptic_nome_05_real64(k) result(q)
        !! Calculate the elliptic nome \( q \) using
        !! the Horner's method and
        !! the following polynomial:
        !! $$
        !! \begin{align*}
        !! \varepsilon &:= \frac{ 1 }{ 2 } \frac{ 1 - \sqrt{ { k }^{ \prime } } }{ 1 + \sqrt{ { k }^{ \prime } } } \\
        !! q(\varepsilon) &:=            \varepsilon
        !!             \\ & +        2 { \varepsilon }^{  5 }
        !! \end{align*}
        !! $$
        !! @note
        !! - The elliptic modulus \( k \) should satisfy \( { k }^{ 2 } \le 1/2 \).
        !! @endnote

        real(real64), intent(in) :: k !! elliptic modulus \( k \)



        real(real64) :: q !! elliptic nome \( q \)



        real(real64) :: comp_k !! \( { k }^{ \prime } \)

        real(real64) :: pw02_k !! \( { k }^{ 2 } \)

        real(real64) :: sqrt_comp_k !! \( \sqrt{ { k }^{ \prime } } \)

        real(real64) :: pw01_eps !! auxiliary parameter \( \varepsilon \)

        real(real64) :: pw04_eps !! \( { \varepsilon }^{ 4 } \)



        call evaluate_modulus(k = k, pw02_k = pw02_k, comp_k = comp_k)

        call calculate_pw04_epsilon( &!
        &        pw02_k   =      pw02_k   , &!
        &        comp_k   =      comp_k   , &!
        &   sqrt_comp_k   = sqrt_comp_k   , &!
        &        pw01_eps =      pw01_eps , &!
        &        pw04_eps =      pw04_eps   &!
        )

        q = &!
            elliptic_nome_by_epsilon_05( &!
                pw01_eps = pw01_eps , &!
                pw04_eps = pw04_eps   &!
            )

    end function elliptic_nome_05_real64



    elemental function elliptic_nome_09_real64(k) result(q)
        !! Calculate the elliptic nome \( q \) using
        !! the Horner's method and
        !! the following polynomial:
        !! $$
        !! \begin{align*}
        !! \varepsilon &:= \frac{ 1 }{ 2 } \frac{ 1 - \sqrt{ { k }^{ \prime } } }{ 1 + \sqrt{ { k }^{ \prime } } } \\
        !! q(\varepsilon) &:=            \varepsilon
        !!             \\ & +        2 { \varepsilon }^{  5 }
        !!             \\ & +       15 { \varepsilon }^{  9 }
        !! \end{align*}
        !! $$
        !! @note
        !! - The elliptic modulus \( k \) should satisfy \( { k }^{ 2 } \le 1/2 \).
        !! @endnote

        real(real64), intent(in) :: k !! elliptic modulus \( k \)



        real(real64) :: q !! elliptic nome \( q \)



        real(real64) :: comp_k !! \( { k }^{ \prime } \)

        real(real64) :: pw02_k !! \( { k }^{ 2 } \)

        real(real64) :: sqrt_comp_k !! \( \sqrt{ { k }^{ \prime } } \)

        real(real64) :: pw01_eps !! auxiliary parameter \( \varepsilon \)

        real(real64) :: pw04_eps !! \( { \varepsilon }^{ 4 } \)



        call evaluate_modulus(k = k, pw02_k = pw02_k, comp_k = comp_k)

        call calculate_pw04_epsilon( &!
        &        pw02_k   =      pw02_k   , &!
        &        comp_k   =      comp_k   , &!
        &   sqrt_comp_k   = sqrt_comp_k   , &!
        &        pw01_eps =      pw01_eps , &!
        &        pw04_eps =      pw04_eps   &!
        )

        q = &!
            elliptic_nome_by_epsilon_09( &!
                pw01_eps = pw01_eps , &!
                pw04_eps = pw04_eps   &!
            )

    end function elliptic_nome_09_real64



    elemental function elliptic_nome_13_real64(k) result(q)
        !! Calculate the elliptic nome \( q \) using
        !! the Horner's method and
        !! the following polynomial:
        !! $$
        !! \begin{align*}
        !! \varepsilon &:= \frac{ 1 }{ 2 } \frac{ 1 - \sqrt{ { k }^{ \prime } } }{ 1 + \sqrt{ { k }^{ \prime } } } \\
        !! q(\varepsilon) &:=            \varepsilon
        !!             \\ & +        2 { \varepsilon }^{  5 }
        !!             \\ & +       15 { \varepsilon }^{  9 }
        !!             \\ & +      150 { \varepsilon }^{ 13 }
        !! \end{align*}
        !! $$
        !! @note
        !! - The elliptic modulus \( k \) should satisfy \( { k }^{ 2 } \le 1/2 \).
        !! @endnote

        real(real64), intent(in) :: k !! elliptic modulus \( k \)



        real(real64) :: q !! elliptic nome \( q \)



        real(real64) :: comp_k !! \( { k }^{ \prime } \)

        real(real64) :: pw02_k !! \( { k }^{ 2 } \)

        real(real64) :: sqrt_comp_k !! \( \sqrt{ { k }^{ \prime } } \)

        real(real64) :: pw01_eps !! auxiliary parameter \( \varepsilon \)

        real(real64) :: pw04_eps !! \( { \varepsilon }^{ 4 } \)



        call evaluate_modulus(k = k, pw02_k = pw02_k, comp_k = comp_k)

        call calculate_pw04_epsilon( &!
        &        pw02_k   =      pw02_k   , &!
        &        comp_k   =      comp_k   , &!
        &   sqrt_comp_k   = sqrt_comp_k   , &!
        &        pw01_eps =      pw01_eps , &!
        &        pw04_eps =      pw04_eps   &!
        )

        q = &!
            elliptic_nome_by_epsilon_13( &!
                pw01_eps = pw01_eps , &!
                pw04_eps = pw04_eps   &!
            )

    end function elliptic_nome_13_real64



    elemental function elliptic_nome_17_real64(k) result(q)
        !! Calculate the elliptic nome \( q \) using
        !! the Horner's method and
        !! the following polynomial:
        !! $$
        !! \begin{align*}
        !! \varepsilon &:= \frac{ 1 }{ 2 } \frac{ 1 - \sqrt{ { k }^{ \prime } } }{ 1 + \sqrt{ { k }^{ \prime } } } \\
        !! q(\varepsilon) &:=            \varepsilon
        !!             \\ & +        2 { \varepsilon }^{  5 }
        !!             \\ & +       15 { \varepsilon }^{  9 }
        !!             \\ & +      150 { \varepsilon }^{ 13 }
        !!             \\ & +     1707 { \varepsilon }^{ 17 }
        !! \end{align*}
        !! $$
        !! @note
        !! - The elliptic modulus \( k \) should satisfy \( { k }^{ 2 } \le 1/2 \).
        !! @endnote

        real(real64), intent(in) :: k !! elliptic modulus \( k \)



        real(real64) :: q !! elliptic nome \( q \)



        real(real64) :: comp_k !! \( { k }^{ \prime } \)

        real(real64) :: pw02_k !! \( { k }^{ 2 } \)

        real(real64) :: sqrt_comp_k !! \( \sqrt{ { k }^{ \prime } } \)

        real(real64) :: pw01_eps !! auxiliary parameter \( \varepsilon \)

        real(real64) :: pw04_eps !! \( { \varepsilon }^{ 4 } \)



        call evaluate_modulus(k = k, pw02_k = pw02_k, comp_k = comp_k)

        call calculate_pw04_epsilon( &!
        &        pw02_k   =      pw02_k   , &!
        &        comp_k   =      comp_k   , &!
        &   sqrt_comp_k   = sqrt_comp_k   , &!
        &        pw01_eps =      pw01_eps , &!
        &        pw04_eps =      pw04_eps   &!
        )

        q = &!
            elliptic_nome_by_epsilon_17( &!
                pw01_eps = pw01_eps , &!
                pw04_eps = pw04_eps   &!
            )

    end function elliptic_nome_17_real64



    elemental function elliptic_nome_21_real64(k) result(q)
        !! Calculate the elliptic nome \( q \) using
        !! the Horner's method and
        !! the following polynomial:
        !! $$
        !! \begin{align*}
        !! \varepsilon &:= \frac{ 1 }{ 2 } \frac{ 1 - \sqrt{ { k }^{ \prime } } }{ 1 + \sqrt{ { k }^{ \prime } } } \\
        !! q(\varepsilon) &:=            \varepsilon
        !!             \\ & +        2 { \varepsilon }^{  5 }
        !!             \\ & +       15 { \varepsilon }^{  9 }
        !!             \\ & +      150 { \varepsilon }^{ 13 }
        !!             \\ & +     1707 { \varepsilon }^{ 17 }
        !!             \\ & +    20910 { \varepsilon }^{ 21 }
        !! \end{align*}
        !! $$
        !! @note
        !! - The elliptic modulus \( k \) should satisfy \( { k }^{ 2 } \le 1/2 \).
        !! @endnote

        real(real64), intent(in) :: k !! elliptic modulus \( k \)



        real(real64) :: q !! elliptic nome \( q \)



        real(real64) :: comp_k !! \( { k }^{ \prime } \)

        real(real64) :: pw02_k !! \( { k }^{ 2 } \)

        real(real64) :: sqrt_comp_k !! \( \sqrt{ { k }^{ \prime } } \)

        real(real64) :: pw01_eps !! auxiliary parameter \( \varepsilon \)

        real(real64) :: pw04_eps !! \( { \varepsilon }^{ 4 } \)



        call evaluate_modulus(k = k, pw02_k = pw02_k, comp_k = comp_k)

        call calculate_pw04_epsilon( &!
        &        pw02_k   =      pw02_k   , &!
        &        comp_k   =      comp_k   , &!
        &   sqrt_comp_k   = sqrt_comp_k   , &!
        &        pw01_eps =      pw01_eps , &!
        &        pw04_eps =      pw04_eps   &!
        )

        q = &!
            elliptic_nome_by_epsilon_21( &!
                pw01_eps = pw01_eps , &!
                pw04_eps = pw04_eps   &!
            )

    end function elliptic_nome_21_real64



    elemental function elliptic_nome_25_real64(k) result(q)
        !! Calculate the elliptic nome \( q \) using
        !! the Horner's method and
        !! the following polynomial:
        !! $$
        !! \begin{align*}
        !! \varepsilon &:= \frac{ 1 }{ 2 } \frac{ 1 - \sqrt{ { k }^{ \prime } } }{ 1 + \sqrt{ { k }^{ \prime } } } \\
        !! q(\varepsilon) &:=            \varepsilon
        !!             \\ & +        2 { \varepsilon }^{  5 }
        !!             \\ & +       15 { \varepsilon }^{  9 }
        !!             \\ & +      150 { \varepsilon }^{ 13 }
        !!             \\ & +     1707 { \varepsilon }^{ 17 }
        !!             \\ & +    20910 { \varepsilon }^{ 21 }
        !!             \\ & +   268616 { \varepsilon }^{ 25 }
        !! \end{align*}
        !! $$
        !! @note
        !! - The elliptic modulus \( k \) should satisfy \( { k }^{ 2 } \le 1/2 \).
        !! @endnote

        real(real64), intent(in) :: k !! elliptic modulus \( k \)



        real(real64) :: q !! elliptic nome \( q \)



        real(real64) :: comp_k !! \( { k }^{ \prime } \)

        real(real64) :: pw02_k !! \( { k }^{ 2 } \)

        real(real64) :: sqrt_comp_k !! \( \sqrt{ { k }^{ \prime } } \)

        real(real64) :: pw01_eps !! auxiliary parameter \( \varepsilon \)

        real(real64) :: pw04_eps !! \( { \varepsilon }^{ 4 } \)



        call evaluate_modulus(k = k, pw02_k = pw02_k, comp_k = comp_k)

        call calculate_pw04_epsilon( &!
        &        pw02_k   =      pw02_k   , &!
        &        comp_k   =      comp_k   , &!
        &   sqrt_comp_k   = sqrt_comp_k   , &!
        &        pw01_eps =      pw01_eps , &!
        &        pw04_eps =      pw04_eps   &!
        )

        q = &!
            elliptic_nome_by_epsilon_25( &!
                pw01_eps = pw01_eps , &!
                pw04_eps = pw04_eps   &!
            )

    end function elliptic_nome_25_real64



    elemental function elliptic_nome_29_real64(k) result(q)
        !! Calculate the elliptic nome \( q \) using
        !! the Horner's method and
        !! the following polynomial:
        !! $$
        !! \begin{align*}
        !! \varepsilon &:= \frac{ 1 }{ 2 } \frac{ 1 - \sqrt{ { k }^{ \prime } } }{ 1 + \sqrt{ { k }^{ \prime } } } \\
        !! q(\varepsilon) &:=            \varepsilon
        !!             \\ & +        2 { \varepsilon }^{  5 }
        !!             \\ & +       15 { \varepsilon }^{  9 }
        !!             \\ & +      150 { \varepsilon }^{ 13 }
        !!             \\ & +     1707 { \varepsilon }^{ 17 }
        !!             \\ & +    20910 { \varepsilon }^{ 21 }
        !!             \\ & +   268616 { \varepsilon }^{ 25 }
        !!             \\ & +  3567400 { \varepsilon }^{ 29 }
        !! \end{align*}
        !! $$
        !! @note
        !! - The elliptic modulus \( k \) should satisfy \( { k }^{ 2 } \le 1/2 \).
        !! @endnote

        real(real64), intent(in) :: k !! elliptic modulus \( k \)



        real(real64) :: q !! elliptic nome \( q \)



        real(real64) :: comp_k !! \( { k }^{ \prime } \)

        real(real64) :: pw02_k !! \( { k }^{ 2 } \)

        real(real64) :: sqrt_comp_k !! \( \sqrt{ { k }^{ \prime } } \)

        real(real64) :: pw01_eps !! auxiliary parameter \( \varepsilon \)

        real(real64) :: pw04_eps !! \( { \varepsilon }^{ 4 } \)



        call evaluate_modulus(k = k, pw02_k = pw02_k, comp_k = comp_k)

        call calculate_pw04_epsilon( &!
        &        pw02_k   =      pw02_k   , &!
        &        comp_k   =      comp_k   , &!
        &   sqrt_comp_k   = sqrt_comp_k   , &!
        &        pw01_eps =      pw01_eps , &!
        &        pw04_eps =      pw04_eps   &!
        )

        q = &!
            elliptic_nome_by_epsilon_29( &!
                pw01_eps = pw01_eps , &!
                pw04_eps = pw04_eps   &!
            )

    end function elliptic_nome_29_real64



    elemental function elliptic_nome_33_real64(k) result(q)
        !! Calculate the elliptic nome \( q \) using
        !! the Horner's method and
        !! the following polynomial:
        !! $$
        !! \begin{align*}
        !! \varepsilon &:= \frac{ 1 }{ 2 } \frac{ 1 - \sqrt{ { k }^{ \prime } } }{ 1 + \sqrt{ { k }^{ \prime } } } \\
        !! q(\varepsilon) &:=            \varepsilon
        !!             \\ & +        2 { \varepsilon }^{  5 }
        !!             \\ & +       15 { \varepsilon }^{  9 }
        !!             \\ & +      150 { \varepsilon }^{ 13 }
        !!             \\ & +     1707 { \varepsilon }^{ 17 }
        !!             \\ & +    20910 { \varepsilon }^{ 21 }
        !!             \\ & +   268616 { \varepsilon }^{ 25 }
        !!             \\ & +  3567400 { \varepsilon }^{ 29 }
        !!             \\ & + 48555069 { \varepsilon }^{ 33 }
        !! \end{align*}
        !! $$
        !! @note
        !! - The elliptic modulus \( k \) should satisfy \( { k }^{ 2 } \le 1/2 \).
        !! @endnote

        real(real64), intent(in) :: k !! elliptic modulus \( k \)



        real(real64) :: q !! elliptic nome \( q \)



        real(real64) :: comp_k !! \( { k }^{ \prime } \)

        real(real64) :: pw02_k !! \( { k }^{ 2 } \)

        real(real64) :: sqrt_comp_k !! \( \sqrt{ { k }^{ \prime } } \)

        real(real64) :: pw01_eps !! auxiliary parameter \( \varepsilon \)

        real(real64) :: pw04_eps !! \( { \varepsilon }^{ 4 } \)



        call evaluate_modulus(k = k, pw02_k = pw02_k, comp_k = comp_k)

        call calculate_pw04_epsilon( &!
        &        pw02_k   =      pw02_k   , &!
        &        comp_k   =      comp_k   , &!
        &   sqrt_comp_k   = sqrt_comp_k   , &!
        &        pw01_eps =      pw01_eps , &!
        &        pw04_eps =      pw04_eps   &!
        )

        q = &!
            elliptic_nome_by_epsilon_33( &!
                pw01_eps = pw01_eps , &!
                pw04_eps = pw04_eps   &!
            )

    end function elliptic_nome_33_real64



    elemental function elliptic_nome_auto_real64(k) result(q)
        !! Calculate the elliptic nome \( q \) 
        !! for the given elliptic modulus \( k \)
        !! @note
        !! - The elliptic modulus \( k \) should satisfy \( { k }^{ 2 } \le 1/2 \).
        !! - If the calculation does not converge, it returns NaN.
        !! @endnote

        real(real64), intent(in) :: k !! elliptic modulus \( k \)



        real(real64) :: q !! elliptic nome \( q \)



        real(real64) :: comp_k !! \( { k }^{ \prime } \)

        real(real64) :: pw01_eps !! auxiliary parameter \( \varepsilon \)

        real(real64) :: pw02_k !! \( { k }^{ 2 } \)

        real(real64) :: pw04_eps !! \( { \varepsilon }^{ 4 } \)

        real(real64) :: q_ref

        real(real64) :: sqrt_comp_k !! \( \sqrt{ { k }^{ \prime } } \)



        call evaluate_modulus(k = k, pw02_k = pw02_k, comp_k = comp_k)

        call calculate_pw04_epsilon( &!
        &        pw02_k   =      pw02_k   , &!
        &        comp_k   =      comp_k   , &!
        &   sqrt_comp_k   = sqrt_comp_k   , &!
        &        pw01_eps =      pw01_eps , &!
        &        pw04_eps =      pw04_eps   &!
        )



        q_ref = pw01_eps

        q = &!
            elliptic_nome_by_epsilon_05( &!
                pw01_eps = pw01_eps , &!
                pw04_eps = pw04_eps   &!
            )

        if ( abs(q - q_ref) .lt. q_err ) return



        q_ref = q

        q = &!
            elliptic_nome_by_epsilon_09( &!
                pw01_eps = pw01_eps , &!
                pw04_eps = pw04_eps   &!
            )

        if ( abs(q - q_ref) .lt. q_err ) return



        q_ref = q

        q = &!
            elliptic_nome_by_epsilon_13( &!
                pw01_eps = pw01_eps , &!
                pw04_eps = pw04_eps   &!
            )

        if ( abs(q - q_ref) .lt. q_err ) return



        q_ref = q

        q = &!
            elliptic_nome_by_epsilon_17( &!
                pw01_eps = pw01_eps , &!
                pw04_eps = pw04_eps   &!
            )

        if ( abs(q - q_ref) .lt. q_err ) return



        q_ref = q

        q = &!
            elliptic_nome_by_epsilon_21( &!
                pw01_eps = pw01_eps , &!
                pw04_eps = pw04_eps   &!
            )

        if ( abs(q - q_ref) .lt. q_err ) return



        q_ref = q

        q = &!
            elliptic_nome_by_epsilon_25( &!
                pw01_eps = pw01_eps , &!
                pw04_eps = pw04_eps   &!
            )

        if ( abs(q - q_ref) .lt. q_err ) return



        q_ref = q

        q = &!
            elliptic_nome_by_epsilon_29( &!
                pw01_eps = pw01_eps , &!
                pw04_eps = pw04_eps   &!
            )

        if ( abs(q - q_ref) .lt. q_err ) return



        q_ref = q

        q = &!
            elliptic_nome_by_epsilon_33( &!
                pw01_eps = pw01_eps , &!
                pw04_eps = pw04_eps   &!
            )

        if ( abs(q - q_ref) .lt. q_err ) return




        q = ieee_value(q, ieee_quiet_nan)

    end function elliptic_nome_auto_real64



    elemental function elliptic_nome_by_epsilon_05(pw01_eps, pw04_eps) result(q)
        !! Calculate the elliptic nome \( q \) using
        !! the Horner's method and
        !! the following polynomial:
        !! $$
        !! \begin{align*}
        !! q(\varepsilon) &:=            \varepsilon
        !!             \\ & +        2 { \varepsilon }^{  5 }
        !! \end{align*}
        !! $$

        real(real64), intent(in) :: pw01_eps !! auxiliary parameter \( \varepsilon \)

        real(real64), intent(in) :: pw04_eps !! \( { \varepsilon }^{ 4 } \)



        real(real64) :: q !! elliptic nome \( q \)



        q = &!
            elliptic_nome_by_epsilon_05_horner( &!
                pw01_eps = pw01_eps , &!
                pw04_eps = pw04_eps , &!
                pre_step = c_05       &!
            )

    end function elliptic_nome_by_epsilon_05



    elemental function elliptic_nome_by_epsilon_09(pw01_eps, pw04_eps) result(q)
        !! Calculate the elliptic nome \( q \) using
        !! the Horner's method and
        !! the following polynomial:
        !! $$
        !! \begin{align*}
        !! q(\varepsilon) &:=            \varepsilon
        !!             \\ & +        2 { \varepsilon }^{  5 }
        !!             \\ & +       15 { \varepsilon }^{  9 }
        !! \end{align*}
        !! $$

        real(real64), intent(in) :: pw01_eps !! auxiliary parameter \( \varepsilon \)

        real(real64), intent(in) :: pw04_eps !! \( { \varepsilon }^{ 4 } \)



        real(real64) :: q !! elliptic nome \( q \)



        q = &!
            elliptic_nome_by_epsilon_09_horner( &!
                pw01_eps = pw01_eps , &!
                pw04_eps = pw04_eps , &!
                pre_step = c_09       &!
            )

    end function elliptic_nome_by_epsilon_09



    elemental function elliptic_nome_by_epsilon_13(pw01_eps, pw04_eps) result(q)
        !! Calculate the elliptic nome \( q \) using
        !! the Horner's method and
        !! the following polynomial:
        !! $$
        !! \begin{align*}
        !! q(\varepsilon) &:=            \varepsilon
        !!             \\ & +        2 { \varepsilon }^{  5 }
        !!             \\ & +       15 { \varepsilon }^{  9 }
        !!             \\ & +      150 { \varepsilon }^{ 13 }
        !! \end{align*}
        !! $$

        real(real64), intent(in) :: pw01_eps !! auxiliary parameter \( \varepsilon \)

        real(real64), intent(in) :: pw04_eps !! \( { \varepsilon }^{ 4 } \)



        real(real64) :: q !! elliptic nome \( q \)



        q = &!
            elliptic_nome_by_epsilon_13_horner( &!
                pw01_eps = pw01_eps , &!
                pw04_eps = pw04_eps , &!
                pre_step = c_13       &!
            )

    end function elliptic_nome_by_epsilon_13



    elemental function elliptic_nome_by_epsilon_17(pw01_eps, pw04_eps) result(q)
        !! Calculate the elliptic nome \( q \) using
        !! the Horner's method and
        !! the following polynomial:
        !! $$
        !! \begin{align*}
        !! q(\varepsilon) &:=            \varepsilon
        !!             \\ & +        2 { \varepsilon }^{  5 }
        !!             \\ & +       15 { \varepsilon }^{  9 }
        !!             \\ & +      150 { \varepsilon }^{ 13 }
        !!             \\ & +     1707 { \varepsilon }^{ 17 }
        !! \end{align*}
        !! $$

        real(real64), intent(in) :: pw01_eps !! auxiliary parameter \( \varepsilon \)

        real(real64), intent(in) :: pw04_eps !! \( { \varepsilon }^{ 4 } \)



        real(real64) :: q !! elliptic nome \( q \)



        q = &!
            elliptic_nome_by_epsilon_17_horner( &!
                pw01_eps = pw01_eps , &!
                pw04_eps = pw04_eps , &!
                pre_step = c_17       &!
            )

    end function elliptic_nome_by_epsilon_17



    elemental function elliptic_nome_by_epsilon_21(pw01_eps, pw04_eps) result(q)
        !! Calculate the elliptic nome \( q \) using
        !! the Horner's method and
        !! the following polynomial:
        !! $$
        !! \begin{align*}
        !! q(\varepsilon) &:=            \varepsilon
        !!             \\ & +        2 { \varepsilon }^{  5 }
        !!             \\ & +       15 { \varepsilon }^{  9 }
        !!             \\ & +      150 { \varepsilon }^{ 13 }
        !!             \\ & +     1707 { \varepsilon }^{ 17 }
        !!             \\ & +    20910 { \varepsilon }^{ 21 }
        !! \end{align*}
        !! $$

        real(real64), intent(in) :: pw01_eps !! auxiliary parameter \( \varepsilon \)

        real(real64), intent(in) :: pw04_eps !! \( { \varepsilon }^{ 4 } \)



        real(real64) :: q !! elliptic nome \( q \)



        q = &!
            elliptic_nome_by_epsilon_21_horner( &!
                pw01_eps = pw01_eps , &!
                pw04_eps = pw04_eps , &!
                pre_step = c_21       &!
            )

    end function elliptic_nome_by_epsilon_21



    elemental function elliptic_nome_by_epsilon_25(pw01_eps, pw04_eps) result(q)
        !! Calculate the elliptic nome \( q \) using
        !! the Horner's method and
        !! the following polynomial:
        !! $$
        !! \begin{align*}
        !! q(\varepsilon) &:=            \varepsilon
        !!             \\ & +        2 { \varepsilon }^{  5 }
        !!             \\ & +       15 { \varepsilon }^{  9 }
        !!             \\ & +      150 { \varepsilon }^{ 13 }
        !!             \\ & +     1707 { \varepsilon }^{ 17 }
        !!             \\ & +    20910 { \varepsilon }^{ 21 }
        !!             \\ & +   268616 { \varepsilon }^{ 25 }
        !! \end{align*}
        !! $$

        real(real64), intent(in) :: pw01_eps !! auxiliary parameter \( \varepsilon \)

        real(real64), intent(in) :: pw04_eps !! \( { \varepsilon }^{ 4 } \)



        real(real64) :: q !! elliptic nome \( q \)



        q = &!
            elliptic_nome_by_epsilon_25_horner( &!
                pw01_eps = pw01_eps , &!
                pw04_eps = pw04_eps , &!
                pre_step = c_25       &!
            )

    end function elliptic_nome_by_epsilon_25



    elemental function elliptic_nome_by_epsilon_29(pw01_eps, pw04_eps) result(q)
        !! Calculate the elliptic nome \( q \) using
        !! the Horner's method and
        !! the following polynomial:
        !! $$
        !! \begin{align*}
        !! q(\varepsilon) &:=            \varepsilon
        !!             \\ & +        2 { \varepsilon }^{  5 }
        !!             \\ & +       15 { \varepsilon }^{  9 }
        !!             \\ & +      150 { \varepsilon }^{ 13 }
        !!             \\ & +     1707 { \varepsilon }^{ 17 }
        !!             \\ & +    20910 { \varepsilon }^{ 21 }
        !!             \\ & +   268616 { \varepsilon }^{ 25 }
        !!             \\ & +  3567400 { \varepsilon }^{ 29 }
        !! \end{align*}
        !! $$

        real(real64), intent(in) :: pw01_eps !! auxiliary parameter \( \varepsilon \)

        real(real64), intent(in) :: pw04_eps !! \( { \varepsilon }^{ 4 } \)



        real(real64) :: q !! elliptic nome \( q \)



        q = &!
            elliptic_nome_by_epsilon_29_horner( &!
                pw01_eps = pw01_eps , &!
                pw04_eps = pw04_eps , &!
                pre_step = c_29       &!
            )

    end function elliptic_nome_by_epsilon_29



    elemental function elliptic_nome_by_epsilon_33(pw01_eps, pw04_eps) result(q)
        !! Calculate the elliptic nome \( q \) using
        !! the Horner's method and
        !! the following polynomial:
        !! $$
        !! \begin{align*}
        !! q(\varepsilon) &:=            \varepsilon
        !!             \\ & +        2 { \varepsilon }^{  5 }
        !!             \\ & +       15 { \varepsilon }^{  9 }
        !!             \\ & +      150 { \varepsilon }^{ 13 }
        !!             \\ & +     1707 { \varepsilon }^{ 17 }
        !!             \\ & +    20910 { \varepsilon }^{ 21 }
        !!             \\ & +   268616 { \varepsilon }^{ 25 }
        !!             \\ & +  3567400 { \varepsilon }^{ 29 }
        !!             \\ & + 48555069 { \varepsilon }^{ 33 }
        !! \end{align*}
        !! $$

        real(real64), intent(in) :: pw01_eps !! auxiliary parameter \( \varepsilon \)

        real(real64), intent(in) :: pw04_eps !! \( { \varepsilon }^{ 4 } \)



        real(real64) :: q !! elliptic nome \( q \)



        q = &!
            elliptic_nome_by_epsilon_33_horner( &!
                pw01_eps = pw01_eps , &!
                pw04_eps = pw04_eps , &!
                pre_step = c_33       &!
            )

    end function elliptic_nome_by_epsilon_33





    elemental function elliptic_nome_by_epsilon_05_horner(pw01_eps, pw04_eps, pre_step) result(q)
        !! Calculate the following for the given \( \varepsilon \) and \( { \varepsilon }^{ 4 } \):
        !! $$ \varepsilon ( 1 + { \varepsilon }^{ 4 } \cdot \texttt{pre_step} ) $$

        real(real64), intent(in) :: pw01_eps !! auxiliary parameter \( \varepsilon \)

        real(real64), intent(in) :: pw04_eps !! \( { \varepsilon }^{ 4 } \)

        real(real64), intent(in) :: pre_step



        real(real64) :: q !! elliptic nome \( q \)



        q = pw01_eps * (c_01 + pre_step * pw04_eps)

    end function elliptic_nome_by_epsilon_05_horner



    elemental function elliptic_nome_by_epsilon_09_horner(pw01_eps, pw04_eps, pre_step) result(q)
        !! Calculate the following for the given \( \varepsilon \) and \( { \varepsilon }^{ 4 } \):
        !! $$
        !! \begin{aligned}
        !! \varepsilon \cdot {}
        !! & ( 1 + { \varepsilon }^{ 4 } \cdot \\
        !! & ( 2 + { \varepsilon }^{ 4 } \cdot \texttt{pre_step} ) \\
        !! & ) \\
        !! \end{aligned}
        !! $$

        real(real64), intent(in) :: pw01_eps !! auxiliary parameter \( \varepsilon \)

        real(real64), intent(in) :: pw04_eps !! \( { \varepsilon }^{ 4 } \)

        real(real64), intent(in) :: pre_step



        real(real64) :: q !! elliptic nome \( q \)



        q = &!
            elliptic_nome_by_epsilon_05_horner(&!
                pw01_eps = pw01_eps                   , &!
                pw04_eps = pw04_eps                   , &!
                pre_step = pw04_eps * pre_step + c_05   &!
            )

    end function elliptic_nome_by_epsilon_09_horner



    elemental function elliptic_nome_by_epsilon_13_horner(pw01_eps, pw04_eps, pre_step) result(q)
        !! Calculate the following for the given \( \varepsilon \) and \( { \varepsilon }^{ 4 } \):
        !! $$
        !! \begin{aligned}
        !! \varepsilon \cdot {}
        !! & (  1 + { \varepsilon }^{ 4 } \cdot \\
        !! & (  2 + { \varepsilon }^{ 4 } \cdot \\
        !! & ( 15 + { \varepsilon }^{ 4 } \cdot \texttt{pre_step} ) \\
        !! & ) \\
        !! & ) \\
        !! \end{aligned}
        !! $$

        real(real64), intent(in) :: pw01_eps !! auxiliary parameter \( \varepsilon \)

        real(real64), intent(in) :: pw04_eps !! \( { \varepsilon }^{ 4 } \)

        real(real64), intent(in) :: pre_step



        real(real64) :: q !! elliptic nome \( q \)



        q = &!
            elliptic_nome_by_epsilon_09_horner(&!
                pw01_eps = pw01_eps                   , &!
                pw04_eps = pw04_eps                   , &!
                pre_step = pw04_eps * pre_step + c_09   &!
            )

    end function elliptic_nome_by_epsilon_13_horner



    elemental function elliptic_nome_by_epsilon_17_horner(pw01_eps, pw04_eps, pre_step) result(q)
        !! Calculate the following for the given \( \varepsilon \) and \( { \varepsilon }^{ 4 } \):
        !! $$
        !! \begin{aligned}
        !! \varepsilon \cdot {}
        !! & (   1 + { \varepsilon }^{ 4 } \cdot \\
        !! & (   2 + { \varepsilon }^{ 4 } \cdot \\
        !! & (  15 + { \varepsilon }^{ 4 } \cdot \\
        !! & ( 150 + { \varepsilon }^{ 4 } \cdot \texttt{pre_step} ) \\
        !! & ) \\
        !! & ) \\
        !! & ) \\
        !! \end{aligned}
        !! $$

        real(real64), intent(in) :: pw01_eps !! auxiliary parameter \( \varepsilon \)

        real(real64), intent(in) :: pw04_eps !! \( { \varepsilon }^{ 4 } \)

        real(real64), intent(in) :: pre_step



        real(real64) :: q !! elliptic nome \( q \)



        q = &!
            elliptic_nome_by_epsilon_13_horner(&!
                pw01_eps = pw01_eps                   , &!
                pw04_eps = pw04_eps                   , &!
                pre_step = pw04_eps * pre_step + c_13   &!
            )

    end function elliptic_nome_by_epsilon_17_horner



    elemental function elliptic_nome_by_epsilon_21_horner(pw01_eps, pw04_eps, pre_step) result(q)
        !! Calculate the following for the given \( \varepsilon \) and \( { \varepsilon }^{ 4 } \):
        !! $$
        !! \begin{aligned}
        !! \varepsilon \cdot {}
        !! & (    1 + { \varepsilon }^{ 4 } \cdot \\
        !! & (    2 + { \varepsilon }^{ 4 } \cdot \\
        !! & (   15 + { \varepsilon }^{ 4 } \cdot \\
        !! & (  150 + { \varepsilon }^{ 4 } \cdot \\
        !! & ( 1707 + { \varepsilon }^{ 4 } \cdot \texttt{pre_step} ) \\
        !! & ) \\
        !! & ) \\
        !! & ) \\
        !! & ) \\
        !! \end{aligned}
        !! $$

        real(real64), intent(in) :: pw01_eps !! auxiliary parameter \( \varepsilon \)

        real(real64), intent(in) :: pw04_eps !! \( { \varepsilon }^{ 4 } \)

        real(real64), intent(in) :: pre_step



        real(real64) :: q !! elliptic nome \( q \)



        q = &!
            elliptic_nome_by_epsilon_17_horner(&!
                pw01_eps = pw01_eps                   , &!
                pw04_eps = pw04_eps                   , &!
                pre_step = pw04_eps * pre_step + c_17   &!
            )

    end function elliptic_nome_by_epsilon_21_horner



    elemental function elliptic_nome_by_epsilon_25_horner(pw01_eps, pw04_eps, pre_step) result(q)
        !! Calculate the following for the given \( \varepsilon \) and \( { \varepsilon }^{ 4 } \):
        !! $$
        !! \begin{aligned}
        !! \varepsilon \cdot {}
        !! & (     1 + { \varepsilon }^{ 4 } \cdot \\
        !! & (     2 + { \varepsilon }^{ 4 } \cdot \\
        !! & (    15 + { \varepsilon }^{ 4 } \cdot \\
        !! & (   150 + { \varepsilon }^{ 4 } \cdot \\
        !! & (  1707 + { \varepsilon }^{ 4 } \cdot \\
        !! & ( 20910 + { \varepsilon }^{ 4 } \cdot \texttt{pre_step} ) \\
        !! & ) \\
        !! & ) \\
        !! & ) \\
        !! & ) \\
        !! & ) \\
        !! \end{aligned}
        !! $$

        real(real64), intent(in) :: pw01_eps !! auxiliary parameter \( \varepsilon \)

        real(real64), intent(in) :: pw04_eps !! \( { \varepsilon }^{ 4 } \)

        real(real64), intent(in) :: pre_step



        real(real64) :: q !! elliptic nome \( q \)



        q = &!
            elliptic_nome_by_epsilon_21_horner(&!
                pw01_eps = pw01_eps                   , &!
                pw04_eps = pw04_eps                   , &!
                pre_step = pw04_eps * pre_step + c_21   &!
            )

    end function elliptic_nome_by_epsilon_25_horner



    elemental function elliptic_nome_by_epsilon_29_horner(pw01_eps, pw04_eps, pre_step) result(q)
        !! Calculate the following for the given \( \varepsilon \) and \( { \varepsilon }^{ 4 } \):
        !! $$
        !! \begin{aligned}
        !! \varepsilon \cdot {}
        !! & (      1 + { \varepsilon }^{ 4 } \cdot \\
        !! & (      2 + { \varepsilon }^{ 4 } \cdot \\
        !! & (     15 + { \varepsilon }^{ 4 } \cdot \\
        !! & (    150 + { \varepsilon }^{ 4 } \cdot \\
        !! & (   1707 + { \varepsilon }^{ 4 } \cdot \\
        !! & (  20910 + { \varepsilon }^{ 4 } \cdot \\
        !! & ( 268616 + { \varepsilon }^{ 4 } \cdot \texttt{pre_step} ) \\
        !! & ) \\
        !! & ) \\
        !! & ) \\
        !! & ) \\
        !! & ) \\
        !! & ) \\
        !! \end{aligned}
        !! $$

        real(real64), intent(in) :: pw01_eps !! auxiliary parameter \( \varepsilon \)

        real(real64), intent(in) :: pw04_eps !! \( { \varepsilon }^{ 4 } \)

        real(real64), intent(in) :: pre_step



        real(real64) :: q !! elliptic nome \( q \)



        q = &!
            elliptic_nome_by_epsilon_25_horner(&!
                pw01_eps = pw01_eps                   , &!
                pw04_eps = pw04_eps                   , &!
                pre_step = pw04_eps * pre_step + c_25   &!
            )

    end function elliptic_nome_by_epsilon_29_horner



    elemental function elliptic_nome_by_epsilon_33_horner(pw01_eps, pw04_eps, pre_step) result(q)
        !! Calculate the following for the given \( \varepsilon \) and \( { \varepsilon }^{ 4 } \):
        !! $$
        !! \begin{aligned}
        !! \varepsilon \cdot {}
        !! & (       1 + { \varepsilon }^{ 4 } \cdot \\
        !! & (       2 + { \varepsilon }^{ 4 } \cdot \\
        !! & (      15 + { \varepsilon }^{ 4 } \cdot \\
        !! & (     150 + { \varepsilon }^{ 4 } \cdot \\
        !! & (    1707 + { \varepsilon }^{ 4 } \cdot \\
        !! & (   20910 + { \varepsilon }^{ 4 } \cdot \\
        !! & (  268616 + { \varepsilon }^{ 4 } \cdot \\
        !! & ( 3567400 + { \varepsilon }^{ 4 } \cdot \texttt{pre_step} ) \\
        !! & ) \\
        !! & ) \\
        !! & ) \\
        !! & ) \\
        !! & ) \\
        !! & ) \\
        !! & ) \\
        !! \end{aligned}
        !! $$

        real(real64), intent(in) :: pw01_eps !! auxiliary parameter \( \varepsilon \)

        real(real64), intent(in) :: pw04_eps !! \( { \varepsilon }^{ 4 } \)

        real(real64), intent(in) :: pre_step



        real(real64) :: q !! elliptic nome \( q \)



        q = &!
            elliptic_nome_by_epsilon_29_horner(&!
                pw01_eps = pw01_eps                   , &!
                pw04_eps = pw04_eps                   , &!
                pre_step = pw04_eps * pre_step + c_29   &!
            )

    end function elliptic_nome_by_epsilon_33_horner



    elemental subroutine calculate_pw01_epsilon(pw02_k, comp_k, sqrt_comp_k, pw01_eps)
        !! calculate the auxiliary parameter \( \varepsilon \)
        !! for the given elliptic modulus \( k \)

        real(real64), intent(in) :: pw02_k !! \( { k }^{ 2 } \)

        real(real64), intent(in) :: comp_k !! \( { k }^{ \prime } \)

        real(real64), intent(out) :: sqrt_comp_k !! \( \sqrt{ { k }^{ \prime } } \)

        real(real64), intent(out) :: pw01_eps !! auxiliary parameter \( \varepsilon \)



        real(real64) :: pls1_comp_k !! \( 1 + { k }^{ \prime } \)



        pls1_comp_k =      comp_k  + 1.0_real64
        sqrt_comp_k = sqrt(comp_k)

        pw01_eps = 0.5_real64 * pw02_k &!
        &        / ( pls1_comp_k * ( pls1_comp_k + sqrt_comp_k + sqrt_comp_k ) )

    end subroutine calculate_pw01_epsilon



    elemental subroutine calculate_pw04_epsilon(pw02_k, comp_k, sqrt_comp_k, pw01_eps, pw04_eps)
        !! calculate the auxiliary parameter \( \varepsilon \) and \( { \varepsilon }^{ 4 } \)
        !! for the given elliptic modulus \( k \)

        real(real64), intent(in) :: pw02_k !! \( { k }^{ 2 } \)

        real(real64), intent(in) :: comp_k !! \( { k }^{ \prime } \)

        real(real64), intent(out) :: sqrt_comp_k !! \( \sqrt{ { k }^{ \prime } } \)

        real(real64), intent(out) :: pw01_eps !! auxiliary parameter \( \varepsilon \)

        real(real64), intent(out) :: pw04_eps !! \( { \varepsilon }^{ 4 } \)



        real(real64) :: pw02_eps



        call calculate_pw01_epsilon( &!
        &        pw02_k   = pw02_k      , &!
        &        comp_k   = comp_k      , &!
        &   sqrt_comp_k   = sqrt_comp_k , &!
        &        pw01_eps = pw01_eps      &!
        )

        pw02_eps = pw01_eps * pw01_eps
        pw04_eps = pw02_eps * pw02_eps

    end subroutine calculate_pw04_epsilon



    elemental subroutine evaluate_modulus(k, pw02_k, comp_k)
        !! Calculate \( { k }^{ 2 } \) 
        !! and \( { k }^{ \prime } := \sqrt{ 1 - { k }^{ 2 } } \) 
        !! for the given elliptic modulus \( k \)

        real(real64), intent(in) :: k !! elliptic modulus \( k \)

        real(real64), intent(out) :: pw02_k !! \( { k }^{ 2 } \)

        real(real64), intent(out) :: comp_k !! \( { k }^{ \prime } \)



        pw02_k = k * k
        comp_k = sqrt(1.0_real64 - pw02_k)

    end subroutine evaluate_modulus

end module elliptic_nome_fortran_real64
