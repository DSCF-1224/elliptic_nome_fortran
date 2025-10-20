module elliptic_nome_fortran
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

    use, non_intrinsic :: elliptic_nome_fortran_real32
    use, non_intrinsic :: elliptic_nome_fortran_real64
    use, non_intrinsic :: elliptic_nome_fortran_real128



    implicit none



    private



    public :: elliptic_nome_01
    public :: elliptic_nome_05
    public :: elliptic_nome_09
    public :: elliptic_nome_13
    public :: elliptic_nome_17
    public :: elliptic_nome_21
    public :: elliptic_nome_25
    public :: elliptic_nome_29
    public :: elliptic_nome_33
    public :: elliptic_nome_auto



    interface elliptic_nome_01
        !! Calculate the elliptic nome \( q \) using
        !! the Horner's method and
        !! the following polynomial:
        !! $$
        !! \begin{align*}
        !! \varepsilon &:= \frac{ 1 }{ 2 } \frac{ 1 - \sqrt{ { k }^{ \prime } } }{ 1 + \sqrt{ { k }^{ \prime } } } \\
        !! q(\varepsilon) &:=     \varepsilon
        !! \end{align*}
        !! $$

        module procedure :: elliptic_nome_01_real32
        module procedure :: elliptic_nome_01_real64
        module procedure :: elliptic_nome_01_real128

    end interface elliptic_nome_01



    interface elliptic_nome_05
        !! Calculate the elliptic nome \( q \) using
        !! the Horner's method and
        !! the following polynomial:
        !! $$
        !! \begin{align*}
        !! \varepsilon &:= \frac{ 1 }{ 2 } \frac{ 1 - \sqrt{ { k }^{ \prime } } }{ 1 + \sqrt{ { k }^{ \prime } } } \\
        !! q(\varepsilon) &:=     \varepsilon
        !!             \\ & + 2 { \varepsilon }^{  5 }
        !! \end{align*}
        !! $$

        module procedure :: elliptic_nome_05_real32
        module procedure :: elliptic_nome_05_real64
        module procedure :: elliptic_nome_05_real128

    end interface elliptic_nome_05



    interface elliptic_nome_09
        !! Calculate the elliptic nome \( q \) using
        !! the Horner's method and
        !! the following polynomial:
        !! $$
        !! \begin{align*}
        !! \varepsilon &:= \frac{ 1 }{ 2 } \frac{ 1 - \sqrt{ { k }^{ \prime } } }{ 1 + \sqrt{ { k }^{ \prime } } } \\
        !! q(\varepsilon) &:=      \varepsilon
        !!             \\ & +  2 { \varepsilon }^{  5 }
        !!             \\ & + 15 { \varepsilon }^{  9 }
        !! \end{align*}
        !! $$

        module procedure :: elliptic_nome_09_real32
        module procedure :: elliptic_nome_09_real64
        module procedure :: elliptic_nome_09_real128

    end interface elliptic_nome_09



    interface elliptic_nome_13
        !! Calculate the elliptic nome \( q \) using
        !! the Horner's method and
        !! the following polynomial:
        !! $$
        !! \begin{align*}
        !! \varepsilon &:= \frac{ 1 }{ 2 } \frac{ 1 - \sqrt{ { k }^{ \prime } } }{ 1 + \sqrt{ { k }^{ \prime } } } \\
        !! q(\varepsilon) &:=       \varepsilon
        !!             \\ & +   2 { \varepsilon }^{  5 }
        !!             \\ & +  15 { \varepsilon }^{  9 }
        !!             \\ & + 150 { \varepsilon }^{ 13 }
        !! \end{align*}
        !! $$

        module procedure :: elliptic_nome_13_real32
        module procedure :: elliptic_nome_13_real64
        module procedure :: elliptic_nome_13_real128

    end interface elliptic_nome_13



    interface elliptic_nome_17
        !! Calculate the elliptic nome \( q \) using
        !! the Horner's method and
        !! the following polynomial:
        !! $$
        !! \begin{align*}
        !! \varepsilon &:= \frac{ 1 }{ 2 } \frac{ 1 - \sqrt{ { k }^{ \prime } } }{ 1 + \sqrt{ { k }^{ \prime } } } \\
        !! q(\varepsilon) &:=        \varepsilon
        !!             \\ & +    2 { \varepsilon }^{  5 }
        !!             \\ & +   15 { \varepsilon }^{  9 }
        !!             \\ & +  150 { \varepsilon }^{ 13 }
        !!             \\ & + 1707 { \varepsilon }^{ 17 }
        !! \end{align*}
        !! $$

        module procedure :: elliptic_nome_17_real32
        module procedure :: elliptic_nome_17_real64
        module procedure :: elliptic_nome_17_real128

    end interface elliptic_nome_17



    interface elliptic_nome_21
        !! Calculate the elliptic nome \( q \) using
        !! the Horner's method and
        !! the following polynomial:
        !! $$
        !! \begin{align*}
        !! \varepsilon &:= \frac{ 1 }{ 2 } \frac{ 1 - \sqrt{ { k }^{ \prime } } }{ 1 + \sqrt{ { k }^{ \prime } } } \\
        !! q(\varepsilon) &:=         \varepsilon
        !!             \\ & +     2 { \varepsilon }^{  5 }
        !!             \\ & +    15 { \varepsilon }^{  9 }
        !!             \\ & +   150 { \varepsilon }^{ 13 }
        !!             \\ & +  1707 { \varepsilon }^{ 17 }
        !!             \\ & + 20910 { \varepsilon }^{ 21 }
        !! \end{align*}
        !! $$

        module procedure :: elliptic_nome_21_real32
        module procedure :: elliptic_nome_21_real64
        module procedure :: elliptic_nome_21_real128

    end interface elliptic_nome_21



    interface elliptic_nome_25
        !! Calculate the elliptic nome \( q \) using
        !! the Horner's method and
        !! the following polynomial:
        !! $$
        !! \begin{align*}
        !! \varepsilon &:= \frac{ 1 }{ 2 } \frac{ 1 - \sqrt{ { k }^{ \prime } } }{ 1 + \sqrt{ { k }^{ \prime } } } \\
        !! q(\varepsilon) &:=          \varepsilon
        !!             \\ & +      2 { \varepsilon }^{  5 }
        !!             \\ & +     15 { \varepsilon }^{  9 }
        !!             \\ & +    150 { \varepsilon }^{ 13 }
        !!             \\ & +   1707 { \varepsilon }^{ 17 }
        !!             \\ & +  20910 { \varepsilon }^{ 21 }
        !!             \\ & + 268616 { \varepsilon }^{ 25 }
        !! \end{align*}
        !! $$

        module procedure :: elliptic_nome_25_real32
        module procedure :: elliptic_nome_25_real64
        module procedure :: elliptic_nome_25_real128

    end interface elliptic_nome_25



    interface elliptic_nome_29
        !! Calculate the elliptic nome \( q \) using
        !! the Horner's method and
        !! the following polynomial:
        !! $$
        !! \begin{align*}
        !! \varepsilon &:= \frac{ 1 }{ 2 } \frac{ 1 - \sqrt{ { k }^{ \prime } } }{ 1 + \sqrt{ { k }^{ \prime } } } \\
        !! q(\varepsilon) &:=           \varepsilon
        !!             \\ & +       2 { \varepsilon }^{  5 }
        !!             \\ & +      15 { \varepsilon }^{  9 }
        !!             \\ & +     150 { \varepsilon }^{ 13 }
        !!             \\ & +    1707 { \varepsilon }^{ 17 }
        !!             \\ & +   20910 { \varepsilon }^{ 21 }
        !!             \\ & +  268616 { \varepsilon }^{ 25 }
        !!             \\ & + 3567400 { \varepsilon }^{ 29 }
        !! \end{align*}
        !! $$

        module procedure :: elliptic_nome_29_real32
        module procedure :: elliptic_nome_29_real64
        module procedure :: elliptic_nome_29_real128

    end interface elliptic_nome_29



    interface elliptic_nome_33
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

        module procedure :: elliptic_nome_33_real32
        module procedure :: elliptic_nome_33_real64
        module procedure :: elliptic_nome_33_real128

    end interface elliptic_nome_33



    interface elliptic_nome_auto
        module procedure :: elliptic_nome_auto_real32
        module procedure :: elliptic_nome_auto_real64
        module procedure :: elliptic_nome_auto_real128
    end interface elliptic_nome_auto

end module elliptic_nome_fortran
