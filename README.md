# elliptic_nome_fortran

Compute nome $ q $ for elliptic integrals corresponding to modulus $ k $

$$
\begin{align*}
K (k) &\coloneqq \int_{ 0 }^{ \pi / 2 } \frac{ d \theta }{ \sqrt{ 1 - { k }^{ 2 } \sin^{ 2 } \theta } } , &
{ k }^{ \prime } &\coloneqq \sqrt{ 1 - { k }^{ 2 } } , &
{ K }^{ \prime } (k) &\coloneqq K ( { k }^{ \prime } ) , &
q &:= e^{ - \pi { K }^{ \prime } (k) / K (k) }
\end{align*}
$$

## Reference

- 山内二郎, 宇野利雄, 一松信 共編  
  電子計算機のための数値計算法 3  
  培風館, 1972.  
  数理科学シリーズ ; 5  
  [NDLサーチ](https://ndlsearch.ndl.go.jp/books/R100000039-I2422322)
