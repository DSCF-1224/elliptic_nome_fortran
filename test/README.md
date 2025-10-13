# test

## How to generate the reference values

### [k.dat](k.dat)

```mathematica
N[1/Sqrt[2],39]
```

### [q.dat](q.dat)

```mathematica
N[EllipticNomeQ[Power[Range[0,99]/100,2]],39]
N[EllipticNomeQ[Power[Sqrt[2],-2]],39]
```

[EllipticNomeQ: Nome of an elliptic function—Wolfram Documentation](https://reference.wolfram.com/language/ref/EllipticNomeQ.html)
