# tune-params
Tune parameters of a zeta zeros algorithm.

Requires [flint](https://github.com/wbhart/flint2)
and [arb](https://github.com/fredrik-johansson/arb).

Here's an example of how it's used. In the output, x is log height.
```
> gcc tune.c
> ./a.out > table.txt
> Rscript fit.R table.txt
logJ = -0.506132895307753 + 0.497565124541755*x + 2.25705259466534e-06*x*x;
K = 103.630048145445 + -0.867460371390699*x + 0.00106534037344869*x*x;
grid = 4085.55616952791 + -1.36695690825918*x + -0.0198932963863186*x*x;
interp = 18.9999999999996 + 1.58654931768833e-14*x + -1.69448761125292e-16*x*x;
h = 142.523364391385 + -0.057841628146848*x + 0.0157997314814483*x*x;
H = 0.727682145349253 + -0.00353076385259656*x + 0.000111471602416957*x*x;
```
