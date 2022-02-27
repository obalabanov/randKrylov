# Randomized Krylov methods 

randKrylov is a MATLAB library featuring randomized Krylov methods such as GMRES
for the solution of large linear systems.

The function `randgmres` can provide solution to a linear system several times faster than 
the built-in `gmres`. A speedup over `gmres` can be observed when the GMRES method takes 
many inner iterations: if *m* > 10 log(*n*), where *n* and *m*, respectively, are the dimension of the 
system and the number of GMRES inner iterations. 

`randgmres` can be more robust than `gmres` in terms of numerical stability.
<!--- The efficiency gains should be greater for larger *m*. --->

An additional distinguishing feature of `randgmres` is the ability to perform the dominant operations
in single precision that can halve the memory usage and speed up the computations. 

## Syntax

  `x = randgmres(A,b)` attempts to solve *n* x *n* linear system *Ax* = *b*.
  By default, we use unrestarted method with *m* = min(*n*/10,500) total
  iterations.

  `x = randgmres(afun,b)` takes a function handle *afun* instead of matrix *A*,
  where *afun*(*t*) outputs the matrix-vector products *At*.
  To incorporate preconditioning with preconditioner *P* = (*M*<sub>1</sub>*M*<sub>2</sub>)<sup>-1</sup>,
  use `afun = @(x) A*(M2\(M1\x))` or `afun = @(x) M2\(M1\(A*x))`.
  In the following syntaxes, the matrix *A* can be replaced by *Afun*.

  `x = randgmres(A,b,restart)` restarts the algorithm every *m = restart*
  iterations. If *restart* is *n* or `[]` then `randgmres` uses unrestarted method.

  `x = randgmres(A,b,restart,tol)` specifies desired accuracy of the solution
  measured by relative residual error. If *tol* is `[]` then we take *tol* = 1e-6.

  `x = randgmres(A,b,restart,tol,maxit)` specifies the maximum number of
  outer iterations. If *maxit* is `[]` then `randgmres` uses the default value
  min(*n*/*restart*,10).

  `x = randgmres(A,b,restart,tol,maxit,x0)` specifies the initial guess.
  If *x*<sub>0</sub> is `[]` then `randgmres` takes *x*<sub>0</sub> as *n* x 1 zero vector.

  `x = randgmres(A,b,restart,tol,maxit,x0,thetagenfun)` specifies
  parameterless function `thetagenfun()` that draws a random sketching
  matrix *Theta*, given as function handle, i.e., *Theta*(*t*) = *Theta* *t*.
  By default, we take *Theta* as SRHT matrix with *k* = 2*m* log(*n*)/log(*m*)
  rows.

  `x = randgmres(A,b,restart,tol,maxit,x0,thetagenfun,lssolver)` specifies
  the method used for solving the sketched least-squares problems in
  `randgmres`. The options include '3reorth', '5reorth' (default),
  '20reorth' corresponding to, repsectively, 3, 5, or 20 Richardson
  iterations, or 'CG' corresponding to 20 iterations of CG, applied to the
  normal equation.

  `x = randgmres(A,b,restart,tol,maxit,x0,thetagenfun,lssolver,lowprecision)`
  specifies a flag whether to compute the Krylov basis vectors in single
  precision, while performing other (minor) operations in double
  precision. By default, lowprecision = 0, i.e., all the operations are
  performed in double precision.

 `[x,output] = randgmres(___)` returns the output containing the estimated
  residual and stability measure at each iteration, and stagnation/convergence
  flags.


## Benchmarks

Below is a runtime/robustness comparison of `randgmres` with built-in `gmres`. 
The tests were performed in MATLAB R2021b on a node with 192GB of RAM and 2x Cacade Lake Intel Xeon 5218 16 cores 2.4GHz processor.
The linear systems were taken from [SuiteSparse matrix collection](https://sparse.tamu.edu/) with random right-hand-side vectors. Linear systems ML_Geer, Ga41As41H72 and vas_stokes_1M were preconditioned with `ilu` factorization.  For SiO2 and atmosmodd, the `randgmres` function was executed in the low-precision mode (i.e., by taking `lowprecision = 1`) requiring half the memory used by `gmres`.

|system       | size |restart | maxit| `gmres` (time) | `gmres` (error) | `randgmres` (time) | `randgmres` (error)|
| :-------------: | :-------------: |:-------------: | :-------------: |:-------------: | :-------------: |:-------------: | :-------------: |
|ML_Geer      | 1.5x10<sup>6</sup>|1500  |5 | 2663s | 3x10<sup>-5</sup>  | 1338s |3.4x10<sup>-10</sup>|
|Ga41As41H72  | 2.7x10<sup>5</sup>| 2000 | 1 | 945s | 2.7x10<sup>-8</sup>| 320s |6.4x10<sup>-8</sup>|
|vas_stokes_1M| 1.1x10<sup>6</sup>|800  | 1 |664s | 2.3x10<sup>-13</sup> | 322s |1.3x10<sup>-13</sup>|
|SiO2         | 1.5x10<sup>5</sup>|400  |5 | 70s | 9.9x10<sup>-11</sup>  | 45s   |2.9x10<sup>-11</sup>|
|atmosmodd    | 1.3x10<sup>6</sup>|300  |5 | 172s | 9.8x10<sup>-11</sup> | 135s |2.5x10<sup>-11</sup>|


## Reference

[O. Balabanov and L. Grigori, "Randomized Gram-Schmidt process  with application to GMRES", *in press, SIAM J. Sci. Comput.* (2022).](https://arxiv.org/abs/2011.05090)
