# Randomized Krylov methods 

randKrylov is a MATLAB library featuring randomized Krylov methods for the solution of large linear systems and eigenvalue problems.
So far it contains functions `randgmres` and `randeigs`, which are counterparts of the built-in `gmres` and `eigs`.

`randgmres` can provide solutions to linear systems several times faster than 
the built-in `gmres`. A speedup over `gmres` can be observed when the GMRES method takes 
many inner iterations: if $m > 10 \log(n)$, where $n$ and $m$, respectively, are the dimension of the 
system and the number of GMRES inner iterations. Moreover, `randgmres` can be not only more efficient 
but also more robust than `gmres` in terms of numerical stability, majorly thanks to the new criterion
for certification of solution. In particular, when `gmres` reports error stagnation, `randgmres` can detect
this stagnation much sooner, or even detect no stagnation and reduce the error to a greater extent (see the experiments).

`randeigs` can provide solutions to unsymmetric eigenvalue problems more than **10x** faster than 
the built-in `eigs`. The speedup is greater when a large Krylov space is used. 
`randeigs` can compute the solution either with the classical Rayleigh–Ritz approximation, as is done in `eigs`, 
or randomized Rayleigh–Ritz approximation, which can allow even more reduction of the computational cost but can be slightly less accurate. 

An additional distinguishing feature of `randgmres` and `randeigs` is the ability to perform the dominant operations
in single precision that can halve the memory usage and speed up the computations. 

## Syntax for randgmres
 `x = randgmres(A,b)` attempts to solve *n* x *n* linear system *Ax* = *b*.
 By default, we use unrestarted method with *m* = min(*n*/10,500) total
 iterations.

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
 lags.

  
## Syntax for randeigs

 `D = randeigs(A)` returns 6 largest eigenvalues of a square matrix *A*.

 `D = randeigs(afun)` takes a function handle *afun* instead of matrix *A*,
  where *afun*(*t*) outputs the matrix-vector products *At*. 
  In the following syntaxes, the matrix A can be replaced by *afun*.

  `D = randeigs(A,B)` solves the generalized eigenvalue problem 
  *Ax* = $\lambda$ *Bx* rather than the standard one. *B* must be a square matrix.
  In the following syntaxes, if the input *B* is missing or is `[]`, then
  randeigs solves the standard eigenvalue problem  *Ax* = $\lambda$ *x*.

  `[V,D] = randeigs(A,B)` returns a diagonal matrix *D* composed of the largest 
  eigenvalues of *A*, and a matrix *V* whose columns are the corresponding 
  eigenvectors. 

  `[V,D,flag] = randeigs(A,B)` also returns a convergence flag indicating 
  the convergence of the eigenvalues.

  `randeigs(A,B,K)` returns *K* largest eigenvalues.

  `randeigs(A,B,K,method)` returns *K* eigenvalues. The options for `method` are 
  same as in built-in eigs function:

      'largestabs' or 'smallestabs' - largest or smallest magnitude
    'largestreal' or 'smallestreal' - largest or smallest real part
                     'bothendsreal' - K/2 values with largest and
                                      smallest real part, respectively
                                      (one more from largest if K is odd)

  For nonsymmetric problems, `method` can also be:
  
    'largestimag' or 'smallestimag' - largest or smallest imaginary part
                     'bothendsimag' - K/2 values with largest and
                                     smallest imaginary part, respectively
                                     (one more from largest if K is odd)

  If `method = sigma` (scalar), then randeigs finds the eigenvalues 
  closest to `sigma`.

  `randeigs(A,B,K,sigma,name,value)` specifies additional options given by the option's name followed by its value. The variable name can be: 

      'Tolerance' - convergence tolerance. Default value: 1e-14.

      'MaxIterations' - maximum number of iterations. Default value: 10

      'SubspaceDimension' - size of Krylov subspace. Default value:
       max(2*K,500).

      'StartVector' - starting vector. Default value: n x 1 random vector.


      'Display'- display diagnostic messages. Default value: 0.

      'InvertOperator' - if method = sigma (scalar) or method = 'smallestabs',
      then we take method = 'largestabs', lambda = sigma+1/lambda 
      and A = (A-sigma*B)^(-1) computed implicitly with Cholesky or lu 
      decomposition. In this case, if used, Afun must output products 
      with (A-sigma*B)^(-1). Default value: 1.

      'SketchingMatrixFun' - specifies parameterless function thetagenfun() 
      that draws a random sketching matrix Theta, given as function handle,
      i.e., Theta(t) = Theta*t. By default, we take Theta as SRHT matrix with 
      k = 2*m*log(n)/log(m) rows, where m is the subspace dimension. 

      'LSsolver' - method used for solving the sketched least-squares problems 
      The options include '3reorth', '5reorth' (default), '20reorth' 
      corresponding to, repsectively, 3, 5, or 20 Richardson iterations, or 'CG' 
      corresponding to 20 iterations of CG, applied to the normal equation.

      'LowPrecision' - compute the Krylov basis vectors in single precision, 
      while performing other (minor) operations in double precision. 
      By default, LowPrecision = 0, i.e., all the operations are performed in 
      double precision.

      'ExactProj' - augment rand.  Arnoldi algorithm with Cholesky QR of the 
      Krylov basis matrix. Performing this step provides a classical 
      Rayleigh-Ritz approximation of the eigenpairs, which can be preferred to 
      its randomized variant due to the strong optimality guarantees. For some 
      problems, however, this can be an unecessary expensive computation. 
      Default value: 1.




## Benchmarks for randgmres

Below is a runtime/robustness comparison of `randgmres` with built-in `gmres`. 
The tests were performed in MATLAB R2021b on a node with 192GB of RAM and 2x Cacade Lake Intel Xeon 5218 16 cores 2.4GHz processor.
The linear systems were taken from [SuiteSparse matrix collection](https://sparse.tamu.edu/) with random right-hand-side vectors. Linear systems ML_Geer, Ga41As41H72 and vas_stokes_1M were preconditioned with `ilu` factorization. It is revealed that for them using a large Krylov space (of dimension $m \geq 400$) is essential. `randgmres` took up to $4$ times lower time to reach the convergence, and is more robust.  For SiO2 and atmosmodd, the `randgmres` function was executed in the low-precision mode (i.e., by taking `lowprecision = 1`) requiring half the memory used by `gmres`.

|system       | size |restart | maxit| `gmres` (time) | `gmres` (error) | `randgmres` (time) | `randgmres` (error)|
| :-------------: | :-------------: |:-------------: | :-------------: |:-------------: | :-------------: |:-------------: | :-------------: |
|ML_Geer      | 1.5x10<sup>6</sup>|1500  |5 | 3486s | 3x**10<sup>-5</sup>**  | 1340s |3.4x10<sup>-10</sup>|
|ML_Geer      | 1.5x10<sup>6</sup>|500  |50 | 4452s | 9.8x10<sup>-10</sup>  | 2512s |5.7x10<sup>-10</sup>|
|ML_Geer      | 1.5x10<sup>6</sup>|200  |150 | 23299s | 9.8x**10<sup>-1</sup>**  | 384s |1.1x**10<sup>**0**</sup>**|
|Ga41As41H72  | 2.7x10<sup>5</sup>| 2000 | 1 | 1205s | 1.1x10<sup>-7</sup>| 304s |1.9x10<sup>-7</sup>|
|Ga41As41H72  | 2.7x10<sup>5</sup>| 500 | 50 | 6590s |2.4x**10<sup>-2</sup>**| 890s | 2.8x**10<sup>-2</sup>**| 
|Ga41As41H72  | 2.7x10<sup>5</sup>| 200 | 150 | 4890s |4.4x**10<sup>-2</sup>**| 257s | 5.5x**10<sup>-2</sup>**| 
|vas_stokes_1M| 1.1x10<sup>6</sup>|800  | 1 |819s | 1.1x10<sup>-12</sup> | 336s |6.13x10<sup>-13</sup>|
|vas_stokes_1M| 1.1x10<sup>6</sup>|400  | 10 |1332s | 6.5x10<sup>-13</sup> | 767s |4.4x10<sup>-13</sup>|
|vas_stokes_1M| 1.1x10<sup>6</sup>|200  | 40 |3441s | 7.7x**10<sup>-3</sup>** | 287s |9.8x**10<sup>-3</sup>**|
|SiO2         | 1.5x10<sup>5</sup>|400  |5 | 90s | 1x10<sup>-10</sup>  | 39s   |2.7x10<sup>-11</sup>|
|atmosmodd    | 1.3x10<sup>6</sup>|300  |5 | 211s | 9.7x10<sup>-11</sup> | 127s |2.4x10<sup>-11</sup>|

## Benchmarks for randeigs

Next we compare `randeigs` with built-in `eigs`. The tests again were executed in MATLAB R2021b on a node with 192GB of RAM and 2x Cacade Lake Intel Xeon 5218 16 cores 2.4GHz processor. The goal was to compute the $K$ smallest eigenvalues, and the corresponding eigenvectors of ML_Geer, Ga41As41H72, vas_stokes_1M and t2em matrices from [SuiteSparse matrix collection](https://sparse.tamu.edu/). The matrices ML_Geer, Ga41As41H72 and vas_stokes_1M could not be efficiently inverted with an `lu` decomposition, therefore they were shifted using `A = speye(n,n)-Ain/(lambdamax+0.1)`, with *A<sub>in</sub>* denoting the input matrix of size *n* and `lambdamax` the magnitude of the largest-magnitute eigenvalue of *A<sub>in</sub>* . Then the eigenpairs of *A* with the largest-magnitude eigenvalues were computed (taking `method = largestabs`). The matrix t2em, on the other hand, could be efficiently inverted with `lu`. Consequently, its eigenpairs with the smallest eigenvalues were computed directly without shifting (taking `A = Ain` and `method = smallestabs`). In the experiments the options `'SubspaceDimension'`, `'MaxIterations'`, `'Tolerance'` and `'ExactProj'` were taken as m, maxiter, tol and 1.


|system       | size |K|m|maxit|tol| `eigs` (time) | `randeigs` (time) | 
| :-------------: | :-------------: |:-------------: | :-------------: |:-------------: | :-------------: |:-------------: | :-------------: |
|ML_Geer      | 1.5x10<sup>6</sup>|20  |1000 |50|10<sup>-10</sup> |-- | 1.5x10<sup>4</sup>s|
|ML_Geer      | 1.5x10<sup>6</sup>|20  |400 |150|10<sup>-10</sup> |6.7x10<sup>4</sup>s| 1.3x10<sup>4</sup>s|
|ML_Geer      | 1.5x10<sup>6</sup>|20  |100 |500|10<sup>-10</sup> |2.6x10<sup>4</sup>s | --|
|vas_stokes_1M| 1.1x10<sup>6</sup>|50  |1500 |25|10<sup>-10</sup> | 1.8x10<sup>**5**</sup>s | 1.3x10<sup>4</sup>s|
|vas_stokes_1M| 1.1x10<sup>6</sup>|50  |200 |2000|10<sup>-10</sup> | 3.1x10<sup>4</sup>s | --|
|t2em      | 9.2x10<sup>5</sup>|300  |800 |5|10<sup>-10</sup> |2.3x10<sup>3</sup>s | 2.8x10<sup>2</sup>s |
|t2em      | 9.2x10<sup>6</sup>|300  |600 |5|10<sup>-10</sup> |1.6x10<sup>3</sup>s | 2.9x10<sup>2</sup>s |


## Reference

[O. Balabanov and L. Grigori, "Randomized Gram-Schmidt process  with application to GMRES", *SIAM J. Sci. Comput.* 44.3 (2022): A1450-A1474.](https://arxiv.org/abs/2011.05090)
