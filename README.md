# Randomized Krylov methods 

randKrylov is a MATLAB library that features randomized Krylov methods for solving large linear systems and eigenvalue problems. The functions `randgmres` and `randeigs` can be **over 4x and over 10x faster**, respectively, than the built-in analogues `gmres` and `eigs`, while providing the same (or even greater) accuracy and robustness.

`randgmres` can provide solutions to linear systems several times faster than the built-in `gmres` when the GMRES method takes many inner iterations:
if $m > 10 \log(n)$, where $n$ and $m$ represent the dimension of the system and the number of inner iterations, respectively.
The speedup of GMRES iterations is attributed to orthogonalizing Krylov basis with the randomized Gram-Schmidt (RGS) process, which requires 4x fewer flops and is more cache/memory efficient than the Householder process used in `gmres`. At the same time, the RGS yields as accurate solution with a very high probability $>0.99999999999$. Moreover, `randgmres`  uses a new randomized criterion for the solution certification, which proves to be more robust than the criterion used in `gmres`. In particular, when `gmres` reports error stagnation, `randgmres` can report the stagnation much sooner, or even detect the absence of stagnation and reduce the error to a greater extent (see experiments).

`randeigs` can solve nonsymmetric eigenvalue problems over 10x faster than the built-in `eigs`. The reduction in computational cost is achieved primarily due to orthogonalizing the Krylov basis with the RGS process. The speedup is greater when a large Krylov basis is used.  By default, `randeigs` re-orthogonalizes the computed Krylov basis with a Cholesky QR step at the end of the computation. Such re-orthogonalization should have only minor cost and is employed to ensure that `randeigs` provides the same solution as `eigs`. For many problems this step may be unnecessary and can be ignored by specifying the corresponding flag. In this case `randeigs` will provide a randomized Rayleigh-Ritz solution which can have a similar quality as the regular Rayleigh-Ritz solution from `eigs`.

An additional distinguishing feature of `randgmres` and `randeigs` is the ability to perform the dominant operations in single precision that can halve the memory usage and speed up the computations.

## Syntax for randgmres

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

Here is a runtime/robustness comparison of `randgmres` with the built-in `gmres`.
The tests were performed in MATLAB R2021b on a node with 192GB of RAM and 2x Cacade Lake Intel Xeon 5218 16 cores 2.4GHz processor.
The linear systems were taken from the [SuiteSparse matrix collection](https://sparse.tamu.edu/) with random right-hand-side vectors. The linear systems ML_Geer, Ga41As41H72 and vas_stokes_1M were preconditioned with `ilu` factorization. For SiO2 and atmosmodd, the `randgmres` function was executed in the low-precision mode (i.e., by taking `lowprecision = 1`) requiring half the memory used by `gmres`. 

The following table reports the observed runtimes and residual errors of `gmres` and `randgmres` for varying inputs. It is revealed that `randgmres` has provided an accurate solution up to 4x faster than `gmres`.

For the ML_Geer, Ga41As41H72, and vas_stokes_1M systems, the use of a large number $m \geq 400$ of inner iterations is crucial to circumvent error stagnation. Noteworthy, when a small number of inner iterations was used, `randgmres` was able to identify the stagnation up to 60x faster than `gmres` due to the new randomized criterion.

For the ML_Geer system with $m = 1500$, `gmres` reported an error stagnation with a residual error of 3x10<sup>-5</sup>, while `randgmres` refined the solution to a much greater extent (in addition to being thrice as fast), which demonstrates the robustness of `randgmres`.

|system       | size |restart | maxit| `gmres` time | `gmres` error | `randgmres` time | `randgmres` error|
| :-------------: | :-------------: |:-------------: | :-------------: |-------------: | -------------: |-------------: | -------------: |
|ML_Geer      | 1.5x10<sup>6</sup>|1500  |5 | 3500s | 3x**10<sup>-5</sup>**  | 1300s |3.4x10<sup>-10</sup>|
|ML_Geer      | 1.5x10<sup>6</sup>|500  |50 | 4500s | 9.8x10<sup>-10</sup>  | 2500s |5.7x10<sup>-10</sup>|
|ML_Geer      | 1.5x10<sup>6</sup>|200  |150 | 23300s | 9.8x**10<sup>-1</sup>**  | 380s |1.1x**10<sup>**0**</sup>**|
|Ga41As41H72  | 2.7x10<sup>5</sup>| 2000 | 1 | 1200s | 1.1x10<sup>-7</sup>| 300s |1.9x10<sup>-7</sup>|
|Ga41As41H72  | 2.7x10<sup>5</sup>| 500 | 50 | 6600s |2.4x**10<sup>-2</sup>**| 890s | 2.8x**10<sup>-2</sup>**| 
|Ga41As41H72  | 2.7x10<sup>5</sup>| 200 | 150 | 4890s |4.4x**10<sup>-2</sup>**| 260s | 5.5x**10<sup>-2</sup>**| 
|vas_stokes_1M| 1.1x10<sup>6</sup>|800  | 1 |820s | 1.1x10<sup>-12</sup> | 340s |6.13x10<sup>-13</sup>|
|vas_stokes_1M| 1.1x10<sup>6</sup>|400  | 10 |1300s | 6.5x10<sup>-13</sup> | 770s |4.4x10<sup>-13</sup>|
|vas_stokes_1M| 1.1x10<sup>6</sup>|200  | 40 |3400s | 7.7x**10<sup>-3</sup>** | 290s |9.8x**10<sup>-3</sup>**|
|SiO2         | 1.5x10<sup>5</sup>|400  |5 | 90s | 1x10<sup>-10</sup>  | 40s   |2.7x10<sup>-11</sup>|
|atmosmodd    | 1.3x10<sup>6</sup>|300  |5 | 210s | 9.7x10<sup>-11</sup> | 130s |2.4x10<sup>-11</sup>|

The table below shows the average runtimes taken by a single inner iteration in `gmres` and `randgmres` for varying inputs.

|system       | size |restart | maxit| `gmres` time | `randgmres` time | 
| :-------------: | :-------------: |:-------------: | :-------------: |-------------: | -------------: |
|ML_Geer      | 1.5x10<sup>6</sup>|1500  |5 | 2.6s per it.| 0.78s per it.|
|ML_Geer      | 1.5x10<sup>6</sup>|500  |50 | 1.42s per it.| 0.72s per it.|
|Ga41As41H72  | 2.7x10<sup>5</sup>| 2000 | 1 | 0.67s per it.| 0.17s per it.|
|vas_stokes_1M| 1.1x10<sup>6</sup>|800  | 1 |1.03s per it.| 0.43s per it.|
|vas_stokes_1M| 1.1x10<sup>6</sup>|400  | 10 |0.57s per it.| 0.36s per it.|
|SiO2         | 1.5x10<sup>5</sup>|400  |5 | 0.11s per it.| 0.04s per it.|
|atmosmodd    | 1.3x10<sup>6</sup>|300  |5 | 0.36s per it.| 0.2s per it.|



## Benchmarks for randeigs

Next we compare `randeigs` with built-in `eigs`. The tests again were executed in MATLAB R2021b on a node with 192GB of RAM and 2x Cacade Lake Intel Xeon 5218 16 cores 2.4GHz processor. The goal was to compute the $K$ smallest eigenvalues, and the corresponding eigenvectors of ML_Geer, vas_stokes_1M and t2em matrices from [SuiteSparse matrix collection](https://sparse.tamu.edu/). The matrices ML_Geer and vas_stokes_1M could not be efficiently inverted with an `lu` decomposition, therefore they were shifted using `A = speye(n,n)-Ain/(lambdamax+0.1)`, with *A<sub>in</sub>* denoting the input matrix of size *n* and `lambdamax` the magnitude of the largest-magnitute eigenvalue of *A<sub>in</sub>* . Then the eigenpairs of *A* with the largest-magnitude eigenvalues were computed (taking `method = largestabs`). The matrix t2em, on the other hand, could be efficiently inverted with `lu`. Consequently, its eigenpairs with the smallest eigenvalues were computed directly without shifting (taking `A = Ain` and `method = smallestabs`). In the experiments the options `'SubspaceDimension'`, `'MaxIterations'`, `'Tolerance'` and `'ExactProj'` were taken as m, maxiter, tol and 1.

In all the tests, `randeigs` and `eigs` provided the same solutions up to machine precision. Their runtimes are shown below.

|system       | size |K|m|maxit|tol| `eigs` time | `randeigs` time | 
| :-------------: | :-------------: |:-------------: | :-------------: |:-------------: | :-------------: |-------------: | -------------: |
|ML_Geer      | 1.5x10<sup>6</sup>|20  |1000 |50|10<sup>-10</sup> |151200s | 13400s|
|ML_Geer      | 1.5x10<sup>6</sup>|20  |400 |150|10<sup>-10</sup> |66800s| 13300s|
|ML_Geer      | 1.5x10<sup>6</sup>|20  |100 |500|10<sup>-10</sup> |26200s | --|
|vas_stokes_1M| 1.1x10<sup>6</sup>|50  |1500 |25|10<sup>-10</sup> | 177500s | 13000s|
|vas_stokes_1M| 1.1x10<sup>6</sup>|50  |200 |2000|10<sup>-10</sup> | 31100s | --|
|t2em      | 9.2x10<sup>5</sup>|300  |800 |5|10<sup>-10</sup> |2300s | 280s |
|t2em      | 9.2x10<sup>5</sup>|300  |600 |5|10<sup>-10</sup> |1600s | 290s |


## Reference

[O. Balabanov and L. Grigori, "Randomized Gram-Schmidt process  with application to GMRES", *SIAM J. Sci. Comput.* 44.3 (2022): A1450-A1474.](https://arxiv.org/abs/2011.05090)
