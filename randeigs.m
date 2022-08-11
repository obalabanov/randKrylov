function [V, D, flag] = randeigs(varargin)
%  Compute several eigenvalues and eigenvectors of a matrix by using  
%  Krylov-Schur algorithm based on rand. Arnoldi process. 
%  
%  D = randeigs(A) returns 6 largest eigenvalues of a square matrix A.
%
%  D = randeigs(afun) takes a function handle afun instead of matrix A,
%  where afun(t) outputs the matrix-vector products A*t. 
%  In the following syntaxes, the matrix A can be replaced by afun.
%
%  D = randeigs(A,B) solves the generalized eigenvalue problem 
%  Ax = lambda Bx rather than the standard one. B must be a square matrix.
%  In the following syntaxes, if the input B is missing or is [], then
%  randeigs solves the standard eigenvalue problem Ax = lambda x.
%
%  [V,D] = randeigs(A,B) returns a diagonal matrix D composed of the largest 
%  eigenvalues of A, and a matrix V whose columns are the corresponding 
%  eigenvectors. 
%
%  [V,D,flag] = randeigs(A,B) also returns a convergence flag indicating 
%  the convergence of the eigenvalues.
%
%  randeigs(A,B,K) returns K largest eigenvalues.
%
%  randeigs(A,B,K,method) returns K eigenvalues. The options for method are 
%  same as in built-in eigs function:
%
%      'largestabs' or 'smallestabs' - largest or smallest magnitude
%    'largestreal' or 'smallestreal' - largest or smallest real part
%                     'bothendsreal' - K/2 values with largest and
%                                      smallest real part, respectively
%                                      (one more from largest if K is odd)
%
%  For nonsymmetric problems, method can also be:
%    'largestimag' or 'smallestimag' - largest or smallest imaginary part
%                     'bothendsimag' - K/2 values with largest and
%                                     smallest imaginary part, respectively
%                                     (one more from largest if K is odd)
%
%  If method = sigma (scalar), then randeigs finds the eigenvalues 
%  closest to sigma.
%
%  randeigs(A,B,K,sigma,name,value) specifies additional options given by 
%  the option's name followed by its value. The variable name can be: 
%
%  'Tolerance' - convergence tolerance. Default value: 1e-14.
%
%  'MaxIterations - maximum number of iterations. Default value: 10
%
%  'SubspaceDimension' - size of Krylov subspace. Default value:
%   max(2*K,500).
%
%  'StartVector' - starting vector. Default value: n x 1 random vector.
%
%  'Display'- display diagnostic messages. Default value: 0.
%
%  'InvertOperator' - if method = sigma (scalar) or method = 'smallestabs',
%  then we take method = 'largestabs', lambda = sigma+1/lambda 
%  and A = (A-sigma*B)^(-1) computed implicitly with chol or lu 
%  decomposition. In this case, if used, Afun must output products 
%  with (A-sigma*B)^(-1). Default value: 1.
%
%  'SketchingMatrixFun' - specifies parameterless function thetagenfun() 
%  that draws a random sketching matrix Theta, given as function handle,
%  i.e., Theta(t) = Theta*t. By default, we take Theta as SRHT matrix with 
%  k = 2*m*log(n)/log(m) rows, where m is the subspace dimension. 
%
%  'LSsolver' - method used for solving the sketched least-squares problems 
%  The options include '3reorth', '5reorth' (default), '20reorth' 
%  corresponding to, repsectively, 3, 5, or 20 Richardson iterations, or 'CG' 
%  corresponding to 20 iterations of CG, applied to the normal equation.
%
%  'LowPrecision' - compute the Krylov basis vectors in single precision, 
%  while performing other (minor) operations in double precision. 
%  By default, LowPrecision = 0, i.e., all the operations are performed in 
%  double precision.
%
%  'ExactProj' - augment rand.  Arnoldi algorithm with Cholesky QR of the 
%  Krylov basis matrix. Performing this step provides a classical 
%  Rayleigh-Ritz approximation of the eigenpairs, which can be preferred to 
%  its randomized variant due to the strong optimality guarantees. For some 
%  problems, however, this can be an unecessary expensive computation. 
%  Default value: 1.
%
% 
%  References:
%  O. Balabanov and L. Grigori, "Randomized Gram-Schmidt process  with 
%  application to GMRES", SIAM J. Sci. Comput. 44.3 (2022): A1450-A1474.
%
%  G.W. Stewert, "A Krylov-Schur Algorithm for Large Eigenproblems."
%  SIAM. J. Matrix Anal. & Appl., 23(3), 601-614.
%
%  Copyright (c) 2022, Oleg Balabanov.
%
%  This program is free software: you can redistribute it and/or modify 
%  it under the terms of the GNU Lesser General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or 
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful, 
%  but WITHOUT ANY WARRANTY; without even the implied warranty of 
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
%  General Public License for more details.
%
%  You should have received a copy of the GNU Lesser General Public License 
%  along with this program. If not, see <https://www.gnu.org/licenses/>. 



% Assign values to problem's parameters.
[A,M,n,K,method,b,tol,m,maxit,display,thetagenfun,Theta,k,lssolver,lowprecision,...
reorthog,postprocesslambda] = checkInputs(varargin{:});

clear varargin options

% Allocate space for Krylov basis matrix.
if lowprecision
    Q = zeros(n,m,'single');  
    V = zeros(n,K,'single');  
else
    Q = zeros(n,m,'double');
    V = zeros(n,K,'double');  
end

K0 = K; % save original K
H = []; % initialize the Hessenberg matrix
nconv = 0;
sizeQ = 1;

% Allocate space for small matrices and vectors.
S = zeros(k,m);
P = zeros(k,m+1);
q = zeros(n,1);
w = zeros(n,1);
s = zeros(k,1);
p = zeros(k,1);

% Perform first inner iteration.
s = Theta(b);
p = s;
P(:,1) = p;
r = [norm(s); zeros(m-1,1)];
s = s/r(1);
q = b/r(1);  

Q(:,1) = q;
S(:,1) = s;
  
for outiter = 1:maxit  
    
    % Measure orthogonality of the basis obtained at previous iteration.
    stability(sizeQ) = norm(eye(sizeQ,sizeQ)-S(:,1:sizeQ)'*S(:,1:sizeQ),'fro');
     
    % Enrich the basis with more Krylov vectors.
    for initer = sizeQ:m
      Q(:,initer) = q;
      S(:,initer) = s;
         
      if lowprecision
         w = A(double(q));
      else
         w = A(q);
      end
      
      % Perform RGS orthogonalization.
      p = Theta(w);
      P(:,initer+1) = p;
      
      r = leastsquares(S,p,lssolver,initer);     
      r(initer+1:m) = zeros(m-initer,1);
       
      if lowprecision
          q = single(w)-Q*single(r);
          s = Theta(double(q));
      else
          q = w -Q*r;
          s = Theta(q);
      end
    
      r(initer+1) = norm(s);
      q = q/r(initer+1);
      s = s/r(initer+1); 
      
      % If numerically invariant subspace is found, then reorthogonalize q.
      if norm(s'*S)>1e-4
        p2 = Theta(q);
        r2 = leastsquares(S,p2,lssolver,initer);  
        r2(initer+1:m) = zeros(m-initer,1);
       
        if lowprecision
          q = single(q)-Q*single(r2);
          s = Theta(double(q));
        else
          q = q -Q*r2;
          s = Theta(q);
        end
        r2(initer+1) = norm(s);
        q = q/r2(initer+1);
        s = s/r2(initer+1); 
        
        % Retrive r.
        if initer<m
            S(:,initer+1) = s;
            r = leastsquares(S,p,lssolver,initer);  
            r(initer+2:m) = zeros(m-initer-1,1);
        else 
            r = leastsquares([S,s],p,lssolver,initer);  
            r(initer+2:m) = zeros(m-initer-1,1);
        end
      end
      
      % Measure orthogonality of the Krylov basis.  
      stability(initer+1) = sqrt(stability(initer)^2+(s'*s-1)^2+norm(s'*S)^2);
      
      % only display if the corresponding flag is 1 
      if display && (mod(initer+1,floor(m/5))== 0 || initer == m-1)
         fprintf('randeigs: Krylov dimension %.2d; ',initer+1)
         fprintf('stability measure %.1d;\n',stability(initer+1))
      end

      H(1:initer+1,initer) = r(1:initer+1);
    end
    
    % Orthogonalize Krylov basis with respect to the M-inner product with
    % Cholesky QR.
    if reorthog
        if isdiag(M) && ~all(abs(diag(M) - 1)) % check whether M = I
           QQ = [[Q'*Q, (q'*Q)']; [(q'*Q), q'*q]];
        else
           Mq = M*q;
           QQ = [[Q'*M*Q, (Mq'*Q)']; [(Mq'*Q), Mq'*q]]; 
        end
        r1 = chol(double(QQ));
        %Q = Q/r1(1:m,1:m); done implicitly
        S = S/r1(1:m,1:m); 
        q = (q - Q*(r1(1:m,1:m)\r1(1:m,end)))/r1(end,end);
        s = Theta(q);
        H = r1*H(1:m+1,1:m)/r1(1:m,1:m);
    end    
    
    % Compute Schur decomposition of H.
    [X, T] = schur(H(1:end-1, :));    
    
    % Compute eigenpairs.
    [U, d] = eig(T, 'vector');
    U = X*U;
    
    % Measure the error of each Ritz pair.
    res = abs(H(end, :)*U);
    
    % Sort eigenvalues and residuals (see built-in eigs function).
    ind = whichEigenvalues(d, method);
    d = d(ind);
    res = res(ind);
    
    % Number of converged eigenpairs (see built-in eigs function):
    nconvold = nconv;
    isNotConverged = ~(res(1:K0)' < tol*max(eps^(2/3), abs(d(1:K0))));
    nconv = nnz(~isNotConverged);
    
    if nconv >= K0 || outiter == maxit
        % Stop the algorithm now (see built-in eigs function)
        
        % only display if the corresponding flag is 1 
        if display && nconv >= K0
            fprintf('\nrandeigs converged at outer(inner) iteration %.2d \n \n', outiter);
        elseif display 
            fprintf('\nrandeigs stopped because the maximum number of iterations was reached. \n');    
            fprintf('number of converged eigenvalues: %.1d. \n \n', nconv);
        end
        break;
    else
        % Adjust K to prevent stagnation (see built-in eigs function).
        K = K0 + min(nconv, floor((m - K0) / 2));
        if K == 1 && m > 3
            K = floor(m / 2);
        end
        if K + 1 < m && nconvold > nconv
            K = K + 1;
        end
    end
    
    % Get original ordering of eigenvalues back (see built-in eigs function).
    d = ordeig(T);
    
    % Choose desired eigenvalues in d to create a Boolean select vector (see built-in eigs function).
    ind = whichEigenvalues(d, method);
    ind = ind(1:K);
    select = false(1, m);
    select(ind) = true;
    
    % If H is real, make sure both parts of a conjugate pair are present. (see built-in eigs function)
    if isreal(H)
        for i = ind'
            if i < m && T(i+1,i) ~= 0 && ~select(i+1)
                select(i+1) = true;
                K = K+1;
            end
            if i > 1 && T(i, i-1) ~= 0 && ~select(i-1)
                select(i-1) = true;
                K = K+1;
            end
        end
    end
    
    
    % Only display if the corresponding flag is 1.
    if display
      fprintf('\nrandeigs: outer iteration %.2d; number of converged eigenvalues: %.1d. \n \n', outiter, nconv);
    end
    
    % Reorder X and T based on select (see built-in eigs function).
    [X, T] = ordschur(X, T, select);
    
    % Compute quantitites for restarting the algorithm (see built-in eigs function).
    H = [T(1:K,1:K); H(end,:)*X(:,1:K)];
    if reorthog
        Q(:,1:K) = Q*(r1(1:m,1:m)\X(:, 1:K)); 
    else
        Q(:,1:K) = Q*X(:, 1:K); 
    end
    Q(:,K+1) = q;
    
    % Regenerate random sketching matrix.
    Theta = thetagenfun();
    
    % Orthogonalize Q(:,1:K+1) with respect to <Theta.,Theta.> using
    % sketched Cholesky QR.
    for i = 1:K+1
       S(:,i) = Theta(Q(:,i));
    end
    [q1,r1] = qr(S(:,1:K+1),0);
    invr1 = inv(r1(1:K,1:K)); % stable since r1 is well-conditioned 
    S(:,1:K+1) =  q1;
    S(:,K+2:m) = zeros(k,m-K-1);
    Q(:,1:K) = Q(:,1:K)*invr1;
    q = (q - Q(:,1:K)*r1(1:K,end))/r1(end,end);
    Q(:,K+1) = q;
    s = Theta(q);
    H = r1*H(1:K+1,1:K)*invr1;   
    
    % Compute Schur decomposition and update H and Q,S,P
    [X, T] = schur(H(1:end-1, :));
    H = [T(1:K, 1:K); H(end, :) * X(:, 1:K)];
    Q(:,1:K) = Q(:,1:K) * X(:, 1:K); 
    S(:,1:K) = S(:,1:K) *  X(:, 1:K); 
    P(:,2:K+1) = P(:,2:K+1)* X(:, 1:K);
    
    sizeQ = K+1;
end

% Compute solution.
if reorthog
    C = r1(1:m,1:m)\U(:,ind(1:K0));
else
    C = U(:,ind(1:K0));
end

if isreal(C) || ~isreal(Q)
    V = Q*C;
else
    V = Q*real(C)+1i*(Q*imag(C));
end

% Normalize.
normV = sqrt(sum(abs(V).^2,1));
V = V*diag(1./normV);

if nargout == 1
    V = postprocesslambda(d(1:K0));
else
    D = diag(postprocesslambda(d(1:K0)));
end
   flag = isNotConverged;
end




    
function [A,M,n,K,method,b,tol,m,maxit,display,thetagenfun,Theta,k,lssolver,...
lowprecision,reorthog,postprocesslambda]= checkInputs(varargin)

if issparse(varargin{1})
    n = size(varargin{1},1);
    A = varargin{1};
    argind = 2;
    if issparse(varargin{2})||isa(varargin{2},'function_handle')
        B = varargin{2};
        argind = 3;
    elseif isempty(varargin{2})
        B = speye(n,n);
        argind = 3;
    else
        B = speye(n,n);
    end
elseif isa(varargin{1},'function_handle')    
    A = varargin{1};
    if isa(varargin{2},'sparse') || isa(varargin{2},'function_handle')
        B = varargin{2};
        n = varargin{3};
        argind = 4;
    elseif isempty(varargin{2})
        n = varargin{3};
        B = speye(n,n);
        argind = 4;
    else
        n = varargin{2};
        B = speye(n,n);
        argind = 3;
    end  
end

if (nargin < argind) || isempty(varargin{argind}) || (varargin{argind} > n) ||  (varargin{argind} < 1)
    K = 6;
else
    K = varargin{argind};
    if K > floor(n/10)
        warning(message('randeigs: too many eigenvalues'));
    end
end

argind = argind+1;

if (nargin < argind) || isempty(varargin{argind})
    method = 'largestabs';
else
    method = varargin{argind};
end

argind = argind+1;

if (nargin >= argind)
    p = inputParser; 
    addOptional(p,'StartVector',[]);
    addOptional(p,'Tolerance',1e-14);
    addOptional(p,'SubspaceDimension',max(2*K,min(floor(n/10),500)));
    addOptional(p,'MaxIterations',[]);
    addOptional(p,'Display',0);
    addOptional(p,'InvertOperator',1);
    addOptional(p,'SketchingMatrixFun',[]);
    addOptional(p,'LSsolver','5reorth');
    addOptional(p,'LowPrecision',0);
    addOptional(p,'ExactProj',1);

    parse(p,varargin{argind:nargin});
    options = p.Results;  
    
    if ~isempty(options.StartVector)
        b = options.StartVector;
    else
        randStr = RandStream('dsfmt19937','Seed',0);
        b = randn(randStr, n,1);
    end
    tol = options.Tolerance;
    m = options.SubspaceDimension;
    if ~isempty(options.MaxIterations)
        maxit = options.MaxIterations;
    else
        maxit = min(ceil(n/m),10);
    end
    display = options.Display;
    
    if ~isempty(options.SketchingMatrixFun)
        thetagenfun = options.SketchingMatrixFun;
        Theta = thetagenfun();
        k = size(Theta(zeros(n,1)),1); 
    else
        k = min(n,ceil(2*m*log(n)/log(m)));
        thetagenfun = @() SRHT(n,k); 
        Theta = thetagenfun();        
    end
    
    lssolver = options.LSsolver; 
    lowprecision = options.LowPrecision; 
    reorthog = options.ExactProj;
    invop = options.InvertOperator;
    postprocesslambda = @(lambda) lambda;
    M = speye(n,n);
    if isscalar(method)
        sigma = method;
        method = 'smallestabs';
        if issparse(A)
            A = A-sigma*B;
        end
        postprocesslambda = @(lambda) postprocesslambda(lambda+sigma);
    end
    
    if ~(strcmp(method,'smallestabs') && invop) % A = B^(-1)*A
       flagB = 1;
       if ishermitian(B)
          [R,flagB,P] = chol(B);
       end
       if ~flagB  %B is SPD
          A = @(x) P*(R\(((A*x)'*P)/R)'); 
       end
       if flagB %B is not SPD
          [L,U,P,Q,D] = lu(B); 
          A = @(x) D*(Q*(U\(L\(P*(A*x))))); 
       end
    end   
    
    if strcmp(method,'smallestabs') && invop %shift and invert      
        method = 'largestabs';
        % Check whether B is SPD
        flagB = 1;
        if ishermitian(B)
              [R,flagB,P] = chol(B);
        end
        if flagB % A = A^{-1}*B
           flag = 1;
            if ~isa(A,'function_handle') && ishermitian(A)
               [R,flag,P] = chol(A);
            end
            if ~isa(A,'function_handle') && ~flag 
               A = @(x) P*(R\(((B*x)*P)'/R)'); 
            end
            if ~isa(A,'function_handle') && flag 
                [L,U,P,Q,D] = lu(A); 
                A = @(x) D*(Q*(U\(L\(P*(B*x))))); 
            end 
            if isa(A,'function_handle') 
                A = @(x) A*(B*x);
            end
        else % A = A^{-1}, M = B, Theta = Theta*R
            M = B; 
            thetagenfun = @() addmetric(thetagenfun,R);
            Theta = thetagenfun();    
            flag = 1;
            if ~isa(A,'function_handle') && ishermitian(A)
               [R,flag,P] = chol(A);
            end
            if ~isa(A,'function_handle') && ~flag 
               A = @(x) P*(R\((x'*P)/R)'); 
            end
            if ~isa(A,'function_handle') && flag 
                [L,U,P,Q,D] = lu(A); 
                A = @(x) (Q*(U\(L\(P*(D\x)))));
            end 
        end        
        postprocesslambda = @(lambda) postprocesslambda(1./lambda);    
    end
end
end


%% Solving nearly orthogonal least-squares problem
function r = leastsquares(S,p,lssolver,initer)
    if strcmp(lssolver,'3reorth') 
        j = 3;
    elseif strcmp(lssolver,'5reorth') 
        j = 5;
    elseif strcmp(lssolver,'20reorth') 
        j = 20;   
    else
        j = 0;
    end
    if j ~= 0
        % Richardson iterations
        r = (p'*S)';
        p = p - S*r;
        for i=1:j-1
            dr = (p'*S)';
            p = p - S*dr;
            r = r+dr;
        end
    else
        % Conjugate Gradient
        Stemp = S(:,1:initer);
        [rtemp,~,~,~,~] = pcg(@(x) ((Stemp*x)'*Stemp)',Stemp'*p,1.0e-14,20);
        r = [rtemp; zeros(size(S,2)-initer,1)];
    end
end


%% Subsampled randomized Hadamard Transform (SRHT).
function Theta = SRHT(n,k)
    D = randi([0 1], n,1)*2-1;
    N = 2^ceil(log(n)/log(2));
    perm =  randperm(N,k);   
    select = @(t,ind) t(ind); 
    Theta = @(t) (1/sqrt(k))*select(myfwht(D.*t),perm);
end

%% Add metric to the Theta generating function
function thetafun = addmetric(thetagenfun,R)
    thetafun = thetagenfun ();
    thetafun = @(x)  thetafun(R*x);
end


%% Fast Walsh Hadamard Transform
function z = myfwht(a)
    if ~isreal(a)
       z = complex(myfwht(real(a)),myfwht(imag(a)));
       return
    end
    h = 1;
    n = length(a);
    N = 2^ceil(log(n)/log(2));
    z = zeros(N,1);
    z(1:n) = a;
    while h < N  
        for i=1:2*h:N
            for j=i:(i+h-1)
               x = z(j);
               y = z(j + h);
               z(j) = x + y;
               z(j + h) = x - y;
            end
        end
        h = 2*h;
    end
end

% Taken from built-in eigs function
function ind = whichEigenvalues(d, method)

switch method
    case 'largestabs'
        [~, ind] = sort(abs(d), 'descend');
    case 'largestreal'
        [~, ind] = sort(real(d), 'descend');
    case 'smallestreal'
        [~, ind] = sort(real(d), 'ascend');
    case 'largestimag'
        [~, ind] = sort(imag(d), 'descend');
    case 'smallestimag'
        [~, ind] = sort(imag(d), 'ascend');
    case 'bothendsreal'
        [~, ind] = sort(real(d), 'descend');
        ind2 = [ind, flip(ind)]';
        ind2 = ind2(:);
        ind = ind2(1:size(d,1));
    case 'bothendsimag'
        [~, ind] = sort(imag(d), 'descend');
        ind2 = [ind, flip(ind)]';
        ind2 = ind2(:);
        ind = ind2(1:size(d,1));
    case 'smallestabs'
        [~,ind] = sort(abs(d), 'ascend');
    case 'smallestimagabs'
        [~,ind] = sort(abs(imag(d)), 'ascend');
end

end
