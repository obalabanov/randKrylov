function [xout,output] = randgmres(A,b,restart,tol,maxit,x,thetagenfun,lssolver,lowprecision)
%  Generalized Minimum Residual Method based on rand. Gram-Schmidt process
%
%  x = randgmres(A,b) attempts to solve n x n linear system A*x = b.
%  By default, we use unrestarted method with m = min(n/10,500) total
%  iterations.
%
%  x = randgmres(afun,b) takes a function handle afun instead of matrix A,
%  where afun(t) outputs the matrix-vector products A*t.
%  To incorporate preconditioning with preconditioner P = inv(M1*M2),
%  use afun = @(x) A*(M2\(M1\x)) or afun = @(x) M2\(M1\(A*x)).
%  In the following syntaxes, the matrix A can be replaced by Afun.
%
%  x = randgmres(A,b,restart) restarts the algorithm every m = restart
%  iterations. If restart is n or [] then randGMRES uses unrestarted method.
%
%  x = randgmres(A,b,restart,tol) specifies desired accuracy of the solution
%  measured by relative residual error. If tol is [] then we take tol = 1e-6.
%
%  x = randgmres(A,b,restart,tol,maxit) specifies the maximum number of
%  outer iterations. If maxit is [] then randGMRES uses the default value
%  min(n/restart,10).
%
%  x = randgmres(A,b,restart,tol,maxit,x0) specifies the initial guess.
%  If x0 is [] then randGMRES takes x0 as n x 1 zero vector.
%
%  x = randgmres(A,b,restart,tol,maxit,x0,thetagenfun) specifies
%  parameterless function thetagenfun() that draws a random sketching
%  matrix Theta, given as function handle, i.e., Theta(t) = Theta*t.
%  By default, we take Theta as SRHT matrix with k = 2*m*log(n)/log(m)
%  rows.
%
%  x = randgmres(A,b,restart,tol,maxit,x0,thetagenfun,lssolver) specifies
%  the method used for solving the sketched least-squares problems in
%  randGMRES. The options include '3reorth', '5reorth' (default),
%  '20reorth' corresponding to, repsectively, 3, 5, or 20 Richardson
%  iterations, or 'CG' corresponding to 20 iterations of CG, applied to the
%  normal equation.
%
%  x = randgmres(A,b,restart,tol,maxit,x0,thetagenfun,lssolver,lowprecision)
%  specifies a flag whether to compute the Krylov basis vectors in single
%  precision, while performing other (minor) operations in double
%  precision. By default, lowprecision = 0, i.e., all the operations are
%  performed in double precision.
%
%  [x,output] = randgmres(___) returns the output containing the estimated
%  residual and stability measure at each iteration, and stagnation/convergence
%  flags.
%
%  Copyright (c) 2022, Oleg Balabanov. 
%
%  This program is free software: you can redistribute it and/or modify 
%  it under the terms of the GNU Lesser General Public License as 
%  published by the Free Software Foundation, either version 3 of 
%  the License, or (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful, 
%  but WITHOUT ANY WARRANTY; without even the implied warranty of 
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
%  See the GNU Lesser General Public License for more details.
%
%  You should have received a copy of the GNU Lesser General Public License 
%  along with this program. If not, see <https://www.gnu.org/licenses/>. 
 
    if (nargin < 2)
        error(message('randgmres: not enough input variables'));
    end
    
    if isa(A,'sparse')||isa(A,'double')
        A = @(x) A*x;
    elseif ~isa(A,'function_handle')
        error(message('randgmres: l.h.s. must be a square matrix or function handle'));
    end
    
    n = length(b);
    
    
    % Assign default values to unspecified parameters.
    if (nargin < 3) || isempty(restart) || (restart == n)
        restart = min(floor(n/10), 500);
        restarted = 0;
    else
        restart = max(restart, 1);
        if restart > floor(n/10)
            warning(message('randgmres: too many inner iterations'));
        end
        restarted = 1;
    end
    
    m = restart;
    
    if (nargin < 4) || isempty(tol)
        tol = 1e-6;
    end
    
    if ((nargin < 5) || isempty(maxit)) && restarted
        maxit = min(ceil(n/restart), 10);
    elseif ~restarted
        maxit = 1;
    end
    
    if (nargin >= 6) && ~isempty(x)
        xout = double(x);
        rhs = A(xout) - b;
    else
        xout = zeros(n,1);
        rhs = b;
    end
    
    if (nargin < 7) || isempty(thetagenfun)
        k = min(n, ceil(2*m*log(n)/log(m)));
        thetagenfun = @() SRHT(n,k);
    else
        Theta = thetagenfun();
        k = size(Theta(zeros(n,1)),1);
    end
    
    if (nargin < 8) || isempty(lssolver)
        lssolver = '5reorth';
    end
    
    if (nargin < 9) || isempty(lowprecision)
        lowprecision = 0;
    end
    
    if (nargin > 9)
        error(message('randgmres: too many input variables'));
    end
    
    % Normalize r.h.s.
    n2b = norm(b);
    
    % if b is zero vector
    if n2b == 0
        xout = b;
        output.stag = 0;
        output.conv = 1;
        output.Hres = 0;
        output.skres = 0;
        output.res = 0;
        output.stability = 0;
        return
    end
    
    b = b/n2b;
    rhs = rhs/n2b;
    
    % Allocate space for Krylov basis matrix.
    if lowprecision
        Q = zeros(n,m,'single');
    else
        Q = zeros(n,m,'double');
    end
    
    %% randomized GMRES
    
    for outiter = 1:maxit
        
        % Allocate space for small matrices and vectors.
        R = zeros(m,m);
        S = zeros(k,m);
        P = zeros(k,m);
        q = zeros(n,1);
        w = zeros(n,1);
        s = zeros(k,1);
        p = zeros(k,1);
        
        % Generate random sketching matrix.
        Theta = thetagenfun();
        
        % Perform first iteration.
        s = Theta(rhs);
        p = s;
        P(:,1) = p;
        r = [norm(s);zeros(m-1,1)];
        s = s/r(1);
        q = rhs/r(1);
        
        % error estimator
        skres = zeros(1,m);
        Hres = zeros(1,m);
        skres(1) = r(1);
        Hres(1) = r(1);
        
        % stability measure |S'S-I|_F
        stability = zeros(1,m);
        stability(1) = s'*s - 1;
        
        % solution's coordinates in the Krylov basis
        yout = zeros(m,1);
        
        % stagnation/convergence flags
        stag = 0;
        conv = 0;
        
        % Givens rotation quantities
        cs(1:n) = zeros(n,1);
        sn(1:n) = zeros(n,1);
        g = zeros(m,1);
        g(1) = r(1);
        
        for initer=1:m-1
            
            Q(:,initer) = q;
            S(:,initer) = s;
            R(:,initer) = r;
            
            % Perform RGS iteration.
            if lowprecision
                w = A(double(q));
            else
                w = A(q);
            end
            
            p = Theta(w);
            P(:,initer+1) = p;
            
            r = leastsquares(S,p,lssolver,initer);
            r(initer+1:m) = zeros(m-initer,1);
            
            if lowprecision
                q = single(w) - Q*single(r);
                s = Theta(double(q));
            else
                q = w -Q*r;
                s = Theta(q);
            end
            
            r(initer+1) = norm(s);
            q = q/r(initer+1);
            s = s/r(initer+1);
            
            % Characterize solution.
            H(1:initer+1,initer) = r(1:initer+1);
            
            % Apply Givens rotation.
            for j = 1:initer-1
                tmp           = cs(j)*H(j,initer) + conj(sn(j))*H(j+1,initer);
                H(j+1,initer) = -sn(j)*H(j,initer) + conj(cs(j))*H(j+1,initer);
                H(j,initer)   = tmp;
            end
            
            % Form initer-th rotation matrix.
            [cs(initer),sn(initer)] = givens(H(initer,initer),H(initer+1,initer));
            
            tmp                 = cs(initer)*g(initer);
            g(initer+1)         = -sn(initer)*g(initer);
            g(initer)           = tmp;
            H(initer,initer)    = cs(initer)*H(initer,initer) + conj(sn(initer))*H(initer+1,initer);
            H(initer+1,initer)  = 0;
            
            % Find solution's representation.
            y = H(1:initer,1:initer) \ g(1:initer);
            
            % Approximate residual error.
            Hres(initer+1) = abs(g(initer+1));
            
            % Measure orthogonality of the Krylov basis.
            stability(initer+1) = sqrt(stability(initer)^2 + (s'*s-1)^2 + 2*norm(s'*S)^2);
            
            % Estimate the residual error with random sketching.
            skres(initer+1) = norm(P(:,1) - P*[0; y; zeros(m-initer-1,1)]);
            
            % Check for error stagnation.
            if initer>2 && skres(initer+1) < skres(initer)
                yout = [y; zeros(m-initer,1)];
            elseif Hres(initer+1) < skres(initer+1)/2
                stag = 1;
            end
            
            % only display if the output flag is not used
            if nargout < 2 && (skres(initer+1) < tol/4  ...
                    || mod(initer+1,floor(m/5)) == 0 || initer == m-1 ...
                    || stability(initer+1) > 1 || stag == 1)
                fprintf('randgmres: inner iteration %.2d; ',initer+1)
                fprintf('estimated residual %.1d; ',skres(initer+1))
                fprintf('stability measure %.1d;\n',stability(initer+1))
            end
            
            if skres(initer+1) < tol/4 || initer == m-1 ...
                    || stability(initer+1) > 1 || stag == 1
                % only display if the output flag is not used
                if nargout < 2 && stag == 1
                    fprintf('randgmres: inner iteration %.2d; ',initer+1)
                    fprintf('stagnated.\n')
                end
                if nargout < 2 && stability(initer+1) > 1
                    fprintf('randgmres: inner iteration %.2d; ',initer+1)
                    fprintf('the orthogonality of Krylov basis was lost.\n')
                end
                break
            end
        end
        
        % Update solution.
        x = Q*yout;
        tmp = xout + double(x);
        rhs = b - A(double(tmp));
        res(outiter) = norm(rhs);
        
        % Check accuracy of error estimator.
        if (skres(initer+1) < tol/4) && res(outiter) > 2*tol
            stag = 1;
        end
        
        % only display if the output flag is not used
        if nargout < 2 && res(outiter) < tol
            xout = tmp;
            conv = 1;
            fprintf('\nrandgmres converged at outer(inner) iteration %.2d(%.2d)', outiter, initer+1);
            fprintf(' to relative residual %.1d. \n \n', res(outiter));
            break
        elseif nargout < 2 &&  outiter>1 && res(outiter) > res(outiter-1)*(1-100*eps)
            stag = 1;
            fprintf('\nrandgmres stopped at outer(inner) iteration %.2d(%.2d)', outiter, initer+1);
            fprintf(' because the method stagnated with relative residual %.1d. \n \n', res(outiter));
            break
        elseif nargout < 2 && outiter == maxit
            xout = tmp;           
            if stag == 1
                fprintf('\nrandgmres stopped at outer(inner) iteration %.2d(%.2d)', outiter, initer+1);
                fprintf(' because the method stagnated with relative residual %.1d. \n \n', res(outiter));
                break
            elseif stability(initer+1) > 1
                fprintf('\nrandgmres stopped with relative residual %.1d', res(outiter));
                fprintf(' because the orthogonality of Krylov basis was lost. \n \n');
                break
            else
                fprintf('\nrandgmres stopped with relative residual %.1d', res(outiter));
                fprintf(' because the maximum number of iterations was reached. \n \n');
                break
            end
        elseif nargout < 2
            fprintf('\nrandgmres: outer iteration %.2d; residual %.1d. \n \n', outiter, res(outiter));
        end
        
        stag = 0;
        xout = tmp;
    end
    
    xout = xout*n2b;
    output.stag = stag;
    output.conv = conv;
    output.Hres = Hres;
    output.skres = skres;
    output.res = res;
    output.stability = stability;
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
            r = r + dr;
        end
    else
        % Conjugate Gradient
        Stemp = S(:,1:initer);
        [rtemp,~,~,~,~] = pcg(@(x) ((Stemp*x)'*Stemp)',Stemp'*p,1.0e-14,20);
        r = [rtemp; zeros(size(S,2) - initer,1)];
    end
end


%% Subsampled randomized Hadamard Transform (SRHT).
function Theta = SRHT(n,k)
    D = randi([0 1], n,1)*2 - 1;
    N = 2^ceil(log(n)/log(2));
    perm = randperm(N,k);
    select = @(t,ind) t(ind);
    Theta = @(t) (1/sqrt(k)) * select(myfwht(D.*t),perm);
end


%% Fast Walsh Hadamard Transform
function z = myfwht(a)
    h = 1;
    n = length(a);
    N = 2^ceil(log(n)/log(2));
    z = zeros(N,1);
    z(1:n) = a;
    while h < N
        for i = 1:2*h:N
            for j = i:(i + h - 1)
                x = z(j);
                y = z(j + h);
                z(j) = x + y;
                z(j + h) = x - y;
            end
        end
        h = 2*h;
    end
end

%% Givens rotation
function [cs,sn] = givens(a,b)
    if(b == 0)
        cs = 1;
        sn = 0;
    elseif (abs(b) > abs(a))
        tmp = a/b;
        sn = 1/sqrt(1 + tmp^2);
        cs = tmp*sn;
    else
        tmp = b/a;
        cs = 1/sqrt(1 + tmp^2);
        sn = tmp*cs;
    end
end