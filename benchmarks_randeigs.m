%   Copyright (c) 2022, Oleg Balabanov. 
%
%   This program is free software: you can redistribute it and/or modify 
%   it under the terms of the GNU Lesser General Public License as 
%   published by the Free Software Foundation, either version 3 of 
%   the License, or (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful, 
%   but WITHOUT ANY WARRANTY; without even the implied warranty of 
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
%   See the GNU Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public License 
%   along with this program. If not, see <https://www.gnu.org/licenses/>. 
%% 

fprintf('First test \n \n');  

% Download ML_Geer.mat problem from https://sparse.tamu.edu/

% Load l.h.s. matrix.
t = open('ML_Geer.mat');

A = t.Problem.A;
clear t

m = 1000;
maxiter = 50;
tol = 1e-10;
K = 20;

% Shift A, since lu(A) is too expensive.
[v,lambdamax] = eigs(A,1);
lambdamax = norm(A*v)/norm(v);
n = size(A,1);
b = randn(n,1);
A = speye(n,n)-A/(lambdamax+0.1); 


fprintf('randeigs with classical Rayleigh–Ritz approximation (classical Galerkin)... \n');  

tic
[V,D,flag] = randeigs(A,[],K,'largestabs','SubspaceDimension',m,'MaxIterations', maxiter,'Tolerance',tol,'ExactProj',1,'Display',1,'InvertOperator',1);
toc
clear V

%
fprintf('randeigs with randomozed Rayleigh–Ritz approximation (sketched Galerkin)... \n');  
tic
[V,D,flag] = randeigs(A,[],K,'largestabs','SubspaceDimension',m,'MaxIterations', maxiter,'Tolerance',tol,'ExactProj',0,'Display',1,'InvertOperator',1);
toc
clear V


%Warning: eigs for unsymmetric matrices can be computationally heavy!
fprintf('built-in eigs... \n');  
tic
[V,D,flag] = eigs(A, K,'largestabs','SubspaceDimension',m,'MaxIterations',maxiter,'Tolerance',tol,'FailureTreatment','keep', 'Display',1);
toc
clear V

%% 

fprintf('Second test \n \n');  

% Download vas_stokes_1M.mat problem from https://sparse.tamu.edu/

% Load l.h.s. matrix.
t = open('vas_stokes_1M.mat');

A = t.Problem.A;
clear t

m = 1500;
maxiter = 25;
tol = 1e-10;
K = 50;

% Shift A, since lu(A) is too expensive.
[v,lambdamax] = eigs(A,1);
lambdamax = norm(A*v)/norm(v);
n = size(A,1);
b = randn(n,1);
A = speye(n,n)-A/(lambdamax+0.1); 


fprintf('randeigs with classical Rayleigh–Ritz approximation (classical Galerkin)... \n');  

tic
[V,D,flag] = randeigs(A,[],K,'largestabs','SubspaceDimension',m,'MaxIterations', maxiter,'Tolerance',tol,'ExactProj',1,'Display',1,'InvertOperator',1);
toc
clear V

%
fprintf('randeigs with randomized Rayleigh–Ritz approximation (sketched Galerkin)... \n');  
tic
[V,D,flag] = randeigs(A,[],K,'largestabs','SubspaceDimension',m,'MaxIterations', maxiter,'Tolerance',tol,'ExactProj',0,'Display',1,'InvertOperator',1);
toc
clear V


%Warning: eigs for unsymmetric matrices can be computationally heavy!
fprintf('built-in eigs... \n');  
tic
[V,D,flag] = eigs(A, K,'largestabs','SubspaceDimension',m,'MaxIterations',maxiter,'Tolerance',tol,'FailureTreatment','keep', 'Display',1);
toc
clear V


%% 

fprintf('Third test \n \n');  

% Download t2em.mat problem from https://sparse.tamu.edu/

% Load l.h.s. matrix.
t = open('t2em.mat');

A = t.Problem.A;
clear t

m = 600;
maxiter = 5;
tol = 1e-10;
K = 300;

% No need to shift A, since lu(A) is cheap.

fprintf('randeigs with classical Rayleigh–Ritz approximation (classical Galerkin)... \n');  

tic
[V,D,flag] = randeigs(A,[],K,'smallestabs','SubspaceDimension',m,'MaxIterations', maxiter,'Tolerance',tol,'ExactProj',1,'Display',1,'InvertOperator',1);
toc
clear V

%
fprintf('randeigs with randomized Rayleigh–Ritz approximation (sketched Galerkin)... \n');  
tic
[V,D,flag] = randeigs(A,[],K,'smallestabs','SubspaceDimension',m,'MaxIterations', maxiter,'Tolerance',tol,'ExactProj',0,'Display',1,'InvertOperator',1);
toc
clear V


%Warning: eigs for unsymmetric matrices can be computationally heavy!
fprintf('built-in eigs... \n');  
tic
[V,D,flag] = eigs(A, K,'smallestabs','SubspaceDimension',m,'MaxIterations',maxiter,'Tolerance',tol,'FailureTreatment','keep', 'Display',1);
toc
clear V
















