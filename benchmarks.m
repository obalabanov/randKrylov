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

% Construct preconditioner.
n = size(A,1);
[L,U] = ilu(A);

% Precondition the initial matrix.
A = @(x) A*(U\(L\x));

% Define r.h.s.
b = randn(n,1);

% rand-GMRES
m = 1500;
fprintf('rand-GMRES... \n');  
tic
sol = randgmres(A,b,m,1e-9,5);
toc

% built-in GMRES
fprintf('built-in-GMRES... \n');  
tic 
sol2 = gmres(A,b,m,1e-9,5);
toc


%%

fprintf('Second test \n \n');  

% Download Ga41As41H72.mat problem from https://sparse.tamu.edu/

% Load l.h.s. matrix.
t = open('Ga41As41H72.mat');
A = t.Problem.A;
clear t

% Construct preconditioner.
n = size(A,1);
[L,U] = ilu(A);

% Precondition the initial matrix.
A = @(x) A*(U\(L\x));

% Define r.h.s.
b = randn(n,1);

% rand-GMRES
m = 2000;
fprintf('rand-GMRES... \n');  
tic
sol = randgmres(A,b,m,1e-10,1);
toc

% built-in GMRES
fprintf('built-in-GMRES... \n');  
tic 
sol2 = gmres(A,b,m,1e-10,1);
toc


%% 

fprintf('Third test  \n \n');  

% Download vas_stokes_1M.mat problem from https://sparse.tamu.edu/

% Load l.h.s. matrix.
t = open('vas_stokes_1M.mat');
A = t.Problem.A;
clear t

% Construct preconditioner.
n = size(A,1);
[L,U] = ilu(A);

% Precondition the initial matrix.
A = @(x) A*(U\(L\x));

% Define r.h.s.
b = randn(n,1);

% rand-GMRES
m = 800;
fprintf('rand-GMRES... \n');  
tic
sol = randgmres(A,b,m,1e-14,1);
toc

% built-in GMRES
fprintf('built-in-GMRES... \n');  
tic 
sol2 = gmres(A,b,m,1e-14,1);
toc


%% 
fprintf('Fourth test \n \n');  

% Download SiO2.mat problem from https://sparse.tamu.edu/

% Load l.h.s. matrix.
t = open('SiO2.mat');
A = t.Problem.A;
clear t

% no preconditioner.
n = size(A,1);

% Define r.h.s.
b = randn(n,1);

% low-precision (reduced memory) rand-GMRES 
m = 400;
fprintf('low-precision (reduced memory) rand-GMRES... \n');  
tic
sol = randgmres(A,b,m,1e-10,5,[],[],[],1);
toc

% built-in GMRES
fprintf('built-in-GMRES... \n');  
tic 
sol2 = gmres(A,b,m,1e-10,5);
toc

%% 
fprintf('Fifth test \n \n');  

% Download atmosmodd.mat problem from https://sparse.tamu.edu/

% Load l.h.s. matrix.
t = open('atmosmodd.mat');
A = t.Problem.A;
clear t

% no preconditioner.
n = size(A,1);

% Define r.h.s.
b = randn(n,1);

% low-precision (reduced memory) rand-GMRES 
m = 300;
fprintf('low-precision (reduced memory) rand-GMRES... \n');  
tic
sol = randgmres(A,b,m,1e-10,5,[],[],[],1);
toc

% built-in GMRES
fprintf('built-in-GMRES... \n');  
tic 
sol2 = gmres(A,b,m,1e-10,5);
toc


