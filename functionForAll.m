function [x_next, u_opt, J_opt] = functionForAll(Np, xk, parameters, u)
%% Expand parameters
xcon_init = xk(1);
xaux_init = xk(2);
x_max = parameters.x_max;
N_max_gr = parameters.N_max_gr;
gamma0 = parameters.gamma0;
gamma1 = parameters.gamma1; 
gamma2 = parameters.gamma2;

x_eff = parameters.x_eff;
psi = parameters.psi;
a1 = parameters.A(1); a2 = parameters.A(2); a3 = parameters.A(3);
b1 = parameters.B(1); b2 = parameters.B(2); b3 = parameters.B(3);

lambda = parameters.lambda;

Mu = 2; mu = 0; Mcon = 70; mcon = 0; Maux = N_max_gr; maux = 0; eps = 10^-3;

Nx = Np*26+2;

%% Inequalities - delta and z definitions

A1 = zeros(63, 26); B1 = zeros(63,1); signs = char(zeros(63,1)); 

i=0;
% delta0
i = i+1; A1(i,3) = 1; A1(i,4) = Mu; B1(i) = Mu; signs(i) = 'U'; % <=
i = i+1; A1(i,3) = 1; A1(i, 4) = -(mu-eps); B1(i) = eps; signs(i) = 'L'; % >=

% delta 0n - NOT
i = i+1; A1(i,4) = 1; A1(i, 5) = 1; B1(i) = 1; signs(i) = 'U';
i = i+1; A1(i,4) = -1; A1(i, 5) = -1; B1(i) = -1; signs(i) = 'U';

% delta0a
i = i+1; A1(i,1) = 1; A1(i, 6) = Mcon-30 + eps; B1(i) = Mcon; signs(i) = 'U';
i = i+1; A1(i,1) = 1; A1(i,6) = -(mcon-30-eps); B1(i) = 30; signs(i) = 'L';

% delta 0an - NOT
i = i+1; A1(i,6) = 1; A1(i, 7) = 1; B1(i) = 1; signs(i) = 'U';
i = i+1; A1(i,6) = -1; A1(i, 7) = -1; B1(i) = -1; signs(i) = 'U';

% delta0b
i = i+1; A1(i,1) = 1; A1(i, 8) = Mcon-50 + eps; B1(i) = Mcon; signs(i) = 'U';
i = i+1; A1(i,1) = 1; A1(i,8) = -(mcon-50-eps); B1(i) = 50; signs(i) = 'L';

% delta 0bn - NOT
i = i+1; A1(i,8) = 1; A1(i, 9) = 1; B1(i) = 1; signs(i) = 'U';
i = i+1; A1(i,8) = -1; A1(i, 9) = -1; B1(i) = -1; signs(i) = 'U';

% delta 0,1 - 2xAND
i = i+1; A1(i, 4) = -1; A1(i,10) = 1; B1(i) = 0; signs(i) = 'U'; 
i = i+1; A1(i, 6) = -1; A1(i,10) = 1; B1(i) = 0; signs(i) = 'U'; 
i = i+1; A1(i,4) = 1; A1(i,6) = 1; A1(i,10) = -1; B1(i) = 1; signs(i) = 'U';

% delta 0,2 - 3xAND
i = i+1; A1(i, 4) = -1; A1(i,11) = 1; B1(i) = 0; signs(i) = 'U'; 
i = i+1; A1(i, 7) = -1; A1(i,11) = 1; B1(i) = 0; signs(i) = 'U'; 
i = i+1; A1(i, 8) = -1; A1(i,11) = 1; B1(i) = 0; signs(i) = 'U'; 
i = i+1; A1(i,4) = 1; A1(i,7) = 1; A1(i,8) = 1; A1(i,11) = -1; B1(i) = 2; signs(i) = 'U';

% delta 0,3 - 3xAND
i = i+1; A1(i, 4) = -1; A1(i,12) = 1; B1(i) = 0; signs(i) = 'U'; 
i = i+1; A1(i, 7) = -1; A1(i,12) = 1; B1(i) = 0; signs(i) = 'U'; 
i = i+1; A1(i, 9) = -1; A1(i,12) = 1; B1(i) = 0; signs(i) = 'U'; 
i = i+1; A1(i,4) = 1; A1(i,7) = 1; A1(i,9) = 1; A1(i,12) = -1; B1(i) = 2; signs(i) = 'U';

% delta M1
i = i+1; A1(i,3) = 1; A1(i,13) = Mu-1; B1(i) = Mu; signs(i) = 'U';
i = i+1; A1(i,3) = 1; A1(i, 13) = -(mu-1-eps); B1(i) = eps + 1; signs(i) = 'L';

% delta 1 - 2xAND
i = i+1; A1(i, 5) = -1; A1(i,14) = 1; B1(i) = 0; signs(i) = 'U'; 
i = i+1; A1(i, 13) = -1; A1(i,14) = 1; B1(i) = 0; signs(i) = 'U'; 
i = i+1; A1(i,5) = 1; A1(i,13) = 1; A1(i,14) = -1; B1(i) = 1; signs(i) = 'U';

% delta 1n - not
i = i+1; A1(i,14) = 1; A1(i, 15) = 1; B1(i) = 1; signs(i) = 'U';
i = i+1; A1(i,14) = -1; A1(i, 15) = -1; B1(i) = -1; signs(i) = 'U';

% delta1a
i = i+1; A1(i,1) = 1; A1(i, 16) = Mcon-x_eff; B1(i) = Mcon; signs(i) = 'U';
i = i+1; A1(i,1) = 1; A1(i,16) = -(mcon-x_eff-eps); B1(i) = eps + x_eff; signs(i) = 'L';

% delta 1an - not
i = i+1; A1(i,16) = 1; A1(i, 17) = 1; B1(i) = 1; signs(i) = 'U';
i = i+1; A1(i,16) = -1; A1(i, 17) = -1; B1(i) = -1; signs(i) = 'U';

% delta1b
i = i+1; A1(i,2) = 1; A1(i, 18) = Maux-N_max_gr + eps; B1(i) = Maux; signs(i) = 'U';
i = i+1; A1(i,2) = 1; A1(i,18) = -(mcon-N_max_gr-eps); B1(i) = N_max_gr; signs(i) = 'L';

% delta 1,1 - 3xAND
i = i+1; A1(i, 14) = -1; A1(i,19) = 1; B1(i) = 0; signs(i) = 'U'; 
i = i+1; A1(i, 16) = -1; A1(i,19) = 1; B1(i) = 0; signs(i) = 'U'; 
i = i+1; A1(i, 18) = -1; A1(i,19) = 1; B1(i) = 0; signs(i) = 'U'; 
i = i+1; A1(i,14) = 1; A1(i,16) = 1; A1(i,18) = 1; A1(i,19) = -1; B1(i) = 2; signs(i) = 'U';

% delta 1,2 - 3xAND
i = i+1; A1(i, 14) = -1; A1(i,20) = 1; B1(i) = 0; signs(i) = 'U'; 
i = i+1; A1(i, 17) = -1; A1(i,20) = 1; B1(i) = 0; signs(i) = 'U'; 
i = i+1; A1(i, 18) = -1; A1(i,20) = 1; B1(i) = 0; signs(i) = 'U'; 
i = i+1; A1(i,14) = 1; A1(i,17) = 1; A1(i,18) = 1; A1(i,20) = -1; B1(i) = 2; signs(i) = 'U';

% delta 1 - 2xAND
i = i+1; A1(i, 5) = -1; A1(i,21) = 1; B1(i) = 0; signs(i) = 'U'; 
i = i+1; A1(i, 15) = -1; A1(i,21) = 1; B1(i) = 0; signs(i) = 'U'; 
i = i+1; A1(i,5) = 1; A1(i,15) = 1; A1(i,21) = -1; B1(i) = 1; signs(i) = 'U';

% start z1
% 48
i = i+1;
A1(i,10) = -Mcon; A1(i,22) = 1; 
B1(i) = 0;
signs(i) = 'U';
% 49
i = i+1;
A1(i,10) = -mcon; A1(i,22) = 1;
B1(i) = 0;
signs(i) = 'L';
% 50
i = i+1;
A1(i,10) = -mcon; A1(i,1) = -1; A1(i,22) = 1; 
B1(i) = -mcon;
signs(i) = 'U';
% 51
i = i+1;
A1(i,10) = -Mcon; A1(i,1) = -1; A1(i,22) = 1;
B1(i) = -Mcon;
signs(i) = 'L';
%---end

% start z2
% 52
i = i+1;
A1(i,11) = -Mcon; A1(i,23) = 1;
B1(i) = 0;
signs(i) = 'U';
% 53 
i = i+1;
A1(i,11) = -mcon; A1(i,23) = 1;
B1(i) = 0;
signs(i) = 'L';
% 54
i = i+1;
A1(i,11) = -mcon; A1(i,1) = -1; A1(i,23) = 1;
B1(i) = -mcon;
signs(i) = 'U';
% 55 
i = i+1;
A1(i,11) = -Mcon; A1(i,1) = -1; A1(i,23) = 1;

B1(i) = -Mcon;
signs(i) = 'L';
%---end

% start z3
% 56
i = i+1;
A1(i,12) = -Mcon; A1(i,24) = 1;
B1(i) = 0;
signs(i) = 'U';
% 57 
i = i+1;
A1(i,12) = -mcon; A1(i,24) = 1;
B1(i) = 0;
signs(i) = 'L';
% 58
i = i+1;
A1(i,12) = -mcon; A1(i,1) = -1; A1(i,24) = 1;
B1(i) = -mcon;
signs(i) = 'U';
% 59 
i = i+1;
A1(i,12) = -Mcon; A1(i,1) = -1; A1(i,24) = 1;

B1(i) = -Mcon;
signs(i) = 'L';
%---end

% start z4
% 60
i = i+1;
A1(i,20) = -Mcon; A1(i,25) = 1;
B1(i) = 0;
signs(i) = 'U';
% 61 
i = i+1;
A1(i,20) = -mcon; A1(i,25) = 1;
B1(i) = 0;
signs(i) = 'L';
% 62
i = i+1;
A1(i,20) = -mcon; A1(i,1) = -1; A1(i,25) = 1;
B1(i) = -mcon;
signs(i) = 'U';
% 63
i = i+1;
A1(i,20) = -Mcon; A1(i,1) = -1; A1(i,25) = 1;
B1(i) = -Mcon;
signs(i) = 'L';
% start z5
% 64
i = i+1;
A1(i,21) = -Maux; A1(i,26) = 1;
B1(i) = 0;
signs(i) = 'U';
% 65
i = i+1;
A1(i,21) = -maux; A1(i,26) = 1;
B1(i) = 0;
signs(i) = 'L';
% 66
i = i+1;
A1(i,2) = -1; A1(i,21) = -maux; A1(i,26) = 1;
B1(i) = -maux;
signs(i) = 'U';
% 67
i = i+1;
A1(i,2) = -1; A1(i,21) = -Maux; A1(i,26) = 1;
B1(i) = -Maux;
signs(i) = 'L';
%---end
%% Inequalities - State solution definition

Acon = zeros(1,26); Aaux = zeros(1,26);
Acon(10) = b1; Acon(11) = b2; Acon(12) = b3; Acon(20) = -psi*x_eff;
Acon(22) = a1; Acon(23) = a2; Acon(24) = a3; Acon(25) = psi;

Aaux(2) = 1; Aaux(19) = 1; Aaux(20) = 1; Aaux(26) = -1;


A3con = zeros(Np, Nx); A3aux = zeros(Np, Nx);

% Define
for i=1:Np
    A3con(i,(i-1)*26+1:i*26+1) = [Acon -1];
    A3aux(i,(i-1)*26+1:i*26+2) = [Aaux 0 -1];
end

A3 = [A3con ; A3aux];

% Add u functionality (Forced action for real system)
B3 = zeros(size(A3, 1),1);

%{
A_uConstraints=[]; B_uConstraints = []; signs_uConstraints = [];
if size(u,1) ~= 1 % is u a vector?
    A_uConstraints = zeros(Np, 2+Np*25);
    B_uConstraints = zeros(Np,1);
    for i = 1:Np % across all time steps in Np, u = Np x 1
        A_uConstraints(i, 25*(i-1)+1:25*i) = [0 0 1 zeros(1,22)];
        B_uConstraints(i) = u(i);
    end
    signs_uConstraints = 'S'; signs_uConstraints = kron(ones(Np,1), signs_uConstraints);
%}
if u ~= -1
    A3 = [A3; 0 0 1 zeros(1,size(A3,2)-3)];
    B3 = [B3; u];
end

signs3=char(zeros(size(A3,1), 1)); signs3(:) = 'S';



%% Cost function

Cset = zeros(1,26);
Cset(24) = 1; % xk+1
Cset(2) = gamma0 * lambda; % mode 0
Cset(12) = gamma1 * lambda; % mode 1
Cset(19) = gamma2 * lambda; % mode 2
C = [0 0 kron(ones(1,Np), Cset)];

%% Initial Conditions

A_xconinit = [1 zeros(1,Nx-1)];
A_xauxinit = [0 1 zeros(1,Nx-2)];
A_init = [A_xconinit; A_xauxinit]; B_init = [xcon_init; xaux_init];
signsinit = ['S'; 'S'];


%% Boundary condition (x_max)
A_xconmax = zeros(1,26); A_xconmax(25) = 1; A_xconmax = [zeros(Np,2) kron(eye(Np), A_xconmax)];
B_xconmax = x_max; 
signs_xconmax = kron(ones(Np,1), 'U');

%% Together

% Constructing total delta/z definitions across Np
A1big = kron(eye(Np), A1);
A1big = [A1big zeros(size(A1big, 1), 2)];
B1 = kron(ones(Np,1), B1);
bigSigns = kron(ones(Np,1), signs);

% Building variable types across Np
variableType = ['C', 'C', 'I', 'B',  'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'C', 'C', 'C', 'C', 'C']';
variableTypeList = [];
for i=1:Np
    variableTypeList = [variableTypeList; variableType];
end
variableTypeList = [variableTypeList; 'C'; 'C'];


% In total then
Atotal = [A1big; A3; A_init; A_xconmax]; 
signsTotal = [bigSigns; signs3; signsinit; signs_xconmax];
Btotal = [B1; B3; B_init; B_xconmax];

% Upper Bounds / Lower Bounds (just for simulation)
UB = 100*ones(Nx, 1); LB = -100*ones(Nx,1);

[xopt, J_opt, MILP_status, MILP_extras] = glpk(C, Atotal, Btotal, LB, UB, signsTotal, variableTypeList, 1);
u_opt = xopt(3); x_next = xopt(27:28);
end