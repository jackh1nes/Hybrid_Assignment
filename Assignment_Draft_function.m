function [u_opt, J_opt, MILP_status, MILP_extras] = Assignment_Draft_function(Np, xcon_init, xaux_init, Parameters)
%% Expand parameters
x_max = Parameters.x_max;
N_max_gr = Parameters.N_max_gr;
gamma0 = Parameters.gamma0;
gamma1 = Parameters.gamma1; 
gamma2 = Parameters.gamma2;

x_eff = Parameters.x_eff;
psi = Parameters.psi;
a1 = Parameters.A(1); a2 = Parameters.A(2); a3 = Parameters.A(3);
b1 = Parameters.B(1); b2 = Parameters.B(2); b3 = Parameters.B(3);

Mu = 2; mu = 0; Mcon = 70; mcon = 0; Maux = N_max_gr; maux = 0; eps = 10^-3;

Nx = Np*25+2;

%% Inequalities - delta and z definitions

A1 = zeros(63, 25); B1 = zeros(63,1); signs = char(zeros(63,1)); 

% delta0
i = 1; A1(i,3) = 1; A1(i,4) = Mu; B1(i) = Mu; signs(i) = 'U'; % <=
i = 2; A1(i,3) = 1; A1(i, 4) = -(mu-eps); B1(i) = eps; signs(i) = 'L'; % >=

% delta 0n - NOT
i = 3; A1(i,4) = 1; A1(i, 5) = 1; B1(i) = 1; signs(i) = 'U';
i = 4; A1(i,4) = -1; A1(i, 5) = -1; B1(i) = -1; signs(i) = 'U';

% delta0a
i = 5; A1(i,1) = 1; A1(i, 6) = Mcon-30 + eps; B1(i) = Mcon; signs(i) = 'U';
i = 6; A1(i,1) = 1; A1(i,6) = -(mcon-30-eps); B1(i) = 30; signs(i) = 'L';

% delta 0an - NOT
i = 7; A1(i,6) = 1; A1(i, 7) = 1; B1(i) = 1; signs(i) = 'U';
i = 8; A1(i,6) = -1; A1(i, 7) = -1; B1(i) = -1; signs(i) = 'U';

% delta0b
i = 9; A1(i,1) = 1; A1(i, 8) = Mcon-50 + eps; B1(i) = Mcon; signs(i) = 'U';
i = 10; A1(i,1) = 1; A1(i,8) = -(mcon-50-eps); B1(i) = 50; signs(i) = 'L';

% delta 0bn - NOT
i = 11; A1(i,8) = 1; A1(i, 9) = 1; B1(i) = 1; signs(i) = 'U';
i = 12; A1(i,8) = -1; A1(i, 9) = -1; B1(i) = -1; signs(i) = 'U';

% delta 0,1 - 2xAND
i = 13; A1(i, 4) = -1; A1(i,10) = 1; B1(i) = 0; signs(i) = 'U'; 
i = 14; A1(i, 6) = -1; A1(i,10) = 1; B1(i) = 0; signs(i) = 'U'; 
i = 15; A1(i,4) = 1; A1(i,6) = 1; A1(i,10) = -1; B1(i) = 1; signs(i) = 'U';

% delta 0,2 - 3xAND
i = 16; A1(i, 4) = -1; A1(i,11) = 1; B1(i) = 0; signs(i) = 'U'; 
i = 17; A1(i, 7) = -1; A1(i,11) = 1; B1(i) = 0; signs(i) = 'U'; 
i = 18; A1(i, 8) = -1; A1(i,11) = 1; B1(i) = 0; signs(i) = 'U'; 
i = 19; A1(i,4) = 1; A1(i,7) = 1; A1(i,8) = 1; A1(i,11) = -1; B1(i) = 2; signs(i) = 'U';

% delta 0,3 - 3xAND
i = 20; A1(i, 4) = -1; A1(i,12) = 1; B1(i) = 0; signs(i) = 'U'; 
i = 21; A1(i, 7) = -1; A1(i,12) = 1; B1(i) = 0; signs(i) = 'U'; 
i = 22; A1(i, 9) = -1; A1(i,12) = 1; B1(i) = 0; signs(i) = 'U'; 
i = 23; A1(i,4) = 1; A1(i,7) = 1; A1(i,9) = 1; A1(i,12) = -1; B1(i) = 2; signs(i) = 'U';

% delta M1
i = 24; A1(i,3) = 1; A1(i,13) = Mu-1; B1(i) = Mu; signs(i) = 'U';
i = 25; A1(i,3) = 1; A1(i, 13) = -(mu-1-eps); B1(i) = eps + 1; signs(i) = 'L';

% delta 1 - 2xAND
i = 26; A1(i, 5) = -1; A1(i,14) = 1; B1(i) = 0; signs(i) = 'U'; 
i = 27; A1(i, 13) = -1; A1(i,14) = 1; B1(i) = 0; signs(i) = 'U'; 
i = 28; A1(i,5) = 1; A1(i,13) = 1; A1(i,14) = -1; B1(i) = 1; signs(i) = 'U';

% delta 1n - not
i = 29; A1(i,14) = 1; A1(i, 15) = 1; B1(i) = 1; signs(i) = 'U';
i = 30; A1(i,14) = -1; A1(i, 15) = -1; B1(i) = -1; signs(i) = 'U';

% delta1a
i = 31; A1(i,1) = 1; A1(i, 16) = Mcon-x_eff; B1(i) = Mcon; signs(i) = 'U';
i = 32; A1(i,1) = 1; A1(i,16) = -(mcon-x_eff-eps); B1(i) = eps + x_eff; signs(i) = 'L';

% delta 1an - not
i = 33; A1(i,16) = 1; A1(i, 17) = 1; B1(i) = 1; signs(i) = 'U';
i = 34; A1(i,16) = -1; A1(i, 17) = -1; B1(i) = -1; signs(i) = 'U';

% delta1b
i = 35; A1(i,2) = 1; A1(i, 18) = Maux-N_max_gr + eps; B1(i) = Maux; signs(i) = 'U';
i = 36; A1(i,2) = 1; A1(i,18) = -(mcon-N_max_gr-eps); B1(i) = N_max_gr; signs(i) = 'L';

% delta 1,1 - 3xAND
i = 37; A1(i, 14) = -1; A1(i,19) = 1; B1(i) = 0; signs(i) = 'U'; 
i = 38; A1(i, 16) = -1; A1(i,19) = 1; B1(i) = 0; signs(i) = 'U'; 
i = 39; A1(i, 18) = -1; A1(i,19) = 1; B1(i) = 0; signs(i) = 'U'; 
i = 40; A1(i,14) = 1; A1(i,16) = 1; A1(i,18) = 1; A1(i,19) = -1; B1(i) = 2; signs(i) = 'U';

% delta 1,2 - 3xAND
i = 41; A1(i, 14) = -1; A1(i,20) = 1; B1(i) = 0; signs(i) = 'U'; 
i = 42; A1(i, 17) = -1; A1(i,20) = 1; B1(i) = 0; signs(i) = 'U'; 
i = 43; A1(i, 18) = -1; A1(i,20) = 1; B1(i) = 0; signs(i) = 'U'; 
i = 44; A1(i,14) = 1; A1(i,17) = 1; A1(i,18) = 1; A1(i,20) = -1; B1(i) = 2; signs(i) = 'U';

% delta 1 - 2xAND
i = 45; A1(i, 5) = -1; A1(i,21) = 1; B1(i) = 0; signs(i) = 'U'; 
i = 46; A1(i, 15) = -1; A1(i,21) = 1; B1(i) = 0; signs(i) = 'U'; 
i = 47; A1(i,5) = 1; A1(i,15) = 1; A1(i,21) = -1; B1(i) = 1; signs(i) = 'U';


% z1
signs(48) = 'U'; signs(49) = 'L'; signs(50) = 'U'; signs(51) = 'L';
% z2
signs(52) = 'U'; signs(53) = 'L'; signs(54) = 'U'; signs(55) = 'L';
% z3
signs(56) = 'U'; signs(57) = 'L'; signs(58) = 'U'; signs(59) = 'L';
% z4
signs(60) = 'U'; signs(61) = 'L'; signs(62) = 'U'; signs(63) = 'L';

% start z1
% 48

A1(48,10) = -Mcon; A1(48,22) = 1;
B1(48) = 0;

% 49

A1(49,10) = -mcon; A1(49,22) = 1;
B1(49) = 0;

% 50

A1(50,10) = -mcon; A1(50,1) = -1; A1(50,22) = 1;
B1(50) = -mcon;

% 51

A1(51,10) = -Mcon; A1(51,1) = -1; A1(51,22) = 1;

B1(51) = -Mcon;
%---end

% start z2
% 52

A1(52,11) = -Mcon; A1(52,23) = 1;
B1(52) = 0;

% 53 

A1(53,11) = -mcon; A1(53,23) = 1;
B1(53) = 0;

% 54

A1(54,11) = -mcon; A1(54,1) = -1; A1(54,23) = 1;
B1(54) = -mcon;

% 55 

A1(55,11) = -Mcon; A1(55,1) = -1; A1(55,23) = 1;

B1(55) = -Mcon;
%---end

% start z3
% 56

A1(56,12) = -Mcon; A1(56,24) = 1;
B1(56) = 0;

% 57 

A1(57,12) = -mcon; A1(57,24) = 1;
B1(57) = 0;

% 58

A1(58,12) = -mcon; A1(58,1) = -1; A1(58,24) = 1;
B1(58) = -mcon;

% 59 

A1(59,12) = -Mcon; A1(59,1) = -1; A1(59,24) = 1;

B1(59) = -Mcon;
%---end

% start z4
% 60

A1(60,20) = -Mcon; A1(60,25) = 1;
B1(60) = 0;

% 61 

A1(61,20) = -mcon; A1(61,25) = 1;
B1(61) = 0;

% 62

A1(62,20) = -mcon; A1(62,1) = -1; A1(62,25) = 1;
B1(62) = -mcon;

% 63

A1(63,20) = -Mcon; A1(63,1) = -1; A1(63,25) = 1;
B1(63) = -Mcon;
%---end
%% Inequalities - State solution definition

Acon = zeros(1,25); Aaux = zeros(1,25);
Acon(10) = b1; Acon(11) = b2; Acon(12) = b3; Acon(20) = -psi*x_eff;
Acon(22) = a1; Acon(23) = a2; Acon(24) = a3; Acon(25) = psi;

Aaux(2) = 1; Aaux(19) = 1; Aaux(20) = 1;

signs3=char(zeros(Np*2, 1)); signs3(:) = 'S';


A3con = zeros(Np, Nx); A3aux = zeros(Np, Nx);

% Define
for i=1:Np
    A3con(i,(i-1)*25+1:i*25+1) = [Acon -1];
    A3aux(i,(i-1)*25+1:i*25+2) = [Aaux 0 -1];
end
A3 = [A3con ; A3aux];
B3 = zeros(size(A3, 1),1);
%% Cost function

% initialize params
lambda = 700;

C = zeros(1,Nx);

for i = 1:7
    for j = 1:25
        if ((i-1)*25+j)<153
            if (i>1) && (j==1 || j==2)
                C((i-1)*25+j) = 1;
            elseif (j == 19) || (j == 20) 
                C((i-1)*25+j) = gamma1*lambda;
            elseif (j == 21)
                C((i-1)*25+j) = gamma2*lambda;
            end
        end
    end
end

%% Initial Conditions

A_xconinit = [1 zeros(1,Nx-1)];
A_xauxinit = [0 1 zeros(1,Nx-2)];
A_init = [A_xconinit; A_xauxinit]; B_init = [xcon_init; xaux_init];
signsinit = ['S'; 'S'];
%% Together

% Constructing total delta/z definitions across Np
A1big = kron(eye(Np), A1);
A1big = [A1big zeros(size(A1big, 1), 2)];
B1 = kron(ones(Np,1), B1);
bigSigns = kron(ones(Np,1), signs);

% Building variable types across Np
variableType = ['C', 'C', 'I', 'B',  'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'B', 'C', 'C', 'C', 'C']';
variableTypeList = [];
for i=1:Np
    variableTypeList = [variableTypeList; variableType];
end
variableTypeList = [variableTypeList; 'C'; 'C'];


% In total then
Atotal = [A1big; A3; A_init]; 
signsTotal = [bigSigns; signs3; signsinit];
Btotal = [B1; B3; B_init];

% Upper Bounds / Lower Bounds (just for simulation)
UB = 100*ones(Nx, 1); LB = -100*ones(Nx,1);

[xopt, J_opt, MILP_status, MILP_extras] = glpk(C, Atotal, Btotal, LB, UB, signsTotal, variableTypeList, 1);
u_opt = xopt(3);

end