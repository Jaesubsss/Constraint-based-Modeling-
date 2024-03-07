%% TASK 1 Duality
%%%%%%%%%%%%%%%%%%
clear; clc;
% initCobraToolbox(false)
%
% a) Write the dual problem for the following optimization problem
%
% primal LP:
%  min C = 0.5v1 + 0.3v2 + 0.2v3
%  s.t.
%  v1 + v3 >= 5
%  v2 + v3 >= 10
%  v3 <= 2
%  v >= 0
%
% dual LP:
%  max D = 5y1 + 10y2 + 2y3
%  s.t.
%  y1 <= 0.5
%  y2 <= 0.3
%  y1 + y2 + y3 <= 0.2
%  y1,y2 >= 0
%  y3 <= 0
%
% b) Use Matlab to check that both problems (primal and dual) have the same optimum.

% primal
c = [0.5 0.3 0.2];
A = [-1  0 -1; % minus one to change from >= to <=
      0 -1 -1; % minus one to change from >= to <=
      0  0  1];
b = [-5; -10; 2]; % -5 and -10 to change from >= to <=
lb = [0 0 0];
ub = [1000 1000 1000];
[~,f_primal] = linprog(c, A, b, [], [], lb, ub)

% dual
c_dual = -[5 10 2];
A_dual = [1 0 0;
          0 1 0;
          1 1 1];
b_dual = [0.5; 0.3; 0.2]; % -5 and -10 to change from >= to <=
lb_dual = [0 0 -1000];
ub_dual = [1000 1000 0];
[sp,f_dual] = linprog(c_dual, A_dual, b_dual, [], [], lb_dual, ub_dual)
f_dual = -f_dual;

disp('optimal value primal and dual LP is the same');...
disp(f_primal == f_dual)

% c) Calculate the shadow prices for the optimization problem above.
% The shadow price is given by the values of the variables in the dual LP
% Using optimizeCbModel you find the shadow price in the solution field .y

% d) Given a metabolic model with all irreversible reactions, the objective is to maximize the flux
% towards the biomass reaction when the model is at the steady state. Compare the number
% of variables and constraints in the primal problem and in the dual problem.
%
% If the model has no reversible reactions the number of variables and constraints in the
% primal LP is given by the number of columns/ rows in the stoichiometric matrix
% N of the model.
% Given that N has dimension m x n, we have m constraints in the primal and
% n variables. The corresponding dual LP has n constraints and m variables.
%

%% Simplified OptKnock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
% initCobraToolbox(false)
% read the E. coli core model
model = readCbModel('e_coli_core.mat');

% Set the lower bound of ‘Fumarase’ (FUM) and ‘ATP maintenance requirement’ reaction to zero.
ngam_idx = find(strcmp(model.rxnNames, 'ATP maintenance requirement'));
model.lb(ngam_idx) = 0;
fum_idx = find(strcmp(model.rxns, 'FUM'));
model.lb(fum_idx)=0;

% Calculate the optimum biomass flux z^*, under the steady state constraints.
fba_sol = optimizeCbModel(model);
bio_opt = fba_sol.f;

% Calculate the optimum flux through fumarase reaction w_FUM, under steady state constraints at the optimum biomass.
biomass_idx = find(strcmp(model.rxns, 'BIOMASS_Ecoli_core_w_GAM'));
model.lb(biomass_idx)=bio_opt;
model.ub(biomass_idx)=bio_opt;
model.c(:)=0;
model.c(fum_idx)=1;
fum_opt = optimizeCbModel(model).f;

% Check if reversible reactions in the model have to be replaced by two irreversible reactions to solve the MILP below.
model = convertToIrreversible(model);
fum_idx = find(strcmp(model.rxns, 'FUM'));
biomass_idx = find(strcmp(model.rxns, 'BIOMASS_Ecoli_core_w_GAM'));

% The problem is to find the minimum number of eliminations such that while
% having at least 90% of the maximum biomass in the E. coli model, flux through
% fumarase reaction should be increased by 50%.

% Next, we construct the mixed-integer linear problem.
RXNS = size(model.S,2);
epsilon = 1e-6;

problem = struct;

% inequality constraints
problem.Aineq = [
%      v            y
    -eye(RXNS)  -epsilon*eye(RXNS); % (1-yi)eps <= vi
    eye(RXNS)   diag(model.ub-epsilon) % vi <= (1-yi)vmax+yi*eps
];

problem.bineq = [repmat(-epsilon,RXNS,1); model.ub];

% equality constraints
problem.Aeq = [
%      v            y
    model.S zeros(size(model.S))
];
problem.beq = model.b;

% lower and upper bounds
problem.lb = [model.lb; zeros(RXNS,1)];
problem.ub = [model.ub; ones(RXNS,1)];

% set lower bounds for biomass and citrate synthase reaction
problem.lb(biomass_idx) = 0.7*bio_opt;
problem.lb(fum_idx) = 1.5*fum_opt;

% define the objective: minimize the sum of y
problem.f = [zeros(RXNS,1); ones(RXNS,1)];

% Finally define integer variables
problem.intcon = RXNS+1:2*RXNS;

% disable display
problem.solver = 'intlinprog';
problem.options = optimoptions('intlinprog');
problem.options.Display = 'off';

[OK_SOL,OK_FVAL,EXIT_FLAG] = intlinprog(problem);

if EXIT_FLAG == 1
    ko_rxn_idx = round(OK_SOL(RXNS+1:end),4) > 0;
    fprintf('%d reactions must be knocked out:\n',sum(ko_rxn_idx))
    disp(model.rxnNames(ko_rxn_idx))
else
    warning('OptKnock could not find a solution')
end


