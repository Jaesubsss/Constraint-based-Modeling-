%% Exercise 10 solution

% read the model
model = readCbModel('e_coli_core.mat');

% predict optimal growth rate
fba_solution = optimizeCbModel(model);
biomass_opt = fba_solution.f;

% determine reaction weights
w = getRxnWeights(model);

% write LP for pFBA
[METS,RXNS] = size(model.S);
problem = struct('options',optimoptions('linprog'),'solver','linprog');
problem.Aeq = [model.S zeros(METS,2*RXNS); -eye(RXNS) eye(RXNS) -eye(RXNS)];
problem.beq = [model.b; zeros(RXNS,1)];
problem.A = [];
problem.b = [];
problem.lb = [model.lb; zeros(2*RXNS,1)];
problem.lb(find(model.c)) = biomass_opt;
problem.ub = [model.ub; 1000*ones(2*RXNS,1)];
problem.ub(find(model.c)) = biomass_opt;
problem.f = [zeros(RXNS,1); repmat(w,2,1)];

[x,f] = linprog(problem);

flux_values = x(1:RXNS);

function w = getRxnWeights(model)
% Determines reaction weights by the number of unique genes associated with
% a reaction in the model.
% INPUT
%   struct model:   metabolic model
% OUTPUT
%   double w:       reaction weights
w = cellfun(@(x)numel(unique(regexp(x,'x\(\d+\)','match'))),model.rules);
end
