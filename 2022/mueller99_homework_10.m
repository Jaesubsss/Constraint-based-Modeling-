clear

load("e_coli_core.mat");
model = e_coli_core;

%Simplified OptKnock
idx_Fum = find(strcmp(model.rxns, "FUM"));
idx_ATP = find(strcmp(model.rxnNames, "ATP maintenance requirement"));
idx_biomass = find(strcmp(model.rxns, "BIOMASS_Ecoli_core_w_GAM"));

model.lb(idx_Fum) = 0;
model.lb(idx_ATP) = 0;
modelr = split_rev_reactions(model);
idx_Fum = find(strcmp(modelr.rxns, "FUM"));
idx_biomass = find(strcmp(modelr.rxns, "BIOMASS_Ecoli_core_w_GAM"));

%calculate optimal biomass flux z*
[Sol.x, Sol.f] = linprog(-modelr.c, [],[], modelr.S, modelr.b, modelr.lb, modelr.ub);
flux_bio = -Sol.f;

% FVA
a = idx_biomass;
    for j = 1:size(modelr.S,2)
        f = zeros(size(modelr.c));
        f(j) = 1;
        ub_bio = modelr.ub;
        lb_bio = modelr.lb;
        ub_bio(a) = flux_bio;
        lb_bio(a) = flux_bio;
        OPTIONS = optimset('linprog');
        OPTIONS.Display = 'off';
        %tilde is used to supress the output of values, so I just get the
        %variable for my desired variable j
        [~, vmax(j,1)] = linprog(-f, [], [], modelr.S, modelr.b, modelr.lb, modelr.ub, OPTIONS); %falsch
        [~, vmin(j,1)] = linprog(f, [], [], modelr.S, modelr.b, modelr.lb, modelr.ub, OPTIONS); %falsch
        [~, v_U(j,1)] = linprog(-f, [], [], modelr.S, modelr.b, lb_bio, ub_bio, OPTIONS); %=sol
        [~, v_L(j,1)] = linprog(f, [], [], modelr.S, modelr.b, lb_bio, ub_bio, OPTIONS);    %=sol
    end
    vmax = -vmax;
    v_U = -v_U;

%MILP

epsilon = 1e-6;
C = 0.001;

% calculate v_U' and v_L'
v_U_prime = v_U*(1-C) + vmax*C;
v_L_prime = v_L*(1-C) + vmin*C;

%find rev reactions
rev_match = zeros(size(modelr.rxns));
for i=1:length(modelr.rxns)
    if contains(modelr.rxns(i), "_forward")
        string = strrep(modelr.rxns{i}, "_forward", "_backward");
        idx = find(strcmp(modelr.rxns, string));
        rev_match(i, i) = -1;
        rev_match(i, idx) = -1;
    end
    if contains(modelr.rxns(i), "_backward")
        string = strrep(modelr.rxns{i}, "_backward", "_forward");
        idx = find(strcmp(modelr.rxns, string));
        rev_match(i, i) = -1;
        rev_match(i, idx) = -1;
    end
end

%remove zeros rows
rev_match(~any(rev_match,2),:) = [];

%matrix KO foward & backward
rev_M = 0; 
for i=1:length(modelr.rxns)
    if contains(modelr.rxns(i), "_forward")
        string = strrep(modelr.rxns{i}, "_forward", "_backward");
        idx = find(strcmp(modelr.rxns, string));
        rev_M(i, i) = 1;
        rev_M(i, idx) = -1;
    end
    if contains(modelr.rxns(i), "_backward")
        string = strrep(modelr.rxns{i}, "_backward", "_forward");
        idx = find(strcmp(modelr.rxns, string));
        rev_M(i, i) = 1;
        rev_M(i, idx) = -1;
    end
end
rev_M(~any(rev_M,2),:) = [];

I = eye(length(modelr.rxns));
Z = zeros(length(modelr.rxns));

A_k = [-I, epsilon*I, Z, Z;
       I, (-vmax + epsilon).*I, Z, Z];
A_d = [-I, Z, (v_L - vmin).*I, Z;
       I, Z, (v_L_prime -vmax).*I, Z];
A_u = [-I, I*0, I*0, (vmin-v_U_prime).*I;
       I, Z, Z, (vmax -v_U).*I];
%include Task 3 and task 4 in A
Z1 = zeros(size(rev_match));
A = [Z, -I, -I, -I;
    Z1, Z1, Z1, rev_match;
    Z1, Z1, rev_match, Z1;
    Z1, rev_M, Z1, Z1];
Aineq = [A_k; A_d; A_u; A];

one = ones(length(modelr.rxns),1);
b_k = [0*one; epsilon.*one];
b_d = [-vmin.*one; v_L_prime.*one];
b_u = [-v_U_prime.*one; vmax.*one];
%include task 3 and 4 in b
b = [-one; -ones(size(rev_match,1),1); -ones(size(rev_match,1),1); zeros(size(rev_M,1),1)];
bineq = [b_k;b_d;b_u; b];


%MILP:

Aeq = [modelr.S zeros(size(modelr.mets,1), 3*size(modelr.rxns,1))];
beq = modelr.b;

r = size(modelr.S,2);
intcon = r+1:4*r;

f = zeros(size(modelr.c));
f(r+1:4*r) = -1;

modelr.lb(idx_biomass) = 0.9*flux_bio;
modelr.lb(idx_Fum) = 1.7*v_U(idx_Fum);
%modelr.ub(idx_Fum) = 1000;
z = zeros(size(modelr.ub));
o = ones(size(modelr.ub));
ub = [modelr.ub; o; o; o];
lb = [modelr.lb; z; z; z];

abc = intlinprog(f,intcon, Aineq, bineq, Aeq, beq, lb, ub,[]);

%Output Table:
T = table('Size', [length(modelr.rxns), 2], 'VariableTypes', {'string', 'string'}, 'VariableNames', {'Reaction', 'Status'});

for i = 1:length(modelr.rxns)
    if abc(2*i) == 0
        T.Reaction(i) = modelr.rxns{i};
        T.Status(i) = "Knocked-Out";
    elseif abc(3*i) == 0
        T.Reaction(i) = modelr.rxns{i};
        T.Status(i) = "Downregulation";
    elseif abc(4*i) == 0
        T.Reaction(i) = modelr.rxns{i};
        T.Status(i) = "Upregulation";
    end
end

% Remove rows where both columns are empty
T(any(ismissing(T), 2), :) = [];

disp(T)

%Task: Find alternative solution using integer cut

function model_reduced = split_rev_reactions(model)

model_reduced.S = model.S;
model_reduced.c = model.c;
model_reduced.b = model.b;
model_reduced.lb = model.lb;
model_reduced.ub = model.ub;
model_reduced.mets = model.mets;
model_reduced.rxns = model.rxns;

 % For your model model_reduced, change the model such that each reversible reaction is splitted into two irreversible reactions. 
rev_rxn_set = find(model_reduced.lb<0 & model_reduced.ub>0);

model_reduced.S(:,end+1:end+length(rev_rxn_set)) = model_reduced.S(:,rev_rxn_set)*-1; % add columns corresponding to backward direction

% Update the other model fields accordingly, such that dimensions of all fields are correct. 
% update other fields
model_reduced.c(end+1:end+length(rev_rxn_set)) = 0;
model_reduced.lb(end+1:end+length(rev_rxn_set)) = 0;
model_reduced.ub(end+1:end+length(rev_rxn_set)) = model_reduced.lb(rev_rxn_set)*-1; % lower bound becomes upper bound of backward rxn
model_reduced.lb(rev_rxn_set) = 0; % forward direction now has lower bound 0
model_reduced.rxns(end+1:end+length(rev_rxn_set)) = strcat(model_reduced.rxns(rev_rxn_set), '_backward'); % add 'backward' and 'forward' to reaction name
model_reduced.rxns(rev_rxn_set) = strcat(model_reduced.rxns(rev_rxn_set), '_forward');

end


