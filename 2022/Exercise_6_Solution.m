clear
load('e_coli_core.mat')

%% 1. Remove blocked reactions and convert to irreversible
% multiple options:
% a) use function from coba model
    blk = findBlockedReaction(e_coli_core);
% b) use FVA and check iof there are reactions i for which min_i=0 and max_i=0
%    instead of FVA we can also use function FlxRange from previous
%    exercise
%     [mini,maxi] = FluxRange(e_coli_core,1:length(e_coli_core.rxns));
%     [mini,maxi] = fluxVariability(e_coli_core);
%     blk = find(mini==0 & maxi==0);

% remove blk reactions
e_coli_core = removeRxns(e_coli_core,blk);

% split reversible reactions
e_coli_core_irrev = convertToIrreversible(e_coli_core);

%% 2.	Manipulating growth conditions
%
% -	Which metabolites are imported?

import_rxns = find(all(e_coli_core_irrev.S>=0));
[import_mets,~] = find(e_coli_core_irrev.S(:,import_rxns)~=0);

disp('Import metabolites (default):');...
    disp(e_coli_core_irrev.mets(import_mets))

% -	Is the model simulating aerobic or anaerobic growth? What is the carbon source?

disp('The model simulates aerobic growth with glucose and CO2 as carbon sources.')

% -	Change the model such that it simulates anaerobic growth with glucose as only carbon source.

idx_O2 = import_rxns(find(strcmp(e_coli_core_irrev.mets(import_mets),'o2_e')));
idx_CO2 = import_rxns(find(strcmp(e_coli_core_irrev.mets(import_mets),'co2_e')));

e_coli_core_irrev_anaerobic = e_coli_core_irrev;

% block import of CO2 and O2
e_coli_core_irrev_anaerobic.lb([idx_O2 idx_CO2]) = 0;
e_coli_core_irrev_anaerobic.ub([idx_O2 idx_CO2]) = 0;
%
%% 3. Calculate the optimal growth rate under aerobic and anaerobic growth

Sol_aerobic = optimizeCbModel(e_coli_core_irrev);
Sol_anaerobic = optimizeCbModel(e_coli_core_irrev_anaerobic);

disp('optimal growth anaerobic:')
Sol_anaerobic.f
disp('optimal growth aerobic:')
Sol_aerobic.f

%
%% 4. Write MATLAB code to classify reactions based on their relation to optimal specific growth rate using FVA under the aerobic and aerobic growth conditions.
% check that biomass is the objective
e_coli_core_irrev.rxns(find(e_coli_core_irrev.c~=0))
e_coli_core_irrev_anaerobic.rxns(find(e_coli_core_irrev_anaerobic.c~=0))

% the lower bound of biomass should be 90%
% fluxVariability can take a second argument that specifies how much of
% optimal biomass must be guaranteed
[mini_aerobic,maxi_aerobic] = fluxVariability(e_coli_core_irrev,0.9);
[mini_anaerobic,maxi_anaerobic] = fluxVariability(e_coli_core_irrev_anaerobic,0.9);

% create vectors indicating the coupling type
% 1: hard-coupled
% 2: partially coupled
% 3: not coupled
% 4: blocked
% set all partially coupled and refill with other coupling types
class_aerobic = ones(length(e_coli_core_irrev_nob.rxns),1)*2;
class_anaerobic = ones(length(e_coli_core_irrev_anaerobic_nob.rxns),1)*2;

% -	Hard-coupled to objective if its flux varies exactly like objective
class_aerobic(intersect(find(mini_aerobic==(Sol_aerobic_nob.f*0.9)),find(maxi_aerobic==Sol_aerobic_nob.f))) = 1;
class_anaerobic(intersect(find(mini_anaerobic==(Sol_anaerobic_nob.f*0.9)),find(maxi_anaerobic==Sol_anaerobic_nob.f))) = 1;

% -	Partially-coupled to objective if its flux is non-zero, but can vary
% = 2

% -	Not coupled to objective if its flux can take a zero value
class_aerobic(find(mini_aerobic==0)) = 3;
class_anaerobic(find(mini_anaerobic==0)) = 3;

% -	Blocked at objective if its flux takes only a value of zero
class_aerobic(intersect(find(mini_aerobic==0),find(maxi_aerobic==0))) = 4;
class_anaerobic(intersect(find(mini_anaerobic==0),find(maxi_anaerobic==0))) = 4;

%% Find reactions whose flux is highly correlated to biomass flux under aerobic and anaerobic conditions.

% sample biomass values from range [0.1*z z], where z is optimum biomass
rand_bio_ae = (0.1*Sol_aerobic.f) + (Sol_aerobic.f - (0.1*Sol_aerobic.f))*rand(100,1);
rand_bio_anae = (0.1*Sol_anaerobic.f) + (Sol_anaerobic.f - (0.1*Sol_anaerobic.f))*rand(100,1);

idx_bio = find(e_coli_core_irrev.c~=0);

% minimize sum over all fluxes
e_coli_core_irrev.c(:) = 1;

for i=1:length(rand_bio)

    e_coli_core_irrev.lb(idx_bio) = rand_bio(i);
    e_coli_core_irrev.ub(idx_bio) = rand_bio(i);

    sample_v(:,i) = optimizeCbModel(e_coli_core_irrev,'min').v;
end

% calculate correlation between reaction i and biomass
for j=1:size(sample_v,1)
    C(j)=corr(sample_v(j,:),rand_bio);
end

[val,rxn_idx]=sort(C,'descend');
