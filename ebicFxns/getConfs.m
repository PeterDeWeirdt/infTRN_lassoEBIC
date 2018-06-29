function confs_rnks_sgn = getConfs(X, Ys, Prior, options)
%% Goal: For each column of Y calculate a confidence and rank
%% Inputs: 
% X -- Samples by predictors matrix
% Ys -- Response matrix, with each column corresponding to a response
% Prior -- Our prior matrix
% options -- Options for GLMNET
%% References:
% (1) Qian et al. (2013) "Glmnet for Matlab."
% http://www.stanford.edu/~hastie/glmnet_matlab/
% (2) Greenfield et al. (2013) "Robust data-driven incorporation of prior
% knowledge into the inference of dynamic regulatory netowrks"
% (3) Chen et al. (2008) "Extended bayesian information criteria for model
% selection with large model spaces"
% (4) Castro, De Veaux, Miraldi, Bonneau "Multitask learning for joint
%   inference of gene regulatory networks form several expression datasets"
%% Output:
% confs_rnk_sgn -- confidences in the first cell, then ranks, then sign of
% interactions. 
%% Author:
% Peter Deweirdt, Summer Intern, Divisions of Immunobiology and Biomedical
%   Informatics, Cincinnati Children's Hospital
%% Date: 6/15/2018
nResp = size(Ys, 1);
nTFs = size(X, 1);
confs = zeros(nResp, nTFs);
rnks = confs; 
sgn = confs;
for i = 1:nResp
    % get (finite) predictor indices for each response // filter
    currWeights = Prior(i,:);
    % limit to predictors with finite lambda penalty (e.g., to exclude TF mRNA self-interaction loops)
    predInds = find(isfinite(currWeights)); 
    currPreds = zscore(X(predInds,:)');
    currResp = zscore(Ys(i,:)');
    currOpts = options;
    currOpts.penalty_factor = Prior(i,predInds)';    
    lsoln = glmnet(currPreds,currResp,'gaussian',currOpts);
    % Use output of GLMNET to check fit
    currBetas = fliplr(lsoln.beta); % flip so that the lambdas are increasing
    currA0 = fliplr(lsoln.a0');
    [confs_i, rnk_i] = calc_confs(currBetas, currA0, currPreds, currResp);
    e_confs_i = zeros(1, nTFs);
    e_rnks_i = repelem(max(rnk_i), nTFs);
    e_signs_i = e_confs_i;
    e_confs_i(predInds) = confs_i;
    e_rnks_i(predInds) = rnk_i;
    e_signs_i(predInds) = sign(currBetas);
    confs(i,:) = e_confs_i;
    rnks(i,:) = e_rnks_i;
    sgn(i,:) = e_signs_i;
end
confs_rnks_sgn = {confs, rnks, sgn};
