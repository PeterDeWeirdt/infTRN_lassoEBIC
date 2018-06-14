%% EBIC + Confidence Score example
% Use the Extended Bayesian Information Criterion (EBIC) 
% To select the best parameter value(s) for one gene.
% Then, obtain a confidence score for the nonzero predictors in the optimal
% model.
%% Author: Peter DeWeirdt 
% Summer Intern Cincinnati Children's Hospital
%% Date: 6/12/2018
clear all
close all
restoredefaultpath

matlabDir = '..';

addpath(fullfile(matlabDir,'glmnet'))
addpath(fullfile(matlabDir,'ebicFxns'))

%%
n_preds = 300; % > 3
n_samps = 100;
lambdaBias = 0.5;

% set up vector of potential values of lambda

lambdaMin = .01;
lambdaMax = 1;
totLogLambdaSteps = 25;
lamLog10step = 1/totLogLambdaSteps;
logLamRange =  log10(lambdaMin):lamLog10step:log10(lambdaMax);
lambdaRange = (10.^logLamRange)';

% Use glmnet to fit a weights coefficient for each parameter
options = glmnetSet;
options.lambda = fliplr(lambdaRange'); % must flip for glmnet
currOptions = options;
prior_weights = ones(1, n_preds) - [0.5 0.5 0.5 zeros(1, n_preds - 3)];
%prior_weights = priorWeightsMat(res,predInds)'
currOptions.penalty_factor = prior_weights; 
fake_beta = [1 2 10 100 zeros(1, n_preds - 4)]';
currPreds = rand(n_samps, n_preds);
currResponses = currPreds*fake_beta + rand(n_samps,1);
lsoln = glmnet(currPreds,currResponses,'gaussian',currOptions);

% Use output of GLMNET to then check EBIC value
currBetas = fliplr(lsoln.beta); % flip so that the lambdas are increasing
currA0 = fliplr(lsoln.a0');
RSSs = sum(bsxfun(@minus,(bsxfun(@plus, currPreds*currBetas, currA0)),...
    currResponses).^2, 1); 
n_samps = size(currResponses, 1);
n_nonzero = fliplr(lsoln.df');
tot_preds = size(currBetas, 1);
gamma = 1;
[EBIC, optInd] = ebic(RSSs, n_samps, n_nonzero, tot_preds, gamma);

% Plot EBIC values
s = scatter(lambdaRange, EBIC, 36,  n_nonzero);
s.MarkerFaceColor = 'flat';
set(gca, 'XScale', 'log');
xlabel('labmda');
ylabel('EBIC');
colormap(jet(max(n_nonzero + 1)));
c = colorbar;
ylabel(c, 'num. nonzero');

% Get a confidence score and ranking for each entry in the optimal beta
optBeta = currBetas(:,optInd);
optA0 = currA0(:,optInd);
[confs, rnk] = calc_confs(optBeta, optA0, currPreds, currResponses);


