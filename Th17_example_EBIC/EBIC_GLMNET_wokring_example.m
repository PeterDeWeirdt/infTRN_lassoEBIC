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
n_preds = 500; % > 3
n_samps = 200;
n_genes = 500;
lambdaBias = 0.5;

% set up vector of potential values of lambda

lambdaMin = .001;
lambdaMax = 1;
totLogLambdaSteps = 5;
lamLog10step = 1/totLogLambdaSteps;
logLamRange =  log10(lambdaMin):lamLog10step:log10(lambdaMax);
lambdaRange = (10.^logLamRange)';

% Use glmnet to fit a weights coefficient for each parameter
options = glmnetSet;
options.lambda = fliplr(lambdaRange'); % must flip for glmnet
Prior = ones(n_genes, n_preds) - [rand(n_genes, 2), zeros(n_genes, n_preds - 2)];
%prior_weights = priorWeightsMat(res,predInds)'
fake_beta = [rand(3,n_genes); zeros(n_preds-3,n_genes)];
Preds = rand(n_samps, n_preds);
Responses = Preds*fake_beta + 0.1*rand(n_samps,n_genes);
method = 'ebic';
[test1, test2] = getFit(Preds', Responses', Prior, options, method);
lambdaFit = outFit{2};
optInds = outFit{3};
histogram(optInds); %in this example, a significant portion are trying to be more sparse 
set(gca,'XTick',1:max(optInds))
set(gca,'XTickLabel',lambdaRange)
set(gca,'XTickLabelRotation',45)
xlabel('lambda')
boxplot(lambdaFit,'PlotStyle','compact')
set(gca,'XTick',1:3:max(optInds), 'XTickLabel',round(lambdaRange(1:3:max(optInds)), 3))
xlabel('lambda')
ylabel(method)

subplot(3,1,1)
method = 'aic';
[test1, test2] = SelectLam(Preds, Responses, Prior, options, method);
lambdaFit = outFit{2};
optInds = outFit{3};
histogram(optInds); %in this example, a significant portion are trying to be more sparse 
set(gca,'XTick',1:max(optInds))
set(gca,'XTickLabel',lambdaRange)
set(gca,'XTickLabelRotation',45)
xlabel('lambda')
ylabel(method)

subplot(3,1,2)
method = 'bic';
outFit = SelectLam(Preds, Responses, Prior, options, method);
lambdaFit = outFit{2};
optInds = outFit{3};
histogram(optInds); %in this example, a significant portion are trying to be more sparse 
set(gca,'XTick',1:max(optInds))
set(gca,'XTickLabel',lambdaRange)
set(gca,'XTickLabelRotation',45)
xlabel('lambda')
ylabel(method)

subplot(3,1,3)
method = 'ebic';
outFit = SelectLam(Preds, Responses, Prior, options, method);
lambdaFit = outFit{2};
optInds = outFit{3};
histogram(optInds); %in this example, a significant portion are trying to be more sparse 
set(gca,'XTick',1:max(optInds))
set(gca,'XTickLabel',lambdaRange)
set(gca,'XTickLabelRotation',45)
xlabel('lambda')
ylabel(method)

subplot(3,1,1)
method = 'aic';
boxplot(lambdaFit,'PlotStyle','compact')
set(gca,'XTick',1:3:max(optInds), 'XTickLabel',round(lambdaRange(1:3:max(optInds)), 3))
xlabel('lambda')
ylabel(method)

subplot(3,1,2)
method = 'bic';
boxplot(lambdaFit,'PlotStyle','compact')
set(gca,'XTick',1:3:max(optInds), 'XTickLabel',round(lambdaRange(1:3:max(optInds)), 3))
xlabel('lambda')
ylabel(method)

subplot(3,1,3)
method = 'bic';
boxplot(lambdaFit,'PlotStyle','compact')
set(gca,'XTick',1:3:max(optInds), 'XTickLabel',round(lambdaRange(1:3:max(optInds)), 3))
xlabel('lambda')
ylabel(method)






