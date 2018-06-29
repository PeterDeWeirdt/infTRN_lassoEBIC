function [geneFit, netFit] = getFit(X, Ys, Prior, Opts, method, parallel)
%% Goal: Get our best fit model for each column of Y and calculate how good our fit
% using a variety of methods, including ebic, aic, and bic
% Note: 'MATLAB:mir_warning_maybe_uninitialized_temporary' is turned off
% for the time being. Warned about assigning temporary variables that might
% be lost 
%% Inputs: 
% X -- Samples by predictors matrix
% Ys -- Response matrix, with each column corresponding to a response
% Prior -- Our prior matrix
% Opts -- Options for GLMNET
% Method -- string for ebic, aic, or bic
% parralel -- boolean whether to use parralel for loop in glmnet
%% References:
% (1) Qian et al. (2013) "Glmnet for Matlab."
% http://www.stanford.edu/~hastie/glmnet_matlab/
% (2) Greenfield et al. (2013) "Robust data-driven incorporation of prior
% knowledge into the inference of dynamic regulatory netowrks"
% (3) Chen et al. (2008) "Extended bayesian information criteria for model
% selection with large model spaces"
%% Output:
% geneFit -  matrix with genes in rows and lambdas in columns
% netFit - median fit for each lambda
%% Author:
% Peter Deweirdt, Summer Intern, Divisions of Immunobiology and Biomedical
%   Informatics, Cincinnati Children's Hospital
%% Date: 6/15/2018
nResp = size(Ys, 1);
lambdaFit = zeros(nResp , size(Opts.lambda, 2));
nonzeroMat = lambdaFit;
optLams = zeros(1, nResp);
optNonzero = optLams;
lambdaRange = fliplr(Opts.lambda);
n_samps = size(Ys, 2);
if parallel
    parforArg = Inf;
else
    parforArg = 0;
end
parfor (i = 1:nResp, parforArg)
    % get (finite) predictor indices for each response // filter
    currPrior = Prior;
    currWeights = currPrior(i,:);
    % limit to predictors with finite lambda penalty (e.g., to exclude TF mRNA self-interaction loops)
    predInds = find(isfinite(currWeights)); 
    currPreds = zscore(X(predInds,:)');
    currResp = zscore(Ys(i,:)');
    currOpts = Opts;
    currOpts.penalty_factor = currPrior(i,predInds)';    
    if i == 1
        display(strcat(['Size of x is: ' num2str(size(currPreds)) ', and size of y is: ' num2str(size(currResp))]))
    end
    lsoln = glmnet(currPreds,currResp,'gaussian',currOpts);
    % Use output of GLMNET to check fit
    currBetas = fliplr(lsoln.beta); % flip so that the lambdas are increasing
    currA0 = fliplr(lsoln.a0');
    RSSs = sum(bsxfun(@minus,(bsxfun(@plus, currPreds*currBetas, currA0)),...
        currResp).^2, 1);
    n_nonzero = fliplr(lsoln.df');
    if strcmp(method, 'ebic')
        gamma = 1; % Parameter for how much the 'extend portion of EBIC is weighted
        tot_preds = size(currBetas, 1);
        [Fit, optInd] = ebic(RSSs, n_samps, n_nonzero, tot_preds, gamma);
    elseif strcmp(method,  'bic')
        [Fit, optInd] = mybic(RSSs, n_samps, n_nonzero);
    elseif strcmp(method,  'aic')
        [Fit, optInd] = myaic(RSSs, n_samps, n_nonzero);
    else
        display(strcat(method, ' not yet supported'))
    end
    lambdaFit(i,:) = Fit; 
    nonzeroMat(i,:) = n_nonzero;
    optLams(i) = lambdaRange(optInd);
    optNonzero(i) = n_nonzero(optInd);
end
geneFit = lambdaFit;
netFit = median(lambdaFit);

%% Plot ranges of solutions
if size(Opts.lambda, 2) > 10
    logx = 1:floor(size(Opts.lambda, 2)/10):size(Opts.lambda, 2);
else
    logx = 1:1:size(Opts.lambda, 2);
end

set(0, 'defaultAxesFontSize', 12);
subplot(4,1,1)
boxplot(geneFit,'PlotStyle','compact')
title('Gene Fit')
set(gca,'XTick',logx, 'XTickLabel',round(lambdaRange(logx), 3))
xlabel('\lambda Range')
ylabel(method)

subplot(4,1,2)
scatter(lambdaRange, netFit, 30, mean(nonzeroMat), 'filled')
cb = colorbar;
ylabel(cb, 'num. pred.')
title('Network Fit')
xlabel('lambda')
ylabel(method)
box on
set(gca,'YTick',min(geneFit(:)):((max(geneFit(:)) - min(geneFit(:)))/4):max(geneFit(:)))
set(gca,'xscale','log')
grid on
axis([lambdaRange(1) lambdaRange(end) min(geneFit(:)) max(geneFit(:))])
xlabel('\lambda Range')

subplot(4,1,3)
binedges = round(lambdaRange(logx), 3);
histogram(optLams, binedges)
set(gca, 'XScale', 'log')
xlabel('\lambda Range')
ylabel('Count')
title('Optimal \lambda per gene')
grid on
axis tight

subplot(4,1,4)
histogram(optNonzero)
xlabel('Num. Nonzero')
ylabel('Count')
title('Optimal edges per gene')
grid on 
axis tight

