function estimateInstabilitiesTRN(geneExprMat,tfaMat,lambdaBias,tfaOpt,...
    totSS,lambdaMin,lambdaMax,totLogLambdaSteps,subsampleFrac,instabOutMat,...
    leaveOutSampleList)
%% estimateInstabilitiesTRN(geneExprMat,tfaMat,lambdaBias,tfaOpt,...
%     totSS,lambdaMin,lambdaMax,totLogLambdaSteps,subsampleFrac,instabOutMat,...
%     leaveOutSampleList)
%% Goal: Estimate mLASSO-StARS instabilities for given input prior, prior 
% reinforcement, TFA methods, and FIXED range of lambda penalties (as
% opposed to closely related function, estimateInstabilitiesTRNbStARS.m,
% which calculates upper and lower bounds on lambda range containing target
% instability, with the goal of speeding computation)
%% References:
% Miraldi et al. "Leveraging chromatin accessibility for 
%   transcriptional regulatory network inference in T Helper 17 Cells"
%  for Matlab (2013) Qian, J., Hastie, T., Friedman, J., Tibshirani, 
%   R. and Simon, N. -- http://www.stanford.edu/~hastie/glmnet_matlab/
% Liu, Roeder, Wasserman (2010) "Stability Approach to Regularization 
%   Selection (StARS) for High Dimensional Graphical Models". Adv. Neural.
%   Inf. Proc.
%% Author: Emily R. Miraldi, Ph.D., Divisions of Immunobiology and Biomedical
%   Informatics, Cincinnati Children's Hospital
%% INPUTS:
% geneExprMat -- a .mat file containing gene expression data, gene lists
%   (target genes, potential regulators,...), e.g., as generated by
%   importGeneExpGeneLists.m
% tfaMat -- a .mat file containing the prior of TF-gene interactions as
%   well as TFA (prior-based and TF mRNA), e.g., as generated by 
%   integratePrior_estTFA.m
% lambdaBias -- a fractional lambda penalty that a TF-gene interaction in 
%   the prior matrix will have (e.g., penalty term for prior-supported edge
%   is reduced to bias * lambda, where lambda corresponds to the "reference" 
%   penalty applied to edges not in the prior), belongs to range [0,1]
% tfaOpt -- two options recognized:
%   '' --> TFA based on target gene and priors is used, 
%   '_TFmRNA' --> TF based on TF mRNA levels are used
% totSS -- total number of subsamples for the final instability estimates
% lambdaMin -- lambda penalty lower bound [0, infinity), note glmnet
%   is very slow for lambda < .01
% lambdaMax -- lambda penalty upper bound [0, infinity)
% totLogLambdaSteps -- number of steps per log10 lambda range
% subsampleFrac -- fraction of samples to use for subsampling. Liu et al.,
%   recommend subsample size = floor(100/sqrt(N)), where N = total samples.
%   Given that some TRN inference datasets have < 100 samples, Miraldi et
%   al, used .63*N.
% instabOutMat -- full file name and path for output .mat file, will also
%   be used as file base name for figures 
% leaveOutSampleList -- a text file, where each line corresponds to a
%   sample condition to be left-out of the inference procedure (e.g., for
%   the purposes of cross-validation)
%% OUTPUTS:
% instabOutMat -- contains network- and gene-level instabilities,
%   lambdaRange, number of nonzero subsamples per edge (used to rank
%   TF-gene interactions by subsequent scripts)
% instabOutMat.fig + .pdf -- showing network- and gene-level instabilities
%   as a function of final lambda range
%% NOTE: Here edges in the prior are treated in a binary manner (present, 
% nonzero, or absent, zero). This code could be modified to include
% real-valued prior edge confidences.

%% load gene expression and TFA
load(geneExprMat)
totSamps = size(targGeneMat,2);
responseMat = targGeneMat;
load(tfaMat)
% have to match prior names with target gene expression and TFA
if tfaOpt
    disp('noTfa option')
    pRegs = pRegsNoTfa;
    pTargs = pTargsNoTfa;
    priorMatrix = priorMatrixNoTfa;
end
[uniNoPriorRegs,uniNoPriorRegInds] = setdiff(potRegs_mRNA,pRegs);
allPredictors = cellstr(strvcat(strvcat(pRegs),strvcat(uniNoPriorRegs)));
totPreds = length(allPredictors);

[vals, targGeneInds, priorGeneInds] = intersect(targGenes,pTargs);
totTargGenes = length(targGenes);
totPRegs = length(strvcat(pRegs));
priorMat = zeros(totTargGenes,totPreds);
priorMat(targGeneInds,1:totPRegs) = priorMatrix(priorGeneInds,:);

%% set input priors and predictors    
predictorMat = [medTfas; potRegMat_mRNA(uniNoPriorRegInds,:)];
if tfaOpt % use the mRNA levels of TFs
    currPredMat = zeros(totPreds,totSamps);
    for prend = 1:totPreds
        prendInd = find(ismember(potRegs_mRNA,allPredictors{prend}));
        currPredMat(prend,:) = potRegMat_mRNA(prendInd,:);
    end        
    predictorMat = currPredMat;
    disp(['TF mRNA used.'])
end
priorWeightsMat = ones(totTargGenes,totPreds) - (1-lambdaBias)*abs(sign(priorMat));

if tfaOpt
    %% set lambda penalty to infinity for positive feedback edges where TF 
    % mRNA levels serves both as gene expression and TFA estimate
    for pr = 1:totPreds
        targInd = find(ismember(targGenes,allPredictors{pr})); 
        if length(targInd) % set lambda penalty to infinity, avoid predicting a TF's mRNA based on its own mRNA level
            priorWeightsMat(targInd,pr) = inf; % i.e., target gene is its own predictor
        end
    end    
else % have to set prior inds to zero for TFs in TFA that don't have prior info
    for pr = 1:totPreds        
        if sum(abs(priorMat(:,pr))) == 0 % we have no target edges to estimate TF's TFA
            targInd = find(ismember(targGenes,allPredictors{pr}));
            if length(targInd) % And TF is in the predictor set
                priorWeightsMat(targInd,pr) = inf;
            end
        end
    end
end

%% Check whether to use full gene expression matrix or exclude leave-out set 
if leaveOutSampleList
    disp(['Leave-out set detected: ' leaveOutSampleList])
    % get leave-out set of samples
    fin = fopen(leaveOutSampleList,'r');
    C = textscan(fin,'%s','HeaderLines',0);
    fclose(fin);
    testSamples = C{1};
    testInds = find(ismember(conditionsc,testSamples));
    trainInds = setdiff(1:totSamps,testInds);
else        
    disp(['Full gene expression matrix used.'])
    trainInds = 1:totSamps; % all training samples used
    testInds = [];
end
subsampleSize = floor(subsampleFrac*length(trainInds));

%% get instability per edge 
% get total steps
lamLog10step = 1/totLogLambdaSteps;
logLamRange =  log10(lambdaMin):lamLog10step:log10(lambdaMax);
lambdaRange = [10.^logLamRange]';

figure(3), clf
% get Network-wise edge instabilities, and gene-wise edge instabilities    
[geneInstabilities,netInstabilities,ssMatrix] = ...
    getMLassoStARSinstabilitiesPerGeneAndNet(predictorMat(:,trainInds),...
    responseMat(:,trainInds),priorWeightsMat,lambdaRange,...
    subsampleSize, totSS);     

figInf = strrep(instabOutMat,'.mat','');
saveas(gcf,figInf,'fig')
fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [6 5]);
print('-painters','-dpdf','-r150',[figInf '.pdf'])
disp(figInf)        
save([instabOutMat '.mat'],'-v7.3',...
    'geneInstabilities','netInstabilities','ssMatrix','predictorMat',...
    'responseMat','priorMat','lambdaBias','lambdaRange','trainInds',...
    'subsampleSize','totSS',...
    'targGenes','priorWeightsMat','allPredictors','conditionsc') 