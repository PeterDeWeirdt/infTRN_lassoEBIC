%% example_workflow_Th17_EBIC
% Use mLASSO-EBIC to build a TRN from gene expression and prior
% information in four steps. Please refer to each function's help
% annotations for descriptions of inputs, outputs and other information.
%% Highlited Changes:
% Steps 1,2, and 5 (-> 4 here) are consistent with "example_workflow_TH17", 
% but step 3 now uses EBIC in order to select model parameters for bootstrapped 
% samples, edges are ranked by confidence scores and rank combined. 
%% References: 
% (1) Miraldi et al. (2018) "Leveraging chromatin accessibility for 
% transcriptional regulatory network inference in T Helper 17 Cells"
% (2) Qian et al. (2013) "Glmnet for Matlab."
% http://www.stanford.edu/~hastie/glmnet_matlab/
% (3) Liu, Roeder, Wasserman (2010) "Stability Approach to Regularization 
%   Selection (StARS) for High Dimensional Graphical Models". Adv. Neural.
%   Inf. Proc.
% (4) Muller, Kurtz, Bonneau. "Generalized Stability Approach for Regularized
%   Graphical Models". 23 May 2016. arXiv.
%% Authors: Emily R. Miraldi, Ph.D., Divisions of Immunobiology and Biomedical
%   Informatics, Cincinnati Children's Hospital
% Peter DeWeirdt, Summer Intern, Divisions of Immunobiology and Biomedical
%   Informatics, Cincinnati Children's Hospital
%% Date: June 14, 2018 -PD

clear all
close all
restoredefaultpath

matlabDir = '..';

addpath(fullfile(matlabDir,'ebicFxns'))
addpath(fullfile(matlabDir,'infLassoStARS'))
addpath(fullfile(matlabDir,'glmnet'))
addpath(fullfile(matlabDir,'customMatlabFxns'))

%% 1. Import gene expression data, list of regulators, list of target genes
% into a Matlab .mat object
geneExprTFAdir = './outputs/processedGeneExpTFA';
mkdir(geneExprTFAdir)
normGeneExprFile = './inputs/geneExpression/th17_RNAseq254_DESeq2_VSDcounts.txt';
targGeneFile = './inputs/targRegLists/short_targetGenes_names.txt';
potRegFile = './inputs/targRegLists/potRegs_names.txt';
tfaGeneFile = './inputs/targRegLists/genesForTFA.txt';
geneExprMat = fullfile(geneExprTFAdir,'geneExprGeneLists.mat');

disp('1. importGeneExpGeneLists.m')
importGeneExpGeneLists(normGeneExprFile,targGeneFile,potRegFile,...
    tfaGeneFile,geneExprMat)

%% 2. Given a prior of TF-gene interactions, estimate transcription factor 
% activities (TFAs) using prior-based TFA and TF mRNA levels
priorName = 'ATAC_Th17';
priorFile = ['./inputs/priors/' priorName '.tsv']; % Th17 ATAC-seq prior
edgeSS = 50;
minTargets = 3;
[xx, priorName, ext] = fileparts(priorFile);
tfaMat = fullfile(geneExprTFAdir,[priorName '_ss' num2str(edgeSS) '.mat']);

disp('2. integratePrior_estTFA.m')
integratePrior_estTFA(geneExprMat,priorFile,edgeSS,...
     minTargets, tfaMat)

%% 3. use EBIC to select model parameters for bootstrapped 
% samples, edges are ranked by confidence scores and rank combined.  
lambdaBias = .5;
tfaOpt = ''; % options are '_TFmRNA' or ''
nBoots = 2;
minFrac = 0.5; % Minimum fraction of bootstraps that an edge must be in to be included in the final network
lambdaMin = .01;
lambdaMax = 1;
totLogLambdaSteps = 25; % will have this many steps per log10 within bStARS lambda range
leaveOutSampleList = '';
leaveOutInf = ''; % leave out information 
TRNweightedDir = fullfile('./outputs',strrep(['TRNweighted_boots' ...
    num2str(nBoots) '_frac' num2str(minFrac) leaveOutInf],'.','p'));
mkdir(TRNweightedDir)
netSummary = [priorName '_bias' strrep(num2str(100*lambdaBias),'.','p') tfaOpt];
trnOutMat = fullfile(TRNweightedDir,netSummary);

disp('3. TRNebic.m')
TRNebic(geneExprMat,tfaMat,lambdaBias,tfaOpt,...
    nboots,minFrac,lambdaMin,lambdaMax,totLogLambdaSteps,...
    trnOutMat,leaveOutSampleList)

%% 4. Calculate precision-recall relative to KO-ChIP G.S.
gsFile = './inputs/priors/KC1p5_sp.tsv';
prNickName = 'KC1p5';
rankColTrn = 3;
prTargGeneFile = './inputs/priors/goldStandardGeneLists/targGenesPR_mm9mm10.txt';
gsRegsFile = '';
prDir = fullfile(networkSubDir,['PR_' prNickName]);
mkdir(prDir)
prMatBase = fullfile(prDir,netSummary);
prFigBase = fullfile(prDir,netSummary);

display('5. calcPRinfTRNs')
calcPRinfTRNs(outNetFileSparse,gsFile,rankColTrn,...
    prTargGeneFile,gsRegsFile,prMatBase,prFigBase)