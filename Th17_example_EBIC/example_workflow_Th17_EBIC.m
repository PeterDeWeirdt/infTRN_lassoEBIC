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
targGeneFile = './inputs/targRegLists/targetGenes_names.txt';
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

%% 3. Select model parameters using a specified method

lambdaBias = .5;
tfaOpt = ''; % options are '_TFmRNA' or ''
lambdaMin = 0.01;
lambdaMax = 1;
totLogLambdaSteps = 5; % will have this many steps per log10 within bStARS lambda range
leaveOutSampleList = '';
leaveOutInf = ''; % leave out information 
method = 'ebic'; % options rn are ebic, aic, bic... to do 'cv', 'StARS'
fitDir = fullfile('./outputs',strrep(['fits_' method leaveOutInf],'.','p'));
mkdir(fitDir)
netSummary = [priorName '_bias' strrep(num2str(100*lambdaBias),'.','p') tfaOpt];
fitOutMat = fullfile(fitDir,netSummary);
parallel = true; 
if parallel
    if isempty(gcp('nocreate'))
        mypool = parpool();
    end
end

disp('3. estimateFitTRN.m')
EstimateFitTRN(geneExprMat,tfaMat,lambdaBias,tfaOpt,...
    method,lambdaMin,lambdaMax,totLogLambdaSteps,...
    fitOutMat,leaveOutSampleList, parallel)

%% 4. For the minimum fit score, rank TF-gene
% interactions, calculate confidences and network file for jp_gene_viz
% visualizations
priorMergedTfsFile = ['./inputs/priors/' priorName '_mergedTfs.txt'];
try % not all priors have merged TFs and merged TF files
    ls(priorMergedTfsFile) 
catch
    priorMergedTfsFile = '';
end
nboot = 50;
bootCut = .5;
rankMethod = 'confidence'; % rank or confidence
networkDir = strrep(fitDir,'fits','networks');
mkdir(networkDir);
networkSubDir = fullfile(networkDir,[num2str(nboot) 'bootstraps_' ...
    strrep(num2str(bootCut), '.', 'p') 'cutoff']);
mkdir(networkSubDir)
trnOutMat = fullfile(networkSubDir,netSummary);
outNetFileSparse = fullfile(networkSubDir,[netSummary '_sp.tsv']);
networkHistDir = fullfile(networkSubDir,'Histograms');
mkdir(networkHistDir)
bootsHistPdf = fullfile(networkHistDir,[netSummary '_bsHist']);

disp('4. buildTRNs_mLassoStARS.m')
buildTRNs_mLassoFit(fitOutMat,tfaMat,priorMergedTfsFile,...
    bootCut, nboot, parallel, rankMethod , bootsHistPdf, trnOutMat,outNetFileSparse)
%% 5. Calculate precision-recall relative to KO-ChIP G.S.
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
