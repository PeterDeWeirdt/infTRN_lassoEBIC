%% visPerformance_th17_multiParams_3way
% looking at CV = 5 results, 5th training set, limited to 50 samples
% enable subplots
% visualize results from performanceTh17_1704
%% DIFFERENT FROM ALL PREVIOUS PERFORMANCE.....M 
% -- _rawPR --> uses aupr_step (not aupr_fastP.m, which is not good when
% dealing with a small set of rank values)
% -- a new version of the performance code.  I didn't realize that I had limited recovery
% to genes that have edges in the gold standard and the model genes

%% Debugging
clear all
%close all
%% Path to required functions
% path2ems_functions = fullfile('~','erm','MATLAB','emily_functions');

%% plotting
% axVals = [0 1 0 1]; % for P-R plots
zoomRecalls = [1 .25];
fontSize = 14; % for subplots
regWidth = 2; % for p-r plots
controlWidth = 3;
fontSizeBar = 18; % for AUPR bar graphs
xSize = 5;       % for P-R plots
ySize = 4.75;
lineWidth = 2; % for random

figPause = 1;  % figPause --> 0, don't set figure for pdf, figPause --> 1, you can drag the figure windows to make the figure larger for pdf

darkBlue = [0 0 170]./255;
mediumBlue = [0 85 255]./255;
lightBlue = [0.5843    0.8157    0.9882];
darkRed = [170 0 0]./255;
mediumRed = [228 26 28]./255;
pink = [ 0.9686    0.5059    0.7490];
mediumGreen = [0 .45 0];
darkGreen = [0 .75 0];
lightGreen = [0 .25 0];

% lightBlue2 = [91, 207, 244] / 255; 
% purple = [103, 2, 94] / 255;

%% Gold standards
gsFolder = './RNAseq_inputs_overlap/priors/';
% gsPlotGroups = { ...
%     'KO75','KO75',{'th17_whole_K_cut_prcnt_20_num_tfs_28_sam_0_deseq_cut_0.25_Aug_8_2012_priorCut0p75'},.3};
%    ...'KC1p5','KC1p5',{'th17_whole_KC_cut_prcnt_20_num_tfs_28_sam_0_deseq_cut_0.25_Aug_8_2012_priorCut1p5'},.2;
%       'KC1','KC1',{'th17_whole_KC_cut_prcnt_20_num_tfs_28_sam_0_deseq_cut_0.25_Aug_8_2012_priorCut1p0'},.3;
%       'ChIP75','ChIP75',{'th17_whole_C_cut_prcnt_0_num_tfs_28_sam_0_deseq_cut_1_Aug_8_2012_priorCut0p75'},.4};

% gsPlotGroups = { ...
%     'KO75','KO75',{'th17_whole_K_cut_prcnt_20_num_tfs_28_sam_0_deseq_cut_0.25_Aug_8_2012_priorCut0p75'},.3;
%     'KC1','KC1',{'th17_whole_KC_cut_prcnt_20_num_tfs_28_sam_0_deseq_cut_0.25_Aug_8_2012_priorCut1p0'},.3;
%     'ChIP75','ChIP75',{'th17_whole_C_cut_prcnt_0_num_tfs_28_sam_0_deseq_cut_1_Aug_8_2012_priorCut0p75'},.4;
%     'KC1p5','KC1p5',{'th17_whole_KC_cut_prcnt_20_num_tfs_28_sam_0_deseq_cut_0.25_Aug_8_2012_priorCut1p5'},.2;
%     'KOrk','KOrk',{'regevTh17_RNAseq_TF_KO_sp.tsv'},.2;
%     'KOall','KOall',{'KO75_KOrk_sp.tsv'},.3;
% };

gsFile = './RNAseq_inputs_overlap/priors/KC1p5_sp.tsv';
gsPlotGroups = {'KC1p5_sp.tsv'};
gsRegsFile = './RNAseq_inputs_overlap/priors/goldStandardGeneLists/targGenesPR_mm9mm10.txt';
targGeneFile = './RNAseq_inputs_overlap/targRegLists/targetGenes_names.txt';
if targGeneFile
    geneIn = fopen(targGeneFile,'r');
    C = textscan(geneIn,'%s');
    fclose(geneIn);
    potTargGenes = C{1};
    totTargGenes = length(potTargGenes);
end

totGsGroups = size(gsPlotGroups,1);
groupInds = 1:totGsGroups;
gsInfs = [];    % file base, text, regs, reg-targs, random P-R
gsInfs.regs = {};      % keep track of the gs TFs across all sets, and target genes, because LASSO networks are big
fid = fopen(gsFile,'r');
% get first line and see what regulators we have
tline = fgetl(fid);    
totCols = length(cellstr(strsplit(tline,'\t')));
fclose(fid);
% get the rest of the data using textscan
fid = fopen(gsFile,'r');
if totCols == 3
    C = textscan(fid,['%s%s%f'],'Delimiter','\t','Headerlines',1);
elseif totCols > 3
    C = textscan(fid,['%s%s%f' repmat('%s',1,totCols-3)],'Delimiter','\t','Headerlines',1);
else
    error('Prior matrix should be in sparse format (TF,gene,edge weight).')
end
fclose(fid);
gsRegsTmp = C{1};
gsTargsTmp = C{2};
gsWeights = C{3};
keepWeights = find(gsWeights); % indices of edges with nonzero weight
% limit to TF-gene interactions considered by the model
if gsRegsFile
    fid = fopen(gsRegsFile,'r');
    C = textscan(fid,'%s','Headerlines',0);
    fclose(fid);
    Regs = C{1};
    gsRegInds = find(ismember(gsRegsTmp,Regs));
else
    gsRegInds = 1:length(gsRegsTmp); % consider all TFs
end
if targGeneFile
    gsTargInds = find(ismember(gsTargsTmp,potTargGenes));
else
    gsTargInds = 1:length(gsRegsTmp); % consider all target genes
end
keepInds = intersect(intersect(gsRegInds,gsTargInds),keepWeights);
gsRegs = cellstr(strvcat(gsRegsTmp{keepInds}));
gsTargs = cellstr(strvcat(gsTargsTmp{keepInds}));
totGsInts = length(keepInds);
uniGsRegs = unique(gsRegs);
totGsRegs = length(uniGsRegs);
gsInfs.regs = uniGsRegs;
if not(length(targGeneFile))  % all targets in the prior are used if no target gene list was supplied
    potTargGenes = unique(gsTargs);
    totTargGenes = length(potTargGenes);
end
gsInfs.totPotInts = totTargGenes*totGsRegs;
randPR = totGsInts/gsInfs.totPotInts;
[dd, currFileBase,ext] = fileparts(gsFile);
gsInfs.fileBase = gsFile;
gsInfs.figText = strrep(strrep(currFileBase,'_sp',''),'_',' '); % use file name 
        
ints = [gsRegs gsTargs];
totInts = length(ints);
gsInfs.edges = cell(totInts,1);
for ii = 1:totInts
    gsInfs.edges{ii} = strjoin({ints{ii,:}},',');
end    
gsInfs.randPR = randPR;
gsInfs.targs = unique(gsTargs);
    
gsInfs.edgesByTf = cell(totGsRegs,1);
gsInfs.randAuprByTf = zeros(totGsRegs,1);
for gind = 1:totGsRegs
    currInds = find(ismember(gsRegs,uniGsRegs(gind)));
    gsInfs.edgesByTf{gind} = {gsInfs.edges{currInds}}';
    gsInfs.randAuprByTf(gind) = length(gsInfs.edgesByTf{gind})/totTargGenes;
end        
% we will fill in these values and then save a version of gsInfs for
% each network
gsInfs.precisions = {};
gsInfs.recalls = {};
gsInfs.fprs = {};
gsInfs.arocs = NaN;
gsInfs.auprs = NaN;
gsInfs.auprsByTf = zeros(totGsRegs,1);
gsInfs.arocsByTf = zeros(totGsRegs,1);
gsInfs.precisionsByTf = cell(totGsRegs,1);
gsInfs.recallsByTf = cell(totGsRegs,1);
gsInfs.fprsByTf = cell(totGsRegs,1);
    
fprintf([num2str(totGsInts) '\t' gsInfs.figText '\n']);


% gsInputs = {'KO75','mm10_merged/th17/KO75_sp.tsv';
% 	'ChIP75','mm10_merged/th17/ChIP75_sp.tsv';
% 	'KC1','mm10_merged/th17/KC1_sp.tsv';
%     'KC1p5','mm10_merged/th17/KC1p5_sp.tsv';
%     'KOrk','mm10_merged/th17/regevTh17_RNAseq_TF_KO_sp.tsv';
%     'KOall','mm10_merged/th17/KO75_KOrk_sp.tsv';
%     };  
  
%% Inf P-R results
% outInf = '_ext3_SS100';
% outInf = '_SSedge_SS50';
% outInf = '_ext3_SS100';
% outInf = '_SS50';
% % outInf = '_ext_SS100';
% % outInf = '_eNet0p5';
% % prInfLS = '1802_sN';%G';
% % prInfLS = '1802_s_ext3_SS100_legacy';
% prInfLS = '1802_SSedge_SS50_legacy';
% prInfLS = '1802_s_ext3_SS100_Gene_noMerge';
% prInfLS = '1802_sN_ext_SS100_legacy';
% prInfLS = '1802_s_ext_SS100_';
% % prInfLS = '1802_s_eNet0p5_legacy3';
% % prInfLS = '1802_s_ext_SS100_legacyStabRange';
% prInfLS = '1802_s_ext3_SS100_legacyStabL80';
% prInfLS = '1802_s_ext3_SS100_setDiff_legacyStabL80';
% prInfLS = '1802_s_ext3_SS100_rmLB_legacy80';    % from working code and works
% prInfLS = '1802_s_ext3_SS100_debugSetDiff_L80'; % rmLB in nonwork code -- doesn't work
% prInfLS = '1802_s_ext3_SS100_totSS_L80';        % change totSS
% prInfLS = '_fullMods/stabNets_sN5_20gene_finerRes';        % change totSS
% 
% prInfBB = '1801';


addOutBit = '_PR';

% params: (1) weight (g-prior or bias), (2) TFA option, (3) BBSR-BIC or
% LASSO-StARS, (4) line style, (5) line color, (6) line width
params = { % NaN --> No Prior
    1,'TFA','mLASSO-cv-net','-', darkRed,regWidth;
    1,'TFA','mLASSO-cv-net-a-0.99','-', pink,regWidth;
    1,'TFA','mLASSO-cv-net-refit','--', darkRed,regWidth;
    1,'TFA','mLASSO-cv-gene','--', mediumRed,regWidth
     1, 'TFA', 'mLASSO-ebic','-',mediumBlue, regWidth;
     1, 'TFA', 'mLASSO-StARS', '-', darkGreen, regWidth};
    
%    1, 'TFA', 'LASSO-StARS-bulk','-', darkRed, regWidth};



%     1.1,'','BB','-',darkBlue,regWidth;
%     1.5,'','BB',':',darkBlue,regWidth;
%     1.1,'_noTFA','BB','-',mediumBlue,regWidth;
%     1.5,'_noTFA','BB',':',mediumBlue,regWidth
%     NaN,'','LS','-',pink,controlWidth; % NaN --> No Prior
%     1,'','LS','--',darkRed,regWidth;
%     .5,'','LS','-',darkRed,regWidth;
%     .25,'','LS',':',darkRed,regWidth;
%     .5,'_noTFA','LS','-',mediumRed,regWidth;
%     .25,'_noTFA','LS',':',mediumRed,regWidth};

% priorList = '/Users/ngurb2/Desktop/prior_miraldi_Th17_48h_cut4_bp10000_sATAC_p1Em5_huA/priorList_pr.txt';
%  = [1 5 6 7 8];% 11];%[1 2 5 8 11];%26];%:13 18:19];  % correspond to rows in the priorList
% priorInds = [1:2];% 8];%[8];% 2];%[6:8]; % 8 --> ChIP, 1 --> sA(Th17) 1, 2, 6-8
% priorInds = [2:7 9 10 12];

%% Performance Analysis  
% add path to performance and other functions
% addpath(path2ems_functions) % add path to functions
%load linecolors.mat         % matrix of plotting colors for figures

controlColor = .66*[1 1 1];

% outInf
% prInfLS

%% read in prior file

priorInds = 1;
for priorInd = priorInds

% infModelDirACTION = '/Users/ngurb2/Desktop/action';
% % noPriorResultsLS = 'Th17_48h_cut4_sA_p5_huA_bias100_noTFA.mat';%'ChBod10kMac10bin_bias100_noTFA.mat'; 
% infModelDirSCENIC = '/Users/ngurb2/Desktop/pyscenic/Th17_scRNA_bryson';
% noPriorResultsBB = 'Th17_48h_cut4_sA_p5_huA_bs50_w1_noTFA.mat';
outDir = fullfile('./matlab_pr/PRFiguresAll_0627');
% gompers
mkdir(outDir)

totInfResults = size(params,1);
legendInf = cell(totInfResults+1,1);
legendInf{1} = 'Random';
totGsGroups = size(gsPlotGroups,1);
% need to get some info on priors to start and so will load first results
% load(fullfile(infModelDirLS,noPriorResultsLS))

% Plot random performance

for groupCount = groupInds
    auprMax = 0.5;
    currGsTitle = 'Th17 KC';
%     currGsOut = [gsPlotGroups{groupCount,2} addOutBit];
%     currGroups = gsPlotGroups{groupCount,3};
%     auprMax = gsPlotGroups{groupCount,4};
%     [xx, gsInd] = intersect({gsInfs.figText},currGsTitle);
%     currGS = gsInfs(gsInd).figText;
    figure(1), clf % P-R
    % plot random P-R ( : grey )
    plot([0 1],gsInfs.randPR*[1 1],':','LineWidth',lineWidth+1,'Color',controlColor)
    hold on, grid on
    set(gca,'FontSize',fontSize,'FontWeight','Bold')
    xlabel('Recall','FontWeight','Bold','FontSize',fontSize)
    ylabel('Precision','FontWeight','Bold','FontSize',fontSize)
    axis([ 0 1 0 1])
    grid on, grid minor
    
    figure(2), clf % AUPR
    plot(gsInfs.randPR*[1 1],[.4 totInfResults+1.6],'-','LineWidth',lineWidth+1,'Color',controlColor)
    hold on, grid on, grid minor
    set(gca,'FontSize',fontSize,'FontWeight','Bold')
    xlabel('AUPR','FontWeight','Bold','FontSize',fontSizeBar)

        
for tind = 1:totInfResults
    weight = params{tind,1};
    tfaOpt = params{tind,2};
    methodInf = params{tind,3};
%     priorInf = params{tind,4};
    currMark = params{tind,4};
    currColor = params{tind,5};
    currWidth = params{tind,6};
    
    if isnan(weight) % No Prior condition
        legendInf{tind+1} = ['No Prior, ' methodInf];
%         prMat = fullfile(infModelDir,noPriorResults);
    else % Model with Prior
        legendInf{tind+1} = [methodInf];
        if (find(ismember({'mLASSO-cv-net'},methodInf)))
            prMat = fullfile('/Users/dewpz7/Documents/infTRN_lassoStARS/Th17_example_EBIC/CV_outputs_TH17_full/networks_10Fold_cv/50bootstraps_0p01cutoff/PR_KC1p5/ATAC_Th17_bias50.mat');
        elseif find(ismember({'mLASSO-ebic'},methodInf))
            prMat = fullfile('/Users/dewpz7/Documents/infTRN_lassoStARS/Th17_example_EBIC/ebic_outputs_TH17_full/networks_ebic/50bootstraps_0p01cutoff/PR_KC1p5/ATAC_Th17_bias50.mat');
        elseif find(ismember({'mLASSO-StARS'},methodInf))
            prMat = fullfile('/Users/dewpz7/Documents/infTRN_lassoStARS-master/Th17_example/outputs/networks_targ0p05_SS50_bS5/Network0p05_15tfsPerGene/PR_KC1p5/ATAC_Th17_bias50.mat');
        elseif find(ismember({'mLASSO-cv-gene'},methodInf))
            prMat = fullfile('/Users/dewpz7/Documents/infTRN_lassoStARS/Th17_example_EBIC/CV_outputs_TH17_full/networks_5Fold_cv/50bootstraps_0p01cutoff_rank_confidenceselection_gene/PR_KC1p5/ATAC_Th17_bias50.mat');
        elseif find(ismember({'mLASSO-cv-net-a-0.99'},methodInf))
            prMat = fullfile('/Users/dewpz7/Documents/infTRN_lassoStARS/Th17_example_EBIC/CV_outputs_TH17_full/networks_5Fold_cv/50bootstraps_0p01cutoff__rank_confidence_selection_network/PR_KC1p5/ATAC_Th17_bias50alpha_0p99.mat');
        elseif find(ismember({'mLASSO-cv-net-refit'},methodInf))
            prMat = '/Users/dewpz7/Documents/infTRN_lassoStARS/Th17_example_EBIC/CV_outputs_TH17_full/networks__10Foldcv/50bootstraps_0p01cutoff__rank_confidence_selection_network/PR_KC1p5/ATAC_Th17_bias50_alpha1.mat';
        else
            disp(methodInf)
        end
    end 
    load(prMat)

    figure(1) % plot P-R 
    plot(gsInfs.recalls,gsInfs.precisions,...
        currMark,'LineWidth',currWidth,'Color',currColor)
    
    figure(2), % AUPR
    % add no Prior AUPR
    barh(totInfResults-tind+1,gsInfs.auprs,...
        'FaceColor',currColor)

end

    figInf = fullfile(outDir);
    figure(1), % P-R
    
    for zoomRecall = zoomRecalls
        figure(1)
        title('Precision-Recall Th17')
        axis([0 zoomRecall 0 1]) 
        set(gca,'XTick',[0:.1:zoomRecall])
        set(gca,'YTick',[0:.1:1])
        grid on, grid minor
        figure(1)
        currFile = [figInf '_prZoom' num2str(zoomRecall*100)];
        fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [xSize ySize]);
        print('-painters','-dpdf','-r150',[currFile '.pdf'])
        saveas(gcf,currFile,'fig')
        disp(currFile)
        
    end
%     gompers
    figure(1)
    legend(legendInf,'Location','East')
    currFile = [figInf '_legend'];
    fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [xSize xSize]);
    print('-painters','-dpdf','-r150',[currFile '.pdf'])
    saveas(gcf,currFile,'fig')
    disp(currFile)
        
%     gompers
    
    figure(2), % AUPR
    title('Precision-Recall Th17')
%     if figPause
%         pause
%     end
%     currFile = [figInf '_aupr'];
%     fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [xSize ySize]);
%     print('-painters','-dpdf','-r150',[currFile '.pdf'])
%     saveas(gcf,currFile,'fig')
%     disp(currFile)
    
        axis([0 auprMax .4 totInfResults+.6]) 
        set(gca,'XTick',[0:.05:auprMax],'YTick',[])
        grid on, grid minor

    currFile = [figInf '_aupr' num2str(auprMax*100)];
    fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [xSize ySize]);
    print('-painters','-dpdf','-r150',[currFile '.pdf'])
    disp(currFile)


end
end
% end
disp('Finished')