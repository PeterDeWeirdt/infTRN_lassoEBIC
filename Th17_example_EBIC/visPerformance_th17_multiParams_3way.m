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
% close all

%% Path to required functions
path2ems_functions = fullfile('~','erm','MATLAB','emily_functions');

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

figPause = 0;  % figPause --> 0, don't set figure for pdf, figPause --> 1, you can drag the figure windows to make the figure larger for pdf

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
gsFolder = '/Users/emiraldi/erm/MariaP/Inferelator/input/GeneralPriors';
% gsPlotGroups = { ...
%     'KO75','KO75',{'th17_whole_K_cut_prcnt_20_num_tfs_28_sam_0_deseq_cut_0.25_Aug_8_2012_priorCut0p75'},.3};
%    ...'KC1p5','KC1p5',{'th17_whole_KC_cut_prcnt_20_num_tfs_28_sam_0_deseq_cut_0.25_Aug_8_2012_priorCut1p5'},.2;
%       'KC1','KC1',{'th17_whole_KC_cut_prcnt_20_num_tfs_28_sam_0_deseq_cut_0.25_Aug_8_2012_priorCut1p0'},.3;
%       'ChIP75','ChIP75',{'th17_whole_C_cut_prcnt_0_num_tfs_28_sam_0_deseq_cut_1_Aug_8_2012_priorCut0p75'},.4};

gsPlotGroups = { ...
    'KO75','KO75',{'th17_whole_K_cut_prcnt_20_num_tfs_28_sam_0_deseq_cut_0.25_Aug_8_2012_priorCut0p75'},.3;
    'KC1','KC1',{'th17_whole_KC_cut_prcnt_20_num_tfs_28_sam_0_deseq_cut_0.25_Aug_8_2012_priorCut1p0'},.3;
    'ChIP75','ChIP75',{'th17_whole_C_cut_prcnt_0_num_tfs_28_sam_0_deseq_cut_1_Aug_8_2012_priorCut0p75'},.4;
    'KC1p5','KC1p5',{'th17_whole_KC_cut_prcnt_20_num_tfs_28_sam_0_deseq_cut_0.25_Aug_8_2012_priorCut1p5'},.2;
    'KOrk','KOrk',{'regevTh17_RNAseq_TF_KO_sp.tsv'},.2;
    'KOall','KOall',{'KO75_KOrk_sp.tsv'},.3;
};

totGsGroups = size(gsPlotGroups,1);
groupInds = 1:totGsGroups;

% gsInputs = {'KO75','mm10_merged/th17/KO75_sp.tsv';
% 	'ChIP75','mm10_merged/th17/ChIP75_sp.tsv';
% 	'KC1','mm10_merged/th17/KC1_sp.tsv';
%     'KC1p5','mm10_merged/th17/KC1p5_sp.tsv';
%     'KOrk','mm10_merged/th17/regevTh17_RNAseq_TF_KO_sp.tsv';
%     'KOall','mm10_merged/th17/KO75_KOrk_sp.tsv';
%     };  
  
%% Inf P-R results
% outInf = '_ext3_SS100';
outInf = '_SSedge_SS50';
outInf = '_ext3_SS100';
outInf = '_SS50';
% outInf = '_ext_SS100';
% outInf = '_eNet0p5';
% prInfLS = '1802_sN';%G';
% prInfLS = '1802_s_ext3_SS100_legacy';
prInfLS = '1802_SSedge_SS50_legacy';
prInfLS = '1802_s_ext3_SS100_Gene_noMerge';
prInfLS = '1802_sN_ext_SS100_legacy';
prInfLS = '1802_s_ext_SS100_';
% prInfLS = '1802_s_eNet0p5_legacy3';
% prInfLS = '1802_s_ext_SS100_legacyStabRange';
prInfLS = '1802_s_ext3_SS100_legacyStabL80';
prInfLS = '1802_s_ext3_SS100_setDiff_legacyStabL80';
prInfLS = '1802_s_ext3_SS100_rmLB_legacy80';    % from working code and works
prInfLS = '1802_s_ext3_SS100_debugSetDiff_L80'; % rmLB in nonwork code -- doesn't work
prInfLS = '1802_s_ext3_SS100_totSS_L80';        % change totSS
prInfLS = '_fullMods/stabNets_sN5_20gene_finerRes';        % change totSS

prInfBB = '1801';


instCut = .1;  % cutoff for instability
outBit = [''];
timeOpt = ''; % added to get older versions of Th17 TRNs
addOutBit = '_PR';
kfoldCV = 5;
bsTot = 60;

% params: (1) weight (g-prior or bias), (2) TFA option, (3) BBSR-BIC or
% LASSO-StARS, (4) line style, (5) line color, (6) line width
params = {NaN,'','BB','-',lightBlue,controlWidth; % NaN --> No Prior
    1,'','BB','--',darkBlue,regWidth;
    1.1,'','BB','-',darkBlue,regWidth;
    1.5,'','BB',':',darkBlue,regWidth;
    1.1,'_noTFA','BB','-',mediumBlue,regWidth;
    1.5,'_noTFA','BB',':',mediumBlue,regWidth
    NaN,'','LS','-',pink,controlWidth; % NaN --> No Prior
    1,'','LS','--',darkRed,regWidth;
    .5,'','LS','-',darkRed,regWidth;
    .25,'','LS',':',darkRed,regWidth;
    .5,'_noTFA','LS','-',mediumRed,regWidth;
    .25,'_noTFA','LS',':',mediumRed,regWidth};

dataID = 'av9TRD';
priorList = '/Users/emiraldi/erm/MariaP/Inferelator/input/priorLists/priorList_1712validation.txt';
priorList = ['/Users/emiraldi/erm/MariaP/Inferelator/input/priorLists/priorList_1709validationTEXT_CV.txt'];
priorList = '/Users/emiraldi/erm/MariaP/Inferelator/input/priorLists/priorList_1712validationTEXT.txt';
%  = [1 5 6 7 8];% 11];%[1 2 5 8 11];%26];%:13 18:19];  % correspond to rows in the priorList
priorInds = [1:5];% 8];%[8];% 2];%[6:8]; % 8 --> ChIP, 1 --> sA(Th17) 1, 2, 6-8
% priorInds = [2:7 9 10 12];

%% Performance Analysis  
% add path to performance and other functions
addpath(path2ems_functions) % add path to functions
load linecolors.mat         % matrix of plotting colors for figures
load cbcolors
controlColor = .66*[1 1 1];

outInf
prInfLS

%% read in prior file
pIn = fopen(priorList,'r');
C = textscan(pIn,'%s%s','Delimiter','\t');
fclose(pIn);

priorFileTexts = C{1};
priorNames = C{2};
totPoi = length(priorFileTexts);

for priorInd = priorInds

infModelDirLS = ['/Users/emiraldi/erm/MariaP/Inferelator/input/lymph36broad_cv/Th17_TARGpap17v0FC0p58FDR10_REGSexp80nPapFC0p58FDR25_ivTh/Results_lassoStARS' outInf '/precisionRecall' prInfLS];
noPriorResultsLS = 'Th17_48h_cut4_sA_p5_huA_bias100_noTFA.mat';%'ChBod10kMac10bin_bias100_noTFA.mat'; 
infModelDirBB = ['/Users/emiraldi/erm/MariaP/Inferelator/output/Th17_TARGpap17v0FC0p58FDR10_REGSexp80nPapFC0p58FDR25_ivTh/precisionRecall' prInfBB];
noPriorResultsBB = 'Th17_48h_cut4_sA_p5_huA_bs50_w1_noTFA.mat';



outDir = fullfile(infModelDirLS,['PRfiguresAll'])
gompers
mkdir(outDir)




totInfResults = size(params,1);
legendInf = cell(totInfResults+1,1);
legendInf{1} = 'Random';
totGsGroups = size(gsPlotGroups,1);
% need to get some info on priors to start and so will load first results
load(fullfile(infModelDirLS,noPriorResultsLS))

% Plot random performance



for groupCount = groupInds
    currGsTitle = gsPlotGroups{groupCount,1};
    currGsOut = [gsPlotGroups{groupCount,2} addOutBit];
    currGroups = gsPlotGroups{groupCount,3};
    auprMax = gsPlotGroups{groupCount,4};
    [xx, gsInd] = intersect({gsInfs.figText},currGsTitle);
    currGS = gsInfs(gsInd).figText;
    figure(1), clf % P-R
    % plot random P-R ( : grey )
    plot([0 1],gsInfs(gsInd).randPR*[1 1],':','LineWidth',lineWidth+1,'Color',controlColor)
    hold on, grid on
    set(gca,'FontSize',fontSize,'FontWeight','Bold')
    xlabel('Recall','FontWeight','Bold','FontSize',fontSize)
    ylabel('Precision','FontWeight','Bold','FontSize',fontSize)
    axis([ 0 1 0 1])
    grid on, grid minor
    
    figure(2), clf % AUPR
    plot(gsInfs(gsInd).randPR*[1 1],[.4 totInfResults+1.6],'-','LineWidth',lineWidth+1,'Color',controlColor)
    hold on, grid on, grid minor
    set(gca,'FontSize',fontSize,'FontWeight','Bold')
    xlabel('AUPR','FontWeight','Bold','FontSize',fontSizeBar)

        
for tind = 1:totInfResults
    weight = params{tind,1};
    tfaOpt = params{tind,2};
    methodInf = params{tind,3};
    currMark = params{tind,4};
    currColor = params{tind,5};
    currWidth = params{tind,6};
    if find(ismember({'LS'},methodInf))
        infModelDir = infModelDirLS;
        noPriorResults = noPriorResultsLS;
        disp(methodInf)
    elseif find(ismember({'BB'},methodInf))
        infModelDir = infModelDirBB;
        noPriorResults = noPriorResultsBB;
        disp(methodInf)
    else
        error('Inference method not recognized as "BB" or "LS".')
    end
    if isnan(weight) % No Prior condition
        legendInf{tind+1} = ['No Prior, ' methodInf];
        prMat = fullfile(infModelDir,noPriorResults);
    else % Model with Prior
        priorName = priorNames{priorInd};
        priorFileText = priorFileTexts{priorInd};
        legendInf{tind+1} = [priorName ', ' methodInf ', w=' num2str(weight)];
        if find(ismember({'LS'},methodInf))
            prMat = fullfile(infModelDir,[priorFileText '_bias' num2str(100*weight) tfaOpt ]);
        else
            prMat = fullfile(infModelDir,[priorFileText '_bs50_w' strrep(num2str(weight),'.','p') tfaOpt]);
        end
    end 
    
    load(prMat)

    figure(1) % plot P-R 
    plot(gsInfs(gsInd).recalls,gsInfs(gsInd).precisions,...
        currMark,'LineWidth',currWidth,'Color',currColor)
    
    figure(2), % AUPR
    % add no Prior AUPR
    barh(totInfResults-tind+1,gsInfs(gsInd).auprs,...
        'FaceColor',currColor)

end

    figInf = fullfile(outDir,[currGsOut '_' priorFileText]);
    figure(1), % P-R
    
    for zoomRecall = zoomRecalls
        figure(1)
        title([priorName ' (' currGsTitle ')'])
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
    title([priorName ' (' currGsTitle ')'])
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