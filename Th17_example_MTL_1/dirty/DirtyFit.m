function [S, B, W, mEBIC, Lambdas] = DirtyFit(Xs,Ys, Slambdas, Blambdas, opts)
%% Adapted from: Castro, De Veaux, Miraldi, Bonneau "Multitask learning for joint
%   inference of gene regulatory networks form several expression datasets"
%% Goal: Find the optimal fit for the dirty multi-task model over a grid of parameters
% 
%   1/2n*sum_(k=1)^r [RSS^(k)] + lam_s * ||S||_(1,1) + lam_b * ||B||_(1,sup)
%               W = S + B  
%   n -- number of samples
%   r -- number of tasks
%   RSS -- residual sum of squares
%   lam_s/lam_b -- tuning parameters
%   S -- sparse matrix (unique to each model)
%   B -- block matrix (similar accross models)
%   ||S||_(1,1) = sum of the absolute values of all the elements of S
%   ||B||_(1,sup) = sum of the absolute values of the row maximum
%% Inputs: 
%   Xs: Cell of predictor matrices
%   Ys: Cell of response matrices
%   Slambdas: sparse lambdas
%   Blambdas: block sparse lambdas
%   opts: options. Currently includes
%       .pior: a vector of prior confidence for each coef. 0 indicates the
%       most confidence. This will be used as the initial beta vector if
%       'warm' is true
%       .warm: True if we use the previous fit to initialize our weight
%       vectors
%       .maxIter: maximum number of iterations
%       .tol: tolerance for RSS
%% Outputs:
%       Ss: Sparse weight matrices for each pair of tuning parameters (tps)
%       Bs: Block weight matrices for tps
%       mRSSs: mean RSS for tps  
%       n_nonzeros: number of nonzero parameters for W for tps
%       lambdas: 2 column matrix of tps
%% Reference: Jalali, Ali, et al. "A dirty model for multi-task learning."
% Advances in neural information processing systems. 2010.
%% Author: Peter DeWeirdt
%% Date:6/22/2018

ntasks = length(Ys);
npreds = size(Xs{1},2);
nlams = length(Slambdas) *length(Blambdas);
Ss = cell(1, nlams);
Bs = Ss;
mRSSs = Ss;
lambdas = Ss;

if ~isfield(opts, 'warm')
    opts.warm = 1;
end

if ~isfield(opts, 'maxIter')
    opts.maxIter = 1000;
end

if ~isfield(opts, 'tol')
    opts.tol = 1e-2;
end

if ~isfield(opts, 'penalty_factor')
    opts.penalty_factor = ones(npreds, 1);
end


%Calculate the covariance update terms
C, D = covarianceUpdateTerms(X,Y, ntasks, npreds); % Need to write

%Note: for warm start, we assume that lambdas are decreasing
for i = 1:length(Slambdas);
    lamS = Slambdas(i);
    if mod(i,2)
        iterVec = 1:length(Blambdas);
    else
        iterVec = length(Blambdas):-1:1;
    end
    for j = iterVec;
        lamB = Blambdas(j);
        lambdas(lami) = [lamS lamB];
        if lami == 1
            CurrB = zeros(npreds,ntasks);
            if opts.warm
                CurrS = repmat(ones(npreds,1) - opts.penalty_factor, 1, ntasks); % To Do: test whether starting with sparse or block first is better
            else
                CurrS = CurrB;
            end
            CurrW = CurrS + CurrB;
            LastW = ones(npreds,1);
        end
        whilei = 1;
        while (whilei < opts.maxIter) && (max(abs(LastW - CurrW)) > opts.tol)
            % Find our weight matrices
            CurrS = updateS(C, D, CurrB, CurrS, lamS, opts.penalty_factor); %to write
            CurrB = updateB(C, D, CurrB, CurrS, lamB, opts.penalty_factor); %to write
            LastW = CurrW;
            CurrW = CurrS + CurrB;
            whilei = whilei + 1;
        end
        B = CurrB;
        S = CurrS;
        W = B + S;
        RSSs = getRSSs(W, Xs, Ys);
        mEBIC = getmEbic(RSS, ,W)
    end
end




    