function [confs, rnk] = calc_confs(B, A0, X, Y)
%% Description:
% Caclulate confidence scores for each entry in B, beta. Confidence
% scores are calculated as:
%           1 - sigma^2/sigma~i^2.
% Where sigma^2 is the variance of the residuals for the whole model, and
% sigma~i^2 is the variance of the residuals of the model without the ith
% predictor. 
%% Inputs:
% B - coefficient matrix
% A0 - intercept
% X - predictor matrix
% Y - response matrix
%% Output:
% confs - the confidence score for each coefficient
% rnk - vector of ranks for each edge. 1 indicates the coefficient we are
% most confident in. Gives all of the 0 confidence scores the lowest
% possible rank. 
%% Author:
% Peter DeWeirdt - Cincinnati Children's Summer Intern
%% Date: 6/13/2018

sigma = var((X*B + A0) - Y);
nzero = find(B)';
confs = zeros(size(B, 1),1);
for i = nzero
    Bno_i = B;
    Bno_i(i) = 0;
    sigma_i =  var((X*Bno_i + A0) - Y);
    confs(i) = 1 - (sigma/sigma_i);
end
[~,~,rev_rnk] = unique(confs);
rnk = rev_rnk + (max(rev_rnk) - rev_rnk) - (rev_rnk - min(rev_rnk));
end
