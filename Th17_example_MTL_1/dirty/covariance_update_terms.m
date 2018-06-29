function [C, D] = covariance_update_terms(X,Y)
%% Adapted from: Castro, De Veaux, Miraldi, Bonneau "Multitask learning for joint
%   inference of gene regulatory networks form several expression datasets"
%% Goal: caclulate terms for covariance update for OLS fit
%% Outputs:
%        C: transpose(X_j)*Y for each feature j
%        D: transpose(X_j)*X_l for each feature j for each feature l
%% Reference: Friedman, Hastie, Tibshirani, 2010 in Journal of Statistical Software
%        Regularization Paths for Generalized Linear Models via Coordinate Descent.
%% Author: Peter DeWeirdt
%% Date:6/22/2018

