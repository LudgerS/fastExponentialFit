function [T2map, M0est] = T2analysis_withM0(data, TEs)
% Implementation of Pei et al., "Algorithm for Fast Monoexponential Fitting Based on
% Auto-Regression on Linear Operations (ARLO) of Data", 2015
% https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.25137
%
% Considers the corresponding erratum
% https://onlinelibrary.wiley.com/doi/10.1002/mrm.27807
%
% Fits the model 'data = M0*exp(-TE/T2)'
%
% Assumes that the last dimension of data is the time dimension. Should
% work for an arbitrary number of spatial dimensions.
%
% Can be applied just as well to any other mono-exponential 
% fit with equally spaced data points besides T2 relaxometry.
%
% Notation as in the publication
% m is the decaying signal
% deltaTE the spacing between echos
%
% Written by Ludger Starke; Max Delbrück Center for Molecular Medicine in
% the Helmholtz Association, Berlin; 21-03-31
%
% License: GNU GPLv3 
 

%% sort data
dim = size(data);

data = reshape(data, [prod(dim(1:(end-1))), dim(end)]);

[TEs, I] = sort(TEs);
data = data(:,I);

deltaTE = mean(TEs(2:end) - TEs(1:(end-1)));


%% compute T2 map
s = deltaTE/3*(data(:,1:(end-2)) + 4*data(:,2:(end-1)) + data(:,3:end));
delta = data(:,1:(end-2)) - data(:,3:end);
T2map = (sum(s.^2, 2) + deltaTE/3*sum(s.*delta, 2))./(deltaTE/3*sum(delta.^2, 2) + sum(s.*delta, 2));


%% restore original dimensions
T2map = reshape(T2map, dim(1:(end-1)));
T2map(T2map < 0) = 0;

data = reshape(data, dim);

%% estimate M0
R = exp(-repmat(permute(TEs(:), [2:length(dim), 1]), [dim(1:(end-1)), 1]) ./ repmat(T2map, [ones(1, length(dim)-1), numel(TEs)]));
M0est = sum(data.*R, length(dim)) ./ sum(R.^2, length(dim));
M0est(isnan(M0est)) = 0; 

