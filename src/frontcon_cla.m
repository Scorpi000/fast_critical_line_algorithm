function [PortRisk, PortReturnOut, PortWts] = frontcon_cla(ExpReturn, ExpCovariance, NumPorts, PortReturn, AssetBounds)
% FRONTCON_CLA  Computes the mean-variance efficient frontier using the CLA.
%   Computes the mean-variance efficient frontier with 0 as a lower bound
%   for asset weight and returns risk, return and asset weights for each
%   point requested on the frontier. Points are requested by either
%   specifying the returns of the points or by specifying the number of
%   points on the frontier that should be computed. In the latter case
%   points are equally spaced between the minimal and the maximal
%   possible return. This function can be used like FRONTCON in the
%   financial toolbox with the difference that the parameters 
%   Groups and GroupBounds cannot be specified and that the Critical Line
%   Algorithm is used.
%
%   [PortRisk, PortReturn, PortWts] = frontcon(ExpReturn, ExpCovariance, ...
%      [, NumPorts[, PortReturn[, AssetBounds]]])
%
%   Inputs: 
%     ExpReturn: 1xNASSETS vector with the expected return of each asset.
%    
%     ExpCovariance: NASSETSxNASSETS covariance matrix
%    
%     NumPorts (optional): number of efficient portfolios, NPORTS.  If entered as the
%     empty matrix [] or not specified the default value is entered.
%
%     PortReturn (optional): NPORTSx1 vector with the target return values on the frontier.
%
%     AssetBounds (optional): 2xNASSETS matrix containing the lower and
%     upper bounds on the portfolio weights. Default values: 0 for lower
%     and +infinity for upper bounds (i.e. no upper bounds)
%
%   Outputs: 
%     PortRisk: NPORTSx1 vector of the standard deviation of return for
%     portfolio on the efficient frontier. 
%    
%     PortReturn: NPORTSx1 vector of the expected return of portfolios on the 
%     efficient frontier.
%    
%     PortWts: NPORTSxNASSETS matrix of asset weights of portfolios on the
%     efficient frontier. The weights of each portfolio sum up to 1.
%            
%
%   Notes: 
%     As in FRONTCON the efficient frontier is plotted if no output arguments
%     are specified.
%
%
%   See also FRONTCON (from the Matlab Financial Toolbox), FIND_POINT_ON_FRONTIER, 
%   TURNINGPOINTS.
%
%   A description of the algorithm used to compute the turning points is
%   given in "Applying Markowitz's Critical Line Algorithm", by Andras Niedermayer 
%   and Daniel Niedermayer, 
%   http://www.vwi.unibe.ch/publikationen/download/dp0701.pdf
%
% (C) 2006 Andras Niedermayer and Daniel Niedermayer


if nargin<2
    error('At least two arguments needed. Type "help frontcon_cla" for usage description.')
end

n_assets = length(ExpReturn);

% if number of portfolios not specified set the default value
if nargin<3 | (nargin<4 & all(size(NumPorts)==[0 0]))
    NumPorts = 10; % default value
end

% compute the turning points
if nargin<5 % if no lower and upper bounds given
    [mu_tps, sig_tps, wts_tps] = turningpoints(ExpReturn', ExpCovariance);
else
    [mu_tps, sig_tps, wts_tps] = turningpoints(ExpReturn', ExpCovariance, AssetBounds(1,:)', AssetBounds(2,:)');
end
    
n_turningpoints = length(mu_tps);

% if portfolio returns not specified, set them
if nargin<4
    max_mu_tps = mu_tps(1);  % turning points are given in descending order in mu
    min_mu_tps = mu_tps(n_turningpoints);
    if NumPorts>1
        stepsize=(max_mu_tps-min_mu_tps)/(NumPorts-1);
        PortReturn = min_mu_tps:stepsize:max_mu_tps;
    else
        PortReturn = max_mu_tps;
    end
else
    NumPorts = length(PortReturn);
end

PortRisk=zeros(size(PortReturn));
PortWts=zeros(NumPorts, n_assets);
% compute the points on the efficient frontier
for i=1:NumPorts
    [PortRisk(i), PortWts(i,:)] = find_point_on_frontier(ExpReturn', ExpCovariance, PortReturn(i), mu_tps, sig_tps, wts_tps);
end

PortReturnOut = PortReturn;

% Without output parameters the efficient frontier is plotted
if nargout < 1
    plot(PortRisk, PortReturn);
    grid on;
    title('Mean-Variance Efficient Frontier', 'Color', 'k');
    xlabel('Risk (Standard Deviation)');
    ylabel('Expected Return');
    clear PortRisk;
end
