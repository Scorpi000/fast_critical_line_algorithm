function [sig_port, wts_port] = find_point_on_frontier(mu, sigma, mu_port, mu_tps, sig_tps, wts_tps)
% FIND_POINT_ON_FRONTIER  Finds a point on the efficient frontier using the CLA.
%   Finds a point on the efficient frontier given the expected returns and
%   the covariance matrix of the NASSETS assets. The function returns the
%   minimal risk (standard deviation of the return) associated with the
%   given return of the portfolio. It also returns the weights of the
%   portfolio which achieves the minimal risk. Information about the
%   NTP turning points is (returned as output by TURNINGPOINTS) is needed.
%
%   [PortRisk, PortWeights] = find_point_on_frontier(ExpReturn, ExpCovariance, ...
%     PortReturn, mu_tps, sig_tps, wts_tps)
%
%   Inputs:
%     ExpReturn: NASSETSx1 vector of returns
%
%     ExpCovariance: NASSETSxNASSETS covariance matrix of assets
%
%     mu_tps: 1xNTP vector of expected returns of turning
%     points. Portfolios are sorted descending in expected returns.
%
%     sig_tps: 1xNTP vector of risks of turning points
%
%     wts_pts: NASSETSxNTP matrix of asset weights for each turning point
%
%   Output:
%     PortRisk: risk of the minimum variance portfolio with return
%     ExpReturn
%
%     PortWeights: weights of assets in the minimum variance portfolio
%
%   See also TURNINGPOINTS, FRONTCON_CLA
%
% (C) 2006 Andras Niedermayer and Daniel Niedermayer


epsilon=1e-11; % maximal tolerated error

ind = find(mu_tps<mu_port); % indices of portfolios with lower returns than mu_port

if length(ind)<1
    if abs(mu_port-mu_tps(length(mu_tps)))>epsilon
        error('Required portfolio return not feasible (too small).')
    end
    % if mu given corresponds to lowest possible return
    wts_port = wts_tps(:,length(mu_tps))';
elseif ind(1)==1
    if abs(mu_tps(1)-mu_port)>epsilon
        error('Required portfolio return not feasible (too large).')
    end
    % if mu given corresponds to highest possible return
    wts_port = wts_tps(:,1)';
else
	i=ind(1); % i: index of portfolio with higher return, i-1: ... with lower return
	alpha = (mu_port-mu_tps(i))/(mu_tps(i-1)-mu_tps(i));
    % interpolate weights according to mu
	wts_port = (alpha*wts_tps(:,i-1) + (1-alpha)*wts_tps(:,i))';
end
sig_port = sqrt(wts_port*sigma*wts_port');

