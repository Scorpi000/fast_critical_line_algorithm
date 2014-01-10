% TURNINGPOINTS  Calculates turning points of the efficient frontier using the CLA.
%   This function returns the NTP turning points of the efficient frontier
%   specified by ExpReturn and ExpCovariance of NASSETS risky assets. 
%   For each turning point the return, the standard deviation, the 
%   weights of the portfolio, and the Lagrange multiplier lambda is given.
%   The weights returned by TURNINGPOINTS can be slightly below the lower or 
%   slightly above the upper bound due to numerical imprecision. The help
%   of CHECK_WEIGHTS describes how to deal with such cases.
%   (See "Applying Markowitz's Critical Line Algorithm", by Andras Niedermayer 
%   and Daniel Niedermayer, 
%   http://www.vwi.unibe.ch/publikationen/download/dp0701.pdf for a
%   description of the algorithm and for an explanation of lambda.)
%
%   [mu_tp, sigma_tp, weights_tp, lambda_tp] = turningpoints(ExpReturn, ExpCovariance ...
%      [, LowerBound, UpperBound][, MaxNTP])
%
%   Inputs:
%     ExpReturn: NASSETSx1 vector with expected returns
%
%     ExpCovariance: NASSETSxNASSETS with the covariance matrix
%
%     LowerBound (optional): NASSETSx1 vector with allowed lower bounds for assets, default: 0
%
%     UpperBound (optional): NASSETSx1 vector with allower upper bounds, default: +infinity
%
%     MaxNTP (optional): maximal number of permissible turningpoints, if NTP exceeds
%     MaxNTP computation is interrupted and an error message is returned,
%     default value: 2*NASSETS
%
%   Outputs:
%     mu_tp: 1xNTP vector with expected return of each turning point
%
%     sigma_tp: 1xNTP vector with standard deviations
%
%     weights_tp: NASSETSxNTP matrix with the weights of the efficient
%     portfolio for each turning point
%
%     lambda_tp: 1xNTP vector with the Lagrange multipliers
%
%   See also FIND_POINT_ON_FRONTIER, FRONTCON_CLA, CHECK_WEIGHTS
%
% (C) 2006 Andras Niedermayer and Daniel Niedermayer

