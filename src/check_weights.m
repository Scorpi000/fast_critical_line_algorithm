function [correct] = check_weights(weights, lowerbounds, upperbounds, epsilon)
% CHECK_WEIGHTS. Checks whether the weights are within their upper and
%   lower bounds. A deviation of max. epsilon is allowed.
%
%   The weights returned by TURNINGPOINTS can be slightly below the lower or 
%   slightly above the upper bound due to numerical imprecision. You can
%   check this by running 
%
%   correct = check_weights(weights_tp, LowerBound, UpperBound, epsilon)
%
%   with e.g. epsilon=1e-10 being the error tolerance and weights_tp, LowerBound, and UpperBound
%   as described in the help of TURNINGPOINTS. If this results in correct = 0
%   you may have to rerun turningpoints with slightly modified ExpReturn or
%   ExpCovariance to get rid of the numerical imprecisions.
%
%   Inputs:
%     weights_tp: NASSETSxNTP matrix with the weights of the efficient
%     portfolio for each turning point (NASSETS and NTP as described in
%     TURNINGPOINTS)
%
%     LowerBound: NASSETSx1 vector with allowed lower bounds for assets, default: 0
%
%     UpperBound: NASSETSx1 vector with allower upper bounds, default: +infinity
%
%     epsilon: maximum allowed violation of bounds
%
%   Outputs:
%     correct: 1 if no weight violates the bounds by more than epsilon, 
%     0 otherwise
%
%   See also TURNINGPOINTS
%
% (C) 2006 Andras Niedermayer and Daniel Niedermayer


correct = 1;

for i=1:size(weights,2)
    if any(weights(:,i)<lowerbounds-epsilon) | any(weights(:,i) > upperbounds+epsilon)
        correct = 0;
    end
end
