% Tests our implementation of the critical line algorithm and compares it
% with the standard implementation in Matlab's financial toolbox. A plot is
% produced showing the efficient frontier calculated with the standard
% toolbox and the critical line algorithm.
% The variables n_assets and n_portfolios can be changed to test with
% different numbers of assets and different numbers of points on the
% frontier.
%
% (C) 2006 Andras Niedermayer and Daniel Niedermayer

n_assets=20;      % number of assets
n_portfolios=10;  % number of portfolios (points) on the efficient frontier


disp('Number of assets: '); disp(n_assets);

% random covariance matrix sigma is created
disp('creating covariance matrix...');
tic;
sigma=zeros(n_assets, n_assets);

for i=1:n_assets+5
    x=rand(n_assets,1);
    sigma=sigma+x*x';
end
creationtime=toc;
disp('Time for creation of random covariance matrix: '); disp(creationtime);

% random vector of expected returns is created
mu=rand(1,n_assets);

% compute efficient frontier with standard toolbox
disp('computing efficient frontier with standard toolbox...');
tic;
[stdrisk, stdreturn, stdwts] = frontcon(mu, sigma, n_portfolios);
standard_computation_time=toc;
disp('Computation time standard toolbox: '); disp(standard_computation_time);

% compute efficient frontier with critical line algorithm
disp('computing efficient frontier with critical line algorithm...');
tic;
[ourrisk, ourreturn, ourwts] = frontcon_cla(mu, sigma, [], stdreturn);
our_computation_time=toc;
disp('Computation time critical line algorithm: '); disp(our_computation_time);

% plot the 2 results
plot(stdrisk, stdreturn, 'b+', ourrisk, ourreturn, 'rO');
title('Comparison of standard (+) and CLA (O) results for efficient frontier');
grid on;