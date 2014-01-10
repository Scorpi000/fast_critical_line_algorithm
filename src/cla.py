"""Python interface to fortran code computing the minimum variance frontier
with the Critical Line Algorithm.

There are three ways of accessing the fortran subroutines:
1) object oriented through the class PortfolioChoice
2) structural through the functions turningpoints, find_points_on_frontier,
   frontcon_cla (in analogy to the Matlab interface provided to our algorithm
   and to the Matlab Financial Toolbox function FRONTCON)
3) direct access to the fortran subroutines calcturningpoints (constraint
   w>=0) and calcturningpointsgen (general constraint lowerbound<=w<=upperbound).

The doc strings and the source code contain more information about the class and
the functions. To see how to access the fortran functions directly, type
 print cla.calcturningpoints.__doc__
and
 print cla.calcturningpointsgen.__doc__

Sample usage:

>>> mu = [1., 2., 3.]
>>> sigma = [[3., 0., 0.], [0., 2., 0.], [0., 0., 1.]]
>>> lower = [0.,0.,0.]
>>> upper = [1.1,1.1,1.1]
>>> # object-oriented
>>> pc = PortfolioChoice(mu, sigma, lower, upper)
>>> print pc.mu_tps
[ 3.          2.8         1.42857143  1.        ]
>>> pc.compute_frontier(NumPorts = 5)
>>> print pc.PortRisk
[ 1.73205081  1.10867789  0.81649658  0.75        1.        ]
>>> print pc.PortReturn
[ 1.   1.5  2.   2.5  3. ]
>>> # non-object-oriented
>>> PortRisk, PortReturnOut, PortWts = frontcon_cla(mu, sigma, NumPorts=5)
>>> print PortRisk
[ 1.73205081  1.10867789  0.81649658  0.75        1.        ]


Also see the "if __name__='__main__'" at the end of the file for
usage examples.

(See "Applying Markowitz's Critical Line Algorithm", by Andras Niedermayer 
and Daniel Niedermayer, 
http://www.vwi.unibe.ch/publikationen/download/dp0701.pdf for a
description of the algorithm and for an explanation of lambda.)
"""

import numpy as N
from cla_fortran import calcturningpoints, calcturningpointsgen

class CLAError(Exception): pass

def turningpoints(ExpReturn, ExpCovariance, LowerBound=None, UpperBound=None, MaxNTP=None):
    """ mu_tp, sigma_tp, weights_tp, lambda_tp = turningpoints(ExpReturn, ExpCovariance,
                           LowerBound=None, UpperBound=None, MaxNTP=None)
        Returns mu, sigma, the weights, and lambda for each turning point on the minimum
        variance frontier. The default value for LowerBound is zeros(nn), for UpperBound ones(nn),
        and for MaxNTP 2*nn, where nn is the number of assets (=len(ExpReturn)).
        
        Parameters:
            ExpReturn: any object convertible to a numpy array of dimension (nn,1)
            ExpCovariance:              -''-                                (nn,nn)
            LowerBound:                 -''-                                (nn,1)
                        vector of lower bounds for asset weights,
                        default value: zeros(nn)
            UpperBound:                 -''-                                (nn,1)
                        vector of upper bounds for asset weights,
                        default value: 2.*ones(nn)
            MaxNTP:     maximum space allocated for saving turning points
                        default value: 2*nn
    """
    nn = len(ExpReturn)
    if MaxNTP is None:
        MaxNTP = 2*nn
    # if simple version of algorithm is to be used (LowerBound=0, UpperBound=1)
    if LowerBound is None and UpperBound is None:
        mu_tp, sigma_tp, weights_tp, lambda_tp, n_tp, error_flag = calcturningpoints(
                                                            ExpReturn, ExpCovariance, MaxNTP)
    else: # if general version of algorithm is to be used
        if LowerBound is None:
            LowerBound = N.zeros(nn)
        if UpperBound is None:
            UpperBound = 2.*N.ones(nn)
        mu_tp, sigma_tp, weights_tp, lambda_tp, n_tp, error_flag = calcturningpointsgen(
                                                            ExpReturn, ExpCovariance,
                                                            LowerBound, UpperBound, MaxNTP)
    if error_flag == 4:
        raise CLAError('turningpoints: Error during execution: '
                        'maximal space allocated for turning points exceeded. '
                        'Try increasing MaxNTP')
    if error_flag != 0:
        raise CLAError('turningpoints: Unknown error. (Error code %d)' % error_flag)
    # only the first n_tp elements of the returned arrays contain the turning points,
    # the rest is filled up with zeros
    weights_tp.resize(nn, n_tp)
    for vector in mu_tp, sigma_tp, lambda_tp:
        vector.resize(n_tp, refcheck=False) # we don't need refcheck since the only
                                           # extra reference is vector
    #~ return mu_tp[:n_tp], sigma_tp[:n_tp], weights_tp[:,:n_tp], lambda_tp[:n_tp] 
    return mu_tp, sigma_tp, weights_tp, lambda_tp 

#function [sig_port, wts_port] = find_point_on_frontier(mu, sigma, mu_port, mu_tps, sig_tps, wts_tps)
def find_points_on_frontier(mu, sigma, mu_port, mu_tps, sig_tps, wts_tps):
    """Finds points on the efficient frontier using the CLA.
    Finds points on the efficient frontier given the expected returns and
    the covariance matrix of the NASSETS assets. The function returns the
    minimal risk (standard deviation of the return) associated with the
    given return of each portfolio. It also returns the weights of the
    portfolios which achieves the minimal risk. Information about the
    NTP turning points is (returned as output by turningpoints) is needed.

    PortRisk, PortWeights = find_point_on_frontier(ExpReturn, ExpCovariance, 
      PortReturn, mu_tps, sig_tps, wts_tps)

    Parameters:
      ExpReturn: NASSETSx1 vector of returns

      ExpCovariance: NASSETSxNASSETS covariance matrix of assets
      
      PortReturn: NPORTSx1 matrix of returns of the portfolios

      mu_tps: 1xNTP vector of expected returns of turning
      points. Portfolios are sorted descending in expected returns.

      sig_tps: 1xNTP vector of risks of turning points

      wts_pts: NASSETSxNTP matrix of asset weights for each turning point

    Returns:
      PortRisk: risk of the minimum variance portfolio with return
      ExpReturn

      PortWeights: weights of assets in the minimum variance portfolio
    """
    epsilon = 1e-11 # maximal tolerated error
    mu_port = N.array(mu_port)
    mu_tps_reverse = mu_tps[::-1].copy()
    i = len(mu_tps)-mu_tps_reverse.searchsorted(mu_port) # find turning point to the right of mu_port
    # i: index of portfolio with lower return, i-1: ... with higher return
    toosmall = i == len(mu_tps)
    toolarge = i == 0
    if N.any(toosmall):
        if N.any(abs(mu_port[toosmall]-mu_tps[-1])>epsilon):
            raise CLAError('Required portfolio return(s) not feasible (too small).\n'
                            'min mu_tps: %s, required portfolio return(s): %s' % (mu_tps[-1], mu_port[toosmall]))
        # if mu given corresponds to lowest possible return
        i[toosmall] = len(mu_tps)-1
    if N.any(toolarge):
        if N.any(abs(mu_tps[0]-mu_port[toolarge])>epsilon):
            raise CLAError('Required portfolio return(s) not feasible (too large).\n'
                            'max mu_tps: %s, required portfolio return(s): %s' % (mu_tps[0], mu_port[toolarge]))
        # if mu given corresponds to highest possible return
        i[toolarge] = 1
    alpha = (mu_port-mu_tps[i])/(mu_tps[i-1]-mu_tps[i])
    # interpolate weights according to mu
    wts_port = (alpha*wts_tps[:, i-1] + (1-alpha)*wts_tps[:, i])
    # take sqrt of diagonal elements of w*sigma*w'
    sig_port = N.sqrt(N.sum(N.dot(wts_port.T, sigma)*wts_port.T, axis=1))
    return sig_port, wts_port

def frontcon_cla(ExpReturn, ExpCovariance, NumPorts=10, PortReturn=None, AssetBounds=None, 
                 Plot=False, PlotTurningpoints=False):
    """Computes the mean-variance efficient frontier using the CLA.
    Computes the mean-variance efficient frontier with lower and upper bounds
    for asset weights. Returns risk, return and asset weights for each
    point requested on the frontier. Points are requested by either
    specifying the returns of the points or by specifying the number of
    points on the frontier that should be computed. In the latter case
    points are equally spaced between the minimal and the maximal
    possible return. This function can be used like FRONTCON in Matlab's
    financial toolbox with the main differences that the parameters 
    Groups and GroupBounds cannot be specified and that the Critical Line
    Algorithm is used.

    PortRisk, PortReturn, PortWts = frontcon(ExpReturn, ExpCovariance, NumPorts=10, 
            PortReturn=None, AssetBounds=None, Plot=False, PlotTurningpoints=False)

    Parameters: 
      ExpReturn: 1xNASSETS vector with the expected return of each asset.
     
      ExpCovariance: NASSETSxNASSETS covariance matrix
     
      NumPorts (optional): number of efficient portfolios, NPORTS.

      PortReturn (optional): NPORTSx1 vector with the target return values on the frontier.
      Default value: NPORTS returns equidistantly distributed between the minimal and the
      maximal achievable return. If contradicting values for PortReturn and NumPorts are
      provided, NPORTS is set to the size of PortReturn.

      AssetBounds (optional): 2xNASSETS matrix containing the lower and
      upper bounds on the portfolio weights. Default values: 0 for lower
      and +infinity for upper bounds (i.e. no upper bounds)

    Returns: 
      PortRisk: NPORTSx1 vector of the standard deviation of return for
      portfolio on the efficient frontier. 
     
      PortReturn: NPORTSx1 vector of the expected return of portfolios on the 
      efficient frontier.
     
      PortWts: NPORTSxNASSETS matrix of asset weights of portfolios on the
      efficient frontier. The weights of each portfolio sum up to 1.
             

    Notes: 
      As in FRONTCON the efficient frontier is plotted if no output arguments
      are specified.


    See also FRONTCON (from the Matlab Financial Toolbox).
"""
    # convert to type array
    ExpReturn = N.array(ExpReturn, copy=False)
    ExpCovariance = N.array(ExpCovariance, copy=False)
    # compute the turning points
    if AssetBounds is None: # if no lower and upper bounds given
        mu_tps, sig_tps, wts_tps, lambda_tps = turningpoints(ExpReturn.T, ExpCovariance)
    else:
        AssetBounds = N.array(AssetBounds, copy=False)
        mu_tps, sig_tps, wts_tps, lambda_tps = turningpoints(ExpReturn.T, ExpCovariance, 
                                                  AssetBounds[0,:].T, AssetBounds[1,:].T)
    # if portfolio returns not specified, set them
    if PortReturn is None:
        PortReturn = _create_port_return(NumPorts, mu_tps)
    # compute the points on the efficient frontier
    PortRisk, PortWts = find_points_on_frontier(ExpReturn.T, ExpCovariance, PortReturn, mu_tps, sig_tps, wts_tps)
    # Plot the efficient frontier if required
    if Plot or PlotTurningpoints:
        _plot_frontier(PortRisk, PortReturn, mu_tps, sig_tps, Plot, PlotTurningpoints)
        
    return PortRisk, PortReturn, PortWts

def _create_port_return(NumPorts, mu_tps):
    max_mu_tps = mu_tps[0]  # turning points are given in descending order in mu
    min_mu_tps = mu_tps[-1]
    if NumPorts > 1:
        stepsize = (max_mu_tps-min_mu_tps)/(NumPorts-1)
        PortReturn = min_mu_tps + stepsize*N.arange(NumPorts, dtype='d')
    else:
        PortReturn = N.array([max_mu_tps])
    return PortReturn

def _plot_frontier(PortRisk, PortReturn, mu_tps, sig_tps, Plot, PlotTurningpoints):
    import pylab as P
    P.title('Mean-Variance Efficient Frontier')
    P.xlabel('Risk (Standard Deviation)')
    P.ylabel('Expected Return')
    if Plot:
        P.plot(PortRisk, PortReturn)
        if PlotTurningpoints:
            P.plot(sig_tps, mu_tps,'o')
    elif PlotTurningpoints:
        P.plot(sig_tps, mu_tps,'o:')
    P.show()


class PortfolioChoice(object):
    """Object representing a portfolio choice problem.
    Provides the turning points of the efficient frontier and also computes
    points on the frontier.
    The first access to any of the properties mu_tps, sigma_tps, weights_tps,
    or lambda_tps triggers the computation of the turningpoints. Points
    on the frontier can be computed with compute_frontier, information about
    the points accessed through the properties PortRisk, PortReturn, PortWts
   """
    __tps_updated = False
    __frontier_updated = False

    __mu_tps, __sigma_tps, __weights_tps, __lambda_tps = 4*(None,)
    
    # read-only properties
    @property
    def mu_tps(self): self._update_tps(); return self.__mu_tps
    @property
    def sigma_tps(self): self._update_tps(); return self.__sigma_tps
    @property
    def weights_tps(self): self._update_tps(); return self.__weights_tps
    @property
    def lambda_tps(self): self._update_tps(); return self.__lambda_tps
    @property
    def PortRisk(self): self._update_frontier(); return self.__PortRisk
    @property
    def PortReturn(self): self._update_frontier(); return self.__PortReturn
    @property
    def PortWts(self): self._update_frontier(); return self.__PortWts

    def __init__(self, mu, sigma, LowerBound=None, UpperBound=None):
        """mu: expected returns of assets  (nn vector)
        sigma: covariance matrix of assets (nn x nn matrix)
        the default values for LowerBound and UpperBound are zeros(nn)
        and 2.*ones(nn), respectively"""
        self.mu = mu
        self.sigma = sigma
        self.lower = LowerBound
        self.upper = UpperBound
        

    def _update_tps(self):
        if not self.__tps_updated:
            self.__mu_tps, self.__sigma_tps, self.__weights_tps, self.__lambda_tps = \
                            turningpoints(self.mu, self.sigma, self.lower, self.upper)
        self.__tps_updated = True
        
    def _update_frontier(self):
        if not self.__frontier_updated:
            self.compute_frontier()
        
    def compute_frontier(self, NumPorts=10, PortReturn=None):
        """Computes the efficient frontier for the points with returns specified
        in the vector PortReturn. If PortReturn is not specified, NumPorts points
        are taken between the highest and the lowest possible returns.
        Results are accessible through the properties PortRisk, PortReturn, and
        PortWts."""
        self._update_tps()
        # if portfolio returns not specified, set them
        if PortReturn is None:
            PortReturn = _create_port_return(NumPorts, self.__mu_tps)
        self.__PortReturn = PortReturn
        # compute the points on the efficient frontier
        self.__PortRisk, self.__PortWts = find_points_on_frontier(self.mu, self.sigma, PortReturn, self.__mu_tps, 
                                self.__sigma_tps, self.__weights_tps)
        self.__frontier_updated = True

    def plot(self, PlotFrontier=True, PlotTurningpoints=False):
        self._update_tps()
        self._update_frontier()
        if PlotFrontier or PlotTurningpoints:
            _plot_frontier(self.__PortRisk, self.__PortReturn, self.__mu_tps, self.__sigma_tps, 
                           PlotFrontier, PlotTurningpoints)

def oo_example():
    from numpy import zeros, ones
    #mu=array([1.,2.,3.])
    mu = [1., 2., 3.]
    #sigma=array([[3.,0.,0.],[0.,2.,0.],[0.,0.,1.]])
    sigma = [[3., 0., 0.], [0., 2., 0.], [0., 0., 1.]]
    nn = len(mu)
    lower = zeros(nn)
    upper = 1.01*ones(nn)
    pc = PortfolioChoice(mu, sigma, lower, upper)

    print pc.mu_tps
    print pc.sigma_tps
    print pc.weights_tps
    print pc.lambda_tps

    pc.compute_frontier(NumPorts=100)

    print pc.PortRisk
    print pc.PortReturn
    print pc.PortWts

    pc.plot(PlotTurningpoints=True)

def non_oo_example():
    #from numpy import zeros, ones
    #mu=array([1.,2.,3.])
    mu = [1., 2., 3.]
    #sigma=array([[3.,0.,0.],[0.,2.,0.],[0.,0.,1.]])
    sigma = [[3., 0., 0.], [0., 2., 0.], [0., 0., 1.]]
    #nn = len(mu)
    #lower=zeros(nn)
    #upper=1.01*ones(nn)

    PortRisk, PortReturn, PortWts = frontcon_cla(mu, sigma, Plot=True, PlotTurningpoints=True)

    print PortRisk
    print PortReturn
    print PortWts

def test_cla():
    import time
    
    n_assets = 100      # number of assets
    n_portfolios = 50  # number of portfolios (points) on the efficient frontier


    print 'Number of assets: ', n_assets

    # random covariance matrix sigma is created
    print 'creating covariance matrix...'
    tic = time.time()
    x = N.random.rand(n_assets, n_assets+5)
    sigma = N.dot(x, x.T)
    toc = time.time()-tic
    print 'Time for creation of random covariance matrix: ', toc

    # random vector of expected returns is created
    mu = N.random.rand(1, n_assets)

    # compute efficient frontier with critical line algorithm
    print 'computing efficient frontier with critical line algorithm...'
    tic = time.time()
    ourrisk, ourreturn, ourwts = frontcon_cla(mu, sigma, NumPorts=n_portfolios)
    toc = time.time()-tic
    print 'Computation time for whole efficient frontier with critical line algorithm: ', toc

    # compute a point (with mu_p = (max(mu)-min(mu))/2)
    # on the efficient frontier with cvxopt if installed
    try:
        # (un)comment following line if necessary
        #~ raise ImportError("fake exception: don't use cvxopt")
        import cvxopt.cvxprog as cvx
        import cvxopt.base as B
        # use mosek solver if installed
        try:
            # (un)comment following line if necessary
            #~ raise ImportError("fake exception: don't use pymosek")
            import pymosek
            print 'computing efficient frontier with cvxopt (mosek solver)'
            cvx.options['MOSEK'] = {pymosek.iparam.log:0}
            solver_used = 'mosek'
        except ImportError:
            print 'computing efficient frontier with cvxopt (default solver)'
            cvx.options['show_progress'] = False
            solver_used = None
        # the solver solves:
        # min (1/2)*x'*PP*x + q'*x
        # s.t. G*x <= h
        #      A*x = b
        PP = sigma
        q = N.zeros((n_assets, 1))
        G = -N.eye(n_assets)
        h = N.zeros((n_assets, 1))
        A = N.vstack((N.ones((1, n_assets)),
                      mu))
        b = N.array([[1.],
                    [(N.max(mu)-N.min(mu))/2]])
        tic = time.time()
        retval = cvx.qp(B.matrix(PP), B.matrix(q), B.matrix(G), B.matrix(h), 
                        B.matrix(A), B.matrix(b), solver=solver_used)
        x = retval['x']
        mucvx = B.matrix(mu)*x
        sigcvx = N.sqrt(x.T*B.matrix(sigma)*x)
        toc = time.time()-tic
        print 'Computation time for one point on efficient frontier with cvxopt: ', toc
    except ImportError:
        print "(you can install the python package 'cvxopt' to make a timing comparison)"
        mucvx = None

    # plot the our results
    try:
        import pylab as P
        P.plot(ourrisk, ourreturn) #, 'r')
        if mucvx is not None:
            P.plot(sigcvx, mucvx,'o')
        P.title('CLA results for efficient frontier')
        P.show()
    except ImportError:
        print "(you can install the python package 'matplotlib' to plot results)"
    
def _test():
    import doctest
    doctest.testmod()

if __name__ == '__main__':
    #oo_example()
    #non_oo_example()
    test_cla()
    #_test()
