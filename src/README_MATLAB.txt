Matlab tool accompanying "Applying Markowitz' Critical Line Algorithm"
======================================================================
Andras Niedermayer and Daniel Niedermayer
May 9, 2006

The provided files offer an interface to the Fortran 90 code implementing
Markowitz' Critical Line Algorithm. The function FRONTCON_CLA computes the
points of the efficient frontier with the constraint that all asset 
weights have to be between lower and upper bounds specified for each
asset individually. It works the same way as FRONTCON in the 
Matlab Financial Toolbox  with the difference that it uses the critical 
line algorithm.  FRONTCON_CLA uses two functions: 
FIND_POINT_ON_FRONTIER and TURNINGPOINTS. FIND_POINT_ON_FRONTIER does what 
its name suggests. TURNINGPOINTS is the Matlab interface to the Fortran 90 
code implementing the algorithm. It finds all turning points (points where 
an asset comes into the portfolio or gets out of the portfolio). The 
algorithm is described in the paper (see below). We further provide the
function CHECK_WEIGHTS to check whether numerical imprecisions have lead
to small violations of asset bounds (see comments in check_weights.m).

Usage
-----

Copy all files from Matlab_CLA.zip to the working directory of Matlab. You
can compute the efficient frontier by calling

>> frontcon_cla(mu, sigma)

where mu is the vector of expected returns and sigma is the covariance 
matrix. Here the default values 0 and infinity are used for the lower and upper bounds,
respectively. You can get further information by typing

>> help frontcon_cla

or looking at the source code. Of course you can also access the functions
FIND_POINT_ON_FRONTIER and TURNINGPOINTS directly.

You may also want to take a look at the M-file test_cla.m which tests the 
implementation of the algorithm and compares the results with the standard
implementation in the financial toolbox. (You need to have the financial 
toolbox installed in order to run the test_cla.)

System Requirements
-------------------

The software has been tested with Matlab versions 6.1, 6.5, and 7.0 on Windows 
XP.  However,  it should work with higher versions of Matlab and other Windows 
versions as well. We also provide a tool for Matlab version 7.3 on 32-bit Linux 
(this seems not to run on lower versions of Matlab on Linux).

Description of Algorithm
------------------------

A description of the algorithm implemented  (a variation of Markowitz' 
Critical Line Algorithm)  is provided in the article "Applying Markowitz' 
Critical Line Algorithm" in the Discussion Paper Series of the University 
of Bern. The article is available for download at
http://www.vwi.unibe.ch/publikationen/download/dp0701.pdf