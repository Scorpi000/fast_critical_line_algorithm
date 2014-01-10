Python module accompanying "Applying Markowitz' Critical Line Algorithm"
========================================================================
Andras Niedermayer and Daniel Niedermayer
January 13, 2008

The provided files offer an interface to the Fortran 90 code implementing
Markowitz' Critical Line Algorithm. Note that this Python module is currently 
only available for Linux. A Matlab interface for both Windows and Linux is
available upon request from the authors.

Prerequisites
-------------

You need the packages "python-numpy" and "libgfortran2". If you want to use the
plotting functions in the module, you need "python-matplotlib". This can be
installed on Ubuntu with
 $ sudo apt-get update
 $ sudo apt-get install python-numpy libgfortran2 python-matplotlib
The universe repository has to be enabled.

You can also install the package "python-cvxopt" to compare our runtimes with 
a convex optimization toolbox. A commercial optimization toolbox that cvxopt 
can plug in to is available at http://www.mosek.com.


Usage
-----

Copy the files from this archive to your current working directory (or to a
a directory PYTHONPATH points to).

To see an example of the usage of the module, type
 $ python cla.py

A description of the functions and classes is provided in the doc strings in
cla.py. Alternatively, type
>>> import cla
>>> help(cla)
>>> help(cla.frontcon_cla)
in the Python interpreter.

System Requirements
-------------------

The 32-bit version of this software has been tested on Ubuntu 7.10, Python 2.5.1,
numpy 1.0.3. The 64-bit version on Gentoo Base System release 1.12.10, Python 2.4.4,
numpy 1.0.4. However, it should work on other Linux/Python/numpy distributions
and versions as well.


Description of Algorithm
------------------------

A description of the algorithm implemented  (a variation of Markowitz' 
Critical Line Algorithm)  is provided in the article "Applying Markowitz' 
Critical Line Algorithm" in the Discussion Paper Series of the University 
of Bern. The article is available for download at
http://www.vwi.unibe.ch/publikationen/download/dp0701.pdf
