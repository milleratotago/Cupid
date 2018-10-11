# Cupid
Object-oriented MATLAB toolbox for working with probability distributions

[1]  Overview
------------------------

This toolbox can be helpful for working with univariate probability distributions.
For any given distribution, it can compute all of the standard quantities: density and
cumulative density functions, hazard function, mean, variance, skewness, kurtosis, etc,
and it can generate random numbers.
It can also estimate parameter values (i.e., by maximum likelihood and
several other methods).
Parameters estimates can be constrained within a single distribution, and it is possible
to estimate the parameters of multiple distributions jointly while specifying constraints
among the parameters of different distributions.
More than 40 standard distributions (e.g., normal, exponential, gamma, ...) have already
been implemented, and it is easy to add new ones.
(Any distribution for which you know either the PDF or the CDF can be added; everything
else can be computed numerically from either of those functions.)

A somewhat novel feature of CUPID is the ability to create what might be called
"derived" distributions.
For example, you can form a new distribution by truncating one of the standard ones
(e.g., a standard normal distribution restricted to the range of -2 to +3).
Within the object-oriented framework, this truncated distribution is simply a new
distribution object, so its functions can be computed and parameters estimated
just like those of the standard distributions.
In addition to truncation, you can form derived distributions by transformations of
standard distributions (e.g., linear, inverse, exponential, power, or log).
You can also form derived distribution that involve independent random variables
from two or more standard distributions.
For example, this can be done by convolution (i.e., the distribution of the sum
of random variables from two or more standard distributions),
by a probability mixture of two or more distributions,
or by an order statistic of two or more random variables
[e.g., the minimum a standard normal and a uniform(0,1)].

Moreover, the object-oriented framework is fully recursive, so you can form derived
distributions from other derived distributions rather than from standard ones (e.g., form
the distribution of 1/sqrt(x), where x is uniform over some range, or form the
convolution of a mixture distribution and a log-transformed distribution).
All of the same functions (e.g., PDF) and operations (e.g., parameter estimation)
that are available for the standard distributions are also available for the
derived distributions, although the required numerical computations might
be much slower and less accurate.

[2]  Standard Distributions
---------------------------

Here are most of the standard probability distributions that have been implemented so far
(several alternative parameterizations are available for some of them):

* Beta
* Cauchy
* Chi
* ChiSq
* ChiSqNoncentral
* DblMon (double monomial)
* ExGaussian (sum of exponential and normal)
* Exponential
* ExpSum (sum of exponentials)
* ExpSumT (truncated sum of exponentials)
* ExtrVal1 (extreme value, type 1)
* ExtrVal2 (extreme value, type 2)
* ExWald (sum of exponential and Wald)
* F
* FNoncentral
* Gamma
* GenNor1 (general error distribution, type 1)
* GenNor2 (general error distribution, type 2)
* HyperbolicTan
* JohnsonSB
* JohnsonSU
* Laplace
* Logistic
* Lognormal
* NakaRush (Naka-Rushton)
* Normal
* Pareto
* Quantal
* Quick
* r
* Rayleigh
* Recinormal (1/normal)
* RNGamma
* rNoncentral
* Rosin
* SkewNor (skew normal)
* t
* tNoncentral
* Triangular
* Uniform
* VonMises
* Wald
* Weibull


[3]  Derived Distributions
---------------------------

Here are most of the types of derived distributions that have been implemented so far:
(RVs are the distributions from which these are derived).

* AddTrans (add a constant to 1 RV)
* ArcsinTrans (arcsin transformation of 1 RV)
* Convolution (sum of 2+ RVs)
* Difference (difference of 2 RVs)
* ExpTrans (exponential transformation of 1 RV)
* InfMix (infinite mixture; i.e., a parameter of one distribution is not fixed but instead varies according to some other distribution)
* InverseTrans (inverse transformation of 1 RV)
* LinearTrans (linear transformation of 1 RV)
* LogTrans (log transformation of 1 RV)
* MinBound (distribution implied by Frechet's bound)
* Mixture (probability mixture of 2+ RVs)
* MonotoneTrans (any monotonic transformation of 1 RV)
* MultTrans (multiplicative transformation of 1 RV)
* Order (k'th order statistic of any 2+ RVs)
* OrderIID (k'th order statistic of K IID RVs)
* PowerTrans (power transformation of 1 RV)
* Product (product of 2 RVs)
* Ratio (ratio of 2 RVs)
* TruncatedX (truncated version of 1 RV)

[4]  Functions and Operations Available for Any Distribution
------------------------------------------------------------

Here are most of the values that can be computed for any distribution, standard or derived:

* Probability density and cumulative distribution functions
* Spline approximations of PDF and CDF, sometimes useful for speed
* Inverse cumulative probability
* Mean, standard deviation, variance, skewness, kurtosis, coefficient of variation
* Generate random numbers
* Log likelihood of a set of observations
* Moments: Raw and central; unconditional and conditional
* Hazard function
* Moment generating function
* Maximum likelihood parameter estimates & standard errors based on Fisher information
* Parameter estimates based on method of moments
* Parameter estimates based on minimum chi-square for bin frequencies
* Parameter estimates based on closest match to estimated percentiles
* Parameter estimates to produce a desired probability within a certain range
* Parameter estimates for probit models with yes/no and m-AFC tasks

[5]  Additional requirements
------------------------------------------------------------

The special-purpose multivariate random number generator class RandGen also requires the
ExtractNameVal package from https://github.com/milleratotago/ExtractNameVal.

Copyright (C) 2018 Jeffrey Owen Miller
  
      This program is free software: you can redistribute it and/or modify
      it under the terms of the GNU General Public License as published by
      the Free Software Foundation, either version 3 of the License, or
      any later version.
  
      This program is distributed in the hope that it will be useful,
      but WITHOUT ANY WARRANTY; without even the implied warranty of
      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      GNU General Public License for more details.
  
      You should have received a copy of the GNU General Public License
      along with this program (00License.txt). If not, see 
      <http://www.gnu.org/licenses/>.
 
