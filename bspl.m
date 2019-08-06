function [fdobj] = bspl(Pi,Xi,nbas,fdn)
% B-spline fits on Multivariate Hydrographic Profiles
% 
% BSPL This function fits B-splines on multivariate hydrographic profiles and return a functional data object.
%
% ARGUMENTS
% Xi        ... Array containing the profiles stored in this order levels x stations x variables
% Pi        ... Vector containing the levels
% NBAS  ... number of Bsplines (coefficients)
%                Default is set to 20.
% FDN    ... List of the variable names
%		   Default is set to list('Temperature','Salinity').
%
% RETURN
% FDOBJ ... fd objects of the splines construction containing coefficients, basis etc... The coefficients are stored in an array of shape nbasis x stations x variables
%
% DEPENDENCIES
% The method uses the fdaM Toolbox by Jim Ramsay.
% http://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/Matlab/
% You will need to install this toolbox and add it to the matlab path to use this software
%
% CONTACT
% This code was written by Etienne Pauthenet, David Nerini and Fabien Roquet. 
% Questions, comments and bugs can be sent to: 
% etienne.pauthenet@gmail.com
% 
% REFERENCES 
% Pauthenet et al. (2017) A linear decomposition of the Southern Ocean thermohaline structure. Journal of Physical Oceanography, http://dx.doi.org/10.1175/JPO-D-16-0083.1
% Ramsay, J. O., and B. W. Silverman, 2005: Functional Data Analysis. 2nd Edition Springer, 426 pp., Isbn : 038740080X.
%
% See also function fpca for functional principal component analysis of multivariate hydrographic profiles.

if ~exist('nbas','var'), nbas = 20; end
if ~exist('fdn','var'), fdn = {'Temperature','Salinity'}; end
fdnames{1} = 'Level';
fdnames{2} = 'Station';
fdnames{3} = fdn;

prange = [Pi(1) Pi(end)];
Breaks=prange(1)+(prange(2)-prange(1))*tan(0:1/(nbas-3):1)'/tan(1);
basis = create_bspline_basis(prange,nbas,4,Breaks);

fdobj = data2fd(Pi,Xi,basis,fdnames);
end
