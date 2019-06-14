function pc = proj(fdobj,pca)
% Projection of multivariate hydrographic profiles on a basis of
% decomposition
%
% PROJ This function projects multivariate hydrographic profiles on a basis
% defined with the function fpca. It returns the Principal Components (PCs)
% of the profiles.
%
% ARGUMENTS
% FDOBJ ... functional object containing the coefficients of the spline in an array
% in that order nbas * nobs * nvar
% PCA ... structure containing the vertical modes to project on
% 
% RETURN
% PC ... the principal components in a matrix of shape nobs x (nbas*ndim)
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
basis = pca.basis;
nbas = getnbasis(basis);
ndim = pca.ndim;
Cm = pca.Cm;
coef  = getcoef(fdobj);
nobs = size(coef,2);

C = zeros(nobs,ndim*nbas);
for kk=1:ndim
    C(:,(kk-1)*nbas+1:kk*nbas) = squeeze(coef(:,:,kk))';
end
Cc = C - repmat(Cm,nobs,1);

pc = Cc * pca.W' * pca.M * pca.vectors;

