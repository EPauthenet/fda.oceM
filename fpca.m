function [pca] = fpca(fdobj)
% Functional Principal Component Analysis of multivariate hydrographic profiles
%
% FPCA This function computes a FPCA on multivariate hydrographic profiles,
% it returns the basis of decomposition, also called the vertical modes.
%
% ARGUMENTS
% FDOBJ ... functional object containing the coefficients of the spline in an array
% in that order nbas * nobs * nvar
%
% RETURN
% PCA ... structure containing the eigenvectors, eigenvalues,...
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
% See also function proj for computing the principal components (PCs) of a dataset.

coef  = getcoef(fdobj);
basis  = getbasis(fdobj);
metric = eval_penalty(basis);
nbas   = getnbasis(basis);
nobs   = size(coef,2);
ndim   = size(coef,3);
prange = getbasisrange(basis);
depth  = (prange(1):prange(2))';

C = zeros(nobs,ndim*nbas);
for kk=1:ndim,
    C(:,(kk-1)*nbas+1:kk*nbas) = squeeze(coef(:,:,kk))';
end

Cm = mean(C,1);
Cc = C - repmat(Cm,nobs,1);
pca.Cm = Cm;

% inertia
inertia=zeros(ndim,1);
for kk=1:ndim,
    V = Cc(:,(kk-1)*nbas+1:kk*nbas)'*Cc(:,(kk-1)*nbas+1:kk*nbas)*metric/nobs;
    inertia(kk) = trace(V);
end
pca.inertia = inertia;

% metric M
M = zeros(ndim*nbas);
Mdeminv = zeros(ndim*nbas);
for kk=1:ndim,
    M((kk-1)*nbas+1:kk*nbas,(kk-1)*nbas+1:kk*nbas)=diag(ones(nbas,1)/inertia(kk));
    Mdeminv((kk-1)*nbas+1:kk*nbas,(kk-1)*nbas+1:kk*nbas)=diag(ones(nbas,1)*sqrt(inertia(kk)));
end
Mdem = sqrt(M);

% metric W
W = [];
for kk=1:ndim,
    W = blkdiag(W,metric);
end
W = (W+W')/2;
pca.W = W;
Wdem = chol(W);
Wdeminv = inv(Wdem);

% Cov matrix and eigen values/vectors
V = 1/nobs * Mdem * Wdem * Cc' * Cc * Wdem' * Mdem;
[pca_vectors,pca_values] = eig(V,'vector');
[pca_values,I] = sort(pca_values,'descend');
pca_vectors = pca_vectors(:,I);

% complete pca structure
pca.values = pca_values;
pca.pval = 100 * pca.values / sum(pca.values);
pca.vecnotWM = pca_vectors;
pca.M = M;
pca.vectors = Mdeminv * Wdeminv * pca_vectors;
pca.axes = pca.vectors .* repmat(sqrt(pca.values),1,ndim*nbas);

pca.basis = basis;
pca.metric = metric;
pca.nbas = nbas;
pca.nobs = nobs;
pca.ndim = ndim;
pca.fdnames = getnames(fdobj);

%% Comments
% Verif0 = pca.vectors(:,1)' * W * M * pca.vectors(:,1) - 1;
% Verif1 = pca.vectors(:,1)' * W * M * pca.vectors(:,2);
% Verif2 = pca.values(1) - 1/nobs*pca.pc(:,1)'*pca.pc(:,1);

