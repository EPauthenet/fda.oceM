function pc = proj3(pca,fdobj,scalars)
% Projection of Temperature and Salinity profiles
%
% Projection of any Temperature and Salinity (T-S) profiles on a chosen basis. It is essential to locate a profile relatively to a climatology. The profiles to project must have the same \code{myb} than the modes to project on (i.e. same length of profile \code{c(dmin,dmax)} and same number of Bsplines \code{mybn}).
%
% @param temp.fd,sal.fd fd objects (list) of the T-S profile(s) to project. This is produced by the function \code{bspl}.
% @param pca list of the modes to project on, produced by the function \code{fpca}.
% @param te how many PCs to use in the reconstruction, default is set to the total number of PC, \code{te = mybn}.
% @return \code{Npc} The principal components of the profiles \code{temp.fd} and \code{sal.fd} projected on the modes contained in \code{pca}.
%
% @references Ramsay J. O., and B. W. Silverman, 2005: Functional Data Analysis. Springer, 426 pp.
%
% @seealso \code{\link{bspl}} for bsplines fit on T-S profiles, \code{\link{fpca}} FPCA of T-S profiles, \code{\link{PCmap}} for plotting a map of PC, \code{\link{kde_pc}} for kernel density estimation of two PCs...

basis = pca.basis;
nbas = getnbasis(basis);
ndim = pca.ndim;
nsca = pca.nsca; % must equal size(scalars,2)

Cm = pca.Cm;
coef  = getcoef(fdobj);
nobs = size(coef,2); % must equal size(scalars,1)

C = zeros(nobs,ndim*nbas+nsca);
for kk=1:ndim,
    C(:,(kk-1)*nbas+1:kk*nbas) = squeeze(coef(:,:,kk))';
end
if nsca,
	C(:,nbas*ndim+1:nbas*ndim+nsca) = scalars';
end
Cc = C - repmat(Cm,nobs,1);

pc = Cc * pca.W' * pca.M * pca.vectors;

