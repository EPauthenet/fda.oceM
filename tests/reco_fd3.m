function [fd_reco,scalars_reco] = reco_fd3(pca,pc,Ntrunc)
% Reconstruction of Temperature and Salinity profiles
%
% Reconstruction of Temperature and Salinity (T-S) profiles with a chosen number of Principal Components (PCs).
%
% @param pca list produced by the function \code{fpca}
% @param te how many PCs to use in the reconstruction, default is set to the total number of PC, \code{te = mybn}.
% @param depth how many depth levels to get back.
% @return \code{recotemp} and \code{recosali} matrix containing the reconstruction of \code{temp} and \code{sali} with the number of PCs \code{te}.
%
% @references Ramsay J. O., and B. W. Silverman, 2005: Functional Data Analysis. Springer, 426 pp.
%
% @seealso \code{\link{bspl}} for bsplines fit on T-S profiles, \code{\link{fpca}} FPCA of T-S profiles, \code{\link{PCmap}} for plotting a map of PC, \code{\link{kde_pc}} for kernel density estimation of two PCs...

basis = pca.basis;
prange = getbasisrange(basis);
nbas = pca.nbas;
ndim = pca.ndim;
nsca = pca.nsca;

if ~exist('Ntrunc','var'), Ntrunc = ndim*nbas+nsca; end

nobs = size(pc,1);
coef = zeros(nbas,nobs,ndim);
for kk=1:ndim,
    coef(:,:,kk) = repmat(pca.Cm((kk-1)*nbas+1:kk*nbas)',1,nobs) + pca.vectors((kk-1)*nbas+1:kk*nbas,1:Ntrunc)*pc(:,1:Ntrunc)';
end

fd_reco = fd(coef,basis,pca.fdnames);

scalars_reco = [];
if nsca,
    scalars_reco = zeros(nobs,nsca);
    for kk=1:nsca,
        scalars_reco(:,kk) = pca.Cm(nbas*ndim+nsca) + pca.vectors(nbas*ndim+nsca,1:Ntrunc)*pc(:,1:Ntrunc)';
    end
end

