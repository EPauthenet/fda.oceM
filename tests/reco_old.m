function fd_reco = reco_fd(pca,pc,Ntrunc)
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
mybn = getnbasis(basis);
nobs = size(pc,1);

if ~exist('Ntrunc','var'), Ntrunc = 2*mybn; end

Tcoef = repmat(pca.Cm(1:mybn)',1,nobs) + pca.vectors(1:mybn,1:Ntrunc)*pc(:,1:Ntrunc)';
Scoef = repmat(pca.Cm(mybn+1:2*mybn)',1,nobs) + pca.vectors(mybn+1:2*mybn,1:Ntrunc)*pc(:,1:Ntrunc)';
fd_reco.T = fd(Tcoef,basis);
fd_reco.S = fd(Scoef,basis);

