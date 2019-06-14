function [pca,pc] = fpca3(fdobj,scalars,weights)
% Functional Principal Component Analysis of Temperature and Salinity profiles
%
% Functional Principal Component Analysis (FPCA) in the multivariate case, applied on Temperature and Salinity (T-S) profiles seen as curves (Bsplines)
%
% @param temp.fd,sal.fd fd objects (list) of the splines construction containing coefficients, etc... This is produced by the function \code{bspl}.
% @param plot,plot3d if TRUE, plot the first two or three PCs.
% @return \code{pca} list containing the eigen (\code{values}), eigen (\code{vectors}), eigen vectors not weighted by WM (\code{vecnotWM}), principal components (\code{pc}), deformation induced by an axis (\code{axe}) percentage of each axis (\code{pval}).
%
% @references Ramsay J. O., and B. W. Silverman, 2005: Functional Data Analysis. Springer, 426 pp.
%
% @seealso \code{\link{bspl}} for bsplines fit on T-S profiles, \code{\link{PCmap}} for plotting a map of PC, \code{\link{kde_pc}} for kernel density estimation of two PCs...

basis  = getbasis(fdobj);
coef  = getcoef(fdobj);
metric = eval_penalty(basis);

nbas   = getnbasis(basis);
nobs   = size(coef,2);
ndim   = size(coef,3);
prange = getbasisrange(basis);
depth  = (prange(1):prange(2))';

nsca = 0;
if exist('scalars','var') & ~isempty(scalars),
    nsca = size(scalars,1);
    if size(scalars,2)~=nobs, error('size(scalars,2) must be equal to number of observations (nobs)'); end
end

if exist('weights','var'),
    if length(weights)~=ndim+nsca, error('weights should be a vector of dimension ndim + nsca'); end
else
    weights = ones(ndim+nsca,1);
end
    
    
pca = [];
pca.basis = basis;
pca.metric = metric;
pca.nbas = nbas;
pca.nobs = nobs;
pca.ndim = ndim;
pca.fdnames = getnames(fdobj);
pca.nsca = nsca;
pca.weights = weights;

C = zeros(nobs,ndim*nbas+nsca);
for kk=1:ndim,
    C(:,(kk-1)*nbas+1:kk*nbas) = squeeze(coef(:,:,kk))';
end
if nsca,
	C(:,nbas*ndim+1:nbas*ndim+nsca) = scalars';
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

inertia_sca = [];
if nsca,
    for kk=1:nsca,
        inertia_sca(kk) = Cc(:,nbas*ndim+kk)'*Cc(:,nbas*ndim+kk)/nobs;
    end
end
pca.inertia_sca = inertia_sca;

% metric M
M = zeros(ndim*nbas+nsca);
Mdeminv = zeros(ndim*nbas+nsca);
for kk=1:ndim,
    M((kk-1)*nbas+1:kk*nbas,(kk-1)*nbas+1:kk*nbas)=diag(ones(nbas,1)*weights(kk)/inertia(kk));
    Mdeminv((kk-1)*nbas+1:kk*nbas,(kk-1)*nbas+1:kk*nbas)=diag(ones(nbas,1)*sqrt(inertia(kk)/weights(kk)));
end
if nsca,
    for kk=1:nsca,
        M(nbas*ndim+kk,nbas*ndim+kk) = weights(ndim+kk)/inertia_sca(kk);
        Mdeminv(nbas*ndim+kk,nbas*ndim+kk) = sqrt(inertia_sca(kk)/weights(ndim+kk));
    end
end
Mdem = sqrt(M);

% metric W
W = [];
for kk=1:ndim,
    W = blkdiag(W,metric);
end
if nsca,
    for kk=1:nsca,
        W(nbas*ndim+kk,nbas*ndim+kk) = 1;
    end
end
W = (W+W')/2;
pca.W = W;
Wdem = chol(W);
Wdeminv = inv(Wdem);

% Cov matrix
V = 1/nobs * Mdem * Wdem * Cc' * Cc * Wdem' * Mdem;
[pca_vectors,pca_values] = eig(V,'vector');
[pca_values,I] = sort(pca_values,'descend');
pca_vectors = pca_vectors(:,I);

% create FPCA
pca.values = pca_values;
pca.pval = 100 * pca.values / sum(pca.values);
pca.vecnotWM = pca_vectors;
pca.M = M;
pca.vectors = Mdeminv * Wdeminv * pca_vectors;
pca.axes = pca.vectors .* repmat(sqrt(pca.values),1,ndim*nbas+nsca);

%% comments
% Verif0 = pca.vectors(:,1)' * W * M * pca.vectors(:,1) - 1;
% Verif1 = pca.vectors(:,1)' * W * M * pca.vectors(:,2);
% Verif2 = pca.values(1) - 1/nobs*pca.pc(:,1)'*pca.pc(:,1);
% 
% plot(pca.pc(:,1),pca.pc(:,2),'.')
% plot3(pca.pc(:,1),pca.pc(:,2),pca.pc(:,3),'.')

