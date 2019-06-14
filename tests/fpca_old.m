function [pca,pc] = fpca(Tfdobj,Sfdobj)
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

basis  = getbasis(Tfdobj); %same basis for T and S
Tcoef  = getcoef(Tfdobj);
Scoef  = getcoef(Sfdobj);
metric = eval_penalty(basis);

mybn   = getnbasis(basis);
nobs   = size(Tcoef,2);
prange = getbasisrange(basis);
depth  = (prange(1):prange(2))';

pca = [];
pca.basis = basis;
pca.metric = metric;
pca.nobs = nobs;

C = [Tcoef' Scoef'];
nobs = size(C,1);
Cm = mean(C,1);
Cc = C - repmat(Cm,nobs,1);
pca.Cm = Cm;

% inertia
VT = 1/nobs*Cc(:,1:mybn)'*Cc(:,1:mybn)*metric;
inerT = trace(VT);
pca.inerT = inerT;
VS = 1/nobs*Cc(:,mybn+1:end)'*Cc(:,mybn+1:end)*metric;
inerS = trace(VS);
pca.inerS = inerS;

% metric M
Mdem = diag([ones(mybn,1)/sqrt(inerT);ones(mybn,1)/sqrt(inerS)]);
Mdeminv = diag([ones(mybn,1)*sqrt(inerT);ones(mybn,1)*sqrt(inerS)]);
M = diag([ones(mybn,1)/inerT;ones(mybn,1)/inerS]);

% metric W
W = blkdiag(metric,metric);
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
pca.axes = pca.vectors .* repmat(sqrt(pca.values),1,2*mybn);

%% comments
% Verif0 = pca.vectors(:,1)' * W * M * pca.vectors(:,1) - 1;
% Verif1 = pca.vectors(:,1)' * W * M * pca.vectors(:,2);
% Verif2 = pca.values(1) - 1/nobs*pca.pc(:,1)'*pca.pc(:,1);
% 
% plot(pca.pc(:,1),pca.pc(:,2),'.')
% plot3(pca.pc(:,1),pca.pc(:,2),pca.pc(:,3),'.')

