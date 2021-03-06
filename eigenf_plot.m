function eigenf_plot(pca,te,sign)
% Plot of the vertical modes
%
% EIGENF_PLOT Plot the effect of each vertical mode (i. e. eigenfunction) on the mean profile.
%
% ARGUMENTS
% PCA  ... list containing the vertical modes used (computed with the function \code{fpca}).
% TE   ... Which PC to plot.
% N    ... vector of the PC to plot, it can be of length 2 or 3. Default is n = c(1,2)
% SIGN ... The PCs are invariant by their sign, so the choice of the sign depends on the feature to represent. The parameter \code{sign} can contain a vector of \code{1} or \code{-1} to inverse the sign of the PCs if wanted. \code{sign} can also be used as a factor to increase the effect of eigenfunctions and see better the small variations.
%
% RETURN
% Plot of the effects of the vertical mode \code{te} on the mean profiles. The curves show the mean profile (solid) and the effects of adding (red) and subtracting (blue) eigenfunctions. The percentages next to the header titles are the amount of variance explained by the mode displayed. The percentages in the horizontal axis label are the variance contained by each variable (T and S) on the mode displayed. \code{sign} can be used as a factor to enhance the deformation of + and - curves.
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

  if ~exist('sign','var'), sign = 1; end
  nbas  = pca.nbas;
  ndim   = pca.ndim;
  basis = pca.basis;
  Cm    = pca.Cm;
  prange = getbasisrange(basis);
  depth  = prange(1):prange(2);

  figure(1),clf
  for k = 1:ndim
    d  = ((k-1)*nbas+1):(k*nbas);                                   % Position of the variable k in the eigenvectors
    pb = round(100*sum(pca.vecnotWM(d,te).^2),0);                   % Percentage of the bloc
    fdobj_vp = fd(Cm(d)' + sign*pca.axes(d,te),pca.basis,pca.fdnames);
    fdobj_vm = fd(Cm(d)' - sign*pca.axes(d,te),pca.basis,pca.fdnames);
    fdobj_m   = fd(Cm(d)',pca.basis,pca.fdnames);
    vp = eval_fd(depth,fdobj_vp);
    vm = eval_fd(depth,fdobj_vm);
    prof_m = eval_fd(depth,fdobj_m);

    subplot(1,2,k)
    plot(prof_m, -depth,'k');
    xlabel(char([pca.fdnames{3}{k}," (",num2str(pb)," %)"]))
    ylabel(pca.fdnames(1))
    xlim([min([vp;vm;prof_m]) max([vp;vm;prof_m])])
    ylim([-max(depth) 0])
    hold on
    plot(vp, -depth, 'r')
    plot(vm, -depth, 'b')
  end
  title(['PC',num2str(te),' (',num2str(round(pca.pval(te),0)),'%)'])
end
