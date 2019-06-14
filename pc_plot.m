function pc_plot(pca,pc,n)
% Graphical representation of the PC on the main principal plans
%
% PC_PLOT This function plot the PC values with the appropriate axis and percentage values on the labels.
%
% ARGUMENTS
% PCA ... list containing the vertical modes used (computed with the function \code{fpca}).
% PC  ... The principal components of the profiles (computed with the function \code{proj}).
% N   ... vector of the PC to plot, it can be of length 2 or 3. Default is n = c(1,2)
%
% RETURN
% Plot the main PC with their corresponding variance on the axis labels.
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

 if ~exist('n','var'), n = [1 2]; end
  if length(n)==2 
    figure(1),clf
    plot(pc(:,n(1)),pc(:,n(2)),'.k');
    xlabel(['PC',num2str(n(1)),' (',num2str(round(pca.pval(n(1)),2)),'%)'])
    ylabel(['PC',num2str(n(2)),' (',num2str(round(pca.pval(n(2)),2)),'%)'])
    yline(0);
    xline(0);
  end
  
  if length(n)==3
    figure(1),clf
    plot3(pc(:,n(1)),pc(:,n(2)),pc(:,n(3)),'.k');
    xlabel(['PC',num2str(n(1)),' (',num2str(round(pca.pval(n(1)),2)),'%)'])
    ylabel(['PC',num2str(n(2)),' (',num2str(round(pca.pval(n(2)),2)),'%)'])
    zlabel(['PC',num2str(n(3)),' (',num2str(round(pca.pval(n(3)),2)),'%)'])
  end


