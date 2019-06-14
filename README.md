# fda.oceM
Functional Data Analysis of oceanographic profiles

**Functional Data Analysis** is a set of tools to study curves or functions. Here we see vertical hydrographic profiles of several variables (temperature, salinity, oxygen,...) as curves and apply a functional principal component analysis (FPCA) in the multivaraite case to reduce the dimensionality of the system. The classical case is done with couples of temperature and salinity. It can be used for front detection, water mass identification, unsupervised or supervised classification, model comparison, data calibration ...
This repository is also available in an R package [fda.oce](https://github.com/EPauthenet/fda.oce).

*Dependencies*:
The method uses the fdaM Toolbox by Jim Ramsay.
http://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/Matlab/
You will need to install this toolbox and add it to the matlab path to use this software

*References*: 
- Pauthenet et al. (2018) Seasonal meandering of the Polar Front upstream of the Kerguelen Plateau. Geophysical Research Letters, [10.1029/2018GL079614](https://doi.org/10.1029/2018GL079614)
- Pauthenet et al. (2017) A linear decomposition of the Southern Ocean thermohaline structure. Journal of Physical Oceanography, [10.1175/JPO-D-16-0083.1](http://dx.doi.org/10.1175/JPO-D-16-0083.1)
- Ramsay, J. O., and B. W. Silverman, 2005: Functional Data Analysis. 2nd Edition Springer, 426 pp., Isbn : 038740080X.


# Demo
Here is an example of how to use these functions. We compute the modes for temperature and salinity profiles of the reanalysis [GLORYS](http://marine.copernicus.eu/services-portfolio/access-to-products/) in the Southern Ocean for December of 2015.

First we load the data and fit the Bsplines on the 1691 profiles of the example. By default the function fit 20 Bsplines. It returns a fd object named 'fdobj' :
``` Matlab
load GLORYS_SO_2015-12.mat
fdobj = bspl(Xi,Pi);
```

Then we apply the FPCA on the fd object :
``` Matlab
fpca(fdobj);
```

The profiles can be projected on the modes defined by the FPCA, to get the principal components (PCs) :
``` Matlab
proj(fdobj,pca);
```

Visualisation of the 2 first PCs :
``` Matlab
pc_plot(pca,pc);
```
<img src="https://github.com/EPauthenet/fda.oceM/blob/master/figures/pc_plot.png" alt="drawing" width="1000px"/>

Visualisation of the 2 first eigenfunctions effect on the mean profile (red (+1) and blue (-1)) :
``` Matlab
eigenf_plot(pca,1);
eigenf_plot(pca,2);
```
<img src="https://github.com/EPauthenet/fda.oceM/blob/master/figures/eigenf1.png" alt="drawing" width="370px"/> <img src="https://github.com/EPauthenet/fda.oceM/blob/master/figures/eigenf2.png" alt="drawing" width="370px"/>

The profiles can then be reconstructed with less PCs than the total number, removing the small variability. For example with only 5 modes :
``` Matlab
te = 5;
fdobj_reco = reco(pca,pc,te);
```

To transform fd objects back in a the variable space, we use the function eval.fd ("fda" package) :
``` Matlab
X = eval_fd(Pi,fdobj);
X_reco = eval_fd(Pi,fdobj_reco);
```

And finally we can represent the profiles reconstructed compared to the original data :
``` Matlab
i = 3  %index of a profile
figure(1),clf
for k = 1:ndim    %Loop for each variable                        
  subplot(1,2,k)
  plot(Xi(:,i,k),-Pi,'ko')              %Plot of the raw data
  xlim([min([Xi(:,i,k);X_reco(:,i,k)]) max([Xi(:,i,k);X_reco(:,i,k)])])
  ylim([-max(depth) 0]);
  xlabel(char([pca.fdnames{3}{k}]))
  ylabel(pca.fdnames(1))
  hold on
  plot(X(:,i,k),-Pi,'r')           %Plot of the B-spline fit
  plot(X_reco(:,i,k),-Pi,'g')     %Plot of the reconstructed profiles
end
legend('raw','spline','reconstructed')
```
<img src="https://github.com/EPauthenet/fda.oceM/blob/master/figures/reco_prof3.png" alt="drawing" width="1000px"/>

