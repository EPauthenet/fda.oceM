cd('~/Documents/MATLAB/fda.oceM')
clear all

%% Load the data
load Argos_imos_profil.mat

plot(data.TEMP(:,3800),-data.PRES(:,3800),data.PSAL(:,3800)-35,-data.PRES(:,3800))
PRES=data.PRES(:,3800);
TEMP=data.TEMP(:,3800);
PSAL=data.PSAL(:,3800);

I=find(~isnan(PRES));
PRES=PRES(I);
TEMP=TEMP(I);
PSAL=PSAL(I);

%% shape the T and S data
prange=[0 500];
nbas = 20;
Pi = (prange(1):prange(2))';
Breaks=prange(1)+(prange(2)-prange(1))*tan(0:1/(nbas-3):1)'/tan(1);

Ti=zeros(length(Pi),data.np);
Si=zeros(length(Pi),data.np);
count = 1;
for kk=1:data.np,
    
    % T 
    P = data.PRES(:,kk);
    T = data.TEMP(:,kk);
    if max(P)<prange(2), continue; end
    I = find(~isnan(P) & ~isnan(T) & T<25 & P>prange(1) & P<prange(2));
    if length(I)<5, continue; end
    P=P(I); T=T(I);
    [P,ia,ib] = unique(P);
    T = T(ia);
    if min(P)>prange(1),
        P = [prange(1);P];
        T = [T(1);T];
    end
    if max(P)<prange(2),
        P = [P;prange(2)];
        T = [T;T(end)];
    end
    Pt = P;
    
    % S
    P = data.PRES(:,kk);
    S = data.PSAL(:,kk);
    I = find(~isnan(P) & ~isnan(S)& S<38 & S>32 & P>prange(1) & P<prange(2));
    if length(I)<5, continue; end
    P=P(I); S=S(I);
    [P,ia,ib] = unique(P);
    S = S(ia);
    if min(P)>prange(1),
        P = [prange(1);P];
        S = [S(1);S];
    end
    if max(P)<prange(2),
        P = [P;prange(2)];
        S = [S;S(end)];
    end
    Ps = P;
    
    % interp
    Ti(:,count) = interp1(Pt,T,Pi);
    Si(:,count) = interp1(Ps,S,Pi);
    count = count + 1;

end 

Ti = Ti(:,1:count-1);
Si = Si(:,1:count-1);

Xi(:,:,1) = Ti;
Xi(:,:,2) = Si;

%% Bsplines
fdn = {'Temperature','Salinity'};
nbas = 20;
fdobj = bspl(Xi,Pi,nbas,fdn);

%% fpca
pca = fpca(fdobj); % compute a pc basis !
pc = proj(fdobj,pca);

%% projection
fdn = {'Temperature','Salinity'};
nbas = 20;
fdobj2 = bspl(Xi(:,3800:3900,:),Pi,nbas,fdn);
pc2 = proj(fdobj2,pca); % project a profile on a pc basis!

%% plots
pc_plot(pca,pc)
hold on
plot(pc2(:,1),pc2(:,2),'ro')

eigenf_plot(pca,2)

%% reconstruction
fd_reco = reco(pca,pc2,3);
Xreco = eval_fd(Pi,fd_reco);
Treco = squeeze(Xreco(:,:,1));
Sreco = squeeze(Xreco(:,:,2));

figure(2),clf
plot(Ti(:,3800),-Pi,'k'); hold on
plot(eval_fd(Pi,data2fd(Pi,Ti(:,3800),basis)),-Pi,':r'); hold on
plot(Treco(:,1),-Pi,'b');
legend('raw','spline','reco3')


