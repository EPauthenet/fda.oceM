%% test fda

cd('/Users/roquet/Mirounga Dropbox/MEOP/FDA_test')
addpath('./fdaM/')
clear all

%%
load Argos_imos_profil.mat

plot(data.TEMP(:,3800),-data.PRES(:,3800),data.PSAL(:,3800)-35,-data.PRES(:,3800))
PRES=data.PRES(:,3800);
TEMP=data.TEMP(:,3800);
PSAL=data.PSAL(:,3800);

I=find(~isnan(PRES));
PRES=PRES(I);
TEMP=TEMP(I);
PSAL=PSAL(I);


%% create bspline basis
prange=[0 500];
I=find(PRES>prange(1)&PRES<prange(2));
t1=TEMP(1);
t2 = interp1(PRES,TEMP,prange(2));
P = [prange(1);PRES(I);prange(2)];
T = [t1;TEMP(I);t2];

Pi = (prange(1):prange(2))';
Ti = interp1(P,T,Pi);


mybn=30;
Breaks=prange(1)+(prange(2)-prange(1))*tan(0:1/(mybn-3):1)'/tan(1);

basis = create_bspline_basis(prange,mybn,4,Breaks);
fdobj = data2fd(Pi,Ti,basis);
Ti1 = eval_fd(Pi,fdobj);

basis2 = create_bspline_basis(prange,mybn,4);
fdobj2 = data2fd(Pi,Ti,basis2);
Ti2 = eval_fd(Pi,fdobj2);

figure(1),clf
subplot(1,2,1)
plot(Ti1-Ti,-Pi,Ti2-Ti,-Pi)
legend(sprintf('sp1: rms=%g',rms(Ti1-Ti)),sprintf('sp2: rms=%g',rms(Ti2-Ti)),'location','SE')
subplot(1,2,2)
plot(Ti1,-Pi,Ti2,-Pi,Ti,-Pi)
legend('sp1','sp2','ref','location','SE')

figure(4),clf,plot(basis)
figure(5),clf,plot(basis2)


%% compute coeffs
prange=[0 500];
mybn=30;
Breaks=prange(1)+(prange(2)-prange(1))*tan(0:1/(mybn-3):1)'/tan(1);
basis = create_bspline_basis(prange,mybn,4,Breaks);
metric = eval_penalty(basis);
Pi = (prange(1):prange(2))';

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

sst = Ti(1,:) + randn(1,size(Ti,2));

%% old method
% fpca
Tfdobj = data2fd(Pi,Ti,basis);
Sfdobj = data2fd(Pi,Si,basis);
pca = fpca(Tfdobj,Sfdobj); % compute a pc basis !
pc = proj(Tfdobj,Sfdobj,pca);

% projection
Tfdobj2 = data2fd(Pi,Ti(:,3800:3900),pca.basis);
Sfdobj2 = data2fd(Pi,Si(:,3800:3900),pca.basis);
pc2 = proj(Tfdobj2,Sfdobj2,pca); % project a profile on a pc basis!
figure(1),clf
plot(pc(:,1),pc(:,2),'.',pc2(:,1),pc2(:,2),'ro')

% reconstruction
figure(2),clf
plot(Ti(:,3800),-Pi,'k'); hold on
plot(eval_fd(Pi,data2fd(Pi,Ti(:,3800),basis)),-Pi,':r'); hold on
fd_reco = reco_fd(pca,pc2,3);
Treco = eval_fd(Pi,fd_reco.T);
plot(Treco(:,1),-Pi,'b');
legend('raw','spline','reco3')


%% new method
% fpca2
Xi(:,:,1) = Ti;
Xi(:,:,2) = Si;
fdnames{1} = 'Pressure';
fdnames{2} = 'Station';
fdnames{3} = {'Temperature','Salinity'};
fdobj = data2fd(Pi,Xi,basis,fdnames);
pca = fpca2(fdobj); % compute a pc basis !
pc = proj2(fdobj,pca);
pc_old = pc;
Npc=1; varT_percent_PC1 = 100*sum(pca.vecnotWM(1:pca.nbas,Npc).^2);
Npc=1; varS_percent_PC1 = 100*sum(pca.vecnotWM(pca.nbas+1:end,Npc).^2);
Npc=2; varT_percent_PC2 = 100*sum(pca.vecnotWM(1:pca.nbas,Npc).^2);
Npc=2; varS_percent_PC2 = 100*sum(pca.vecnotWM(pca.nbas+1:end,Npc).^2);

% projection2
fdobj2 = data2fd(Pi,Xi(:,3800:3900,:),basis,fdnames);
pc2 = proj2(fdobj2,pca); % project a profile on a pc basis!
figure(2),clf
plot(pc(:,1),pc(:,2),'.',pc2(:,1),pc2(:,2),'ro')
figure(4),clf
plot3(pc(:,1),pc(:,2),pc(:,3),'.',pc2(:,1),pc2(:,2),pc2(:,3),'ro')

% reconstruction2
figure(2),clf
plot(Ti(:,3800),-Pi,'k'); hold on
plot(eval_fd(Pi,data2fd(Pi,Ti(:,3800),basis)),-Pi,':r'); hold on
fd_reco = reco_fd2(pca,pc2,3);
Xreco = eval_fd(Pi,fd_reco);
Treco = squeeze(Xreco(:,:,1));
Sreco = squeeze(Xreco(:,:,2));
plot(Treco(:,1),-Pi,'b');
legend('raw','spline','reco3')


%% new method + sst scalar
% fpca2
Xi(:,:,1) = Ti;
Xi(:,:,2) = Si;
fdnames{1} = 'Pressure';
fdnames{2} = 'Station';
fdnames{3} = {'Temperature','Salinity'};
fdobj = data2fd(Pi,Xi,basis,fdnames);
pca = fpca3(fdobj,sst,[1/2 1 1/2]); % compute a pc basis !
pc = proj3(pca,fdobj,sst);
pc_new = pc;
pca_new=pca;

Npc=1; varT_percent_PC1 = 100*sum(pca.vecnotWM(1:pca.nbas,Npc).^2)
Npc=1; varS_percent_PC1 = 100*sum(pca.vecnotWM(pca.nbas+1:end-1,Npc).^2)
Npc=1; varScal_percent_PC1 = 100*sum(pca.vecnotWM(end,Npc).^2)
Npc=2; varT_percent_PC2 = 100*sum(pca.vecnotWM(1:pca.nbas,Npc).^2);
Npc=2; varS_percent_PC2 = 100*sum(pca.vecnotWM(pca.nbas+1:end,Npc).^2);

% projection2
fdobj2 = data2fd(Pi,Xi(:,3800:3900,:),basis,fdnames);
sst2 = sst(3800:3900);
pc2 = proj3(pca,fdobj2,sst2); % project a profile on a pc basis!
figure(2),clf
plot(pc(:,1),pc(:,2),'.',pc2(:,1),pc2(:,2),'ro')
figure(3),clf
plot3(pc(:,1),pc(:,2),pc(:,3),'.',pc2(:,1),pc2(:,2),pc2(:,3),'ro')

% reconstruction2
[fd_reco,sst_reco] = reco_fd3(pca,pc2,1);
[fd_reco2,sst_reco2] = reco_fd3(pca,pc2,2);

figure(2),clf
plot(Ti(:,3800),-Pi,'k'); hold on
plot(eval_fd(Pi,data2fd(Pi,Ti(:,3800),basis)),-Pi,':r'); hold on
Xreco = eval_fd(Pi,fd_reco);
Treco = squeeze(Xreco(:,:,1));
Sreco = squeeze(Xreco(:,:,2));
plot(Treco(:,1),-Pi,'b');
legend('raw','spline','reco3')

figure(3),clf
plot(1:length(sst2),sst2,1:length(sst2),sst_reco,1:length(sst2),sst_reco2)


%% canonical analysis: pc_500 vs pc_200
prange=[0 300];
mybn=30;
Breaks=prange(1)+(prange(2)-prange(1))*tan(0:1/(mybn-3):1)'/tan(1);
basis2 = create_bspline_basis(prange,mybn,4,Breaks);
Pi2 = (prange(1):prange(2))';
Xi2 = [];
Xi2(:,:,1) = Ti(1:301,:);
Xi2(:,:,2) = Si(1:301,:);
fdobj = data2fd(Pi2,Xi2,basis2,fdnames);
pca = fpca3(fdobj,sst); % compute a pc basis !
pc = proj3(pca,fdobj,sst);
pc_200 = pc;


X=pc_200;
Y=pc_new;
N=pca.nobs;
[A,B,r,U,V,stats] = canoncorr(pc_200,pc_new);
U = (X-repmat(mean(X),N,1))*A;
V = (Y-repmat(mean(Y),N,1))*B;
Y_reco = repmat(mean(Y),N,1) + (X-repmat(mean(X),N,1)) * A * inv(B);

figure(1),clf,plot(Y(:,1),Y_reco(:,1),'.')
figure(2),clf,plot(Y(:,2),Y_reco(:,2),'.')

[fd_reco,sst_reco] = reco_fd3(pca_new,Y,40);
[fd_reco2,sst_reco] = reco_fd3(pca_new,Y_reco,40);

Xreco = eval_fd(Pi,fd_reco);
Treco = squeeze(Xreco(:,:,1));
Sreco = squeeze(Xreco(:,:,2));
Xreco2 = eval_fd(Pi,fd_reco2);
Treco2 = squeeze(Xreco2(:,:,1));
Sreco2 = squeeze(Xreco2(:,:,2));

figure(3),clf
plot(Ti(:,3800),-Pi,'k'); hold on
plot(eval_fd(Pi,data2fd(Pi,Ti(:,3800),basis)),-Pi,':r'); hold on
plot(Treco(:,3800),-Pi,Treco2(:,3800),-Pi,'b');
legend('raw','spline','reco','reco2')




