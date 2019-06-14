%% test fda

cd('~/Documents/MATLAB/fda.oceM')
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


%% create bspline basis (and test the breaks change)
%The increase of number of breaks near the surface improve the good fit of
%the data.
prange=[0 500];
I=find(PRES>prange(1)&PRES<prange(2));
t1=TEMP(1);
t2 = interp1(PRES,TEMP,prange(2));
P = [prange(1);PRES(I);prange(2)];
T = [t1;TEMP(I);t2];

Pi = (prange(1):prange(2))';
Ti = interp1(P,T,Pi);


mybn=50;
Breaks=prange(1)+(prange(2)-prange(1))*tan(0:1/(mybn-3):1)'/tan(1); %Function for the breaks spacing.

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


%% shape the T and S data
prange=[0 500];
mybn=65;
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

Tfdobj = data2fd(Pi,Ti,basis);
Sfdobj = data2fd(Pi,Si,basis);

%% fpca
pca = fpca_old(Tfdobj,Sfdobj); % compute a pc basis !
pc = proj_old(Tfdobj,Sfdobj,pca);

%% projection
Tfdobj2 = data2fd(Pi,Ti(:,3800:3900),pca.basis);
Sfdobj2 = data2fd(Pi,Si(:,3800:3900),pca.basis);
pc2 = proj_old(Tfdobj2,Sfdobj2,pca); % project a profile on a pc basis!
figure(1),clf
plot(pc(:,1),pc(:,2),'.',pc2(:,1),pc2(:,2),'ro')

%% reconstruction
figure(2),clf
plot(Ti(:,3800),-Pi,'k'); hold on
plot(eval_fd(Pi,data2fd(Pi,Ti(:,3800),basis)),-Pi,':r'); hold on
fd_reco = reco_old(pca,pc2,3);
Treco = eval_fd(Pi,fd_reco.T);
plot(Treco(:,1),-Pi,'b');
legend('raw','spline','reco3')

%% fpca2
Xi(:,:,1) = Ti;
Xi(:,:,2) = Si;
fdnames{1} = 'Pressure';
fdnames{2} = 'Station';
fdnames{3} = {'Temperature','Salinity'};
fdobj = data2fd(Pi,Xi,basis,fdnames);
pca = fpca(fdobj); % compute a pc basis !
pc = proj(fdobj,pca);
Npc=1; varT_percent_PC1 = 100*sum(pca.vecnotWM(1:pca.nbas,Npc).^2);
Npc=1; varS_percent_PC1 = 100*sum(pca.vecnotWM(pca.nbas+1:end,Npc).^2);
Npc=2; varT_percent_PC2 = 100*sum(pca.vecnotWM(1:pca.nbas,Npc).^2);
Npc=2; varS_percent_PC2 = 100*sum(pca.vecnotWM(pca.nbas+1:end,Npc).^2);

%% projection2
fdobj2 = data2fd(Pi,Xi(:,3800:3900,:),basis,fdnames);
pc2 = proj(fdobj2,pca); % project a profile on a pc basis!
figure(2),clf
plot(pc(:,1),pc(:,2),'.',pc2(:,1),pc2(:,2),'ro')

%% reconstruction2
figure(2),clf
plot(Ti(:,3800),-Pi,'k'); hold on
plot(eval_fd(Pi,data2fd(Pi,Ti(:,3800),basis)),-Pi,':r'); hold on
fd_reco = reco(pca,pc2,3);
Xreco = eval_fd(Pi,fd_reco);
Treco = squeeze(Xreco(:,:,1));
Sreco = squeeze(Xreco(:,:,2));
plot(Treco(:,1),-Pi,'b');
legend('raw','spline','reco3')






