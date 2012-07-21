% Make some plots for MCQMC 2012 paper
%% Garbage collection
set(0,'defaultaxesfontsize',20,'defaulttextfontsize',20)
clear all, close all
format compact

%% Set fixed constants
alpha=0.01;
talpha=1-sqrt(1-alpha);
beta=0.01;
A=0.56;

%% Set eps/sigma and kappa range
tolovsigvec=-norminv(talpha/2)./sqrt(10.^(4:0.05:9)');
ntol=length(tolovsigvec);
kappavec=[10 100 1000];
nkappa=length(kappavec);

%% Define functions
NCLT = @(tolovsig) ceil((norminv(alpha/2)./tolovsig).^2); %CLT sample
NCheb = @(tolovsig) ceil(1./((tolovsig.^2).*talpha)); %Chebychev sample size
NBEnonfunzero=@(logsqrtn,tolovsig,rho) normcdf(-exp(logsqrtn).*tolovsig) + ...
    A*rho.*exp(-logsqrtn)./(1+exp(logsqrtn).*tolovsig).^3 - talpha/2;
    %solve for non-uniform Berry-Esseen sample size
NBE = @(tolovsig,rho) ...
    ceil(exp(2*fzero(@(x) NBEnonfunzero(x,tolovsig,rho),log(sqrt(NCLT(tolovsig))))));
NChebBE = @(tolovsig,rho) min(NCheb(tolovsig),NBE(tolovsig,rho));
v2weight = @(alpha,beta,fudge2) ...
    (fudge2)+(fudge2-1).*sqrt(talpha.*(1-beta)/(beta*(1-alpha)));
fudge2fun = @(kappa,nsigma) 1./(1 - sqrt((kappa-(nsigma-3)./(nsigma-1)).* ...
    ((1-alpha)./(alpha*nsigma))));
Nmubound = @(tolovsig,kappa,nsigma) ...
    max(nsigma,NChebBE(tolovsig/v2weight(talpha,beta,fudge2fun(kappa,nsigma)),kappa.^(3/4)));
Ntot = @(tolovsig,kappa,nsigma) nsigma+Nmubound(tolovsig,kappa,nsigma);

%% Initialize sample sizes
NCLTvec=NCLT(tolovsigvec); %CLT sample size
NChebBEvec=zeros(ntol,nkappa); %Berry-Esseen sample size
Ntotvec0=NChebBEvec; %Upper bound on total sample ssize
nsigoptvec=NChebBEvec; %Optimal nsigma
Ntotoptvec=NChebBEvec;
fudgeoptvec=NChebBEvec;
nsigma0vec=zeros(1,nkappa);
fudge0vec=nsigma0vec;

%% Main loop
for k=1:nkappa
    
    tic
    %% Set kurtosis
    kappa=kappavec(k);
    rho=kappa^(3/4);

    %% Plot some stuff to look
%     nsigmavec=10.^(3:0.1:8)';
%     nsig=length(nsigmavec);
%     tolovsig=1e-2;
%     Ntotvec=zeros(size(nsigmavec));
%     for i=1:nsig
%         Ntotvec(i)=Ntot(tolovsig,kappa,nsigmavec(i));
%     end
%     figure
%     loglog(nsigmavec,Ntotvec,'b-','linewidth',2)

    %% Compute sample sizes
    nsigma0=kappa*1e3;
    nsigma0vec(k)=nsigma0;
    fudge0=sqrt(fudge2fun(kappa,nsigma0));
    fudge0vec=fudge0;
    for i=1:ntol
        tolovsig=tolovsigvec(i);
        NChebBEvec(i,k)=NChebBE(tolovsig,rho);
        Ntotvec0(i,k)=Ntot(tolovsig,kappa,nsigma0);
        minlog10n=max(3,log10(kappa*(1-talpha)/talpha));
        nsigopt=10.^(fminbnd(@(x) Ntot(tolovsig,kappa,10.^x),minlog10n,10));
        nsigoptvec(i,k)=round(nsigopt);
        Ntotoptvec(i,k)=Ntot(tolovsig,kappa,nsigoptvec(i,k));
    end
    fudgeoptvec(:,k)=sqrt(fudge2fun(kappa,nsigoptvec(:,k)));

    toc
end

%% Plot ratios of sample sizes
figure;
h=loglog( ... %NCLTvec,NChebBEvec./repmat(NCLTvec,1,nkappa),'k-.',...
    NCLTvec,Ntotvec0./repmat(NCLTvec,1,nkappa),'k--',...
    NCLTvec,Ntotoptvec./repmat(NCLTvec,1,nkappa),'k-',...
    'linewidth',2);
set(gca,'Linewidth',2);
xlabel('{\it N}_{CLT}')
ylabel('Cost Ratios')
%ylabel('{\it N_{\mu}}/{\it N}_{CLT}')
axis([1e4 1e9 1 100])
print -deps 'MCSampleSizes.eps'

%% Plot optima nsigma and fudge
figure;
line(NCLTvec,fudgeoptvec,'color','k','linestyle','--','linewidth',2)
ax1 = gca;
set(ax1,'XColor','k','YColor','k',...
    'XLim',[1e4 1e9],...
    'Xscale','log','Xtick',10.^(4:9),...
    'Yscale','linear','Ytick',1:0.1:1.7, ...
    'YLim',[1 1.7],'Linewidth',2);
xlabel('{\it N}_{CLT}')
ylabel('{\it C}')
axes('Position',get(ax1,'Position'),...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none',...
    'XColor','k','YColor','k',...
    'XLim',[1e4 1e9],...
    'YLim',[1e4 1e9],...
    'LineWidth',2, ...
    'XScale','log','Yscale','log', ...
    'XTick',[],'YTick',10.^(4:9));
ylabel('{\it n_{\sigma}}')
line(NCLTvec,nsigoptvec,'color','k','linestyle','-','linewidth',2)
print -deps 'MCnsigmafudge.eps'



