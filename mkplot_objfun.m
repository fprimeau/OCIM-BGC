%plot fminunc_iterations
% plot objective function as a function of parameter values (using fminunc
% saved iterations)
%mkplot_objfun

if ismac
    load('/Users/megansullivan/Documents/UC Irvine/GitHub/TraitModel_output/Tv4_PCCellv5_DOC0.25_DOP0_historysave.mat')
    xhat = fmin_history.xhat;
    figDir = '/Users/megansullivan/Documents/UC Irvine/GitHub/TraitModel_output/FIGS_v5_07Jun21/' ;
elseif isunix
    load('/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/MSK90/Tv4_PCCellv4_DOC0.25_DOP0_history.mat')
    load('/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/MSK90/Tv4_PCCellv4_DOC0.25_DOP0_xhat.mat')
    figDir = '/DFS-L/DATA/primeau/meganrs/OCIM_BGC_OUTPUT/MSK90/FIGS_PCCellv4_DOC0.25_DOP0/';
end

f = fmin_history.fval;

pnames = fieldnames(xhat);
pnames = pnames(1:end-3)
npara = numel(pnames);

pindx = struct;
for nField = 1:npara
    pField = pnames{nField} ;
    pindx.(pField) = nField ;
end
%% sigma
if isfield(xhat,'sigma')
    isigma        = pindx.sigma;
    sigma         = exp(fmin_history.x(:,isigma));
    dsigma        = fmin_history.gradient(:,isigma);
    
    figure; set(gcf,'Position',[0 0.2 0.35 0.4]);
    t = tiledlayout(1,1,'Padding','compact');
    nexttile;
    yyaxis left
        plot(sigma,f,'-^'); hold on
        plot(xhat.sigma,xhat.f,'b^','MarkerFaceColor','b','MarkerSize',8)
        ylabel('objective function value');
    yyaxis right
        plot(sigma,dsigma,'-o'); hold on
        plot(sigma,zeros(length(sigma),1),'k--')
        plot(xhat.sigma,xhat.fx(isigma),'ro','MarkerFaceColor','r','MarkerSize',8)
        ylabel('df/dlogsigma')
        xlabel('sigma')
        grid on
    title('OCIM-BGC objective function as a function of sigma')
    
    figTitle = 'objfun_sigma';
    print(gcf,[figDir 'FIG_' figTitle '.png'],'-dpng')

end
%% bC_T
if isfield(xhat,'bC_T')
    ibC_T         = pindx.bC_T;
    bC_T          = fmin_history.x(:,ibC_T);
    dbC_T         = fmin_history.gradient(:,ibC_T);
    
    figure; set(gcf,'Position',[0 0.2 0.35 0.4]);
    t = tiledlayout(1,1,'Padding','compact');
    nexttile;
    yyaxis left
        plot(bC_T,f,'-^'); hold on
        plot(xhat.bC_T,xhat.f,'b^','MarkerFaceColor','b','MarkerSize',10)
        ylabel('objective function value')
    yyaxis right
        plot(bC_T,dbC_T,'-o'); hold on
        ylabel('df/dbC_T')
        xlabel('bC_T')
        plot(bC_T,zeros(length(bC_T),1),'k--')
        plot(xhat.bC_T,xhat.fx(ibC_T),'ro','MarkerFaceColor','r','MarkerSize',10)
        grid on
        %axis tight
    title('OCIM-BGC objective function as a function of bC_T')
    
    figTitle = 'objfun_bC_T';
    print(gcf,[figDir 'FIG_' figTitle '.png'],'-dpng')
end

%% --------- Cell Model ----------
notedim = [0.15 0.65 0.1 0.2];

%% Q10Photo

if isfield(xhat,'Q10Photo')
    iQ10Photo     = pindx.Q10Photo;
    Q10Photo      = exp(fmin_history.x(:,iQ10Photo))  ;
    dQ10Photo     = fmin_history.gradient(:,iQ10Photo);
    
    figure; set(gcf,'Position',[0 0.2 0.35 0.4]);
    t = tiledlayout(1,1,'Padding','compact');
    nexttile;
    yyaxis left
        plot(Q10Photo,f,'-^'); hold on
        ylabel('objective function value')
        plot(xhat.Q10Photo,xhat.f,'b^','MarkerFaceColor','b','MarkerSize',10)
    yyaxis right
        plot(Q10Photo,dQ10Photo,'-o'); hold on
        ylabel('df/dlogQ10Photo')
        xlabel('Q10Photo')
        grid on
        plot(Q10Photo,zeros(length(Q10Photo),1),'k--')
        plot(xhat.Q10Photo,xhat.fx(iQ10Photo),'ro','MarkerFaceColor','r','MarkerSize',10)
    title('OCIM-BGC objective function as a function of Q10Photo')
    
    notes{1,1} = ['xhat obj fun value   = ' num2str(xhat.f)];
    notes{2,1} = ['xhat Q10Photo        = ' num2str(xhat.Q10Photo)]; 
    notes{3,1} = ['xhat df/dlogQ10Photo = ' num2str(xhat.fx(iQ10Photo),'%5.4e')];
    annotation('textbox',notedim,'String',notes,'FitBoxToText','on');
    
    figTitle = 'objfun_Q10Photo';
    print(gcf,[figDir 'FIG_' figTitle '.png'],'-dpng')
end

%% Q10Photo

if isfield(xhat,'Q10Photo')
    iQ10Photo     = pindx.Q10Photo;
    Q10Photo      = exp(fmin_history.x(:,iQ10Photo))  ;
    dQ10Photo     = fmin_history.gradient(:,iQ10Photo);
    
    figure; set(gcf,'Position',[0 0.2 0.35 0.4]);
    t = tiledlayout(1,1,'Padding','compact');
    nexttile;
    yyaxis left
        plot(Q10Photo,f,'-^'); hold on
        ylabel('objective function value')
        plot(xhat.Q10Photo,xhat.f,'b^','MarkerFaceColor','b','MarkerSize',10)
    yyaxis right
        plot(Q10Photo,dQ10Photo,'-o'); hold on
        ylabel('df/dlogQ10Photo')
        xlabel('Q10Photo')
        grid on
        plot(Q10Photo,zeros(length(Q10Photo),1),'k--')
        plot(xhat.Q10Photo,xhat.fx(iQ10Photo),'ro','MarkerFaceColor','r','MarkerSize',10)
    title('OCIM-BGC objective function as a function of Q10Photo')
    
    notes{1,1} = ['xhat obj fun value   = ' num2str(xhat.f)];
    notes{2,1} = ['xhat Q10Photo        = ' num2str(xhat.Q10Photo)]; 
    notes{3,1} = ['xhat df/dlogQ10Photo = ' num2str(xhat.fx(iQ10Photo),'%5.4e')];
    annotation('textbox',notedim,'String',notes,'FitBoxToText','on');
    
    figTitle = 'objfun_Q10Photo';
    %print(gcf,[figDir 'FIG_' figTitle '.png'],'-dpng')
end
%% fStorage
if isfield(xhat,'fStorage')
    ifStorage       = pindx.fStorage;
    fStorage  = exp(fmin_history.x(:,ifStorage));
    dfStorage = fmin_history.gradient(:,ifStorage);
    
    figure; set(gcf,'Position',[0 0.2 0.35 0.4]);
    t = tiledlayout(1,1,'Padding','compact');
    nexttile;
    yyaxis left
        plot(fStorage,f,'-^'); hold on
        ylabel('objective function value')
        plot(xhat.fStorage,xhat.f,'b^','MarkerFaceColor','b','MarkerSize',10)
    yyaxis right
        plot(fStorage,dfStorage,'-o'); hold on
        ylabel('df/dlog(fStorage)')
        xlabel('fStorage')
        grid on
        plot(fStorage,zeros(length(fStorage),1),'k--')
        plot(xhat.fStorage,xhat.fx(ifStorage),'ro','MarkerFaceColor','r','MarkerSize',10)
    title('OCIM-BGC objective function as a function of fStorage')
    
    notes{1,1} = ['xhat obj fun value          = ' num2str(xhat.f)];
    notes{2,1} = ['xhat fStorage               = ' num2str(xhat.fStorage,'%5.4e')]; 
    notes{3,1} = ['xhat df/dlog(fStorage)      = ' num2str(xhat.fx(ifStorage),'%5.4e')];
    annotation('textbox',notedim,'String',notes,'FitBoxToText','on');
    
    figTitle = 'objfun_fStorage';
    print(gcf,[figDir 'FIG_' figTitle '.png'],'-dpng')
end

