%% Plot time course decoding of each feature, as a subset of each level of another feature

addpath('~/Dropbox (Personal)/MATLAB/matplotlib/')

scalemax = [.28 .45 .45 .4];

durations = {'soa150','soa50'};h=[];
targetlabels = {'ori','sf','color','contrast'};
targetnames = {'Orientation','SF','Colour','Contrast'};

co = tab10;co3=[];
for c = 1:4
    co3(:,:,c) = [linspace(1,co(c,1),5);...
        linspace(1,co(c,2),5);...
        linspace(1,co(c,3),5)]';
end
co3(1,:,:) = [];
co3 = flipud(co3);

%% plot by feature
load('results/stats_decoding_interactions.mat')

for d=1:2
    f=figure(d);clf
    f.Position = [f.Position(1:2) 1200 700];
    
    for fplot = 1:length(targetlabels) % decoding contrast
        for fcont = 1:length(targetlabels) % for each level of 
            
            if fplot == fcont
                continue
            end
            
            subplot(4,4,(fplot-1)*4+fcont)
            hold on
            for lev = 1:4
                
                dat = stats.(durations{d}).(targetlabels{fplot}).(targetlabels{fcont}).(sprintf('level%d',lev)).mu;
                
                plot(stats.timevect,dat,'Color',co3(lev,:,fplot),'LineWidth',1.5,'DisplayName',sprintf('%s %d',targetnames{fcont},lev));
                                
            end
            
            legend('Box','off')
%             title(sprintf('%s by %s',targetlabels{fplot},targetlabels{fcont}))
            xlim([-100 600]);
            ylim([.242 scalemax(fplot)])
            
        end
    end
end

%% plot mean decoding on diagonal
load('results/stats_decoding.mat')

for d=1:2
    f=figure(d);
    
    for fplot = 1:length(targetlabels)
        subplot(4,4,(fplot-1)*4+fplot)
        
        dat = stats.(durations{d}).(targetlabels{fplot}).mu;
        
        plot(stats.timevect,dat,'Color',co(fplot,:),'LineWidth',3);
        
%         title(sprintf('%s',targetlabels{fplot}))
        xlim([-100 600]);
        ylim([.242 scalemax(fplot)])
        box off
        
    end
    
    %% labels
    
    for l = 1:4
        annotation('textbox',[.03 0.82-.22*(l-1) .05 .05],'String',targetnames{l},...
            'FontSize',24,'Color',co(l,:),...
            'HorizontalAlignment','center','LineStyle','none');
    end
    
    %%
    fn = sprintf('figures/Supp_interactions_%s',durations{d});
    print(gcf,'-dpng','-r500',fn)
    im=imread([fn '.png']);
    [i,j]=find(mean(im,3)<255);margin=2;
    imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');
    
end

