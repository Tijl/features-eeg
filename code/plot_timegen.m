%% add paths
addpath('~/Repository/CommonFunctions/matplotlib/')
addpath('~/CoSMoMVPA/mvpa/')

%% timegens

load('../results/stats_timegen.mat')

cmap = cividis();

timv = stats.timevect;
xlab = 0:100:800;

ts = [0.24 0.26; .22 .28; .22 .28; .22 .28]; % scales
durations = {'soa150','soa50'};h=[];
durationsplot = {'6.67Hz','20Hz'};
targetlabels = {'ori','sf','color','contrast'};
targetlabelsplot = {'Orientation','SF','Colour','Contrast'};

f=figure(1);clf
f.Position(3:4)=[500 800];
plotnr=0;
for t = 1:length(targetlabels)
    for d=1:2
        plotnr=plotnr+1;
        a=subplot(4,2,plotnr);hold on

        dat = stats.(durations{d}).(targetlabels{t});
        mu = reshape(dat.mu,length(timv),length(timv));
        bf = reshape(dat.bf,length(timv),length(timv));
        mu(bf<10) = .25;

        %%% plot time gen decoding
        ax=imagesc(timv,timv,(mu),ts(t,:));
        colormap(gca,cmap)
        set(gca,'YDir','normal')
        cb=colorbar();
        cb.Ticks = [ts(t,1) 0.25 ts(t,2)];
        xlim(minmax(timv'))
        ylim(minmax(timv'))
        set(gca,'XTick',xlab,'XTickLabel',num2cell(xlab));
        set(gca,'YTick',xlab,'YTickLabel',num2cell(xlab));

        %set(gca,'YTickLabel',[])

        if plotnr>6
            xlabel('Training time (ms)');
        end
        if mod(plotnr,2)
            ylabel('Testing time (ms)');
        end
        axis square
        title(sprintf('%s %s',durationsplot{d},targetlabelsplot{t}),'Color',co(t,:))
        set(gca,'Fontsize',12)

        %grid
        g = -100:100:600;
        gridcol = 'w';
        for gg=g
            plot([g*0+gg;g]',[g;g*0+gg]',':','Color',gridcol)
        end
        %diagonal
        plot(g,g,'r-','LineWidth',1)
    end
end

%% save
fn = '../figures/Figure_timegen';
print(gcf,'-dpng','-r500',fn)
im=imread([fn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');
