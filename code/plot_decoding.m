%%
addpath('~/Repository/CommonFunctions/Matplotlibcolors2/')
addpath('~/Repository/CommonFunctions/distributionPlot/')
load('../results/stats_decoding.mat')

f=figure(1);clf
f.Position = [f.Position(1:2) 800 1000];
durations = {'soa150','soa50'};h=[];
targetlabels = {'ori','sf','color','contrast'};
for d=1:2
    a=subplot(3,2,d);hold on
    a.YLim=[.2 .4];
    a.FontSize=12;
    co = vega10;
    for t=1:4
        s = stats.(durations{d}).(targetlabels{t});
        fill([s.tv fliplr(s.tv)],[s.mu-s.se fliplr(s.mu+s.se)],co(t,:),...
            'FaceAlpha',.2,'LineStyle','none');
        h(t) = plot(s.tv,s.mu,...
            'DisplayName',s.targetname,'Color',co(t,:),...
            'LineWidth',2);
        [~,x] = max(s.mu);
    end
    legend(h)
    ylabel('accuracy')
    xlabel('time (ms)')
    xlim(minmax(s.tv))
    title(sprintf('%.2f Hz',1000/str2double(durations{d}(4:end))))

    %topo
    a=subplot(15,2,11+ (d==2));hold on

    %bf
    for t=1:4
        a=subplot(15,2,11 + 2*t + (d==2));hold on
        s = stats.(durations{d}).(targetlabels{t});
        plot(s.tv,s.bf,'o','MarkerFaceColor',co(t,:),'MarkerEdgeColor',co(t,:),'Clipping','off');
        a.YScale='log';
        ylim(10.^(5*[-1 1]))
        xlim(minmax(s.tv))
    end
    
    %peak
    a=subplot(3,2,6);hold on
    a.FontSize=12;
    co = vega10;
    for t=1:4
        s = stats.(durations{d}).(targetlabels{t});
        [~,peak] = max(s.mu);
        peak_boot = [];
        rng(10);
        for b=1:1000
            bidx = randsample(1:s.n,s.n,1);
            [~,peak_boot(b)] = max(mean(s.mu_all(bidx,:)));
        end
        ci = prctile(peak_boot,[5 95]);
%         plot(s.tv(peak_boot),t+linspace(.2,.3,numel(peak_boot)),'.',...
%             'MarkerSize',1,'Color',co(t,:));
%         hd=distributionPlot(s.tv(peak_boot)','Color',co(t,:),'xValues',t,...
%             'xyOri','flipped','distWidth',.95,'showMM',0);
        %plot(s.tv(peak),t,'.','Color',co(t,:),'MarkerSize',20);
        line(s.tv(round(ci))+[-.5 .5],[t t]+.5*(d==2),'Color',co(t,:),'LineWidth',10)
        drawnow
    end
    tt = repelem(targetlabels,2,1);
    a.YTickLabel = arrayfun(@(x) sprintf('%s %.2f Hz',tt{x},1000/str2double(durations{1 + ~mod(x,2)}(4:end))), 1:8,'UniformOutput',false);
    a.YTick=[1:.5:4.5];a.YDir='reverse';
    ylim([0 5])
    xlim([-100 600])
    legend(h)
    xlabel('time of peak (ms)')
    %title(sprintf('%.2f Hz',1000/str2double(durations{d}(4:end))))
    title('peak decoding (95% CI)')
end

%%
fn = '../figures/decoding';
print(gcf,'-dpng','-r500',fn)
im=imread([fn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');