%%
addpath('~/Repository/CommonFunctions/Matplotlibcolors2/')
addpath('~/Repository/CommonFunctions/distributionPlot/')
addpath('~/CoSMoMVPA/mvpa/')
addpath('~/fieldtrip');ft_defaults;
load('../results/stats_decoding.mat')

%%
res_cell = {};
for s=1:16
    x=load(sprintf('../results/sub-%02i_decoding_searchlight.mat',s));
    res_cell{s} = x.res;
end
res_searchlight = cosmo_fx(cosmo_stack(res_cell),@mean,{'soaduration','targetfeature'});
res_searchlight.sa.d = 1+ (res_searchlight.sa.soaduration>0.1);

layout=cosmo_meeg_find_layout(res_searchlight);
ft = ft_timelockanalysis(cfg,cosmo_map2meeg(res_searchlight));

%%
f=figure(1);clf
f.Position = [f.Position(1:2) 800 1000];
durations = {'soa150','soa50'};h=[];
targetlabels = {'ori','sf','color','contrast'};
targetlabelsplot = {'ori','sf','colour','contrast'};
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
            'DisplayName',targetlabelsplot{t},'Color',co(t,:),...
            'LineWidth',2);
        [~,x] = max(s.mu);
    end
    legend(h)
    ylabel('accuracy')
    xlabel('time (ms)')
    xlim(minmax(s.tv))
    title(sprintf('%.2f Hz',1000/str2double(durations{d}(4:end))))

    %topo
    a1=subplot(15,2,11+(d==2));hold on
    axis off
    for t=1:4
        s = stats.(durations{d}).(targetlabels{t});
        [~,tidx] = max(s.mu);
        x = res_searchlight.samples(res_searchlight.sa.targetfeature==t & res_searchlight.sa.d==d,res_searchlight.fa.time==tidx);
        
        ft.avg(:,:) = repmat(x',1,size(ft.avg,2));
        
        mrange = prctile(x,[5 95]);

        a1pos = a1.Position;
        aw = a1pos(3);
        aw2 = aw/4;
        a=axes('Position',a1pos.*[1 1 0 0]+[t*aw2-aw2 -.35*a1pos(4) aw2 aw2]);
        hold on
        co3 = [linspace(1,co(t,1),100);...
               linspace(1,co(t,2),100);...
               linspace(1,co(t,3),100)]';
        % show figure with plots for each sensor
        cfg = [];
        cfg.zlim = mrange;
        cfg.xlim = s.tv([tidx tidx+1]);
        cfg.layout = layout;
        cfg.showscale = 'no';
        cfg.comment = 'no';
        cfg.markersymbol = '.';
        cfg.style = 'straight';
        cfg.gridscale = 128;
        ft_info off verbose
        ft_topoplotER(cfg, ft);
        a.FontSize = 12;
        a.Colormap = co3;
        set(a.Children,'LineWidth',.5)
        drawnow
    end

    %bf
    for t=1:4
        a=subplot(15,2,11 + 2*t + (d==2));hold on
        s = stats.(durations{d}).(targetlabels{t});
        plot(s.tv,1+0*s.tv,'k-');
        idx = s.bf>6;
        plot(s.tv(idx),s.bf(idx),'o','MarkerFaceColor',co(t,:),'MarkerEdgeColor',co(t,:),'Clipping','off');
        idx = s.bf<1/6;
        plot(s.tv(idx),s.bf(idx),'o','MarkerFaceColor',.5*[1 1 1],'MarkerEdgeColor',.5*[1 1 1],'Clipping','off');
        idx = s.bf>1/6 & s.bf<6;
        plot(s.tv(idx),s.bf(idx),'o','MarkerEdgeColor','k','Clipping','off');
        a.YScale='log';
        ylim(10.^(6*[-1 1]))
        xlim(minmax(s.tv))
        ylabel('BF')
        if t==4
            xlabel('time (ms)')
        else
            a.XTick = [];
        end
    end
    
    %onset
    a=subplot(3,2,5);hold on


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

        text(s.tv(max(ci)),t+.5*(d==2),sprintf(' %s %.2f Hz',targetlabelsplot{t},1000/str2double(durations{d}(4:end))),'Color',co(t,:))
        drawnow
    end
    tt = repelem(targetlabelsplot,2,1);
    %a.YTickLabel = arrayfun(@(x) sprintf('%s %.2f Hz',tt{x},1000/str2double(durations{1 + ~mod(x,2)}(4:end))), 1:8,'UniformOutput',false);
    %a.YTick=[1:.5:4.5];a.YDir='reverse';   
    a.YTick=[];a.YDir='reverse';   
    ylim([0.5 5])
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