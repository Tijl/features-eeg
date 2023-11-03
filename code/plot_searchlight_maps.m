%% add paths
addpath('~/Repository/CommonFunctions/Matplotlibcolors2/')
addpath('~/CoSMoMVPA/mvpa/')
addpath('~/fieldtrip');

ft_defaults;

%% collate searchlight results

res_cell = {};
for s=1:16
    x=load(sprintf('results/sub-%02i_decoding_searchlight.mat',s));
    res_cell{s} = x.res;
end
res_searchlight = cosmo_fx(cosmo_stack(res_cell),@mean,{'soaduration','targetfeature'});
res_searchlight.sa.d = 1+ (res_searchlight.sa.soaduration>0.1);

%% get layout
layout=cosmo_meeg_find_layout(res_searchlight,'label_threshold',.80);

cfg={};
ft = ft_timelockanalysis(cfg,cosmo_map2meeg(res_searchlight));

%% plot
load('results/stats_decoding.mat')

durations = {'soa150','soa50'};h=[];
durnum = [.15 .05];
targetlabels = {'ori','sf','color','contrast'};
co = tab10;
timestoplot = [0:60:300];
tv = stats.timevect;
for d=1:2
    f=figure(d);clf
    f.Position = [690 82 1600 621];

    for t=1:4
        s = stats.(durations{d}).(targetlabels{t});
        for tp = 1:length(timestoplot)

            tidx=find(tv==timestoplot(tp));
            %topo
            a1=subplot(4,length(timestoplot),(t-1)*length(timestoplot)+tp);hold on
            axis off

            x = res_searchlight.samples(res_searchlight.sa.targetfeature==t & res_searchlight.sa.soaduration==durnum(d)&res_searchlight.fa.time==tidx);

            ft.avg(:,:) = repmat(x,1,size(ft.avg,2));

            mrange = [.25 .275];

            a1pos = a1.Position;
            aw = a1pos(1)-(tp-1)*.05;
            a=axes('Position',[aw a1pos(2:4)]);
            hold on
            co3 = [linspace(1,co(t,1),100);...
                linspace(1,co(t,2),100);...
                linspace(1,co(t,3),100)]';
            % show figure with plots for each sensor
            cfg = [];
            cfg.zlim = mrange;
            cfg.xlim = tv([tidx tidx+1]);
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
            %             colorbar
            %             title(sprintf('%s, %d ms',targetlabelsplot{t},tv(tidx)))
            drawnow

            if tp == length(timestoplot) % add colorbar at end of row
                origsize = get(gca,'Position');
                c=colorbar;
                c.Position = [c.Position(1)+.03 c.Position(2:4)];
                set(a,'Position',origsize)
            end
        end


    end

    % labels
    targetlabelsplot = {'Orientation','SF','Colour','Contrast'};

    for t = 1:length(targetlabelsplot)
        annotation('textbox',[0.07 0.663-(t-1)*.22 .07 .2],...
            'String',targetlabelsplot{t},'FontSize',20,'LineStyle','none','Color',co(t,:),'HorizontalAlignment','right')
    end

    for p = 1:length(timestoplot)
        annotation('textbox',[0.17+(p-1)*.082 0.8 .2 .2],...
            'String',sprintf('%d ms',timestoplot(p)),'FontSize',20,'LineStyle','none','Color','k')
    end

    % save
    fn = sprintf('./figures/SuppFigure_searchlight_%s',durations{d});
    print(gcf,'-dpng','-r500',fn)
    im=imread([fn '.png']);
    [i,j]=find(mean(im,3)<255);margin=2;
    imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');

end



