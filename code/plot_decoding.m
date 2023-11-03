%% add paths
addpath('~/Repository/CommonFunctions/matplotlib/')
addpath('~/CoSMoMVPA/mvpa/')
addpath('~/fieldtrip');

ft_defaults;

%% collate searchlight results

res_cell = {};
for s=1:16
    x=load(sprintf('./results/sub-%02i_decoding_searchlight.mat',s));
    res_cell{s} = x.res;
end
res_searchlight = cosmo_fx(cosmo_stack(res_cell),@mean,{'soaduration','targetfeature'});
res_searchlight.sa.d = 1+ (res_searchlight.sa.soaduration>0.1);

%% get layout
layout=cosmo_meeg_find_layout(res_searchlight,'label_threshold',.80);

cfg={};
ft = ft_timelockanalysis(cfg,cosmo_map2meeg(res_searchlight));

%% plot
load('./results/stats_decoding.mat')

f=figure(1);clf
f.Position = [f.Position(1:2) 800 1000];
durations = {'soa150','soa50'};h=[];
targetlabels = {'ori','sf','color','contrast'};
targetlabelsplot = {'Orientation','SF','Colour','Contrast'};

tv = stats.timevect;
ha=[];
for d=1:2
    clear a
    a=subplot(3,2,d);hold on
    a.YLim=[.23 .5];
    a.FontSize=12;
    % chance line
    plot([min(tv) max(tv)],[.25 .25],'Color',[.5 .5 .5])
    co = tab10; %vega10;
    for t=1:4

        % plot decoding and standard error
        s = stats.(durations{d}).(targetlabels{t});
        fill([tv fliplr(tv)],[s.mu-s.se fliplr(s.mu+s.se)],co(t,:),...
            'FaceAlpha',.2,'LineStyle','none');
        h(t) = plot(tv,s.mu,...
            'DisplayName',targetlabelsplot{t},'Color',co(t,:),...
            'LineWidth',2);

        % plot onset and peak CIs
        yval = 0.47-t*.015;
        plot([s.onset s.peak],[yval yval],'.','Color',co(t,:),'MarkerSize',12)
        line(s.onsetci+[-.5 .5],[yval yval],'Color',co(t,:),'LineWidth',2)
        line(s.peakci+[-.5 .5],[yval yval],'Color',co(t,:),'LineWidth',2)
        drawnow

    end
    text(min(s.onsetci)-100,0.47-.015,'onset')
    text(max(s.peakci)+40,0.47-.015,'peak')

    legend(h)
    ylabel('Accuracy')
    xlabel('Time (ms)')
    xlim(minmax(tv))
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
        cfg.xlim = tv([tidx tidx+1]);
        cfg.layout = layout;
        cfg.showscale = 'no';
        cfg.figure = a;
        cfg.interactive ='no';
        cfg.comment = 'no';
        cfg.markersymbol = '.';
        cfg.style = 'straight';
        cfg.gridscale = 128;
        ft_info off verbose
        ft_topoplotER(cfg, ft);
        a.FontSize = 12;
        a.Colormap = co3;
        set(a.Children,'LineWidth',.5)
        ha(t,d)=a;
    end

    %bf
    for t=1:4
        a=subplot(15,2,11 + 2*t + (d==2));hold on
        a.FontSize=12;
        s = stats.(durations{d}).(targetlabels{t});
        plot(tv,1+0*tv,'k-');
        co3 = [.5 .5 .5;1 1 1;co(t,:)];
        idx = [s.bf<1/10,1/10<s.bf & s.bf<10,s.bf>10]';
        for i=1:3
            x = tv(idx(i,:));
            y = s.bf(idx(i,:));
            if ~isempty(x)
                stem(x,y,'Marker','o','Color',.6*[1 1 1],'BaseValue',1,'MarkerSize',5,'MarkerFaceColor',co3(i,:),'Clipping','off');
                plot(x,y,'o','Color',.6*[1 1 1],'MarkerSize',5,'MarkerFaceColor',co3(i,:),'Clipping','off');
            end
        end
        a.YScale='log';
        ylim(10.^(6*[-1 1]))
        a.YTick = 10.^([-5 0 5]);
        xlim(minmax(tv))
        ylabel('BF')
        if t==4
            xlabel('Time (ms)')
        else
            a.XTick = [];
        end
    end
end
for d=1:2
    for t=1:4
        co3 = [linspace(1,co(t,1),100);...
            linspace(1,co(t,2),100);...
            linspace(1,co(t,3),100)]';
        set(ha(t,d),'Colormap',co3);
    end
end

%% timegens

load('./results/stats_timegen.mat')

%cmap = brewermap(1000,'PiYG');
cmap = cividis();

timv = stats.timevect;
timerv = fliplr(timv);
xlab = 0:200:600;
ylab = 600:-200:0;

ts = [0.24 0.26; .22 .28; .22 .28; .22 .28]; % scales

for d= 1
    for f = 1:length(targetlabels)

        a=subplot(7,4,20+f);hold on

        dat = stats.(durations{d}).(targetlabels{f});
        mu = reshape(dat.mu,length(timv),length(timv));
        bf = reshape(dat.bf,length(timv),length(timv));
        mu(bf<10) = .25;

        %%% plot time gen decoding
        ax=imagesc(timv,timv,(mu),ts(f,:));
        colormap(gca,cmap)
        set(gca,'YDir','normal')
        cb=colorbar();
        cb.Ticks = [ts(f,1) 0.25 ts(f,2)];
        xlim([0 400])
        ylim([0 400])
        set(gca,'XTick',[0:200:400 max(timv)],'XTickLabel',num2cell(xlab));
        set(gca,'YTick',[0:200:400 max(timv)],'YTickLabel',num2cell(xlab));
        if f>1
            set(gca,'YTickLabel',[])
        end

        xlabel('Training time (ms)');
        if f ==1
            ylabel('Testing time (ms)');
        end
        axis square
        title(targetlabelsplot{f},'Color',co(f,:))
        set(gca,'Fontsize',12)

    end
end

%% annotate

annotation('textbox',[0.05 .86 .1 .1],...
    'String','A','FontSize',20,'LineStyle','none')

annotation('textbox',[0.05 .255 .1 .1],...
    'String','B','FontSize',20,'LineStyle','none')

%% save
fn = './figures/Figure2';
print(gcf,'-dpng','-r500',fn)
im=imread([fn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');

