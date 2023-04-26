%%
addpath('../MATLAB/matplotlib/')

%% get data
load('results/stats_decoding_conjunction.mat')
conj = stats;
load('results/stats_decoding.mat')

combos = fieldnames(conj.soa150);

cnames = {'Colour x Contrast','SF x Contrast','SF x Colour',...
    'Orientation x Contrast','Orientation x Colour','Orientation x SF'};

feats = cellfun(@(x) strsplit(x, '_'), combos, 'UniformOutput', false);
feats = vertcat(feats{:}); % To remove nesting of cell array newA
targetlabels = {'Orientation', 'SF','Colour','Contrast'};

durations = {'soa150','soa50'};h=[];

tv = stats.timevect;

%% set up plotting

co = tab10;
featco = [3 2 2 1 1 1; 4 4 3 4 3 2]'; % map pairs of feats to individual targetlabels

onsetys = [.475 .475 .475 .49 .49 .49];
peakys = [.522 .527 .547 .507 .507 .508];

yl = [.47 .56; .47 .56; .47 .565; .488 .52; .488 .52; .488 .52];

%% plot interaction per plot

for d=1:2 % only plot soa150
    f=figure(d);clf
    f.Position = [f.Position(1:2) 1200 900];

    for fplot = 1:length(combos)

        %% plot decoding accuracy
        plotnr = (ceil(fplot/3)-1)*6+fplot;
        sp = subplot(6,3,[plotnr plotnr+3]);
        sp.Position = [sp.Position(1) sp.Position(2)+.04 sp.Position(3) sp.Position(4) - .04];
        hold on
        set(gca,'FontSize',18)

        % get data
        dat = conj.(durations{d}).(combos{fplot});
        mu = dat.mu;
        se = dat.se;
        onsetci = dat.onsetci;
        peakci = dat.peakci;

        % chance line
        plot([min(tv) max(tv)],[.5 .5],'Color',[.5 .5 .5])

        % plot decoding accuracy
        fill([stats.timevect fliplr(stats.timevect)],[mu+se fliplr(mu-se)],'k','FaceAlpha',.2,'LineStyle','none','HandleVisibility','off');
        plot(stats.timevect,mu,'LineWidth',2,'Color','k');

        % plot onsets
        if ~isnan(dat.onset) % if reliable onset time
            % plot onset and peak CIs
            plot(dat.onset,onsetys(fplot),'.k','MarkerSize',12)
            line(onsetci+[-.5 .5],[onsetys(fplot) onsetys(fplot)],'Color','k','LineWidth',2)
            text(max(onsetci)+10,onsetys(fplot),sprintf(' %s',cnames{fplot}),'Color','k','Fontsize',14);
        end
        drawnow

        % plot peaks
        plot(dat.peak,peakys(fplot),'.k','MarkerSize',12)
        line(peakci+[-.5 .5],[peakys(fplot) peakys(fplot)],'Color','k','LineWidth',2)
        text(max(peakci)+10,peakys(fplot),sprintf(' %s',cnames{fplot}),'Color','k','Fontsize',14);

        % now do onsets and peaks for individual features
        if fplot<4
            offs = .007;
        else
            offs = .0025;
        end
        for fe=1:2 % each feature
            fname = feats{fplot,fe};
            fnamelabel = targetlabels{featco(fplot,fe)};

            fdat = stats.(durations{d}).(fname);

            % onset
            plot(fdat.onset,onsetys(fplot)+offs*fe,'.','Color',co(featco(fplot,fe),:),'MarkerSize',12)
            line(fdat.onsetci+[-.5 .5],[onsetys(fplot)+offs*fe onsetys(fplot)+offs*fe],'Color',co(featco(fplot,fe),:),'LineWidth',2)
            text(max(fdat.onsetci)+10,onsetys(fplot)+offs*fe,sprintf(' %s',fnamelabel),'Color',co(featco(fplot,fe),:),'Fontsize',14);

            % peak
            plot(fdat.peak,peakys(fplot)+offs*fe,'.','Color',co(featco(fplot,fe),:),'MarkerSize',12)
            line(fdat.peakci+[-.5 .5],[peakys(fplot)+offs*fe peakys(fplot)+offs*fe],'Color',co(featco(fplot,fe),:),'LineWidth',2)
            text(max(fdat.peakci)+10,peakys(fplot)+offs*fe,sprintf(' %s',fnamelabel),'Color',co(featco(fplot,fe),:),'Fontsize',14);
        end

        text(min(fdat.onsetci)-120,onsetys(fplot)+offs,'onset','Color','k','Fontsize',14);
        text(min(fdat.peakci)-110,peakys(fplot)+offs,'peak','Color','k','Fontsize',14);

        title(cnames{fplot})
        xlim([-100 600]);
        ylim(yl(fplot,:))

        if ismember(fplot,[1,4])
            ylabel('Accuracy');
        end

        %% now bfs

        supb = subplot(6,3,(ceil(fplot/3)-1)*6+6+fplot);
        supb.Position = [supb.Position(1) supb.Position(2)+.07 supb.Position(3) supb.Position(4)-.04];
        supb.FontSize=18;
        hold on

        s=dat;

        %bf
        plot(tv,1+0*tv,'k-');
        co3 = [.5 .5 .5;1 1 1;0 0 0];
        idx = [s.bf<1/10,1/10<s.bf & s.bf<10,s.bf>10]';
        for i=1:3
            x = tv(idx(i,:));
            y = s.bf(idx(i,:));
            if ~isempty(x)
                stem(x,y,'Marker','o','Color',.6*[1 1 1],'BaseValue',1,'MarkerSize',5,'MarkerFaceColor',co3(i,:),'Clipping','off');
                plot(x,y,'o','Color',.6*[1 1 1],'MarkerSize',5,'MarkerFaceColor',co3(i,:),'Clipping','off');
            end
        end
        supb.YScale='log';
        ylim(10.^(6*[-1 1]))
        supb.YTick = 10.^([-5 0 5]);
        xlim([-100 600]);
        if ismember(fplot,[1,4])
            ylabel('BF')
        end
        xlabel('Time (ms)')

    end

    %% save
    if d ==1 
        fn = sprintf('figures/Figure4_%s',durations{d});
    else
        fn = sprintf('figures/Supp_conjunction_%s',durations{d});
    end
    
    print(gcf,'-dpng','-r500',fn)
    im=imread([fn '.png']);
    [i,j]=find(mean(im,3)<255);margin=2;
    imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');

end
