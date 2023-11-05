
%% get data
load('results/stats_decoding_conjunction_eachcontrast.mat')
conj = stats;
load('results/stats_decoding_conjunction.mat')
conj_coll = stats;
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

onsetys = [.465 .465 .465 .48 .48 .48];
peakys = [.55 .57 .615 .515 .515 .515];

yl = [.46 .64; .46 .64; .46 .64; .475 .53; .475 .53; .475 .53];

%% plot interaction per plot
cols = viridis(36);

for d=1:2 % soas
    f=figure(d);clf
    f.Position = [f.Position(1:2) 1200 900];

    for fplot = 1:length(combos)

        %% plot decoding accuracy
        sp = subplot(2,3,fplot);
        sp.Position = [sp.Position(1) sp.Position(2)+.04 sp.Position(3) sp.Position(4) - .04];
        hold on
        set(gca,'FontSize',18)
        
        % chance line
        plot([min(tv) max(tv)],[.5 .5],'Color',[.5 .5 .5])
        
        for n = 1:36
            % get data
            dat = conj.(durations{d}).(combos{fplot}).(sprintf('comb%d',n));

            allmu(n,:) = dat.mu;
%             allse(n,:) = dat.se;
            peak(n) = max(dat.mu);

        end
        
        [~,j]=sort(peak,'descend'); % get order of max strength
        peakorder.(durations{d}).(combos{fplot}) = j;

        for n = 1:36
            % plot decoding accuracy in order of max strength
            plot(stats.timevect,allmu(j(n),:),'LineWidth',1,'Color',cols(n,:));
        end

        % plot onsets
        datall = conj_coll.(durations{d}).(combos{fplot});
        onsetci = datall.onsetci;
        peakci = datall.peakci;

        if ~isnan(datall.onset) % if reliable onset time
            % plot onset and peak CIs
            plot(datall.onset,onsetys(fplot),'.k','MarkerSize',12)
            line(onsetci+[-.5 .5],[onsetys(fplot) onsetys(fplot)],'Color','k','LineWidth',2)
            text(max(onsetci)+10,onsetys(fplot),sprintf(' %s',cnames{fplot}),'Color','k','Fontsize',14);
        end
        drawnow

        % plot peaks
        plot(datall.peak,peakys(fplot),'.k','MarkerSize',12)
        line(peakci+[-.5 .5],[peakys(fplot) peakys(fplot)],'Color','k','LineWidth',2)
        text(max(peakci)+10,peakys(fplot),sprintf(' %s',cnames{fplot}),'Color','k','Fontsize',14);
        
        % now do onsets and peaks for individual features
        if fplot<4
            offs = .01;
        else
            offs = .0035;
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
        if ismember(fplot,4:6)
            xlabel('Time (ms)')
        end

    end

    %% save
    fn = sprintf('figures/Supp_conjlevels_%s',durations{d});
    print(gcf,'-dpng','-r500',fn)
    im=imread([fn '.png']);
    [i,j]=find(mean(im,3)<255);margin=2;
    imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');

end
