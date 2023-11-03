%% plot feature x feature interaction peaks with stats

durations = {'soa150','soa50'};
targetlabels = {'ori','sf','color','contrast'};
targetnames = {'Orientation','SF','Colour','Contrast'};

% set colour gradients, one per feature
co = tab10;co3=[];
for c = 1:4
    co3(:,:,c) = [linspace(1,co(c,1),5);...
        linspace(1,co(c,2),5);...
        linspace(1,co(c,3),5)]';
end
co3(1,:,:) = [];
% co3 = flipud(co3);

% feature-specific details
scales = [.25 .28; .25 .45; .25 .45; .25 .4]; % ylim for each feature
peaktimes = [120 112 112 96; 120 112 112 100];

load('results/stats_decoding_interactions.mat')

%% plot

for d=1:2
    f=figure(d);clf
    f.Position = [f.Position(1:2) 1000 860];
    x=0;
    for fplot = 1:length(targetlabels) % decoding contrast

        timestoplot = find(stats.timevect==peaktimes(d,fplot));
        timestoplot = (timestoplot-2):(timestoplot+2);

        for fcont = 1:length(targetlabels) % for each level of

            if fplot == fcont % can't be done

                continue

            else % plot decoding by feature
                
                x=x+1;
                subplot(4,3,x)
                set(gca,'FontSize',18)
                hold on

                for lev = 1:4

                    dat_all = stats.(durations{d}).(targetlabels{fplot}).(targetlabels{fcont}).(sprintf('level%d',lev)).mu_all;

                    % calculate mean and se for that time window
                    dat_mu = mean(dat_all(:,timestoplot),2);
                    se(lev) = std(dat_mu)/sqrt(length(dat_mu));

                    peaks(lev) = mean(dat_mu);
                    names{lev} = sprintf('%s %d',targetlabels{fcont},lev);

                    h(lev) = bar(lev,peaks(lev));
                    set(h(lev),'FaceColor', co3(lev,:,fplot));

                end

                errorbar(1:4,peaks,se,se,'.k','HandleVisibility','off','Marker','none')
            end


            % add bits and pieces
            set(gca,'XTick',1:4)
            xlim([.4 4.6])
            if fplot==1
                set(gca,'YTick',0.25:.01:.45);
            else
                set(gca,'YTick',0.25:.05:.45)
            end

            if mod(x,3)==1
                ylabel({'Decoding'; sprintf('%s',targetnames{fplot})},'Color',co(fplot,:))
            else
                set(gca,'YTickLabel','')
            end
            xlabel(sprintf('%s per %s',targetnames{fplot},targetnames{fcont}))
            ylim(scales(fplot,:))
        end
    end

    %% saving
    if d == 1
        fn = sprintf('figures/Figure4_decoding_interaction_peak%s',durations{d});
    else
        fn = sprintf('figures/Supp_decoding_interaction_peak%s',durations{d});        
    end
    print(gcf,'-dpng','-r500',fn)
    im=imread([fn '.png']);
    [i,j]=find(mean(im,3)<255);margin=2;
    imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');
    
    exportgraphics(gcf,[fn '.pdf'],'ContentType','vector')
end

