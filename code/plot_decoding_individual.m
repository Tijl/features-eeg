%% plot_decoding_indi
addpath('~/Repository/CommonFunctions/Matplotlibcolors2/')

load('../results/stats_decoding.mat')

durations = {'soa150','soa50'};
targetlabels = {'ori','sf','color','contrast'};
tv = stats.timevect;
co = tab10;
f=figure(1);clf
f.Position(3:4) = [1000 700];
for s=1:16
    a=subplot(4,4,s);a.FontSize=12;hold on
    plot(tv,x*0+.25,'k-');
    h=[];
    for feat=1:4
        x = stats.soa150.(targetlabels{feat}).mu_all(s,:);
        h(feat) = plot(tv,x,'Color',co(feat,:),'LineWidth',2);

        d = x(tv<=50); %-100 to 50ms to get as many baseline points as we can
        ps = mean(x<=d');
        ps = ps<fdr(ps);
        plot(tv(ps),tv(ps)*0+.235-.0075*feat,'.','Color',co(feat,:));
    end
    ylim([.2 .5])
    xlim(minmax(tv))
    if s==1
        leg=legend(h,{'Orientation','SF','Colour','Contrast'});
        leg.Position = leg.Position+[0.025 0 0 0];
    end
    title(sprintf('Subject %02i',s))
    if mod(s,4)==1
        ylabel('accuracy')
    end
    if s>12
        xlabel('time (ms)')
    end
end

%% save
fn = '../figures/Figure_individual_decoding';
tn = tempname;
print(gcf,'-dpng','-r500',tn)
im=imread([tn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');
