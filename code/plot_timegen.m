%% cosmo
if isempty(which('cosmo_wtf'))
    addpath(genpath('~/Dropbox/MATLAB/CoSMoMVPA'))
end

%% set up details

resfiles = dir('results/*decoding_timegen.mat');
features = {'ori' 'sf' 'color' 'contrast'};
duration = [0.05 0.15];
fname = 'results/timegen.mat';
collate=0; %rerun collation of data from individual participants

% get data
if ~collate
    load(fname)
else
    % concatenate different participants' data
    for p = 1:length(resfiles)
        fprintf('Loading %s\n',resfiles(p).name);
        m1 = matfile(sprintf('%s/%s',resfiles(p).folder,resfiles(p).name));
        dat = m1.res_cell;
        res = cosmo_stack(dat);

        for d = 1:length(duration)
            for f = 1:length(features)

                fdat = cosmo_slice(res,res.sa.targetfeature==f&res.sa.soaduration==duration(d));

                alltimegen.(sprintf('soa%d',duration(d)*1000)).(features{f})(:,:,p) = cosmo_unflatten(fdat,1);
            end
        end
    end
    timevect = res.a.sdim.values{1};
    save(fname,'alltimegen','timevect');
end

%% plotting details
% colour stuff
addpath('~/Dropbox (Personal)/MATLAB/Colormaps/Colormaps (5)/Colormaps')
addpath('~/Dropbox (Personal)/MATLAB/mseb')
addpath('~/Dropbox (Personal)/MATLAB/altmany-export_fig')

% font sizes
fs=14; % axes
tfs=18; % titles
sfs=10; % stim labels

% plotting stuff
colormap viridis
cmap = viridis(5);
% cmap(3,:) = [0 0 0];

bfcutoff = [inf 10 3 1 .33 0];
bfcols = [1 1 0 0 1]; % whether to plot filled or open circles for each one

%% timegen figures

for d = 1:length(duration)

    % Visualize results
    f=figure(d+1);clf;
    set(f,'Resize','off')
    f.Position=[-1500 373 1600 1050];

    ts = [0.23 0.3]; % scale

    timv = timevect;
    timerv = fliplr(timv);
    xlab = -100:100:600;
    ylab = 600:-100:-100;

    for f = 1:length(features)

        dat = alltimegen.(sprintf('soa%d',duration(d)*1000)).(features{f});
        mu = mean(dat,3);

        supA=subplot(2,2,f);

        %%% plot time gen decoding
        ax=imagesc(flipud(mu),ts);
        colormap(gca,viridis())
        cb=colorbar();
        cb.Ticks = [ts(1) 0.25 ts(2)];

        dlabels = {sprintf('%.2f below chance',ts(1)) '0.25 chance' sprintf('%.2f above chance',ts(2))};
        dlabels = cellfun(@(x) strrep(x,' ','\newline'), dlabels,'UniformOutput',false);
        cb.TickLabels = dlabels;

        set(gca,'XTick',[1:25:length(timv) length(timv)],'XTickLabel',num2cell(xlab));
        set(gca,'YTick',[1:25:length(timv) length(timv)],'YTickLabel',num2cell(ylab));
        ylabel('Training time (ms)');
        xlabel(['Testing time (ms)']);
        set(gca,'FontSize',fs)
        axis square

        hold on

        % plot diagonal and on/off times per training/testing time
        plot(supA.XLim,fliplr(supA.YLim),'Color',[.3 .3 .3],'LineWidth',2);
        plot(supA.XLim,[find(flipud(timerv) == 0) find(flipud(timerv) == 0)],'--k','LineWidth',3);%,'Color','k');%[.3 .3 .3]);
        plot([find(timv == 0) find(timv == 0)],supA.YLim,'--k','LineWidth',3);%,'Color','k');%[.3 .3 .3]);

        title(features{f})

    end
    export_fig(sprintf('figures/time_generalisation_soa%d.png',duration(d)*1000),'-r300')

end
%%
% export_fig('figures/time_generalisation.png','-transparent','-r300')


%% plot diags to check

figure(4);clf;hold on
for f = 1:length(features)

    mu = alltimegen.soa150.(features{f});
    mudiag = diag(mean(mu,3));

    plot(timevect,mudiag,'Color',cmap(f,:),'LineWidth',2)
end