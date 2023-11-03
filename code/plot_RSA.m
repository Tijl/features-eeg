%%
% addpath('~/Repository/CommonFunctions/matplotlib/')
% addpath('~/Repository/CommonFunctions/tSNE_matlab/')
% addpath('~/Repository/CommonFunctions/CircStat/')

%% load neural data
durations = [150 50];
if exist('results/EEG_RDMs.mat')
    fprintf('loading results/EEG_RDMs.mat')
    load('results/EEG_RDMs.mat')
else
    rdm_cell_all = {};
    for d = 1:2
        for s=1:16
            x=load(sprintf('results/sub-%02i_decoding_rdm_%ims.mat',s,durations(d)));
            rdm = cosmo_average_samples(x.rdm,'split_by',{'target1','target2'});
            stimtable = x.stimtable;
            stimtable.f_colour = stimtable.f_color;
            tv = rdm.a.fdim.values{1};
            
            X = nan(175,256,256);
            for r = 1:length(rdm.sa.target1)
                t1 = 1+rdm.sa.target1(r);
                t2 = 1+rdm.sa.target2(r);
                X(:,t1,t2) = rdm.samples(r,:);
                X(:,t2,t1) = rdm.samples(r,:);
            end
            rdm_cell_all{d,s} = X
        end
    end
    groupRDMsoa150 = permute(cat(4,rdm_cell_all{1,:}),[4 1 2 3]);
    groupRDMsoa50 = permute(cat(4,rdm_cell_all{2,:}),[4 1 2 3]);
    
    meanRDMsoa150 = squeeze(mean(groupRDMsoa150));
    meanRDMsoa50 = squeeze(mean(groupRDMsoa50));
    
    % save RDM (so loading is faster on second run)
    save('results/EEG_RDMs.mat','groupRDMsoa150','groupRDMsoa50','meanRDMsoa150','meanRDMsoa50','tv','stimtable','-v7.3');
end

%% read stimuli
stimfns = dir('stimuli_small/stim*');
stims={};
IM={};
for s=1:256
    stims{s,1} = fullfile('stimuli_small',stimfns(s).name);
    [IM{s,1}] = imread(fullfile('stimuli',stimfns(s).name));
end

%% load behavioural data
if exist('results/stats_behaviour.mat')
    load('results/stats_behaviour.mat','RDMmean','T')
else
    fns = dir('data_online_triplets/jatos_results_*.txt');
    T = [];subnr = 0;
    for f=1:numel(fns)
        fn = fullfile(fns(f).folder,fns(f).name)
        lines = strsplit(fileread(fn),'\n');
        for i=1:(numel(lines)-1)
            txt = lines{i}(2:end-1);
            lines2 = strsplit(txt(2:end-1),'},{');
            TS = struct();idx = 0;
            for j=1:numel(lines2)
                X = jsondecode(['{',lines2{j},'}']);
                if isfield(X,'test_part')
                    idx = idx+1;
                    F = fieldnames(X);
                    for jj = 1:numel(F)
                        TS(idx).(F{jj}) = X.(F{jj});
                    end
                end
            end
            subnr = subnr+1;
            TS = struct2table(TS);
            TS.subjectnr(:) = subnr;
            TS.stimulus = [];
            T = [T; TS];
            if ~mod(i,10)
                fprintf('line %i/%i\n',i,numel(lines));
            end
        end
    end
    
    % re-map stimuli
    fprintf('re-map\n')
    [~,T.stim0number]=ismember(T.stim0,stims);
    [~,T.stim1number]=ismember(T.stim1,stims);
    [~,T.stim2number]=ismember(T.stim2,stims);
    fprintf('Done\n')
    
    % make RDM
    fprintf('create RDM\n')
    allcombs = [T.stim0number T.stim1number T.stim2number];
    choice = T.button_pressed;
    RDMsum = zeros(256);
    RDMcounts = zeros(256);
    for i=1:numel(choice)
        % for every choice, add 1 to the similarity of the two items that were
        % not the odd-one-out (three choices were coded as 0, 1, 2)
        v = allcombs(i,(0:2)~=choice(i));
        RDMsum(v(1),v(2)) = RDMsum(v(1),v(2))+1;
        RDMsum(v(2),v(1)) = RDMsum(v(2),v(1))+1;
        % add 1 to the counts of all items compared, to compute the mean later
        for v = combnk(allcombs(i,:),2)'
            RDMcounts(v(1),v(2)) = RDMcounts(v(1),v(2))+1;
            RDMcounts(v(2),v(1)) = RDMcounts(v(2),v(1))+1;
        end
    end
    RDMmean = 1 - RDMsum./RDMcounts;
    fprintf('Done\n')
    
    % save RDM
    save('results/stats_behaviour.mat','RDMmean','T')
end

%% create embedding
fprintf('create embedding\n')
rng(1)
Y = RDMmean;
Y(eye(size(Y))==1)=0;
X1 = mdscale(Y,2,'Start','random','Criterion','metricstress');

%% models
modelnames = {'ori','sf','colour','contrast','behaviour'};
modelvec = [];models={};
loweridx = find(tril(ones(256),-1)==1);
for m=1:5
    if m==5
        models{m} = RDMmean;
    elseif m==1
        ori = stimtable.feature_ori;
        %make ori RDM circular
        r = min(cat(3,abs(ori - ori'),abs(ori - (ori-180)'),abs(ori - (ori+180)')),[],3);
        models{m} = r;
    else
        models{m} = squareform(pdist(1+stimtable.(['f_' modelnames{m}]),'euclidean'));
    end
    modelvec(:,m) = models{m}(loweridx);
end

%% compute model correlations with EEG and their BF
if exist('results/stats_rsa.mat')
    load('results/stats_rsa.mat','rmodel','rmodel_bf','beh')
else
    % models
    fprintf('Computing model correlations/n')
    rdm_cell = {groupRDMsoa150,groupRDMsoa50};
    rmodel = [];
    rmodel_bf = [];
    beh = struct();
    
    for d=1:2
        rdm = rdm_cell{d};
        loweridx = find(tril(ones(256),-1)==1);
        for m=1:5
            for s=1:size(rdm,1)
                rmodel(s,:,m,d) = corr(squeeze(rdm(s,:,loweridx))', modelvec(:,m),'type','Spearman');
            end
            rmodel_bf(:,m,d) = bayesfactor_R_wrapper(rmodel(:,:,m,d)','returnindex',2,'verbose',false,'args','mu=0,rscale="medium",nullInterval=c(-Inf,0.5)');
            fprintf('.')
        end
        
        % calculate onsets and peaks for behavioural model
        r_beh = squeeze(mean(rmodel(:,:,5,d),1));
        r_bf = squeeze(rmodel_bf(:,5,d));
        
        n = size(rmodel,1);
        
        s=struct();
        s.mu = r_beh;
        s.mu_all = rmodel(:,:,5,d);
        s.bf = r_bf;
        
        % group mean onset & peak
        [~,ot] = max(movmean(s.bf>10,[0 2])==1); % onset 3 consecutive tp BF>10
        [~,peak] = max(s.mu); % peak
        
        % use jackknife method and calculate bayesfactors to get onsets
        xx = [];
        combs = combnk(1:n,n-2); %leave-2-out
        for i=1:size(combs,1)
            xx = [xx s.mu_all(combs(i,:),:)];
        end
        fprintf('Calculating jackknifed BFs for %d ms duration...\n',durations(d))
        b = bayesfactor_R_wrapper(xx',...
            'returnindex',2,'verbose',false,'args','mu=0,rscale="medium",nullInterval=c(-Inf,0.5)');
        s.bf_jackknife = reshape(b,[],size(combs,1));
        
        % get onsets and peaks per jackknife
        onset_jack = [];peak_jack=[];
        for b=1:size(s.bf_jackknife,2)
            % get onset (first three consecutive time points BF>10)
            [~,onset_jack(b)] = max(movmean(s.bf_jackknife(:,b)>10,[0 2])==1);
            
            % get peak
            [~,peak_jack(b)] = max(mean(s.mu_all(combs(b,:),:)));
        end
        
        % convert onsets and peaks to time in ms (rather than time point number)
        s.onset = tv(ot);
        s.onsetjk = tv(onset_jack);
        if isempty(onset_jack)
            s.onsetci=[];
        else
            ci = prctile(onset_jack,[5 95]);
            s.onsetci = tv([floor(ci(1)) ceil(ci(2))]);
        end
        s.peak = tv(peak);
        s.peakjk = tv(peak_jack);
        s.peakci = tv(round(prctile(peak_jack,[5 95])));
        
        beh.(sprintf('soa%d',durations(d))) = s;
        
    end
    fprintf('done\n')
    % save RDM
    save('results/stats_rsa.mat','rmodel','rmodel_bf','beh')
end


%% make models to link to behaviour

x=load(sprintf('results/sub-%02i_decoding_rdm_%ims.mat',1,durations(1)));
stimtable = x.stimtable;
stimtable = stimtable(:,3:6);

% two way conjunction models
modcomb = combnk(1:4,2);

mnum = 5;
for mc = 1:length(modcomb)
    
    m1 = table2array(stimtable(:,modcomb(mc,1)));
    m2 = table2array(stimtable(:,modcomb(mc,2)));
    m3 = ones(256,256);
    
    for s = 1:size(m1,1) % stim1 to be compared
        for s2 = 1:size(m1,1) % stim2 to be compared
            if (m1(s) == m1(s2)) && (m2(s) == m2(s2))
                m3(s,s2) = 0;
            end
        end
    end
    
    mnum = mnum+1;
    modelvec(:,mnum) = squareform(m3);
    modelnames{mnum} = sprintf('%sX%s',modelnames{modcomb(mc,1)},modelnames{modcomb(mc,2)});
end

mnames = {'Orientation' 'Spatial frequency' 'Colour' 'Contrast' 'Perceptual' ...
    'Colour x Contrast' 'SF x Contrast' 'SF x Colour'...
    'Orientation x Contrast' 'Orientation x Colour' 'Orientation x SF'};

%% make figure
f=figure(1);clf
f.Position = [f.Position(1:2) 800 800];
f.Resize='off';

co = tab10;%vega10;

for m = 1:5 % feature models
    s=subplot(7,5,m)
    s.Position(1) = s.Position(1)-.03;
    s.Position(4) = s.Position(4)+.03;
    imagesc(squareform(modelvec(:,m)))
    
    % colour
    co3 = [linspace(1,co(m,1),4);...
        linspace(1,co(m,2),4);...
        linspace(1,co(m,3),4)]';
    colormap(gca,co3)
    colormap(gray)
    axis square
    set(gca,'Fontsize',12)
    title(mnames{m},'Color',co(m,:))
    axis off
end

%% plot behaviour mds
a=subplot(7,2,[5 5+2 5+4 5+6 5+8]);hold on
a.FontSize=12;
a.Position(1) = a.Position(1)-.03;
a.Position(2) = a.Position(2)-.1-.02;
a.Position(4) = a.Position(4) + .03;
X=X1;
aw=range(X(:))*.02;
axis equal
axis square
imrange = 340:450;
plot(X(:,1),X(:,2),'.')
for i=1:256
    image('XData',X(i,1)+aw*[-1 1],'YData',X(i,2)+aw*[-1 1],'CData',IM{i}(imrange,imrange,1:3))
end
a.XTick = [];
a.YTick = [];
% title({'Perceptual similarity'})
axis off

%% plot correlations with behaviour

[r p] = corr(modelvec,'type','Spearman');

a=subplot(7,4,5)
a.Position = [a.Position(1)-.03 a.Position(2)-.09-.02 a.Position(3) a.Position(4)+.05]
set(gca, 'Fontsize',12)
b=bar(r(5,1:4));
for f = 1:4
    b.FaceColor = 'flat' ; b.CData (f,:) = co(f,:);
end
set(gca,'XTickLabel',mnames([1:4 6:end]))
ylabel('Correlation')
box off
ylim([0 0.6])
title('Features')

a=subplot(7,4,6)
a.Position = [a.Position(1)-.03 a.Position(2)-.09-.02 a.Position(3) a.Position(4)+.05]
set(gca, 'Fontsize',12)
b=bar(r(5,6:11));
for f = 1:6
    b.FaceColor = 'flat' ; b.CData (f,:) = [0 0 0];
end
set(gca,'XTickLabel',mnames(6:end))
% ylabel('Correlation')
title('Conjunctions')
box off
ylim([0 0.6])

%% plot neural-behaviour correlations
plotnr = 3;
co = tab10(5);

for d=1:2
    plotnr=plotnr+(d-1)*5+1;
    a = subplot(7,2,[plotnr plotnr+2]);hold on
    a.Position(2) = a.Position(2)+0.02*(d-1)-.02
    a.Position(4) = a.Position(4)-.02;
    a.FontSize=12;a.ColorOrder=co;
    
    h=[];
    i=5;
    s= beh.(sprintf('soa%d',durations(d)));
    
    mu = s.mu;
    se = std(s.mu_all)./sqrt(size(s.mu_all,1));
    
    fill([tv fliplr(tv)],[mu-se fliplr(mu+se)],co(i,:),...
        'FaceAlpha',.2,'LineStyle','none','HandleVisibility','off');
    h(i) = plot(tv,mu,'Color',co(i,:),...
        'LineWidth',2);
    
    % plot onset and peak CIs
    yval = 0.24;
    plot([s.onset s.peak],[yval yval],'.','Color',co(i,:),'MarkerSize',12)
    line(s.onsetci+[-.5 .5],[yval yval],'Color',co(i,:),'LineWidth',2)
    line(s.peakci+[-.5 .5],[yval yval],'Color',co(i,:),'LineWidth',2)
    drawnow
    
    text(min(s.onsetci)-100,yval,'onset')
    text(max(s.peakci)+30,yval,'peak')
    
    a.YTick = 0:.05:.25;
    ylim([-.02 0.25])
    xlim([-100 600])
    xlabel('Time (ms)')
    
    legend(sprintf('%.2f Hz',1000/durations(d)))
    %     t=title(sprintf('%.2f Hz EEG-behaviour correlation',1000/durations(d)));
    %     t.Position(2) = t.Position(2)+.005;
    ylabel('Correlation')
    drawnow
end

%% plot BF
plotnr = 7;
for d=1:2
    plotnr=plotnr+1+(d-1)*5;
    a = subplot(7,2,plotnr); %make a plot, then split it in 5 plots later
    a.Position(4) = a.Position(4)+.07;
    a.Position(2) = a.Position(2)-.07+.03*(d-1)-.02;
    a.FontSize=12;
    axis off
    apos = a.Position;
    apos(2) = apos(2)-.015;
    XX=[];
    i=5
    a2 = axes('Position',[apos(1),apos(2)+apos(4)-apos(4)./3,apos(3),apos(4)./5]);
    a2.FontSize = 12;hold on
    
    bf = rmodel_bf(:,i,d);
    
    plot(tv,1+0*tv,'k-');
    
    co3 = [.5 .5 .5;1 1 1;co(i,:)];
    idx = [bf<1/10,1/10<bf & bf<10,bf>10]';
    for c=1:3
        x = tv(idx(c,:));
        y = bf(idx(c,:));
        if ~isempty(x)
            stem(x,y,'Marker','o','Color',.6*[1 1 1],'BaseValue',1,'MarkerSize',5,'MarkerFaceColor',co3(c,:),'Clipping','off');
            plot(x,y,'o','Color',.6*[1 1 1],'MarkerSize',5,'MarkerFaceColor',co3(c,:),'Clipping','off');
        end
    end
    a2.YScale='log';
    a2.YLim = 10.^(6*[-1 1]);
    a2.YTick = 10.^(5*[-1  1]);
    a2.Color = .9*[1 1 1];
    a2.XLim = minmax(tv);
    ylabel('Bayes Factor')
    xlabel('Time (ms)')
end


%% add plot labels

annotation('textbox',[0.05 .905 .1 .1],...
    'String','A. Stimulus models','FontSize',18,'LineStyle','none')

annotation('textbox',[0.05 .71 .1 .1],...
    'String','B. Correlation with perceptual model','FontSize',18,'LineStyle','none')

annotation('textbox',[0.5 .71 .1 .1],...
    'String','D. EEG-perceptual correlation','FontSize',18,'LineStyle','none')

annotation('textbox',[0.05 .4 .1 .1],...
    'String','C. Perceptual similarity','FontSize',18,'LineStyle','none')

%% save
fn = 'figures/Figure6';
tn = tempname;
print(gcf,'-dpng','-r500',tn)
im=imread([tn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');
