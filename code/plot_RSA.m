%%
addpath('~/Repository/CommonFunctions/matplotlib/')
addpath('~/Repository/CommonFunctions/tSNE_matlab/')
addpath('~/Repository/CommonFunctions/CircStat/')

%% load neural data
durations = [150 50];
if exist('../results/EEG_RDMs.mat')
    fprintf('loading ../results/EEG_RDMs.mat')
    load('../results/EEG_RDMs.mat')
else
    rdm_cell_all = {};
    for d = 1:2
        for s=1:16
            x=load(sprintf('../results/sub-%02i_decoding_rdm_%ims.mat',s,durations(d)));
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
    save('../results/EEG_RDMs.mat','groupRDMsoa150','groupRDMsoa50','meanRDMsoa150','meanRDMsoa50','tv','stimtable','-v7.3');
end

%% read stimuli
stimfns = dir('../stimuli_small/stim*');
stims={};
IM={};
for s=1:256
    stims{s,1} = fullfile('stimuli_small',stimfns(s).name);
    [IM{s,1}] = imread(fullfile('../stimuli',stimfns(s).name));
end

%% load behavioural data
if exist('../results/stats_behaviour.mat')
    load('../results/stats_behaviour.mat','RDMmean','T')
else
    fns = dir('../data_online_triplets/jatos_results_*.txt');
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
    save('../results/stats_behaviour.mat','RDMmean','T')
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
if exist('../results/stats_rsa.mat')
    load('../results/stats_rsa.mat','rmodel','rmodel_bf')
else
    % models
    fprintf('Computing model correlations')
    rdm_cell = {groupRDMsoa150,groupRDMsoa50};
    rmodel = [];
    rmodel_bf = [];
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
    end
    fprintf('done\n')
    % save RDM
    save('../results/stats_rsa.mat','rmodel','rmodel_bf')
end

%% make figure
f=figure(1);clf
f.Position = [f.Position(1:2) 800 700];
f.Resize='off';
nsub = numel(unique(T.subjectnr));
tv = -100:4:596;

% mds
a=subplot(2,2,2);hold on
a.FontSize=12;
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
title({'2-dimensional embedding of similarities','(mdscale-metricstress)'})
axis off

%model mds
a=subplot(2,2,1);
a.FontSize=16;
rng(10)
R = corr(modelvec,'type','spearman');
Y = mdscale(squareform(pdist(modelvec','spearman')),2);
aw = .9*[-1 1]*max(range(minmax(Y')));

for i=5
    for j=1:4
        line(Y([i j],1)',Y([i j],2)','LineWidth',20*abs(R(i,j)),'Color','k');hold on    
    end
end

modelnames = {'orientation','SF','colour','contrast','behaviour'};
for i=1:5
    x = Y(i,1);
    y = Y(i,2);
    plot(x,y,'o','Color',co(i,:));hold on
    axis equal
    image('XData',x+aw,'YData',y+aw,'CData',squareform(255*normalize(modelvec(:,i),'range')))
    text(x,max(y+aw),modelnames{i},'Color',co(i,:),'FontSize',16,'VerticalAlignment','bottom','HorizontalAlignment','center')
end
axis off

% plot
plotnr = 4;
for d=1:2
    plotnr=plotnr+1;
    a = subplot(4,2,plotnr);hold on
    a.FontSize=12;a.ColorOrder=co;

    h=[];
    for i=1:5

        mu = squeeze(mean(rmodel(:,:,i,d)));
        se = squeeze(std(rmodel(:,:,i,d)))./sqrt(size(rmodel,1)-1);

        fill([tv fliplr(tv)],[mu-se fliplr(mu+se)],co(i,:),...
            'FaceAlpha',.2,'LineStyle','none');
        h(i) = plot(tv,mu,'Color',co(i,:),...
            'LineWidth',2);
    end
    legend(h,modelnames)
    ylim([-.02 0.25])
    xlim([-100 600])
    xlabel('Time (ms)')

    title(sprintf('%.2f Hz',1000/durations(d)))
    ylabel('model-EEG correlation')
    drawnow
end

% plot BF
plotnr = 6;
for d=1:2
    plotnr=plotnr+1;
    a = subplot(4,2,plotnr); %make a plot, then split it in 5 plots later
    a.FontSize=12;
    axis off
    apos = a.Position;
    co = tab10(5);
    XX=[];
    for i=1:5
        a2 = axes('Position',[apos(1),apos(2)+apos(4)-i*apos(4)./5,apos(3),apos(4)./5]);
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
        if mod(i,2)
            a2.YTick = 10.^(5*[-1  1]);
            a2.Color = .9*[1 1 1];
        else
            a2.YTick =[];
        end
        a2.XLim = minmax(tv);
        if i==3
            ylabel('Bayes Factor')
        end
        if i==5
            xlabel('Time (ms)')
        else
            a2.XTick = [];
        end

    end
end


%% save
fn = '../figures/rsa';
tn = tempname;
print(gcf,'-dpng','-r500',tn)
im=imread([tn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');

