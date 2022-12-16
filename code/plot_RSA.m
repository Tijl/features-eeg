%%
addpath('~/Repository/CommonFunctions/matplotlib/')
addpath('~/Repository/CommonFunctions/tSNE_matlab/')
addpath('~/Repository/CommonFunctions/CircStat/')

%% load neural data
durations = [150 50];
for d = 1:2
    rdm={};
    for s=1:16
        x=load(sprintf('../results/sub-%02i_decoding_rdm_%ims.mat',s,durations(d)));
        rdm{s} = cosmo_average_samples(x.rdm,'split_by',{'target1','target2'});
        stimtable = x.stimtable;
        stimtable.f_colour = stimtable.f_color;
    end
    rdm_cell{d} = cosmo_average_samples(cosmo_stack(rdm),'split_by',{'target1','target2'});
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
if exists('../results/stats_behaviour.mat')
    load('../results/stats_behaviour.mat','RDMmean')
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
    save('../results/stats_behaviour.mat','RDMmean')
end

%% create embedding
fprintf('create embedding\n')
Y = RDMmean;
Y(eye(size(Y))==1)=0;
X1 = mdscale(Y,2,'Start','random','Criterion','metricstress','Replicates',10);

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

%% make figure
f=figure(1);clf
f.Position = [f.Position(1:2) 800 1000];
f.Resize='off';
nsub = numel(unique(T.subjectnr));

% behavioural results
a=subplot(3,2,1);
imagesc(RDMmean);axis square
a.YDir='normal';
a.FontSize=12;
p = a.Position;
c=colorbar;
a.Position = p;
c.Label.String='dissimilarity';
title(sprintf('Behavioural RDM (n=%i)',nsub))
colormap inferno
xlabel('stimulus')
ylabel('stimulus')

% mds 
a=subplot(3,2,2);hold on
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

%model correlations
a=subplot(3,2,3);
R = corr(modelvec,'type','spearman');
R(eye(size(R))==1)=nan;
imagesc(R)
axis square
hold on
a.FontSize=14;
colormap inferno;
a.YDir = 'normal';
%c=colorbar;c.Label.String='model correlation';
a.XTick = 1:5;
a.YTick = 1:5;

co = tab10;
for m=1:5
    cmodelnames{m} = ['\bf{\',sprintf('color[rgb]{%f %f %f}%s}',co(m,1),co(m,2),co(m,3),modelnames{m})];
end

%cmodelnames = {'\bf{\color[rgb]{1 1 1}ori}','\color{red}sf','\color{red}colour','\color{red}contrast','behaviour'};
a.XTickLabel = cmodelnames;
a.YTickLabel = cmodelnames;
for x=1:5
    for y=1:5
        if (x==5 && y<5) | (y==5 && x<5) | 1
            r = R(x,y);
            if r>0 || r<1
                cc = flipud(inferno(200));
                c = cc(find(r<=linspace(a.CLim(1),a.CLim(2),size(cc,1)),1),:);
                text(x,y,sprintf('%.3f',r),'Color',c,'HorizontalAlignment','center')
            end
        end
    end
end

%model mds
a=subplot(3,2,4);
a.FontSize=16;
rng(1)
Y = mdscale(squareform(pdist(modelvec','spearman')),2);
aw = .9*[-1 1]*max(range(minmax(Y')));
for i=1:5
    x = Y(i,1);
    y = Y(i,2);
    plot(x,y,'o','Color',co(i,:));hold on
    axis equal
    image('XData',x+aw,'YData',y+aw,'CData',squareform(255*normalize(modelvec(:,i),'range')))
    text(x,max(y+aw),modelnames{i},'Color',co(i,:),'FontSize',16,'VerticalAlignment','bottom','HorizontalAlignment','center')
end
axis off

% corr with neural data
a=subplot(3,2,5);

% models
rmodel = [];
for d=1:2
    rdm = rdm_cell{d};
    tv = rdm.a.fdim.values{1};
    X = nan(175,256,256);
    for r = 1:length(rdm.sa.target1)
        t1 = 1+rdm.sa.target1(r);
        t2 = 1+rdm.sa.target2(r);
        X(:,t1,t2) = rdm.samples(r,:);
        X(:,t2,t1) = rdm.samples(r,:);
    end
    loweridx = find(tril(ones(256),-1)==1);
    for m=1:5
        rmodel(:,m,d) = corr(X(:,loweridx)', modelvec(:,m),'type','Spearman');
    end
end

% plot
plotnr = 4;
for d=1:2
    plotnr=plotnr+1;
    a = subplot(3,2,plotnr);hold on
    a.FontSize=12;a.ColorOrder=co;
    rdm = rdm_cell{d};
    tv = rdm.a.fdim.values{1};
    X = nan(175,256,256);
    for r = 1:length(rdm.sa.target1)
        t1 = 1+rdm.sa.target1(r);
        t2 = 1+rdm.sa.target2(r);
        X(:,t1,t2) = rdm.samples(r,:);
        X(:,t2,t1) = rdm.samples(r,:);
    end
    loweridx = find(tril(ones(256),-1)==1);    
    plot(tv,rmodel(:,:,d),'LineWidth',2)
    legend(modelnames)
    ylim([-.05 0.4])
    xlim([-100 600])
    xlabel('time (ms)')

    title(sprintf('%.2f Hz',1000/durations(d)))
    ylabel('model correlation with EEG')
    drawnow
end



%% save
fn = '../figures/rsa';
tn = tempname;
print(gcf,'-dpng','-r500',tn)
im=imread([tn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');

