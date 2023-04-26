function run_decoding_conjunctions(subjectnr)

%%
if ismac
    if isempty(which('cosmo_wtf'))
        addpath('~/CoSMoMVPA/mvpa')
    end
    nproc = 2;
else %on HPC
    addpath('../CoSMoMVPA/mvpa');
    % start cluster, give it a unique directory
    % starting a pool can fail when 2 procs are requesting simultaneous
    % thus try again after a second until success
    pool=[];
    while isempty(pool)
        try
            pc = parcluster('local');
            pc.JobStorageLocation=tempdir;
            pool=parpool(pc);
        catch err
            disp(err)
            delete(gcp('nocreate'));
            pause(1)
        end
    end
    nproc=cosmo_parallel_get_nproc_available();
end

%% set up output file
outfn = sprintf('results/sub-%02i_decoding_2wayconjunction.mat',subjectnr);

if exist(outfn,'file')==2 % already a file
    return
end

%% load data
fn = sprintf('data/derivatives/cosmomvpa/sub-%02i_task-rsvp_cosmomvpa.mat',subjectnr);
fprintf('loading %s\n',fn);tic
ds = matfile(fn);
ds = ds.ds;
fprintf('loading data finished in %i seconds\n',ceil(toc))
dsbackup = ds;

%% decoding details
measure = @cosmo_crossvalidation_measure;
ma = {};
ma.classifier = @cosmo_classify_lda;
ma.nproc = nproc;
ma.progress=0;

%% get feature combinations
targetnames = {'ori','sf','color','contrast'};
fc = combnk(1:4,2);

% 2 way combinations of levels eg orientation 1-4 combined with contrast
% 1-4 as class 1 for decoding, and reverse as class 2

%% two way combinations
% for two way feature combinations, use 2 stimulus types per class. 
% use 2 levels per each of 2 features

% THESE ARE THE SPECIFIC FEATURE LEVELS TO USE PER FEATURE
% 2 levels x 2 features x all possible combinations
featlevels = zeros(2,2,36);% 2 feature levels x 2 features x 36 combinations 
x=0;
combs = combnk(0:3,2); % choice of two from 4 feature levels
for f1 = 1:length(combs) 
    for f2 = 1:length(combs)
        x=x+1;
        featlevels(:,1,x) = combs(f1,:); % feature1
        featlevels(:,2,x) = combs(f2,:); % feature2
    end
end

% THIS IS THE GENERAL COMBINATION STRATEGY
% get indices for which stimulus feature goes into each decoding class
% combine features for two classes that are equal on single features
% but differ on the 2-way interactions
lfeats = zeros(2,2,2); % stimnumber (2 per class) x feature type (2 of ori,sf,col,con) x class (decoding class 1 vs 2)
lfeats(:,1,1) = [1 2]; % feature1, class 1
lfeats(:,1,2) = [2 1]; % feature1, class 2
lfeats(:,2,1) = [1 2]; % feature2, class 1
lfeats(:,2,2) = [1 2]; % feature2, class 2

%% decode
durations = [0.15 0.05];
res_cell={};
for d = 1:length(durations)
    ds = cosmo_slice(dsbackup,dsbackup.sa.soaduration==durations(d));
    ds.sa.chunks = ds.sa.sequencenumber;
    nh = cosmo_interval_neighborhood(ds,'time','radius',0);

    for fcnum=1:length(fc) % which two features to combine
        feats = targetnames(fc(fcnum,:));
        fprintf('Decoding: %d: %s and %s conjunction\n',1000*durations(d),targetnames{fc(fcnum,1)},targetnames{fc(fcnum,2)})
        
        for levcomb = 1:size(featlevels,3) % which combination of features
            fprintf('feature combination %d out of %d\n',levcomb,size(featlevels,3));
            
            fl = squeeze(featlevels(:,:,levcomb)); % level x feature type. which 2 specific ori/sf/col/contrast to use
            ds.sa.targets = zeros(size(ds.sa.f_ori));

            % set class targets
            for i = 1:size(lfeats,1) % for each stim in that class
                for n = 1:size(lfeats,3) % two classes
                    % get indices for stim each feature
                    classidx = lfeats(i,:,n); % which level number for each feature according to decoding class
                    idx = ds.sa.(sprintf('f_%s',feats{1}))==fl(classidx(1),1)&...
                        ds.sa.(sprintf('f_%s',feats{2}))==fl(classidx(2),2); % subtypeA - level for that feature
                    ds.sa.targets(idx) = n; % class num
                end
            end
            
            dsd = cosmo_slice(ds,ds.sa.targets>0);
            dsd.sa.combos = (dsd.sa.(sprintf('f_%s',feats{1}))+1)*10+dsd.sa.(sprintf('f_%s',feats{2}))+1;
            combo1 = unique(dsd.sa.combos);
            
            ma.partitions = cosmo_nfold_partitioner(dsd);
            
            res = cosmo_searchlight(dsd,nh,measure,ma);
            res.sa.subject = subjectnr;
            res.sa.soaduration = dsd.sa.soaduration(1);
            res.sa.featurecombonum = fcnum;
            res.sa.featurecombo = fc(fcnum,1)*10+fc(fcnum,2);
            res.sa.feature1 = feats(1);
            res.sa.feature2 = feats(2);
            res.sa.featurenames = {sprintf('%s_%s',feats{1},feats{2})};
            res.sa.levcomb = levcomb;
            res.sa.levelcombos = {sprintf('%d_%d_%d_%d',combo1(1),combo1(2),combo1(3),combo1(4))};

            res_cell{end+1} = res;
        end
        save(outfn,'res_cell','-v7.3')
    end
end

res = cosmo_stack(res_cell);
save(outfn,'res','-v7.3')

