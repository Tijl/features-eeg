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
% ma.progress=0;

%% resume stopped output files
outfn = sprintf('results/sub-%02i_decoding_2wayconjunction.mat',subjectnr);

if exist(outfn,'file')==2 % already a file
    load(outfn,'res')
    n = size(res.samples,1); % completed decoding contrasts
    if n == 72*12 % all done
        return
    end
    ncombos = floor(n/72); % completed number of combinations
    startdur = floor(ncombos/6)+1;  % 1 or 2
    starttf =  mod(ncombos,6)+1; % 1-6
    startlev = mod(n,72)+1;
else
    res = {};
    startdur = 1;
    starttf = 1;
    startlev = 1;
end
origres = res;

%% get feature combinations
targetnames = {'ori','sf','color','contrast'};
fc = combnk(1:4,2);

% 2 way combinations of levels eg orientation 1-4 combined with contrast
% 1-4 as class 1 for decoding, and reverse as class 2
x=0;lcombs=zeros(4,1);lfeats=[];
for f = 1:4
    for b = 1:4

        class1A = f*10+b;
        
        for f2 = 1:4
            for b2 = 1:4
                class1B = f2*10+b2;
        
                class2A = f*10+b2;
                class2B = f2*10+b;
                
                if f == f2 || b == b2 % need to vary over both features
                    continue
                elseif sum(ismember(find(lcombs(1,:)==class2A),find(lcombs(2,:)==class2B)))>0 % combination already in array
                    continue
                end
                
                x=x+1;
                lcombs(:,x) = [class1A class1B class2A class2B];
                lfeats(:,:,x) = [f f2; b b2]; % lfeats(1,:,x) = feature 1 levels, lfeats(2,:,x) = feature 2 levels
            end
        end
    end
end

%% decode
durations = [0.15 0.05];
res_cell{1}=origres;
for d = startdur:length(durations)
    ds = cosmo_slice(dsbackup,dsbackup.sa.soaduration==durations(d));
    ds.sa.chunks = ds.sa.sequencenumber;
    nh = cosmo_interval_neighborhood(ds,'time','radius',0);

    for tf=starttf:length(fc) % combo of feat1 and feat2

        fprintf('Decoding: %d: %s and %s conjunction\n',1000*durations(d),targetnames{fc(tf,1)},targetnames{fc(tf,2)})
        
        
        for levcomb = startlev:size(lfeats,3) % all combinations of levels for feat1 and feat2
            
            fprintf('feature combination %d out of %d\n',levcomb,size(lfeats,3));
            
            % class 1
            clear idx1a idx1b
            for fnum = 1:2 % for each feature
                idx1a(:,fnum) = ds.sa.(sprintf('f_%s',targetnames{fc(tf,fnum)}))==lfeats(fnum,1,levcomb)-1; % subtypeA - level for that feature
                idx1b(:,fnum) = ds.sa.(sprintf('f_%s',targetnames{fc(tf,fnum)}))==lfeats(fnum,2,levcomb)-1; % subtypeB - other level for that feature
            end
            
            dsd1 = cosmo_slice(ds,(idx1a(:,1)&idx1a(:,2)) | (idx1b(:,1)&idx1b(:,2)));
            dsd1.sa.targets = ones(size(dsd1.sa.f_color));
            dsd1.sa.combos = (dsd1.sa.(sprintf('f_%s',targetnames{fc(tf,1)}))+1)*10+dsd1.sa.(sprintf('f_%s',targetnames{fc(tf,2)}))+1;
            combo1 = unique(dsd1.sa.combos);
            
            % class 2 = other feat level combinations
            dsd2 = cosmo_slice(ds,(idx1a(:,1)&idx1b(:,2)) | (idx1b(:,1)&idx1a(:,2)));
            dsd2.sa.targets = ones(size(dsd1.sa.f_color))*2;
            dsd2.sa.combos = (dsd2.sa.(sprintf('f_%s',targetnames{fc(tf,1)}))+1)*10+dsd2.sa.(sprintf('f_%s',targetnames{fc(tf,2)}))+1;
            combo2 = unique(dsd2.sa.combos);
            
            dsd = cosmo_stack({dsd1,dsd2});
            
            ma.partitions = cosmo_nfold_partitioner(dsd);
            
            res = cosmo_searchlight(dsd,nh,measure,ma);
            res.sa.subject = subjectnr;
            res.sa.soaduration = dsd.sa.soaduration(1);
            res.sa.featurecombonum = tf;
            res.sa.featurecombo = fc(tf,1)*10+fc(tf,2);
            res.sa.feature1 = targetnames(fc(tf,1));
            res.sa.feature2 = targetnames(fc(tf,2));
            res.sa.featurenames = {sprintf('%s_%s',res.sa.feature1{:},res.sa.feature2{:})};
            res.sa.levcomb = levcomb;
            res.sa.levelcombos = {sprintf('%dand%dvs%dand%d',combo1(1),combo1(2),combo2(1),combo2(2))};

            res_cell{end+1} = res;
        end
        res = cosmo_stack(res_cell);
        save(outfn,'res','-v7.3')
    end
end

