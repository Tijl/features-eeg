function run_searchlight_decoding(subjectnr)

%%
if ismac
    if isempty(which('cosmo_wtf'))
        addpath('~/CoSMoMVPA/mvpa')
    end
    addpath('~/Dropbox (Personal)/MATLAB/fieldtrip-20171104/')
else %on HPC
    addpath('../CoSMoMVPA/mvpa');
    addpath('../fieldtrip')
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
end
nproc=cosmo_parallel_get_nproc_available();
ft_defaults;

cosmo_warning('off')

%% load files
fn = sprintf('data/derivatives/cosmomvpa/sub-%02i_task-rsvp_cosmomvpa.mat',subjectnr);
fprintf('loading %s\n',fn);tic
load(fn,'ds')
fprintf('loading data finished in %i seconds\n',ceil(toc))
outfn = sprintf('results/sub-%02i_decoding_searchlight.mat',subjectnr);
dsbackup = ds;

%% set up decoding
ma={};
ma.classifier = @cosmo_classify_lda;
ma.nproc = nproc;
ma.progress = 0;

measure = @cosmo_crossvalidation_measure;

% get channel neighbours
nh1 = cosmo_meeg_chan_neighborhood(ds,'count',4,'label','dataset','label_threshold',.99);
nh2 = cosmo_interval_neighborhood(ds,'time','radius',0);
nh = cosmo_cross_neighborhood(ds,{nh1,nh2});

%% decode
durations = [0.15 0.05];
res_cell={};
for d = 1:2
    ds = cosmo_slice(dsbackup,dsbackup.sa.soaduration==durations(d));
    ds.sa.chunks = ds.sa.sequencenumber;
    
    ma.partitions = cosmo_nfold_partitioner(ds);
    
    targetnames = {'ori','sf','color','contrast'};
    for tf=1:4
        fprintf('Decoding: %ims %s\n',1000*durations(d),targetnames{tf})
        ds.sa.targets = ds.sa.(['f_' targetnames{tf}]);
        res = cosmo_searchlight(ds,nh,measure,ma);
        res.sa.subject = subjectnr;
        res.sa.soaduration = ds.sa.soaduration(1);
        res.sa.targetname = targetnames(tf);
        res.sa.targetfeature = tf;
        res_cell{end+1} = res;
    end
end
res = cosmo_stack(res_cell);
save(outfn,'res','-v7.3')
