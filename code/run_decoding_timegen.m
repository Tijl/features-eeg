function run_decoding_timegen(subjectnr)

%% temporal generalisation for each feature

if ismac
    if isempty(which('cosmo_wtf'))
        addpath('~/CoSMoMVPA/mvpa')
    end
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
end
nproc=cosmo_parallel_get_nproc_available();

%% decoding pilot
fn = sprintf('data/derivatives/cosmomvpa/sub-%02i_task-rsvp_cosmomvpa.mat',subjectnr);
fprintf('loading %s\n',fn);tic
load(fn,'ds')
fprintf('loading data finished in %i seconds\n',ceil(toc))
outfn = sprintf('results/sub-%02i_decoding_timegen.mat',subjectnr);
dsbackup = ds;

%% decoding arguments
measure=@cosmo_dim_generalization_measure;

ma=struct();
ma.measure=@cosmo_crossvalidation_measure;
ma.classifier=@cosmo_classify_lda;
ma.dimension='time';
ma.nproc=nproc;
if ~ismac
    ma.progress=0;
end

%% decode
durations = [0.15 0.05];
targetnames = {'ori','sf','color','contrast'};
res_cell={};
for d = 1:2 % soa duration
    ds = cosmo_slice(dsbackup,dsbackup.sa.soaduration==durations(d));
    ds.sa.chunks = ds.sa.sequencenumber;
    
    for tf=1:4
        fprintf('Decoding: %ims %s\n',1000*durations(d),targetnames{tf})
        ds.sa.targets = ds.sa.(['f_' targetnames{tf}]);
        nchunks=2; % two chunks are required for this analysis
        ds.sa.chunks=cosmo_chunkize(ds,nchunks);

        % do timegen
        ds_tr=cosmo_dim_transpose(ds,'time',1);
        res=measure(ds_tr,ma);

        res.sa.subject = repmat(subjectnr,length(res.sa.labels),1);
        res.sa.soaduration = repmat(ds.sa.soaduration(1),length(res.sa.labels),1);
        res.sa.targetname = repmat(targetnames(tf),length(res.sa.labels),1);
        res.sa.targetfeature = repmat(tf,length(res.sa.labels),1);
        res_cell{end+1} = res;
    end
end
% res = cosmo_stack(res_cell);
save(outfn,'res_cell','-v7.3')
