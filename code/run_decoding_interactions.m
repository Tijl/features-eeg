function run_decoding_interactions(subjectnr)

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

%% decoding pilot
fn = sprintf('data/derivatives/cosmomvpa/sub-%02i_task-rsvp_cosmomvpa.mat',subjectnr);
fprintf('loading %s\n',fn);tic
load(fn,'ds')
fprintf('loading data finished in %i seconds\n',ceil(toc))
outfn = sprintf('results/sub-%02i_decoding_interactions.mat',subjectnr);
dsbackup = ds;

%% decode
durations = [0.15 0.05];
res_cell={};
for d = 1:2
    ds = cosmo_slice(dsbackup,dsbackup.sa.soaduration==durations(d));
    ds.sa.chunks = ds.sa.sequencenumber;
    nh = cosmo_interval_neighborhood(ds,'time','radius',0);
    measure = @cosmo_crossvalidation_measure;
    ma = {};
    ma.classifier = @cosmo_classify_lda;
    ma.nproc = nproc;

    targetnames = {'ori','sf','color','contrast'};

    for tf=1:4 % target to decode
        ds.sa.targets = ds.sa.(['f_' targetnames{tf}]);


        for fu = 1:4 % feature under which to decode
            if tf == fu
                continue
            end
            fprintf('Decoding: %ims %s by %s\n',1000*durations(d),targetnames{tf},targetnames{fu})

            for l = 0:3 %level of fu to decode under
                dsd = cosmo_slice(ds,ds.sa.(sprintf('f_%s',targetnames{fu}))==l);

                ma.partitions = cosmo_nfold_partitioner(dsd);

                res = cosmo_searchlight(dsd,nh,measure,ma);
                res.sa.subject = subjectnr;
                res.sa.soaduration = dsd.sa.soaduration(1);
                res.sa.targetname = targetnames(tf);
                res.sa.targetfeature = tf;
                res.sa.constantname = targetnames(fu);
                res.sa.constantfeature = fu;
                res.sa.constantfeatlevel = l;
                res_cell{end+1} = res;
            end
        end
    end
end
res = cosmo_stack(res_cell);
save(outfn,'res','-v7.3')
