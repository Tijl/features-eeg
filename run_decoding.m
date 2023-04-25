function run_decoding(subjectnr)
    
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
    load(fn,'ds')
    fprintf('loading data finished in %i seconds\n',ceil(toc))
    outfn = sprintf('results/sub-%02i_decoding.mat',subjectnr);
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
        ma.partitions = cosmo_nfold_partitioner(ds);
        ma.nproc = nproc;
        
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
    