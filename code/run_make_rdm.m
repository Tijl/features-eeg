function run_make_rdm(subjectnr)
    
    %%
    if ismac
        if isempty(which('cosmo_wtf'))
            addpath('~/CoSMoMVPA/mvpa')
        end
        nproc = 1;
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
    
    %%
    fn = sprintf('data/derivatives/cosmomvpa/sub-%02i_task-rsvp_cosmomvpa.mat',subjectnr);
    fprintf('loading %s\n',fn);tic
    load(fn,'ds')
    fprintf('loading data finished in %i seconds\n',ceil(toc))
    ds_splits = cosmo_split(ds,'soaduration');
    ds=[];    
    
    %% 
    for split = 1:numel(ds_splits)
        ds = ds_splits{split};        
        fprintf('Making %.2fHz RDM\n',1/ds.sa.soaduration(1))
        outfn = sprintf('results/sub-%02i_decoding_rdm_%ims.mat',subjectnr,1000*ds.sa.soaduration(1));
        ds.sa.chunks = mod(ds.sa.sequencenumber,2);
        ds.sa.targets = ds.sa.stimnumber;
        % store info about the stimuli
        x = cosmo_average_samples(cosmo_slice(ds,ds.sa.blocksequencenumber==1),'split_by',{'stimnumber'});
        stimtable = struct2table(x.sa);
        % find locations of stimuli and chunks
        ut = unique(ds.sa.targets);
        uc = unique(ds.sa.chunks);
        combs = combnk(1:length(ut),2);
        idx_k = [];
        for k = 1:length(uc)
            idx_k(:,k) = ds.sa.chunks==uc(k);
        end
        idx_e = [];
        for e = 1:length(ut)
            idx_e(:,e) = ds.sa.targets == ut(e);
        end
        % make partitions
        partitions = struct('train_indices',[],'test_indices',[]);
        sa = struct('target1',[],'target2',[],'lefoutchunk',[]);
        cc = clock(); mm = '';
        for c = 1:length(combs)
            idx_c = any(idx_e(:,combs(c,:)),2);
            for k = 1:length(uc)
                partitions.train_indices{1,end+1} = find(idx_c & ~idx_k(:,k));
                partitions.test_indices{1,end+1} = find(idx_c & idx_k(:,k));
                sa.target1(end+1,1) = ut(combs(c,1));
                sa.target2(end+1,1) = ut(combs(c,2));
                sa.lefoutchunk(end+1,1) = k;
            end
            if ~mod(c,100)
                mm = cosmo_show_progress(cc,c/length(combs),sprintf('%i/%i',c,length(combs)),mm);
            end
        end
        % set up decoding
        fprintf('\nDecoding all pairs\n')
        ma = {};
        ma.partitions = partitions;
        ma.nproc = nproc;
        ma.classifier = @cosmo_classify_lda;
        ma.check_partitions = 0;
        ma.output = 'fold_accuracy';
        measure = @cosmo_crossvalidation_measure;
        nh = cosmo_interval_neighborhood(ds,'time','radius',0);
        
        % run searchlight
        rdm = cosmo_searchlight(ds,nh,measure,ma);
        rdm.sa = cosmo_structjoin(rdm.sa,sa);
        rdm.sa.split = split * ones(size(rdm.samples,1),1);
        rdm.sa.soaduration = ds.sa.soaduration(1) * ones(size(rdm.samples,1),1);
        
        % save
        fprintf('saving %s\n',outfn);tic
        save(outfn,'rdm','stimtable','-v7.3')
        fprintf('saving data finished in %i seconds\n',ceil(toc))
    end
    