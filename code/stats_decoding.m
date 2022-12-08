function stats_decoding()

    if isempty(which('cosmo_wtf'))
        addpath('~/CoSMoMVPA/mvpa')
    end
     
    %% stack results
    fprintf('Loading data\n')
    files = dir('results/sub-*_decoding.mat');
    res_cell={};
    cc = clock();mm='';
    for f=1:length(files)
        fn = sprintf('results/%s',files(f).name);
        load(fn,'res');
        res.sa.contrast = (1:8)';
        res_cell{f} = res;
        mm = cosmo_show_progress(cc,f/length(files),sprintf('%i/%i',f,length(files)),mm);
    end
    res_all = cosmo_stack(res_cell);
    
    %% stats
    fprintf('Computing stats\n')
    stats = struct();
    for c = 1:8
        idx = res_all.sa.contrast==c;
        x = res_all.samples(idx,:);
        s = struct();
        s.soaduration = res_all.sa.soaduration(find(idx,1));
        s.targetname = res_all.sa.targetname{find(idx,1)};
        s.n = size(x,1);
        s.mu = mean(x);
        s.mu_all = x;
        s.se = std(x)./sqrt(s.n);
        h0mean = .25;
        s.tstat = (s.mu-h0mean)./s.se;
        % calculate bayesfactors
        s.bf = bayesfactor_R_wrapper(x',...
            'returnindex',2,'verbose',false,'args','mu=0.25,rscale="medium",nullInterval=c(-Inf,0.5)');

        s.p_uncor = arrayfun(@(y) signrank(x(:,y),h0mean,'tail','right'),1:size(x,2));
        [~,~,s.fdr_adj_p]=fdr(s.p_uncor);
        s.tv = res_all.a.fdim.values{1};
        
        stats.(sprintf('soa%i',round(1000*s.soaduration))).(s.targetname) = s;
    end
    
    stats.format = 'stats.soa.targetfeature.x';
    stats.timevect = res_all.a.fdim.values{1};
    
    %%
    fprintf('Saving\n')
    save('results/stats_decoding.mat','stats');
    fprintf('Done\n')