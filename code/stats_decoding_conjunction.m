function stats_decoding_conjunction()

    if isempty(which('cosmo_wtf'))
        addpath('~/CoSMoMVPA/mvpa')
    end
        
    %% stack results
    fprintf('Loading data\n')
    files = dir('results/sub-*_decoding_2wayconjunction.mat');
    res_cell={};
    cc = clock();mm='';
    for f=1:length(files)
        fn = sprintf('results/%s',files(f).name);
        load(fn,'res');
        res.sa.contrast = (1:length(res.sa.subject))';
        res_cell{f} = res;
        mm = cosmo_show_progress(cc,f/length(files),sprintf('%i/%i',f,length(files)),mm);
    end
    res_all = cosmo_stack(res_cell);
    res_all = cosmo_average_samples(res_all,'split_by',{'subject', 'soaduration','featurecombo'});
    
    %% stats
    fprintf('Computing stats\n')
    stats = struct();
    dur = [.15 .05];

    for d = 1:length(dur)
        for c = 1:length(unique(res_all.sa.featurecombonum))
            idx = res_all.sa.featurecombonum==c&res_all.sa.soaduration==dur(d);
            x = res_all.samples(idx,:);
            s = struct();
            s.soaduration = res_all.sa.soaduration(find(idx,1));
            s.targetname = res_all.sa.featurenames{find(idx,1)};
            s.n = size(x,1);
            s.mu = mean(x);
            s.mu_all = x;
            s.se = std(x)./sqrt(s.n);
            h0mean = .5;

            % calculate bayesfactors
            s.bf = bayesfactor_R_wrapper(x',...
                'returnindex',2,'verbose',false,'args','mu=0.5,rscale="medium",nullInterval=c(-Inf,0.5)');

            s.p_uncor = arrayfun(@(y) signrank(x(:,y),h0mean,'tail','right'),1:size(x,2));
            [~,~,s.fdr_adj_p]=fdr(s.p_uncor);
            s.tv = res_all.a.fdim.values{1};

            stats.(sprintf('soa%i',round(1000*s.soaduration))).(s.targetname) = s;
        end
    end
    stats.format = 'stats.soa.targetcombo.x';
    stats.timevect = res_all.a.fdim.values{1};
    
    %%
    fprintf('Saving\n')
    save('results/stats_decoding_conjunction.mat','stats');
    fprintf('Done\n')