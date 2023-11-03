function stats_decoding_conjunction()

% analyses for conjunction of features in fast features experiment:
% - ori x SF
% - ori x contrast
% - ori x col
% - sf x contrast
% - sf x col
% - contrast x col
% calculate group (N=16) means, se, bayes factors for each conjunction
% above chance (50%)
% calculate onsets and peaks of decoding using leave-two-participants-out jackknife

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
stats.timevect = res_all.a.fdim.values{1};

dur = [.15 .05];
h0mean = .5;

for d = 1:length(dur)
    for c = 1:length(unique(res_all.sa.featurecombonum))
        
        %% get group means and stats

        idx = res_all.sa.featurecombonum==c&res_all.sa.soaduration==dur(d);
        x = res_all.samples(idx,:);
        s = struct();
        s.soaduration = res_all.sa.soaduration(find(idx,1));
        s.targetname = res_all.sa.featurenames{find(idx,1)};
        s.n = size(x,1);
        s.mu = mean(x);
        s.mu_all = x;
        s.se = std(x)./sqrt(s.n);
        
        % calculate bayesfactors
        s.bf = bayesfactor_R_wrapper(x',...
            'returnindex',2,'verbose',false,'args','mu=0.5,rscale="medium",nullInterval=c(-Inf,0.5)');

        s.p_uncor = arrayfun(@(y) signrank(x(:,y),h0mean,'tail','right'),1:size(x,2));
        [~,~,s.fdr_adj_p]=fdr(s.p_uncor);

        stats.(sprintf('soa%i',round(1000*s.soaduration))).(s.targetname) = s;

        %% calculate onsets and peaks
        fprintf('Calculating onsets for %s\n',s.targetname)
        
        % group mean onset & peak
        [xo,ot] = max(movmean(s.bf>10,[0 2])==1); % onset 3 consecutive tp BF>10
        [~,peak] = max(s.mu); % peak

        % jackknife method on bayesfactors to get distribution of BFs
        x = [];
        combs = combnk(1:16,14); %leave-2-out
        for i=1:size(combs,1) % concatenate all data for jackknife
            x = [x s.mu_all(combs(i,:),:)];
        end
        b = bayesfactor_R_wrapper(x',...
            'returnindex',2,'verbose',false,'args','mu=0.5,rscale="medium",nullInterval=c(-Inf,0.5)');
        s.bf_jackknife = reshape(b,[],size(combs,1)); 
        
        % get onsets and peaks per jackknife
        onset_jack = [];peak_jack=[];
        xx=0;
        for b=1:size(s.bf_jackknife,2) % for each jackknife permutation
            
            bfj = s.bf_jackknife(:,b); % bfs for that perm
            
            % get onsets (first three consecutive time points BF>10)
            exc = movmean(bfj>10,[0 2])==1;
            if sum(exc)>0 % is there is an onset
                [~,t] = max(exc);
                xx=xx+1;
                onset_jack(xx) = t;
            end

            % calculate peaks per permutation
            [~,peak_jack(b)] = max(mean(s.mu_all(combs(b,:),:)));

        end

        % convert onsets and peaks to time (rather than tp)
        if xo == 0 % no reliable onset
            s.onset = NaN;
        else
            s.onset = stats.timevect(ot);
        end
        s.onsetjk = stats.timevect(onset_jack);
        if isempty(onset_jack)
            s.onsetci=[];
        else
            ci = prctile(onset_jack,[5 95]);
            s.onsetci = stats.timevect([floor(ci(1)) ceil(ci(2))]);
        end
        s.peak = stats.timevect(peak);
        s.peakjk = stats.timevect(peak_jack);
        s.peakci = stats.timevect(round(prctile(peak_jack,[5 95])));
        
        stats.(sprintf('soa%i',round(1000*s.soaduration))).(s.targetname) = s;
        fprintf('Saving\n')
        save('results/stats_decoding_conjunction.mat','stats');
        
    end
end
stats.format = 'stats.soa.targetcombo.x';

%%
fprintf('Saving\n')
save('results/stats_decoding_conjunction.mat','stats');
fprintf('Done\n')