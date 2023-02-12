function stats_decoding()

% analyses for individual features in fast features experiment:
% features: orientation, spatial frequency, contrast and colour
% calculate group (N=16) means, se, bayes factors for each feature above chance
% calculate onsets and peaks of decoding using leave-two-participants-out jackknife

if isempty(which('cosmo_wtf'))
    addpath('~/CoSMoMVPA/mvpa')
end

%% stack results
fprintf('Loading data\n')
files = dir('../results/sub-*_decoding.mat');
res_cell={};
cc = clock();mm='';
for f=1:length(files)
    fn = sprintf('../results/%s',files(f).name);
    load(fn,'res');
    res.sa.contrast = (1:8)';
    res_cell{f} = res;
    mm = cosmo_show_progress(cc,f/length(files),sprintf('%i/%i',f,length(files)),mm);
end
res_all = cosmo_stack(res_cell);

%% stats
fprintf('Computing stats\n')
stats = struct();
stats.timevect = res_all.a.fdim.values{1};

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

    %% calculate onsets and peaks
    fprintf('Calculating onsets for %s\n',s.targetname)

    % group mean onset & peak
    [~,ot] = max(movmean(s.bf>6,[0 2])==1); % onset 3 consecutive tp BF>6
    [~,peak] = max(s.mu); % peak

    % use jackknife method and calculate bayesfactors to get onsets
    x = [];
    combs = combnk(1:16,14); %leave-2-out
    for i=1:size(combs,1)
        x = [x s.mu_all(combs(i,:),:)];
    end
    b = bayesfactor_R_wrapper(x',...
        'returnindex',2,'verbose',false,'args','mu=0.25,rscale="medium",nullInterval=c(-Inf,0.5)');
    s.bf_jackknife = reshape(b,[],size(combs,1)); %[~,onsets] = max(movmean(bf_jackknife>6,[0 2])==1);

    % get onsets and peaks per jackknife
    onset_jack = [];peak_jack=[];
    for b=1:size(s.bf_jackknife,2)

        % get onset (first three consecutive time points BF>6)
        [~,onset_jack(b)] = max(movmean(s.bf_jackknife(:,b)>6,[0 2])==1);

        % get peak
        [~,peak_jack(b)] = max(mean(s.mu_all(combs(b,:),:)));

    end

    % convert onsets and peaks to time (rather than tp)
    s.onset = stats.timevect(ot);
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
    fprintf('%i/8\n',c)
end
stats.format = 'stats.soa.targetfeature.x';
stats.timevect = res_all.a.fdim.values{1};

%%
fprintf('Saving\n')
save('../results/stats_decoding.mat','stats');
fprintf('Done\n')