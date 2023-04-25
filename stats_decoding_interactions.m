function stats_decoding_interactions()

if isempty(which('cosmo_wtf'))
    addpath('~/CoSMoMVPA/mvpa')
end

%% stack results
fprintf('Loading data\n')
files = dir('results/sub-*_decoding_interactions.mat');
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

%% stats above chance per contrast
fprintf('Computing stats\n')
stats = struct();
for c = 1:length(unique(res.sa.contrast))
    idx = res_all.sa.contrast==c;
    x = res_all.samples(idx,:);
    s = struct();
    s.soaduration = res_all.sa.soaduration(find(idx,1));
    s.targetname = res_all.sa.targetname{find(idx,1)};
    s.constantname = res_all.sa.constantname{find(idx,1)};
    s.constantlevel = res_all.sa.constantfeatlevel(find(idx,1));
    s.n = size(x,1);
    s.mu = mean(x);
    s.mu_all = x;
    s.se = std(x)./sqrt(s.n);
    h0mean = .25;
    % calculate bayesfactors
    s.bf = bayesfactor_R_wrapper(x',...
        'returnindex',2,'verbose',false,'args','mu=0.5,rscale="medium",nullInterval=c(-Inf,0.5)');

    s.p_uncor = arrayfun(@(y) signrank(x(:,y),h0mean,'tail','right'),1:size(x,2));
    [~,~,s.fdr_adj_p]=fdr(s.p_uncor);
    s.tv = res_all.a.fdim.values{1};

    stats.(sprintf('soa%i',round(1000*s.soaduration))).(s.targetname).(s.constantname).(sprintf('level%d',s.constantlevel+1)) = s;
end

stats.format = 'stats.soa.targetfeature.constantfeature.constantlevel.x';
stats.timevect = res_all.a.fdim.values{1};

%% stats for comparison across levels at peak
durations = {'soa150','soa50'};
targetlabels = {'ori','sf','color','contrast'};
peaktimes = [120 112 112 96; 120 112 112 100]; % overall peak for ori, sf, col and con - use for calculating diffs per level
comblev = combnk(1:4,2);

for d=1:2 % each SOA

    for feat_dec = 1:length(targetlabels) % decoding feature

        timestoplot = find(stats.timevect==peaktimes(d,feat_dec));
        timestoplot = (timestoplot-2):(timestoplot+2);

        for feat_by = 1:length(targetlabels) % for each level of other feature

            if feat_dec == feat_by % no such thing

                continue

            else % statistically compare different levels of feat_by
                fprintf('Computing stats for %s by level of %s\n',targetlabels{feat_dec},targetlabels{feat_by})
                stats.(durations{d}).(targetlabels{feat_dec}).(targetlabels{feat_by}).diffbf = zeros(4,4);
                
                for c = 1:size(comblev,1)

                    dat1 = stats.(durations{d}).(targetlabels{feat_dec}).(targetlabels{feat_by}).(sprintf('level%d',comblev(c,1))).mu_all;
                    dat2 = stats.(durations{d}).(targetlabels{feat_dec}).(targetlabels{feat_by}).(sprintf('level%d',comblev(c,2))).mu_all;

                    % calculate mean and se for that time window
                    dat_mu1 = mean(dat1(:,timestoplot),2);
                    dat_mu2 = mean(dat2(:,timestoplot),2);
                    dat_diff = dat_mu1-dat_mu2;

                    % calculate bayesfactors
                    s.bf = bayesfactor_R_wrapper(dat_diff',...
                        'returnindex',2,'verbose',false,'args','mu=0,rscale="medium",nullInterval=c(-0.5,0.5)');

                    stats.(durations{d}).(targetlabels{feat_dec}).(targetlabels{feat_by}).diffbf(comblev(c,1),comblev(c,2)) = s.bf;

                end
            end
        end
    end
end

%%
fprintf('Saving\n')
save('results/stats_decoding_interactions.mat','stats');
fprintf('Done\n')