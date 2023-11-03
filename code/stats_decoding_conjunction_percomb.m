function stats_decoding_conjunction_percomb()

% analyses for conjunction of features in fast features experiment:
% for each of 36 different contrasts
% - ori x SF
% - ori x contrast
% - ori x col
% - sf x contrast
% - sf x col
% - contrast x col
% calculate group (N=16) means, se, bayes factors for each conjunction
% contrast above chance (50%)

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

%% stats
fprintf('Computing stats\n')
stats = struct();
stats.timevect = res_all.a.fdim.values{1};

dur = [.15 .05];
h0mean = .5;

for d = 1:length(dur)
    for c = 1:length(unique(res_all.sa.featurecombonum))
        for n = 1:length(unique(res_all.sa.levcomb))

            %% get group means and stats

            idx = res_all.sa.featurecombonum==c&res_all.sa.soaduration==dur(d)&res_all.sa.levcomb==n;
            x = res_all.samples(idx,:);
            s = struct();
            s.soaduration = res_all.sa.soaduration(find(idx,1));
            s.targetname = res_all.sa.featurenames{find(idx,1)};
            s.levcomb = n;
            s.n = size(x,1);
            s.mu = mean(x);
            s.mu_all = x;
            s.se = std(x)./sqrt(s.n);

            % calculate bayesfactors
            s.bf = bayesfactor_R_wrapper(x',...
                'returnindex',2,'verbose',false,'args','mu=0.5,rscale="medium",nullInterval=c(-Inf,0.5)');

            stats.(sprintf('soa%i',round(1000*s.soaduration))).(s.targetname).(sprintf('comb%d',s.levcomb)) = s;

            fprintf('Saving\n')
            save('results/stats_decoding_conjunction_eachcontrast.mat','stats');
        end
    end
end
stats.format = 'stats.soa.targetcombo.x';

%%
fprintf('Saving\n')
save('results/stats_decoding_conjunction_eachcontrast.mat','stats');
fprintf('Done\n')