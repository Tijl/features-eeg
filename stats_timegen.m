function stats_timegen()

% analyses for individual features in fast features experiment for time x time generalisation analyses:
% features: orientation, spatial frequency, contrast and colour
% calculate group (N=16) means, se, bayes factors for each feature 
% test for above/below chance for every train/test timepoint combination

if isempty(which('cosmo_wtf'))
    addpath('~/CoSMoMVPA/mvpa')
end

%% stack results

resfiles = dir('results/*decoding_timegen.mat');
features = {'ori' 'sf' 'color' 'contrast'};
duration = [0.05 0.15];
fname = 'results/stats_timegen.mat';

% concatenate different participants' data
fprintf('Loading data\n')
cc = clock();mm='';
res_cell={};
for p = 1:length(resfiles)
    m1 = matfile(sprintf('%s/%s',resfiles(p).folder,resfiles(p).name));
    dat = m1.res_cell;
    res_cell{end+1} = cosmo_stack(dat);
    mm = cosmo_show_progress(cc,p/length(resfiles),sprintf('%i/%i',p,length(resfiles)),mm);
end
res_all = cosmo_stack(res_cell);
timevect = dat{1}.a.sdim.values{1};
n=length(resfiles);

%% stats

fprintf('Computing stats\n')
stats = struct();
stats.timevect = timevect;
cc = clock();mm='';

for d = 1:length(duration)
    for f = 1:length(features)
        
        fdat = cosmo_slice(res_all,res_all.sa.targetfeature==f&res_all.sa.soaduration==duration(d));

        x = reshape(fdat.samples,[],n)';
        s = struct();
        s.soaduration = unique(fdat.sa.soaduration);
        s.targetname = unique(fdat.sa.targetname);
        s.traintime = fdat.sa.train_time(1:size(x,2));
        s.testtime = fdat.sa.test_time(1:size(x,2));
        s.n = n;
        s.mu = mean(x)';
        s.mu_all = x';
        s.se = (std(x)./sqrt(s.n))';
        
        % calculate bayesfactors
        s.bf = bayesfactor_R_wrapper(x',...
            'returnindex',2,'verbose',false,'args','mu=0.25,rscale="medium",nullInterval=c(-.5,0.5)');
        
        stats.(sprintf('soa%d',duration(d)*1000)).(features{f}) = s;
        mm = cosmo_show_progress(cc,((d-1)*4+f)/8,sprintf('%i/8',((d-1)*4+f)),mm);

    end
end

stats.format = 'stats.soa.targetfeature.x';

%%
fprintf('Saving\n')
save(fname,'stats');
fprintf('Done\n')