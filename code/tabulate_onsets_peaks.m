%% tabulate onsets and peaks

durations = {'soa150','soa50'};h=[];
targetlabels = {'ori','sf','color','contrast'};
combs = combnk(1:16,14); %leave-2-out

allinfo=struct();

%% get onsets and peaks for individual feature decoding
load('features-eeg/results/stats_decoding.mat')


x=0;
for d = 1:2
    for t=1:4
        s = stats.(durations{d}).(targetlabels{t});
        
        x=x+1;
        allinfo.SOA{x,1} = durations{d};
        allinfo.Condition{x,1} = targetlabels{t};
        allinfo.Onset(x,1) = s.onset;
        if ~isempty(s.onsetci)
            allinfo.Onset_lowCI(x,1) = s.onsetci(1);
            allinfo.Onset_highCI(x,1) = s.onsetci(2);
        end
        allinfo.Peak(x,1) = s.peak;
        allinfo.Peak_lowCI(x,1) = s.peakci(1);
        allinfo.Peak_highCI(x,1) = s.peakci(2);
        
    end
end

%% get onsets and peaks for feature conjunction decoding
load('results/stats_decoding_conjunction.mat')

feats = fieldnames(stats.soa150);

for d = 1:2

    for t=1:length(feats)
        s = stats.(durations{d}).(feats{t});
        
        x=x+1;
        allinfo.SOA{x,1} = durations{d};
        allinfo.Condition{x,1} = feats{t};
        allinfo.Onset(x,1) = s.onset;
        if ~isempty(s.onsetci)
            allinfo.Onset_lowCI(x,1) = s.onsetci(1);
            allinfo.Onset_highCI(x,1) = s.onsetci(2);
        end
        allinfo.Peak(x,1) = s.peak;
        allinfo.Peak_lowCI(x,1) = s.peakci(1);
        allinfo.Peak_highCI(x,1) = s.peakci(2);
    end
end

allinfo = struct2table(allinfo)
