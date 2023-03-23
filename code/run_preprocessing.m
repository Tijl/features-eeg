function run_preprocessing(partid)

    %% eeglab
    if ~ismac
        addpath('../CoSMoMVPA/mvpa')
        addpath('../eeglab')
    else
        addpath('~/CoSMoMVPA/mvpa')
        addpath('~/PHD/Matlabtoolboxes/eeglab')
    end
    eeglab
    
    %% get files

    datapath = 'data';
    
    contfn = sprintf('%s/derivatives/eeglab/sub-%02i_task-rsvp_continuous.set',datapath,partid);
    if isfile(contfn)
        fprintf('Using %s\n',contfn)
    	EEG_cont = pop_loadset(contfn);
    else
        % load EEG file
        EEG_raw = pop_loadbv(sprintf('%s/sub-%02i/eeg/',datapath,partid), sprintf('sub-%02i_task-rsvp_eeg.vhdr',partid));
        EEG_raw = eeg_checkset(EEG_raw);
        EEG_raw.setname = partid;
        EEG_raw = eeg_checkset(EEG_raw);

        % re-reference
        EEG_raw = pop_chanedit(EEG_raw, 'append',1,'changefield',{2 'labels' 'Cz'},'setref',{'' 'Cz'});
        EEG_raw = pop_reref(EEG_raw, [],'refloc',struct('labels',{'Cz'},'sph_radius',{[]},'sph_theta',{[]},'sph_phi',{[]},'theta',{[]},'radius',{[]},'X',{[]},'Y',{[]},'Z',{[]},'type',{''},'ref',{''},'urchan',{[]},'datachan',{0}));
        
        % high pass filter
        EEG_raw = pop_eegfiltnew(EEG_raw, 0.1,[]);

        % low pass filter
        EEG_raw = pop_eegfiltnew(EEG_raw, [],100);

        % downsample
        EEG_cont = pop_resample(EEG_raw, 250);
        EEG_cont = eeg_checkset(EEG_cont);
        
        % save
        pop_saveset(EEG_cont,contfn);
    end
    
    %% add eventinfo to events
    eventsfncsv = sprintf('data/sub-%02i/eeg/sub-%02i_task-rsvp_events.csv',partid,partid);
    eventsfntsv = sprintf('data/sub-%02i/eeg/sub-%02i_task-rsvp_events.tsv',partid,partid);
    eventlist = readtable(eventsfncsv);
        
    idx = find(strcmp({EEG_cont.event.type},'E  1'));
    nrem=0;
    if length(idx) < size(eventlist,1)
        idx_seq = strcmp({EEG_cont.event.type},'E  3');
        nseq = sum(idx_seq);
        nseqr = length(unique(eventlist.sequencenumber));
        warning(sprintf('sub-%02i: found only %i sequences instead of %i. Removing first %i sequences',partid,nseq,nseqr,nseqr-nseq))
        eventlist = eventlist(eventlist.sequencenumber>=(nseqr-nseq),:);
        nrem = length(idx)-size(eventlist,1);
        idx = idx((nrem+1):end);
    end
    onset = vertcat(EEG_cont.event(idx).latency);
    duration = ones(size(onset));
    neweventlist = [table(onset,duration,'VariableNames',{'onset','duration'}) eventlist];
    
    writetable(neweventlist,eventsfntsv,'filetype','text','Delimiter','\t')

    %% create epochs
    EEG_epoch = pop_epoch(EEG_cont, {'E  1'}, [-0.100 0.600]);
    EEG_epoch = eeg_checkset(EEG_epoch);
    
    %% convert to cosmo
    ds = cosmo_flatten(permute(EEG_epoch.data(:,:,(nrem+1):end),[3 1 2]),{'chan','time'},{{EEG_epoch.chanlocs.labels},EEG_epoch.times},2);
    ds.a.meeg=struct(); %or cosmo thinks it's not a meeg ds 
    ds.sa = table2struct(eventlist,'ToScalar',true);
    cosmo_check_dataset(ds,'meeg');
    
    %% save epochs
    save(sprintf('%s/derivatives/cosmomvpa/sub-%02i_task-rsvp_cosmomvpa.mat',datapath,partid),'ds','-v7.3')
    
end
