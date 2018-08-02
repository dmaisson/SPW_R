%% bring in data
cd('X:\08. Lab personnel\Current\David\Projects\Ephys\Re Suppression - HC Modulation\2. Output\Ephys\Data Analysis\Hermann\Sample Day');
load 'Sample Day 1 - 12-11-2017_Processed.mat';
for i = 1:16;
    for j = 1:24;
    lfp{i,j} = DDelay{i,j};
    end
end
cfg = [];

%% filter in SWR band
cd 'X:\03. Lab Procedures and Protocols\MATLABToolbox\chronux\spectral_analysis\continuous';
cfg = [];
cfg.f = [140 220];
cfg.display_filter = 0;
for i = 1:16;
    for j = 1:24;
SWRf{i,j} = skaggs_filter_var(lfp{i,j},140,220,params.Fs);
    end
end
 
%% obtain power and z-score it
SWRp = LFPpower([],SWRf);
SWRp_z = zscore_tsd(SWRp);
 
%% detect events
cfg = [];
cfg.method = 'raw';
cfg.threshold = 3;
cfg.operation =  '>'; % return intervals where threshold is exceeded
cfg.merge_thr = 0.05; % merge events closer than this
cfg.minlen = 0.05; % minimum interval length
 
SWR_evt = TSDtoIV(cfg,SWRp_z);
 
%% to each event, add a field with the max z-scored power (for later selection)
cfg = [];
cfg.method = 'max'; % 'min', 'mean'
cfg.label = 'maxSWRp'; % what to call this in iv, i.e. usr.label
 
SWR_evt = AddTSDtoIV(cfg,SWR_evt,SWRp_z);
 
%% select only those events of >5 z-scored power
cfg = [];
cfg.operation = '>';
cfg.threshold = 5;
 
SWR_evt = SelectIV(cfg,SWR_evt,'maxSWRp');
 
%% plot events in highlighted on top of full lfp
PlotTSDfromIV([],SWR_evt,lfp);
 
%% ..or the events alone (fixed 200ms window centered at event time)
close all;
 
cfg = [];
cfg.display = 'iv';
cfg.mode = 'center';
cfg.fgcol = 'k';
 
PlotTSDfromIV(cfg,SWR_evt,lfp);
%% ..hold on (highlight edges of event on top of previous plot)
cfg = [];
cfg.display = 'iv';
cfg.fgcol = 'r';
 
PlotTSDfromIV(cfg,SWR_evt,lfp);