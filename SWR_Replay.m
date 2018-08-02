%% load txt files
cd 'H:\Classwork Files\2017\Spring\Spatial Cognition\Project\1101-16'
% cd 'Z:\01.Experiments\Completed Studies\DualTask_CDAlternation_HippocampusRecording\1101\1101-16' %%amy's comp.
load 'TT1_SS_01.txt';
load 'TT1_SS_02.txt';
load 'TT6hh_SS_01.txt';
load 'TT6hh_SS_02.txt';
load 'TT6hh_SS_03.txt';
load 'TT6hh_SS_04.txt';
load 'TT13hh_SS_01.txt';
load 'TT15hh_SS_01.txt';
load 'TT15hh_SS_02.txt';
load 'TT15hh_SS_03.txt';
load 'TT15hh_SS_04.txt';
load('Intervals.mat')

%% organize spikes into Cell Array
Trial10Spikes{1,1}=(TT1_SS_01>Int1(10,1) & TT1_SS_01<Int1(10,8));
totalspikes{1,1}=TT1_SS_01;
totalspikes{1,2}=TT1_SS_02;
totalspikes{1,3}=TT6hh_SS_01;
totalspikes{1,4}=TT6hh_SS_02;
totalspikes{1,5}=TT6hh_SS_03;
totalspikes{1,6}=TT6hh_SS_04;
totalspikes{1,7}=TT13hh_SS_01;
totalspikes{1,8}=TT15hh_SS_01;
totalspikes{1,9}=TT15hh_SS_02;
totalspikes{1,10}=TT15hh_SS_03;
totalspikes{1,11}=TT15hh_SS_04;
spiketimes_11cells=totalspikes;
clear totalspikes

%% identify spikes during trial 10
for i=1:size(spiketimes_11cells,2);
    Trial10Spikes{1,i}=find(spiketimes_11cells{1,i}>Int1(10,1) & spiketimes_11cells{1,i}<Int1(10,8));
end;
for i=1:size(spiketimes_11cells,2);
    Trial10Spikes_log{1,i}=log10(find(spiketimes_11cells{1,i}>Int1(10,1) & spiketimes_11cells{1,i}<Int1(10,8)));
end;

%% plot spike sequence
figure; subplot 212; hold on;
plot(Trial10Spikes_log{1,1},1,'r.');
plot(Trial10Spikes_log{1,3},2,'g.');
plot(Trial10Spikes_log{1,2},3,'c.');
plot(Trial10Spikes_log{1,8},4,'k.');
plot(Trial10Spikes_log{1,6},5,'.');
plot(Trial10Spikes_log{1,7},6,'m.');

%% collect sequenced spikes
seq_cells{1,1}=spiketimes_11cells{1,1}(:)';
seq_cells{1,2}=spiketimes_11cells{1,3}(:)';
seq_cells{1,3}=spiketimes_11cells{1,2}(:)';
seq_cells{1,4}=spiketimes_11cells{1,8}(:)';
seq_cells{1,5}=spiketimes_11cells{1,6}(:)';
seq_cells{1,6}=spiketimes_11cells{1,7}(:)';
%% LFP
load 'H:\Classwork Files\2017\Spring\Spatial Cognition\Project\1101-16\CSC6.mat'
CSC6=Samples(:)';
CSC_Timestamps=linspace(Timestamps(1,1),Timestamps(1,end),length(CSC6));

%% break out trial 10 LFP
CSC6_T10=CSC6(find(CSC_Timestamps>Int1(10,1) & CSC_Timestamps<Int1(10,8)));
CSC_Timestamps_10=find(CSC_Timestamps>Int1(10,1) & CSC_Timestamps<Int1(10,8));
%% adjust time stamps to second scale
CSC_Timestamps_10_real=CSC_Timestamps_10 ./ 32544;
Time_T10=CSC_Timestamps_10_real - 359.5205;
subplot 211; plot(Time_T10,CSC6_T10);
%% end of traversal
CSC6_T10_Timestamps_end=find(Time_T10>10 & Time_T10<20);
CSC6_T10_end=CSC6_T10(find(Time_T10>10 & Time_T10<20));
Time_T10_end=CSC6_T10_Timestamps_end ./ 32544;
figure; subplot 211; plot(Time_T10_end, CSC6_T10_end);
%% BPass Filter 150-250Hz
% cd 'x:\03. Lab Procedures and Protocols\MATLABToolbox\chronux\spectral_analysis\continuous';
% CSC6_filtered=skaggs_filter_var(CSC6,150,250,32544);
% 
% %% Identify +/- 3stdv from mean abs. amplitude
CSC6_T10_end_amp=abs(CSC6_T10_end);
CSC6_T10_end_avg=mean(CSC6_T10_end_amp);
CSC6_T10_end_max=CSC6_T10_end_avg + 4*std(CSC6_T10_end_amp);

%% SWR detection
SWR_events=CSC6_T10_end_amp>CSC6_T10_end_max;
SWR_events_key=find(SWR_events(1,:)>0);

%% isolate segments of trial 10 signal that have amplitudes in range