%% Script for Detecting SPW-R events that last for at least 15ms
gaus=
%% Directory
cd 'x:\03. Lab Procedures and Protocols\MATLABToolbox\chronux\spectral_analysis\continuous';

%% LFP Load
% load '';
% 

%% Build the Gaussian Kernel (4-ms s.d - per Karlsson, 2009)
pts_ms=params.Fs/1000; %points/ms in ms scale
pts_env=(1:(ceil(pts_ms*15))); pts_env=ones(1,length(pts_env)); %num data points expected within a 15ms window
sig=pts_ms*4; %4ms worth of data (per Karlsson, 2009 s.d.)
max=ceil(3*sig);
gaus=fspecial('gaussian',2*max+1,sig); %generate 2D gaussian
n=ceil(length(gaus)/2);%round to the half-way row
gaus=gaus(n,:); %make vector equal to middle-row
gaus=gaus./sum(gaus); %correct gaussian kernel to sum to 1
%% Pull out LFP from given channel - NO LIGHT
for k=1:12;
data{1,k}=DelayNoLight{1,k};

%% Bandpass filter

data_filt{1,k}=skaggs_filter_var(data{1,k},140,220,params.Fs); %filters for 150-250Hz

%% Hilbert Transform
data_hil{1,k}=abs(hilbert(data_filt{1,k}));

%% Convolve with Gaussian Kernel - smoothing
data_smooth{1,k}=conv(data_hil{1,k},gaus,'same');

%% Identify + 3stdv from mean abs. amplitude
data_amp_avg{1,k}=mean(data_smooth{1,k}); %calculate the avg amplitude
data_std{1,k}=std(data_smooth{1,k}); %standard dev for absolute amplitudes
data_normax{1,k}=data_amp_avg{1,k}+((data_std{1,k})*3); %maximum expected value in normal distribution defined as 3stds above mean


%% SPW-R detection
data_flag{1,k}=data_smooth{1,k}>data_normax{1,k};% returns boolean for all cells in the absoluted signal that are greater than normal max
data_env{1,k}=strfind(data_flag{1,k},pts_env);% cases in which minimum length of TRUE is found in boolean signal
[data_flagCol{1,k}, data_reps{1,k}, data_ind{1,k}]=RunLength(data_flag{1,k});% Collapsed unbroken sequence of same TRUE/FALSE value; number of items in collapse; cell column value in original signal in which TRUE/FALSE value change
data_multiples{1,k} = find(data_reps{1,k}>length(pts_env));%cases in which the number of collapsed values is greater than minimum
data_SWRchk{1,k}=data_ind{1,k}(data_multiples{1,k}); %the first of cell in which the collapsed sequence started

for q=1:length(data_SWRchk{1,k});
    if data_flag{1,k}(data_SWRchk{1,k}(q))==1; %if the sequence met the duration criteria, is that a sequence of 1's
        data_SWRstart{1,k}(q)=data_SWRchk{1,k}(q); %if so, enter here
    else
        data_SWRstart{1,k}(q)=NaN; %if not enter as NaN
    end
end

data_SWRstart{1,k}(isnan(data_SWRstart{1,k})) = []; %clear NaN's

data_SWRcount_NLt{1,k}=length(data_SWRstart{1,k}); %total number, for each trial, of SPW-R events

%% index events from original signal
for d = 1:data_SWRcount{1,k}; %for the number of detected events
    SWRs{1,k}(d,:)=data{1,k}(1,data_SWRstart{1,k}(1,d)-250:data_SWRstart{1,k}(1,d)+250); %make a matrix with rows representing the event, and columns populated with signal values from start-250 to start+250
end

%filter for ripple band
for d=1:data_SWRcount{1,k};
    SWRs_f{1,k}(d,:)=skaggs_filter_var(SWRs{1,k}(d,1:end),140,220,params.Fs);
end

%filter for theta band
for d=1:data_SWRcount{1,k};
    SWRs_f_the{1,k}(d,:)=skaggs_filter_var(SWRs{1,k}(d,1:end),6,12,params.Fs);
end

end



%% Pull out LFP from given channel - LIGHT
for k=1:12;
data{1,k}=DelayLight{1,k};

%% Bandpass filter

data_filt{1,k}=skaggs_filter_var(data{1,k},150,250,params.Fs); %filters for 150-250Hz

%% Hilbert Transform
data_hil{1,k}=abs(hilbert(data_filt{1,k}));

%% Convolve with Gaussian Kernel - smoothing
data_smooth{1,k}=conv(data_hil{1,k},gaus,'same');

%% Identify + 3stdv from mean abs. amplitude
data_amp_avg{1,k}=mean(data_smooth{1,k}); %calculate the avg amplitude
data_std{1,k}=std(data_smooth{1,k}); %standard dev for absolute amplitudes
data_normax{1,k}=data_amp_avg{1,k}+((data_std{1,k})*3); %maximum expected value in normal distribution defined as 3stds above mean


%% SPW-R detection
data_flag{1,k}=data_smooth{1,k}>data_normax{1,k};% returns boolean for all cells in the absoluted signal that are greater than normal max
data_env{1,k}=strfind(data_flag{1,k},pts_env);% cases in which minimum length of TRUE is found in boolean signal
[data_flagCol{1,k}, data_reps{1,k}, data_ind{1,k}]=RunLength(data_flag{1,k});% Collapsed unbroken sequence of same TRUE/FALSE value; number of items in collapse; cell column value in original signal in which TRUE/FALSE value change
data_multiples{1,k} = find(data_reps{1,k}>length(pts_env));%cases in which the number of collapsed values is greater than minimum
data_SWRchk{1,k}=data_ind{1,k}(data_multiples{1,k}); %the first of cell in which the collapsed sequence started

for q=1:length(data_SWRchk{1,k});
    if data_flag{1,k}(data_SWRchk{1,k}(q))==1; %if the sequence met the duration criteria, is that a sequence of 1's
        data_SWRstart{1,k}(q)=data_SWRchk{1,k}(q); %if so, enter here
    else
        data_SWRstart{1,k}(q)=NaN; %if not enter as NaN
    end
end

data_SWRstart{1,k}(isnan(data_SWRstart{1,k})) = []; %clear NaN's

data_SWRcount_Lt{1,k}=length(data_SWRstart{1,k}); %total number, for each trial, of SPW-R events
end
%% pull time stamps of those event starts from the searched boolean
% format long g;
% HC_NLt_SWRevent{1,k} = []; %open empty var
% 
% for h = 1:length(HC_NLt_SWRstart{1,k});
%     HC_NLt_SWRevent{1,k}(h) = DelayTimesNoLight{1,k}(HC_NLt_SWRstart{1,k}(h)); %retrieve timestamps from cells matching those identified as having SWR amps
% end

%% Save output
save ('X:\08. Lab personnel\Current\David\Projects\Re Suppression - HC Modulation\2. Output\Ephys\Data Analysis\Santiago\DelayPhaseLtNLt_SWRs.mat','-v7.3');
