function [spindle_det] = fun_spindle_detection_double_run(Data, ch_names, Fs, tthresh, ampl_factor,sq_opt,params)

 %% Input:
 % Data       : double (channels x timpepoints)
 % ch_names   : cell  {1xn} with the channel labels
 % Fs         : double,  sampling frequency
 %-------------- NEW 3/9/19 time is [s], not [ms] -------------------------
 % tthresh    : double,  time (s) that the signal has to esceed a threshold (default := 0.400)
 %-------------------------------------------------------------------------
 % ampl_facot : double amplification factor for defining baseline (default := 9)
 % sq_opt     : string, 'Square'|'NoSquare' Option to double square the signal (default := 'NoSquare')
 %-------------- NEW 3/9/19 -----------------------------------------------
 % params     : structure to specify fc, fb, center_freq_desired [Hz]
 %-------------------------------------------------------------------------
 
% Mother wavelet center frequency
fc = params.fc;
% Mother wavelet bandwidth
fb = params.fb;
%-------------- NEW 3/9/19 fb and fc were switched ------------------------
% Define the wavelet
wname = ['cmor' num2str(fb) '-' num2str(fc)];
%--------------------------------------------------------------------------

% Seconds before and after detection to determine spindle parameters
segmDurBef = 2; % in secs
segmDurAft = 2;

% Choose the scale based on the desired center frequency.
center_freq_desired = params.center_freq_desired;
scale = scal2frq(center_freq_desired, wname, 1/Fs);
f     = scal2frq(scale,               wname, 1/Fs);
fprintf(['Scale is ' num2str(scale) ', Center freq is ' num2str(f) '\n'])

%-------------- NEW 3/9/19 visualize the wavelet --------------------------
% Visualize spindle detector
% figure(10)
% [psi,xval] = wavefun(wname, 20);
% plot(xval, real(psi));
% title(['fc ' num2str(fc) ', fb ' num2str(fb) ', Center freq ' num2str(f)])
%--------------------------------------------------------------------------

% Create bandpass filter for spike detection.
bpFilt = designfilt('bandpassfir', ... %'FilterOrder',200, ...
    'StopbandFrequency1', 1, ...
    'PassbandFrequency1', 3, ...
    'PassbandFrequency2', 25, ...
    'StopbandFrequency2', 30, ...
    'PassbandRipple', 0.1, ...
    'StopbandAttenuation1', 40, ...
    'StopbandAttenuation2', 20, ...
    'SampleRate', Fs);

nchannel=size(Data,1); % define number of channels

%% Calculate the continuous wavelet transform for each channel
for x=1:nchannel
    current_data  = Data(x,:);
    EEGWave       = cwt(current_data,scale,wname);  % conducts the wavelet tranformation on workspace variable "EEG"
    EEGChannel    = real(EEGWave.^2);               % Takes only the real component of the coefficients
    EEGChannel    = EEGChannel';
    Data_new(:,x) = EEGChannel;
    DataFilt(x,:) = filtfilt(bpFilt, Data(x,:));
end

%% Set the Sampling Frequency and take Moving Average
switch sq_opt
    case 'Square'
        Data_new = Data_new.^2;
    case 'NoSquare'
        Data_new =  abs(Data_new);
end

window   = ones(round(Fs/10),1)/(Fs/10);    % create 100ms window to convolve with
Data2    = filter(window,1,Data_new);       % take the moving average using the above window

%% Compute Threshhold (Separately for Each Channel)

spindle_det = [];

for ch=1:nchannel
%ch=1;

    msg = ['Working on Channel ',ch_names{ch} '\n'];
    fprintf(msg);
    
    %% Get data from a channel.
    current_data = Data2(:,ch);
    data_filt    = DataFilt(ch,:);
    keep_detection = 0;
    
    while any(keep_detection == 0)
        
        %% Compute the threshold.
        switch params.threshold_type
            case 'mean'
                signalmean = nanmean(current_data);     % compute mean amplitude of rectified signal
                
            case 'median'
                signalmean = nanmedian(current_data);   % compute mean amplitude of rectified signal
        end
        threshold  = signalmean.*ampl_factor;
    
        %% Detect spindles
        [spindle_det0, keep_detection] = detect_spindles(current_data, data_filt, threshold, Fs, tthresh, segmDurBef, segmDurAft, params.slowest_spindle_period);
        
        fprintf(['Num of bad detections is ' num2str(sum(keep_detection==0)) ', threshold = ' num2str(threshold, 4) '\n'])
            
        for k=1:length(keep_detection)                  % if there are detections to ignore, then set data at those times to NaN, and repeat.
            if keep_detection(k) == 0 
                current_data(spindle_det0.startSample(k):spindle_det0.endSample(k)) = NaN;
            end
        end
        
    end
    
    if spindle_det0.spindle_count > 0
    
        spindle_det(ch).label       = ch_names{ch};
        spindle_det(ch).sample      = spindle_det0.sample;
        spindle_det(ch).energy      = spindle_det0.energy;
        spindle_det(ch).peakAmp     = spindle_det0.peakAmp;
        spindle_det(ch).peakLoc     = spindle_det0.peakLoc;
        spindle_det(ch).startSample = spindle_det0.startSample;
        spindle_det(ch).endSample   = spindle_det0.endSample;
        spindle_det(ch).duration    = spindle_det0.duration;
        spindle_det(ch).peakFreq    = spindle_det0.peakFreq;
        spindle_det(ch).spindle_count = length(spindle_det0.duration);
        spindle_det(ch).threshold   = threshold;
        
    else
        spindle_det(ch).label       = ch_names{ch};
        spindle_det(ch).spindle_count = 0;
        spindle_det(ch).threshold   = threshold;
        
    end
        
    
end

%fprintf(repmat('\b',1,nmsg)); % Delete previous msg
disp('Finished detecting spindles.');
disp('**************');

