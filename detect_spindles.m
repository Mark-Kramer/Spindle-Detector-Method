function [spindle_det, keep_detection] = detect_spindles(current_data, Data, data_filt, threshold, Fs, tthresh, segmDurBef, segmDurAft, params)

  current_data(isnan(current_data)) = 0;

  spindle_det = [];

  if any(isnan(current_data))
      spindle_det(1).bads = 1;
  else
      spindle_det(1).bads = 0;
  end

  over = current_data>threshold; % Mark all points over threshold as '1'
  locs = (zeros(1,length(current_data)))'; % Create a vector of zeros the length of the MS signal

  min_duration = round(Fs*tthresh);
  for i=1:(length(current_data)-min_duration) % for the length of the signal, if the sum of 40 concurrent points = Fs*0.4, mark a spindle
      i_start = i;
      i_stop  = i_start + min_duration - 1;
      if sum(over(i_start:i_stop)) == min_duration
          locs(i,1)=1;
      end
  end

  spin=zeros((length(locs)),1);  % only mark a spindle in vector 'spin' at the end of a 400ms duration peak
  for i=1:length(locs)
      if locs(i,1)==1 && locs(i+1,1)==0
          spin(i,1)=1;
      end
  end

  % added 9/17/2012: for every spindle marked in 'spin', delete the spindle if there is also a spindle within the second preceeding it.
  for i=(Fs+1):length(spin)
      if spin(i,1)==1 && sum(spin((i-Fs):(i-1)))>0
          spin(i,1)=0;
      end
  end

  %spindle_det(1).label = ch_names;
  spindle_det(1).sample = find(spin);
  spindle_det(1).spindle_count = sum(spin);
  %spindle_det(1).backgr_mean = signalmean(1);

  %%%% Discard spindles that are closer than 2 sec to the boundaries
  bound_idx = find(spindle_det(1).sample + Fs*segmDurAft > length(current_data) | spindle_det(1).sample - Fs*segmDurAft < 1);
  if ~isempty(bound_idx)
      spindle_det(1).sample(bound_idx) = [];
      spindle_det(1).spindle_count = length(spindle_det(1).sample);
  end

  %%%% Calculate Spindle Parameters
  numSpindles = spindle_det(1).spindle_count;
  numSpindles(isnan(numSpindles)) = 0;
  
  if numSpindles == 0 % no spindles on this channel
      keep_detection = NaN;
      return;
  else
    
    segms = cell(numSpindles,1);
    
    % Make 4 sec segments around each spindle event
    for jj = 1:numSpindles
        startSmp = max(spindle_det(1).sample(jj)-Fs*segmDurBef+1,1);
        endSmp = min(spindle_det(1).sample(jj)+Fs*segmDurAft,length(current_data));
        segms{jj} = Data(startSmp:endSmp)';
    end
    
    % get rid of truncated segments
    tmp = cellfun(@length, segms); % Find shorter segments
    segms(tmp ~= Fs*(segmDurBef+segmDurAft)) = {NaN*ones(1,Fs*(segmDurBef+segmDurAft))}; %Turn their data to nans
    segms = cell2mat(segms);
    
    [spindle_det(1).duration, spindle_det(1).energy, spindle_det(1).peakAmp, spindle_det(1).peakLoc, spindle_det(1).startSample, spindle_det(1).endSample, spindle_det(1).peakFreq, spindle_det(1).SigmaPower] = fun_Calculate_Spindle_Features(segms, Fs); % ******
    
    % Convert to whole-signal timeline (samples given relative to 1 sec before detection point)
    spindle_det(1).startSample  = spindle_det(1).startSample + spindle_det(1).sample' - Fs*(segmDurBef-1)*ones(size([spindle_det(1).startSample]));
    spindle_det(1).endSample    = spindle_det(1).endSample   + spindle_det(1).sample' - Fs*(segmDurBef-1)*ones(size([spindle_det(1).endSample]));
    spindle_det(1).peakLoc      = spindle_det(1).peakLoc     + spindle_det(1).sample' - Fs*(segmDurBef-1)*ones(size([spindle_det(1).peakLoc]));
    
    %% Discard spindles that last less than x samples
    spindle_det(1).sample([spindle_det(1).duration]<tthresh)      = [];
    spindle_det(1).energy([spindle_det(1).duration]<tthresh)      = [];
    spindle_det(1).peakAmp([spindle_det(1).duration]<tthresh)     = [];
    spindle_det(1).peakLoc([spindle_det(1).duration]<tthresh)     = [];
    spindle_det(1).startSample([spindle_det(1).duration]<tthresh) = [];
    spindle_det(1).endSample([spindle_det(1).duration]<tthresh)   = [];
    spindle_det(1).peakFreq([spindle_det(1).duration]<tthresh)    = [];
    spindle_det(1).SigmaPower([spindle_det(1).duration]<tthresh)  = [];
    spindle_det(1).duration([spindle_det(1).duration]<tthresh)    = [];
    spindle_det(1).spindle_count = length(spindle_det(1).duration);
    
    %% Discard spindles that overlap
    if length(spindle_det(1).startSample)>=2
        overl_idx = find(spindle_det(1).startSample(2:end)-spindle_det(1).endSample(1:end-1)<0);
        
        spindle_det(1).sample(overl_idx)      = [];
        spindle_det(1).energy(overl_idx)      = [];
        spindle_det(1).peakAmp(overl_idx)     = [];
        spindle_det(1).peakLoc(overl_idx)     = [];
        spindle_det(1).startSample(overl_idx) = [];
        spindle_det(1).endSample(overl_idx)   = [];
        spindle_det(1).duration(overl_idx)    = [];
        spindle_det(1).peakFreq(overl_idx)    = [];
        spindle_det(1).SigmaPower(overl_idx)    = [];
        spindle_det(1).spindle_count = length(spindle_det(1).duration);
        
    end
    
    fprintf(['Spindle count ' num2str(spindle_det(1).spindle_count) '\n'])
    
    %% Discard spindles that co-occur with a spike.
    %spike_threshold = 80*1e-6;
    %slowest_spindle_period  = 1/9.5;
    spike_threshold        = params.spike_threhsold;
    slowest_spindle_period = params.slowest_spindle_period;
    
    fastest_spindle_period  = 1/35;
    MinPeakDistance         = ceil(fastest_spindle_period * Fs);
    
    keep_detection = zeros(spindle_det(1).spindle_count,1);
    for k=1:spindle_det(1).spindle_count
        d0 = data_filt(spindle_det(1).startSample(k):spindle_det(1).endSample(k));
        [pos_pks, pos_locs] = findpeaks( d0, 'MinPeakDistance', MinPeakDistance, 'MinPeakProminence',2e-6);
        [neg_pks, neg_locs] = findpeaks(-d0, 'MinPeakDistance', MinPeakDistance, 'MinPeakProminence',2e-6);
        neg_pks = -neg_pks;
        % Compare each peak to subsequent trough.
        ht_forward = NaN(size(pos_locs));
        for n=1:length(pos_locs)
            i_pos = pos_locs(n);
            i_neg = neg_locs( neg_locs > i_pos );
            p_neg = neg_pks ( neg_locs > i_pos );
            [~, imn] = min(abs(i_neg - i_pos));
            if ~isempty(imn)
                ht_forward(n) = pos_pks(n) - p_neg(imn);
            end
        end
        % Compare each peak to previous trough.
        ht_backward = NaN(size(pos_locs));
        for n=1:length(pos_locs)
            i_pos = pos_locs(n);
            i_neg = neg_locs( neg_locs < i_pos );
            p_neg = neg_pks ( neg_locs < i_pos );
            [~, imn] = min(abs(i_neg - i_pos));
            if ~isempty(imn)
                ht_backward(n) = pos_pks(n) - p_neg(imn);
            end
        end
        
        ISI  = diff(pos_locs);
        fano = var(ISI)/mean(ISI);
        
        min_num_cycles = floor(spindle_det(1).duration(k) / slowest_spindle_period);
        max_num_cycles = floor(spindle_det(1).duration(k) / fastest_spindle_period);
        %if sum([ht_forward, ht_backward] > spike_threshold) > 2 || sum(isfinite(ht_forward)) < min_num_cycles || sum(isfinite(ht_forward)) > max_num_cycles
        if length(pos_pks) < min_num_cycles || length(pos_pks) > max_num_cycles || fano > 5 || length(pos_pks) <= 2 || any([ht_forward, ht_backward] > spike_threshold)
            keep_detection(k) = 0;
        else
            keep_detection(k) = 1;
        end
        
%                 t0 = (1:length(d0))/Fs;
%                 plot(t0,d0)
%                 hold on
%                 plot(pos_locs/Fs, pos_pks, 'o')
%                 plot(neg_locs/Fs, neg_pks, 'o')
%                 hold off
%                 axis tight
%                 xlabel(['Time [s]'])
%                 title([  'Too few cycles? ' num2str(length(pos_pks) < min_num_cycles) ...
%                     ', Too big? ' num2str(any([ht_forward, ht_backward] > spike_threshold)) ...
%                     ', Fano? ' num2str(fano) ...
%                     ', Keep it? ' num2str(keep_detection(k))])
    end
  end
end