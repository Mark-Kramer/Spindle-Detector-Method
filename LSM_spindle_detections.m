% Give the spindle probabilities, compute the spindle detections (start and end times of each spindle).
%
% INPUT:
%  spindle_probabilities = [channels] structure output by LSM_spindle_probabilities.m
%
% OUTPUT:
%  spindle_det = [channels] structure with ...
%             .startSample = vector of spindle start times [s].
%             .endSample   = vector of spindle end times [s].
%     	      .label       = channel label.
%             .Fs          = sampling frequency [Hz].

function [spindle_det] = LSM_spindle_detections(spindle_probabilities, varargin)

  % Don't alter this unless you know what you're doing. -------------------
  prob_threshold = 0.95;
  
  if nargin > 1
      
      alter_threshold = find(strcmp(varargin, 'prob_threshold'));
      if ~isempty(alter_threshold); prob_threshold = varargin{alter_threshold+1}; end
      
  end

  min_spindle_duration         = 0.5;           % [s]
  spindle_separation_threshold = 1;             % [s]
  % -----------------------------------------------------------------------
  
  spindle_det = [];
  for k=1:length(spindle_probabilities)                                 % For each channel,
  
      %%%% 1) Get initial set of detections. ------------------------------
      prob = spindle_probabilities(k).prob;                             % Get probability for this channel.
      t    = spindle_probabilities(k).t;                                % Get time for this channel.
      good = find(prob >= prob_threshold);                              % Where probability > threshold,
      detections = zeros(size(prob));   
      detections(good) = 1;                                             % ... set detections = 1.
      
      if all(detections)
          startTimes = t(1);                                            % Catch case of all data = detection
          endTimes   = t(end);
      else
          detections(1) = 0;                                            % Confirm no detections at first index,
          detections(end)=0;                                            % ... or last index.
          startTimes = t(find(diff(detections)==1)+1);                  % Get startTimes,
          endTimes   = t(find(diff(detections)==-1)+1);                 % ... and endTimes from detections.
      end
      
      %%%% 2) Keep detections that are long enough. -----------------------
      detection_duration = endTimes - startTimes;
      long_enough        = find(detection_duration >= min_spindle_duration);
      startTimes         = startTimes(long_enough);
      endTimes           = endTimes(long_enough);
      
      %%%% 3) Merge detections too close. ----------------------------- % InterSpindle Interval
      ISI       = [nan, startTimes(2:end) - endTimes(1:end-1)];         % ... Start time of next spindle - End time of this spindle.
      too_close = find(ISI < spindle_separation_threshold);             % Find spindles too close
      
      for j=length(too_close):-1:1                                      % ... and merge them.
          endTimes(too_close(j)-1) = endTimes(too_close(j));            % End time of previous spindle = end time of next spindle.
          startTimes(too_close(j)) = nan;                               % ... and remove the next spindle,
          endTimes(too_close(j))   = nan;                               % ... since it's now overlapped by previous spindle.
      end
      startTimes = startTimes(isfinite(startTimes));                    %  Keep only spindles not too close.
      endTimes   = endTimes(isfinite(endTimes));
              
      %%%% 4) Return detections for this elec. ----------------------------
      Fs = spindle_probabilities(k).Fs;
      spindle_det(k).startSample = round(startTimes * Fs);
      spindle_det(k).endSample   = round(endTimes * Fs);
      spindle_det(k).label       = spindle_probabilities(k).label;
      spindle_det(k).Fs          = Fs;
      
  end
  
end