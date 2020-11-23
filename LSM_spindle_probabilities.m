% Compute the probability of a spindle in 0.5 s intervals.
%
% INPUTS:
%  data = [channels, time]
%  hdr  = structure that MUST HAVE ...
%       hdr.info.sfreq      = sampling frequency [Hz]
%       hdr.info.ch_names   = cell of channel names.
%  options
%       options.MinPeakProminence  = sets MinPeakPromience.
%       options.StartFrequency     = sets spindle analysis to frequencies [X,Y].
%              .StopFrequency        NOTE: must set both Start and Stop frequencies.
%
% OUTPUT:
%  spindle_probabilities = [channels] structure with ...
%                       .label       = channel name
%                       .prob        = probabiltiy of a spindle at time t.
%                       .t           = time t [s]
%                       .Fs          = sampling frequency [Hz]
%
% The output is typically an input to function LSM_spindle_detections.m

function spindle_probabilities = LSM_spindle_probabilities(data, hdr, options)

  MinPeakProminence = 2e-6;                                         % Default value for HD scalp EEG.
  start_frequency = [];                                             % 9-15 Hz analysis
  stop_frequency  = [];
  feature = 'broadband';                                            % ... is broadband.
  
  if nargin>2                                                       % ---- Adjust default settings. ----
      if isfield(options,'MinPeakProminence')                       % Set MinPeakProminence for Fano step.
          MinPeakProminence = options.MinPeakProminence;
      else
          MinPeakProminence = 2e-6;
      end
          
      if isfield(options,'StartFrequency')                          % Set frequencies for narrowband analsysis.
          start_frequency = options.StartFrequency;
          stop_frequency  = options.StopFrequency;
          feature = 'narrowband';
      end
  end

  fprintf(['Detecting spindles with spectral features: ' feature ' ' num2str(start_frequency) ' ' num2str(stop_frequency) '\n'])                                 
                                             
  % Build filter 3-25 Hz.
  Fs = hdr.info.sfreq;
  bpFilt = designfilt('bandpassfir', ...
      'StopbandFrequency1', 1, ...
      'PassbandFrequency1', 3, ...
      'PassbandFrequency2', 25, ...
      'StopbandFrequency2', 30, ...
      'PassbandRipple', 0.1, ...
      'StopbandAttenuation1', 40, ...
      'StopbandAttenuation2', 20, ...
      'SampleRate', Fs);
  
  % Make sure to use the MATLAB version of findpeaks.
  current_dir = pwd;                                                % Get current directory.
  s = regexp(path, pathsep, 'split');                               % Create list of all paths.
  i0 = find(endsWith(s, 'toolbox/signal/signal'));                  % Find path location of 'findpeaks'.
  cd(s{i0})                                                         % Go there,
  findpeaks_vMAT = str2func('findpeaks');                           % ... and point to 'findpeaks' function.
  cd(current_dir);                                                  % Return to current directory.

  electrodes_to_analyze = hdr.info.ch_names;
  K = length(electrodes_to_analyze);  
  spindle_probabilities = struct('label',cell(K,1),'prob',cell(K,1), 't',cell(K,1), 'Fs',cell(K,1));
  
  for i_channel = 1:K                                               % NOTE: This can be parfor
    
      [likelihood, mu, params, sigma, transition_matrix] = load_inputs();  
      channel = electrodes_to_analyze{i_channel};
      fprintf(['... ' num2str(channel) '(' num2str(i_channel) ' of ' num2str(length(electrodes_to_analyze)) ') \n'])
      
      i0 = find(strcmp(hdr.info.ch_names, channel));
        
      if ~isempty(i0)
          
          d = data(i0,:);                                           % Get channel to analyze.  

          if any(isnan(d))
              fprintf(['... detected NaNs in data, replacing with 0s for filter only. \n'])
              d0 = d; d0(isnan(d0))=0;
              dfilt = filtfilt(bpFilt, d0);
          else
              dfilt  = filtfilt(bpFilt, d);                         % Filter the data.
          end
          
          extent = max(d) - min(d);
          if (MinPeakProminence == 2e-6) && (extent > 300e-6 || extent < 10e-6)
              fprintf(['Are your data in microvolts? If not, set options.MinPeakProminence \n'])
              break
          end

          % Initialize
          p = [0.5; 0.5];
          [prob, t] = deal([]);
          window_size = round(params.window_duration*Fs);
          step_size   = round(params.step_duration*Fs);
          time_to_analyze = length(d)-window_size;

          i=1;
          while i < time_to_analyze                                  % For each time interval, ...
              
              % Get the power features in this interval.
              interval  = i:i+window_size-1;
              d0        = d(interval);
              d0        = detrend(d0);
              pow       = abs(fft(hann(length(d0)).*d0', round(Fs)));
              pow       = pow / sum(pow);
              
              pow_theta = pow(params.theta_index);
              pow_9_15  = mean(pow(params.nine_15_index));
              
              % Analyze fano factor in this interval.
              d0 = dfilt(interval);
              d0 = detrend(d0);
              fastest_spindle_period  = 1/35;
              MinPeakDistance         = ceil(fastest_spindle_period * Fs);
              [~, pos_locs] = findpeaks_vMAT( d0, 'MinPeakDistance', MinPeakDistance, 'MinPeakProminence',MinPeakProminence);
              [~, neg_locs] = findpeaks_vMAT(-d0, 'MinPeakDistance', MinPeakDistance, 'MinPeakProminence',MinPeakProminence);
              
              if length(pos_locs)>1 && length(neg_locs)>1           % if you have more than 1 peak, and more than 1 trough,
                  ISI  = [diff(pos_locs) diff(neg_locs)];           % ... then compute ISI for each, and average.
              else
                  ISI  = nan;                                       % otherwise, not enough points, so ISI = nan.
              end
              
              if length(ISI)>1                                      % if you have more than 1 ISI,
                  fano = var(ISI)/mean(ISI);                        % ... then compute the fano factor.
              else
                  fano=nan;                                         % otherwise, not enough ISI, so fano=nan.
              end
              instant_freq = 1/( mean(ISI)/Fs );                    % Compute instantaneous freq.
              
              % Get the 1 step prediction.
              p = transition_matrix * p;
              
              % Normalize 1 step prediction.
              p = p / sum(p);
              
              % Select feature and compute posterior.
              switch feature
                  
                  case 'broadband'
                      if isnan(log(fano)) || ~isfinite(log(fano))
                          posterior = [likelihood(log(pow_9_15),  mu.log_P_9_15_1, sigma.log_P_9_15_1, 1)...
                                      *likelihood(log(pow_theta), mu.log_P_theta1, sigma.log_P_theta1, 1)...
                                  ... *likelihood(fano,      mu.F1,       sigma.F1,       scale.F1)...
                                      *p(1),...
                                       likelihood(log(pow_9_15),  mu.log_P_9_15_0, sigma.log_P_9_15_0, 1)...
                                      *likelihood(log(pow_theta), mu.log_P_theta0, sigma.log_P_theta0, 1)...
                                   ...*likelihood(fano,      mu.F0,       sigma.F0,       scale.F0)...
                                      *p(2)];
                      else
                          posterior = [likelihood(log(pow_9_15),  mu.log_P_9_15_1, sigma.log_P_9_15_1, 1)...
                                      *likelihood(log(pow_theta), mu.log_P_theta1, sigma.log_P_theta1, 1)...
                                      *likelihood(log(fano),      mu.F1,           sigma.F1,           1)...
                                      *p(1),...
                                       likelihood(log(pow_9_15),  mu.log_P_9_15_0, sigma.log_P_9_15_0, 1)...
                                      *likelihood(log(pow_theta), mu.log_P_theta0, sigma.log_P_theta0, 1)...
                                      *likelihood(log(fano),      mu.F0,           sigma.F0,           1)...
                                      *p(2)];
                      end
                      
                  %%%% Narrowband analysis ---------------------------------------------------------------------------------------------------------------
                      
                  case 'narrowband'
                      if isnan(log(fano)) || ~isfinite(log(fano))
                          posterior = [likelihood(log(pow_9_15),  mu.log_P_9_15_1, sigma.log_P_9_15_1, 1)...
                                      *likelihood(log(pow_theta), mu.log_P_theta1, sigma.log_P_theta1, 1)...
                                      *narrowband_likelihood(instant_freq, [start_frequency,stop_frequency]) ...
                                  ... *likelihood(fano,      mu.F1,       sigma.F1,       scale.F1)...
                                      *p(1),...
                                       likelihood(log(pow_9_15),  mu.log_P_9_15_0, sigma.log_P_9_15_0, 1)...
                                      *likelihood(log(pow_theta), mu.log_P_theta0, sigma.log_P_theta0, 1)...
                                   ...*likelihood(fano,      mu.F0,       sigma.F0,       scale.F0)...
                                      *p(2)];
                      else
                          posterior = [likelihood(log(pow_9_15),  mu.log_P_9_15_1, sigma.log_P_9_15_1, 1)...
                                      *likelihood(log(pow_theta), mu.log_P_theta1, sigma.log_P_theta1, 1)...
                                      *likelihood(log(fano),      mu.F1,           sigma.F1,           1)...
                                      *narrowband_likelihood(instant_freq, [start_frequency,stop_frequency]) ...
                                      *p(1),...
                                       likelihood(log(pow_9_15),  mu.log_P_9_15_0, sigma.log_P_9_15_0, 1)...
                                      *likelihood(log(pow_theta), mu.log_P_theta0, sigma.log_P_theta0, 1)...
                                      *likelihood(log(fano),      mu.F0,           sigma.F0,           1)...
                                      *p(2)];
                      end
                      
              end
              
              % Normalize posterior;
              p = posterior / sum(posterior);
              p = transpose(p);
              
              if all(isnan(p))                                      % if the data are omitted,
                  p = [0.5; 0.5];                                   % ... then set p to unknown state.
              end
              
              prob = [prob, p];                                     % save the probabilities,
              t    = [t, i/Fs];                                     % ... and the time
              i = i + step_size;
          end
          
          spindle_probabilities(i_channel).label       = channel;   % save results for this channel
          spindle_probabilities(i_channel).prob        = prob(1,:);
          spindle_probabilities(i_channel).t           = t;
          spindle_probabilities(i_channel).Fs          = Fs;
      end
            
  end
  
end 

function [likelihood, mu, params, sigma, transition_matrix] = load_inputs()
  % Use default likelihood file (trained on Chu-lab CECTS data). Don't alter this unless you know what you're doing.
  load('likelihood_and_transition_matrix_LSM_LOO_min_manual_duration_0_0.5_0.1_pt999_PDF_use_all_data_compute_all_features_chu.mat')
  fprintf(['... using default likelihood file. \n'])
  
  % Use alternative likelihood file trained on Manoach-lab data.
  %load('likelihood_and_transition_matrix_LSM_LOO_min_manual_duration_0_0.5_0.1_pt999_PDF_use_all_data_compute_all_features_manoach.mat')
  %fprintf(['... using alternative likelihood file. \n'])
end

function [L0] = narrowband_likelihood(f, range)
  faxis = (0:0.1:100);
  L     = zeros(size(faxis));
  i0    = find(faxis >= range(1) & faxis <= range(2));
  L(i0) = 1;
  [~, imn] = min(abs(faxis - f));
  L0    = L(imn);
end 