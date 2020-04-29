% Visualize spindle detections at one channel.
%
% INPUTS:
%  data        = [channels, time]
%  hdr         = structure that MUST HAVE ...
%                 hdr.info.sfreq      = sampling frequency [Hz]
%                 hdr.info.ch_names   = cell of channel names.
%  spindle_det = output structure from LSM_spindle_detections.m
%  channel     = string with channel name to analyze.
%
% OUTPUT:
%  None
%
% The functions plots the data at "channel" and indicates spindles with red
% boxes.

function LSM_spindle_visualizer(data, hdr, spindle_det, channel)

  i0 = find(strcmp(hdr.info.ch_names, channel));        % Channel to visualize.
  
  Fs = hdr.info.sfreq;
  t  = (1:size(data,2))/Fs;                             % Time axis for plotting.
  
  plot(t, data(i0,:), 'k')
  xlabel('Time [s]')
  title(channel)
  
  spindle_det = spindle_det(i0);

  ax = axis;                                            % Indicate each detection with a red box.
  h0 = (ax(4)-ax(3));
  for k=1:length(spindle_det.startSample)
      %k=1;
      x = spindle_det.startSample(k)/Fs;
      y = ax(3);
      w = spindle_det.endSample(k)/Fs - x;
      h = h0;
      rectangle('Position', [x,y,w,h], 'FaceColor', [1 0 0 0.1], 'LineStyle', 'none')
  end
  

