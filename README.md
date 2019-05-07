# Spindle-Detector-Method

To run the code:

`[spindle_det] = fun_spindle_detection_double_run(d, ch, Fs, tthresh, amp_factor,'NoSquare',params);`

where

`d` = the data, [channels x timpepoints]

`ch` = cell (1xn) with the channel labels

`Fs` = sampling frequency [Hz]

`tthresh` = time [s] that the signal has to exceed a threshold

`amp_factor` = amplification factor for defining baseline

`'NoSquare'` = do not change this

`params` = structure to specify analysis parameters with these fields,

	`params.fc` = Morlet wavelet center frequency
	`params.fb` = Morelet wavelet bandwidth
	`params.center_freq_desired` = the center frequency desired [Hz]
	`params.threshold_type` = set to `median`
	`params.slowest_spindle_period = the lowest period we allow for a spindle
	`params.spike_threhsold        = the largest peak-to-trough changes we allow for a spindle.

Our recommended settings, which cover spindles from [9.5 to 16.5] Hz, are,

	tthresh                    = 0.4;
	amp_factor                 = 3;
	params.fc                  = 0.5;
	params.fb                  = 5;
	params.center_freq_desired = 13;
	params.threshold_type      = 'median';
	params.spike_threhsold        = 80*1e-6;
	params.slowest_spindle_period = 1/9.5;





