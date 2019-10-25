# Spindle-Detector-Method

The latent state model (LSM) detector runs in two steps:

1. `spindle_prob = LSM_spindle_probabilities(d, hdr);`
2. `spindle_det  = LSM_spindle_detections(spindle_prob);`

where

`d` = the data, [channels x timpepoints]

`hdr` = structure that MUST HAVE:

	`hdr.info.sfreq`      = sampling frequency (Hz).
  
	`hdr.info.ch_names`   = cell of channel names.

Step (1) is slow, and step (2) is fast.

**The method is configured to run on Chu-lab BECTS data.**

----

## When the code crashes

If your code stops, and you encounter this returned text:

`Are your data in microvolts? If not, set options.MinPeakProminence`

Then you must add a **third input** to `LSM_spindle_probabilities`:

1. `spindle_prob = LSM_spindle_probabilities(d, hdr, options);`

where

`options.minPeakProminence` indicates how much the peak must "stand out" to be identified. See [here](https://www.mathworks.com/help/signal/ref/findpeaks.html#buff2uu).
