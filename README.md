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

To run other types of data, talk to Cat or Mark first. 



