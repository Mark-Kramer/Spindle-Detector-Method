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
## Preprocessing steps

	Downsample the data:

		d0 = data(k,:);					Get data from one sensor.

		d0 = decimate(d0,d_factor);			Apply `decimate` to downsample the data.
		
Repeat the steps above for all sensors in `data`.

		hdr.info.sfreq = hdr.info.sfreq / d_factor;	Remember to adjust the sampling rate.	

For example, for an original sampling rate of 1024 Hz, set `d_factor = 4`. The new sampling rate is 256 Hz.

----

## To visualize spindles

`LSM_spindle_visualizer(data, hdr, spindle_det, channel)`

![alt text](https://github.com/Mark-Kramer/Spindle-Detector-Method/blob/master/example_spindles.png)

----

## When the code crashes (`Are your data in microvolts?`)

If your code stops, and you encounter this returned text:

`Are your data in microvolts? If not, set options.MinPeakProminence`

Then you must add a **third input** to `LSM_spindle_probabilities`:

1. `spindle_prob = LSM_spindle_probabilities(d, hdr, options);`

where

`options.minPeakProminence` indicates how much the peak must "stand out" to be identified. See [here](https://www.mathworks.com/help/signal/ref/findpeaks.html#buff2uu).
