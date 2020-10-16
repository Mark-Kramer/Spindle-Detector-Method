# Spindle-Detector-Method

The latent state (LS) detector runs in two steps:

1. `spindle_prob = LSM_spindle_probabilities(d, hdr);`
2. `spindle_det  = LSM_spindle_detections(spindle_prob);`

where

`d` = the data, [channels x time points]

`hdr` = structure that **must have**:

`hdr.info.sfreq`      = sampling frequency (Hz).
  
`hdr.info.ch_names`   = cell of channel names.

Step (1) is slow, step (2) is fast.

----

## To visualize spindles

`LSM_spindle_visualizer(data, hdr, spindle_det, channel)`

![alt text](https://github.com/Mark-Kramer/Spindle-Detector-Method/blob/master/example_spindles.png)

----

## When the code crashes (`Are your data in microvolts?`)

If your code stops, and you encounter this returned text:

`Are your data in microvolts? If not, set options.MinPeakProminence`

Then you must add a **third and fourth input** to `LSM_spindle_probabilities`:

`spindle_prob = LSM_spindle_probabilities(data, hdr, 'MinPeakPromience', X)`

where `X` is a number indicating how much the peak must "stand out" to be identified. See [here](https://www.mathworks.com/help/signal/ref/findpeaks.html#buff2uu).
