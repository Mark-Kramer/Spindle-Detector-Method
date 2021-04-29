# Spindle-Detector-Method

This is the ORIGINAL VERSION of the code was used for analysis in the publication [Kramer MA, Stoyell SM, Chinappen D, Ostrowski LM, Spencer ER, Morgan AK, Emerton BC, Jing J, Westover MB, Eden UT, Stickgold R, Manoach DS, Chu CJ *Focal Sleep Spindle Deficits Reveal Focal Thalamocortical Dysfunction and Predict Cognitive Deficits in Sleep Activated Developmental Epilepsy.* The Journal of Neuroscience 41, no. 8 (February 24, 2021): 1816–29](https://www.jneurosci.org/content/41/8/1816).

**Please find bug fixes and improvements on other branches.**

----

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

Then you must add a **third input** to `LSM_spindle_probabilities`:

`spindle_prob = LSM_spindle_probabilities(data, hdr, options)`

where `options.MinPeakProminence` is a number indicating how much the peak must "stand out" to be identified. See [here](https://www.mathworks.com/help/signal/ref/findpeaks.html#buff2uu).
