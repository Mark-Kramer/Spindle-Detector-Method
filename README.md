# Spindle-Detector-Method

This repository contains the spindle detector method used in [Kramer MA, Stoyell SM, Chinappen D, Ostrowski LM, Spencer ER, Morgan AK, Emerton BC, Jing J, Westover MB, Eden UT, Stickgold R, Manoach DS, Chu CJ *Focal Sleep Spindle Deficits Reveal Focal Thalamocortical Dysfunction and Predict Cognitive Deficits in Sleep Activated Developmental Epilepsy.* The Journal of Neuroscience 41, no. 8 (February 24, 2021): 1816–29](https://www.jneurosci.org/content/41/8/1816).

# Versions

- The branch [`master`](https://github.com/Mark-Kramer/Spindle-Detector-Method) contains the most up-to-date, working version of the spinlde detector.

- For the original method applied in 2021 publication, see this branch [`j-neurosci-2021-publication-code`](https://github.com/Mark-Kramer/Spindle-Detector-Method/tree/j-neurosci-2021-publication-code).
 
# Basic use

The latent state (LS) detector runs in two steps:

1. `spindle_prob = LSM_spindle_probabilities(d, hdr);`
2. `spindle_det  = LSM_spindle_detections(spindle_prob);`

where

`d` = the data, [channels x time points]

`hdr` = structure that **must have**:

`hdr.info.sfreq`      = sampling frequency (Hz).
  
`hdr.info.ch_names`   = cell of channel names.

Step (1) is slow, step (2) is fast.

---

# To visualize spindles

`LSM_spindle_visualizer(data, hdr, spindle_det, channel)`

![alt text](https://github.com/Mark-Kramer/Spindle-Detector-Method/blob/master/example_spindles.png)

----

# Example
| Code |  Note |
| --- | --- |
|`load('example_data.mat')`  |  Simulated `data` and `hdr`, available [here](https://github.com/Mark-Kramer/Spindle-Detector-Method/blob/master/example_data.mat).
|`spindle_prob = LSM_spindle_probabilities(data, hdr);`| Compute probabilities of spindles.
|`spindle_det  = LSM_spindle_detections(spindle_prob);`| Convert probabilities to detections.
|`LSM_spindle_visualizer(data, hdr, spindle_det, 'Example')` | Visualize the results.

----

# Advanced use

For narrowband analysis, include a **third input** to `LSM_spindle_probabilities`:

`spindle_prob = LSM_spindle_probabilities(data, hdr, options)`

where

`options.StartFrequency`     = low frequency for narrowband analysis [Hz].
`options.StopFrequency`      = high frequency for narrowband analysis [Hz].

Both options must be specified.  For example, to focus analysis on 9-12 Hz, 

```
options                = []
options.StartFrequency = 9;
options.StopFrequency  = 12;
spindle_prob = LSM_spindle_probabilities(data, hdr, options)
```

----

# When the code warns (`Are your data in microvolts?`)

If the code produces the warning:

`Are your data in microvolts? If not, set options.MinPeakProminence`

Then consider adding a **third input** to `LSM_spindle_probabilities`:

`spindle_prob = LSM_spindle_probabilities(data, hdr, options)`

where `options.MinPeakProminence` is a number indicating how much the peak must "stand out" to be identified. See [here](https://www.mathworks.com/help/signal/ref/findpeaks.html#buff2uu).

The code will run, but results might be meaningless.
