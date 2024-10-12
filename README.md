# Spindle-Detector-Method

- This repository contains the spindle detector methods used in different publications.
- While the method is the same, the parameter files differ. You must choose a parameter file to use.
- An [example](#example) below indicates how to apply the detector.

---

# Parameter files

| Publication | Parameter File |
| --- | --- |
| [Kramer et al., 2021](https://www.jneurosci.org/content/41/8/1816) | `Original_Spindle_Detector.mat` |
| [Kwon et al., 2023](https://doi.org/10.1093/sleep/zsad017) | `AllAges_Spindle_Detector.mat` |
| [McLaren et al., 2023](https://doi.org/10.1002/acn3.51840) | `ESES_Spindle_Detector.mat` |
| [Berja et al., 2024](https://doi.org/10.1016/j.clinph.2024.08.017) | `Infant_Spindle_Detector.mat` |


----

# Example

- This example applies the spindle detector to 60 s of simulated data, available [here](https://github.com/Mark-Kramer/Spindle-Detector-Method/blob/master/example_data.mat).
- There are 5 spindles near times (10 s, 20 s, 30 s, 40 s, 50 s).
- We use the spindle detector with parameters from [Kwon et al., 2023](https://doi.org/10.1093/sleep/zsad017).


| Code |  Note |
| --- | --- |
|`load('example_data.mat')`  |  Simulated `data` and `hdr`.
| `parameter_file = 'AllAges_Spindle_Detector'` | Choose the parameter file.
|`spindle_prob = LSM_spindle_probabilities(data, hdr, parameter_file);`| Compute probabilities of spindles.
|`spindle_det  = LSM_spindle_detections(spindle_prob);`| Convert probabilities to detections.
|`LSM_spindle_visualizer(data, hdr, spindle_det, 'Example')` | Visualize the results; see [here](#to-visualize-spindles).

---

# Basic use

The latent state model (LSM) spindle detector runs in two steps:

1. `spindle_prob = LSM_spindle_probabilities(data, hdr, parameter_file);`
2. `spindle_det  = LSM_spindle_detections(spindle_prob);`

where

`data` = the data, [channels x time points]

`hdr` = structure that **must have**:

`hdr.info.sfreq`      = sampling frequency (Hz).
  
`hdr.info.ch_names`   = cell of channel names.

`parameter_file`  = specifiy which [parameter file](#parameter-files) to use.

Step (1) is slow, step (2) is fast.

---

# To visualize spindles

`LSM_spindle_visualizer(data, hdr, spindle_det, channel)`

![alt text](https://github.com/Mark-Kramer/Spindle-Detector-Method/blob/master/example_spindles.png)


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
