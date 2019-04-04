# ScatteringHhtPy

ScatteringHhtPy is a tool based on the Hilbert-Huang transform to identify scattering noise in LIGO. The main function is located in the _CheckScatCorr.py_ file and the list of channels in the file named _Channel_List.txt_.

It needs two input variables:
* `observatory` that can be `L1` or `H1` 
* `gps` time

Additionally, the user can enter the optional parameters: 
* `-p primary-channel` which by default is `GDS-CALIB_STRAIN`
* `-d duration` which by default is `100` _note: this feature needs to be fixed_ 
* `-f channel-list-file` which by default is `Channel_List.txt` _note: it changes depending on the epoch_

Example:

```$ ./CheckScatCorr.py L1 1129835218 -p GDS-CALIB_STRAIN -d 100 -f Channel_List.txt```

returns

```Maximum correlation of 0.50 with channel L1:SUS-SRM_M1_DAMP_L_IN1_DQ```
