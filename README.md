# electric_signal_analysis

A collection of scripts that analyze recordings of electrical signals from electric fish.

## Signal acquisition
The signal is taken with an electrode, amplified, and recorded in a portable MP3 recorder.

## gymnotSignalAnalysis.R
This is a collection of functions in R built to automatically analyze sound recordings of the electric signal of the pulse-type electric fish <i>Brachyhypopomus</i>. It requires the R package tuneR (by Uwe Ligges) to read in sound recordings.

### Function breakdown: identifyPeaks()
This function takes in a sound recording and identifies all the pulses in it, and tabulates their location (in samples) and amplitude (for both the high amplitude peak and the negative amplitude peak).

### Function breakdown: signalSummaryTable()
For a single recording, the signalSummaryTable() function will get a user-defined amount of pulses in a recording, and identify: their high-amplitude peak ('P1', location and amplitude value), their low amplitude peak ('P2', location and amplitude), the location at which the signal amplitude crosses from positive to negative in between those two peaks, and the boundaries of the pulse (when it begins and when it ends). These values are then used to calculate: the P1:P2 ratio, the total time of a pulse (and the time of all its stages - 1. activation and reaching P1, 2. polarity switch from P1 to P2, and 3. return to 0 from P2), and to determine the decay constant lambda at which the signal goes from P2 (negative amplitude peak) back to 0.
