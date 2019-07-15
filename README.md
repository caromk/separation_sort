# Separation Sort

Blind source separation based spike sorting for ephysiological recordings using large, densely packed multi-electrode arrays. Traces containing multiple units are first "unmixed" using a blind source separation method, Independent Component Analysis (ICA). Then the resulting components are evaluated on how successful the unmixing, or separation, has been. Components determined to have been "well-separated" are 

## NOTE

Code and documention in process of being made more user-friendly for publication and use. Some extraneous scripts remain. A port of the code to Python is planned.

Paper will be posted on bioRxiv soon, a previous version is available in my dissertation: https://hdl.handle.net/2144/19751

## Directories

### SpikeSort

Spike sorting code. 

### Opt

Code for input options. Functionality has been replaced by Matlab's input parser, but code has not all been updated. 
