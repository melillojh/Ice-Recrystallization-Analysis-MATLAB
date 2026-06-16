# Ice-Recrystallization-Analysis-MATLAB

MATLAB workflow for quantitative analysis of ice recrystallization from microscopy image sequences.



This software was developed to quantify the evolution of ice crystal populations during recrystallization experiments. The workflow combines Cellpose-based segmentation, crystal tracking, and quantitative morphometric analysis to measure crystal growth, size distributions, and ice coverage over time.



\---



\## Features



The software performs:



\* Interactive region-of-interest (ROI) selection.

\* Automated ice crystal segmentation using Cellpose.

\* Crystal tracking between consecutive frames.

\* Measurement of:



&#x20; \* Crystal area

&#x20; \* Crystal perimeter

&#x20; \* Equivalent radius

&#x20; \* Total ice-covered area

&#x20; \* Mean crystal area

&#x20; \* Median crystal area

&#x20; \* Maximum crystal area

&#x20; \* Mean crystal radius

&#x20; \* Ice area fraction

\* Export of quantitative results as CSV files.

\* Generation of annotated images with scale bars.



\---



\## Spatial Calibration



Image calibration is calculated from the camera pixel size and microscope objective:



```matlab

cameraPixelSize\_um = 4.4;

objectiveMagnification = 10;



um\_per\_pixel = cameraPixelSize\_um / objectiveMagnification;

```



All measurements are reported in micrometers (µm) and square micrometers (µm²).



\---



\## Ice Area Fraction



The software calculates the fraction of the analyzed region occupied by ice:



```text

IceAreaFraction = TotalIceArea / ROIArea

```



For experiments performed with approximately constant sample thickness, the area fraction can be interpreted as an approximation of the ice volume fraction.



\---



\## Crystal Tracking Strategy



The analysis starts from the \*\*last image of the experiment\*\*, which is typically the image containing the largest crystals after recrystallization.



First, all crystals in the last image are detected and assigned a unique identification number (ID).



The software then moves \*\*backwards in time\*\*, analyzing the previous image in the sequence. For each newly detected crystal, the software compares its position (centroid) with the crystal positions identified in the later image. If two crystals are sufficiently close, they are considered part of the same crystal lineage and receive the same tracking ID.



This process is repeated until the first image is reached:



```text

Last Image

&#x20;   ↓

Previous Image

&#x20;   ↓

Previous Image

&#x20;   ↓

...

&#x20;   ↓

First Image

```



Importantly, crystal contours are \*\*not copied\*\* between images. Each image is segmented independently using Cellpose. The information transferred between frames is only the crystal identity, determined from the centroid positions.



This backward-tracking strategy is particularly useful for recrystallization experiments because the final images typically contain fewer and larger crystals that are easier to identify and track throughout the sequence.



\---



\## Workflow



```text

Microscopy Images

&#x20;       ↓

ROI Selection

&#x20;       ↓

Cellpose Segmentation

&#x20;       ↓

Crystal Identification

&#x20;       ↓

Backward Crystal Tracking

&#x20;       ↓

Morphological Analysis

&#x20;       ↓

CSV Export

```



\---



\## Outputs



\### Per-Crystal Measurements



For every detected crystal, the software exports:



\* Area (µm²)

\* Perimeter (µm)

\* Radius (µm)

\* Centroid position

\* Tracking ID



Output file:



```text

all\_crystals\_all\_frames.csv

```



\### Frame Summary



For each frame, the software exports:



\* Total ice area

\* ROI area

\* Ice area fraction

\* Mean crystal area

\* Median crystal area

\* Maximum crystal area

\* Mean crystal radius

\* Radius standard deviation

\* Area standard deviation



Output file:



```text

summary\_per\_frame.csv

```



\### Annotated Images



The software generates annotated images showing:



\* Detected crystals

\* Segmentation overlays

\* Scale bars



These images can be used for visual inspection and quality control.



\---



\## Requirements



The software requires:



\* MATLAB

\* Image Processing Toolbox

\* Medical Imaging Toolbox

\* Cellpose support for MATLAB



\### Installing Cellpose in MATLAB



This workflow uses the MATLAB implementation of Cellpose available through the Medical Imaging Toolbox.



Before running the code:



1\. Open MATLAB.

2\. Go to \*\*Apps → Get More Apps\*\*.

3\. Search for and install:



```text

Medical Imaging Toolbox

```



4\. Verify that the following functions are available:



```matlab

cellpose

segmentCells2D

```



These functions are used for automated crystal segmentation.



More information can be found in the MATLAB documentation:



https://www.mathworks.com/help/medical-imaging/



\### Cellpose Model



The current implementation uses:



```matlab

Model = "cyto"

```



Users may modify the segmentation model and parameters according to their own datasets.



\---



\## Repository Structure



```text

Ice-Recrystallization-Analysis-MATLAB/

│

├── README.md

├── LICENSE

├── CITATION.cff

│

├── src/

│   └── ice\_recrystallization\_analysis.m

│

├── examples/

│

└── results/

```



\---



## Citation

If you use this software or adapt parts of the workflow for your own research, please cite:

Melillo, J.H. (2026).
*Ice Recrystallization Analysis MATLAB* (Version 1.0.0).

https://doi.org/10.5281/zenodo.20717587


\---



\## Author



\*\*Jorge H. Melillo\*\*



Experimental Physicist | Cryophysics Researcher | Image Analysis



ORCID:

https://orcid.org/0000-0001-7642-0368



GitHub:

https://github.com/melillojh



\---



\## License



This project is distributed under the MIT License.



