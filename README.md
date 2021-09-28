# SunPy Tutorial at the Hinode-14/IRIS-11 Meeting

Notebooks for the SunPy tutorial at the Hinode-14/IRIS-11 Meeting

## Outline

* Emphasize three areas of core functionality of sunpy:
  - Download
  - Data structures
  - Coordinates
* Choose a particular (recent) AR that was observed by EIS
  - Is there a recent IHOP that we could use to have synergy with the IRIS tutorial?
  - Ideally this observation would be used in the `eispac` tutorial 
* Pull in the EIS observations, demonstrate capabilities of Map
* Use date from the EIS map to construct query for the following data
  - AIA
  - HMI
  - STEREO-A EUVI
  - SolO EUI
  - GOES (If this AR coincided with a flare?)
* Emphasize how to construct complex, multi-instrument queries
* Highlight use of `sunpy-soar` package to query EUI observations and note that it is a client instantiated outside of `sunpy`
* Show configuration of all of the observatories
* Crop all images to the EIS FOV
* Emphasize that all images are contained in the same data structure
* Do field extrapolation
  - Reproject HMI image into synoptic map
  - Perform field extrapolation
  - Trace field lines
* Emphasize benefits of `SkyCoord`
* Show reprojection of HMI fieldlines onto AIA, EIS, EUVI, and EUI
  - This again emphasizes ease of use of EIS data in `Map` format
  - Also emphasizes use of multiple instruments
* (Optionally: plot GOES time series)
