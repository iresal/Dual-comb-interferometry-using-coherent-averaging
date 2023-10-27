# Dual-comb-interferometry-using-coherent-averaging
Script dedicated to post-processing the data measured in the DCI setup.
The system requires a calibration step using a first measurement as a reference, this data is processed
and the complex amplitude is retrieved. The subsequent measurements are called SIGNAL, the same 
processing is performed but with the additional step of subtracting the reference complex amplitude. 
Import MAT files for both cases
    
    1. REFERENCE. Output of the balanced photodetector (BPD) and 25 MHz signal without a DUT.
    2. SIGNAL. BPD output and 25MHz signal with the DUT (WS, fiber, etc).  
    
Both files were recorded using electro-optic combs, although the same approach can be applied using other platforms such as microcombs (see https://doi.org/10.1038/s42005-023-01424-5)

Created on 2021-06-26
@autor: Israel Rebolledo (israel.rebolledo@ri.se)
