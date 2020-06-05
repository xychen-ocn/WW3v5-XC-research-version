### WaveWatchIII (v5.16) source code used in the following research:
_Impact of Shoaling Ocean Surface Waves on Wind stress and Drag Coefficient_

#### 1. About this Repository
The **Diagnostic Flux modules (FLD1/2)** originally created by Dr. Reichl have been 
modified for the use in shallower water.


Data in the following two publications are generated with the code uploaded here.
 _1. Chen, X., T. Hara, and I. Ginis, 2020. Impact of Shoaling Ocean Surface Waves on Wind Stress and Drag Coefficient in Coastal Waters: Part I Uniform Wind. J. Geophys. Res., In review._
 _2. Chen, X., I. Ginis, and T. Hara, 2020. Impact of Shoaling Ocean Surface Waves on Wind Stress and Drag Coefficient in Coastal Waters: Part II Tropical Cyclones. J. Geophys. Res., In press._
 
#### 2. Folders in this Repository 
I. **_ftn_** folder contains the modified WW3v5 code used in the above research. Main changes are made to w3fld1md.ftn, w3fld2md.ftn; other files have only been modified to assist output extra User-defined quantities (flags) associated with the diagnostic flux. 
 

II. **_inp_** folder contains the ww3_grid.inp, ww3_shel.inp which are required to be modified to be used properly by the modified version of the WaveWatch code. Details of the modifications can be find therein.
 

 



