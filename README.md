## WaveWatchIII (v5.16) source code used in the following research:
_Impact of Shoaling Ocean Surface Waves on Wind stress and Drag Coefficient_

### 1. About this Repository
The **Diagnostic Flux modules (FLD1/2)** originally created by Dr. Reichl have been 
modified for the use in shallower water.


Data in the following two publications are generated with the code uploaded here.
[_1. Chen, X., T. Hara, and I. Ginis, 2020. Impact of Shoaling Ocean Surface Waves on Wind Stress and Drag Coefficient in Coastal Waters: Part I Uniform Wind. J. Geophys. Res._][paper1]

 [_2. Chen, X., I. Ginis, and T. Hara, 2020. Impact of Shoaling Ocean Surface Waves on Wind Stress and Drag Coefficient in Coastal Waters: Part II Tropical Cyclones. J. Geophys. Res._][paper2]
 
### 2. Folders in this Repository 
I. **_ftn_** folder contains the modified WW3v5 code used in the above research. Main changes are made to w3fld1md.ftn, w3fld2md.ftn; other files have only been modified to assist output extra user-defined quantities (flags) associated with the diagnostic flux. 

II. **_inp_** folder contains the ww3_grid.inp, ww3_shel.inp which are required to be modified to be used properly by the modified version of the WaveWatch code. Details of the modifications can be find therein.


### 3. How to use the source code posted here

Download WAVEWATCH III (v5.16) first, then replace the official **_ftn_** with the source code package here for code compilation. 
 
[paper1]: https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2020JC016222
[paper2]: https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2020JC016223



