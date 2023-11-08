# SIGMA
Real-Time Peak Flow Prediction Based on Signal Matching

This package contains the R codes and data used for the EMS paper entitled "Real-Time Peak Flow Prediction Based on Signal Matching".

Developed by:

      Dr. Xiuquan (Xander) Wang, P.Eng. 
      Canadian Centre for Climate Change and Adaptation
      University of Prince Edward Island
      Website: http://climatesmartlab.ca
      Email: xiuquan.wang@gmail.com or xxwang@upei.ca


1. This folder contains the following items:

   - SIGMA.r
     ==> this is the R codes for the proposed SIGMA model

   - events_data (folder)
     ==> this subfolder contains all the data used for the case study in the manuscript

   - SIGMA_9years_triangular_comparison.png
     ==> this is the figure showing the comparison between model simulation (based on triangular 
         sub-hourly precipitation pattern and observations)

   - SIGMA_9years_uniform_comparison.png
     ==> this is the figure showing the comparison between model simulation (based on uniform
         sub-hourly precipitation pattern and observations)

   - ReadMe.txt	
     ==> the file you are reading now.

2. To test the R codes, please simply unzip the downloaded zip file. Open R and change the work directory to this unzipped folder
   Then in R, run the code by typing:
   
   - source("SIGMA.r")

   You can open the SIGMA.r to change the sub-hourly precipitation pattern generation method (i.e., uniform or triangular). You should be able to see the results (i.e., those two png files and some outputs files under the events_data folder).

4. Enjoy!


