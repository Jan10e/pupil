## Post hoc pupil analysis with JPEG images in MatLab ##


**Goal**
----

This is a pupil analysis script to extract pupil axis and location, which is used to look at pupil size fluctuations as an indication of state. 


**Experimental set-up**
---
Data is collected from head-fixed mice that run on a wheel while neuronal recordings are collected. The frame rate of the camera is 30Hz and the pupil measurements are normalized before analysis of the pupil axis.



**Process**
-------

The idea is to use edge detection to extract features to estimate the coefficients for the ellipse fitting formula. This ellipse fitting formula is then used to extract the pupil axis and location parameters, which will be used to compare cortical states. 

----------


*author:* <br> &nbsp; Jantine Broek <br><br>
*e-mail:* <br> &nbsp; jantine dot broek at yale dot edu <br><br>
*date:* <br> &nbsp; &nbsp; May 2017 <br><br>
*lab:* <br> &nbsp; &nbsp; McCormick lab, Yale University


----------
