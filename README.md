# MacrophageRecruitmentAnalysis
Macrophage recruitment index analysis

An Image Analyzer algorithm used in a closed loop controller setup.
The controller gives a current (A) value that gets applied by a current source. The current source
applies the current creating an electric field (EF) to the cells in the chamber. An image is taken
and is evaluated by the image analyzer. In the image analyzer, the cells from images
are identified and tracked for each time step using the Trackpy Python package. The x, y position
values of the cells at every time step are used to calculate the speed and directedness of the
cells. Cells with less than 25 percentile of migration speed are filtered to remove the immobile
cells. Cell directedness is calculated for every 30 min i.e., 6 frames ahead as 

![Screen Shot 2023-02-27 at 11 50 10 AM](https://user-images.githubusercontent.com/89179388/221667924-80eabd25-0c29-4c10-9286-797336f70598.png)

The recruitment index (RI) value is feedback to determine the error closing the loop which is calculated using the directedness values as

![Screen Shot 2023-02-27 at 11 51 08 AM](https://user-images.githubusercontent.com/89179388/221668120-7fd6fcf1-1ce9-4b34-93ae-e42cc137c66d.png)

where

. Cells to the anode are cells with directedness > 0.01

. Cells to the cathode are cells with directedness < -0.01

. Total cells are the sum of the cells to anode, cathode and also the cells with -0.01< directedness < 0.01
