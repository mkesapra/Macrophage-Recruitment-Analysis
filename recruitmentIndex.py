#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 16:29:35 2021

@author: Manasa Kesapragada
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas import DataFrame, Series 

import pims
import trackpy as tp

### Step-1: Read the data
#frames = pims.open('/Users/kesaprm/FY19_20/Spring2020/Project/FS-Tracer_Data/2_2_2022/M2_S3/*.tif')

#frames = pims.open('/Users/kesaprm/FY19_20/Spring2020/Project/Images_Darpa/Switch_Polarity/S2/*.tif')

frames = pims.open('/Users/kesaprm/FY19_20/Spring2020/Project/Images_Darpa/Macrophage galvanotaxis_08_02_22/4v/3/*.tif')
### Uncomment the below line to plot the image
plt.imshow(frames[0]);

### Step-2: Locate Features
f=tp.locate(frames[0],51)
#f.head() #to check the first 5 rows of f
#tp.annotate(f,frames[0])

### Step-3: Refine parameters
### Uncomment the below lines to check the mass histogram
# fig,ax = plt.subplots()
# ax.hist(f['mass'],bins=20)
# # Optionally, label the axes.
# ax.set(xlabel = 'mass', ylabel = 'count')

# Here I am trying to change the minmass looking at the mass-count plot above. 
# Decreasing the minmass value includes more cells, we could find this from the f dataframe
# 51 acts like a threshold 

# Lets assign the finalized threshold and minmass in variables so that they could be used in 
# multiple places
## M1_S4 : thresh = 75, min_mass =1.4*1e4  --121
## M1_S5 : thresh = 85, min_mass =1*1e4 --159
## M2_S3 : thresh = 65, min_mass =1*1e4 --101
## M2_S8 : thresh = 75, min_mass =2.5*1e4 --125
## S2: thresh = 81, min_mass =2.8*1e4, t= tp.link(f,141)
thresh = 43
min_mass =1e4
f = tp.locate(frames[0],thresh,minmass = min_mass)
### Uncomment the below line to plot the segmented image
tp.annotate(f, frames[0], plot_style={'markersize': 9})

#checking for subpixel accuracy
#tp.subpx_bias(f)


### Step-4: Locate features in all frames
f = tp.batch(frames[:], thresh, minmass=min_mass);


###Step-5,6: Link features into particle trajectories
# 60 is the no. of pixels moved and memory keeps track of disappeared particles and
# maintains their ID for up to some number of frames after their last appearance
t= tp.link(f,189)#,memory=1)

###Step-7: Filter trajectories that exist in all the frames
t1 = tp.filter_stubs(t, len(frames))
print('No. of tracks before filtering:', t['particle'].nunique())
print('No. of tracks after filtering:', t1['particle'].nunique())




###Step-8: Results - Trajectories plotting
background = np.mean(frames, axis=0)

### Uncomment the below lines to plot the trajectories
plt.figure()
#superimpose will help to plot the trajectories on the background image
tp.plot_traj(t1,superimpose=background,ax=plt.gca(),  plot_style={'linewidth': 2}); #label=True,
plt.show()

#### Below code is for reference. Commented as not used currently
#plt.rcParams['font.size'] = '8'
# plt.imshow(frames[0])
# plt.quiver(t1.x, t1.y, t1.x, -t1.y, pivot='middle', headwidth=4, headlength=6, color='red')

#plt.xlim(1000, 2000)
#plt.ylim(1000, 2000);

###Step-9: Results - Convert x,y to relative values

#Convert the x, y values from pixels to micrometers -- 1 Pixels= 264.5833 Micrometers
#Correction- Refer Yao-Hui email Dt:Sep-14, the resolution in all images is 4.31 pix per micron, so 1 pixel = 1/4.31 microns
##generate relative x,y-values
all_rel_x = []
all_rel_y = []
for k in t1.particle.unique():
    curr_cell_x = t1.x[t1.particle == k]*(1/4.31) #2.155 for 10x, 4.31 for 20x
    curr_cell_y = t1.y[t1.particle == k]*(1/4.31) #2.155 for 10x ,4.31 for 20x
    rel_x = []; rel_y =[];
    for i in range(0 ,len(frames)): 
        rel_x.append(curr_cell_x[i] - curr_cell_x[0])
        rel_y.append(curr_cell_y[i] - curr_cell_y[0])
    all_rel_x.append(rel_x) 
    all_rel_y.append(rel_y)




###Step-10: Results - Plot the relative x, y values
# Uncomment the below lines when needed to plot
plt.figure()      

for j in range(0,  len(all_rel_x)):
    if(all(i <= 0 for i in all_rel_x[j])):
        plt.plot(all_rel_x[j],all_rel_y[j],'k-', linewidth=0.5)
    else:
        plt.plot(all_rel_x[j],all_rel_y[j],'r-' ,linewidth=0.5)
plt.xlabel(r'x [µm]',fontSize="14" )
plt.ylabel('y [µm]',fontSize="14")
plt.title('Relative trajectories',fontSize="16")
plt.show()

##Directedness calculation
directedness = []
for k in range(0, len(all_rel_x)):
    d_in_loop = []
    for i in range(1, len(all_rel_x[k])): #first val is 0, hence the range starts from 1
        x_val = all_rel_x[k][i]
        y_val = all_rel_y[k][i]
        euc_d = np.sqrt(x_val**2 + y_val**2)
        d_in_loop.append(x_val/ euc_d)
    directedness.append(d_in_loop)
    
###Step-11: Results - Plot the Directedness
avgD = []
for m in range(0, len(directedness)):
    avgD.append(np.sum(directedness[m])/len(directedness[m]))

# Bar plot for directedness
for k in range(0, len(all_rel_x)): 
    if (avgD[k] < 0):
        colors = 'k'
    else:
        colors = 'r'
    plt.bar(k,avgD[k], color = colors)   
plt.xlabel("Cells",fontSize="14");
plt.ylabel("Directedness",fontSize="14")
plt.title("Avg Directedness of the cells",fontSize="16")
plt.show()

## Changed directedness calculation as per previous point i.e. cos(theta) = Sum_i=0^n-1 (x_i+1 - xi / sqrt((x_i+1 - xi)^2 + (y_i+1 - yi)^2)) 
directedness_ch = []
for k in range(0, len(all_rel_x)):
    d_in_loop_ch = []
    for i in range(0, len(all_rel_x[k])-1): #first val is 0, hence the range starts from 1
        x_i = all_rel_x[k][i]
        y_i = all_rel_y[k][i]
        x_i1 = all_rel_x[k][i+1]
        y_i1 = all_rel_y[k][i+1]
        euc_dd = np.sqrt((x_i1 - x_i)**2 + (y_i1 - y_i)**2)
        d_in_loop_ch.append((x_i1 - x_i)/ euc_dd)
    directedness_ch.append(d_in_loop_ch)


#find avg directedness of all the cells at eacgh point
def Extract(lst,i):
    return [item[i] for item in lst]

print(np.mean(Extract(directedness_ch,0)))

mean_direct_ch = []
for n in range(0,  26):
    mean_direct_ch.append(np.mean(Extract(directedness_ch,n)))


plt.figure()
for j in range(0,  len(all_rel_x)):
        plt.plot(range(0,  len(all_rel_x[k])-1),directedness_ch[j])
plt.plot(range(0,  26),mean_direct_ch , color='k', linewidth=3)
plt.xlabel(r'Time in frames',fontSize="14")
plt.ylabel('Cell Directedness',fontSize="14")
plt.title('Cell Directedness',fontSize="16")
plt.show()

##Normal directedness-previous measure

mean_direct = []
for n in range(0,  26):
    mean_direct.append(np.mean(Extract(directedness,n)))


plt.figure()
for j in range(0,  len(all_rel_x)):
        plt.plot(range(0,  len(all_rel_x[k])-1),directedness[j])
plt.plot(range(0,  26),mean_direct , color='k', linewidth=3)
plt.xlabel(r'Time in frames',fontSize="14")
plt.ylabel('Cell Directedness',fontSize="14")
plt.title('Cell Directedness',fontSize="16")
plt.show()


###Step-11: Results - Plot the Directedness change
avgD_ch = []
for m in range(0, len(directedness_ch)):
    avgD_ch.append(np.sum(directedness_ch[m])/len(directedness_ch[m]))

# Bar plot for directedness
for k in range(0, len(all_rel_x)): 
    if (avgD_ch[k] < 0):
        colors = 'k'
    else:
        colors = 'r'
    plt.bar(k,avgD_ch[k], color = colors)   
plt.xlabel("Cells",fontSize="14");
plt.ylabel("Directedness",fontSize="14")
plt.title("Avg Directedness of the cells -Changed Measure",fontSize="16")
plt.show()



###Distance calculation

distance = []
for k in range(0, len(all_rel_x)):
    dist_in_loop = []
    for i in range(0, len(all_rel_x[k])-1): #first val is 0, hence the range starts from 1
        x0 = all_rel_x[k][i]
        x1 = all_rel_x[k][i+1]
        y0 = all_rel_y[k][i]
        y1 = all_rel_y[k][i+1]
        dist = np.sqrt((x1-x0)**2 + (y1-y0)**2)
        dist_in_loop.append(dist)
    distance.append(dist_in_loop)
 
plt.figure()
for j in range(0,  len(all_rel_x)):
        plt.plot(range(0,  len(all_rel_x[k])-1),distance[j])
plt.xlabel(r'time in frames')
plt.ylabel('Distance travelled by a cell[microns]')
plt.title('Distance')
plt.show()

##Displacement calculation

displacement =[]
for k in range(0, len(all_rel_x)):
    x0 = all_rel_x[k][0]
    xN = all_rel_x[k][len(all_rel_x[k])-1]
    y0 = all_rel_y[k][0]
    yN = all_rel_y[k][len(all_rel_x[k])-1]
    disp = np.sqrt((xN-x0)**2 + (yN-y0)**2)
    displacement.append(disp)

plt.figure()

plt.bar(range(0,  len(all_rel_x)),displacement)
plt.xlabel(r'No. of cells')
plt.ylabel('Displacement[µm]')
plt.title('Cell Displacement')
plt.show()


plt.boxplot(displacement)
plt.xticks([1], ['Avg. displacement of the cells'])
##Save all the values to export in a dataframe
d = {'CellNo': range(1,len(all_rel_x)+1),'X': all_rel_x , 'Y': all_rel_y , 'Directedness': directedness, 'Displacement': displacement, 'AvgDirect': avgD}
df_final =  pd.DataFrame(data=d)

##Step-12: Export the x, y relative values and directedness to a csv
df_final.to_csv(r'M2_S8.csv',index=None)
print('X,Y values in micrometers and Directedness values are written to the CSV File successfully.')


