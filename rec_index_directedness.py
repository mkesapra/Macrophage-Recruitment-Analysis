#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 20:00:49 2022

@author: kesaprm
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas import DataFrame, Series 

import pims
import trackpy as tp

### Step-1: Read the data
frames = pims.open('/Users/kesaprm/FY19_20/Spring2020/Project/Images_Darpa/Macrophage galvanotaxis_08_02_22/4V/*.tif')


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
##0v - thresh = 29, min_mass = 1.5*1e3, tplink = 55 --81.87
##0.5v - thresh = 29, min_mass = 0.8*1e3, tplink = 55, --80.23
##1v - thresh = 29, min_mass = 5*1e3, tplink = 49  --78.28%
##2v - thresh = 29, min_mass = 5*1e3, tplink = 59 --81.2%
##4v - thresh = 25, min_mass = 5*1e3, tplink = 65 --85.5%

thresh = 25
min_mass = 5*1e3
f = tp.locate(frames[0],thresh,minmass = min_mass)
### Uncomment the below line to plot the segmented image
tp.annotate(f, frames[0], plot_style={'markersize': 8})


#checking for subpixel accuracy
#tp.subpx_bias(f)


### Step-4: Locate features in all frames
f = tp.batch(frames[:], thresh,minmass=min_mass);


###Step-5,6: Link features into particle trajectories
# 60 is the no. of pixels moved and memory keeps track of disappeared particles and
# maintains their ID for up to some number of frames after their last appearance
t= tp.link(f,65)

###Step-7: Filter trajectories that exist in all the frames
t1 = tp.filter_stubs(t, len(frames))
#t1 = t
print('No. of tracks before filtering:', t['particle'].nunique())
print('No. of tracks after filtering:', t1['particle'].nunique())

#t1 = t


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
    for i in range(1 , len(frames)-1): #len(curr_cell_x)-1)
        rel_x.append(curr_cell_x.iloc[i] - curr_cell_x.iloc[0])
        rel_y.append(curr_cell_y.iloc[i] - curr_cell_y.iloc[0])
    all_rel_x.append(rel_x) 
    all_rel_y.append(rel_y)


#all_rel_x  = [i for i in all_rel_x if i != []]

###Step-10: Results - Plot the relative x, y values
# Uncomment the below lines when needed to plot
plt.figure()
num = 0;
for j in range(0,  len(all_rel_x)):
    if(all(i <= 0 for i in all_rel_x[j])):
        plt.plot(all_rel_x[j],all_rel_y[j],'k-', linewidth=0.5)
    else:
        plt.plot(all_rel_x[j],all_rel_y[j],'r-' ,linewidth=0.5)
        num+=1
plt.xlabel(r'x [microns]')
plt.ylabel('y [microns]')
plt.title('Relative trajectories')
plt.show()

##Directedness at each pt calculation
directedness = []
for k in range(0, len(all_rel_x)):
    d_in_loop = []
    for i in range(0, len(all_rel_x[k])): #first val is 0, hence the range starts from 1
        x_val = all_rel_x[k][i]
        y_val = all_rel_y[k][i]
        euc_d = np.sqrt(x_val**2 + y_val**2)
        d_in_loop.append(x_val/ euc_d)
    directedness.append(d_in_loop)
    
###Step-11: Results - Plot the Directedness
avgD = []
for m in range(0, len(directedness)):
    if (len(directedness[m]) != 0):
        avgD.append(np.sum(directedness[m])/len(directedness[m]))

# Bar plot for directedness
for k in range(0, len(all_rel_x)): 
    if (avgD[k] < 0):
        colors = 'k'
    else:
        colors = 'r'
    plt.bar(k,avgD[k], color = colors)   
plt.xlabel("Cells");
plt.ylabel("Directedness")
plt.title("Avg Directedness of the cells")
plt.show()


# Changed directedness calculation as per previous point i.e. cos(theta) = Sum_i=0^n-1 (x_i+1 - xi / sqrt((x_i+1 - xi)^2 + (y_i+1 - yi)^2)) 
directedness_ch = []
for k in range(0, len(all_rel_x)):
    d_in_loop_ch = [0]
    for i in range(0, len(all_rel_x[k])-1): #first val is 0, hence the range starts from 1
        x_i = all_rel_x[k][i]
        y_i = all_rel_y[k][i]
        x_i1 = all_rel_x[k][i+1]
        y_i1 = all_rel_y[k][i+1]
        euc_dd = np.sqrt((x_i1 - x_i)**2 + (y_i1 - y_i)**2)
        d_in_loop_ch.append((x_i1 - x_i)/ euc_dd)
    directedness_ch.append(d_in_loop_ch)

# ## Changed directedness calculation as per previous 5th point i.e. cos(theta) = Sum_i=0^n-1 (x_i+5 - xi / sqrt((x_i+5 - xi)^2 + (y_i+5 - yi)^2)) 
# directedness_ch = []
# for k in range(0, len(all_rel_x)):
#     d_in_loop_ch = []
#     for i in range(0, len(all_rel_x[k])-5,5): #first val is 0, hence the range starts from 1
#         x_i = all_rel_x[k][i]
#         y_i = all_rel_y[k][i]
#         x_i1 = all_rel_x[k][i+5]
#         y_i1 = all_rel_y[k][i+5]
#         euc_dd = np.sqrt((x_i1 - x_i)**2 + (y_i1 - y_i)**2)
#         d_in_loop_ch.append((x_i1 - x_i)/ euc_dd)
#     directedness_ch.append(d_in_loop_ch)



#find avg directedness of all the cells at eacgh point
def Extract(lst,i):
    return [item[i] for item in lst]

print(np.mean(Extract(directedness_ch,0)))

mean_direct_ch = []
for n in range(0,  len(directedness_ch[0])):
    mean_direct_ch.append(np.mean(Extract(directedness_ch,n)))

plt.figure()
for j in range(0,  len(all_rel_x)):
        plt.plot(range(0,  len(all_rel_x[k])),directedness_ch[j])
plt.plot(range(0, 27),mean_direct_ch , color='k', linewidth=3)
plt.xlabel(r'Time in frames',fontSize="14")
plt.ylabel('Cell Directedness',fontSize="14")
plt.title('Cell Directedness',fontSize="16")
plt.show()


##find the directedness>0 cells
###Filling the cell values with 0
for j in range(1,  len(all_rel_x)):
    if (len(directedness[j]) < 27):
         directedness[j].extend([0 for i in range(27 - len(directedness[j]))])
cnt = [0]*27

for j in range(1,  len(all_rel_x)):
    for i in range(0,27):
            if (directedness[j][i] >=0):
                cnt[i] += 1



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



# ##Save all the values to export in a dataframe
# d = {'track': range(1,len(all_rel_x)+1),'x': all_rel_x , 'y': all_rel_y , 'cum_dir': directedness_ch, 'ef': 0, 'set': 2}
# df_final =  pd.DataFrame(data=d)

#### to csv

csv_out =[]
csv_out.append(',track,slice,x,y,cum_dir,ef,set')
ef_val = 400
set_val = 2
col_count =0
for indx1 in range(len(all_rel_x)):
    for indx2 in range(len(all_rel_x[indx1])) :
        csv_out.append(str(col_count)+','+str(indx1+1)+','+str(indx2+1)+','+\
                       str(all_rel_x[indx1][indx2])+','+\
                       str(all_rel_y[indx1][indx2])+','+\
                       str(directedness_ch[indx1][indx2])+','+\
                       str(ef_val)+','+str(set_val))
        col_count = col_count+1
        
csv_out_f = open("csv_out4V.csv","w")
for line in csv_out:
    csv_out_f.write(line+'\n')
    print (line)
            
        
##merge csvs

import os
import glob
import pandas as pd
os.chdir("/Users/kesaprm/Learning/galv")

#Use glob to match the pattern ‘csv’
extension = 'csv'
all_filenames = [i for i in glob.glob('*.{}'.format(extension))]


#combine all files in the list
combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames ])
#export to csv
combined_csv.to_csv( "combined_csv.csv", index=False, encoding='utf-8-sig')




















