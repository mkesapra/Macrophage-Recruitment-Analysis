#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 10:11:10 2022

@author: kesaprm
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas import DataFrame, Series 

import pims
import trackpy as tp



### Step-2: Locate Features
#f=tp.locate(frames[0],51)
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
##0.5v - thresh = 29, min_mass = 0.8*1e3, tplink = 55, --80.23
##0v - thresh = 29, min_mass = 1.5*1e3, tplink = 55 --81.87
##1v - thresh = 29, min_mass = 5*1e3, tplink = 49  --78.28%
##2v - thresh = 29, min_mass = 5*1e3, tplink = 55 --81.2%
##4v - thresh = 25, min_mass = 5*1e3, tplink = 65 --85.5%




#checking for subpixel accuracy
#tp.subpx_bias(f)

#threshV = 29; min_massV = 1.5*1e3; tplinkV = 55; minDist = 10 ###0V
#threshV = 29; min_massV = 0.8*1e3; tplinkV = 55; minDist = 10 ###0.5V
#threshV = 29; min_massV = 5*1e3; tplinkV = 47; minDist = 10 ###1V
#threshV = 29; min_massV = 5*1e3; tplinkV = 55; minDist = 10 ###2V
#threshV = 25; min_massV = 5*1e3; tplinkV = 65; minDist = 10 ###4V






def cal_dist(threshV,min_massV,tplinkV,minDist,frames) :
    ### Uncomment the below line to plot the image
    plt.imshow(frames[0]);
    
    thresh = threshV
    min_mass = min_massV
    f = tp.locate(frames[0],thresh,minmass = min_mass)
    
    ### Uncomment the below line to plot the segmented image
    tp.annotate(f, frames[0], plot_style={'markersize': 8})
    
    ### Step-4: Locate features in all frames
    f = tp.batch(frames[:], thresh,minmass=min_mass);
    
    ###Step-5,6: Link features into particle trajectories
    # 60 is the no. of pixels moved and memory keeps track of disappeared particles and
    # maintains their ID for up to some number of frames after their last appearance
    t= tp.link(f,tplinkV,memory=1)
    
    print('No. of tracks before filtering:', t['particle'].nunique())
    
    frameDs = {}
    finalxDs = {}
    finaldDs = {}
    tfList = []
    fList = []
    for k in range(0,len(t.frame)):
        if not t.frame[k] in frameDs.keys():
            frameDs[t.frame[k]] = {}
        frameDs[t.frame[k]][t.particle[k]] = [t.x[k],t.y[k]]
    
    for f in range(0,len(frameDs.keys())-1) : 
        if not f in finalxDs.keys():
            finaldDs[f] = {}
            finalxDs[f] = {}
        for p in frameDs[f].keys():
            if p in frameDs[int(f+1)] :
                finaldDs[f][p] = np.sqrt((frameDs[f+1][p][0]-frameDs[f][p][0])**2 + (frameDs[f+1][p][1]-frameDs[f][p][1])**2)
                if (frameDs[f+1][p][0] - frameDs[f][p][0]) >= 0 :
                    finalxDs[f][p] = 1
                else :
                    finalxDs[f][p] = 0
                    
    for f in range(0,len(finalxDs.keys())):
        count = 0
        countNum = 0
        for k in finalxDs[f].keys():
            if finaldDs[f][k] > minDist:
                count = count+finalxDs[f][k]
            countNum = countNum+1
        if (countNum > 0):
            fList.append(int((count/countNum)*100))
        
    return fList
    
mDist = 20
    
### Step-1: Read the data
frames = pims.open('/Users/kesaprm/FY19_20/Spring2020/Project/Images_Darpa/Macrophage galvanotaxis_08_02_22/0V/*.tif')
fList0 = cal_dist(29,1.5*1e3,55,mDist,frames) ###0V
frames = pims.open('/Users/kesaprm/FY19_20/Spring2020/Project/Images_Darpa/Macrophage galvanotaxis_08_02_22/05V/*.tif')
fList05 = cal_dist(29,0.8*1e3,55,mDist,frames) ###05V
frames = pims.open('/Users/kesaprm/FY19_20/Spring2020/Project/Images_Darpa/Macrophage galvanotaxis_08_02_22/1V/*.tif')
fList1 = cal_dist(29,5*1e3,47,mDist,frames) ###1V
frames = pims.open('/Users/kesaprm/FY19_20/Spring2020/Project/Images_Darpa/Macrophage galvanotaxis_08_02_22/2V/*.tif')
fList2 = cal_dist(29,5*1e3,55,mDist,frames) ###2V
frames = pims.open('/Users/kesaprm/FY19_20/Spring2020/Project/Images_Darpa/Macrophage galvanotaxis_08_02_22/4V/*.tif')
fList4 = cal_dist(25,5*1e3,65,mDist,frames) ###4V

fullSnap = [fList0,fList05,fList1,fList2,fList4]
###Step-7: Filter trajectories that exist in all the frames
#t1 = tp.filter_stubs(t, len(frames))
#print('No. of tracks before filtering:', t['particle'].nunique())
#print('No. of tracks after filtering:', t1['particle'].nunique())

# def cal_flow(threshV,min_massV,tplinkV,minDist):
#     frameDs = {}
#     finalxDs = {}
#     finaldDs = {}
#     fList = []
    
#     thresh = threshV
#     min_mass = min_massV
#     f = tp.locate(frames[0],thresh,minmass = min_mass)
    
#     ### Uncomment the below line to plot the segmented image
#     tp.annotate(f, frames[0], plot_style={'markersize': 8})
    
#     ### Step-4: Locate features in all frames
#     f = tp.batch(frames[:], thresh,minmass=min_mass);
    
#     ###Step-5,6: Link features into particle trajectories
#     # 60 is the no. of pixels moved and memory keeps track of disappeared particles and
#     # maintains their ID for up to some number of frames after their last appearance
#     t= tp.link(f,tplinkV,memory=1)
    
#     print('No. of tracks before filtering:', t['particle'].nunique())
    
#     for k in range(0,len(t.frame)):
#         if not t.frame[k] in frameDs.keys():
#             frameDs[t.frame[k]] = {}
#         frameDs[t.frame[k]][t.particle[k]] = [t.x[k],t.y[k]]
    
#     for f in range(0,len(frameDs.keys())-1) : 
#         if not f in finalxDs.keys():
#             finaldDs[f] = {}
#             finalxDs[f] = {}
#         for p in frameDs[f].keys():
#             finaldDs[f][p] = np.sqrt((frameDs[f+1][p][0]-frameDs[f][p][0])**2 + (frameDs[f+1][p][1]-frameDs[f][p][1])**2)
#             if p in frameDs[int(f+1)] :
#                 if (frameDs[f+1][p][0] - frameDs[f][p][0]) >= 0 :
#                     finalxDs[f][p] = 1
#                 else :
#                     finalxDs[f][p] = 0
                    
                    
#     # for f in range(0,len(finalDs.keys())):
#     #     count = 0
#     #     for k in finalDs[f].keys():
#     #         count = count+finalDs[f][k]
#     #     fList.append(count/len(finalDs[f].keys()))
    
#     # return fList
                        
    
# x = cal_flow(29,1.5*1e3,55,1)

































t1 = t

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
for k in t1.particle:#.unique():
    curr_cell_x = t1.x[t1.particle == k]*(1/4.31) #2.155 for 10x, 4.31 for 20x
    curr_cell_y = t1.y[t1.particle == k]*(1/4.31) #2.155 for 10x ,4.31 for 20x
    rel_x = []; rel_y =[];
    for i in range(1 , len(curr_cell_x)-1):#len(frames)-1): 
        rel_x.append(curr_cell_x.iloc[i] - curr_cell_x.iloc[0])
        rel_y.append(curr_cell_y.iloc[i] - curr_cell_y.iloc[0])
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
plt.xlabel(r'x [microns]')
plt.ylabel('y [microns]')
plt.title('Relative trajectories')
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
plt.xlabel("Cells");
plt.ylabel("Directedness")
plt.title("Avg Directedness of the cells")
plt.show()

#find num of cells Directedness>0 in all the frames
# cnt = [0]* (len(directedness)+1);
# num = [0]*(len(directedness)+1);
# #cnt_dir = [];
# for m in range(0, len(directedness)):
#     for n in range(0,len(directedness[m])):
#         num[n]+=1
#         if (directedness[m][n] > 0):
#            cnt[n]+=1 
# cnt_dir = [i / j for i, j in zip(cnt, num)]

        


###Distance calculation

###Distance and speed at each pt calculation
distance = []
speed_eachpt = []
abc = len(all_rel_x[0])-1

for k in range(1, len(all_rel_x)):
    dist_in_loop = []
    for i in range(0, len(all_rel_x[k])-1): #first val is 0, hence the range starts from 1
        x0 = all_rel_x[k][i]
        x1 = all_rel_x[k][i+1]
        y0 = all_rel_y[k][i]
        y1 = all_rel_y[k][i+1]
        dist = np.sqrt((x1-x0)**2 + (y1-y0)**2)
        dist_in_loop.append(dist)
    distance.append(dist_in_loop)
    speed_eachpt.append([number/((k)*5) for number in dist_in_loop])
 
plt.figure()
for j in range(0,  len(all_rel_x)-1):
        plt.plot(range(0,  len(all_rel_x[j+1])-1),speed_eachpt[j])
plt.xlabel(r'time in frames')
plt.ylabel('Speed(um/min)')
plt.ylim(0, 1)
#plt.ylabel('Distance travelled by a cell[microns]')
plt.title('Distance')
plt.show()


plt.figure()
for j in range(0,  len(all_rel_x)-1):
        plt.plot(range(0,  len(all_rel_x[j])),directedness[j])
        print(directedness[j])
plt.xlabel(r'time in frames')
plt.ylabel('directedness(cos(theta))')


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
plt.ylabel('Cell Displacement[microns]')
plt.title('Cell Displacement')
plt.show()


plt.boxplot(displacement)
plt.xticks([1], ['Avg. displacement of the cells'])
# ##Save all the values to export in a dataframe
# d = {'CellNo': range(1,len(all_rel_x)+1),'X': all_rel_x , 'Y': all_rel_y , 'Directedness': directedness}
# df_final =  pd.DataFrame(data=d)

###Step-12: Export the x, y relative values and directedness to a csv
# df_final.to_csv(r'Image_Analysis_vals.csv',index=None)
# print('X,Y values in micrometers and Directedness values are written to the CSV File successfully.')


