#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 20:42:27 2022

@author: kesaprm
"""

from matplotlib import pyplot as plt
import pandas as pd
import numpy as np


#xlsx = pd.ExcelFile('/Users/kesaprm/FY19_20/Spring2020/Project/Images_Darpa/Macrophage galvanotaxis_08_02_22/Brandon_Manual_track/Macrophage_galvanotaxis_Summary.xlsx')

xlsx = pd.ExcelFile('/Users/kesaprm/FY19_20/Spring2020/Project/Images_Darpa/Macrophage galvanotaxis_08_02_22/dist_20.xlsx')

df_recIndex = xlsx.parse('Sheet1')
df_recIndex.fillna(0, inplace=True)

norm_0 = []
norm_05 = []
norm_1 = []
norm_2 = []
norm_4 = []


for j in range(0,25):
    norm_0.append((df_recIndex['0v'][j] - df_recIndex['0v'][0])/df_recIndex['0v'][0])
    norm_05.append((df_recIndex['0.5v'][j] - df_recIndex['0.5v'][0])/df_recIndex['0.5v'][0])
    norm_1.append((df_recIndex['1v'][j] - df_recIndex['1v'][0])/df_recIndex['1v'][0])
    norm_2.append((df_recIndex['2v'][j] - df_recIndex['2v'][0])/df_recIndex['2v'][0])
    norm_4.append((df_recIndex['4v'][j] - df_recIndex['4v'][0])/df_recIndex['4v'][0])

norm_4New = []

for j in range(0,36):
    norm_4New.append((df_recIndex['4v_new'][j] - df_recIndex['4v_new'][0])/df_recIndex['4v_new'][0])
    
plt.plot(range(0,6),norm_0[0:26:5], label='0v')
plt.plot(range(0,6),norm_05[0:26:5], label='0.5v')
plt.plot(range(0,6),norm_1[0:26:5], label='1v')
plt.plot(range(0,6),norm_2[0:26:5], label='2v')
plt.plot(range(0,6),norm_4[0:26:5], label='4v')
plt.legend()
plt.xlabel('Time in frames')
plt.ylabel('Number of cells')

    
plt.plot(range(0,25),norm_0, label='0v')
plt.plot(range(0,25),norm_05, label='0.5v')
plt.plot(range(0,25),norm_1, label='1v')
plt.plot(range(0,25),norm_2, label='2v')
plt.plot(range(0,25),norm_4, label='4v')
plt.legend()
plt.xlabel('Time in frames')
plt.ylabel('Number of cells')


plt.plot(range(0,28),np.ceil(df_recIndex['0v']), label='0v')
plt.plot(range(0,28),np.ceil(df_recIndex['0.5v']), label='0.5v')
plt.plot(range(0,28),np.ceil(df_recIndex['1v']), label='1v')
plt.plot(range(0,28),np.ceil(df_recIndex['2v']), label='2v')
plt.plot(range(0,28),np.ceil(df_recIndex['4v']), label='4v')
plt.legend()
plt.xlabel('Time in frames')
plt.ylabel('Number of cells')



avg_n0_list = [0]; avg_n05_list = [0]; avg_n1_list = [0]; avg_n2_list = [0]; avg_n4_list = [0]; 
for k in range(0,28,5):
    avg_n0_list.append(np.mean(df_recIndex['0v'][k:k+5]))
    avg_n05_list.append(np.mean(df_recIndex['0.5v'][k:k+5]))
    avg_n1_list.append(np.mean(df_recIndex['1v'][k:k+5]))
    avg_n2_list.append(np.mean(df_recIndex['2v'][k:k+5]))
    avg_n4_list.append(np.mean(df_recIndex['4v'][k:k+5]))
    
plt.plot(range(0,35,5),[i  for i in avg_n0_list], label='0v')
plt.plot(range(0,35,5),[i  for i in avg_n05_list], label='0.5v')
plt.plot(range(0,35,5),[i  for i in avg_n1_list], label='1v')
plt.plot(range(0,35,5),[i  for i in avg_n2_list], label='2v')
plt.plot(range(0,35,5),[i  for i in avg_n4_list], label='4v')
plt.legend()
plt.axhline(y=30, color='k', linestyle='dashdot')
plt.xlabel('Time in frames')
plt.ylabel('Relative % change in cells')

     

    
plt.plot(range(0,6),df_recIndex['0v'][0:27:5], label='0v')
plt.plot(range(0,6),df_recIndex['0.5v'][0:27:5], label='0.5v')
plt.plot(range(0,6),df_recIndex['1v'][:27:5], label='1v')
plt.plot(range(0,6),df_recIndex['2v'][:27:5], label='2v')
plt.plot(range(0,6),df_recIndex['4v'][:27:5], label='4v')
plt.legend()
plt.xlabel('Time in frames for every 5 frames')
plt.ylabel('Number of cells')


##normalizing instead to cell count in each image
tot_0 = 54; tot_05 = 50; tot_1 = 27; tot_2 = 34; tot_4 = 28; 

n_0 = []
n_05 = []
n_1 = []
n_2 = []
n_4 = []

for j in range(0,27):
    n_0.append((tot_0 - df_recIndex['0v'][j] )/tot_0)
    n_05.append((tot_05 - df_recIndex['0.5v'][j])/tot_05)
    n_1.append((tot_1 - df_recIndex['1v'][j])/tot_1)
    n_2.append((tot_2 - df_recIndex['2v'][j])/tot_2)
    n_4.append((tot_4 - df_recIndex['4v'][j])/tot_4)

plt.plot(range(0,6),n_0[0:27:5], label='0v')
plt.plot(range(0,6),n_05[0:27:5], label='0.5v')
plt.plot(range(0,6),n_1[0:27:5], label='1v')
plt.plot(range(0,6),n_2[0:27:5], label='2v')
plt.plot(range(0,6),n_4[0:27:5], label='4v')
plt.legend()
plt.xlabel('Time in frames')
plt.ylabel('Number of cells')



df_Vel = xlsx.parse('Velocity')
df_Dir = xlsx.parse('Directedness')

#inplace=True means the operation would work on the original object. axis=1 means we are dropping the column, not the row.
df_Vel.drop('Unnamed: 0', inplace=True, axis=1)
df_Dir.drop('Unnamed: 0', inplace=True, axis=1)

#axis=0 - Marks the rows in the dataframe to be deleted
df_Vel.drop(df_Vel.index[218:226], axis=0, inplace=True)
df_Dir.drop(df_Dir.index[218:226], axis=0, inplace=True)

df_Vel.fillna(0, inplace=True)
df_Dir.fillna(0, inplace=True)



df_Dir['0v']*df_Vel['0v']

# initialize data of lists.
data = {'0v': df_Dir['0v']*df_Vel['0v'],
        '0.5v': df_Dir['0.5v']*df_Vel['0.5v'],
        '1v': df_Dir['1v']*df_Vel['1v'],
        '2v': df_Dir['2v']*df_Vel['2v'],
        '4v': df_Dir['4v']*df_Vel['4v'],
        '4v_new': df_Dir['4v_new']*df_Vel['4v_new']}
 
# Create DataFrame
df = pd.DataFrame(data)

# mean of all the cells in each EF
avg_xDist = df[['0v','0.5v','1v', '2v', '4v','4v_new']].mean()


#4v_new compared to rest of the EF
avg4_01 = (avg_xDist['4v_new'] - avg_xDist['0v'])/avg_xDist['0v']
avg4_051 = (avg_xDist['4v_new'] - avg_xDist['0.5v'])/avg_xDist['0.5v']
avg4_11 = (avg_xDist['4v_new'] - avg_xDist['1v'])/avg_xDist['1v']
avg4_21 = (avg_xDist['4v_new'] - avg_xDist['2v'])/avg_xDist['2v']
avg4_4 = (avg_xDist['4v_new'] - avg_xDist['4v'])/avg_xDist['4v']

#4v compared to rest of the EF
avg4_0 = (avg_xDist['4v'] - avg_xDist['0v'])/avg_xDist['0v']
avg4_05 = (avg_xDist['4v'] - avg_xDist['0.5v'])/avg_xDist['0.5v']
avg4_1 = (avg_xDist['4v'] - avg_xDist['1v'])/avg_xDist['1v']
avg4_2 = (avg_xDist['4v'] - avg_xDist['2v'])/avg_xDist['2v']

#2v compared to rest of the EF
avg2_0 = (avg_xDist['2v'] - avg_xDist['0v'])/avg_xDist['0v']
avg2_05 = (avg_xDist['2v'] - avg_xDist['0.5v'])/avg_xDist['0.5v']
avg2_1 = (avg_xDist['2v'] - avg_xDist['1v'])/avg_xDist['1v']

#1v compared to rest of the EF
avg1_0 = (avg_xDist['1v'] - avg_xDist['0v'])/avg_xDist['0v']
avg1_05 = (avg_xDist['1v'] - avg_xDist['0.5v'])/avg_xDist['0.5v']

#0.5v compared to rest of the EF
avg05_0 = (avg_xDist['0.5v'] - avg_xDist['0v'])/avg_xDist['0v']


