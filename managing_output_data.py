#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 10:33:47 2018

@author: paul
"""
# importing libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from os import listdir

parent_folder = "/home/paul/Desktop/MSc DSBA/7. Geometric Methods for Data Analysis/Github/GMDA-CS/"


# loading the data
def save_plots(parent_folder_path):
    output_data_frames = []
    for type_of_tree in ['regular_tree', 'rotated_tree']:
        filenames = listdir(parent_folder_path+ type_of_tree)
        print(filenames)
        
        
        list_of_dataframes = []
        
        for filename in filenames:
            list_of_dataframes.append(pd.read_csv(parent_folder_path +type_of_tree+"/"+filename,
                                                  sep = "|", skiprows = 1,
                                                  names = ["id", "depth", "diameter", "size"]))
        
        
        df = pd.concat(list_of_dataframes, axis = 0)
        
        df = df.drop("id", axis = 1).groupby("depth").agg({"diameter":[np.mean, np.std], "size":[np.mean, np.std]})[:100]
        
        output_data_frames.append(df)
        
        for feature in ['diameter', 'size']:
            fig, ax = plt.subplots(1, figsize = (16,8))
            ax.plot(df.index, df[feature]['mean'], lw=2, label=feature, color='blue')
            ax.fill_between(df.index, df[feature]['mean']+df[feature]['std'], df[feature]['mean']-df[feature]['std'], facecolor='yellow', alpha=0.5)
            
            ax.set_title('evolution of the ' + feature +" - 100 first levels - "+type_of_tree.replace("_", " "))
            ax.legend(loc='upper left')
            ax.set_xlabel('depth')
            ax.set_ylabel(feature)
            ax.grid()
            fig.savefig(parent_folder_path + type_of_tree + "_" + feature + "_evolution.png")
            del fig, ax
    return output_data_frames


output_auslan = save_plots(parent_folder + "results_auslan/")


regular = output_auslan[0]
rotated = output_auslan[1]


plt.figure(figsize = (16,8))
plt.plot(regular.index, regular['diameter']['mean'])
plt.plot(regular.index, rotated['diameter']['mean'])
plt.fill_between(regular.index, regular['diameter']['mean'], rotated['diameter']['mean'], color = 'green', alpha = .25)
plt.legend(['regular tree', 'rotated tree'], loc = "upper right")
plt.grid()
plt.title('Auslan Dataset - diameter evolution comparison')
plt.savefig(parent_folder + "results_auslan/auslan_diff.png")
plt.show()


output_swissroll = save_plots(parent_folder + "/results_swissroll/")

# output_mnist = save_plots(parent_folder + "/results_mnist/")
