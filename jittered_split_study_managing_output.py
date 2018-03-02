#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 11:44:48 2018

@author: paul
"""

from os import listdir
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

parent_folder = "/home/paul/Desktop/MSc DSBA/7. Geometric Methods for Data Analysis/Github/GMDA-CS/results_swissroll_jittered/"

dataframes = {}

for folder in listdir(parent_folder):
    print(folder)
    tree_specific_dataframes = []
    for filename in listdir(parent_folder +'/'+ folder):
        print(filename)
        tree_specific_dataframes.append(pd.read_csv(parent_folder+'/'+folder+'/'+filename,
                                                  sep = "|", skiprows = 1,
                                                  names = ["id", "depth", "diameter", "size"]))
        df = pd.concat(tree_specific_dataframes, axis = 0)
        df = df.drop("id", axis = 1).groupby("depth").agg({"diameter":[np.mean, np.std], "size":[np.mean, np.std]})[:100]

        dataframes[folder] = df


rotated_dfs = []
regular_dfs = []
for key in dataframes.keys():
    if 'rotated' in key:
        df = pd.DataFrame(dataframes[key]['diameter']['mean'])
        df.columns = [key]
        rotated_dfs.append(df.reset_index())
    else:
        df = pd.DataFrame(dataframes[key]['diameter']['mean'])
        df.columns = [key]
        regular_dfs.append(df.reset_index())

df_rotated = rotated_dfs[0]
for df in rotated_dfs[1:]:
    df_rotated = pd.merge(left = df_rotated, right = df, on = ["depth"], how = 'outer')
    df_rotated['depth'] = df_rotated.index
df_rotated = df_rotated.set_index('depth', drop = True)


df_rotated.drop(["rotated_tree_4"], axis = 1).plot(figsize = (16,8))
plt.show()



df_regular = regular_dfs[0]
for df in regular_dfs[1:]:
    df_regular = pd.merge(left = df_regular, right = df, on = ["depth"], how = 'outer')
    df_regular['depth'] = df_regular.index
df_regular = df_regular.set_index('depth', drop = True)


df_regular.drop(["regular_tree_4"], axis = 1).plot(figsize = (16,8))
plt.show()


df_total = pd.merge(left  = df_regular.drop(['regular_tree_4'], axis = 1).reset_index(), 
                    right = df_rotated.drop(['rotated_tree_4'], axis = 1).reset_index(),
                    on = "depth",
                    how = "outer")
df_total = df_total[sorted(df_total.columns)]

df_total.set_index('depth', drop = True).plot(figsize = (16,8), style = ['-', '-', '-', '--', '--', '--'], 
                   color = ['r', 'orange', 'b', 'r', 'orange', 'b'])
plt.legend(['regular tree - jittered factor = 1',
            'regular tree - jittered factor = 0.5',
            'regular tree - jittered factor = 0.25',
            'rotated tree - jittered factor = 1',
            'rotated tree - jittered factor = 0.5',
            'rotated tree - jittered factor = 0.25'], loc = "upper right")
plt.savefig(parent_folder + 'jittered_results.png')
plt.show()