#!/bin/python3
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import os

filelist=[]
global mean_scores
mean_scores=[]
global score_std
score_std=[]

#loop through each animal file in directory and save in alpha order. calculate mean and st. dev. for each animal's diversity scores and save to variable in order
def get_mean_sd (filepath):
    for file in sorted(glob.glob('{}'.format(filepath))):
        filelist.append(file)
    for file in filelist:
        with open (file,'r') as text:
            scores=text.read().split() #removes \n and makes a list of lines
            flt_scores=[float(score)for score in scores]
            score_array=np.array(flt_scores)
            mean_score=np.mean(score_array)
            scores_std=np.std(score_array)
            mean_scores.append(mean_score)
            score_std.append(scores_std)
    return mean_scores,score_std
get_mean_sd('./inputfiles/diversityScores/*.diversity.txt')

global animals
animals=[]

#read in clinical data file as df, add new averages/sd columns with values from get_mean_sd function, identify animals with two highest averages,1 lowest and save to list
def add_data(filename):
    clinical_data_frame=pd.read_csv('{}'.format(filename), sep='\t')
    clinical_data_frame['Averages']=mean_scores
    clinical_data_frame['Std']=score_std
    clinical_data_frame.to_csv('clinical_data.stats.txt',sep='\t')
    sorted_frame=clinical_data_frame.sort_values('Averages', ascending=False)
    highest_div_average=sorted_frame.iloc[0]['code_name']
    animals.append(highest_div_average)
    second_highest_average=sorted_frame.iloc[1]['code_name']
    animals.append(second_highest_average)
    lowest_div_average=sorted_frame.iloc[-1]['code_name']
    animals.append(lowest_div_average)   
add_data('./inputfiles/clinical_data.txt')

animal_file_list=[]
#find the matching distance files for animals in list, make numpy arrays, plot the distances
def plot_animals(filepath):
    for file in glob.glob('{}'.format(filepath)):
        animal_name=(os.path.basename(file).split('.')[0])
        if animal_name in animals:
            animal_file_list.append(file)
    for animal_file in animal_file_list:
        animal_name=(os.path.basename(animal_file).split('.')[0])
        with open (animal_file,'r') as distance:
            distance_values=distance.read().split()
        x_values=[float(i.split(",")[0])for i in distance_values]
        y_values=[float(i.split(",")[1]) for i in distance_values]
        x_arr=np.array(x_values)
        y_arr=np.array(y_values)
        print(x_arr,y_arr)
        plt.scatter(x_arr,y_arr, c='c', marker='o', s=10)
        plt.xlabel("MDS1")
        plt.ylabel("MDS2")
        plt.suptitle('{} distances'.format(animal_name).title())
        plt.show()
        plt.savefig('{}.pdf'.format(animal_name))
       
        
plot_animals('./inputfiles/distanceFiles/*.distance.txt')
