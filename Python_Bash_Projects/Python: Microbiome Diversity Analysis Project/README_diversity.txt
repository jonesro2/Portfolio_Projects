#RBIF100 Assignment 4 - Interpreting the diversity of the microbiome
#Author: Robin Jones

This purpose of this script is to identify the organisms with the highest and lowest average diversity scores by calculating the mean diversity score with standard deviation for each organism and updating the data table. Distance scores from multidimensional scaling are then plotted for the two organisms with highest diversity scores and the organism with the lowest score.

Note: This script can be customized to any organism data table.  You will need to update the function calls for get_mean_sd, add_data, and plot_animals functions.

Execute the following script inside the week10 directory by copying or typing the following into the command line:
    python3 diversity.py

This script takes the clinical_data.txt,*.diversity.txt and *.distance.txt files as input. * represents any one of 50 animal names.

It will:
  a. calculate mean and standard deviation of the diversity scores for each animal
  b. update the current data table with columns for mean diversity score and it's standard deviation and save as clinical_data.stats.txt
  c. identify the animals with the two highest average diversity scores and the one with the lowest  score
  d. plot the distance matrix scores on a scatterplot for the animals from step c and save them as PDF files

Input files for script are found:
    a. *.diversity.txt -- within diversityScores directory which is found within the inputfiles directory inside the week10 folder
    b. clinical_data.txt -- within inputfiles directory inside the week10 folder
    c. *.distance.txt -- within distanceFiles directory which is found within the inputfiles directory inside the week10 folder

Output locations: 
    a. clinical_data.stats.txt -- week10 folder
    b. 3 distance plots saved as PDF -- week10 folder
