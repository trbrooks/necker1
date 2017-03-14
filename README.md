# Compatibility effects in bistable perception of images

**PYTHON FILES (PsychoPy2)**

cong = congruent condition

incong = incongruent condition

tr = training phase

te = test phase

Thus, "incong_te.py" is the incongruent condition of the test phase of the experiment.

**DATA FILES**

In necker_csv_data

Or you can use dataframe.RDA, an R object, and avoid having to use py-to-R7.R

**R Files**

py-to-R7.R uses .csv files to create an R data object. Creates dataframe.RDA from necker_csv_data. Or you can just use dataframe.RDA

necker_analyses_figures.R uses the resulting data object to generate all plots and analyses.

dataframe.R is a dataframe of all data. Load it into an R object called "df1" to run anaylses and generate figures.

**PNG Files**

These images are used by the PsychoPy2 program to create the display.
