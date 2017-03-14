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

"py-to-R7.R" uses .csv files to create an R data object. Creates an R object called "df1" from necker_csv_data. Or you can just use dataframe.RDA (see below).

dataframe.Rda is a dataframe of all data. Load it into an R object called "df1" to run anaylses and generate figures.

necker_analyses_figures.R uses the dataframe df1 to generate all plots and analyses.


**PNG Files**

These images are used by the PsychoPy2 program to create the display.

"downcube.png" is the disambiguated cube, which can be displayed in either direction.
