# GME_Movement
Movement analysis for Mara/Grumeti data


1. Data Merge
Takes cleaned tracking datasets and merges into a single dataframe
NOTE: Grumeti data is at 30-min intervals!!
#TODO: Add the grumeti metadata in the GR cleaner code
#TODO: Add region names to collar metadata - should be possible with DAS

2. Raster Extraction
Raster extraction for the full dataset. Limited to ecosystem-wide layers. 
#TODO: Parallelize raster extraction -- requires fixing extents to create a raster stack

3. HWC clustering
Applies GMM models to the full data set
Option to filter and apply to only a single data set
#TODO: Subsample Grumeti data to 1 hour before extracting ag values?

4. HWC temporal use
Generates cummulative ag use graphs (full or single dataset)
Applies smoothing function over data for rolling average of ag use (full or single dataset)
#TODO: Subsample Grumeti data to 1 hour before extracting ag values?

5. HMM
#TODO: Create HMM code to run on the full dataset
#TODO: Run HMM on Grumeti only at 30-min intervals 
