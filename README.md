# GME_Movement
Movement analysis for Mara/Grumeti data


1. Data Merge
Takes cleaned tracking datasets and merges into a single dataframe
NOTE: Grumeti data is at 30-min intervals!!
#TODO: Add the grumeti metadata in the GR cleaner code
#TODO: Add region names to collar metadata - should be possible with DAS

2. Raster Extraction
Raster extraction for the full dataset. Limited to ecosystem-wide layers. 

3. MovData Filter
Filter out individuals based on overlap with spatial layer extents and time (must overlap at least one crop season)

4. Tactic Clustering

GMM/Tactic Clustering
Calculates ag use stats (mean, 90-day max, etc.)
Applies GMM models to the full data set
Cross validation of model
Option to filter and apply to only a single data set

GMM/GMM_cluster_results
Get cluster cutpoints and classify individual-year tactics
Figures and tables summarizing aggregate and individual-year tactic cluster results

GMM/Ag_Regression_Models
Regression models of ag use in relation to movement and sex/age characteristics

GMM/Tactic_Change_Models
Regression models on tactic change


5. HMM
Fit HMM, evaluate, simulation, and activity budgets

HMM_prepData_pop_GME
Preps data and calculates log speed.
Filters out individuals with high % missing fixes

HMM_population_fit_GME
Fit HMM candidate models

HMM model assess plots & Simulation
Assess the model outputs

HMM_ActivityBudgets_GME
Apply viterbi algorithm to estimate latent states
Calculate activity time budgets and plot (Fig S11)
Calculate activity density bugets and plot (Fig 4)
Conduct overlap tests and apply regression models




