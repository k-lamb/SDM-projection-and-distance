
### Species Distribution Modeling

# Code repository for:
 - Data preparation for modeling
 - Random Forest (RF) and Maximum Entropy (ME) modeling of the complete species range
 - RF and ME sliding window distribution modeling
 - LGM forward- and backward-projection
 - Effective geographic distance

# Important (non-code) files include:
 - iNaturalist_coords.csv -> coordinates of all research-grade iNaturalist observations to May 2022.
 - 'pops rear edge.txt' -> coordinates of sequenced populations from the southern monophyletic cluster.
 - '*_2m30s.csv' -> climate and coordinate data for presence and pseudo-absence to run model code without going through dat_prep.R scripts. 
 - West_seq.txt -> coordinates of sequences populations which effective geographic distances are to be calculated between (between all possible pairs).


## The Pipeline

DATA PREPARATION - dat_prep.R 

The script takes iNaturalist coordinate data and removes outlier points by determining the number of neighboring iNaturalist observations within X distance. The script then down-samples presence data to one point per lat/long grid cell in (X^2 km size). Then, a minimum concave polygon is generating and used to mask the creation of pseudo-absence points far beyong the range margin. Pseudo-absence data is then generated. An equal number of psudo-absence data as presence data exists are generated. The script accounts for points removed in spatial down-sampling. Finally, WorldClim bioclimatic data is downloaded and extracted for all presence and pseudo-absence data points.

Variables to declare:
 - grid size -> kilometer resolution for down-sampling of presence data.
 - cluster size -> presence points can be filtered using the number of neighbors N within X km distance to determine if the point is real or not.
 - buffer distance -> how far to extend the creation of pseudo-absences beyond the current range limit.
 - resolution -> WorldClim bioclimate variable resolution (2.5 recommended).
 - nearest neighbor (T/F) -> whether to use cluster size to discriminate probability that a population is real.
 - presence/absence replication -> if set to 1, the number of pseudo-absence points created will be equal to the number of presence points. Can scale up for smaller subsections of the range or for finer resolution of where observations are known to exist.
 - sanity plot (T/F) -> whether to plot sanity check maps in the course of running the script.

# FINAL - unified (complete) modeling of the species range (/SDM/final)

MAXIMUM ENTROPY - maxent_model.R 

The script data formatted and prepared by the dat_prep.R script, then computes a PCA using the bioclimate data. Then, the number of PC which explain a declared percentage of the data are retained. Next, a training model with 70% data split is run and the type I and II errors are evaluated using the remainder 30%. Finally, the full data set is run and predictions calculated for the species range.

Variables to declare:
 - pca.perc -> cumulative percentage of data explained by PCA used to determine how many PCA to retain as variables in the final model.

RANDOM FOREST - rf_model.R

Similar to the maxent script but, in order to account for variance between random forest runs (which cause ±1% change in accuracy), X number of RF models are averaged for final predictions.

Variables to declare:
 - pca.perc -> cumulative percentage of data explained by PCA used to determine how many PCA to retain as variables in the final model.
 - n.model -> number of RF models to average to account for variance due to seed. 10-100 recommended.


# GRID SEARCH - sliding window approach to SDM to account for regional local adaptation

MODELING - maxent_model.R & rf_model.R 

For descriptions of how the modeling functions in these scripts, see FINAL - *_model.R above. To perform the grid search/sliding window, the data set is subset by kilometer gridding and models performed on each window. Windows are slid along the longitude, then latitude. Each point is sampled at least N^2 times. The highest prediction for each grid cell are then retained. Regions of the range which do not vary in their presence or absence are not modeled (i.e., only presences or only absences exist). 

Variables to declare:
 - pca.perc -> cumulative percentage of data explained by PCA used to determine how many PCA to retain as variables in the final model.
 - n.model -> number of RF models to average to account for variance due to seed. 10-100 recommended.
 - window.fac -> number of discrete windows along the latitudinal and longitudinal plane.
 - pts.rep -> number of times which a point should be replicated in windows. essentially determines how much the window is slid.


# LGM Projection - forward- and backward-projecting using contemporary and ancient climate data

These scripts are intended to map current presence/absence data onto LGM climate data and forward-project where suitability is now using contemporary LGM climate data. The inverse can also be performed, where a contemporary SDM is used to project suitable habitat in the LGM. This script is distinct from others in that it focuses on a monophyletic group at the rear-edge.

DATA PREPARATION - dat_prep_RearEdge_contemporary.R & dat_prep_RearEdge_LGM.R

This script performs data preparation in a similar fashion as described above, but uses a mixed approach to declare buffers. First, a concave buffer is generated around sequenced populations known to be a part of the rear-edge monophyly. iNaturalist presence points are then retained if they are within a given buffer distance of the sequenced populations and are north of a predetermined cutoff where all populations below the cutoff (33.5°N) are members of the monophyly. South of this cut off, all iNaturalist points are retained. The contemporary script downloads contemporary bioclimate data while the LGM script downloads LGM bioclimate data.

Variables to declare:
 - grid size -> kilometer resolution for down-sampling of presence data.
 - cluster size -> presence points can be filtered using the number of neighbors N within X km distance to determine if the point is real or not.
 - buffer distance N -> northern buffer to determine how far iNaturalist presence points should be retained from known sequenced populations.
 - buffer distance PA -> how far to extend the creation of pseudo-absences beyond the current range limit.
 - resolution -> WorldClim bioclimate variable resolution (2.5 recommended).
 - nearest neighbor (T/F) -> whether to use cluster size to discriminate probability that a population is real.
 - presence/absence replication -> if set to 1, the number of pseudo-absence points created will be equal to the number of presence points. Can scale up for smaller subsections of the range or for finer resolution of where observations are known to exist.
 - sanity plot (T/F) -> whether to plot sanity check maps in the course of running the script.
 - size.t -> font size for plots.
 - lat.cutoff -> latitude cut off for where all iNaturalist points south of the cutoff will be retained.
 - lat.cutoff.TF -> whether to use the latitude cutoff at all

MODELING - maxent_model.R & rf_model.R

These scripts follow identical logic as those in FINAL, but forward- and backward-project the SDM. 


# SDM Distance - calculates the least cost path in kilometers through an SDM prediction raster

LEAST COST PATH / EFFECTIVE GEOGRAPHIC DISTANCE - leastcostpath.R 

This script requires a prediction raster for the species range (r), which can be provided as an XYZ data frame stored as a '.csv' (automatically generated by *_model.R scripts above), and a file of points which you want to determine the effective distance between on the SDM ('West_seq.txt'). This script can calculate these distances, and provides a method for mapping the route as a sanity check.
