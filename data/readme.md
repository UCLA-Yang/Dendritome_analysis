- 2466CPneurons_with_standardized_z.csv:<br>the standardized Z (soma-to-surface) value of all 2466 reconstructed MSNs in the P56 dataset
- TME_morpho_box_May-20-2024_2466neurons.csv: <br>the morphometrics and assigned box for all 2466 MSNs in the P56 dataset
- all_htme_brains_with_registration_all.csv: <br>the morphologies and region information of 1296 reconstructed MSNs for the 12m dataset (12m WT MSNs and 12m HD MSNs)
- morpho_2466CPneurons_with_community.csv: <br>the morphologies and region information of 2466 reconstructed MSNs for the P56 dataset (WT)
- Other datasets that are used in the scripts: <br>
-   12m WT and 12m Q140 (HD) MSN morphologies: <br>
all_htme_brains_with_registration_1168CPneurons_onlyCP.csv can be created by all_htme_brains_with_registration_all.csv by limiting Sriatal.Subregion to CPr, CPi, CPc

-   P56 MSN morphologies: <br>
-     morpho_2248CPneurons_with_community.csv can be created by morpho_2466CPneurons_with_community.csv by limiting Striatal.Subregion to CPr, CPi, CPc
-     morpho_1871CPneurons_excluded30.csv can be created by further limiting the standardized_z to >= 0.3
-     morpho_brain_adjusted_box_cluster_xyz_D1D2.csv can be created by running box_analysis.R
