# HybridMicrobiomeFramework
Data and Code for Camper, Benjamin, et al. "A Conceptual Framework for Host-Associated Microbiomes of Hybrid Organisms." bioRxiv (2023): 2023-05.

https://zenodo.org/records/10265622

The OTU table and associated metadata files for the Kikihia cicadas is found in the Cicadas folder
The OTU table and associated metadata files for the Aspidoscelis lizards is found in the Lizards folder
The OTU table and associated metadata files for the Neotoma woodrats is found in the Woodrats folder
The OTU table and associated metadata fils for the Maize systems is found in the Maize folder

To recreate the figures and data from the cicada,lizard and woodrat systems, use code makefigure2abc.R and makefigure3abc.R
To recreate the figures and data from the maize systems, use code makefigure2def.R and makefigure3def.R
These four scripts will produce quaternary plots similar to the plots in the citataion above. They will also print out summary information consistent with the data in Tables 1-4 in the citation above. Note that the results will not be 'identical' due to differences in the random sampling processes during both OTU table rarefaction and bootstrapping. Both scripts rely on the associated 'HybridMicrobiomes' CRAN package available at: https://cran.r-project.org/web/packages/HybridMicrobiomes/index.html
