#Metagenomic Analysis folder

The six original datasets may be found in the **Metagenomic-Analysis** repository on this Github page.

Because the ENNB analysis package takes into account the ordering of the samples themselves within each dataset, the outputs in the _resorted outputs_ folder represent the analysis performed after data samples have been rearranged in the dataset. 
All six original datasets have been analyzed using both the EENB with a normalization method (*sigtest* files) and a common dispersion method (*sigtest-method2* files).

##By Watered-Drought status

These datasets were sorted and compared by treatment: either watered conditions or drought conditions. 
Samples in all datasets were resorted based on treatment status (i.e. watered samples are listed first, followed by drought samples).
The outputs from these analyses are CSV files with names following the format *X*.resorted_sigtest(method2).csv

##By City

These datasets were sorted and compared by city: HF, DE, or CA. However, data was then subset by watered/drought treatment; 
therefore, the data was resorted by city first, then treatment within city. For example, all of the HF-watered samples were listed first, then HF-drought, then DE-watered, and so on.

