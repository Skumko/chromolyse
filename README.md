# chromolyse

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

chromolyse is an R package designed for the analysis of mass parallel sequencing data and the detection of chromothripsis and chromoplexy events. 

## Installation

You can install chromolyse from GitHub using the `devtools` package:

```R
devtools::install_github("Skumko/chromolyse")
```

## Usage

Here's a basic example of how to use chromolyse:

```R
library(chromolyse)
library(dplyr)

cnv_file_path <- "[the_path_to_your_file]" #should contain CNV data in the Nexus or Battenberg format
sv_dataset <- chromolyse::readExcel("./P12.WES.Somatic.Variants.xlsx", sheetname="POST") # an example of Structural Somatic Variation dataset
cnv_dataset <- read.table(cnv_file_path, header = TRUE, sep = '\t')

# modify CNv datset to expected format based on the format
modified_cnv_dataset <- chromolyse::preprocessRawNexusFile(cnv_dataset)
modified_cnv_dataset <- chromolyse::preprocessRawBattenbergFile(cnv_dataset)

# fill missing CNV segments if there are any
filled_cnv_dataset <- chromolyse::fillMissingSegments(modified_cnv_dataset)

#analyse the SV data using the main function
result <- chromolyse::analyseGeneChains(sv_dataset)
# the function creates a circos visualisation by default if any events are detected and puts them all into one graph.

# extract all parts of the result 
graph <- result$graph
identifiedEvents <- result$events
clusteredDataset <- result$clusteredData
mutationData <- result$mutations
allPaths <- result$allPaths

# check number of identified paths/chains
print(length(identifiedEvents))

#visualise the dataset
chromolyse::visualiseCircos(clusteredDataset, filled_cnv_dataset)
# only show ctx links
chromolyse::visualiseCircos(clusteredDataset, filled_cnv_dataset, sv_focus = "ctx")

#visualise some of the retrieved paths, the first one for example
visualiseCircos(clusteredDataset, filled_cnv_dataset, path = identifiedEvents[[1]]$path, sv_focus = "ctx")

# only show chromosomes that the path goes through by extracting chromosome labels from the path itself
chromosomes_to_visualise <- as.numeric(unique(lapply(identifiedEvents[[1]]$path, function(x) strsplit(x, "_")[[1]][1])))
# or use the initial circos visualisation to manually select chromosomes
chromosomes_to_visualise <- c(3,6,10,17)

chromolyse::visualiseCircos(clusteredDataset, filled_cnv_dataset, chromosome_selection = chromosomes_to_visualise)
# or only the CTX edges
chromolyse::visualiseCircos(clusteredDataset, filled_cnv_dataset, chromosome_selection = chromosomes_to_visualise, sv_focus = "ctx")
# we can even combine the path visualisation with chromosome filtering
visualiseCircos(clusteredDataset, filled_cnv_dataset, path = identifiedEvents[[1]]$path, chromosome_selection = chromosomes_to_visualise, sv_focus = "ctx")

#from the path, extract the translocation data
pathInfo <- chromolyse::createDataframeFromEvent(graph, identifiedEvents[[1]])
View(pathInfo)

# or, optimally, the mutation dataset contains all translocations in all identified paths ordered by the rating for analysis
View(mutationData)

# based on this data, we might want to look at a specific path, which is possible using the pathIndex
# let's say that the path 2 is most interesting
pathIndexToExtract <- 2
mutationDataSelection <- mutationData |> filter(pathIndex == pathIndexToExtract)
# then visualise the path, in whatever fashion we like
visualiseCircos(clusteredDataset, filled_cnv_dataset, path = identifiedEvents[[2]]$path, sv_focus = "ctx")

# now we can show only the chromosomes that the path hits - 6, 17, 20
visualiseCircos(clusteredDataset, filled_cnv_dataset, path = identifiedEvents[[2]]$path, chromosome_selection = c(6,17,20), sv_focus = "ctx")

View(mutationDataSelection)

```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

This package was developed as a diploma thesis. Further development from the creator is not expected. For inquiries contact the Faculty of Informatics and Information Technologies of Slovak University of Technology in Bratislava. 