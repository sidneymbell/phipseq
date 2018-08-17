## Installation  
_I'm assuming you already have [python 2.7](https://www.python.org/download/releases/2.7/) and [gitkraken](https://www.gitkraken.com/download) installed. If that is not the case, do that first._
1. Install [nteract](http://nteract.io)  
2. Install the python kernel for nteract by running these two commands in your terminal:  
`python -m pip install ipykernel`  
`python -m ipykernel install --user`  
3. Clone the repository:  
In gitkraken: file > clone repo > clone with url  
Choose where you want the repo to live (e.g., `~/Documents/`).  
Put `https://github.com/sidneymbell/phipseq` in the `URL` box.  
Hit "clone the repo!"  
4. Install required packages  
In the terminal, navigate to wherever you put the repo.  
`cd ~/Documents/phipseq/`  
Run `pip install -r requirements.txt --user` to install the packages we need.  
5. For online/interactive binding maps, set up your plot.ly credentials following the instructions [here](https://plot.ly/python/getting-started/#initialization-for-online-plotting).  
*You only have to do this whole section once and you should be good to go :) *

## Analysis overview  
This set of analyses was designed to enable robust interactive exploration of datasets, rather than automated processing. 
#### Input    
An annotated_counts.xlsx file is required. An example of the correct input format can be found [here](https://github.com/sidneymbell/phipseq/blob/master/experiments/2018-03-23/raw/2018.03.23.annotatedCounts.xlsx). Be careful with your column header format -- this is how downstream analyses identify replicates.  
  
### Analysis steps  
The template analyses can be found in `./scripts/`. Simply duplicate the notebook, move the copy to `./experiments/date/`, and you're ready to get started!  
  
#### Quality control  
**Input:** annotated counts spreadsheet in `./experiments/date/raw/annotated_counts.xlsx`  
**Output:**    
* Visualizations of basic QC metrics in the interactive notebook  
* Basic metadata cleaning  
* A tidy version of the annotated counts in `counts.csv`  
* Dataset with each column normalized to sum to 1, found in `proportions.csv`  
* A list of samples to drop in `drop.txt`  
    
This is the first step in the analytical pipeline, intended to sanity-check replicates and controls (and flag any problems along the way). 
See notebook [here](./scripts/qc-template.ipynb) for step-by-step explanations.  
  
#### Enrichment scores  
**Input:** `drop.txt`, `proportions.csv`  
**Output:**    
* Visualizations of how enrichment scores vary across samples and oligos  
* Each [optionally aggregated by replicate] column converted to fold enrichment over [specified] input, in `enrichment.csv`  
 
This is a short step in the pipeline that mostly handles the aggregation of replicates and the generation of enrichment scores.  
  
#### Exploratory  
**Input:** `enrichment.csv`, [optional: `samples.tsv`]    
**Output:**    
* Calculates sitewise enrichment scores for each virus genome x serum sample in `sitewise-enrichment/`
* Generates static and interactive maps of binding footprints  
* If the optional `samples.tsv` key (which lists `sample-name\tserum` for each sample) is provided, also generates comparisons of enrichment scores for autologous/heterologous epitopes for each serum.  
  
This is where the majority of actual analysis lives, and where most of the plots come from. See notebook for detailed explanations.  
