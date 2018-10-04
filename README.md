# pyCNV

Python copy number variant (pyCNV) detector.

## Installation

The config.py file defines a output directory. For each new capture a CAPv00 directory will be created during the first analysis for that capture. After installation a custom python config file (.py extension) can be given with the --configfile option.

## Preparation
The sample_interval_summary data from GATK's Depth of coverage walker will be kept in a database, so the first step is to create the databases needed to store the information.
To create the necessary tables the panel needs to be in the [diagostic test repository](https://github.com/martinhaagmans/ngstargets)

* Create the tables:

    ```bash
    CNV --create --capture CAPv00
    ```

## Analysis

* Add the  data for every sample in a serie to the database
    ```bash
    CNV --addonly --capture CAPv00 --serie Serie00 --new DoC.sample_interval_summary 
    ```

* Perform QC and analysis for the entire serie

    ```bash
    CNV --capture CAPv00 --serie Serie00 --output /path/to/outputfolder
    ```
