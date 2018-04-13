## Preparation
The Depth of coverage data of a panel will be kept in a database, so the first
step is to create the databases needed to store the information.
***
### Targets of interest
To create the databases you need a BED file with the exact same regions used to run GATK's
DepthOfCoverage tool. To make sure the regions are the same you can create the BED file from
the DoC file.

```bash
$ cut -f 1 Serie00CAPv00.DoC.sample_interval_summary | sed 's/[:-]/\t/gi'  > CAPv00.bed
```

Next annotate the BED file (or combine this step and the previous like [this](https://gist.github.com/zaag/2a327b2fb94c07a6a17992b7ace72017)).

```bash
$ awk -F '\t' '{printf("select DISTINCT \"%s\",\"%s\",\"%s\", name2,strand from refGene  \
where chrom=\"%s\" and  (%s<txEnd and %s>txStart);\n",$1,$2,$3,$1,$2,$3);}' CAPv00.bed |\
mysql --user=genome --host=genome-mysql.cse.ucsc.edu  -A  -N -D hg19 > CAPv00.annotated.bed
```

### Create databases
Use the annotated BED file to create the database.

```bash
$ pyCNV.py --todo create --capture CAPv00 --bedfile CAPv00.annotated.bed
```

Optionally a general database can be created. This is going to be used to store positive controls, samples that should be excluded and genes of particular interest (a-genes).

```bash
$ pyCNV.py --todo create --poscontroles poscons.txt --ingestuurd samples_to_exclude.txt
```

## Analysis

### Fill databases
Obviously the next step is adding the available data to the database.

```bash
$ pyCNV.py --todo analyze --capture CAPv00 --serie Serie00 --new DoC.sample_interval_summary --addonly
```
When adding a DoC file to the database the analysis for that serie will be performed by default.
The addonly flag is optional and adds the data to the database without analysing the added serie.

The analysis can also be started with data already present in the database.

```bash
$ pyCNV --todo analyze --capture CAPv00 --serie Serie00
```
The config.py file can be used to define a output directory. For each new capture a CAPv00 directory will be created during the first analysis for that capture.
