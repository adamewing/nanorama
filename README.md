# ont-lrm
ONT Live Read Monitor: monitor your adaptive sampling run in real time (or afterwards)

### Install:
```
git clone https://github.com/adamewing/ont-lrm.git
cd ont-lrm
mamba env create -f ont-lrm.yml
conda activate ont-lrm
pip install -e $PWD
```
### Usage:
```
$ lrm_mapper.py -h
usage: lrm_mapper.py [-h] -o OUTDIR [OUTDIR ...] -r REF -t TARGETS -d PLOTDIR [--maponly]

fastq live coverage monitor

options:
  -h, --help            show this help message and exit
  -o OUTDIR [OUTDIR ...], --outdir OUTDIR [OUTDIR ...] minknow output directory containing fastqs, can be called more than once
  -r REF, --ref REF      reference genome .fasta
  -t TARGETS, --targets  TARGETS target .bed
  -d PLOTDIR, --plotdir  PLOTDIR output directory
```
Example:
```
lrm_mapper.py -o path/to/fastq_pass -o path/to/fastq_fail -r hg38.fasta -t RRMS.bed -d example_output
```
Where `-o` points to the directory where fastq files are being deposited by MinKNOW (you'll need basecalling on, fast model is fine). Once it has started processing .fastq files it outputs statistics to the terminal using [plotext](https://github.com/piccolomo/plotext):

![terminal plots](https://github.com/adamewing/ont-lrm/assets/1037202/847241ba-7438-4708-af1e-2b242ec69e0e)


Once at least one fastq has been processed, you can start the web UI:
```
$ lrm_viewer.py -h
usage: lrm_viewer.py [-h] -d PLOTDIR -t TARGETS [--tx TX] [--external]

live coverage monitor aimed at fastqs generated by ONT devices

options:
  -h, --help            show this help message and exit
  -d PLOTDIR, --plotdir PLOTDIR output directory
  -t TARGETS, --targets TARGETS target .bed
  --tx TX               transcripts (or any other annotation) .bed
  --external
  --port
```
The `-d` and `-t` options should match those given to `lrm_mapper.py`. The `--tx` option can be used to show an additional track from a .bed file (optionally BED3+1 with column 4 being an annotation).

A running example of the UI may be accessed [here](http://genome.coffee:8055), if that is inaccessible (down or your local IT policy blocks the port) here's a screenshot:

![screencapture-genome-coffee-8055-2023-06-12-13_48_09](https://github.com/adamewing/ont-lrm/assets/1037202/e5b84649-63a9-4590-8fc1-cc4b060c13ab)


