# SQUICH Data Processing
Created processed data from SQUICH experiments from the raw fastqs.

## Files included

`20180829_SQUISH14`: Folder containing raw fastqs from SQUICH experiment

`generate_counts.py`: Script to process the raw reads in `20180829_SQUISH14` to get the code count table used for downstream SQUICH analysis

`20180829_SQUISH14_collapsed_codes_align.csv`: Parsed output file of `generate_counts.py` (created with the `--one_tube` flag)

`20180829_SQUISH14_collapsed_codes_align_newcodelabel.csv`: Parsed output file of `generate_counts.py` (created without the `--one_tube` flag)

`all_libraries.csv`: Contains metadata for each experiment corresponding to a fastq file in `20180829_SQUISH14`

`all_codes.csv`: Contains concentrations of codes used in the experiments

`all_degs.csv`: Contains concentrations of degraders used in the experiments

`all_spikes.csv`: Contains concentrations of degraders used in the experiments

`sq_utils.py`: Contains some functions for working with sequence data


## Running the Processing Script

To create the processed output for yourself using the script, clone this repo, cd inside it, and run the `generate_counts.py` script:

```
git clone https://github.com/salzmanlab/SQUICH_Data_Processing.git
cd SQUICH_Data_Processing
python generate_counts.py
```