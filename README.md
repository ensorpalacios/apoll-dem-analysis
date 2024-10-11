# Analysis of effects of air pollution on dementia in UK Biobank cohort

This repository contains the code for reproducing the analysis in **paper**. The analysis was conducted using the UK Biobank Research Analysis Platform (RAP), and in particular was run from a DNAnexus job within a docker container, containing all the packages and dependencies needed for the analysis.

Steps to setup the analysis from UKB RAP

From within UKB Rap run ttyd command to open a web-based terminal as a DNAnexux job; set as output location the project containing the .dataset file. Advised instance settings are mem1_ssd1_v2_x8 or higher and priority=high. From within the DNAnexus job, donwload the ukbb_analysis folder (including analysis scripts and data) and the docker image:

```dx download ukbb_analysis docker_image.tar.gz```

then load the docker image:

``` docker load -i docker_image.tar.gz```

and finally run the docker container:

``` chmod +x ukbb_analysis/init_container.sh ```

You are now within a container and are ready to run the analysis.

Analysis scripts (folder "scr")

- data_prepare.R: script to clean and recode phenotipic data contained in the /data/raw/data_participant.csv; the output is saved in the /data/processed/df_ukbb.rds.

- data_display.R: script to generate descriptive tables and figures; required output of data_prepare.R, and saves output the /output/descriptive/ folder.

- cox_analysis.R: script to run main analysis; requires output of data_prepare.R; saves output in the /output/results/cox/savedata/ folder.

- cox_display.R: script to visualise results of main analysis; requires output of cox_analysis.R; saves output in the /output/results/cox/display/ folder.

Use ``` ./run.sh``` to run the whole analysis. Alternatively, run each script at a time.


