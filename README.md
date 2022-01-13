# primed-imputation
**Development version**

This repository is for WDL workflows to submit jobs to the NHLBI BioData Catalyst TOPMed and Michigan imputation servers. This is a work in progress and not yet intended for general use.

## Overview

The [Michigan imputation server](https://imputationserver.sph.umich.edu/) and the [TOPMed imputation server](https://imputation.biodatacatalyst.nhlbi.nih.gov/) are cloud instances of the same imputation software with different reference panels. The [imputationbot](https://github.com/lukfor/imputationbot) software allows a user to submit files to either server using the command line. Users must create an account on the server they wish to use, then download an authentication token. VCF files are uploaded to the server, imputed, and the results are available for download for 7 days.

## Docker

The Dockerfile creates a docker image containing the
[imputationbot](https://github.com/lukfor/imputationbot) software. The
image is available on Docker Hub as
[uwgac/primed-imputation](https://hub.docker.com/r/uwgac/primed-imputation).

## WDL

The workflows are written in the Workflow Description Language ([WDL](https://docs.dockstore.org/en/stable/getting-started/getting-started-with-wdl.html)). This GitHub repository contains the Dockerfile, the WDL code, and JSON files containing inputs to each workflow, both for testing and to serve as examples. The bash scripts were used for testing imputationbot and are not used in the workflows.

### Authenticating

imputationbot has a command to register a token interactively; however, this command writes the token to a file in the user’s home directory. The WDLs take the token as an input string and write the configuration file directly, before running any other commands.

### Submitting a job

The user must specify the following inputs:

input | description
--- | ---
token | string with the authentication token (note the example in the JSON file is not a real token)
hostname | URL for either the TOPMed or Michigan server
refpanel | “topmed-r2” is the only option for TOPMed, but there are multiple options for Michigan
population | “all” for TOPMed, multiple options for Michigan
vcf_files | files to impute
password | string that must also be supplied to the results workflow for download. Specifying the password during job submission means the user doesn’t have to rely on receiving the password by email.

When VCF files are submitted to the imputation server, a job_id is
assigned. The submit workflow returns this job_id as an output, and it
must be provided to the results workflow.

### Retrieving results

The user must specify the following inputs:

input | description
--- | ---
token | string with the authentication token (note the example in the JSON file is not a real token)
hostname | URL for either the TOPMed or Michigan server
job_id | string returned by the submission workflow
password | string that must also be supplied to the results workflow for download. Specifying the password during job submission means the user doesn’t have to rely on receiving the password by email.

The imputed genotypes and accompanying files (log, QC report, statistics, md5) are downloaded to the user’s workspace.

![data flow diagram](data_flow_diagram.png)
