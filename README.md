Prediction Covariate Lookback
=============

<img src="https://img.shields.io/badge/Study%20Status-Repo%20Created-lightgray.svg" alt="Study Status: Repo Created">

- Analytics use case(s): **Patient-Level Prediction**
- Study type: **Methods Research**
- Tags: **-**
- Study lead: **Jill Hardin**
- Study lead forums tag: **[jill_hardin](https://forums.ohdsi.org/u/jill_hardin)**
- Study start date: **June 2019**
- Study end date: **-**
- Protocol: **-**
- Publications: **-**
- Results explorer: **-**

This study aims to empirically investigate the impact of covariate lookback on internal and external discriminative performance when developing patient-level prediction models.


Instructions To Install and Run Package From Github
===================

- Make sure you have PatientLevelPrediction installed (this requires having Java and the OHDSI FeatureExtraction R package installed):

```r
  # get the latest PatientLevelPrediction
  install.packages("devtools")
  devtools::install_github("OHDSI/PatientLevelPrediction", ref = 'development')
  # check the package
  PatientLevelPrediction::checkPlpInstallation()
```

- Then install the study package:
```r
  # install the network package
  devtools::install_github("ohdsi-studies/PredictionCovariateLookback")
```

- Execute the study by running the code below but make sure to edit the settings:
```r
library(PredictionCovariateLookback)
# USER INPUTS
#=======================

# The folder where the study intermediate and result files will be written:
outputFolder <- "./OhdsiResults"

# how many test/train splits to repeat for model training -it will train this number
# of models per lookback per prediction analysis (default 10)
repeats <- 10

# Details for connecting to the server:
dbms <- "you dbms"
user <- 'your username'
pw <- 'your password'
server <- 'your server'
port <- 'your port'

connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = dbms,
                                                                server = server,
                                                                user = user,
                                                                password = pw,
                                                                port = port)

# Model development settings
# ======================
# Add the database containing the OMOP CDM data
cdmDatabaseSchema <- 'cdm database schema'
cdmDatabaseName <- 'A friendly name for the database name'
# Add a database with read/write access as this is where the cohorts will be generated
cohortDatabaseSchema <- 'work database schema'

oracleTempSchema <- NULL

# table name where the cohorts will be generated
cohortTable <- 'PredictionCovariateLookbackCohort'

# If you need to sample the data for speed (not useful when doing validation as model application is quick)
sampleSize <- NULL


# Model validation settings if running (need another database)
# ======================
runValidation <- F  # make this T to run validation 
cdmDatabaseSchemaV <- 'cdm database schema for validation'
cdmDatabaseNameV <- 'A friendly name for the database name for validation'
cohortDatabaseSchemaV <- 'work database schema for validation'
oracleTempSchemaV <- NULL
cohortTableV <- 'PredictionCovariateLookbackCohortVal'
sampleSizeV <- NULL


# RECOMMENDED TO NOT EDIT BELOW THIS:

# 1) Model developement:
#=======================
execute(connectionDetails = connectionDetails,
        cdmDatabaseSchema = cdmDatabaseSchema,
        cdmDatabaseName = cdmDatabaseName,
        cohortDatabaseSchema = cohortDatabaseSchema,
	oracleTempSchema = oracleTempSchema,
        cohortTable = cohortTable,
        outputFolder = outputFolder,
        sampleSize = sampleSize,
        repeats = repeats
        studyStartDate = "20160101",
        studyEndDate = "",
        createCohorts = T,
        runAnalyses = T,
        externalValidate = F,
        minCellCount= 5,
        verbosity = "TRACE",
        cdmVersion = 5)
        
        
# 2) Model validation:
#=======================  
if(runValidation){
execute(connectionDetails = connectionDetails,
        cdmDatabaseSchema = cdmDatabaseSchemaV,
        cdmDatabaseName = cdmDatabaseNameV,
        cohortDatabaseSchema = cohortDatabaseSchemaV,
        oracleTempSchema = oracleTempSchemaV,
        cohortTable = cohortTableV,
        outputFolder = outputFolder,
        sampleSize = sampleSizeV,
        repeats = repeats,
        studyStartDate = "20160101",
        studyEndDate = "",
        createCohorts = T,
        runAnalyses = F,
        externalValidate = T,
        modelDatabaseName = cdmDatabaseName,      
        minCellCount= 5,
        verbosity = "TRACE",
        cdmVersion = 5)
        }
        
        
# 3) Summarise results:
#=======================  
execute(connectionDetails = connectionDetails,
        cdmDatabaseSchema = cdmDatabaseSchema,
        cdmDatabaseName = cdmDatabaseName,
        cohortDatabaseSchema = cohortDatabaseSchema,
        oracleTempSchema = oracleTempSchema,
        cohortTable = cohortTable,
        outputFolder = outputFolder,
        sampleSize = sampleSize,
        repeats = repeats,
        studyStartDate = "20160101",
        studyEndDate = "",
        createCohorts = F,
        runAnalyses = F,
        externalValidate = F, 
        summarizeResults = T,
        modelDatabaseName = cdmDatabaseName,      
        minCellCount= 5,
        verbosity = "TRACE",
        cdmVersion = 5)


```


# Output

After running the code go to the location you specified as outputFolder - there will be a directory in there called 'CovariateLookBack' and insid that directory you will find a directory with the cdmDatabaseName you specified. Inside this location you should see two summary csv files: 'aucSummary.csv' and 'covCount.csv' in addition to the full data/results.

# Development status
Under development.
