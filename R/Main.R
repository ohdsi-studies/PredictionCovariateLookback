# Copyright 2020 Observational Health Data Sciences and Informatics
#
# This file is part of PredictionCovariateLookback
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Execute the Study
#'
#' @details
#' This function executes the PredictionCovariateLookback Study.
#' 
#' @param connectionDetails    An object of type \code{connectionDetails} as created using the
#'                             \code{\link[DatabaseConnector]{createConnectionDetails}} function in the
#'                             DatabaseConnector package.
#' @param cdmDatabaseSchema    Schema name where your patient-level data in OMOP CDM format resides.
#'                             Note that for SQL Server, this should include both the database and
#'                             schema name, for example 'cdm_data.dbo'.
#' @param cdmDatabaseName      Shareable name of the database 
#' @param cohortDatabaseSchema Schema name where intermediate data can be stored. You will need to have
#'                             write priviliges in this schema. Note that for SQL Server, this should
#'                             include both the database and schema name, for example 'cdm_data.dbo'.
#' @param cohortTable          The name of the table that will be created in the work database schema.
#'                             This table will hold the target population cohorts used in this
#'                             study.
#' @param oracleTempSchema     Should be used in Oracle to specify a schema where the user has write
#'                             priviliges for storing temporary tables.
#' @param outputFolder         Name of local folder to place results; make sure to use forward slashes
#'                             (/). Do not use a folder on a network drive since this greatly impacts
#'                             performance.
#' @param sampleSize           The number of patients in the target population to sample
#' @param repeats              How many test/train splits to run
#' @param studyStartDate       Only cohort start dates after this date are included
#' @param studyEndDate         Only cohort end dates before this date are included    
#' @param createCohorts        Create the cohortTable table with the target population and outcome cohorts?
#' @param runAnalyses          Run the model development
#' @param externalValdiate     Externally validate the models you developed
#' @param summarizeResults     Summarise the results in the outputFolder for the cdmDatabaseName models  
#' @param minCellCount         The minimum number of subjects contributing to a count before it can be included 
#'                             in packaged results.
#' @param verbosity            Sets the level of the verbosity. If the log level is at or higher in priority than the logger threshold, a message will print. The levels are:
#'                                         \itemize{
#'                                         \item{DEBUG}{Highest verbosity showing all debug statements}
#'                                         \item{TRACE}{Showing information about start and end of steps}
#'                                         \item{INFO}{Show informative information (Default)}
#'                                         \item{WARN}{Show warning messages}
#'                                         \item{ERROR}{Show error messages}
#'                                         \item{FATAL}{Be silent except for fatal errors}
#'                                         }                              
#' @param cdmVersion           The version of the common data model                             
#'
#' @examples
#' \dontrun{
#' connectionDetails <- createConnectionDetails(dbms = "postgresql",
#'                                              user = "joe",
#'                                              password = "secret",
#'                                              server = "myserver")
#'
#' execute(connectionDetails,
#'         cdmDatabaseSchema = "cdm_data",
#'         cdmDatabaseName = 'shareable name of the database'
#'         cohortDatabaseSchema = "study_results",
#'         cohortTable = "cohort",
#'         oracleTempSchema = NULL,
#'         outputFolder = "c:/temp/study_results", 
#'         createCohorts = T,
#'         runAnalyses = T,
#'         minCellCount = 5,
#'         verbosity = "INFO",
#'         cdmVersion = 5)
#' }
#'
#' @export
execute <- function(connectionDetails,
                    cdmDatabaseSchema,
                    cdmDatabaseName = 'friendly database name',
                    cohortDatabaseSchema = cdmDatabaseSchema,
                    cohortTable = "cohort",
                    oracleTempSchema = cohortDatabaseSchema,
                    outputFolder,
                    sampleSize = NULL,
                    repeats = 5,
                    studyStartDate = "20160101",
                    studyEndDate = "",
                    createCohorts = T,
                    runAnalyses = T,
                    externalValidate = F,
                    modelDatabaseName,
                    summarizeResults,
                    minCellCount= 5,
                    verbosity = "INFO",
                    cdmVersion = 5) {
  
  
  if(runAnalyses && externalValidate){
    stop('Must pick development or validation - cant pick both together')
  }
  
  if(externalValidate){
    if(missing(modelDatabaseName)){
      stop('Must specify database to validate: modelDatabaseName')
    }
    outputFolder <- file.path(outputFolder, 'CovariateLookBack', modelDatabaseName)
  } else{
    outputFolder <- file.path(outputFolder, 'CovariateLookBack', cdmDatabaseName)
  }
  
  
  if (!file.exists(outputFolder))
    dir.create(outputFolder, recursive = TRUE)
  
  ParallelLogger::addDefaultFileLogger(file.path(outputFolder, "log.txt"))
  
  
  if (createCohorts) {
    ParallelLogger::logInfo("Creating cohorts")
    createCohorts(connectionDetails = connectionDetails,
                  cdmDatabaseSchema = cdmDatabaseSchema,
                  cohortDatabaseSchema = cohortDatabaseSchema,
                  cohortTable = cohortTable,
                  oracleTempSchema = oracleTempSchema,
                  outputFolder = outputFolder)
  }
  
  if(runAnalyses){
    
    settings <- utils::read.csv(system.file("settings", "settings.csv", package = "PredictionCovariateLookback"))
    
    for(k in 1:nrow(settings)){
    
    ParallelLogger::logInfo(paste0("Extracting Data for setting ", k))
    plpData <- getData(outputFolder = outputFolder,
                       connectionDetails = connectionDetails,
                       cdmDatabaseSchema = cdmDatabaseSchema,
                       cohortId = settings$cohortId[k], 
                       outcomeIds = settings$outcomeId[k], 
                       studyStartDate = studyStartDate, 
                       studyEndDate = studyEndDate, 
                       cohortDatabaseSchema = cohortDatabaseSchema,
                       cohortTable = cohortTable,
                       cdmVersion = cdmVersion,
                       firstExposureOnly = FALSE,
                       washoutPeriod = 0,
                       sampleSize = sampleSize,
                       analysisId = paste0('Analysis_',k)
                       )
    
    ParallelLogger::logInfo(paste0("Developing Models for setting ", k))
    internalResult <- runModels(outputFolder = outputFolder,
                         outcomeId = settings$outcomeId[k],
                         repeats = repeats,
                         analysisId = paste0('Analysis_',k),
                         verbosity = verbosity
    )
    }
    
 
  }
  
  if (externalValidate) {
    
    settings <- utils::read.csv(system.file("settings", "settings.csv", package = "PredictionCovariateLookback"))
    
    for(k in 1:nrow(settings)){
      
      if(dir.exists(file.path(outputFolder,'results', paste0('Analysis_',k)))){
        ParallelLogger::logInfo(paste0("Externally Validating Models in ", paste0('Analysis_',k)))
        
        externalResult <- externalValidateModel(connectionDetails = connectionDetails,
                                                outputFolder = outputFolder,
                                                validationSchemaTarget = cohortDatabaseSchema,
                                                validationSchemaOutcome = cohortDatabaseSchema, 
                                                validationSchemaCdm = cdmDatabaseSchema, 
                                                databaseNames = cdmDatabaseName, 
                                                validationTableTarget = cohortTable, 
                                                validationTableOutcome = cohortTable, 
                                                validationIdTarget = settings$cohortId[1],
                                                validationIdOutcome = settings$outcomeId[1],
                                                analysisId = paste0('Analysis_',k),
                                                oracleTempSchema = oracleTempSchema,
                                                sampleSize = sampleSize
        )
      } else{
        ParallelLogger::logInfo(paste0('Analysis ',k, ' has no results'))
      }
    }
  }
  
  if(summarizeResults){
    getSummary(outputFolder)
  }
  
  
  invisible(NULL)
}




