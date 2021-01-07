getData <- function(outputFolder,
                    connectionDetails,
                    cdmDatabaseSchema,
                    cohortId,
                    outcomeIds,
                    studyStartDate = "20160101",
                    studyEndDate = "",
                    cohortDatabaseSchema,
                    cohortTable = "cohort",
                    cdmVersion = "5",
                    firstExposureOnly = FALSE,
                    washoutPeriod = 0,
                    sampleSize = NULL,
                    analysisId
){
  
  if(!dir.exists(file.path(outputFolder,'data',analysisId))){
    dir.create(file.path(outputFolder,'data',analysisId), recursive = T)
  }
  
  for (value in c(-99999, -730, -365, -180, -90, -30, -14)) {
    
    if(!file.exists( file.path(outputFolder,'data',analysisId, paste0('plpData_lookback_', abs(value)),'cohorts.rds'))){
    covset <- FeatureExtraction::createCovariateSettings(useDemographicsAgeGroup =  TRUE, useDemographicsGender = TRUE,
                                                        useConditionGroupEraLongTerm = TRUE, useDrugGroupEraLongTerm = TRUE,
                                                        useProcedureOccurrenceLongTerm = TRUE, longTermStartDays = value) 
    
    plpdata <- PatientLevelPrediction::getPlpData(connectionDetails = connectionDetails,
                         cdmDatabaseSchema = cdmDatabaseSchema,
                         oracleTempSchema = cdmDatabaseSchema,
                         cohortId =  cohortId,
                         outcomeIds = outcomeIds,
                         studyStartDate = studyStartDate,
                         studyEndDate = studyEndDate,
                         cohortDatabaseSchema = cohortDatabaseSchema,
                         cohortTable = cohortTable,
                         outcomeDatabaseSchema = cohortDatabaseSchema,
                         outcomeTable = cohortTable,
                         cdmVersion = cdmVersion,
                         firstExposureOnly = firstExposureOnly,
                         washoutPeriod = washoutPeriod,
                         sampleSize = sampleSize,
                         covariateSettings = covset)
    
    # add analysis_id to save 
    PatientLevelPrediction::savePlpData(plpdata, file.path(outputFolder,'data',analysisId, paste0('plpData_lookback_', abs(value))))
    }
  }
}

#============================
#step 2 - create the AUCS for each lookback time, eventually will create loop for executing 10 times to get auc distributions
#============================
runModels <- function(outputFolder,
                      outcomeId,
                      repeats = 5,
                      analysisId){
  #============================
  #step 2.a - create the repeats population splits  
  #============================
  # create the same split for all the time periods (picking one randomly)  
  # pick any plpData as they should all have the same patients
  indexFolder <- file.path(outputFolder,'data',analysisId, paste0('plpData_lookback_', abs(730)))
  plpDataNew <- PatientLevelPrediction::loadPlpData(indexFolder)
  
  if(!dir.exists(file.path(outputFolder,'populations',analysisId))){
    dir.create(file.path(outputFolder,'populations',analysisId), recursive = T)
  }
  
  for( i in 1:repeats){
    pops <- PatientLevelPrediction::createStudyPopulation(plpData=plpDataNew, 
                                                          outcomeId = outcomeId, 
                                                          riskWindowStart =  1,
                                                          riskWindowEnd = 365, # is this correct? 
                                                          requireTimeAtRisk = F, 
                                                          removeSubjectsWithPriorOutcome = TRUE, #add additional setting for time for actues
                                                          #for acute we want to remove outcome in past 60 days; requires a boolean?
                                                          binary = TRUE)
    split <- PatientLevelPrediction::subjectSplitter(pops, 
                                                     test=0.2, 
                                                     train=NULL, 
                                                     nfold=3, 
                                                     seed=i)
    trow <- nrow(split)
    allpops <- merge(split, pops, by= 'rowId')  
    if(nrow(allpops) != trow){
      print(paste0('Issue in the code as row counts should be the same ',nrow(allpops),' - ',trow))
    }
    
    saveRDS(allpops[, c('subjectId','cohortStartDate', 'index')], file.path(outputFolder,'populations',analysisId, paste0('pop',i,'.rds')))
  }
  
  #============================
  #step 2.b - Train and evaluate all the models usign the same test/train/fold split
  #============================
  
  runPlpI <- function(x){
    
    
    outputFolder <- x$outputFolder
    n <- x$n
    value <- x$value
    outcomeId <- x$outcomeId
    aId <- x$aId
    
    #plpDataNew <- PatientLevelPrediction::loadPlpData(indexFolder)
    plpDataNew <- PatientLevelPrediction::loadPlpData(file.path(outputFolder,'data',aId, paste0('plpData_lookback_',abs(value))))
    
    # create the population
    pops <- PatientLevelPrediction::createStudyPopulation(plpData = plpDataNew, 
                                                          outcomeId = outcomeId, 
                                                          riskWindowStart =  1,
                                                          riskWindowEnd = 365, # is this correct? 
                                                          requireTimeAtRisk = F, 
                                                          removeSubjectsWithPriorOutcome = TRUE, #for acute we want to remove outcome in past 60 days?
                                                          binary = TRUE)
    
    # load the test/train split so it is the same for each covariate setting 
    splitInfo <- readRDS(file.path(outputFolder,'populations',aId, paste0('pop',n,'.rds')))
    indexes <- merge(pops, splitInfo, by = c('subjectId','cohortStartDate'))[, c('rowId','index')]
    md <- attr(pops, 'metaData')
    pops <- merge(pops, splitInfo[,c('subjectId','cohortStartDate')], 
                  by = c('subjectId','cohortStartDate'))
    attr(pops, 'metaData') <- md
    
    # train the model and save everything
    mod.lr = PatientLevelPrediction::runPlp(population = pops, 
                    plpData = plpDataNew, 
                    minCovariateFraction = 0.001,
                    normalizeData = T,
                    modelSettings = PatientLevelPrediction::setLassoLogisticRegression(seed = 453),
                    indexes = indexes, 
                    saveDirectory = file.path(outputFolder,'results',aId),
                    savePlpData = F,
                    savePlpResult = T,
                    savePlpPlots = T,
                    saveEvaluation = T,
                    verbosity = "INFO",
                    timeStamp = FALSE,
                    analysisId = paste0('plpData_lookback_', abs(value),'_', n)
                    )
    
  } # end runPlpI
  
  
  # create the cluster
  ParallelLogger::logInfo(paste0('Number of cores not specified'))
  cores <- 3 #parallel::detectCores()
  ParallelLogger::logInfo(paste0('Using this many cores ', cores))
  ParallelLogger::logInfo(paste0('Set cores input to use fewer...'))
  
  
  cluster <- ParallelLogger::makeCluster(numberOfThreads = cores)
  ParallelLogger::clusterRequire(cluster, c("PatientLevelPrediction", "Andromeda"))
  
  
  logger <- ParallelLogger::createLogger(name = "PAR",
                                         threshold = 'INFO', 
                                         appenders = list(ParallelLogger::createFileAppender(layout = ParallelLogger::layoutParallel, 
                                                                                             fileName = file.path(outputFolder,'parlog.txt'))))
  ParallelLogger::registerLogger(logger)
  
  for (n in c(1:repeats)){
    
    values <- c(-99999, -730, -365, -180, -90, -30, -14)
    
    getRunSettings <- function(i){
      result <- list(outputFolder = outputFolder,
                     n = n,
                     value = values[i],
                     outcomeId = outcomeId,
                     aId = analysisId)
      return(result)
    }
    
    runSettings <- lapply(1:length(values), getRunSettings)
    
    allResults <- ParallelLogger::clusterApply(cluster = cluster, 
                                               x = runSettings, 
                                               fun = runPlpI, 
                                               stopOnError = FALSE,
                                               progressBar = TRUE)
    
  }
  ParallelLogger::stopCluster(cluster)
}

#============================
#step 3 - Externally validate a model per lookback
#============================
externalValidateModel <- function(connectionDetails,
                                  outputFolder,
                                  validationSchemaTarget,
                                  validationSchemaOutcome,
                                  validationSchemaCdm,
                                  databaseNames,
                                  validationTableTarget,
                                  validationTableOutcome,
                                  validationIdTarget = 16631,
                                  validationIdOutcome = 16474
){
  
  for (value in c(-730, -365, -180, -90, -30, -14))
  {
    plpModelLoc <- file.path(...)
    result <- PatientLevelPrediction::loadPlpResult(plpModelLoc)
    valResults <- PatientLevelPrediction::externalValidatePlp(
      plpResult = result, 
      connectionDetails = connectionDetails,
      validationSchemaTarget = validationSchemaTarget,
      validationSchemaOutcome = validationSchemaOutcome,
      validationSchemaCdm = validationSchemaCdm,
      databaseNames = databaseNames,
      validationTableTarget = validationTableTarget,
      validationTableOutcome = validationTableOutcome ,
      validationIdTarget = validationIdTarget,
      validationIdOutcome = validationIdOutcome
    )
    
    
    saveRDS(valResults, paste0(outputFolder,'external',  '/externalval_', abs(value),'_',n, '.rds'))
  }
}