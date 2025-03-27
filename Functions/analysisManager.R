# Analysis manager

# Analysis type. Note: could add individual meta-analysis as a separate type as well.

am.standardAnalysisTypes = function(index = NULL)
{
  out = c("preprocessing", "individual analysis", "WGCNA", "consensus WGCNA", "causal analysis",
          "external analysis");
  if (!is.null(index)) out[index] else out;
}

am.standardObjectTypes = function(analysisType = NULL)
{
  out = list(preprocessing = c(
                  "analysis info",
                  "input data", 
                  "input sample data",
                  "input sample information",
                  "input feature data",
                  "input data mapped to Entrez",
                  "feature data for individual analysis",
                  "principal components of filtered data",
                  "RUV factors of filtered data",
                  "SVA factors of filtered data",

                  "OR data for individual analysis",
                  "OR data for WGCNA",
                  "OR data for consensus individual analysis",
                  "OR data for consensus WGCNA",
                  "OR sample data",
                  "OR plot sample data",
                  "OR sample information",
                  "OR feature data for WGCNA",
                  "OR feature data for individual analysis",
                  "OR principal components of filtered input data",
                  "principal components of filtered and OR data",
                  "RUV factors of filtered and OR data",
                  "SVA factors of filtered and OR data",
                  "OR and QN data for individual analysis",
                  "OR and QN data for WGCNA",

                  "IO removed data for individual analysis",
                  "IO removed data for WGCNA",
                  "IO replaced data for individual analysis",
                  "IO replaced data for WGCNA",
                  "IO replacement weights for data for individual analysis",
                  "IO replacement weight factors for data for individual analysis",
                  
                  "IO removed data for consensus WGCNA",
                  "IO replaced data for consensus WGCNA",
                  "IO removed data for consensus individual analysis",
                  "IO replaced data for consensus individual analysis",

                  "IO removed feature data for individual analysis",
                  "IO removed feature data for WGCNA",
                  "IO replaced feature data for individual analysis",
                  "IO replaced feature data for WGCNA",
                  # Note: IO replaced feature data are the same as the corresponding IO removed feature data
                  # But keep it here so the user doesn't have to remeber this.

                  "principal components of IO replaced data for WGCNA",
                  "principal components of IO removed data for WGCNA",

                  "adjusted IO removed data for individual analysis",
                  "adjusted IO removed data for WGCNA",
                  "adjusted IO removed feature data for individual analysis",
                  "adjusted IO removed feature data for WGCNA",
                  "adjusted IO replaced data for individual analysis",
                  "adjusted IO replaced data for WGCNA",
                  "adjusted IO replaced feature data for individual analysis",
                  "adjusted IO replaced feature data for WGCNA",

                  "OR collapsed data for individual analysis",
                  "OR collapsed data for WGCNA",
                  "OR collapsed feature data for WGCNA",
                  "OR collapsed feature data for individual analysis",
                  "OR and QN collapsed data for individual analysis",
                  "OR and QN collapsed data for WGCNA",
                  "collapsed IO replacement weights for data for individual analysis",
                  "collapsed IO replacement weight factors for data for individual analysis",
                  "IO removed collapsed data for individual analysis",
                  "IO removed collapsed data for WGCNA",
                  "IO removed collapsed feature data for individual analysis",
                  "IO removed collapsed feature data for WGCNA",
                  "IO replaced collapsed data for individual analysis",
                  "IO replaced collapsed data for WGCNA",
                  "IO replaced collapsed feature data for individual analysis",
                  "IO replaced collapsed feature data for WGCNA",

                  "adjusted IO removed collapsed data for individual analysis",
                  "adjusted IO removed collapsed data for WGCNA",
                  "adjusted IO removed collapsed feature data for individual analysis",
                  "adjusted IO removed collapsed feature data for WGCNA",
                  "adjusted IO replaced collapsed data for individual analysis",
                  "adjusted IO replaced collapsed data for WGCNA",
                  "adjusted IO replaced collapsed feature data for individual analysis",
                  "adjusted IO replaced collapsed feature data for WGCNA",

                  "final data for individual analysis",
                  "final data for WGCNA",
                  "final feature data for individual analysis",
                  "final feature data for WGCNA",

                  "final collapsed data for individual analysis",
                  "final collapsed data for WGCNA",
                  "final collapsed feature data for individual analysis",
                  "final collapsed feature data for WGCNA"),
             `individual analysis` = c(
                  "analysis info",
                  "input data",
                  "input data weights",
                  "sample data",
                  "feature data",
                  "design info",
                  "individual results",
                  "combined individual results",
                  "extended combined individual results",
                  "selected genes for individual tests",
                  "selected genes for combined tests",
                  "numbers of significant genes",
                  "numbers of significant genes with top genes",
                  "enrichment analysis",
                  "anRichment collection",
                  "prettify list"),
             WGCNA = c(
                  "analysis info",
                  "input data",
                  "input data weights",
                  "sample data",
                  "feature data",
                  "scale-fee topoly analysis",
                  "individual TOM info",
                  "resampling labels",
                  "blockwise modules",
                  "module labels",
                  "module enrichment analysis",
                  "module enrichment labels",
                  "module eigengenes",
                  "module association",
                  "selected modules, stringent",
                  "selected modules, loose"),
             `consensus WGCNA` = c(
                  "analysis info",
                  "input data",
                  "sample data",
                  "feature data",
                  "scale-fee topoly analysis",
                  "individual TOM info",
                  "consensus TOM info",
                  "resampling labels",
                  "blockwise modules",
                  "module labels",
                  "module enrichment analysis",
                  "module enrichment labels",
                  "module eigengenes",
                  "module association",
                  "selected modules, stringent",
                  "selected modules, loose"),
             `causal analysis` = c(
                  "analysis info",
                  "input data",
                  "sample data",
                  "feature data",
                  "results",
                  "selected results, stringent",
                  "selected results, loose"),
              `external analysis` = c("analysis info",
                  "input data",
                  "sample data",
                  "feature data",
                  "design info",
                  "individual results",
                  "combined individual results",
                  "selected genes for individual tests",
                  "selected genes for combined tests",
                  "numbers of significant genes",
                  "numbers of significant genes with top genes",
                  "enrichment analysis",
                  "anRichment collection",
                  "prettify list",
                  "module labels",
                  "module features"));
  if (!is.null(analysisType))
  {
    index = am.matchAnalysisType(analysisType);
    out = out[[index]];
  }
  unlist(out);
}

am.matchAnalysisType = function(s, punctuation = "-._= ", throw = TRUE)
{
  s2 = gsub(spaste("[", punctuation, "]"), "", tolower(s));
  s3 = gsub("analysis", "", s2);
  out = match(s3, multiGSub(c(spaste("[", punctuation, "]"), "analysis"), c("", ""),
                             tolower( am.standardAnalysisTypes())));
  if (is.na(out) && throw) stop("Argument s: ", s, " matches no known analysis type.");
  out;
}

am.standardizeAnalysisType = function(s, punctuation = "-._= ", throw = TRUE, allowEmpty = TRUE)
{
  if (length(s)==0) 
  {
    if (allowEmpty) 
    {
      return(s)
    } else
      stop("Input analysis type cannot be empty.")
  }
  i = am.matchAnalysisType(s, punctuation = punctuation, throw = throw)
  if (is.finite(i)) am.standardAnalysisTypes(i) else s;
}

am.matchObjectType = function(s, analysisType = NULL, punctuation = "-._= ", throw = TRUE)
{
  if (length(s)!=1) stop("'s' must be a character scalar.");
  s2 = gsub(spaste("[", punctuation, "]"), "", tolower(s));
  
  out = match(s2, gsub(spaste("[", punctuation, "]"), "", 
              tolower(am.standardObjectTypes(analysisType = analysisType))));
  if (is.na(out) && throw) stop("Argument s: ", s, " matches no known standard object type.");
  out;
}

am.standardizeObjectType = function(s, analysisType=NULL, punctuation = "-._= ", throw = TRUE,
                                    allowEmpty = TRUE)
{
  if (length(s)==0) 
  {
    if (allowEmpty)
    {
      return(s)
    } else
      stop("Input object type cannot be empty.")
  }
  if (length(s)>1)
  {
    if (throw) {
      stop("Object type (name) must be a single character string.")
    } else
      return(s);
  }
  index = am.matchObjectType(s, analysisType = analysisType, punctuation = punctuation, throw = throw);
  if (is.na(index)) s else am.standardObjectTypes(analysisType)[index];
}
                   
.am.analysisIDPrefix = "am.obj.";

am.initialize = function(
  storageDir = file.path(om.getBaseDir(), "ObjectManager/AnalysisManager"),
  storageFile = "analysisManager.RData",
  backupDir = file.path(storageDir, "Backups"),
  maxBackups = 1000,
  backupFileBase = "objectManagerBackup-",
  backupFileExtension = ".RData",
  mustExist = FALSE)

{
  if (file.exists(file.path(storageDir, storageFile)))
  {
    amd = loadAsList(file.path(storageDir, storageFile))$.analysisManagerData;
    if (is.null(amd)) 
       stop("File ", storageFile, " does not contain expected object '.analysisManagerData'.");
    if (is.null(amd$analyses))
       stop("Analysis '.analysisManagerData' in file ", storageFile, " does not contain the expected components.");
    amd$storageDir = storageDir;
    amd$storageFile = storageFile;
    amd$backupDir = backupDir;
    amd$backupFileBase = backupFileBase;
    amd$backupFileExtension = backupFileExtension;
    if (!is.null(maxBackups)) amd$maxBackups = maxBackups;

    amd$haveLock = FALSE;
    .analysisManagerData <<- amd;
  } else {
    if (mustExist) stop("File ", storageFile, " does not exist");
    suppressWarnings(dir.create(storageDir, recursive = TRUE));
    .analysisManagerData <<- list(
         analyses = list(), 
         storageDir = storageDir,
         storageFile = storageFile,
         backupDir = backupDir,
         storageFile = storageFile,
         backupFileBase = backupFileBase,
         backupFileExtension = backupFileExtension,
         maxBackups = if (!is.null(maxBackups)) maxBackups else 1000,
         lockFile = "00-analysisManagerLock.RData",
         haveLock = FALSE);
    am.saveManagerData(FORCE = TRUE);
  }
}

am.backupDir = function() .analysisManagerData$backupDir;
am.maxBackups = function() .analysisManagerData$maxBackups;
am.backupFileBase = function() .analysisManagerData$backupFileBase;
am.backupFileExtension = function() .analysisManagerData$backupFileExtension;
am.storageDir = function() .analysisManagerData$storageDir;
am.storageFile = function() .analysisManagerData$storageFile;
am.lockFile = function() .analysisManagerData$lockFile;
am.haveLock = function() .analysisManagerData$haveLock;

#==========================================================================================================
#
# Disk operations
#
#==========================================================================================================

am.reread = function(haveLock = am.haveLock())
{
  amd = loadAsList(file.path(am.storageDir(), am.storageFile()))$.analysisManagerData;
  amd$haveLock = haveLock;
  .analysisManagerData <<- amd;
  amd;
}

am.haveLock = function()
{
  .analysisManagerData$haveLock;
}

am.lock = function(timeout = 30, sleepTime = 0.25)
{
  maxAttempts = timeout/sleepTime
  att = 1;
  pid = Sys.getpid()
  lockFile = file.path(am.storageDir(), am.lockFile());
  while (att<=maxAttempts && file.exists(lockFile)) {att = att + 1; Sys.sleep(sleepTime)};
  if (file.exists(lockFile))
    stop("Could not get a lock in ", am.storageDir(), ". Please check the file ", am.lockFile(), ".");
  on.exit(unlink(lockFile));
  lock = file(lockFile, open = "wt");
  writeLines(as.character(pid), con = lock);
  close(lock);
  on.exit(NULL);
  .analysisManagerData$haveLock<<-TRUE
}

am.lockAndRead = function(timeout = 10, sleepTime = 0.5)
{
  am.lock(timeout, sleepTime);
  am.reread();
}

am.checkLock = function()
{
  if (!am.haveLock()) return(FALSE);
  lockFile = file.path(am.storageDir(), am.lockFile());
  if (!file.exists(lockFile))
    stop("Could not find a lock file, aborting.");

  lock = file(lockFile, open = "rt");
  pid = readLines(con = lock, n=1, ok = FALSE, warn = FALSE);
  close(lock);

  if (pid!=as.character(Sys.getpid()))
    stop("Lock file exists but belongs to a different process, aborting.");

  return(TRUE);
}


am.unlock = function()
{
  if (am.checkLock())
  {
    unlink(file.path(am.storageDir(), am.lockFile()));
    .analysisManagerData$haveLock <<- FALSE;
  }
}

am.write = function(minBackupTimeDifference = 12*3600)
{
  if (am.checkLock()) {
     tempF = tempfile(tmpdir = am.storageDir(), fileext = ".RData");
     on.exit(try(suppressWarnings(file.remove(tempF)), silent = TRUE));
     saveFile = file.path(am.storageDir(), am.storageFile());
     save(.analysisManagerData, file = file.path(tempF))
     if (filesDiffer(tempF, saveFile))
     {
        am.backup(minTimeDifference = minBackupTimeDifference);
        file.rename(tempF, saveFile);
     }
  } else
    stop("Cannot write: do not have lock.");
}

am.backup = function(minTimeDifference = 12*3600)
{
  if (am.checkLock())
  {
    if (!file.exists(file.path(am.storageDir(), am.storageFile()))) return(NULL);
    suppressWarnings(dir.create(am.backupDir(), recursive = TRUE));
    bf = list.files(am.backupDir(),
               pattern = spaste("^", am.backupFileBase(), "[0-9]*", am.backupFileExtension()));
    if (length(bf)> 0)
    {
      mtimes = file.mtime(file.path(am.backupDir(), bf));
      lastBackupTime = max(mtimes);
      currentTime = Sys.time();
      if (currentTime - minTimeDifference < lastBackupTime ) return(NULL);
      numbers = replaceMissing(as.numeric(multiSub(
                     c(am.backupFileBase(), am.backupFileExtension()), c("", ""), bf, fixed = TRUE)));
      newNumber = max(numbers) + 1;
    } else
      newNumber = 1;

    newFile = spaste(am.backupFileBase(), prependZeros(newNumber, 8), am.backupFileExtension());
    file.copy(from = file.path(am.storageDir(), am.storageFile()),
              to = file.path(am.backupDir(), newFile), copy.mode = TRUE, copy.date = FALSE);
    if (length(bf)>=am.maxBackups())
    {
      firstInd = which.min(mtimes)
      file.remove(file.path(am.backupDir(), bf[firstInd]));
    }
  } else
    stop("Cannot backup: do not have the lock.");
}
    
am.writeAndUnlock = function(minBackupTimeDifference = 12*3600)
{
  tempF = tempfile(tmpdir = am.storageDir(), fileext = ".RData");
  saveFile = file.path(am.storageDir(), am.storageFile());
  on.exit(try(suppressWarnings(file.remove(tempF)), silent = TRUE));
  if (am.checkLock()) 
  {
    save(.analysisManagerData, file = file.path(tempF))
    if (filesDiffer(tempF, saveFile))
    {
      am.backup(minTimeDifference = minBackupTimeDifference);
      file.rename(tempF, saveFile);
    } 
    am.unlock();
    invisible(.analysisManagerData);
  } else
    stop("Do not have lock: cannot write and unlock.");
}

am.clearAll = function()
{
  printFlush(spaste(
      "Warning: this function de-synchronizes object manager from saved data.\n",
      "Data on disk are not removed."));
  .analysisManagerData$analyses <<- list(); 
  invisible(.analysisManagerData);
}

am.saveManagerData = function(
      storageDir = .analysisManagerData$storageDir,
      storageFile = .analysisManagerData$storageFile,
      ...)
{
 args = list(...);
  if ("FORCE" %in% names(args))
  {
    suppressWarnings(dir.create(storageDir, recursive = TRUE));
    save(.analysisManagerData, file = file.path(storageDir, storageFile));
  } else
    stop("This function is for emergencies only.");
}

#=====================================================================================================
#
# Existing name vectors and their corresponding IDs.
#
#=====================================================================================================

am.existingAnalysisIDs = function(reread = !am.haveLock())
{
  if (reread) am.reread();
  sapply(.analysisManagerData$analyses, getElement, "ID");
}

am.existingAnalysisNames = function(reread = !am.haveLock())
{
  if (reread) am.reread();
  sapply(.analysisManagerData$analyses, getElement, "name");
}

am.nAnalyses = function(reread = !am.haveLock())
{
  if (reread) am.reread();
  length(.analysisManagerData$analyses);
}

am.getNewAnalysisID = function(reread = !am.haveLock())
{
  if (reread) am.reread();
  ex.n = as.numeric(sub(.am.analysisIDPrefix, "", am.existingAnalysisIDs(), fixed= TRUE));
  out = suppressWarnings(max(ex.n, na.rm = TRUE) + 1);
  if (!is.finite(out)) out = 1;
  spaste(.am.analysisIDPrefix, out);
}

#========================================================================================================
#
# newAnalysis
#
#========================================================================================================

# For now make it a flat list. In the future could have a separate list for each sub-analysis.
am.newAnalysis = function(
   name,
   description,
   analysisType = character(0),
   organism,
   tissue,
   shortName = name,
   prettyName = name,
   metaData = list(),
   dataSource,
   allowNonstandardType = FALSE,
   updateAnalysisManager = FALSE)
{

  analysisType = am.standardizeAnalysisType(analysisType, throw = !allowNonstandardType)
  
  out = list(ID = character(0),
      name = name,
      description = description,
      analysisType = analysisType,
      organism = organism,
      tissue = tissue,
      dataSource = dataSource,
      shortName = shortName,
      prettyName = prettyName,
      metaData = metaData,
      objectIDs = character(0));
  class(out) = c("AM.Analysis", class(out));

  if (updateAnalysisManager)
    out$ID = am.insertAnalysis(out, stopOnDuplicate = TRUE);

  out;
}

am.insertAnalysis = function(analysis, stopOnDuplicate = FALSE, verbose = 1, indent = 0)
{
  if (!am.haveLock())
  {
    am.lockAndRead();
    on.exit(am.unlock());
    doWrite = TRUE;
  } else 
    doWrite = FALSE

  spaces = indentSpaces(indent);
  amd = .analysisManagerData;

  if (!inherits(analysis, "AM.Analysis"))
    stop("'analysis' must be of class 'AM.Analysis'.");

  # try finding the name among existing analyses
  pos = am.analysisIndex(analysis$name, analysis$ID);
  if (is.na(pos))
  {
    if (length(analysis$ID)==1) stop("'analysis' contains invalid ID. Remove the invalid ID and re-run.");
    if (verbose) printFlush(paste(spaces, "Adding new analysis to analysis manager."));
    analysis$ID = am.getNewAnalysisID();
    amd$analyses = c(.analysisManagerData$analyses, list(analysis));
  } else {
    if (stopOnDuplicate) 
      stop("Analysis of the same name or ID already exists.");
    if (length(analysis$ID) > 0  && !all(analysis$nameVector==amd$analyses[[pos]]$nameVector))
       stop("ID in given 'analysis' exists in analysis manager but has different names.");
    if (verbose) printFlush(paste(spaces, "Replacing existing analysis in analysis manager."));
    analysis$ID = amd$analyses[[pos]]$ID;
    amd$analyses[[pos]] = analysis;
  }
  .analysisManagerData <<- amd
  if (doWrite) am.write();
  analysis$ID;
}

## Note: this is almost equivalent to am.newAnalysis(..., updateAnalysisManager = TRUE) except this returns the analysis
## ID instead of the actual analysis.
am.insertNewAnalysis = function(...) am.insertAnalysis(am.newAnalysis(...));

#==========================================================================================================
#
# Removing an analysis from analysis (and object) manager
#
#==========================================================================================================

am.renameAnalysis = function(name, newName, renameObjects = TRUE, dryRun = FALSE)
{
  if (!am.haveLock())
  {
    am.lockAndRead();
    on.exit(am.unlock());
    doWrite = TRUE;
  } else
    doWrite = FALSE

  ind = am.analysisIndex(name = name, throw = TRUE);

  ind2 = am.analysisIndex(name = newName, throw = FALSE);
  if (is.finite(ind2)) 
    stop("'newName' already exists as a name of an analysis in the analysis manager. ");

  amd = .analysisManagerData;
  amd$analyses[[ind]]$name = newName;
  if (!dryRun) .analysisManagerData <<- amd

  if (renameObjects)
    om.renameObjects(oldName = name, newName = newName, matchComponents = 1, dryRun = dryRun);

  if (doWrite) am.write();

  TRUE;
  
}


am.removeAnalysisFromManager = function(ID, removeObjectsFromManager = TRUE)
{
  if (!am.haveLock())
  {
    am.lockAndRead();
    on.exit(am.unlock());
    doWrite = TRUE;
  } else
    doWrite = FALSE
  index = am.analysisIndex(ID = ID, throw = TRUE)
  analysis = am.getAnalysis(ID);
  if (removeObjectsFromManager && length(analysis$objectIDs) > 0)
  {
    # Check that all object IDs exist...
    if (!om.haveLock())
    {
       om.lockAndRead();
       on.exit(om.writeAndUnlock(), add = TRUE);
    }
    inManager = analysis$objectIDs %in% om.existingObjectIDs();
    if (any(inManager))
    {
      xx = om.removeObjectsFromManager(analysis$objectIDs[inManager], dryRun = TRUE);
      if (all(xx))
        om.removeObjectsFromManager(analysis$objectIDs[inManager], dryRun = FALSE);
    }
  }
  .analysisManagerData$analyses <<- .analysisManagerData$analyses[-index];
  if (doWrite) am.write();
  return(TRUE);
}

am.matchObjectsToAnalyses = function(IDs)
{
  existingIDs = lapply(.analysisManagerData$analyses, getElement, "objectIDs");
  analysisIDs = sapply(.analysisManagerData$analyses, getElement, "ID");
  existingIDs.x = do.call(rbind, mymapply(function(.objIDs, .anaID, .index) 
                        data.frame(objectIDs = .objIDs, 
                                   analysisID = rep(.anaID, length(.objIDs)), 
                                   analysisIndex = rep(.index, length(.objIDs))),
                                   existingIDs, analysisIDs, 1:length(existingIDs)));
  data.frame(inputID = IDs, existingIDs.x[match(IDs, existingIDs.x$objectIDs), ]);
}
  
am.removeObjectsFromManagers = function(IDs, removeFromObjectManager = TRUE, dryRun = TRUE)
{
  if (!am.haveLock())
  {
    am.lockAndRead();
    on.exit(am.unlock());
    doWrite = TRUE;
  } else
    doWrite = FALSE
  
  if (removeFromObjectManager)
  {
    rmd.om = om.removeObjectsFromManager(ID = IDs, dryRun = dryRun);
    if (any(!rmd.om))
      stop("Some given object IDs are not valied object manager IDs. Please check and try again.");

  }

  anaInfo = am.matchObjectsToAnalyses(IDs = IDs);
  anaInfo = anaInfo[!is.na(anaInfo$analysisIndex), ];
  modifyAnalyses = unique(anaInfo$analysisIndex);
  for (ana in modifyAnalyses)
  {
    objIDs = .analysisManagerData$analyses[[ana]]$objectIDs
    objIDs = setdiff(objIDs, anaInfo$objectIDs[ anaInfo$analysisIndex==ana]);
    if (!dryRun)
        .analysisManagerData$analyses[[ana]]$objectIDs <<- objIDs;
  }

  if (doWrite) am.write();
  invisible(anaInfo);
}

# clean a component of an analysis. Returns the name vectors of removed objects.

am.cleanAnalysisComponent = function(name, ID = am.analysisID(name), component, removeFromObjectManager = TRUE, dryRun = TRUE)
{
  objIDs = am.allObjectIDsInAnalysis(analysisID = ID, analysisType = component);
  if (length(objIDs)==0)
  {
    warning("am.cleanAnalysisComponent: no objects in specified component ", component);
    return(NULL);
  }
  nameVectors = om.existingNameVectors(objIDs);
  am.removeObjectsFromManagers(objIDs, removeFromObjectManager = removeFromObjectManager, dryRun = dryRun);
  invisible(nameVectors);
}

am.existingAnalysisComponents = function(name, ID = am.analysisID(name))
{
  objIDs = am.allObjectIDsInAnalysis(analysisID = ID);
  nameVectors = om.existingNameVectors(objIDs);
  unique(sapply(nameVectors, `[[`, 2));
}

    
#==========================================================================================================
#
# Searching for analyses
#
#==========================================================================================================

am.analysisIndex = function(name = NULL, ID = NULL, noMatch = NA, throw = FALSE)
{
  if (!am.haveLock())
  {
    am.lockAndRead();
    on.exit(am.unlock());
  }
  if (length(ID) > 0)
  {
    ids = am.existingAnalysisIDs();
    out = match(ID, ids, nomatch = noMatch)
    if (throw && any(is.na(out)))
      stop("Some entries in 'ID' are not among existing IDs.");
    return (out);
  }
  if (length(name) > 0)
  {
    exNames = am.existingAnalysisNames();
    out = match(name, exNames, nomatch = noMatch);
    if (throw && any(is.na(out)))
      stop("Some entries in 'name' are not among existing names.");
    return (out);
  }
  stop("'name' or 'ID' must be specified.");
}

am.analysisID = function(name, noMatch = NA, throw = FALSE)
{
  if (!am.haveLock())
  {
    am.lockAndRead();
    on.exit(am.unlock());
  }

  n = length(name);
  ind = am.analysisIndex(name = name, noMatch = NA, throw = throw);
  out = rep(noMatch, n);
  exist = is.finite(ind)
  if (any(exist))
    out[exist] = sapply(.analysisManagerData$analyses[ind[exist]], getElement, "ID");
  out;
}

am.analysisIndex.fromTags = function(tags,
      searchType = c("any", "all"),
      searchComponents = c("name", "description", "organism", "tissue", "dataSource"),
      invert = FALSE,
      exact = TRUE,
      fixed = TRUE,
      ignoreCase = TRUE,
      ignorePunctuation = exact,
      punctuation = "-._= ")
{
  if (!am.haveLock())
  {
    am.lockAndRead();
    on.exit(am.unlock());
  }

  searchType = match.arg(searchType);

  validComponents = c("name", "description", "organism", "tissue", "dataSource");

  if (any(!searchComponents %in% validComponents))
    stop(formatStrings(spaste("'searchComponents' must contain one or more of ", 
            paste(validComponents, collapse = ", ")), maxCharPerLine = 70, capitalMultiplier =1 ));

  if (ignoreCase)
  {
     tags = tolower(tags);
  }

  if (ignorePunctuation)
  {
    pattern = spaste("[", punctuation, "]");
    tags = gsub(pattern, "", tags);
  }

  tagMatch = lapply(tags, function(t)
  {
    out = do.call(cbind, lapply(searchComponents, function(comp)
    {
      searchText = lapply(.analysisManagerData$analyses, getElement, comp);
      if (ignoreCase)
        searchText = tolower(searchText);
      if (ignorePunctuation)
        searchText = gsub(pattern, "", searchText);
      if (exact)
      {
        match = sapply(searchText, function(x) any(x==t))
      } else
        match = sapply(searchText, function(x) length( grep(t, x, fixed = fixed) ) > 0);
    }));
    which(rowSums(out + 0) > 0);
  });

  if (searchType=="any") {
    match = multiUnion(tagMatch)
  } else
    match = multiIntersect(tagMatch);

  if (invert) match = setdiff(1:am.nAnalyses(), match);
  return(match)
}


am.analysisID.fromTags = function(...)  # Same arguments as am.analysisIndex.fromTags
{
  am.existingAnalysisIDs()[am.analysisIndex.fromTags(...)];
}

am.analysisNames.fromTags = function(...)  # Same arguments as am.analysisIndex.fromTags
{
  am.existingAnalysisNames()[am.analysisIndex.fromTags(...)];
}

am.analysisNames.fromIDs = function(IDs, throw = TRUE)
{
  am.existingAnalysisNames()[am.analysisIndex(ID = IDs, throw = TRUE, noMatch = NA)];
}


am.getAnalysis = function(analysisID)
{
  if (!am.haveLock())
  {
    am.lockAndRead();
    on.exit(am.unlock());
  }
  .analysisManagerData$analyses[[ am.analysisIndex(ID = analysisID, throw = TRUE)]];
}

#=====================================================================================================
#
# Adding objects to analyses
#
#=====================================================================================================

# This function can operate on given 'analysis' or can retrieve the analysis first from the manager if analysisID is given
# instead. 

am.addExistingObjectToAnalysis = function(objectID, 
                       analysisID = NULL,
                       analysis = am.getAnalysis(analysisID),
                       anaType = character(0),
                       organism = character(0), 
                       tissue = character(0),
                       dataSource = character(0),
                       allowNonstandardType = FALSE,
                       updateAnalysisManager = !is.null(analysisID))
{
  if (!inherits(analysis, "AM.Analysis"))
    stop("'analysis' must be of class AM.Analysis.");

  if (! objectID %in% analysis$objectIDs)
    analysis$objectIDs = c(analysis$objectIDs, objectID);

  if (is.na(om.objectIndex(ID = objectID)))
    stop("Given 'objectID' does not identify an object in the object manager.");

  anaType = am.standardizeAnalysisType(anaType, throw = !allowNonstandardType)
  analysis$analysisType = unique(c(analysis$analysisType, anaType));
  analysis$organism = unique(c(analysis$organism, organism));
  analysis$tissue = unique(c(analysis$tissue, tissue));
  analysis$dataSource = unique(c(analysis$dataSource, dataSource));
  if (updateAnalysisManager) am.insertAnalysis(analysis);
  analysis;
}

am.addExistingObjectToAnalysisManager = function(objectID, analysisID, ...) 
# ...same rest of arguments as am.addObjectToAnalysis
{
  if (!am.haveLock())
  {
    am.lockAndRead();
    on.exit(am.unlock());
  }

  out = am.insertAnalysis(
    am.addExistingObjectToAnalysis(objectID, am.getAnalysis(analysisID), ...)); 
  am.write();
  out;
}


am.addNewObjectToAnalysis = function(
  # object information
   name,  # This could be a single string or a vector.
   fileName,
   fileType = typeFromExtension(fileName),
   fileReadArgs = defaultReadArguments(fileType),
   fileObject = character(0),
   objectSubsetting = numeric(0),
   metaData = list(),
   fileIsRelative = pathIsRelative(fileName),
  # analysis information
   analysisID = NULL,
   analysis = am.getAnalysis(analysisID),
   anaType,
   organism = character(0),
   tissue = character(0),
   dataSource = character(0),
   allowNonstandardAnalysisType = FALSE,
   allowNonstandardObjectType = FALSE,

   returnBoth = TRUE,
   updateAnalysisManager = !is.null(analysisID),
   verbose = 1, indent = 0)
{
  if (!inherits(analysis, "AM.Analysis"))
    stop("'analysis' must be of class AM.Analysis.");
  anaType = am.standardizeAnalysisType(anaType, throw = !allowNonstandardAnalysisType, 
                                            allowEmpty = FALSE);
  name = am.standardizeObjectType(name, throw = !allowNonstandardObjectType,
                                  allowEmpty = FALSE);
  nameVector = c(analysis$name, anaType, name)

  objectID = om.insertNewObject(nameVector = nameVector,
      fileName = fileName,
      fileType = fileType,
      fileReadArgs = fileReadArgs,
      fileObject = fileObject,
      objectSubsetting = objectSubsetting,
      metaData = metaData,
      fileIsRelative = fileIsRelative, verbose = verbose, indent = indent)
  analysis = am.addExistingObjectToAnalysis(
      objectID,
      analysis = analysis,
      anaType = anaType,
      organism= organism,
      tissue = tissue,
      dataSource = dataSource,
      allowNonstandardType = allowNonstandardAnalysisType);
  if (updateAnalysisManager)
     am.insertAnalysis(analysis, verbose = verbose, indent = indent);

  if (returnBoth) list(objectID = objectID, analysis = analysis) else analysis;
}

# Convenience function for multiple objects within a single file and with single names...
# Assume "names" are a vector with one name per new object to be added.
# fileObjects must be of the same length as names

.subsetOrNothing = function(x, i) 
  if (length(x)==0) x else x[[i]];

am.addMultipleObjectsToAnalysis = function(
   names, 
   fileName,
   fileObjects = character(0),
   fileType = typeFromExtension(fileName),
   fileReadArgs = defaultReadArguments(fileType),
   objectSubsetting = numeric(0),
   metaData = list(),
   fileIsRelative = pathIsRelative(fileName),
  # analysis information
   analysisID = NULL,
   analysis = am.getAnalysis(analysisID),
   anaType,
   organism = character(0),
   tissue = character(0),
   dataSource = character(0),
   allowNonstandardAnalysisType = FALSE,
   allowNonstandardObjectType = FALSE,
   updateAnalysisManager = !is.null(analysisID),
   verbose = 1, indent = 0)
{
  n = length(names)
  if (n>1 && fileType!="RData")
  {
    stop("Multiple objects can only be stored in an RData file.");
  }
  if (fileType=="RData" && n!=length(fileObjects)) stop("length of 'names' and 'fileObjects' must be the same.");

  if (n==0) return(analysis);

  anaType = checkOrExtend(anaType, n);
  organism = checkOrExtend(organism, n);
  tissue = checkOrExtend(tissue, n);
  dataSource = checkOrExtend(dataSource, n);

  for (i in 1:n)
    analysis = am.addNewObjectToAnalysis(
      name = names[i],
      fileName = fileName,
      fileObject = subsetOrNothing(fileObjects, i),
      fileReadArgs = fileReadArgs,
      objectSubsetting = subsetOrNothing(objectSubsetting, i),
      metaData = subsetOrNothing(metaData, i),
      fileIsRelative = fileIsRelative,
      analysis = analysis,
      anaType = subsetOrNothing(anaType, i),
      organism = subsetOrNothing(organism, i),
      tissue = subsetOrNothing(tissue, i),
      dataSource = subsetOrNothing(dataSource, i),
      allowNonstandardAnalysisType = allowNonstandardAnalysisType,
      allowNonstandardObjectType = allowNonstandardObjectType, verbose = verbose, indent = indent,
      updateAnalysisManager = FALSE)$analysis;

  if (updateAnalysisManager)
    am.insertAnalysis(analysis, verbose = verbose, indent = indent);

  analysis;
}

am.addNewObjectToAnalysisManager = function(analysisID, ...)
# ...same rest of arguments as am.addNewObjectToAnalysis
{
  if (!am.haveLock())
  {
    am.lockAndRead();
    on.exit(am.unlock());
    doWrite = TRUE;
  } else
    doWrite = FALSE;

  out = am.addNewObjectToAnalysis(analysisID = analysisID, ..., updateAnalysisManager = TRUE)$analysis;

  if (doWrite) am.write()
  out;
}

#=====================================================================================================
#
# Retrieving objects from analyses
#
#=====================================================================================================

am.allObjectIDsInAnalysis = function(analysisID, analysisType = NULL)
{
  analysis = am.getAnalysis(analysisID);
  allIDs = analysis$objectIDs;
  if (!is.null(analysisType))
  {
    analysisType = am.standardizeAnalysisType(analysisType, throw = FALSE);
    nameVecs = om.existingNameVectors(allIDs);
    keep = sapply(nameVecs, `[`, 2) %in% analysisType;
    allIDs = allIDs[keep];
  } 
  allIDs;
}

  
am.allObjectNamesInAnalysis = function(analysisID, analysisType = NULL)
{
  om.existingNameVectors(am.allObjectIDsInAnalysis(analysisID, analysisType));
}
 

