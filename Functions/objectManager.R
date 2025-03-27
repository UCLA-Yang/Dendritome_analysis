# General object manager.

#===================================================================================================
#
# Basic setup
#
#===================================================================================================

.om.objectIDPrefix = "om.obj.";
.om.crossReferenceIDPrefix = "om.cr."
.om.tempFile = "tmp.om.RData";

om.getBaseDir = function()
{ 
 file.path(sub("/Work/.*", "", getwd()), "Work");
}

om.startupOptions = function(
    useLocalCopy = FALSE,
    copyDir = "ObjectManagerLocalData",
    stopOnFileCollision = TRUE)
{
  list(useLocalCopy = useLocalCopy,
       copyDir = copyDir, 
       stopOnFileCollision = stopOnFileCollision);
}

om.useLocalCopy = function() { om.options("useLocalCopy") }
om.copyDir = function() { om.options("copyDir") }


om.initialize = function(
  storageDir = file.path(om.getBaseDir(), "ObjectManager/Storage"),
  storageFile = "objectManager.RData",
  backupDir = file.path(storageDir, "Backups"),
  maxBackups = 1000,
  backupFileBase = "objectManagerBackup-",
  backupFileExtension = ".RData",
  mustExist = FALSE,
  options = om.startupOptions())
{
  storageFile.full = file.path(storageDir, storageFile);
  if (file.exists(storageFile.full))
  {
    omd = loadAsList(storageFile.full)$.objectManagerData;
    if (is.null(omd)) 
       stop("File ", storageFile.full, " does not contain expected object '.objectManagerData'.");
    if (is.null(omd$objects) || is.null(omd$crossReferences))
       stop("Object '.objectManagerData' in file ", storageFile.full, " does not contain the expected components.");
    omd$storageDir = storageDir;
    omd$storageFile = storageFile;
    omd$backupDir = backupDir;
    omd$backupFileBase = backupFileBase;
    omd$backupFileExtension = backupFileExtension;
    if (!is.null(maxBackups)) omd$maxBackups = maxBackups;
    omd$haveLock = FALSE;
    finfo = file.info(storageFile.full, extra_cols = FALSE);
    omd$fileInfo = finfo;
    omd$options = options;
    .objectManagerData <<- omd;
  } else {
    if (mustExist) stop("File ", storageFile, " does not exist");
    suppressWarnings(dir.create(storageDir, recursive = TRUE));
    suppressWarnings(dir.create(backupDir, recursive = TRUE));
    .objectManagerData <<- list(
         objects = list(),
         crossReferences = list(),
         storageDir = storageDir,
         backupDir = backupDir,
         storageFile = storageFile,
         backupFileBase = backupFileBase,
         backupFileExtension = backupFileExtension,
         maxBackups = if (!is.null(maxBackups)) maxBackups else 1000,
         lockFile = "00-objectManagerLock.RData",
         haveLock = FALSE,
         options = options);
    om.saveManagerData(FORCE = TRUE);
  }
}

om.storageDir = function() .objectManagerData$storageDir;
om.backupDir = function() .objectManagerData$backupDir;
om.maxBackups = function() .objectManagerData$maxBackups;
om.backupFileBase = function() .objectManagerData$backupFileBase;
om.backupFileExtension = function() .objectManagerData$backupFileExtension;
om.storageFile = function() .objectManagerData$storageFile;
om.lockFile = function() .objectManagerData$lockFile;
om.haveLock = function() .objectManagerData$haveLock;

om.options = function(...)
{
  args = list(...);
  if (length(args) > 1) stop("om.options can take at most one argument");
  options = .objectManagerData$options;
  if (length(args)==0) return(options);
  if (!is.null(names(args)))
  {
    comp = names(args)[1];
    options[[comp]] = args[[comp]];
    .objectManagerData$options <<- options;
  } else {
    comp = args[[1]];
    if (!comp %in% names(options)) return(NULL) else return(options[[comp]]);
  }
}
 
#=================================================================================================
#
# Disk operations
#
#=================================================================================================

om.reread = function(haveLock = om.haveLock())
{
  f = file.path(om.storageDir(), om.storageFile());
  if (!file.exists(f)) return(invisible(.objectManagerData));
  finfo = file.info(f);
  finfo.old = .objectManagerData$fileInfo;
  if ( (is.null(finfo.old)) || finfo$size[1]!=finfo.old$size[1] || finfo$mtime!=finfo.old$mtime)
  {
    #printFlush("Rereading...");
    omd = loadAsList(file.path(om.storageDir(), om.storageFile()))$.objectManagerData;
    omd$haveLock = haveLock;
    omd$fileInfo = finfo;
    .objectManagerData <<- omd;
  }
  invisible(.objectManagerData);
}

om.lock = function(timeout = 180, sleepTime = 0.25)
{
  if (!om.haveLock())
  {
    maxAttempts = timeout/sleepTime
    att = 1;
    pid = Sys.getpid()
    lockFile = file.path(om.storageDir(), om.lockFile());
    while (att<=maxAttempts && file.exists(lockFile)) {att = att + 1; Sys.sleep(sleepTime)};
    if (file.exists(lockFile))
      stop("Could not get a lock in ", om.storageDir(), ". Please check the file ", om.lockFile(), ".");
    on.exit(unlink(lockFile));
    lock = file(lockFile, open = "wt");
    writeLines(as.character(pid), con = lock);
    close(lock);
    on.exit(NULL);
    .objectManagerData$haveLock<<-TRUE
  }
}

om.lockAndRead = function(timeout = 180, sleepTime = 0.25)
{
  om.lock(timeout, sleepTime);
  om.reread();
}

om.checkLock = function()
{
  if (!om.haveLock()) return(FALSE);
  lockFile = file.path(om.storageDir(), om.lockFile());
  if (!file.exists(lockFile))
    stop("Could not find a lock file, aborting.");

  lock = file(lockFile, open = "rt");
  pid = readLines(con = lock, n=1, ok = FALSE, warn = FALSE);
  close(lock);

  if (pid!=as.character(Sys.getpid()))
    stop("Lock file exists but belongs to a different process, aborting.");
  TRUE;
}

om.unlock = function()
{
  if (om.checkLock())
  {
    unlink(file.path(om.storageDir(), om.lockFile()));
    .objectManagerData$haveLock <<- FALSE;
  }
}

om.backup = function(minTimeDifference = 12*3600)
{
  if (om.checkLock())
  {
    if (!file.exists(file.path(om.storageDir(), om.storageFile()))) return(NULL);
    suppressWarnings(dir.create(om.backupDir(), recursive = TRUE));
    bf = list.files(om.backupDir(), 
               pattern = spaste("^", om.backupFileBase(), "[0-9]*", om.backupFileExtension()));
    if (length(bf)> 0)
    {
      mtimes = file.mtime(file.path(om.backupDir(), bf));
      lastBackupTime = max(mtimes);
      currentTime = Sys.time();
      if (currentTime - minTimeDifference < lastBackupTime ) return(NULL);
      numbers = replaceMissing(as.numeric(multiSub(
                     c(om.backupFileBase(), om.backupFileExtension()), c("", ""), bf, fixed = TRUE)));
      newNumber = max(numbers) + 1;
    } else
      newNumber = 1;

    newFile = spaste(om.backupFileBase(), prependZeros(newNumber, 8), om.backupFileExtension());
    file.copy(from = file.path(om.storageDir(), om.storageFile()),
              to = file.path(om.backupDir(), newFile), copy.mode = TRUE, copy.date = FALSE);
    if (length(bf)>=om.maxBackups()) 
    {
      firstInd = which.min(mtimes)
      file.remove(file.path(om.backupDir(), bf[firstInd]));
    }
  } else
    stop("Cannot backup: do not have the lock.");
}

om.writeAndUnlock = function(minBackupTimeDifference = 12*3600)
{
  om.write(minBackupTimeDifference = minBackupTimeDifference);
  om.unlock();
}

om.write = function(minBackupTimeDifference = 12*3600)
{
  if (om.checkLock())
  {
    f = file.path(om.storageDir(), om.storageFile());
    save(.objectManagerData, file = file.path(om.storageDir(), .om.tempFile));
    if (filesDiffer(file.path(om.storageDir(), .om.tempFile), f))
    {
      om.backup(minTimeDifference = minBackupTimeDifference);
      file.rename(file.path(om.storageDir(), .om.tempFile), f);
    } else
      file.remove(file.path(om.storageDir(), .om.tempFile));
    .objectManagerData$fileInfo <<- file.info(f);
    invisible(.objectManagerData);
  } else
    stop("Cannot write: do not have the lock.");
}

om.fileRename = function(from, to, pattern = NULL, objectIDs = om.existingObjectIDs(), fixed = TRUE,
                         checkValidity = TRUE)
{
  if (!om.haveLock())
  {
    on.exit(om.unlock());
    om.lockAndRead();
  }
  index = match(objectIDs, om.existingObjectIDs());
  if (any(is.na(index)))
    stop("Some given 'objectIDs' are invalid.");
  objects = .objectManagerData$objects;
  objects[index] = lapply(objects[index], function(obj)
  {
    if (length(pattern) ==0 || grepl(pattern, obj$fileName, fixed = fixed)) 
    {
       newFileName = gsub(from, to, obj$fileName, fixed = fixed);
       if (checkValidity)
       {
         if (file.exists(newFileName)) {
              obj$fileName = newFileName;
         } else 
            warning("om.fileRename: new file does not exist: \n", formatLabels(newFileName, maxCharPerLine = 80, split = "/"));
       } else  obj$fileName = newFileName;
    }
    obj;
  });
  .objectManagerData$objects <<- objects;
  om.writeAndUnlock()
  invisible(.objectManagerData$objects);
}

om.checkFileExistence = function()
{
  fileNames = sapply(.objectManagerData$objects, getElement, "fileName");
  IDs = sapply(.objectManagerData$objects, getElement, "ID");
  dirs = dirname(fileNames);

  dirs2 = unique(dirs);
  out = list();
  checkDirs = file.exists(dirs2)
  if (!all(checkDirs))
  {
    out$missingDirs = dirs2[!checkDirs];
    out$objectIDs.missingDirs = IDs[ dirs %in% out$missingDirs];
    printFlush("Warning: ", sum(!checkDirs), " directories are missing.");
  }

  checkIndex = which(!dirs %in% dirs2[!checkDirs]);
  checkFiles = fileNames[ checkIndex];
  fe = file.exists(checkFiles);
  if (any(!fe)) 
  {
    out$missingFiles = checkFiles[!fe];
    out$objectIDs.missingFiles = IDs[checkIndex][!fe];
    printFlush("Warning: ", sum(!fe), " files are missing (in existing directories).");
  }
  invisible(out);
}

om.renameObjects = function(oldName, newName, matchComponents = NULL, dryRun = FALSE)
{
  if (!om.haveLock())
  {
    on.exit(om.unlock());
    om.lockAndRead();
  }
  objects = .objectManagerData$objects;
  count = 0;
  for (obj in 1:length(objects))
  {
    if (is.null(matchComponents)) matchComponents1 = 1:length(objects[[obj]]$nameVector) else 
      matchComponents1 = intersect(1:length(objects[[obj]]$nameVector), matchComponents);
    fix = match(oldName, objects[[obj]]$nameVector);
    if (!is.na(fix)) {
       printFlush(formatLabels(spaste("Renaming object with name vector ", paste(objects[[obj]]$nameVector, collapse = " | ")),
                   maxCharPerLine = 100));
       objects[[obj]]$nameVector[fix] = newName;
       count = count + 1;
    }
  }
  printFlush(spaste("Renamed ", count, " objects.", if (dryRun) " (DRY RUN)" else ""));
  if (!dryRun) .objectManagerData$objects <<- objects;
  om.writeAndUnlock()
  invisible(.objectManagerData$objects);
}


# Check for file collision for a given object, i.e., whether any other object points to the same file.
# Returns character(0) if there are no collisions or the IDs with which the file collides
om.fileCollision = function(object, fileName = object$fileName, throw = om.options("stopOnFileCollision"))
{
  f.n = normalizePath(fileName, mustWork = TRUE);
  id = object$ID;
  ## Try not locking here: this is a read-only function.
  #if (!om.haveLock())
  #{
  #  on.exit(om.unlock());
  #  om.lockAndRead();
  #}
  om.reread();
  objects = .objectManagerData$objects;
  files = sapply(objects, getElement, "fileName");
  files.norm = normalizePath(files);
  allIDs = sapply(objects, getElement, "ID");

  collisionIDs = setdiff(allIDs[ files.norm == f.n ], id);
  if (length(collisionIDs)==0) return(NULL);

  collisionIndex = match(collisionIDs, allIDs);
   
  if (length(id)==0)
  {
    # must check name vector equality. This is because one can replace an existing object by creating a new one (which is
    # created without ID) and the ID will only be assigned upon insering into OM.
    collisionNameVectors = lapply(objects[ collisionIndex ], getElement, "nameVector");
    match1 = .om.matchNameVector(object$nameVector, collisionNameVectors)
    if (!is.na(match1)) {
      collisionIDs = collisionIDs[-match1];
      collisionIndex = collisionIndex[-match1];
    }
  }

  if (length(collisionIDs)==0) return(NULL)

  # For RData files, the check is more complicated since several objects can be saved in a single file.
  # We want to drop the object that has the same file type, file object and object subsetting.
  if (object$fileType=="RData")
  {
    collisionObjects = objects[collisionIndex];
    collisionFileTypes = sapply(collisionObjects, getElement, "fileType");
    collisionFileObjects = lapply(collisionObjects, getElement, "fileObject");
    collisionFObjSubset = lapply(collisionObjects, getElement, "objectSubsetting");
    if (!all(collisionFileTypes=="RData"))
    {
      printFlush("Apparent file type error: same file name is marked as type 'RData' in one object but not another.\n", 
                 "Dropping into browser.");
      browser();
    }
    # Check file object. 
    # if object$fileObject is empty, all objects automatically collide since empty file objects mean a collision 
    # and non-empty mean the files represent different sets of objects and hence collide 
    if (length(object$fileObject)>0)
    {
      keep.fo = sapply(collisionFileObjects, function(fo) isTRUE(all.equal(object$fileObject, fo)));
    } else keep.fo = rep(TRUE, length(collisionFileObjects));
    # Check file object subsetting. Same rules as above apply.
    if (length(object$objectSubsetting)>0)
    {
      keep.os = sapply(collisionFObjSubset, function(os) isTRUE(all.equal(object$objectSubsetting, os)));
    } else keep.os = rep(TRUE, length(collisionObjects));
    keep = keep.fo & keep.os;
    collisionIDs = collisionIDs[keep];
    collisionIndex = collisionIndex[keep];
  }

  #out = list(ID = character(0),
  #           nameVector = as.character(nameVector),
  #           fileName = fileName,
  #           fileType = fileType,
  #           fileReadArgs = fileReadArgs,
  #           fileObject = fileObject,
  #           objectSubsetting = objectSubsetting,
  #           metaData = metaData);

  if (length(collisionIDs) ==0) return(NULL)

  if (throw)
  {
    .om.reportCollisions(object, collisionIDs);
    browser("file name collision.")
    stop("file name collision, see details above.");
  }
  return(collisionIDs);
}

om.fileCollision.fromName = function(nameVector, fileName, throw = om.options("stopOnFileCollision"))
{
  ind = om.objectIndex(nameVector = nameVector);  ## This includes a re-read
  if (is.na(ind)) return(NULL);
  om.fileCollision(.objectManagerData$objects[[ind]], fileName = fileName, throw = throw);
}

om.clearAll = function()
{
  printFlush(spaste(
      "Warning: this function de-synchronizes object manager from saved data.\n",
      "Data on disk are not removed."));
  #om.lockAndRead();
  omd = .objectManagerData;
  omd$objects = list();
  omd$crossReferences = list();
  .objectManagerData <<-omd;
  omd;
  #om.writeAndUnlock();
}

om.saveManagerData = function(
      storageDir = .objectManagerData$storageDir,
      storageFile = .objectManagerData$storageFile, ...)
{
  args = list(...);
  if ("FORCE" %in% names(args))
  {
    suppressWarnings(dir.create(storageDir, recursive = TRUE));
    save(.objectManagerData, file = file.path(storageDir, storageFile));
  } else
    stop("This function is for emergencies only.");
}

# Existing name vectors and their corresponding IDs.

om.existingObjectIDs = function(reread = !om.haveLock())
{
  if (reread) om.reread();
  sapply(.objectManagerData$objects, getElement, "ID");
}

om.existingCrossRefIDs = function(reread = !om.haveLock())
{
  if (reread) om.reread();
  sapply(.objectManagerData$crossReferences, getElement, "ID");
}

om.existingNameVectors = function(objectIDs = NULL, reread = !om.haveLock())
{
  if (reread) om.reread();
  if (is.null(objectIDs)) {
    lapply(.objectManagerData$objects, getElement, "nameVector");
  } else 
    lapply(.objectManagerData$objects[om.existingObjectIDs() %in% objectIDs],
             getElement, "nameVector");
}

om.nObjects = function(reread = !om.haveLock())
{
  if (reread) om.reread();
  length(.objectManagerData$objects);
}

om.getObject = function(index = NULL, ID = NULL, nameVector = NULL,
                        reread = !om.haveLock(),
                        useLocalCopy = om.options("useLocalCopy"),
                        copyDir = om.options("copyDir"))
{
  if (reread) om.reread();
  if (length(index)==0) 
  {
    index = om.objectIndex(nameVector = nameVector, ID = ID, throw = !useLocalCopy && is.null(ID))
    if (is.na(index))
    {
       file = file.path(copyDir, om.fileNameFromNameVector(nameVector));
       if (!file.exists(file)) stop("Local object copy file ", file, " does not exist.");
       return(readRDS(file))
    }
  }
  if (length(index) > 1)
    stop("om.getObject returns a single object.");
  .objectManagerData$objects[[index]];
}

om.getNewObjectID = function()
{
  if (!om.haveLock())
    stop("Must lock the disk image of object manager before running this function.");
  ex.n = as.numeric(sub(.om.objectIDPrefix, "", om.existingObjectIDs(), fixed= TRUE));
  out = suppressWarnings(max(ex.n, na.rm = TRUE) + 1);
  if (!is.finite(out)) out = 1;
  spaste(.om.objectIDPrefix, out);
}

om.getNewCrossReferenceID = function()
{
  if (!om.haveLock())
    stop("Must lock the disk image of object manager before running this function.");
  ex.n = as.numeric(sub(.om.crossReferenceIDPrefix, "", om.existingCrossRefIDs(), fixed= TRUE));
  out = suppressWarnings(max(ex.n, na.rm = TRUE) + 1);
  if (!is.finite(out)) out = 1;
  spaste(.om.crossReferenceIDPrefix, out);
}

#=====================================================================================================
#
# definitions of classes via new... functions
#
#=====================================================================================================

om.newCrossReferenceIndex = function(
   type = c("names", "dimnames", "component", "value", "attributes"),
   reference = character(0))
{
  type = match.arg(type);
  if (type=="dimnames")
  {
    if (!is.numeric(reference)) 
      stop("When 'type' is 'dimnames', 'reference' must be a single positive integer.");
    if (length(reference)!=1)
      stop("When 'type' is 'dimnames', 'reference' must be a single positive integer.");
    if (reference < 1 || as.integer(reference)!=reference)
      stop("When 'type' is 'dimnames', 'reference' must be a single positive integer.");
  }
  if (type=="value")
  {
    if (!is.numeric(reference))
      stop("When 'type' is 'value', 'reference' must be a vector of integers with exactly one zero.");
    if (length(reference)<2)
      stop("When 'type' is 'value', 'reference' must be a vector of integers with exactly one zero.");
    if (any(reference < 0 || as.integer(reference)!=reference))
      stop("When 'type' is 'value', 'reference' must be a vector of integers with exactly one zero.");
    if (sum(reference==0)!=1)
      stop("When 'type' is 'value', 'reference' must be a vector of integers with exactly one zero.");
  }
  if (type=="attributes" || type=="component")
  {   
    if (!is.character(reference))
      stop("When 'type' is 'attributes', 'reference' must be a single character string.");
    if (length(reference)!=1)
      stop("When 'type' is 'attributes', 'reference' must be a single character string.");
  }
  out = list(type = type, reference = reference)
  class(out) = c("OM.CrossReferenceIndex", class(out));
  out;
}
         

#=====================================================================================================
#
# CrossReference
#
#=====================================================================================================

om.newCrossReference = function(
  objectID1, objectID2,
  matchingReference1, matchingReference2, 
  nameToward1 = "", nameToward2 = "",
  metaData = list())
{
  if (!inherits(matchingReference1, "OM.CrossReferenceIndex"))
    stop("'matchingReference1' must be of type OM.CrossReferenceIndex.");
  if (!inherits(matchingReference2, "OM.CrossReferenceIndex"))
    stop("'matchingReference2' must be of type OM.CrossReferenceIndex.");
  out = list(ID = character(0), 
             objectID1 = objectID1, objectID2 = objectID2,
             matchingReference1 = matchingReference1,
             matchingReference2 = matchingReference2,
             nameToward1 = nameToward1, nameToward2 = nameToward2,
             metaData = metaData);
  class(out) = c("OM.CrossReference", class(out));
  out;
}

om.insertCrossReference = function(crossRef)
{
  if (!om.haveLock())
  {
    om.lockAndRead();
    on.exit(om.unlock());
  }
  if (!inherits(crossRef, "OM.CrossReference")) 
    stop("'crossRef' must be of class 'OM.CrossReference'.");
  omd = .objectManagerData;
  crid = om.crossReferenceID(crossRef$objectID1, crossRef$objectID2, noMatch = NA, throw = FALSE)
  if (length(crossRef$ID)==0) 
  {
    if (is.na(crid)) 
    {
       crossRef$ID = om.getNewCrossReferenceID();
       omd$crossReferences = c(omd$crossReferences, list(crossRef));
       .objectManagerData <<- omd;
       return(crossRef$ID);
    } else 
       crossRef$ID = crid;
  } else
    if (!replaceMissing(crid==crossRef$ID))
      stop("'crossRef' ID does not match that of a cross reference connecting the same objects in data.");

  pos = match(crossRef$ID, om.existingCrossRefIDs());
  if (is.na(pos)) stop("'crossRef' contains invalid ID. Remove the invalid ID and re-run.");
  omd$crossReferences[[pos]] = crossRef;
  .objectManagerData <<- omd;
  om.write();
  crossRef$ID;
}

om.addNewCrossReference = function(...) om.insertCrossReference(om.newCrossReference(...));

#=====================================================================================================
#
# Object
#
#=====================================================================================================


om.newObject = function(nameVector, fileName, 
                        fileType = typeFromExtension(fileName),
                        fileReadArgs = defaultReadArguments(fileType),
                        fileObject = character(0),
                        objectSubsetting = numeric(0),
                        metaData = list(),
                        fileIsRelative = pathIsRelative(fileName))
{
  if (length(nameVector)==0) stop("'nameVector' cannot be empty.");
  if (length(fileName)!=1) stop("'fileName' must be a single character string.");
  if (!file.exists(fileName)) stop("The file referenced by 'fileName' must exist.");

  if (fileType=="RData" && length(fileObject)==0)
    stop("When 'fileType' is 'RData', 'fileObject' must contain the name of an R object saved in the file.");

  if (fileIsRelative) fileName = file.path(getwd(), fileName);

  if (length(objectSubsetting) > 0)
  {
    l = sapply(objectSubsetting, length);
    if (any(l!=1)) stop("All elements in 'objectSubsetting' must be of length 1.");
  }

  out = list(ID = character(0),
             nameVector = as.character(nameVector),
             fileName = fileName,
             fileType = fileType,
             fileReadArgs = fileReadArgs,
             fileObject = fileObject,
             objectSubsetting = objectSubsetting,
             metaData = metaData);

  class(out) = c("OM.Object", class(out));
  out;
}

.om.reportCollisions = function(object, IDs, doReport = om.options("stopOnFileCollision"))
{
  if (length(IDs) > 0 && doReport)
  {
    printFlush(
     spaste("ObjectManager: detected file name collision. Details follow. \n",
            "Inserting object ID ", object$ID, ", name vector \n",
            paste(object$nameVector, collapse = "\n"),
            "\nfile name ", object$fileName,
            "\ncollides with the following objects:"));
    for (id in IDs)
    {
       obj = om.getObject(ID = id);
       printFlush("object ID ", obj$ID, ", name vector \n",
            paste(obj$nameVector, collapse = "\n"),
            "\nfile name ", obj$fileName,
            "\n===============================================================================")
    }
  }
}

om.insertObject = function(object, stopOnDuplicate = FALSE, stopOnFileChange = FALSE, 
                             stopOnFileCollision = om.options("stopOnFileCollision"), verbose = 1, indent = 0)
{
  if (!om.haveLock())
  {
    om.lockAndRead();
    on.exit(om.unlock());
  }
  spaces = indentSpaces(indent);
  omd = .objectManagerData;

  if (!inherits(object, "OM.Object"))
    stop("'object' must be of class 'OM.Object'."); 

  collisionIDs = om.fileCollision(object, throw = stopOnFileCollision);
  if (length(collisionIDs) > 0) .om.reportCollisions(object, collisionIDs);

  # try finding the name vector among existing objects
  pos = om.objectIndex(object$nameVector, object$ID);
  if (is.na(pos))
  {
    if (length(object$ID)==1) stop("'object' contains invalid ID. Remove the invalid ID and re-run.");
    if (verbose>0) printFlush(paste(spaces, "Adding new object to object manager."));
    object$ID = om.getNewObjectID();
    omd$objects = c(.objectManagerData$objects, list(object));
  } else {
    if (stopOnDuplicate) 
      stop("Object of the same name or ID already exists.");
    if (length(object$ID) > 0  && !all(object$nameVector==omd$objects[[pos]]$nameVector))
       stop("ID in given 'object' exists in object manager but has different names.");
    oldObject = omd$objects[[pos]];
    if (stopOnFileChange)
    {
      if (file.exists(oldObject$fileName) &&  (oldObject$fileName != object$fileName))
      {
        s1 = oldObject$fileName;
        s2 = object$fileName;
        i = 1; ind = 1;
        while (s1!="" && s2!="" && substring(s1, 1, 1)==substring(s2, 1, 1)) 
        { 
          i = i + 1;
          if (substring(s1, 1, 1)=="/") ind = i;
          s1 = substring(s1, 2); s2 = substring(s2, 2)
        };
        printFlush("om.insertObject: object has new file name while old file still exists.\n",
             "Old file name (dif part): ", s1, "\n",
             "New file name (dif part): ", s2, "\n",
             "Old file name (dif part): ", substring(oldObject$fileName, ind), "\n",
             "New file name (dif part): ", substring(object$fileName, ind), "\n",
             "Full old file name: ", oldObject$fileName, "\n",
             "Fill new file name: ", object$fileName);
        stop("See above.")
      }
    }
    if (verbose>0) printFlush(paste(spaces, "Replacing existing object in object manager."));
    object$ID = omd$objects[[pos]]$ID;
    omd$objects[[pos]] = object;
  }
  .objectManagerData <<- omd
  om.write();
  object$ID;
}

om.insertNewObject = function(..., verbose = 1, indent = 0) 
   om.insertObject(om.newObject(...), verbose = verbose, indent = indent);

om.removeObjectsFromManager = function(ID = NULL, nameVector = NULL, dryRun = TRUE)
{
  if (is.null(ID) && is.null(nameVector))
  {
    warning("No 'ID' or 'nameVector' specified.");
    return(FALSE);
  }

  if (is.atomic(nameVector)) nameVector = list(nameVector);

  if (!om.haveLock())
  {
    om.lockAndRead();
    on.exit(om.unlock());
  }

  pos = om.objectIndex(nameVector = nameVector, ID = ID, throw = FALSE);

  if (any(is.na(pos)))
  {
    warning("Some of the specified ID or nameVector do not match any objects in objectManager.");
  }

  pos.use = pos[is.finite(pos)];
  
  IDs = sapply(.objectManagerData$objects[pos.use], getElement, "ID");
  crIndex = om.objectCrossReferenceIndices(IDs[is.finite(pos)]);
  if (!dryRun)
  {
    .objectManagerData$objects <<- .objectManagerData$objects[-pos.use] 
    if (length(crIndex) > 0) .objectManagerData$crossReferences  <<- .objectManagerData$crossReferences[-crIndex];
    om.write();
  }
  out = is.finite(pos);
  names(out) = if (length(ID)==0) sapply(nameVector, base::paste, collapse = "|") else ID;
  return(out);
}

..changedCompression = 5;

om.stringMatchScore = function(ref, candidate, minPatternLen = 4, matchStart = 5)
{
  if (ref==candidate) return(0);
  if (multiSubr(c("\\.gz$", "\\.bz2$"), ref)==multiSubr(c("\\.gz$", "\\.bz2$"), candidate)) return(..changedCompression);

  # Change in case
  ref1 = tolower(ref);
  candidate1 = tolower(candidate)
  if (tolower(ref)==tolower(candidate)) return(10);

  # Change in "punctuation"
  ref2 = gsub("[^[:alnum:]]", "_", ref1);
  candidate2 = gsub("[^[:alnum:]]", "_", candidate1)

  if (ref2==candidate2) return(20);

  # repeated pattern in ref that only appears once in candidate...

  nch1 = nchar(ref2);
  nch2 = nchar(candidate2);

  if (nch1-nch2 >= minPatternLen && substring(ref2, 1, matchStart)==substring(candidate2, 1, matchStart))
  {
    patternLen = nch1-nch2;
    if (patternLen <= nch2) 
    {
      st = 1;
      en = nch2-patternLen + 1;
      for (i in st:en)
      {
        pat = substring(candidate, i, i+patternLen-1);
        pos1 = gregexpr(pat, ref, fixed = TRUE)[[1]]
        pos2 = gregexpr(pat, candidate, fixed = TRUE)[[1]];
        if (length(pos1)==2 && length(pos2)==1 && all(pos1>matchStart) && all(pos2 > matchStart) && pos1[1]==pos2[1])
        {
          ref3 = spaste(substring(ref, 1, pos1[2]-1),
                        if (pos1[2] + patternLen <= nch1) substring(ref,  pos1[2] + patternLen, nch1) else "");
          if (ref3==candidate) return(30)
          score1 = 30 + om.stringMatchScore(ref3, candidate);
          if (score1 < 1e6) return(score1);
        }
      }
    }
  }
  return(1e6);
}

om.matchStrings = function(ref, candidates, ...)
{
  printFlush("matching", ref);
  scores = sapply(candidates, om.stringMatchScore, ref = ref, ...);
  if (any(scores < 1e6)) return(candidates[which.min(scores)]) else return(NA);
}
  

om.searchForFile = function(objectIDs, ...)
{
  objectIndex = om.objectIndex(ID = objectIDs);

  files = sapply(.objectManagerData$objects[objectIndex], getElement, "fileName");

  dirs = dirname(files);
  bases = basename(files);

  dirs.unique = unique(dirs);
  fileLists = lapply(dirs.unique, list.files);
  candidates = as.vector(mapply(function(file, candidates)
          om.matchStrings(ref = file, candidates = candidates, ...),
          bases, fileLists[match(dirs, dirs.unique)]));
  data.frame(objectID = objectIDs, recordedFile = bases, candidate = candidates, dir = dirs);
}

#========================================================================================================
#
# R Object retrieval
#
#========================================================================================================

om.subsetRObject = function(x, subsetVector)
{
  if (length(subsetVector)==0) return(x);
  om.subsetRObject(x[[ subsetVector[[1]] ]], subsetVector[-1]);
}

.recognizedExtensions = function()
{
  list(RData = c("RData", "rdata", "rda"),
       RDS = c("RDS", "rds"),
       csv = "csv",
       csv.gz = "csv.gz",
       csv.bz2 = "csv.bz2",
       txt = "txt",
       txt.gz = "txt.gz",
       txt.bz2 = "txt.bz2")
}

om.fileNameFromNameVector = function(nameVector)
{
    require(digest)
    s = paste(nameVector, collapse = " | omomomom | "); 
        ## Something that's not likely to appear in actual name vectors and hence be more resistant to potential
        ## collisions
    spaste(digest(s, algo = "sha1"), ".omObject.RDS");
}
  
om.fileNameFromObject = function(object, type = c("data", "metadata"))
{  
  require(digest);
  type = match.arg(type);
  if (type=="data")
  {
    fn = object$fileName;
    exts = unlist(.recognizedExtensions());
    ind = partialMatch(fn, spaste("\\.", exts, "$"));
    hash = digest(fn, algo = "sha1");
    out = spaste(hash, ".", exts[ind]);
  } else if (type=="metadata") {
    out = om.fileNameFromNameVector(object$nameVector);
  }
  out;
}

om.retrieveRObject = function(object, overrideFileReadArgs = FALSE,
                              fileReadArgs = NULL,
                              useLocalCopy = om.options("useLocalCopy"),
                              copyDir = om.options("copyDir"))
{
  fileName = object$fileName;
  if (!file.exists(fileName))
  {
    if (!useLocalCopy)
      stop(spaste("File\n\n  ", fileName, "\n\n  does not exist."));
    fileName = file.path(copyDir, om.fileNameFromObject(object, type = "data"))
    if (!file.exists(fileName)) 
      stop(spaste("File\n\n  ", object$fileName, "\n\n  does not exist (nor does a local copy)."));
  }
  if (overrideFileReadArgs) object$fileReadArgs = fileReadArgs;
  out = switch(object$fileType, 
    RData = {
      lst = loadAsList(fileName);
      if (!object$fileObject %in% names(lst))
        stop("File ", object$fileName, " does not contain object named ", object$fileObject);
      om.subsetRObject(lst[[object$fileObject]], object$objectSubsetting);
    }, 
    RDS = do.call(readRDS, c(list(file = fileName), object$fileReadArgs)),
    csv = do.call(read.csv, c(list(file = fileName), object$fileReadArgs)),
    csv.bz2 = do.call(read.csv, c(list(file = bzfile(fileName)), object$fileReadArgs)),
    csv.gz = do.call(read.csv, c(list(file = gzfile(fileName)), object$fileReadArgs)),
    txt = do.call(read.table, c(list(file = fileName), object$fileReadArgs)),
    txt.gz = do.call(read.table, c(list(file = gzfile(fileName)), object$fileReadArgs)),
    txt.bz2 = do.call(read.table, c(list(file = bzfile(fileName)), object$fileReadArgs)),
    stop("Unsupported 'fileType' ", object$fileType));
  if (file.exists(object$fileName) && useLocalCopy)
  {
     modTime = file.mtime(object$fileName);
     objFile = file.path(copyDir, om.fileNameFromObject(object, type = "metadata"))
     dataFile = file.path(copyDir, om.fileNameFromObject(object, type = "data"));
     if (file.exists(copyDir))
     {
       if (!file_test("-d", copyDir)) 
          stop("Unable to create 'copyDir'", copyDir, ".\n",
               ". A non-directory file of that name already exists.");
     } else dir.create(copyDir, recursive = TRUE);
     ## Always write the object out. The object could have changed (e.g., changed name) without changing the file.
     saveRDS(object, file = objFile)
     ## Copy the file into dataFile only if data file does not exist or is older.
     if (!file.exists(dataFile) || (file.mtime(object$fileName) > file.mtime(dataFile)))
       file.copy(object$fileName, dataFile);
  }
  out;
}

om.retrieveRObject.fromName = function(nameVector = NULL,
                        reread = !om.haveLock(),
                        overrideFileReadArgs = FALSE, fileReadArgs = NULL,
                        useLocalCopy = om.options("useLocalCopy"),
                        copyDir = om.options("copyDir"))
{
  om.retrieveRObject(om.getObject(nameVector = nameVector, reread = reread, useLocalCopy = useLocalCopy, copyDir = copyDir),
                     overrideFileReadArgs = overrideFileReadArgs, fileReadArgs = fileReadArgs);
}

#========================================================================================================
#
# Object index / search
#
#========================================================================================================


.om.matchNameVector = function(needle, haystack, noMatch = NA)
{
  l1 = length(needle)
  if (l1==0) stop("'needle' is empty.");
  lens = sapply(haystack, length);
  candidates = which(lens==l1);
  index = 0;
  while (length(candidates) > 0 && index < l1)
  {
     index = index + 1;
     candidates = candidates[ which( sapply(haystack[candidates], `[`, index)==needle[index])];
  }
  if (index==l1 && length(candidates)>0)
  {
     return(candidates)
  } else {
     return(noMatch);
  }
}

  

# Match a name vector to an ID
om.objectIndex = function(nameVector = NULL, ID = NULL, noMatch = NA, throw = FALSE, uniqueOnly = TRUE,
                          reread = TRUE, lock = FALSE)
{
  if (lock && !om.haveLock()) 
  {
    om.lockAndRead();
    on.exit(om.unlock());
  } else if (reread) om.reread();

  if (length(ID) > 0)
  {
    ids = om.existingObjectIDs(reread = FALSE);
    if (throw && any(!ID %in% ids))
      stop("Some entries in 'ID' are not among existing IDs.");

    return (match(ID, ids, nomatch = noMatch));
  }
  if (length(nameVector) > 0)
  {
    if (!is.atomic(nameVector))
    {
      # Recursively call itself...
      return(sapply(nameVector, om.objectIndex, ID = NULL, noMatch = noMatch, throw = throw, 
                    uniqueOnly = uniqueOnly, reread = FALSE, lock = FALSE));
    }
    candidates = .om.matchNameVector(nameVector, om.existingNameVectors(reread = FALSE), noMatch = noMatch);
    if (length(candidates) > 1 && uniqueOnly) 
         stop("More than one match for name vector ", paste(nameVector, collapse = " | "));
    if (length(candidates)==1 && is.na(candidates) && throw) 
         stop("No match found for name vector ", paste(nameVector, collapse = " | "));
    return(candidates);
  }
  stop("'nameVector' or 'ID' must be specified.");
}

om.objectID = function(nameVector, noMatch = NA, throw = FALSE, uniqueOnly = TRUE)
{
  if (!om.haveLock()) 
  {
    om.lockAndRead();
    on.exit(om.unlock());
  }
  n = if (is.atomic(nameVector)) 1 else length(nameVector);
  ind = om.objectIndex(nameVector = nameVector, noMatch = NA, throw = throw, uniqueOnly = uniqueOnly,
           reread = FALSE, lock = FALSE);
  out = rep(noMatch, n);
  exist = is.finite(ind)
  out[exist] = sapply(.objectManagerData$objects[ind[exist]], getElement, "ID");
  out;
}

om.crossReferenceIndex = function(objectID1, objectID2, noMatch = NA, throw = FALSE)
{
  if (!om.haveLock()) 
  {
    om.lockAndRead();
    on.exit(om.unlock());
  }
  ids1 = sapply(.objectManagerData$crossReferences, getElement, "objectID1");
  ids2 = sapply(.objectManagerData$crossReferences, getElement, "objectID2");

  match = which( objectID1==ids1 & objectID2==ids2 | objectID1==ids2 & objectID2==ids1 );
  if (length(match)==0) 
  {
    if (throw) stop("crossReference of objects ", objectID1, " and ", objectID2, " not found.");
    return(noMatch);
  }
  if (length(match) > 1)
    stop("Internal error: 2 or more crossReferences link objects ", objectID1, " and ", objectID2, ".");
  return(match);
}

om.crossReferenceID = function(objectID1, objectID2, noMatch = NA, throw = FALSE)
{
  if (!om.haveLock()) 
  {
    om.lockAndRead();
    on.exit(om.unlock());
  }
  ind = om.crossReferenceIndex(objectID1, objectID2, noMatch = NA, throw = throw);
  if (is.na(ind)) noMatch else .objectManagerData$crossReferences[[ind]]$ID;
}

om.objectIndex.fromNames = function(nameVector, 
      searchType = c("any", "all"),
      invert = FALSE,
      partial = FALSE,
      exact = TRUE,
      fixed = TRUE,
      ignoreCase = TRUE,
      ignorePunctuation = exact,
      punctuation = "-._= ",
      reread = TRUE, lock = FALSE)
{
  if (lock && !om.haveLock()) 
  {
    om.lockAndRead();
    on.exit(om.unlock());
  } else if (reread) om.reread();

  searchType = match.arg(searchType);
  if (!is.atomic(nameVector))
  {
    # Recursively call itself...
    lst = sapply(nameVector, om.objectIndex.fromNames, 
                      searchType = "any", invert = FALSE, partial = partial, exact = exact,
                      fixed = fixed, ignoreCase = ignoreCase, ignorePunctuation = ignorePunctuation,
                      punctuation = punctuation);
    if (searchType=="any") {
      out = multiUnion(lst)
    } else
      out = multiIntersect(lst)
    if (invert) out = setdiff(1:om.nObjects(), out);
    return(out);
  }

  l1 = length(nameVector);
  if (l1==0) stop("'nameVector' is empty.");

  exNames = om.existingNameVectors();
  if (ignoreCase) 
  {
     nameVector = tolower(nameVector);
     exNames = lapply(exNames, tolower);
  }

  if (ignorePunctuation)
  {
    pattern = spaste("[", punctuation, "]");
    nameVector = gsub(pattern, "", nameVector);
    exNames = lapply(exNames, function(x) gsub(pattern, "", x));
  }

  lens = sapply(exNames, length);
  candidates = if (partial) which(lens >= l1) else which(lens==l1);
  index = 0;
  while (length(candidates) > 0 && index < l1)
  {
     index = index + 1;
     if (exact) {
       match = which( sapply(exNames[candidates], `[`, index)==nameVector[index]);
     } else {
       match = grep(nameVector[index], sapply(exNames[candidates], `[`, index), fixed = fixed);
     }
     candidates = candidates[ match ];
  }
  if (invert) candidates = setdiff(1:om.nObjects(), candidates);
  return(candidates)
}


om.objectIndex.fromTags = function(tags,
      searchType = c("any", "all"),
      searchComponents = c("nameVector", "metaData"),
      searchMetaDataComponents = character(0),
      invert = FALSE,
      exact = TRUE,
      fixed = TRUE,
      ignoreCase = TRUE,
      ignorePunctuation = exact,
      punctuation = "-._= ")
{
  if (!om.haveLock()) 
  {
    om.lockAndRead();
    on.exit(om.unlock());
  }
  searchType = match.arg(searchType);

  if ("nameVector" %in% searchComponents) 
  {
    searchText = om.existingNameVectors();
  } else searchText = listRep(character(0), om.nObjects());

  if ("metaData" %in% searchComponents && length(searchMetaDataComponents)>0)
  {
    searchText = mymapply(match.fun("c"), searchText, lapply(.objectManagerData$objects, function(obj)
                    lapply(searchMetaDataComponents, function(comp) getElement(obj$metaData, comp))));
  }

  if (ignoreCase) 
  {
     tags = tolower(tags);
     searchText = lapply(searchText, tolower);
  }

  if (ignorePunctuation)
  {
    pattern = spaste("[", punctuation, "]");
    tags = gsub(pattern, "", tags);
    searchText = lapply(searchText, function(x) gsub(pattern, "", x));
  }

  tagMatch = lapply(tags, function(t)
  {
    if (exact)
    {
      match = sapply(searchText, function(x) any(x==t))
    } else
      match = sapply(searchText, function(x) length( grep(t, x, fixed = fixed) ) > 0);
    which(match);
  });


  if (searchType=="any") {
    match = multiUnion(tagMatch)
  } else 
    match = multiIntersect(tagMatch);

  if (invert) match = setdiff(1:om.nObjects(), match);
  return(match)
}

om.objectIDs.fromTags = function(
      tags,
      searchType = c("any", "all"),
      searchComponents = c("nameVector", "metaData"),
      searchMetaDataComponents = character(0),
      invert = FALSE,
      exact = TRUE,
      fixed = TRUE,
      ignoreCase = TRUE,
      ignorePunctuation = exact,
      punctuation = "-._= ")
{
  if (!om.haveLock())
  {
    om.lockAndRead();
    on.exit(om.unlock());
  }
  om.existingObjectIDs()[ om.objectIndex.fromTags(
    tags = tags,
    searchType = searchType, 
    searchComponents = searchComponents,
    searchMetaDataComponents = searchMetaDataComponents,
    invert = invert,
    exact = exact,
    fixed = fixed,
    ignoreCase = ignoreCase,
    ignorePunctuation = ignorePunctuation,
    punctuation = punctuation)];
}


om.hint = function(tags, exact = FALSE, fixed = TRUE, ignoreCase = TRUE,
      ignorePunctuation = exact,
      punctuation = "-._= ") 
{
  om.nameVectors.fromTags(tags, searchType = "all", exact = exact, fixed = fixed, ignoreCase = ignoreCase,
      ignorePunctuation = ignorePunctuation, punctuation = punctuation);
}

om.nameVectors.fromTags = function(
      tags,
      searchType = c("any", "all"),
      searchComponents = c("nameVector", "metaData"),
      searchMetaDataComponents = character(0),
      invert = FALSE,
      exact = TRUE,
      fixed = TRUE,
      ignoreCase = TRUE,
      ignorePunctuation = exact,
      punctuation = "-._= ")
{
  if (!om.haveLock())
  {
    om.lockAndRead();
    on.exit(om.unlock());
  }
  om.existingNameVectors()[ om.objectIndex.fromTags(
    tags = tags,
    searchType = searchType,
    searchComponents = searchComponents,
    searchMetaDataComponents = searchMetaDataComponents,
    invert = invert,
    exact = exact,
    fixed = fixed,
    ignoreCase = ignoreCase,
    ignorePunctuation = ignorePunctuation,
    punctuation = punctuation)];
}

om.hint = function(tags, searchType = "all", searchComponents = "nameVector", 
                   exact = FALSE, fixed = FALSE)
{
  om.nameVectors.fromTags(tags, searchType = searchType, searchComponents = searchComponents,
           exact = exact, fixed = fixed);
}






#=====================================================================================================
#
# Search for linked objects
#
#=====================================================================================================  

om.objectCrossReferenceIndices = function(objectIDs)
{
  if (!om.haveLock()) 
  {
    om.lockAndRead();
    on.exit(om.unlock());
  }
  ids1 = sapply(.objectManagerData$crossReferences, getElement, "objectID1");
  ids2 = sapply(.objectManagerData$crossReferences, getElement, "objectID2");

  which(ids1 %in% objectIDs | ids2 %in% objectIDs);
}

om.objectCrossReferenceIDs = function(objectIDs)
{
  if (!om.haveLock()) 
  {
    om.lockAndRead();
    on.exit(om.unlock());
  }
  indices = om.objectCrossReferenceIndices(objectIDs);
  sapply(.objectManagerData$crossReferences[indices], getElement, "ID");
}

om.referencedObjectIDs = function(objectIDs)
{
  if (!om.haveLock()) 
  {
    om.lockAndRead();
    on.exit(om.unlock());
  }
  indices = om.objectCrossReferenceIndices(objectIDs);
  setdiff(unlist(lapply(.objectManagerData$crossReferences[indices], function(x) c(x$objectID1, x$objectID2))),
          objectIDs);
}

om.referencedObjectIndices = function(objectIDs)
{
  om.objectIndex(ID = om.referencedObjectIDs(objectIDs));
}

#=====================================================================================================
#
# Utility functions
#
#=====================================================================================================  
  
pathIsRelative = function(path)
{
  if (path=="") stop("invalid 'path': cannot be an empty string.");
  !(substr(path, 1, 1) %in% c("/", "~"))
}

typeFromExtension = function(name)
{
  recognized = list(RData = c("RData", "rdata", "rda"),
                   RDS = c("RDS", "rds"),
                   csv = "csv",
                   csv.gz = "csv.gz",
                   csv.bz2 = "csv.bz2",
                   txt = "txt",
                   txt.gz = "txt.gz",
                   txt.bz2 = "txt.bz2")
  type = sapply(recognized, function(exts)
  {
    any(sapply(exts, function(ext) length(grep(spaste("\\.", ext,"$"), name))>0 ));
  });
  if (any(type)) return(names(recognized)[type]);
  return(NA);
}

defaultReadArguments = function(type)
{
  switch(type, 
         txt = list(sep = "\t", quote = "\"", header = TRUE),
         txt.gz = list(sep = "\t", quote = "\"", header = TRUE),
         txt.bz2 = list(sep = "\t", quote = "\"", header = TRUE),
         csv = list(quote = "\""),
         csv.gz = list(quote = "\""),
         csv.bz2 = list(quote = "\""),
         list());
}
  
#on.getNewIndex = function(nameVector, sep = ".", nZeros = 4)
#{
#  l1 = length(nameVector);
#  if (l1==0) stop("'nameVector' is empty.");
#
#
#  if (om.nObjects() == 0) 
#  {
#    idVec = rep(1, l1);
#  } else {
#    exNames = om.existingNameVectors();
#    exIDs = om.existingObjectIDs();
#    exIDs.split = lapply(strsplit(exIDs, split = sep), as.numeric);
#
#    lens = sapply(exNames, length);
#    idVec = numeric(l1);
#    candidates = 1:om.nObjects();
#    i = 1;
#    while (length(candidates) > 0 && i <= l1)
#    {
#      candidates.old = candidates;
#      candidates = candidates[ lens[candidates] >=i];
#      candidates = candidates[ sapply(exNames, `[`, i)==nameVector[i] ];
#      i = i+1;
#    }
