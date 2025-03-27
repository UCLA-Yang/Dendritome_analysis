
# Individual analysis for general data, either counts (DESEq2) or transformed data (limma)

# version 001: based on RNAseq analysis version 2
#  . Adding code calling limma
#  . Adding analysis manager calls to store the results.
#
# Version 001-10: clean up in file names. For clean file naming, combineNames should not contain analysis
# name.

# Version 002: Adding weights for DESeq2 analysis, supported by DESeq from version 1.16.

# Version 003: Removing calculation of weights since they are best calculated in preprocessing.

# Version 004: Changing arguments sent to combResultsPPFnc, compatible with GNVFunctions-006, plus adding
# arguments regarding log-transformation of data.

# Version 004-02: default for enrichment.compress changed to TRUE and for interleaveResults to FALSE

# Version 004-03: adding test names/pretty to pretifyList is optional

# Version 005: Reorganizing the enrichment calculation. Calculation is done on a single combined collection;
# the results are then split according to given collections. The enrichment table component of the 
# programmatic output of standardAnalysis will have component names consistent with output of
# enrichmentAnalysis 

# bugfix: some files were saved with .gz suffix but were not compressed.
# Important bugfix: it appears that gene IDs for rescued/exacerbated genes were recorded incorrectly.

# Version 006: Cleaning up and changing the names of output lists that contain gene IDs and gene symbols.

# Version 007: By default, only enrichment in the combined collection is returned; this is done so that the
# analysis does not take as much memory (and disk space upon saving). By default, no columns from the
# enrichment table are dropped to allow the split.
# . designInfo can now contain additional columns combineSignificantHits which signals whether significant
# hits in both directions should be combined into a (third) list of genes, and splitSignificantHitsBy which,
# if not NA, should be a test name by whose sign the significant hits should be split. 

# Version 007-01: 
#   . Fix the default in combineResults
#   . Add nSignif.byCombine, nSignif.byCombine.combinedThresholds, nSignif.all.combinedThresholds to output
#   . Adding ability to run parallel threads for DESeq
#   . Argument amComponent can be used to specify analysisType other than "individual analysis", e.g., "individual
#     analysis with 1 SVA covariate" etc.
#   . Column referenced in splitSignificantHitsBy can also be external and contained in any of the geneAnnotation tables.
#   . Added argument dropRepeated.exceptions that allows one to specify patterns of columns names in combined results that
#     should be protected from removal because of duplication. Default value is "FDR.for". 
#     changePAdjToIt turns out that FDR columns were being removed for being
#     duplicated (all FDRs equal a constant value near 1 that is the same between two columns).
#   . Added ability to include in results columns from the mcols accessor of the DESeq2 object which contain actual fit
#     coefficients, deviance, variances etc.

# 007-02:
#   . Bugfix in a call to geneSetFnc may change the output geneSetFnc where backupGeneIDs are used.
#   . geneSetFnc is also supplied designInfo as one of the arguments.

# 008:
#   . Instead of combining all thresholds into one pile when doing enrichment analysis, split the gene sets by thresholds
#     since it does not really make sense to run Bonferroni correction (not to mention FDR) over different thresholds.

# 008-02:
#   . Optionally and by default, combined individual results are returned and saved in two versions, one containing
#     columns in geneAnnotation.end and one without them. The one without is saved as combined individual results in
#     analysis manager, the extended as "extended combined individual results".
#   . infinite values in limma-voom Z statistics are replaced by appropriately signed 1.4*maximum abs. finite Z.

# 008-03:
#   . New defaults for combineResults etc that should produce more sensible results
#   . Fixed what looks a strange bug in splitEnrichmentResultsBySubcollections: apparently the matrices 
#     containing pvalues and overlap counts do not contain rows for all sets in the collection; this was 
#     assumed by the function. The bug fix makes the function work without making this assumption.
#   . Fixed a bug that gave duplicate gene set IDs in createCollectionOfDEGenes. This looks like a long standing bug that
#     affected most of the collections of DE genes.
#   . prettifyList is ordered in decreasing order of nchar of the from component.

# 009: 
#   . createCollectionOfDEGenes now by default drops empty sets; it accepts argument dropEmptySets that can make
#     it retain the empty sets. A corresponding argument has been added to the main function.
#   . Add checking for zero-variance genes after weighing which sometimes happens.
#   . geneSetFnc now gets the extended results
#   . Adding ability to set size factor calculation arguments for DESeq2 (DESeqSizeFactorArgs)
#   . Adding ability to exclude specific genes (features) from the results and lists of significant genes
#   . Adding argument includedFeatureIndex to geneSetFnc
#   . Changing the way amAnalysis is handled. Although this should be completely transparent to the user, the whole code
#     now carries just the analysis ID, not the actual analysis; upon each addition to the analysis, the analysis manager
#     is updated immediately.
#   . Slight chage in the name of enrichment folders when tag %threshold is used: the underscores are stripped because
#     Windows seems to have a problem with the folder names containing either ._ or _.
#   . Adding arguments exportNSignificant.perThreshold (default FALSE) and exportNSignificant.combinedThresholds (default
#     TRUE)
#   . Switching standardized group loading from a file to internal anRichmentMethods. The file is not distributed with
#     code so code won't work on other people's machines.
#
# 010:
#   . enrichmentAnalysis is switched to enrichmentAnalysis.Entrez
#   . Generalizing the association function

#   Note: for now, the combined result post-processing function combResultsPPFnc is given the full results, including
#   excluded features.


# For designInfo for DESeq2: 
#      . if contrast.name is NA, coefficientName is used;
#      . if contrast.name is "", the numerator and denominator are split using ", " as the separator and used
#           as numerator and denominator.
#      . For limma: simple contrasts (name, numerator, denominator) should work but are not tested yet. Note
#      that limma will not work with a contrast containing a level that is not explicitly included in in the
#      model matrix, necessitating a design without intercept (the intercept is replaced by the last level of
#      a factor). 

# Although DESeq converts : in interaction terms to ., supply the interaction names with :. The code will
# convert them automatically, and the proper names are necessary for removal of dependent covariates to work.

# version 001-01: 
#  . cleaned up interpretation of designInfo
#       if contrast.* is non-missing, it is used, otherwise coefficientName is used.
#  . designInfo can have two additional columns: testType and reducedDesign
#  . New argument testColName for labeling the test column in output tables.
#  . DESeq is run automatically such that same designs are carried out sequentially

# version 002-00:
#  . create a collection of DE genes
#  . record information in analysis manager


pairwiseContrastInfo = function(levelOrder, checkNames = TRUE)
{
  nLevels = length(levelOrder);
  levelMatrix = matrix(levelOrder, nrow = length(levelOrder), ncol = length(levelOrder));
  # manual ordering:
  lev1 = t(levelMatrix)[lower.tri(levelMatrix)];  # fold change denominators
  lev2 = levelMatrix[lower.tri(levelMatrix)];  # fold change numerators
  lev1.n = if (checkNames) make.names(lev1) else lev1;
  lev2.n = if (checkNames) make.names(lev2) else lev2;
  nContrasts = length(lev1)
  contrasts = do.call(rbind, mymapply(c, lev2.n, lev1.n))
  setColnames(contrasts, c("contrast.numerator", "contrast.denominator"));
}

checkModelRegularityForContrasts = function(sampleAnnot, varName, levels, covars, tol = 1e-10)
{
  sampleAnnot = as.data.frame(sampleAnnot);
  keepSamples = sampleAnnot[[varName]] %in% levels;
  dataMat = sampleAnnot[ keepSamples, union(varName, covars) ]
  mm = model.matrix(as.formula(spaste("~", paste(union(varName, covars), collapse = "+"))),
                    data = dataMat);
  means = colMeans(mm);
  means[means==0] = 1;
  mm = mm/matrix(means, nrow(mm), ncol(mm), byrow = TRUE);
  prod = t(mm) %*% mm;
  rcond(prod) > tol;
}
  
removeDependentCovariates = function(modelMat, termsOfInterest)
{
  nPred = ncol(modelMat);
  termIndex = match(termsOfInterest, colnames(modelMat));
  depCols = dependentColumns(modelMat, testOrder = c(termIndex, setdiff(1:nPred, termIndex)),
                             returnColumnNames = TRUE);
  if (any(depCols %in% termsOfInterest))
    stop("Some of the 'termsOfInterest' are among the auto-removed predictors.\n",
         " Please check and correct the design or data.");

  if (any(! (depCols %in% colnames(predictors))))
      stop("Some of the dependent columns are interaction terms. These cannot\n",
           "be auto-removed at this time; please remove manually:\n",
           depCols[! (depCols %in% colnames(predictors))]);
  modelMat[, ! (colnames(modelMat) %in% depCols), drop = FALSE];
}

removeDependentTermsFromDesign = function(design, modelMat, coeffName)
{
  if (is.na(coeffName)) return(design);
  termLabels = attr(terms(formula(design)), "term.labels")
  modelMat = modelMat[, ! colnames(modelMat) %in% "(Intercept)", drop = FALSE]; 
  if (!isTRUE(all.equal(termLabels, colnames(modelMat))))
  {
    warning(immediate. = TRUE, 
            "removeDependentTermsFromDesign: term.labels and colnames of modelMat do not agree, \n",
            "will use design as is.");
    return(design);
  }
  termIndex = match(make.names(coeffName), make.names(colnames(modelMat)));
  if (is.na(termIndex)) 
     stop("Cannot find ", coeffName, " among column names of modelMat:\n",
          paste(colnames(modelMat), collapse = ", "));
      
  
  nPred = ncol(modelMat);
  depCols = dependentColumns(modelMat, testOrder = c(termIndex, setdiff(1:nPred, termIndex)),
                             returnColumnNames = FALSE);
  out = paste("~", paste(colnames(modelMat)[!depCols], collapse = " + "));
  out;
}


runDESeq = function(data, pheno, 
      dataWeights = NULL, designInfo,
      sizeFactorArgs = list(type = "ratio", locfunc = stats::median),
      DESeqArgs,
      goodSamplesGenes, parallel = FALSE, nWorkers = 2,
      dropColumns, addColumns, changePAdjTo, changeCoeffTo, testNameSep,
      resultsArgs, 
      keepFullModels,
      verbose, indent,
      ...)
{
  require(DESeq2);
  extraArgs = list(...)
  if (length(extraArgs) > 0)
    warning(.immediate = TRUE,
       spaste("runDESeq: ignoring extra arguments ", paste(names(extraArgs), collapse = ", ")));

  designs = designInfo$design;
  effectiveDesigns = designs;
  spaces = indentSpaces(indent);
  nDesigns = length(designs);
  deseq = results = vector(mode = "list", length = nDesigns);
  oldDesign = "";
  oldTestType = "";
  designOrder = order(designs);
  testTypes = designInfo$testType;
  if (is.null(testTypes)) testTypes = rep("Wald", nDesigns);
  reducedDesigns = designInfo$reducedDesign;
  if (is.null(reducedDesigns)) reducedDesigns = rep("", nDesigns);
  testNames = designInfo$testName;
  if (parallel) DESeqArgs = c(DESeqArgs, list(parallel = parallel, BPPARAM = MulticoreParam(workers = nWorkers)));
  if (is.null(rownames(pheno))) rownames(pheno) = spaste("Sample.", 1:nrow(pheno));
  if (length(dataWeights) > 0)
    dimnames(dataWeights) = dimnames(data);
  if (!goodSamplesGenes$allOK)
  {
    data = data[goodSamplesGenes$goodSamples, goodSamplesGenes$goodGenes, drop = FALSE];
    if (length(dataWeights) > 0) 
      dataWeights = dataWeights[goodSamplesGenes$goodSamples, goodSamplesGenes$goodGenes, drop = FALSE];
    pheno = pheno[goodSamplesGenes$goodSamples, drop = FALSE];
  }
  contrasts = mymapply(c, designInfo$contrast.name, designInfo$contrast.numerator, 
                        designInfo$contrast.denominator);
  for (dd in 1:nDesigns)
  {
    d = designOrder[dd];
    if (verbose > 1) printFlush(spaste(spaces, " ..working on design ", designs[d], "\n", 
                                spaces, "   and reduced design ", designInfo$reducedDesign[d]));
    mm = model.matrix(as.formula(designs[d]), data = pheno);
      ## Note to self: model.matrix drops rows that have missing data. perhaps options(na.action) could
      ## change the behaviour but for now go with it as is.
    coeffName = designInfo$contrast.name[d];
    if (is.na(coeffName)) coeffName = designInfo$coefficientName[d];
    effectiveDesigns[d] = removeDependentTermsFromDesign(designs[d], mm, coeffName);
    if (verbose > 1) 
        printFlush(spaste(spaces), " ..using effective design ", effectiveDesigns[d]);
    keepSamples = rownames(pheno) %in% rownames(mm);
    print(table(keepSamples));
    if (!all(keepSamples))
    {
      pheno1 = pheno[keepSamples, , drop = FALSE];
      data1 = data[keepSamples, , drop = FALSE];
      dataWeights1 = if (!is.null(dataWeights)) dataWeights[keepSamples, , drop = FALSE] else NULL;
    } else {
      pheno1 = pheno;
      data1 = data;
      dataWeights1 = dataWeights;
    }
    # Check that the data are suitable for analysis... DESeq2 sometimes crashes in situations where only a few (or none)
    # of the samples of a particular gene have both non-zero values and relatively high weight
    dataWeights2 = dataWeights1;
    dataWeights2[dataWeights1 < 0.5] = 0;
    gs = goodSamplesGenes(data1*dataWeights2, weights = dataWeights2, minFraction = 0.25, minNSamples = 1, minNGenes = 1);
    if (!gs$allOK)
    {
       data1 = data1[gs$goodSamples, gs$goodGenes];
       pheno1 = pheno1[gs$goodSamples,  ];
       if (!is.null(dataWeights1)) dataWeights1 = dataWeights1[gs$goodSamples, gs$goodGenes];
    }
    ver = unlist(packageVersion("DESeq2"))
    if (ver[1] == 1 && ver[2] < 28)
    {
      if (!is.null(dataWeights1)) dimnames(dataWeights1) = NULL;  ## This seems to be required for assigning weights to 
                                                          ## assays(object)[["weights"]]
      ####### !!!!!! In new DESeq2 this seems to cause errors (!!)
    }
    if (designs[d]==oldDesign && testTypes[d]=="Wald" && oldTestType=="Wald")
    {
      deseq1 = deseq[[designOrder[dd-1]]];
    } else {
      object = DESeqDataSetFromMatrix.PL(t(data1),
                                sizeFactorArgs = sizeFactorArgs, 
                                colData = pheno1, design = as.formula(effectiveDesigns[d]));
      if (!is.null(dataWeights)) assays(object, withDimnames  = FALSE)[["weights"]] = t(dataWeights1);
      runArgs = c(list(object, test = testTypes[d]), 
                  DESeqArgs);
      if (testTypes[d]=="LRT") 
      {
          runArgs$reduced = as.formula(designInfo$reducedDesign[d]);
          runArgs$betaPrior = FALSE;
      }
      deseq1 = try(do.call(DESeq, runArgs));
      if (inherits(deseq1, "try-error"))
      {
        printFlush("DESeq returned an error, dropping into browser.");
        browser();
        design1 = effectiveDesigns[d];
        data1 = t(data1);
        rownames(pheno1) = colnames(data1)
        weights = t(dataWeights1);
        save(data1, pheno1, design1, weights, file = "~/errorData.RData");
        ds = DESeqDataSetFromMatrix(data1, colData = pheno1, design = as.formula(design1));
        assays(ds, withDimnames  = FALSE)[["weights"]] = weights
      }
      oldDesign = designs[d];
      oldTestType = testTypes[d];
    }

    deseq[[d]] = deseq1;
    # Get results
    if (verbose > 0)
       printFlush(spaste(spaces, " ..getting results ", testNames[d]));
    haveCoeff = TRUE;
    if (is.na(contrasts[[d]] [1])) 
    {
      if (is.na(designInfo$coefficientName[d]) && testTypes[d]=="LRT")
      {
        results[[d]] = try(do.call(DESeq2::results, c(list(object = deseq1), resultsArgs)));
        haveCoeff = FALSE;
      } else {
        results[[d]] = try(do.call(DESeq2::results,
                     c(list(object = deseq1, name = make.names(designInfo$coefficientName[d])),
                       resultsArgs)));
      }
      if (inherits(results[[d]], "try-error")) {
        printFlush("DESeq2::results returned an error (place 1). Dropping into browser.");
        browser();
      }
    } else {
      if (contrasts[[d]] [1]=="")
      {
        contrast.arg = strsplit(contrasts[[d]] [2:3], split = ", ", fixed = TRUE);
      } else 
        contrast.arg = contrasts[[d]];
      results[[d]] = try(do.call(DESeq2::results,
                   c(list(object = deseq1, contrast = contrast.arg),
                     resultsArgs)))
      if (inherits(results[[d]], "try-error")){
        printFlush("DESeq2::results returned an error (place 2). Dropping into browser.");
        browser();
      }

    }
    if (testTypes[d]=="LRT")
    {
      # I am assuming a two-sided test here.
      if (haveCoeff) 
      {
        results[[d]]$stat = sign(results[[d]]$log2FoldChange) * qnorm(results[[d]]$pvalue/2, lower.tail = FALSE);
      } else {
        results[[d]]$stat = qnorm(results[[d]]$pvalue/2, lower.tail = FALSE)
        ### Note: in this case all Z statistics will be positive; keep the half in p/2 so that formulas for calculating p
        ### from Z used in getTopGenes are valid.
      }
    }
    if (length(addColumns)>0)
    {
      mc = as.data.frame(mcols(ds, use.names = TRUE));
      inMC = addColumns %in% names(mc);
      if (!all(inMC)) stop("Not all entries in addColumns are valid column names in 'mcols(ds)'.\n",
         "Valid choices are: ", formatLabels(paste(names(mc), collapse = ", "), maxCharPerLine = 80), "\n",
         "\nInvalid entries are: ", formatLabels(paste(addColumns[!inMC], collapse = ", "), maxCharPerLine = 80));
      results[[d]] = cbind(results[[d]], mc[match(rownames(out), rownames(mc)), addColumns, drop = FALSE]);
    }
    results[[d]] = results[[d]][, ! colnames(results[[d]]) %in% dropColumns];
    results[[d]] = results[[d]][ match(goodSamplesGenes$colNames, rownames(results[[d]])), , drop = FALSE];
    rownames(results[[d]]) = goodSamplesGenes$colNames;
    colnames(results[[d]]) = spaste(multiSub(c("baseMean", "^log2FoldChange", "^stat", "^pvalue", "^padj"), 
                                     c("meanExpr", changeCoeffTo, "Z", "p", changePAdjTo), colnames(results[[d]])), 
                            testNameSep, testNames[d]);
  }
  names(deseq) = testNames;
  names(results) = designInfo$testName;
  out = list(results = results);
  if (keepFullModels) out$fullModels = deseq;
  out;
}

getNormalizedCounts = function(deseqObjects, designs)
{
  oldDesign = "";
  currentCounts = NULL;
  nDesigns = length(designs);
  counts = list();
  for (d in 1:nDesigns)
  {
    if (d>1 && designs[d]==oldDesign)
    {
       counts[[d]] = counts[[d-1]]
    } else {
       counts[[d]] = t(counts(deseqObjects[[d]], normalized = TRUE));
       oldDesign = designs[d];
    }
  }
  names(counts) = names(deseqObjects);
  counts;
}

# Code to run limma functions

runLimmaVoom = function(...)
{
  rumLimma(voom = TRUE, ...);
}

runLimma = function(data, pheno, dataWeights, designInfo, lmFitArgs,
       normalizationMethod = "TMM",
       voom, 
       dataNameInDesign,
       coxphArgs,
       dropColumns, changePAdjTo, changeCoeffTo, testNameSep,
       topTableArgs, 
       keepFullModels, 
       coxph.checkRobustResults = TRUE,  
          ## This option causes the resulting significance to be the lower of the standard and robust significances
       gcInterval = 5000,
       verbose, indent, ...)
{
  require(limma); 
  spaces = indentSpaces(indent);
  extraArgs = list(...)
  if (length(extraArgs) > 0)
    printFlush(spaste(spaces, " runLimma: ignoring extra arguments ", paste(names(extraArgs), collapse = ", ")));
  designs = designInfo$design;
  nDesigns = length(designs);
  lmFits= vector(mode = "list", length = nDesigns);
  oldDesign = "";
  oldTestType = "";
  designOrder = order(designs);
  testNames = designInfo$testName;
  testTypes = designInfo$testType;
  testTypes = multiSub(c("Wald", "LRT"), c("t", "F"), testTypes)
  #if (!is.null(dataWeights)) dataWeights = t(dataWeights);
  rownames(pheno) = rownames(data);
  if (is.null(rownames(pheno))) rownames(pheno) = rownames(data) = spaste("Sample.", 1:nrow(pheno));
  lmFitResults = list();
  for (dd in 1:nDesigns)
  {
    d = designOrder[dd];
    if (verbose > 1) printFlush(spaste(spaces, " ..working on design ", designs[d], "\n",
                                spaces, "   and reduced design ", designInfo$reducedDesign[d]));
    #if (designs[d]==oldDesign)
    #{
    #  lmFits[[d]] = lmFits[[designOrder[dd-1]]];
    #} else {
      if (testTypes[d]=="CoxPH")
      {
        modelData0 = data.frame(lkqejhlkxejhfqoiesjhxlwkhflwkx = 1:nrow(pheno), pheno);
        names(modelData0)[1] = dataNameInDesign;
      } else modelData0 = pheno;
      mm = model.matrix(as.formula(designs[d]), data = modelData0);
        ## Note to self: model.matrix drops rows that have missing data. perhaps options(na.action) could
        ## change the behaviour but for now go with it as is.
      keepSamples = rownames(pheno) %in% rownames(mm);
      #print(table(keepSamples));
      if (!all(keepSamples))
      {
        pheno1 = pheno[keepSamples, , drop = FALSE];
        data1 = data[keepSamples, , drop = FALSE];
        if (!is.null(dataWeights)) 
            dataWeights1 = dataWeights[keepSamples, , drop = FALSE] else dataWeights1 = NULL
        modelData0 = modelData0[keepSamples, , drop = FALSE];
      } else {
        pheno1 = pheno;
        data1 = data;
        dataWeights1 = dataWeights;
      }
      if (!is.null(dataWeights))
         tDataWeights1 = t(dataWeights1) else tDataWeights1 = NULL;
      if (voom)
      {
         dge = DGEList(counts = t(data1));
         dge <- calcNormFactors(dge, method = normalizationMethod);
         v <- do.call(limma::voom, 
                 c(list(counts = dge, design = mm), 
                    if (!is.null(dataWeights)) list(weights = tDataWeights1) else list(),
                   lmFitArgs));
         if (!is.null(dataWeights)) v$weights = v$weights * tDataWeights1;
      }
      if (testTypes[d] %in% c("F", "t"))
      {
         if (voom) 
         {
            lmFits[[d]] = do.call(lmFit,
                 c(list(object = v,
                        design = mm),
                        lmFitArgs))
         } else 
            lmFits[[d]] = do.call(lmFit,
                 c(list(object = t(data1),
                        weights = tDataWeights1,
                        design = mm),
                        lmFitArgs));
         # Get the results right away
         fitCoeffNames = colnames(lmFits[[d]]@.Data[[1]]);
         if (!is.na(designInfo$contrast.name[d]))
         {
           conCoeffNames = spaste(designInfo$contrast.name[d], 
                                     c(designInfo$contrast.numerator[d], designInfo$contrast.denominator[d]));
           conInFit = conCoeffNames %in% fitCoeffNames;
           if (sum(conInFit)==2)
           {
             contrast.matrix = makeContrasts(contrasts = paste(conCoeffNames, collapse = " - "),
                                             levels = fitCoeffNames)
             lmFits[[d]] = contrasts.fit(lmFits[[d]], contrast.matrix)
             name = 1; sign = 1;
           } else if (sum(conInFit)==1) {
             name = conCoeffNames[conInFit];
             sign = conInFit[1] - conInFit[2];
           } else  {
             cat("Something is wrong: neither of the two contrast coefficients in among fit coefficients.\n");
             browser();
           }
         } else {
           sign = 1;
           name = designInfo$coefficientName[d];
         }
         name1 = name;
         if (testTypes[d]=="F") name1 = NULL;
         if (is.na(name)) name1 = NULL; ## topTable will return the F statistic if coef is NULL
         if (testTypes[d]=="F" || is.null(name1)) {
           dropColumns1 = c(setdiff(fitCoeffNames, name), setdiff(dropColumns, "F"))
         } else dropColumns1 = dropColumns;
         eb = eBayes(lmFits[[d]]);
         tt = try(do.call(topTable, c(list(fit = eb, coef = name), topTableArgs)));
         if (inherits(tt, "try-error"))
         {
           printFlush("Error in lmFitResults. Dropping into browser.");
           browser()
         }
         if (testTypes[d] =="F")
         {
           if (!name %in% names(tt))
           {
             printFlush(spaste("coefficientName '", name, "' not in column names of topTable returned by limma. Column names:\n",
                  formatLabels(paste(colnames(tt), collapse = ", "), maxCharPerLine = 80, split = ", ")));
             browser();
           }
           names(tt)[ names(tt)==name ] = "logFC";
         } else {
           tt$t = sign * tt$t;  ## t is not among the columns of top table for F test
         }
         tt$logFC = sign*tt$logFC;
         Z = ZfromP(tt$logFC, tt$P.Value);
         if (any(Z %in% c(-Inf, Inf)))
         {
           maxZ = max(abs(Z[is.finite(Z)]), na.rm = TRUE);
           if (!is.finite(maxZ)) maxZ = 100;
           Z[Z==-Inf] = -1.4 * maxZ;
           Z[Z==Inf] = 1.4 * maxZ;
         }
         out = cbind(tt, Z = Z);
         out = out[, !(colnames(out) %in% dropColumns)];
         if ((!"F" %in% dropColumns) && (!"F" %in% colnames(out)))
         {
           F = classifyTestsF(eb, fstat.only = TRUE)
           out$F = as.numeric(F);
           out$F.df1 = attr(F, "df1");
           out$F.df2 = attr(F, "df2");
         }
         lmFitResults[[d]] = setColnames(out, spaste(multiSub(c("AveExpr", "^logFC", "^stat", "^P.Value", "^adj.P.Val"),
                                          c("meanExpr", changeCoeffTo, "Z", "p", changePAdjTo), colnames(out)),
                                 testNameSep, designInfo$testName[d]));
      } else if (testTypes[d]=="CoxPH")
      {
        if (voom)
        {
          fitData = t(v$E)
          fitWeights = v$weights;
          if (!is.null(dataWeights)) fitWeights = fitWeights * tDataWeights1;
          printFlush(spaste(spaces, " runLimma: Using voom in Cox PH testing."));
        } else  {
          fitData = data1;
          fitWeights = tDataWeights1;
        }
        # If a non-robust model is explictly requested, we have to remove the weights if they contain non-integer
        # weights.
        if ("robust" %in% names(coxphArgs))
        {
          if (!coxphArgs$robust & any(fitWeights!=round(fitWeights))) 
          {
            fitWeights = NULL;
            printFlush(spaste(spaces, " runLimma warning: Non-robust Cox PH model requested but weights contain ",
                                      "non-integer entries. Removing weights."));
          }
        }
        # Fit Cox proportional hazard models.
        nGenes = ncol(data1);
        res = matrix(NA, nGenes, 5)
        colnames(res) = spaste(c("meanExpr", changeCoeffTo, "Z", "p", changePAdjTo), testNameSep, designInfo$testName[d]);
        rownames(res) = colnames(data1);
        if (is.null(fitWeights)) 
        {
           res[, 1] = colMeans(fitData);
        } else {
           tFitWeights = t(fitWeights);
           res[, 1] = colMeans(fitData * tFitWeights)/colMeans(tFitWeights);
        }
        updateInterval = max(round(nGenes/100), 1);
        if (verbose > 0) pind = initProgInd();
        if ("control" %in% names(coxphArgs))
        {
           maxIter = coxphArgs$control$iter.max;
        } else
           maxIter = coxph.control()$iter.max;
        nFail = 0;
        nExceedIter = 0;
        for (g in 1:nGenes)
        {
          modelData = modelData0;
          #if (inherits(try( {modelData[, 1] = fitData[, g]}), "try-error")) browser();
          modelData[, 1] = fitData[, g];
          if (!is.null(fitWeights))
          {
             if (!is.null(dataWeights1))
             {
                keepSamples2 = tFitWeights[, g] > 0;
                modelData = modelData[keepSamples2, , drop = FALSE];
                w1 = tFitWeights[keepSamples2, g]
             } else w1 = tFitWeights[, g]
          } else 
             w1 = NULL
          coxModel1 = try(do.call(coxph, c(list(formula = as.formula(designs[d]), data = modelData,
                                                weights = w1),
                                           coxphArgs)), silent = TRUE);
          if (!inherits(coxModel1, "try-error") && coxModel1$iter <= maxIter)
          {
            sum1 = summary(coxModel1);
            coeff1 = sum1$coefficients;
            if ("robust se" %in% colnames(coeff1))
            {
               if (coxph.checkRobustResults)
               {
                 se = max(coeff1[dataNameInDesign, c("se(coef)", "robust se")])
                 coeff1[dataNameInDesign, "z"] = coeff1[dataNameInDesign, "coef"]/se
                 coeff1[dataNameInDesign, "Pr(>|z|)"] = 2*pnorm(abs(coeff1[dataNameInDesign, "z"]), lower.tail = FALSE);
               }
            } 
            #if (replaceMissing(abs(coeff1[dataNameInDesign, "z"])) > 10) browser()
            res[g, 2:4] = coeff1[dataNameInDesign, c("coef", "z", "Pr(>|z|)")];
          } else {
            if (inherits(coxModel1, "try-error")) nFail = nFail + 1 
               else if (coxModel1$iter> maxIter) nExceedIter = nExceedIter +1;
          }
          if (verbose > 0 && g %% updateInterval == 0) pind = updateProgInd(g/nGenes, pind);
        } 
        if (verbose > 0) printFlush();
        if (verbose > 0) printFlush(spaste("Number of failed regressions: ", nFail, 
                                           "\nNumber of exceeded iterations: ", nExceedIter));
        res[, 5] = p.adjust(res[, 4], method = "fdr");
        lmFitResults[[d]] = data.frame.ncn(res)
      }
      oldDesign = designs[d];
    #}
  }
  names(lmFits) = testNames;
  names(lmFitResults) = designInfo$testName;
  out = list(results = lmFitResults);
  if (keepFullModels)
    out$fullModels = lmFits;
  out;
}

#==========================================================================================================
#
# get top genes. This function should work for both limma and deseq results.
#
#==========================================================================================================

getTopGenes.simple = function(
  Z, 
  p = 2*pnorm(abs(Z), lower.tail = FALSE),
  pAdj = p.adjust(p, method = "fdr"),
  logFC = NULL,
  selectAtLeast = 0,
  geneIDs = NULL,
  pThreshold = 0.05,
  pAdjThreshold = 0.1,
  logFCThreshold = NULL,
  direction = c("Downregulated", "Upregulated", "Differentially expressed"),
  order = FALSE)
{
  direction = match.arg(direction);
  sign = c(-1, 1, 0)[ match(direction, c("Downregulated", "Upregulated", "Differentially expressed"))];
  if (is.null(logFC) | is.null(logFCThreshold))
  {
    keepByFC = rep(TRUE, length(Z));
  } else
    keepByFC = abs(logFC) > logFCThreshold;
  index = which(Z*sign >= 0 & p < pThreshold & pAdj < pAdjThreshold & keepByFC);
  if (length(index) < selectAtLeast)
    index = which(base::rank(-Z*sign, ties.method = "first") <=selectAtLeast);
  if (order) index = index[ order(p[index])];
  if (!is.null(geneIDs)) geneIDs[index] else index;
}
                         

getTopGenes = function(testResults, 
   testNames = names(testResults), 
   thresholdName = "",
   pThreshold = 0.05,
   pAdjThreshold = 0.05,
   logFCThreshold = 0,
   # Minimum number of genes. If the number of significantly associated genes is less than this, this number
   # of top genes will be included on the list.
   topList.minN = 0,
   # Gene IDs to use. If not given, rownames of results will be used.
   geneIDs = NULL,
   backupGeneIDs = NULL,
   # Column match patterns...
   FCColPattern = "^log2FC\\.for",
   ZColPattern = "^Z\\.for",
   pColPattern = "^p\\.for",
   pAdjColPattern = "^FDR\\.for",
   DEDirectionNames = c("Downregulated", "Upregulated"),
   addDirectional = NULL,
   addCombined = NULL,
   combinedName = "Differentially.expressed",
   addSplitBy = NULL,
   geneAnnot = NULL,
   testNameSep = ".for.",
   organizeIntoPairs = FALSE,
   returnFlatList = FALSE
   )
{
   #if (FCThreshold != 0) stop("Fold change threshold is currently ignored, do not use."); ## Should work now.

   if (is.null(logFCThreshold)) logFCThreshold = 0;
   if (is.null(geneIDs)) geneIDs = rownames(testResults[[1]]);
   ZCol = grepSingleOrError(ZColPattern, colnames(testResults[[1]]));
   pCol = grepSingleOrError(pColPattern, colnames(testResults[[1]]));
   FCCol = grepSingleOrError(FCColPattern, colnames(testResults[[1]]));
   pAdjCol = grepSingleOrError(pAdjColPattern, colnames(testResults[[1]]));
   if (length(addSplitBy) == 0) addSplitBy = rep(NA, length(testResults));
   if (length(addSplitBy) != length(testResults)) 
     stop("Length of 'addSplitBy' must be 0 or agree with length of 'testResults'.");

   if (length(addDirectional) == 0) addDirectional = rep(TRUE, length(testResults));
   if (length(addDirectional) != length(testResults)) 
     stop("Length of 'addDirectional' must be 0 or agree with length of 'testResults'.");

   if (length(addCombined) == 0) addCombined = rep(FALSE, length(testResults));
   if (length(addCombined) != length(testResults)) 
     stop("Length of 'addCombined' must be 0 or agree with length of 'testResults'.");
   if (organizeIntoPairs && any(addCombined))
     stop("'organizeIntoPairs' cannot be TRUE when any of 'addCombined' is TRUE.");

   out = mymapply(function(res, addSplitBy1, addDirectional1, addCombined1, testName1)
   { 
     if (!addDirectional1 && !addCombined1) return(NULL);
     if (addDirectional1) 
     {
       Zs = list(Downregulated = res[, ZCol], Upregulated = -res[, ZCol]);
     } else Zs = list() 
     logFCs = res[, FCCol];
     if (addCombined1) Zs = c(Zs, list(DE = -abs(res[, ZCol])));
     names(Zs) = spaste(c(if (addDirectional1) DEDirectionNames else NULL, 
                    if (addCombined1) combinedName else NULL), testNameSep, 
                    testName1, thresholdName);
     directions1 = c(DEDirectionNames, if (addCombined1) combinedName else NULL);
     Zs = mymapply(function(Z1, dirWord) 
     {
       attr(Z1, "testInfo") = list(
         testName = testName1,
         thresholdName = thresholdName,
         direction = dirWord,
         splitBy = NA,
         splitDirection = NA);
       Z1;
     }, Zs, directions1);
     if (!is.na(addSplitBy1))
     {
       if (!addSplitBy1 %in% c(names(testResults), names(geneAnnot)))
          stop("When 'addSplitBy' is given, it must be one of 'names(testResults)' or 'names(geneAnnot)'. \n",
               "  addSplitBy[1]: ", addSplitBy1,
             "\n  names(testResults): ", formatLabels(paste(names(testResults), collapse = ", "),
                                                      maxCharPerLine = options("width")$width, split = ", "),
             "\n  names(geneAnnot): ", formatLabels(paste(names(geneAnnot), collapse = ", "),
                                                      maxCharPerLine = options("width")$width, split = ", "));
       splitByZ = if (addSplitBy1 %in% names(testResults)) testResults[[addSplitBy1]] [[ZCol]] else
                    geneAnnot[[addSplitBy1]];
       Zs.1 = mymapply(maskVectorByReferenceDirection,
              x = removeListNames(Zs),
              name.x = names(Zs),
              MoreArgs = list(ref = splitByZ, name.ref = addSplitBy1, includeAll = TRUE));
       Zs.1 = mymapply(function(Z1, dirWord)
         mymapply(function(Z2, splitBy, dirWord.split)
         {
           attr(Z2, "testInfo") = list(testName = testName1,
             thresholdName = thresholdName,
             direction = dirWord,
             splitBy = addSplitBy1,
             splitDirection = dirWord.split);
           Z2;
         }, Z1, splitBy = c(NA, addSplitBy1, addSplitBy1), dirWord.split = c(NA, "down", "up")), 
         Zs.1, directions1);
       if (organizeIntoPairs) {
         Zs = lapply(1:length(Zs.1[[1]]), function(i) c(Zs.1[[1]] [i], Zs.1[[2]] [i]));
       } else 
         Zs = list(unlist(Zs.1, recursive = FALSE))
     } else Zs = list(Zs);

     out = lapply(Zs, lapply, function(Z1) 
      {
        x1 = significantOrTop(geneIDs, backupIDs = backupGeneIDs, Z = Z1, 
                p = res[, pCol], pThreshold = pThreshold,
                p.adjusted = res[, pAdjCol], pAdjThreshold = pAdjThreshold, 
                zThreshold = 0, 
                extraStat = abs(logFCs), extraThreshold = logFCThreshold,
                nTop = topList.minN, warn = FALSE);
        attr(x1, "testInfo") = attr(Z1, "testInfo");
        x1;
      });
     if (!organizeIntoPairs) out = unlist(out, recursive = FALSE);
     out;
   }, testResults, addSplitBy, addDirectional, addCombined, testNames);
   out = out[sapply(out, length) > 0];
   if (organizeIntoPairs) out = unlist(removeListNames(out), recursive = FALSE);
   if (returnFlatList) out = unlist(removeListNames(out), recursive = FALSE);
   out;
}

#======================================================================================================
#
# create a collection of DE genes
#
#======================================================================================================

# input: topGenes from the main function: a list with one component per threshold; each is a list
# with one component per test/direction/split.

# geneSetNamePattern: character string with tags %t and %d, %p, %m for test, direction, significance and
# modifier
# geneSetDescriptionPattern: same 

standardGroupFile = function()
{
  home = sub("Work/.*", "", getwd());
  spaste(home, "Work/RLibs/inHouseGeneAnnotation/StandardizedGroups/standardizedGroups-general.csv");
}

modifiedTestName = function(
  testInfo, 
  pretty = FALSE, 
  modifierSeparator = if (pretty) " and " else ".and.", 
  modifierTestSeparator = if (pretty) " for " else ".for.",
  designInfo = NULL,
  prettifyList = NULL)
{
  if (pretty)
  {
    if (!is.null(designInfo))
    {
       trtab = designInfo[ c("testName", "testName.pretty")];
       testInfo$splitBy = translateUsingTable(testInfo$splitBy, trtab);
       testInfo$testName = translateUsingTable(testInfo$testName, trtab);
    } else if (!is.null(prettifyList)) {
       testInfo$splitBy = prettifyStrings(testInfo$splitBy, prettifyList);
       testInfo$testName = prettifyStrings(testInfo$testName, prettifyList);
    }
  }
  
  modifierString = if (is.na(testInfo$splitBy) || is.na(testInfo$splitDirection)) "" else 
         spaste(modifierSeparator, testInfo$splitDirection, modifierTestSeparator, 
                testInfo$splitBy);
  spaste(testInfo$testName, modifierString);
}


createCollectionOfDEGenes = function(
  topGenes, 
  keepThresholds,
  thresholdNames.pretty,
  keepTests,
  designInfo,  # For the mapping of testName and testName.pretty
  geneEvidence,
  geneSource,
  IDBase,
  geneSetNamePattern,
  geneSetDescriptionPattern,
  organism,
  internalClassification,
  geneSetGroups,
  groupList, 
  collectionFile = NULL,  # Tags: %t for threshold, %i for threhsold index
  separateCollectionsByThreshold = TRUE,
  modifierSeparator = " and ",
  modifierTestSeparator = " for ",
  dropEmptySets = TRUE,
  addToAnalysisManager = FALSE,
  analysisID = NULL,
  amComponent)
{
  stdTags =c("%t", "%d", "%p");
  #nTopThresholds = length(topGenes);
  geneSets = list();
  indexStart = 1;
  geneListInfo = lapply(topGenes, lapply, attr, "testInfo");
  testNames = sapply(geneListInfo[[1]], getElement, "testName");
  if (any(!testNames %in% designInfo$testName)) 
    stop("Some elements of 'testNames' are not among 'designInfo$testName':\n",
         paste(testNames[!testNames %in% designInfo$testName], collapse = ", "));
  if (is.null(keepTests)) keepTests = unique(testNames);
  assocCollections = list();
  collIndex = 1;
  stdGroups = anRichmentMethods::standardizedGroups();
  collectionFiles = character(0);
  for (tt in which(names(topGenes) %in% keepThresholds))
  {
    nGenes = sapply(topGenes[[tt]], length);
    testIndex = which(testNames %in% keepTests & (nGenes > 0 | !dropEmptySets));
    n1 = length(testIndex);
    if (n1==0) next;
    geneSets = mymapply(function(entrez, testInfo, setIndex, prettyTest)
    {
      prettyTest.x = modifiedTestName(testInfo, pretty = TRUE, designInfo = designInfo);
      newGeneSet(geneEntrez = entrez,
                 geneEvidence = geneEvidence,
                 geneSource = geneSource,,
                 ID = spaste(IDBase, ".", prependZeros(setIndex, 4)),
                 name = WGCNA:::.substituteTags(geneSetNamePattern, stdTags, 
                                      c(prettyTest.x, testInfo$direction, thresholdNames.pretty[tt])),
                 description = WGCNA:::.substituteTags(geneSetDescriptionPattern, stdTags,
                                    c(prettyTest.x, testInfo$direction, thresholdNames.pretty[tt])),
                 source = geneSource,
                 internalClassification = internalClassification,
                 groups = geneSetGroups,
                 organism = organism,
                 lastModified = Sys.Date())
      
    }, topGenes[[tt]] [testIndex], testInfo = geneListInfo[[tt]] [testIndex],
       setIndex = indexStart:(indexStart + n1 -1),
       prettyTest = translateUsingTable(testNames[testIndex], designInfo[ c("testName", "testName.pretty")]));
    
    assocCollection = dropUnreferencedGroups(newCollection(geneSets, groups = c(stdGroups, groupList)),
                                             verbose = 0);
    if (separateCollectionsByThreshold && !is.null(collectionFile)) 
    {
      file1 = multiSub(c("%i", "%t"), c(collIndex, gsub("_", "", names(topGenes)[tt])), collectionFile);
      save(assocCollection, file = file1);
      collectionFiles[[collIndex]] = file1;
      if (addToAnalysisManager)
         am.addNewObjectToAnalysis(
             name = c("anRichment collection", thresholdNames.pretty[tt]),
             fileName = file1,
             fileObject = "assocCollection",
             analysisID = analysisID,
             anaType = amComponent,
             allowNonstandardAnalysisType = TRUE,
             allowNonstandardObjectType = TRUE,
             returnBoth = FALSE, updateAnalysisManager = TRUE);
    }
    assocCollections[[collIndex]] = assocCollection;
    indexStart = indexStart + n1;
    collIndex = collIndex + 1;
  }
  if (separateCollectionsByThreshold) {
    list(collections = assocCollections, files = collectionFiles);
  } else {
    assocCollection = do.call(mergeCollections, assocCollections);
    save(assocCollection, file = collectionFile);
    if (addToAnalysisManager)
       am.addNewObjectToAnalysis(
           name = c("anRichment collection"),
           fileName = collectionFile,
           fileObject = "assocCollection",
           analysisID = analysisID,
           anaType = amComponent,
           allowNonstandardAnalysisType = TRUE,
           allowNonstandardObjectType = TRUE,
           returnBoth = FALSE, updateAnalysisManager = TRUE);
    list(collections = list(assocCollection), files = collectionFile);
  }
    
}

defaultOverallNormFunction = function(pipeline)
{
  switch(pipeline,
    DESeq2 = "getOverallNormalizedData.DESeq2",
    `limma-voom` = "getOverallNormalizedData.limma",
    stop("No default normalization for pipeline ", pipeline));
}

defaultOverallNormArgs = function(overallNormalizedDataFnc)
{
  if (length(overallNormalizedDataFnc)==0) return(list());
  if (length(overallNormalizedDataFnc) > 1) 
     stop("defaultOverallNormalizedCountArgs: 'overallNormalizedDataFnc' must be of length no more than 1.");
  switch(overallNormalizedDataFnc,
   getOverallNormalizedData.DESeq2 = list(design = "~1"),
   getOverallNormalizedData.limma = list(normalizationMethod = "TMM"),
   stop("defaultOverallNormalizedCountArgs: Do not have arguments for function ", overallNormalizedDataFnc));
}



defaultDEAnalysisFnc = function(analysisPipeline)
{
  switch(analysisPipeline,
    DESeq2 = "runDESeq",
    `limma-voom` = "runLimmaVoom",
    `limma-lm` = "runLimma",
    `custom` = stop("'custom' pipeline requires user-supplied 'DEAnalysisFnc' and 'DEAnalysisArgs'."),
    stop("unrecognized 'analysisPipeline': ", analysisPipeline))
}


defaultDEAnalysisArgs = function(DEAnalysisFnc)
{
  if (!is.character(DEAnalysisFnc))
    stop("defaultDEAnalysisArgs: argument must be a character string.");

  switch(DEAnalysisFnc,
    runDESeq = list(
            sizeFactorArgs = list(type = "ratio", locfunc = stats::median),
            DESeqArgs = list(minReplicatesForReplace = Inf, fitType = "parametric", quiet = TRUE),
            resultsArgs = list(cooksCutoff = Inf, independentFiltering = FALSE),
            parallel = FALSE,
            dropColumns = "lfcSE",
            addColumns = character(0)),
    runLimma = list(
            lmFitArgs = list(),
            topTableArgs = list(number = Inf, sort.by = "none"),
            dropColumns = c("t", "F", "B", "genelist"),
            voom = FALSE),
    runLimmaVoom = list(
            lmFitArgs = list(),
            topTableArgs = list(number = Inf, sort.by = "none"),
            dropColumns = c("t", "F", "B", "genelist"),
            normalizationMethod = "TMM",
            voom = TRUE),
     custom = stop("'custom' pipeline requires user-supplied 'DEAnalysisFnc' and 'DEAnalysisArgs'."),
     stop("unrecognized 'DEAnalysisFnc': ", DEAnalysisFnc))
}


getOverallNormalizedData.DESeq2 = function(data, pheno, design, ...)
{
  dsData = estimateSizeFactors(DESeqDataSetFromMatrix.PL(t(data),
                   colData = pheno, design = as.formula(design)));
  t(counts(dsData, normalized = TRUE));
}

getOverallNormalizedData.limma = function(data, normalizationMethod = "TMM", ...)
{
  dge = DGEList(counts = t(data));
  dge <- calcNormFactors(dge, method = normalizationMethod);
  f = dge$samples$lib.size * dge$samples$norm.factors;
  normCounts = data/f;
  normCounts;
}

#======================================================================================================
#
# Main function
#
#======================================================================================================


individualAssociation.general = function(
   data,
   pheno,

   dataWeights = NULL,

   # choice of analysis pipeline
   # Choices: DESeq2, limma-lm, limma-voom, custom
   # Default: DESeq2 for data that are integer.

   analysisPipeline = NULL,

    
   # Flow control, for some measure of interactivity
   intermediateResults = list(),
   stopAt = c("End", "10-individual association", "30-Combine results",
              "40-Gene lists", "50-Enrichment"),
  
   dataIsLogTransformed = NULL,
   logBase = 2,
   logShift = 1,

   # Analysis name
   analysisName = "",
   organism,
   tissue = NULL,

   # analysis manager options
   addToAnalysisManager = TRUE,
   amAnalysisName,
   amAnalysisDescription,
   amMetaData = list(),
   amDataSource,
   amAnalysisName.short = analysisName,
   amAnalysisName.pretty = amAnalysisName,
   amComponent = "individual analysis",

   # A data frame that contains the following information
   designInfo,
   # For general referring to the tests. The following columns:
   # testName,
   # design,
   #   A character vector with formulas in character forms. Each test should have its own design even if
   #   designs are repeated. designs can have names that will be used in combining and contrasting results.
   # testType:
   #   for DESeq2: one of "Wald" or "LRT"; for limma: one of "F" or "t"
   # reducedDesign:
   #   the design of the reduced model for LRT, otherwise ignored.
   # contrast.name, contrast.numerator, contrast.denominator
   #   name, numerator, denominator for each contrast. If NA, the corresponding entry of coefficientName will
   #   be used
   # coefficientName
   #   Names of coefficients (from resultsNames).
   # testName.pretty
   #   Pretty test names.

   getOverallNormalizedData = TRUE,
   overallNormalizedDataFnc = defaultOverallNormFunction(analysisPipeline),
   overallNormalizedDataArgs = defaultOverallNormArgs(overallNormalizedDataFnc),
   #overallCountDesign = "~1",
   #getNormalizedCountsPerDesign = TRUE,


   DEAnalysisFnc = defaultDEAnalysisFnc(analysisPipeline),
   DEAnalysisArgs = defaultDEAnalysisArgs(DEAnalysisFnc),

   keepFullModels = FALSE,
   testNameSep = ".for.",
   testNameSep.pretty = " for ",

   levelSep = ".vs.",
   levelSep.pretty = " vs. ",

   ## Arguments for running DESeq
   #DESeqSizeFactorArgs = list(type = "ratio", locfunc = stats::median),
   #DESeqArgs = list(minReplicatesForReplace = Inf, 
   #                 fitType = "parametric", quiet = TRUE),
   #keepFullDESeqObjects = FALSE,
#
#   resultsArgs = list(cooksCutoff = Inf, independentFiltering = FALSE),
#   dropDESeqResultColumns = c("lfcSE"),
   #addDESeqColumns = character(0),

   # Arguments for running limma's lmFit
   #lmFitArgs = list(),
   #keepFullLmFitObjects = FALSE,
   #topTableArgs = list(number = Inf, sort.by = "none"),
   #dropColumns = c("t", "F", "B", "genelist"),

   changePAdjTo = "FDR",
   changeCoeffTo = "log2FC",

   # Additional info to include in each result table
   geneAnnotation.start,
   geneAnnotation.end = NULL,
   # Additional gene information potentially useful for splitSignificantHitsBy. 
   # This table is not included in the results.
   geneAnnotation.extra = NULL,

   # features that are to be excluded from results
   excludeFeatures = NULL,

   # Creating a collection of DE genes  
   createCollectionOfSignificantGenes = TRUE,
   collectionGeneIDcolName = "Entrez",
   DEDirectionNames = c("Downregulated", "Upregulated"),
   collection.keepTests = NULL,
   collection.keepThresholds = topThresholdNames[1],
   separateCollectionsByThreshold = length(collection.keepThresholds) > 1,
   collection.dropEmptySets = TRUE,
   geneEvidence = "IEP",
   geneSource,
   IDBase = analysisName,
   # character string with tags %t and %d, %p for test, direction and significance.
   # If it exists in designInfo, testName.pretty.ext will be used.
   geneSetNamePattern,
   geneSetDescriptionPattern,
   internalClassification,
   geneSetGroups,
   groupList,
   collectionFile = spaste("RData/collectionOfDEGenes", prependDash(analysisName), ".RData"),


   # A list of vectors (corresponding to design$testName) specifying how to combine the results. 
   combineResults = list(all = designInfo$testName),
   combineNames = names(combineResults),
   combineNames.pretty = combineNames,
   combineColNameExtensions = if (length(combineResults)==1) "" else combineNames,
   combineNameSep = if (length(combineColNameExtensions)==1) "" else ".in.",
   combineNameSep.pretty = if (length(combineColNameExtensions)==1) "" else " in ",
   interleaveResults = FALSE,
   dropRepeatedColumns = TRUE,
   dropRepeated.exceptions = c(changePAdjTo),
   # Combined results postprocessing fnc. Results will replace the combined ones, so the function needs to
   # return its results argument, possibly modified or with extra columns. Arguments should include the deseq
   # objects and one data frame of combined results (res), as well as any other arguments. Also supplied are
   # topThresholds and topThresholdNames but the function may ignore them.
   combResultsPPFnc = NULL,
   combResultsPPArgs = list(),

   # Creating list of top genes for each test
   topThresholds = list(._FDR.lessThan.0.05._ = list(pAdj = 0.05, p = 1, logFC = 0),
                        ._p.lessThan.0.05._ = list(pAdj = 2, p = 0.05, logFC = 0)),

   topThresholdNames = names(topThresholds),
   thresholdNames.pretty = NULL,  # if not given, will be generated using the complete prettifyList

   # Minimum number of genes. If the number of significantly associated genes is less than this, this number
   # of top genes will be included on the list.
   topList.minN = 0,
   # The column in geneAnnotation.start to use for the list
   topList.machineColName = "Entrez",
   topList.humanColName = "Symbol",
   topList.backupColName = topList.humanColName,

   # Add at most this number of genes to tables of numbers of significant genes
   maxGenesInTables = 15,
   maxmaxGenesInTables = 20,

   # Additional gene lists can be generated from the (extended) combined results. This function should take 
   # the full combined and  extended results, and return a named list of gene lists. A pretty name should be
   # specified in "prettyName" attribute.
   # The function is also expected to take argument topThresholds which is a list of significance thesholds
   # for identifying top genes (but may ignore them if they are not necessary).
   geneSetFnc = NULL,
   geneSetFncArgs = list(),

   # Enrichment analysis
   enrichment.useResultsCombinations = NULL,
   enrichment.useThresholds = topThresholdNames,
   enrichment.collections = NULL,
   enrichment.combinedCollection = do.call(mergeCollections, 
              c(enrichment.collections, list(stopOnDuplicates = TRUE))),
   enrichment.annotIDCol = "Entrez",
   enrichment.minN = 100,
   # Optional specification of DE tests that should not be included in enrichment results. 
   enrichment.dropTests = NULL,
   enrichment.otherArgs = list(getOverlapEntrez = FALSE, 
                               getOverlapSymbols = TRUE,
                               maxReportedOverlapGenes = 10000,
                               nBestDataSets = 5,
                               getDataSetDetails = FALSE,
                               useBackground = "given"),
   enrichment.dropCols = character(0), # "dataSetID",
   enrichment.compress = TRUE,
   enrichment.combinedCollectionName = "CombinedCollection",

   minimizeEnrichmentSize = FALSE,
   
   # general result export options
   resultDir = "Results",
   exportDigits = 3,

   # Export options for result tables
   exportIndividualResults = FALSE,
   addGeneAnnotationToExportedIndividualResults = TRUE,

   exportCombinedResults = TRUE,
   separateExtendedResults = TRUE,
   associationDir = file.path(resultDir, "01-Association"),
   resultFileBase = spaste("associationWithPhenotypes", prependDash(analysisName)),
   extendedResultFileBase = spaste("associationWithPhenotypes-extended", prependDash(analysisName)),
   topGeneFileBase = spaste("topGenes", prependDash(analysisName)),
   testColName = "Test",

   exportExcludedFeatures = TRUE,
   excludedFeaturesResultFileBase = spaste("associationWithPhenotypes-technicalFeatures", prependDash(analysisName)),

   exportNSignificant.perThreshold = FALSE,
   exportNSignificant.combinedThresholds = TRUE,

   # Export options for enrichment results
   exportEnrichmentTables = TRUE,
   enrichmentDir = file.path(resultDir, "02-Enrichment"),
     # Use %combine and %threshold to add names of results combinations and threshold name to enrichmentDir
     # analysisName can be added as an argument here.
   enrichmentFileBase = spaste("enrichmentOfDEGenes-", analysisName),

   # Translation table between programmatic and human-readable names
   prettifyList = data.frame(character(0), character(0)),
   extendPrettifyList = FALSE,

   # Options for generation of plots.
   plotDir = "Plots",
   plotFnc = "pdf",
   plotExt = ".pdf",

   # behavior
   collectGarbage = TRUE,
   forceCalculation = FALSE,

   # forceExport is useful if a change in prettifyList/pretty names is the only reason to re-run calculation.
   forceExport = FALSE,

   # Parallel calculation options
   runParallel = FALSE,
   nWorkers = 2,

   verbose = 1,
   indent = 0
)
{
   prettifyNames.loc = function(df) prettifyNames(df, prettifyList);

   stopAt = match.arg(stopAt);
   spaces = indentSpaces(indent);
   plotFnc = match.fun(plotFnc);

   if (is.null(analysisPipeline))
   {
     if (all(round(data)==data, na.rm = TRUE))
     {
       analysisPipeline = "DESeq2"
     } else
       analysisPipeline = "limma-lm";
   }
   recognizedPipelines = c("DESeq2", "limma-lm", "limma-voom", "custom");
   anaPipe = recognizedPipelines[charmatch(tolower(analysisPipeline), tolower(recognizedPipelines))]
   printFlush(spaste(spaces, "Using analysis pipeline ", anaPipe));
   countPipeline = anaPipe %in% c("DESeq2", "limma-voom");

   if (is.na(anaPipe))
       stop("If given, 'analysisPipeline' must be one of\n   ", paste(recognizedPipelines, collapse = ", "));

   if (anaPipe == "DESeq2") storage.mode(data) = "integer"

   if (is.null(dataIsLogTransformed))
   {
     if (countPipeline) dataIsLogTransformed = FALSE;
     if (anaPipe=="lm-Limma") dataIsLogTransformed = max(data, na.rm = TRUE) < 100;
     printFlush(spaste(spaces, "Assuming data is", if (dataIsLogTransformed) "" else " not", " log-transformed."));
   }

   designInfo = as.data.frame(designInfo);
   designs = designInfo$design;
   nDesigns = length(designs);

   # Fill in missing designInfo entries with default values.
   if (is.null(designInfo$coefficientName)) designInfo$coefficientName = rep(NA, nrow(designInfo));
   if (is.null(designInfo$testType)) designInfo$testType = rep("Wald", nrow(designInfo));
   if (is.null(designInfo$combineSignificantHits)) designInfo$combineSignificantHits = rep(FALSE, nrow(designInfo));
   if (is.null(designInfo$splitSignificantHitsBy)) designInfo$splitSignificantHitsBy = rep(NA, nrow(designInfo));
   if (is.null(designInfo$getDirectionalHits)) 
     designInfo$getDirectionalHits = with(designInfo, 
         testType %in% c("Wald", "CoxPH", "t") | !is.na(coefficientName) | !is.na(contrast.name));

   nTopThresholds = length(topThresholds);

   testNames = designInfo$testName;
   testNames.pretty = designInfo$testName.pretty;
   testNames.pretty.ext = 
      if (is.null(designInfo$testName.pretty.ext)) designInfo$testName.pretty else designInfo$testName.pretty.ext;

   if (is.null(combineNames))
   {
     if (isTRUE(all.equal(combineResults, testNames))) {
       names(combineResults) = testNames
     } else stop("When 'combineResults' is given and different from 'testNames', it must have names.");
   } 

   if (length(combineNames)!=length(combineResults))
      stop("Lengths of 'combineResults' and 'combineNames' must be the same.");

   if (length(combineColNameExtensions)!=length(combineResults))
      stop("Lengths of 'combineResults' and 'combineColNameExtensions' must be the same.");
   names(combineResults) = combineNames;

   if (any(sapply(combineResults, length)==0)) 
   {
     printFlush("'combineResults' contains empty components. These will be removed.");
     keep = sapply(combineResults, length) > 0;
     if (sum(keep)==0) stop("'combineResults' does not have non-empty components.");
     combineResults = combineResults[keep];
     combineNames = combineNames[keep];
     combineColNameExtensions = combineColNameExtensions[keep];
   }
     

   if (length(testNames)!=length(designs))
     stop("Lengths of 'testNames' and 'designs' must be the same.");

   if (is.null(names(designs))) names(designs) = testNames;

   rownames(pheno) = rownames(data);

   if (!is.null(dataWeights))
   {
     if (!is.null(dim(dataWeights)))
     {
       if (!isTRUE(all.equal(dim(dataWeights), dim(data))))
          stop("When 'dataWeights' is an array, it must have the same dimensions as 'data'.");
     } else if (length(dataWeights)!=nrow(data) && length(dataWeights)!=ncol(data))
          stop("When 'dataWeights' is a dimensionless vector, its length must equal the number\n",
               "of genes or samples in 'data'.");
     if (any(!is.finite(dataWeights)))
     {
       warning(immediate. = TRUE,
               "Some 'dataWeights' are missing or not finite; these will be set to zero.");
       dataWeights[!is.finite(dataWeights)] = 0;
     }
   }

   goodSG = goodSamplesGenes(data, weights = dataWeights, verbose = verbose -1, indent = indent + 1);
   goodSG$colNames = colnames(data);

   if (is.null(enrichment.useResultsCombinations))
      enrichment.useResultsCombinations = combineNames;

   if (any(!enrichment.useResultsCombinations %in% combineNames))
     stop("If given, entries in 'enrichment.useResultsCombinations' must also be present in 'combineNames'.\n",
          "The following are not: \n    ", 
            paste(enrichment.useResultsCombinations[!enrichment.useResultsCombinations %in% combineNames],
                  collapse = "\n    "));
 
   enrichment.useCRIndex = match(enrichment.useResultsCombinations, combineNames);
  
   enrThresholdIndex = match(enrichment.useThresholds, topThresholdNames);
   if (any(is.na(enrThresholdIndex)))
     stop("'enrichment.useThresholds' must contain valid top threshold names.");
   nEnrichmentThresholds = length(enrichment.useThresholds);

   testNames = make.names(testNames);
   combineNames = make.names(combineNames);
   analysisName = if (analysisName!="") substring(make.names(spaste("X", analysisName)), 2) else "";
      ### This construct should get rid of all non-standard characters in the middle of analysisName but leave a leading
      ### number if that's present (e.g., 5xFAD)

   if (length(prettifyList)>0)
   {
     prettifyList = as.data.frame(prettifyList)
     names(prettifyList) = c("from", "to");
   }
   if (extendPrettifyList)
   {
     prettifyList = data.frame(from=c(prettifyList[[1]], testNames, combineNames),
                         to=c(prettifyList[[2]], testNames.pretty, combineNames.pretty));
   }

   prettifyList = as.data.frame(removeDuplicatesFromPrettifyList(
                     rbind.ncn(prettifyList, 
                                    data.frame( c(testNameSep, combineNameSep, levelSep),
                                          c(testNameSep.pretty, combineNameSep.pretty, levelSep.pretty)))))
   # Clean prettifyList
   keep = prettifyList[[1]]!="" & prettifyList[[1]] != prettifyList[[2]];
   prettifyList = prettifyList[keep, ];

   order1 = order(-nchar(prettifyList[[1]]));
   prettifyList = prettifyList[order1, ];
    
   prettifyList.plots = prettifyList;
   prettifyList.plots[[2]] = fix.dN17(prettifyList[[2]]);
   prettifyList.plots = rbind.ncn(prettifyList.plots, c("dN17", "\U0394N17"));

   intermediateResults$analysisName = analysisName;
   intermediateResults$testNames = testNames;
   intermediateResults$testNames.pretty = testNames.pretty;
   intermediateResults$combineNames.pretty = combineNames.pretty;
   intermediateResults$combineNames = combineNames;
   intermediateResults$prettifyList = prettifyList;
   intermediateResults$prettifyList.plots = prettifyList.plots;

   intermediateResults$topThresholds = topThresholds;
   intermediateResults$topThresholdNames = topThresholdNames;
   intermediateResults$goodSamplesGenes = goodSG;
   
   if (is.null(thresholdNames.pretty)) 
     thresholdNames.pretty = prettifyStrings(topThresholdNames, prettifyList);

   if (!is.null(intermediateResults$designs))
   {
     if (!isTRUE(all.equal(as.character(designs), as.character(intermediateResults$designs))))
     {
       printFlush(spaste(spaces, "'designs' have changed from ones saved in 'intermediateResults'.\n",
                         spaces, "Will start calculations from the beginning."));
       forceCalculation = TRUE;
     }
   } else {
     intermediateResults$designs = designs;
   }

   if (is.null(geneAnnotation.start)) {
      geneAnnotation.start = data.frame(ID = colnames(data));
   } else {
      if (nrow(geneAnnotation.start)!=ncol(data)) 
        stop("Number of rows in 'geneAnnotation.start' must equal number of columns in 'data'.");
      if (!isTRUE(all.equal(as.character(geneAnnotation.start[[1]]), as.character(colnames(data)))))
        warning(immediate. = TRUE,
                "**** Possible discrepancy between 'geneAnnotation.start' and columns of 'data':\n",
                "****   First column of 'geneAnnotation.start' does not agree with 'colnames(data)'.");
   }

   if (!is.null(geneAnnotation.end)) {
      if (nrow(geneAnnotation.start) != nrow(geneAnnotation.end) )
        stop("Number of rows in 'geneAnnotation.start' and 'geneAnnotation.end' must be the same.");
   }

   dir.create(associationDir, recursive = TRUE, showWarnings = FALSE);

   if (addToAnalysisManager)
   {
      am.analysisID = am.analysisID(name = amAnalysisName, throw = FALSE);
      if (is.na(am.analysisID))
      {
        am.analysisID = am.insertNewAnalysis(name = amAnalysisName,
           description = amAnalysisDescription,
           analysisType = amComponent,
           allowNonstandardType = TRUE,
           organism = organism,
           tissue = tissue,
           shortName = amAnalysisName.short,
           prettyName = amAnalysisName.pretty,
           metaData = amMetaData,
           dataSource = amDataSource);
      } 
   }

   if (!is.null(excludeFeatures))
   {
     if (is.logical(excludeFeatures) || is.numeric(excludeFeatures)) excludeFeatures = colnames(data)[excludeFeatures];
     if (any(!excludeFeatures %in% colnames(data)))
     {
        warning(.immediate = TRUE, "Some elements of 'excludeFeatures' are not among 'colnames(data)'.");
     }
     excludedFeatureIndex = match(excludeFeatures, colnames(data));
     includedFeatureIndex = setdiff(1:ncol(data), excludedFeatureIndex);
   } else {
     excludedFeatureIndex = numeric(0);
     includedFeatureIndex = 1:ncol(data);
   }
   geneAnnotation.start.included = geneAnnotation.start[includedFeatureIndex, , drop = FALSE];
   if (length(geneAnnotation.end) > 0) 
      geneAnnotation.end.included = geneAnnotation.end[includedFeatureIndex, , drop = FALSE];
   if (length(geneAnnotation.extra) > 0) 
      geneAnnotation.extra.included = geneAnnotation.extra[includedFeatureIndex, , drop = FALSE];

   #-------------------------------------------------------------------------------------------------------
   # Individual association screening
   #-------------------------------------------------------------------------------------------------------

   if (forceCalculation || is.null(intermediateResults$indivResults))
   {
     if (verbose > 0) printFlush(spaste(spaces, "Running ", DEAnalysisFnc))
     indivResults00 = do.call(DEAnalysisFnc, c(
            list(data = data, pheno = pheno, dataWeights = dataWeights, designInfo = designInfo,
                 keepFullModels = keepFullModels, 
                 changePAdjTo = changePAdjTo, changeCoeffTo = changeCoeffTo,
                 testNameSep = testNameSep,
                 goodSamplesGenes = goodSG,
                 verbose = verbose, indent = indent),
            DEAnalysisArgs));

     if (collectGarbage) gc();
     if (keepFullModels) intermediateResults$fullModels = indivResults0$fullModels;
     #if (getNormalizedCountsPerDesign) {
     #     normalizedCounts.perDesign = getNormalizedCounts(deseqObjects, designs);
     #   } else 
            normalizedCounts.perDesign = NULL;
     indivResults0 = indivResults00$results;

     if (getOverallNormalizedData)
     {
        overallNormalizedData = do.call(overallNormalizedDataFnc, 
                 c(list(data = data, pheno = pheno), overallNormalizedDataArgs));
     } else overallNormalizedData = NULL
     indivResults.excludedFeatures = lapply(indivResults0, function(x) x[excludedFeatureIndex, ]);
     indivResults = if (length(excludedFeatureIndex) > 0) lapply(indivResults0, function(x) x[includedFeatureIndex, ]) else
                     indivResults0;
     intermediateResults$indivResults0 = indivResults0;
     intermediateResults$indivResults = indivResults;
     intermediateResults$indivResults.excludedFeatures = indivResults.excludedFeatures;
     intermediateResults$designInfo = designInfo;
     intermediateResults$normalizedCounts.perDesign = normalizedCounts.perDesign; 
     intermediateResults$overallNormalizedData = 
        if (getOverallNormalizedData) overallNormalizedData else NULL;
     if (stopAt=="10-individual association") return(intermediateResults)
     forceCalculation = TRUE;
   } else {
     if (anaPipe == "DESeq2")
     {
       deseqObjects = intermediateResults$deseqObjects;
     } else {
       lmFits = intermediateResults$lmFitObjects;
     }
     normalizedCounts.perDesign = intermediateResults$normalizedCounts.perDesign;
     overallNormalizedData = intermediateResults$overallNormalizedData;
     if (is.null(overallNormalizedData) && getOverallNormalizedData)
     { 
        overallNormalizedData = do.call(overallNormalizedDataFnc,
                 c(list(data = data, pheno = pheno), overallNormalizedDataArgs));
     }
     indivResults = intermediateResults$indivResults;
     indivResults0 = intermediateResults$indivResults0;
     indivResults.excludedFeatures = intermediateResults$indivResults.excludedFeatures;
   }

   if (exportIndividualResults)
     for (d in 1:nDesigns)
     {
       f = file.path(associationDir, spaste(resultFileBase, prependDash(testNames[d]), ".csv.gz"));
       write.csv(prettifyNames.loc(
                   dataForTable(signifNumeric(indivResults[[d]], exportDigits), 
                                transpose = FALSE, IDcolName = "ID")),
                 file = gzfile(f),
                 quote = TRUE, row.names = FALSE);
       if (addToAnalysisManager) 
         am.addNewObjectToAnalysis(
            name = c("individual results", designInfo$testName[d]),
            fileName = f,
            fileReadArgs = c(defaultReadArguments("csv.gz"), list(check.names = FALSE)),
            analysisID = am.analysisID,
            anaType = amComponent,
            allowNonstandardAnalysisType = TRUE,
            allowNonstandardObjectType = TRUE,
            returnBoth = FALSE, updateAnalysisManager = TRUE);
       if (length(excludedFeatureIndex) > 0 && exportExcludedFeatures)
       {
         f = file.path(associationDir, spaste(excludedFeatureResultFileBase, prependDash(testNames[d]), ".csv.gz"));
         write.csv(prettifyNames.loc(
                     dataForTable(signifNumeric(indivResults.excludedFeatures[[d]], exportDigits),
                                  transpose = FALSE, IDcolName = "ID")),
                   file = gzfile(f),
                   quote = TRUE, row.names = FALSE);
         if (addToAnalysisManager) am.addNewObjectToAnalysis(
            name = c("individual results for technical features", designInfo$testName[d]),
            fileName = f,
            fileReadArgs = c(defaultReadArguments("csv.gz"), list(check.names = FALSE)),
            analysis = am.analysisID,
            anaType = amComponent,
            allowNonstandardAnalysisType = TRUE,
            allowNonstandardObjectType = TRUE,
            returnBoth = FALSE, updateAnalysisManager = TRUE);
       }
     }
              
   if (stopAt=="20-individual results") 
   {
     return(intermediateResults);
   }
   if (collectGarbage) gc();

   #-------------------------------------------------------------------------------------------------------
   # Combining (organizing) individual results 
   #-------------------------------------------------------------------------------------------------------

   if (forceCalculation || (!"combinedDESeqResults" %in% names(intermediateResults)) || 
           length(intermediateResults[["combinedDESeqResults"]])==0)
     ## Note to self: intermediateResults$combinedDESeqResults can match intermediateResults$combinedDESeqResults or
     ## intermediateResults$combinedDESeqResults.ext, therefore the more elaborate testing for presence of valid
     ## combinedDESeqResults
   {
     if (verbose > 0) printFlush(spaste(spaces, "Combining results..."));
     combinedDESeqResults0 = mymapply(function(index, name, colNameExt)
     {
       comb = removeListNames(indivResults0[index]);
       comb = lapply(comb, as.data.frame)
       if (interleaveResults) {
	 out = interleave(comb, sep = "");
       } else out = do.call(cbind, comb);
       if (dropRepeatedColumns)
         out = dropRepeatedColumns(out, exceptions = dropRepeated.exceptions);
       nc1 = ncol(out);
       if (length(combResultsPPFnc) > 0)
       {
         out = do.call(match.fun(combResultsPPFnc),
                           c(list(normalizedData = overallNormalizedData,
                                  dataIsLogTransformed = dataIsLogTransformed,
                                  logBase = logBase,
                                  logShift = logShift,
                                  pheno = pheno,     
                                  designInfo = designInfo,
                                  res = out, testNames = index,
                                  combName = name,
                                  topThresholds = topThresholds,
                                  topThresholdNames = topThresholdNames),
                             combResultsPPArgs));
       }
       nc2 = ncol(out);
       if (colNameExt!="") colnames(out) = spaste(colnames(out), combineNameSep, colNameExt); 
       out = data.frame(geneAnnotation.start, out, check.names = FALSE);
       attr(out, "combName") = name; 
       attr(out, "colNameExtension") = colNameExt; 
       attr(out, "constituentTests") = index
       attr(out, "columnsAddedInPP") = nc1+ ncol(geneAnnotation.start) + seq_len(nc2-nc1);
       out;
     }, combineResults, combineNames, combineColNameExtensions)
     
     if (length(excludeFeatures)>0)
     {
       combinedDESeqResults = lapply(combinedDESeqResults0, function(x) x[includedFeatureIndex, ]);
       combinedDESeqResults.excludedFeatures = lapply(combinedDESeqResults0, function(x) x[excludedFeatureIndex, ])
     } else {
       combinedDESeqResults = combinedDESeqResults0;
       combinedDESeqResults.excludedFeatures = NULL;
     }
     if (!is.null(geneAnnotation.end)) {
       combinedDESeqResults.ext = lapply(combinedDESeqResults, data.frame, geneAnnotation.end.included, check.names = FALSE);
     } else 
       combinedDESeqResults.ext = combinedDESeqResults;
     if (separateExtendedResults) {
       intermediateResults$combinedDESeqResults.ext = combinedDESeqResults.ext;
     } else {
       combinedDESeqResults = combinedDESeqResults.ext;
     }
     nCombineResults = length(combineResults);
     intermediateResults$combinedDESeqResults = combinedDESeqResults;
     intermediateResults$combinedDESeqResults.excludedFeatures = combinedDESeqResults.excludedFeatures;
     intermediateResults$nCombineResults = nCombineResults;
     if (collectGarbage) gc();
     forceCalculation = TRUE;
   } else {
     combinedDESeqResults = intermediateResults$combinedDESeqResults;
     combinedDESeqResults.ext = intermediateResults$combinedDESeqResults.ext;
     combinedDESeqResults.excludedFeatures = intermediateResults$combinedDESeqResults.excludedFeatures;
     nCombineResults = intermediateResults$nCombineResults;
   }

   if (exportCombinedResults && (forceCalculation || forceExport))
       for (cr in 1:nCombineResults)
       {
         # resultFileBase already contains analysis name in it...
         f = file.path(associationDir, spaste(resultFileBase, prependDash(combineNames[cr]), ".csv.gz"));
	 write.csv(prettifyNames.loc(signifNumeric(combinedDESeqResults[[cr]], exportDigits)),
		   file = gzfile(f), quote = TRUE, row.names = FALSE);

         if (addToAnalysisManager) am.addNewObjectToAnalysis(
            name = c("combined individual results", combineNames[cr]),
            fileName = f,
            fileReadArgs = c(defaultReadArguments("csv.gz"), list(check.names = FALSE)),
            analysisID = am.analysisID,
            anaType = amComponent,
            allowNonstandardAnalysisType = TRUE,
            allowNonstandardObjectType = TRUE,
            returnBoth = FALSE, updateAnalysisManager = TRUE);

         if (separateExtendedResults)
         {
           f = file.path(associationDir, spaste(extendedResultFileBase, prependDash(combineNames[cr]), ".csv.gz"));
           write.csv(prettifyNames.loc(signifNumeric(combinedDESeqResults.ext[[cr]], exportDigits)),
                     file = gzfile(f), quote = TRUE, row.names = FALSE);

           if (addToAnalysisManager) am.addNewObjectToAnalysis(
              name = c("extended combined individual results", combineNames[cr]),
              fileName = f,
              fileReadArgs = c(defaultReadArguments("csv.gz"), list(check.names = FALSE)),
              analysisID = am.analysisID,
              anaType = amComponent,
              allowNonstandardAnalysisType = TRUE,
              allowNonstandardObjectType = TRUE,
              returnBoth = FALSE, updateAnalysisManager = TRUE);
         }
         if (length(excludeFeatures) > 0 && exportExcludedFeatures)
         {
           f = file.path(associationDir, spaste(excludedFeaturesResultFileBase, prependDash(combineNames[cr]), ".csv.gz"));
           write.csv(prettifyNames.loc(signifNumeric(combinedDESeqResults.excludedFeatures[[cr]], exportDigits)),
                     file = gzfile(f), quote = TRUE, row.names = FALSE);
         }
       }


   if (stopAt=="30-Combine results") return(intermediateResults)

   #-------------------------------------------------------------------------------------------------------
   # Generate and optionally export lists of DE genes
   #-------------------------------------------------------------------------------------------------------

   geneAnnotation.all = data.frame.nonNULL(geneAnnotation.start, geneAnnotation.end, geneAnnotation.extra, 
                            check.names = FALSE);
   geneAnnotation.all.included = geneAnnotation.all[includedFeatureIndex, ];
   if (forceCalculation || is.null(intermediateResults$topOrNSymbols.combined)
       || is.null(intermediateResults$topSymbols.individual))
   {
     topSymbols.individual = mymapply(function(tt, tName)
       getTopGenes(indivResults, 
        geneAnnot = geneAnnotation.all.included,
        testNames = designInfo$testName,
	pThreshold = tt$p,
	pAdjThreshold = tt$pAdj,
        logFCThreshold = tt$logFC,
        FCColPattern = spaste("^", changeCoeffTo, "\\.for"),
	topList.minN = topList.minN,
	geneIDs = getElement(geneAnnotation.start.included, topList.humanColName),
        backupGeneIDs = getElement(geneAnnotation.start.included, topList.backupColName),
        DEDirectionNames = DEDirectionNames,
        addDirectional = designInfo$getDirectionalHits,
        addSplitBy = designInfo$splitSignificantHitsBy,
        addCombined = designInfo$combineSignificantHits,
        testNameSep = testNameSep,
        thresholdName = tName,
        returnFlatList = TRUE),
      topThresholds, topThresholdNames);

     topGeneIDs.individual = mymapply(function(tt, tName)
       getTopGenes(indivResults,
        geneAnnot = geneAnnotation.all.included,
        testNames = designInfo$testName,
        thresholdName = tName,
        pThreshold = tt$p,
        pAdjThreshold = tt$pAdj,
        logFCThreshold = tt$logFC,
        FCColPattern = spaste("^", changeCoeffTo, "\\.for"),
        topList.minN = topList.minN,
        geneIDs = getElement(geneAnnotation.start.included, topList.machineColName),
        DEDirectionNames = DEDirectionNames,
        addDirectional = designInfo$getDirectionalHits,
        addSplitBy = designInfo$splitSignificantHitsBy,
        addCombined = designInfo$combineSignificantHits,
        testNameSep = testNameSep, returnFlatList = TRUE),
      topThresholds, topThresholdNames);

     topSymbols.testNames = sapply(topSymbols.individual[[1]], function(x) attr(x, "testInfo")$testName);

     intermediateResults$topSymbols.individual = topSymbols.individual;
     intermediateResults$topGeneIDs.individual = topGeneIDs.individual;

     # Combine top genes according to combineResults: these will have pretty names.

     combineTestIndex = lapply(combineResults, function(x) which(topSymbols.testNames %in% x));

     topOrNSymbols.combined = lapply(combineTestIndex, function(i)
       interleaveLists.withoutNames(lapply(topSymbols.individual, function(ts) ts[i])));

     if (!is.null(geneSetFnc))
     {
       # Get the geneID-filled extra lists
       topIDs.extraLists = lapply(combinedDESeqResults.ext, function(res)
       {
          out = do.call(geneSetFnc, c(list(results = res, topList.minN = 0, 
                                topThresholds = topThresholds, 
                                geneIDs = getElement(geneAnnotation.start.included, topList.machineColName),
                                designInfo = designInfo, 
                                includedFeatureIndex = includedFeatureIndex),
                                geneSetFncArgs));
          out;
       });

       topOrNIDs.extraLists = lapply(combinedDESeqResults.ext, function(res)
       {
          out = do.call(geneSetFnc, c(list(results = res, topList.minN = topList.minN,
                                topThresholds = topThresholds, 
                                geneIDs = getElement(geneAnnotation.start.included, topList.machineColName),
                                designInfo = designInfo, 
                                includedFeatureIndex = includedFeatureIndex),
                                geneSetFncArgs));
          out;
       });

       # Get the geneID-filled extra lists
       topSymbols.extraLists = lapply(combinedDESeqResults.ext, function(res)
       {
          out = do.call(geneSetFnc, c(list(results = res, topList.minN = 0,
                                topThresholds = topThresholds,
                                geneIDs = getElement(geneAnnotation.start.included, topList.humanColName),
                                backupGeneIDs = getElement(geneAnnotation.start.included, topList.backupColName),
                                designInfo = designInfo, 
                                includedFeatureIndex = includedFeatureIndex),
                                geneSetFncArgs));
          out;
       });

       topOrNSymbols.extraLists = lapply(combinedDESeqResults.ext, function(res)
       {
          out = do.call(geneSetFnc, c(list(results = res, topList.minN = topList.minN,
                                topThresholds = topThresholds,
                                geneIDs = getElement(geneAnnotation.start.included, topList.humanColName),
                                backupGeneIDs = getElement(geneAnnotation.start.included, topList.backupColName),
                                designInfo = designInfo, 
                                includedFeatureIndex = includedFeatureIndex),
                                geneSetFncArgs));
          out;
       });


       # Save both in intermediate results...
       intermediateResults$topSymbols.extraLists = topSymbols.extraLists;
       intermediateResults$topIDs.extraLists = topIDs.extraLists;
       intermediateResults$topOrNSymbols.extraLists = topOrNSymbols.extraLists;
       intermediateResults$topOrNIDs.extraLists = topOrNIDs.extraLists;
       topOrNSymbols.combined = mymapply(c, topOrNSymbols.combined, topOrNSymbols.extraLists);
     }

     intermediateResults$topOrNSymbols.combined = topOrNSymbols.combined;
     forceCalculation = TRUE;
   } else {
     topOrNSymbols.combined = intermediateResults$topOrNSymbols.combined;
     topSymbols.extraLists = intermediateResults$topSymbols.extraLists;
     topIDs.extraLists = intermediateResults$topIDs.extraLists;
     topOrNSymbols.extraLists = intermediateResults$topOrNSymbols.extraLists;
     topOrNIDs.extraLists = intermediateResults$topOrNIDs.extraLists;
     topSymbols.individual = intermediateResults$topSymbols.individual
     topGeneIDs.individual = intermediateResults$topGeneIDs.individual
   }

   # Add test names to "Down"/"Upregulated"

   if (exportIndividualResults && (forceCalculation || forceExport))
     for (th in 1:nTopThresholds) for (d in 1:nDesigns)
     {
       f = file.path(associationDir,
              spaste(topGeneFileBase, prependDash(topThresholdNames[th]), prependDash(testNames[d]), ".csv.gz"));
       write.csv(prettifyNames.loc(irregularList2DataFrame(
                       topSymbols.individual[[th]] [topSymbols.testNames==testNames[g]])),
                 file = gzfile(f), quote = TRUE, row.names = FALSE);
       if (addToAnalysisManager) am.addNewObjectToAnalysis(
          name = c("selected genes for individual tests", thresholdNames.pretty[th], designInfo$testName.pretty[d]),
          fileName = f,
          fileReadArgs = c(defaultReadArguments("csv"), list(check.names = FALSE)),
          analysisID = am.analysisID,
          anaType = amComponent,
          allowNonstandardAnalysisType = TRUE,
          allowNonstandardObjectType = TRUE,
          returnBoth = FALSE, updateAnalysisManager = TRUE);
     }

   if (exportCombinedResults && (forceCalculation || forceExport))
     for (cr in 1:nCombineResults)
     {
       f = file.path(associationDir, spaste(topGeneFileBase, prependDash(combineNames[cr]), ".csv.gz"));
       write.csv(prettifyNames.loc(irregularList2DataFrame(topOrNSymbols.combined[[cr]])),
                 file = gzfile(f), quote = TRUE, row.names = FALSE);
       if (addToAnalysisManager) am.addNewObjectToAnalysis(
          name = c("selected genes for combined tests", combineNames[cr]),
          fileName = f,
          fileReadArgs = c(defaultReadArguments("csv"), list(check.names = FALSE)),
          analysisID = am.analysisID,
          anaType = amComponent,
          allowNonstandardAnalysisType = TRUE,
          allowNonstandardObjectType = TRUE,
          returnBoth = FALSE, updateAnalysisManager = TRUE);

     }

   #-------------------------------------------------------------------------------------------------------
   # Collection of DE genes for anRichment
   #-------------------------------------------------------------------------------------------------------
   if (createCollectionOfSignificantGenes)
   {
     topGeneIDs.collection = mymapply(function(tt, tName)
       suppressWarnings(getTopGenes(indivResults,
        geneAnnot = geneAnnotation.all.included,
        testNames = designInfo$testName,
        thresholdName = tName,
        pThreshold = tt$p,
        pAdjThreshold = tt$pAdj,
        logFCThreshold = tt$logFC,
        FCColPattern = spaste("^", changeCoeffTo, "\\.for"),
        topList.minN = topList.minN,
        geneIDs = getElement(geneAnnotation.start.included, collectionGeneIDcolName),
        DEDirectionNames = DEDirectionNames,
        addDirectional = designInfo$getDirectionalHits,
        addSplitBy = designInfo$splitSignificantHitsBy,
        addCombined = designInfo$combineSignificantHits,
        testNameSep = testNameSep,
        returnFlatList = TRUE)),
      topThresholds, topThresholdNames);

      collectionsOfDEgenes = createCollectionOfDEGenes(
            topGenes = topGeneIDs.collection,
            keepThresholds = collection.keepThresholds,
            thresholdNames.pretty = thresholdNames.pretty,
            keepTests = collection.keepTests,
            designInfo = designInfo,
            geneEvidence = geneEvidence,
            geneSource = geneSource,
            IDBase = IDBase,
            geneSetNamePattern = geneSetNamePattern,
            geneSetDescriptionPattern = geneSetDescriptionPattern,
            organism = organism,
            internalClassification = internalClassification,
            geneSetGroups = geneSetGroups,
            groupList = groupList,
            collectionFile = collectionFile,
            separateCollectionsByThreshold = separateCollectionsByThreshold,
            dropEmptySets = collection.dropEmptySets,
            addToAnalysisManager = addToAnalysisManager,
            analysisID = am.analysisID,
            amComponent = amComponent)
   }

   #-------------------------------------------------------------------------------------------------------
   # Numbers of significant genes
   #-------------------------------------------------------------------------------------------------------

   topIDs.0 = mymapply(function(tt, tName)
     getTopGenes(indivResults,
        geneAnnot = geneAnnotation.all.included,
        testNames = designInfo$testName,
        thresholdName = tName,
        pThreshold = tt$p,
        pAdjThreshold = tt$pAdj,
        logFCThreshold = tt$logFC,
        FCColPattern = spaste("^", changeCoeffTo, "\\.for"),
        topList.minN = 0, geneIDs = geneAnnotation.start.included[, topList.machineColName],
        addDirectional = designInfo$getDirectionalHits,
        addSplitBy = designInfo$splitSignificantHitsBy,
        addCombined = rep(FALSE, nrow(designInfo)),
        DEDirectionNames = DEDirectionNames,
        organizeIntoPairs = TRUE,
        returnFlatList = FALSE),
     topThresholds, topThresholdNames);

   topSymbols.0 = mymapply(function(tt, tName)
     getTopGenes(indivResults,
        geneAnnot = geneAnnotation.all.included,
        testNames = designInfo$testName,
        thresholdName = tName,
        pThreshold = tt$p,
        pAdjThreshold = tt$pAdj,
        logFCThreshold = tt$logFC,
        FCColPattern = spaste("^", changeCoeffTo, "\\.for"),
        topList.minN = 0, geneIDs = geneAnnotation.start.included[, topList.humanColName],
        backupGeneIDs = geneAnnotation.start.included[, topList.backupColName],
        addDirectional = designInfo$getDirectionalHits,
        addSplitBy = designInfo$splitSignificantHitsBy,
        addCombined = rep(FALSE, nrow(designInfo)),
        DEDirectionNames = DEDirectionNames,
        organizeIntoPairs = TRUE,
        returnFlatList = FALSE),
     topThresholds, topThresholdNames);

   nSignif.testNames = sapply(topIDs.0[[1]], function(x) 
      attr(x[[1]], "testInfo")$testName);

   nSignif.all = mymapply(function(te0, thrName)
   {
      out = do.call(rbind, lapply(te0, sapply, length))
      rownames(out) = sapply( lapply(te0, function(x) attr(x[[1]], "testInfo")),
                              modifiedTestName, pretty = FALSE)
      setColnames(out, spaste(DEDirectionNames, thrName));
   }, topIDs.0, topThresholdNames);

   nSignif.all.comb = interleave(nSignif.all, nameBase = rep("", length(nSignif.all)),
                              sep = "");

   if (exportIndividualResults)
   {
     if (exportNSignificant.perThreshold) for (tt in 1:nTopThresholds) 
     {
       f = file.path(associationDir, spaste("numbersOfSignificantHits", prependDash(analysisName),
                               "-allIndividual",
                               prependDash(topThresholdNames[tt]), ".csv"));
       write.csv(prettifyNames(
          prettifyColumns(dataForTable(nSignif.all[[tt]], transpose = FALSE, IDcolName = testColName),
                 column = testColName, prettifyList = prettifyList),
          prettifyList = prettifyList),
          file = f, row.names = FALSE);
       if (addToAnalysisManager) am.addNewObjectToAnalysis(
            name = c("numbers of significant genes", thresholdNames.pretty[tt]),
            fileName = f,
            fileReadArgs = c(defaultReadArguments("csv"), list(check.names = FALSE)),
            analysisID = am.analysisID,
            anaType = amComponent,
            allowNonstandardAnalysisType = TRUE,
            allowNonstandardObjectType = TRUE,
            returnBoth = FALSE, updateAnalysisManager = TRUE);
     }

     if (exportNSignificant.combinedThresholds)
     {
       f = file.path(associationDir, spaste("numbersOfSignificantHits", prependDash(analysisName),
                                 "-allIndividual",
                                 "-allThresholds.csv"));
       write.csv(prettifyNames(
          prettifyColumns(dataForTable(nSignif.all.comb, transpose = FALSE, IDcolName = testColName),
                 column = testColName, prettifyList = prettifyList),
          prettifyList = prettifyList),
          file = f, row.names = FALSE);

       if (addToAnalysisManager) am.addNewObjectToAnalysis(
            name = c("numbers of significant genes", "all thresholds"),
            fileName = f,
            fileReadArgs = c(defaultReadArguments("csv"), list(check.names = FALSE)),
            analysisID = am.analysisID,
            anaType = amComponent,
            allowNonstandardAnalysisType = TRUE,
            allowNonstandardObjectType = TRUE,
            returnBoth = FALSE, updateAnalysisManager = TRUE);
     }
   }
   intermediateResults$nSignif.all = nSignif.all;
   intermediateResults$nSignif.all.combinedThresholds = nSignif.all.comb;

   # Add top genes to the tables
   nAdd = function(n) if (n<=maxmaxGenesInTables) n else maxGenesInTables;
   char.top = function(x) 
   {
     if (length(x) > 0) spaste(" (",
              paste(x[1:nAdd(length(x))], collapse = ", "), ")") else "";
   }

   nSignif.all.ext = mymapply(function(tg0, nSig1) 
   {
     char.top = do.call(rbind, lapply(tg0, sapply, char.top));
     out = nSig1;
     out[, 1] = spaste(nSig1[, 1], char.top[, 1])
     out[, 2] = spaste(nSig1[, 2], char.top[, 2])
     out;
   }, topSymbols.0, nSignif.all)

   nSignif.all.ext.comb = interleave(nSignif.all.ext, nameBase = rep("", length(nSignif.all.ext)),
                              sep = "");

   if (exportIndividualResults)
   {
     if (exportCombinedResults) allName1 = "-0-allIndividual" else allName1 = "";

     if (exportNSignificant.perThreshold) for (tt in 1:nTopThresholds) 
     {
       f = file.path(associationDir, spaste("numbersOfSignificantHits-withTopGenes", prependDash(analysisName), 
                               allName1, prependDash(topThresholdNames[tt]), ".csv"));
       write.csv( prettifyNames(
             prettifyColumns(dataForTable(nSignif.all.ext[[tt]], transpose = FALSE, IDcolName = testColName),
                 column = testColName, prettifyList = prettifyList),
             prettifyList),
         file = f, row.names = FALSE);
       if (addToAnalysisManager) am.addNewObjectToAnalysis(
            name = c("numbers of significant genes with top genes", thresholdNames.pretty[tt]),
            fileName = f,
            fileReadArgs = c(defaultReadArguments("csv"), list(check.names = FALSE)),
            analysisID = am.analysisID,
            anaType = amComponent,
            allowNonstandardAnalysisType = TRUE,
            allowNonstandardObjectType = TRUE,
            returnBoth = FALSE, updateAnalysisManager = TRUE);
     }

     if (exportNSignificant.combinedThresholds)
     {
       f = file.path(associationDir, spaste("numbersOfSignificantHits-withTopGenes", prependDash(analysisName),
                                 allName1, "-allThresholds.csv"));
       write.csv( prettifyNames(
             prettifyColumns(dataForTable(nSignif.all.ext.comb, transpose = FALSE, IDcolName = testColName),
                 column = testColName, prettifyList = prettifyList),
             prettifyList),
         file = f, row.names = FALSE);
       if (addToAnalysisManager) am.addNewObjectToAnalysis(
            name = c("numbers of significant genes with top genes", "all thresholds"),
            fileName = f,
            fileReadArgs = c(defaultReadArguments("csv"), list(check.names = FALSE)),
            analysisID = am.analysisID,
            anaType = amComponent,
            allowNonstandardAnalysisType = TRUE,
            allowNonstandardObjectType = TRUE,
            returnBoth = FALSE, updateAnalysisManager = TRUE);
     }
   }

   combineIndex = lapply(combineResults, function(cr) which(nSignif.testNames %in% cr));
   nSignif.byCombine = lapply(combineIndex, function(cr)
      lapply(nSignif.all, function(ns1) ns1[cr, , drop = FALSE]))

   nSignif.byCombine.comb = lapply(nSignif.byCombine, interleave.withoutNames);

   nSignif.byCombine.ext = lapply(combineResults, function(cr)
      lapply(nSignif.all.ext, function(ns1) ns1[cr, , drop = FALSE]))

   nSignif.byCombine.ext.comb = lapply(nSignif.byCombine.ext, interleave.withoutNames);

   intermediateResults$nSignif.byCombine = nSignif.byCombine;
   intermediateResults$nSignif.byCombine.combinedThresholds = nSignif.byCombine.comb;
 
   if (exportCombinedResults) for (cr in 1:nCombineResults) 
   {
     if (exportNSignificant.perThreshold) for (tt in 1:nTopThresholds) 
     {
       f= file.path(associationDir, spaste("numbersOfSignificantHits", prependDash(analysisName),
                               prependDash(prependZeros(cr, nchar(nCombineResults))),
                               prependDash(combineNames[cr]), prependDash(topThresholdNames[tt]), ".csv"));
       write.csv( prettifyNames(
             prettifyColumns(dataForTable(nSignif.byCombine[[cr]] [[tt]], transpose = FALSE, 
                                          IDcolName = testColName),
                   column = testColName, prettifyList = prettifyList),
             prettifyList = prettifyList),
          file = f, row.names = FALSE);

       if (addToAnalysisManager) am.addNewObjectToAnalysis(
            name = c("numbers of significant genes for combined tests", thresholdNames.pretty[tt],
                     combineNames[cr]),
            fileName = f,
            fileReadArgs = c(defaultReadArguments("csv"), list(check.names = FALSE)),
            analysisID = am.analysisID,
            anaType = amComponent,
            allowNonstandardAnalysisType = TRUE,
            allowNonstandardObjectType = TRUE,
            returnBoth = FALSE, updateAnalysisManager = TRUE);

       f = file.path(associationDir, spaste("numbersOfSignificantHits-withTopGenes", prependDash(analysisName),
                               prependDash(prependZeros(cr, nchar(nCombineResults))),
                               prependDash(combineNames[cr]), prependDash(topThresholdNames[tt]), ".csv"));
       write.csv(
          prettifyNames(
             prettifyColumns(dataForTable(nSignif.byCombine.ext[[cr]] [[tt]], transpose = FALSE, 
                                          IDcolName = testColName),
                   column = testColName, prettifyList = prettifyList),
             prettifyList = prettifyList),
          file = f, row.names = FALSE);
       if (addToAnalysisManager) am.addNewObjectToAnalysis(
            name = c("numbers of significant genes with top genes for combined tests", thresholdNames.pretty[tt],
                     combineNames[cr]),
            fileName = f,
            analysisID = am.analysisID,
            anaType = amComponent,
            allowNonstandardAnalysisType = TRUE,
            allowNonstandardObjectType = TRUE,
            fileReadArgs = c(defaultReadArguments("csv"), list(check.names = FALSE)),
            returnBoth = FALSE, updateAnalysisManager = TRUE);
     }

     if (exportNSignificant.combinedThresholds)
     {
       f= file.path(associationDir, spaste("numbersOfSignificantHits", prependDash(analysisName),
                               prependDash(prependZeros(cr, nchar(nCombineResults))),
                               prependDash(combineNames[cr]), "-allThreshold.csv"));
       write.csv( prettifyNames(
             prettifyColumns(dataForTable(nSignif.byCombine.comb[[cr]], transpose = FALSE,
                                          IDcolName = testColName),
                   column = testColName, prettifyList = prettifyList),
             prettifyList = prettifyList),
          file = f, row.names = FALSE);

       if (addToAnalysisManager) am.addNewObjectToAnalysis(
            name = c("numbers of significant genes for combined tests", "all thresholds",
                     combineNames[cr]),
            fileName = f,
            fileReadArgs = c(defaultReadArguments("csv"), list(check.names = FALSE)),
            analysisID = am.analysisID,
            anaType = amComponent,
            allowNonstandardAnalysisType = TRUE,
            allowNonstandardObjectType = TRUE,
            returnBoth = FALSE, updateAnalysisManager = TRUE);

       f = file.path(associationDir, spaste("numbersOfSignificantHits-withTopGenes", prependDash(analysisName),
                               prependDash(prependZeros(cr, nchar(nCombineResults))),
                               prependDash(combineNames[cr]), "-allThresholds.csv"));
       write.csv(
          prettifyNames(
             prettifyColumns(dataForTable(nSignif.byCombine.ext.comb[[cr]], transpose = FALSE,
                                          IDcolName = testColName),
                   column = testColName, prettifyList = prettifyList),
             prettifyList = prettifyList),
          file = f, row.names = FALSE);
       if (addToAnalysisManager) am.addNewObjectToAnalysis(
            name = c("numbers of significant genes with top genes for combined tests", "all thresholds",
                     combineNames[cr]),
            fileName = f,
            analysisID = am.analysisID,
            anaType = amComponent,
            allowNonstandardAnalysisType = TRUE,
            allowNonstandardObjectType = TRUE,
            fileReadArgs = c(defaultReadArguments("csv"), list(check.names = FALSE)),
            returnBoth = FALSE, updateAnalysisManager = TRUE);
     }
   }

   if (stopAt=="40-Gene lists") 
   {
     return(intermediateResults);
   }

   #-------------------------------------------------------------------------------------------------------
   # Enrichment analysis
   #-------------------------------------------------------------------------------------------------------

   if ( length(enrichment.collections) > 0 && nDataSets(enrichment.combinedCollection)>0 && 
        length(enrichment.useCRIndex) > 0 && length(enrichment.useThresholds) > 0)
   {
     if (forceCalculation || is.null(intermediateResults$enrichment))
     {
       nCollections = length(enrichment.collections) + 1;
       collectionNames = c(enrichment.combinedCollectionName, names(enrichment.collections));

       if (!enrichment.annotIDCol %in% c(names(geneAnnotation.start), names(geneAnnotation.end),
               names(geneAnnotation.extra)))
         stop("'enrichment.annotIDCol' must specify a valid column name in 'geneAnnotation.start', 'geneAnnotation.end' ",
              "or 'geneAnnotation.extra'.");
       if (enrichment.annotIDCol!="Entez") printFlush("Using enrichment annotation ID column", enrichment.annotIDCol);
       geneIDsForEnrichment = geneAnnotation.start.included[[enrichment.annotIDCol]];
       if (is.null(geneIDsForEnrichment)) geneIDsForEnrichment = geneAnnotation.end.included[[enrichment.annotIDCol]];
       if (is.null(geneIDsForEnrichment)) geneIDsForEnrichment = geneAnnotation.extra.included[[enrichment.annotIDCol]];
       if (is.null(geneIDsForEnrichment)) stop("Something's wrong, maybe a NULL component?");
       
       topEntrez.enrichment = mymapply(function(tt, tName)
          suppressWarnings(getTopGenes(indivResults,
            geneAnnot = geneAnnotation.all.included,
            testNames = designInfo$testName,
            thresholdName = tName,
            pThreshold = tt$p,
            pAdjThreshold = tt$pAdj,
            logFCThreshold = tt$logFC,
            FCColPattern = spaste("^", changeCoeffTo, "\\.for"),
            topList.minN = enrichment.minN, geneIDs = geneIDsForEnrichment,
            DEDirectionNames = DEDirectionNames,
            addDirectional = designInfo$getDirectionalHits,
            addSplitBy = designInfo$splitSignificantHitsBy,
            addCombined = designInfo$combineSignificantHits,
            testNameSep = testNameSep, returnFlatList = TRUE)), 
        topThresholds[enrichment.useThresholds], enrichment.useThresholds);

       topEntrez.testNames = as.character(sapply(topEntrez.enrichment[[1]], function(x)
        attr(x, "testInfo")$testName));

       combineResults.enr = combineResults;

       if (length(enrichment.dropTests) > 0) 
         combineResults.enr = lapply(combineResults, setdiff, enrichment.dropTests)

       combineTestIndex.enr = lapply(combineResults.enr, function(x) which(topEntrez.testNames %in% x));

       topEntrez.comb.enrichment =  lapply(combineTestIndex.enr, function(i)
          lapply(topEntrez.enrichment, function(ts) ts[i]));
       # The result here is a 3-level list, top level is combineIndex, 2nd level threshold, 3rd level is test 

       if (!is.null(geneSetFnc)) {
         # Get the Entrez-filled extra lists
         topOrNIDs.extraLists.enrichment = lapply(combinedDESeqResults, function(res)
         {
            out = lapply(enrichment.useThresholds, function(tt)
                    do.call(geneSetFnc, c(list(results = res, topList.minN = enrichment.minN,
                                        topThresholds = topThresholds,
                                        topThreshold1 = tt,
                                        geneIDs = geneIDsForEnrichment),
                                  geneSetFncArgs)));
                  ## Note to self: note the argument topThreshold1, so named to remind that it should carry a single
                  ## threshold name to be selected from topThresholds
            out;
         });

         topEntrez.ext.enrichment = mymapply(function(x, y) mymapply(c, x, y), 
                                              topEntrez.comb.enrichment, topOrNIDs.extraLists.enrichment);
       } else 
          topEntrez.ext.enrichment = topEntrez.comb.enrichment;

      # Drop length-zero gene lists
       topEntrez.enr = lapply(topEntrez.ext.enrichment, lapply, function(lst) lst[sapply(lst, length)>0]);

       enrichment.otherArgs$verbose = verbose-1;
       enrichment.otherArgs$indent = indent + 1;
       enrichment = list();
       nUseCombineResults = length(enrichment.useCRIndex);
       for (cri in 1:nUseCombineResults)
       {
         if (verbose > 1) printFlush(spaste("  Enrichment calculation on combine index ", cri));
         cr = enrichment.useCRIndex[cri];
         enrichment[[cri]] = list();
         for (tt in 1:nEnrichmentThresholds)
         {
           if (length(topEntrez.enr[[cri]][[tt]]) > 0)
           {
             if (verbose > 1) printFlush(spaste("  Enrichment calculation on significance treshold ", tt));
             enr1 = do.call(enrichmentAnalysis.Entrez,
                         c(list(active = topEntrez.enr[[cr]][[tt]], inactive = geneIDsForEnrichment,
                                 refCollection = enrichment.combinedCollection,
                                 combineEnrichmentTables = TRUE),
                           enrichment.otherArgs));
             enrTable = enr1$enrichmentTable[, !colnames(enr1$enrichmentTable) %in% enrichment.dropCols];
             enrichment[[cri]][[tt]] = list(CombinedCollection = c(list( enrichmentTable = enrTable),
                      enr1[if (minimizeEnrichmentSize) numeric(0) else c("countsInDataSet", "pValues")]));
           } else 
             enrichment[[cri]][[tt]] = list(CombinedCollection = list());
         }
         names(enrichment[[cri]]) = enrichment.useThresholds;
       }
       names(enrichment) = combineNames[enrichment.useCRIndex];
       intermediateResults$enrichment = enrichment;
       intermediateResults$topEntrez.ext.enrichment = topEntrez.ext.enrichment;
       intermediateResults$topEntrezForEnrichment = topEntrez.enr;
       intermediateResults$enrichment.useCRIndex = enrichment.useCRIndex;
       intermediateResults$enrichment.useResultsCombinations = enrichment.useResultsCombinations;
       forceCalculation = TRUE;
     } else {
       enrichment = intermediateResults$enrichment;
       #if (length(enrichment)!= length(enrichment.useCRIndex)) 
       #  stop("Length of 'enrichment.useResultCombinations' is not consistent with length of 'enrichment'\n",
       #       "  in 'intermediateResults'.")
       topEntrez.ext.enrichment = intermediateResults$topEntrez.ext.enrichment;
       nCollections = length(enrichment.collections) + 1;
       collectionNames = c(enrichment.combinedCollectionName, names(enrichment.collections));
       enrichment.useResultsCombinations = intermediateResults$enrichment.useResultsCombinations;
       enrichment.useCRIndex = intermediateResults$enrichment.useCRIndex;
     }
     if (length(intermediateResults$enrichment)>0 && exportEnrichmentTables && (forceCalculation || forceExport))
     {
       nUseCombineResults = length(enrichment.useCRIndex);
       for (cri in 1:nUseCombineResults) for (tt in 1:nEnrichmentThresholds) 
          if (length(enrichment[[cri]][[tt]]$CombinedCollection) > 0)
          {
            enrDir1 = multiGSub(c("%combine", "%threshold"), 
                               c(enrichment.useResultsCombinations[cri], 
                                 gsub("_", "", enrichment.useThresholds[tt])), enrichmentDir);
            dir.create(enrDir1, recursive = TRUE, showWarnings = FALSE);
            cr = enrichment.useCRIndex[cri];
            enr.subColl = splitEnrichmentTableBySubcollections(enrichment[[cri]][[tt]]$CombinedCollection$enrichmentTable, 
                            enrichment.combinedCollection, enrichment.collections, dropColumns = enrichment.dropCols);
            enrichment.temp = c(list(enrichment[[cri]][[tt]]$CombinedCollection$enrichmentTable), enr.subColl);
            for (col in 1:nCollections) if (nrow(enrichment.temp[[col]]) > 0)
            {
              f = file.path(enrDir1, spaste(enrichmentFileBase, "-", combineNames[cr],
                                       "-", enrichment.useThresholds[tt], "-", prependZeros(col-1, 2), "-", 
                                       collectionNames[col], ".csv"));
              if (enrichment.compress) f = spaste(f, ".gz");
              write.csv(prettifyColumns(
                          signifNumeric(enrichment.temp[[col]], exportDigits),
                          columns = "class", prettifyList = prettifyList),
                    file = if (enrichment.compress) gzfile(f) else f, quote = TRUE, row.names = FALSE);
              if (addToAnalysisManager) am.addNewObjectToAnalysis(
                 name = c("enrichment analysis", combineNames[cr], enrichment.useThresholds[tt], collectionNames[col]),
                 fileName = f,
                 analysisID = am.analysisID,
                 anaType = amComponent,
                 allowNonstandardAnalysisType = TRUE,
                 allowNonstandardObjectType = TRUE,
                 returnBoth = FALSE, updateAnalysisManager = TRUE);
            }
          }
     }
   }

    
   #if (stopAt=="50-Enrichment") 
   #{
   #  return(intermediateResults);
   #}

   return(intermediateResults);   
}
