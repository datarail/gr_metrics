.GRcalculate = function(inputData, groupingVariables, cap = FALSE, case = "A"){
  # declaring values NULL to avoid note on package check
  cell_count = NULL
  cell_count__time0 = NULL
  cell_count__ctrl = NULL
  log2nn = with(inputData, log2(cell_count/cell_count__time0))
  log2nn_ctrl = with(inputData, log2(cell_count__ctrl/cell_count__time0))
  GR = 2^(log2nn/log2nn_ctrl) - 1
  input_edited = inputData
  input_edited$log10_concentration = log10(input_edited$concentration)
  input_edited$GR = GR
  tmp<-input_edited[,groupingVariables, drop = FALSE]
  experimentNew = (apply(tmp,1, function(x) (paste(x,collapse=" "))))
  if(cap == TRUE) {
    input_edited$GR[input_edited$GR > 1] = 1
    input_edited$GR[input_edited$GR < -1] = -1
  }
  if(length(groupingVariables) > 0) {
    input_edited$experiment = as.factor(experimentNew)
  } else {
    input_edited$experiment = as.factor("All Data")
  }
  return(input_edited)
}

.GRlogisticFit = function(inputData, groupingVariables, force = FALSE,
                          cap = FALSE) {
  # declaring values NULL to avoid note on package check
  experiment = NULL
  GEC50 = NULL
  GRinf = NULL
  Hill = NULL
  experiments = levels(inputData$experiment)
  parameters = matrix(data = NA, ncol = 3, nrow = length(experiments))
  parameters = as.data.frame(parameters)
  colnames(parameters) = c('Hill','GRinf','GEC50')
  if(length(groupingVariables) > 0) {
    metadata = matrix(data = NA, ncol = length(groupingVariables),
                      nrow = length(experiments))
    metadata = as.data.frame(metadata)
    colnames(metadata) = groupingVariables
  } else {
    metadata = NULL
  }
  pval = NULL
  GRmax = NULL
  GR_mean = NULL
  AOC = NULL
  R_square = NULL
  for(i in 1:length(experiments)) {
    # print(i)
    data_exp = inputData[inputData$experiment == experiments[i], ]
    concs = sort(unique(data_exp$concentration))
    l = length(concs)
    max_concs = data_exp[data_exp$concentration %in% concs[c(l,l-1)],]
    GRmax[i] = min(max_concs$GR, na.rm = TRUE)
    #     metadata[i,] = data_exp[1,1:5]
    if(!is.null(metadata)) {
      metadata[i,] = data_exp[1,groupingVariables, drop = FALSE]
    }
    GR_mean[i] = mean(data_exp$GR, na.rm = TRUE)
    #===== constrained fit ============
    c = unique(data_exp$concentration)
    priors = c(2, 0.1, stats::median(c))
    lower = c(.1, -1, min(c)*1e-2)
    upper = c(5, 1, max(c)*1e2)
    if(dim(data_exp)[1] > 1) {
      controls = drc::drmc()
      controls$relTol = 1e-06
      controls$errorm = FALSE
      controls$noMessage = TRUE
      controls$rmNA = TRUE
      output_model_new = try(drc::drm(
        GR~concentration, experiment, data=data_exp,
        fct=drc::LL.3u(names=c('Hill','GRinf','GEC50')), start = priors,
        lowerl = lower, upperl = upper, control = controls,
        na.action = na.omit))
      if(class(output_model_new)!="try-error") {
        parameters[i,] = c(as.numeric(stats::coef(output_model_new)))
        # F-test for the significance of the sigmoidal fit
        Npara = 3 # N of parameters in the growth curve
        Npara_flat = 1 # F-test for the models
        RSS2 = sum(stats::residuals(output_model_new)^2, na.rm = TRUE)
        RSS1 = sum((data_exp$GR - mean(data_exp$GR, na.rm = TRUE))^2,
                   na.rm = TRUE)
        df1 = (Npara - Npara_flat)
        df2 = (length(na.omit(data_exp$GR)) - Npara + 1)
        f_value = ((RSS1-RSS2)/df1)/(RSS2/df2)
        f_pval = stats::pf(f_value, df1, df2, lower.tail = FALSE)
        pval[i] = f_pval
        R_square[i] = 1 - RSS2/RSS1
      }
    }
    #==================================

    # Trapezoid rule for integration of GR_AOC
    GRavg = NULL
    for(j in 1:length(concs)) {
      data_trapz = data_exp[data_exp$concentration == concs[j],]
      GRavg[j] = mean(data_trapz$GR, na.rm = TRUE)
    }
    AOC[i] = sum((1 - (GRavg[1:(length(GRavg)-1)]+GRavg[2:length(GRavg)])/2)*
                   diff(log10(concs), lag = 1), na.rm = TRUE)/
      (log10(concs[length(concs)]) - log10(concs[1]))
  }

  # Calculate GR50 from parameters
  parameters$GR50 = with(parameters,GEC50*((1-GRinf)/(0.5-GRinf) - 1)^(1/Hill))
  parameters$GRmax = GRmax
  parameters$GR_AOC = AOC
  parameters$r2 = R_square
  if(is.null(pval)) {pval = NA}
  parameters$pval = pval
  # Re-order rows to match reference_output
  parameters$experiment = experiments
  # Threshold for F-test pval
  pcutoff = ifelse(force == FALSE, .05 , 1)
  for(i in 1:dim(parameters)[1]) {
    if(!is.na(parameters$pval[i])) {
      parameters$fit[i] = ifelse(parameters$pval[i] >= pcutoff |
                                 is.na(parameters$GEC50[i]), "flat","sigmoid")
    } else {
      parameters$fit[i] = ifelse(is.na(parameters$GEC50[i]), "flat", "sigmoid")
    }
  }
  # changed to above code to deal with NAs
  #parameters$fit = with(parameters, ifelse(pval >= pcutoff | is.na(GEC50),
  #"flat","sigmoid"))
  # Add values for flat fits: GEC50 = 0, Hill = 0.01 and GR50 = +/- Inf
  parameters$flat_fit = GR_mean
  for(i in 1:dim(parameters)[1]) {
    if(parameters$fit[i] == "flat") {
      parameters$GEC50[i] = 0
      parameters$Hill[i] = 0.01
      parameters$GR50[i] = ifelse(parameters$flat_fit[i] > .5, Inf, -Inf)
      parameters$GRinf[i] = parameters$flat_fit[i]
    }
  }
  # Add GR50 = +/-Inf for any curves that don't reach GR = 0.5
  for(i in 1:dim(parameters)[1]) {
    if(is.na(parameters$GR50[i])) {
      parameters$GR50[i] = ifelse(parameters$flat_fit[i] > .5, Inf, -Inf)
    }
  }
  for(i in 1:dim(parameters)[1]) {
    if(parameters$fit[i] == "sigmoid") {
      parameters$flat_fit[i] = NA
    }
  }
  parameters = parameters[,c('GR50','GRmax','GR_AOC','GEC50','GRinf','Hill',
                             'r2','pval','experiment', 'fit','flat_fit')]
  if(!is.null(metadata)) {
    parameters = cbind(metadata, parameters)
  }
  return(parameters)
}

.GRlogistic_3u = function(c, GRinf, GEC50, Hill){
  GRinf + (1 - GRinf)/(1 + (c/GEC50)^Hill)
  }

.trim_mean = function(x, percent) {
  x = x[!is.na(x)]
  n = length(x)
  k = n*(percent/100)/2
  # round down if k is half an integer
  if(round(k) != k & round(k*2) == k*2) {
    lo = floor(k) + 1
    hi = n - lo + 1
  } else {
    lo = round(k) + 1
    hi = n - lo + 1
  }
  x = sort(x)[lo:hi]
  return(mean(x))
}

.convert = function(inputData, case) {
  if(case == "A") {
    if(length(intersect(colnames(inputData), c('concentration', 'cell_count',
                                               'cell_count__ctrl',
                                               'cell_count__time0'))) == 4) {
      return(inputData)
    } else {
      stop("There must be columns named 'concentration', 'cell_count',
            'cell_count__ctrl', and 'cell_count__time0' in inputData")
    }
  } else if(case == "C") {
    if(length(intersect(colnames(inputData), c('concentration', 'cell_count',
                                               'time'))) != 3) {
      stop("There must be columns named 'concentration', 'cell_count',
           and 'time' in inputData")
    }
    delete_cols = which(colnames(inputData) %in% c('concentration',
                                                   'cell_count'))
    keys = colnames(inputData)[-delete_cols]
    time0 = inputData[inputData$time == 0, c(keys, 'cell_count')]
    ctrl = inputData[inputData$concentration == 0 & inputData$time > 0,
                     c(keys, 'cell_count')]
    data = inputData[inputData$concentration != 0 & inputData$time > 0, ]
    time0_keys = NULL
    ctrl_keys = NULL
    for(i in 1:length(keys)) {
      time0_keys[i] = length(intersect(time0[[ keys[i] ]],
                                       data[[ keys[i] ]])) > 0
      ctrl_keys[i] = length(intersect(ctrl[[ keys[i] ]],
                                      data[[ keys[i] ]])) > 0
    }
    ctrl_keys = keys[ctrl_keys]
    time0_keys = keys[time0_keys]

    temp = ctrl[, ctrl_keys]
    ctrl$key = apply(temp, 1, function(x) paste(x, collapse = ' '))

    temp = time0[, time0_keys]
    time0$key = apply(temp, 1, function(x) paste(x, collapse = ' '))

    temp = data[, ctrl_keys]
    data$key_ctrl = apply(temp, 1, function(x) paste(x, collapse = ' '))

    temp = data[, time0_keys]
    data$key_time0 = apply(temp, 1, function(x) paste(x, collapse = ' '))

    data$cell_count__ctrl = NA
    data$cell_count__time0 = NA
    for(key in unique(ctrl$key)) {
      trimmed_mean = .trim_mean(ctrl[ctrl$key == key,]$cell_count, 50)
      data[data$key_ctrl == key, 'cell_count__ctrl'] = trimmed_mean
    }

    for(key in unique(time0$key)) {
      trimmed_mean = .trim_mean(time0[time0$key == key,]$cell_count, 50)
      data[data$key_time0 == key, 'cell_count__time0'] = trimmed_mean
    }

    delete_cols = which(colnames(data) %in% c('key_ctrl', 'key_time0'))
    data = data[, -delete_cols]

    row.names(data) = 1:dim(data)[1]
    inputData = data
    return(inputData)
  }
}

#' Extract GR parameters from a dataset
#'
#' This function takes in a dataset with information about concentration,
#' cell counts over time, and additional grouping variables for a dose-response
#' assay and calculates growth-rate inhibition (GR) metrics for each experiment
#' in the dataset. The data must be in a specific format: either that specified
#' by case "A" or case "C" described in the details below.
#'
#' @param inputData a data table in one of the specified formats (Case A or
#' Case C). See details below for description. See \code{data(inputCaseA)} or
#' \code{data(inputCaseC)} for example input data frames. See help files for
#' \code{\link{inputCaseA}} and \code{\link{inputCaseC}} for description of
#' these examples.
#' @param groupingVariables a vector of column names from inputData. All of the
#' columns in inputData except for those identified here will be averaged over.
#' @param case either "A" or "C", indicating the format of the input data. See
#' below for descriptions of these formats.
#' @param force a logical value indicating whether to attempt to "force" a
#' sigmoidal fit, i.e. whether to allow fits with F-test p-values greater
#' than .05
#' @param cap a logical value indicating whether to cap GR values in between
#'  -1 and 1 (the range of the GR dose-response curve). If true, all GR values
#'  greater than 1 will be set to 1 and all values less than -1 will be set
#'  to -1.
#' @return A SummarizedExperiment object containing GR metrics parameters
#' (GR50, GRmax, etc.) as well as goodness of fit measures is returned. The
#' object also contains, in its metadata, a table of the original data
#' converted to the style of "Case A" (with calculated GR values for each row)
#' and a vector of the grouping variables used for the calculation.
#' @author Nicholas Clark
#' @details
#' Calculation of GR values is performed by the function \code{.GRcalculate}
#' according to the "Online Methods" section of Hafner and Niepel et al.
#' (\url{http://dx.doi.org/10.1038/nmeth.3853}).
#'
#' The fitting of the logistic curve is performed by the \code{.GRlogisticFit}
#' function, which calls the \code{drm} function from the \code{drc} package
#' to solve for the curve parameters. The GR curve fit function is
#' given by f(c) = GRinf + (1 - GRinf)/(1 + (c/GEC50)^Hill) where c is
#' concentration. The fit is performed under following constraints: Hill slope
#' in [.1, 5], GRinf in [-1, 1], and GEC50 in [min(c)*1e-2, max(c)*1e2] (c is
#' concentration). The initial conditions for the fitting algorithm are Hill
#' slope = 2, GRinf = 0.1 and GEC50 = median(c).
#'
#' The parameters of the curve for each experiment are fitted separately. An
#' F-test is used to compare the sigmoidal fit to a flat line fit. If the
#' p-value of the F-test is less than .05, the sigmoidal fit is accepted. If
#' the p-value is greater than or equal to .05, a flat horizontal line fit is
#' given, with y equal to the mean of the GR values. In the parameters data
#' table, for each flat fit, GEC50 is set to 0, Hill is set to 0.01, GRinf is
#' set to the y value of the flat fit (the mean of the GR values), and GR50 is
#' set to +/-Inf depending on whether GRinf is greater or less than .5.
#'
#' The mandatory columns for inputData for Case "A" are the following as
#' well as other grouping columns.
#'
#' 1. concentration - column with concentration values (not log transformed)
#' of the perturbagen on which dose-response curves will be evaluated
#'
#' 2. cell_count - column with the measure of cell number (or a surrogate of
#' cell number) after treatment
#'
#' 3. cell_count__time0 - column with initial (Time 0) cell counts - the
#' measure of cell number in untreated wells grown in parallel until the
#' time of treatment
#'
#' 4. cell_count__ctrl - column with the Control cell count: the measure of
#' cell number in control (e.g. untreated or DMSO-treated) wells from the
#' same plate
#'
#' All other columns will be treated as additional keys on which the data
#' will be grouped (e.g. cell_line, drug, time, replicate)
#'
#' The mandatory columns for inputData for Case "C" are the following as
#' well as other grouping columns.
#'
#' 1. concentration - column with concentration values (not log transformed)
#' of the perturbagen on which dose-response curves will be evaluated
#'
#' 2. cell_count - column with the measure of cell number (or a surrogate of
#' cell number)
#'
#' 3. time - column with the time at which a cell count is observed
#'
#' All other columns will be treated as additional keys on which the data
#' will be grouped (e.g. *cell_line*, *drug*, *replicate*)
#' @note
#' To see the underlying code, use (\code{getAnywhere(.GRlogistic_3u)}),
#' (\code{getAnywhere(.GRcalculate)}), and (\code{getAnywhere(.GRlogisticFit)})
#' @seealso See \code{\link{drm}} for the general logistic fit function that
#' solves for the parameters GRinf, GEC50, and Hill slope. See
#' \code{\link{drmc}} for
#' options of this function. Use the functions \code{\link{GRdrawDRC}},
#' \code{\link{GRbox}}, and \code{\link{GRscatter}} to create visualizations
#' using the output from this function. For online GR calculator and browser,
#' see \url{http://www.grcalculator.org}.
#' @references Hafner, M., Niepel, M. Chung, M. and Sorger, P.K.,
#' "Growth Rate Inhibition Metrics Correct For Confounders In Measuring
#' Sensitivity To Cancer Drugs". \emph{Nature Methods} 13.6 (2016): 521-527.
#' \url{http://dx.doi.org/10.1038/nmeth.3853}
#'
#' @examples
#' # Load Case A (example 1) input
#' data("inputCaseA")
#' head(inputCaseA)
#' # Run GRfit function with case = "A"
#' output1 = GRfit(inputData = inputCaseA,
#' groupingVariables = c('cell_line','agent', 'perturbation','replicate',
#' 'time'))
#' # Overview of SummarizedExperiment output data
#' output1
#' # View GR metrics table
#' assay(output1)
#' # View details of each experiment
#' colData(output1)
#' # View descriptions of each GR metric (or goodness of fit measure)
#' \dontrun{
#' View(rowData(output1))
#' }
#' # View table of original data (converted to style of Case A) with GR values
#' metadata(output1)[[1]]
#' # View vector of grouping variables used for calculation
#' metadata(output1)[[2]]
#' # Load Case C (example 4) input
#' # Same data, different format
#' data("inputCaseC")
#' head(inputCaseC)
#' output4 = GRfit(inputData = inputCaseC,
#' groupingVariables = c('cell_line','agent', 'perturbation','replicate',
#' 'time'),
#' case = "C")
#' # Extract data tables and export to .tsv or .csv
#' \dontrun{
#' # Write GR metrics parameter table to tab-separated text file
#' write.table(assay(output1), file = "filename.tsv", quote = FALSE,
#' sep = "\t", row.names = FALSE)
#' # Write original data plus GR values to comma-separated file
#' write.table(metadata(output1)[[1]], file = "filename.csv", quote = FALSE,
#' sep = ",", row.names = FALSE)
#' }
#' @export

GRfit = function(inputData, groupingVariables, case = "A",
                 force = FALSE, cap = FALSE) {
  if('experiment' %in% colnames(inputData)) {
    stop("Change name of 'experiment' column.")
  }
  inputData = .convert(inputData, case)
  gr_table = .GRcalculate(inputData, groupingVariables, cap, case)
  parameter_table = .GRlogisticFit(gr_table, groupingVariables, force, cap)

  colData = parameter_table[ ,c(groupingVariables, 'fit', 'experiment')]
  rownames(colData) = colData$experiment
  colData = S4Vectors::DataFrame(colData)

  Metrics = c("GR50","GRmax","GR_AOC","GEC50","GRinf","Hill",
                "r2","pval","flat_fit")
  assays = parameter_table[ , Metrics]
  rownames(assays) = parameter_table$experiment
  assays = t(assays)

  Description = c("The concentration at which GR(c) = 0.5",
                  "The maximal effect of the drug (minimal GR value",
                  "The 'Area Over the Curve' - The area between the line
                  GR = 1 and the curve, similar to traditional AUC",
                  "The concentration at half-maximal effect",
                  "The Hill coefficient of the fitted curve, which reflects
                  how steep the dose-response curve is",
                  "The coefficient of determination - essentially how well
                  the curve fits to the data points",
                  "The p-value of the F-test comparing the fit of the curve
                  to a horizontal line fit",
                  "For data that doesn't significantly fit better to a curve
                  than a horizontal line fit, the GR value of the flat line")
  rowData = cbind(Metrics, Description)
  rownames(rowData) = Metrics
  rowData = S4Vectors::DataFrame(rowData)
  rowData$Metrics = as.character(rowData$Metrics)
  rowData$Description = as.character(rowData$Description)

  output = SummarizedExperiment::SummarizedExperiment(assays = assays,
                                                      colData = colData,
            rowData = rowData, metadata = list(gr_table, groupingVariables))
  return(output)
}

