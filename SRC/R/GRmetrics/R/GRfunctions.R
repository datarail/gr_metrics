.GRcalculate = function(inputData, groupingVariables, cap = F, case = "A") {
  log2nn = with(inputData, log2(cell_count/cell_count__time0))
  log2nn_ctrl = with(inputData, log2(cell_count__ctrl/cell_count__time0))
  GR = 2^(log2nn/log2nn_ctrl) - 1
  input_edited = inputData
  input_edited$log10_concentration = log10(input_edited$concentration)
  input_edited$GR = GR
  tmp<-input_edited[,groupingVariables, drop = F]
  experimentNew = (apply(tmp,1, function(x) (paste(x,collapse=" "))))
  if(cap == T) {
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

.GRlogisticFit = function(inputData, groupingVariables, force = F, cap = F) {
  experiments = levels(inputData$experiment)
  parameters = matrix(data = NA, ncol = 3, nrow = length(experiments))
  parameters = as.data.frame(parameters)
  colnames(parameters) = c('Hill','GRinf','GEC50')
  if(length(groupingVariables) > 0) {
    metadata = matrix(data = NA, ncol = length(groupingVariables), nrow = length(experiments))
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
    GRmax[i] = min(max_concs$GR, na.rm = T)
    #     metadata[i,] = data_exp[1,1:5]
    if(!is.null(metadata)) {
      metadata[i,] = data_exp[1,groupingVariables, drop = F]
    }
    GR_mean[i] = mean(data_exp$GR, na.rm = T)
    #===== constrained fit ============
    c = unique(data_exp$concentration)
    priors = c(2, 0.1, median(c))
    lower = c(.1, -1, min(c)*1e-2)
    upper = c(5, 1, max(c)*1e2)
    if(dim(data_exp)[1] > 1) {
      controls = drc::drmc()
      controls$relTol = 1e-06 #set to match MATLAB code
      controls$errorm = F
      controls$noMessage = T
      controls$rmNA = T
      output_model_new = try(drc::drm(GR~concentration, experiment, data=data_exp, fct=drc::LL.3u(names=c('Hill','GRinf','GEC50')), start = priors, lowerl = lower, upperl = upper, control = controls, na.action = na.omit))
      if(class(output_model_new)!="try-error") {
        parameters[i,] = c(as.numeric(coef(output_model_new)))
        # F-test for the significance of the sigmoidal fit
        Npara = 3 # N of parameters in the growth curve
        Npara_flat = 1 # F-test for the models
        RSS2 = sum(residuals(output_model_new)^2, na.rm = T)
        RSS1 = sum((data_exp$GR - mean(data_exp$GR, na.rm = T))^2, na.rm = T)
        df1 = (Npara - Npara_flat)
        df2 = (length(na.omit(data_exp$GR)) - Npara + 1)
        f_value = ((RSS1-RSS2)/df1)/(RSS2/df2)
        f_pval = pf(f_value, df1, df2, lower.tail = F)
        pval[i] = f_pval
        R_square[i] = 1 - RSS2/RSS1
      }
    }
    #==================================

    # Trapezoid rule for integration of GR_AOC
    GRavg = NULL
    for(j in 1:length(concs)) {
      data_trapz = data_exp[data_exp$concentration == concs[j],]
      GRavg[j] = mean(data_trapz$GR, na.rm = T)
    }
    AOC[i] = sum((1 - (GRavg[1:(length(GRavg)-1)]+GRavg[2:length(GRavg)])/2)*diff(log10(concs), lag = 1), na.rm = T)/(log10(concs[length(concs)]) - log10(concs[1]))
  }

  # Calculate GR50 from parameters
  parameters$GR50 = with(parameters, GEC50*((1-GRinf)/(0.5-GRinf) - 1)^(1/Hill))
  parameters$GRmax = GRmax
  parameters$GR_AOC = AOC
  parameters$r2 = R_square
  if(is.null(pval)) {pval = NA}
  parameters$pval = pval
  # Re-order rows to match reference_output
  parameters$experiment = experiments
  # Threshold for F-test pval
  pcutoff = ifelse(force == F, .05 , 1)
  for(i in 1:dim(parameters)[1]) {
    if(!is.na(parameters$pval[i])) {
      parameters$fit[i] = ifelse(parameters$pval[i] >= pcutoff | is.na(parameters$GEC50[i]), "flat","sigmoid")
    } else {
      parameters$fit[i] = ifelse(is.na(parameters$GEC50[i]), "flat", "sigmoid")
    }
  }
  # changed to above code to deal with NAs
  #parameters$fit = with(parameters, ifelse(pval >= pcutoff | is.na(GEC50), "flat","sigmoid"))
  # Add specified values for flat fits: GEC50 = 0, Hill = 0.01 and GR50 = +/- Inf
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
  parameters = parameters[,c('GR50','GRmax','GR_AOC','GEC50','GRinf','Hill', 'r2','pval','experiment', 'fit','flat_fit')]
  if(!is.null(metadata)) {
    parameters = cbind(metadata, parameters)
  }
  return(parameters)
}

.GRlogistic_3u = function(c){GRinf + (1 - GRinf)/(1 + (c/GEC50)^Hill)}

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

.convert = function(inputData, case, concentration = 'concentration', cell_count = 'cell_count', time = 'time', cell_count__time0 = 'cell_count__time0', cell_count__ctrl = 'cell_count__ctrl') {
  if(case == "A") {
    if(length(intersect(colnames(inputData), c('concentration', 'cell_count', 'cell_count__ctrl', 'cell_count__time0'))) == 4) {
      return(inputData)
    } else {
      inputData$concentration = inputData[[concentration]]
      inputData$cell_count = inputData[[cell_count]]
      inputData$cell_count__ctrl = inputData[[cell_count__ctrl]]
      inputData$cell_count__time0 = inputData[[cell_count__time0]]
      delete_cols = setdiff(c(concentration, cell_count, cell_count__time0, cell_count__ctrl), c('concentration', 'cell_count', 'cell_count__ctrl', 'cell_count__time0'))
      delete_col_nums = which(colnames(inputData) %in% delete_cols)
      inputData = inputData[,-delete_col_nums]
      return(inputData)
    }
  } else if(case == "C") {
    delete_cols = which(colnames(inputData) %in% c(concentration, cell_count))
    keys = colnames(inputData)[-delete_cols]
    time0 = inputData[inputData[[time]] == 0, c(keys, cell_count)]
    ctrl = inputData[inputData[[concentration]] == 0 & inputData[[time]] > 0, c(keys, cell_count)]
    data = inputData[inputData[[concentration]] != 0 & inputData[[time]] > 0, ]
    time0_keys = NULL
    ctrl_keys = NULL
    for(i in 1:length(keys)) {
      time0_keys[i] = length(intersect(time0[[ keys[i] ]], data[[ keys[i] ]])) > 0
      ctrl_keys[i] = length(intersect(ctrl[[ keys[i] ]], data[[ keys[i] ]])) > 0
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
      trimmed_mean = .trim_mean(ctrl[ctrl$key == key,][[cell_count]], 50)
      data[data$key_ctrl == key, 'cell_count__ctrl'] = trimmed_mean
    }

    for(key in unique(time0$key)) {
      trimmed_mean = .trim_mean(time0[time0$key == key,][[cell_count]], 50)
      data[data$key_time0 == key, 'cell_count__time0'] = trimmed_mean
    }

    delete_cols = which(colnames(data) %in% c('key_ctrl', 'key_time0'))
    data = data[, -delete_cols]

    row.names(data) = 1:dim(data)[1]
    inputData = data
    if(length(intersect(colnames(inputData), c('concentration', 'cell_count', 'time'))) == 3) {
      return(inputData)
    } else {
      inputData$concentration = inputData[[concentration]]
      inputData$cell_count = inputData[[cell_count]]
      inputData$time = inputData[[time]]
      delete_cols = setdiff(c(concentration, cell_count, time), c('concentration', 'cell_count', 'time'))
      delete_col_nums = which(colnames(inputData) %in% delete_cols)
      inputData = inputData[,-delete_col_nums]
      return(inputData)
    }
  }
}

#' Extract GR parameters from a dataset
#'
#' This function takes in a dataset with information about concentration,
#' cell counts over time, and additional grouping variables for a dose-response assay
#' and calculates growth-rate inhibition (GR) metrics for each experiment in
#' the dataset. The data must be in a specific format: either that specified
#' by case "A" or case "C" described here
#' \url{https://github.com/sorgerlab/gr50_tools/blob/master/README.md} and in
#' the details below, although the column names of the key variables
#' (concentration, cell_count, etc.) can be specified in the command.
#'
#' @param inputData a data table in one of the specified formats (Case A or
#' Case C). See details below or
#' \url{https://github.com/sorgerlab/gr50_tools/blob/master/README.md}
#' @param groupingVariables a vector of column names from inputData. All of the
#'  columns in inputData except for those identified here will be averaged over.
#' @param GRtable a character string ("param", "GR", or "both") identifying
#'  whether to return, respectively, only the table of parameters, only the
#'  table of GR values, or both tables. Default is "both".
#' @param force a logical value indicating whether to attempt to "force" a
#' sigmoidal fit, i.e. whether to allow fits with F-test p-values greater
#' than .05
#' @param cap a logical value indicating whether to cap GR values in between
#'  -1 and 1 (the range of the GR dose-response curve). If true, all GR values
#'  greater than 1 will be set to 1 and all values less than -1 will be set
#'  to -1.
#' @param case either "A" or "C", indicating the format of the input data. See
#' here \url{https://github.com/sorgerlab/gr50_tools/blob/master/README.md}
#' or details below for descriptions of these formats.
#' @param concentration (Case A and Case C) the name of the column with
#' concentration values (not log transformed) of the perturbagen on which
#' dose-response curves will be evaluated
#' @param cell_count (Case A and Case C) the name of the column with the
#' measure of cell number (or a surrogate of cell number such as CellTiter-Glo®
#' staining) after treatment
#' @param time (Case C) the name of the column with the time at which a cell
#' count is observed
#' @param cell_count__time0 (Case A) the name of the column with Time 0 cell
#' counts - the measure of cell number in untreated wells grown in parallel
#' until the time of treatment
#' @param cell_count__ctrl (Case A) the name of the column with the Control
#' cell count - the measure of cell number in control (e.g. untreated or
#' DMSO-treated) wells from the same plate
#' @return By default, a list of three elements is returned: 1) a data table of
#' the original data, converted to the style of Case A, with GR values for each
#' experiment 2) a data table of GR metrics parameters (GR50, GRmax, etc.) as
#' well as goodness of fit measures. 3) the vector of grouping variables (this
#' is used for other functions such as \code{\link{GRdrawDRC}}). If GRtable is
#' equal to "GR", then only (1) and (3) are returned. If GRtable is set to
#' "param", then only (2) and (3) are returned.
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
#' @note
#' To see the underlying code, use (\code{getAnywhere(.GRlogistic_3u)}),
#' (\code{getAnywhere(.GRcalculate)}), and (\code{getAnywhere(.GRlogisticFit)})
#' @seealso See \code{\link{drm}} for the general logistic fit function that
#' solves for the parameters GRinf, GEC50, and Hill slope. See \code{\link{drmc}} for
#' options of this function. Use the functions \code{\link{GRdrawDRC}},
#' \code{\link{GRbox}}, and \code{\link{GRscatter}} to create visualizations
#' using the output from this function. For online GR calculator and browser, see
#' \url{http://www.grcalculator.org}.
#' @references Hafner, M., Niepel, M. Chung, M. and Sorger, P.K., "Growth Rate Inhibition Metrics Correct For Confounders In Measuring Sensitivity To Cancer Drugs". \emph{Nature Methods} 13.6 (2016): 521-527.
#' \link{http://dx.doi.org/10.1038/nmeth.3853}
#' @examples
#' # Load Case A (example 1) input
#' data("inputCaseA")
#' head(inputCaseA)
#' # Run GRfit function with case = "A"
#' output1 = GRfit(inputData = inputCaseA, groupingVariables = c('cell_line','agent', 'perturbation','replicate', 'time'))
#' # See parameter table
#' head(output1$parameter_table)
#' # See GR values table
#' head(output1$gr_table)
#' # Load Case C (example 4) input
#' # Same data, different format
#' data("inputCaseC")
#' head(inputCaseC)
#' output4 = GRfit(inputData = inputCaseC, groupingVariables = c('cell_line','agent', 'perturbation','replicate', 'time'), case = "C")
#'
#' @export

GRfit = function(inputData, groupingVariables, GRtable = 'both', force = F, cap = F, case = "A", concentration = 'concentration', cell_count = 'cell_count', time = 'time', cell_count__time0 = 'cell_count__time0', cell_count__ctrl = 'cell_count__ctrl') {
  inputData = .convert(inputData, case, concentration, cell_count, time, cell_count__time0, cell_count__ctrl)
  gr_table = .GRcalculate(inputData, groupingVariables, cap, case)
  if(GRtable == 'GR') {
    returnList = list(gr_table, groupingVariables)
    names(returnList) = c(gr_table, groupingVariables)
    class(returnList) = 'GRfit'
    return(returnList)
  }
  parameter_table = .GRlogisticFit(gr_table, groupingVariables, force, cap)
  if(GRtable == 'param'){
    returnList = list(parameter_table, groupingVariables)
    names(returnList) = c('parameter_table', 'groupingVariables')
    class(returnList) = 'GRfit'
    return(returnList)
  } else if(GRtable == 'both') {
    returnList = list(gr_table,parameter_table, groupingVariables)
    names(returnList) = c('gr_table', 'parameter_table', 'groupingVariables')
    class(returnList) = 'GRfit'
    return(returnList)
  }
}

