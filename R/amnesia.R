#' @title
#' Amnesia and other drug reactions in the MHRA pharmacovigilance spontaneous
#' reporting system
#' 
#' @description
#' `r lifecycle::badge('deprecated')`
#' 
#' For each of 2,446 drugs in the MHRA database (column 1), 
#' the number of cases with amnesia as an adverse event (column 2), 
#' and the number of cases with other adverse event for this drug (column 3).
#' In total, 682,648 adverse drug reactions were reported, among them 2,044
#' cases of amnesia.
#' 
#' **Note**: In future versions, this dataset will be removed. Please use the
#' [`amnesia`][DiscreteDatasets::amnesia] dataset from package
#' [`DiscreteDatasets`][DiscreteDatasets::DiscreteDatasets-package].
#' 
#' @usage data(amnesia)
#'
#' @details 
#' The data was collected from the Drug Analysis Prints published 
#' by the Medicines and Healthcare products Regulatory Agency (MHRA),
#' by Heller & Gur. See references for more details.
#'
#' @format A data frame with 2,446 rows representing drugs with the following
#' three columns:
#' \describe{
#'   \item{DrugName}{The name of the drug.}
#'   \item{AmnesiaCases}{Number of the amnesia cases reported for the drug.}
#'   \item{OtherAdverseCases}{Number of other adverse drug reactions reported
#'                            for the drug.}
#' }
#' @source
#' [Drug Analysis Prints on MHRA site](https://yellowcard.mhra.gov.uk/idaps)
#' 
#' @section References:
#' R. Heller and H. Gur (2011). False discovery rate controlling procedures
#'   for discrete tests. 
#'   [arXiv:1112.4627v2](https://arxiv.org/abs/1112.4627v2) (preprint).
"amnesia"
