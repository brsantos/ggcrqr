#' Data about male breast cancer in Brazil
#'
#' A dataset containing information about men with breast cancer in different
#'  stages with their time to death, collected between 2000 and 2019 in the
#'  state of São Paulo, in Brazil.
#'
#' @format A data frame with 872 rows and 5 variables:
#' \describe{
#'   \item{time_to_d}{time to death, in years.}
#'   \item{age}{age of each person in years.}
#'   \item{cens}{0, if the time is not censored and 1 otherwise.}
#'   \item{age_group}{age group, 1: age < 55 years; 2: 55 years <= age <= 65
#'    years; 3: age > 65 years.}
#'   \item{stage_c}{stage of the cancer, with possibilities: I, II, III e IV.}
#'   \item{surgery}{if the person had surgery as treatment, yes or no.}
#'   \item{radio}{if the person had radiation therapy as treatment, yes or no.}
#'   \item{chemo}{if the person had chemotherapy as treatment, yes or no.}
#'   \item{hormone}{if the person had hormone therapy as treatment, yes or no.}
#' }
#' @source \url{http://www.fosp.saude.sp.gov.br/publicacoes/downloadarquivos}
"m_breast_cancer"
