#' Plot distribution AIC fit values
#'
#' This function evaluates and compares the goodness of fit for various statistical distributions
#' applied to a given dataset using the Akaike Information Criterion (AIC). It includes transformations
#' and several distribution fittings to identify the best model by AIC score.
#'
#' @param values A numeric vector containing the data for which the distribution fitting should be assessed.
#' @return A ggplot object displaying a bar chart of AIC scores for each distribution tested, arranged from the lowest to the highest AIC.
#' @export
#' @examples
#' # Generating some random data
#' set.seed(123)
#' data <- rnorm(100)
#' # Testing distributions on the generated data
#' testDistributions(data)
#'
#' @importFrom tibble tibble
#' @importFrom dplyr arrange
#' @importFrom ggplot2 aes ggplot geom_bar coord_flip xlab ylab ggtitle
#' @importFrom fitdistrplus fitdist
#' @importFrom brms theme_black
testDistributions <- function(values){
  suppressWarnings(
    tibble::tibble(
      aic = c(
        fitdistrplus::fitdist(values, "norm")$aic,
        fitdistrplus::fitdist(as.numeric(na.omit(log(abs(values) + 1))), "norm")$aic,
        fitdistrplus::fitdist(as.numeric(na.omit(sqrt(abs(values) + 1))), "norm")$aic,
        fitdistrplus::fitdist(abs(values) + 1, "gamma")$aic,
        fitdistrplus::fitdist(abs(values) + 1, "weibull")$aic,
        fitdistrplus::fitdist(abs(values) + 1, "lnorm")$aic,
        fitdistrplus::fitdist(abs(values) + 1, "t", start = list(df = 2))$aic,
        fitdistrplus::fitdist(log(abs(values) + 1), "t", start = list(df = 2))$aic,
        fitdistrplus::fitdist(abs(values) + 1, "logis")$aic,
        fitdistrplus::fitdist(abs(values) + 1, "chisq", start = list(df = 82))$aic
      ),
      distribution = c(
        "Normal",
        "Normal, log transform",
        "Normal, sqrt transform",
        "Gamma",
        "Weibull",
        "Lognormal",
        "t, df=2",
        "t, df=82, with log transform",
        "Logistic",
        "Chi-square, df=82"
      )
    ) %>%
      dplyr::arrange(aic) %>%
      ggplot2::ggplot(mapping = ggplot2::aes(x = reorder(distribution, -aic), y = aic, color = distribution, fill = distribution)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::coord_flip() +
      ggplot2::xlab("Distribution & Transformation") +
      ggplot2::ggtitle("Comparison of fit metric (AIC) for Distributions and Transformations") +
      ggplot2::ylab("AIC") +
      brms::theme_black()
  )
}