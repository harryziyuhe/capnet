#' Normalize user input into a stats family object
#' 
#' @param family optional character scalar (e.g. "binomial"), 
#' a function (e.g. \code{stats::binomial}), or a family object 
#' (e.g. \code{stats::binomial()}).
#' 
#' @return A stats family object

normalize_family <- function(family = NULL) {
  if (is.null(family)) {
    return (stats::gaussian())
  }
  
  if (is.character(family) && length(family) == 1L) {
    fam_fun <- get(family, envir = asNamespace("stats"), inherits = FALSE)
    fam <- fam_fun()
    .assert_is_family_object(fam)
    return(fam)
  }
  
  if (is.function(family)) {
    fam <- family()
    .assert_is_family_object(fam)
    return(fam)
  }
  
  if (.is_family_object(family)) {
    return(family)
  }
  
  stop(
    "family must be a stats family object (e.g., binomial()), a stats family function ",
    "(e.g., stats::binomial), a character string (e.g., \"binomial\"), or NULL.",
    call. = FALSE
  )
}

#' Validate the family/link combination is supported
#' 
#' @param fam stats family object
#' 
#' @return invisibly \code{TRUE}, errors otherwise

validate_family_supported <- function(fam){
  .assert_is_family_object(fam)
  
  f <- tolower(fam$family %||% "")
  link <- tolower(fam$link %||% "")
  
  if (identical(f, "gaussian")) {
    if (!identical(link, "identity")) {
      stop(sprintf("Unsupported gaussian link '%s'. Currently supported: identity.", fam$link), call. = FALSE)
    }
    return(invisible(TRUE))
  }
  
  if (identical(f, "binomial")) {
    if (!identical(link, "logit")) {
      stop(sprintf("Unsupported binomial link '%s'. Currently supported: logit.", fam$link), call. = FALSE)
    }
    return(invisible(TRUE))
  }
  
  if (identical(f, "poisson")) {
    if (!identical(link, "log")) {
      stop(sprintf("Unsupported poisson link '%s'. Currently supported: log.", fam$link), call. = FALSE)
    }
    return(invisible(TRUE))
  }
  
  if (identical(f, "gamma")) {
    if (!identical(link, "inverse")) {
      stop(sprintf("Unsupported Gamma link '%s'. Currently supported: inverse.", fam$link), call. = FALSE)
    }
    return(invisible(TRUE))
  }
  
  stop(sprintf("Unsupported family '%s'. Supported families: gaussian, binomial, poisson.", fam$family), call. = FALSE)
}

# ------------------------------------------------------------------------------
# Helpers (internal)
# ------------------------------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

.is_family_object <- function(x) {
  is.list(x) &&
    !is.null(x$family) &&
    !is.null(x$link) &&
    is.function(x$linkinv) &&
    is.function(x$linkfun) &&
    is.function(x$mu.eta)
}

.assert_is_family_object <- function(x) {
  if (!.is_family_object(x)) {
    stop("Invalid family object. Expected a stats family (e.g., stats::binomial()).", call. = FALSE)
  }
  invisible(TRUE)
}


























