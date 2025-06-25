#' @title Summary method for POM objects
#' @description This S3 method provides a minimal summary object for custom null models...
#' @param object An object of class `POM`.
#' @param ... Additional arguments passed to other methods (unused).
#' @return A list containing a single element: `dispersion`.
#' @method summary POM
#' @export
summary.POM <- function(object, ...) {
  list(dispersion = object$dispersion)
}