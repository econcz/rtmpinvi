#' Solve an interactive tabular matrix estimation problem via Convex Least
#' Squares Programming (CLSP).
#'
#' @param ival NULL, numeric matrix, or data.frame.
#'   Prior information on known cell values. If supplied and not entirely
#'   missing, \code{ival} is flattened and used to construct
#'   \code{b_val} and the corresponding identity-subset model matrix
#'   \code{M} internally. Missing entries (\code{NA}) are ignored.
#'   If all entries of \code{ival} are \code{NA}, no prior information
#'   is used and \code{b_val} and \code{M} are set to \code{NULL}.
#'   When \code{ival} is provided, it overrides any \code{b_val} or
#'   \code{M} arguments supplied through \code{...}.
#'
#' @param ibounds NULL, \code{numeric(2)}, or list of \code{numeric(2)}.
#'   Dynamic cell-value bounds passed to \code{tmpinv(bounds = ...)}.
#'   The object supplied to \code{ibounds} may be created or modified
#'   programmatically (for example within \code{preestimation()}).
#'   If a single pair \code{c(low, high)} is provided, it is applied
#'   uniformly to all cells. Alternatively, a list of pairs may be
#'   supplied to specify cell-specific bounds with others set to NA.
#'   When \code{ibounds} is non-NULL, it overrides any \code{bounds}
#'   argument supplied through \code{...}.
#'
#' @param preestimation NULL or function.
#'   A function executed prior to model estimation. If supplied,
#'   it is called as \code{preestimation(ival)} and may perform
#'   arbitrary preparatory steps, such as constructing dynamic
#'   bounds or modifying objects in the calling environment. The return
#'   value is ignored.
#'
#' @param postestimation NULL or function.
#'   A function executed after model estimation. For a full model,
#'   it is called as \code{postestimation(model)}. For reduced
#'   (block-wise) models, it is called as
#'   \code{postestimation(model[[i]], i)} for each block index
#'   \code{i}. The return value is ignored.
#'
#' @param update logical scalar, default = \code{FALSE}.
#'   If \code{TRUE} and \code{ival} is supplied, missing entries
#'   (\code{NA}) in \code{ival} are replaced by the corresponding
#'   fitted values from \code{tmpinvi$result$x}. The updated matrix is
#'   returned in the \code{tmpinvi$data} component.
#'   If \code{FALSE}, the \code{data} component contains the fitted
#'   solution matrix \code{tmpinvi$result$x}.
#'
#' @param ... Additional arguments passed to \code{rtmpinv::tmpinv()}.
#'
#' @return
#' An object of class \code{"tmpinvi"} with components:
#' \itemize{
#'   \item \code{result}: a fitted object of class \code{"tmpinv"}.
#'   \item \code{data}: the processed matrix (either the fitted solution
#'     \code{x} or the updated \code{ival}, depending on \code{update}).
#' }
#'
#' @seealso \link[rtmpinv]{tmpinv}
#'
#' @examples
#' \donttest{
#'   RNGkind("L'Ecuyer-CMRG")
#'   set.seed(123456789)
#'
#'   iso2 <- c("CN", "DE", "JP", "NL", "US")
#'   T    <- 10L
#'   year <- (as.integer(format(Sys.Date(), "%Y")) - T) + seq_len(T)
#'   m    <- length(iso2)
#'
#'   df <- expand.grid(year = year, iso2 = iso2, KEEP.OUT.ATTRS = FALSE)
#'   df <- df[order(df$year, df$iso2), ]
#'
#'   ex_cols <- paste0("EX_", iso2)
#'   df[ex_cols] <- NA_real_
#'   df$EX <- NA_real_
#'   df$IM <- NA_real_
#'
#'   X_true <- vector("list", length(year))
#'   names(X_true) <- as.character(year)
#'
#'   for (t in seq_along(year)) {
#'     scale <- 1000 * (1.05^(t - 1L))
#'     X <- matrix(runif(m * m, 0, scale), m, m)
#'     diag(X) <- 0
#'     X_true[[t]] <- X
#'
#'     rows <- ((t - 1L) * m + 1L):((t - 1L) * m + m)
#'     df$EX[rows] <- rowSums(X)
#'     df$IM[rows] <- colSums(X)
#'
#'     miss <- matrix(runif(m * m) > 0.5, m, m)
#'     X[miss] <- NA_real_
#'     df[rows, ex_cols] <- X
#'   }
#'
#'   cv <- qnorm(0.975)
#'
#'   for (nm in ex_cols) {
#'     fit <- lm(df[[nm]] ~ year * iso2, data = df, na.action = na.exclude)
#'     pr  <- predict(fit, df, se.fit = TRUE)
#'     ub  <- pr$fit + cv * pr$se.fit
#'     ub[ub < 0] <- NA_real_
#'     df[[paste0("_", nm, "_lb")]] <- 0
#'     df[[paste0("_", nm, "_ub")]] <- ub
#'   }
#'
#'   make_bounds <- function(lb, ub)
#'     Map(function(a, b) c(a, b), lb, ub)
#'
#'   df_out <- df
#'
#'   for (step in 1:2) {
#'     for (y in year) {
#'       idx  <- df_out$year == y
#'       d    <- df_out[idx, ]
#'       ival <- as.matrix(d[ex_cols])
#'
#'       lb <- as.vector(t(as.matrix(d[paste0("_EX_", iso2, "_lb")])))
#'       ub <- as.vector(t(as.matrix(d[paste0("_EX_", iso2, "_ub")])))
#'
#'       fit <- tmpinvi(
#'         ival = ival,
#'         ibounds = make_bounds(lb, ub),
#'         b_row = d$EX,
#'         b_col = d$IM,
#'         alpha = 1.0,
#'         update = TRUE
#'       )
#'
#'       df_out[idx, ex_cols] <- fit$data
#'     }
#'   }
#'
#'   drop_cols <- grep("^_EX_.*_(lb|ub)$", names(df_out), value = TRUE)
#'   df_out[drop_cols] <- NULL
#'   df_out
#' }
#' @export
tmpinvi <- function(ival=NULL, ibounds=NULL, preestimation=NULL,
                    postestimation=NULL, update=FALSE, ...) {
  dots     <- list(...)
  # adjust and preprocess options
  if (!is.null(ival)) {
    if (is.data.frame(ival)) ival <- as.matrix(ival)
    if (!is.matrix(ival))    stop("ival must be a data.frame or a 2D matrix.")
    if (!is.null(names(dots)))
      dots <- dots[!(names(dots) %in% c("b_val", "M"))]
  }
  if (!is.null(ibounds) && !is.null(names(dots)))
    dots <- dots[setdiff(names(dots), "bounds")]
  if (!is.null(preestimation) && !is.function(preestimation))
    stop("preestimation must be a function.")
  if (!is.null(postestimation) && !is.function(postestimation))
    stop("postestimation must be a function.")
  
  # run preestimation function (= multiple commands)
  if (!is.null(preestimation))
    preestimation(ival)
  # perform estimation
  M     <- NULL
  b_val <- NULL
  if (!is.null(ival)) {
    if (!all(is.na(ival))) {
      b_val <- as.vector(t(ival))
      M     <- diag(length(b_val))[!is.na(b_val), , drop = FALSE]
      b_val <- b_val[!is.na(b_val)]
    }
  }
  bounds        <- if (!is.null(ibounds)) ibounds else NULL
  result        <- list()
  result$result <- do.call(rtmpinv::tmpinv,
                           c(list(b_val=b_val, M=M, bounds=bounds),
                             dots))
  # run postestimation command/program (= multiple commands)
  if (!is.null(postestimation)) {
    if (isTRUE(result$result$full)) postestimation(result$result$model)
    else for (i in seq_along(result$result$model)) {
      postestimation(result$result$model[[i]], i)
    }
  }
  # generate, update, or replace data from result$x
  if (isTRUE(update)) {
    if (!is.null(ival)) {
      if (nrow(ival) != nrow(result$result$x)                                ||
          ncol(ival) != ncol(result$result$x))
        stop("Dimensions of ival and result$x do not match.")
      ival[is.na(ival)] <- result$result$x[is.na(ival)]
      result$data      <- ival
    } else result$data <- result$result$x
  } else   result$data <- result$result$x
  
  class(result) <- "tmpinvi"
  result
}
#' @export
print.tmpinvi <- function(x, ...) {
  tmp <- if (!is.null(x$result)) x$result else x
  cat("Call:\n")
  if (!is.null(tmp$call)) print(tmp$call) else cat("tmpinvi(...)\n")
}
#' @export
summary.tmpinvi <- function(object, ...) {
  tmp <- if (!is.null(object$result)) object$result else object
  if (inherits(tmp$model, "clsp"))
    return(summary(tmp$model))
  dots <- list(...)
  idx  <- as.integer(dots$i)
  if        (is.null(dots$i))
    stop("Reduced model: please supply the block index using i=#.")
  if (idx < 1L || idx > length(tmp$model))
    stop(sprintf("i must be in 1..%d for reduced model.", length(tmp$model)))
  return(summary(tmp$model[[idx]]))
}
#' @export
print.summary.tmpinvi <- function(x, ...) {
  tmp <- if (!is.null(x$result)) x$result else x
  if (inherits(tmp$model, "clsp"))
    return(print(tmp$model))
  dots <- list(...)
  idx  <- as.integer(dots$i)
  if        (is.null(dots$i))
    stop("Reduced model: please supply the block index using i=#.")
  if (idx < 1L || idx > length(tmp$model))
    stop(sprintf("i must be in 1..%d for reduced model.", length(tmp$model)))
  return(print(tmp$model[[idx]]))
}
