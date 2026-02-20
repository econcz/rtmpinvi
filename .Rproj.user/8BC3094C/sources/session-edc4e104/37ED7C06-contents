#' Solve a tabular matrix estimation problem via Convex Least Squares
#' Programming (CLSP).
#'
#' @param S numeric matrix of size \eqn{(m + p) \times (m + p)}, optional.
#'   A diagonal sign-slack (surplus) matrix with entries in
#'   \eqn{\{0, \pm 1\}}.
#'   * \code{0} enforces equality (== \code{b_row} or \code{b_col}),
#'   * \code{1} enforces a lower-than-or-equal (\eqn{\le}) condition,
#'   * \code{-1} enforces a greater-than-or-equal (\eqn{\ge}) condition.
#'   The first \code{m} diagonal entries correspond to row constraints,
#'   and the remaining \code{p} correspond to column constraints.
#'
#' @param M numeric matrix of size \eqn{k \times (m p)}, optional.
#'   A model matrix, typically with entries in \eqn{\{0,1\}}. Each row
#'   defines a linear restriction on the flattened solution matrix.
#'   The corresponding right-hand-side values must be supplied in
#'   \code{b_val}. This block encodes known cell values.
#'
#' @param b_row numeric vector of length \code{m}.
#'   Right-hand-side vector of row totals.
#'
#' @param b_col numeric vector of length \code{p}.
#'   Right-hand-side vector of column totals.
#'
#' @param b_val numeric vector of length \code{k}.
#'   Right-hand-side vector of known cell values.
#'
#' @param i integer, default = \code{1}.
#'   Number of row groups.
#'
#' @param j integer, default = \code{1}.
#'   Number of column groups.
#'
#' @param zero_diagonal logical scalar, default = \code{FALSE}.
#'   If \code{TRUE}, enforces a structural zero diagonal.
#'
#' @param reduced integer vector of length \code{2}, optional.
#'   Dimensions of the reduced problem. If supplied, estimation is
#'   performed block-wise on contiguous submatrices. For example,
#'   \code{reduced = c(6,6)} yields \eqn{5 \times 5} blocks with one
#'   slack row and one slack column (edge blocks may be smaller).
#'
#' @param symmetric logical scalar, default = \code{FALSE}.
#'   If TRUE, enforces symmetry of the estimated matrix via
#'   \code{x <- 0.5 * (x + t(x))}. This applies to \code{tmpinv$x} only.
#'   For symmetry in the model, add explicit symmetry rows to
#'   \code{M} instead of using this flag.
#'
#' @param bounds NULL, \code{numeric(2)}, or list of \code{numeric(2)}.
#'   Bounds on cell values. If a single pair \code{c(low, high)} is
#'   given, it is applied to all \eqn{m p} cells.
#'   Example: \code{c(0, NA)}.
#'
#' @param replace_value numeric scalar or \code{NA}, default = \code{NA}.
#'   Final replacement value for any cell that violates the bounds by
#'   more than the given tolerance.
#'
#' @param tolerance numeric scalar, default = \code{sqrt(.Machine$double.eps)}.
#'   Convergence tolerance for bounds.
#'
#' @param iteration_limit integer, default = \code{50}.
#'   Maximum number of iterations allowed in the refinement loop.
#'
#' @param r integer scalar, default = \code{1}  
#'   Number of refinement iterations for the first step of the CLSP estimator.
#'
#' @param final logical scalar, default = \code{TRUE}  
#'   If \code{FALSE}, only the first step of the CLSP estimator is performed.
#'
#' @param alpha numeric scalar, numeric vector, or \code{NULL},
#'   Regularization parameter for the second step of the CLSP estimator.
#'
#' @param ... Additional arguments passed to the \pkg{rclsp} solver.
#'
#' @return
#' An object of class \code{"tmpinv"} containing the fitted CLSP
#' model (\code{tmpinv$model}) and solution matrix (\code{tmpinv$x}).
#'
#' @note
#'   1. In the reduced model, \code{S} is ignored. Slack behaviour is
#'      inferred from block-wise marginal totals.  
#'      Likewise, \code{M} must be a unique row subset of an identity
#'      matrix (diagonal-only). Non-diagonal model matrices cannot be
#'      mapped into reduced blocks.
#'   2. Internal keyword arguments \code{b_lim} and \code{C_lim} are
#'      passed to \code{.tmpinv.instance()} and contain cell-value
#'      bounds. These arguments are ignored in the reduced model.
#'
#' @seealso \link[rclsp]{clsp}
#' @seealso \link[CVXR]{CVXR-package}
#'
#' @examples
#' \donttest{
#'   ## Example 1: AP/TM reconstruction on a symmetric 20x20 matrix
#'   ## (10 percent known entries)
#'
#'   RNGkind("L'Ecuyer-CMRG")
#'   set.seed(123456789)
#'
#'   m <- 20L
#'   p <- 20L
#'
#'   # sample (dataset)
#'   X_true <- abs(matrix(rnorm(m * p), nrow = m, ncol = p))
#'   X_true <- 0.5 * (X_true + t(X_true))              # symmetric
#'
#'   idx <- sample.int(
#'       m * p,
#'       size = max(1L, floor(0.1 * (m * p))),         # 10 percent known
#'       replace = FALSE
#'   )
#'
#'   M     <- diag(m * p)[idx, , drop = FALSE]
#'   b_row <- rowSums(X_true)
#'   b_col <- colSums(X_true)
#'   b_val <- matrix(as.numeric(X_true)[idx], ncol = 1L)
#'
#'   # model (unique MNBLUE estimator)
#'   result <- tmpinv(
#'       M = M,
#'       b_row = b_row,
#'       b_col = b_col,
#'       b_val = b_val,
#'       bounds = c(0, NA),                            # non-negativity
#'       symmetric = TRUE,
#'       r = 1L,
#'       alpha = 1.0
#'   )
#'
#'   # coefficients
#'   print("true X:")
#'   print(round(X_true, 4))
#'
#'   print("X_hat:")
#'   print(round(result$x, 4))
#'
#'   # numerical stability
#'   print("\nNumerical stability:")
#'   print(paste("  kappaC :", result$model$kappaC))
#'   print(paste("  kappaB :", result$model$kappaB))
#'   print(paste("  kappaA :", result$model$kappaA))
#'
#'   # diagnostics
#'   print("\nGoodness-of-fit:")
#'   print(paste("  NRMSE :", result$model$nrmse))
#'   print(paste("  Diagnostic band (min):", min(result$model$x_lower)))
#'   print(paste("  Diagnostic band (max):", max(result$model$x_upper)))
#'
#'   # bootstrap NRMSE t-test
#'   tt <- rclsp::ttest(
#'       result$model,
#'       sample_size = 30L,
#'       seed = 123456789L,
#'       distribution = rnorm,
#'       partial = TRUE
#'   )
#'   print("\nBootstrap t-test:")
#'   print(tt)
#'
#'   ## Example 2: AP/TM reconstruction on a 40x40 matrix
#'   ## with zero diagonal and reduced (20,20) submodels
#'   ## (20 percent known entries)
#'
#'   RNGkind("L'Ecuyer-CMRG")
#'   set.seed(123456789)
#'
#'   m <- 40L
#'   p <- 40L
#'
#'   # sample (dataset)
#'   X_true <- abs(matrix(rnorm(m * p), nrow = m, ncol = p))
#'   diag(X_true) <- 0                                 # zero diagonal
#'
#'   idx <- sample.int(
#'       m * p,
#'       size = max(1L, floor(0.2 * (m * p))),         # 20 percent known
#'       replace = FALSE
#'   )
#'
#'   M     <- diag(m * p)[idx, , drop = FALSE]
#'   b_row <- rowSums(X_true)
#'   b_col <- colSums(X_true)
#'   b_val <- matrix(as.numeric(X_true)[idx], ncol = 1L)
#'
#'   # model (reduced models of size 20x20)
#'   result <- tmpinv(
#'       M = M,
#'       b_row = b_row,
#'       b_col = b_col,
#'       b_val = b_val,
#'       zero_diagonal = TRUE,
#'       reduced = c(20L, 20L),
#'       bounds = c(0, NA),
#'       r = 1L,
#'       alpha = 1.0
#'   )
#'
#'   print("true X:")
#'   print(round(X_true, 4))
#'
#'   print("X_hat:")
#'   print(round(result$x, 4))
#'
#'   # numerical stability across submodels
#'   kC <- sapply(result$model, function(CLSP) CLSP$kappaC)
#'   kB <- sapply(result$model, function(CLSP) CLSP$kappaB)
#'   kA <- sapply(result$model, function(CLSP) CLSP$kappaA)
#'
#'   print("\nNumerical stability (min-max across models):")
#'   print(paste("  kappaC :", range(kC)))
#'   print(paste("  kappaB :", range(kB)))
#'   print(paste("  kappaA :", range(kA)))
#'
#'   # diagnostics (min-max)
#'   nrmse <- sapply(result$model, function(CLSP) CLSP$nrmse)
#'   x_low <- unlist(lapply(result$model, function(CLSP) CLSP$x_lower))
#'   x_up  <- unlist(lapply(result$model, function(CLSP) CLSP$x_upper))
#'
#'   print("\nGoodness-of-fit (min-max across models):")
#'   print(paste("  NRMSE :", range(nrmse)))
#'   print(paste("  Diagnostic band (min):", range(x_low)))
#'   print(paste("  Diagnostic band (max):", range(x_up)))
#'
#'   # bootstrap t-tests across all block models
#'   print("\nBootstrap t-tests:")
#'   tests <- lapply(
#'       result$model,
#'       function(CLSP) rclsp::ttest(
#'           CLSP,
#'           sample_size = 30L,
#'           seed = 123456789L,
#'           distribution = rnorm,
#'           partial = TRUE
#'       )
#'   )
#'   print(tests)
#' }
#'
#' @export
tmpinv <- function(S=NULL, M=NULL, b_row=NULL, b_col=NULL, b_val=NULL,
                   i=1L, j=1L, zero_diagonal=FALSE, reduced=NULL,
                   symmetric=FALSE, bounds=NULL, replace_value=NA_real_,
                   tolerance=sqrt(.Machine$double.eps), iteration_limit=50L,
                   r=1L, final=TRUE, alpha=NULL, ...) {
  dots     <- list(...)
  # (n_cells) Perform initial estimation and get cell count
  result   <- do.call(.tmpinv.instance,
                      c(list(b_row=b_row, b_col=b_col, b_val=b_val, i=i, j=j,
                             S=S, M=M, zero_diagonal=zero_diagonal,
                             reduced=reduced, symmetric=symmetric,
                             tolerance=tolerance,
                             iteration_limit=iteration_limit, r=r, final=final,
                             alpha=alpha),
                        dots))
  n_cells  <- nrow(result$x) * ncol(result$x)
  if        (is.null(bounds))      {
    bounds <- list(c(NA_real_, NA_real_))
  } else if (is.numeric(bounds) && length(bounds) == 2L)
    bounds <- list(bounds)
  if                              (length(bounds) == 1L) {
    bounds <- rep(bounds, n_cells)                     # replicate (low, high)
  } else if (length(bounds) !=   n_cells)
    stop(sprintf("Bounds length %d does not match number of variables %d.",
                 length(bounds), n_cells))
  bounds   <- lapply(bounds, function(v)   {
    if (length(v) != 2L)     stop("Each bounds element must have length 2.")
    vapply(v, function(x)    if (is.null(x) || is.na(x) || length(x) != 1L )
      NA_real_                 else           as.numeric(x), numeric(1L))
  })
  if (all(vapply(bounds,     function(v)    all(is.na(v)), logical(1L))))
    return(result)                                     # finish if unbounded
  
  # (result) Perform bound-constrained iterative refinement
  b_lim      <- matrix(rbind(vapply(bounds, function(v)  if (is.na(v[2]))  Inf
                                    else                           v[2],
                                    numeric(1L)),
                             vapply(bounds, function(v)  if (is.na(v[1])) -Inf
                                    else                           v[1],
                                    numeric(1L))),          ncol=     1L)
  C_lim      <- matrix(rbind(diag(n_cells), diag(n_cells)), ncol=n_cells)
  S          <- if (!is.null(S)) S  else    matrix(0,  nrow=length(b_row) +
                                                     length(b_col),
                                                   ncol=0L)
  S          <- rbind(cbind(S, matrix(0, nrow=nrow(S), ncol=      2 * n_cells)),
                      cbind(   matrix(0, nrow=n_cells, ncol=          ncol(S)),
                               diag(n_cells),
                               matrix(0, nrow=n_cells, ncol=          n_cells)),
                      -cbind(  matrix(0, nrow=n_cells, ncol=ncol(S) + n_cells),
                               diag(n_cells)))
  finite_rows   <- is.finite(b_lim[, 1L])              # drop rows with +-np.inf
  nonzero_cols  <- colSums(abs(S[c(rep(TRUE, nrow(S) - length(finite_rows)),
                                   finite_rows), ,
                                 drop=FALSE]))   > 0   # reduce S width
  for   (iter in seq_len(iteration_limit)) {
    M_idx  <- integer(0)
    b_val  <- numeric(0)
    x      <- matrix(result$x, ncol=1L)
    for (k in seq_len(n_cells)) {
      if ((!is.na(bounds[[k]][1]) && x[k] < bounds[[k]][1] - tolerance) ||
          (!is.na(bounds[[k]][2]) && x[k] > bounds[[k]][2] + tolerance)) {
        next                                           # skip out-of-bounds
      }
      M_idx  <- c(M_idx, k)
      b_val  <- c(b_val, x[k])
    }
    if (length(M_idx) < n_cells) {
      M      <- diag(n_cells)[M_idx, , drop=FALSE]
      result <- suppressWarnings(do.call(.tmpinv.instance,
                                         c(list(b_row=b_row, b_col=b_col,  b_val=b_val,
                                                b_lim=b_lim[finite_rows, , drop=FALSE], i=i, j=j,
                                                C_lim=C_lim[finite_rows, , drop=FALSE],
                                                S=S[c(rep(TRUE, nrow(S) -  nrow(b_lim)),
                                                      finite_rows),        nonzero_cols,
                                                    drop=FALSE],           M=M,
                                                zero_diagonal=zero_diagonal, reduced=reduced,
                                                symmetric=symmetric, tolerance=tolerance,
                                                iteration_limit=iteration_limit, r=r,
                                                final=final, alpha=alpha),
                                           dots)))
    } else {
      break
    }
  }
  
  # (result) Replace out-of-bound values with replace_value
  x        <- matrix(result$x, ncol=1L)
  x_lb     <- vapply(bounds, function(v) if (is.na(v[1])) -Inf
                     else                          v[1], numeric(1L))
  x_ub     <- vapply(bounds, function(v) if (is.na(v[2]))  Inf
                     else                          v[2], numeric(1L))
  x[(x < x_lb - tolerance) | (x > x_ub + tolerance)] <- replace_value
  result$x <- matrix(x, nrow=nrow(result$x),  ncol=ncol(result$x),
                     byrow=TRUE)
  
  # enforce (final) symmetry
  if (isTRUE(symmetric)) if (nrow(result$x) == ncol(result$x))
    result$x <-          0.5 *   (result$x  +     t(result$x))         else
      stop("symmetric=True requires a square matrix (m == p).")
  
  result
}
#' @export
print.tmpinv <- function(x, ...) {
  cat("Call:\n")
  if (!is.null(x$call)) print(x$call) else cat("tmpinv(...)\n")
}
#' @export
summary.tmpinv <- function(object, ...) {
  if (inherits(object$model, "clsp"))
    return(summary(object$model))
  dots <- list(...)
  idx  <- as.integer(dots$i)
  if        (is.null(dots$i))
    stop("Reduced model: please supply the block index using i=#.")
  if (idx < 1L || idx > length(object$model))
    stop(sprintf("i must be in 1..%d for reduced model.", length(object$model)))
  return(summary(object$model[[idx]]))
}
#' @export
print.summary.tmpinv <- function(x, ...) {
  if (inherits(x$model, "clsp"))
    return(print(x$model))
  dots <- list(...)
  idx  <- as.integer(dots$i)
  if        (is.null(dots$i))
    stop("Reduced model: please supply the block index using i=#.")
  if (idx < 1L || idx > length(     x$model))
    stop(sprintf("i must be in 1..%d for reduced model.", length(     x$model)))
  return(print(x$model[[idx]]))
}
################################################################################
# Ancillary functions
################################################################################
.tmpinv.instance <- function(S=NULL, M=NULL, b_row=NULL, b_col=NULL, b_val=NULL,
                             i=1L, j=1L, zero_diagonal=FALSE, reduced=NULL,
                             symmetric=FALSE, ...) {
  dots  <- list(...)
  # (m), (p) Process the parameters, assert conformity, and get dimensions
  if (is.null(b_row)         || is.null(b_col)        )
    stop("Both b_row and b_col must be provided.")
  if (length(b_row) < 2L     || length(b_col) < 2L    )
    stop("Minimum length for b_row and b_col is 2.")
  if (!all(is.finite(b_row)) || !all(is.finite(b_col)))
    stop("b_row and b_col must not contain +-Inf or NA")
  b_row <- matrix(as.numeric(b_row), ncol=1L, byrow=TRUE)
  b_col <- matrix(as.numeric(b_col), ncol=1L, byrow=TRUE)
  m     <- nrow(b_row) * i
  p     <- nrow(b_col) * j
  if (!is.null(S))                    {
    n_rows <- m + p  + (if (!is.null(dots$C_lim)) nrow(dots$C_lim)     else 0L)
    S      <- as.matrix(S)
    if (nrow(S) != n_rows)        stop(sprintf("S must have %d rows.", n_rows))
    if (!all(S == -1 | S == 0 | S == 1) ||
        max(colSums(abs(S))) > 1        ||
        max(rowSums(abs(S))) > 1) stop("S must be a zero-padded subset of +-I.")
  }
  if (!is.null(M))                    {
    if (is.null(b_val))           stop("Both M and b_val must be defined.")
    if (is.null(dim(M)) || length(dim(M)) == 1L)    M <- matrix(M, nrow=1L)
    if (ncol(M) != m * p)         stop(sprintf(paste0("M must have exactly %d ",
                                                      "columns."),  m * p))
  }
  if (!is.null(b_val))                {
    if (is.null(M))               stop("Both M and b_val must be defined." )
    if (!all(is.finite(b_val)))   stop("b_val must not contain +-Inf or NA.")
    b_val <- matrix(b_val, ncol=1L, byrow=TRUE)
  }
  if (!is.null(M) && !is.null(b_val)) {
    if (nrow(M) != nrow(b_val))   stop(sprintf(paste0("M and b_val must have ",
                                                      "the same number of ",
                                                      "rows: %d vs %d"),
                                               nrow(M), nrow(b_val)))
  }
  
  # perform full estimation and return the result
  if (is.null(reduced)) {
    result       <- list(full=TRUE,                   model=NULL,
                         x   =matrix(NA_real_,        nrow=m,    ncol=p))
    b_blocks     <- rbind(b_row,  b_col )
    if (!is.null(dots$b_lim)) b_blocks <- rbind(b_blocks, matrix(dots$b_lim,
                                                                 ncol=1L))
    if (!is.null(b_val))      b_blocks <- rbind(b_blocks, matrix(b_val,
                                                                 ncol=1L))
    b            <- matrix(b_blocks,      ncol=1L)
    result$model <- do.call(rclsp::clsp, c(list(problem="ap", b=b, C=dots$C_lim,
                                                S=S, M=M, m=m, p=p, i=i, j=j,
                                                zero_diagonal=zero_diagonal),
                                           dots[setdiff(names(dots),
                                                        c("C_lim", "b_lim"))]))
    result$x     <- result$model$x
    
    # perform reduced estimation and return the result
  } else {
    reduced      <- as.integer(reduced)
    if (length(reduced)  != 2L                                               ||
        isTRUE(reduced[1] < 3L)                                              ||
        isTRUE(reduced[2] < 3L))  stop(paste0("Each reduced block must be at ",
                                              "least (3, 3) to allow a ",
                                              "solvable CLSP submatrix with a ",
                                              "slack (surplus) structure."))
    result       <- list(full=FALSE,                  model=list(),
                         x   =matrix(NA_real_,        nrow=m,       ncol=p))
    m_subset     <- reduced[1] - 1L
    p_subset     <- reduced[2] - 1L
    if (isTRUE(zero_diagonal)) {
      b_diag <-  matrix(0, nrow=min(m, p), ncol=        1L)
      M_diag <-  matrix(0, nrow=min(m, p), ncol=     m * p)
      for (k in  seq_len(min(m, p)))   M_diag[k,   (k - 1L) * p + k]   <- 1
      b_val  <-  if (is.null(b_val) || length(b_val) == 0L) b_diag     else
        rbind(b_val, b_diag)
      M      <-  if (is.null(M)     || length(M)     == 0L) M_diag     else
        rbind(M,      M_diag)
      ord    <- do.call(order, as.data.frame(M))
      idx    <- ord[!duplicated(M[ord, , drop=FALSE])]
      M      <- M[idx, ,                 drop=FALSE]
      b_val  <- b_val[idx,               drop=FALSE]
    }
    if (!is.null(M))           {
      if (!(all(        abs(M    ) <   sqrt(.Machine$double.eps)  |
                        abs(M - 1) <   sqrt(.Machine$double.eps))        &&
            all(rowSums(abs(M - 1) <   sqrt(.Machine$double.eps)) == 1)  &&
            all(colSums(abs(M - 1) <   sqrt(.Machine$double.eps)) <= 1)))
        stop(paste0("M must be a unique row subset ",
                    "of the identity matrix in the ",
                    "reduced model."))
      X_true   <- matrix(NA_real_,     nrow=1L*m, ncol=p)
      for (idx in seq_len(nrow(M))) {
        col          <- which.max(M[idx, ])
        r            <- (col - 1L) %/% p  + 1L         # M has m * p columns
        c            <- (col - 1L) %%  p  + 1L
        X_true[r, c] <- b_val[idx]
      }
    }
    if (!is.null(dots$b_lim))  {
      X_lim    <- matrix(NA_real_,     nrow=2L*m, ncol=p)
      l_idx    <- if (nrow(pos<-which(
        S[(nrow(S)-nrow(dots$b_lim)+1L):nrow(S), ] == -1,
        arr.ind=TRUE)) > 0) pos[1,1] else nrow(dots$b_lim)
      if (l_idx > 1L               ) {
        for (idx in 1:(l_idx - 1L        )) {
          col  <- which.max(dots$C_lim[idx, ])
          r    <- floor((col - 1L) /  p) + 1L
          c    <-      ((col - 1L) %% p) + 1L
          X_lim[  r, c] <-  dots$b_lim[idx  ]
        }
      }
      if (l_idx <= nrow(dots$b_lim)) {
        for (idx in l_idx:nrow(dots$C_lim)) {
          col  <- which.max(dots$C_lim[idx, ])
          r    <- floor((col - 1L) /  p) + 1L
          c    <-      ((col - 1L) %% p) + 1L
          X_lim[m+r, c] <-  dots$b_lim[idx  ]
        }
      }
    }
    if (!is.null(S) && any(S[1:(if (is.null(dots$b_lim)) nrow(S) else nrow(S)-
                                nrow(dots$b_lim)), ] != 0)) {
      warning("User-provided S is ignored in the reduced model.", call.=FALSE)
    }
    for   (row_block in seq_len(ceiling(m / m_subset))) {
      for (col_block in seq_len(ceiling(p / p_subset))) {
        m_start  <- (row_block  - 1L) * m_subset + 1L
        m_end    <- min(m_start + m_subset - 1L,   m)
        p_start  <- (col_block  - 1L) * p_subset + 1L
        p_end    <- min(p_start + p_subset - 1L,   p)
        C_subset <-  NULL
        S_subset <- diag((m_end - m_start  + 1L) +    (p_end - p_start  + 1L))
        b_subset <- matrix(c(b_row[m_start:m_end],
                             b_col[p_start:p_end]),                  ncol=1L)
        M_subset <-  NULL
        if (!is.null(M))      {
          subset     <- matrix(X_true[m_start:m_end, p_start:p_end], ncol=1L)
          non_empty  <- !is.na(subset) & is.finite(subset)
          if (any(non_empty)) {
            M_subset <- diag(length(subset))[non_empty, , drop=FALSE]
            b_subset <- matrix(c(b_subset, subset[non_empty]),       ncol=1L)
          }
        }
        if (!is.null(dots$b_lim))  {
          subset     <- as.vector(t(rbind(
            X_lim[(m_start  ):(m_end  ), p_start:p_end],
            X_lim[(m+m_start):(m+m_end), p_start:p_end]
          )))
          non_empty  <- !is.na(subset) & is.finite(subset)
          if (any(non_empty)) {
            b_subset <- matrix(c(b_subset, subset[non_empty]),       ncol=1L)
            C_subset <- rbind(diag(h <- length(subset) %/% 2),
                              diag(h)      )[non_empty, , drop=FALSE]
            S        <- cbind(rbind(diag(h),
                                    matrix(0, nrow=h, ncol=h))[non_empty, ,
                                                               drop=FALSE],
                              rbind(matrix(0, nrow=h, ncol=h),
                                    -diag(h)                  )[non_empty, ,
                                                                drop=FALSE])
            S_subset <- rbind(cbind(S_subset, matrix(0, nrow=nrow(S_subset),
                                                     ncol=ncol(S) )),
                              cbind(matrix(0, nrow=nrow(S),
                                           ncol=ncol(S_subset)), S))
          }
        }
        tmp      <-  do.call(rclsp::clsp, c(list(problem="ap", b=b_subset,
                                                 C=C_subset,   S=S_subset,
                                                 M=M_subset,
                                                 m=m_end - m_start + 1L,
                                                 p=p_end - p_start + 1L,
                                                 i=i, j=j, zero_diagonal=FALSE),
                                            dots[setdiff(names(dots),
                                                         c("C_lim", "b_lim"))]))
        result$model[[length(result$model) + 1L]] <-  tmp
        result$x[m_start:m_end, p_start:p_end]    <-  tmp$x
      }
    }
  }
  
  # enforce symmetry
  if (isTRUE(symmetric)) if (nrow(result$x) == ncol(result$x))
    result$x    <-       0.5 *   (result$x  +     t(result$x))         else
      stop("symmetric=True requires a square matrix (m == p).")
  
  class(result) <-   "tmpinv"
  result
}
