test_that("tmpinvi works for an interactive AP/TM example", {
  skip_if_not_installed("tmpinv")
  
  set.seed(123)
  
  iso2 <- c("CN", "DE", "JP", "NL", "US")
  T    <- 3L                                           # small for CRAN
  year <- 1:T
  m    <- length(iso2)
  
  df <- expand.grid(year = year, iso2 = iso2, KEEP.OUT.ATTRS = FALSE)
  df <- df[order(df$year, df$iso2), ]
  
  ex_cols <- paste0("EX_", iso2)
  df[ex_cols] <- NA_real_
  df$EX <- NA_real_
  df$IM <- NA_real_
  
  # deterministic stable matrices
  for (t in seq_len(T)) {
    base <- matrix(seq_len(m * m), m, m)
    X <- base * (1 + 0.1 * t)
    diag(X) <- 0
    
    rows <- ((t - 1L) * m + 1L):((t - 1L) * m + m)
    df$EX[rows] <- rowSums(X)
    df$IM[rows] <- colSums(X)
    
    # 50 percent missing (deterministic)
    mask <- matrix((seq_len(m * m) %% 2) == 0, m, m)
    X[mask] <- NA_real_
    
    df[rows, ex_cols] <- X
  }
  
  make_bounds <- function(lb, ub)
    Map(function(a, b) c(a, b), lb, ub)
  
  df_out <- df
  
  for (y in year) {
    idx  <- df_out$year == y
    d    <- df_out[idx, ]
    ival <- as.matrix(d[ex_cols])
    
    # simple deterministic bounds
    lb <- rep(0, m * m)
    ub <- rep(max(d$EX), m * m)
    
    fit <- tmpinvi(
      ival = ival,
      ibounds = make_bounds(lb, ub),
      b_row = d$EX,
      b_col = d$IM,
      alpha = 1.0,
      update = TRUE
    )
    
    expect_s3_class(fit, "tmpinvi")
    expect_true(is.list(fit$result))
    expect_true(is.matrix(fit$data))
    expect_equal(dim(fit$data), c(m, m))
    
    df_out[idx, ex_cols] <- fit$data
  }
})
