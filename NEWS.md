# rtmpinvi 2.0.0

## Changes
* Synchronized the package with `rtmpinv` 2.0.0.

# rtmpinvi 1.1.0

## Changes
* Updated the minimum R version to R 4.3 to match CVXR 1.8.x requirements.
* Synchronized the package with `rtmpinv` 1.1.0.

## Bug fixes
* Fixed handling of all-missing `ival` inputs so that `b_val` and `M` are not
  passed to `rtmpinv::tmpinv()` when no usable prior cell information is
  available. `b_val`, `M`, and `bounds` are also no longer passed as explicit
  `NULL` values to `rtmpinv::tmpinv()`.

# rtmpinvi 1.0.0

## Changes
* Synchronized the package with `rtmpinv` 1.0.0.

# rtmpinvi 0.2.0

## Bug fixes
* Updated CVXR integration for CVXR 1.8.x compatibility.
