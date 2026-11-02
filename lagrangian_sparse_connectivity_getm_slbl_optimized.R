#!/usr/bin/env Rscript

################################################################################
# Lagrangian particle-tracking sparse connectivity for GETM 2D surface/bottom
# files: getm_ref_phy_dm_2d_slbl_YYYYMM.nc4
#
# This script replaces the one-day shifted-cell overlap graph with a particle
# tracker plus sparse transition matrices.
#
# Main idea:
#   1. Read daily velocity slices from NetCDF: ssu/ssv, sbu/sbv.
#   2. Release particles from a pipe cell, specified cell IDs, or all wet cells.
#   3. Move particles with sub-daily numerical steps using bilinear interpolation.
#   4. Convert endpoints back to grid cells.
#   5. Build sparse transition matrices P where P[i, j] is the probability that
#      water released from source cell i reaches destination cell j.
#   6. Write daily sparse edges, a daily cell water-budget CSV, optional matrix
#      RDS files, and PNG maps of source removal and transported destination gain.
#
# Important limitation:
#   The NetCDF has only surface and bottom velocities. The "mixed" mode is a
#   scenario assumption, not resolved vertical physics. Run surface, bottom, and
#   mixed ensembles to represent vertical uncertainty.
################################################################################

# ==============================================================================
# Performance and reproducibility
#   - Default "exact" preserves the original seeded numerical results and output
#     order. Existing scientific parameters below are unchanged.
#   - Each monthly NetCDF is copied ONCE from the input drive to local staging.
#     Required daily fields (all four velocity variables together) are then read
#     once into a shared, on-disk daily store. Workers never read the input drive.
#   - Scenario/month jobs use Windows-compatible PSOCK workers. A single bounded
#     worker pool avoids nested parallelism. Only the coordinator writes outputs.
#   - All jobs use the ORIGINAL global time axis and horizon-specific valid starts.
#     Trajectories may cross months; required following slices are included.
#     There is no change_days statistic in this script. No artificial boundary
#     is introduced, and daily results are reduced in original order (not means
#     of monthly means). Existing treatment of gaps in the input is unchanged.
#   - Velocity memory is an LRU cache of at most velocity_cache_days per worker.
#     Staging and temporary results use local DISK; disk use grows with the run.
#   - Exact horizon reuse is possible only without release jitter and diffusion.
#     Original jitter uses new random positions for EACH horizon. For Kh > 0,
#     draw counts depend on surviving particles, so exact legacy runs remain
#     serial; nondiffusive scenarios on either side still run in parallel.
#   - Opt-in "shared" uses one reproducible stream per scenario/release and one
#     trajectory to the longest valid horizon, with intermediate checkpoints.
#     It parallelizes every Kh/month, but CHANGES legacy stochastic samples and
#     the coupling between horizons. It is NOT an exact-output optimization.
#   - Timing/log/metadata fields necessarily differ; CSV schemas, matrix format,
#     filenames, physical equations and map definitions are preserved.
#
# 0. USER PARAMETERS - edit this block first
# ==============================================================================

# Main folder containing getm_ref_phy_dm_2d_slbl_YYYYMM.nc4 files.
input_dir <- "Z:/Data/MSI_models/EMARIS/rita_clim_test/getm/ref/phy/dm"

# Folder where CSV, RDS, log, and PNG map outputs are written.
output_dir <- "Z:/Data/MSI_models/connectivity model/temp/lagrangian_connectivity_output"

# Normal monthly-file discovery. The supplied inventory covered 1980-01 to 1981-12.
years  <- 1980
months <- 1
file_name_template <- "getm_ref_phy_dm_2d_slbl_%04d%02d.nc4"
product_id <- "getm_ref_phy_dm_2d_slbl"

# Optional explicit file list. Leave character() to use input_dir + template above.
# Example for a single local test file:
# explicit_input_files <- c("getm_ref_phy_dm_2d_slbl_198001.nc4")
explicit_input_files <- character()

# NetCDF variable names.
time_var <- "time"
surface_u_var <- "ssu"
surface_v_var <- "ssv"
bottom_u_var  <- "sbu"
bottom_v_var  <- "sbv"

# Depth/scenario handling. Recommended: at least surface, bottom, and mixed.
# mixed_w = 1 means surface; mixed_w = 0 means bottom; mixed_w = 0.5 is simple average.
velocity_modes <- c("surface", "bottom", "mixed")
mixed_w_values <- c(0.5)

# Horizontal random-walk diffusion coefficients in m2/s.
# If turbulence is unknown, run sensitivity tests, e.g. c(0, 1, 5, 10).
Kh_values <- c(0)

# Particle tracking settings.
dt_seconds <- 3600                 # hourly substeps recommended for daily fields
tracking_horizons_days <- c(1, 3, 5, 7)
N_particles_per_source <- 100L
release_jitter_within_cell <- TRUE # TRUE releases uniformly inside the source cell
random_seed <- 1234L
R_earth_m <- 6371000

# Source selection.
#   "pipe"          = nearest wet cell to pipe_lon/pipe_lat.
#   "cell_ids"      = use source_cell_ids below.
#   "all_wet_cells" = build transitions from every wet cell; expensive.
source_selection <- "all_wet_cells"
pipe_lon <- 24.75
pipe_lat <- 59.45
source_cell_ids <- integer()

# Start-day controls. NULL means all possible release days for each horizon.
# For a quick pilot, set max_start_days <- 2L or 3L.
start_day_stride <- 1L
max_start_days <- 2L
start_global_indices <- NULL       # optional integer vector overriding stride/max_start_days

# Optional spatial subset in lon/lat degrees, applied to source cells only.
# Leave as NA to use all selected source cells.
source_subset_lonlat <- list(lon_min = NA_real_, lon_max = NA_real_,
                             lat_min = NA_real_, lat_max = NA_real_)

# Output controls.
# TRUE lets you rerun after a failed/partial run without manually deleting outputs.
# Set to FALSE when you want to protect completed outputs from replacement.
overwrite_existing <- TRUE
output_prefix <- NULL
save_daily_sparse_edges_csv <- TRUE
save_daily_cell_budget_csv <- TRUE
save_daily_transition_matrices <- TRUE
save_mean_transition_matrices <- TRUE
# TRUE keeps matrix subfolder/file names short. This avoids Windows path-length
# failures such as gzfile()/saveRDS() saying "No such file or directory"
# when the full RDS path exceeds the legacy 260-character limit.
short_matrix_file_names <- TRUE
save_mean_sparse_edges_csv <- TRUE
save_particle_endpoints <- FALSE   # can become large
save_cell_lookup <- TRUE
save_metadata <- TRUE
save_day_qc_csv <- TRUE
save_maps <- TRUE

# The cell-budget CSV normally includes only source cells and cells that receive
# transported particles. TRUE writes all wet cells for every day/scenario/horizon.
budget_include_all_wet_cells <- FALSE

# Sparse edge CSV threshold. Matrix RDS files are not thresholded.
min_edge_probability_to_write <- 0

# Maps: "all" or a numeric vector matching tracking_horizons_days, e.g. c(1, 7).
map_horizons_days <- "all"
map_width_px <- 2200L
map_height_px <- 1800L
map_res_dpi <- 220L
map_point_cex <- 0.20
map_na_colour <- "#EEEEEE"
map_low_high_colours <- c("#F5F5F5", "#111111")

# Coordinate and dimension detection candidates.
lon_coord_candidates <- c("lonc", "lon", "longitude", "clon", "nav_lon", "xlon")
lat_coord_candidates <- c("latc", "lat", "latitude", "clat", "nav_lat", "ylat")
lon_dim_candidates <- c("lonc", "lon", "longitude", "x")
lat_dim_candidates <- c("latc", "lat", "latitude", "y")

# Numerical controls.
fill_value_cutoff <- -9000
interp_requires_all_four_corners <- TRUE
inactive_particles_are_lost <- TRUE
compact_edge_accumulator_every <- 10L

# Performance controls. No extra R packages are needed (parallel ships with R).
trajectory_reuse <- "shared"        # "exact" (same output) or opt-in "shared"
parallel_workers <- NULL            # NULL: up to 4 cores; 1: serial execution
parallel_split <- "month"          # "month" or "none" (scenario-only jobs)
velocity_cache_days <- NULL         # NULL: max horizon + 1; minimum 2
# Must be on a LOCAL disk with room for inputs, daily fields and pending results.
# Example on Windows: staging_dir <- "E:/RActions/connectivity_temp"
staging_dir <- tempdir()
keep_temporary_files <- FALSE       # TRUE retains this run's unique temp folder

# Thread control for data.table. NULL leaves the coordinator at package default
# and uses one data.table thread per worker to avoid nested parallelism.
data_table_threads <- NULL

# ==============================================================================
# 1. PACKAGE SETUP
# ==============================================================================

required_packages <- c("ncdf4", "data.table", "Matrix")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                              logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "),
       ". Install them before running this script.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(ncdf4)
  library(data.table)
  library(Matrix)
})

if (!is.null(data_table_threads)) {
  data.table::setDTthreads(data_table_threads)
}

# ==============================================================================
# 2. SMALL UTILITIES
# ==============================================================================

log_lines <- character()

log_message <- function(..., level = "INFO") {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " [", level, "] ",
                paste(..., collapse = ""))
  message(msg)
  log_lines <<- c(log_lines, msg)
  invisible(msg)
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

sanitize_label <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

ensure_dir <- function(path) {
  if (is.null(path) || length(path) == 0L || is.na(path) || !nzchar(path)) {
    stop("Invalid directory path.", call. = FALSE)
  }
  if (!dir.exists(path)) {
    ok <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!isTRUE(ok) && !dir.exists(path)) {
      stop("Could not create directory: ", path, call. = FALSE)
    }
  }
  if (!dir.exists(path)) {
    stop("Directory does not exist after creation attempt: ", path, call. = FALSE)
  }
  invisible(path)
}

ensure_parent_dir <- function(path) {
  parent <- dirname(path)
  ensure_dir(parent)
  invisible(parent)
}

assert_can_write <- function(path, overwrite) {
  ensure_parent_dir(path)
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop("Output file already exists and overwrite_existing is FALSE: ", path,
         call. = FALSE)
  }
  invisible(TRUE)
}

prepare_output_file <- function(path, overwrite) {
  assert_can_write(path, overwrite)
  if (file.exists(path) && isTRUE(overwrite)) unlink(path)
  invisible(path)
}

safe_save_rds <- function(object, path, overwrite, compress = TRUE) {
  prepare_output_file(path, overwrite)
  tryCatch(
    {
      saveRDS(object, path, compress = compress)
    },
    error = function(e) {
      msg <- paste0("Could not write RDS file:\n  ", path,
                    "\nParent directory exists: ", dir.exists(dirname(path)),
                    "\nPath length: ", nchar(path, type = "chars"),
                    " characters. Original error: ", conditionMessage(e))
      stop(msg, call. = FALSE)
    }
  )
  if (!file.exists(path)) {
    stop("saveRDS completed but the file was not found afterwards: ", path,
         call. = FALSE)
  }
  invisible(path)
}


append_fwrite <- function(dt, path) {
  ensure_parent_dir(path)
  data.table::fwrite(dt, path, append = file.exists(path),
                     col.names = !file.exists(path))
  invisible(path)
}

short_run_label <- function(params) {
  yr_label <- if (length(params$years) == 1L) {
    as.character(params$years)
  } else {
    paste0(min(params$years), "-", max(params$years))
  }
  mo_label <- paste(sprintf("%02d", sort(unique(params$months))), collapse = "-")
  src_label <- switch(params$source_selection,
                      all_wet_cells = "allwet",
                      cell_ids = "cells",
                      pipe = "pipe",
                      params$source_selection)
  sanitize_label(paste0("y", yr_label, "_m", mo_label, "_", src_label))
}

make_daily_matrix_dir <- function(params) {
  if (isTRUE(params$short_matrix_file_names)) {
    file.path(params$output_dir, paste0("P_daily_", short_run_label(params)))
  } else {
    file.path(params$output_dir, paste0(params$output_prefix, "_daily_transition_matrices"))
  }
}

make_mean_matrix_dir <- function(params) {
  if (isTRUE(params$short_matrix_file_names)) {
    file.path(params$output_dir, paste0("P_mean_", short_run_label(params)))
  } else {
    file.path(params$output_dir, paste0(params$output_prefix, "_mean_transition_matrices"))
  }
}

make_daily_matrix_path <- function(daily_matrix_dir, matrix_key, start_global_index, params) {
  if (isTRUE(params$short_matrix_file_names)) {
    file.path(daily_matrix_dir, paste0("P_", matrix_key, "_g", start_global_index, ".rds"))
  } else {
    file.path(daily_matrix_dir,
              paste0(params$output_prefix, "_", matrix_key,
                     "_g", start_global_index, "_P.rds"))
  }
}

make_mean_matrix_path <- function(mean_matrix_dir, key, params) {
  if (isTRUE(params$short_matrix_file_names)) {
    file.path(mean_matrix_dir, paste0("Pmean_", key, ".rds"))
  } else {
    file.path(mean_matrix_dir, paste0(params$output_prefix, "_", key, "_P_mean.rds"))
  }
}

coord_edges <- function(x) {
  x <- as.numeric(x)
  if (length(x) < 2L) stop("At least two coordinate values are needed.", call. = FALSE)
  if (!all(is.finite(x))) stop("Coordinate vector contains non-finite values.", call. = FALSE)
  if (!all(diff(x) > 0)) stop("Coordinate vector must be strictly increasing.", call. = FALSE)
  mids <- (x[-1L] + x[-length(x)]) / 2
  c(x[1L] - (mids[1L] - x[1L]), mids,
    x[length(x)] + (x[length(x)] - mids[length(mids)]))
}

plausible_lonlat <- function(lon, lat) {
  lon <- as.numeric(lon)
  lat <- as.numeric(lat)
  if (!any(is.finite(lon)) || !any(is.finite(lat))) return(FALSE)
  lon_rng <- range(lon, finite = TRUE)
  lat_rng <- range(lat, finite = TRUE)
  isTRUE(lon_rng[1L] >= -360 && lon_rng[2L] <= 360 &&
         lat_rng[1L] >= -90 && lat_rng[2L] <= 90)
}

parse_cf_time <- function(values, units) {
  values <- as.numeric(values)
  out <- rep(as.POSIXct(NA_real_, origin = "1970-01-01", tz = "UTC"), length(values))
  if (is.null(units) || length(units) == 0L || is.na(units)) return(out)
  units <- as.character(units)[1L]
  m <- regexec("^([A-Za-z]+) since ([0-9]{4}-[0-9]{2}-[0-9]{2})(?:[ T]([0-9:]+))?", units)
  parts <- regmatches(units, m)[[1L]]
  if (length(parts) < 3L) return(out)

  unit_name <- tolower(parts[2L])
  origin_time <- if (length(parts) >= 4L && nzchar(parts[4L])) parts[4L] else "00:00:00"
  origin <- as.POSIXct(paste(parts[3L], origin_time), tz = "UTC")
  if (is.na(origin)) return(out)

  multiplier <- switch(unit_name,
                       second = 1, seconds = 1, sec = 1, secs = 1,
                       minute = 60, minutes = 60, min = 60, mins = 60,
                       hour = 3600, hours = 3600,
                       day = 86400, days = 86400,
                       NA_real_)
  if (!is.finite(multiplier)) return(out)
  origin + values * multiplier
}

format_time_safe <- function(x) {
  if (inherits(x, "POSIXct")) {
    out <- format(x, "%Y-%m-%d %H:%M:%S", tz = "UTC")
    out[is.na(x)] <- NA_character_
    return(out)
  }
  as.character(x)
}

make_colour_vector <- function(values, params) {
  cols <- rep(params$map_na_colour, length(values))
  finite <- is.finite(values)
  if (!any(finite)) return(cols)

  vals <- values[finite]
  vmin <- min(vals, na.rm = TRUE)
  vmax <- max(vals, na.rm = TRUE)
  ramp <- grDevices::colorRampPalette(params$map_low_high_colours)(100L)
  if (!is.finite(vmin) || !is.finite(vmax) || vmax <= vmin) {
    idx <- if (is.finite(vmax) && vmax > 0) rep(100L, length(vals)) else rep(1L, length(vals))
  } else {
    idx <- floor(1 + 99 * (vals - vmin) / (vmax - vmin))
    idx <- pmin(100L, pmax(1L, idx))
  }
  cols[finite] <- ramp[idx]
  cols
}

# ==============================================================================
# 3. PARAMETER COLLECTION AND FILE DISCOVERY
# ==============================================================================

collect_params <- function() {
  params <- list(
    input_dir = input_dir,
    output_dir = output_dir,
    years = as.integer(years),
    months = as.integer(months),
    file_name_template = file_name_template,
    product_id = product_id,
    explicit_input_files = explicit_input_files,
    time_var = time_var,
    surface_u_var = surface_u_var,
    surface_v_var = surface_v_var,
    bottom_u_var = bottom_u_var,
    bottom_v_var = bottom_v_var,
    velocity_modes = velocity_modes,
    mixed_w_values = as.numeric(mixed_w_values),
    Kh_values = as.numeric(Kh_values),
    dt_seconds = as.numeric(dt_seconds),
    tracking_horizons_days = as.integer(tracking_horizons_days),
    N_particles_per_source = as.integer(N_particles_per_source),
    release_jitter_within_cell = release_jitter_within_cell,
    random_seed = random_seed,
    R_earth_m = as.numeric(R_earth_m),
    source_selection = source_selection,
    pipe_lon = as.numeric(pipe_lon),
    pipe_lat = as.numeric(pipe_lat),
    source_cell_ids = as.integer(source_cell_ids),
    start_day_stride = as.integer(start_day_stride),
    max_start_days = max_start_days,
    start_global_indices = start_global_indices,
    source_subset_lonlat = source_subset_lonlat,
    overwrite_existing = overwrite_existing,
    output_prefix = output_prefix,
    save_daily_sparse_edges_csv = save_daily_sparse_edges_csv,
    save_daily_cell_budget_csv = save_daily_cell_budget_csv,
    save_daily_transition_matrices = save_daily_transition_matrices,
    save_mean_transition_matrices = save_mean_transition_matrices,
    short_matrix_file_names = short_matrix_file_names,
    save_mean_sparse_edges_csv = save_mean_sparse_edges_csv,
    save_particle_endpoints = save_particle_endpoints,
    save_cell_lookup = save_cell_lookup,
    save_metadata = save_metadata,
    save_day_qc_csv = save_day_qc_csv,
    save_maps = save_maps,
    budget_include_all_wet_cells = budget_include_all_wet_cells,
    min_edge_probability_to_write = as.numeric(min_edge_probability_to_write),
    map_horizons_days = map_horizons_days,
    map_width_px = as.integer(map_width_px),
    map_height_px = as.integer(map_height_px),
    map_res_dpi = as.integer(map_res_dpi),
    map_point_cex = map_point_cex,
    map_na_colour = map_na_colour,
    map_low_high_colours = map_low_high_colours,
    lon_coord_candidates = lon_coord_candidates,
    lat_coord_candidates = lat_coord_candidates,
    lon_dim_candidates = lon_dim_candidates,
    lat_dim_candidates = lat_dim_candidates,
    fill_value_cutoff = as.numeric(fill_value_cutoff),
    interp_requires_all_four_corners = interp_requires_all_four_corners,
    inactive_particles_are_lost = inactive_particles_are_lost,
    compact_edge_accumulator_every = as.integer(compact_edge_accumulator_every),
    data_table_threads = data_table_threads,
    trajectory_reuse = trajectory_reuse,
    parallel_workers = parallel_workers,
    parallel_split = parallel_split,
    velocity_cache_days = velocity_cache_days,
    staging_dir = staging_dir,
    keep_temporary_files = keep_temporary_files
  )

  if (length(params$velocity_modes) == 0L ||
      any(!params$velocity_modes %in% c("surface", "bottom", "mixed"))) {
    stop("velocity_modes must contain any of: surface, bottom, mixed.", call. = FALSE)
  }
  if (length(params$tracking_horizons_days) == 0L ||
      any(!is.finite(params$tracking_horizons_days)) || any(params$tracking_horizons_days < 1L)) {
    stop("tracking_horizons_days must contain positive integers.", call. = FALSE)
  }
  if (!is.finite(params$dt_seconds) || params$dt_seconds <= 0) {
    stop("dt_seconds must be positive.", call. = FALSE)
  }
  steps_per_day <- 86400 / params$dt_seconds
  if (abs(steps_per_day - round(steps_per_day)) > 1e-9) {
    stop("dt_seconds must divide 86400 exactly for this script. Example: 1800 or 3600.",
         call. = FALSE)
  }
  if (!is.finite(params$N_particles_per_source) || params$N_particles_per_source < 1L) {
    stop("N_particles_per_source must be a positive integer.", call. = FALSE)
  }
  if (!params$source_selection %in% c("pipe", "cell_ids", "all_wet_cells")) {
    stop("source_selection must be 'pipe', 'cell_ids', or 'all_wet_cells'.", call. = FALSE)
  }
  if (params$source_selection == "cell_ids" && length(params$source_cell_ids) == 0L) {
    stop("source_selection is 'cell_ids' but source_cell_ids is empty.", call. = FALSE)
  }
  if (!is.finite(params$start_day_stride) || params$start_day_stride < 1L) {
    stop("start_day_stride must be a positive integer.", call. = FALSE)
  }
  if (!is.finite(params$min_edge_probability_to_write) || params$min_edge_probability_to_write < 0) {
    stop("min_edge_probability_to_write must be non-negative.", call. = FALSE)
  }

  if (!params$trajectory_reuse %in% c("exact", "shared")) {
    stop("trajectory_reuse must be 'exact' or 'shared'.", call. = FALSE)
  }
  if (!params$parallel_split %in% c("month", "none")) {
    stop("parallel_split must be 'month' or 'none'.", call. = FALSE)
  }
  positive_integer <- function(x) length(x) == 1L && is.finite(x) && x >= 1 && x == floor(x)
  if (is.null(params$parallel_workers)) {
    cores <- parallel::detectCores(logical = FALSE)
    if (!is.finite(cores)) cores <- 1L
    params$parallel_workers <- as.integer(max(1L, min(4L, cores - 1L)))
  }
  if (!positive_integer(params$parallel_workers)) {
    stop("parallel_workers must be NULL or a positive integer.", call. = FALSE)
  }
  if (is.null(params$velocity_cache_days)) {
    params$velocity_cache_days <- max(params$tracking_horizons_days) + 1L
  }
  if (!positive_integer(params$velocity_cache_days) || params$velocity_cache_days < 2L) {
    stop("velocity_cache_days must be NULL or an integer >= 2.", call. = FALSE)
  }
  if (anyDuplicated(params$tracking_horizons_days)) {
    stop("tracking_horizons_days must be unique to avoid duplicate output paths.", call. = FALSE)
  }

  if (is.null(params$output_prefix)) {
    yr_label <- if (length(params$years) == 1L) as.character(params$years) else paste0(min(params$years), "-", max(params$years))
    mo_label <- paste(sprintf("%02d", sort(unique(params$months))), collapse = "-")
    src_label <- params$source_selection
    params$output_prefix <- sanitize_label(paste0("lagrangian_connectivity_", params$product_id,
                                                  "_", src_label, "_y", yr_label,
                                                  "_m", mo_label))
  } else {
    params$output_prefix <- sanitize_label(params$output_prefix)
  }

  params
}

discover_files <- function(params) {
  if (length(params$explicit_input_files) > 0L) {
    dt <- data.table(file_order = seq_along(params$explicit_input_files),
                     full_path = normalizePath(params$explicit_input_files, mustWork = FALSE))
    dt[, exists := file.exists(full_path)]
    if (any(!dt$exists)) {
      missing <- dt[exists == FALSE, full_path]
      stop("Missing explicit NetCDF files:\n", paste(missing, collapse = "\n"), call. = FALSE)
    }
    dt[, `:=`(year = NA_integer_, month = NA_integer_, file_name = basename(full_path))]
    return(dt[])
  }

  request <- data.table(year = rep(params$years, each = length(params$months)),
                        month = rep(params$months, times = length(params$years)))
  setorder(request, year, month)
  request[, file_order := .I]
  request[, file_name := sprintf(params$file_name_template, year, month)]
  request[, full_path := file.path(params$input_dir, file_name)]
  request[, exists := file.exists(full_path)]

  if (any(!request$exists)) {
    missing <- request[exists == FALSE, paste0(year, "-", sprintf("%02d", month),
                                               " -> ", full_path)]
    stop("Missing requested NetCDF files:\n", paste(missing, collapse = "\n"),
         call. = FALSE)
  }

  request[]
}

# ==============================================================================
# 4. NETCDF GRID AND VARIABLE READERS
# ==============================================================================

get_var_dim_names <- function(nc, var_name) {
  if (!var_name %in% names(nc$var)) {
    stop("Variable not found in NetCDF: ", var_name, call. = FALSE)
  }
  vapply(nc$var[[var_name]]$dim, function(d) d$name, character(1))
}

get_time_length <- function(nc, time_name, fallback_var = NULL) {
  if (time_name %in% names(nc$dim)) return(nc$dim[[time_name]]$len)
  if (!is.null(fallback_var) && fallback_var %in% names(nc$var)) {
    dnames <- get_var_dim_names(nc, fallback_var)
    pos <- match(time_name, dnames)
    if (!is.na(pos)) return(nc$var[[fallback_var]]$dim[[pos]]$len)
    pos <- grep("time", dnames, ignore.case = TRUE)[1L]
    if (!is.na(pos)) return(nc$var[[fallback_var]]$dim[[pos]]$len)
  }
  stop("Could not determine time dimension length.", call. = FALSE)
}

choose_spatial_dims <- function(nc, params) {
  u_dims <- get_var_dim_names(nc, params$surface_u_var)
  time_dim <- if (params$time_var %in% u_dims) {
    params$time_var
  } else {
    td <- grep("time", u_dims, ignore.case = TRUE, value = TRUE)
    if (length(td) == 0L) stop("Could not identify time dimension in ", params$surface_u_var,
                               call. = FALSE)
    td[1L]
  }

  spatial_dims <- setdiff(u_dims, time_dim)
  lon_dim <- intersect(params$lon_dim_candidates, spatial_dims)[1L]
  lat_dim <- intersect(params$lat_dim_candidates, spatial_dims)[1L]

  if (is.na(lon_dim) || length(lon_dim) == 0L) {
    lon_dim <- grep("lon|x", spatial_dims, ignore.case = TRUE, value = TRUE)[1L]
  }
  if (is.na(lat_dim) || length(lat_dim) == 0L) {
    lat_dim <- grep("lat|y", spatial_dims, ignore.case = TRUE, value = TRUE)[1L]
  }
  if (is.na(lon_dim) || is.na(lat_dim)) {
    stop("Could not identify longitude and latitude dimensions from ",
         params$surface_u_var, " dimensions: ", paste(u_dims, collapse = ", "),
         call. = FALSE)
  }

  list(lon_dim = lon_dim, lat_dim = lat_dim, time_dim = time_dim,
       nlon = nc$dim[[lon_dim]]$len, nlat = nc$dim[[lat_dim]]$len)
}

read_coord_1d <- function(nc, dim_name, candidates, expected_len) {
  dim_vals <- nc$dim[[dim_name]]$vals
  if (length(dim_vals) == expected_len && any(is.finite(dim_vals))) return(as.numeric(dim_vals))

  for (nm in unique(c(dim_name, candidates))) {
    if (!nm %in% names(nc$var)) next
    vals <- drop(ncdf4::ncvar_get(nc, nm))
    if (length(dim(vals)) <= 1L && length(vals) == expected_len && any(is.finite(vals))) {
      return(as.numeric(vals))
    }
  }
  stop("Could not read 1D coordinate values for dimension: ", dim_name, call. = FALSE)
}

read_daily_var <- function(nc, var_name, time_index, grid_def, params) {
  if (!var_name %in% names(nc$var)) stop("Variable not found: ", var_name, call. = FALSE)

  var <- nc$var[[var_name]]
  dnames <- vapply(var$dim, function(d) d$name, character(1))
  time_pos <- match(grid_def$time_dim, dnames)
  if (is.na(time_pos)) time_pos <- match(params$time_var, dnames)
  if (is.na(time_pos)) time_pos <- grep("time", dnames, ignore.case = TRUE)[1L]
  if (is.na(time_pos)) stop("Could not identify time dimension for variable ", var_name,
                            call. = FALSE)

  start <- rep(1L, length(dnames))
  count <- vapply(var$dim, function(d) as.integer(d$len), integer(1))
  start[time_pos] <- as.integer(time_index)
  count[time_pos] <- 1L

  arr <- ncdf4::ncvar_get(nc, var_name, start = start, count = count)
  arr <- drop(arr)
  remaining_dims <- dnames[-time_pos]

  if (length(dim(arr)) != 2L) {
    stop("Expected a 2D daily slice for ", var_name, ", got dimensions: ",
         paste(dim(arr), collapse = " x "), call. = FALSE)
  }

  if (identical(remaining_dims, c(grid_def$lat_dim, grid_def$lon_dim))) {
    mat <- arr
  } else if (identical(remaining_dims, c(grid_def$lon_dim, grid_def$lat_dim))) {
    mat <- t(arr)
  } else if (all(dim(arr) == c(grid_def$nlat, grid_def$nlon))) {
    mat <- arr
  } else if (all(dim(arr) == c(grid_def$nlon, grid_def$nlat))) {
    mat <- t(arr)
  } else {
    stop("Could not orient daily slice for ", var_name, ". Variable dimensions are ",
         paste(dnames, collapse = ", "), ".", call. = FALSE)
  }

  storage.mode(mat) <- "double"
  mat[mat <= params$fill_value_cutoff] <- NA_real_
  mat
}

prepare_grid_info <- function(first_file, params) {
  log_message("Preparing grid from ", basename(first_file), ".")
  nc <- ncdf4::nc_open(first_file)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  required_vars <- c(params$surface_u_var, params$surface_v_var,
                     params$bottom_u_var, params$bottom_v_var)
  missing_vars <- setdiff(required_vars, names(nc$var))
  if (length(missing_vars) > 0L) {
    stop("Required variables missing from NetCDF: ", paste(missing_vars, collapse = ", "),
         call. = FALSE)
  }

  grid_def <- choose_spatial_dims(nc, params)
  lon <- read_coord_1d(nc, grid_def$lon_dim, params$lon_coord_candidates, grid_def$nlon)
  lat <- read_coord_1d(nc, grid_def$lat_dim, params$lat_coord_candidates, grid_def$nlat)

  if (!plausible_lonlat(lon, lat)) {
    stop("Detected coordinates are not plausible WGS84 lon/lat values.", call. = FALSE)
  }
  if (!all(diff(lon) > 0) || !all(diff(lat) > 0)) {
    stop("This script expects increasing 1D lon/lat coordinates.", call. = FALSE)
  }

  lon_edges <- coord_edges(lon)
  lat_edges <- coord_edges(lat)

  n_time <- get_time_length(nc, params$time_var, fallback_var = params$surface_u_var)
  wet_mask <- NULL
  wet_slice <- NA_integer_
  for (tt in seq_len(n_time)) {
    su <- read_daily_var(nc, params$surface_u_var, tt, grid_def, params)
    sv <- read_daily_var(nc, params$surface_v_var, tt, grid_def, params)
    bu <- read_daily_var(nc, params$bottom_u_var,  tt, grid_def, params)
    bv <- read_daily_var(nc, params$bottom_v_var,  tt, grid_def, params)
    candidate <- (is.finite(su) & is.finite(sv)) | (is.finite(bu) & is.finite(bv))
    if (any(candidate)) {
      wet_mask <- candidate
      wet_slice <- tt
      break
    }
  }
  if (is.null(wet_mask) || !any(wet_mask)) {
    stop("Could not infer any wet cells from finite velocity values.", call. = FALSE)
  }

  grid <- data.table(lon_index = rep(seq_len(grid_def$nlon), times = grid_def$nlat),
                     lat_index = rep(seq_len(grid_def$nlat), each = grid_def$nlon))
  grid[, cell_id := (lat_index - 1L) * grid_def$nlon + lon_index]
  grid[, `:=`(lon = lon[lon_index], lat = lat[lat_index],
              lon_min = lon_edges[lon_index], lon_max = lon_edges[lon_index + 1L],
              lat_min = lat_edges[lat_index], lat_max = lat_edges[lat_index + 1L])]
  grid[, wet := wet_mask[cbind(lat_index, lon_index)]]
  wet_grid <- grid[wet == TRUE]
  setkey(wet_grid, cell_id)

  wet_by_cell_id <- rep(FALSE, grid_def$nlat * grid_def$nlon)
  wet_by_cell_id[grid$cell_id] <- grid$wet

  log_message("Detected dimensions: time x ", grid_def$lat_dim, "(", grid_def$nlat,
              ") x ", grid_def$lon_dim, "(", grid_def$nlon, ").")
  log_message("Wet cells: ", format(nrow(wet_grid), big.mark = ","),
              " of ", format(nrow(grid), big.mark = ","),
              "; wet-mask source time index: ", wet_slice, ".")
  log_message("Lon range: ", paste(round(range(lon), 5), collapse = " to "),
              "; lat range: ", paste(round(range(lat), 5), collapse = " to "), ".")

  list(grid_def = grid_def,
       lon = lon,
       lat = lat,
       lon_edges = lon_edges,
       lat_edges = lat_edges,
       grid = grid,
       wet_grid = wet_grid,
       wet_mask = wet_mask,
       wet_by_cell_id = wet_by_cell_id,
       wet_mask_time_index = wet_slice,
       n_cells = grid_def$nlat * grid_def$nlon)
}

build_slice_table <- function(files, grid_def, params) {
  rows <- vector("list", nrow(files))
  global_counter <- 0L

  for (ii in seq_len(nrow(files))) {
    f <- files$full_path[ii]
    nc <- ncdf4::nc_open(f)
    n_time <- get_time_length(nc, params$time_var, fallback_var = params$surface_u_var)

    time_values <- rep(NA_real_, n_time)
    time_units <- NA_character_
    if (params$time_var %in% names(nc$var)) {
      time_values <- as.numeric(ncdf4::ncvar_get(nc, params$time_var))
      time_units <- ncdf4::ncatt_get(nc, params$time_var, "units")$value %||% NA_character_
    } else if (params$time_var %in% names(nc$dim)) {
      time_values <- as.numeric(nc$dim[[params$time_var]]$vals)
      time_units <- ncdf4::ncatt_get(nc, params$time_var, "units")$value %||% NA_character_
    }
    parsed_time <- parse_cf_time(time_values, time_units)

    rows[[ii]] <- data.table(global_index = global_counter + seq_len(n_time),
                             file_order = files$file_order[ii],
                             file_name = basename(f),
                             full_path = f,
                             year = files$year[ii],
                             month = files$month[ii],
                             time_index = seq_len(n_time),
                             time_value = time_values,
                             time_units = time_units,
                             velocity_time_utc = format_time_safe(parsed_time))
    global_counter <- global_counter + n_time
    ncdf4::nc_close(nc)
  }

  out <- rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (nrow(out) < 2L) {
    stop("At least two time slices are needed for temporal interpolation.", call. = FALSE)
  }
  out[]
}

# ==============================================================================
# 5. VELOCITY READING, INTERPOLATION, AND PARTICLE TRACKING
# ==============================================================================

# The cache stores composed u/v for one job; raw fields stay on local disk.
new_velocity_cache <- function(params) {
  cache <- new.env(parent = emptyenv())
  cache$values <- new.env(parent = emptyenv())
  cache$order <- character()
  cache$limit <- as.integer(params$velocity_cache_days)
  cache$peak_days <- 0L
  cache$reads <- 0L
  cache
}

read_velocity_slice <- function(global_idx, scenario, slice_table, grid_def, params, cache_env) {
  key <- paste(global_idx, scenario$mode, format(scenario$w, digits = 17), sep = "_")
  if (exists(key, envir = cache_env$values, inherits = FALSE)) {
    cache_env$order <- c(cache_env$order[cache_env$order != key], key)
    return(get(key, envir = cache_env$values, inherits = FALSE))
  }
  while (length(cache_env$order) >= cache_env$limit) {
    rm(list = cache_env$order[1L], envir = cache_env$values)
    cache_env$order <- cache_env$order[-1L]
  }
  path <- file.path(params$velocity_store, paste0("g", global_idx, ".rds"))
  if (!file.exists(path)) stop("Missing staged velocity slice: ", global_idx, call. = FALSE)
  raw <- readRDS(path)
  if (scenario$mode == "surface") {
    out <- list(u = raw$su, v = raw$sv)
  } else if (scenario$mode == "bottom") {
    out <- list(u = raw$bu, v = raw$bv)
  } else if (scenario$mode == "mixed") {
    out <- list(u = scenario$w * raw$su + (1 - scenario$w) * raw$bu,
                v = scenario$w * raw$sv + (1 - scenario$w) * raw$bv)
  } else {
    stop("Unknown velocity mode: ", scenario$mode, call. = FALSE)
  }
  assign(key, out, envir = cache_env$values)
  cache_env$order <- c(cache_env$order, key)
  cache_env$peak_days <- max(cache_env$peak_days, length(cache_env$order))
  cache_env$reads <- cache_env$reads + 1L
  out
}

interp_grid <- function(mat, xlon, ylat, lon, lat, require_all_four = TRUE) {
  j <- findInterval(xlon, lon)
  i <- findInterval(ylat, lat)

  ok <- is.finite(xlon) & is.finite(ylat) &
    j >= 1L & j < length(lon) & i >= 1L & i < length(lat)

  out <- rep(NA_real_, length(xlon))
  if (!any(ok)) return(out)

  jo <- j[ok]
  io <- i[ok]

  x1 <- lon[jo]
  x2 <- lon[jo + 1L]
  y1 <- lat[io]
  y2 <- lat[io + 1L]

  wx <- (xlon[ok] - x1) / (x2 - x1)
  wy <- (ylat[ok] - y1) / (y2 - y1)

  q11 <- mat[cbind(io,      jo)]
  q21 <- mat[cbind(io,      jo + 1L)]
  q12 <- mat[cbind(io + 1L, jo)]
  q22 <- mat[cbind(io + 1L, jo + 1L)]

  if (isTRUE(require_all_four)) {
    good <- is.finite(q11) & is.finite(q21) & is.finite(q12) & is.finite(q22)
    val <- rep(NA_real_, length(q11))
    val[good] <- (1 - wx[good]) * (1 - wy[good]) * q11[good] +
      wx[good]       * (1 - wy[good]) * q21[good] +
      (1 - wx[good]) * wy[good]       * q12[good] +
      wx[good]       * wy[good]       * q22[good]
    out[ok] <- val
  } else {
    # Renormalized finite-corner interpolation. Useful near coasts, but less conservative.
    w11 <- (1 - wx) * (1 - wy)
    w21 <- wx       * (1 - wy)
    w12 <- (1 - wx) * wy
    w22 <- wx       * wy
    weights <- cbind(w11, w21, w12, w22)
    vals <- cbind(q11, q21, q12, q22)
    finite <- is.finite(vals)
    denom <- rowSums(weights * finite)
    num <- rowSums(weights * ifelse(finite, vals, 0))
    val <- ifelse(denom > 0, num / denom, NA_real_)
    out[ok] <- val
  }

  out
}

cell_id_from_lonlat <- function(xlon, ylat, grid_info) {
  j <- findInterval(xlon, grid_info$lon_edges, rightmost.closed = TRUE)
  i <- findInterval(ylat, grid_info$lat_edges, rightmost.closed = TRUE)

  ok <- is.finite(xlon) & is.finite(ylat) &
    j >= 1L & j <= grid_info$grid_def$nlon &
    i >= 1L & i <= grid_info$grid_def$nlat

  id <- rep(NA_integer_, length(xlon))
  id[ok] <- (i[ok] - 1L) * grid_info$grid_def$nlon + j[ok]

  in_wet <- rep(FALSE, length(id))
  valid_id <- !is.na(id)
  in_wet[valid_id] <- grid_info$wet_by_cell_id[id[valid_id]]
  id[!in_wet] <- NA_integer_
  id
}

release_particles <- function(source_dt, grid_info, params) {
  N <- params$N_particles_per_source
  n_sources <- nrow(source_dt)
  source_ids <- rep(source_dt$cell_id, each = N)

  if (isTRUE(params$release_jitter_within_cell)) {
    lon0 <- stats::runif(n_sources * N,
                         min = rep(source_dt$lon_min, each = N),
                         max = rep(source_dt$lon_max, each = N))
    lat0 <- stats::runif(n_sources * N,
                         min = rep(source_dt$lat_min, each = N),
                         max = rep(source_dt$lat_max, each = N))
  } else {
    lon0 <- rep(source_dt$lon, each = N)
    lat0 <- rep(source_dt$lat, each = N)
  }

  list(source_cell_id = source_ids,
       lon = lon0,
       lat = lat0,
       active = rep(TRUE, n_sources * N))
}

advect_particle_checkpoints <- function(source_dt, start_global_index, horizons,
                                        scenario, slice_table, grid_info, params,
                                        cache_env, on_snapshot) {
  horizons <- sort(unique(as.integer(horizons)))
  particles <- release_particles(source_dt, grid_info, params)
  lonp <- particles$lon
  latp <- particles$lat
  active <- particles$active
  source_ids <- particles$source_cell_id
  steps_per_day <- as.integer(round(86400 / params$dt_seconds))
  checkpoint_steps <- horizons * steps_per_day
  n_steps <- max(checkpoint_steps)
  next_checkpoint <- 1L
  emit <- function(horizon) {
    dest_id <- cell_id_from_lonlat(lonp, latp, grid_info)
    if (isTRUE(params$inactive_particles_are_lost)) {
      dest_id[!active] <- NA_integer_
    }

    result <- data.table(source_cell_id = source_ids,
               lon = lonp,
               lat = latp,
               active = active,
               sink_cell_id = dest_id)
    on_snapshot(horizon, result)
    invisible(NULL)
  }
  for (step in seq_len(n_steps)) {
    if (!any(active)) {
      # Later horizons have the same frozen state, including inactive/lost flags.
      for (jj in seq.int(next_checkpoint, length(horizons))) emit(horizons[jj])
      return(invisible(NULL))
    }
    elapsed_days <- (step - 1L) / steps_per_day
    base_day <- floor(elapsed_days)
    a <- elapsed_days - base_day
    idx0 <- start_global_index + base_day
    idx1 <- idx0 + 1L

    vel0 <- read_velocity_slice(idx0, scenario, slice_table, grid_info$grid_def, params, cache_env)
    vel1 <- read_velocity_slice(idx1, scenario, slice_table, grid_info$grid_def, params, cache_env)

    active_idx <- which(active)
    u0 <- interp_grid(vel0$u, lonp[active_idx], latp[active_idx],
                      grid_info$lon, grid_info$lat,
                      require_all_four = params$interp_requires_all_four_corners)
    v0 <- interp_grid(vel0$v, lonp[active_idx], latp[active_idx],
                      grid_info$lon, grid_info$lat,
                      require_all_four = params$interp_requires_all_four_corners)
    u1 <- interp_grid(vel1$u, lonp[active_idx], latp[active_idx],
                      grid_info$lon, grid_info$lat,
                      require_all_four = params$interp_requires_all_four_corners)
    v1 <- interp_grid(vel1$v, lonp[active_idx], latp[active_idx],
                      grid_info$lon, grid_info$lat,
                      require_all_four = params$interp_requires_all_four_corners)

    u <- (1 - a) * u0 + a * u1
    v <- (1 - a) * v0 + a * v1

    good <- is.finite(u) & is.finite(v) & is.finite(lonp[active_idx]) & is.finite(latp[active_idx])
    if (any(!good)) active[active_idx[!good]] <- FALSE
    if (any(good)) {

      ok_idx <- active_idx[good]
      lat_rad <- latp[ok_idx] * pi / 180
      cos_lat <- cos(lat_rad)
      cos_lat[abs(cos_lat) < 1e-12] <- NA_real_

      dlon <- (u[good] * params$dt_seconds / (params$R_earth_m * cos_lat)) * 180 / pi
      dlat <- (v[good] * params$dt_seconds / params$R_earth_m) * 180 / pi

      if (is.finite(scenario$Kh) && scenario$Kh > 0) {
        dx <- stats::rnorm(length(ok_idx), mean = 0, sd = sqrt(2 * scenario$Kh * params$dt_seconds))
        dy <- stats::rnorm(length(ok_idx), mean = 0, sd = sqrt(2 * scenario$Kh * params$dt_seconds))
        dlon <- dlon + (dx / (params$R_earth_m * cos_lat)) * 180 / pi
        dlat <- dlat + (dy / params$R_earth_m) * 180 / pi
      }

      finite_step <- is.finite(dlon) & is.finite(dlat)
      if (any(!finite_step)) active[ok_idx[!finite_step]] <- FALSE
      move_idx <- ok_idx[finite_step]
      lonp[move_idx] <- lonp[move_idx] + dlon[finite_step]
      latp[move_idx] <- latp[move_idx] + dlat[finite_step]
    }
    if (step == checkpoint_steps[next_checkpoint]) {
      emit(horizons[next_checkpoint])
      next_checkpoint <- next_checkpoint + 1L
    }
  }
  invisible(NULL)
}

advect_particles <- function(source_dt, start_global_index, horizon_days,
                             scenario, slice_table, grid_info, params, cache_env) {
  result <- NULL
  advect_particle_checkpoints(source_dt, start_global_index, horizon_days,
    scenario, slice_table, grid_info, params, cache_env,
    function(h, particles) result <<- particles)
  result
}

# ==============================================================================
# 6. SOURCES, SCENARIOS, SPARSE MATRICES, AND BUDGET TABLES
# ==============================================================================

apply_source_subset <- function(source_dt, params) {
  ss <- params$source_subset_lonlat
  if (all(vapply(ss, function(z) !is.finite(z), logical(1)))) return(source_dt)
  keep <- rep(TRUE, nrow(source_dt))
  if (is.finite(ss$lon_min)) keep <- keep & source_dt$lon >= ss$lon_min
  if (is.finite(ss$lon_max)) keep <- keep & source_dt$lon <= ss$lon_max
  if (is.finite(ss$lat_min)) keep <- keep & source_dt$lat >= ss$lat_min
  if (is.finite(ss$lat_max)) keep <- keep & source_dt$lat <= ss$lat_max
  source_dt[keep]
}

select_source_cells <- function(grid_info, params) {
  wet <- copy(grid_info$wet_grid)

  if (params$source_selection == "pipe") {
    if (!is.finite(params$pipe_lon) || !is.finite(params$pipe_lat)) {
      stop("pipe_lon and pipe_lat must be finite when source_selection = 'pipe'.",
           call. = FALSE)
    }
    lat0 <- params$pipe_lat * pi / 180
    approx_dx <- (wet$lon - params$pipe_lon) * cos(lat0)
    approx_dy <- wet$lat - params$pipe_lat
    nearest <- wet[which.min(approx_dx^2 + approx_dy^2)]
    nearest[, source_label := "pipe_nearest_wet_cell"]
    log_message("Pipe coordinate ", params$pipe_lon, ", ", params$pipe_lat,
                " mapped to wet cell_id ", nearest$cell_id,
                " at lon/lat ", round(nearest$lon, 6), ", ", round(nearest$lat, 6), ".")
    return(nearest)
  }

  if (params$source_selection == "cell_ids") {
    source_dt <- wet[cell_id %in% params$source_cell_ids]
    missing <- setdiff(params$source_cell_ids, source_dt$cell_id)
    if (length(missing) > 0L) {
      log_message("Ignoring source_cell_ids that are absent or dry: ",
                  paste(missing, collapse = ", "), ".", level = "WARN")
    }
    if (nrow(source_dt) == 0L) stop("No requested source_cell_ids are wet.", call. = FALSE)
    source_dt[, source_label := "requested_cell_id"]
    source_dt <- apply_source_subset(source_dt, params)
    if (nrow(source_dt) == 0L) stop("Source lon/lat subset removed all source cells.", call. = FALSE)
    return(source_dt)
  }

  if (params$source_selection == "all_wet_cells") {
    source_dt <- apply_source_subset(wet, params)
    if (nrow(source_dt) == 0L) stop("Source lon/lat subset removed all wet cells.", call. = FALSE)
    source_dt[, source_label := "all_wet_cells"]
    log_message("Using all selected wet cells as sources: ",
                format(nrow(source_dt), big.mark = ","),
                ". Particles per release day = ",
                format(nrow(source_dt) * params$N_particles_per_source, big.mark = ","), ".")
    return(source_dt)
  }

  stop("Unsupported source_selection: ", params$source_selection, call. = FALSE)
}

make_scenarios <- function(params) {
  rows <- list()
  for (mode in params$velocity_modes) {
    w_values <- if (mode == "mixed") params$mixed_w_values else if (mode == "surface") 1 else 0
    for (w in w_values) {
      for (Kh in params$Kh_values) {
        rows[[length(rows) + 1L]] <- data.table(mode = mode, w = w, Kh = Kh)
      }
    }
  }
  out <- rbindlist(rows)
  out[, scenario_id := sanitize_label(paste0(mode, "_w", sprintf("%.2f", w),
                                             "_Kh", sprintf("%.3g", Kh)))]
  out[]
}

choose_start_indices <- function(slice_table, params) {
  if (!is.null(params$start_global_indices) && length(params$start_global_indices) > 0L) {
    idx <- sort(unique(as.integer(params$start_global_indices)))
  } else {
    idx <- seq.int(1L, nrow(slice_table) - 1L, by = params$start_day_stride)
    if (!is.null(params$max_start_days)) {
      idx <- head(idx, as.integer(params$max_start_days))
    }
  }
  idx <- idx[idx >= 1L & idx < nrow(slice_table)]
  if (length(idx) == 0L) stop("No valid start days selected.", call. = FALSE)
  idx
}

build_transition <- function(particles, source_dt, grid_info, params) {
  N <- params$N_particles_per_source
  valid <- !is.na(particles$sink_cell_id)

  if (any(valid)) {
    edge_counts <- particles[valid, .(particle_count = .N),
                             by = .(source_cell_id, sink_cell_id)]
    edge_counts[, transfer_probability := particle_count / N]
  } else {
    edge_counts <- data.table(source_cell_id = integer(), sink_cell_id = integer(),
                              particle_count = integer(), transfer_probability = numeric())
  }

  P <- Matrix::sparseMatrix(i = edge_counts$source_cell_id,
                            j = edge_counts$sink_cell_id,
                            x = edge_counts$transfer_probability,
                            dims = c(grid_info$n_cells, grid_info$n_cells))

  list(edge_counts = edge_counts, P = P)
}

make_daily_budget <- function(edge_counts, source_dt, grid_info, params) {
  N <- params$N_particles_per_source
  n_sources <- nrow(source_dt)
  total_released <- N * n_sources

  source_stats <- source_dt[, .(cell_id, released_particles = N)]

  wet_arrived <- if (nrow(edge_counts) > 0L) {
    edge_counts[, .(wet_arrived_particles = sum(particle_count)), by = .(cell_id = source_cell_id)]
  } else {
    data.table(cell_id = integer(), wet_arrived_particles = integer())
  }

  retained <- if (nrow(edge_counts) > 0L) {
    edge_counts[source_cell_id == sink_cell_id,
                .(retained_particles = sum(particle_count)), by = .(cell_id = source_cell_id)]
  } else {
    data.table(cell_id = integer(), retained_particles = integer())
  }

  transported_out <- if (nrow(edge_counts) > 0L) {
    edge_counts[source_cell_id != sink_cell_id,
                .(transported_out_particles = sum(particle_count)), by = .(cell_id = source_cell_id)]
  } else {
    data.table(cell_id = integer(), transported_out_particles = integer())
  }

  incoming <- if (nrow(edge_counts) > 0L) {
    edge_counts[source_cell_id != sink_cell_id,
                .(incoming_added_particles = sum(particle_count),
                  incoming_added_source_equiv_fraction = sum(transfer_probability)),
                by = .(cell_id = sink_cell_id)]
  } else {
    data.table(cell_id = integer(), incoming_added_particles = integer(),
               incoming_added_source_equiv_fraction = numeric())
  }

  if (isTRUE(params$budget_include_all_wet_cells)) {
    cells <- grid_info$wet_grid$cell_id
  } else {
    cells <- sort(unique(c(source_stats$cell_id, incoming$cell_id)))
  }
  budget <- data.table(cell_id = cells)
  budget <- merge(budget, grid_info$grid[, .(cell_id, lon_index, lat_index, lon, lat, wet)],
                  by = "cell_id", all.x = TRUE)
  budget <- merge(budget, source_stats, by = "cell_id", all.x = TRUE)
  budget <- merge(budget, wet_arrived, by = "cell_id", all.x = TRUE)
  budget <- merge(budget, retained, by = "cell_id", all.x = TRUE)
  budget <- merge(budget, transported_out, by = "cell_id", all.x = TRUE)
  budget <- merge(budget, incoming, by = "cell_id", all.x = TRUE)

  numeric_zero_cols <- c("released_particles", "wet_arrived_particles",
                         "retained_particles", "transported_out_particles",
                         "incoming_added_particles", "incoming_added_source_equiv_fraction")
  for (cc in numeric_zero_cols) {
    if (!cc %in% names(budget)) budget[, (cc) := 0]
    set(budget, which(is.na(budget[[cc]])), cc, 0)
  }

  budget[, external_loss_particles := pmax(0, released_particles - wet_arrived_particles)]
  budget[, water_removed_particles := pmax(0, released_particles - retained_particles)]

  budget[, retained_fraction_of_released := fifelse(released_particles > 0,
                                                    retained_particles / released_particles,
                                                    NA_real_)]
  budget[, transported_out_fraction_of_released := fifelse(released_particles > 0,
                                                           transported_out_particles / released_particles,
                                                           NA_real_)]
  budget[, external_loss_fraction_of_released := fifelse(released_particles > 0,
                                                         external_loss_particles / released_particles,
                                                         NA_real_)]
  budget[, water_removed_fraction_of_released := fifelse(released_particles > 0,
                                                         water_removed_particles / released_particles,
                                                         NA_real_)]
  # total_released is a scalar, while incoming_added_particles is one value per
  # budget row. data.table::fifelse() requires the test to be the same length as
  # the yes/no vectors, so use a normal scalar if/else here.
  if (total_released > 0) {
    budget[, incoming_added_fraction_of_total_released := incoming_added_particles / total_released]
  } else {
    budget[, incoming_added_fraction_of_total_released := NA_real_]
  }
  budget[]
}

annotate_edges <- function(edge_counts, source_dt, grid_info, scenario, horizon_days,
                           start_row, params) {
  if (nrow(edge_counts) == 0L) {
    return(data.table())
  }
  out <- copy(edge_counts)
  if (params$min_edge_probability_to_write > 0) {
    out <- out[transfer_probability >= params$min_edge_probability_to_write]
  }
  if (nrow(out) == 0L) return(out)

  src_xy <- grid_info$grid[, .(source_cell_id = cell_id,
                               source_lon = lon, source_lat = lat,
                               source_lon_index = lon_index, source_lat_index = lat_index)]
  sink_xy <- grid_info$grid[, .(sink_cell_id = cell_id,
                                sink_lon = lon, sink_lat = lat,
                                sink_lon_index = lon_index, sink_lat_index = lat_index)]
  out <- merge(out, src_xy, by = "source_cell_id", all.x = TRUE)
  out <- merge(out, sink_xy, by = "sink_cell_id", all.x = TRUE)
  out[, `:=`(scenario_id = scenario$scenario_id,
             mode = scenario$mode,
             w = scenario$w,
             Kh = scenario$Kh,
             horizon_days = as.integer(horizon_days),
             start_global_index = start_row$global_index,
             start_file = start_row$file_name,
             start_time_index = start_row$time_index,
             start_velocity_time_utc = start_row$velocity_time_utc,
             N_particles_per_source = params$N_particles_per_source)]
  setcolorder(out, c("scenario_id", "mode", "w", "Kh", "horizon_days",
                    "start_global_index", "start_file", "start_time_index",
                    "start_velocity_time_utc", "source_cell_id", "sink_cell_id",
                    "source_lon", "source_lat", "sink_lon", "sink_lat",
                    "source_lon_index", "source_lat_index", "sink_lon_index", "sink_lat_index",
                    "particle_count", "N_particles_per_source", "transfer_probability"))
  out[]
}

annotate_budget <- function(budget, scenario, horizon_days, start_row, params) {
  out <- copy(budget)
  out[, `:=`(scenario_id = scenario$scenario_id,
             mode = scenario$mode,
             w = scenario$w,
             Kh = scenario$Kh,
             horizon_days = as.integer(horizon_days),
             start_global_index = start_row$global_index,
             start_file = start_row$file_name,
             start_time_index = start_row$time_index,
             start_velocity_time_utc = start_row$velocity_time_utc,
             N_particles_per_source = params$N_particles_per_source)]
  setcolorder(out, c("scenario_id", "mode", "w", "Kh", "horizon_days",
                    "start_global_index", "start_file", "start_time_index", "start_velocity_time_utc",
                    "cell_id", "lon_index", "lat_index", "lon", "lat", "wet",
                    "released_particles", "wet_arrived_particles", "retained_particles",
                    "transported_out_particles", "external_loss_particles", "water_removed_particles",
                    "incoming_added_particles", "N_particles_per_source",
                    "retained_fraction_of_released", "transported_out_fraction_of_released",
                    "external_loss_fraction_of_released", "water_removed_fraction_of_released",
                    "incoming_added_fraction_of_total_released",
                    "incoming_added_source_equiv_fraction"))
  out[]
}

make_qc_row <- function(particles, edge_counts, source_dt, scenario, horizon_days, start_row, params) {
  N <- params$N_particles_per_source
  n_released <- nrow(source_dt) * N
  n_wet_end <- sum(!is.na(particles$sink_cell_id))
  n_lost <- n_released - n_wet_end
  row_sums <- if (nrow(edge_counts) > 0L) {
    edge_counts[, .(row_sum = sum(transfer_probability)), by = source_cell_id]
  } else {
    data.table(source_cell_id = integer(), row_sum = numeric())
  }
  data.table(scenario_id = scenario$scenario_id,
             mode = scenario$mode,
             w = scenario$w,
             Kh = scenario$Kh,
             horizon_days = as.integer(horizon_days),
             start_global_index = start_row$global_index,
             start_file = start_row$file_name,
             start_time_index = start_row$time_index,
             start_velocity_time_utc = start_row$velocity_time_utc,
             n_source_cells = nrow(source_dt),
             N_particles_per_source = N,
             n_particles_released = n_released,
             n_particles_wet_endpoint = n_wet_end,
             n_particles_lost = n_lost,
             lost_fraction = ifelse(n_released > 0, n_lost / n_released, NA_real_),
             n_nonzero_edges = nrow(edge_counts),
             min_sparse_row_sum = ifelse(nrow(row_sums) > 0L, min(row_sums$row_sum), NA_real_),
             mean_sparse_row_sum = ifelse(nrow(row_sums) > 0L, mean(row_sums$row_sum), NA_real_),
             max_sparse_row_sum = ifelse(nrow(row_sums) > 0L, max(row_sums$row_sum), NA_real_))
}

# ==============================================================================
# 7. ACCUMULATORS FOR MEAN MATRICES AND MAPS
# ==============================================================================

compact_edge_prob_sum <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0L) {
    return(data.table(source_cell_id = integer(), sink_cell_id = integer(), probability_sum = numeric()))
  }
  dt[, .(probability_sum = sum(probability_sum)), by = .(source_cell_id, sink_cell_id)]
}

update_edge_accumulator <- function(edge_env, key, edge_counts) {
  if (!exists(key, envir = edge_env, inherits = FALSE)) {
    assign(key, list(edge_sum = data.table(source_cell_id = integer(), sink_cell_id = integer(),
                                           probability_sum = numeric()),
                     n_days = 0L, n_updates_since_compact = 0L),
           envir = edge_env)
  }
  acc <- get(key, envir = edge_env, inherits = FALSE)
  if (nrow(edge_counts) > 0L) {
    tmp <- edge_counts[, .(source_cell_id, sink_cell_id,
                           probability_sum = transfer_probability)]
    acc$edge_sum <- rbindlist(list(acc$edge_sum, tmp), use.names = TRUE)
  }
  acc$n_days <- acc$n_days + 1L
  acc$n_updates_since_compact <- acc$n_updates_since_compact + 1L
  assign(key, acc, envir = edge_env)
}

maybe_compact_edge_accumulator <- function(edge_env, key, params) {
  acc <- get(key, envir = edge_env, inherits = FALSE)
  if (acc$n_updates_since_compact >= params$compact_edge_accumulator_every) {
    acc$edge_sum <- compact_edge_prob_sum(acc$edge_sum)
    acc$n_updates_since_compact <- 0L
    assign(key, acc, envir = edge_env)
    gc(verbose = FALSE)
  }
  invisible(NULL)
}

update_map_accumulator <- function(map_env, key, budget, grid_info) {
  if (!exists(key, envir = map_env, inherits = FALSE)) {
    assign(key, list(removed_sum = numeric(grid_info$n_cells),
                     added_sum = numeric(grid_info$n_cells),
                     n_days = 0L),
           envir = map_env)
  }
  acc <- get(key, envir = map_env, inherits = FALSE)

  rem <- budget[released_particles > 0 & is.finite(water_removed_fraction_of_released),
                .(cell_id, value = water_removed_fraction_of_released)]
  if (nrow(rem) > 0L) {
    acc$removed_sum[rem$cell_id] <- acc$removed_sum[rem$cell_id] + rem$value
  }
  add <- budget[incoming_added_source_equiv_fraction > 0,
                .(cell_id, value = incoming_added_source_equiv_fraction)]
  if (nrow(add) > 0L) {
    acc$added_sum[add$cell_id] <- acc$added_sum[add$cell_id] + add$value
  }
  acc$n_days <- acc$n_days + 1L
  assign(key, acc, envir = map_env)
  invisible(NULL)
}

# ==============================================================================
# 8. MAPS AND FINAL OUTPUT WRITERS
# ==============================================================================

plot_point_map <- function(dt, value_col, title, output_file, params) {
  assert_can_write(output_file, params$overwrite_existing)
  values <- dt[[value_col]]
  cols <- make_colour_vector(values, params)
  finite <- is.finite(values)
  rng <- if (any(finite)) range(values[finite], finite = TRUE) else c(NA_real_, NA_real_)

  grDevices::png(output_file, width = params$map_width_px,
                 height = params$map_height_px, res = params$map_res_dpi)
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(4.2, 4.5, 4.2, 1.0))
  graphics::plot(dt$lon, dt$lat, pch = 15, cex = params$map_point_cex,
                 col = cols, asp = 1,
                 xlab = "Longitude (degrees east)",
                 ylab = "Latitude (degrees north)",
                 main = title)
  graphics::legend("topright",
                   legend = c(paste0("low: ", signif(rng[1L], 4)),
                              paste0("high: ", signif(rng[2L], 4))),
                   pch = 15,
                   col = c(params$map_low_high_colours[1L], params$map_low_high_colours[2L]),
                   bty = "n")
  invisible(output_file)
}

map_horizon_selected <- function(horizon_days, params) {
  mh <- params$map_horizons_days
  if (identical(mh, "all")) return(TRUE)
  horizon_days %in% as.integer(mh)
}

write_maps <- function(map_env, grid_info, params) {
  if (!isTRUE(params$save_maps)) return(character())
  keys <- ls(envir = map_env)
  if (length(keys) == 0L) return(character())

  written <- character()
  map_grid <- grid_info$wet_grid[, .(cell_id, lon, lat)]
  for (key in keys) {
    acc <- get(key, envir = map_env, inherits = FALSE)
    if (acc$n_days <= 0L) next

    dt <- copy(map_grid)
    dt[, mean_water_removed_fraction := acc$removed_sum[cell_id] / acc$n_days]
    dt[, mean_incoming_added_source_equiv_fraction := acc$added_sum[cell_id] / acc$n_days]

    removal_file <- file.path(params$output_dir,
                              paste0(params$output_prefix, "_", key, "_source_removed.png"))
    added_file <- file.path(params$output_dir,
                            paste0(params$output_prefix, "_", key, "_transported_to.png"))

    plot_point_map(dt, "mean_water_removed_fraction",
                   paste0("Mean water removed from source cells: ", key),
                   removal_file, params)
    plot_point_map(dt, "mean_incoming_added_source_equiv_fraction",
                   paste0("Mean transported water arriving in cells: ", key),
                   added_file, params)
    written <- c(written, removal_file, added_file)
    log_message("Wrote maps for ", key, ".")
  }
  written
}

write_mean_outputs <- function(edge_env, grid_info, params) {
  keys <- ls(envir = edge_env)
  if (length(keys) == 0L) {
    return(list(mean_matrix_files = character(), mean_edges_file = NA_character_))
  }

  mean_matrix_dir <- make_mean_matrix_dir(params)
  if (isTRUE(params$save_mean_transition_matrices)) {
    ensure_dir(mean_matrix_dir)
  }

  mean_edges_path <- file.path(params$output_dir, paste0(params$output_prefix, "_mean_sparse_edges.csv"))
  if (isTRUE(params$save_mean_sparse_edges_csv)) prepare_output_file(mean_edges_path, params$overwrite_existing)

  mean_matrix_files <- character()
  for (key in keys) {
    acc <- get(key, envir = edge_env, inherits = FALSE)
    edge_sum <- compact_edge_prob_sum(acc$edge_sum)
    if (nrow(edge_sum) > 0L) {
      edge_sum[, mean_transfer_probability := probability_sum / acc$n_days]
      edge_sum[, probability_sum := NULL]
    } else {
      edge_sum[, mean_transfer_probability := numeric()]
    }

    if (isTRUE(params$save_mean_sparse_edges_csv) && nrow(edge_sum) > 0L) {
      tmp <- copy(edge_sum)
      tmp[, `:=`(matrix_key = key, n_release_days = acc$n_days)]
      setcolorder(tmp, c("matrix_key", "n_release_days", "source_cell_id", "sink_cell_id",
                         "mean_transfer_probability"))
      append_fwrite(tmp, mean_edges_path)
    }

    if (isTRUE(params$save_mean_transition_matrices)) {
      P_mean <- Matrix::sparseMatrix(i = edge_sum$source_cell_id,
                                     j = edge_sum$sink_cell_id,
                                     x = edge_sum$mean_transfer_probability,
                                     dims = c(grid_info$n_cells, grid_info$n_cells))
      matrix_path <- make_mean_matrix_path(mean_matrix_dir, key, params)
      safe_save_rds(list(P_mean = P_mean,
                         matrix_key = key,
                         n_release_days = acc$n_days,
                         note = "Rows are source cells, columns are sink cells. Row sums can be < 1 when particles leave wet grid/domain."),
                    matrix_path,
                    params$overwrite_existing,
                    compress = FALSE)
      mean_matrix_files <- c(mean_matrix_files, matrix_path)
    }
  }

  list(mean_matrix_files = mean_matrix_files,
       mean_edges_file = if (isTRUE(params$save_mean_sparse_edges_csv)) mean_edges_path else NA_character_)
}

# ==============================================================================
# 9. LOCAL INPUT STORE AND PARALLEL EXECUTION
# ==============================================================================

stage_input_files <- function(files, work_dir) {
  out <- copy(files)
  for (ii in seq_len(nrow(files))) {
    destination <- file.path(work_dir, paste0("input_", ii, ".nc4"))
    log_message("Staging input once: ", files$file_name[ii], ".")
    ok <- file.copy(files$full_path[ii], destination, overwrite = FALSE,
                    copy.mode = FALSE, copy.date = FALSE)
    if (!isTRUE(ok) || file.info(destination)$size != file.info(files$full_path[ii])$size) {
      stop("Could not stage complete input on local disk: ", files$full_path[ii],
           ". Check available space and staging_dir.", call. = FALSE)
    }
    out$full_path[ii] <- destination
  }
  out
}

make_run_plan <- function(scenarios, start_indices, slice_table, params) {
  rows <- list()
  for (sc_i in seq_len(nrow(scenarios))) {
    for (h_i in seq_along(params$tracking_horizons_days)) {
      h <- params$tracking_horizons_days[h_i]
      starts <- start_indices[start_indices + h <= nrow(slice_table)]
      if (!length(starts)) {
        log_message("No valid start days for horizon ", h,
                    " because selected files do not contain enough following velocity slices.",
                    level = "WARN")
        next
      }
      rows[[length(rows) + 1L]] <- data.table(sc_i = sc_i, h_i = h_i,
        horizon = h, start = starts)
    }
  }
  if (!length(rows)) return(data.table())
  out <- rbindlist(rows)
  out[, run_id := .I]
  # Month labels are only scheduling labels. Never reset global indices here.
  month_key <- substr(slice_table$velocity_time_utc, 1L, 7L)
  missing <- is.na(month_key) | !nzchar(month_key)
  month_key[missing] <- paste0("file_", slice_table$file_order[missing])
  # Contiguous groups preserve input order even if a month label recurs.
  out[, part := if (params$parallel_split == "month")
    data.table::rleid(month_key)[start] else 1L]
  out[]
}

stage_velocity_fields <- function(staged_files, slice_table, run_plan, grid_info, params) {
  required <- rep(FALSE, nrow(slice_table))
  if (nrow(run_plan)) {
    # Interpolation at the last substep needs the endpoint day as well.
    windows <- run_plan[, .(last = max(start + horizon)), by = start]
    for (ii in seq_len(nrow(windows))) required[seq.int(windows$start[ii], windows$last[ii])] <- TRUE
  }
  ensure_dir(params$velocity_store)
  bytes <- sum(required) * grid_info$n_cells * 4 * 8
  log_message("Preparing ", sum(required), " required daily bundles (~",
              round(bytes / 1024^3, 2), " GiB of raw fields on local disk).")
  for (ii in seq_len(nrow(staged_files))) {
    rows <- slice_table[required & file_order == staged_files$file_order[ii]]
    if (nrow(rows)) {
      local({
        nc <- ncdf4::nc_open(staged_files$full_path[ii])
        on.exit(ncdf4::nc_close(nc))
        for (jj in seq_len(nrow(rows))) {
          tt <- rows$time_index[jj]
          raw <- list(
            su = read_daily_var(nc, params$surface_u_var, tt, grid_info$grid_def, params),
            sv = read_daily_var(nc, params$surface_v_var, tt, grid_info$grid_def, params),
            bu = read_daily_var(nc, params$bottom_u_var, tt, grid_info$grid_def, params),
            bv = read_daily_var(nc, params$bottom_v_var, tt, grid_info$grid_def, params))
          saveRDS(raw, file.path(params$velocity_store, paste0("g", rows$global_index[jj], ".rds")),
                  compress = FALSE)
        }
      })
    }
    # The shared daily store is now sufficient; do not retain duplicate NetCDFs.
    unlink(staged_files$full_path[ii])
  }
  invisible(sum(required))
}

# Advancing uniform draws does not need the arrays of release coordinates.
# Chunking bounds memory and consumes the same stream as the original two runif
# calls. This is valid only for Kh <= 0 (or non-finite Kh, as in the old tracker).
advance_release_rng <- function(n_particles, jitter) {
  if (!isTRUE(jitter)) return(invisible(NULL))
  remaining <- 2 * as.double(n_particles)
  while (remaining > 0) {
    n <- min(remaining, 1000000)
    invisible(stats::runif(n))
    remaining <- remaining - n
  }
  invisible(NULL)
}

result_path <- function(run_id, params) file.path(params$result_store, paste0("r", run_id, ".rds"))

process_tracking_job <- function(job, context) {
  params <- context$params
  grid_info <- context$grid_info
  source_dt <- context$source_dt
  slice_table <- context$slice_table
  scenario <- as.list(context$scenarios[job$runs$sc_i[1L]])
  runs <- job$runs
  cache <- new_velocity_cache(params)
  finished <- 0L
  save_result <- function(row, particles) {
    start_row <- slice_table[global_index == row$start][1L]
    transition <- build_transition(particles, source_dt, grid_info, params)
    budget <- make_daily_budget(transition$edge_counts, source_dt, grid_info, params)
    result <- list(run_id = row$run_id, edge_counts = transition$edge_counts,
      P = if (isTRUE(params$save_daily_transition_matrices)) transition$P else NULL,
      budget = budget,
      qc = make_qc_row(particles, transition$edge_counts, source_dt,
                      scenario, row$horizon, start_row, params),
      endpoints = if (isTRUE(params$save_particle_endpoints)) particles else NULL)
    # Unique private paths: no two workers append to the same output file.
    saveRDS(result, result_path(row$run_id, params), compress = FALSE)
    finished <<- finished + 1L
    invisible(NULL)
  }
  if (isTRUE(job$reuse)) {
    for (g in unique(runs$start)) {
      selected <- runs[start == g]
      if (!is.null(job$seeds)) {
        assign(".Random.seed", job$seeds[[as.character(selected$run_id[1L])]], envir = .GlobalEnv)
      }
      advect_particle_checkpoints(source_dt, g, selected$horizon, scenario,
        slice_table, grid_info, params, cache,
        function(h, particles) save_result(selected[horizon == h][1L], particles))
    }
  } else {
    for (ii in seq_len(nrow(runs))) {
      row <- runs[ii]
      if (!is.null(job$seeds)) {
        assign(".Random.seed", job$seeds[[as.character(row$run_id)]], envir = .GlobalEnv)
      }
      particles <- advect_particles(source_dt, row$start, row$horizon, scenario,
                                    slice_table, grid_info, params, cache)
      save_result(row, particles)
    }
  }
  list(n_results = finished, peak_velocity_days = cache$peak_days,
       local_velocity_reads = cache$reads)
}

make_tracking_jobs <- function(runs, seeds, context, reuse) {
  groups <- unique(runs[, .(sc_i, part)])
  lapply(seq_len(nrow(groups)), function(ii) {
    selected <- runs[sc_i == groups$sc_i[ii] & part == groups$part[ii]]
    list(runs = selected,
         seeds = if (is.null(seeds)) NULL else seeds[as.character(selected$run_id)],
         reuse = reuse)
  })
}

tracking_worker_entry <- function(job) process_tracking_job(job, tracking_context)

execute_tracking_jobs <- function(run_plan, context) {
  params <- context$params
  if (!nrow(run_plan)) return(list(n_results = 0L, peak_velocity_days = 0L, jobs = 0L))
  needed_workers <- min(params$parallel_workers, nrow(unique(run_plan[, .(sc_i, part)])))
  cluster <- NULL
  # Preserve the caller's RNG through worker creation. Workers never run main().
  if (needed_workers > 1L) {
    seed_before <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      get(".Random.seed", envir = .GlobalEnv) else NULL
    cluster <- tryCatch(parallel::makePSOCKcluster(needed_workers), error = function(e) {
      log_message("Could not start PSOCK workers; using serial execution: ", conditionMessage(e),
                  level = "WARN")
      NULL
    })
    if (!is.null(seed_before)) assign(".Random.seed", seed_before, envir = .GlobalEnv)
    if (!is.null(cluster)) {
      on.exit(parallel::stopCluster(cluster), add = TRUE)
      exports <- c("new_velocity_cache", "read_velocity_slice", "interp_grid",
        "cell_id_from_lonlat", "release_particles", "advect_particle_checkpoints",
        "advect_particles", "build_transition", "make_daily_budget", "make_qc_row",
        "result_path", "process_tracking_job", "tracking_worker_entry")
      parallel::clusterExport(cluster, exports, envir = environment(execute_tracking_jobs))
      parallel::clusterCall(cluster, function(worker_threads) {
        suppressPackageStartupMessages({ library(data.table); library(Matrix) })
        data.table::setDTthreads(worker_threads)
        NULL
      }, params$data_table_threads %||% 1L)
      tracking_context <- context
      parallel::clusterExport(cluster, "tracking_context", envir = environment())
    }
  }
  log_message("Tracking workers: ", if (is.null(cluster)) 1L else length(cluster),
              "; trajectory_reuse = ", params$trajectory_reuse, ".")
  stats <- list()
  dispatch <- function(jobs) {
    if (!length(jobs)) return(invisible(NULL))
    # A bounded queue: one wave per worker, with progress between waves.
    width <- if (is.null(cluster)) 1L else length(cluster)
    for (first in seq.int(1L, length(jobs), by = width)) {
      wave <- jobs[seq.int(first, min(length(jobs), first + width - 1L))]
      if (is.null(cluster)) {
        # Predetermined seeds must not alter the coordinator's future stream.
        seed <- get(".Random.seed", envir = .GlobalEnv)
        ans <- lapply(wave, process_tracking_job, context = context)
        assign(".Random.seed", seed, envir = .GlobalEnv)
      } else {
        ans <- parallel::parLapplyLB(cluster, wave, tracking_worker_entry)
      }
      stats <<- c(stats, ans)
      log_message("Finished ", sum(vapply(ans, `[[`, integer(1), "n_results")),
                  " scenario/horizon/release results in ", length(wave), " job(s).")
    }
    invisible(NULL)
  }
  # Ensure an RNG state also exists for random_seed = NULL, like the first legacy
  # release would do. An unseeded run intentionally has no legacy fixed stream.
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) set.seed(NULL)
  if (params$trajectory_reuse == "shared") {
    log_message("Shared trajectories selected: stochastic samples differ from the original script.", level = "WARN")
    # Stable stream allocation before scheduling; workers/chunking do not affect it.
    RNGkind("L'Ecuyer-CMRG")
    if (!is.null(params$random_seed) && is.finite(params$random_seed)) set.seed(as.integer(params$random_seed))
    stream <- get(".Random.seed", envir = .GlobalEnv)
    seeds <- vector("list", nrow(run_plan)); names(seeds) <- as.character(run_plan$run_id)
    releases <- unique(run_plan[, .(sc_i, start)])
    setorder(releases, sc_i, start)
    for (ii in seq_len(nrow(releases))) {
      ids <- run_plan[sc_i == releases$sc_i[ii] & start == releases$start[ii], run_id]
      for (id in ids) seeds[[as.character(id)]] <- stream
      stream <- parallel::nextRNGStream(stream)
    }
    assign(".Random.seed", stream, envir = .GlobalEnv)
    dispatch(make_tracking_jobs(run_plan, seeds, context, reuse = TRUE))
  } else {
    # Positive-Kh runs are barriers in the original global RNG sequence. Between
    # barriers, preallocate exact release seeds cheaply and dispatch in parallel.
    positive <- with(context$scenarios, is.finite(Kh) & Kh > 0)
    pending <- list()
    for (sc in unique(run_plan$sc_i)) {
      runs <- run_plan[sc_i == sc]
      if (positive[sc]) {
        dispatch(pending); pending <- list()
        log_message("Exact RNG: running diffusive scenario ", context$scenarios$scenario_id[sc],
                    " serially in original horizon/start order.")
        # No seed reset: the actual, survival-dependent rnorm calls advance the
        # coordinator stream for subsequent scenarios exactly as in the original.
        stats[[length(stats) + 1L]] <- process_tracking_job(
          list(runs = runs, seeds = NULL, reuse = FALSE), context)
      } else {
        seeds <- vector("list", nrow(runs)); names(seeds) <- as.character(runs$run_id)
        for (ii in seq_len(nrow(runs))) {
          seeds[[ii]] <- get(".Random.seed", envir = .GlobalEnv)
          advance_release_rng(nrow(context$source_dt) * as.double(params$N_particles_per_source),
                              params$release_jitter_within_cell)
        }
        pending <- c(pending, make_tracking_jobs(runs, seeds, context,
                         reuse = !isTRUE(params$release_jitter_within_cell)))
      }
    }
    dispatch(pending)
  }
  list(n_results = sum(vapply(stats, `[[`, integer(1), "n_results")),
       peak_velocity_days = max(vapply(stats, `[[`, integer(1), "peak_velocity_days")),
       local_velocity_reads = sum(vapply(stats, `[[`, integer(1), "local_velocity_reads")),
       jobs = length(stats), workers = if (is.null(cluster)) 1L else length(cluster))
}

# ==============================================================================
# 10. MAIN WORKFLOW
# ==============================================================================

main <- function() {
  start_time <- Sys.time()
  params <- collect_params()

  if (!is.null(params$random_seed) && is.finite(params$random_seed)) {
    set.seed(as.integer(params$random_seed))
  }

  ensure_dir(params$output_dir)
  ensure_dir(params$staging_dir)
  work_dir <- tempfile("lagrangian_", tmpdir = params$staging_dir)
  ensure_dir(work_dir)
  on.exit({
    if (!isTRUE(params$keep_temporary_files)) unlink(work_dir, recursive = TRUE)
  }, add = TRUE)
  params$velocity_store <- file.path(work_dir, "velocity")
  params$result_store <- file.path(work_dir, "results")
  ensure_dir(params$result_store)
  log_message("Temporary working directory: ", work_dir, ".")

  edge_csv_path <- file.path(params$output_dir, paste0(params$output_prefix, "_daily_sparse_edges.csv"))
  budget_csv_path <- file.path(params$output_dir, paste0(params$output_prefix, "_daily_cell_water_budget.csv"))
  qc_csv_path <- file.path(params$output_dir, paste0(params$output_prefix, "_day_qc.csv"))
  endpoints_csv_path <- file.path(params$output_dir, paste0(params$output_prefix, "_particle_endpoints.csv"))
  lookup_path <- file.path(params$output_dir, paste0(params$output_prefix, "_cell_lookup.csv"))
  metadata_path <- file.path(params$output_dir, paste0(params$output_prefix, "_metadata.rds"))
  log_path <- file.path(params$output_dir, paste0(params$output_prefix, "_run_log.txt"))

  for (path in c(if (isTRUE(params$save_daily_sparse_edges_csv)) edge_csv_path,
                 if (isTRUE(params$save_daily_cell_budget_csv)) budget_csv_path,
                 if (isTRUE(params$save_day_qc_csv)) qc_csv_path,
                 if (isTRUE(params$save_particle_endpoints)) endpoints_csv_path,
                 if (isTRUE(params$save_cell_lookup)) lookup_path,
                 if (isTRUE(params$save_metadata)) metadata_path,
                 log_path)) {
    if (!is.null(path) && !is.na(path)) prepare_output_file(path, params$overwrite_existing)
  }

  daily_matrix_dir <- make_daily_matrix_dir(params)
  matrix_manifest_path <- file.path(params$output_dir, paste0(params$output_prefix, "_daily_matrix_manifest.csv"))
  if (isTRUE(params$save_daily_transition_matrices)) {
    ensure_dir(daily_matrix_dir)
    prepare_output_file(matrix_manifest_path, params$overwrite_existing)
    log_message("Daily transition matrix RDS directory: ", daily_matrix_dir, ".")
  }

  log_message("Starting Lagrangian sparse connectivity workflow.")
  log_message("Modes: ", paste(params$velocity_modes, collapse = ", "),
              "; Kh: ", paste(params$Kh_values, collapse = ", "),
              "; horizons: ", paste(params$tracking_horizons_days, collapse = ", "), " days.")
  log_message("Particles per source cell: ", params$N_particles_per_source,
              "; dt_seconds: ", params$dt_seconds, ".")

  files <- discover_files(params)
  log_message("Found ", nrow(files), " NetCDF file(s).")

  staged_files <- stage_input_files(files, work_dir)
  grid_info <- prepare_grid_info(staged_files$full_path[1L], params)
  slice_table <- build_slice_table(staged_files, grid_info$grid_def, params)
  # Restore original provenance in every public output and metadata table.
  file_match <- match(slice_table$file_order, files$file_order)
  slice_table[, full_path := files$full_path[file_match]]
  slice_table[, file_name := files$file_name[file_match]]
  source_dt <- select_source_cells(grid_info, params)

  if (isTRUE(params$save_cell_lookup)) {
    lookup <- grid_info$grid[, .(cell_id, lon_index, lat_index, lon, lat,
                                 lon_min, lon_max, lat_min, lat_max, wet)]
    data.table::fwrite(lookup, lookup_path)
    log_message("Wrote cell lookup: ", lookup_path, ".")
  }

  scenarios <- make_scenarios(params)
  start_indices <- choose_start_indices(slice_table, params)
  if (anyDuplicated(scenarios$scenario_id)) {
    stop("Scenario labels collide. Use distinct mode/w/Kh settings with distinct labels.", call. = FALSE)
  }
  run_plan <- make_run_plan(scenarios, start_indices, slice_table, params)
  stage_velocity_fields(staged_files, slice_table, run_plan, grid_info, params)
  context <- list(params = params, grid_info = grid_info, source_dt = source_dt,
                  slice_table = slice_table, scenarios = scenarios)
  execution <- execute_tracking_jobs(run_plan, context)
  if (execution$n_results != nrow(run_plan)) stop("Incomplete tracking results.", call. = FALSE)
  log_message("Peak cached velocity days per worker: ", execution$peak_velocity_days,
              " (limit ", params$velocity_cache_days, ").")
  edge_accumulators <- new.env(parent = emptyenv())
  map_accumulators <- new.env(parent = emptyenv())

  day_qc_rows <- list()
  qc_counter <- 0L
  matrix_counter <- 0L

  # Reduce in the EXACT old scenario -> horizon -> release order. This keeps
  # row order, release-day weights and floating-point summation order unchanged.
  log_message("Merging daily results in original order.")
  for (rr in seq_len(nrow(run_plan))) {
    row <- run_plan[rr]
    scenario <- as.list(scenarios[row$sc_i])
    horizon_days <- params$tracking_horizons_days[row$h_i]
    start_global_index <- row$start
    start_row <- slice_table[global_index == start_global_index][1L]
    matrix_key <- sanitize_label(paste0(scenario$scenario_id, "_h", horizon_days, "d"))
    do_maps <- isTRUE(params$save_maps) && map_horizon_selected(horizon_days, params)
    private_path <- result_path(row$run_id, params)
    if (!file.exists(private_path)) stop("Missing result for run ", row$run_id, call. = FALSE)
    result <- readRDS(private_path)
    if (!identical(result$run_id, row$run_id)) stop("Result identity mismatch.", call. = FALSE)
    edge_counts <- result$edge_counts
    P <- result$P
    edge_dt <- annotate_edges(edge_counts, source_dt, grid_info, scenario,
                              horizon_days, start_row, params)
    if (isTRUE(params$save_daily_sparse_edges_csv) && nrow(edge_dt) > 0L) {
      append_fwrite(edge_dt, edge_csv_path)
    }

    budget <- result$budget
    budget_dt <- annotate_budget(budget, scenario, horizon_days, start_row, params)
    if (isTRUE(params$save_daily_cell_budget_csv) && nrow(budget_dt) > 0L) {
      append_fwrite(budget_dt, budget_csv_path)
    }

    if (isTRUE(params$save_particle_endpoints)) {
      endpoints <- result$endpoints[, .(scenario_id = scenario$scenario_id,
                                 mode = scenario$mode, w = scenario$w, Kh = scenario$Kh,
                                 horizon_days = as.integer(horizon_days),
                                 start_global_index = start_global_index,
                                 source_cell_id, sink_cell_id, lon, lat, active)]
      append_fwrite(endpoints, endpoints_csv_path)
    }

    qc_counter <- qc_counter + 1L
    day_qc_rows[[qc_counter]] <- result$qc

    if (isTRUE(params$save_daily_transition_matrices)) {
      matrix_counter <- matrix_counter + 1L
      matrix_path <- make_daily_matrix_path(daily_matrix_dir, matrix_key,
                                            start_global_index, params)
      safe_save_rds(list(P = P,
                         scenario = scenario,
                         horizon_days = horizon_days,
                         start_row = start_row,
                         source_cell_ids = source_dt$cell_id,
                         note = "Rows are source cells, columns are sink cells. Row sums can be < 1 when particles leave wet grid/domain."),
                    matrix_path,
                    params$overwrite_existing,
                    compress = FALSE)
      manifest <- data.table(matrix_key = matrix_key,
                             scenario_id = scenario$scenario_id,
                             mode = scenario$mode,
                             w = scenario$w,
                             Kh = scenario$Kh,
                             horizon_days = as.integer(horizon_days),
                             start_global_index = start_global_index,
                             start_file = start_row$file_name,
                             start_time_index = start_row$time_index,
                             start_velocity_time_utc = start_row$velocity_time_utc,
                             matrix_rds = matrix_path)
      append_fwrite(manifest, matrix_manifest_path)
    }

    update_edge_accumulator(edge_accumulators, matrix_key, edge_counts)
    maybe_compact_edge_accumulator(edge_accumulators, matrix_key, params)
    if (do_maps) update_map_accumulator(map_accumulators, matrix_key, budget, grid_info)
    if (!isTRUE(params$keep_temporary_files)) unlink(private_path)
  }

  day_qc <- if (length(day_qc_rows) > 0L) {
    rbindlist(day_qc_rows, use.names = TRUE, fill = TRUE)
  } else {
    data.table()
  }
  if (isTRUE(params$save_day_qc_csv)) {
    data.table::fwrite(day_qc, qc_csv_path)
    log_message("Wrote day QC: ", qc_csv_path, ".")
  }

  mean_outputs <- write_mean_outputs(edge_accumulators, grid_info, params)
  map_paths <- write_maps(map_accumulators, grid_info, params)

  runtime <- difftime(Sys.time(), start_time, units = "secs")
  metadata <- list(
    product_id = params$product_id,
    created_at = Sys.time(),
    runtime_seconds = as.numeric(runtime),
    parameters = params,
    file_list = files,
    grid_definition = grid_info$grid_def,
    n_cells = grid_info$n_cells,
    wet_cells = nrow(grid_info$wet_grid),
    source_cells = source_dt,
    slice_table = slice_table,
    scenarios = scenarios,
    execution = execution,
    output_files = list(
      daily_sparse_edges_csv = if (isTRUE(params$save_daily_sparse_edges_csv)) edge_csv_path else NA_character_,
      daily_cell_water_budget_csv = if (isTRUE(params$save_daily_cell_budget_csv)) budget_csv_path else NA_character_,
      day_qc_csv = if (isTRUE(params$save_day_qc_csv)) qc_csv_path else NA_character_,
      particle_endpoints_csv = if (isTRUE(params$save_particle_endpoints)) endpoints_csv_path else NA_character_,
      cell_lookup_csv = if (isTRUE(params$save_cell_lookup)) lookup_path else NA_character_,
      daily_matrix_manifest_csv = if (isTRUE(params$save_daily_transition_matrices)) matrix_manifest_path else NA_character_,
      mean_sparse_edges_csv = mean_outputs$mean_edges_file,
      mean_transition_matrix_rds = mean_outputs$mean_matrix_files,
      maps = map_paths,
      run_log = log_path
    ),
    package_versions = list(ncdf4 = as.character(utils::packageVersion("ncdf4")),
                            data.table = as.character(utils::packageVersion("data.table")),
                            Matrix = as.character(utils::packageVersion("Matrix")))
  )

  if (isTRUE(params$save_metadata)) {
    safe_save_rds(metadata, metadata_path, params$overwrite_existing, compress = TRUE)
    log_message("Wrote metadata: ", metadata_path, ".")
  }

  summary_lines <- c(
    "Lagrangian sparse connectivity run summary",
    paste0("Created at: ", Sys.time()),
    paste0("Input directory: ", params$input_dir),
    paste0("Product: ", params$product_id),
    paste0("Files: ", paste(basename(files$full_path), collapse = "; ")),
    paste0("Velocity modes: ", paste(params$velocity_modes, collapse = ", ")),
    paste0("Kh values: ", paste(params$Kh_values, collapse = ", ")),
    paste0("Tracking horizons days: ", paste(params$tracking_horizons_days, collapse = ", ")),
    paste0("dt_seconds: ", params$dt_seconds),
    paste0("Particles per source cell: ", params$N_particles_per_source),
    paste0("Source selection: ", params$source_selection),
    paste0("Wet cells: ", nrow(grid_info$wet_grid)),
    paste0("Source cells: ", nrow(source_dt)),
    paste0("Daily sparse edge CSV: ", if (isTRUE(params$save_daily_sparse_edges_csv)) edge_csv_path else NA_character_),
    paste0("Daily cell water-budget CSV: ", if (isTRUE(params$save_daily_cell_budget_csv)) budget_csv_path else NA_character_),
    paste0("Mean sparse edge CSV: ", mean_outputs$mean_edges_file),
    paste0("Daily matrix manifest: ", if (isTRUE(params$save_daily_transition_matrices)) matrix_manifest_path else NA_character_),
    paste0("Map PNG files: ", paste(map_paths, collapse = "; ")),
    paste0("Runtime seconds: ", round(as.numeric(runtime), 1)),
    paste0("Trajectory reuse: ", params$trajectory_reuse),
    paste0("Parallel split: ", params$parallel_split),
    paste0("Peak velocity cache days: ", execution$peak_velocity_days),
    "",
    "Daily cell water-budget CSV columns to use:",
    "  water_removed_fraction_of_released = fraction of water released from this cell that is no longer retained in this cell after the tracking horizon.",
    "  transported_out_fraction_of_released = fraction that moved from this source cell to another wet cell.",
    "  external_loss_fraction_of_released = fraction that left the wet grid/domain or hit invalid velocity.",
    "  incoming_added_source_equiv_fraction = sum of incoming source-cell probabilities arriving in this cell from other cells.",
    "  incoming_added_fraction_of_total_released = incoming particles divided by all particles released in that simulation.",
    "",
    "Sparse matrix convention:",
    "  P[i, j] = probability that particles released from source cell i end in sink cell j.",
    "  Row sums can be < 1 because particles can leave the wet grid/domain.",
    "",
    "Detailed log:",
    log_lines
  )
  writeLines(summary_lines, log_path)
  log_message("Workflow complete. Runtime: ", round(as.numeric(runtime), 1), " seconds.")

  invisible(list(metadata = metadata,
                 day_qc = day_qc,
                 source_cells = source_dt))
}

# Set options(lagrangian.skip_main = TRUE) before source() for function-only use.
if (!isTRUE(getOption("lagrangian.skip_main", FALSE))) main()
