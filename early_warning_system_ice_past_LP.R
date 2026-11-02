#!/usr/bin/env Rscript

# ------------------------------- Packages -------------------------------
REQUIRED_PACKAGES <- c("terra", "ncdf4", "readxl", "data.table")
OPTIONAL_PACKAGES <- c("sf", "maps")

missing_required_packages <- REQUIRED_PACKAGES[
  !vapply(REQUIRED_PACKAGES, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)
]
if (length(missing_required_packages) > 0L) {
  stop(
    paste0(
      "Missing required packages: ",
      paste(missing_required_packages, collapse = ", "),
      ". Install them before running this script."
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(terra)
  library(ncdf4)
  library(readxl)
  library(data.table)
})

OPTIONAL_PACKAGE_AVAILABLE <- setNames(
  vapply(OPTIONAL_PACKAGES, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE),
  OPTIONAL_PACKAGES
)

# ------------------------------- Parameters -------------------------------
options(stringsAsFactors = FALSE)
Sys.setenv(TZ = "UTC")

# Default parameter values.
# Relative paths below are resolved against the directory where this script is stored.
# Use forward slashes in Windows paths, for example "Z:/Data/MSI_models/early warning system ice".
# nc_dir is selected from the hc/ref or hc/owfp folder according to use_windpark_influence unless overridden.
# Set nc_years to a year such as "2011", or to multiple years such as "2011;2012". Blank means all parseable years.
# Optional nc_start_date and nc_end_date use YYYY-MM-DD and include both boundary dates.
# They filter the existing filename dates together with nc_years; leave nc_years blank for a date-only filter.
# Either date may be blank. Example: --nc_years="" --nc_start_date="2010-11-01" --nc_end_date="2011-04-30".
DEFAULT_PARAMETERS <- list(
  use_windpark_influence = FALSE,
  nc_dir = "Z:/Data/MSI_models/MARTE/hc/ref",
  nc_years = "",
  nc_start_date = "2011-01-01",
  nc_end_date = "2010-02-01",
  inventory_path = "nc_inventory_with_metadata.xlsx",
  config_path = "alert_config_EST05_ice_LP.csv",
  output_dir = "output",
  cache_path = "product_structure_cache_EST05_ice_LP.rds",
  pw4b_grid_path = "//export.hpc.ut.ee/mifiles/Data/planwise4blue/PW4B_EU/GRID/eu_1km_with_est_country_id_72.csv",
  pw4b_lookup_path = "lookup_nc_to_pw4b_grid_EST05_ice_LP.rds",
  write_pw4b = TRUE,
  refresh_pw4b_lookup = FALSE,
  pw4b_crs = "EPSG:3035",
  pw4b_buffer_m = 20000,
  pw4b_nn_distance_m = 500,
  pw4b_idw_k = 1,
  write_png = TRUE,
  refresh_structure = TRUE,
  characterize_only = FALSE,
  overwrite = TRUE,
  window_hours = 24,
  max_window_memory_mb = 2048,
  terra_tempdir = "terra_tmp",
  allow_curvilinear_resample = TRUE,
  include_past_files = TRUE,
  gdal_write_options = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=4"),
  preferred_crs = "EPSG:4326",
  run_main_when_sourced = TRUE
)

# This modified pipeline intentionally accepts only the EST05 sea-ice
# Lagrangian-particle diagnostic products represented by files such as
# EST05_ice_LP_1h_run1.20110101.nc4. The required diagnostic fields are
# validated before processing so older native ice-variable files are skipped.
SUPPORTED_LP_REQUIRED_DATA_VARS <- c("P0", "P1", "P2", "P3", "Pt", "Pc", "Pc1")
SUPPORTED_LP_FILENAME_REGEX <- "^EST05_ice_LP_.*\\.[0-9]{8}([0-9]{2})?\\.nc(4|\\.dap\\.nc4)?$"

# operational_alert_pipeline_EST05_ice_LP_only.R
# Operational alert analysis pipeline for hourly EST05 sea-ice Lagrangian-particle
# diagnostic NetCDF files.
#
# This version is restricted to the linked LP product type, for example:
#   EST05_ice_LP_1h_run1.20110101.nc4
#
# It validates the product before processing and expects these data variables:
#   P0, P1, P2, P3, Pt, Pc, and Pc1
#
# The script keeps the same output families as the original pipeline:
#   - threshold severity GeoTIFFs
#   - change severity GeoTIFFs
#   - combined severity GeoTIFFs
#   - max absolute change GeoTIFFs
#   - optional PNG debug maps
#   - optional PW4B rasters and combined-severity CSVs
#   - output_manifest.csv
#
# Default behavior:
#   * Look for alert_config_EST05_ice_LP.csv next to this script
#   * Look for nc_inventory_with_metadata.xlsx next to this script, if available
#   * Write outputs to ./output
#   * Cache structure in ./product_structure_cache_EST05_ice_LP.rds
#
# Example:
#   Rscript operational_alert_pipeline_EST05_ice_LP_only.R \
#     --nc="Z:/Data/MSI_models/MARTE/hc/ref" \
#     --nc_years="2011" \
#     --config="alert_config_EST05_ice_LP.csv" \
#     --out="output" \
#     --png=TRUE \
#     --refresh_structure=TRUE
#
# Notes:
#   * P0 = stopped Lagrangian particle count per grid cell.
#   * P1 = active Lagrangian particle count per grid cell.
#   * P2 = reserved Lagrangian particle diagnostic field.
#   * P3 = out-of-domain Lagrangian particle count per grid cell.
#   * Pt = integrated active-particle residence time per grid cell.
#   * Pc = integrated active-particle impact proxy per grid cell.
#   * Pc1 = active-particle sample count used for Pt/Pc accumulation.
#   * These LP files are handled as no-depth products.
#   * Thresholds in the default config are starter values and should be adjusted
#     for the intended operational warning logic.


# ---------------------------- Utility helpers ----------------------------

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) {
    y
  } else {
    x
  }
}

is_blank_scalar <- function(x) {
  length(x) == 0L || is.null(x) || all(is.na(x)) || all(trimws(as.character(x)) == "")
}

first_non_blank <- function(...) {
  vals <- list(...)
  for (v in vals) {
    if (!is_blank_scalar(v)) {
      return(v[[1]])
    }
  }
  NA
}

scalar_chr <- function(x, default = NA_character_) {
  if (is_blank_scalar(x)) {
    return(default)
  }
  as.character(x[[1]])
}

scalar_num <- function(x, default = NA_real_) {
  if (is_blank_scalar(x)) {
    return(default)
  }
  suppressWarnings({
    out <- as.numeric(x[[1]])
  })
  if (length(out) == 0L || is.na(out)) default else out
}

scalar_int <- function(x, default = NA_integer_) {
  out <- scalar_num(x, default = default)
  if (is.na(out)) default else as.integer(round(out))
}

scalar_lgl <- function(x, default = FALSE) {
  if (is.logical(x) && length(x) == 1L && !is.na(x)) {
    return(x)
  }
  if (is_blank_scalar(x)) {
    return(default)
  }
  x_chr <- tolower(trimws(as.character(x[[1]])))
  if (x_chr %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (x_chr %in% c("false", "f", "0", "no", "n")) return(FALSE)
  default
}


blank_to_na_chr <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

safe_dir_create <- function(path) {
  if (is_blank_scalar(path)) {
    stop("Directory path is blank and cannot be created.", call. = FALSE)
  }

  path <- normalizePath(path.expand(path), winslash = "/", mustWork = FALSE)

  if (!dir.exists(path)) {
    ok <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!isTRUE(ok) && !dir.exists(path)) {
      stop(
        paste0("Failed to create directory: ", path, ". Check the path and your write permissions."),
        call. = FALSE
      )
    }
  }

  invisible(path)
}

sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nzchar(x), x, "unnamed")
}

split_semicolon_values <- function(x) {
  if (is_blank_scalar(x)) {
    return(character(0))
  }
  out <- trimws(unlist(strsplit(as.character(x), ";", fixed = TRUE)))
  out[nzchar(out)]
}

split_semicolon_int <- function(x) {
  vals <- split_semicolon_values(x)
  if (length(vals) == 0L) {
    return(integer(0))
  }
  suppressWarnings(as.integer(vals))
}

split_year_values <- function(x) {
  if (is_blank_scalar(x)) {
    return(integer(0))
  }

  vals_chr <- trimws(as.character(x))
  vals_chr <- vals_chr[nzchar(vals_chr)]
  if (length(vals_chr) == 0L) {
    return(integer(0))
  }

  vals <- trimws(unlist(strsplit(paste(vals_chr, collapse = ";"), "[,;[:space:]]+", perl = TRUE)))
  vals <- vals[nzchar(vals)]
  if (length(vals) == 0L) {
    return(integer(0))
  }

  suppressWarnings(out <- as.integer(vals))
  out <- out[!is.na(out)]
  out <- out[out >= 1000L & out <= 9999L]
  sort(unique(out))
}

safe_attr <- function(nc, varid, attname) {
  out <- tryCatch(
    ncdf4::ncatt_get(nc, varid = varid, attname = attname),
    error = function(e) NULL
  )
  if (is.null(out) || isTRUE(out$hasatt) == FALSE || is.null(out$value)) {
    return(NA)
  }
  out$value
}


safe_ncvar_get_by_name <- function(nc, var_name) {
  if (is_blank_scalar(var_name)) {
    return(NULL)
  }

  tryCatch(
    ncdf4::ncvar_get(nc, varid = as.character(var_name[[1]])),
    error = function(e) NULL
  )
}

safe_coord_attr <- function(nc, name, attname) {
  if (is_blank_scalar(name)) {
    return(NA)
  }

  name <- as.character(name[[1]])
  out <- safe_attr(nc, name, attname)
  if (!is_blank_scalar(out)) {
    return(out)
  }

  if (name %in% names(nc$dim) && identical(attname, "units")) {
    dim_units <- nc$dim[[name]]$units
    if (!is_blank_scalar(dim_units)) {
      return(dim_units)
    }
  }

  NA
}

coord_units_from_name <- function(nc, name) {
  scalar_chr(safe_coord_attr(nc, name, "units"), default = NA_character_)
}

read_nc_name_values <- function(nc, name, fallback = NULL) {
  if (is_blank_scalar(name)) {
    return(fallback)
  }

  name <- as.character(name[[1]])
  vals <- safe_ncvar_get_by_name(nc, name)
  if (!is.null(vals)) {
    return(vals)
  }

  if (name %in% names(nc$dim)) {
    dim_vals <- nc$dim[[name]]$vals
    if (!is.null(dim_vals)) {
      return(dim_vals)
    }
  }

  fallback
}

coordinate_tokens_from_var_attrs <- function(nc, var_names = names(nc$var)) {
  tokens <- character(0)
  var_names <- var_names[!is.na(var_names) & nzchar(as.character(var_names))]
  if (length(var_names) == 0L) {
    return(tokens)
  }

  for (vn in var_names) {
    coords <- scalar_chr(safe_attr(nc, vn, "coordinates"), default = "")
    if (!is.na(coords) && nzchar(trimws(coords))) {
      tokens <- c(tokens, strsplit(trimws(coords), "[[:space:]]+", perl = TRUE)[[1]])
    }
  }

  unique(tokens[nzchar(tokens)])
}

coordinate_name_matches_role <- function(name, standard_name = "", units = "", axis = "", wanted = c("latitude", "longitude")) {
  wanted <- match.arg(wanted)
  name_l <- tolower(trimws(scalar_chr(name, default = "")))
  std_l <- tolower(trimws(scalar_chr(standard_name, default = "")))
  units_l <- tolower(trimws(scalar_chr(units, default = "")))
  axis_l <- toupper(trimws(scalar_chr(axis, default = "")))

  if (wanted == "latitude") {
    return(
      name_l %in% c("lat", "latitude", "nav_lat", "_nc4_non_coord_lat", "y") ||
        identical(std_l, "latitude") ||
        grepl("degrees_north", units_l, fixed = TRUE) ||
        axis_l %in% c("Y", "LAT")
    )
  }

  name_l %in% c("lon", "longitude", "nav_lon", "_nc4_non_coord_lon", "x") ||
    identical(std_l, "longitude") ||
    grepl("degrees_east", units_l, fixed = TRUE) ||
    axis_l %in% c("X", "LON")
}




units_from_coord_info <- function(nc, info) {
  if (is.null(info) || is.null(info$name)) {
    return(NA_character_)
  }

  units <- coord_units_from_name(nc, info$name)
  if (!is.na(units) && nzchar(units)) {
    return(units)
  }

  if (identical(info$kind, "dim") && info$name %in% names(nc$dim)) {
    return(scalar_chr(nc$dim[[info$name]]$units))
  }

  NA_character_
}

extract_first_numeric <- function(x) {
  txt <- paste(x, collapse = " ")
  hit <- regexpr("[0-9]+(?:\\.[0-9]+)?", txt, perl = TRUE)
  if (length(hit) == 0L || is.na(hit[1]) || hit[1] < 1L) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(regmatches(txt, hit)))
}

detect_total_available_memory_mb <- function() {
  candidates_mb <- numeric(0)

  if (.Platform$OS.type == "windows") {
    ps_out <- tryCatch(
      system2(
        "powershell",
        c(
          "-NoProfile",
          "-Command",
          "$os = Get-CimInstance Win32_OperatingSystem; [double]$os.TotalVisibleMemorySize / 1024"
        ),
        stdout = TRUE,
        stderr = FALSE
      ),
      error = function(e) character(0)
    )
    ps_mb <- extract_first_numeric(ps_out)
    if (is.finite(ps_mb) && ps_mb > 0) {
      candidates_mb <- c(candidates_mb, ps_mb)
    }

    wmic_out <- tryCatch(
      system2(
        "wmic",
        c("computersystem", "get", "TotalPhysicalMemory", "/value"),
        stdout = TRUE,
        stderr = FALSE
      ),
      error = function(e) character(0)
    )
    wmic_bytes <- extract_first_numeric(wmic_out)
    if (is.finite(wmic_bytes) && wmic_bytes > 0) {
      candidates_mb <- c(candidates_mb, wmic_bytes / 1024^2)
    }
  }

  if (file.exists("/proc/meminfo")) {
    meminfo <- tryCatch(readLines("/proc/meminfo", warn = FALSE), error = function(e) character(0))
    memtotal_line <- grep("^MemTotal:", meminfo, value = TRUE)
    memtotal_kb <- extract_first_numeric(memtotal_line)
    if (is.finite(memtotal_kb) && memtotal_kb > 0) {
      candidates_mb <- c(candidates_mb, memtotal_kb / 1024)
    }
  }

  if (length(candidates_mb) == 0L && nzchar(Sys.which("sysctl"))) {
    sysctl_out <- tryCatch(
      system2("sysctl", c("-n", "hw.memsize"), stdout = TRUE, stderr = FALSE),
      error = function(e) character(0)
    )
    sysctl_bytes <- extract_first_numeric(sysctl_out)
    if (is.finite(sysctl_bytes) && sysctl_bytes > 0) {
      candidates_mb <- c(candidates_mb, sysctl_bytes / 1024^2)
    }
  }

  candidates_mb <- candidates_mb[is.finite(candidates_mb) & candidates_mb > 0]
  if (length(candidates_mb) == 0L) {
    return(NA_real_)
  }
  min(candidates_mb)
}

cap_requested_memory_guard_mb <- function(requested_mb, total_available_memory_mb, max_fraction = 0.25, log_file = NULL) {
  if (!is.finite(total_available_memory_mb) || total_available_memory_mb <= 0) {
    stop_with_log(
      "Unable to determine total available memory, so the 25% memory safety cap cannot be enforced.",
      log_file
    )
  }

  cap_mb <- max(1, floor(total_available_memory_mb * max_fraction))
  if (!is.finite(requested_mb) || is.na(requested_mb) || requested_mb <= 0) {
    return(cap_mb)
  }

  min(requested_mb, cap_mb)
}

log_message <- function(level, message, log_file = NULL) {
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC")
  line <- sprintf("[%s] [%s] %s", stamp, level, message)
  message(line)
  if (!is.null(log_file)) {
    cat(line, file = log_file, sep = "\n", append = TRUE)
  }
  invisible(line)
}

stop_with_log <- function(message, log_file = NULL) {
  log_message("ERROR", message, log_file)
  stop(message, call. = FALSE)
}


is_supported_lagrangian_particle_nc_filename <- function(path) {
  filename <- basename(path)
  grepl(SUPPORTED_LP_FILENAME_REGEX, filename, ignore.case = TRUE, perl = TRUE)
}


validate_supported_lagrangian_particle_nc <- function(nc_path, log_file = NULL, stop_on_error = TRUE, nc = NULL) {
  validation_failed <- function(reason) {
    if (isTRUE(stop_on_error)) {
      stop_with_log(
        paste0(
          "Unsupported NetCDF product for this LP-only script: ",
          normalizePath(nc_path, winslash = "/", mustWork = FALSE),
          ". Reason: ", reason
        ),
        log_file
      )
    }

    log_message(
      "WARN",
      paste0(
        "Skipping unsupported NetCDF file: ",
        normalizePath(nc_path, winslash = "/", mustWork = FALSE),
        ". Reason: ", reason
      ),
      log_file
    )
    FALSE
  }

  if (!is_supported_lagrangian_particle_nc_filename(nc_path)) {
    return(validation_failed(
      paste0(
        "filename does not match the supported EST05_ice_LP pattern '",
        SUPPORTED_LP_FILENAME_REGEX,
        "'."
      )
    ))
  }

  close_nc <- FALSE
  if (is.null(nc)) {
    nc <- tryCatch(
      ncdf4::nc_open(nc_path, readunlim = FALSE),
      error = function(e) e
    )
    if (inherits(nc, "error")) {
      return(validation_failed(paste0("unable to open NetCDF file. ", conditionMessage(nc))))
    }
    close_nc <- TRUE
  }
  if (isTRUE(close_nc)) {
    on.exit(ncdf4::nc_close(nc), add = TRUE)
  }

  nc_vars <- names(nc$var)
  if (is.null(nc_vars)) nc_vars <- character(0)
  nc_vars_lower <- tolower(trimws(nc_vars))
  required_lower <- tolower(SUPPORTED_LP_REQUIRED_DATA_VARS)
  missing_required <- SUPPORTED_LP_REQUIRED_DATA_VARS[!required_lower %in% nc_vars_lower]
  if (length(missing_required) > 0L) {
    return(validation_failed(
      paste0(
        "missing required Lagrangian-particle variables: ",
        paste(missing_required, collapse = ", "),
        "."
      )
    ))
  }

  get_actual_var_name <- function(requested_name) {
    idx <- match(tolower(requested_name), nc_vars_lower)
    if (is.na(idx)) NA_character_ else nc_vars[idx]
  }

  reference_dim_names <- NULL
  reference_dim_lengths <- NULL
  for (requested_var in SUPPORTED_LP_REQUIRED_DATA_VARS) {
    actual_var <- get_actual_var_name(requested_var)
    var_obj <- nc$var[[actual_var]]
    dim_names <- vapply(var_obj$dim, function(d) d$name, FUN.VALUE = character(1))
    dim_lengths <- vapply(var_obj$dim, function(d) as.integer(d$len), FUN.VALUE = integer(1))

    if (length(dim_names) != 3L) {
      return(validation_failed(
        paste0(
          "variable '", requested_var,
          "' must have exactly 3 dimensions: time, lat, lon. Detected dimensions: ",
          paste(dim_names, collapse = " x "),
          "."
        )
      ))
    }

    if (!any(tolower(dim_names) %in% c("time", "time_counter"))) {
      return(validation_failed(
        paste0("variable '", requested_var, "' does not have a time dimension.")
      ))
    }

    if (!all(c("lat", "lon") %in% tolower(dim_names))) {
      return(validation_failed(
        paste0(
          "variable '", requested_var,
          "' does not use the expected lat/lon horizontal dimensions. Detected dimensions: ",
          paste(dim_names, collapse = " x "),
          "."
        )
      ))
    }

    if (is.null(reference_dim_names)) {
      reference_dim_names <- tolower(dim_names)
      reference_dim_lengths <- dim_lengths
    } else if (!identical(reference_dim_names, tolower(dim_names)) || !identical(reference_dim_lengths, dim_lengths)) {
      return(validation_failed(
        paste0(
          "variable '", requested_var,
          "' does not share the same dimension order and lengths as the other required LP variables."
        )
      ))
    }
  }

  has_coord_standard_name <- function(wanted) {
    wanted <- match.arg(wanted, c("latitude", "longitude"))

    var_match <- any(vapply(
      nc_vars,
      function(vn) {
        coordinate_name_matches_role(
          name = vn,
          standard_name = safe_coord_attr(nc, vn, "standard_name"),
          units = coord_units_from_name(nc, vn),
          axis = safe_coord_attr(nc, vn, "axis"),
          wanted = wanted
        )
      },
      FUN.VALUE = logical(1)
    ))

    dim_names_nc <- names(nc$dim)
    if (is.null(dim_names_nc)) dim_names_nc <- character(0)
    dim_match <- any(vapply(
      dim_names_nc,
      function(dn) {
        dim_obj <- nc$dim[[dn]]
        coordinate_name_matches_role(
          name = dn,
          standard_name = safe_coord_attr(nc, dn, "standard_name"),
          units = first_non_blank(dim_obj$units, coord_units_from_name(nc, dn)),
          axis = safe_coord_attr(nc, dn, "axis"),
          wanted = wanted
        )
      },
      FUN.VALUE = logical(1)
    ))

    actual_required_vars <- vapply(
      SUPPORTED_LP_REQUIRED_DATA_VARS,
      get_actual_var_name,
      FUN.VALUE = character(1)
    )
    actual_required_vars <- actual_required_vars[!is.na(actual_required_vars)]
    coord_tokens <- coordinate_tokens_from_var_attrs(nc, var_names = actual_required_vars)
    token_match <- any(vapply(
      coord_tokens,
      coordinate_name_matches_role,
      FUN.VALUE = logical(1),
      wanted = wanted
    ))

    var_match || dim_match || token_match
  }

  if (!has_coord_standard_name("latitude")) {
    return(validation_failed("no latitude coordinate variable, coordinate dimension, or coordinates='lat ...' reference was found."))
  }
  if (!has_coord_standard_name("longitude")) {
    return(validation_failed("no longitude coordinate variable, coordinate dimension, _nc4_non_coord_lon variable, or coordinates='... lon' reference was found."))
  }

  TRUE
}

safe_pmax <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  out <- pmax(a, b, na.rm = TRUE)
  if (is.null(dim(out))) {
    template_dims <- dim(a)
    if (is.null(template_dims)) template_dims <- dim(b)
    if (!is.null(template_dims) && length(out) == prod(template_dims)) {
      dim(out) <- template_dims
    }
  }
  both_na <- is.na(a) & is.na(b)
  out[both_na] <- NA_real_
  out
}

safe_pmin <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  out <- pmin(a, b, na.rm = TRUE)
  if (is.null(dim(out))) {
    template_dims <- dim(a)
    if (is.null(template_dims)) template_dims <- dim(b)
    if (!is.null(template_dims) && length(out) == prod(template_dims)) {
      dim(out) <- template_dims
    }
  }
  both_na <- is.na(a) & is.na(b)
  out[both_na] <- NA_real_
  out
}

coerce_xy_matrix <- function(x, nx = NULL, ny = NULL, label = "x") {
  dx <- dim(x)

  if (is.matrix(x)) {
    return(x)
  }

  if (!is.null(dx) && length(dx) == 2L) {
    return(x)
  }

  if (!is.null(dx) && length(dx) > 2L) {
    extra_dims <- dx[-c(1, 2)]
    if (all(extra_dims == 1L, na.rm = TRUE)) {
      return(matrix(as.vector(x), nrow = dx[1], ncol = dx[2]))
    }
  }

  if (is.null(dx) && !is.null(nx) && !is.null(ny) && length(x) == nx * ny) {
    return(matrix(as.vector(x), nrow = nx, ncol = ny))
  }

  stop(
    paste0(
      "Object '", label, "' could not be interpreted as an x/y matrix. ",
      "Returned dimensions: ",
      if (is.null(dx)) "<none>" else paste(dx, collapse = " x "),
      "; length=",
      length(x),
      if (!is.null(nx) && !is.null(ny)) paste0("; expected nx x ny = ", nx, " x ", ny) else ""
    ),
    call. = FALSE
  )
}

nonmissing_mean_pair <- function(sum_mat, count_mat, new_mat) {
  keep <- !is.na(new_mat)
  sum_mat[keep] <- ifelse(is.na(sum_mat[keep]), 0, sum_mat[keep]) + new_mat[keep]
  count_mat[keep] <- ifelse(is.na(count_mat[keep]), 0, count_mat[keep]) + 1
  list(sum = sum_mat, count = count_mat)
}

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(normalizePath(dirname(sub("^--file=", "", file_arg[1])), winslash = "/", mustWork = FALSE))
  }

  frame_ofiles <- vapply(
    sys.frames(),
    function(fr) {
      ofile <- fr$ofile
      if (is.null(ofile) || !nzchar(ofile)) "" else as.character(ofile)
    },
    FUN.VALUE = character(1)
  )
  frame_ofiles <- frame_ofiles[nzchar(frame_ofiles)]
  if (length(frame_ofiles) > 0L) {
    return(normalizePath(dirname(frame_ofiles[[length(frame_ofiles)]]), winslash = "/", mustWork = FALSE))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    editor_path <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(editor_path)) {
      return(normalizePath(dirname(editor_path), winslash = "/", mustWork = FALSE))
    }
  }

  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

is_absolute_path <- function(path) {
  grepl("^([A-Za-z]:[/\\\\]|[/\\\\]{2}|~|/)", path.expand(path))
}

resolve_path_against_base <- function(path, base_dir) {
  if (is_blank_scalar(path)) {
    return(path)
  }

  path <- path.expand(path)
  if (is_absolute_path(path)) {
    return(path)
  }

  file.path(base_dir, path)
}


extract_forecast_date_from_nc_filename <- function(path) {
  filename <- basename(path)
  # EST05 LP filenames may use either YYYYMMDD, e.g. EST05_ice_LP_1h_run1.20260101.nc4,
  # or YYYYMMDDHH. Prefer the last date-like stamp in the filename.
  matches <- regmatches(filename, gregexpr("[0-9]{8}([0-9]{2})?", filename, perl = TRUE))[[1]]
  if (length(matches) == 0L) {
    return(as.Date(NA_character_))
  }

  stamp <- matches[length(matches)]
  date_format <- if (nchar(stamp) == 8L) "%Y%m%d" else "%Y%m%d%H"
  parsed <- as.POSIXct(stamp, format = date_format, tz = "UTC")
  if (is.na(parsed)) {
    return(as.Date(NA_character_))
  }

  as.Date(parsed, tz = "UTC")
}

parse_nc_date_filter <- function(x, argument_name, log_file = NULL) {
  if (is_blank_scalar(x)) {
    return(as.Date(NA_character_))
  }
  value <- trimws(as.character(x))
  if (length(value) != 1L || !grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", value)) {
    stop_with_log(paste0("--", argument_name, " must be a valid date in YYYY-MM-DD format, or blank."), log_file)
  }
  parsed <- suppressWarnings(as.Date(value, format = "%Y-%m-%d"))
  if (is.na(parsed) || format(parsed, "%Y-%m-%d") != value) {
    stop_with_log(paste0("--", argument_name, " must be a valid date in YYYY-MM-DD format, or blank."), log_file)
  }
  parsed
}

list_current_or_future_nc_files <- function(nc_dir, include_past_files = FALSE, nc_years = integer(0), log_file = NULL,
                                            nc_start_date = "", nc_end_date = "") {
  if (is_blank_scalar(nc_dir) || !dir.exists(nc_dir)) {
    stop_with_log(
      "No NetCDF input directory was supplied via --nc_dir, or the directory does not exist.",
      log_file
    )
  }

  start_date <- parse_nc_date_filter(nc_start_date, "nc_start_date", log_file)
  end_date <- parse_nc_date_filter(nc_end_date, "nc_end_date", log_file)
  if (!is.na(start_date) && !is.na(end_date) && start_date > end_date) {
    stop_with_log("--nc_start_date must be on or before --nc_end_date.", log_file)
  }
  date_filter_label <- NULL
  if (!is.na(start_date) || !is.na(end_date)) {
    date_filter_label <- paste0(
      if (is.na(start_date)) "no lower limit" else format(start_date, "%Y-%m-%d"),
      " to ",
      if (is.na(end_date)) "no upper limit" else format(end_date, "%Y-%m-%d"),
      " (inclusive)"
    )
  }

  nc_files <- list.files(
    nc_dir,
    pattern = "\\.nc4?$",
    full.names = TRUE,
    recursive = FALSE,
    ignore.case = TRUE
  )

  nc_files <- nc_files[!endsWith(basename(nc_files), "nc.nc4")]

  filename_ok <- vapply(nc_files, is_supported_lagrangian_particle_nc_filename, FUN.VALUE = logical(1))
  skipped_by_name <- nc_files[!filename_ok]
  if (length(skipped_by_name) > 0L) {
    log_message(
      "WARN",
      paste0(
        "Skipping ", length(skipped_by_name),
        " NetCDF file(s) because their names do not match the supported EST05_ice_LP product pattern: ",
        paste(basename(skipped_by_name), collapse = ", ")
      ),
      log_file
    )
  }
  nc_files <- nc_files[filename_ok]

  if (length(nc_files) == 0L) {
    stop_with_log(
      paste0(
        "No supported EST05_ice_LP .nc or .nc4 files were found in input directory: ",
        nc_dir
      ),
      log_file
    )
  }

  forecast_dates <- as.Date(
    vapply(
      nc_files,
      function(path) {
        d <- extract_forecast_date_from_nc_filename(path)
        if (is.na(d)) NA_character_ else format(d, "%Y-%m-%d")
      },
      FUN.VALUE = character(1)
    ),
    format = "%Y-%m-%d"
  )

  today_utc <- as.Date(Sys.time(), tz = "UTC")
  if (isTRUE(include_past_files)) {
    keep <- !is.na(forecast_dates)
  } else {
    keep <- !is.na(forecast_dates) & forecast_dates >= today_utc
  }

  year_filter <- split_year_values(nc_years)
  if (length(year_filter) > 0L) {
    forecast_years <- suppressWarnings(as.integer(format(forecast_dates, "%Y")))
    keep <- keep & !is.na(forecast_years) & forecast_years %in% year_filter

    log_message(
      "INFO",
      paste0(
        "Restricting NetCDF file selection to forecast year(s): ",
        paste(year_filter, collapse = ", "),
        "."
      ),
      log_file
    )
  }

  if (!is.na(start_date)) {
    keep <- keep & !is.na(forecast_dates) & forecast_dates >= start_date
  }
  if (!is.na(end_date)) {
    keep <- keep & !is.na(forecast_dates) & forecast_dates <= end_date
  }
  if (!is.null(date_filter_label)) {
    log_message(
      "INFO",
      paste0("Restricting NetCDF file selection to filename dates: ", date_filter_label, "."),
      log_file
    )
  }

  skipped_unparseable <- nc_files[is.na(forecast_dates)]
  if (length(skipped_unparseable) > 0L) {
    log_message(
      "WARN",
      paste0(
        "Skipping NetCDF files whose forecast date could not be parsed from a YYYYMMDD or YYYYMMDDHH filename stamp: ",
        paste(basename(skipped_unparseable), collapse = ", ")
      ),
      log_file
    )
  }

  skipped_past <- nc_files[!is.na(forecast_dates) & forecast_dates < today_utc]
  if (length(skipped_past) > 0L) {
    if (isTRUE(include_past_files)) {
      log_message(
        "INFO",
        paste0(
          "include_past_files=TRUE; retaining ", length(skipped_past),
          " past NetCDF file(s), including dated test/reference files."
        ),
        log_file
      )
    } else {
      log_message(
        "INFO",
        paste0(
          "Skipping ", length(skipped_past),
          " past NetCDF file(s) with forecast dates before ", format(today_utc, "%Y-%m-%d"), "."
        ),
        log_file
      )
    }
  }

  if (length(year_filter) > 0L) {
    skipped_year <- nc_files[
      !is.na(forecast_dates) &
        !(suppressWarnings(as.integer(format(forecast_dates, "%Y"))) %in% year_filter)
    ]
    if (length(skipped_year) > 0L) {
      log_message(
        "INFO",
        paste0(
          "Skipping ", length(skipped_year),
          " NetCDF file(s) outside the selected year filter."
        ),
        log_file
      )
    }
  }

  selected_files <- nc_files[keep]
  selected_dates <- forecast_dates[keep]
  if (length(selected_files) == 0L) {
    selection_label <- if (isTRUE(include_past_files)) "parseable" else "current-day or future"
    if (length(year_filter) > 0L) {
      selection_label <- paste0(selection_label, " in selected year(s) ", paste(year_filter, collapse = ", "))
    }
    if (!is.null(date_filter_label)) {
      selection_label <- paste0(selection_label, " within filename date range ", date_filter_label)
    }
    stop_with_log(
      paste0(
        "No ", selection_label, " supported EST05_ice_LP NetCDF files were found in input directory: ",
        nc_dir,
        ". Current day is ", format(today_utc, "%Y-%m-%d"), " UTC."
      ),
      log_file
    )
  }

  type_ok <- vapply(
    selected_files,
    validate_supported_lagrangian_particle_nc,
    FUN.VALUE = logical(1),
    log_file = log_file,
    stop_on_error = FALSE
  )
  selected_files <- selected_files[type_ok]
  selected_dates <- selected_dates[type_ok]

  if (length(selected_files) == 0L) {
    stop_with_log(
      paste0(
        "All date-selected NetCDF files were rejected because they do not match the required EST05_ice_LP particle-diagnostic structure: ",
        paste(SUPPORTED_LP_REQUIRED_DATA_VARS, collapse = ", ")
      ),
      log_file
    )
  }

  ord <- order(selected_dates, basename(selected_files))
  selected_files <- selected_files[ord]

  selection_label <- if (isTRUE(include_past_files)) "date-parseable" else "current-day/future"
  if (length(year_filter) > 0L) {
    selection_label <- paste0(selection_label, " year-filtered")
  }
  if (!is.null(date_filter_label)) {
    selection_label <- paste0(selection_label, " date-filtered")
  }
  log_message(
    "INFO",
    paste0(
      "Selected ", length(selected_files), " ", selection_label,
      " NetCDF file(s) from input directory."
    ),
    log_file
  )

  normalizePath(selected_files, winslash = "/", mustWork = FALSE)
}

parse_cli_args <- function(args) {
  out <- list()
  if (length(args) == 0L) return(out)
  for (arg in args) {
    if (grepl("^--", arg)) {
      key_val <- sub("^--", "", arg)
      if (grepl("=", key_val, fixed = TRUE)) {
        key <- sub("=.*$", "", key_val)
        val <- sub("^[^=]*=", "", key_val)
      } else {
        key <- key_val
        val <- "TRUE"
      }
      out[[key]] <- val
    }
  }
  out
}

warn_missing_optional_packages <- function(write_png = TRUE) {
  if (!isTRUE(write_png)) {
    return(invisible(TRUE))
  }

  missing_optional <- names(OPTIONAL_PACKAGE_AVAILABLE)[!OPTIONAL_PACKAGE_AVAILABLE]
  if (length(missing_optional) > 0L) {
    warning(
      paste0(
        "Optional packages missing for coastline overlays / PNGs: ",
        paste(missing_optional, collapse = ", "),
        ". GeoTIFF outputs will still work."
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# ------------------------------ Defaults ------------------------------

default_settings <- function(script_dir) {
  out <- DEFAULT_PARAMETERS
  out$nc_dir <- if (isTRUE(out$use_windpark_influence)) "Z:/Data/MSI_models/MARTE/hc/owfp" else "Z:/Data/MSI_models/MARTE/hc/ref"
  out$inventory_path <- file.path(script_dir, out$inventory_path)
  out$config_path <- file.path(script_dir, out$config_path)
  out$output_dir <- file.path(script_dir, out$output_dir)
  out$cache_path <- file.path(script_dir, out$cache_path)
  out$pw4b_lookup_path <- file.path(script_dir, out$pw4b_lookup_path)
  out$terra_tempdir <- file.path(script_dir, out$terra_tempdir)
  out
}

merge_settings <- function(defaults, cli_args) {
  out <- defaults
  if (length(cli_args) == 0L) return(out)

  if (!is.null(cli_args$use_windpark_influence)) {
    out$use_windpark_influence <- scalar_lgl(cli_args$use_windpark_influence, out$use_windpark_influence)
    out$nc_dir <- if (isTRUE(out$use_windpark_influence)) "Z:/Data/MSI_models/MARTE/hc/owfp" else "Z:/Data/MSI_models/MARTE/hc/ref"
  }

  if (!is.null(cli_args$nc)) out$nc_dir <- cli_args$nc
  if (!is.null(cli_args$nc_dir)) out$nc_dir <- cli_args$nc_dir
  if (!is.null(cli_args$years)) out$nc_years <- cli_args$years
  if (!is.null(cli_args$nc_years)) out$nc_years <- cli_args$nc_years
  if (!is.null(cli_args$prediction_years)) out$nc_years <- cli_args$prediction_years
  if (!is.null(cli_args$selected_years)) out$nc_years <- cli_args$selected_years
  if (!is.null(cli_args$nc_start_date)) out$nc_start_date <- cli_args$nc_start_date
  if (!is.null(cli_args$nc_end_date)) out$nc_end_date <- cli_args$nc_end_date
  if (!is.null(cli_args$inventory)) out$inventory_path <- cli_args$inventory
  if (!is.null(cli_args$config)) out$config_path <- cli_args$config
  if (!is.null(cli_args$out)) out$output_dir <- cli_args$out
  if (!is.null(cli_args$cache)) out$cache_path <- cli_args$cache
  if (!is.null(cli_args$pw4b_grid)) out$pw4b_grid_path <- cli_args$pw4b_grid
  if (!is.null(cli_args$pw4b_grid_path)) out$pw4b_grid_path <- cli_args$pw4b_grid_path
  if (!is.null(cli_args$pw4b_lookup)) out$pw4b_lookup_path <- cli_args$pw4b_lookup
  if (!is.null(cli_args$pw4b_lookup_path)) out$pw4b_lookup_path <- cli_args$pw4b_lookup_path
  if (!is.null(cli_args$pw4b)) out$write_pw4b <- scalar_lgl(cli_args$pw4b, out$write_pw4b)
  if (!is.null(cli_args$write_pw4b)) out$write_pw4b <- scalar_lgl(cli_args$write_pw4b, out$write_pw4b)
  if (!is.null(cli_args$refresh_pw4b_lookup)) out$refresh_pw4b_lookup <- scalar_lgl(cli_args$refresh_pw4b_lookup, out$refresh_pw4b_lookup)
  if (!is.null(cli_args$pw4b_crs)) out$pw4b_crs <- cli_args$pw4b_crs
  if (!is.null(cli_args$pw4b_buffer_m)) out$pw4b_buffer_m <- scalar_num(cli_args$pw4b_buffer_m, out$pw4b_buffer_m)
  if (!is.null(cli_args$pw4b_nn_distance_m)) out$pw4b_nn_distance_m <- scalar_num(cli_args$pw4b_nn_distance_m, out$pw4b_nn_distance_m)
  if (!is.null(cli_args$pw4b_idw_k)) out$pw4b_idw_k <- scalar_int(cli_args$pw4b_idw_k, out$pw4b_idw_k)
  if (!is.null(cli_args$png)) out$write_png <- scalar_lgl(cli_args$png, out$write_png)
  if (!is.null(cli_args$refresh_structure)) out$refresh_structure <- scalar_lgl(cli_args$refresh_structure, out$refresh_structure)
  if (!is.null(cli_args$characterize_only)) out$characterize_only <- scalar_lgl(cli_args$characterize_only, out$characterize_only)
  if (!is.null(cli_args$overwrite)) out$overwrite <- scalar_lgl(cli_args$overwrite, out$overwrite)
  if (!is.null(cli_args$window_hours)) out$window_hours <- scalar_int(cli_args$window_hours, out$window_hours)
  if (!is.null(cli_args$max_window_memory_mb)) out$max_window_memory_mb <- scalar_num(cli_args$max_window_memory_mb, out$max_window_memory_mb)
  if (!is.null(cli_args$terra_tempdir)) out$terra_tempdir <- cli_args$terra_tempdir
  if (!is.null(cli_args$allow_curvilinear_resample)) out$allow_curvilinear_resample <- scalar_lgl(cli_args$allow_curvilinear_resample, out$allow_curvilinear_resample)
  if (!is.null(cli_args$include_past_files)) out$include_past_files <- scalar_lgl(cli_args$include_past_files, out$include_past_files)
  if (!is.null(cli_args$preferred_crs)) out$preferred_crs <- cli_args$preferred_crs

  out
}

# --------------------------- Configuration ---------------------------

write_default_alert_config <- function(config_path) {
  safe_dir_create(dirname(config_path))
  dt <- data.table::data.table(
    enabled = rep(TRUE, 7L),
    alert_id = c(
      "P0_stopped_particles",
      "P1_active_particles",
      "P2_reserved_particles",
      "P3_out_of_domain_particles",
      "Pt_residence_time",
      "Pc_impact_proxy",
      "Pc1_sample_count"
    ),
    product_var = c("P0", "P1", "P2", "P3", "Pt", "Pc", "Pc1"),
    source_type = rep("native", 7L),
    source_var_1 = c("P0", "P1", "P2", "P3", "Pt", "Pc", "Pc1"),
    source_var_2 = rep("", 7L),
    low_threshold = rep(NA_real_, 7L),
    high_threshold = c(0, 0, 0, 0, 3600, 1, 0),
    change_threshold = c(1, 1, 1, 1, 3600, 1, 1),
    change_days = rep(1, 7L),
    change_direction = rep("abs", 7L),
    depth_mode = rep("surface", 7L),
    depth_indices = rep("", 7L),
    summary_fun = rep("", 7L),
    write_png = rep(TRUE, 7L),
    notes = c(
      "Stopped Lagrangian particle count per grid cell. Starter high_threshold=0 flags any stopped particles; adjust operationally.",
      "Active Lagrangian particle count per grid cell. Starter high_threshold=0 flags any active particles; adjust operationally.",
      "Reserved Lagrangian particle diagnostic field. Usually zero in the linked product; adjust or disable if not needed.",
      "Out-of-domain Lagrangian particle count per grid cell. Starter high_threshold=0 flags any out-of-domain particles; adjust operationally.",
      "Integrated active-particle residence time per grid cell. Starter threshold is 3600 seconds; adjust operationally.",
      "Integrated active-particle impact proxy per grid cell. Starter threshold is 1 proxy unit; adjust operationally.",
      "Active-particle sample count used for Pt and Pc accumulation. Starter high_threshold=0 flags any samples; adjust operationally."
    )
  )
  data.table::fwrite(dt, config_path, na = "")
  invisible(config_path)
}

read_alert_config <- function(config_path, log_file = NULL) {
  if (!file.exists(config_path)) {
    write_default_alert_config(config_path)
    stop_with_log(
      paste0(
        "No alert configuration file was found. A template has been created at:\n  ",
        config_path,
        "\nFill in thresholds / depth rules and rerun."
      ),
      log_file
    )
  }

  cfg <- data.table::fread(config_path, na.strings = c("", "NA", "NaN", "null"), showProgress = FALSE)

  required_cols <- c(
    "enabled", "alert_id", "product_var", "source_type", "source_var_1", "source_var_2",
    "low_threshold", "high_threshold", "change_threshold", "change_days",
    "change_direction", "depth_mode", "depth_indices", "summary_fun", "write_png", "notes"
  )
  missing_cols <- setdiff(required_cols, names(cfg))
  if (length(missing_cols) > 0L) {
    stop_with_log(
      paste0("alert_config_EST05_ice_LP.csv is missing required columns: ", paste(missing_cols, collapse = ", ")),
      log_file
    )
  }

  cfg[, enabled := vapply(enabled, scalar_lgl, FUN.VALUE = logical(1), default = FALSE)]
  cfg[, write_png := vapply(write_png, scalar_lgl, FUN.VALUE = logical(1), default = TRUE)]
  cfg[, source_type := blank_to_na_chr(tolower(trimws(source_type)))]
  cfg[, product_var := blank_to_na_chr(product_var)]
  cfg[, source_var_1 := blank_to_na_chr(source_var_1)]
  cfg[, source_var_2 := blank_to_na_chr(source_var_2)]
  cfg[, depth_mode := blank_to_na_chr(tolower(trimws(depth_mode)))]
  cfg[, depth_indices := blank_to_na_chr(depth_indices)]
  cfg[, summary_fun := blank_to_na_chr(tolower(trimws(summary_fun)))]
  cfg[, change_direction := blank_to_na_chr(tolower(trimws(change_direction)))]
  cfg[, alert_id := blank_to_na_chr(alert_id)]

  cfg <- cfg[enabled == TRUE]

  if (nrow(cfg) == 0L) {
    stop_with_log("alert_config_EST05_ice_LP.csv contains no enabled alert definitions.", log_file)
  }

  bad_source_type <- cfg[!source_type %in% c("native", "vector_magnitude")]
  if (nrow(bad_source_type) > 0L) {
    stop_with_log(
      paste0(
        "Invalid source_type values in alert_config_EST05_ice_LP.csv: ",
        paste(unique(bad_source_type$source_type), collapse = ", "),
        ". Allowed values: native, vector_magnitude."
      ),
      log_file
    )
  }

  bad_depth_mode <- cfg[!depth_mode %in% c("surface", "bottom", "selected", "summary", "max_severity")]
  if (nrow(bad_depth_mode) > 0L) {
    stop_with_log(
      paste0(
        "Invalid depth_mode values in alert_config_EST05_ice_LP.csv: ",
        paste(unique(bad_depth_mode$depth_mode), collapse = ", "),
        ". Allowed values: surface, bottom, selected, summary, max_severity."
      ),
      log_file
    )
  }

  cfg
}

# ---------------------- Inventory / metadata mapping ----------------------

guess_inventory_key_column <- function(dt) {
  if (nrow(dt) == 0L) return(NA_character_)
  nn <- tolower(names(dt))
  preferred <- c("variable", "var", "varname", "variable_name", "short_name", "name", "standard_name", "parameter")
  hit <- match(preferred, nn, nomatch = 0L)
  hit <- hit[hit > 0L]
  if (length(hit) > 0L) return(names(dt)[hit[1]])
  NA_character_
}

read_inventory_metadata <- function(inventory_path, log_file = NULL) {
  if (!file.exists(inventory_path)) {
    log_message(
      "WARN",
      paste0(
        "Inventory workbook not found at ", inventory_path,
        ". Continuing with NetCDF-only metadata. Supply nc_inventory_with_metadata.xlsx for the full intended workflow."
      ),
      log_file
    )
    return(data.table::data.table())
  }

  sheets <- readxl::excel_sheets(inventory_path)
  if (length(sheets) == 0L) {
    log_message("WARN", "Inventory workbook has no readable sheets.", log_file)
    return(data.table::data.table())
  }

  parts <- lapply(
    sheets,
    function(sh) {
      dt <- data.table::as.data.table(readxl::read_excel(inventory_path, sheet = sh, guess_max = 10000))
      dt[, inventory_sheet := sh]
      dt
    }
  )
  inv <- data.table::rbindlist(parts, fill = TRUE, use.names = TRUE)
  if (nrow(inv) == 0L) {
    log_message("WARN", "Inventory workbook was read successfully but contains no rows.", log_file)
    return(inv)
  }
  inv
}

read_inventory_file_summary <- function(inventory_path, nc_path, log_file = NULL) {
  if (!file.exists(inventory_path)) {
    return(NULL)
  }

  sheets <- tryCatch(readxl::excel_sheets(inventory_path), error = function(e) character(0))
  if (!"selected_nc_files" %in% sheets) {
    return(NULL)
  }

  dt <- data.table::as.data.table(
    readxl::read_excel(inventory_path, sheet = "selected_nc_files", guess_max = 10000)
  )
  if (nrow(dt) == 0L) {
    return(NULL)
  }

  target_file <- tolower(basename(nc_path))

  hit <- rep(FALSE, nrow(dt))
  if ("file_name" %in% names(dt)) {
    file_vals <- tolower(trimws(as.character(dt$file_name)))
    hit <- hit | (!is.na(file_vals) & nzchar(file_vals) & file_vals == target_file)
  }
  if ("full_path" %in% names(dt)) {
    full_vals <- tolower(trimws(basename(as.character(dt$full_path))))
    hit <- hit | (!is.na(full_vals) & nzchar(full_vals) & full_vals == target_file)
  }

  out <- dt[hit]
  if (nrow(out) == 0L) {
    return(NULL)
  }

  log_message(
    "INFO",
    paste0(
      "Using file-level inventory summary from sheet 'selected_nc_files' for '",
      basename(nc_path),
      "'."
    ),
    log_file
  )

  out[1]
}

read_inventory_variable_metadata <- function(inventory_path, nc_path = NULL, log_file = NULL) {
  if (!file.exists(inventory_path)) {
    return(data.table::data.table())
  }

  sheets <- tryCatch(readxl::excel_sheets(inventory_path), error = function(e) character(0))
  if (!"variables" %in% sheets) {
    return(data.table::data.table())
  }

  dt <- data.table::as.data.table(
    readxl::read_excel(inventory_path, sheet = "variables", guess_max = 10000)
  )
  dt[, inventory_sheet := "variables"]

  if (nrow(dt) == 0L || is.null(nc_path)) {
    return(dt)
  }

  target_file <- tolower(basename(nc_path))
  hit <- rep(FALSE, nrow(dt))

  if ("file_name" %in% names(dt)) {
    file_vals <- tolower(trimws(as.character(dt$file_name)))
    hit <- hit | (!is.na(file_vals) & nzchar(file_vals) & file_vals == target_file)
  }
  if ("full_path" %in% names(dt)) {
    full_vals <- tolower(trimws(basename(as.character(dt$full_path))))
    hit <- hit | (!is.na(full_vals) & nzchar(full_vals) & full_vals == target_file)
  }

  out <- dt[hit]
  if (nrow(out) > 0L) {
    log_message(
      "INFO",
      paste0(
        "Using variable-level inventory metadata from sheet 'variables' for '",
        basename(nc_path),
        "': ",
        nrow(out),
        " rows retained."
      ),
      log_file
    )
    return(out)
  }

  dt
}


filter_inventory_to_nc_file <- function(inventory_dt, nc_path, log_file = NULL) {
  if (nrow(inventory_dt) == 0L) {
    return(inventory_dt)
  }

  target_file <- tolower(basename(nc_path))
  inv <- data.table::copy(inventory_dt)

  match_mask <- rep(FALSE, nrow(inv))

  if ("file_name" %in% names(inv)) {
    file_vals <- tolower(trimws(as.character(inv$file_name)))
    match_mask <- match_mask | (!is.na(file_vals) & nzchar(file_vals) & file_vals == target_file)
  }

  if ("full_path" %in% names(inv)) {
    full_vals <- tolower(trimws(basename(as.character(inv$full_path))))
    match_mask <- match_mask | (!is.na(full_vals) & nzchar(full_vals) & full_vals == target_file)
  }

  if (any(match_mask, na.rm = TRUE)) {
    filtered <- inv[match_mask]
    log_message(
      "INFO",
      paste0(
        "Filtered inventory workbook to rows relevant to file '",
        basename(nc_path),
        "': ",
        nrow(filtered),
        " rows retained."
      ),
      log_file
    )
    return(filtered)
  }

  log_message(
    "WARN",
    paste0(
      "No file-specific rows were found in the inventory workbook for '",
      basename(nc_path),
      "'. Continuing with the full workbook."
    ),
    log_file
  )
  inv
}

harmonize_inventory <- function(variable_dt, inventory_dt) {
  if (nrow(inventory_dt) == 0L) {
    variable_dt[, inventory_match := FALSE]
    variable_dt[, inventory_sheet := NA_character_]
    variable_dt[, inventory_key := NA_character_]
    return(variable_dt)
  }

  key_col <- guess_inventory_key_column(inventory_dt)
  if (is.na(key_col)) {
    variable_dt[, inventory_match := FALSE]
    variable_dt[, inventory_sheet := NA_character_]
    variable_dt[, inventory_key := NA_character_]
    variable_dt[, inventory_note := "No obvious variable-name column found in inventory workbook."]
    return(variable_dt)
  }

  inv <- data.table::copy(inventory_dt)
  vars <- data.table::copy(variable_dt)
  inv[, join_key := tolower(trimws(as.character(get(key_col))))]
  inv <- inv[!is.na(join_key) & nzchar(join_key)]

  inv_first <- inv[
    , .(
      inventory_match = TRUE,
      inventory_sheet = first_non_blank(inventory_sheet),
      inventory_key = first_non_blank(get(key_col))
    ),
    by = join_key
  ]

  vars[, join_key := tolower(trimws(as.character(var_name)))]
  merged <- merge(vars, inv_first, by = "join_key", all.x = TRUE)
  merged[is.na(inventory_match), inventory_match := FALSE]
  merged[is.na(inventory_sheet), inventory_sheet := NA_character_]
  merged[is.na(inventory_key), inventory_key := NA_character_]
  merged[]
}

# ---------------------- NetCDF structure characterization ----------------------

name_matches <- function(x, candidates) {
  x <- tolower(x)
  any(x %in% tolower(candidates))
}

looks_like_role <- function(name, standard_name = NA, long_name = NA, units = NA, axis = NA, role = "x") {
  name_all <- tolower(paste(name, standard_name, long_name, units, axis))
  switch(
    role,
    x = grepl("(^|[^a-z])(lon|longitude|nav_lon|x)([^a-z]|$)", name_all),
    y = grepl("(^|[^a-z])(lat|latitude|nav_lat|y)([^a-z]|$)", name_all),
    z = grepl("(^|[^a-z])(depth|lev|level|z|sigma)([^a-z]|$)", name_all),
    t = grepl("(^|[^a-z])(time|time_counter|t)([^a-z]|$)", name_all),
    FALSE
  )
}


get_reference_horizontal_dims <- function(nc) {
  var_names <- names(nc$var)
  if (is.null(var_names) || length(var_names) == 0L) {
    return(NULL)
  }

  preferred <- intersect(SUPPORTED_LP_REQUIRED_DATA_VARS, var_names)
  search_vars <- unique(c(preferred, var_names))

  for (vn in search_vars) {
    v <- nc$var[[vn]]
    if (is.null(v) || length(v$dim) < 2L) {
      next
    }

    dim_names <- vapply(v$dim, function(d) d$name, FUN.VALUE = character(1))
    dim_lengths <- vapply(v$dim, function(d) as.integer(d$len), FUN.VALUE = integer(1))
    roles <- vapply(dim_names, infer_dim_role_from_name, FUN.VALUE = character(1))

    x_pos <- match("x", roles)
    y_pos <- match("y", roles)
    if (!is.na(x_pos) && !is.na(y_pos)) {
      return(list(
        var_name = vn,
        dim_names = dim_names,
        dim_lengths = dim_lengths,
        x_dim_name = dim_names[x_pos],
        y_dim_name = dim_names[y_pos],
        x_len = dim_lengths[x_pos],
        y_len = dim_lengths[y_pos]
      ))
    }
  }

  NULL
}

infer_coord_dim_names_from_values <- function(nc, values) {
  value_dims <- dim(values)
  if (is.null(value_dims) || length(value_dims) < 2L) {
    return(character(0))
  }

  ref <- get_reference_horizontal_dims(nc)
  if (is.null(ref)) {
    return(character(0))
  }

  first_two <- as.integer(value_dims[1:2])
  yx <- as.integer(c(ref$y_len, ref$x_len))
  xy <- as.integer(c(ref$x_len, ref$y_len))

  if (identical(first_two, yx)) {
    return(c(ref$y_dim_name, ref$x_dim_name))
  }
  if (identical(first_two, xy)) {
    return(c(ref$x_dim_name, ref$y_dim_name))
  }

  character(0)
}


extract_coord_or_dim <- function(nc, role = c("x", "y", "z", "t"), preferred_name = NULL) {
  role <- match.arg(role)
  var_names <- names(nc$var)
  dim_names <- names(nc$dim)
  if (is.null(var_names)) var_names <- character(0)
  if (is.null(dim_names)) dim_names <- character(0)
  available_names <- unique(c(var_names, dim_names))

  case_match <- function(candidate, names_vec) {
    if (is_blank_scalar(candidate) || length(names_vec) == 0L) {
      return(NA_character_)
    }
    idx <- match(tolower(candidate), tolower(names_vec))
    if (is.na(idx)) NA_character_ else names_vec[idx]
  }

  coord_info_for_name <- function(name) {
    name <- as.character(name[[1]])
    matched_var <- case_match(name, var_names)
    matched_dim <- case_match(name, dim_names)
    actual_name <- if (!is.na(matched_var)) matched_var else if (!is.na(matched_dim)) matched_dim else name

    vals <- read_nc_name_values(nc, actual_name)
    kind <- if (!is.na(matched_var)) "var" else "dim"
    dim_names_for_value <- character(0)

    if (!is.na(matched_var)) {
      dim_names_for_value <- vapply(nc$var[[matched_var]]$dim, function(d) d$name, FUN.VALUE = character(1))
    } else if (!is.null(vals) && !is.null(dim(vals)) && length(dim(vals)) >= 2L) {
      dim_names_for_value <- infer_coord_dim_names_from_values(nc, vals)
    }

    list(kind = kind, name = actual_name, values = vals, dim_names = dim_names_for_value)
  }

  if (!is_blank_scalar(preferred_name)) {
    preferred_name <- as.character(preferred_name[[1]])
    preferred_l <- tolower(preferred_name)
    preferred_candidates <- unique(c(
      # For this LP product a generic "lon" inventory entry may refer to a 1D
      # placeholder dimension; the true 2D longitude coordinate is commonly
      # stored as _nc4_non_coord_lon. Try that first.
      if (role == "x" && preferred_l %in% c("lon", "longitude", "x")) "_nc4_non_coord_lon" else character(0),
      if (role == "y" && preferred_l %in% c("lat", "latitude", "y")) "lat" else character(0),
      preferred_name
    ))

    for (candidate in preferred_candidates) {
      matched <- case_match(candidate, available_names)
      if (!is.na(matched)) {
        info <- coord_info_for_name(matched)
        if (!is.null(info$values) && length(info$values) > 0L) {
          return(info)
        }
      }
      vals <- safe_ncvar_get_by_name(nc, candidate)
      if (!is.null(vals)) {
        return(list(kind = "var", name = candidate, values = vals, dim_names = infer_coord_dim_names_from_values(nc, vals)))
      }
    }
  }

  # Prefer explicit EST05 / CF coordinate names. The linked LP product can expose
  # longitude as _nc4_non_coord_lon and latitude as a dimension-scale named lat.
  preferred_coord_names <- switch(
    role,
    x = c("_nc4_non_coord_lon", "lon", "longitude", "nav_lon"),
    y = c("lat", "_nc4_non_coord_lat", "latitude", "nav_lat"),
    t = c("time", "time_counter"),
    z = c("depth", "deptht", "lev", "level", "z")
  )
  for (candidate in preferred_coord_names) {
    matched <- case_match(candidate, available_names)
    if (!is.na(matched)) {
      return(coord_info_for_name(matched))
    }
    vals <- safe_ncvar_get_by_name(nc, candidate)
    if (!is.null(vals)) {
      return(list(kind = "var", name = candidate, values = vals, dim_names = infer_coord_dim_names_from_values(nc, vals)))
    }
  }

  # 1. Search variables first. For horizontal coordinates, ignore 3D/4D
  # data variables such as P0/P1/Pt whose attributes are diagnostics.
  if (length(var_names) > 0L) {
    score <- vapply(
      var_names,
      function(vn) {
        v <- nc$var[[vn]]
        nd <- length(v$dim)
        if (role %in% c("x", "y") && nd > 2L) return(0L)
        if (role %in% c("t", "z") && nd > 1L) return(0L)

        long_name <- scalar_chr(safe_coord_attr(nc, vn, "long_name"))
        standard_name <- scalar_chr(safe_coord_attr(nc, vn, "standard_name"))
        units <- coord_units_from_name(nc, vn)
        axis <- scalar_chr(safe_coord_attr(nc, vn, "axis"))

        score_i <- if (looks_like_role(vn, standard_name, long_name, units, axis, role)) 1L else 0L
        vn_l <- tolower(vn)
        std_l <- tolower(scalar_chr(standard_name, default = ""))
        units_l <- tolower(scalar_chr(units, default = ""))
        axis_l <- toupper(scalar_chr(axis, default = ""))

        if (role == "x" && (
          vn_l %in% c("lon", "longitude", "nav_lon", "_nc4_non_coord_lon") ||
          std_l == "longitude" || grepl("degrees_east", units_l, fixed = TRUE)
        )) score_i <- score_i + 10L
        if (role == "y" && (
          vn_l %in% c("lat", "latitude", "nav_lat", "_nc4_non_coord_lat") ||
          std_l == "latitude" || grepl("degrees_north", units_l, fixed = TRUE)
        )) score_i <- score_i + 10L
        if (role == "t" && (
          vn_l %in% c("time", "time_counter") ||
          std_l == "time" || axis_l == "T"
        )) score_i <- score_i + 10L
        if (role == "z" && (
          vn_l %in% c("depth", "deptht", "lev", "level", "z") ||
          grepl("depth|vertical", std_l) || axis_l == "Z"
        )) score_i <- score_i + 10L

        score_i
      },
      FUN.VALUE = integer(1)
    )
    if (any(score > 0L)) {
      chosen <- var_names[which.max(score)]
      return(coord_info_for_name(chosen))
    }
  }

  # 2. Fall back to dimensions / dimension-scale coordinate entries. This is
  # needed when ncdf4 classifies lat/lon as dimensions rather than variables.
  if (length(dim_names) > 0L) {
    score <- vapply(
      dim_names,
      function(dn) {
        dim_obj <- nc$dim[[dn]]
        long_name <- scalar_chr(safe_coord_attr(nc, dn, "long_name"))
        standard_name <- scalar_chr(safe_coord_attr(nc, dn, "standard_name"))
        units <- first_non_blank(coord_units_from_name(nc, dn), dim_obj$units)
        axis <- scalar_chr(safe_coord_attr(nc, dn, "axis"))

        score_i <- if (looks_like_role(dn, standard_name, long_name, units, axis, role)) 1L else 0L
        dn_l <- tolower(dn)
        std_l <- tolower(scalar_chr(standard_name, default = ""))
        units_l <- tolower(scalar_chr(units, default = ""))
        axis_l <- toupper(scalar_chr(axis, default = ""))

        if (role == "x" && (dn_l %in% c("lon", "longitude", "nav_lon", "x") ||
            std_l == "longitude" || grepl("degrees_east", units_l, fixed = TRUE))) score_i <- score_i + 10L
        if (role == "y" && (dn_l %in% c("lat", "latitude", "nav_lat", "y") ||
            std_l == "latitude" || grepl("degrees_north", units_l, fixed = TRUE))) score_i <- score_i + 10L
        if (role == "t" && (dn_l %in% c("time", "time_counter", "t") ||
            std_l == "time" || axis_l == "T")) score_i <- score_i + 10L
        if (role == "z" && (dn_l %in% c("depth", "deptht", "lev", "level", "z") ||
            grepl("depth|vertical", std_l) || axis_l == "Z")) score_i <- score_i + 10L

        score_i
      },
      FUN.VALUE = integer(1)
    )
    if (any(score > 0L)) {
      chosen <- dim_names[which.max(score)]
      return(coord_info_for_name(chosen))
    }
  }

  NULL
}

infer_dim_role_from_name <- function(dim_name, time_name = NULL, depth_name = NULL, lon_name = NULL, lat_name = NULL) {
  dn <- tolower(dim_name)

  if (!is.null(time_name) && tolower(time_name) == dn) return("t")
  if (!is.null(depth_name) && tolower(depth_name) == dn) return("z")
  if (!is.null(lon_name) && tolower(lon_name) == dn) return("x")
  if (!is.null(lat_name) && tolower(lat_name) == dn) return("y")

  if (looks_like_role(dim_name, role = "t")) return("t")
  if (looks_like_role(dim_name, role = "z")) return("z")
  if (looks_like_role(dim_name, role = "x")) return("x")
  if (looks_like_role(dim_name, role = "y")) return("y")
  NA_character_
}

align_array_to_xy <- function(arr, dim_names, time_name = NULL, depth_name = NULL, lon_name = NULL, lat_name = NULL) {
  if (is.null(dim(arr))) {
    return(arr)
  }

  roles <- vapply(
    dim_names,
    infer_dim_role_from_name,
    FUN.VALUE = character(1),
    time_name = time_name,
    depth_name = depth_name,
    lon_name = lon_name,
    lat_name = lat_name
  )

  x_pos <- match("x", roles)
  y_pos <- match("y", roles)

  if (is.na(x_pos) || is.na(y_pos)) {
    return(arr)
  }

  perm <- c(x_pos, y_pos, setdiff(seq_along(dim(arr)), c(x_pos, y_pos)))
  aperm(arr, perm = perm)
}

detect_crs <- function(nc, lon_info = NULL, lat_info = NULL, preferred_crs = "EPSG:4326") {
  global_candidates <- c(
    scalar_chr(safe_attr(nc, 0, "crs")),
    scalar_chr(safe_attr(nc, 0, "spatial_ref")),
    scalar_chr(safe_attr(nc, 0, "proj4")),
    scalar_chr(safe_attr(nc, 0, "projection"))
  )
  global_candidates <- global_candidates[!is.na(global_candidates) & nzchar(global_candidates)]

  if (length(global_candidates) > 0L) {
    return(global_candidates[1])
  }

  lon_units <- units_from_coord_info(nc, lon_info)
  lat_units <- units_from_coord_info(nc, lat_info)
  unit_blob <- tolower(paste(lon_units, lat_units))

  if (grepl("degree", unit_blob)) {
    return(preferred_crs)
  }

  ""
}

parse_time_values <- function(time_vals, time_units) {
  if (length(time_vals) == 0L) {
    return(list(times = as.POSIXct(character(0), tz = "UTC"), origin = NA_character_, unit = NA_character_))
  }

  if (is_blank_scalar(time_units)) {
    origin <- as.POSIXct("1970-01-01 00:00:00", tz = "UTC")
    return(list(times = origin + as.numeric(time_vals), origin = "1970-01-01 00:00:00 UTC", unit = "seconds"))
  }

  pattern <- "^\\s*(seconds|second|secs|sec|minutes|minute|mins|min|hours|hour|hrs|hr|days|day)\\s+since\\s+(.+)\\s*$"
  match <- regexec(pattern, time_units, ignore.case = TRUE)
  pieces <- regmatches(time_units, match)[[1]]

  if (length(pieces) < 3L) {
    origin <- as.POSIXct("1970-01-01 00:00:00", tz = "UTC")
    return(list(times = origin + as.numeric(time_vals), origin = "1970-01-01 00:00:00 UTC", unit = "seconds"))
  }

  unit_raw <- tolower(pieces[2])
  origin_str <- gsub("T", " ", pieces[3], fixed = TRUE)
  origin_str <- sub("Z$", "", origin_str)

  origin <- suppressWarnings(as.POSIXct(origin_str, tz = "UTC", tryFormats = c(
    "%Y-%m-%d %H:%M:%S",
    "%Y-%m-%d %H:%M",
    "%Y-%m-%d"
  )))

  if (is.na(origin)) {
    origin <- as.POSIXct("1970-01-01 00:00:00", tz = "UTC")
  }

  multiplier <- if (unit_raw %in% c("seconds", "second", "secs", "sec")) {
    1
  } else if (unit_raw %in% c("minutes", "minute", "mins", "min")) {
    60
  } else if (unit_raw %in% c("hours", "hour", "hrs", "hr")) {
    3600
  } else if (unit_raw %in% c("days", "day")) {
    86400
  } else {
    1
  }

  list(
    times = origin + as.numeric(time_vals) * multiplier,
    origin = format(origin, "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"),
    unit = unit_raw
  )
}

compute_extent_from_coords <- function(x_vals, y_vals) {
  x_sorted <- sort(as.numeric(x_vals))
  y_sorted <- sort(as.numeric(y_vals))

  dx <- if (length(x_sorted) > 1L) stats::median(diff(x_sorted), na.rm = TRUE) else 1
  dy <- if (length(y_sorted) > 1L) stats::median(diff(y_sorted), na.rm = TRUE) else 1

  list(
    xmin = min(x_sorted, na.rm = TRUE) - dx / 2,
    xmax = max(x_sorted, na.rm = TRUE) + dx / 2,
    ymin = min(y_sorted, na.rm = TRUE) - dy / 2,
    ymax = max(y_sorted, na.rm = TRUE) + dy / 2,
    dx = dx,
    dy = dy,
    x_sorted = x_sorted,
    y_sorted = y_sorted
  )
}

compute_curvilinear_template <- function(lon_xy, lat_xy) {
  lon_vec <- as.numeric(lon_xy)
  lat_vec <- as.numeric(lat_xy)
  lon_vec <- lon_vec[is.finite(lon_vec)]
  lat_vec <- lat_vec[is.finite(lat_vec)]

  lon_u <- sort(unique(round(lon_vec, 6)))
  lat_u <- sort(unique(round(lat_vec, 6)))
  dx <- suppressWarnings(stats::median(diff(lon_u), na.rm = TRUE))
  dy <- suppressWarnings(stats::median(diff(lat_u), na.rm = TRUE))

  if (!is.finite(dx) || dx <= 0) {
    dx <- (max(lon_vec, na.rm = TRUE) - min(lon_vec, na.rm = TRUE)) / max(1, nrow(lon_xy) - 1)
  }
  if (!is.finite(dy) || dy <= 0) {
    dy <- (max(lat_vec, na.rm = TRUE) - min(lat_vec, na.rm = TRUE)) / max(1, ncol(lat_xy) - 1)
  }

  xmin <- min(lon_vec, na.rm = TRUE) - dx / 2
  xmax <- max(lon_vec, na.rm = TRUE) + dx / 2
  ymin <- min(lat_vec, na.rm = TRUE) - dy / 2
  ymax <- max(lat_vec, na.rm = TRUE) + dy / 2

  list(
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax,
    dx = dx,
    dy = dy,
    ncols = max(1L, as.integer(round((xmax - xmin) / dx))),
    nrows = max(1L, as.integer(round((ymax - ymin) / dy)))
  )
}

characterize_product <- function(nc_path, inventory_path, preferred_crs = "EPSG:4326", log_file = NULL) {
  log_message("INFO", paste0("Characterizing NetCDF structure from: ", nc_path), log_file)

  inventory_summary <- read_inventory_file_summary(inventory_path, nc_path = nc_path, log_file = log_file)
  inventory_var_dt <- read_inventory_variable_metadata(inventory_path, nc_path = nc_path, log_file = log_file)

  preferred_lon_name <- if (!is.null(inventory_summary)) scalar_chr(inventory_summary$lon_name) else NA_character_
  preferred_lat_name <- if (!is.null(inventory_summary)) scalar_chr(inventory_summary$lat_name) else NA_character_
  preferred_depth_name <- if (!is.null(inventory_summary)) scalar_chr(inventory_summary$depth_name) else NA_character_
  preferred_time_name <- if (!is.null(inventory_summary)) scalar_chr(inventory_summary$time_name) else NA_character_
  inventory_grid_type <- if (!is.null(inventory_summary)) {
    sgt <- tolower(scalar_chr(inventory_summary$spatial_grid_type))
    if (!is.na(sgt) && grepl("regular", sgt)) {
      "regular"
    } else if (!is.na(sgt) && grepl("curvilinear", sgt)) {
      "curvilinear"
    } else {
      NA_character_
    }
  } else {
    NA_character_
  }

  nc <- ncdf4::nc_open(nc_path, readunlim = FALSE)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  validate_supported_lagrangian_particle_nc(
    nc_path,
    log_file = log_file,
    stop_on_error = TRUE,
    nc = nc
  )

  lon_info <- extract_coord_or_dim(nc, "x", preferred_name = preferred_lon_name)
  lat_info <- extract_coord_or_dim(nc, "y", preferred_name = preferred_lat_name)
  time_info <- extract_coord_or_dim(nc, "t", preferred_name = preferred_time_name)
  depth_info <- extract_coord_or_dim(nc, "z", preferred_name = preferred_depth_name)

  if (is.null(lon_info) || is.null(lat_info)) {
    stop_with_log("Unable to identify horizontal coordinates in the NetCDF file.", log_file)
  }

  if (!is.na(inventory_grid_type) && inventory_grid_type == "regular" &&
      (!is.null(dim(lon_info$values)) || !is.null(dim(lat_info$values)))) {
    log_message(
      "WARN",
      "Inventory marks this product as regular, but the LP NetCDF exposes 2D lon/lat coordinates; treating it as curvilinear.",
      log_file
    )
    inventory_grid_type <- NA_character_
  }

  if (!is.na(inventory_grid_type) && inventory_grid_type == "regular") {
    if (!is.null(preferred_lon_name) && preferred_lon_name %in% names(nc$dim)) {
      lon_info <- list(kind = "dim", name = preferred_lon_name, values = nc$dim[[preferred_lon_name]]$vals)
    }
    if (!is.null(preferred_lat_name) && preferred_lat_name %in% names(nc$dim)) {
      lat_info <- list(kind = "dim", name = preferred_lat_name, values = nc$dim[[preferred_lat_name]]$vals)
    }
    if (!is.null(preferred_depth_name) && preferred_depth_name %in% names(nc$dim)) {
      depth_info <- list(kind = "dim", name = preferred_depth_name, values = nc$dim[[preferred_depth_name]]$vals)
    }
    if (!is.null(preferred_time_name) && preferred_time_name %in% names(nc$dim)) {
      time_info <- list(kind = "dim", name = preferred_time_name, values = nc$dim[[preferred_time_name]]$vals)
    }
  }

  grid_type <- if (!is.na(inventory_grid_type)) {
    inventory_grid_type
  } else if (is.null(dim(lon_info$values)) && is.null(dim(lat_info$values))) {
    "regular"
  } else {
    "curvilinear"
  }

  crs_value <- detect_crs(nc, lon_info, lat_info, preferred_crs = preferred_crs)

  var_dt <- data.table::rbindlist(
    lapply(
      names(nc$var),
      function(vn) {
        v <- nc$var[[vn]]
        dim_names <- vapply(v$dim, function(d) d$name, FUN.VALUE = character(1))
        dim_lengths <- vapply(v$dim, function(d) as.integer(d$len), FUN.VALUE = integer(1))

        roles <- vapply(
          dim_names,
          infer_dim_role_from_name,
          FUN.VALUE = character(1),
          time_name = if (is.null(time_info)) NULL else time_info$name,
          depth_name = if (is.null(depth_info)) NULL else depth_info$name,
          lon_name = if (is.null(lon_info)) NULL else lon_info$name,
          lat_name = if (is.null(lat_info)) NULL else lat_info$name
        )

        data.table::data.table(
          var_name = vn,
          ndims = length(v$dim),
          dim_names = paste(dim_names, collapse = "|"),
          dim_lengths = paste(dim_lengths, collapse = "|"),
          x_dim_pos = match("x", roles),
          y_dim_pos = match("y", roles),
          z_dim_pos = match("z", roles),
          t_dim_pos = match("t", roles),
          has_depth = "z" %in% roles,
          has_time = "t" %in% roles,
          long_name = scalar_chr(safe_attr(nc, vn, "long_name")),
          standard_name = scalar_chr(safe_attr(nc, vn, "standard_name")),
          units = scalar_chr(safe_attr(nc, vn, "units")),
          fill_value = first_non_blank(
            safe_attr(nc, vn, "_FillValue"),
            safe_attr(nc, vn, "missing_value")
          )
        )
      }
    ),
    fill = TRUE,
    use.names = TRUE
  )

  if (nrow(inventory_var_dt) == 0L) {
    inventory_dt <- read_inventory_metadata(inventory_path, log_file = log_file)
    inventory_dt <- filter_inventory_to_nc_file(inventory_dt, nc_path = nc_path, log_file = log_file)
    var_dt <- harmonize_inventory(var_dt, inventory_dt)
  } else {
    var_dt <- harmonize_inventory(var_dt, inventory_var_dt)
  }

  time_units <- units_from_coord_info(nc, time_info)
  if ((is.na(time_units) || !nzchar(time_units)) && !is.null(inventory_summary)) {
    time_units <- scalar_chr(inventory_summary$time_units)
  }
  time_raw_values <- if (is.null(time_info)) numeric(0) else (time_info$values %||% numeric(0))
  time_parsed <- parse_time_values(time_raw_values, time_units)
  time_values <- time_parsed$times
  time_step_hours <- if (length(time_values) > 1L) {
    stats::median(as.numeric(diff(time_values), units = "hours"), na.rm = TRUE)
  } else {
    1
  }

  if (grid_type == "regular") {
    grid_extent <- compute_extent_from_coords(lon_info$values, lat_info$values)
    grid <- list(
      type = "regular",
      crs = crs_value,
      x = as.numeric(lon_info$values),
      y = as.numeric(lat_info$values),
      extent = grid_extent,
      nx = length(lon_info$values),
      ny = length(lat_info$values),
      lon_name = lon_info$name,
      lat_name = lat_info$name
    )
  } else {

    read_coord_xy <- function(coord_info, role_label) {
      vals <- coord_info$values

      if (identical(coord_info$kind, "var")) {
        coord_arr <- read_nc_name_values(nc, coord_info$name, fallback = vals)
        coord_dim_names <- coord_info$dim_names

        if ((is.null(coord_dim_names) || length(coord_dim_names) == 0L) &&
            coord_info$name %in% names(nc$var)) {
          coord_var <- nc$var[[coord_info$name]]
          coord_dim_names <- vapply(coord_var$dim, function(d) d$name, FUN.VALUE = character(1))
        }

        if ((is.null(coord_dim_names) || length(coord_dim_names) == 0L) &&
            !is.null(dim(coord_arr)) && length(dim(coord_arr)) >= 2L) {
          coord_dim_names <- infer_coord_dim_names_from_values(nc, coord_arr)
        }

        if (!is.null(dim(coord_arr)) && length(dim(coord_arr)) >= 2L &&
            !is.null(coord_dim_names) && length(coord_dim_names) >= length(dim(coord_arr))) {
          return(align_array_to_xy(
            coord_arr,
            dim_names = coord_dim_names[seq_along(dim(coord_arr))],
            time_name = if (is.null(time_info)) NULL else time_info$name,
            depth_name = if (is.null(depth_info)) NULL else depth_info$name,
            lon_name = lon_info$name,
            lat_name = lat_info$name
          ))
        }

        return(coord_arr)
      }

      if (!is.null(dim(vals)) && length(dim(vals)) >= 2L) {
        coord_dim_names <- coord_info$dim_names
        if (is.null(coord_dim_names) || length(coord_dim_names) == 0L) {
          coord_dim_names <- infer_coord_dim_names_from_values(nc, vals)
        }

        if (!is.null(coord_dim_names) && length(coord_dim_names) >= length(dim(vals))) {
          return(align_array_to_xy(
            vals,
            dim_names = coord_dim_names[seq_along(dim(vals))],
            time_name = if (is.null(time_info)) NULL else time_info$name,
            depth_name = if (is.null(depth_info)) NULL else depth_info$name,
            lon_name = lon_info$name,
            lat_name = lat_info$name
          ))
        }
      }

      vals
    }

    expand_coord_vector <- function(vals, nx, ny, role_label) {
      vals <- as.numeric(vals)
      if (length(vals) == nx && role_label == "x") {
        return(matrix(rep(vals, times = ny), nrow = nx, ncol = ny))
      }
      if (length(vals) == ny && role_label == "x") {
        return(matrix(rep(vals, each = nx), nrow = nx, ncol = ny))
      }
      if (length(vals) == ny && role_label == "y") {
        return(matrix(rep(vals, each = nx), nrow = nx, ncol = ny))
      }
      if (length(vals) == nx && role_label == "y") {
        return(matrix(rep(vals, times = ny), nrow = nx, ncol = ny))
      }
      stop_with_log(
        paste0(
          "Could not expand 1D ", role_label,
          " coordinate to the detected curvilinear grid size ", nx, " x ", ny, "."
        ),
        log_file
      )
    }

    lon_xy <- read_coord_xy(lon_info, "x")
    lat_xy <- read_coord_xy(lat_info, "y")

    if (is.null(dim(lon_xy)) && is.null(dim(lat_xy))) {
      nx0 <- length(lon_xy)
      ny0 <- length(lat_xy)
      lon_xy <- expand_coord_vector(lon_xy, nx0, ny0, "x")
      lat_xy <- expand_coord_vector(lat_xy, nx0, ny0, "y")
    } else if (is.null(dim(lon_xy)) && !is.null(dim(lat_xy))) {
      lon_xy <- expand_coord_vector(lon_xy, dim(lat_xy)[1], dim(lat_xy)[2], "x")
    } else if (!is.null(dim(lon_xy)) && is.null(dim(lat_xy))) {
      lat_xy <- expand_coord_vector(lat_xy, dim(lon_xy)[1], dim(lon_xy)[2], "y")
    }

    if (is.null(dim(lon_xy)) || is.null(dim(lat_xy)) ||
        length(dim(lon_xy)) != 2L || length(dim(lat_xy)) != 2L) {
      stop_with_log("Curvilinear coordinate arrays could not be resolved to 2D lon/lat matrices.", log_file)
    }

    template <- compute_curvilinear_template(lon_xy, lat_xy)
    grid <- list(
      type = "curvilinear",
      crs = crs_value,
      lon_xy = lon_xy,
      lat_xy = lat_xy,
      template = template,
      nx = dim(lon_xy)[1],
      ny = dim(lon_xy)[2],
      lon_name = lon_info$name,
      lat_name = lat_info$name
    )
  }

  structure <- list(
    created_utc = format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"),
    nc_path_used_for_setup = normalizePath(nc_path, winslash = "/", mustWork = FALSE),
    inventory_path = normalizePath(inventory_path, winslash = "/", mustWork = FALSE),
    global_title = scalar_chr(safe_attr(nc, 0, "title")),
    grid = grid,
    time = list(
      name = if (is.null(time_info)) NA_character_ else time_info$name,
      values_raw = time_raw_values,
      values = time_values,
      units = time_units,
      origin = time_parsed$origin,
      unit = time_parsed$unit,
      step_hours = time_step_hours,
      horizon_hours = if (length(time_values) > 0L) {
        as.numeric(max(time_values) - min(time_values), units = "hours")
      } else {
        NA_real_
      }
    ),
    depth = list(
      name = if (is.null(depth_info)) NA_character_ else depth_info$name,
      values = if (is.null(depth_info)) numeric(0) else as.numeric(depth_info$values),
      units = {
        depth_units_value <- units_from_coord_info(nc, depth_info)
        if ((is.na(depth_units_value) || !nzchar(depth_units_value)) && !is.null(inventory_summary)) {
          depth_units_value <- scalar_chr(inventory_summary$depth_units)
        }
        if (is.na(depth_units_value) || !nzchar(depth_units_value)) "m" else depth_units_value
      }
    ),
    variables = var_dt
  )

  if (!is.null(inventory_summary)) {
    inv_depth_length <- scalar_int(inventory_summary$depth_length)
    if (!is.na(inv_depth_length) && length(structure$depth$values) != inv_depth_length) {
      log_message(
        "WARN",
        paste0(
          "Detected depth axis length ",
          length(structure$depth$values),
          " differs from inventory depth_length=",
          inv_depth_length,
          "."
        ),
        log_file
      )
    }

    inv_lon_length <- scalar_int(inventory_summary$lon_length)
    inv_lat_length <- scalar_int(inventory_summary$lat_length)
    if (!is.na(inv_lon_length) && !is.na(inv_lat_length) &&
        (structure$grid$nx != inv_lon_length || structure$grid$ny != inv_lat_length)) {
      log_message(
        "WARN",
        paste0(
          "Detected grid size ",
          structure$grid$nx, " x ", structure$grid$ny,
          " differs from inventory grid size ",
          inv_lon_length, " x ", inv_lat_length,
          "."
        ),
        log_file
      )
    }
  }

  log_message(
    "INFO",
    paste0(
      "Characterization complete. Grid type: ", structure$grid$type,
      "; forecast horizon (hours): ", format(structure$time$horizon_hours, digits = 6),
      "; depth levels: ", length(structure$depth$values),
      "; variables found: ", nrow(structure$variables), "."
    ),
    log_file
  )

  structure
}

quick_validate_against_cache <- function(nc_path, cache, log_file = NULL, stop_on_error = TRUE, nc = NULL) {
  log_message("INFO", "Validating current NetCDF against cached structure assumptions.", log_file)

  validation_failed <- function(reason) {
    if (isTRUE(stop_on_error)) {
      stop_with_log(paste0(reason, " Re-run with --refresh_structure=TRUE."), log_file)
    }

    log_message(
      "WARN",
      paste0(
        "Skipping NetCDF file because it is not compatible with the cached structure: ",
        normalizePath(nc_path, winslash = "/", mustWork = FALSE),
        ". Reason: ", reason
      ),
      log_file
    )
    FALSE
  }

  close_nc <- FALSE
  if (is.null(nc)) {
    nc <- tryCatch(
      ncdf4::nc_open(nc_path, readunlim = FALSE),
      error = function(e) e
    )
    if (inherits(nc, "error")) {
      return(validation_failed(paste0("Unable to open NetCDF file. ", conditionMessage(nc))))
    }
    close_nc <- TRUE
  }
  if (isTRUE(close_nc)) {
    on.exit(ncdf4::nc_close(nc), add = TRUE)
  }

  nc_vars <- names(nc$var)
  required_vars <- unique(cache$variables$var_name)
  missing_vars <- setdiff(required_vars, nc_vars)

  if (length(missing_vars) > 0L) {
    return(validation_failed(
      paste0(
        "Cached structure is not compatible with this NetCDF file. Missing variables: ",
        paste(missing_vars, collapse = ", "),
        "."
      )
    ))
  }

  for (i in seq_len(nrow(cache$variables))) {
    row <- cache$variables[i]
    current_var <- nc$var[[row$var_name]]
    current_dim_names <- vapply(current_var$dim, function(d) d$name, FUN.VALUE = character(1))
    cached_dim_names <- unlist(strsplit(row$dim_names, "\\|", fixed = FALSE))
    if (!identical(tolower(current_dim_names), tolower(cached_dim_names))) {
      return(validation_failed(
        paste0(
          "Cached structure does not match current dimension order for variable '",
          row$var_name, "'."
        )
      ))
    }

    current_dim_lengths <- vapply(current_var$dim, function(d) as.integer(d$len), FUN.VALUE = integer(1))
    cached_dim_lengths <- as.integer(unlist(strsplit(row$dim_lengths, "\\|", fixed = FALSE)))
    if (!identical(current_dim_lengths, cached_dim_lengths)) {
      return(validation_failed(
        paste0(
          "Cached structure does not match current dimension lengths for variable '",
          row$var_name, "'."
        )
      ))
    }
  }

  log_message("INFO", "Cached structure validation succeeded.", log_file)
  TRUE
}

write_structure_summary <- function(structure, output_dir) {
  safe_dir_create(output_dir)

  var_path <- file.path(output_dir, "product_structure_variables.csv")
  data.table::fwrite(structure$variables, var_path, na = "")

  summary_txt <- c(
    paste0("Created UTC: ", structure$created_utc),
    paste0("Setup file: ", structure$nc_path_used_for_setup),
    paste0("Inventory path: ", structure$inventory_path),
    paste0("Grid type: ", structure$grid$type),
    paste0("CRS: ", structure$grid$crs %||% ""),
    paste0("Grid size: ", structure$grid$nx, " x ", structure$grid$ny),
    paste0("Time variable: ", structure$time$name),
    paste0("Time units: ", structure$time$units),
    paste0("Forecast horizon (hours): ", structure$time$horizon_hours),
    paste0("Depth variable: ", structure$depth$name),
    paste0("Depth levels: ", length(structure$depth$values))
  )
  writeLines(summary_txt, con = file.path(output_dir, "product_structure_summary.txt"))
  invisible(var_path)
}

cache_needs_rebuild <- function(cache, inventory_path, nc_path, log_file = NULL) {
  inventory_summary <- read_inventory_file_summary(inventory_path, nc_path = nc_path, log_file = NULL)
  if (is.null(inventory_summary)) {
    return(FALSE)
  }

  inv_depth_length <- scalar_int(inventory_summary$depth_length)
  inv_lon_length <- scalar_int(inventory_summary$lon_length)
  inv_lat_length <- scalar_int(inventory_summary$lat_length)
  inv_grid_type <- tolower(scalar_chr(inventory_summary$spatial_grid_type))

  if (!is.na(inv_depth_length) && length(cache$depth$values) != inv_depth_length) {
    log_message(
      "WARN",
      paste0(
        "Existing cache depth length ",
        length(cache$depth$values),
        " does not match inventory depth_length=",
        inv_depth_length,
        ". Rebuilding cache."
      ),
      log_file
    )
    return(TRUE)
  }

  if (!is.na(inv_lon_length) && !is.na(inv_lat_length) &&
      (cache$grid$nx != inv_lon_length || cache$grid$ny != inv_lat_length)) {
    log_message(
      "WARN",
      paste0(
        "Existing cache grid size ",
        cache$grid$nx, " x ", cache$grid$ny,
        " does not match inventory grid size ",
        inv_lon_length, " x ", inv_lat_length,
        ". Rebuilding cache."
      ),
      log_file
    )
    return(TRUE)
  }

  if (!is.na(inv_grid_type) && grepl("regular", inv_grid_type) && !identical(cache$grid$type, "regular")) {
    log_message(
      "WARN",
      paste0(
        "Existing cache grid type '",
        cache$grid$type,
        "' does not match inventory grid type 'regular'. Rebuilding cache."
      ),
      log_file
    )
    return(TRUE)
  }

  FALSE
}


refresh_structure_time_from_nc <- function(structure, nc_path, log_file = NULL, nc = NULL) {
  time_name <- scalar_chr(structure$time$name)
  if (is_blank_scalar(time_name)) {
    return(structure)
  }

  close_nc <- FALSE
  if (is.null(nc)) {
    nc <- ncdf4::nc_open(nc_path, readunlim = FALSE)
    close_nc <- TRUE
  }
  if (isTRUE(close_nc)) {
    on.exit(ncdf4::nc_close(nc), add = TRUE)
  }

  time_raw_values <- structure$time$values_raw %||% numeric(0)
  time_units <- scalar_chr(structure$time$units)

  if (time_name %in% names(nc$dim)) {
    time_raw_values <- nc$dim[[time_name]]$vals %||% time_raw_values
    dim_units <- scalar_chr(nc$dim[[time_name]]$units)
    if (!is.na(dim_units) && nzchar(dim_units)) {
      time_units <- dim_units
    }
  }

  if (time_name %in% names(nc$var)) {
    var_values <- tryCatch(
      ncdf4::ncvar_get(nc, time_name),
      error = function(e) NULL
    )
    if (!is.null(var_values)) {
      time_raw_values <- var_values
    }
    var_units <- scalar_chr(safe_attr(nc, time_name, "units"))
    if (!is.na(var_units) && nzchar(var_units)) {
      time_units <- var_units
    }
  }

  time_parsed <- parse_time_values(time_raw_values, time_units)
  time_values <- time_parsed$times
  time_step_hours <- if (length(time_values) > 1L) {
    stats::median(as.numeric(diff(time_values), units = "hours"), na.rm = TRUE)
  } else {
    structure$time$step_hours %||% 1
  }

  structure$time$values_raw <- time_raw_values
  structure$time$values <- time_values
  structure$time$units <- time_units
  structure$time$origin <- time_parsed$origin
  structure$time$unit <- time_parsed$unit
  structure$time$step_hours <- time_step_hours
  structure$time$horizon_hours <- if (length(time_values) > 0L) {
    as.numeric(max(time_values) - min(time_values), units = "hours")
  } else {
    NA_real_
  }

  structure
}

# ----------------------- Depth / window definitions -----------------------

nearest_surface_index <- function(depth_vals) {
  if (length(depth_vals) == 0L) return(NA_integer_)
  if (all(depth_vals >= 0, na.rm = TRUE)) {
    return(which.min(depth_vals))
  }
  which.min(abs(depth_vals))
}

nearest_bottom_index <- function(depth_vals) {
  if (length(depth_vals) == 0L) return(NA_integer_)
  if (all(depth_vals >= 0, na.rm = TRUE)) {
    return(which.max(depth_vals))
  }
  which.max(abs(depth_vals))
}

depth_label_from_index <- function(depth_vals, idx, depth_units = "m") {
  if (is.na(idx) || length(depth_vals) == 0L) {
    return("nodepth")
  }
  val <- depth_vals[idx]
  paste0("z_", sprintf("%.3f", val), depth_units)
}

normalize_summary_fun <- function(x) {
  x <- tolower(trimws(x))
  if (x %in% c("mean", "avg", "average")) return("mean")
  if (x %in% c("min", "minimum")) return("min")
  if (x %in% c("max", "maximum")) return("max")
  x
}

resolve_depth_targets <- function(structure, product_row, log_file = NULL) {
  depth_vals <- structure$depth$values
  has_depth <- FALSE

  if (product_row$source_type == "vector_magnitude") {
    src1 <- structure$variables[var_name == product_row$source_var_1]
    src2 <- structure$variables[var_name == product_row$source_var_2]
    has_depth <- (nrow(src1) == 1L && isTRUE(src1$has_depth)) || (nrow(src2) == 1L && isTRUE(src2$has_depth))
  } else {
    src1 <- structure$variables[var_name == product_row$source_var_1]
    has_depth <- nrow(src1) == 1L && isTRUE(src1$has_depth)
  }

  if (!has_depth || length(depth_vals) == 0L) {
    return(list(list(kind = "single", label = "nodepth", depth_index = NA_integer_)))
  }

  depth_mode <- product_row$depth_mode

  if (depth_mode == "surface") {
    idx <- nearest_surface_index(depth_vals)
    return(list(list(kind = "single", label = depth_label_from_index(depth_vals, idx, structure$depth$units), depth_index = idx)))
  }

  if (depth_mode == "bottom") {
    idx <- nearest_bottom_index(depth_vals)
    return(list(list(kind = "single", label = depth_label_from_index(depth_vals, idx, structure$depth$units), depth_index = idx)))
  }

  if (depth_mode == "selected") {
    idx <- split_semicolon_int(product_row$depth_indices)
    idx <- idx[!is.na(idx)]
    if (length(idx) == 0L) {
      stop_with_log(
        paste0("alert_id='", product_row$alert_id, "' uses depth_mode='selected' but depth_indices is empty."),
        log_file
      )
    }
    idx <- unique(idx)
    bad <- idx[idx < 1L | idx > length(depth_vals)]
    if (length(bad) > 0L) {
      stop_with_log(
        paste0(
          "alert_id='", product_row$alert_id, "' requested out-of-range depth_indices: ",
          paste(bad, collapse = ", "), "."
        ),
        log_file
      )
    }
    return(lapply(idx, function(i) {
      list(kind = "single", label = depth_label_from_index(depth_vals, i, structure$depth$units), depth_index = i)
    }))
  }

  if (depth_mode == "summary") {
    summary_fun <- normalize_summary_fun(product_row$summary_fun)
    if (!summary_fun %in% c("mean", "min", "max")) {
      stop_with_log(
        paste0(
          "alert_id='", product_row$alert_id,
          "' uses depth_mode='summary' but summary_fun must be one of mean/min/max."
        ),
        log_file
      )
    }
    idx <- split_semicolon_int(product_row$depth_indices)
    if (length(idx) == 0L) idx <- seq_along(depth_vals)
    idx <- idx[!is.na(idx)]
    bad <- idx[idx < 1L | idx > length(depth_vals)]
    if (length(bad) > 0L) {
      stop_with_log(
        paste0(
          "alert_id='", product_row$alert_id, "' requested out-of-range depth_indices: ",
          paste(bad, collapse = ", "), "."
        ),
        log_file
      )
    }
    return(list(list(kind = "summary", label = paste0("summary_", summary_fun), depth_indices = idx, summary_fun = summary_fun)))
  }

  if (depth_mode == "max_severity") {
    idx <- split_semicolon_int(product_row$depth_indices)
    if (length(idx) == 0L) idx <- seq_along(depth_vals)
    idx <- idx[!is.na(idx)]
    bad <- idx[idx < 1L | idx > length(depth_vals)]
    if (length(bad) > 0L) {
      stop_with_log(
        paste0(
          "alert_id='", product_row$alert_id, "' requested out-of-range depth_indices: ",
          paste(bad, collapse = ", "), "."
        ),
        log_file
      )
    }
    return(list(list(kind = "max_severity", label = "max_severity", depth_indices = idx)))
  }

  stop_with_log(
    paste0("Unsupported depth_mode for alert_id='", product_row$alert_id, "': ", depth_mode),
    log_file
  )
}

build_forecast_windows <- function(time_values, window_hours = 6L) {
  if (length(time_values) == 0L) {
    return(data.table::data.table())
  }

  time_values <- as.POSIXct(time_values, tz = "UTC")
  start_time <- min(time_values)
  offset_hours <- as.numeric(difftime(time_values, start_time, units = "hours"))
  block_id <- floor(offset_hours / window_hours)

  dt <- data.table::data.table(
    time_index = seq_along(time_values),
    forecast_time = time_values,
    offset_hours = offset_hours,
    block_id = block_id
  )

  wins <- dt[, .(
    time_indices = list(time_index),
    start_time = min(forecast_time),
    end_time = max(forecast_time),
    start_hour = min(offset_hours),
    end_hour = max(offset_hours)
  ), by = block_id][order(block_id)]

  wins[, window_label := sprintf("t%03d_to_t%03d", as.integer(floor(start_hour)), as.integer(ceiling(start_hour + window_hours)))]
  wins
}

# ------------------------- Slice reading helpers -------------------------

get_var_meta <- function(structure, var_name, log_file = NULL) {
  target_name <- trimws(as.character(var_name))
  row <- structure$variables[tolower(var_name) == tolower(target_name)]
  if (nrow(row) == 0L) {
    stop_with_log(paste0("Variable not found in characterized structure: ", target_name), log_file)
  }
  row[1]
}

dim_lengths_from_meta <- function(var_meta) {
  as.integer(unlist(strsplit(var_meta$dim_lengths, "\\|", fixed = FALSE)))
}

build_start_count <- function(var_meta, time_index = NULL, depth_index = NULL) {
  lens <- dim_lengths_from_meta(var_meta)
  start <- rep(1L, length(lens))
  count <- lens

  x_pos <- scalar_int(var_meta$x_dim_pos)
  y_pos <- scalar_int(var_meta$y_dim_pos)
  z_pos <- scalar_int(var_meta$z_dim_pos)
  t_pos <- scalar_int(var_meta$t_dim_pos)

  if (!is.na(t_pos) && !is.null(time_index) && length(time_index) == 1L) {
    start[t_pos] <- as.integer(time_index)
    count[t_pos] <- 1L
  }

  if (!is.na(z_pos) && !is.null(depth_index) && length(depth_index) == 1L) {
    start[z_pos] <- as.integer(depth_index)
    count[z_pos] <- 1L
  }

  frozen_dims <- setdiff(seq_along(lens), c(x_pos, y_pos, z_pos, t_pos))
  if (length(frozen_dims) > 0L) {
    start[frozen_dims] <- 1L
    count[frozen_dims] <- 1L
  }

  list(start = start, count = count)
}

read_native_slice_xy <- function(nc, var_meta, time_index = NULL, depth_index = NULL, log_file = NULL) {
  sc <- build_start_count(var_meta, time_index = time_index, depth_index = depth_index)

  arr <- ncdf4::ncvar_get(
    nc,
    varid = var_meta$var_name,
    start = sc$start,
    count = sc$count,
    collapse_degen = FALSE
  )

  arr_dims <- dim(arr)
  if (is.null(arr_dims) || length(arr_dims) < 2L) {
    stop_with_log(
      paste0(
        "Expected at least a 2D slice for variable '", var_meta$var_name,
        "', but the returned object has dimensions: ",
        if (is.null(arr_dims)) "<none>" else paste(arr_dims, collapse = " x "),
        "."
      ),
      log_file
    )
  }

  x_pos <- scalar_int(var_meta$x_dim_pos)
  y_pos <- scalar_int(var_meta$y_dim_pos)

  if (is.na(x_pos) || is.na(y_pos)) {
    stop_with_log(
      paste0("Unable to identify x/y dimensions for variable slice: ", var_meta$var_name),
      log_file
    )
  }

  if (x_pos > length(arr_dims) || y_pos > length(arr_dims)) {
    stop_with_log(
      paste0(
        "x/y dimension positions for variable '", var_meta$var_name,
        "' exceed the returned array rank. Returned dims: ",
        paste(arr_dims, collapse = " x "),
        "; x_pos=", x_pos,
        "; y_pos=", y_pos,
        "."
      ),
      log_file
    )
  }

  perm <- c(x_pos, y_pos, setdiff(seq_along(arr_dims), c(x_pos, y_pos)))
  arr <- aperm(arr, perm = perm)
  arr_dims <- dim(arr)

  if (length(arr_dims) == 2L) {
    mat <- arr
  } else {
    extra_dims <- arr_dims[-c(1, 2)]
    if (any(extra_dims > 1L, na.rm = TRUE)) {
      stop_with_log(
        paste0(
          "Variable slice for '", var_meta$var_name,
          "' still has non-singleton dimensions beyond x/y after selecting time/depth. ",
          "Returned dims after reordering: ", paste(arr_dims, collapse = " x "),
          "; start=", paste(sc$start, collapse = ","),
          "; count=", paste(sc$count, collapse = ","),
          "."
        ),
        log_file
      )
    }
    mat <- matrix(as.vector(arr), nrow = arr_dims[1], ncol = arr_dims[2])
  }

  storage.mode(mat) <- "double"
  mat
}

read_product_slice_xy <- function(nc, structure, product_row, time_index = NULL, depth_index = NULL, log_file = NULL) {
  if (product_row$source_type == "vector_magnitude") {
    u_meta <- get_var_meta(structure, product_row$source_var_1, log_file = log_file)
    v_meta <- get_var_meta(structure, product_row$source_var_2, log_file = log_file)
    u <- read_native_slice_xy(nc, u_meta, time_index = time_index, depth_index = depth_index, log_file = log_file)
    v <- read_native_slice_xy(nc, v_meta, time_index = time_index, depth_index = depth_index, log_file = log_file)
    out <- sqrt(u^2 + v^2)
    return(out)
  }

  src_meta <- get_var_meta(structure, product_row$source_var_1, log_file = log_file)
  read_native_slice_xy(nc, src_meta, time_index = time_index, depth_index = depth_index, log_file = log_file)
}

compose_summary_slice_xy <- function(nc, structure, product_row, time_index, depth_indices, summary_fun, log_file = NULL) {
  summary_fun <- normalize_summary_fun(summary_fun)

  if (summary_fun == "mean") {
    sum_mat <- NULL
    count_mat <- NULL
    for (z in depth_indices) {
      slice <- read_product_slice_xy(nc, structure, product_row, time_index = time_index, depth_index = z, log_file = log_file)
      if (is.null(sum_mat)) {
        sum_mat <- matrix(NA_real_, nrow = nrow(slice), ncol = ncol(slice))
        count_mat <- matrix(NA_real_, nrow = nrow(slice), ncol = ncol(slice))
      }
      tmp <- nonmissing_mean_pair(sum_mat, count_mat, slice)
      sum_mat <- tmp$sum
      count_mat <- tmp$count
    }
    out <- sum_mat / count_mat
    out[count_mat == 0 | is.na(count_mat)] <- NA_real_
    return(out)
  }

  if (summary_fun == "min") {
    agg <- NULL
    for (z in depth_indices) {
      slice <- read_product_slice_xy(nc, structure, product_row, time_index = time_index, depth_index = z, log_file = log_file)
      agg <- safe_pmin(agg, slice)
    }
    return(agg)
  }

  if (summary_fun == "max") {
    agg <- NULL
    for (z in depth_indices) {
      slice <- read_product_slice_xy(nc, structure, product_row, time_index = time_index, depth_index = z, log_file = log_file)
      agg <- safe_pmax(agg, slice)
    }
    return(agg)
  }

  stop_with_log(paste0("Unsupported summary_fun: ", summary_fun), log_file)
}

read_depth_target_slice_xy <- function(nc, structure, product_row, time_index, depth_target, log_file = NULL) {
  if (depth_target$kind == "single") {
    return(read_product_slice_xy(
      nc, structure, product_row,
      time_index = time_index,
      depth_index = depth_target$depth_index,
      log_file = log_file
    ))
  }

  if (depth_target$kind == "summary") {
    return(compose_summary_slice_xy(
      nc, structure, product_row,
      time_index = time_index,
      depth_indices = depth_target$depth_indices,
      summary_fun = depth_target$summary_fun,
      log_file = log_file
    ))
  }

  stop_with_log(paste0("Depth target kind not supported for direct slice reads: ", depth_target$kind), log_file)
}

build_block_start_count <- function(var_meta, time_indices = NULL, depth_indices = NULL) {
  lens <- dim_lengths_from_meta(var_meta)
  start <- rep(1L, length(lens))
  count <- lens

  x_pos <- scalar_int(var_meta$x_dim_pos)
  y_pos <- scalar_int(var_meta$y_dim_pos)
  z_pos <- scalar_int(var_meta$z_dim_pos)
  t_pos <- scalar_int(var_meta$t_dim_pos)

  time_indices <- sort(unique(as.integer(time_indices[!is.na(time_indices)])))
  depth_indices <- sort(unique(as.integer(depth_indices[!is.na(depth_indices)])))

  if (!is.na(t_pos)) {
    if (length(time_indices) > 0L) {
      start[t_pos] <- min(time_indices)
      count[t_pos] <- max(time_indices) - min(time_indices) + 1L
    } else {
      start[t_pos] <- 1L
      count[t_pos] <- 1L
    }
  }

  if (!is.na(z_pos)) {
    if (length(depth_indices) > 0L) {
      start[z_pos] <- min(depth_indices)
      count[z_pos] <- max(depth_indices) - min(depth_indices) + 1L
    } else {
      start[z_pos] <- 1L
      count[z_pos] <- 1L
    }
  }

  frozen_dims <- setdiff(seq_along(lens), c(x_pos, y_pos, z_pos, t_pos))
  if (length(frozen_dims) > 0L) {
    start[frozen_dims] <- 1L
    count[frozen_dims] <- 1L
  }

  list(
    start = start,
    count = count,
    time_indices = time_indices,
    depth_indices = depth_indices,
    read_time = !is.na(t_pos) && length(time_indices) > 0L,
    read_depth = !is.na(z_pos) && length(depth_indices) > 0L
  )
}

read_native_block_xyzt <- function(nc, var_meta, time_indices = NULL, depth_indices = NULL, log_file = NULL) {
  bc <- build_block_start_count(var_meta, time_indices = time_indices, depth_indices = depth_indices)

  arr <- ncdf4::ncvar_get(
    nc,
    varid = var_meta$var_name,
    start = bc$start,
    count = bc$count,
    collapse_degen = FALSE
  )

  arr_dims <- dim(arr)
  if (is.null(arr_dims) || length(arr_dims) < 2L) {
    stop_with_log(
      paste0(
        "Expected at least a 2D block for variable '", var_meta$var_name,
        "', but the returned object has dimensions: ",
        if (is.null(arr_dims)) "<none>" else paste(arr_dims, collapse = " x "),
        "."
      ),
      log_file
    )
  }

  x_pos <- scalar_int(var_meta$x_dim_pos)
  y_pos <- scalar_int(var_meta$y_dim_pos)
  z_pos <- scalar_int(var_meta$z_dim_pos)
  t_pos <- scalar_int(var_meta$t_dim_pos)

  if (is.na(x_pos) || is.na(y_pos)) {
    stop_with_log(
      paste0("Unable to identify x/y dimensions for variable block: ", var_meta$var_name),
      log_file
    )
  }

  core_positions <- c(x_pos, y_pos)
  if (isTRUE(bc$read_depth)) core_positions <- c(core_positions, z_pos)
  if (isTRUE(bc$read_time)) core_positions <- c(core_positions, t_pos)

  if (any(core_positions > length(arr_dims))) {
    stop_with_log(
      paste0(
        "Requested x/y/z/t positions for variable '", var_meta$var_name,
        "' exceed the returned array rank. Returned dims: ",
        paste(arr_dims, collapse = " x "),
        "."
      ),
      log_file
    )
  }

  perm <- c(core_positions, setdiff(seq_along(arr_dims), core_positions))
  arr <- aperm(arr, perm = perm)
  arr_dims <- dim(arr)
  core_ndim <- length(core_positions)

  if (length(arr_dims) > core_ndim) {
    extra_dims <- arr_dims[(core_ndim + 1L):length(arr_dims)]
    if (any(extra_dims > 1L, na.rm = TRUE)) {
      stop_with_log(
        paste0(
          "Variable block for '", var_meta$var_name,
          "' still has non-singleton dimensions beyond the requested x/y/z/t block. ",
          "Returned dims after reordering: ", paste(arr_dims, collapse = " x "),
          "; start=", paste(bc$start, collapse = ","),
          "; count=", paste(bc$count, collapse = ","),
          "."
        ),
        log_file
      )
    }
    dim(arr) <- arr_dims[seq_len(core_ndim)]
  }

  if (isTRUE(bc$read_depth) || isTRUE(bc$read_time)) {
    subscripts <- lapply(dim(arr), seq_len)

    if (isTRUE(bc$read_depth)) {
      depth_rel <- bc$depth_indices - min(bc$depth_indices) + 1L
      subscripts[[3L]] <- depth_rel
    }

    if (isTRUE(bc$read_time)) {
      time_rel <- bc$time_indices - min(bc$time_indices) + 1L
      time_pos_out <- if (isTRUE(bc$read_depth)) 4L else 3L
      subscripts[[time_pos_out]] <- time_rel
    }

    arr <- do.call(`[`, c(list(arr), subscripts, list(drop = FALSE)))
  }

  storage.mode(arr) <- "double"
  arr
}

read_product_block_xyzt <- function(nc, structure, product_row, time_indices = NULL, depth_indices = NULL, log_file = NULL) {
  if (product_row$source_type == "vector_magnitude") {
    u_meta <- get_var_meta(structure, product_row$source_var_1, log_file = log_file)
    v_meta <- get_var_meta(structure, product_row$source_var_2, log_file = log_file)
    u_block <- read_native_block_xyzt(
      nc, u_meta,
      time_indices = time_indices,
      depth_indices = depth_indices,
      log_file = log_file
    )
    v_block <- read_native_block_xyzt(
      nc, v_meta,
      time_indices = time_indices,
      depth_indices = depth_indices,
      log_file = log_file
    )

    if (!identical(dim(u_block), dim(v_block))) {
      stop_with_log(
        paste0(
          "Vector magnitude source blocks do not share the same dimensions for alert_id='",
          product_row$alert_id,
          "'."
        ),
        log_file
      )
    }

    out <- sqrt(u_block^2 + v_block^2)
    storage.mode(out) <- "double"
    return(out)
  }

  src_meta <- get_var_meta(structure, product_row$source_var_1, log_file = log_file)
  read_native_block_xyzt(
    nc, src_meta,
    time_indices = time_indices,
    depth_indices = depth_indices,
    log_file = log_file
  )
}


native_block_cache_key <- function(var_meta, time_indices = NULL, depth_indices = NULL) {
  time_key <- if (is.null(time_indices) || length(time_indices) == 0L) {
    "notime"
  } else {
    paste(as.integer(time_indices), collapse = ",")
  }
  depth_key <- if (is.null(depth_indices) || length(depth_indices) == 0L) {
    "nodepth"
  } else {
    paste(as.integer(depth_indices), collapse = ",")
  }
  paste(scalar_chr(var_meta$var_name), time_key, depth_key, sep = "|")
}

read_native_block_xyzt_cached <- function(
  nc,
  var_meta,
  time_indices = NULL,
  depth_indices = NULL,
  block_cache = NULL,
  log_file = NULL
) {
  if (is.null(block_cache)) {
    return(read_native_block_xyzt(
      nc,
      var_meta,
      time_indices = time_indices,
      depth_indices = depth_indices,
      log_file = log_file
    ))
  }

  key <- native_block_cache_key(var_meta, time_indices = time_indices, depth_indices = depth_indices)
  if (exists(key, envir = block_cache, inherits = FALSE)) {
    return(get(key, envir = block_cache, inherits = FALSE))
  }

  block <- read_native_block_xyzt(
    nc,
    var_meta,
    time_indices = time_indices,
    depth_indices = depth_indices,
    log_file = log_file
  )
  assign(key, block, envir = block_cache)
  block
}

read_product_block_xyzt_cached <- function(
  nc,
  structure,
  product_row,
  time_indices = NULL,
  depth_indices = NULL,
  block_cache = NULL,
  log_file = NULL
) {
  if (product_row$source_type == "vector_magnitude") {
    u_meta <- get_var_meta(structure, product_row$source_var_1, log_file = log_file)
    v_meta <- get_var_meta(structure, product_row$source_var_2, log_file = log_file)
    u_block <- read_native_block_xyzt_cached(
      nc,
      u_meta,
      time_indices = time_indices,
      depth_indices = depth_indices,
      block_cache = block_cache,
      log_file = log_file
    )
    v_block <- read_native_block_xyzt_cached(
      nc,
      v_meta,
      time_indices = time_indices,
      depth_indices = depth_indices,
      block_cache = block_cache,
      log_file = log_file
    )

    if (!identical(dim(u_block), dim(v_block))) {
      stop_with_log(
        paste0(
          "Vector magnitude source blocks do not share the same dimensions for alert_id='",
          product_row$alert_id,
          "'."
        ),
        log_file
      )
    }

    out <- sqrt(u_block^2 + v_block^2)
    storage.mode(out) <- "double"
    return(out)
  }

  src_meta <- get_var_meta(structure, product_row$source_var_1, log_file = log_file)
  read_native_block_xyzt_cached(
    nc,
    src_meta,
    time_indices = time_indices,
    depth_indices = depth_indices,
    block_cache = block_cache,
    log_file = log_file
  )
}

read_depth_target_block_xyt_cached <- function(
  nc,
  structure,
  product_row,
  depth_target,
  time_indices,
  block_cache = NULL,
  log_file = NULL
) {
  if (depth_target$kind == "single") {
    depth_indices_for_read <- NULL
    if (!is.na(depth_target$depth_index)) {
      depth_indices_for_read <- depth_target$depth_index
    }

    block <- read_product_block_xyzt_cached(
      nc,
      structure,
      product_row,
      time_indices = time_indices,
      depth_indices = depth_indices_for_read,
      block_cache = block_cache,
      log_file = log_file
    )
    return(ensure_xyt_block(block, label = paste0(product_row$alert_id, "_block"), log_file = log_file))
  }

  if (depth_target$kind == "summary") {
    block <- read_product_block_xyzt_cached(
      nc,
      structure,
      product_row,
      time_indices = time_indices,
      depth_indices = depth_target$depth_indices,
      block_cache = block_cache,
      log_file = log_file
    )
    return(summarize_block_over_depth_to_xyt(block, depth_target$summary_fun, log_file = log_file))
  }

  stop_with_log(
    paste0("Depth target kind not supported for combined processing: ", depth_target$kind),
    log_file
  )
}

ensure_xyt_block <- function(block, label = "block", log_file = NULL) {
  block_dims <- dim(block)

  if (is.null(block_dims) || length(block_dims) < 2L) {
    stop_with_log(
      paste0(
        "Object '", label,
        "' could not be interpreted as an x/y/t block. Returned dimensions: ",
        if (is.null(block_dims)) "<none>" else paste(block_dims, collapse = " x "),
        "."
      ),
      log_file
    )
  }

  if (length(block_dims) == 2L) {
    out <- array(as.vector(block), dim = c(block_dims[1], block_dims[2], 1L))
  } else if (length(block_dims) == 3L) {
    out <- array(as.vector(block), dim = block_dims)
  } else if (length(block_dims) == 4L && block_dims[3] == 1L) {
    out <- array(as.vector(block), dim = c(block_dims[1], block_dims[2], block_dims[4]))
  } else {
    stop_with_log(
      paste0(
        "Object '", label,
        "' could not be reduced to an x/y/t block. Returned dimensions: ",
        paste(block_dims, collapse = " x "),
        "."
      ),
      log_file
    )
  }

  storage.mode(out) <- "double"
  out
}

summarize_block_over_depth_to_xyt <- function(block_xyzt, summary_fun, log_file = NULL) {
  summary_fun <- normalize_summary_fun(summary_fun)
  block_dims <- dim(block_xyzt)

  if (is.null(block_dims) || length(block_dims) < 3L) {
    return(ensure_xyt_block(block_xyzt, label = "summary_block", log_file = log_file))
  }

  if (length(block_dims) != 4L) {
    stop_with_log(
      paste0(
        "Summary depth block could not be interpreted as x/y/z/t. Returned dimensions: ",
        paste(block_dims, collapse = " x "),
        "."
      ),
      log_file
    )
  }

  nx <- block_dims[1]
  ny <- block_dims[2]
  nz <- block_dims[3]
  nt <- block_dims[4]

  if (nz == 1L) {
    out <- array(as.vector(block_xyzt), dim = c(nx, ny, nt))
    storage.mode(out) <- "double"
    return(out)
  }

  if (summary_fun == "mean") {
    sum_arr <- array(NA_real_, dim = c(nx, ny, nt))
    count_arr <- array(NA_real_, dim = c(nx, ny, nt))

    for (zi in seq_len(nz)) {
      slice_arr <- array(block_xyzt[, , zi, , drop = FALSE], dim = c(nx, ny, nt))
      keep <- !is.na(slice_arr)
      sum_arr[keep] <- ifelse(is.na(sum_arr[keep]), 0, sum_arr[keep]) + slice_arr[keep]
      count_arr[keep] <- ifelse(is.na(count_arr[keep]), 0, count_arr[keep]) + 1L
    }

    out <- sum_arr / count_arr
    out[count_arr == 0 | is.na(count_arr)] <- NA_real_
    return(out)
  }

  if (summary_fun == "min") {
    agg <- NULL
    for (zi in seq_len(nz)) {
      slice_arr <- array(block_xyzt[, , zi, , drop = FALSE], dim = c(nx, ny, nt))
      agg <- safe_pmin(agg, slice_arr)
    }
    return(agg)
  }

  if (summary_fun == "max") {
    agg <- NULL
    for (zi in seq_len(nz)) {
      slice_arr <- array(block_xyzt[, , zi, , drop = FALSE], dim = c(nx, ny, nt))
      agg <- safe_pmax(agg, slice_arr)
    }
    return(agg)
  }

  stop_with_log(paste0("Unsupported summary_fun: ", summary_fun), log_file)
}

compute_window_severity_from_xyt <- function(xyt, current_rel_idx, product_row, change_steps, log_file = NULL) {
  xyt <- ensure_xyt_block(xyt, label = "xyt", log_file = log_file)
  xyt_dims <- dim(xyt)
  nx <- xyt_dims[1]
  ny <- xyt_dims[2]
  nt <- xyt_dims[3]

  current_rel_idx <- as.integer(current_rel_idx)
  current_rel_idx <- current_rel_idx[!is.na(current_rel_idx) & current_rel_idx >= 1L & current_rel_idx <= nt]

  low_threshold <- scalar_num(product_row$low_threshold)
  high_threshold <- scalar_num(product_row$high_threshold)
  change_threshold <- scalar_num(product_row$change_threshold)
  change_direction <- scalar_chr(product_row$change_direction, default = "abs")

  threshold_max <- NULL
  change_max <- NULL
  max_abs_change <- NULL

  for (rel_t in current_rel_idx) {
    slice_xy <- xyt[, , rel_t, drop = TRUE]

    thr_xy <- threshold_severity_xy(
      slice_xy,
      low_threshold = low_threshold,
      high_threshold = high_threshold
    )
    threshold_max <- safe_pmax(threshold_max, thr_xy)

    if (change_steps > 0L) {
      prev_start <- max(1L, rel_t - change_steps)
      prev_end <- rel_t - 1L

      if (prev_end >= prev_start) {
        raw_prev_min <- NULL
        raw_prev_max <- NULL

        for (pj in prev_start:prev_end) {
          prev_xy <- xyt[, , pj, drop = TRUE]
          raw_prev_min <- safe_pmin(raw_prev_min, prev_xy)
          raw_prev_max <- safe_pmax(raw_prev_max, prev_xy)
        }

        raw_inc <- if (!is.null(raw_prev_min)) {
          absolute_change_xy(
            new_xy = slice_xy,
            old_xy = raw_prev_min
          )
        } else {
          NULL
        }
        raw_dec <- if (!is.null(raw_prev_max)) {
          absolute_change_xy(
            new_xy = slice_xy,
            old_xy = raw_prev_max
          )
        } else {
          NULL
        }

        max_abs_change <- safe_pmax(max_abs_change, safe_pmax(raw_inc, raw_dec))
      }
    }

    if (!is.na(change_threshold) && change_steps > 0L) {
      prev_start <- max(1L, rel_t - change_steps)
      prev_end <- rel_t - 1L

      if (prev_end >= prev_start) {
        prev_min <- NULL
        prev_max <- NULL

        if (change_direction %in% c("increase", "abs")) {
          for (pj in prev_start:prev_end) {
            prev_min <- safe_pmin(prev_min, xyt[, , pj, drop = TRUE])
          }
        }

        if (change_direction %in% c("decrease", "abs")) {
          for (pj in prev_start:prev_end) {
            prev_max <- safe_pmax(prev_max, xyt[, , pj, drop = TRUE])
          }
        }

        ch_xy <- NULL
        if (change_direction == "increase" && !is.null(prev_min)) {
          ch_xy <- change_severity_xy(
            new_xy = slice_xy,
            old_xy = prev_min,
            change_threshold = change_threshold,
            direction = "increase"
          )
        } else if (change_direction == "decrease" && !is.null(prev_max)) {
          ch_xy <- change_severity_xy(
            new_xy = slice_xy,
            old_xy = prev_max,
            change_threshold = change_threshold,
            direction = "decrease"
          )
        } else if (change_direction == "abs") {
          ch_inc <- if (!is.null(prev_min)) {
            change_severity_xy(
              new_xy = slice_xy,
              old_xy = prev_min,
              change_threshold = change_threshold,
              direction = "increase"
            )
          } else {
            NULL
          }
          ch_dec <- if (!is.null(prev_max)) {
            change_severity_xy(
              new_xy = slice_xy,
              old_xy = prev_max,
              change_threshold = change_threshold,
              direction = "decrease"
            )
          } else {
            NULL
          }
          ch_xy <- safe_pmax(ch_inc, ch_dec)
        }

        change_max <- safe_pmax(change_max, ch_xy)
      }
    }
  }

  if (is.null(threshold_max)) {
    threshold_max <- matrix(0, nrow = nx, ncol = ny)
  }
  if (is.null(change_max)) {
    change_max <- matrix(0, nrow = nx, ncol = ny)
  }
  if (is.null(max_abs_change)) {
    max_abs_change <- matrix(0, nrow = nx, ncol = ny)
  }

  list(
    threshold = threshold_max,
    change = change_max,
    combined = safe_pmax(threshold_max, change_max),
    max_abs_change = max_abs_change
  )
}

# -------------------------- Alert calculations --------------------------

scale_for_threshold <- function(thr_low, thr_high) {
  refs <- c(abs(thr_low), abs(thr_high))
  refs <- refs[is.finite(refs) & refs > 0]
  if (length(refs) == 0L) return(1)
  refs[1]
}

threshold_severity_xy <- function(slice_xy, low_threshold = NA_real_, high_threshold = NA_real_) {
  slice_xy <- coerce_xy_matrix(slice_xy, nx = nrow(slice_xy), ny = ncol(slice_xy), label = "slice_xy")

  if (all(is.na(c(low_threshold, high_threshold)))) {
    out <- matrix(0, nrow = nrow(slice_xy), ncol = ncol(slice_xy))
    out[is.na(slice_xy)] <- NA_real_
    return(out)
  }

  scale_ref <- scale_for_threshold(low_threshold, high_threshold)
  sev <- matrix(0, nrow = nrow(slice_xy), ncol = ncol(slice_xy))

  if (!is.na(low_threshold)) {
    sev_low <- (low_threshold - slice_xy) / scale_ref
    sev_low[sev_low < 0] <- 0
    sev <- safe_pmax(sev, sev_low)
  }

  if (!is.na(high_threshold)) {
    sev_high <- (slice_xy - high_threshold) / scale_ref
    sev_high[sev_high < 0] <- 0
    sev <- safe_pmax(sev, sev_high)
  }

  sev <- coerce_xy_matrix(sev, nx = nrow(slice_xy), ny = ncol(slice_xy), label = "threshold_severity")
  sev[is.na(slice_xy)] <- NA_real_
  sev
}

absolute_change_xy <- function(new_xy, old_xy) {
  new_xy <- coerce_xy_matrix(new_xy, nx = nrow(new_xy), ny = ncol(new_xy), label = "new_xy")
  old_xy <- coerce_xy_matrix(old_xy, nx = nrow(old_xy), ny = ncol(old_xy), label = "old_xy")

  out <- abs(new_xy - old_xy)
  out <- coerce_xy_matrix(out, nx = nrow(new_xy), ny = ncol(new_xy), label = "absolute_change")
  out[is.na(new_xy) | is.na(old_xy)] <- NA_real_
  out
}

change_severity_xy <- function(new_xy, old_xy, change_threshold, direction = "abs") {
  new_xy <- coerce_xy_matrix(new_xy, nx = nrow(new_xy), ny = ncol(new_xy), label = "new_xy")
  old_xy <- coerce_xy_matrix(old_xy, nx = nrow(old_xy), ny = ncol(old_xy), label = "old_xy")

  if (is.na(change_threshold) || change_threshold <= 0) {
    out <- matrix(0, nrow = nrow(new_xy), ncol = ncol(new_xy))
    out[is.na(new_xy) | is.na(old_xy)] <- NA_real_
    return(out)
  }

  direction <- tolower(direction)
  delta <- new_xy - old_xy

  mag <- switch(
    direction,
    abs = abs(delta),
    increase = delta,
    decrease = -delta,
    abs(delta)
  )
  mag[mag < 0] <- 0

  sev <- (mag - change_threshold) / change_threshold
  sev[sev < 0] <- 0
  sev <- coerce_xy_matrix(sev, nx = nrow(new_xy), ny = ncol(new_xy), label = "change_severity")
  sev[is.na(new_xy) | is.na(old_xy)] <- NA_real_
  sev
}

change_steps_from_days <- function(change_days, steps_per_day) {
  if (is.na(change_days) || change_days <= 0 || is.na(steps_per_day) || steps_per_day <= 0) {
    return(0L)
  }
  max(1L, as.integer(floor(change_days * steps_per_day + 1e-9)))
}

steps_per_day_from_structure <- function(structure) {
  n_time <- length(structure$time$values)
  if (!is.na(n_time) && n_time > 0L) {
    return(as.integer(n_time))
  }

  time_step_hours <- structure$time$step_hours %||% NA_real_
  if (!is.na(time_step_hours) && time_step_hours > 0) {
    return(max(1L, as.integer(round(24 / time_step_hours))))
  }

  24L
}

estimate_peak_window_memory_mb <- function(nx, ny, n_history_slices) {
  # Approximate working set:
  #   - cached history slices
  #   - threshold matrix
  #   - change matrix
  #   - combined matrix
  #   - some overhead / intermediate copies
  bytes <- as.numeric(nx) * as.numeric(ny) * 8 * (n_history_slices + 5)
  bytes / 1024^2
}

guard_memory <- function(structure, n_history_slices, settings, log_file = NULL) {
  peak_mb <- estimate_peak_window_memory_mb(structure$grid$nx, structure$grid$ny, n_history_slices)
  if (is.finite(peak_mb) && peak_mb > settings$max_window_memory_mb) {
    stop_with_log(
      paste0(
        "Memory guard triggered. Estimated peak window demand is about ",
        sprintf("%.1f", peak_mb),
        " MB, which exceeds max_window_memory_mb=",
        settings$max_window_memory_mb,
        ". Lower the change window, reduce depth scope, or increase the guard deliberately."
      ),
      log_file
    )
  }
  invisible(peak_mb)
}

# -------------------- Matrix -> raster / output writing --------------------

matrix_xy_to_raster <- function(mat_xy, grid, allow_curvilinear_resample = TRUE) {
  mat_xy <- coerce_xy_matrix(mat_xy, nx = grid$nx, ny = grid$ny, label = "mat_xy")

  if (grid$type == "regular") {
    x_ord <- order(grid$x)
    y_ord <- order(grid$y)
    mat_ord <- mat_xy[x_ord, y_ord, drop = FALSE]
    m <- t(mat_ord)
    m <- m[nrow(m):1, , drop = FALSE]  # top row should be highest y

    r <- terra::rast(m)
    terra::ext(r) <- c(grid$extent$xmin, grid$extent$xmax, grid$extent$ymin, grid$extent$ymax)
    if (!is_blank_scalar(grid$crs)) {
      terra::crs(r) <- grid$crs
    }
    return(r)
  }

  if (!allow_curvilinear_resample) {
    stop("Curvilinear grid detected but allow_curvilinear_resample is FALSE.", call. = FALSE)
  }

  lon_xy <- grid$lon_xy
  lat_xy <- grid$lat_xy

  if (!all(dim(mat_xy)[1:2] == dim(lon_xy)[1:2]) || !all(dim(mat_xy)[1:2] == dim(lat_xy)[1:2])) {
    stop("Curvilinear coordinate matrices do not match data slice dimensions.", call. = FALSE)
  }

  df <- data.frame(
    lon = as.numeric(lon_xy),
    lat = as.numeric(lat_xy),
    value = as.numeric(mat_xy)
  )
  df <- df[is.finite(df$lon) & is.finite(df$lat) & !is.na(df$value), , drop = FALSE]

  template <- terra::rast(
    nrows = grid$template$nrows,
    ncols = grid$template$ncols,
    xmin = grid$template$xmin,
    xmax = grid$template$xmax,
    ymin = grid$template$ymin,
    ymax = grid$template$ymax,
    crs = grid$crs
  )

  pts <- terra::vect(df, geom = c("lon", "lat"), crs = grid$crs)
  terra::rasterize(pts, template, field = "value", fun = "mean", background = NA_real_)
}

build_output_basename <- function(product_row, depth_label, window_label) {
  paste(
    sanitize_filename(product_row$product_var),
    sanitize_filename(product_row$alert_id),
    sanitize_filename(depth_label),
    sanitize_filename(window_label),
    sep = "__"
  )
}

write_geotiff <- function(r, filename, overwrite = TRUE, gdal_opts = c("COMPRESS=DEFLATE")) {
  terra::writeRaster(
    r,
    filename = filename,
    overwrite = overwrite,
    wopt = list(gdal = gdal_opts, datatype = "FLT4S")
  )
  invisible(filename)
}

# -------------------- PW4B grid integration outputs --------------------
# This block keeps the existing native-grid calculations and GeoTIFF outputs intact,
# then converts the finished native-grid matrices to the PW4B 1 km grid by using a
# cached nearest-neighbour / optional IDW lookup table.

require_pw4b_dependencies <- function(log_file = NULL) {
  needed <- c("sf", "RANN")
  missing <- needed[!vapply(needed, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop_with_log(
      paste0(
        "PW4B output generation requires missing package(s): ",
        paste(missing, collapse = ", "),
        ". Install them or run with --pw4b=FALSE."
      ),
      log_file
    )
  }
  invisible(TRUE)
}

first_matching_column <- function(dt, candidates) {
  nn <- names(dt)
  nn_lower <- tolower(trimws(nn))
  hit <- match(tolower(candidates), nn_lower, nomatch = 0L)
  hit <- hit[hit > 0L]
  if (length(hit) == 0L) NA_character_ else nn[hit[1L]]
}

detect_pw4b_coordinate_columns <- function(pw4b_dt, log_file = NULL) {
  if (ncol(pw4b_dt) >= 5L && all(names(pw4b_dt)[1:5] == paste0("V", 1:5))) {
    cellcode_vals <- as.character(pw4b_dt[[2L]])
    x_vals <- suppressWarnings(as.numeric(pw4b_dt[[3L]]))
    y_vals <- suppressWarnings(as.numeric(pw4b_dt[[4L]]))
    x_finite <- x_vals[is.finite(x_vals)]
    y_finite <- y_vals[is.finite(y_vals)]

    if (length(x_finite) > 0L && length(y_finite) > 0L &&
        any(grepl("^1kmE[0-9]+N[0-9]+$", cellcode_vals)) &&
        stats::median(abs(x_finite), na.rm = TRUE) > 100000 &&
        stats::median(abs(y_finite), na.rm = TRUE) > 100000) {
      log_message(
        "INFO",
        "Detected headerless PW4B grid layout; using V3 and V4 as ETRS LAEA x/y coordinates.",
        log_file
      )
      return(list(x_col = names(pw4b_dt)[3L], y_col = names(pw4b_dt)[4L]))
    }
  }

  x_col <- first_matching_column(
    pw4b_dt,
    c("pw4b_x", "point_x", "x", "xcoord", "x_coord", "x_coordinate", "easting", "etrs_laea_x", "laea_x")
  )
  y_col <- first_matching_column(
    pw4b_dt,
    c("pw4b_y", "point_y", "y", "ycoord", "y_coord", "y_coordinate", "northing", "etrs_laea_y", "laea_y")
  )

  if (is.na(x_col) || is.na(y_col) || identical(x_col, y_col)) {
    numeric_cols <- names(pw4b_dt)[vapply(
      pw4b_dt,
      function(v) is.numeric(v) || is.integer(v),
      FUN.VALUE = logical(1)
    )]
    numeric_cols <- setdiff(numeric_cols, c(x_col, y_col, NA_character_))
    if (length(numeric_cols) >= 2L) {
      x_col <- numeric_cols[1L]
      y_col <- numeric_cols[2L]
      log_message(
        "WARN",
        paste0(
          "PW4B coordinate columns were inferred as '", x_col,
          "' and '", y_col, "'. Verify these are ETRS LAEA x/y coordinates."
        ),
        log_file
      )
    }
  }

  if (is.na(x_col) || is.na(y_col) || identical(x_col, y_col)) {
    stop_with_log("Unable to identify PW4B x/y coordinate columns in the PW4B grid CSV.", log_file)
  }

  list(x_col = x_col, y_col = y_col)
}

detect_pw4b_cellcode_column <- function(pw4b_dt, log_file = NULL) {
  if (ncol(pw4b_dt) >= 2L && all(names(pw4b_dt)[1:2] == paste0("V", 1:2))) {
    cellcode_vals <- as.character(pw4b_dt[[2L]])
    if (any(grepl("^1kmE[0-9]+N[0-9]+$", cellcode_vals))) {
      log_message(
        "INFO",
        "Detected headerless PW4B grid layout; using V2 as the PW4B cellcode column.",
        log_file
      )
      return(names(pw4b_dt)[2L])
    }
  }

  cell_col <- first_matching_column(
    pw4b_dt,
    c("cellcode", "cell_code", "cellcode_1km", "gridcode", "grid_code", "cell_id", "cellid")
  )
  if (!is.na(cell_col)) {
    return(cell_col)
  }

  log_message(
    "WARN",
    "No PW4B cellcode column was detected; creating a row-number cellcode for PW4B CSV output.",
    log_file
  )
  NA_character_
}

detect_pw4b_landsea_column <- function(pw4b_dt, x_col, y_col, log_file = NULL) {
  coord_pos <- match(c(x_col, y_col), names(pw4b_dt))
  candidate_pos <- max(coord_pos, na.rm = TRUE) + 1L
  candidate <- if (candidate_pos <= ncol(pw4b_dt)) names(pw4b_dt)[candidate_pos] else NA_character_

  is_binary_col <- function(col_name) {
    if (is.na(col_name) || !col_name %in% names(pw4b_dt)) return(FALSE)
    vals <- suppressWarnings(as.numeric(pw4b_dt[[col_name]]))
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0L) return(FALSE)
    all(unique(vals) %in% c(0, 1))
  }

  if (!is.na(candidate) && is_binary_col(candidate)) {
    landsea_col <- candidate
  } else {
    landsea_col <- first_matching_column(
      pw4b_dt,
      c("land", "land_sea", "landsea", "is_land", "land_flag", "sea_land", "land_or_sea")
    )
    if (is.na(landsea_col) || !is_binary_col(landsea_col)) {
      binary_cols <- names(pw4b_dt)[vapply(names(pw4b_dt), is_binary_col, FUN.VALUE = logical(1))]
      binary_cols <- setdiff(binary_cols, c(x_col, y_col))
      landsea_col <- if (length(binary_cols) > 0L) binary_cols[1L] else NA_character_
    }
  }

  if (is.na(landsea_col)) {
    stop_with_log(
      "Unable to identify the PW4B land/sea field. The first field after the coordinate fields should be binary with land=1 and sea=0.",
      log_file
    )
  }

  vals <- suppressWarnings(as.numeric(pw4b_dt[[landsea_col]]))
  vals <- vals[is.finite(vals)]
  if (!all(unique(vals) %in% c(0, 1))) {
    stop_with_log(
      paste0(
        "PW4B land/sea field '", landsea_col,
        "' is not coded only as 0/1. Expected land=1 and sea=0."
      ),
      log_file
    )
  }

  log_message(
    "INFO",
    paste0("Using PW4B land/sea field '", landsea_col, "' with land=1 and sea=0."),
    log_file
  )
  landsea_col
}

pw4b_coordinate_spacing <- function(vals, default = 1000) {
  vals <- sort(unique(as.numeric(vals)))
  vals <- vals[is.finite(vals)]
  if (length(vals) < 2L) {
    return(default)
  }
  diffs <- diff(vals)
  diffs <- diffs[is.finite(diffs) & diffs > 0]
  if (length(diffs) == 0L) {
    return(default)
  }
  out <- stats::median(diffs, na.rm = TRUE)
  if (is.finite(out) && out > 0) out else default
}

pw4b_key <- function(x, y) {
  paste0(
    format(round(as.numeric(x), 3), scientific = FALSE, trim = TRUE),
    "|",
    format(round(as.numeric(y), 3), scientific = FALSE, trim = TRUE)
  )
}

select_pw4b_sea_and_coastal_cells <- function(pw4b_dt, x_col, y_col, landsea_col, log_file = NULL) {
  if (nrow(pw4b_dt) == 0L) {
    return(pw4b_dt)
  }

  landsea <- suppressWarnings(as.numeric(pw4b_dt[[landsea_col]]))
  sea_mask <- is.finite(landsea) & landsea == 0
  land_mask <- is.finite(landsea) & landsea == 1

  if (!any(sea_mask, na.rm = TRUE)) {
    log_message("WARN", "No sea cells were found in the selected PW4B extent.", log_file)
    return(pw4b_dt[FALSE])
  }

  dx <- pw4b_coordinate_spacing(pw4b_dt[[x_col]], default = 1000)
  dy <- pw4b_coordinate_spacing(pw4b_dt[[y_col]], default = 1000)

  sea_x <- as.numeric(pw4b_dt[[x_col]][sea_mask])
  sea_y <- as.numeric(pw4b_dt[[y_col]][sea_mask])
  offsets <- data.table::data.table(
    ox = c(-1, 0, 1, -1, 1, -1, 0, 1),
    oy = c(-1, -1, -1, 0, 0, 1, 1, 1)
  )
  neighbour_keys <- unique(unlist(lapply(seq_len(nrow(offsets)), function(i) {
    pw4b_key(sea_x + offsets$ox[i] * dx, sea_y + offsets$oy[i] * dy)
  }), use.names = FALSE))

  cell_keys <- pw4b_key(pw4b_dt[[x_col]], pw4b_dt[[y_col]])
  coastal_land_mask <- land_mask & cell_keys %in% neighbour_keys
  keep <- sea_mask | coastal_land_mask

  log_message(
    "INFO",
    paste0(
      "Selected ", sum(sea_mask, na.rm = TRUE), " PW4B sea cell(s) and ",
      sum(coastal_land_mask, na.rm = TRUE), " coastal land cell(s) inside the buffered NetCDF extent."
    ),
    log_file
  )

  pw4b_dt[keep]
}

build_native_grid_point_table <- function(structure, pw4b_crs = "EPSG:3035", log_file = NULL) {
  require_pw4b_dependencies(log_file = log_file)

  nx <- as.integer(structure$grid$nx)
  ny <- as.integer(structure$grid$ny)
  idx_dt <- data.table::data.table(
    x_index = rep(seq_len(nx), times = ny),
    y_index = rep(seq_len(ny), each = nx)
  )
  idx_dt[, nc_point_id := seq_len(.N)]

  if (structure$grid$type == "regular") {
    idx_dt[, nc_lon := as.numeric(structure$grid$x[x_index])]
    idx_dt[, nc_lat := as.numeric(structure$grid$y[y_index])]
  } else {
    idx_dt[, nc_lon := as.numeric(structure$grid$lon_xy)]
    idx_dt[, nc_lat := as.numeric(structure$grid$lat_xy)]
  }

  idx_dt <- idx_dt[is.finite(nc_lon) & is.finite(nc_lat)]
  if (nrow(idx_dt) == 0L) {
    stop_with_log("No finite native NetCDF longitude/latitude grid points were available for PW4B lookup creation.", log_file)
  }

  nc_sf <- sf::st_as_sf(idx_dt, coords = c("nc_lon", "nc_lat"), crs = 4326, remove = FALSE)
  nc_sf <- sf::st_transform(nc_sf, crs = pw4b_crs)
  xy <- sf::st_coordinates(nc_sf)
  idx_dt[, nc_x := as.numeric(xy[, 1])]
  idx_dt[, nc_y := as.numeric(xy[, 2])]
  idx_dt <- idx_dt[is.finite(nc_x) & is.finite(nc_y)]

  idx_dt
}

create_nc_to_pw4b_lookup <- function(
  structure,
  pw4b_grid_path,
  lookup_path,
  pw4b_crs = "EPSG:3035",
  buffer_m = 20000,
  idw_k = 4L,
  refresh = FALSE,
  log_file = NULL
) {
  require_pw4b_dependencies(log_file = log_file)

  lookup_path <- normalizePath(path.expand(lookup_path), winslash = "/", mustWork = FALSE)
  if (file.exists(lookup_path) && !isTRUE(refresh)) {
    log_message("INFO", paste0("Using existing PW4B lookup table: ", lookup_path), log_file)
    return(lookup_path)
  }

  if (is_blank_scalar(pw4b_grid_path) || !file.exists(pw4b_grid_path)) {
    stop_with_log(
      paste0(
        "PW4B grid CSV was not found at: ", pw4b_grid_path,
        ". Supply --pw4b_grid=... or run with --pw4b=FALSE."
      ),
      log_file
    )
  }

  nc_points <- build_native_grid_point_table(structure, pw4b_crs = pw4b_crs, log_file = log_file)
  bbox <- list(
    xmin = min(nc_points$nc_x, na.rm = TRUE) - buffer_m,
    xmax = max(nc_points$nc_x, na.rm = TRUE) + buffer_m,
    ymin = min(nc_points$nc_y, na.rm = TRUE) - buffer_m,
    ymax = max(nc_points$nc_y, na.rm = TRUE) + buffer_m
  )

  log_message("INFO", paste0("Reading PW4B grid CSV: ", pw4b_grid_path), log_file)
  pw4b_dt <- data.table::fread(pw4b_grid_path, showProgress = FALSE)
  if (nrow(pw4b_dt) == 0L) {
    stop_with_log("PW4B grid CSV was read successfully but contains no rows.", log_file)
  }

  coord_cols <- detect_pw4b_coordinate_columns(pw4b_dt, log_file = log_file)
  x_col <- coord_cols$x_col
  y_col <- coord_cols$y_col
  cellcode_col <- detect_pw4b_cellcode_column(pw4b_dt, log_file = log_file)
  if (is.na(cellcode_col)) {
    pw4b_dt[, cellcode := as.character(seq_len(.N))]
    cellcode_col <- "cellcode"
  }
  landsea_col <- detect_pw4b_landsea_column(pw4b_dt, x_col = x_col, y_col = y_col, log_file = log_file)

  pw4b_dt[, .pw4b_x_tmp := suppressWarnings(as.numeric(get(x_col)))]
  pw4b_dt[, .pw4b_y_tmp := suppressWarnings(as.numeric(get(y_col)))]
  pw4b_region <- pw4b_dt[
    is.finite(.pw4b_x_tmp) & is.finite(.pw4b_y_tmp) &
      .pw4b_x_tmp >= bbox$xmin & .pw4b_x_tmp <= bbox$xmax &
      .pw4b_y_tmp >= bbox$ymin & .pw4b_y_tmp <= bbox$ymax
  ]

  if (nrow(pw4b_region) == 0L) {
    stop_with_log("No PW4B grid cells fell inside the buffered NetCDF extent.", log_file)
  }

  selected_pw4b <- select_pw4b_sea_and_coastal_cells(
    pw4b_region,
    x_col = x_col,
    y_col = y_col,
    landsea_col = landsea_col,
    log_file = log_file
  )
  if (nrow(selected_pw4b) == 0L) {
    stop_with_log("No PW4B sea/coastal cells remained after applying the land/sea filter.", log_file)
  }

  idw_k <- max(1L, scalar_int(idw_k, default = 1L))
  idw_k <- min(idw_k, nrow(nc_points))

  nc_xy_matrix <- as.matrix(nc_points[, .(nc_x, nc_y)])
  pw4b_xy_matrix <- as.matrix(selected_pw4b[, .(.pw4b_x_tmp, .pw4b_y_tmp)])

  nearest <- RANN::nn2(
    data = nc_xy_matrix,
    query = pw4b_xy_matrix,
    k = idw_k
  )

  selected_pw4b[, pw4b_x := .pw4b_x_tmp]
  selected_pw4b[, pw4b_y := .pw4b_y_tmp]
  selected_pw4b[, nc_point_id := nc_points$nc_point_id[nearest$nn.idx[, 1L]]]
  selected_pw4b[, nearest_nc_x := nc_points$nc_x[nearest$nn.idx[, 1L]]]
  selected_pw4b[, nearest_nc_y := nc_points$nc_y[nearest$nn.idx[, 1L]]]
  selected_pw4b[, nearest_nc_lon := nc_points$nc_lon[nearest$nn.idx[, 1L]]]
  selected_pw4b[, nearest_nc_lat := nc_points$nc_lat[nearest$nn.idx[, 1L]]]
  selected_pw4b[, distance_to_nearest_nc_point_m := as.numeric(nearest$nn.dists[, 1L])]

  if (!"cellcode" %in% names(selected_pw4b)) {
    selected_pw4b[, cellcode := as.character(get(cellcode_col))]
  } else {
    selected_pw4b[, cellcode := as.character(cellcode)]
  }

  if (idw_k > 1L) {
    selected_pw4b[, idw_nc_point_ids := apply(nearest$nn.idx, 1L, function(ids) {
      paste(nc_points$nc_point_id[ids], collapse = ";")
    })]
    selected_pw4b[, idw_distances_m := apply(nearest$nn.dists, 1L, function(d) {
      paste(format(as.numeric(d), scientific = FALSE, trim = TRUE), collapse = ";")
    })]
  } else {
    selected_pw4b[, idw_nc_point_ids := as.character(nc_point_id)]
    selected_pw4b[, idw_distances_m := format(distance_to_nearest_nc_point_m, scientific = FALSE, trim = TRUE)]
  }

  selected_pw4b[, .pw4b_x_tmp := NULL]
  selected_pw4b[, .pw4b_y_tmp := NULL]

  attr(selected_pw4b, "pw4b_x_col") <- x_col
  attr(selected_pw4b, "pw4b_y_col") <- y_col
  attr(selected_pw4b, "pw4b_cellcode_col") <- cellcode_col
  attr(selected_pw4b, "pw4b_landsea_col") <- landsea_col
  attr(selected_pw4b, "pw4b_crs") <- pw4b_crs
  attr(selected_pw4b, "created_utc") <- format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC")
  attr(selected_pw4b, "buffer_m") <- buffer_m
  attr(selected_pw4b, "idw_k") <- idw_k

  safe_dir_create(dirname(lookup_path))
  saveRDS(selected_pw4b, lookup_path)

  log_message(
    "INFO",
    paste0(
      "Saved NetCDF-to-PW4B lookup table with ", nrow(selected_pw4b),
      " PW4B sea/coastal cell(s): ", lookup_path
    ),
    log_file
  )

  lookup_path
}

parse_semicolon_numeric <- function(x) {
  if (is_blank_scalar(x)) {
    return(numeric(0))
  }
  suppressWarnings(as.numeric(strsplit(as.character(x), ";", fixed = TRUE)[[1]]))
}

parse_semicolon_integer <- function(x) {
  if (is_blank_scalar(x)) {
    return(integer(0))
  }
  suppressWarnings(as.integer(strsplit(as.character(x), ";", fixed = TRUE)[[1]]))
}

convert_nc_output_to_pw4b <- function(
  mat_xy,
  structure,
  lookup_path,
  value_col = "combined_severity",
  depth = NULL,
  distance_threshold_m = 500,
  use_idw = TRUE,
  log_file = NULL,
  lookup = NULL
) {
  mat_xy <- coerce_xy_matrix(mat_xy, nx = structure$grid$nx, ny = structure$grid$ny, label = value_col)
  if (is.null(lookup)) {
    if (is.null(lookup_path) || !file.exists(lookup_path)) {
      stop_with_log(paste0("PW4B lookup table does not exist: ", lookup_path), log_file)
    }
    lookup <- data.table::as.data.table(readRDS(lookup_path))
  } else {
    lookup <- data.table::as.data.table(lookup)
  }
  required_cols <- c("cellcode", "pw4b_x", "pw4b_y", "nc_point_id", "distance_to_nearest_nc_point_m")
  missing_cols <- setdiff(required_cols, names(lookup))
  if (length(missing_cols) > 0L) {
    stop_with_log(
      paste0("PW4B lookup table is missing required columns: ", paste(missing_cols, collapse = ", ")),
      log_file
    )
  }

  nc_values <- as.numeric(mat_xy)
  nearest_values <- nc_values[as.integer(lookup$nc_point_id)]
  pw4b_values <- nearest_values

  can_idw <- isTRUE(use_idw) && all(c("idw_nc_point_ids", "idw_distances_m") %in% names(lookup))
  far_rows <- which(
    can_idw &
      is.finite(lookup$distance_to_nearest_nc_point_m) &
      lookup$distance_to_nearest_nc_point_m > distance_threshold_m
  )

  if (length(far_rows) > 0L) {
    for (row_idx in far_rows) {
      ids <- parse_semicolon_integer(lookup$idw_nc_point_ids[row_idx])
      distances <- parse_semicolon_numeric(lookup$idw_distances_m[row_idx])
      n <- min(length(ids), length(distances))
      if (n == 0L) next

      ids <- ids[seq_len(n)]
      distances <- distances[seq_len(n)]
      vals <- nc_values[ids]
      valid <- is.finite(vals) & is.finite(distances) & distances >= 0
      if (!any(valid)) {
        pw4b_values[row_idx] <- NA_real_
        next
      }

      vals <- vals[valid]
      distances <- distances[valid]
      zero_distance <- distances == 0
      if (any(zero_distance)) {
        pw4b_values[row_idx] <- vals[which(zero_distance)[1L]]
      } else {
        weights <- 1 / distances^2
        weights <- weights / sum(weights, na.rm = TRUE)
        pw4b_values[row_idx] <- sum(weights * vals, na.rm = TRUE)
      }
    }
  }

  out <- data.table::data.table(
    cellcode = as.character(lookup$cellcode),
    pw4b_x = as.numeric(lookup$pw4b_x),
    pw4b_y = as.numeric(lookup$pw4b_y)
  )
  out[, (value_col) := pw4b_values]

  if (!is.null(depth) && !is_blank_scalar(depth)) {
    out[, depth := as.character(depth)]
    data.table::setcolorder(out, c("cellcode", "pw4b_x", "pw4b_y", "depth", value_col))
  }

  out
}

convert_combined_severity_to_pw4b_csv <- function(
  mat_xy,
  structure,
  lookup_path,
  depth = NULL,
  log_file = NULL,
  lookup = NULL
) {
  value_col <- "combined_severity"
  fill_method_col <- "combined_severity_fill_method"
  mat_xy <- coerce_xy_matrix(mat_xy, nx = structure$grid$nx, ny = structure$grid$ny, label = value_col)
  if (is.null(lookup)) {
    if (is.null(lookup_path) || !file.exists(lookup_path)) {
      stop_with_log(paste0("PW4B lookup table does not exist: ", lookup_path), log_file)
    }
    lookup <- data.table::as.data.table(readRDS(lookup_path))
  } else {
    lookup <- data.table::as.data.table(lookup)
  }
  required_cols <- c("cellcode", "pw4b_x", "pw4b_y", "nc_point_id", "distance_to_nearest_nc_point_m")
  missing_cols <- setdiff(required_cols, names(lookup))
  if (length(missing_cols) > 0L) {
    stop_with_log(
      paste0("PW4B lookup table is missing required columns: ", paste(missing_cols, collapse = ", ")),
      log_file
    )
  }

  nc_values <- as.numeric(mat_xy)
  pw4b_values <- nc_values[as.integer(lookup$nc_point_id)]
  fill_method <- rep("", nrow(lookup))
  has_nearby_cols <- all(c("idw_nc_point_ids", "idw_distances_m") %in% names(lookup))

  rows_needing_valid_nc_value <- which(!is.finite(pw4b_values))
  if (length(rows_needing_valid_nc_value) > 0L) {
    for (row_idx in rows_needing_valid_nc_value) {
      candidate_ids <- as.integer(lookup$nc_point_id[row_idx])
      candidate_distances <- as.numeric(lookup$distance_to_nearest_nc_point_m[row_idx])

      if (has_nearby_cols) {
        nearby_ids <- parse_semicolon_integer(lookup$idw_nc_point_ids[row_idx])
        nearby_distances <- parse_semicolon_numeric(lookup$idw_distances_m[row_idx])
        n <- min(length(nearby_ids), length(nearby_distances))
        if (n > 0L) {
          candidate_ids <- c(candidate_ids, nearby_ids[seq_len(n)])
          candidate_distances <- c(candidate_distances, nearby_distances[seq_len(n)])
        }
      }

      candidates <- data.table::data.table(
        nc_point_id = candidate_ids,
        distance_m = candidate_distances
      )
      candidates <- candidates[
        is.finite(nc_point_id) & nc_point_id >= 1L & nc_point_id <= length(nc_values) &
          is.finite(distance_m) & distance_m >= 0
      ]

      if (nrow(candidates) > 0L) {
        candidates[, nc_value := nc_values[as.integer(nc_point_id)]]
        candidates <- candidates[is.finite(nc_value)]
      }

      if (nrow(candidates) > 0L) {
        data.table::setorder(candidates, distance_m)
        pw4b_values[row_idx] <- candidates$nc_value[1L]
      } else {
        pw4b_values[row_idx] <- 0
        fill_method[row_idx] <- "forced_zero_no_valid_nc_value"
      }
    }
  }

  out <- data.table::data.table(
    cellcode = as.character(lookup$cellcode),
    pw4b_x = as.numeric(lookup$pw4b_x),
    pw4b_y = as.numeric(lookup$pw4b_y)
  )
  out[, (value_col) := pw4b_values]
  out[, (fill_method_col) := fill_method]

  if (!is.null(depth) && !is_blank_scalar(depth)) {
    out[, depth := as.character(depth)]
    data.table::setcolorder(out, c("cellcode", "pw4b_x", "pw4b_y", "depth", value_col, fill_method_col))
  }

  out
}

pw4b_table_to_raster <- function(pw4b_dt, value_col, crs = "EPSG:3035") {
  if (nrow(pw4b_dt) == 0L) {
    stop("Cannot build a PW4B raster from an empty table.", call. = FALSE)
  }

  x <- as.numeric(pw4b_dt$pw4b_x)
  y <- as.numeric(pw4b_dt$pw4b_y)
  values <- as.numeric(pw4b_dt[[value_col]])
  keep <- is.finite(x) & is.finite(y)
  if (!any(keep)) {
    stop("Cannot build a PW4B raster because no finite PW4B x/y coordinates are available.", call. = FALSE)
  }

  x <- x[keep]
  y <- y[keep]
  values <- values[keep]

  dx <- pw4b_coordinate_spacing(x, default = 1000)
  dy <- pw4b_coordinate_spacing(y, default = 1000)

  xmin <- min(x, na.rm = TRUE) - dx / 2
  xmax <- max(x, na.rm = TRUE) + dx / 2
  ymin <- min(y, na.rm = TRUE) - dy / 2
  ymax <- max(y, na.rm = TRUE) + dy / 2

  template <- terra::rast(
    nrows = max(1L, as.integer(round((ymax - ymin) / dy))),
    ncols = max(1L, as.integer(round((xmax - xmin) / dx))),
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax,
    crs = crs
  )

  pts_df <- data.frame(pw4b_x = x, pw4b_y = y, value = values)
  pts <- terra::vect(pts_df, geom = c("pw4b_x", "pw4b_y"), crs = crs)
  r <- terra::rasterize(pts, template, field = "value", fun = "mean", background = NA_real_)
  names(r) <- value_col
  r
}

write_pw4b_outputs <- function(
  matrices,
  structure,
  settings,
  base_name,
  output_dir,
  depth_label = NULL,
  log_file = NULL
) {
  empty <- list(
    pw4b_threshold_file = NA_character_,
    pw4b_change_file = NA_character_,
    pw4b_combined_file = NA_character_,
    pw4b_max_abs_change_file = NA_character_,
    pw4b_combined_csv_file = NA_character_,
    pw4b_lookup_file = settings$pw4b_lookup_path %||% NA_character_
  )

  if (!isTRUE(settings$write_pw4b)) {
    return(empty)
  }

  lookup_path <- settings$pw4b_lookup_path
  pw4b_lookup <- NULL
  if ("pw4b_lookup" %in% names(settings) && !is.null(settings$pw4b_lookup)) {
    pw4b_lookup <- data.table::as.data.table(settings$pw4b_lookup)
  }
  if (!file.exists(lookup_path)) {
    lookup_path <- create_nc_to_pw4b_lookup(
      structure = structure,
      pw4b_grid_path = settings$pw4b_grid_path,
      lookup_path = settings$pw4b_lookup_path,
      pw4b_crs = settings$pw4b_crs,
      buffer_m = settings$pw4b_buffer_m,
      idw_k = settings$pw4b_idw_k,
      refresh = FALSE,
      log_file = log_file
    )
    pw4b_lookup <- NULL
  }
  if (is.null(pw4b_lookup)) {
    pw4b_lookup <- data.table::as.data.table(readRDS(lookup_path))
  }

  depth_for_csv <- if (!is_blank_scalar(depth_label) && !identical(depth_label, "nodepth")) {
    depth_label
  } else {
    NULL
  }

  metric_map <- list(
    threshold_severity = "pw4b_threshold_file",
    change_severity = "pw4b_change_file",
    combined_severity = "pw4b_combined_file",
    max_abs_change = "pw4b_max_abs_change_file"
  )

  out_files <- empty
  out_files$pw4b_lookup_file <- lookup_path

  for (metric_name in names(metric_map)) {
    if (!metric_name %in% names(matrices) || is.null(matrices[[metric_name]])) {
      next
    }

    pw4b_dt <- convert_nc_output_to_pw4b(
      matrices[[metric_name]],
      structure = structure,
      lookup_path = lookup_path,
      value_col = metric_name,
      depth = depth_for_csv,
      distance_threshold_m = settings$pw4b_nn_distance_m,
      use_idw = settings$pw4b_idw_k > 1L,
      log_file = log_file,
      lookup = pw4b_lookup
    )

    pw4b_r <- pw4b_table_to_raster(pw4b_dt, value_col = metric_name, crs = settings$pw4b_crs)
    pw4b_tif <- file.path(output_dir, paste0(base_name, "__pw4b_", metric_name, ".tif"))
    write_geotiff(pw4b_r, pw4b_tif, overwrite = settings$overwrite, gdal_opts = settings$gdal_write_options)
    out_files[[metric_map[[metric_name]]]] <- pw4b_tif

    if (identical(metric_name, "combined_severity")) {
      pw4b_csv_dt <- convert_combined_severity_to_pw4b_csv(
        matrices[[metric_name]],
        structure = structure,
        lookup_path = lookup_path,
        depth = depth_for_csv,
        log_file = log_file,
        lookup = pw4b_lookup
      )
      pw4b_csv <- file.path(output_dir, paste0(base_name, "__pw4b_combined_severity.csv"))
      csv_cols <- c("cellcode", "pw4b_x", "pw4b_y")
      if ("depth" %in% names(pw4b_csv_dt)) {
        csv_cols <- c(csv_cols, "depth")
      }
      csv_cols <- c(csv_cols, "combined_severity", "combined_severity_fill_method")
      data.table::fwrite(pw4b_csv_dt[, ..csv_cols], pw4b_csv, na = "")
      out_files$pw4b_combined_csv_file <- pw4b_csv
    }
  }

  out_files
}

# -------------------- OPTIONAL DEBUG PNG MAPS (separable) --------------------
# This block is intentionally isolated so it can be deleted or disabled later
# without affecting the core analytical workflow.

make_coastline_sf <- function(crs_value = "EPSG:4326", bbox = NULL, log_file = NULL) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    return(NULL)
  }
  if (is_blank_scalar(crs_value) || !grepl("4326", crs_value)) {
    return(NULL)
  }

  bbox_values <- function(bb) {
    if (is.null(bb)) {
      return(NULL)
    }

    vals <- tryCatch(
      c(
        xmin = as.numeric(bb$xmin),
        xmax = as.numeric(bb$xmax),
        ymin = as.numeric(bb$ymin),
        ymax = as.numeric(bb$ymax)
      ),
      error = function(e) rep(NA_real_, 4L)
    )

    if (any(!is.finite(vals))) {
      vals <- tryCatch(
        c(
          xmin = as.numeric(terra::xmin(bb)),
          xmax = as.numeric(terra::xmax(bb)),
          ymin = as.numeric(terra::ymin(bb)),
          ymax = as.numeric(terra::ymax(bb))
        ),
        error = function(e) rep(NA_real_, 4L)
      )
    }

    if (any(!is.finite(vals))) {
      return(NULL)
    }
    vals
  }

  bb <- bbox_values(bbox)
  if (is.null(bb)) {
    return(NULL)
  }

  bb_poly <- sf::st_as_sfc(sf::st_bbox(c(
    xmin = bb["xmin"], xmax = bb["xmax"],
    ymin = bb["ymin"], ymax = bb["ymax"]
  ), crs = sf::st_crs(4326)))

  clip_land_to_bbox <- function(land_sf, source_name) {
    if (is.null(land_sf) || !inherits(land_sf, "sf") || nrow(land_sf) == 0L) {
      return(NULL)
    }

    if (is.na(sf::st_crs(land_sf))) {
      land_sf <- sf::st_set_crs(land_sf, 4326)
    } else if (!identical(sf::st_crs(land_sf), sf::st_crs(4326))) {
      land_sf <- sf::st_transform(land_sf, 4326)
    }

    land_sf <- tryCatch(
      suppressWarnings(sf::st_make_valid(land_sf)),
      error = function(e) land_sf
    )
    land_sf <- land_sf[!sf::st_is_empty(land_sf), , drop = FALSE]
    if (nrow(land_sf) == 0L) {
      return(NULL)
    }

    geom_type <- as.character(sf::st_geometry_type(land_sf, by_geometry = TRUE))
    keep <- geom_type %in% c("POLYGON", "MULTIPOLYGON")
    if (!any(keep)) {
      land_sf <- tryCatch(
        suppressWarnings(sf::st_collection_extract(land_sf, "POLYGON")),
        error = function(e) land_sf
      )
      geom_type <- as.character(sf::st_geometry_type(land_sf, by_geometry = TRUE))
      keep <- geom_type %in% c("POLYGON", "MULTIPOLYGON")
    }
    land_sf <- land_sf[keep, , drop = FALSE]
    if (nrow(land_sf) == 0L) {
      return(NULL)
    }

    land_sf_full <- land_sf
    land_sf <- tryCatch(
      suppressWarnings(sf::st_crop(land_sf, sf::st_bbox(bb_poly))),
      error = function(e) {
        tryCatch(
          suppressWarnings(sf::st_intersection(land_sf_full, bb_poly)),
          error = function(e2) {
            log_message(
              "WARN",
              paste0(
                "Could not clip land overlay to the forecast bounding box; using unclipped land overlay instead. ",
                conditionMessage(e2)
              ),
              log_file
            )
            land_sf_full
          }
        )
      }
    )

    land_sf <- land_sf[!sf::st_is_empty(land_sf), , drop = FALSE]
    if (nrow(land_sf) == 0L) {
      land_sf <- land_sf_full
    }
    if (nrow(land_sf) == 0L) {
      return(NULL)
    }

    attr(land_sf, "debug_land_source") <- source_name
    land_sf
  }

  read_etak_land_layer <- function(layer_name) {
    bbox_param <- paste(
      format(c(bb["xmin"], bb["ymin"], bb["xmax"], bb["ymax"]), scientific = FALSE, trim = TRUE),
      collapse = ","
    )
    page_size <- 5000L
    max_pages <- 10L
    layer_pages <- list()

    for (page_idx in seq_len(max_pages)) {
      start_index <- (page_idx - 1L) * page_size
      url <- paste0(
        "https://gsavalik.envir.ee/geoserver/etak/ows?",
        "service=WFS&version=2.0.0&request=GetFeature",
        "&typeName=", utils::URLencode(layer_name, reserved = TRUE),
        "&outputFormat=application/json",
        "&srsName=EPSG:4326",
        "&count=", page_size,
        "&startIndex=", start_index,
        "&bbox=", bbox_param, ",EPSG:4326"
      )

      page_sf <- tryCatch(
        suppressWarnings(sf::st_read(url, quiet = TRUE, stringsAsFactors = FALSE)),
        error = function(e) NULL
      )
      if (is.null(page_sf) || !inherits(page_sf, "sf") || nrow(page_sf) == 0L) {
        break
      }

      layer_pages[[length(layer_pages) + 1L]] <- page_sf
      if (nrow(page_sf) < page_size) {
        break
      }

      if (page_idx == max_pages) {
        log_message(
          "WARN",
          paste0("Reached the WFS paging limit while reading ETAK layer '", layer_name, "'. Land overlay may be incomplete."),
          log_file
        )
      }
    }

    if (length(layer_pages) == 0L) {
      return(NULL)
    }

    layer_sf <- tryCatch(
      do.call(rbind, layer_pages),
      error = function(e) layer_pages[[1L]]
    )
    if (is.null(layer_sf) || !inherits(layer_sf, "sf") || nrow(layer_sf) == 0L) {
      return(NULL)
    }

    if (is.na(sf::st_crs(layer_sf))) {
      layer_sf <- sf::st_set_crs(layer_sf, 4326)
    } else if (!identical(sf::st_crs(layer_sf), sf::st_crs(4326))) {
      layer_sf <- sf::st_transform(layer_sf, 4326)
    }

    layer_sf <- layer_sf[!sf::st_is_empty(layer_sf), , drop = FALSE]
    if (nrow(layer_sf) == 0L) {
      return(NULL)
    }

    sf::st_sf(
      source_layer = rep(layer_name, nrow(layer_sf)),
      geometry = sf::st_geometry(layer_sf),
      crs = sf::st_crs(layer_sf)
    )
  }

  # Estonian Topographic Database land-cover polygon layers from the public
  # Maa- ja Ruumiamet / ETAK WFS service. These give a much more detailed
  # coastline and island outline around Estonia than maps::map("world").
  etak_land_layers <- c(
    "etak:e_301_muu_kolvik_a",
    "etak:e_302_ou_a",
    "etak:e_303_haritav_maa_a",
    "etak:e_304_lage_a",
    "etak:e_305_puittaimestik_a",
    "etak:e_306_margala_a",
    "etak:e_307_turbavali_a",
    "etak:e_501_tee_a"
  )

  etak_parts <- lapply(etak_land_layers, read_etak_land_layer)
  etak_parts <- etak_parts[!vapply(etak_parts, is.null, FUN.VALUE = logical(1))]

  if (length(etak_parts) > 0L) {
    etak_land <- tryCatch(
      do.call(rbind, etak_parts),
      error = function(e) NULL
    )
    etak_land <- clip_land_to_bbox(
      etak_land,
      "Land overlay: Estonian Topographic Database, Maa- ja Ruumiamet"
    )
    if (!is.null(etak_land)) {
      log_message("INFO", "Using Estonian Topographic Database land polygons for PNG land overlay.", log_file)
      return(etak_land)
    }
  }

  if (!requireNamespace("maps", quietly = TRUE)) {
    log_message("WARN", "Could not build detailed Estonian land overlay and maps package is unavailable for fallback.", log_file)
    return(NULL)
  }

  world <- tryCatch(
    sf::st_as_sf(maps::map("world", plot = FALSE, fill = TRUE)),
    error = function(e) {
      log_message("WARN", paste0("Could not build fallback coastline layer: ", conditionMessage(e)), log_file)
      NULL
    }
  )
  world <- clip_land_to_bbox(world, "Land overlay: fallback maps::map('world')")
  if (is.null(world)) {
    log_message("WARN", "No land polygons were available for PNG overlay.", log_file)
  } else {
    log_message("WARN", "Detailed Estonian land overlay unavailable; using fallback maps::map('world') polygons.", log_file)
  }
  world
}

debug_format_setting_value <- function(x, suffix = "") {
  if (is_blank_scalar(x)) {
    return("unavailable")
  }
  x1 <- x[[1]]
  if (is.na(x1)) {
    return("unavailable")
  }
  num <- suppressWarnings(as.numeric(x1))
  if (length(num) > 0L && is.finite(num)) {
    return(paste0(format(num, trim = TRUE, scientific = FALSE), suffix))
  }
  out <- trimws(as.character(x1))
  if (!nzchar(out)) "unavailable" else paste0(out, suffix)
}

debug_lookup_units <- function(var_name, structure = NULL) {
  if (is.null(structure) || is.null(structure$variables) || is_blank_scalar(var_name)) {
    return(NA_character_)
  }
  vars <- structure$variables
  if (!all(c("var_name", "units") %in% names(vars))) {
    return(NA_character_)
  }
  vn <- trimws(as.character(var_name[[1]]))
  idx <- which(tolower(trimws(as.character(vars$var_name))) == tolower(vn))
  if (length(idx) == 0L) {
    return(NA_character_)
  }
  units <- scalar_chr(vars$units[idx[1]], default = NA_character_)
  if (is_blank_scalar(units)) NA_character_ else units
}

debug_product_units <- function(product_row = NULL, structure = NULL) {
  if (is.null(product_row)) {
    return("unavailable")
  }

  for (unit_col in c("units", "unit", "product_units")) {
    if (unit_col %in% names(product_row) && !is_blank_scalar(product_row[[unit_col]])) {
      return(scalar_chr(product_row[[unit_col]], default = "unavailable"))
    }
  }

  source_units <- unique(stats::na.omit(c(
    debug_lookup_units(product_row$source_var_1, structure),
    debug_lookup_units(product_row$source_var_2, structure),
    debug_lookup_units(product_row$product_var, structure)
  )))
  source_units <- source_units[nzchar(trimws(source_units))]
  if (length(source_units) == 1L) {
    if (tolower(source_units[1]) %in% c("1", "fraction", "unitless") &&
        tolower(scalar_chr(product_row$product_var, default = "")) %in% c("siconc", "icefrac")) {
      return("fraction")
    }
    return(source_units[1])
  }
  if (length(source_units) > 1L) {
    return(paste(source_units, collapse = " / "))
  }

  product_var <- tolower(scalar_chr(product_row$product_var, default = ""))
  fallback_units <- c(
    sla = "m",
    zos = "m",
    siconc = "fraction",
    sithick = "m",
    icefrac = "fraction",
    icethic = "m",
    so = "PSU",
    thetao = "°C",
    current_speed = "m s⁻¹",
    ice_speed = "m s⁻¹",
    p0 = "count",
    p1 = "count",
    p2 = "count",
    p3 = "count",
    pt = "seconds",
    pc = "impact proxy units",
    pc1 = "count",
    uo = "m s⁻¹",
    vo = "m s⁻¹",
    siu = "m s⁻¹",
    siv = "m s⁻¹",
    uice = "m s⁻¹",
    vice = "m s⁻¹"
  )
  if (product_var %in% names(fallback_units)) fallback_units[[product_var]] else "unavailable"
}

debug_alert_settings_lines <- function(product_row = NULL, structure = NULL) {
  if (is.null(product_row)) {
    return(character(0))
  }

  change_threshold <- debug_format_setting_value(product_row$change_threshold)
  change_days <- debug_format_setting_value(product_row$change_days)
  change_suffix <- if (!identical(change_threshold, "unavailable") && !identical(change_days, "unavailable")) {
    paste0(" over ", change_days, " d")
  } else {
    ""
  }

  c(
    paste0("High threshold: ", debug_format_setting_value(product_row$high_threshold)),
    paste0("Low threshold: ", debug_format_setting_value(product_row$low_threshold)),
    paste0("Change threshold: ", change_threshold, change_suffix),
    paste0("Change direction: ", debug_format_setting_value(product_row$change_direction)),
    paste0("Units: ", debug_product_units(product_row, structure))
  )
}

debug_trigger_info <- function(combined_r, threshold_r = NULL, change_r = NULL) {
  combined_vals <- terra::values(combined_r, mat = FALSE)
  n <- length(combined_vals)

  if (is.null(threshold_r) || is.null(change_r)) {
    return(list(r = NULL, labels = character(0), cols = character(0)))
  }

  threshold_vals <- terra::values(threshold_r, mat = FALSE)
  change_vals <- terra::values(change_r, mat = FALSE)

  if (length(threshold_vals) != n || length(change_vals) != n) {
    return(list(r = NULL, labels = character(0), cols = character(0)))
  }

  threshold_hit <- is.finite(threshold_vals) & threshold_vals > 0
  change_hit <- is.finite(change_vals) & change_vals > 0

  trigger_vals <- rep(NA_real_, n)
  trigger_vals[is.finite(combined_vals)] <- 0
  trigger_vals[threshold_hit & !change_hit] <- 1
  trigger_vals[!threshold_hit & change_hit] <- 2
  trigger_vals[threshold_hit & change_hit] <- 3

  trigger_r <- combined_r
  terra::values(trigger_r) <- trigger_vals

  trigger_cols <- c("#FFFFFF", "#6A51A3", "#238B45", "#08519C")
  trigger_counts <- c(
    no_alert = sum(trigger_vals == 0, na.rm = TRUE),
    threshold_only = sum(trigger_vals == 1, na.rm = TRUE),
    change_only = sum(trigger_vals == 2, na.rm = TRUE),
    both = sum(trigger_vals == 3, na.rm = TRUE)
  )
  trigger_labels <- c(
    paste0("No alert: ", trigger_counts["no_alert"], " cells"),
    paste0("Threshold only: ", trigger_counts["threshold_only"], " cells"),
    paste0("Change only: ", trigger_counts["change_only"], " cells"),
    paste0("Both: ", trigger_counts["both"], " cells")
  )

  list(
    r = trigger_r,
    labels = trigger_labels,
    cols = trigger_cols,
    breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5)
  )
}

debug_plot_land_overlay <- function(land_sf, plot_extent = NULL) {
  if (!is.null(land_sf) && inherits(land_sf, "sf")) {
    tryCatch(
      plot(sf::st_geometry(land_sf), add = TRUE, col = "grey90", border = "grey40", lwd = 0.7),
      error = function(e) NULL
    )
    return(invisible(TRUE))
  }

  if (requireNamespace("maps", quietly = TRUE)) {
    ext_vals <- tryCatch(
      c(
        xmin = as.numeric(terra::xmin(plot_extent)),
        xmax = as.numeric(terra::xmax(plot_extent)),
        ymin = as.numeric(terra::ymin(plot_extent)),
        ymax = as.numeric(terra::ymax(plot_extent))
      ),
      error = function(e) rep(NA_real_, 4L)
    )
    xlim <- if (all(is.finite(ext_vals[c("xmin", "xmax")]))) ext_vals[c("xmin", "xmax")] else NULL
    ylim <- if (all(is.finite(ext_vals[c("ymin", "ymax")]))) ext_vals[c("ymin", "ymax")] else NULL

    tryCatch(
      maps::map("world", add = TRUE, fill = TRUE, col = "grey90", border = "grey40", lwd = 0.7, xlim = xlim, ylim = ylim),
      error = function(e) NULL
    )
  }

  invisible(NULL)
}

plot_debug_png <- function(
  combined_r,
  png_file,
  title_lines,
  coast_sf = NULL,
  threshold_r = NULL,
  change_r = NULL,
  product_row = NULL,
  structure = NULL
) {
  safe_dir_create(dirname(png_file))

  vals <- terra::values(combined_r, mat = FALSE)
  vals_bad <- vals[is.finite(vals) & vals > 0]
  has_bad <- length(vals_bad) > 0L
  max_severity <- if (has_bad) max(vals_bad, na.rm = TRUE) else 0

  png(filename = png_file, width = 2400, height = 1800, res = 220, bg = "white")
  on.exit(dev.off(), add = TRUE)

  wrap_title_lines <- function(lines, width = 90L) {
    wrapped <- unlist(lapply(lines, function(x) strwrap(as.character(x), width = width)))
    if (length(wrapped) == 0L) wrapped <- ""
    wrapped
  }

  class_r <- combined_r
  combined_vals_for_png <- terra::values(combined_r, mat = FALSE)
  class_vals <- rep(NA_real_, length(combined_vals_for_png))
  finite_vals <- is.finite(combined_vals_for_png)
  class_vals[finite_vals & combined_vals_for_png <= 0] <- 0
  class_vals[finite_vals & combined_vals_for_png > 0 & combined_vals_for_png <= 0.5] <- 1
  class_vals[finite_vals & combined_vals_for_png > 0.5 & combined_vals_for_png <= 1] <- 2
  class_vals[finite_vals & combined_vals_for_png > 1] <- 3
  terra::values(class_r) <- class_vals

  class_breaks <- c(-0.5, 0.5, 1.5, 2.5, 3.5)
  class_cols <- c("#FFFFFF", "#FEE5D9", "#FC9272", "#CB181D")
  class_labels <- c(
    "No alert: warning level not exceeded",
    "Low: up to 50% above warning level",
    "Medium: 50–100% above warning level",
    "High: >100% above warning level"
  )

  trigger_info <- debug_trigger_info(combined_r, threshold_r = threshold_r, change_r = change_r)
  settings_lines <- debug_alert_settings_lines(product_row = product_row, structure = structure)
  settings_lines <- c(
    settings_lines,
    paste0("Max combined severity: ", format(max_severity, digits = 4, trim = TRUE))
  )

  land_source <- if (!is.null(coast_sf) && inherits(coast_sf, "sf")) attr(coast_sf, "debug_land_source") else NA_character_
  if (is_blank_scalar(land_source) && requireNamespace("maps", quietly = TRUE)) {
    land_source <- "Land overlay: fallback maps::map('world')"
  }
  if (is_blank_scalar(land_source)) {
    land_source <- "Land overlay: unavailable"
  }

  graphics::layout(
    matrix(c(1, 2,
             1, 3,
             1, 4,
             1, 5), nrow = 4, byrow = TRUE),
    widths = c(5.3, 2.05),
    heights = c(1.05, 1.05, 1.15, 1.55)
  )
  par(oma = c(3.4, 2.0, 9.2, 2.0), xaxs = "i", yaxs = "i")

  par(mar = c(1.2, 1.2, 1.2, 0.8))
  terra::plot(
    class_r,
    col = class_cols,
    breaks = class_breaks,
    legend = FALSE,
    axes = FALSE
  )
  debug_plot_land_overlay(coast_sf, plot_extent = class_r)
  box(which = "plot", lwd = 1.2, col = "black")

  par(mar = c(0.4, 0.2, 0.8, 0.2))
  plot.new()
  legend(
    "topleft",
    inset = c(0.02, 0.02),
    legend = class_labels,
    fill = class_cols,
    border = "black",
    bg = "white",
    bty = "o",
    cex = 0.8,
    y.intersp = 1.15,
    title = "Combined severity"
  )
  box(which = "figure", lwd = 1.0, col = "black")

  par(mar = c(0.4, 0.2, 0.8, 0.2))
  plot.new()
  if (length(trigger_info$labels) > 0L) {
    legend(
      "topleft",
      inset = c(0.02, 0.02),
      legend = trigger_info$labels,
      fill = trigger_info$cols,
      border = "black",
      bg = "white",
      bty = "o",
      cex = 0.8,
      y.intersp = 1.15,
      title = "Trigger type"
    )
  } else {
    text(0.5, 0.5, labels = "Trigger type unavailable", cex = 0.85)
  }
  box(which = "figure", lwd = 1.0, col = "black")

  par(mar = c(0.4, 0.2, 0.8, 0.2))
  plot.new()
  if (length(settings_lines) > 0L) {
    legend(
      "topleft",
      inset = c(0.02, 0.02),
      legend = settings_lines,
      bg = "white",
      bty = "o",
      cex = 0.76,
      y.intersp = 1.08,
      title = "Alert settings"
    )
  } else {
    text(0.5, 0.5, labels = "Alert settings unavailable", cex = 0.85)
  }
  box(which = "figure", lwd = 1.0, col = "black")

  par(mar = c(0.7, 0.25, 1.1, 0.25))
  if (!is.null(trigger_info$r)) {
    terra::plot(
      trigger_info$r,
      col = trigger_info$cols,
      breaks = trigger_info$breaks,
      legend = FALSE,
      axes = FALSE
    )
    debug_plot_land_overlay(coast_sf, plot_extent = trigger_info$r)
    box(which = "plot", lwd = 0.9, col = "black")
    title(main = "Trigger type map", line = 0.15, cex.main = 0.82, font.main = 2)
  } else {
    plot.new()
    text(0.5, 0.5, labels = "Trigger type map unavailable", cex = 0.85)
    box(which = "figure", lwd = 1.0, col = "black")
  }

  title_text <- paste(wrap_title_lines(title_lines, width = 88L), collapse = "\n")
  mtext(title_text, side = 3, outer = TRUE, line = 2.2, cex = 0.95, font = 2)
  mtext(
    paste0(
      "Combined severity = max(threshold severity, change severity). Filled land is drawn over alert pixels. ",
      land_source, "."
    ),
    side = 1,
    outer = TRUE,
    line = 0.8,
    cex = 0.58
  )
}

# ----------------------- Window-level processing -----------------------

compute_window_severity_matrices <- function(
  nc,
  structure,
  product_row,
  depth_target,
  window_row,
  settings,
  log_file = NULL
) {
  time_idx <- unlist(window_row$time_indices)
  if (length(time_idx) == 0L) {
    return(NULL)
  }

  change_steps <- change_steps_from_days(product_row$change_days, steps_per_day_from_structure(structure))
  needed_idx <- seq(max(1L, min(time_idx) - change_steps), max(time_idx))

  depth_indices_for_read <- NULL
  if (depth_target$kind == "single" && !is.na(depth_target$depth_index)) {
    depth_indices_for_read <- depth_target$depth_index
  } else if (depth_target$kind == "summary") {
    depth_indices_for_read <- depth_target$depth_indices
  }

  depth_block_count <- if (!is.null(depth_indices_for_read) && length(depth_indices_for_read) > 0L) {
    as.integer(max(depth_indices_for_read) - min(depth_indices_for_read) + 1L)
  } else {
    1L
  }
  source_factor <- if (product_row$source_type == "vector_magnitude") 2L else 1L
  guard_memory(
    structure,
    length(needed_idx) * depth_block_count * source_factor,
    settings,
    log_file = log_file
  )

  if (depth_target$kind == "single") {
    block <- read_product_block_xyzt(
      nc,
      structure,
      product_row,
      time_indices = needed_idx,
      depth_indices = depth_indices_for_read,
      log_file = log_file
    )
    xyt <- ensure_xyt_block(block, label = paste0(product_row$alert_id, "_block"), log_file = log_file)
  } else if (depth_target$kind == "summary") {
    block <- read_product_block_xyzt(
      nc,
      structure,
      product_row,
      time_indices = needed_idx,
      depth_indices = depth_target$depth_indices,
      log_file = log_file
    )
    xyt <- summarize_block_over_depth_to_xyt(block, depth_target$summary_fun, log_file = log_file)
  } else {
    stop_with_log(
      paste0("Depth target kind not supported for block processing: ", depth_target$kind),
      log_file
    )
  }

  current_rel_idx <- match(time_idx, needed_idx)
  sev <- compute_window_severity_from_xyt(
    xyt,
    current_rel_idx = current_rel_idx,
    product_row = product_row,
    change_steps = change_steps,
    log_file = log_file
  )

  rm(block, xyt)
  invisible(gc(verbose = FALSE))

  sev
}

process_window_for_target <- function(
  nc,
  structure,
  product_row,
  depth_target,
  window_row,
  settings,
  coast_sf = NULL,
  output_dir,
  log_file = NULL
) {
  safe_dir_create(output_dir)

  sev <- compute_window_severity_matrices(
    nc = nc,
    structure = structure,
    product_row = product_row,
    depth_target = depth_target,
    window_row = window_row,
    settings = settings,
    log_file = log_file
  )
  if (is.null(sev)) return(NULL)

  base_name <- build_output_basename(product_row, depth_target$label, window_row$window_label)
  threshold_file <- file.path(output_dir, paste0(base_name, "__threshold_severity.tif"))
  change_file <- file.path(output_dir, paste0(base_name, "__change_severity.tif"))
  combined_file <- file.path(output_dir, paste0(base_name, "__combined_severity.tif"))
  max_abs_change_file <- file.path(output_dir, paste0(base_name, "__max_abs_change.tif"))

  threshold_r <- matrix_xy_to_raster(
    sev$threshold,
    grid = structure$grid,
    allow_curvilinear_resample = settings$allow_curvilinear_resample
  )
  change_r <- matrix_xy_to_raster(
    sev$change,
    grid = structure$grid,
    allow_curvilinear_resample = settings$allow_curvilinear_resample
  )
  combined_r <- matrix_xy_to_raster(
    sev$combined,
    grid = structure$grid,
    allow_curvilinear_resample = settings$allow_curvilinear_resample
  )
  max_abs_change_r <- matrix_xy_to_raster(
    sev$max_abs_change,
    grid = structure$grid,
    allow_curvilinear_resample = settings$allow_curvilinear_resample
  )

  names(threshold_r) <- paste0(base_name, "__threshold")
  names(change_r) <- paste0(base_name, "__change")
  names(combined_r) <- paste0(base_name, "__combined")
  names(max_abs_change_r) <- paste0(base_name, "__max_abs_change")

  write_geotiff(threshold_r, threshold_file, overwrite = settings$overwrite, gdal_opts = settings$gdal_write_options)
  write_geotiff(change_r, change_file, overwrite = settings$overwrite, gdal_opts = settings$gdal_write_options)
  write_geotiff(combined_r, combined_file, overwrite = settings$overwrite, gdal_opts = settings$gdal_write_options)
  write_geotiff(max_abs_change_r, max_abs_change_file, overwrite = settings$overwrite, gdal_opts = settings$gdal_write_options)

  png_file <- NA_character_
  if (isTRUE(settings$write_png) && isTRUE(product_row$write_png)) {
    png_file <- file.path(output_dir, paste0(base_name, "__debug.png"))
    title_lines <- c(
      paste0("Variable: ", product_row$product_var, " (", product_row$alert_id, ")"),
      paste0("Depth handling: ", depth_target$label),
      paste0(
        "Forecast window: ",
        format(window_row$start_time, "%Y-%m-%d %H:%M UTC", tz = "UTC"),
        " to ",
        format(window_row$end_time, "%Y-%m-%d %H:%M UTC", tz = "UTC")
      )
    )
    plot_debug_png(
      combined_r,
      png_file,
      title_lines = title_lines,
      coast_sf = coast_sf,
      threshold_r = threshold_r,
      change_r = change_r,
      product_row = product_row,
      structure = structure
    )
  }

  pw4b_files <- write_pw4b_outputs(
    matrices = list(
      threshold_severity = sev$threshold,
      change_severity = sev$change,
      combined_severity = sev$combined,
      max_abs_change = sev$max_abs_change
    ),
    structure = structure,
    settings = settings,
    base_name = base_name,
    output_dir = output_dir,
    depth_label = depth_target$label,
    log_file = log_file
  )

  rm(sev, threshold_r, change_r, combined_r, max_abs_change_r)
  invisible(gc(verbose = FALSE))

  data.table::data.table(
    alert_id = product_row$alert_id,
    product_var = product_row$product_var,
    depth_label = depth_target$label,
    window_label = window_row$window_label,
    window_start_utc = format(window_row$start_time, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    window_end_utc = format(window_row$end_time, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    threshold_file = threshold_file,
    change_file = change_file,
    combined_file = combined_file,
    max_abs_change_file = max_abs_change_file,
    png_file = png_file,
    pw4b_threshold_file = pw4b_files$pw4b_threshold_file,
    pw4b_change_file = pw4b_files$pw4b_change_file,
    pw4b_combined_file = pw4b_files$pw4b_combined_file,
    pw4b_max_abs_change_file = pw4b_files$pw4b_max_abs_change_file,
    pw4b_combined_csv_file = pw4b_files$pw4b_combined_csv_file,
    pw4b_lookup_file = pw4b_files$pw4b_lookup_file
  )
}

process_max_severity_target <- function(
  nc,
  structure,
  product_row,
  depth_target,
  window_row,
  settings,
  coast_sf = NULL,
  output_dir,
  log_file = NULL
) {
  safe_dir_create(output_dir)

  time_idx <- unlist(window_row$time_indices)
  if (length(time_idx) == 0L) {
    return(NULL)
  }

  change_steps <- change_steps_from_days(product_row$change_days, steps_per_day_from_structure(structure))
  needed_idx <- seq(max(1L, min(time_idx) - change_steps), max(time_idx))

  depth_indices_all <- sort(unique(as.integer(depth_target$depth_indices)))
  depth_indices_all <- depth_indices_all[!is.na(depth_indices_all)]
  if (length(depth_indices_all) == 0L) {
    return(NULL)
  }

  source_factor <- if (product_row$source_type == "vector_magnitude") 2L else 1L

  split_depth_batches <- function(depth_indices, max_contiguous_span) {
    depth_indices <- sort(unique(as.integer(depth_indices)))
    depth_indices <- depth_indices[!is.na(depth_indices)]
    if (length(depth_indices) == 0L) {
      return(list(integer(0)))
    }

    run_breaks <- cumsum(c(1L, as.integer(diff(depth_indices) != 1L)))
    runs <- unname(split(depth_indices, run_breaks))

    batches <- list()
    for (run in runs) {
      if (length(run) <= max_contiguous_span) {
        batches[[length(batches) + 1L]] <- run
      } else {
        starts <- seq(1L, length(run), by = max_contiguous_span)
        for (st in starts) {
          en <- min(length(run), st + max_contiguous_span - 1L)
          batches[[length(batches) + 1L]] <- run[st:en]
        }
      }
    }

    batches
  }

  full_depth_block_count <- as.integer(max(depth_indices_all) - min(depth_indices_all) + 1L)
  full_peak_mb <- estimate_peak_window_memory_mb(
    structure$grid$nx,
    structure$grid$ny,
    length(needed_idx) * full_depth_block_count * source_factor
  )

  if (is.finite(full_peak_mb) && full_peak_mb > settings$max_window_memory_mb) {
    bytes_per_xy <- as.numeric(structure$grid$nx) * as.numeric(structure$grid$ny) * 8
    slice_budget <- floor((settings$max_window_memory_mb * 1024^2 / bytes_per_xy) - 5)
    max_depth_span <- floor(slice_budget / (length(needed_idx) * source_factor))
    max_depth_span <- max(1L, as.integer(floor(max_depth_span * 0.9)))

    if (!is.finite(max_depth_span) || max_depth_span < 1L) {
      stop_with_log(
        paste0(
          "Unable to split max_severity processing into a memory-safe depth batch for alert_id='",
          product_row$alert_id,
          "'."
        ),
        log_file
      )
    }

    log_message(
      "WARN",
      paste0(
        "max_severity processing for alert_id='", product_row$alert_id,
        "' would exceed the memory guard as a single full-depth read (estimated ",
        sprintf("%.1f", full_peak_mb),
        " MB). Processing in depth batches of up to ",
        max_depth_span,
        " contiguous levels instead."
      ),
      log_file
    )

    depth_batches <- split_depth_batches(depth_indices_all, max_contiguous_span = max_depth_span)
  } else {
    depth_batches <- list(depth_indices_all)
  }

  current_rel_idx <- match(time_idx, needed_idx)
  agg_threshold <- NULL
  agg_change <- NULL
  agg_max_abs_change <- NULL

  for (depth_subset in depth_batches) {
    if (length(depth_subset) == 0L) next

    depth_block_count <- as.integer(max(depth_subset) - min(depth_subset) + 1L)
    guard_memory(
      structure,
      length(needed_idx) * depth_block_count * source_factor,
      settings,
      log_file = log_file
    )

    block_xyzt <- read_product_block_xyzt(
      nc,
      structure,
      product_row,
      time_indices = needed_idx,
      depth_indices = depth_subset,
      log_file = log_file
    )

    block_dims <- dim(block_xyzt)
    if (is.null(block_dims) || length(block_dims) != 4L) {
      stop_with_log(
        paste0(
          "max_severity processing expected an x/y/z/t block for alert_id='",
          product_row$alert_id,
          "', but received dimensions: ",
          if (is.null(block_dims)) "<none>" else paste(block_dims, collapse = " x "),
          "."
        ),
        log_file
      )
    }

    for (z_pos in seq_len(block_dims[3])) {
      xyt <- array(
        block_xyzt[, , z_pos, , drop = FALSE],
        dim = c(block_dims[1], block_dims[2], block_dims[4])
      )

      sev <- compute_window_severity_from_xyt(
        xyt,
        current_rel_idx = current_rel_idx,
        product_row = product_row,
        change_steps = change_steps,
        log_file = log_file
      )

      agg_threshold <- safe_pmax(agg_threshold, sev$threshold)
      agg_change <- safe_pmax(agg_change, sev$change)
      agg_max_abs_change <- safe_pmax(agg_max_abs_change, sev$max_abs_change)
    }

    rm(block_xyzt)
    invisible(gc(verbose = FALSE))
  }

  if (is.null(agg_threshold) || is.null(agg_change) || is.null(agg_max_abs_change)) {
    return(NULL)
  }

  combined_xy <- safe_pmax(agg_threshold, agg_change)

  threshold_r <- matrix_xy_to_raster(
    agg_threshold,
    grid = structure$grid,
    allow_curvilinear_resample = settings$allow_curvilinear_resample
  )
  change_r <- matrix_xy_to_raster(
    agg_change,
    grid = structure$grid,
    allow_curvilinear_resample = settings$allow_curvilinear_resample
  )
  combined_r <- matrix_xy_to_raster(
    combined_xy,
    grid = structure$grid,
    allow_curvilinear_resample = settings$allow_curvilinear_resample
  )
  max_abs_change_r <- matrix_xy_to_raster(
    agg_max_abs_change,
    grid = structure$grid,
    allow_curvilinear_resample = settings$allow_curvilinear_resample
  )

  base_name <- build_output_basename(product_row, depth_target$label, window_row$window_label)
  threshold_file <- file.path(output_dir, paste0(base_name, "__threshold_severity.tif"))
  change_file <- file.path(output_dir, paste0(base_name, "__change_severity.tif"))
  combined_file <- file.path(output_dir, paste0(base_name, "__combined_severity.tif"))
  max_abs_change_file <- file.path(output_dir, paste0(base_name, "__max_abs_change.tif"))

  names(threshold_r) <- paste0(base_name, "__threshold")
  names(change_r) <- paste0(base_name, "__change")
  names(combined_r) <- paste0(base_name, "__combined")
  names(max_abs_change_r) <- paste0(base_name, "__max_abs_change")

  write_geotiff(threshold_r, threshold_file, overwrite = settings$overwrite, gdal_opts = settings$gdal_write_options)
  write_geotiff(change_r, change_file, overwrite = settings$overwrite, gdal_opts = settings$gdal_write_options)
  write_geotiff(combined_r, combined_file, overwrite = settings$overwrite, gdal_opts = settings$gdal_write_options)
  write_geotiff(max_abs_change_r, max_abs_change_file, overwrite = settings$overwrite, gdal_opts = settings$gdal_write_options)

  png_file <- NA_character_
  if (isTRUE(settings$write_png) && isTRUE(product_row$write_png)) {
    png_file <- file.path(output_dir, paste0(base_name, "__debug.png"))
    title_lines <- c(
      paste0("Variable: ", product_row$product_var, " (", product_row$alert_id, ")"),
      "Depth handling: max_severity across selected/full water column",
      paste0(
        "Forecast window: ",
        format(window_row$start_time, "%Y-%m-%d %H:%M UTC", tz = "UTC"),
        " to ",
        format(window_row$end_time, "%Y-%m-%d %H:%M UTC", tz = "UTC")
      )
    )
    plot_debug_png(
      combined_r,
      png_file,
      title_lines = title_lines,
      coast_sf = coast_sf,
      threshold_r = threshold_r,
      change_r = change_r,
      product_row = product_row,
      structure = structure
    )
  }

  pw4b_files <- write_pw4b_outputs(
    matrices = list(
      threshold_severity = agg_threshold,
      change_severity = agg_change,
      combined_severity = combined_xy,
      max_abs_change = agg_max_abs_change
    ),
    structure = structure,
    settings = settings,
    base_name = base_name,
    output_dir = output_dir,
    depth_label = depth_target$label,
    log_file = log_file
  )

  rm(agg_threshold, agg_change, agg_max_abs_change, combined_xy, threshold_r, change_r, combined_r, max_abs_change_r)
  invisible(gc(verbose = FALSE))

  data.table::data.table(
    alert_id = product_row$alert_id,
    product_var = product_row$product_var,
    depth_label = depth_target$label,
    window_label = window_row$window_label,
    window_start_utc = format(window_row$start_time, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    window_end_utc = format(window_row$end_time, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    threshold_file = threshold_file,
    change_file = change_file,
    combined_file = combined_file,
    max_abs_change_file = max_abs_change_file,
    png_file = png_file,
    pw4b_threshold_file = pw4b_files$pw4b_threshold_file,
    pw4b_change_file = pw4b_files$pw4b_change_file,
    pw4b_combined_file = pw4b_files$pw4b_combined_file,
    pw4b_max_abs_change_file = pw4b_files$pw4b_max_abs_change_file,
    pw4b_combined_csv_file = pw4b_files$pw4b_combined_csv_file,
    pw4b_lookup_file = pw4b_files$pw4b_lookup_file
  )
}

append_xyt_history <- function(history_xyt, current_xyt, log_file = NULL) {
  current_xyt <- ensure_xyt_block(current_xyt, label = "current_xyt", log_file = log_file)
  current_dims <- dim(current_xyt)

  if (is.null(history_xyt)) {
    return(current_xyt)
  }

  history_xyt <- ensure_xyt_block(history_xyt, label = "history_xyt", log_file = log_file)
  history_dims <- dim(history_xyt)

  if (!identical(history_dims[1:2], current_dims[1:2])) {
    stop_with_log(
      paste0(
        "Cannot append historical and current x/y/t blocks because their grids differ. History dims: ",
        paste(history_dims, collapse = " x "),
        "; current dims: ", paste(current_dims, collapse = " x "), "."
      ),
      log_file
    )
  }

  out <- array(NA_real_, dim = c(current_dims[1], current_dims[2], history_dims[3] + current_dims[3]))
  out[, , seq_len(history_dims[3])] <- history_xyt
  out[, , history_dims[3] + seq_len(current_dims[3])] <- current_xyt
  storage.mode(out) <- "double"
  out
}

trim_xyt_history <- function(xyt, max_slices, log_file = NULL) {
  if (is.null(xyt) || is.na(max_slices) || max_slices <= 0L) {
    return(NULL)
  }

  xyt <- ensure_xyt_block(xyt, label = "history_candidate", log_file = log_file)
  xyt_dims <- dim(xyt)
  nt <- xyt_dims[3]
  keep_n <- min(as.integer(max_slices), nt)
  keep_idx <- seq.int(nt - keep_n + 1L, nt)
  out <- array(xyt[, , keep_idx, drop = FALSE], dim = c(xyt_dims[1], xyt_dims[2], keep_n))
  storage.mode(out) <- "double"
  out
}

read_depth_target_block_xyt <- function(nc, structure, product_row, depth_target, time_indices, log_file = NULL) {
  if (depth_target$kind == "single") {
    depth_indices_for_read <- NULL
    if (!is.na(depth_target$depth_index)) {
      depth_indices_for_read <- depth_target$depth_index
    }

    block <- read_product_block_xyzt(
      nc,
      structure,
      product_row,
      time_indices = time_indices,
      depth_indices = depth_indices_for_read,
      log_file = log_file
    )
    return(ensure_xyt_block(block, label = paste0(product_row$alert_id, "_block"), log_file = log_file))
  }

  if (depth_target$kind == "summary") {
    block <- read_product_block_xyzt(
      nc,
      structure,
      product_row,
      time_indices = time_indices,
      depth_indices = depth_target$depth_indices,
      log_file = log_file
    )
    return(summarize_block_over_depth_to_xyt(block, depth_target$summary_fun, log_file = log_file))
  }

  stop_with_log(
    paste0("Depth target kind not supported for combined processing: ", depth_target$kind),
    log_file
  )
}

compute_file_severity_from_current_history <- function(current_xyt, history_xyt, product_row, log_file = NULL) {
  current_xyt <- ensure_xyt_block(current_xyt, label = "current_xyt", log_file = log_file)
  current_dims <- dim(current_xyt)
  nx <- current_dims[1]
  ny <- current_dims[2]
  nt <- current_dims[3]

  low_threshold <- scalar_num(product_row$low_threshold)
  high_threshold <- scalar_num(product_row$high_threshold)
  change_threshold <- scalar_num(product_row$change_threshold)
  change_direction <- scalar_chr(product_row$change_direction, default = "abs")

  threshold_max <- NULL
  change_max <- NULL
  max_abs_change <- NULL

  for (rel_t in seq_len(nt)) {
    slice_xy <- current_xyt[, , rel_t, drop = TRUE]
    thr_xy <- threshold_severity_xy(
      slice_xy,
      low_threshold = low_threshold,
      high_threshold = high_threshold
    )
    threshold_max <- safe_pmax(threshold_max, thr_xy)
  }

  if (!is.null(history_xyt)) {
    history_xyt <- ensure_xyt_block(history_xyt, label = "history_xyt", log_file = log_file)
    history_dims <- dim(history_xyt)

    if (!identical(history_dims[1:2], current_dims[1:2])) {
      stop_with_log(
        paste0(
          "Cannot compute change severity because historical and current grids differ. History dims: ",
          paste(history_dims, collapse = " x "),
          "; current dims: ", paste(current_dims, collapse = " x "), "."
        ),
        log_file
      )
    }

    raw_prev_min <- NULL
    raw_prev_max <- NULL
    for (pj in seq_len(history_dims[3])) {
      prev_xy <- history_xyt[, , pj, drop = TRUE]
      raw_prev_min <- safe_pmin(raw_prev_min, prev_xy)
      raw_prev_max <- safe_pmax(raw_prev_max, prev_xy)
    }

    prev_min <- NULL
    prev_max <- NULL
    if (change_direction %in% c("increase", "abs")) {
      prev_min <- raw_prev_min
    }
    if (change_direction %in% c("decrease", "abs")) {
      prev_max <- raw_prev_max
    }

    for (rel_t in seq_len(nt)) {
      slice_xy <- current_xyt[, , rel_t, drop = TRUE]

      raw_inc <- if (!is.null(raw_prev_min)) {
        absolute_change_xy(
          new_xy = slice_xy,
          old_xy = raw_prev_min
        )
      } else {
        NULL
      }
      raw_dec <- if (!is.null(raw_prev_max)) {
        absolute_change_xy(
          new_xy = slice_xy,
          old_xy = raw_prev_max
        )
      } else {
        NULL
      }
      max_abs_change <- safe_pmax(max_abs_change, safe_pmax(raw_inc, raw_dec))

      if (!is.na(change_threshold)) {
        ch_xy <- NULL
        if (change_direction == "increase" && !is.null(prev_min)) {
          ch_xy <- change_severity_xy(
            new_xy = slice_xy,
            old_xy = prev_min,
            change_threshold = change_threshold,
            direction = "increase"
          )
        } else if (change_direction == "decrease" && !is.null(prev_max)) {
          ch_xy <- change_severity_xy(
            new_xy = slice_xy,
            old_xy = prev_max,
            change_threshold = change_threshold,
            direction = "decrease"
          )
        } else if (change_direction == "abs") {
          ch_inc <- if (!is.null(prev_min)) {
            change_severity_xy(
              new_xy = slice_xy,
              old_xy = prev_min,
              change_threshold = change_threshold,
              direction = "increase"
            )
          } else {
            NULL
          }
          ch_dec <- if (!is.null(prev_max)) {
            change_severity_xy(
              new_xy = slice_xy,
              old_xy = prev_max,
              change_threshold = change_threshold,
              direction = "decrease"
            )
          } else {
            NULL
          }
          ch_xy <- safe_pmax(ch_inc, ch_dec)
        }

        change_max <- safe_pmax(change_max, ch_xy)
      }
    }
  }

  if (is.null(threshold_max)) {
    threshold_max <- matrix(0, nrow = nx, ncol = ny)
  }
  if (is.null(change_max)) {
    change_max <- matrix(0, nrow = nx, ncol = ny)
  }
  if (is.null(max_abs_change)) {
    max_abs_change <- matrix(0, nrow = nx, ncol = ny)
  }

  list(
    threshold = threshold_max,
    change = change_max,
    combined = safe_pmax(threshold_max, change_max),
    max_abs_change = max_abs_change
  )
}

build_combined_window_label <- function(settings) {
  year_filter <- split_year_values(settings$nc_years)
  if (length(year_filter) > 0L) {
    return(paste0("combined_years_", paste(year_filter, collapse = "_")))
  }
  "combined_selected_nc_files"
}

format_combined_period_line <- function(period_start, period_end, processed_files) {
  if (!is.null(period_start) && !is.null(period_end)) {
    return(paste0(
      "Combined selected NetCDF files: ",
      format(period_start, "%Y-%m-%d %H:%M UTC", tz = "UTC"),
      " to ",
      format(period_end, "%Y-%m-%d %H:%M UTC", tz = "UTC"),
      " (", length(processed_files), " file(s))"
    ))
  }

  paste0("Combined selected NetCDF files: ", length(processed_files), " file(s)")
}

process_combined_regular_target <- function(
  nc_files,
  structure,
  product_row,
  depth_target,
  settings,
  coast_sf = NULL,
  output_dir,
  log_file = NULL
) {
  safe_dir_create(output_dir)

  history_xyt <- NULL
  agg_threshold <- NULL
  agg_change <- NULL
  agg_max_abs_change <- NULL
  period_start <- NULL
  period_end <- NULL
  processed_files <- character(0)
  combined_label <- build_combined_window_label(settings)

  for (nc_path in nc_files) {
    current_nc_path <- normalizePath(nc_path, winslash = "/", mustWork = FALSE)
    current_structure <- tryCatch(
      refresh_structure_time_from_nc(structure, current_nc_path, log_file = log_file),
      error = function(e) {
        log_message(
          "WARN",
          paste0(
            "Skipping NetCDF file because its structure/time metadata could not be read: ",
            current_nc_path,
            ". Reason: ", conditionMessage(e)
          ),
          log_file
        )
        NULL
      }
    )
    if (is.null(current_structure)) {
      next
    }
    if (!isTRUE(quick_validate_against_cache(current_nc_path, current_structure, log_file = log_file, stop_on_error = FALSE))) {
      next
    }

    time_indices <- seq_along(current_structure$time$values)
    if (length(time_indices) == 0L) {
      log_message("WARN", paste0("Skipping NetCDF file with no readable time axis: ", current_nc_path), log_file)
      next
    }

    steps_per_day <- steps_per_day_from_structure(current_structure)
    change_steps <- change_steps_from_days(product_row$change_days, steps_per_day)
    history_xyt <- trim_xyt_history(history_xyt, change_steps, log_file = log_file)

    depth_indices_for_read <- NULL
    if (depth_target$kind == "single" && !is.na(depth_target$depth_index)) {
      depth_indices_for_read <- depth_target$depth_index
    } else if (depth_target$kind == "summary") {
      depth_indices_for_read <- depth_target$depth_indices
    }

    depth_block_count <- if (!is.null(depth_indices_for_read) && length(depth_indices_for_read) > 0L) {
      as.integer(max(depth_indices_for_read) - min(depth_indices_for_read) + 1L)
    } else {
      1L
    }
    source_factor <- if (product_row$source_type == "vector_magnitude") 2L else 1L
    history_slices <- if (is.null(history_xyt)) 0L else dim(history_xyt)[3]
    guard_memory(
      current_structure,
      history_slices + length(time_indices) * depth_block_count * source_factor,
      settings,
      log_file = log_file
    )

    log_message("INFO", paste0("Adding NetCDF file to combined map: ", current_nc_path), log_file)

    nc <- ncdf4::nc_open(current_nc_path, readunlim = FALSE)
    current_xyt <- tryCatch(
      {
        read_depth_target_block_xyt(
          nc = nc,
          structure = current_structure,
          product_row = product_row,
          depth_target = depth_target,
          time_indices = time_indices,
          log_file = log_file
        )
      },
      finally = {
        ncdf4::nc_close(nc)
      }
    )

    sev <- compute_file_severity_from_current_history(
      current_xyt = current_xyt,
      history_xyt = history_xyt,
      product_row = product_row,
      log_file = log_file
    )

    agg_threshold <- safe_pmax(agg_threshold, sev$threshold)
    agg_change <- safe_pmax(agg_change, sev$change)
    agg_max_abs_change <- safe_pmax(agg_max_abs_change, sev$max_abs_change)

    if (change_steps > 0L) {
      history_xyt <- trim_xyt_history(
        append_xyt_history(history_xyt, current_xyt, log_file = log_file),
        max_slices = change_steps,
        log_file = log_file
      )
    } else {
      history_xyt <- NULL
    }

    file_times <- current_structure$time$values
    if (length(file_times) > 0L) {
      file_start <- min(as.POSIXct(file_times, tz = "UTC"))
      file_end <- max(as.POSIXct(file_times, tz = "UTC"))
      period_start <- if (is.null(period_start)) file_start else min(period_start, file_start)
      period_end <- if (is.null(period_end)) file_end else max(period_end, file_end)
    }
    processed_files <- c(processed_files, current_nc_path)

    rm(current_xyt, sev)
    invisible(gc(verbose = FALSE))
  }

  if (is.null(agg_threshold) || is.null(agg_change) || is.null(agg_max_abs_change)) {
    return(NULL)
  }

  combined_xy <- safe_pmax(agg_threshold, agg_change)

  threshold_r <- matrix_xy_to_raster(
    agg_threshold,
    grid = structure$grid,
    allow_curvilinear_resample = settings$allow_curvilinear_resample
  )
  change_r <- matrix_xy_to_raster(
    agg_change,
    grid = structure$grid,
    allow_curvilinear_resample = settings$allow_curvilinear_resample
  )
  combined_r <- matrix_xy_to_raster(
    combined_xy,
    grid = structure$grid,
    allow_curvilinear_resample = settings$allow_curvilinear_resample
  )
  max_abs_change_r <- matrix_xy_to_raster(
    agg_max_abs_change,
    grid = structure$grid,
    allow_curvilinear_resample = settings$allow_curvilinear_resample
  )

  base_name <- build_output_basename(product_row, depth_target$label, combined_label)
  threshold_file <- file.path(output_dir, paste0(base_name, "__threshold_severity.tif"))
  change_file <- file.path(output_dir, paste0(base_name, "__change_severity.tif"))
  combined_file <- file.path(output_dir, paste0(base_name, "__combined_severity.tif"))
  max_abs_change_file <- file.path(output_dir, paste0(base_name, "__max_abs_change.tif"))

  names(threshold_r) <- paste0(base_name, "__threshold")
  names(change_r) <- paste0(base_name, "__change")
  names(combined_r) <- paste0(base_name, "__combined")
  names(max_abs_change_r) <- paste0(base_name, "__max_abs_change")

  write_geotiff(threshold_r, threshold_file, overwrite = settings$overwrite, gdal_opts = settings$gdal_write_options)
  write_geotiff(change_r, change_file, overwrite = settings$overwrite, gdal_opts = settings$gdal_write_options)
  write_geotiff(combined_r, combined_file, overwrite = settings$overwrite, gdal_opts = settings$gdal_write_options)
  write_geotiff(max_abs_change_r, max_abs_change_file, overwrite = settings$overwrite, gdal_opts = settings$gdal_write_options)

  png_file <- NA_character_
  if (isTRUE(settings$write_png) && isTRUE(product_row$write_png)) {
    png_file <- file.path(output_dir, paste0(base_name, "__debug.png"))
    title_lines <- c(
      paste0("Variable: ", product_row$product_var, " (", product_row$alert_id, ")"),
      paste0("Depth handling: ", depth_target$label),
      format_combined_period_line(period_start, period_end, processed_files)
    )
    plot_debug_png(
      combined_r,
      png_file,
      title_lines = title_lines,
      coast_sf = coast_sf,
      threshold_r = threshold_r,
      change_r = change_r,
      product_row = product_row,
      structure = structure
    )
  }

  pw4b_files <- write_pw4b_outputs(
    matrices = list(
      threshold_severity = agg_threshold,
      change_severity = agg_change,
      combined_severity = combined_xy,
      max_abs_change = agg_max_abs_change
    ),
    structure = structure,
    settings = settings,
    base_name = base_name,
    output_dir = output_dir,
    depth_label = depth_target$label,
    log_file = log_file
  )

  out <- data.table::data.table(
    alert_id = product_row$alert_id,
    product_var = product_row$product_var,
    depth_label = depth_target$label,
    window_label = combined_label,
    window_start_utc = if (!is.null(period_start)) format(period_start, "%Y-%m-%d %H:%M:%S", tz = "UTC") else "",
    window_end_utc = if (!is.null(period_end)) format(period_end, "%Y-%m-%d %H:%M:%S", tz = "UTC") else "",
    threshold_file = threshold_file,
    change_file = change_file,
    combined_file = combined_file,
    max_abs_change_file = max_abs_change_file,
    png_file = png_file,
    pw4b_threshold_file = pw4b_files$pw4b_threshold_file,
    pw4b_change_file = pw4b_files$pw4b_change_file,
    pw4b_combined_file = pw4b_files$pw4b_combined_file,
    pw4b_max_abs_change_file = pw4b_files$pw4b_max_abs_change_file,
    pw4b_combined_csv_file = pw4b_files$pw4b_combined_csv_file,
    pw4b_lookup_file = pw4b_files$pw4b_lookup_file,
    input_nc_count = length(processed_files),
    input_nc_files = paste(processed_files, collapse = ";")
  )

  rm(agg_threshold, agg_change, agg_max_abs_change, combined_xy, threshold_r, change_r, combined_r, max_abs_change_r, history_xyt)
  invisible(gc(verbose = FALSE))

  out
}

process_combined_max_severity_target <- function(
  nc_files,
  structure,
  product_row,
  depth_target,
  settings,
  coast_sf = NULL,
  output_dir,
  log_file = NULL
) {
  safe_dir_create(output_dir)

  depth_indices_all <- sort(unique(as.integer(depth_target$depth_indices)))
  depth_indices_all <- depth_indices_all[!is.na(depth_indices_all)]
  if (length(depth_indices_all) == 0L) {
    return(NULL)
  }

  history_by_depth <- list()
  agg_threshold <- NULL
  agg_change <- NULL
  agg_max_abs_change <- NULL
  period_start <- NULL
  period_end <- NULL
  processed_files <- character(0)
  combined_label <- build_combined_window_label(settings)
  source_factor <- if (product_row$source_type == "vector_magnitude") 2L else 1L

  for (nc_path in nc_files) {
    current_nc_path <- normalizePath(nc_path, winslash = "/", mustWork = FALSE)
    current_structure <- tryCatch(
      refresh_structure_time_from_nc(structure, current_nc_path, log_file = log_file),
      error = function(e) {
        log_message(
          "WARN",
          paste0(
            "Skipping NetCDF file because its structure/time metadata could not be read: ",
            current_nc_path,
            ". Reason: ", conditionMessage(e)
          ),
          log_file
        )
        NULL
      }
    )
    if (is.null(current_structure)) {
      next
    }
    if (!isTRUE(quick_validate_against_cache(current_nc_path, current_structure, log_file = log_file, stop_on_error = FALSE))) {
      next
    }

    time_indices <- seq_along(current_structure$time$values)
    if (length(time_indices) == 0L) {
      log_message("WARN", paste0("Skipping NetCDF file with no readable time axis: ", current_nc_path), log_file)
      next
    }

    steps_per_day <- steps_per_day_from_structure(current_structure)
    change_steps <- change_steps_from_days(product_row$change_days, steps_per_day)

    log_message("INFO", paste0("Adding NetCDF file to combined max-severity map: ", current_nc_path), log_file)

    nc <- ncdf4::nc_open(current_nc_path, readunlim = FALSE)
    tryCatch(
      {
        for (depth_idx in depth_indices_all) {
          depth_key <- as.character(depth_idx)
          depth_history <- history_by_depth[[depth_key]]
          depth_history <- trim_xyt_history(depth_history, change_steps, log_file = log_file)
          history_slices <- if (is.null(depth_history)) 0L else dim(depth_history)[3]
          guard_memory(
            current_structure,
            history_slices + length(time_indices) * source_factor,
            settings,
            log_file = log_file
          )

          block <- read_product_block_xyzt(
            nc,
            current_structure,
            product_row,
            time_indices = time_indices,
            depth_indices = depth_idx,
            log_file = log_file
          )
          current_xyt <- ensure_xyt_block(block, label = paste0(product_row$alert_id, "_depth_", depth_idx), log_file = log_file)

          sev <- compute_file_severity_from_current_history(
            current_xyt = current_xyt,
            history_xyt = depth_history,
            product_row = product_row,
            log_file = log_file
          )

          agg_threshold <- safe_pmax(agg_threshold, sev$threshold)
          agg_change <- safe_pmax(agg_change, sev$change)
          agg_max_abs_change <- safe_pmax(agg_max_abs_change, sev$max_abs_change)

          if (change_steps > 0L) {
            history_by_depth[[depth_key]] <- trim_xyt_history(
              append_xyt_history(depth_history, current_xyt, log_file = log_file),
              max_slices = change_steps,
              log_file = log_file
            )
          } else {
            history_by_depth[[depth_key]] <- NULL
          }

          rm(block, current_xyt, sev)
          invisible(gc(verbose = FALSE))
        }
      },
      finally = {
        ncdf4::nc_close(nc)
      }
    )

    file_times <- current_structure$time$values
    if (length(file_times) > 0L) {
      file_start <- min(as.POSIXct(file_times, tz = "UTC"))
      file_end <- max(as.POSIXct(file_times, tz = "UTC"))
      period_start <- if (is.null(period_start)) file_start else min(period_start, file_start)
      period_end <- if (is.null(period_end)) file_end else max(period_end, file_end)
    }
    processed_files <- c(processed_files, current_nc_path)
  }

  if (is.null(agg_threshold) || is.null(agg_change) || is.null(agg_max_abs_change)) {
    return(NULL)
  }

  combined_xy <- safe_pmax(agg_threshold, agg_change)

  threshold_r <- matrix_xy_to_raster(
    agg_threshold,
    grid = structure$grid,
    allow_curvilinear_resample = settings$allow_curvilinear_resample
  )
  change_r <- matrix_xy_to_raster(
    agg_change,
    grid = structure$grid,
    allow_curvilinear_resample = settings$allow_curvilinear_resample
  )
  combined_r <- matrix_xy_to_raster(
    combined_xy,
    grid = structure$grid,
    allow_curvilinear_resample = settings$allow_curvilinear_resample
  )
  max_abs_change_r <- matrix_xy_to_raster(
    agg_max_abs_change,
    grid = structure$grid,
    allow_curvilinear_resample = settings$allow_curvilinear_resample
  )

  base_name <- build_output_basename(product_row, depth_target$label, combined_label)
  threshold_file <- file.path(output_dir, paste0(base_name, "__threshold_severity.tif"))
  change_file <- file.path(output_dir, paste0(base_name, "__change_severity.tif"))
  combined_file <- file.path(output_dir, paste0(base_name, "__combined_severity.tif"))
  max_abs_change_file <- file.path(output_dir, paste0(base_name, "__max_abs_change.tif"))

  names(threshold_r) <- paste0(base_name, "__threshold")
  names(change_r) <- paste0(base_name, "__change")
  names(combined_r) <- paste0(base_name, "__combined")
  names(max_abs_change_r) <- paste0(base_name, "__max_abs_change")

  write_geotiff(threshold_r, threshold_file, overwrite = settings$overwrite, gdal_opts = settings$gdal_write_options)
  write_geotiff(change_r, change_file, overwrite = settings$overwrite, gdal_opts = settings$gdal_write_options)
  write_geotiff(combined_r, combined_file, overwrite = settings$overwrite, gdal_opts = settings$gdal_write_options)
  write_geotiff(max_abs_change_r, max_abs_change_file, overwrite = settings$overwrite, gdal_opts = settings$gdal_write_options)

  png_file <- NA_character_
  if (isTRUE(settings$write_png) && isTRUE(product_row$write_png)) {
    png_file <- file.path(output_dir, paste0(base_name, "__debug.png"))
    title_lines <- c(
      paste0("Variable: ", product_row$product_var, " (", product_row$alert_id, ")"),
      "Depth handling: max_severity across selected/full water column",
      format_combined_period_line(period_start, period_end, processed_files)
    )
    plot_debug_png(
      combined_r,
      png_file,
      title_lines = title_lines,
      coast_sf = coast_sf,
      threshold_r = threshold_r,
      change_r = change_r,
      product_row = product_row,
      structure = structure
    )
  }

  pw4b_files <- write_pw4b_outputs(
    matrices = list(
      threshold_severity = agg_threshold,
      change_severity = agg_change,
      combined_severity = combined_xy,
      max_abs_change = agg_max_abs_change
    ),
    structure = structure,
    settings = settings,
    base_name = base_name,
    output_dir = output_dir,
    depth_label = depth_target$label,
    log_file = log_file
  )

  out <- data.table::data.table(
    alert_id = product_row$alert_id,
    product_var = product_row$product_var,
    depth_label = depth_target$label,
    window_label = combined_label,
    window_start_utc = if (!is.null(period_start)) format(period_start, "%Y-%m-%d %H:%M:%S", tz = "UTC") else "",
    window_end_utc = if (!is.null(period_end)) format(period_end, "%Y-%m-%d %H:%M:%S", tz = "UTC") else "",
    threshold_file = threshold_file,
    change_file = change_file,
    combined_file = combined_file,
    max_abs_change_file = max_abs_change_file,
    png_file = png_file,
    pw4b_threshold_file = pw4b_files$pw4b_threshold_file,
    pw4b_change_file = pw4b_files$pw4b_change_file,
    pw4b_combined_file = pw4b_files$pw4b_combined_file,
    pw4b_max_abs_change_file = pw4b_files$pw4b_max_abs_change_file,
    pw4b_combined_csv_file = pw4b_files$pw4b_combined_csv_file,
    pw4b_lookup_file = pw4b_files$pw4b_lookup_file,
    input_nc_count = length(processed_files),
    input_nc_files = paste(processed_files, collapse = ";")
  )

  rm(agg_threshold, agg_change, agg_max_abs_change, combined_xy, threshold_r, change_r, combined_r, max_abs_change_r, history_by_depth)
  invisible(gc(verbose = FALSE))

  out
}


make_combined_processing_state <- function(product_row, depth_target, settings) {
  list(
    product_row = data.table::copy(product_row),
    depth_target = depth_target,
    history_xyt = NULL,
    history_by_depth = list(),
    agg_threshold = NULL,
    agg_change = NULL,
    agg_max_abs_change = NULL,
    period_start = NULL,
    period_end = NULL,
    processed_files = character(0),
    combined_label = build_combined_window_label(settings)
  )
}

initialize_combined_processing_states <- function(cfg, structure, settings, log_file = NULL) {
  states <- list()

  for (i in seq_len(nrow(cfg))) {
    row <- cfg[i]

    log_message(
      "INFO",
      paste0(
        "Processing combined alert_id='", row$alert_id,
        "', product_var='", row$product_var,
        "', depth_mode='", row$depth_mode, "'."
      ),
      log_file
    )

    depth_targets <- resolve_depth_targets(structure, row, log_file = log_file)
    for (depth_target in depth_targets) {
      states[[length(states) + 1L]] <- make_combined_processing_state(
        product_row = row,
        depth_target = depth_target,
        settings = settings
      )
    }
  }

  states
}

update_combined_period <- function(state, current_structure, current_nc_path) {
  file_times <- current_structure$time$values
  if (length(file_times) > 0L) {
    file_start <- min(as.POSIXct(file_times, tz = "UTC"))
    file_end <- max(as.POSIXct(file_times, tz = "UTC"))
    state$period_start <- if (is.null(state$period_start)) file_start else min(state$period_start, file_start)
    state$period_end <- if (is.null(state$period_end)) file_end else max(state$period_end, file_end)
  }
  state$processed_files <- c(state$processed_files, current_nc_path)
  state
}

update_combined_regular_state_for_file <- function(
  state,
  nc,
  current_structure,
  current_nc_path,
  settings,
  block_cache = NULL,
  log_file = NULL
) {
  product_row <- state$product_row
  depth_target <- state$depth_target

  time_indices <- seq_along(current_structure$time$values)
  if (length(time_indices) == 0L) {
    log_message("WARN", paste0("Skipping NetCDF file with no readable time axis: ", current_nc_path), log_file)
    return(state)
  }

  steps_per_day <- steps_per_day_from_structure(current_structure)
  change_steps <- change_steps_from_days(product_row$change_days, steps_per_day)
  state$history_xyt <- trim_xyt_history(state$history_xyt, change_steps, log_file = log_file)

  depth_indices_for_read <- NULL
  if (depth_target$kind == "single" && !is.na(depth_target$depth_index)) {
    depth_indices_for_read <- depth_target$depth_index
  } else if (depth_target$kind == "summary") {
    depth_indices_for_read <- depth_target$depth_indices
  }

  depth_block_count <- if (!is.null(depth_indices_for_read) && length(depth_indices_for_read) > 0L) {
    as.integer(max(depth_indices_for_read) - min(depth_indices_for_read) + 1L)
  } else {
    1L
  }
  source_factor <- if (product_row$source_type == "vector_magnitude") 2L else 1L
  history_slices <- if (is.null(state$history_xyt)) 0L else dim(state$history_xyt)[3]
  guard_memory(
    current_structure,
    history_slices + length(time_indices) * depth_block_count * source_factor,
    settings,
    log_file = log_file
  )

  log_message("INFO", paste0("Adding NetCDF file to combined map: ", current_nc_path), log_file)

  current_xyt <- read_depth_target_block_xyt_cached(
    nc = nc,
    structure = current_structure,
    product_row = product_row,
    depth_target = depth_target,
    time_indices = time_indices,
    block_cache = block_cache,
    log_file = log_file
  )

  sev <- compute_file_severity_from_current_history(
    current_xyt = current_xyt,
    history_xyt = state$history_xyt,
    product_row = product_row,
    log_file = log_file
  )

  state$agg_threshold <- safe_pmax(state$agg_threshold, sev$threshold)
  state$agg_change <- safe_pmax(state$agg_change, sev$change)
  state$agg_max_abs_change <- safe_pmax(state$agg_max_abs_change, sev$max_abs_change)

  if (change_steps > 0L) {
    state$history_xyt <- trim_xyt_history(
      append_xyt_history(state$history_xyt, current_xyt, log_file = log_file),
      max_slices = change_steps,
      log_file = log_file
    )
  } else {
    state$history_xyt <- NULL
  }

  state <- update_combined_period(state, current_structure, current_nc_path)

  rm(current_xyt, sev)
  invisible(gc(verbose = FALSE))

  state
}

update_combined_max_severity_state_for_file <- function(
  state,
  nc,
  current_structure,
  current_nc_path,
  settings,
  block_cache = NULL,
  log_file = NULL
) {
  product_row <- state$product_row
  depth_target <- state$depth_target

  depth_indices_all <- sort(unique(as.integer(depth_target$depth_indices)))
  depth_indices_all <- depth_indices_all[!is.na(depth_indices_all)]
  if (length(depth_indices_all) == 0L) {
    return(state)
  }

  time_indices <- seq_along(current_structure$time$values)
  if (length(time_indices) == 0L) {
    log_message("WARN", paste0("Skipping NetCDF file with no readable time axis: ", current_nc_path), log_file)
    return(state)
  }

  steps_per_day <- steps_per_day_from_structure(current_structure)
  change_steps <- change_steps_from_days(product_row$change_days, steps_per_day)
  source_factor <- if (product_row$source_type == "vector_magnitude") 2L else 1L

  log_message("INFO", paste0("Adding NetCDF file to combined max-severity map: ", current_nc_path), log_file)

  for (depth_idx in depth_indices_all) {
    depth_key <- as.character(depth_idx)
    depth_history <- state$history_by_depth[[depth_key]]
    depth_history <- trim_xyt_history(depth_history, change_steps, log_file = log_file)
    history_slices <- if (is.null(depth_history)) 0L else dim(depth_history)[3]
    guard_memory(
      current_structure,
      history_slices + length(time_indices) * source_factor,
      settings,
      log_file = log_file
    )

    block <- read_product_block_xyzt_cached(
      nc,
      current_structure,
      product_row,
      time_indices = time_indices,
      depth_indices = depth_idx,
      block_cache = block_cache,
      log_file = log_file
    )
    current_xyt <- ensure_xyt_block(block, label = paste0(product_row$alert_id, "_depth_", depth_idx), log_file = log_file)

    sev <- compute_file_severity_from_current_history(
      current_xyt = current_xyt,
      history_xyt = depth_history,
      product_row = product_row,
      log_file = log_file
    )

    state$agg_threshold <- safe_pmax(state$agg_threshold, sev$threshold)
    state$agg_change <- safe_pmax(state$agg_change, sev$change)
    state$agg_max_abs_change <- safe_pmax(state$agg_max_abs_change, sev$max_abs_change)

    if (change_steps > 0L) {
      state$history_by_depth[[depth_key]] <- trim_xyt_history(
        append_xyt_history(depth_history, current_xyt, log_file = log_file),
        max_slices = change_steps,
        log_file = log_file
      )
    } else {
      state$history_by_depth[[depth_key]] <- NULL
    }

    rm(block, current_xyt, sev)
    invisible(gc(verbose = FALSE))
  }

  state <- update_combined_period(state, current_structure, current_nc_path)
  state
}

write_combined_state_outputs <- function(
  state,
  structure,
  settings,
  coast_sf = NULL,
  output_dir,
  log_file = NULL
) {
  safe_dir_create(output_dir)

  if (is.null(state$agg_threshold) || is.null(state$agg_change) || is.null(state$agg_max_abs_change)) {
    return(NULL)
  }

  product_row <- state$product_row
  depth_target <- state$depth_target
  combined_xy <- safe_pmax(state$agg_threshold, state$agg_change)

  threshold_r <- matrix_xy_to_raster(
    state$agg_threshold,
    grid = structure$grid,
    allow_curvilinear_resample = settings$allow_curvilinear_resample
  )
  change_r <- matrix_xy_to_raster(
    state$agg_change,
    grid = structure$grid,
    allow_curvilinear_resample = settings$allow_curvilinear_resample
  )
  combined_r <- matrix_xy_to_raster(
    combined_xy,
    grid = structure$grid,
    allow_curvilinear_resample = settings$allow_curvilinear_resample
  )
  max_abs_change_r <- matrix_xy_to_raster(
    state$agg_max_abs_change,
    grid = structure$grid,
    allow_curvilinear_resample = settings$allow_curvilinear_resample
  )

  base_name <- build_output_basename(product_row, depth_target$label, state$combined_label)
  threshold_file <- file.path(output_dir, paste0(base_name, "__threshold_severity.tif"))
  change_file <- file.path(output_dir, paste0(base_name, "__change_severity.tif"))
  combined_file <- file.path(output_dir, paste0(base_name, "__combined_severity.tif"))
  max_abs_change_file <- file.path(output_dir, paste0(base_name, "__max_abs_change.tif"))

  names(threshold_r) <- paste0(base_name, "__threshold")
  names(change_r) <- paste0(base_name, "__change")
  names(combined_r) <- paste0(base_name, "__combined")
  names(max_abs_change_r) <- paste0(base_name, "__max_abs_change")

  write_geotiff(threshold_r, threshold_file, overwrite = settings$overwrite, gdal_opts = settings$gdal_write_options)
  write_geotiff(change_r, change_file, overwrite = settings$overwrite, gdal_opts = settings$gdal_write_options)
  write_geotiff(combined_r, combined_file, overwrite = settings$overwrite, gdal_opts = settings$gdal_write_options)
  write_geotiff(max_abs_change_r, max_abs_change_file, overwrite = settings$overwrite, gdal_opts = settings$gdal_write_options)

  png_file <- NA_character_
  if (isTRUE(settings$write_png) && isTRUE(product_row$write_png)) {
    png_file <- file.path(output_dir, paste0(base_name, "__debug.png"))
    depth_line <- if (depth_target$kind == "max_severity") {
      "Depth handling: max_severity across selected/full water column"
    } else {
      paste0("Depth handling: ", depth_target$label)
    }
    title_lines <- c(
      paste0("Variable: ", product_row$product_var, " (", product_row$alert_id, ")"),
      depth_line,
      format_combined_period_line(state$period_start, state$period_end, state$processed_files)
    )
    plot_debug_png(
      combined_r,
      png_file,
      title_lines = title_lines,
      coast_sf = coast_sf,
      threshold_r = threshold_r,
      change_r = change_r,
      product_row = product_row,
      structure = structure
    )
  }

  pw4b_files <- write_pw4b_outputs(
    matrices = list(
      threshold_severity = state$agg_threshold,
      change_severity = state$agg_change,
      combined_severity = combined_xy,
      max_abs_change = state$agg_max_abs_change
    ),
    structure = structure,
    settings = settings,
    base_name = base_name,
    output_dir = output_dir,
    depth_label = depth_target$label,
    log_file = log_file
  )

  out <- data.table::data.table(
    alert_id = product_row$alert_id,
    product_var = product_row$product_var,
    depth_label = depth_target$label,
    window_label = state$combined_label,
    window_start_utc = if (!is.null(state$period_start)) format(state$period_start, "%Y-%m-%d %H:%M:%S", tz = "UTC") else "",
    window_end_utc = if (!is.null(state$period_end)) format(state$period_end, "%Y-%m-%d %H:%M:%S", tz = "UTC") else "",
    threshold_file = threshold_file,
    change_file = change_file,
    combined_file = combined_file,
    max_abs_change_file = max_abs_change_file,
    png_file = png_file,
    pw4b_threshold_file = pw4b_files$pw4b_threshold_file,
    pw4b_change_file = pw4b_files$pw4b_change_file,
    pw4b_combined_file = pw4b_files$pw4b_combined_file,
    pw4b_max_abs_change_file = pw4b_files$pw4b_max_abs_change_file,
    pw4b_combined_csv_file = pw4b_files$pw4b_combined_csv_file,
    pw4b_lookup_file = pw4b_files$pw4b_lookup_file,
    input_nc_count = length(state$processed_files),
    input_nc_files = paste(state$processed_files, collapse = ";")
  )

  rm(combined_xy, threshold_r, change_r, combined_r, max_abs_change_r)
  invisible(gc(verbose = FALSE))

  out
}

process_alerts_by_file <- function(
  nc_files,
  structure,
  cfg,
  settings,
  coast_sf = NULL,
  output_dir,
  log_file = NULL
) {
  safe_dir_create(output_dir)

  states <- initialize_combined_processing_states(
    cfg = cfg,
    structure = structure,
    settings = settings,
    log_file = log_file
  )

  if (length(states) == 0L) {
    return(data.table::data.table())
  }

  for (nc_path in nc_files) {
    current_nc_path <- normalizePath(nc_path, winslash = "/", mustWork = FALSE)
    nc <- tryCatch(
      ncdf4::nc_open(current_nc_path, readunlim = FALSE),
      error = function(e) e
    )
    if (inherits(nc, "error")) {
      log_message(
        "WARN",
        paste0(
          "Skipping NetCDF file because it could not be opened: ",
          current_nc_path,
          ". Reason: ", conditionMessage(nc)
        ),
        log_file
      )
      next
    }

    tryCatch(
      {
        if (!isTRUE(validate_supported_lagrangian_particle_nc(
          current_nc_path,
          log_file = log_file,
          stop_on_error = FALSE,
          nc = nc
        ))) {
          next
        }

        current_structure <- tryCatch(
          refresh_structure_time_from_nc(structure, current_nc_path, log_file = log_file, nc = nc),
          error = function(e) {
            log_message(
              "WARN",
              paste0(
                "Skipping NetCDF file because its structure/time metadata could not be read: ",
                current_nc_path,
                ". Reason: ", conditionMessage(e)
              ),
              log_file
            )
            NULL
          }
        )
        if (!is.null(current_structure) &&
            isTRUE(quick_validate_against_cache(current_nc_path, current_structure, log_file = log_file, stop_on_error = FALSE, nc = nc))) {
          block_cache <- new.env(parent = emptyenv())
          for (state_idx in seq_along(states)) {
            if (states[[state_idx]]$depth_target$kind == "max_severity") {
              states[[state_idx]] <- update_combined_max_severity_state_for_file(
                state = states[[state_idx]],
                nc = nc,
                current_structure = current_structure,
                current_nc_path = current_nc_path,
                settings = settings,
                block_cache = block_cache,
                log_file = log_file
              )
            } else {
              states[[state_idx]] <- update_combined_regular_state_for_file(
                state = states[[state_idx]],
                nc = nc,
                current_structure = current_structure,
                current_nc_path = current_nc_path,
                settings = settings,
                block_cache = block_cache,
                log_file = log_file
              )
            }
          }
        }
      },
      finally = {
        ncdf4::nc_close(nc)
      }
    )

    invisible(gc(verbose = FALSE))
  }

  manifest_parts <- lapply(
    states,
    write_combined_state_outputs,
    structure = structure,
    settings = settings,
    coast_sf = coast_sf,
    output_dir = output_dir,
    log_file = log_file
  )

  data.table::rbindlist(manifest_parts, fill = TRUE, use.names = TRUE)
}

process_alert_row <- function(
  nc_files,
  structure,
  product_row,
  settings,
  coast_sf = NULL,
  output_dir,
  log_file = NULL
) {
  log_message(
    "INFO",
    paste0(
      "Processing combined alert_id='", product_row$alert_id,
      "', product_var='", product_row$product_var,
      "', depth_mode='", product_row$depth_mode, "'."
    ),
    log_file
  )

  depth_targets <- resolve_depth_targets(structure, product_row, log_file = log_file)
  manifest_parts <- list()

  for (depth_target in depth_targets) {
    if (depth_target$kind == "max_severity") {
      manifest_parts[[length(manifest_parts) + 1L]] <- process_combined_max_severity_target(
        nc_files = nc_files,
        structure = structure,
        product_row = product_row,
        depth_target = depth_target,
        settings = settings,
        coast_sf = coast_sf,
        output_dir = output_dir,
        log_file = log_file
      )
    } else {
      manifest_parts[[length(manifest_parts) + 1L]] <- process_combined_regular_target(
        nc_files = nc_files,
        structure = structure,
        product_row = product_row,
        depth_target = depth_target,
        settings = settings,
        coast_sf = coast_sf,
        output_dir = output_dir,
        log_file = log_file
      )
    }
  }

  data.table::rbindlist(manifest_parts, fill = TRUE, use.names = TRUE)
}

# ------------------------------- Main --------------------------------

main <- function() {
  script_dir <- get_script_dir()
  cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
  settings <- merge_settings(default_settings(script_dir), cli_args)

  settings$nc_dir <- resolve_path_against_base(settings$nc_dir, script_dir)
  settings$inventory_path <- resolve_path_against_base(settings$inventory_path, script_dir)
  settings$config_path <- resolve_path_against_base(settings$config_path, script_dir)
  settings$output_dir <- resolve_path_against_base(settings$output_dir, script_dir)
  settings$cache_path <- resolve_path_against_base(settings$cache_path, script_dir)
  settings$pw4b_grid_path <- resolve_path_against_base(settings$pw4b_grid_path, script_dir)
  settings$pw4b_lookup_path <- resolve_path_against_base(settings$pw4b_lookup_path, script_dir)
  settings$terra_tempdir <- resolve_path_against_base(settings$terra_tempdir, script_dir)

  settings$nc_dir <- normalizePath(settings$nc_dir, winslash = "/", mustWork = FALSE)
  settings$inventory_path <- normalizePath(settings$inventory_path, winslash = "/", mustWork = FALSE)
  settings$config_path <- normalizePath(settings$config_path, winslash = "/", mustWork = FALSE)
  settings$output_dir <- normalizePath(settings$output_dir, winslash = "/", mustWork = FALSE)
  settings$cache_path <- normalizePath(settings$cache_path, winslash = "/", mustWork = FALSE)
  settings$pw4b_grid_path <- normalizePath(settings$pw4b_grid_path, winslash = "/", mustWork = FALSE)
  settings$pw4b_lookup_path <- normalizePath(settings$pw4b_lookup_path, winslash = "/", mustWork = FALSE)
  settings$terra_tempdir <- normalizePath(settings$terra_tempdir, winslash = "/", mustWork = FALSE)

  safe_dir_create(settings$output_dir)
  safe_dir_create(settings$terra_tempdir)

  log_file <- file.path(
    settings$output_dir,
    paste0("processing_log_", format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC"), ".txt")
  )

  settings$total_available_memory_mb <- detect_total_available_memory_mb()
  settings$max_window_memory_mb <- cap_requested_memory_guard_mb(
    settings$max_window_memory_mb,
    total_available_memory_mb = settings$total_available_memory_mb,
    max_fraction = 0.25,
    log_file = log_file
  )

  log_message("INFO", paste0("Confirmed output directory exists: ", settings$output_dir), log_file)

  warn_missing_optional_packages(write_png = settings$write_png)

  terra::terraOptions(
    tempdir = settings$terra_tempdir,
    memfrac = min(0.25, settings$max_window_memory_mb / settings$total_available_memory_mb),
    progress = 1
  )

  log_message("INFO", paste0("Script directory: ", script_dir), log_file)
  log_message("INFO", paste0("Output directory: ", settings$output_dir), log_file)
  log_message("INFO", paste0("Input NetCDF directory: ", settings$nc_dir), log_file)
  log_message("INFO", paste0("Include past NetCDF files: ", settings$include_past_files), log_file)
  log_message("INFO", paste0("Selected NetCDF year filter: ", ifelse(length(split_year_values(settings$nc_years)) == 0L, "all parseable years", paste(split_year_values(settings$nc_years), collapse = ", "))), log_file)

  nc_files <- list_current_or_future_nc_files(
    settings$nc_dir,
    include_past_files = settings$include_past_files,
    nc_years = settings$nc_years,
    nc_start_date = settings$nc_start_date,
    nc_end_date = settings$nc_end_date,
    log_file = log_file
  )
  settings$nc_path <- nc_files[1]

  structure <- NULL
  if (!file.exists(settings$cache_path) || isTRUE(settings$refresh_structure)) {
    structure <- characterize_product(
      nc_path = settings$nc_path,
      inventory_path = settings$inventory_path,
      preferred_crs = settings$preferred_crs,
      log_file = log_file
    )
    saveRDS(structure, settings$cache_path)
    write_structure_summary(structure, settings$output_dir)
    log_message("INFO", paste0("Wrote structure cache to: ", settings$cache_path), log_file)
  } else {
    structure <- readRDS(settings$cache_path)

    if (cache_needs_rebuild(structure, settings$inventory_path, settings$nc_path, log_file = log_file)) {
      structure <- characterize_product(
        nc_path = settings$nc_path,
        inventory_path = settings$inventory_path,
        preferred_crs = settings$preferred_crs,
        log_file = log_file
      )
      saveRDS(structure, settings$cache_path)
      log_message("INFO", paste0("Refreshed structure cache at: ", settings$cache_path), log_file)
    } else {
      quick_validate_against_cache(settings$nc_path, structure, log_file = log_file)
    }

    write_structure_summary(structure, settings$output_dir)
  }

  if (isTRUE(settings$characterize_only)) {
    log_message("INFO", "characterize_only=TRUE, exiting after structure setup.", log_file)
    return(invisible(NULL))
  }

  cfg <- read_alert_config(settings$config_path, log_file = log_file)

  # Basic row validation
  required_sources <- unique(c(cfg$source_var_1, cfg$source_var_2))
  required_sources <- required_sources[!is.na(required_sources) & nzchar(trimws(required_sources))]
  nc_var_names <- trimws(as.character(structure$variables$var_name))
  missing_sources <- required_sources[!tolower(required_sources) %in% tolower(nc_var_names)]
  if (length(missing_sources) > 0L) {
    stop_with_log(
      paste0(
        "The following source variables referenced in alert_config_EST05_ice_LP.csv are missing from the NetCDF product: ",
        paste(missing_sources, collapse = ", ")
      ),
      log_file
    )
  }

  no_logic_rows <- cfg[
    is.na(low_threshold) &
      is.na(high_threshold) &
      (is.na(change_threshold) | is.na(change_days))
  ]
  if (nrow(no_logic_rows) > 0L) {
    log_message(
      "WARN",
      paste0(
        "Some enabled alert definitions contain no usable alert logic and will be skipped: ",
        paste(no_logic_rows$alert_id, collapse = ", ")
      ),
      log_file
    )
    cfg <- cfg[!alert_id %in% no_logic_rows$alert_id]
  }

  if (nrow(cfg) == 0L) {
    stop_with_log("No enabled alert definitions remain after validation.", log_file)
  }

  partial_change_rows <- cfg[
    xor(is.na(change_threshold), is.na(change_days))
  ]
  if (nrow(partial_change_rows) > 0L) {
    log_message(
      "WARN",
      paste0(
        "Some alert definitions have incomplete change settings and their change alert component will be ignored: ",
        paste(partial_change_rows$alert_id, collapse = ", ")
      ),
      log_file
    )
  }

  if (isTRUE(settings$write_pw4b)) {
    settings$pw4b_lookup_path <- create_nc_to_pw4b_lookup(
      structure = structure,
      pw4b_grid_path = settings$pw4b_grid_path,
      lookup_path = settings$pw4b_lookup_path,
      pw4b_crs = settings$pw4b_crs,
      buffer_m = settings$pw4b_buffer_m,
      idw_k = settings$pw4b_idw_k,
      refresh = settings$refresh_pw4b_lookup,
      log_file = log_file
    )
    settings$refresh_pw4b_lookup <- FALSE
    settings$pw4b_lookup <- data.table::as.data.table(readRDS(settings$pw4b_lookup_path))
  }

  bbox <- if (structure$grid$type == "regular") {
    structure$grid$extent
  } else {
    structure$grid$template
  }
  coast_sf <- NULL
  if (isTRUE(settings$write_png)) {
    coast_sf <- tryCatch(
      make_coastline_sf(structure$grid$crs, bbox = bbox, log_file = log_file),
      error = function(e) {
        log_message(
          "WARN",
          paste0("PNG coastline overlay disabled for this run: ", conditionMessage(e)),
          log_file
        )
        NULL
      }
    )
  }

  manifest_list <- list()
  combined_output_dir <- file.path(settings$output_dir, build_combined_window_label(settings))
  safe_dir_create(combined_output_dir)

  log_message(
    "INFO",
    paste0(
      "Creating combined LP diagnostic alert maps from ", length(nc_files),
      " selected NetCDF file(s). Each selected NetCDF file is treated as one day for change_days."
    ),
    log_file
  )

  res <- process_alerts_by_file(
    nc_files = nc_files,
    structure = structure,
    cfg = cfg,
    settings = settings,
    coast_sf = coast_sf,
    output_dir = combined_output_dir,
    log_file = log_file
  )

  if (!is.null(res) && nrow(res) > 0L) {
    manifest_list[[length(manifest_list) + 1L]] <- res
  }
  invisible(gc(verbose = FALSE))

  manifest <- if (length(manifest_list) > 0L) {
    data.table::rbindlist(manifest_list, fill = TRUE, use.names = TRUE)
  } else {
    data.table::data.table()
  }
  if (nrow(manifest) > 0L) {
    manifest[, created_utc := format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC")]
  }

  manifest_path <- file.path(settings$output_dir, "output_manifest.csv")
  data.table::fwrite(manifest, manifest_path, na = "")

  log_message("INFO", paste0("Finished successfully. Output manifest: ", manifest_path), log_file)
  invisible(manifest)
}

should_run_main <- function() {
  if (!interactive()) {
    return(TRUE)
  }
  scalar_lgl(DEFAULT_PARAMETERS$run_main_when_sourced, default = TRUE)
}

if (should_run_main()) {
  main()
} else {
  message("Script loaded without executing main(). Set run_main_when_sourced = TRUE or call main() manually.")
}
