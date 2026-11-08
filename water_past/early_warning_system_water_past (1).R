#!/usr/bin/env Rscript
# Historical water alerts from monthly or multi-month RDS data.frames.
# Inputs: EST05 physics daily min/max/mean; ERGOM state daily; ERGOM sediment/rates.
# Edit SETTINGS below. Keep alert_config_water_past.csv beside this script.
# install.packages(c("data.table", "callr", "jsonlite", "processx", "terra"))
# Install 64-bit Python 3.9+ for streaming monthly preparation (standard library only).
# No NetCDF files, inventory workbook, forecast dates or external PW4B grid required.
#
# Memory: selected inputs are prepared on DISK sequentially under a 4 GiB Python
# guard. Modelling uses up to PARALLEL_WORKERS fresh R processes. Their vector-heap
# limits share PARALLEL_VECTOR_HEAP_GB (not a total-process RAM cap). Each saves
# its result and exits; summaries merge chronologically. Each month must fit.
# Preparation uses temporary uncompressed disk space, not R RAM. Prepared monthly
# copies are retained for restart; KEEP_PREPARED controls their removal on success.
# A forced OS kill can leave temporary preparation or .pending_* output folders.
# Only folders with a valid complete.rds receipt are considered completed.
#
# Dates: use the input Date column unchanged. ERGOM time_utc may be an interval
# end label on the next day. Never infer represented dates from that timestamp.
# change_days=N compares each day with available observations in [date-N,date-1]
# for the SAME grid cell, including previous SELECTED files. Calendar gaps are not bridged
# as if adjacent rows were adjacent days. First-day/no-history change is NA.
# Filenames must end _YYYY-MM-DD_to_YYYY-MM-DD.rds when date limits are set.
# Only files whose ranges overlap the inclusive requested interval are opened.
# Earlier excluded files are never read, including for change history. An
# overlapping file is read in full and can supply earlier rows inside that file.
# Parallel and sequential results use this SAME selected input set. Parallel
# workers may reread selected predecessors for history; no unselected file is read.
#
# Physics: threshold_stat=extrema uses daily min for low and max for high;
# mean/min/max select a single statistic. change_stat selects mean/min/max.
# ERGOM: threshold_stat=change_stat=value. No unit conversion is performed.
# Vector magnitude from daily mean components is the magnitude of the DAILY MEAN
# VECTOR, not daily mean/max instantaneous speed. Component extrema cannot recover
# speed extrema. vector_magnitude therefore accepts physics mean/mean only.
#
# Severity matches the supplied water/ice scripts: threshold exceedance divided
# by first positive abs(low,high), or 1 when both zero; change exceedance divided
# by change_threshold; combined=max(available active severities). Equality gives 0.
# Missing inputs stay NA. assessment_complete records whether all active tests
# were evaluable; alert_triggered is NA if incomplete without a known exceedance.
# Monthly and period summaries contain counts and maxima per cell, not a sum of
# dimensionless severities across variables. Each alert is reported separately.
#
# Original cellcode IDs are preserved. CSV summaries include cellcode for PW4B.
# GeoTIFFs use the exact 1 km EPSG:3035 cells encoded by 1kmE...N... IDs, after
# verifying their centres against the supplied grid. No nearest-value filling.
# Only surface/bottom exports are available; full-column/depth-index alerts fail.
# Config examples are carried over from the supplied scripts, not calibrated
# ecological thresholds. Fill/enable ERGOM rules in native units before using them.

SETTINGS <- list(
  INPUT_DIR = "Z:/Data/MSI_models/early warning system water/final_script_past/rds",                         # Folder containing the three RDS types
  GRID_RDS = "MSI_coastal_grid.rds",        # Relative to this script
  CONFIG_CSV = "alert_config_water_past.csv",
  OUTPUT_DIR = "water_past_output_new",
  CACHE_DIR = "water_past_monthly_cache",  # On a disk with ample free space
  TEMP_DIR = "",                          # Blank = CACHE_DIR; temporary decompression
  RECURSIVE = FALSE,
  FILE_PATTERNS = c(physics="^EST05_physics_daily.*[.]rds$",
                   state="^ERGOM_state_daily.*[.]rds$",
                   sediment_rates="^ERGOM_sediment_rates_daily.*[.]rds$"),
  START_DATE = "2023-01-01",                        # Optional YYYY-MM-DD, inclusive
  END_DATE = "2023-06-30",                          # Inclusive filename/row date filter
  KEEP_PREPARED = TRUE,
  PYTHON_EXE = "",                        # Optional full path to python.exe
  SPLIT_MEMORY_GB = 4,
  RESERVE_DISK_GB = 5,
  MAX_STAGING_GB = 1000,
  R_VECTOR_HEAP_GB = 16,                  # Vector heap only, not total RAM
  MAX_MONTH_FILE_GIB = 4,                 # Guard against unintended huge monthly files
  THREADS = 2L,
  PARALLEL_WORKERS = 2L,                  # 1 = original sequential execution
  PARALLEL_VECTOR_HEAP_GB = 16,           # Shared vector-heap budget: 2 workers => 8 GiB each
  MAX_CHANGE_DAYS = 366L,
  WRITE_DAILY_RDS = TRUE,
  WRITE_GEOTIFF = TRUE,                   # Requires terra
  WRITE_PNG = TRUE,                       # Respects each rule's write_png too
  MAX_RASTER_CELLS = 5000000L,
  WORKER_ARCH = "same"                    # Leave unchanged on ordinary R installations
)
RUN_MODEL <- TRUE
MODEL_VERSION <- "water-past-1.1"

script_path <- function() {
  files <- vapply(sys.frames(),function(f) if(is.null(f$ofile)) "" else as.character(f$ofile)[1],character(1))
  files <- files[nzchar(files)]
  if(length(files)) return(normalizePath(tail(files,1),winslash="/",mustWork=TRUE))
  a <- commandArgs(FALSE); a <- sub("^--file=","",a[grepl("^--file=",a)])
  if(length(a)) return(normalizePath(a[1],winslash="/",mustWork=TRUE))
  stop("Source the saved script, or run it with Rscript. Do not paste selected lines.")
}
MODEL_SCRIPT <- tryCatch(script_path(),error=function(e) "")

fail <- function(...) stop(...,call.=FALSE)
say <- function(...) message(format(Sys.time(),"%H:%M:%S"),"  ",...)
mkdir <- function(p) {dir.create(p,recursive=TRUE,showWarnings=FALSE);if(!dir.exists(p)) fail("Cannot create directory: ",p);invisible(p)}
absolute_path <- function(p,base) {
  if(!grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)",p)) p<-file.path(base,p)
  normalizePath(p,winslash="/",mustWork=FALSE)
}
object_hash <- function(x) {p<-tempfile();on.exit(unlink(p));saveRDS(x,p,compress=FALSE,version=2);unname(tools::md5sum(p))}
file_signature <- function(p) {i<-file.info(p);list(path=normalizePath(p,winslash="/",mustWork=TRUE),size=i$size,mtime=as.numeric(i$mtime))}
save_rds <- function(x,p) {con<-gzfile(p,"wb",compression=1);tryCatch(saveRDS(x,con,version=3),finally=close(con));invisible(p)}
finite_max <- function(x) {x<-x[is.finite(x)];if(length(x)) max(x) else NA_real_}
max_available <- function(a,b) {x<-pmax(a,b,na.rm=TRUE);x[!is.finite(x)]<-NA_real_;x}
min_available <- function(a,b) {x<-pmin(a,b,na.rm=TRUE);x[!is.finite(x)]<-NA_real_;x}
parse_bool <- function(x,what) {
  s<-tolower(trimws(as.character(x)));if(anyNA(s)||any(!s %in% c("true","false","1","0","yes","no"))) fail("Invalid TRUE/FALSE in ",what)
  s %in% c("true","1","yes")
}

read_config <- function(path,settings) {
  if(!file.exists(path)) fail("Missing config: ",path)
  cfg<-data.table::fread(path,colClasses="character",na.strings=NULL,showProgress=FALSE)
  required<-c("enabled","alert_id","product_var","source_type","source_var_1","source_var_2",
              "low_threshold","high_threshold","change_threshold","change_days","change_direction","depth_mode")
  if(length(setdiff(required,names(cfg)))) fail("Missing config columns: ",paste(setdiff(required,names(cfg)),collapse=", "))
  defaults<-list(dataset="physics",threshold_stat="mean",change_stat="mean",allowed_qc="0;1",
                 allow_hold="FALSE",min_valid_hours="1",write_png="TRUE",notes="",depth_indices="",summary_fun="")
  for(nm in names(defaults)) if(!nm %in% names(cfg)) cfg[,(nm):=defaults[[nm]]]
  for(nm in names(cfg)) data.table::set(cfg,j=nm,value=trimws(cfg[[nm]]))
  for(nm in c("enabled","allow_hold","write_png")) data.table::set(cfg,j=nm,value=parse_bool(cfg[[nm]],nm))
  if(anyDuplicated(cfg$alert_id)||any(!grepl("^[A-Za-z][A-Za-z0-9_-]*$",cfg$alert_id))) fail("Use unique filename-safe alert_id values.")
  cfg<-cfg[enabled==TRUE]
  if(!nrow(cfg)) fail("No enabled rules. Fill thresholds and set enabled=TRUE in the config.")
  for(nm in c("low_threshold","high_threshold","change_threshold","change_days","min_valid_hours")) {
    raw<-cfg[[nm]];v<-suppressWarnings(as.numeric(raw));blank<-raw %in% c("","NA")
    if(any(!blank & !is.finite(v))) fail("Invalid numeric value in ",nm)
    v[blank]<-NA_real_;data.table::set(cfg,j=nm,value=v)
  }
  if(any(!cfg$dataset %in% names(settings$FILE_PATTERNS))) fail("dataset must be physics, state or sediment_rates.")
  if(any(!cfg$source_type %in% c("native","vector_magnitude"))) fail("source_type must be native or vector_magnitude.")
  if(any(!cfg$depth_mode %in% c("surface","bottom"))||any(nzchar(cfg$depth_indices))||any(nzchar(cfg$summary_fun))) fail("RDS exports only have surface/bottom. Leave depth_indices and summary_fun blank.")
  if(any(!cfg$change_direction %in% c("abs","increase","decrease"))) fail("Invalid change_direction.")
  if(any(is.na(cfg$low_threshold)&is.na(cfg$high_threshold)&is.na(cfg$change_threshold))) fail("An enabled rule has no threshold: disable it or fill a threshold.")
  if(any(!is.na(cfg$low_threshold)&!is.na(cfg$high_threshold)&cfg$low_threshold>cfg$high_threshold)) fail("low_threshold exceeds high_threshold.")
  active<-!is.na(cfg$change_threshold)
  if(any(active & (is.na(cfg$change_days)|cfg$change_days<1|cfg$change_days!=floor(cfg$change_days)|cfg$change_days>settings$MAX_CHANGE_DAYS))||any(cfg$change_threshold[active]<=0)) fail("Change rules need a positive threshold and integer change_days within MAX_CHANGE_DAYS.")
  if(any(!active & !is.na(cfg$change_days))) fail("Blank change_threshold requires blank change_days.")
  if(any(is.na(cfg$min_valid_hours)|cfg$min_valid_hours<1|cfg$min_valid_hours>24|cfg$min_valid_hours!=floor(cfg$min_valid_hours))) fail("min_valid_hours must be 1..24.")
  for(i in seq_len(nrow(cfg))) {
    r<-cfg[i];phy<-r$dataset=="physics"
    if(!r$threshold_stat %in% (if(phy)c("extrema","mean","min","max") else "value")||!r$change_stat %in% (if(phy)c("mean","min","max") else "value")) fail("Invalid statistic for ",r$alert_id)
    if(r$source_type=="vector_magnitude" && (!phy||r$threshold_stat!="mean"||r$change_stat!="mean"||!nzchar(r$source_var_2))) fail("Vector rules require physics mean/mean and two components; speed extrema cannot be recovered.")
    if(!nzchar(r$source_var_1)) fail("Blank source_var_1 for ",r$alert_id)
    if(r$source_type=="native"&&nzchar(r$source_var_2)) fail("source_var_2 must be blank for native rules.")
    q<-strsplit(r$allowed_qc,";",fixed=TRUE)[[1]]
    if(!length(q)||any(!q %in% c("0","1"))) fail("allowed_qc must be 0, 1 or 0;1. Other codes indicate unavailable/rejected data.")
    if(phy && !setequal(q,c("0","1"))) fail("Physics daily summaries do not retain row-level QC flags. Use allowed_qc=0;1 and min_valid_hours; observed-only means require re-aggregating hourly input.")
  }
  as.data.frame(cfg)
}

read_grid <- function(path) {
  g<-readRDS(path)
  if(!is.data.frame(g)||!all(c("cell_id","cell_x","cell_y") %in% names(g))) fail("Grid must contain cell_id, cell_x, cell_y.")
  crs<-attr(g,"crs",exact=TRUE)
  if(is.null(crs)||!grepl("4326",crs,fixed=TRUE)) fail("Expected supplied EPSG:4326 coastal grid.")
  g<-data.frame(cell_id=as.character(g$cell_id),cell_x=g$cell_x,cell_y=g$cell_y)
  if(anyNA(g)||anyDuplicated(g$cell_id)||any(!nzchar(g$cell_id))||any(!is.finite(g$cell_x)|!is.finite(g$cell_y))) fail("Invalid/duplicate coastal grid cells.")
  attr(g,"crs")<-"EPSG:4326";g
}

resolve_rule <- function(rule,dt,meta) {
  suffix<-if(rule$depth_mode=="surface") "surf" else "bot"
  aliases<-c(icefrac="siconc",icethic="sithick",uice="siu",vice="siv")
  resolve<-function(x) {
    if(x %in% names(aliases)) x<-unname(aliases[x])
    if(grepl("_(min|max|mean|n_valid|__qc)$",x)) fail("Use base variable names and statistic columns in config, not statistic suffixes: ",x)
    b<-if(grepl("_(surf|bot)$",x)) x else paste0(x,"_",suffix)
    if(!endsWith(b,paste0("_",suffix))) fail("Variable layer disagrees with depth_mode: ",b)
    b
  }
  bases<-resolve(rule$source_var_1)
  if(rule$source_type=="vector_magnitude") bases<-c(bases,resolve(rule$source_var_2))
  if(anyDuplicated(bases)) fail("Vector components must be different columns.")
  cols<-character();units<-character();statuses<-character()
  for(b in bases) {
    need<-if(rule$dataset=="physics") paste0(b,"_",unique(c(if(rule$threshold_stat=="extrema") c("min","max") else rule$threshold_stat,rule$change_stat))) else b
    if(!all(need %in% names(dt))) fail("Missing columns for ",rule$alert_id,": ",paste(setdiff(need,names(dt)),collapse=", "))
    if(!is.data.frame(meta)||!all(c("column","units","status") %in% names(meta))) fail("Input variable_metadata is missing required column/units/status fields.")
    mr<-meta[as.character(meta$column) %in% need,,drop=FALSE]
    if(anyDuplicated(mr$column)||!all(need %in% mr$column)) fail("Incomplete or duplicate variable metadata for ",b)
    if(anyNA(mr$status)||anyNA(mr$units)||any(!nzchar(mr$units))) fail("Missing status/units for ",b)
    if(any(grepl("^hold",mr$status))&&!rule$allow_hold) fail("Variable is on HOLD in source metadata: ",b,". Review it before explicitly setting allow_hold=TRUE.")
    if(any(!grepl("^hold",mr$status)&!mr$status %in% c("accepted","present","available","ok"))) fail("Unrecognized metadata status for ",b,": ",paste(unique(mr$status),collapse=","))
    units<-c(units,unique(as.character(mr$units)));statuses<-c(statuses,unique(as.character(mr$status)))
    if(rule$dataset=="physics") {
      if(!paste0(b,"_n_valid") %in% names(dt)) fail("Missing hourly validity counts for ",b)
      cols<-c(cols,need,paste0(b,"_n_valid"))
    } else {
      if(!paste0(b,"__qc") %in% names(dt)) fail("Missing QC column for ",b)
      cols<-c(cols,need,paste0(b,"__qc"))
    }
  }
  if(length(unique(units))!=1L) fail("Components/statistics have incompatible units for ",rule$alert_id)
  list(bases=bases,columns=unique(cols),units=unique(units),status=paste(unique(statuses),collapse=";"))
}

rule_series <- function(dt,rule,mapping) {
  get_value<-function(b,stat) {
    nm<-if(rule$dataset=="physics") paste0(b,"_",stat) else b
    if(!is.numeric(dt[[nm]])||is.object(dt[[nm]])) fail("Expected ordinary numeric variable: ",nm)
    x<-as.numeric(dt[[nm]])
    if(rule$dataset=="physics") {
      n<-dt[[paste0(b,"_n_valid")]]
      if(!is.numeric(n)||any(!is.na(n)&(!is.finite(n)|n<0|n>24|n!=floor(n)))) fail("Invalid hourly counts for ",b)
      good<-!is.na(n)&n>=rule$min_valid_hours
    } else {
      q<-dt[[paste0(b,"__qc")]];allowed<-as.integer(strsplit(rule$allowed_qc,";",fixed=TRUE)[[1]])
      if(!is.numeric(q)||any(!is.na(q)&(!is.finite(q)|q!=floor(q)|!q %in% 0:5))) fail("Unexpected QC codes for ",b)
      good<-!is.na(q)&q %in% allowed
    }
    x[!is.finite(x)|!good]<-NA_real_;x
  }
  b<-mapping$bases
  if(rule$source_type=="vector_magnitude") {
    a<-get_value(b[1],"mean");v<-get_value(b[2],"mean");x<-sqrt(a*a+v*v)
    x[!is.finite(x)]<-NA_real_;return(list(low=x,high=x,change=x))
  }
  lo<-if(rule$threshold_stat=="extrema") "min" else rule$threshold_stat
  hi<-if(rule$threshold_stat=="extrema") "max" else rule$threshold_stat
  list(low=get_value(b,lo),high=get_value(b,hi),change=get_value(b,rule$change_stat))
}

threshold_values <- function(low,high,rule) {
  haslo<-!is.na(rule$low_threshold);hashi<-!is.na(rule$high_threshold)
  if(!haslo&&!hashi) return(list(severity=rep(NA_real_,length(low)),complete=rep(TRUE,length(low))))
  refs<-abs(c(rule$low_threshold,rule$high_threshold));refs<-refs[is.finite(refs)&refs>0];scale<-if(length(refs)) refs[1] else 1
  a<-if(haslo) pmax(0,(rule$low_threshold-low)/scale) else rep(NA_real_,length(low))
  b<-if(hashi) pmax(0,(high-rule$high_threshold)/scale) else rep(NA_real_,length(low))
  list(severity=max_available(a,b),complete=(!haslo|is.finite(low))&(!hashi|is.finite(high)))
}

# A small dated matrix (grid cells x at most change_days days) crosses month
# boundaries on DISK. Only one alert's history is loaded at a time.
compute_daily <- function(dt,grid_index,grid,rule,mapping,history=NULL) {
  n<-nrow(dt);ng<-nrow(grid);s<-rule_series(dt,rule,mapping)
  thr<-threshold_values(s$low,s$high,rule)
  change<-raw<-rep(NA_real_,n);pairs<-integer(n)
  change_active<-!is.na(rule$change_threshold)
  if(is.null(history)) history<-list(dates=numeric(),values=matrix(numeric(),nrow=ng,ncol=0))
  if(nrow(history$values)!=ng) fail("History grid mismatch.")
  days<-sort(unique(as.numeric(dt$date)))
  bounds<-dt[,list(first=.I[1],last=.I[.N]),by="date"]
  for(j in seq_along(days)) {
    day<-days[j];ix<-seq.int(bounds$first[j],bounds$last[j]);gi<-grid_index[ix]
    if(change_active) {
      keep<-which(history$dates>=day-rule$change_days & history$dates<day)
      history$dates<-history$dates[keep];history$values<-history$values[,keep,drop=FALSE]
      if(length(keep)) {
        prev<-history$values[gi,,drop=FALSE]
        pminv<-rep(NA_real_,length(ix));pmaxv<-pminv
        for(k in seq_len(ncol(prev))) {pminv<-min_available(pminv,prev[,k]);pmaxv<-max_available(pmaxv,prev[,k])}
        cur<-s$change[ix];ok<-is.finite(cur);pairs[ix]<-rowSums(is.finite(prev))*as.integer(ok)
        delta<-switch(rule$change_direction,increase=cur-pminv,decrease=pmaxv-cur,abs=max_available(abs(cur-pminv),abs(cur-pmaxv)))
        change[ix]<-pmax(0,(delta-rule$change_threshold)/rule$change_threshold)
        raw[ix]<-max_available(abs(cur-pminv),abs(cur-pmaxv))
        rm(prev)
      }
      v<-rep(NA_real_,ng);v[gi]<-s$change[ix]
      history$dates<-c(history$dates,day);history$values<-cbind(history$values,v)
    }
  }
  if(change_active) {
    keep<-which(history$dates>max(days)-rule$change_days)
    history$dates<-history$dates[keep];history$values<-history$values[,keep,drop=FALSE]
  }
  complete<-thr$complete & (!change_active|pairs>0)
  combined<-max_available(thr$severity,change)
  trigger<-rep(NA, n);trigger[complete]<-combined[complete]>0;trigger[!is.na(combined)&combined>0]<-TRUE
  result<-data.table::data.table(cell_id=dt$cell_id,date=dt$date,
      low_test_value=s$low,high_test_value=s$high,change_test_value=s$change,
      threshold_severity=thr$severity,change_severity=change,combined_severity=combined,
      max_abs_change=raw,history_observations=pairs,assessment_complete=complete,alert_triggered=trigger)
  list(daily=result,history=history)
}

summarize_daily <- function(daily,grid) {
  s<-daily[,list(n_days_present=.N,n_days_complete=sum(assessment_complete),
    n_days_evaluable=sum(is.finite(combined_severity)),n_days_alert=sum(alert_triggered %in% TRUE),
    n_days_unknown=sum(is.na(alert_triggered)),
    threshold_severity=finite_max(threshold_severity),change_severity=finite_max(change_severity),
    combined_severity=finite_max(combined_severity),max_abs_change=finite_max(max_abs_change),
    first_alert=if(any(alert_triggered %in% TRUE)) min(date[alert_triggered %in% TRUE]) else as.Date(NA),
    last_alert=if(any(alert_triggered %in% TRUE)) max(date[alert_triggered %in% TRUE]) else as.Date(NA)),by="cell_id"]
  out<-merge(data.table::as.data.table(grid),s,by="cell_id",all.x=TRUE,sort=FALSE)
  out<-out[match(grid$cell_id,out$cell_id)]
  for(nm in grep("^n_days_",names(out),value=TRUE)) data.table::set(out,i=which(is.na(out[[nm]])),j=nm,value=0L)
  out[,cellcode:=cell_id];data.table::setDF(out);out
}
merge_summaries <- function(old,new) {
  if(is.null(old)) return(new)
  if(!identical(old$cell_id,new$cell_id)) fail("Period summary grid mismatch.")
  for(nm in grep("^n_days_",names(new),value=TRUE)) new[[nm]]<-as.double(old[[nm]])+new[[nm]]
  for(nm in c("threshold_severity","change_severity","combined_severity","max_abs_change")) new[[nm]]<-max_available(old[[nm]],new[[nm]])
  new$first_alert<-as.Date(min_available(as.numeric(old$first_alert),as.numeric(new$first_alert)),origin="1970-01-01")
  new$last_alert<-as.Date(max_available(as.numeric(old$last_alert),as.numeric(new$last_alert)),origin="1970-01-01")
  new
}

map_template <- function(grid,settings) {
  if(!settings$WRITE_GEOTIFF&&!settings$WRITE_PNG) return(NULL)
  if(!requireNamespace("terra",quietly=TRUE)) fail("Install terra for GeoTIFF/PNG exports, or turn both switches off.")
  if(any(!grepl("^1kmE[0-9]+N[0-9]+$",grid$cell_id))) fail("Map export requires the supplied 1kmE...N... grid IDs.")
  x<-as.numeric(sub("^1kmE([0-9]+)N[0-9]+$","\\1",grid$cell_id))*1000+500
  y<-as.numeric(sub("^1kmE[0-9]+N([0-9]+)$","\\1",grid$cell_id))*1000+500
  pts<-terra::vect(data.frame(x=x,y=y),geom=c("x","y"),crs="EPSG:3035")
  ll<-terra::crds(terra::project(pts,"EPSG:4326"))
  if(any(abs(ll[,1]-grid$cell_x)>1e-5|abs(ll[,2]-grid$cell_y)>1e-5)) fail("Encoded grid centres disagree with supplied coordinates; no shifted raster was written.")
  nr<-(max(y)-min(y))/1000+1;nc<-(max(x)-min(x))/1000+1
  if(nr*nc>settings$MAX_RASTER_CELLS) fail("Map extent exceeds MAX_RASTER_CELLS.")
  r<-terra::rast(nrows=nr,ncols=nc,xmin=min(x)-500,xmax=max(x)+500,ymin=min(y)-500,ymax=max(y)+500,crs="EPSG:3035")
  cells<-terra::cellFromXY(r,cbind(x,y));if(anyDuplicated(cells)) fail("Duplicate raster cells.")
  list(raster=r,cells=cells)
}
write_maps <- function(summary,prefix,rule,settings,template,title) {
  if(is.null(template)) return(invisible(NULL))
  for(nm in c("threshold_severity","change_severity","combined_severity","max_abs_change")) {
    r<-template$raster;v<-rep(NA_real_,terra::ncell(r));v[template$cells]<-summary[[nm]];terra::values(r)<-v
    names(r)<-nm
    if(settings$WRITE_GEOTIFF) terra::writeRaster(r,paste0(prefix,"__",nm,".tif"),overwrite=FALSE,wopt=list(datatype="FLT8S",gdal=c("COMPRESS=DEFLATE","PREDICTOR=3")))
    if(nm=="combined_severity"&&settings$WRITE_PNG&&rule$write_png) {
      grDevices::png(paste0(prefix,"__combined_severity.png"),width=2400,height=1800,res=300)
      tryCatch({
        if(any(is.finite(v))) terra::plot(r,main=paste(rule$alert_id,title,sep="\n"),col=hcl.colors(60,"YlOrRd",rev=TRUE),axes=FALSE)
        else {graphics::plot.new();graphics::title(main=paste(rule$alert_id,title,sep="\n"));graphics::text(.5,.5,"No evaluable data")}
      },finally=grDevices::dev.off())
    }
    rm(r,v)
  }
  invisible(NULL)
}

month_worker <- function(job,previous,config,settings,run_id,defer_merge=FALSE) {
  suppressPackageStartupMessages(library(data.table));setDTthreads(settings$THREADS);setNumericRounding(0L)
  mem.maxVSize(settings$R_VECTOR_HEAP_GB*1024)
  if(!is.finite(mem.maxVSize())||mem.maxVSize()>settings$R_VECTOR_HEAP_GB*1024) fail("Cannot establish R vector-heap limit.")
  if(file.info(job$path)$size/1024^3>settings$MAX_MONTH_FILE_GIB) fail("Monthly input exceeds MAX_MONTH_FILE_GIB: ",job$path)
  grid<-read_grid(settings$GRID_RDS);template<-map_template(grid,settings)
  seed_dir<-""
  if(defer_merge) {
    seed_dir<-tempfile("water_selected_history_");mkdir(seed_dir)
    on.exit(unlink(seed_dir,recursive=TRUE),add=TRUE)
    seed_selected_history(job,config,settings,grid,seed_dir)
  }
  say("Reading one monthly RDS: ",basename(job$path));dt<-readRDS(job$path)
  if(!is.data.frame(dt)||!all(c("cell_id","cell_x","cell_y","date") %in% names(dt))||!inherits(dt$date,"Date")) fail("Expected daily data.frame with cell_id, cell_x, cell_y, Date.")
  if(!nrow(dt)||anyDuplicated(names(dt))||anyNA(dt$cell_id)||anyNA(dt$date)) fail("Empty data or missing/ambiguous IDs/dates.")
  if(any(!is.finite(as.numeric(dt$date))|as.numeric(dt$date)!=floor(as.numeric(dt$date)))) fail("Invalid dates.")
  dates<-sort(unique(dt$date))
  if(any(format(dates,"%Y-%m")!=job$month)||min(dates)!=as.Date(job$first_date)||max(dates)!=as.Date(job$last_date)) fail("Input dates differ from monthly preparation manifest.")
  meta<-attr(dt,"variable_metadata",exact=TRUE)
  source_meta<-list(variable_metadata=meta,time_semantics=attr(dt,"time_semantics",exact=TRUE),
                    qc_codes=attr(dt,"qc_codes",exact=TRUE),repair=attr(dt,"repair",exact=TRUE),
                    source_files=attr(dt,"source_files",exact=TRUE),export_signature=attr(dt,"export_signature",exact=TRUE))
  setDT(dt)
  if(anyDuplicated(dt,by=c("cell_id","date"))) fail("Duplicate cell/date rows. Resolve them before modelling; no silent averaging or double counting.")
  gi<-match(dt$cell_id,grid$cell_id)
  if(anyNA(gi)) fail("RDS contains cell IDs absent from the supplied grid.")
  if(any(!is.finite(dt$cell_x)|!is.finite(dt$cell_y)|abs(dt$cell_x-grid$cell_x[gi])>1e-6|abs(dt$cell_y-grid$cell_y[gi])>1e-6)) fail("RDS coordinates disagree with coastal grid.")
  maps<-lapply(seq_len(nrow(config)),function(i) resolve_rule(config[i,,drop=FALSE],dt,meta))
  needed<-unique(c("cell_id","date",unlist(lapply(maps,`[[`,"columns"))))
  # Drop unused variables and metadata by reference before creating results.
  for(nm in setdiff(names(attributes(dt)),c("names","row.names","class",".internal.selfref"))) setattr(dt,nm,NULL)
  for(nm in setdiff(names(dt),needed)) set(dt,j=nm,value=NULL)
  setorderv(dt,c("date","cell_id"));gi<-match(dt$cell_id,grid$cell_id);gc(FALSE)
  mkdir(dirname(job$output));stage<-tempfile(".pending_",tmpdir=dirname(job$output));mkdir(stage)
  on.exit(unlink(stage,recursive=TRUE),add=TRUE)
  mkdir(file.path(stage,"history"));mkdir(file.path(stage,"period"))
  export_rows<-dt$date>=settings$START_DATE & dt$date<=settings$END_DATE
  for(i in seq_len(nrow(config))) {
    rule<-config[i,,drop=FALSE];id<-rule$alert_id;say("Evaluating ",id," for ",job$month)
    hp<-if(defer_merge) file.path(seed_dir,paste0(id,".rds")) else if(is.null(previous)) "" else file.path(previous,"history",paste0(id,".rds"))
    history<-if(nzchar(hp)&&file.exists(hp)) readRDS(hp) else NULL
    got<-compute_daily(dt,gi,grid,rule,maps[[i]],history);rm(history)
    save_rds(got$history,file.path(stage,"history",paste0(id,".rds")))
    daily<-got$daily[which(export_rows)];rm(got)
    metadata<-list(model_version=MODEL_VERSION,rule=rule,resolved_sources=maps[[i]],run_id=run_id,
      dataset=job$dataset,source_rds=job$source,monthly_input=job$path,grid_rds=settings$GRID_RDS,
      represented_dates=sort(unique(daily$date)),date_semantics="Original represented Date; UTC end labels do not shift it.",
      history_semantics="Available observations within preceding change_days calendar days; current day excluded.",
      missing_semantics="No interpolation. Unknown is not zero. See assessment_complete and alert_triggered.",
      velocity_semantics=if(rule$source_type=="vector_magnitude") "Magnitude of daily mean components, not mean or maximum instantaneous speed." else NULL,
      source_metadata=source_meta)
    if(nrow(daily)) {
      if(settings$WRITE_DAILY_RDS) {d<-as.data.frame(daily);attr(d,"model_metadata")<-metadata;save_rds(d,file.path(stage,paste0(id,"__daily.rds")));rm(d)}
      sm<-summarize_daily(daily,grid);attr(sm,"model_metadata")<-metadata
      save_rds(sm,file.path(stage,paste0(id,"__summary.rds")));fwrite(sm,file.path(stage,paste0(id,"__summary.csv")),na="NA")
      write_maps(sm,file.path(stage,id),rule,settings,template,paste(min(daily$date),max(daily$date),sep=" to "))
    } else sm<-NULL
    pp<-if(defer_merge||is.null(previous)) "" else file.path(previous,"period",paste0(id,".rds"))
    old<-if(nzchar(pp)&&file.exists(pp)) readRDS(pp) else NULL
    running<-if(is.null(sm)) old else merge_summaries(old,sm)
    if(!is.null(running)) {
      # Period metadata describes the entire emitted period rather than only the latest month.
      pm<-metadata;pm$represented_dates<-NULL;pm$period_start<-if(is.null(old)) as.character(min(daily$date)) else attr(old,"model_metadata")$period_start
      pm$period_end<-if(nrow(daily)) as.character(max(daily$date)) else attr(old,"model_metadata")$period_end
      pm$source_metadata<-NULL;pm$monthly_input<-NULL;pm$source_rds<-NULL
      attr(running,"model_metadata")<-pm;save_rds(running,file.path(stage,"period",paste0(id,".rds")))
    }
    rm(daily,sm,old,running,metadata);gc(FALSE)
  }
  receipt<-list(run_id=run_id,job=job,previous=previous,status=if(defer_merge) "computed" else "complete",worker_pid=Sys.getpid(),
                input_rows=nrow(dt),emitted_rows=sum(export_rows),input_days=as.character(dates),completed_utc=format(Sys.time(),tz="UTC"))
  receipt$files<-list.files(stage,recursive=TRUE,all.files=FALSE)
  receipt$sizes<-unname(file.info(file.path(stage,receipt$files))$size)
  save_rds(receipt,file.path(stage,"complete.rds"));rm(dt);gc(FALSE)
  if(dir.exists(job$output)) fail("Output appeared during processing; refusing overwrite.")
  destination<-if(defer_merge) job$computed_output else job$output
  if(dir.exists(destination)||!file.rename(stage,destination)) fail("Could not publish monthly result folder.")
  say("Saved month; worker exiting and releasing memory.")
  list(status=receipt$status,output=job$output,worker_pid=Sys.getpid(),input_rows=receipt$input_rows)
}

final_worker <- function(last,config,settings,run_id) {
  suppressPackageStartupMessages(library(data.table));setDTthreads(settings$THREADS)
  mem.maxVSize(settings$R_VECTOR_HEAP_GB*1024)
  final<-file.path(settings$OUTPUT_DIR,"combined_period",names(last)[1])
  if(dir.exists(final)) {
    r<-file.path(final,"complete.rds")
    if(file.exists(r)) {
      rr<-readRDS(r);fs<-file.path(final,rr$files)
      if(identical(rr$run_id,run_id)&&!is.null(rr$files)&&all(file.exists(fs))&&identical(unname(file.info(fs)$size),rr$sizes)) return(list(status="already_complete",output=final))
    }
    fail("Existing final output is incomplete or belongs to another run: ",final)
  }
  mkdir(dirname(final));stage<-tempfile(".pending_",tmpdir=dirname(final));mkdir(stage);on.exit(unlink(stage,recursive=TRUE),add=TRUE)
  g<-read_grid(settings$GRID_RDS);template<-map_template(g,settings)
  for(i in seq_len(nrow(config))) {
    rule<-config[i,,drop=FALSE];id<-rule$alert_id;p<-file.path(last[[1]],"period",paste0(id,".rds"))
    if(!file.exists(p)) next
    sm<-readRDS(p);save_rds(sm,file.path(stage,paste0(id,"__summary.rds")))
    fwrite(sm,file.path(stage,paste0(id,"__summary.csv")),na="NA")
    meta<-attr(sm,"model_metadata");write_maps(sm,file.path(stage,id),rule,settings,template,paste(meta$period_start,meta$period_end,sep=" to "))
    rm(sm);gc(FALSE)
  }
  files<-list.files(stage,recursive=TRUE)
  save_rds(list(run_id=run_id,status="complete",files=files,sizes=unname(file.info(file.path(stage,files))$size)),file.path(stage,"complete.rds"))
  if(!file.rename(stage,final)) fail("Could not publish period summaries.")
  list(status="complete",output=final)
}

start_child <- function(fun,args,settings) {
  env<-callr::rcmd_safe_env();env[c("R_MAX_VSIZE","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS")]<-c(paste0(floor(settings$R_VECTOR_HEAP_GB*1024),"M"),"1","1")
  p<-callr::r_bg(function(script,fun,args) {
    options(water_past.run=FALSE);e<-new.env(parent=globalenv());sys.source(script,envir=e);do.call(e[[fun]],args)
  },args=list(script=settings$SCRIPT,fun=fun,args=args),libpath=.libPaths(),env=env,
    stdout="|",stderr="|",supervise=TRUE,user_profile=FALSE,system_profile=FALSE,arch=settings$WORKER_ARCH)
  p
}

run_child <- function(fun,args,settings) {
  p<-start_child(fun,args,settings)
  on.exit(try(if(p$is_alive()) p$kill(),silent=TRUE),add=TRUE)
  last<-Sys.time()
  drain<-function() {lines<-c(p$read_output_lines(),p$read_error_lines());if(length(lines)){for(line in lines) message(line);last<<-Sys.time()}}
  while(p$is_alive()) {p$poll_io(250);drain();if(as.numeric(difftime(Sys.time(),last,units="secs"))>15){say("Worker is running; only one month is in memory.");last<-Sys.time()};Sys.sleep(.05)}
  drain();result<-p$get_result();rm(p);gc(FALSE);result
}

prepare_source <- function(path,dataset,settings) {
  filename_dates<-source_date_range(path,require_dates=FALSE)
  sig<-file_signature(path);key<-substr(object_hash(sig),1,16)
  folder<-file.path(settings$CACHE_DIR,dataset,key);mkdir(dirname(folder))
  receipt<-file.path(folder,"preparation_receipt.rds")
  if(dir.exists(folder)) {
    if(!file.exists(receipt)||!identical(readRDS(receipt),sig)) fail("Incomplete/stale preparation folder: ",folder,". Rename/remove that generated folder and rerun.")
  } else {
    split_large_rds_monthly(input_rds=path,output_dir=folder,date_column="date",temp_dir=settings$TEMP_DIR,
      python_exe=settings$PYTHON_EXE,memory_gb=settings$SPLIT_MEMORY_GB,compress=TRUE,
      reserve_disk_gb=settings$RESERVE_DISK_GB,max_staging_gb=settings$MAX_STAGING_GB,group_timezone="UTC")
    save_rds(sig,receipt)
  }
  m<-jsonlite::fromJSON(file.path(folder,"split_manifest.json"))
  if(m$status!="complete"||m$date_column!="date") fail("Unverified monthly preparation manifest.")
  jobs<-data.table::as.data.table(m$outputs)
  if(!is.null(filename_dates) && (min(as.Date(jobs$first_date))!=filename_dates$first_date||max(as.Date(jobs$last_date))!=filename_dates$last_date)) fail("Actual data dates disagree with the source filename: ",basename(path))
  jobs[,`:=`(path=file.path(folder,filename),dataset=dataset,source=normalizePath(path,winslash="/",mustWork=TRUE),cache_folder=folder)]
  if(any(!file.exists(jobs$path))||any(file.info(jobs$path)$size!=jobs$bytes)) fail("Prepared monthly files missing/changed in ",folder)
  jobs[,output:=file.path(settings$OUTPUT_DIR,"months",dataset,month,paste0(first_date,"_to_",last_date,"_",key))]
  jobs
}

# Read only the filename here: excluded source RDS files are never opened.
source_date_range <- function(path,require_dates=TRUE) {
  hit<-regmatches(basename(path),regexec("_([0-9]{4}-[0-9]{2}-[0-9]{2})_to_([0-9]{4}-[0-9]{2}-[0-9]{2})[.]rds$",basename(path),ignore.case=TRUE))[[1]]
  if(length(hit)!=3L) {
    if(require_dates) fail("Cannot filter an undated filename without reading it: ",basename(path),". Expected ..._YYYY-MM-DD_to_YYYY-MM-DD.rds.")
    return(NULL)
  }
  d<-as.Date(hit[2:3],format="%Y-%m-%d")
  if(anyNA(d)||any(format(d,"%Y-%m-%d")!=hit[2:3])||d[1]>d[2]) fail("Invalid filename date range: ",basename(path))
  list(first_date=d[1],last_date=d[2])
}

# Build the same incoming change history as the sequential checkpoints, using
# ONLY already-selected input files. Read previous months one at a time and
# spill each rule's small history to disk before loading the current month.
# This permits adjacent months to compute concurrently without using unfinished
# outputs, dropping boundary changes, or holding two complete input tables.
seed_selected_history <- function(job,config,settings,grid,directory) {
  active<-which(!is.na(config$change_threshold))
  if(!length(active)||!length(job$history_inputs)) return(invisible(NULL))
  first<-as.numeric(as.Date(job$first_date));ng<-nrow(grid)
  for(hjob in job$history_inputs) {
    if(file.info(hjob$path)$size/1024^3>settings$MAX_MONTH_FILE_GIB) fail("History input exceeds MAX_MONTH_FILE_GIB: ",hjob$path)
    say("Reading selected history input: ",basename(hjob$path))
    h<-readRDS(hjob$path)
    if(!is.data.frame(h)||!all(c("cell_id","cell_x","cell_y","date") %in% names(h))||!inherits(h$date,"Date")) fail("Invalid history table.")
    if(!nrow(h)||anyDuplicated(names(h))||anyNA(h$cell_id)||anyNA(h$date)) fail("Empty or ambiguous history table.")
    d<-as.numeric(h$date)
    if(any(!is.finite(d)|d!=floor(d))||any(format(h$date,"%Y-%m")!=hjob$month)||min(h$date)!=as.Date(hjob$first_date)||max(h$date)!=as.Date(hjob$last_date)) fail("History dates differ from monthly preparation manifest.")
    meta<-attr(h,"variable_metadata",exact=TRUE);data.table::setDT(h)
    if(anyDuplicated(h,by=c("cell_id","date"))) fail("Duplicate cell/date rows in selected history.")
    gi<-match(h$cell_id,grid$cell_id)
    if(anyNA(gi)) fail("History cell IDs absent from coastal grid.")
    if(any(!is.finite(h$cell_x)|!is.finite(h$cell_y)|abs(h$cell_x-grid$cell_x[gi])>1e-6|abs(h$cell_y-grid$cell_y[gi])>1e-6)) fail("History coordinates disagree with coastal grid.")
    mappings<-lapply(active,function(i) resolve_rule(config[i,,drop=FALSE],h,meta))
    needed<-unique(c("cell_id","date",unlist(lapply(mappings,`[[`,"columns"))))
    for(nm in setdiff(names(attributes(h)),c("names","row.names","class",".internal.selfref"))) data.table::setattr(h,nm,NULL)
    for(nm in setdiff(names(h),needed)) data.table::set(h,j=nm,value=NULL)
    # Release the full table before allocating rule-history matrices.
    keep<-which(d>=first-max(config$change_days[active]) & d<first)
    h<-h[keep];rm(d,gi,keep,meta);gc(FALSE)
    if(!nrow(h)) {rm(h);next}
    data.table::setorderv(h,c("date","cell_id"));gi<-match(h$cell_id,grid$cell_id)
    bounds<-h[,list(first=.I[1],last=.I[.N]),by="date"]
    for(k in seq_along(active)) {
      rule<-config[active[k],,drop=FALSE];path<-file.path(directory,paste0(rule$alert_id,".rds"))
      old<-if(file.exists(path)) readRDS(path) else list(dates=numeric(),values=matrix(numeric(),nrow=ng,ncol=0))
      values<-rule_series(h,rule,mappings[[k]])$change
      for(j in which(as.numeric(bounds$date)>=first-rule$change_days)) {
        ix<-seq.int(bounds$first[j],bounds$last[j]);day<-as.numeric(bounds$date[j])
        if(day %in% old$dates) fail("Overlapping selected history dates.")
        v<-rep(NA_real_,ng);v[gi[ix]]<-values[ix]
        old$dates<-c(old$dates,day);old$values<-cbind(old$values,v)
      }
      save_rds(old,path);rm(old,values);gc(FALSE)
    }
    rm(h,gi,mappings,bounds);gc(FALSE)
  }
  invisible(NULL)
}

validate_saved_month <- function(job,previous,run_id) {
  if(!dir.exists(job$output)) return(NULL)
  p<-file.path(job$output,"complete.rds")
  if(!file.exists(p)) fail("Incomplete output directory: ",job$output)
  rr<-readRDS(p);fs<-file.path(job$output,rr$files)
  if(!identical(rr$status,"complete")||!identical(rr$run_id,run_id)||!identical(rr$previous,previous)||any(!file.exists(fs))||any(file.info(fs)$size!=rr$sizes)) fail("Saved month failed resume checks: ",job$output)
  list(status="resumed",output=job$output,worker_pid=rr$worker_pid,input_rows=rr$input_rows)
}

# Merge the small cell summaries in the original order. The computed folder is
# unpublished until its cumulative checkpoint and receipt are both complete.
commit_parallel_month <- function(job,previous,config,run_id) {
  stage<-job$computed_output;p<-file.path(stage,"complete.rds")
  if(!file.exists(p)) fail("Parallel worker did not complete its staging folder.")
  rr<-readRDS(p)
  if(!identical(rr$run_id,run_id)||!identical(rr$status,"computed")||!identical(rr$previous,previous)) fail("Invalid parallel worker receipt.")
  fs<-file.path(stage,rr$files)
  if(any(!file.exists(fs))||any(file.info(fs)$size!=rr$sizes)) fail("Incomplete parallel worker files.")
  for(id in config$alert_id) {
    np<-file.path(stage,"period",paste0(id,".rds"))
    op<-if(is.null(previous)) "" else file.path(previous,"period",paste0(id,".rds"))
    new<-if(file.exists(np)) readRDS(np) else NULL
    old<-if(nzchar(op)&&file.exists(op)) readRDS(op) else NULL
    if(!is.null(old)) {
      if(is.null(new)) new<-old else {
        meta<-attr(new,"model_metadata");meta$period_start<-attr(old,"model_metadata")$period_start
        new<-merge_summaries(old,new);attr(new,"model_metadata")<-meta
      }
      save_rds(new,np)
    }
    rm(new,old);gc(FALSE)
  }
  rr$status<-"complete";rr$job$computed_output<-NULL;rr$job$history_inputs<-NULL
  rr$files<-setdiff(list.files(stage,recursive=TRUE,all.files=FALSE),"complete.rds")
  rr$sizes<-unname(file.info(file.path(stage,rr$files))$size);save_rds(rr,p)
  if(dir.exists(job$output)||!file.rename(stage,job$output)) fail("Could not publish merged monthly result: ",job$output)
  list(status="complete",output=job$output,worker_pid=rr$worker_pid,input_rows=rr$input_rows)
}

run_month_jobs <- function(jobs,cfg,settings,run_id) {
  n<-nrow(jobs);workers<-min(as.integer(settings$PARALLEL_WORKERS),n)
  worker_settings<-settings
  worker_settings$R_VECTOR_HEAP_GB<-min(settings$R_VECTOR_HEAP_GB,floor(settings$PARALLEL_VECTOR_HEAP_GB*1024/workers)/1024)
  say("Modelling workers: ",workers,"; vector-heap limit per worker: ",worker_settings$R_VECTOR_HEAP_GB," GiB (not a total RAM limit).")
  last<-list();log<-list()
  record<-function(i,job,result) {
    last[[job$dataset]]<<-job$output
    log[[i]]<<-data.frame(dataset=job$dataset,month=job$month,first_date=job$first_date,last_date=job$last_date,status=result$status,worker_pid=result$worker_pid,output=result$output)
    data.table::fwrite(data.table::rbindlist(log),file.path(settings$OUTPUT_DIR,"processing_log.csv"));gc(FALSE)
  }
  error_at<-function(job,e) {
    writeLines(conditionMessage(e),file.path(settings$OUTPUT_DIR,"last_error.txt"))
    fail("Stopped at ",job$dataset," ",job$month,". Completed months are retained. ",conditionMessage(e))
  }
  if(workers==1L) {
    # Original chronological processing and history checkpoints.
    for(i in seq_len(n)) {
      job<-as.list(jobs[i]);previous<-last[[job$dataset]]
      say("MONTH ",i,"/",n,": ",job$dataset," ",job$first_date," to ",job$last_date)
      result<-validate_saved_month(job,previous,run_id)
      if(is.null(result)) result<-tryCatch(run_child("month_worker",list(job=job,previous=previous,config=cfg[cfg$dataset==job$dataset,,drop=FALSE],settings=worker_settings,run_id=run_id),worker_settings),error=function(e) error_at(job,e))
      record(i,job,result)
    }
    return(list(last=last,log=log))
  }
  previous<-vector("list",n);paths<-list();ready<-vector("list",n);work_jobs<-vector("list",n)
  for(i in seq_len(n)) {
    job<-as.list(jobs[i]);previous[i]<-list(paths[[job$dataset]]);paths[[job$dataset]]<-job$output
    ready[i]<-list(validate_saved_month(job,previous[[i]],run_id))
    work_jobs[[i]]<-job
  }
  active<-list();temporary<-character();launch<-1L;commit<-1L;last_message<-Sys.time()
  on.exit({
    for(a in active) try({if(a$p$is_alive()) a$p$kill();a$p$wait(5000)},silent=TRUE)
    for(folder in temporary) unlink(folder,recursive=TRUE)
  },add=TRUE)
  drain<-function(a) {
    lines<-c(a$p$read_output_lines(),a$p$read_error_lines())
    if(length(lines)) {for(line in lines) message("[",a$index,"/",n,"] ",line);last_message<<-Sys.time()}
  }
  while(commit<=n) {
    while(launch<=n && length(active)<workers) {
      i<-launch;launch<-launch+1L
      if(!is.null(ready[[i]])) next
      job<-work_jobs[[i]];rules<-cfg[cfg$dataset==job$dataset,,drop=FALSE]
      lookback<-max(c(0,rules$change_days),na.rm=TRUE)
      hi<-which(seq_len(n)<i & jobs$dataset==job$dataset & as.Date(jobs$last_date)>=as.Date(job$first_date)-lookback)
      job$history_inputs<-lapply(hi,function(k) as.list(jobs[k]))
      mkdir(dirname(job$output));job$computed_output<-tempfile(".pending_parallel_",tmpdir=dirname(job$output))
      temporary<-c(temporary,job$computed_output);work_jobs[[i]]<-job
      p<-tryCatch(start_child("month_worker",list(job=job,previous=previous[[i]],config=rules,settings=worker_settings,run_id=run_id,defer_merge=TRUE),worker_settings),error=function(e) error_at(job,e))
      active[[as.character(i)]]<-list(p=p,index=i)
      say("Started month ",i,"/",n,": ",job$dataset," ",job$month,"; active workers ",length(active),"/",workers)
    }
    for(key in names(active)) {
      a<-active[[key]];a$p$poll_io(0);drain(a)
      if(!a$p$is_alive()) {
        drain(a);i<-a$index
        ready[[i]]<-tryCatch(a$p$get_result(),error=function(e) error_at(work_jobs[[i]],e))
        active[[key]]<-NULL
      }
    }
    while(commit<=n && !is.null(ready[[commit]])) {
      job<-work_jobs[[commit]];result<-ready[[commit]]
      if(result$status=="computed") result<-tryCatch(commit_parallel_month(job,previous[[commit]],cfg[cfg$dataset==job$dataset,,drop=FALSE],run_id),error=function(e) error_at(job,e))
      record(commit,job,result);say("Saved/merged month ",commit,"/",n,": ",job$dataset," ",job$month)
      ready[commit]<-list(NULL);commit<-commit+1L
    }
    if(commit<=n) {
      if(!length(active)&&launch>n) fail("Parallel queue ended before all months were committed.")
      if(as.numeric(difftime(Sys.time(),last_message,units="secs"))>=15) {say("Parallel workers running: ",length(active),"; completed results merge in date order.");last_message<-Sys.time()}
      Sys.sleep(.05)
    }
  }
  list(last=last,log=log)
}

run_model <- function(settings=SETTINGS,script=MODEL_SCRIPT) {
  need<-c("data.table","callr","jsonlite","processx")
  missing<-need[!vapply(need,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) fail("Install packages: install.packages(c(",paste(sprintf('"%s"',missing),collapse=","),"))")
  if(.Machine$sizeof.pointer!=8) fail("Use 64-bit R.")
  if(!nzchar(script)||!file.exists(script)) fail("Cannot locate the saved model script.")
  settings$SCRIPT<-normalizePath(script,winslash="/",mustWork=TRUE);base<-dirname(settings$SCRIPT)
  if(!nzchar(settings$INPUT_DIR)) {
    if(interactive()&&.Platform$OS.type=="windows") settings$INPUT_DIR<-utils::choose.dir(caption="Choose folder containing historical daily RDS files")
    else fail("Set INPUT_DIR near the top of the script.")
  }
  if(is.na(settings$INPUT_DIR)||!nzchar(settings$INPUT_DIR)) return(invisible(NULL))
  for(nm in c("INPUT_DIR","GRID_RDS","CONFIG_CSV","OUTPUT_DIR","CACHE_DIR")) settings[[nm]]<-absolute_path(settings[[nm]],base)
  if(!dir.exists(settings$INPUT_DIR)||!file.exists(settings$GRID_RDS)) fail("Input directory or coastal grid does not exist.")
  if(settings$OUTPUT_DIR==settings$INPUT_DIR||settings$CACHE_DIR==settings$INPUT_DIR||settings$CACHE_DIR==settings$OUTPUT_DIR) fail("Input, output and cache directories must be distinct.")
  if(!is.finite(settings$R_VECTOR_HEAP_GB)||settings$R_VECTOR_HEAP_GB<1||settings$R_VECTOR_HEAP_GB>20) fail("Use R_VECTOR_HEAP_GB between 1 and 20; default 16 for 32 GB RAM.")
  if(!is.finite(settings$MAX_MONTH_FILE_GIB)||settings$MAX_MONTH_FILE_GIB<=0||!settings$THREADS %in% 1:64) fail("Invalid file/thread guard.")
  if(is.null(settings$PARALLEL_WORKERS)) settings$PARALLEL_WORKERS<-1L
  if(is.null(settings$PARALLEL_VECTOR_HEAP_GB)) settings$PARALLEL_VECTOR_HEAP_GB<-settings$R_VECTOR_HEAP_GB
  if(length(settings$PARALLEL_WORKERS)!=1L||!is.finite(settings$PARALLEL_WORKERS)||!settings$PARALLEL_WORKERS %in% 1:8) fail("PARALLEL_WORKERS must be an integer from 1 to 8.")
  if(length(settings$PARALLEL_VECTOR_HEAP_GB)!=1L||!is.finite(settings$PARALLEL_VECTOR_HEAP_GB)||settings$PARALLEL_VECTOR_HEAP_GB<settings$PARALLEL_WORKERS||settings$PARALLEL_VECTOR_HEAP_GB>20) fail("PARALLEL_VECTOR_HEAP_GB must be at least PARALLEL_WORKERS and at most 20; default 16 on 32 GB RAM.")
  date_limited<-nzchar(settings$START_DATE)||nzchar(settings$END_DATE)
  mkdir(settings$OUTPUT_DIR);mkdir(settings$CACHE_DIR)
  settings$TEMP_DIR<-if(nzchar(settings$TEMP_DIR)) absolute_path(settings$TEMP_DIR,base) else settings$CACHE_DIR;mkdir(settings$TEMP_DIR)
  settings$START_DATE<-if(nzchar(settings$START_DATE)) as.Date(settings$START_DATE) else as.Date("0001-01-01")
  settings$END_DATE<-if(nzchar(settings$END_DATE)) as.Date(settings$END_DATE) else as.Date("9999-12-31")
  if(anyNA(c(settings$START_DATE,settings$END_DATE))||settings$START_DATE>settings$END_DATE) fail("Invalid START_DATE/END_DATE.")
  cfg<-read_config(settings$CONFIG_CSV,settings)
  grid<-read_grid(settings$GRID_RDS);invisible(map_template(grid,settings));rm(grid);gc(FALSE)
  files<-list.files(settings$INPUT_DIR,pattern="[.]rds$",full.names=TRUE,recursive=settings$RECURSIVE,ignore.case=TRUE)
  files<-normalizePath(files,winslash="/",mustWork=TRUE)
  files<-files[files!=settings$GRID_RDS & !startsWith(files,paste0(settings$OUTPUT_DIR,"/")) & !startsWith(files,paste0(settings$CACHE_DIR,"/"))]
  sources<-list()
  for(ds in names(settings$FILE_PATTERNS)) {
    f<-files[grepl(settings$FILE_PATTERNS[[ds]],basename(files),ignore.case=TRUE)]
    if(!any(cfg$dataset==ds)) {say("No enabled rules for ",ds,"; enable them in the config to process this dataset.");next}
    if(!length(f)) fail("No input files matched enabled dataset ",ds)
    ranges<-lapply(f,source_date_range,require_dates=date_limited)
    keep<-vapply(ranges,function(r) is.null(r)||(r$last_date>=settings$START_DATE && r$first_date<=settings$END_DATE),logical(1))
    say(ds,": selected ",sum(keep),"/",length(f)," source files by inclusive filename date overlap; excluded files will not be opened.")
    f<-f[keep]
    if(!length(f)) {say("No ",ds," files overlap the requested interval.");next}
    for(p in sort(f)) sources[[length(sources)+1L]]<-list(dataset=ds,path=p,signature=file_signature(p))
  }
  if(!length(sources)) fail("No input filename date ranges overlap the requested period.")
  if(anyDuplicated(vapply(sources,`[[`,character(1),"path"))) fail("Input filename matches multiple dataset patterns.")
  signature_settings<-settings;signature_settings$WORKER_ARCH<-NULL
  run_id<-object_hash(list(version=MODEL_VERSION,script=unname(tools::md5sum(settings$SCRIPT)),config=cfg,
      grid=unname(tools::md5sum(settings$GRID_RDS)),settings=signature_settings,sources=sources))
  run_receipt<-file.path(settings$OUTPUT_DIR,"run_definition.rds")
  if(file.exists(run_receipt)) {
    if(!identical(readRDS(run_receipt)$run_id,run_id)) fail("Output folder belongs to different data/config/settings. Choose a new OUTPUT_DIR to keep results consistent.")
  } else {
    save_rds(list(run_id=run_id,settings=settings,config=cfg,sources=sources),run_receipt)
    data.table::fwrite(cfg,file.path(settings$OUTPUT_DIR,"effective_alert_config.csv"),na="")
  }
  # Stream/split sources sequentially. No original data.frame is read here.
  jobs<-data.table::rbindlist(lapply(sources,function(s) prepare_source(s$path,s$dataset,settings)),fill=TRUE)
  data.table::setorderv(jobs,c("dataset","first_date","last_date"))
  for(ds in unique(jobs$dataset)) {
    z<-jobs[jobs$dataset==ds,]
    if(nrow(z)>1&&any(as.Date(z$first_date[-1])<=as.Date(head(z$last_date,-1)))) fail("Overlapping input date ranges for ",ds,". Keep only one version of each cell/day dataset; do not include both original and manually split copies.")
  }
  jobs<-jobs[as.Date(first_date)<=settings$END_DATE & as.Date(last_date)>=settings$START_DATE]
  if(!nrow(jobs)||!any(as.Date(jobs$last_date)>=settings$START_DATE)) fail("No input dates overlap requested period.")
  data.table::fwrite(jobs,file.path(settings$OUTPUT_DIR,"input_manifest.csv"))
  processed<-run_month_jobs(jobs,cfg,settings,run_id)
  last<-processed$last;log<-processed$log
  for(ds in names(last)) run_child("final_worker",list(last=last[ds],config=cfg[cfg$dataset==ds,,drop=FALSE],settings=settings,run_id=run_id),settings)
  if(!settings$KEEP_PREPARED) for(folder in unique(jobs$cache_folder)) unlink(folder,recursive=TRUE)
  say("Finished. Monthly results and combined-period summaries: ",settings$OUTPUT_DIR)
  invisible(data.table::rbindlist(log))
}

# ---- Embedded bounded-memory monthly preparation engine ----
MONTHLY_PYTHON_WORKER <- c(
  "\"\"\"Streaming, lossless row partitioner for ordinary atomic-column R data.frames.",
  "",
  "Copies original serialized payload bytes. No pandas, NumPy or RDS full reader.",
  "Format reference: R Internals, Serialization Formats; binary versions 2 and 3.",
  "Unsupported reference objects and layouts fail before any outputs are published.",
  "\"\"\"",
  "import os, sys, json, struct, math, time, gzip, bz2, lzma",
  "import hashlib, sqlite3, tempfile, shutil, csv, datetime, codecs",
  "from contextlib import contextmanager",
  "",
  "def limit_memory(gib):",
  "    \"\"\"Apply the worker's OS limit BEFORE opening the input or importing extras.\"\"\"",
  "    size = int(gib * 1024**3)",
  "    if struct.calcsize('P') != 8:",
  "        raise RuntimeError('64-bit Python is required.')",
  "    if os.name == 'nt':",
  "        import ctypes as C",
  "        from ctypes import wintypes as W",
  "        class Basic(C.Structure):",
  "            _fields_ = [('PerProcessUserTimeLimit', C.c_longlong),",
  "                        ('PerJobUserTimeLimit', C.c_longlong), ('LimitFlags', W.DWORD),",
  "                        ('MinimumWorkingSetSize', C.c_size_t), ('MaximumWorkingSetSize', C.c_size_t),",
  "                        ('ActiveProcessLimit', W.DWORD), ('Affinity', C.c_size_t),",
  "                        ('PriorityClass', W.DWORD), ('SchedulingClass', W.DWORD)]",
  "        class IO(C.Structure):",
  "            _fields_ = [(n, C.c_ulonglong) for n in ('ReadOperationCount', 'WriteOperationCount',",
  "                         'OtherOperationCount', 'ReadTransferCount', 'WriteTransferCount', 'OtherTransferCount')]",
  "        class Extended(C.Structure):",
  "            _fields_ = [('BasicLimitInformation', Basic), ('IoInfo', IO),",
  "                        ('ProcessMemoryLimit', C.c_size_t), ('JobMemoryLimit', C.c_size_t),",
  "                        ('PeakProcessMemoryUsed', C.c_size_t), ('PeakJobMemoryUsed', C.c_size_t)]",
  "        k = C.WinDLL('kernel32', use_last_error=True)",
  "        k.CreateJobObjectW.argtypes = [C.c_void_p, W.LPCWSTR]",
  "        k.CreateJobObjectW.restype = W.HANDLE",
  "        k.GetCurrentProcess.restype = W.HANDLE",
  "        k.SetInformationJobObject.argtypes = [W.HANDLE, C.c_int, C.c_void_p, W.DWORD]",
  "        k.SetInformationJobObject.restype = W.BOOL",
  "        k.AssignProcessToJobObject.argtypes = [W.HANDLE, W.HANDLE]",
  "        k.AssignProcessToJobObject.restype = W.BOOL",
  "        job = k.CreateJobObjectW(None, None)",
  "        if not job:",
  "            raise C.WinError(C.get_last_error())",
  "        limits = Extended()",
  "        limits.BasicLimitInformation.LimitFlags = 0x100  # PROCESS_MEMORY (commit limit)",
  "        limits.ProcessMemoryLimit = size",
  "        if not k.SetInformationJobObject(job, 9, C.byref(limits), C.sizeof(limits)):",
  "            raise C.WinError(C.get_last_error())",
  "        if not k.AssignProcessToJobObject(job, k.GetCurrentProcess()):",
  "            raise C.WinError(C.get_last_error())",
  "        # Retain the job handle until process exit.",
  "        globals()['_memory_job'] = (k, job)",
  "        return 'Windows per-process committed-memory limit'",
  "    if sys.platform.startswith('linux'):",
  "        import resource",
  "        current = resource.getrlimit(resource.RLIMIT_AS)",
  "        if current[1] != resource.RLIM_INFINITY:",
  "            size = min(size, current[1])",
  "        resource.setrlimit(resource.RLIMIT_AS, (size, size))",
  "        return 'Linux virtual-address-space limit'",
  "    raise RuntimeError('This guarded version supports Windows and Linux. No scan started on this platform.')",
  "",
  "BUF = 4 * 1024**2",
  "WIDTH = {10:4, 13:4, 14:8, 15:16, 24:1}",
  "ATOMIC = set(WIDTH) | {16}",
  "LAST_PROGRESS = 0",
  "",
  "def progress(message, force=False):",
  "    global LAST_PROGRESS",
  "    now = time.monotonic()",
  "    if force or now - LAST_PROGRESS >= 10:",
  "        print(message, flush=True)",
  "        LAST_PROGRESS = now",
  "",
  "def require(ok, message):",
  "    if not ok: raise ValueError(message)",
  "",
  "def room(folder, cfg):",
  "    require(shutil.disk_usage(folder).free >= cfg['reserve_disk_gb'] * 1024**3,",
  "            'Free disk space fell below RESERVE_DISK_GB; stopped without publishing output.')",
  "",
  "@contextmanager",
  "def decoded(path):",
  "    with open(path, 'rb') as raw:",
  "        magic = raw.read(8); raw.seek(0)",
  "        if magic.startswith(b'\\x1f\\x8b'):",
  "            with gzip.GzipFile(fileobj=raw, mode='rb') as f: yield f, 'gzip'",
  "        elif magic.startswith(b'BZh'):",
  "            with bz2.BZ2File(raw,'rb') as f: yield f, 'bzip2'",
  "        elif magic.startswith(b'\\xfd7zXZ\\0'):",
  "            with lzma.LZMAFile(raw,'rb') as f: yield f, 'xz'",
  "        elif magic.startswith(b'\\x28\\xb5\\x2f\\xfd'):",
  "            try: import zstandard",
  "            except ImportError: raise ValueError('zstd input requires: python -m pip install zstandard')",
  "            with zstandard.ZstdDecompressor().stream_reader(raw) as f: yield f, 'zstd'",
  "        else: yield raw, 'none'",
  "",
  "def seekable_source(path, temp, cfg):",
  "    with decoded(path) as (src, compression):",
  "        if compression == 'none': return path, compression",
  "        out = os.path.join(temp,'source_uncompressed.rds')",
  "        progress('Decompressing to temporary DISK storage (not RAM).', True)",
  "        total = 0",
  "        with open(out,'wb') as dst:",
  "            while True:",
  "                b = src.read(BUF)",
  "                if not b: break",
  "                total += len(b)",
  "                require(total <= cfg['max_staging_gb']*1024**3,",
  "                        'Decompressed source exceeds MAX_STAGING_GB; increase it only with sufficient disk space.')",
  "                room(temp,cfg)",
  "                dst.write(b)",
  "                progress('Staged %.2f GiB of uncompressed RDS.' % (total/1024**3))",
  "        return out, compression",
  "",
  "class Node:",
  "    def __init__(self, flags):",
  "        self.flags = flags",
  "        self.typ = flags & 255",
  "        self.n = None",
  "        self.attrs = []",
  "        self.children = []",
  "        self.data = self.end = None",
  "        self.sha = None",
  "        self.text = None",
  "        self.seq = None",
  "",
  "class Parser:",
  "    def __init__(self, f, size=None, verify=False):",
  "        self.f, self.size, self.verify = f, size, verify",
  "        self.pos = 0",
  "        self.endian = '>'",
  "        self.refs = []",
  "        self.nodes = 0",
  "        self.encoding = 'utf-8'",
  "",
  "    def read(self,n):",
  "        require(0 <= n <= BUF, 'Read buffer guard reached.')",
  "        b = self.f.read(n)",
  "        require(len(b)==n, 'Truncated RDS at byte %d.' % self.pos)",
  "        self.pos += n",
  "        return b",
  "",
  "    def i32(self): return struct.unpack(self.endian+'i',self.read(4))[0]",
  "",
  "    def length(self):",
  "        n = self.i32()",
  "        if n == -1:",
  "            hi,lo = struct.unpack(self.endian+'II',self.read(8)); n = (hi<<32)|lo",
  "        require(0 <= n <= 2**52,'Invalid serialized vector length.')",
  "        return n",
  "",
  "    def skip(self,n,h=None):",
  "        require(n>=0,'Negative payload length.')",
  "        if self.size is not None and not self.verify and h is None:",
  "            require(self.pos+n <= self.size,'Vector extends beyond source file.')",
  "            self.f.seek(n,1); self.pos += n",
  "        else:",
  "            while n:",
  "                b=self.read(min(n,BUF)); n-=len(b)",
  "                if h is not None: h.update(b)",
  "                progress('Reading serialized stream: %.2f GiB.' % (self.pos/1024**3))",
  "",
  "    def parse(self):",
  "        magic=self.read(2); require(magic in (b'X\\n',b'B\\n'),'Only binary RDS v2/v3 is supported; no full-file fallback.')",
  "        v=self.read(4)",
  "        self.endian='>' if magic==b'X\\n' else ('<' if struct.unpack('<i',v)[0] in (2,3) else '>')",
  "        self.version=struct.unpack(self.endian+'i',v)[0]",
  "        require(self.version in (2,3),'Unsupported serialization version.')",
  "        extra=self.read(8)",
  "        self.header=magic+v+extra",
  "        if self.version==3:",
  "            b=self.read(4); n=struct.unpack(self.endian+'i',b)[0]",
  "            require(0<=n<=1024,'Invalid native encoding header.')",
  "            enc=self.read(n); self.header+=b+enc",
  "            self.encoding=enc.decode('ascii') or 'utf-8'; codecs.lookup(self.encoding)",
  "        root=self.item()",
  "        require(not self.f.read(1),'Unexpected trailing content after RDS object; refusing to omit it.')",
  "        return root",
  "",
  "    def attributes(self,depth):",
  "        out=[]",
  "        while True:",
  "            flags=self.i32(); typ=flags&255",
  "            if typ in (0,254): return out",
  "            require(typ==2 and flags&1024 and not flags&512,'Unsupported attribute-pair encoding.')",
  "            key=self.item(depth+1)",
  "            require(key.typ==1,'Attribute tag is not a symbol.')",
  "            out.append((key,self.item(depth+1)))",
  "            require(len(out)<=10000,'Too many attributes.')",
  "",
  "    def item(self,depth=0):",
  "        require(depth<=100,'Nested structure exceeds depth limit.')",
  "        self.nodes+=1; require(self.nodes<=100000,'Too many metadata objects for bounded splitting.')",
  "        flags=self.i32(); t=flags&255",
  "        if t==255:",
  "            i=((flags&0xffffffff)>>8) or self.i32()",
  "            require(1<=i<=len(self.refs),'Invalid reference index.')",
  "            return self.refs[i-1]",
  "        z=Node(flags)",
  "        if t in (0,254): return z",
  "        if t==1:",
  "            s=self.item(depth+1); require(s.typ==9,'Unsupported symbol encoding.')",
  "            z.children=[s]; z.text=s.text",
  "            require(len(self.refs)<100000,'Too many symbols.'); self.refs.append(z)",
  "            return z",
  "        if t==9:",
  "            require(not flags&512,'Legacy string attributes are unsupported.')",
  "            z.n=self.i32(); require(z.n>=-1,'Invalid string length.')",
  "            z.data=self.pos",
  "            # Only standalone metadata strings/symbol names use this branch.",
  "            require(z.n<=65536,'A metadata string exceeds 64 KiB.')",
  "            b=self.read(max(0,z.n)); z.end=self.pos",
  "            z.sha=hashlib.sha256(b).hexdigest()",
  "            if z.n<0: z.text=None",
  "            else:",
  "                gp=(flags&0xffffffff)>>12",
  "                enc='utf-8' if gp&(8|64) else ('latin-1' if gp&(2|4) else self.encoding)",
  "                z.text=b.decode(enc,errors='strict')",
  "            return z",
  "        if t==238:",
  "            info=self.item(depth+1); state=self.item(depth+1); z.attrs=self.attributes(depth+1)",
  "            name=info.children[0].text if info.typ==2 and info.children else None",
  "            require(name in ('compact_intseq','compact_realseq'),",
  "                    'Unsupported ALTREP class %r. This file needs a structure-specific splitter.' % name)",
  "            vals=self.values(state)",
  "            require(len(vals)==3 and all(isinstance(v,(int,float)) and math.isfinite(v) for v in vals),'Invalid compact sequence.')",
  "            n,start,step=vals; require(n==int(n) and 0<=n<=2**52,'Invalid sequence length.')",
  "            z.typ=13 if name=='compact_intseq' else 14",
  "            # Convert storage representation only; values and attributes remain unchanged.",
  "            z.flags=z.typ | (256 if flags&256 else 0) | (512 if z.attrs else 0)",
  "            z.n=int(n); z.seq=(start,step)",
  "            return z",
  "        require(t in ATOMIC|{2,19,20},",
  "                'Unsupported R object type %d. Expected a plain data.frame with atomic columns; no data written.' % t)",
  "        if t==2:",
  "            if flags&512: z.attrs=self.attributes(depth+1)",
  "            tag=self.item(depth+1) if flags&1024 else None",
  "            z.children=[self.item(depth+1),self.item(depth+1)]",
  "            z.tag=tag",
  "            return z",
  "        z.n=self.length(); z.data=self.pos",
  "        if t in WIDTH:",
  "            h=hashlib.sha256() if self.verify else None",
  "            self.skip(z.n*WIDTH[t],h)",
  "            z.sha=h.hexdigest() if h is not None else None",
  "        elif t==16:",
  "            h=hashlib.sha256() if self.verify else None",
  "            for i in range(z.n):",
  "                b=self.read(8); f,n=struct.unpack(self.endian+'ii',b)",
  "                require(f&255==9 and not f&512 and n>=-1,'Unsupported character-vector encoding.')",
  "                if h is not None: h.update(b)",
  "                self.skip(max(0,n),h)",
  "                if i%100000==0: progress('Indexing character values: %s / %s.' % (i,z.n))",
  "            z.sha=h.hexdigest() if h is not None else None",
  "        else:",
  "            require(z.n<=100000-self.nodes,'List is too large for the supported flat-table layout.')",
  "            z.children=[self.item(depth+1) for _ in range(z.n)]",
  "        z.end=self.pos",
  "        if flags&512: z.attrs=self.attributes(depth+1)",
  "        return z",
  "",
  "    def values(self,z):",
  "        require(z.n is not None and z.n<=10000,'Structural metadata is too large to interpret safely.')",
  "        if z.seq is not None: return [z.seq[0]+i*z.seq[1] for i in range(z.n)]",
  "        old=self.f.tell(); self.f.seek(z.data)",
  "        try:",
  "            if z.typ in (10,13,14):",
  "                fmt={10:'i',13:'i',14:'d'}[z.typ]; w=WIDTH[z.typ]",
  "                b=self.f.read(z.n*w); require(len(b)==z.n*w,'Truncated metadata.')",
  "                return list(struct.unpack(self.endian+str(z.n)+fmt,b))",
  "            if z.typ==16:",
  "                out=[]; total=0",
  "                for _ in range(z.n):",
  "                    f,n=struct.unpack(self.endian+'ii',self.f.read(8)); total+=max(0,n)",
  "                    require(total<=8*1024**2,'Structural names exceed the metadata limit.')",
  "                    if n==-1: out.append(None); continue",
  "                    gp=(f&0xffffffff)>>12",
  "                    enc='utf-8' if gp&(8|64) else ('latin-1' if gp&(2|4) else self.encoding)",
  "                    out.append(self.f.read(n).decode(enc))",
  "                return out",
  "            raise ValueError('Unsupported structural metadata value type.')",
  "        finally: self.f.seek(old)",
  "",
  "def attr(z,name):",
  "    got=[v for k,v in z.attrs if k.text==name]",
  "    require(len(got)<=1,'Duplicate attribute: '+name)",
  "    return got[0] if got else None",
  "",
  "def signature(z):",
  "    return dict(flags=z.flags,typ=z.typ,n=z.n,sha=z.sha,",
  "                children=[signature(c) for c in z.children],",
  "                attrs=[(k.text,signature(v)) for k,v in z.attrs],",
  "                tag=signature(z.tag) if getattr(z,'tag',None) is not None else None)",
  "",
  "def validate_table(p,root,cfg):",
  "    cl=attr(root,'class')",
  "    classes=p.values(cl) if cl else []",
  "    require(root.typ==19 and 'data.frame' in classes and 'data.table' not in classes,",
  "            'Expected an ordinary data.frame (including an atomic-column tibble). Arrays, nested datasets and data.table reference objects require a tailored splitter.')",
  "    nm=attr(root,'names'); require(nm is not None,'Table has no column names.')",
  "    names=p.values(nm)",
  "    require(len(names)==len(root.children) and all(isinstance(n,str) for n in names),'Invalid column names.')",
  "    require(root.children,'Zero-column tables cannot be split by a date column.')",
  "    n=root.children[0].n",
  "    require(n is not None and 0<n<=2147483647,'This splitter supports 1 to 2,147,483,647 rows.')",
  "    for name,col in zip(names,root.children):",
  "        require(col.typ in ATOMIC and col.n==n and attr(col,'dim') is None,",
  "                'Column %r is not a simple atomic vector of the table row count.' % name)",
  "        cn=attr(col,'names')",
  "        require(cn is None or cn.n==n,'Column names attribute has invalid length.')",
  "    rn=attr(root,'row.names'); require(rn is not None,'Missing row.names.')",
  "    compact=False",
  "    if rn.typ==13 and rn.n==2:",
  "        vals=p.values(rn); compact=(vals[0]==-2147483648 and abs(vals[1])==n)",
  "    require(compact or (rn.typ in (13,16) and rn.n==n),'Unsupported row.names encoding.')",
  "    name=cfg['date_column']",
  "    if not name:",
  "        candidates=[x for x in ('date','time_utc','datetime','time') if names.count(x)==1]",
  "        require(candidates,'Set DATE_COLUMN to the exact Date/POSIXct column name.')",
  "        name=candidates[0]",
  "    require(names.count(name)==1,'DATE_COLUMN is missing or duplicated: '+name)",
  "    date_col=root.children[names.index(name)]",
  "    dc=attr(date_col,'class'); dc=p.values(dc) if dc else []",
  "    require(date_col.typ in (13,14) and ('Date' in dc or 'POSIXct' in dc),",
  "            'DATE_COLUMN must have R class Date or POSIXct. Text dates and unlabelled numeric dates are rejected to avoid guessing.')",
  "    require(cfg['group_timezone']=='UTC','This version groups POSIXct timestamps by UTC; GROUP_TIMEZONE must be UTC. Original timezone attributes are preserved.')",
  "    return names,n,date_col,('days' if 'Date' in dc else 'seconds'),compact",
  "",
  "def values_blocks(p,z,rows=131072):",
  "    for first in range(0,z.n,rows):",
  "        n=min(rows,z.n-first)",
  "        if z.seq is not None:",
  "            yield first,(z.seq[0]+(first+j)*z.seq[1] for j in range(n))",
  "        else:",
  "            p.f.seek(z.data+first*WIDTH[z.typ]); b=p.f.read(n*WIDTH[z.typ])",
  "            require(len(b)==n*WIDTH[z.typ],'Truncated date vector.')",
  "            yield first,(a[0] for a in struct.iter_unpack(p.endian+('d' if z.typ==14 else 'i'),b))",
  "",
  "def build_index(p,date,unit,path,cfg):",
  "    db=sqlite3.connect(path)",
  "    db.execute('PRAGMA cache_size=-32768'); db.execute('PRAGMA temp_store=FILE')",
  "    db.execute('CREATE TABLE runs (month TEXT, start INTEGER, n INTEGER)')",
  "    months={}; previous=None; start=0; count=0; row=0; cached_day=None; cached=None",
  "    epoch=datetime.date(1970,1,1); pending=[]",
  "    for first,vals in values_blocks(p,date):",
  "        for value in vals:",
  "            require(math.isfinite(value) and not (date.typ==13 and value==-2147483648),",
  "                    'Missing/invalid date at original row %d. Nothing is dropped; correct the date before splitting.' % (row+1))",
  "            day=math.floor(value/(86400 if unit=='seconds' else 1))",
  "            if day!=cached_day:",
  "                dt=epoch+datetime.timedelta(days=day)",
  "                iso_date=dt.isoformat()",
  "                cached=(iso_date[:7],iso_date); cached_day=day",
  "            month,iso=cached",
  "            if month not in months:",
  "                require(len(months)<2400,'More than 2,400 months; refusing unexpected date range.')",
  "                months[month]=dict(rows=0,first_date=iso,last_date=iso)",
  "            meta=months[month]; meta['rows']+=1",
  "            meta['first_date']=min(meta['first_date'],iso); meta['last_date']=max(meta['last_date'],iso)",
  "            if month==previous: count+=1",
  "            else:",
  "                if previous is not None: pending.append((previous,start,count))",
  "                previous=month; start=row; count=1",
  "            row+=1",
  "        if pending: db.executemany('INSERT INTO runs VALUES (?,?,?)',pending); pending=[]",
  "        db.commit(); room(os.path.dirname(path),cfg)",
  "        progress('Indexed dates: %s / %s rows; %s months.' % (row,date.n,len(months)))",
  "    if previous is not None: db.execute('INSERT INTO runs VALUES (?,?,?)',(previous,start,count))",
  "    db.execute('CREATE INDEX month_runs ON runs(month,start)'); db.commit()",
  "    require(row==date.n and sum(m['rows'] for m in months.values())==date.n,'Row-accounting validation failed.')",
  "    return db,months",
  "",
  "def ranges(db,month):",
  "    for start,n in db.execute('SELECT start,n FROM runs WHERE month=? ORDER BY start',(month,)):",
  "        yield start,n",
  "",
  "class Writer:",
  "    def __init__(self,p,f,db,month,n,cfg,folder):",
  "        self.p,self.f,self.db,self.month,self.n,self.cfg,self.folder=p,f,db,month,n,cfg,folder",
  "        self.since_check=0",
  "",
  "    def write(self,b,h=None):",
  "        self.f.write(b)",
  "        if h is not None: h.update(b)",
  "        self.since_check+=len(b)",
  "        if self.since_check>=64*1024**2:",
  "            room(self.folder,self.cfg); self.since_check=0",
  "            progress('Writing month %s in bounded blocks.' % self.month)",
  "",
  "    def i32(self,n): self.write(struct.pack(self.p.endian+'i',n))",
  "    def length(self,n):",
  "        if n<=2147483647: self.i32(n)",
  "        else: self.i32(-1); self.write(struct.pack(self.p.endian+'II',n>>32,n&0xffffffff))",
  "",
  "    def copy(self,offset,n,h):",
  "        self.p.f.seek(offset)",
  "        while n:",
  "            b=self.p.f.read(min(n,BUF)); require(b,'Truncated payload during copy.')",
  "            self.write(b,h); n-=len(b)",
  "",
  "    def attrs(self,attributes,subset_names=False,root=False,compact=False):",
  "        out=[]",
  "        for key,val in attributes:",
  "            self.i32(1026); self.emit(key)",
  "            if root and key.text=='row.names':",
  "                sig=self.row_names(val) if compact else self.emit(val,True)",
  "            else: sig=self.emit(val,subset_names and key.text=='names')",
  "            out.append((key.text,sig))",
  "        self.i32(254)",
  "        return out",
  "",
  "    def selected(self,subset,z):",
  "        return ranges(self.db,self.month) if subset else iter([(0,z.n)])",
  "",
  "    def sequence(self,z,subset,h):",
  "        fmt='i' if z.typ==13 else 'd'; start,step=z.seq",
  "        for a,n in self.selected(subset,z):",
  "            for i in range(a,a+n,65536):",
  "                count=min(65536,a+n-i)",
  "                vals=[start+(i+j)*step for j in range(count)]",
  "                if z.typ==13: vals=[int(v) for v in vals]",
  "                self.write(struct.pack(self.p.endian+str(count)+fmt,*vals),h)",
  "",
  "    def characters(self,z,subset,h):",
  "        # Character payloads are scanned sequentially; their original encodings",
  "        # and NA markers are copied byte for byte. No strings are decoded here.",
  "        selected=iter(self.selected(subset,z)); current=next(selected,None)",
  "        self.p.f.seek(z.data)",
  "        for i in range(z.n):",
  "            while current is not None and i>=current[0]+current[1]: current=next(selected,None)",
  "            if current is None: break",
  "            b=self.p.f.read(8); require(len(b)==8,'Truncated character header.')",
  "            flags,n=struct.unpack(self.p.endian+'ii',b)",
  "            require(flags&255==9 and n>=-1,'Invalid character payload.')",
  "            keep=i>=current[0]",
  "            if keep: self.write(b,h)",
  "            remain=max(0,n)",
  "            if not keep: self.p.f.seek(remain,1)",
  "            else:",
  "                while remain:",
  "                    block=self.p.f.read(min(BUF,remain)); require(block,'Truncated string.')",
  "                    self.write(block,h); remain-=len(block)",
  "            if i%100000==0: progress('Copying character values for %s: %d / %d.' % (self.month,i,z.n))",
  "",
  "    def row_names(self,z):",
  "        # Preserve original automatic row identities, including interleaved months.",
  "        self.i32(13); self.length(self.n); h=hashlib.sha256()",
  "        for a,n in ranges(self.db,self.month):",
  "            for i in range(a,a+n,65536):",
  "                k=min(65536,a+n-i)",
  "                self.write(struct.pack(self.p.endian+str(k)+'i',*range(i+1,i+k+1)),h)",
  "        return dict(flags=13,typ=13,n=self.n,sha=h.hexdigest(),children=[],attrs=[],tag=None)",
  "",
  "    def emit(self,z,subset=False,root=False,compact=False):",
  "        flags=z.flags; n=self.n if subset else z.n",
  "        self.i32(flags)",
  "        sig=dict(flags=flags,typ=z.typ,n=n,sha=None,children=[],attrs=[],tag=None)",
  "        if z.typ in (0,254): return sig",
  "        if z.typ==1:",
  "            sig['children']=[self.emit(z.children[0])]; return sig",
  "        if z.typ==9:",
  "            self.i32(z.n); h=hashlib.sha256(); self.copy(z.data,max(0,z.n),h)",
  "            sig['sha']=h.hexdigest(); return sig",
  "        if z.typ==2:",
  "            if flags&512: sig['attrs']=self.attrs(z.attrs)",
  "            if flags&1024: sig['tag']=self.emit(z.tag)",
  "            sig['children']=[self.emit(c) for c in z.children]; return sig",
  "        self.length(n)",
  "        if z.typ in (19,20):",
  "            require(not subset,'Cannot subset a list column.')",
  "            sig['children']=[self.emit(c,subset=root) for c in z.children]",
  "        else:",
  "            h=hashlib.sha256()",
  "            if z.seq is not None: self.sequence(z,subset,h)",
  "            elif z.typ==16: self.characters(z,subset,h)",
  "            else:",
  "                for a,k in self.selected(subset,z): self.copy(z.data+a*WIDTH[z.typ],k*WIDTH[z.typ],h)",
  "            sig['sha']=h.hexdigest()",
  "        if flags&512: sig['attrs']=self.attrs(z.attrs,subset_names=subset,root=root,compact=compact)",
  "        return sig",
  "",
  "def verify(path,expected):",
  "    progress('Verifying every serialized value and attribute in '+os.path.basename(path),True)",
  "    with decoded(path) as (f,kind):",
  "        p=Parser(f,verify=True); got=p.parse()",
  "    require(signature(got)==expected,'Output verification failed; the batch will not be published.')",
  "",
  "def file_sha(path):",
  "    h=hashlib.sha256()",
  "    with open(path,'rb') as f:",
  "        while True:",
  "            b=f.read(BUF)",
  "            if not b: break",
  "            h.update(b)",
  "            progress('Computing file checksum: '+os.path.basename(path))",
  "    return h.hexdigest()",
  "",
  "def split(cfg):",
  "    source=os.path.abspath(cfg['input_rds']); destination=os.path.abspath(cfg['output_dir'])",
  "    require(os.path.isfile(source),'Input RDS does not exist.')",
  "    require(not os.path.exists(destination),'OUTPUT_DIR already exists. Choose a new folder; existing files are never replaced.')",
  "    parent=os.path.dirname(destination); require(os.path.isdir(parent),'Parent of OUTPUT_DIR must exist.')",
  "    temporary=os.path.abspath(cfg['temp_dir']) if cfg['temp_dir'] else parent",
  "    require(os.path.isdir(temporary),'TEMP_DIR does not exist.')",
  "    original=os.stat(source)",
  "    stage=tempfile.mkdtemp(prefix='.monthly_output_',dir=parent)",
  "    work=tempfile.mkdtemp(prefix='rds_monthly_work_',dir=temporary)",
  "    committed=False; db=None",
  "    try:",
  "        local,compression=seekable_source(source,work,cfg)",
  "        with open(local,'rb') as f:",
  "            p=Parser(f,os.path.getsize(local)); root=p.parse()",
  "            names,n,date,unit,compact=validate_table(p,root,cfg)",
  "            chosen=names[root.children.index(date)]",
  "            progress('Validated %d rows, %d columns. Grouping by %s (%s).' % (n,len(names),chosen,unit),True)",
  "            db,months=build_index(p,date,unit,os.path.join(work,'row_ranges.sqlite'),cfg)",
  "            room(stage,cfg)",
  "            stem=os.path.splitext(os.path.basename(source))[0]",
  "            entries=[]",
  "            for month,meta in sorted(months.items()):",
  "                filename='%s_%s_to_%s.rds' % (stem,meta['first_date'],meta['last_date'])",
  "                path=os.path.join(stage,filename)",
  "                progress('Writing %s: %d rows.' % (filename,meta['rows']),True)",
  "                with open(path,'xb') as raw:",
  "                    if cfg['compress']:",
  "                        out=gzip.GzipFile(fileobj=raw,mode='wb',compresslevel=1,mtime=0)",
  "                    else: out=raw",
  "                    try:",
  "                        out.write(p.header)",
  "                        writer=Writer(p,out,db,month,meta['rows'],cfg,stage)",
  "                        expected=writer.emit(root,root=True,compact=compact)",
  "                    finally:",
  "                        if out is not raw: out.close()",
  "                verify(path,expected)",
  "                entries.append(dict(month=month,filename=filename,**meta,bytes=os.path.getsize(path),sha256=file_sha(path),verified=True))",
  "            now=os.stat(source)",
  "            require((now.st_size,now.st_mtime_ns)==(original.st_size,original.st_mtime_ns),'Source file changed during splitting; batch not published.')",
  "            require(sum(e['rows'] for e in entries)==n,'Output row counts do not sum to source row count.')",
  "            with open(os.path.join(stage,'original_row_ranges.csv'),'w',newline='',encoding='utf-8') as out:",
  "                w=csv.writer(out); w.writerow(['month','original_first_row','original_last_row','rows'])",
  "                for month,a,k in db.execute('SELECT month,start,n FROM runs ORDER BY start'):",
  "                    w.writerow([month,a+1,a+k,k])",
  "            manifest=dict(status='complete',source=source,source_bytes=original.st_size,source_compression=compression,",
  "                          source_rows=n,columns=names,date_column=chosen,group_timezone='UTC',memory_limit_gb=cfg['memory_gb'],",
  "                          outputs=entries,notes=[",
  "                              'All source rows occur once across the monthly data tables; original order is retained within each month.',",
  "                              'Atomic data payloads and attributes are copied without numeric or text conversion; compact sequences are materialized exactly.',",
  "                              'Automatic row.names are expanded to original row identities. The row-range CSV records how to restore global row order.',",
  "                              'Source-wide metadata is copied unchanged into every file; embedded source counts/date extents remain source-wide.',",
  "                              'Column names, types, factor levels, timezone attributes and missing-value payloads are preserved.',",
  "                              'Each output is re-read as a bounded stream and its complete parsed structure and payload SHA-256 digests are checked.',",
  "                              'This does not guarantee that an individual month will fit in RAM when later opened with readRDS().'])",
  "            with open(os.path.join(stage,'split_manifest.json'),'w',encoding='utf-8') as out: json.dump(manifest,out,ensure_ascii=False,indent=2)",
  "            require(not os.path.exists(destination),'Output folder was created by another process; refusing to replace it.')",
  "            os.rename(stage,destination); committed=True",
  "            progress('COMPLETE: %d monthly RDS files; %d / %d rows retained.' % (len(entries),sum(e['rows'] for e in entries),n),True)",
  "            return dict(status='complete',output_dir=destination,months=len(entries),source_rows=n)",
  "    finally:",
  "        if db is not None: db.close()",
  "        shutil.rmtree(work,ignore_errors=True)",
  "        if not committed: shutil.rmtree(stage,ignore_errors=True)",
  "",
  "def main():",
  "    # Keep progress readable when Windows redirects Python output to RStudio.",
  "    if hasattr(sys.stdout,'reconfigure'): sys.stdout.reconfigure(encoding='utf-8',errors='replace')",
  "    if hasattr(sys.stderr,'reconfigure'): sys.stderr.reconfigure(encoding='utf-8',errors='replace')",
  "    with open(sys.argv[1],encoding='utf-8-sig') as f: cfg=json.load(f)",
  "    try:",
  "        require(0.25<=cfg['memory_gb']<=8,'Memory guard must be 0.25 to 8 GiB.')",
  "        guard=limit_memory(cfg['memory_gb'])",
  "        progress('Worker memory guard: %s GiB (%s).' % (cfg['memory_gb'],guard),True)",
  "        result=split(cfg)",
  "    except MemoryError: result=dict(status='failed',error='Worker memory limit reached. No completed batch was published.')",
  "    except Exception as exc: result=dict(status='failed',error=type(exc).__name__+': '+str(exc))",
  "    with open(cfg['result_json'],'w',encoding='utf-8') as f: json.dump(result,f)",
  "    if result['status']!='complete': print(result['error'],flush=True); return 2",
  "    return 0",
  "",
  "if __name__=='__main__': sys.exit(main())"
)

monthly_find_python <- function(exe="") {
  candidates <- if (nzchar(exe)) exe else unique(unname(Sys.which(c("python","python3","py"))))
  for (candidate in candidates[nzchar(candidates)]) {
    prefix <- if(tolower(basename(candidate)) %in% c("py","py.exe")) "-3" else character(0)
    probe <- tryCatch(processx::run(candidate,c(prefix,"-c",
      "import sys,struct; assert sys.version_info >= (3,9) and struct.calcsize('P')==8; print(sys.executable)"),
      timeout=10000,error_on_status=FALSE),error=function(e) NULL)
    if (!is.null(probe) && identical(probe$status,0L) && file.exists(trimws(probe$stdout))) return(trimws(probe$stdout))
  }
  stop("64-bit Python 3.9+ was not found. Install it, restart RStudio, or set PYTHON_EXE. No source data was loaded.",call.=FALSE)
}

split_large_rds_monthly <- function(input_rds,output_dir="",date_column="",temp_dir="",
                                   python_exe=PYTHON_EXE,memory_gb=WORKER_MEMORY_GB,
                                   compress=COMPRESS_OUTPUT,reserve_disk_gb=RESERVE_DISK_GB,
                                   max_staging_gb=MAX_STAGING_GB,group_timezone=GROUP_TIMEZONE) {
  need <- c("jsonlite","processx")
  missing <- need[!vapply(need,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) stop("Install required R packages: install.packages(c(\"jsonlite\", \"processx\"))",call.=FALSE)
  if(length(input_rds)!=1L || !is.character(input_rds) || is.na(input_rds) || !file.exists(input_rds) || dir.exists(input_rds)) stop("Select one existing RDS file.")
  input_rds <- normalizePath(input_rds,winslash="/",mustWork=TRUE)
  for(x in list(output_dir,date_column,temp_dir,python_exe)) {
    if(!is.character(x) || length(x)!=1L || is.na(x)) stop("Paths and date_column must each be a single text string.")
  }
  if(!nzchar(output_dir)) output_dir <- paste0(tools::file_path_sans_ext(input_rds),"_monthly")
  output_dir <- path.expand(output_dir)
  if(file.exists(output_dir) || dir.exists(output_dir)) stop("OUTPUT_DIR already exists. Choose a NEW folder; existing outputs are never replaced.")
  if(!dir.exists(dirname(output_dir))) stop("The parent of OUTPUT_DIR must already exist.")
  output_dir <- file.path(normalizePath(dirname(output_dir),winslash="/",mustWork=TRUE),basename(output_dir))
  if(nzchar(temp_dir)) temp_dir <- normalizePath(temp_dir,winslash="/",mustWork=TRUE)
  scalar <- function(x,low,high) is.numeric(x) && length(x)==1L && is.finite(x) && x>=low && x<=high
  if(!scalar(memory_gb,0.25,8)) stop("memory_gb must be between 0.25 and 8; keep 4 on a 32 GB computer.")
  if(!scalar(reserve_disk_gb,0.1,1e5) || !scalar(max_staging_gb,1,1e5)) stop("Invalid disk guard settings.")
  if(!is.logical(compress) || length(compress)!=1L || is.na(compress)) stop("compress must be TRUE or FALSE.")
  if(!identical(group_timezone,"UTC")) stop("This version groups POSIXct dates by UTC; keep GROUP_TIMEZONE = 'UTC'.")
  python <- monthly_find_python(python_exe)
  jobdir <- tempfile("monthly_rds_job_"); dir.create(jobdir)
  proc <- NULL
  on.exit({
    if(!is.null(proc)) try(if(proc$is_alive()) proc$kill(),silent=TRUE)
    unlink(jobdir,recursive=TRUE)
  },add=TRUE)
  script <- file.path(jobdir,"split_worker.py")
  job <- file.path(jobdir,"job.json"); result_path <- file.path(jobdir,"result.json")
  writeLines(MONTHLY_PYTHON_WORKER,script,useBytes=TRUE)
  config <- list(input_rds=input_rds,output_dir=output_dir,date_column=date_column,
                 temp_dir=temp_dir,memory_gb=memory_gb,compress=compress,
                 reserve_disk_gb=reserve_disk_gb,max_staging_gb=max_staging_gb,
                 group_timezone=group_timezone,result_json=result_path)
  jsonlite::write_json(config,job,auto_unbox=TRUE,pretty=FALSE,digits=NA)
  message("Splitting: ",input_rds,"\nNew output folder: ",output_dir,
          "\nMemory limit: ",memory_gb," GiB. Compressed input requires temporary DISK space.")
  proc <- processx::process$new(python,c("-u",script,job),stdout="|",stderr="|",cleanup=TRUE,windows_hide_window=TRUE)
  last <- Sys.time(); started <- last; recent <- character(0)
  drain <- function() {
    lines <- c(proc$read_output_lines(),proc$read_error_lines())
    if(length(lines)) {
      for(line in lines) message(line)
      recent <<- tail(c(recent,lines),20L); last <<- Sys.time()
    }
  }
  while(proc$is_alive()) {
    proc$poll_io(250); drain()
    if(as.numeric(difftime(Sys.time(),last,units="secs"))>=10) {
      message("Working... ",round(as.numeric(difftime(Sys.time(),started,units="secs")))," seconds elapsed; worker limit ",memory_gb," GiB.")
      last <- Sys.time()
    }
    Sys.sleep(0.05)
  }
  drain()
  if(!file.exists(result_path)) stop("Worker stopped without a completion report. The original is unchanged. Last output: ",paste(tail(recent,3),collapse=" "),call.=FALSE)
  result <- jsonlite::fromJSON(result_path)
  if(!identical(proc$get_exit_status(),0L) || !identical(result$status,"complete")) stop(result$error,call.=FALSE)
  message("Saved ",result$months," verified monthly RDS files containing all ",result$source_rows," source rows.\nFolder: ",result$output_dir,
          "\nSee split_manifest.json for counts/checksums and original_row_ranges.csv for original row order.")
  invisible(result)
}

# Run by sourcing this saved script, or by Rscript with optional named paths:
# Rscript early_warning_system_water_past.R --input="D:/past" --grid="D:/past/MSI_coastal_grid.rds" --out="D:/past_results"
if(isTRUE(RUN_MODEL) && !identical(getOption("water_past.run"),FALSE)) {
  a<-if(interactive()) character() else commandArgs(trailingOnly=TRUE)
  aliases<-c(input="INPUT_DIR",grid="GRID_RDS",config="CONFIG_CSV",out="OUTPUT_DIR",cache="CACHE_DIR",temp="TEMP_DIR",python="PYTHON_EXE",start="START_DATE",end="END_DATE")
  for(arg in a) {
    m<-regexec("^--([^=]+)=(.*)$",arg);v<-regmatches(arg,m)[[1]]
    if(length(v)!=3||!v[2] %in% names(aliases)) fail("Use --input=..., --grid=..., --config=..., --out=..., --cache=..., --temp=..., --python=..., --start=..., --end=...")
    SETTINGS[[aliases[[v[2]]]]]<-v[3]
  }
  run_model()
}
