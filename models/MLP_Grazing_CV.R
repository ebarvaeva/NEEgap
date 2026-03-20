# =============================================================================
# MLP_Grazing_CV.R — MLP: BASE + Grazing_days_since (PI disabled)
# =============================================================================
#
# PURPOSE
#   Variant of MLP_CV.R for the management variable ablation study.
#   Identical in all respects to MLP_CV.R except that:
#
#     PI_ENABLED <- FALSE   # Hardcoded for this ablation variant
#
#   is hardcoded, preventing the Phytomass Index regardless of the
#   predictor list.  Grazing_days_since is passed as a direct numeric predictor.
#
#   See MLP_CV.R for full section documentation.
#
# OUTPUTS  (written to RESULTS_DIR)
#   df_cv_all_predictions.rds, training_loss_curves/, progress.log, run_info.txt
#
# =============================================================================


# =============================================================================
# SECTION 1 — Reproducibility
# =============================================================================
set.seed(42)
Sys.setenv(PYTHONHASHSEED="0", CUDA_VISIBLE_DEVICES="-1",
           OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           TF_NUM_INTRAOP_THREADS="1", TF_INTEROP_THREADS="1",
           TF_DETERMINISTIC_OPS="1")


# =============================================================================
# SECTION 2 — Progress logger
# =============================================================================
if (!exists("RESULTS_DIR", inherits=TRUE) || is.null(RESULTS_DIR))
  RESULTS_DIR <- file.path(tempdir(), "mlp_cv_fallback")
dir.create(RESULTS_DIR, recursive=TRUE, showWarnings=FALSE)
.run_log <- list()
log_msg <- function(...) {
  txt <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse=""))
  message(txt); .run_log <<- append(.run_log, list(txt))
  try(silent=TRUE, {writeLines(unlist(.run_log), file.path(RESULTS_DIR,"progress.log"))
    saveRDS(.run_log, file.path(RESULTS_DIR,"progress_log.rds"))}); invisible(txt)
}
log_msg("MLP (NEE gap-fill + post-hoc Reco/GPP) started.  RESULTS_DIR = ", RESULTS_DIR)


# =============================================================================
# SECTION 3 — Helper functions
# =============================================================================
partition_timeseries_into_blocks <- function(n, base_size, tol) {
  stopifnot(n>0,base_size>=1,tol>=0)
  k<-max(1L,round(n/base_size)); sizes<-rep(base_size,k); delta<-n-sum(sizes)
  if(delta!=0){m<-min(k,max(1L,ceiling(abs(delta)/max(1,tol)))); step<-sign(delta)*floor(abs(delta)/m); rem<-sign(delta)*(abs(delta)-abs(step)*m)
  sizes[(k-m+1):k]<-sizes[(k-m+1):k]+step; if(rem!=0) sizes[(k-abs(rem)+1):k]<-sizes[(k-abs(rem)+1):k]+sign(rem)}
  for(i in k:2){dev<-sizes[i]-base_size; if(abs(dev)>tol){ex<-sign(dev)*(abs(dev)-tol); sizes[i]<-base_size+sign(dev)*tol; sizes[i-1]<-sizes[i-1]+ex}}
  sizes[sizes<1]<-1L; sizes[k]<-sizes[k]+(n-sum(sizes)); starts<-c(1L,head(cumsum(sizes)+1L,-1L)); ends<-cumsum(sizes)
  stopifnot(ends[length(ends)]==n); list(sizes=sizes,indices=purrr::map2(starts,ends,seq.int))}

create_cv_gaps <- function(flux_data, base_size, tol, prefix) {
  n<-nrow(flux_data); bl<-partition_timeseries_into_blocks(n,as.integer(round(base_size)),as.integer(round(tol)))
  fn<-paste0(prefix,seq_along(bl$indices)); nn<-paste0("NEE_",fn)
  fm<-matrix(FALSE,nrow=n,ncol=length(fn),dimnames=list(NULL,fn)); nm<-matrix(NA_real_,nrow=n,ncol=length(fn),dimnames=list(NULL,nn))
  for(i in seq_along(bl$indices)){fm[bl$indices[[i]],i]<-TRUE; nm[,i]<-flux_data$NEE_orig; nm[bl$indices[[i]],i]<-NA_real_}
  list(data=dplyr::bind_cols(flux_data,tibble::as_tibble(as.data.frame(fm)),tibble::as_tibble(as.data.frame(nm))),sizes=bl$sizes)}

.vp<-function(o,p) is.finite(as.numeric(o))&is.finite(as.numeric(p))
calc_mae  <-function(o,p){ok<-.vp(o,p); if(!any(ok)) return(NA_real_); mean(abs(as.numeric(p)[ok]-as.numeric(o)[ok]))}
calc_rmse <-function(o,p){ok<-.vp(o,p); if(!any(ok)) return(NA_real_); sqrt(mean((as.numeric(p)[ok]-as.numeric(o)[ok])^2))}
calc_r2   <-function(o,p){ok<-.vp(o,p); if(sum(ok)<2) return(NA_real_); ov<-as.numeric(o)[ok]; pv<-as.numeric(p)[ok]; den<-sum((ov-mean(ov))^2)*sum((pv-mean(pv))^2); if(den<=0) return(NA_real_); (sum((ov-mean(ov))*(pv-mean(pv)))^2)/den}

timestamp_to_season<-function(x){m<-as.integer(format(x,"%m")); dplyr::case_when(m%in%c(11,12,1)~"Winter",m%in%c(2,3,4)~"Spring",m%in%c(5,6,7)~"Summer",m%in%c(8,9,10)~"Autumn",TRUE~NA_character_)}


# =============================================================================
# SECTION 4 — Phytomass Index  (managed grassland only)
# =============================================================================
PI_ENABLED <- FALSE
if (!PI_ENABLED) log_msg("Phytomass Index disabled.")

.find_ppfd_col<-function(df){p<-c("PPFD","PPFD_IN","PPFD_IN_F","PPFD_F","PPFD_1_1_1_gapfilled","PPFD_1_1_1"); h<-p[p%in%names(df)]; if(length(h)) return(h[1]); h2<-grep("PPFD",names(df),ignore.case=TRUE,value=TRUE); if(length(h2)) return(h2[1]); NULL}

add_phytomass_index<-function(flux_data,gap_prefixes=c("VL","L","M","S"),ppfd_col=NULL,night_ppfd=1,day_ppfd=700,window_days=21L){
  if(!requireNamespace("zoo",quietly=TRUE)) stop("Install 'zoo' first.")
  if(!inherits(flux_data$timestamp,"POSIXt")) flux_data$timestamp<-suppressWarnings(lubridate::parse_date_time(flux_data$timestamp,orders=c("Ymd HMS","Ymd HM","Ymd","Y/m/d HMS","Y/m/d")))
  if(is.null(ppfd_col)) ppfd_col<-.find_ppfd_col(flux_data); if(is.null(ppfd_col)||!ppfd_col%in%names(flux_data)){log_msg("PI: no PPFD — skip."); return(flux_data)}
  window_days<-as.integer(window_days); if(!is.finite(window_days)||window_days<1) window_days<-21L; if(window_days%%2==0) window_days<-window_days+1L
  flux_data$.cd<-as.Date(flux_data$timestamp); nv<-as.numeric(flux_data$NEE_orig); pv<-suppressWarnings(as.numeric(flux_data[[ppfd_col]]))
  ok<-is.finite(nv)&is.finite(pv); isn<-ok&pv<night_ppfd; isd<-ok&pv>day_ppfd
  dt<-tibble::tibble(date=flux_data$.cd,nn=dplyr::if_else(isn,nv,0),cn=as.integer(isn),nd=dplyr::if_else(isd,nv,0),cd=as.integer(isd))%>%
    dplyr::group_by(date)%>%dplyr::summarise(nn=sum(nn,na.rm=T),cn=sum(cn,na.rm=T),nd=sum(nd,na.rm=T),cd=sum(cd,na.rm=T),.groups="drop")%>%dplyr::arrange(date)
  rol<-function(x) zoo::rollapply(x,window_days,function(v)sum(v,na.rm=T),align="center",fill=NA_real_)
  build_pi<-function(gl){gi<-which(flux_data[[gl]]%in%c(TRUE,1)); if(!length(gi)) return(NULL)
  gd<-tibble::tibble(date=flux_data$.cd[gi],nn=dplyr::if_else(isn[gi],nv[gi],0),cn=as.integer(isn[gi]),nd=dplyr::if_else(isd[gi],nv[gi],0),cd=as.integer(isd[gi]))%>%
    dplyr::group_by(date)%>%dplyr::summarise(nn=sum(nn,na.rm=T),cn=sum(cn,na.rm=T),nd=sum(nd,na.rm=T),cd=sum(cd,na.rm=T),.groups="drop")
  dd<-dplyr::left_join(dt,gd,by="date",suffix=c("","_g"))%>%dplyr::mutate(across(ends_with("_g"),~dplyr::coalesce(.,0)),an=pmax(nn-nn_g,0),acn=pmax(cn-cn_g,0L),ad=pmax(nd-nd_g,0),acd=pmax(cd-cd_g,0L))%>%dplyr::arrange(date)
  mn<-dplyr::if_else(rol(dd$acn)>0,rol(dd$an)/rol(dd$acn),NA_real_); md<-dplyr::if_else(rol(dd$acd)>0,rol(dd$ad)/rol(dd$acd),NA_real_)
  PR<-dplyr::if_else(is.finite(mn)&is.finite(md),mn-md,NA_real_); mx<-suppressWarnings(max(PR,na.rm=T))
  PN<-if(!is.finite(mx)||mx<=0) rep(NA_real_,length(PR)) else pmax(0,pmin(PR/mx,1))
  ph<-PN[match(flux_data$.cd,dd$date)]; ph[gi]<-NA_real_; ord<-order(flux_data$timestamp); po<-ph[ord]; xo<-as.numeric(flux_data$timestamp[ord])
  rl<-c(if(is.na(po[1]))2 else 1,if(is.na(po[length(po)]))2 else 1)
  if(sum(is.finite(po))>=2) po<-zoo::na.approx(po,x=xo,na.rm=FALSE,rule=rl)
  ph[ord]<-po; pf<-PN[match(flux_data$.cd,dd$date)]; pf[gi]<-ph[gi]; pf}
  for(pfx in gap_prefixes){gls<-names(flux_data)[grepl(paste0("^",pfx,"\\d+$"),names(flux_data))]; gls<-gls[order(as.integer(sub(pfx,"",gls)))]; if(!length(gls)) next
  log_msg("PI: ",length(gls)," columns for '",pfx,"'."); for(gl in gls){pv2<-build_pi(gl); if(!is.null(pv2)) flux_data[[paste0("PI_",gl)]]<-pv2}}
  flux_data$.cd<-NULL; flux_data}


# =============================================================================
# SECTION 5 — Regrowth period detection
# =============================================================================
build_regrowth_periods<-function(days_since,threshold=1.0){
  x<-suppressWarnings(as.numeric(days_since)); is_zero<-is.finite(x)&(x==0); x_prev<-dplyr::lag(x)
  cumsum(ifelse(is_zero&(is.na(x_prev)|!is.finite(x_prev)|(x_prev>threshold)),1L,0L))}


# SECTION 7 — Flux partitioning helpers
# =============================================================================
# Three functions shared by both ground-truth computation and gap prediction:
#
#  compute_reco_gpp_from_nee_orig()
#    Partitions NEE_orig into Reco_orig and GPP_orig per regrowth period.
#    These two columns serve as GROUND TRUTH when computing Reco/GPP metrics.
#    Called ONCE after regrowth periods are assigned.
#
#  derive_reco_gpp_from_filled_nee()
#    For each gap label Si, reconstructs a complete NEE series:
#      NEE_Si_with_predicted = NEE_Si  (observed outside gap)
#                             + NEE_{SIZE}_{model}_predicted  (at Si rows)
#    Then partitions this full series per regrowth period → Reco_Si, GPP_Si.
#    Writes only the Si rows into Reco_{SIZE}_{model}_predicted and
#    GPP_{SIZE}_{model}_predicted.
#    This ensures partitioning sees a complete, gap-free time series — the
#    same way standard post-hoc partitioning is applied in real pipelines.
#
# Night threshold : PPFD < 10 µmol m⁻² s⁻¹
# Day   threshold : PPFD > 10 µmol m⁻² s⁻¹
# Reco model      : Lloyd-Taylor Arrhenius  R10 estimated by OLS
# GPP  model      : Thornley non-rectangular hyperbola  fitted by BFGS

arrhenius_temp_scaling <- function(temp_celsius) {
  exp(309 * ((1 / (283.2 - 230)) - (1 / ((temp_celsius + 273.2) - 230))))
}

light_response_gpp <- function(params, ppfd_vec) {
  a <- params[1]; b <- params[2]; c <- params[3]
  num  <- a * ppfd_vec + b
  disc <- pmax(num^2 - 4 * c * (ppfd_vec * a * b), 0)
  (num - sqrt(disc)) / (2 * c)
}

gpp_ssr <- function(params, ppfd_vec, gpp_obs)
  sum((light_response_gpp(params, ppfd_vec) - gpp_obs)^2)

# Internal helper: fit R10 + GPP params and apply to a full period
.partition_one_period <- function(rows_in_period, nee_vec, ppfd_vec, temp_vec,
                                  is_night, is_day, gpp_init, ppfd_night_thr) {
  train_rows <- rows_in_period[is.finite(nee_vec[rows_in_period])]
  if (length(train_rows) < 10)
    return(list(reco = rep(NA_real_, length(rows_in_period)),
                gpp  = rep(NA_real_, length(rows_in_period))))
  
  night_rows <- train_rows[is_night[train_rows] & nee_vec[train_rows] > 0]
  R10 <- if (length(night_rows) > 0) {
    arr <- arrhenius_temp_scaling(temp_vec[night_rows])
    d   <- sum(arr^2, na.rm = TRUE)
    if (is.finite(d) && d > 0) sum(nee_vec[night_rows] * arr, na.rm = TRUE) / d
    else 1.0
  } else 1.0
  
  reco_period <- R10 * arrhenius_temp_scaling(temp_vec[rows_in_period])
  
  day_rows <- train_rows[is_day[train_rows] & nee_vec[train_rows] < 0]
  if (length(day_rows) > 10) {
    reco_day <- R10 * arrhenius_temp_scaling(temp_vec[day_rows])
    gpp_obs  <- reco_day - nee_vec[day_rows]
    opt <- tryCatch(
      optim(par = gpp_init, fn = gpp_ssr,
            ppfd_vec = ppfd_vec[day_rows], gpp_obs = gpp_obs, method = "BFGS"),
      error = function(e) list(par = gpp_init))
    gpp_params <- if (!is.null(opt$par)) opt$par else gpp_init
  } else {
    gpp_params <- gpp_init
  }
  
  gpp_raw <- light_response_gpp(gpp_params, ppfd_vec[rows_in_period])
  gpp_raw[!is_day[rows_in_period]]         <- 0
  gpp_raw[is.na(ppfd_vec[rows_in_period])] <- NA_real_
  
  list(reco = reco_period, gpp = gpp_raw)
}

# Internal helper: apply .partition_one_period across all regrowth periods
# and handle the fallback for period 0 (before first grazing event)
.partition_full_series <- function(nee_vec, ppfd_vec, temp_vec, regrowth,
                                   rg_ids, ppfd_night_thr, gpp_init) {
  n        <- length(nee_vec)
  is_night <- is.finite(ppfd_vec) & ppfd_vec < ppfd_night_thr
  is_day   <- is.finite(ppfd_vec) & ppfd_vec > ppfd_night_thr
  reco_out <- rep(NA_real_, n)
  gpp_out  <- rep(NA_real_, n)
  
  for (rg_id in rg_ids) {
    rip <- which(regrowth == rg_id)
    res <- .partition_one_period(rip, nee_vec, ppfd_vec, temp_vec,
                                 is_night, is_day, gpp_init, ppfd_night_thr)
    reco_out[rip] <- res$reco
    gpp_out[rip]  <- res$gpp
  }
  
  # Fallback: global fit for rows not covered by any named period
  fb <- which(is.na(reco_out) & is.finite(temp_vec))
  if (length(fb) > 0) {
    all_rows <- seq_len(n)
    res_fb <- .partition_one_period(all_rows, nee_vec, ppfd_vec, temp_vec,
                                    is_night, is_day, gpp_init, ppfd_night_thr)
    reco_out[fb] <- res_fb$reco[fb]
    gpp_out[fb]  <- res_fb$gpp[fb]
  }
  
  list(reco = reco_out, gpp = gpp_out)
}

# ---- Ground truth: Reco_orig and GPP_orig from NEE_orig --------------------
compute_reco_gpp_from_nee_orig <- function(flux_data, ppfd_night_thr = 10,
                                           gpp_init = c(0.08, 15, 0.5)) {
  if (!all(c("PPFD","Temp","regrowth_id","NEE_orig") %in% names(flux_data))) {
    log_msg("compute_reco_gpp_from_nee_orig: required columns missing — skipping.")
    return(flux_data)
  }
  res <- .partition_full_series(
    nee_vec  = as.numeric(flux_data$NEE_orig),
    ppfd_vec = as.numeric(flux_data$PPFD),
    temp_vec = as.numeric(flux_data$Temp),
    regrowth = flux_data$regrowth_id,
    rg_ids   = sort(unique(stats::na.omit(flux_data$regrowth_id))),
    ppfd_night_thr = ppfd_night_thr,
    gpp_init       = gpp_init
  )
  flux_data$Reco_orig <- res$reco
  flux_data$GPP_orig  <- res$gpp
  log_msg("Reco_orig and GPP_orig computed from NEE_orig (ground truth, per regrowth period).")
  flux_data
}

# ---- Per-gap: reconstruct NEE_Si_with_predicted, partition, write Si rows --
derive_reco_gpp_from_filled_nee <- function(flux_data,
                                            gap_size_cat,
                                            model_key,
                                            ppfd_night_thr = 10,
                                            gpp_init       = c(0.08, 15, 0.5)) {
  nee_pred_col  <- paste0("NEE_",  gap_size_cat, "_", model_key, "_predicted")
  reco_pred_col <- paste0("Reco_", gap_size_cat, "_", model_key, "_predicted")
  gpp_pred_col  <- paste0("GPP_",  gap_size_cat, "_", model_key, "_predicted")
  
  if (!nee_pred_col %in% names(flux_data)) {
    log_msg("derive_reco_gpp [", gap_size_cat, " / ", model_key,
            "]: NEE prediction column missing — skipping.")
    return(flux_data)
  }
  
  gap_labels <- names(flux_data)[grepl(paste0("^", gap_size_cat, "\\d+$"), names(flux_data))]
  gap_labels <- gap_labels[order(as.integer(sub(gap_size_cat, "", gap_labels)))]
  if (!length(gap_labels)) return(flux_data)
  
  ppfd_vec <- as.numeric(flux_data[["PPFD"]])
  temp_vec <- as.numeric(flux_data[["Temp"]])
  regrowth <- flux_data[["regrowth_id"]]
  rg_ids   <- sort(unique(stats::na.omit(regrowth)))
  nee_pred_cat <- as.numeric(flux_data[[nee_pred_col]])
  
  flux_data[[reco_pred_col]] <- NA_real_
  flux_data[[gpp_pred_col]]  <- NA_real_
  
  for (gap_lbl in gap_labels) {
    gap_rows <- which(flux_data[[gap_lbl]] %in% c(TRUE, 1))
    if (!length(gap_rows)) next
    
    # Step 1 — Reconstruct the complete NEE series for gap label Si:
    #   start from NEE_Si (observed outside gap, NA inside gap)
    #   fill the gap rows with model predictions from NEE_{SIZE}_{model}_predicted
    masked_col <- paste0("NEE_", gap_lbl)
    nee_complete <- if (masked_col %in% names(flux_data))
      as.numeric(flux_data[[masked_col]])          # NA at Si rows
    else
      flux_data$NEE_orig                           # fallback if column absent
    nee_complete[gap_rows] <- nee_pred_cat[gap_rows]   # fill Si with predictions
    
    # Step 2 — Partition the full reconstructed series per regrowth period
    res <- .partition_full_series(
      nee_vec  = nee_complete,
      ppfd_vec = ppfd_vec,
      temp_vec = temp_vec,
      regrowth = regrowth,
      rg_ids   = rg_ids,
      ppfd_night_thr = ppfd_night_thr,
      gpp_init       = gpp_init
    )
    
    # Step 3 — Write only Si rows into category-level prediction columns
    flux_data[[reco_pred_col]][gap_rows] <- res$reco[gap_rows]
    flux_data[[gpp_pred_col]][gap_rows]  <- res$gpp[gap_rows]
  }
  
  log_msg("Reco/GPP from reconstructed NEE_Si_with_predicted [",
          gap_size_cat, " / ", model_key, "].")
  flux_data
}

# =============================================================================
# SECTION 7 — MLP-specific helpers (scaling, train/val split, builder)
# =============================================================================

fit_feature_scaler <- function(df, train_rows, feature_cols) {
  mins <- vapply(feature_cols, function(c) suppressWarnings(min(as.numeric(df[[c]][train_rows]),na.rm=T)), numeric(1))
  maxs <- vapply(feature_cols, function(c) suppressWarnings(max(as.numeric(df[[c]][train_rows]),na.rm=T)), numeric(1))
  meds <- vapply(feature_cols, function(c) suppressWarnings(median(as.numeric(df[[c]][train_rows]),na.rm=T)), numeric(1))
  dens <- pmax(maxs-mins,0); dens[!is.finite(dens)|dens<=0] <- 1
  list(mins=mins, dens=dens, medians=meds, feature_cols=feature_cols)
}
apply_feature_scaler <- function(df, row_idx, scaler) {
  X <- as.matrix(df[row_idx, scaler$feature_cols, drop=FALSE])
  for(j in seq_along(scaler$feature_cols)){
    x <- as.numeric(X[,j]); x[!is.finite(x)] <- scaler$medians[j]
    X[,j] <- (x-scaler$mins[j])/scaler$dens[j]}
  X[X<0]<-0; X[X>1]<-1; storage.mode(X)<-"double"; colnames(X)<-scaler$feature_cols; X
}
fit_target_scaler   <- function(y){y<-as.numeric(y); mn<-suppressWarnings(min(y,na.rm=T)); mx<-suppressWarnings(max(y,na.rm=T)); if(!is.finite(mn)) mn<-0; if(!is.finite(mx)) mx<-1; d<-mx-mn; if(!is.finite(d)||d<=0) d<-1; list(min=mn,den=d)}
scale_target        <- function(y,s) (as.numeric(y)-s$min)/s$den
unscale_target      <- function(ys,s) as.numeric(ys)*s$den+s$min

stratified_train_val_split<-function(flux_data,obs_rows,val_frac=0.10,seed=42L){
  n_total<-length(obs_rows); if(!n_total) return(list(train_rows=integer(0),val_rows=integer(0)))
  n_val<-max(1L,round(val_frac*n_total)); ts<-flux_data$timestamp[obs_rows]
  if(inherits(ts,"Date")) ts<-as.POSIXct(ts)
  if(!inherits(ts,"POSIXt")) ts<-suppressWarnings(lubridate::parse_date_time(ts,orders=c("Ymd HMS","Ymd HM","Ymd")))
  seas<-timestamp_to_season(ts); mon<-as.integer(format(ts,"%m")); ngt<-suppressWarnings(as.integer(flux_data$night[obs_rows]))
  ok<-!(is.na(seas)|is.na(mon)|is.na(ngt)); idx_ok<-obs_rows[ok]; sk<-seas[ok]; mk<-mon[ok]; nk<-ngt[ok]
  set.seed(seed); n0<-floor(n_val/2L); n1<-n_val-n0; pick<-integer(0)
  rr<-function(pools,target){if(!length(pools)||target<=0) return(integer(0)); szs<-vapply(pools,length,integer(1)); q<-integer(length(pools)); left<-target; ord<-order(names(pools))
  while(left>0){prog<-FALSE; for(j in ord){if(left<=0) break; if(q[j]<szs[j]){q[j]<-q[j]+1L; left<-left-1L; prog<-TRUE}}; if(!prog) break}
  unlist(Map(function(ids,k) if(k>0L) sample(ids,k) else integer(0),pools,as.list(q)),use.names=FALSE)}
  if(any(ok)&&any(nk==0L)) pick<-c(pick,rr(split(idx_ok[nk==0L],paste(sk[nk==0L],mk[nk==0L],sep="|")),n0))
  if(any(ok)&&any(nk==1L)) pick<-c(pick,rr(split(idx_ok[nk==1L],paste(sk[nk==1L],mk[nk==1L],sep="|")),n1))
  sf<-n_val-length(pick); if(sf>0){rem<-setdiff(idx_ok,pick); if(length(rem)) pick<-c(pick,sample(rem,min(sf,length(rem))))}
  sf<-n_val-length(pick); if(sf>0){rem<-setdiff(obs_rows,pick); if(length(rem)) pick<-c(pick,sample(rem,min(sf,length(rem))))}
  list(train_rows=sort(setdiff(obs_rows,pick)),val_rows=sort(unique(pick)))}

build_mlp <- function(n_features, hidden=c(128L,64L,32L)) {
  keras_model_sequential() |>
    layer_dense(units=hidden[1], activation="relu", input_shape=n_features) |>
    layer_dense(units=hidden[2], activation="relu") |>
    layer_dense(units=hidden[3], activation="relu") |>
    layer_dense(units=1) |>
    compile(optimizer=optimizer_adam(), loss="mse", metrics=list("mae"))
}

extract_epoch_history <- function(h) {
  if(is.data.frame(h)&&nrow(h)){if(!"epoch"%in%names(h)) h$epoch<-seq_len(nrow(h)); return(h)}
  hh<-tryCatch(h$history,error=function(e)NULL); if(is.list(hh)&&length(hh)){
    out<-tibble::tibble(epoch=seq_len(max(vapply(hh,length,0L)))); for(k in names(hh)) out[[k]]<-suppressWarnings(as.numeric(hh[[k]])); return(out)}
  tryCatch({df2<-as.data.frame(h); if(!nrow(df2)) return(NULL); if(!"epoch"%in%names(df2)) df2$epoch<-seq_len(nrow(df2)); df2},error=function(e)NULL)
}

save_mlp_loss_curves <- function(history_df, gap_size_cat, gap_label) {
  if (is.null(history_df)||!nrow(history_df)) return(invisible(NULL))
  names(history_df) <- gsub("val_mean_absolute_error","val_mae",
                            gsub("mean_absolute_error","mae",names(history_df)))
  out_dir <- file.path(RESULTS_DIR, "training_loss_curves")
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  do_plot <- function(cols, ylabel, suffix) {
    present <- intersect(cols, names(history_df)); if(!length(present)) return(invisible(NULL))
    d <- history_df[,c("epoch",present),drop=FALSE] |> tidyr::pivot_longer(-epoch,names_to="dataset",values_to="value")
    p <- ggplot2::ggplot(d,ggplot2::aes(epoch,value,linetype=dataset)) +
      ggplot2::geom_line() + ggplot2::geom_point(size=0.6) +
      ggplot2::labs(title=paste0("MLP NEE ",gap_size_cat," '",gap_label,"' — ",ylabel),
                    subtitle="Scaled units [0,1]",x="Epoch",y=ylabel,linetype=NULL) +
      ggplot2::theme_minimal(base_size=11) +
      ggplot2::theme(plot.background=ggplot2::element_rect(fill="white"),legend.position="top")
    ggplot2::ggsave(file.path(out_dir,sprintf("MLP_NEE_%s_%s_%s.png",gap_size_cat,gap_label,suffix)),
                    p,width=7,height=4,dpi=200)
  }
  do_plot(c("loss","val_loss"),"MSE Loss (scaled)","loss")
  do_plot(c("mae","val_mae"),  "MAE (scaled)",     "mae")
  invisible(NULL)
}


# =============================================================================
# SECTION 8 — MLP cross-validation  (NEE only)
# =============================================================================

run_mlp_for_gap_size <- function(flux_data, gap_size_cat, feature_cols, val_frac=0.10) {
  log_msg("MLP [", gap_size_cat, "]: starting.")
  gap_labels <- names(flux_data)[grepl(paste0("^",gap_size_cat,"\\d+$"),names(flux_data))]
  gap_labels <- gap_labels[order(as.integer(sub(gap_size_cat,"",gap_labels)))]
  if (!length(gap_labels)) { log_msg("MLP [",gap_size_cat,"]: no gaps — skip."); return(flux_data) }
  log_msg("MLP [",gap_size_cat,"]: ",length(gap_labels)," gaps.")
  
  pred_col             <- paste0("NEE_", gap_size_cat, "_mlp_predicted")
  flux_data[[pred_col]] <- NA_real_
  
  for (gap_lbl in gap_labels) {
    masked_nee <- paste0("NEE_", gap_lbl)
    if (!masked_nee %in% names(flux_data)) next
    
    feats <- feature_cols
    if (PI_ENABLED) { pc<-paste0("PI_",gap_lbl); if(pc%in%names(flux_data)) feats<-unique(c(feats,pc)) }
    
    obs_rows <- which(is.finite(flux_data[[masked_nee]]))
    gap_rows <- which(flux_data[[gap_lbl]] %in% c(TRUE,1))
    log_msg("MLP [",gap_size_cat,"] ",gap_lbl,": n_obs=",length(obs_rows)," n_gap=",length(gap_rows))
    if (length(obs_rows)<10 || !length(gap_rows)) next
    
    spl     <- stratified_train_val_split(flux_data, obs_rows, val_frac, seed=42L)
    tr_rows <- spl$train_rows; val_rows <- spl$val_rows
    if (length(tr_rows)<3 || !length(val_rows)) next
    
    fs   <- fit_feature_scaler(flux_data, tr_rows, feats)
    Xtr  <- apply_feature_scaler(flux_data, tr_rows,  fs)
    Xval <- apply_feature_scaler(flux_data, val_rows, fs)
    Xgap <- apply_feature_scaler(flux_data, gap_rows, fs)
    
    ytr_raw  <- as.numeric(flux_data[[masked_nee]][tr_rows])
    yval_raw <- as.numeric(flux_data[[masked_nee]][val_rows])
    ts       <- fit_target_scaler(ytr_raw)
    ytr_s    <- scale_target(ytr_raw, ts)
    yval_s   <- scale_target(yval_raw, ts)
    
    mlp        <- build_mlp(ncol(Xtr))
    early_stop <- callback_early_stopping(monitor="val_loss", patience=10,
                                          restore_best_weights=TRUE)
    ep_log <- list()
    recorder <- callback_lambda(on_epoch_end=function(epoch,logs){
      row<-as.list(logs); row$epoch<-epoch+1L
      ep_log[[length(ep_log)+1L]] <<- as.data.frame(row,check.names=FALSE)})
    
    raw_hist <- fit(mlp, Xtr, ytr_s, epochs=20, batch_size=32,
                    validation_data=list(Xval,yval_s),
                    callbacks=list(early_stop,recorder), verbose=0, shuffle=TRUE)
    
    hist_df <- if (length(ep_log)) {
      h <- dplyr::bind_rows(ep_log)
      names(h) <- gsub("val_mean_absolute_error","val_mae",
                       gsub("mean_absolute_error","mae",names(h))); h
    } else extract_epoch_history(raw_hist)
    save_mlp_loss_curves(hist_df, gap_size_cat, gap_lbl)
    
    flux_data[[pred_col]][gap_rows] <- unscale_target(
      as.numeric(predict(mlp, Xgap)), ts)
  }
  
  log_msg("MLP [",gap_size_cat,"] complete — MAE = ",
          signif(calc_mae(flux_data$NEE_orig,flux_data[[pred_col]]),3),
          "  RMSE = ",signif(calc_rmse(flux_data$NEE_orig,flux_data[[pred_col]]),3),
          "  R² = ",  signif(calc_r2(  flux_data$NEE_orig,flux_data[[pred_col]]),3))
  flux_data
}


# =============================================================================
# SECTION 9 — Data preparation
# =============================================================================
flux_data <- df %>% dplyr::filter(is.finite(NEE_orig))
log_msg("Rows with finite NEE_orig retained: ", nrow(flux_data))

RECO_GPP_ENABLED <- all(c("PPFD","Temp") %in% names(flux_data))
if (!RECO_GPP_ENABLED) log_msg("WARNING: PPFD or Temp missing — Reco/GPP will be skipped.")

if (RECO_GPP_ENABLED) {
  if ("Grazing_days_since" %in% names(flux_data)) {
    flux_data <- flux_data %>%
      dplyr::mutate(regrowth_id=build_regrowth_periods(Grazing_days_since,threshold=1.0))
    log_msg("Regrowth periods: ",dplyr::n_distinct(flux_data$regrowth_id)," unique periods.")
  } else {
    flux_data$regrowth_id <- 1L
    log_msg("Grazing_days_since absent — using one global regrowth period.")
  }
  # Compute ground-truth Reco_orig / GPP_orig from NEE_orig per regrowth period
  flux_data <- compute_reco_gpp_from_nee_orig(flux_data)
}


# =============================================================================
# SECTION 10 — Construct artificial gaps
# =============================================================================
gap_cfg<-list(VL=list(base=30*48,tol=250),L=list(base=14*48,tol=120),
              M=list(base=7*48,tol=60),S=list(base=3*48,tol=10))
log_msg("Constructing artificial gaps ...")
cv_VL<-create_cv_gaps(flux_data,  gap_cfg$VL$base,gap_cfg$VL$tol,"VL")
cv_L <-create_cv_gaps(cv_VL$data, gap_cfg$L$base, gap_cfg$L$tol, "L")
cv_M <-create_cv_gaps(cv_L$data,  gap_cfg$M$base, gap_cfg$M$tol, "M")
cv_S <-create_cv_gaps(cv_M$data,  gap_cfg$S$base, gap_cfg$S$tol, "S")
flux_data<-cv_S$data
log_msg("Gaps ready — VL:",length(cv_VL$sizes)," L:",length(cv_L$sizes),
        " M:",length(cv_M$sizes)," S:",length(cv_S$sizes))


# =============================================================================
# SECTION 11 — Phytomass Index
# =============================================================================
if (PI_ENABLED) {
  flux_data <- add_phytomass_index(flux_data, gap_prefixes=c("VL","L","M","S"),
                                   window_days=21L)
  log_msg("PI columns added.")
}


# =============================================================================
# SECTION 12 — Keras / TensorFlow initialisation
# =============================================================================
log_msg("Loading Keras / TensorFlow ...")
library(reticulate); use_virtualenv("r-tensorflow", required=TRUE)
Sys.setenv(CUDA_VISIBLE_DEVICES="-1")
library(tensorflow); library(keras3)
set.seed(42)
import("numpy",    convert=TRUE)$random$seed(42L)
import("tensorflow",convert=TRUE)$random$set_seed(42L)
import("random",   convert=TRUE)$seed(42L)
log_msg("TensorFlow ready.")


# =============================================================================
# SECTION 13 — Run MLP cross-validation  (NEE only)
# =============================================================================
log_msg("=== Running MLP cross-validation — NEE ===")
for (cat in c("VL","L","M","S"))
  flux_data <- run_mlp_for_gap_size(flux_data, cat, predictors)
log_msg("All MLP NEE gap sizes complete.")


# =============================================================================
# SECTION 14 — Derive Reco and GPP from gap-filled NEE
# =============================================================================
if (RECO_GPP_ENABLED) {
  log_msg("=== Deriving Reco and GPP from MLP gap-filled NEE ===")
  for (cat in c("VL","L","M","S"))
    flux_data <- derive_reco_gpp_from_filled_nee(flux_data, cat, model_key="mlp")
  log_msg("Reco/GPP derivation complete.")
}


# =============================================================================
# SECTION 15 — Save outputs
# =============================================================================
dir.create(RESULTS_DIR, recursive=TRUE, showWarnings=FALSE)
saveRDS(flux_data, file.path(RESULTS_DIR,"df_cv_all_predictions.rds"))
log_msg("Saved: df_cv_all_predictions.rds")
writeLines(c(
  paste("Run finished    :", as.character(Sys.time())),
  paste("R version       :", R.version.string),
  paste("Source RDS      :", rds_name),
  paste("RESULTS_DIR     :", RESULTS_DIR),
  paste("Targets         :", if(RECO_GPP_ENABLED) "NEE (MLP), Reco+GPP (post-hoc)" else "NEE only"),
  paste("Predictors      :", paste(predictors,collapse=", "))
), file.path(RESULTS_DIR,"run_info.txt"))
log_msg("Saved: run_info.txt — MLP script finished.")
# ============================= end MLP_CV.R ===================================
