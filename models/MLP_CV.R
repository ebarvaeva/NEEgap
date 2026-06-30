# MLP_CV.R — MLP Cross-Validation Gap-Filling for NEE
#
# WHAT THIS SCRIPT DOES
#   Trains a Multi-Layer Perceptron (keras3 / TensorFlow) to gap-fill NEE.
#
# PRE-REQUISITES
#   This script expects df to be the output of Script 08 + Script 09, i.e.
#   JCi_cv.rds loaded by run_model.R.  That file already contains:
#     - Only rows where NEE_orig is observed  (Script 08 filtered NA rows)
#     - Gap flag columns  VL1…VLk, L1…Lk, M1…Mk, S1…Sk  (Script 08)
#     - Masked-NEE columns  NEE_VL1…, NEE_L1…, NEE_M1…, NEE_S1…  (Script 08)
#     - PI columns  PI_VL1…, PI_L1…, PI_M1…, PI_S1…  (Script 09)
#   The PI column for each gap label is selected automatically inside
#   run_mlp_for_gap_size(): for gap S1 it uses PI_S1, for M3 it uses PI_M3.
#
# HOW TO USE
#   Source from run_model.R, which must define:
#     df          — loaded from data/data_prepared/JCi_cv.rds
#     predictors  — character vector of base predictor column names
#     RESULTS_DIR — output directory path
#     PI_ENABLED  — logical; passed from run_model.R user settings
#
# OUTPUTS  (written to RESULTS_DIR)
#   df_cv_all_predictions.rds
#   training_loss_curves/  — MLP_NEE_{SIZE}_{LABEL}_{loss|mae}.png


# Reproducibility
set.seed(42)
Sys.setenv(PYTHONHASHSEED="0", CUDA_VISIBLE_DEVICES="-1",
           OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
           TF_NUM_INTRAOP_THREADS="1", TF_INTEROP_THREADS="1",
           TF_DETERMINISTIC_OPS="1")


# Output directory (provided by run_model.R; temp-dir fallback if sourced directly)
if (!exists("RESULTS_DIR", inherits=TRUE) || is.null(RESULTS_DIR))
  RESULTS_DIR <- file.path(tempdir(), "mlp_cv_fallback")
dir.create(RESULTS_DIR, recursive=TRUE, showWarnings=FALSE)


# Season labels
timestamp_to_season<-function(x){m<-as.integer(format(x,"%m")); dplyr::case_when(m%in%c(11,12,1)~"Winter",m%in%c(2,3,4)~"Spring",m%in%c(5,6,7)~"Summer",m%in%c(8,9,10)~"Autumn",TRUE~NA_character_)}


# Phytomass Index flag (from run_model.R); if TRUE, the gap-specific PI_{label} column is appended per gap
if (!exists("PI_ENABLED")) PI_ENABLED <- FALSE   # safe fallback if sourced directly


# MLP-specific helpers (scaling, train/val split, builder)

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


# MLP cross-validation: one model per gap (train on rows outside the gap, predict the gap rows; PI_{gap} added when PI_ENABLED)

run_mlp_for_gap_size <- function(flux_data, gap_size_cat, feature_cols, val_frac=0.10) {
  gap_labels <- names(flux_data)[grepl(paste0("^",gap_size_cat,"\\d+$"),names(flux_data))]
  gap_labels <- gap_labels[order(as.integer(sub(gap_size_cat,"",gap_labels)))]
  if (!length(gap_labels)) return(flux_data)
  
  pred_col             <- paste0("NEE_", gap_size_cat, "_mlp_predicted")
  flux_data[[pred_col]] <- NA_real_
  
  for (gap_lbl in gap_labels) {
    masked_nee <- paste0("NEE_", gap_lbl)
    if (!masked_nee %in% names(flux_data)) next
    
    feats <- feature_cols
    if (PI_ENABLED) { pc<-paste0("PI_",gap_lbl); if(pc%in%names(flux_data)) feats<-unique(c(feats,pc)) }
    
    obs_rows <- which(is.finite(flux_data[[masked_nee]]))
    gap_rows <- which(flux_data[[gap_lbl]] %in% c(TRUE,1))
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
  
  flux_data
}


# Data preparation: df is pre-filtered and pre-built by Scripts 08-09
flux_data <- df


# Keras / TensorFlow initialisation
library(reticulate); use_virtualenv("r-tensorflow", required=TRUE)
Sys.setenv(CUDA_VISIBLE_DEVICES="-1")
library(tensorflow); library(keras3)
set.seed(42)
import("numpy",    convert=TRUE)$random$seed(42L)
import("tensorflow",convert=TRUE)$random$set_seed(42L)
import("random",   convert=TRUE)$seed(42L)


# Run MLP cross-validation
for (cat in c("VL","L","M","S"))
  flux_data <- run_mlp_for_gap_size(flux_data, cat, predictors)


# Save outputs
saveRDS(flux_data, file.path(RESULTS_DIR,"df_cv_all_predictions.rds"))
