# =============================================================================
# XGBoost_CV.R — XGBoost Cross-Validation Gap-Filling for NEE
# =============================================================================
#
# WHAT THIS SCRIPT DOES
#   Trains an XGBoost model to gap-fill NEE only.
#   Reco and GPP are derived separately by Script 10 (flux_partitioning.R),
#   which reads the NEE predictions saved here.
#
# PRE-REQUISITES
#   This script expects df to be the output of Script 08 + Script 09, i.e.
#   JCi_cv.rds loaded by run_model.R.  That file already contains:
#     - Only rows where NEE_orig is observed  (Script 08 filtered NA rows)
#     - Gap flag columns  VL1…VLk, L1…Lk, M1…Mk, S1…Sk  (Script 08)
#     - Masked-NEE columns  NEE_VL1…, NEE_L1…, NEE_M1…, NEE_S1…  (Script 08)
#     - PI columns  PI_VL1…, PI_L1…, PI_M1…, PI_S1…  (Script 09)
#   Sections for NA filtering, gap construction, PI computation, and flux
#   partitioning are therefore no longer needed and have been removed.
#   The PI column for each gap label is selected automatically inside
#   run_xgb_for_gap_size(): for gap S1 it uses PI_S1, for M3 it uses PI_M3.
#
# HOW TO USE
#   Source from run_model.R, which must define:
#     df          — loaded from data/data_prepared/JCi_cv.rds
#     predictors  — character vector of base predictor column names
#     RESULTS_DIR — output directory path
#     rds_name    — source file name (written to run_info.txt)
#     PI_ENABLED  — logical; passed from run_model.R user settings
#
# OUTPUTS  (written to RESULTS_DIR)
#   df_cv_all_predictions.rds        ← NEE predictions only; Reco/GPP added by Script 10
#   boosting_evaluation_curves/  — XGB_NEE_{SIZE}_{LABEL}_rmse.png
#   progress.log | run_info.txt
#
# =============================================================================


# =============================================================================
# SECTION 1 — Reproducibility
# =============================================================================
set.seed(42)
Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")


# =============================================================================
# SECTION 2 — Progress logger
# =============================================================================
if (!exists("RESULTS_DIR",inherits=TRUE)||is.null(RESULTS_DIR))
  RESULTS_DIR <- file.path(tempdir(),"xgb_cv_fallback")
dir.create(RESULTS_DIR, recursive=TRUE, showWarnings=FALSE)
.run_log <- list()
log_msg <- function(...) {
  txt<-paste0("[",format(Sys.time(),"%Y-%m-%d %H:%M:%S"),"] ",paste(...,collapse=""))
  message(txt); .run_log <<- append(.run_log,list(txt))
  try(silent=TRUE,{writeLines(unlist(.run_log),file.path(RESULTS_DIR,"progress.log"))
    saveRDS(.run_log,file.path(RESULTS_DIR,"progress_log.rds"))}); invisible(txt)}
log_msg("XGBoost (NEE gap-fill) started.  RESULTS_DIR = ", RESULTS_DIR)


# =============================================================================
# SECTION 3 — Helper functions
# =============================================================================

# 3a  Metrics
.vp<-function(o,p) is.finite(as.numeric(o))&is.finite(as.numeric(p))
calc_mae  <-function(o,p){ok<-.vp(o,p); if(!any(ok)) return(NA_real_); mean(abs(as.numeric(p)[ok]-as.numeric(o)[ok]))}
calc_rmse <-function(o,p){ok<-.vp(o,p); if(!any(ok)) return(NA_real_); sqrt(mean((as.numeric(p)[ok]-as.numeric(o)[ok])^2))}
calc_r2   <-function(o,p){ok<-.vp(o,p); if(sum(ok)<2) return(NA_real_); ov<-as.numeric(o)[ok]; pv<-as.numeric(p)[ok]; den<-sum((ov-mean(ov))^2)*sum((pv-mean(pv))^2); if(den<=0) return(NA_real_); (sum((ov-mean(ov))*(pv-mean(pv)))^2)/den}

# 3b  Season labels
timestamp_to_season<-function(x){m<-as.integer(format(x,"%m")); dplyr::case_when(m%in%c(11,12,1)~"Winter",m%in%c(2,3,4)~"Spring",m%in%c(5,6,7)~"Summer",m%in%c(8,9,10)~"Autumn",TRUE~NA_character_)}

# 3c  Stratified train/val split
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

# 3d  XGBoost matrix builder
make_xgb_matrix<-function(df,row_idx,cols){X<-as.matrix(df[row_idx,cols,drop=FALSE]); for(j in seq_along(cols)) X[,j]<-suppressWarnings(as.numeric(X[,j])); storage.mode(X)<-"double"; colnames(X)<-cols; X}


# =============================================================================
# SECTION 4 — Phytomass Index flag
# =============================================================================
# PI_ENABLED is set by run_model.R and passed into this script via run_env.
# When TRUE, PI_{gap_label} columns (pre-computed by Script 09) are appended
# to the feature set for each gap: PI_S1 for gap S1, PI_M3 for gap M3, etc.
# When FALSE, no PI column is used and the base predictor set is unchanged.

if (!exists("PI_ENABLED")) PI_ENABLED <- FALSE   # safe fallback if sourced directly
log_msg("Phytomass Index: ", if (PI_ENABLED) "ENABLED" else "DISABLED")


# =============================================================================
# SECTION 5 — XGBoost boosting curve saver
# =============================================================================
xgb_params <- list(objective="reg:squarederror", eval_metric="rmse", nthread=1L)

save_xgb_boosting_curves<-function(trained_xgb, gap_size_cat, gap_label){
  elog<-trained_xgb$evaluation_log
  if(is.null(elog)||!nrow(elog)){log_msg("XGB: no eval log for [",gap_size_cat,"] ",gap_label); return(invisible(NULL))}
  ct<-grep("^train_.*rmse$",names(elog),value=TRUE); ce<-grep("^eval_.*rmse$",names(elog),value=TRUE)
  if(length(ct)!=1||length(ce)!=1){log_msg("XGB: RMSE cols not found — skip curve."); return(invisible(NULL))}
  d<-tibble::tibble(round=seq_len(nrow(elog)),`Training set`=as.numeric(elog[[ct]]),
                    `Validation set`=as.numeric(elog[[ce]]))%>%
    tidyr::pivot_longer(-round,names_to="dataset",values_to="rmse")
  out_dir<-file.path(RESULTS_DIR,"boosting_evaluation_curves")
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  p<-ggplot2::ggplot(d,ggplot2::aes(round,rmse,linetype=dataset))+
    ggplot2::geom_line()+ggplot2::geom_point(size=0.5)+
    ggplot2::labs(title=paste0("XGBoost NEE ",gap_size_cat," '",gap_label,"' — boosting curve"),
                  subtitle="Each round adds ONE new tree. RMSE = cumulative ensemble error.",
                  x="Boosting Round",y=expression(RMSE~(mu*mol~m^{-2}~s^{-1})),linetype=NULL)+
    ggplot2::theme_minimal(base_size=11)+
    ggplot2::theme(plot.background=ggplot2::element_rect(fill="white"),legend.position="top")
  ggplot2::ggsave(file.path(out_dir,sprintf("XGB_NEE_%s_%s_rmse.png",gap_size_cat,gap_label)),
                  p,width=7.5,height=4.5,dpi=200)
  log_msg("Saved boosting curve: XGB_NEE_",gap_size_cat,"_",gap_label,"_rmse.png")
  invisible(NULL)}


# =============================================================================
# SECTION 6 — XGBoost cross-validation  (NEE only)
# =============================================================================
# For each gap label gl within gap_size_cat:
#   - Training rows = rows where NEE_{gl} is finite (i.e. outside the gap)
#   - If PI_ENABLED, PI_{gl} is appended to feats — the gap-specific PI
#     from Script 09, computed from outside-gl rows only
#   - One XGBoost model trained per gap; predictions written back to gap rows only

run_xgb_for_gap_size<-function(flux_data, gap_size_cat, feature_cols,
                               val_frac=0.10, n_rounds=100L, seed=42L){
  log_msg("XGB [",gap_size_cat,"]: starting.")
  if(!length(feature_cols)) stop("feature_cols is empty.")
  absent<-setdiff(feature_cols,names(flux_data)); if(length(absent)) stop("Missing predictors: ",paste(absent,collapse=", "))
  gap_labels<-names(flux_data)[grepl(paste0("^",gap_size_cat,"\\d+$"),names(flux_data))]
  gap_labels<-gap_labels[order(as.integer(sub(gap_size_cat,"",gap_labels)))]
  if(!length(gap_labels)){log_msg("XGB [",gap_size_cat,"]: no gaps — skip."); return(flux_data)}
  log_msg("XGB [",gap_size_cat,"]: ",length(gap_labels)," gaps.")
  pred_col             <- paste0("NEE_",gap_size_cat,"_xgb_predicted")
  flux_data[[pred_col]] <- NA_real_
  for(gap_lbl in gap_labels){
    masked_nee<-paste0("NEE_",gap_lbl); if(!masked_nee%in%names(flux_data)) next
    feats<-feature_cols
    if(PI_ENABLED){pc<-paste0("PI_",gap_lbl); if(pc%in%names(flux_data)) feats<-unique(c(feats,pc))}
    obs_rows<-which(is.finite(flux_data[[masked_nee]])); gap_rows<-which(flux_data[[gap_lbl]]%in%c(TRUE,1))
    log_msg("XGB [",gap_size_cat,"] ",gap_lbl,": n_obs=",length(obs_rows)," n_gap=",length(gap_rows))
    if(length(obs_rows)<10||!length(gap_rows)) next
    spl<-stratified_train_val_split(flux_data,obs_rows,val_frac,seed); tr<-spl$train_rows; vr<-spl$val_rows
    if(length(tr)<5||!length(vr)) next
    dtrain<-xgboost::xgb.DMatrix(data=make_xgb_matrix(flux_data,tr,feats),  label=as.numeric(flux_data[[masked_nee]][tr]),missing=NA_real_)
    dval  <-xgboost::xgb.DMatrix(data=make_xgb_matrix(flux_data,vr,feats),  label=as.numeric(flux_data[[masked_nee]][vr]),missing=NA_real_)
    dtest <-xgboost::xgb.DMatrix(data=make_xgb_matrix(flux_data,gap_rows,feats),missing=NA_real_)
    set.seed(seed)
    xgb_model<-xgboost::xgb.train(params=xgb_params,data=dtrain,nrounds=as.integer(n_rounds),
                                  watchlist=list(train=dtrain,eval=dval),verbose=0)
    save_xgb_boosting_curves(xgb_model, gap_size_cat, gap_lbl)
    flux_data[[pred_col]][gap_rows]<-as.numeric(predict(xgb_model,dtest))}
  log_msg("XGB [",gap_size_cat,"] complete — MAE = ",
          signif(calc_mae(flux_data$NEE_orig,flux_data[[pred_col]]),3),
          "  RMSE = ",signif(calc_rmse(flux_data$NEE_orig,flux_data[[pred_col]]),3),
          "  R² = ",  signif(calc_r2(  flux_data$NEE_orig,flux_data[[pred_col]]),3))
  flux_data}


# =============================================================================
# SECTION 7 — Data preparation
# =============================================================================
# df arrives pre-filtered (no NA NEE_orig rows) and pre-built (all gap flag,
# masked-NEE, and PI columns already present from Scripts 08–09).

flux_data <- df
log_msg("Rows in pre-built CV data: ", nrow(flux_data))

# Report gap structure found in the pre-built data
for (pfx in c("VL","L","M","S")) {
  gcols <- names(flux_data)[grepl(paste0("^", pfx, "\\d+$"), names(flux_data))]
  if (length(gcols))
    log_msg("Gap columns found — ", pfx, ": ", length(gcols),
            " (", gcols[1], " … ", gcols[length(gcols)], ")")
}


# =============================================================================
# SECTION 8 — Run XGBoost cross-validation  (NEE only)
# =============================================================================
library(xgboost)
log_msg("=== Running XGBoost cross-validation — NEE ===")
for (cat in c("VL","L","M","S"))
  flux_data <- run_xgb_for_gap_size(flux_data, cat, predictors, n_rounds=100L)
log_msg("All XGBoost NEE gap sizes complete.")


# =============================================================================
# SECTION 9 — Save outputs
# =============================================================================
# Reco and GPP are NOT derived here.  Run Script 10 (flux_partitioning.R)
# after this script to add Reco_{SIZE}_xgb_predicted and GPP_{SIZE}_xgb_predicted
# columns to df_cv_all_predictions.rds.

dir.create(RESULTS_DIR,recursive=TRUE,showWarnings=FALSE)
saveRDS(flux_data, file.path(RESULTS_DIR,"df_cv_all_predictions.rds"))
log_msg("Saved: df_cv_all_predictions.rds  (NEE predictions only — run Script 10 for Reco/GPP)")

writeLines(c(
  paste("Run finished    :", as.character(Sys.time())),
  paste("R version       :", R.version.string),
  paste("Source RDS      :", rds_name),
  paste("RESULTS_DIR     :", RESULTS_DIR),
  paste("Targets         :", "NEE only — Reco/GPP via Script 10"),
  paste("Predictors      :", paste(predictors,collapse=", "))
), file.path(RESULTS_DIR,"run_info.txt"))
log_msg("Saved: run_info.txt — XGBoost script finished.")
# =========================== end XGBoost_CV.R =================================