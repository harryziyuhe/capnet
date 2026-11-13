# Regenerate every manuscript numeric table and application figure from saved fits.
out <- "manuscript/results"; tables <- "manuscript/tables"
dir.create(tables, showWarnings=FALSE, recursive=TRUE)
write_csv <- function(x, name) write.csv(as.data.frame(x),file.path(out,paste0(name,".csv")),row.names=FALSE)
tex_table <- function(headers, rows, name, align=NULL) {
  if (is.null(align)) align <- paste0("l",paste(rep("r",length(headers)-1),collapse=""))
  lines <- c(paste0("\\begin{tabular}{",align,"}"),"\\toprule",
    paste0(paste(headers,collapse=" & ")," \\\\"),"\\midrule",
    vapply(rows,function(row) paste0(paste(row,collapse=" & ")," \\\\"),character(1)),
    "\\bottomrule","\\end{tabular}")
  writeLines(lines,file.path(tables,paste0(name,".tex")))
}
fmt <- function(x,digits=3) {
  x[abs(x) < 0.5*10^(-digits)] <- 0
  formatC(x,format="f",digits=digits)
}
sim <- readRDS(file.path(out,"simulation.rds"))
write_csv(sim$raw,"simulation-replications"); write_csv(sim$summary,"simulation-summary")
for (evaluation in c("joint","walk")) {
  d <- sim$summary[sim$summary$scenario=="mixture" & sim$summary$evaluation==evaluation,]
  b <- d[d$gamma==0,]; rows <- list(c("0 (uncapped)",paste0("\\multicolumn{4}{c}{",fmt(b$mse),"}")))
  best <- min(d$mse[d$gamma>0])
  for (g in sim$config$gammas) {
    cells <- vapply(sim$config$kappas,function(k) {
      r <- d[d$gamma==g & !is.na(d$kappa) & d$kappa==k,]
      value <- fmt(r$mse)
      if (isTRUE(all.equal(r$mse,best,tolerance=1e-12))) paste0("\\textbf{",value,"}") else value
    },character(1))
    rows[[length(rows)+1]] <- c(as.character(g),cells)
  }
  tex_table(c("$\\gamma$",paste0("$\\kappa=",sim$config$kappas,"$")),rows,paste0("simulation-",evaluation))
}
finance <- readRDS(file.path(out,"finance.rds")); wiki <- readRDS(file.path(out,"wiki.rds"))
for (name in c("finance","wiki")) {
  fit <- if(name=="finance") finance else wiki; date_col <- if(name=="finance") "X" else "date"
  score <- paper_scores(fit,date_col)
  score <- score[order(match(score$method,c("Uncapped","Post-Fit Cap","Walk-forward Capnet")),score$gamma),]
  write_csv(score,paste0(name,"-scores")); write_csv(fit$diagnostics,paste0(name,"-diagnostics"))
  write_csv(fit$lambdas,paste0(name,"-lambdas"))
  digits <- if(name=="finance") 5 else 3
  rows <- lapply(seq_len(nrow(score)),function(i) c(
    if(score$method[i]=="Walk-forward Capnet") "Walk-forward \\pkg{capnet}" else score$method[i],
    if(score$method[i]=="Walk-forward Capnet") format(score$gamma[i],scientific=FALSE,trim=TRUE) else "--",
    fmt(score$oos_r2[i],digits)))
  tex_table(c("Method","$\\gamma$","OOS $R^2$"),rows,paste0(name,"-scores"))
  chosen <- if(name=="finance") .1 else {
    d <- score[score$method=="Walk-forward Capnet",]; d$gamma[which.max(d$oos_r2)]
  }
  plots <- paper_plots(fit,date_col,chosen,if(name=="finance")252 else 180,wiki=name=="wiki")
  prefix <- if(name=="finance") "" else "wiki_"
  if(name=="finance") ggplot2::ggsave(file.path("manuscript/figs",paste0(prefix,"oosrsq.png")),plots$rolling,width=9,height=5.6,dpi=180)
  ggplot2::ggsave(file.path("manuscript/figs",paste0(prefix,"max_contributions.png")),plots$contributions,width=9,height=5.6,dpi=180)
  # CSV of the exact plotted statistic provides a check independent of a PNG.
  write_csv(plots$maximum_data[,c("date","method","gamma","maximum")],paste0(name,"-maximum-contributions"))
  uncapped <- fit$predictions[fit$predictions$method=="Uncapped",]
  worst <- uncapped[which.max((uncapped$actual-uncapped$predicted)^2),]
  plot_summary <- data.frame(selected_gamma=chosen,rolling_window=if(name=="finance")252 else 180,
    min_rolling_uncapped=min(plots$rolling_data$r2_roll[plots$rolling_data$method=="Uncapped"],na.rm=TRUE),
    worst_date=as.character(worst[[date_col]]),n_predictions=nrow(uncapped))
  write_csv(plot_summary,paste0(name,"-plot-summary"))
  if(name=="wiki") {
    write_csv(paper_scores(fit,date_col,2024),"wiki-2024-scores")
    # Conditional Pearson dispersion is descriptive, not a Poisson assumption test.
    diagnostics <- fit$predictions |> dplyr::group_by(method,gamma) |>
      dplyr::summarise(pearson_mean=mean((actual-predicted)^2/pmax(predicted,1e-10)),
        mean_actual=mean(actual),mean_predicted=mean(predicted),.groups="drop")
    write_csv(diagnostics,"wiki-mean-diagnostics")
  }
  features <- if(name=="finance")fit$feature_cols else c("google","microsoft","tesla","nvidia","chatgpt","netflix","sony","elon_musk","sam_altman","youtube")
  snapshot <- paper_snapshot(fit,date_col,if(name=="finance")"2023-05-24" else "2023-03-14",chosen,features,if(name=="finance")1e4 else 1)
  write_csv(snapshot,paste0(name,"-snapshot"))
  snapshot <- rbind(snapshot,data.frame(feature="Sum (displayed)",as.list(colSums(snapshot[,-1])),check.names=FALSE))
  pretty <- c(google="Google",microsoft="Microsoft",tesla="Tesla",nvidia="Nvidia",chatgpt="ChatGPT",netflix="Netflix",sony="Sony",elon_musk="Elon Musk",sam_altman="Sam Altman",youtube="YouTube")
  rows <- lapply(seq_len(nrow(snapshot)),function(i) c(if(name=="wiki" && snapshot$feature[i] %in% names(pretty))pretty[[snapshot$feature[i]]] else snapshot$feature[i],fmt(as.numeric(snapshot[i,-1]))))
  tex_table(c("Feature","Uncapped","Post-fit cap","Walk-forward"),rows,paste0(name,"-snapshot"))
}
if (file.exists(file.path(out,"runtime.rds"))) {
runtime <- readRDS(file.path(out,"runtime.rds"))
write_csv(runtime,"runtime-replications")
runtime_summary <- runtime |> dplyr::group_by(model) |>
  dplyr::summarise(median=median(elapsed),min=min(elapsed),max=max(elapsed),repetitions=dplyr::n(),.groups="drop")
write_csv(runtime_summary,"runtime-summary")
tex_table(c("Model","Median (s)","Range (s)"),lapply(seq_len(nrow(runtime_summary)),function(i){
  r<-runtime_summary[i,];c(paste0("\\code{",gsub("_","\\\\_",r$model),"}"),fmt(r$median,2),paste0(fmt(r$min,2),"--",fmt(r$max,2)))
}),"runtime")
}
# Regenerate descriptive statistics from exactly the covariates used in fitting.
for(name in c("finance","wiki")) {
  d <- if(name=="finance")paper_finance_data() else paper_wiki_data()
  moments <- calc_moments(d[,setdiff(names(d),c("X","date","target"))])
  write_csv(moments,paste0(name,"-moments"))
  chosen_features <- if(name=="finance")finance$feature_cols else
    c("google","tesla","nvidia","chatgpt","sony","reddit","elon_musk","sam_altman","microsoft","youtube")
  displayed <- moments[match(chosen_features,moments$variable),]
  titles <- c(google="Google",tesla="Tesla",nvidia="Nvidia",chatgpt="ChatGPT",sony="Sony",
    reddit="Reddit",elon_musk="Elon Musk",sam_altman="Sam Altman",microsoft="Microsoft",youtube="YouTube")
  tex_table(c("Feature","Mean","Variance","Excess kurtosis"),
    lapply(seq_len(nrow(displayed)),function(i) {
      r<-displayed[i,]; c(if(name=="wiki")titles[[r$variable]] else r$variable,
        fmt(r$mean,2),fmt(r$variance,2),fmt(r$excess_kurtosis,2))
    }),paste0(name,"-moments"))
}
message("All manuscript tables and application figures exported.")

source("tools/write_paper_manifest.R")
