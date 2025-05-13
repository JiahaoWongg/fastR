open_r_porfile <- function(){
    file_path = paste0(.libPaths(), "/base/R/Rprofile")
    shell.exec(file_path)
}

rm_all <- function(x) rm(list = ls(envir = parent.frame()), envir = parent.frame())

READ <- function(file, ...){
    format = subString(file, 1, "\\.", rev = TRUE)
    if(format == "csv"){
        mydata = read.csv(file, ...)
    } else if(format == "txt"){
        mydata = read.table(file, ...)
    } else if(format == "RData"){
        mydata = sapply(load2(file, ...), function(x) get(x), simplify = FALSE)
        if(length(mydata) == 1){
            mydata = mydata[[1]]
        } else {
            attach(mydata)            
        }
    } else if(format == "RDS"){
        mydata = readRDS(file, ...)
    } else if(format == "xlsx"){
        mydata = openxlsx::read.xlsx(file, ...)
    } else if(format == "xls"){
        mydata = suppressMessages(readxl::read_excel(file, ...))
    } else if(format == "bedgraph"){
        mydata = read.faster(file = file, header = FALSE, sep = "\t", skip = 1, ...)
    } else if(format == "gmt"){
        geneSetDB = readLines(file)
        geneSetDB = strsplit(geneSetDB, "\t")
        names(geneSetDB) = sapply(geneSetDB, "[", 1)
        geneSetDB = lapply(geneSetDB, "[", -1:-2)
        geneSetDB = lapply(geneSetDB, function(x) x[which(x != "")])
        mydata = geneSetDB
    } else {
        cat("Error: Unrecognized file type:", format, "\n")
        invisible(return())
    }
    return(mydata)
}

case_when <- function(...) dplyr::case_when(...)

gc2 <- function() invisible(gc())

split_10x_3file <- function(input_dir, gene_name = "features"){
    files = lf("rf", input_dir, "barcodes")
    SIDs = subString(basename(files), 1, "_barcodes")
    single_dir = paste0(input_dir, "_single")
    for(SID in SIDs){
        xsingle_dir = paste0(single_dir, "/", SID)
        mkdir(xsingle_dir)
        features_file = ifelse(gene_name == "features", "features.tsv.gz", "genes.tsv.gz")
        from = paste0(input_dir, "/", SID, "_", c("barcodes.tsv.gz", features_file, "matrix.mtx.gz"))
        to = paste0(xsingle_dir, "/",  c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"))
        invisible(file.copy(from, to))
    }
}

convertId <- function (genes = NULL, from, to, OrgDb = "org.Hs.eg.db", drop = TRUE, remove = TRUE){
    if (is.null(genes)) {
        cat("\nAvailable gene types:\n")
        cat("\n\tACCNUM, ALIAS, ENSEMBL, ENSEMBLPROT, ENSEMBLTRANS, ENTREZID, ENZYME\n")
        cat("\n\tEVIDENCE, EVIDENCEALL, GENENAME, GO, GOALL, IPI, MAP, OMIM, ONTOLOGY\n")
        cat("\n\tONTOLOGYALL, PATH, PFAM, PMID, PROSITE, REFSEQ, SYMBOL, UCSCKG, UNIGENE, UNIPROT\n\n")
    }
    else {
        suppressPackageStartupMessages(library("clusterProfiler"))
        res = suppressMessages(suppressWarnings(bitr(genes, fromType = from, toType = to, OrgDb = OrgDb)))
        if(length(to) == 1 & drop){
            if(remove){
               res = unique(res[, 2])
               res = res[res != ""] 
            } else {
               res = res[, 2][match(genes, res[, 1])]
           }
        } else {
           res = res[res != "", ] 
        }
        return(res)
    }
}

unload_all_pkgs <- function(){
    loaded_packages = search()
    loaded_packages = loaded_packages[grepl("^package:", loaded_packages)]
    loaded_packages = sub("^package:", "", loaded_packages)
    base_packages = c("base", "stats", "utils", "grDevices", "graphics", "methods", "datasets")
    packages_to_unload = setdiff(loaded_packages, base_packages)

    for(pkg in packages_to_unload){
        try(detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE, force = TRUE), silent = TRUE)
    }
}

unload_pkgs <- function(...){
    pkgs = as.character(substitute(list(...)))[-1]
    for(pkg in pkgs){
        try(detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE, force = TRUE), silent = TRUE)
        message("已卸载包: ", pkg)
    }
}

p_value_to_stars <- function(p_value) {
  if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("")
  }
}

lf <- function (type = "x", folder, pattern = "", remove = NULL, include.dirs = FALSE,
                ignore.case = FALSE, list.files = FALSE, ...){
    if (!type %in% c("rf", "f", "r", "x"))
        stop("\nOnly support below types: \n
            \t\t\t\t- rf: recursive = TRUE, full.names = TRUE \n
            \t\t\t\t-  f: recursive = TRUE, full.names = FALSE \n
            \t\t\t\t-  r: recursive = FALSE, full.names = TRUE \n
            \t\t\t\t-  x: recursive = FALSE, full.names = FALSE\n")
    res = switch(type, 
                 rf = list.files(folder, pattern, full.names = TRUE, recursive = TRUE, 
                                 include.dirs = include.dirs, ignore.case = ignore.case, ...),
                 f  = list.files(folder, pattern, full.names = TRUE, recursive = FALSE,
                                 include.dirs = include.dirs, ignore.case = ignore.case, ...),
                 r  = list.files(folder, pattern, full.names = FALSE, recursive = TRUE,
                                 include.dirs = include.dirs, ignore.case = ignore.case, ...),
                 x  = list.files(folder, pattern, full.names = FALSE, recursive = FALSE,
                                 include.dirs = include.dirs, ignore.case = ignore.case, ...))
    if (!is.null(remove)) {
        for (x in remove) {
            res = res[!grepl(x, res)]
        }
    }
    return(res)
}

# 生存分析
run_unicox_survival <- function(tdata, genes, out_name, pvalue_cut = 0.05){  
    loadp(survival, survminer)
    colnames(tdata)[1:2] = c("futime", "fustat")
    genes_rm = setdiff(genes, colnames(tdata))    
    if(length(genes_rm) > 0){
        cat2("Not fonud genes", paste0(genes_rm, collapse = " "), "in data!\n")
    }
    tdata$futime = convert_time_to_year(tdata$futime)
    genes_inter = intersect(genes, colnames(tdata))
    for(gene in genes_inter){
        cox = coxph(Surv(futime, fustat) ~ tdata[, gene], data = tdata)
        coxSummary = summary(cox)
        coxP = coxSummary$coefficients[,"Pr(>|z|)"]
        HR = coxSummary$conf.int[,"exp(coef)"]
        HR.95L = coxSummary$conf.int[,"lower .95"]
        HR.95H = coxSummary$conf.int[,"upper .95"]
        pvalue = coxSummary$coefficients[,"Pr(>|z|)"]
        xresT = data.frame(gene, HR, HR.95L, HR.95H, pvalue)
        if(gene == candidate_genes2[1]){
            resT = xresT
        } else {
            resT = rbind(resT, xresT)
        }
    }
    resT_sig = resT[resT$pvalue < pvalue_cut, ]
    resT_sig = resT_sig[order(resT_sig$pvalue), ]
    rownames(resT_sig) = NULL
    write.csv(resT_sig, file = "03_TCGA_uniCox.csv", row.names = F)
    return(resT_sig)
}

convert_time_to_year <- function(times, from = NULL){
    times = as.numeric(times)
    if(is.null(from)){
        if(max(times) > 50){
            if(max(times) > 500){
                res = signif(times / 365, 6)
            } else{
                res = signif(times / 12, 6)
            }
        } else {
            res = signif(times, 6)
        }
    } else {
        if(from == "day"){
            res = signif(times / 365, 6)
        } else if(from == "month"){
            res = signif(times / 12, 6)
        } else {
            res = signif(times, 6)
        }
    }
}

plot_survival_forest <- function(pdata, 
    out_name = "Survival_forest", 
    cols = c("#F38329", "#1C79B7"), h = NULL, top = 30
    ){
    pdata = pdata[1:ifelse(nrow(pdata) > 30, 30, nrow(pdata)), ]
    gene = pdata[, 1]
    HR = sprintf("%.3f",pdata$"HR")
    HR_Low = sprintf("%.3f",pdata$"HR.95L")
    HR_High = sprintf("%.3f",pdata$"HR.95H")
    Hazard.ratio = paste0(HR, "(",HR_Low, "-", HR_High,")")
    pvalue = signif(pdata$pvalue, 4)
    
    if(is.null(h)) height = nrow(pdata) / 12 + 5
    pdf(paste0(out_name, ".pdf"), width = 7.5, height = height)
        n <- nrow(pdata)
        n_row <- n + 1
        ylim <- c(1, n_row)
        layout(matrix(c(1, 2), nc = 2), width = c(3, 2.5))
        
        xlim = c(0,3)
        par(mar=c(4,2.5,2,1))
        plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
        text.cex=0.8
        text(0,n:1,gene,adj=0,cex=text.cex)
        text(1.5-0.5*0.2,n:1,pvalue,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
        text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
        
        par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
        xlim = c(0,max(as.numeric(HR_Low),as.numeric(HR_High)))
        plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
        arrows(as.numeric(HR_Low),n:1,as.numeric(HR_High),n:1,angle=90,code=3,length=0.05,col="black",lwd=2.5)
        abline(v=1,col="black",lty=2,lwd=2)
        boxcolor = ifelse(as.numeric(HR) > 1, cols[1], cols[2])
        points(as.numeric(HR), n:1, pch = 15, col = boxcolor, cex=1.6)
        axis(1)
    dev.off()
}

run_survival_lasso <- function(tdata, genes, out_name = "Survial_LASSO", seed = 1234){
    loadp(glmnet, survival)
    tdata = tdata[!is.na(tdata$futime) & !is.na(tdata$fustat), ]
    x = as.matrix(tdata[, genes])
    y = data.matrix(Surv(tdata$futime, tdata$fustat))
    set.seed(seed)
    fit = glmnet(x, y, family = "cox", alpha = 1, nlambda = 100)

    set.seed(seed)
    cvfit = cv.glmnet(x, y, family = "cox", alpha = 1,nfolds = 10)
    best_lamda = signif(log(cvfit$lambda.min), 6)

    pdata = broom::tidy(fit)
    pdata = pdata[pdata$term != "(Intercept)", ]
    if(length(unique(pdata$term)) > 10) legend =  FALSE
    text_y = calculate_percentile_point(min(pdata$estimate), max(pdata$estimate), 4)
    p <- ggplot(pdata, aes(log(lambda), estimate, group = term, color = term))
    p <- p + geom_line(size = 0.5)
    p <- p + geom_vline(xintercept = best_lamda, color = "grey60", alpha = 0.8, linetype = 2)
    p <- p + geom_hline(yintercept = 0)
    p <- p + ylab("Coefficients")
    p <- p + scale_x_continuous(expand = c(0.02, 0.02))
    p <- p + setText(20, theme = "bw")
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    if(legend){
        p <- p + theme(legend.position = "top")
    } else {
        p <- p + theme(legend.position = "none")
    }
    p <- p + annotate("text", x = best_lamda, y = text_y, label = paste("Optimal Lambda =", signif(cvfit$lambda.min, 4)))
    save.graph(p, file = paste0(out_name, "_01_lambda"), 7, 4.5)

    pdata_cv = broom::tidy(cvfit)
    text_y = calculate_percentile_point(min(pdata_cv$estimate), max(pdata_cv$estimate), 4)
    p <- ggplot()
    p <- p + geom_point(data = pdata_cv, aes(log(lambda), estimate))
    p <- p + geom_errorbar(data = pdata_cv, aes(x = log(lambda), ymin = conf.low, ymax = conf.high))
    p <- p + ylab("Coefficients")
    p <- p + scale_x_continuous(expand = c(0.02, 0.02))
    p <- p + setText(20, theme = "bw")
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + annotate("text", x = best_lamda, y = text_y, label = paste("Optimal Lambda =", signif(cvfit$lambda.min, 4)))
    p <- p + geom_vline(xintercept = best_lamda, color = "grey60", alpha = 0.8, linetype = 2)   
    save.graph(p, file = paste0(out_name, "_02_cvfit"), 7, 4.5)

    coef = coef(fit, s = cvfit$lambda.min)
    cat2("lambda.min:", cvfit$lambda.min, "\n")

    idx = coef[, 1] != 0
    actCoef = coef[idx]
    lasso_genes = row.names(coef)[idx]
    cat2("Get", length(lasso_genes), "genes\n")

    lasso_res = data.frame(gene = lasso_genes, coef = actCoef)
    write.csv(lasso_res, file = paste0(out_name, "_03_result.csv"), row.names = FALSE)
    return(lasso_res)
}

run_multicox_survival <- function(tdata, genes, out_name = "multiCox", pvalue_cut = 0.05){
    loadp(survival, survminer)
    colnames(tdata)[1:2] = c("futime", "fustat")
    genes_rm = setdiff(genes, colnames(tdata))    
    if(length(genes_rm) > 0){
        cat2("Not fonud genes", paste0(genes_rm, collapse = " "), "in data!\n")
    }
    tdata$futime = convert_time_to_year(tdata$futime)
    genes_inter = intersect(genes, colnames(tdata))

    data_multi = tdata[, genes]
    multiCox = coxph(Surv(tdata$futime, tdata$fustat) ~ ., data = data_multi)
    multi_summary = summary(multiCox)
    pvalue = multi_summary$coefficients[,"Pr(>|z|)"]
    HR=multi_summary$conf.int[,"exp(coef)"]
    HR.95L=multi_summary$conf.int[,"lower .95"]
    HR.95H=multi_summary$conf.int[,"upper .95"]
    multi_res = data.frame(gene = genes, HR, HR.95L, HR.95H, pvalue)
    write.csv(multi_res, file = paste0(out_name, "_multiCox.csv"), row.names = FALSE)

    idx = multi_res$pvalue < pvalue_cut
    multi_sig = data.frame(gene = multi_res$gene, coef = multi_res$HR - 1)[idx, ]
    rownames(multi_sig) = NULL
    return(multi_sig)
}

plot_coefficient_barplot <- function(multi_res, pvalue_cut = 0.05, out_name = "Coefficient_barplot"){
    loadp(ggpubr)
    multi_res$group = ifelse(multi_res$coef > 0, "pos", "neg")
    p <- ggbarplot(multi_res, x="gene", y="coef", fill = "group", color = "white",
                   palette = c("#00AFBB", "#E7B800"), sort.val = "desc", sort.by.groups = FALSE,
                   x.text.angle=90, ylab = "Coefficient", xlab = "Genes", rotate=TRUE)
    p <- p + setText(20, theme = "bw", y.title.vjust = 2) + theme(legend.position = "none")
    save.graph(p, file = paste0(out_name, "_Coefficient_barplot"), 5, 6)
    write.csv(multi_res, )
    return(p)
}

# build_survial_model <- function(fit_res, pValue_cut = 0.05){}

# assessment_survial_model(tdata, genes, coefs, out_name = "00_TCGA")
assessment_survial_model <- function(
    tdata, genes, coefs, 
    cut_method = "best",
    out_name = "Survial",
    label.x = 4,
    label.y = 0.95,
    km_cols = c("#cd0000", "#118aff"),
    roc_cols = c("#1C79B7", "#F38329", "#2DA248"),
    time_point = c(1, 3, 5)
    ){
    loadp(survival, survminer, ggplot2, pheatmap, patchwork, timeROC)
    colnames(tdata)[1:2] = c("futime", "fustat")
    tdata = tdata[!is.na(tdata$futime) & !is.na(tdata$fustat), ]
    tdata$futime = convert_time_to_year(tdata$futime)
    genes_rm = setdiff(genes, colnames(tdata))    
    if(length(genes_rm) > 0){
        cat2("Not fonud genes", paste0(genes_rm, collapse = " "), "in data!\n")
    }

    genes_inter = intersect(genes, colnames(tdata))
    coefs_inter = coefs[genes %in% genes_inter]
    gene_exp = tdata[, genes_inter]
    my_fun = function(x) {crossprod(as.numeric(x), coefs_inter)}
    tdata$risk_score = apply(gene_exp, 1, my_fun)

    if(cut_method == "best"){
        sur.cut = surv_cutpoint(tdata, "futime", "fustat", "risk_score")
        sur.group = surv_categorize(sur.cut)
        tdata$risk_group = ifelse(sur.group$risk_score == "low", "Low risk", "High risk")

        # Cutoff detect
        p <- plot(sur.cut, "risk_score", palette = "npg")
        p1 <- p$risk_score[[1]]
        p1 <- p1 + gg_style(15)
        p1 <- p1 + suppressMessages(scale_color_manual(values = c("#ff3c41", "#0ebeff"), name = "Group", breaks = c("high", "low"), labels = c("High risk", "Low risk")))
        save.graph(p1, file = paste0(out_name, "_01_Cutoff_detect_Dot"), 6, 4)

        p2 <- p$risk_score[[2]]
        p2 <- p2 + gg_style(15)
        p2 <- p2 + scale_fill_manual(values = c("#ff3c41", "#0ebeff"), name = "Group", breaks = c("high", "low"), labels = c("High risk", "Low risk"))
        p2 <- p2 + guides(color = "none")
        save.graph(p2, file = paste0(out_name, "_01_Cutoff_detect_Line"), 6, 4)
    } else {
        tdata$risk_group = ifelse(tdata$risk_score < median(tdata$risk_score), "Low risk", "High risk")
    }
    write.csv(tdata, file = paste0(out_name, "_01_Risk_score.csv"))


    fit = survfit(Surv(futime, fustat) ~ risk_group, data = tdata)
    p <- ggsurvplot(fit, data = tdata, pval = FALSE, risk.table = TRUE,
                    legend.title = "Risk",
                    legend.labs = c("High risk", "Low risk"),
                    xlab = " Time (Years)", palette = km_cols)
    p$plot <- p$plot + theme(legend.title = element_text(size = 15))

    res_cox <- coxph(Surv(futime, fustat) ~ risk_group, data = tdata)
    text_HR = paste("HR =", round(1 / summary(res_cox)$conf.int[1],2))
    text_CI = paste0("95CI%:(", round(1 / summary(res_cox)$conf.int[4],2), "-", round(1 / summary(res_cox)$conf.int[3],2),")")
    text_P  = paste(" =", formatC(summary(res_cox)$coef[5], format = "e", digits = 2))
    p$plot <- p$plot + annotate("text", x = label.x, y = label.y, size = 6, hjust = 0, label = bquote(italic("P") * .(text_P)))
    p$plot <- p$plot + annotate("text", x = label.x, y = label.y - 0.1, size = 4, hjust = 0, label = paste0(text_HR, ", ", text_CI))
    p$plot <- p$plot + setText(20, y.title.vjust = 3) + theme(legend.title = element_text(size = 15), legend.position = "top")

    p$table <- p$table + setText(20, y.title.vjust = 3) + theme(plot.title = element_text(hjust = 0, face = "plain"))
    p$table <- p$table + theme(axis.text.y = element_text(color = rev(km_cols)))

    loadp(patchwork)
    layout <- "\nA\nA\nA\nB"
    p_merge <- p$plot / p$table + plot_layout(design = layout)
    save.graph(p_merge, file = paste0(out_name, "_02_KM_curve"), 8, 6)

    cat2(paste0("Pvalue", text_P, "\n"))
    cat2(text_HR, text_CI, "\n")

    # ROC
    ROC_rt = timeROC(tdata$futime, tdata$fustat, marker = tdata$risk_score, 
                   cause = 1, weighting = 'aalen', times = time_point, ROC = TRUE)
    cat2("AUC at", time_point[1], "years: ", sprintf("%.03f",ROC_rt$AUC[1]), "\n")
    cat2("AUC at", time_point[2], "years: ", sprintf("%.03f",ROC_rt$AUC[2]), "\n")
    cat2("AUC at", time_point[3], "years: ", sprintf("%.03f",ROC_rt$AUC[3]), "\n")
    pdf(file = paste0(out_name, "_03_ROC.pdf"), 5, 5.7)
        plot(ROC_rt, time = time_point[1], col = roc_cols[1], title = FALSE, lwd = 2)
        plot(ROC_rt, time = time_point[2], col = roc_cols[2], add = TRUE, title = FALSE, lwd =2)
        plot(ROC_rt, time = time_point[3], col = roc_cols[3], add = TRUE, title = FALSE, lwd =2)
        legend("bottomright", c(paste("AUC at", time_point[1], "years: ", sprintf("%.03f",ROC_rt$AUC[1])),
                                paste("AUC at", time_point[2], "years: ", sprintf("%.03f",ROC_rt$AUC[2])),
                                paste("AUC at", time_point[3], "years: ", sprintf("%.03f",ROC_rt$AUC[3]))),
               col = roc_cols, lwd = 2, bty = 'n')
    dev.off()

    # 散点图
    pdata = tdata
    pdata$Sample = pdata$SID
    pdata = pdata[order(pdata$risk_score), ]

    risk_class = pdata$risk_group
    low_length = length(risk_class[risk_class=="Low risk"])
    high_length = length(risk_class[risk_class=="High risk"])
    line = pdata$risk_score

    pdf(file = paste0(out_name, "_04_Risk_score_curve.pdf"), 8, 5)
        plot(line, type = "p", pch = 20,
             xlab = "Patients (increasing risk socre)", 
             ylab = "Risk score",
             col = c(rep("#0ebeff",low_length),
                     rep("#ff3c41",high_length)))
        legend("topleft", pch = 21, 
               legend = c("High risk","Low risk"),
               col = c("#ff3c41", "#0ebeff"), 
               pt.bg = c("#ff3c41", "#0ebeff"))
        abline(h = sur.cut$cutpoint$cutpoint,v=low_length,lty=2)
    dev.off()

    #绘制生存状态图
    color = as.vector(tdata$fustat)
    color[color == 1] ="#ff3c41"
    color[color == 0] ="#0ebeff"
    pdf(file = paste0(out_name, "_05_Survival_stat.pdf"), 8, 5)
        plot(tdata$futime,
             pch=19,
             xlab="Patients (increasing risk socre)",
             ylab="Survival time (Years)",
             col=color)
        legend("topleft", legend = c("Dead","Alive"), col = c("#ff3c41","#0ebeff"), pch = 21, pt.bg = c("#F39B7FFF","#91D1C2FF"))  
        abline(v = low_length, lty = 2)
    dev.off()

    # heatmap
    SID_High = rownames(tdata)[tdata$risk_group == "High risk"]
    SID_Low  = rownames(tdata)[tdata$risk_group == "Low risk"]
    pdata = t(tdata[c(SID_Low, SID_High), genes_inter])
    annotation_col = data.frame(Risk = rep_by_len(c("Low", "High"), list(SID_Low, SID_High)), row.names = colnames(pdata))
    annotation_colors = list(Risk = c("High" = km_cols[1], "Low" = km_cols[2]))

    p <- pheatmap(pdata, annotation_col = annotation_col, cluster_cols = FALSE, cluster_rows = FALSE,
                  scale = "column", annotation_colors = annotation_colors, show_colnames = FALSE, silent = TRUE)
    p <- ggplotify::as.ggplot(p) + theme(plot.margin = margin(l = 5, r = 10))
    save.graph(p, file = paste0(out_name, "_06_Gene_pheatmap"), 6, 3)
}

get_GEO_group <- function(obj, group_name, case_text, case_name, control_text, control_name){
    sample_anno = obj$sample_anno
    SID_case = sample_anno$geo_accession[grepl(case_text, sample_anno[, group_name])]
    SID_control = sample_anno$geo_accession[grepl(control_text, sample_anno[, group_name])]
    SID = c(SID_case, SID_control)
    group = rep_by_len(c(case_name, control_name), list(SID_case, SID_control))
    expM = obj$expM[, SID]
    sample_anno = sample_anno[SID, ]
    res = list2(expM, sample_anno, SID, SID_case, SID_control, group)
    return(res)
}

convert_pdf_to_png <- function(input, dpi = 300) {
    loadp(pdftools)
    if (dir.exists(input)) {
        pdf_files <- list.files(input, pattern = "\\.pdf$", recursive = TRUE, full.names = TRUE)
        for (pdf_file in pdf_files) {
            png_file <- sub("\\.pdf$", ".png", pdf_file)
            pdf_convert(pdf_file, format = "png", pages = NULL, dpi = dpi, filenames = png_file)
        }
    } else if (file.exists(input) && grepl("\\.pdf$", input)) {
        png_file <- sub("\\.pdf$", ".png", input)
        suppressMessages(pdf_convert(input, format = "png", pages = NULL, dpi = dpi, filenames = png_file))
    } else {
        stop("Input is not a valid PDF file or directory.")
    }
}

set_fanyi <- function(){
    loadp(fanyi)
    set_translate_option(appid = "20250118002256089", key = "II_MRM6yQMjMJPPHYX_Y", source = "baidu")
}

report_pathway <- function(pathways, out_name, top = 5, CN = TRUE, prepfix = ""){
    if(CN){
        loadp(fanyi)
        set_fanyi()
        text_trans = translate2(pathways[1:(ifelse(length(pathways) > top, top, length(pathways)))])
        text = prepfix
        text = paste0(text, "富集的top", top, "通路分别为：")
        text = paste0(text, paste(text_trans[1:(top-1)], "（", pathways[1:(top-1)], "）", collapse = "、"))
        text = paste0(text, "和", text_trans[top], "（", pathways[top], "）。")
    } else {
        text = prepfix
        text = paste0(text, "The top", top, " of significant enriched pathways include: ")
        text = paste0(text, paste(pathways[1:(top-1)], collapse = ", "))
        text = paste0(text, " and ", pathways[top], ".")
    }
    cat(text, file = paste0(out_name, "_top", top, "_text.txt"))
    invisible(text)
}

get_package_versions <- function(..., output_file = "package_versions.csv") {
  all_packages <- installed.packages()
  packages <- as.character(substitute(list(...)))[-1]
  
  if (length(packages) == 0) {
    package_table <- data.frame(
      Package = all_packages[, "Package"],
      Version = all_packages[, "Version"],
      stringsAsFactors = FALSE
    )
  } else {
    package_table <- data.frame(
      Package = packages,
      Version = sapply(packages, function(pkg) {
        if (pkg %in% all_packages[, "Package"]) {
          all_packages[pkg, "Version"]
        } else {
          NA
        }
      }),
      stringsAsFactors = FALSE
    )
  }
  rownames(package_table) = NULL

  if (length(packages) == 0) {
    write.csv(package_table, file = output_file, row.names = FALSE)
    message("Package versions saved to: ", output_file)
  } else {
      return(package_table)
  }  
}

numToStr <- function(num, len){
    paste0(paste(rep("0", len - nchar(num)), collapse = ""), num)
}

read.faster <- function (file, header = FALSE, sep = "\t", showProgress = FALSE, skip = 0, nrows = -1, ...){
    suppressPackageStartupMessages(library("data.table"))
    Tx = data.frame(fread(file = file, header = header, sep = sep, 
        showProgress = showProgress, skip = skip, nrows = nrows, ...))
    return(Tx)
}


dotplot3 <- function(enrich_tb, paths = NULL, top = 10, trans = FALSE, ...){
    if(is.null(paths)){
        pdata = enrich_tb
        pdata = pdata[1:ifelse(nrow(pdata) > top, top, nrow(pdata)), ] 
    } else {
        pdata = enrich_tb[match(paths, enrich_tb$Description), ]    
    }
    if(trans){
        loadp(extrafont)
        # font_import()
        loadfonts()
        pdata$Description = paste0(translate3(pdata$Description), "\n", pdata$Description)
    }
    invisible(dotplot2(pdata, PDF = TRUE, ...))
}

count_to_exp <- function(count, group){
    loadp(DESeq2)
    design = factor(group, levels = levels(factor(group)))
    dds = DESeqDataSetFromMatrix(count, DataFrame(design), design = ~ design)
    dds = DESeq(dds, fitType = "local")
    expM = data.frame(assay(vst(dds)), check.names = FALSE)
    return(expM)
}

tocap <- function(string) Hmisc::capitalize(tolower(string))

save.graph <- function(p, file, w = 4, h = 4, mar = NULL, PDF = TRUE, PNG = TRUE, open = FALSE, ...){
    if(!is.null(mar)){
        p <- p + theme(plot.margin = margin(mar[1], mar[2], mar[3], 
            mar[4]))
    }
    if(PDF){
        ggsave(paste0(file, ".pdf"), p, device = cairo_pdf, width = w, height = h, ...)
    }
    if(PNG){
        ggsave(paste0(file, ".png"), p, width = w, height = h)
    }
    if(open) shell.exec(paste0(file, ".pdf"))
}

# sc_mini_axis(p, obj, axis.len = 0.3)

sc_mini_axis <- function(p, 
                         obj, 
                         reduction = "umap", 
                         pos.just = 1,
                         label.just = 0.93,
                         axis.len = 0.2, 
                         axis.label.cex = 3, 
                         arrow.len = 0.3, 
                         arrow.size = 0.5){
    p <- p + theme_blank() + coord_fixed(ratio = 1)
    reduM = Embeddings(obj, reduction = reduction)
    pdata = data.frame(Dim_1 = reduM[, 1], Dim_2 = reduM[, 2])
    x.vector = pdata[, 1]
    y.vector = pdata[, 2]
    range_x = max(x.vector) - min(x.vector)
    range_y = max(y.vector) - min(y.vector)
    position.just = range_x * (pos.just - 1) * 0.1

    segTable = data.frame(
                      x = c(min(x.vector) * 1.1 * 1.0025 + position.just, min(x.vector) * 1.1 + position.just),
                      y = c(min(y.vector) * 1.1 / 1.0025 + position.just, min(y.vector) * 1.1 + position.just),
                      vx = c(range_x * axis.len, 0),
                      vy = c(0, range_x * axis.len))
    p <- p + geom_segment(data = segTable, 
                          mapping = aes(x = x, 
                                        y = y, 
                                        xend = x + vx, 
                                        yend = y + vy),
                          arrow = arrow(length = unit(arrow.len, "cm")), 
                          size = arrow.size)
    axis_title = ifelse(reduction == "umap", "UMAP", "tSNE")
    p <- p + annotate("text", 
                      x = min(x.vector) * label.just + position.just, 
                      y = min(y.vector) * 1.18 + position.just, 
                      label = paste0(axis_title, "_1"), 
                      size = axis.label.cex, hjust = 0.5, parse = TRUE)
    p <- p + annotate("text", 
                      x = min(x.vector) * 1.18 + position.just, 
                      y = min(y.vector) * label.just + position.just, 
                      label = paste0(axis_title, "_2"), 
                      size = axis.label.cex, hjust = 0.5, parse = TRUE, angle = 90)
    return(p)
}

DimPlot2 <- function(obj, 
                     group.by, 
                     legend.title = NULL,
                     title = "All", 
                     cols = col.cluster2.1,
                     reduction = "umap", 
                     cex = 20, 
                     size = 0.3, 
                     label = TRUE,
                     label.fill = "auto",
                     label.cex = 3, 
                     mini.axis = FALSE, 
                     legend = TRUE, 
                     pos.manual = NULL, ...){
    loadp(Seurat, dplyr, purrr, ggplot2)
    if(is.null(legend.title)) legend.title = group.by
    reduM = Embeddings(obj, reduction = reduction)
    pdata = data.frame(Dim_1 = reduM[, 1], Dim_2 = reduM[, 2], Cluster = obj@meta.data[, group.by])
    pdata %>% 
            group_by(Cluster) %>% 
            do(model = kmeans(.[c('Dim_1', 'Dim_2')], 1)) %>%
            ungroup() %>% group_by(Cluster) %>% 
            do(map_df(.$model, broom::tidy)) %>%
            ungroup() %>% select(Cluster, Dim_1, Dim_2) %>% data.frame() %>% 
            dplyr::rename(x.center = Dim_1, y.center = Dim_2, Cluster = Cluster) -> label.data
    if(!is.null(pos.manual)){
        for(celltype in names(pos.manual)){
            label.data[match(celltype, label.data$Cluster), "x.center"] = pos.manual[[celltype]][1]
            label.data[match(celltype, label.data$Cluster), "y.center"] = pos.manual[[celltype]][2]
        }
    }
    axis_title = ifelse(reduction == "umap", "UMAP", "tSNE")
    p <- ggplot(pdata) + geom_point(aes(Dim_1, Dim_2, color = Cluster), size = size, shape = 16, fill = "red")
    p <- p + scale_color_manual(values = cols, name = legend.title)
    p <- p + labs(x = paste0(axis_title, "_1"),
                  y = paste0(axis_title, "_2"),
                  title = paste0(title, " (n = ", format(ncol(obj), big.mark = ","), ")"))
    p <- p + theme(plot.title = element_text(size = cex / 20 * 18, face = "plain"))
    p <- p + setText(cex)

    if(label){
        if(label.fill == "auto"){
            p <- p + geom_label(data = label.data, aes(label = Cluster, x = x.center, y = y.center, fill = Cluster), size = label.cex)    
        } else if(label.fill == "none"){
            p <- p + geom_label(data = label.data, aes(label = Cluster, x = x.center, y = y.center), size = label.cex, fill = alpha(c("red"), 0))
        } else {
            p <- p + geom_label(data = label.data, aes(label = Cluster, x = x.center, y = y.center), size = label.cex, fill = label.fill)
        }
    }        

    p <- p + scale_fill_manual(values = col.cluster2.1, name = legend.title)
    p <- p + guides(fill = "none")
    p <- p + guides(color = guide_legend(override.aes = list(size = 5)))
    if(!legend) p <- p + theme(legend.position = "none") 
    if(mini.axis) p <- sc_mini_axis(p, obj, reduction = reduction, ...)
    return(p)
}

gg_style <- function(title.size   = 10,
                    title.face = "bold",
                    title.hjust = 0.5,
                    title.vjust = 0,
                    x.title.size = NA,
                    x.title.face = NA,
                    x.title.vjust = 0.5,
                    x.text.size  = NA,
                    x.text.angle = 0,
                    x.text.hjust = 0.5,
                    x.text.vjust = 1,
                    y.title.size = NA,
                    y.title.vjust = 2,
                    y.text.size  = NA,
                    face         = "plain",
                    y.title.face = NA,
                    y.text.angle = 0,
                    y.text.hjust = 1,
                    legend.position = "right",
                    plot.margin = c(5.5, 5.5, 5.5, 15),
                    color = NULL,
                    fill  = NULL,
                    theme = "classic"
){
    
    checkNA <- function(x, y, ratio = NA){
        
        if(is.na(x))
            if(is.na(ratio)){
                x = y
            } else{
                x = y * ratio
            }
        return(x)
    }
    
    x.title.size = checkNA(x.title.size, title.size, 0.8)
    y.title.size = checkNA(y.title.size, title.size, 0.8)
    x.text.size  = checkNA(x.text.size, title.size, 0.6)
    y.text.size  = checkNA(y.text.size, title.size, 0.6)
    x.title.face = checkNA(x.title.face, face)
    y.title.face = checkNA(y.title.face, face)
    if(x.text.angle){
        if(x.text.hjust != 45){
            x.text.hjust = 1
        }
    }
    if(y.text.angle)
        y.text.hjust = 1
    if(x.text.angle == 90){
        x.text.vjust = 0.5
    }

    mytheme = switch(theme,
            "bw" = theme_bw(),
            "classic"  = theme_classic(),
            "box" = theme_bw() + theme(panel.grid =element_blank())
            )
    p <- mytheme
    p <- p + theme(
               plot.title   = element_text(color = "black", size = title.size, face = title.face, hjust = title.hjust, vjust = title.vjust),
               axis.title.x = element_text(color = "black", size = x.title.size, face = x.title.face, vjust = x.title.vjust),
               axis.title.y = element_text(color = "black", size = y.title.size, face = y.title.face, vjust = y.title.vjust),
               axis.text.x  = element_text(color = "black", size = x.text.size, angle = x.text.angle, hjust = x.text.hjust, vjust = x.text.vjust),
               axis.text.y  = element_text(color = "black", size = y.text.size, angle = y.text.angle, hjust = y.text.hjust),
               plot.margin  = margin(t = 5.5, r = 5.5, b = 5.5, l = 15, unit = "pt"),
               legend.title    = element_text(size = unit(x.title.size, "pt")),
               legend.key.size = unit(x.text.size, "pt"),
               legend.text     = element_text(size = unit(x.text.size, "pt")),
               legend.position = legend.position)
    if(theme == "blank"){
        p <- p + theme(axis.line = element_blank(), 
                       axis.text.x = element_blank(),
                       axis.text.y = element_blank(), 
                       axis.ticks = element_blank(),
                       axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       panel.background = element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), plot.background = element_blank())
    }
    if(!is.null(color)) p <- p + scale_color_manual(values = color)
    if(!is.null(fill))  p <- p + scale_fill_manual(values = fill)
    return(p)
} 
setText <- gg_style

# prepare_gdc_download <- function
sort_by2 <- function(tb, by, group.by, ntop = NULL, decreasing = TRUE, select = NULL){
    tb_list = split(tb, f = tb[, group.by])
    order_list = lapply(tb_list, function(x) x[order(x[, by], decreasing = decreasing), ])
    if(!is.null(ntop)){
        order_list = lapply(order_list, function(x) x[1:(ifelse(ntop > nrow(x), nrow(x), ntop)), ])
    }
    res = Reduce(rbind, order_list)
    if(!is.null(select)){
        res = res[, c(group.by, by, select)]
    }
    return(res)
}

add_url <- function(string, type){
    url = switch(type,
        "pubchem" = "https://pubchem.ncbi.nlm.nih.gov/#query=",
        "geo_id"  = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",
        "pmid"    = "https://pubmed.ncbi.nlm.nih.gov/?term=",
        "genecard"    = "https://www.genecards.org/Search/Keyword?queryString=",
        "google"  = "https://www.google.com.hk/search?q="
        )
    res = paste0(url, string)
    return(res)
}


plot_enrich_circle_GO <- function(enrich_tb, gene_anno, 
    terms = 10, 
    n_gene = 20,
    space = 0.02, 
    gene.space = 0.2, 
    gene.size = 4,
    gene.col = c('firebrick3', 'white','steelblue')
    ){
    loadp(GOplot)
    if(class(terms) == "numeric"){
        tb = enrich_tb[1:terms, ]
    } else {
        tb = enrich_tb[match(terms, enrich_tb$Description), ]
    }

    ribbon.col = brewer.pal(nrow(tb), "Set3")
    enrich_data <- data.frame(
        Category = "ONTOLOGY", 
        ID = tb$ID,
        term = tb$Description, 
        Genes = gsub("/", ", ", tb$geneID), 
        adj_pval = tb$pvalue)
    colnames(gene_anno)[1] = "ID"
    circ = circle_dat(enrich_data, gene_anno)
    
    # check gene
    genes_select = names(sort(table(circ$genes), decreasing = TRUE))[1:n_gene]
    if(sum(unlist(strsplit(gene_anno$ID[1:10], "")) %in% letters) > 1){
      genes_select = tocap(genes_select)
      circ$genes = tocap(circ$genes)
    }

    chord = chord_dat(circ, gene_anno[match(genes_select, gene_anno$ID), ], tb$Description)
    chord = chord[, apply(chord, 2, function(y) !all(y == 0))]
    chord = na.omit(chord)
    ribbon.col = ribbon.col[1:(ncol(chord)-1)]
    p <- GOChord(chord, space = space, gene.order = colnames(gene_anno)[2], 
        gene.space = gene.space, 
        gene.size = gene.size,
        lfc.col = gene.col,
        ribbon.col = ribbon.col)
    p <- p + guides(size = guide_legend(title = "", order = 2, ncol = 1, byrow = F, override.aes = list(shape = 22, fill = ribbon.col, size = 8)),
                    fill = guide_colorbar(title.position = "top", title.hjust = 0.5, order = 1),
                    shape = "none")
    p <- p + theme(legend.position = "right", 
        legend.background = element_rect(fill = "transparent"), 
        legend.box = "vertical", legend.direction = "horizontal")
    p <- p + coord_fixed(ratio = 1)
    N = nrow(p$layers[[4]]$data)
    angle = (seq(180 / (N * 2), 180, 180 / N))
    p$layers[[4]]$data$angle = angle
    return(p)
}

plot_enrich_circle_KEGG <- function(enrich_tb, gene_anno, 
    terms = 10, 
    n_gene = 20,
    space = 0.02, 
    gene.space = 0.2, 
    gene.size = 4,
    gene.col = c('firebrick3', 'white','steelblue')
    ){
    loadp(GOplot)
    if(class(terms) == "numeric"){
        tb = enrich_tb[1:terms, ]
    } else {
        tb = enrich_tb[match(terms, enrich_tb$Description), ]
    }
    ribbon.col = brewer.pal(nrow(tb), "Set3")

    enrich_data <- data.frame(
        Category = "KEGG", 
        ID = tb$ID,
        term = tb$Description, 
        Genes = gsub("/", ", ", tb$Symbol), 
        adj_pval = tb$pvalue)
    colnames(gene_anno)[1] = "ID"
    circ = circle_dat(enrich_data, gene_anno)
    genes_select = names(sort(table(circ$genes), decreasing = TRUE))[1:n_gene]

    # check gene
    genes_select = names(sort(table(circ$genes), decreasing = TRUE))[1:n_gene]
    if(sum(unlist(strsplit(gene_anno$ID[1:10], "")) %in% letters) > 1){
      genes_select = tocap(genes_select)
      circ$genes = tocap(circ$genes)
    }

    chord = chord_dat(circ, gene_anno[match(genes_select, gene_anno$ID), ], tb$Description)
    chord = chord[, apply(chord, 2, function(y) !all(y == 0))]
    chord = na.omit(chord)
    ribbon.col = ribbon.col[1:(ncol(chord)-1)]
    p <- GOChord(chord, space = space, gene.order = colnames(gene_anno)[2], 
        gene.space = 0.25, 
        gene.size = 3,
        lfc.col = gene.col,
        ribbon.col = ribbon.col)
    p <- p + guides(size = guide_legend(title = "", order = 2, ncol = 1, byrow = F, override.aes = list(shape = 22, fill = ribbon.col, size = 8)),
                    fill = guide_colorbar(title.position = "top", title.hjust = 0.5, order = 1),
                    shape = "none")
    p <- p + theme(legend.position = "right", 
        legend.background = element_rect(fill = "transparent"), 
        legend.box = "vertical", legend.direction = "horizontal")
    p <- p + coord_fixed(ratio = 1)
    N = nrow(p$layers[[4]]$data)
    angle = (seq(180 / (N * 2), 180, 180 / N))
    p$layers[[4]]$data$angle = angle
    return(p)
}

runEGO <- function (sigG, pvalueCutoff = 0.05, qvalueCutoff = 0.05, universe = NULL, 
    out = NULL, keyType = "ENSEMBL", ont = "ALL", sig = TRUE, species = "human",
    raw = FALSE) 
{
    if (sum(c("clusterProfiler", "org.Hs.eg.db") %in% 
        .packages()) != 3) {
        message("Loading packages ...")
        loadp(clusterProfiler, org.Hs.eg.db)
    }

    if(species == "human"){
        OrgDb = "org.Hs.eg.db"
    } else if(species == "mouse"){
        OrgDb = "org.Mm.eg.db"
    }
    TYPE = ifelse(subString(sigG[1], 1:4) == "ENSG", 
                  "ENSEMBL", 
                  ifelse(subString(sigG[1], 1) %in% LETTERS, "SYMBOL", "ENTREZID"))
    if (TYPE != "ENSEMBL") 
        sigG = convertId(sigG, TYPE, "ENSEMBL", OrgDb = OrgDb)
    message("Running ONTOLOGY: ", ont, " ...")
    suppressMessages(ego <- enrichGO(gene = sigG, keyType = keyType, 
        OrgDb = OrgDb, ont = ont, pAdjustMethod = "fdr", 
        pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, 
        readable = TRUE))
    if (sig) {
        ego = ego[ego$pvalue < 0.05, ]
    }
    if (raw) 
        return(ego)
    egoT = convertEnrichT("go", data.frame(ego))
    egoT$Count = sapply(strsplit(egoT$geneID, "/"), length)
    if (!is.null(out)) 
        write.csv(egoT, file = out)
    if (is.null(nrow(egoT))) {
        message("No terms enriched")
        egoT[1, ] = rep(NA, ncol(egoT))
    }
    else {
        message(paste(nrow(egoT), " terms enriched"))
    }
    return(egoT)
}

runDO <- function (sigG, pvalueCutoff = 0.05, qvalueCutoff = 0.05, universe = NULL, 
    out = NULL, ont = "DO", sig = TRUE, species = "human",
    raw = FALSE) 
{
    suppressPackageStartupMessages(library("clusterProfiler"))
    suppressPackageStartupMessages(library("org.Hs.eg.db"))
    if(species == "human"){
        OrgDb = "org.Hs.eg.db"
    } else if(species == "mouse"){
        OrgDb = "org.Mm.eg.db"
    }
    TYPE = ifelse(subString(sigG[1], 1:4) == "ENSG", "ENSEMBL", 
        ifelse(subString(sigG[1], 1) %in% LETTERS, "SYMBOL", 
            "ENTREZID"))
    if (TYPE != "ENSEMBL") 
        sigG = convertId(sigG, TYPE, "ENTREZID", OrgDb = OrgDb)                    
    organism = ifelse(species == "human", "hsa", "mmu")

    suppressMessages(edo <- enrichDO(gene = sigG, 
        ont = ont, pAdjustMethod = "fdr", 
        pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, 
        readable = TRUE))
    if (sig) {
        edo = edo[edo$pvalue < 0.05, ]
    }
    if (raw) 
        return(edo)
    # egoT = convertEnrichT("go", data.frame(edo))
    edoT = edo
    edoT$Count = sapply(strsplit(edoT$geneID, "/"), length)
    if (!is.null(out)) 
        write.csv(edoT, file = out)
    if (is.null(nrow(edoT))) {
        message("No terms enriched")
        edoT[1, ] = rep(NA, ncol(edoT))
    }
    else {
        message(paste(nrow(edoT), " terms enriched"))
    }
    return(edoT)
}

runKEGG <- function(sigG, pvalueCutoff = 0.05, qvalueCutoff = 0.05, species = "human", out = NULL){
    suppressPackageStartupMessages(library("clusterProfiler"))
    suppressPackageStartupMessages(library("org.Hs.eg.db"))
    if(species == "human"){
        OrgDb = "org.Hs.eg.db"
    } else if(species == "mouse"){
        OrgDb = "org.Mm.eg.db"
    }
    TYPE = ifelse(subString(sigG[1], 1:4) == "ENSG", "ENSEMBL", 
        ifelse(subString(sigG[1], 1) %in% LETTERS, "SYMBOL", 
            "ENTREZID"))
    if (TYPE != "ENTREZID") 
        sigG = convertId(sigG, TYPE, "ENTREZID", OrgDb = OrgDb)                    
    organism = ifelse(species == "human", "hsa", "mmu")
    suppressMessages(ekegg <- enrichKEGG(gene = sigG, organism = organism, 
        keyType = "kegg", pAdjustMethod = "fdr", pvalueCutoff = pvalueCutoff, 
        qvalueCutoff = qvalueCutoff, use_internal_data = FALSE))
    ekeggT = convertEnrichT("kegg", ekegg, species = species)
    ekeggT = ekeggT[ekeggT$pvalue < 0.05, ]
    if (!is.null(out)) 
        write.csv(ekeggT, file = out)
    if (is.null(nrow(ekeggT))) {
        cat("No terms enriched\n")
        ekeggT[1, ] = rep(NA, ncol(ekeggT))
    }
    else {
        cat(paste(nrow(ekeggT), "terms enriched\n"))
    }
    return(ekeggT)
}

convertEnrichT <- function(type, enrichT, species = "human"){
    if(species == "human"){
        OrgDb = "org.Hs.eg.db"
    } else if(species == "mouse"){
        OrgDb = "org.Mm.eg.db"
    }
    Tx = data.frame(enrichT)
    if (type == "go") {
        Tx$geneID = as.character(unlist(sapply(Tx$geneID, function(x) paste0(unique(strsplit(x, 
            "/")[[1]]), collapse = "/"))))
    }
    else {
        Tx$Symbol = as.character(unlist(sapply(Tx$geneID, function(x) paste0(unique(sort(convertId(strsplit(x, 
            "/")[[1]], "ENTREZID", "SYMBOL", OrgDb = OrgDb))), collapse = "/"))))
    }
    return(Tx)
}



FindAllMarkers_genesort <- function(obj, ntop = 100, group.by = NULL){
    check_pkgs(genesorteR)
    if(!is.null(group.by)) Idents(obj) = group.by
    gs = suppressWarnings(sortGenes(obj@assays$RNA@data, Idents(obj)))
    gs_marker = getMarkers(gs, quant = 0.99)
    marker_tb = apply(spec_scores, 2, function(x) rownames(spec_scores)[order(spec_scores[, 1], decreasing = TRUE)[1:ntop]])
    marker_tb = reshape2::melt(a)[, c(2, 3)]
    colnames(marker_tb) = c("cluster", "gene")
    return(marker_tb)
}

read_10x_dir <- function(path){
    counts = Read10X(path)
    obj = CreateSeuratObject(counts = counts)
    return(obj)
}
    
read_10x_dir_raw <- function(path){
    # rename
    files = lf("rf", path)
    new_names = gsub("genes.tsv.gz", "features.tsv.gz", files)
    file.rename(files, new_names)

    files = lf("rf", path, "barcodes")
    SIDs = subString(basename(files), 1, "_barcodes")
    single_dir = paste0(path, "_single")
    for(SID in SIDs){
        xsingle_dir = paste0(single_dir, "/", SID)
        mkdir(xsingle_dir)
        from = paste0(path, "/", SID, "_", c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"))
        to = paste0(xsingle_dir, "/",  c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"))
        file.copy(from, to)
    }
    
    loadp(Seurat)
    input_dir = paste0(single_dir, "/", SIDs)
    obj_list = sapply(input_dir, read_10x_dir)
    obj_all <- merge(obj_list[[1]], y = obj_list[-1], add.cell.ids =  SIDs) 
    return(obj_all)
}


merge_strings <- function(strings, max.len = 5000, collapse = "|"){
    
    combined_str <- paste(strings, collapse = collapse)
    pos <- c(0)
    current_len <- 0
    for (i in 1:length(strings)) {
        current_len <- current_len + nchar(strings[i]) + 1  # 包括一个"|"的长度
        if (current_len > max.len) {
            pos <- c(pos, i - 1)
            current_len <- nchar(strings[i]) + 1
        }
    }
    result <- c()
    for(i in 1:length(pos)){
        start = pos[i] + 1
        end = pos[i+1]
        if(i == length(pos))  end = length(strings)
        part = paste(strings[start:end], collapse = collapse)
        result = c(result, part)
    }
    return(result)
}

translate2 <- function(strings, max.len = 5000){
    loadp(fanyi)
    strings_merge = merge_strings(strings, max.len = max.len, collapse = "=")
    trans_merge = translate(strings_merge)
    results = unlist(strsplit(trans_merge, split = "\\="))
    if(is.null(results)) results = translate2(strings, max.len)
    return(results)
}

translate3 <- function(strings, sleep = 1, ...){
    Sys.sleep(sleep)
    loadp(fanyi)
    translate(strings, ...)
}


gene_summary2 <- function(genes){
    loadp(fanyi)
    id = convertId(genes, "SYMBOL", "ENTREZID")
    tb = gene_summary(id)
    result = data.frame(
        gene_name = tb$name, 
        gene_name_cn = paste0(tb$name, "（", translate2(tb$description), "）"),
        description = tb$description,
        description_cn = translate2(tb$description),
        summary = tb$summary,
        summary_cn = translate2(tb$summary)
        )
    return(result)
}

rep_by_len <- function(x, y){
    if(length(x) == 1){
        rep(x, length(y))
    } else {
        Reduce(c, sapply(1:length(x), function(i) rep(x[i], length(y[[i]]))))
    }
}

adjust_range <- function(numbers, min_num, max_num) {
    min_input <- min(numbers)
    max_input <- max(numbers)    
    scaled_numbers <- ((numbers - min_input) / (max_input - min_input)) * (max_num - min_num) + min_num    
    return(scaled_numbers)
}

check_pkgs <- function(..., load.pkg = TRUE){
    pkgs = as.character(substitute(list(...)))[-1]
    not_installed = setdiff(pkgs, rownames(installed.packages()))
    if(length(not_installed) > 0){
        cat(paste0("Please install package(s) '", paste(not_installed, collapse = "' '"), "' first!\n"))
        stop()
    }
    if(load.pkg)
        suppressMessages(for (pkg in pkgs) loadp(pkg, character.only = TRUE))
}

check_colnames <- function(df, Names){
    not_exists = setdiff(Names, colnames(df))
    if(length(not_exists) > 0){
        cat(paste0("Colnames '", paste(not_exists, collapse = "' '"), "' not exists in data!\n"))
        stop()
    }
}

seurat2cds <- function(obj, group.by, out_dir, ordering_genes = NULL, method = "marker", add.genes = NULL){
    choices = c("marker", "HDG", "HVG")
    match.arg(method, choices)

    loadp(Seurat, monocle)
    cds = as.CellDataSet(obj)
    cds = estimateSizeFactors(cds)
    cds = estimateDispersions(cds)
    cds = detectGenes(cds, min_expr = 1)

    # 三种方法: https://lishensuo.github.io/posts/bioinfo/008%E5%8D%95%E7%BB%86%E8%83%9E%E5%88%86%E6%9E%90%E5%B7%A5%E5%85%B7--monocle%E8%BD%A8%E8%BF%B9%E5%88%86%E6%9E%90/
    if(is.null(ordering_genes)){
        if(method == "marker"){
            # marker gene by Seurat, 首选, 其他的都太慢
            Idents(obj) = group.by
            marker_all = FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25)
            ordering_genes = rownames(marker_all)
        } else if(method == "HDG"){
            gene_Disp = dispersionTable(cds)
            ordering_genes = gene_Disp %>% 
              dplyr::filter(mean_expression >= 0.1,
                            dispersion_empirical >= dispersion_fit) %>% 
              pull(gene_id) %>% unique()
        } else if(method == "HVG"){
            ordering_genes = VariableFeatures(obj)
        }
    }

    if(!is.null(add.genes)) ordering_genes = unique(c(ordering_genes, add.genes))

    gbm_cds = setOrderingFilter(cds, ordering_genes = ordering_genes)
    gbm_cds = reduceDimension(gbm_cds, max_components = 2, verbose = T,check_duplicates = F,num_dim = 10)
    # save(gbm_cds, file = paste0(out_dir, "/01_before_orderCells.RData"))

    gbm_cds = orderCells(gbm_cds, reverse = F)
    # save(gbm_cds, file = paste0(out_dir, "/02_orderCells.RData"))
    return(gbm_cds)
}


seurat2cellchat <- function(obj, group.by){
    meta <- obj@meta.data
    cellchat <- createCellChat(object = obj, meta = meta, group.by = group.by)
    cellchat <- addMeta(cellchat, meta = meta)
    cellchat <- setIdent(cellchat, ident.use = group.by)
    groupSize <- as.numeric(table(cellchat@idents))

    CellChatDB <- CellChatDB.human
    CellChatDB.use <- CellChatDB
    cellchat@DB <- CellChatDB.use
    cellchat <- subsetData(cellchat)

    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    # cellchat <- projectData(cellchat, PPI.human)
    cellchat <- computeCommunProb(cellchat)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    return(cellchat)
}

sc_pipeline <- function(obj_raw, npc = 30, ncluster = c(14, 16), resolution = 0.2, 
                              sample_name = "SID", out_name = "sc_pipeline",
                              cell.min = 3, mt.pct = 10, hb.pct = 3,
                              feature.min = 200, feature.max = 8000, tSNE = TRUE,
                              count.min = 0, count.max = Inf, species = "human",
                              verbose = FALSE
                             ){
    # # test
    # obj_raw = obj
    # npc = 30
    # ncluster = c(14, 16)
    # resolution = 0.2
    # sample_name = "SID"
    # out_dir = "."
    # cell.min = 3
    # mt.pct = 10
    # hb.pct = 3
    # feature.min = 200
    # feature.max = 6000
    # count.min = 0
    # count.max = Inf
    # species = "human"

    loadp(ggplot2, Seurat)

    mkdir(dirname(out_name))

    # 确保Idents唯一
    if(length(unique(Idents(obj_raw))) != 1){
        Idents(obj_raw) = "All cells"
    }

    report("Checking sample name ...")
    if(! sample_name %in% colnames(obj_raw@meta.data)){
        stop(paste0("Sample name '", sample_name, "' not in meta data"))
    }

    report("Plotting MT and HB percent ...")
    obj_raw[["percent.MT"]] = PercentageFeatureSet(obj_raw, pattern = ifelse(species == "human", "^MT-", "^mt-"))
    HB.genes = c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
    if(species == "mouse"){
        HB.genes = homologene::homologene(HB.genes, inTax = 9606, outTax = 10090)[, 2]
    }

    HB.genes = na.omit(rownames(obj_raw@assays$RNA)[match(HB.genes, rownames(obj_raw@assays$RNA))])
    obj_raw[["percent.HB"]] = PercentageFeatureSet(obj_raw, features = HB.genes) 

    features = c("nFeature_RNA", "nCount_RNA", "percent.MT", "percent.HB")
    col.num = length(levels(obj_raw@active.ident))

    plotL = list()
    for(feature in features){
        p <- suppressWarnings(VlnPlot(obj_raw, features = feature, cols = col.cluster2.1, pt.size = 0.1, ncol = 1))
        p <- p + setText(25, x.text.angle = 90)
        p <- p + theme(legend.position = "none") + labs(x = "")
        plotL[[feature]] <- p
    }
    set_image(12, 7)
    p <- gridExtra::grid.arrange(grobs = plotL, nrow = 1) 
    save.graph(p, file = paste0(out_name, "_01_QC_before"), w = 15, h = 7)
    flush.console()

    # qc
    obj_qc = obj_raw
    obj_qc = obj_qc[, obj_qc$nFeature_RNA > feature.min & obj_qc$nFeature_RNA < feature.max]
    obj_qc = obj_qc[, obj_qc$nCount_RNA > count.min & obj_qc$nCount_RNA < count.max]
    obj_qc = obj_qc[, obj_qc$percent.MT < mt.pct & obj_qc$percent.HB < hb.pct]

    plotL = list()
    for(feature in features){
        p <- suppressWarnings(VlnPlot(obj_qc, features = feature, cols = col.cluster2.1, pt.size = 0.1, ncol = 1))
        p <- p + setText(25, x.text.angle = 90)
        p <- p + theme(legend.position = "none") + labs(x = "")
        plotL[[feature]] <- p
    }
    p <- gridExtra::grid.arrange(grobs = plotL, nrow = 1) 
    save.graph(p, file = paste0(out_name, "_02_QC_after"), w = 15, h = 7)
    flush.console()

    report("Plotting high variable genes")
    obj = obj_qc
    obj = NormalizeData(obj, verbose = verbose)
    obj = FindVariableFeatures(obj, selection.method = "vst", verbose = verbose)
    obj = ScaleData(obj, features = VariableFeatures(obj), verbose = verbose)
    top10 = head(VariableFeatures(obj), 10)
    plot1 <- VariableFeaturePlot(obj) + theme(legend.position = "bottom")
    plot2 <- suppressMessages(LabelPoints(plot = plot1, points = top10, repel = TRUE))
    p <- plot1 + plot2
    suppressWarnings(save.graph(p, file = paste0(out_name, "_03_High_variable_genes"), w = 10, h = 5))

    gc2()
    report("Running PCA")
    obj <- RunPCA(obj, features = VariableFeatures(obj), verbose = verbose)
    p <- DimPlot(obj, reduction = "pca", group.by = sample_name, cols = rep(col.cluster2.1, 3))
    save.graph(p, paste0(out_name, "_04_PCA_by_samples"), 7, 6)

    pdf(paste0(out_name, "_05_PC_Heatmap.pdf"), w = 10, h = 15)
        DimHeatmap(obj, dims = 1:15, cells = 100, balanced = TRUE)
    dev.off()

    p <- ElbowPlot(obj, ndims = 50, reduction = "pca")
    save.graph(p, paste0(out_name, "_06_PCA_importance"), 7, 6)

    report("Clustering ...")
    loadp(Seurat, harmony)
    pc.num = 1:npc
    obj = RunHarmony(obj, group.by.vars = sample_name, verbose = FALSE, project.dim = FALSE)
    obj = suppressWarnings(RunUMAP(obj, reduction = "harmony", dims = pc.num, verbose = FALSE))
    obj = FindNeighbors(obj, reduction = "harmony", dims = pc.num, verbose = FALSE)

    gc2()
    if(tSNE) obj = RunTSNE(obj, reduction = "harmony", dims = pc.num, check_duplicates = FALSE)
    obj = FindClusters2(obj, cluster.range = ncluster, res = resolution, verbose = TRUE)

    p <- DimPlot(obj, reduction = "umap", label = FALSE, group.by = sample_name) 
    p <- p + labs(title = "UMAP by Sample")
    p <- p + setText(20, y.title.vjust = 3)
    p <- p + scale_color_manual(values = col.cluster2.1, name = "Sample")
    save.graph(p, file = paste0(out_name, "_07_UMAP_by_Sample"), 7, 6)

    if(tSNE){
        p <- DimPlot(obj, reduction = "tsne", label = FALSE, group.by = sample_name) 
        p <- p + labs(title = "tSNE by Sample")
        p <- p + setText(20, y.title.vjust = 3)
        p <- p + scale_color_manual(values = col.cluster2.1, name = "Sample")
        save.graph(p, file = paste0(out_name, "_07_tSNE_by_Sample"), 7, 6)
    }

    obj$Cluster = obj$seurat_clusters
    p <- DimPlot(obj, reduction = "umap", label = FALSE, group.by = "Cluster") 
    p <- p + labs(title = "UMAP by Cluster")
    p <- p + setText(20, y.title.vjust = 3)
    p <- p + scale_color_manual(values = col.cluster2.1, name = "Cluster")
    save.graph(p, file = paste0(out_name, "_08_UMAP_by_Cluster"), 7, 6)

    if(tSNE){
        obj$Cluster = obj$seurat_clusters
        p <- DimPlot(obj, reduction = "tsne", label = FALSE, group.by = "Cluster") 
        p <- p + labs(title = "tSNE by Cluster")
        p <- p + setText(20, y.title.vjust = 3)
        p <- p + scale_color_manual(values = col.cluster2.1, name = "Cluster")
        save.graph(p, file = paste0(out_name, "_08_tSNE_by_Cluster"), 7, 6)
    }
    sc_dotplot_general(obj, group = "Cluster", file = paste0(out_name, "_09_DotPlot_General"), species = species)
    return(obj)
}

add_mini_axis <- function(p, obj, width = 6){

    umapM = Embeddings(obj, reduction = "umap")
    pdata = data.frame(UMAP_1 = umapM[, 1], UMAP_2 = umapM[, 2])

    x.vector = pdata[, 1]
    y.vector = pdata[, 2]
    segTable = data.frame(x = c(min(x.vector) - 0.8 - width / 120, min(x.vector) - 0.8),
                          y = c(min(y.vector) - 0.8, min(y.vector) - 0.8),
                          vx = c(width - 0.5, 0),
                          vy = c(0, width))
    p <- p + theme_blank() + coord_fixed()
    p <- p + geom_segment(data = segTable, mapping = aes(x = x, y = y, xend = x + vx, yend = y + vy),
                                                         arrow = arrow(length=unit(0.25, "cm")), 
                                                         size = 0.5)
    p <- p + annotate("text", x = segTable[1, "x"] + width / 8, y = segTable[1, "y"] - width / 6, label = "UMAP_1", size = 3, hjust = 0, parse = TRUE)
    p <- p + annotate("text", x = segTable[2, "x"] - width / 6, y = segTable[2, "y"] + width / 6, label = "UMAP_2", size = 3, hjust = 0, parse = TRUE, angle = 90)
    return(p)
}


set_image <- function(width = 8, height = 8){
    options(repr.plot.width = width, repr.plot.height = height)
}

options(Seurat.object.assay.version = "v3")

add_clonotype <- function(tcr_prefix, seurat_obj, type="t"){
    tcr <- read.csv(paste0(tcr_prefix, "/filtered_contig_annotations.csv"))
    tcr <- tcr[!duplicated(tcr$barcode), ]
    tcr <- tcr[,c("barcode", "raw_clonotype_id")]
    names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"

    clono <- read.csv(paste0(tcr_prefix, "/clonotypes.csv"))
    tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")])
    names(tcr)[names(tcr) == "cdr3s_aa"] <- "cdr3s_aa"

    tcr <- tcr[, c(2,1,3)]
    rownames(tcr) <- tcr[,1]
    tcr[,1] <- NULL
    colnames(tcr) <- paste(type, colnames(tcr), sep="_")
    clono_seurat <- AddMetaData(object=seurat_obj, metadata=tcr)
    return(clono_seurat)
}

extractVDJ <- function(vdjOutDir, type = NULL){

    if(is.null(type)){
        stop("Please specific 'type', t/b.")
    }
    contigAnno = read.csv(paste0(vdjOutDir, "/filtered_contig_annotations.csv"))
    contigAnno = contigAnno[!duplicated(contigAnno$barcode), ]
    contigAnno = contigAnno[,c("barcode", "raw_clonotype_id")]
    names(contigAnno)[names(contigAnno) == "raw_clonotype_id"] <- "clonotype_id"

    clonoTypes = read.csv(paste0(vdjOutDir, "/clonotypes.csv"))
    contigAnno = merge(contigAnno, clonoTypes[, c("clonotype_id", "cdr3s_aa")])

    contigAnno = contigAnno[, c(2, 1, 3)]
    rownames(contigAnno) = contigAnno[, 1]
    contigAnno[, 1] = NULL
    colnames(contigAnno) <- paste0(type, "_", colnames(contigAnno))
    return(contigAnno)
}

readr <- function(file_path){
    res = get(load2(file_path))
    return(res)
}

readRData <- function(file_path, idx = NULL, rm = FALSE, verbose = FALSE){
    obj_name = load(file_path)

    if(is.null(idx) & length(obj_name) > 1){
            cat2("This RData contains multiple objects: ", paste0(obj_name, collapse = ", "), "")
            cat2("Default return the first object, use 'idx' to select manully.")           
    }

    if(is.null(idx)){
        idx = 1
    }

    if(length(idx) > 1){
        if(verbose){
            cat2("Loadding object: ", paste0(obj_name[idx], collapse = " "))    
        }
        res = sapply(idx, function(x) get(obj_name[x]), simplify = FALSE)
        names(res) = obj_name[idx]
    } else if(idx == "all"){
        idx = 1:length(obj_name)
        res = sapply(idx, function(x) get(obj_name[x]), simplify = FALSE)
        names(res) = obj_name[idx]
    } else {        
        if(verbose){
            cat2("Loadding object: ", paste0(obj_name[idx], collapse = " "))    
        }
        res = get(obj_name[idx])
    }
    return(res)
}

# to do: cluster可指定多个
sc_subset_obj <- function(obj, cluster, group = "seurat_clusters", verbose = TRUE){
    if(!group %in% colnames(obj@meta.data)){
        info = paste0("Given group '", group, "' can not be found in meta data, please check it!")
        report(info, "E")
        stop()
    }

    if(!cluster %in% obj@meta.data[, group]){
        info = paste0("Given cluster '", cluster, "' can not be found in specific group '", group, "', please check it!")
        report(info, "E")
        stop()       
    }
    xobj = obj[, obj@meta.data[, group] == cluster]
    if(verbose){
        report(paste0("Subset ", ncol(xobj), "/", ncol(obj), " cells from cluster '", cluster, "' in group '", group, "'"))
    }    

    return(xobj)
}

sc_get_count <- function(obj, cluster = NULL, group = "seurat_clusters", verbose = TRUE){
    if(!is.null(cluster)){
        cluster %in% obj@meta.data[, group]
        obj = sc_subset_obj(obj, cluster, group, verbose = verbose)
    }
    count_tb = obj@assays$RNA@counts
    return(count_tb)
}

sc_VolcanoPlot <- function(diff, 
                           log2FC = log2(1.5), 
                           padj = 0.05, 
                           label.symbols = NULL, 
                           label.max = 30,
                           cols = c("#F93463", "#27A2EF"), 
                           title = ""){
    if( all( !c("log2FoldChange", "padj", "symbol") %in% colnames(diff) )){
    stop("Colnames must include: log2FoldChange, padj, symbol")
    }
    rownames(diff) = diff$symbol

    # (1) define up and down
    diff$threshold="ns";
    if(length(which(diff$log2FoldChange > log2FC & diff$padj < padj)) != 0){
        diff[which(diff$log2FoldChange > log2FC & diff$padj <padj),]$threshold="up"
    }

    if(length(which(diff$log2FoldChange < (-log2FC) & diff$padj < padj)) != 0){
        diff[which(diff$log2FoldChange < (-log2FC) & diff$padj < padj),]$threshold="down";
    }

    diff$threshold=factor(diff$threshold, levels=c('down','ns','up'))
    #head(diff)
    #
    tb2=table(diff$threshold)# ; print(tb2)
    library(ggplot2)
    # (2) plot
    g1 = ggplot(data=diff, aes(x=log2FoldChange, y=-log10(padj), color=threshold)) +
    geom_point(alpha=0.8, size=1.5) +
    geom_vline(xintercept = c(-log2FC, log2FC), linetype=2, color="grey")+
    geom_hline(yintercept = -log10(padj), linetype=2, color="grey")+
    labs(title= ifelse(""==title, "", paste(title)))+
    xlab(bquote(Log[2]*FoldChange))+
    ylab(bquote(-Log[10]*italic(P.adj)) )+
    theme_classic(base_size = 14) +
    theme(legend.box = "horizontal",
          legend.position="top",
          legend.spacing.x = unit(0, 'pt'),
          legend.text = element_text( margin = margin(r = 20) ),
          legend.margin=margin(b= -10, unit = "pt"),
          plot.title = element_text(hjust = 0.5, size=10)
          ) +
    coord_flip() + setText(20)

    no_up = length(which(diff$log2FoldChange > log2FC & diff$padj < padj)) == 0
    no_dn = length(which(diff$log2FoldChange < (-log2FC) & diff$padj < padj)) == 0

    if(no_up){
        g1 <- g1 + scale_color_manual('',labels=c(paste0("down(",tb2[[1]],')'),'ns'),
                           values=c(cols[1], "grey") )
    }

    if(no_dn){
        g1 <- g1 + scale_color_manual(labels=c('ns',
                                           paste0("up(",tb2[[3]],')' )),
                               values=c("grey", cols[2]) )
    }

    if(!no_up & !no_dn){
        g1 <- g1 + scale_color_manual('',labels=c(paste0("down(",tb2[[1]],')'),'ns',
                                       paste0("up(",tb2[[3]],')' )),
                           values=c(cols[1], "grey", cols[2]) )
    }        

    g1 <- g1 + guides(color=guide_legend(override.aes = list(size=3, alpha=1)))

    # (3)label genes
    if(is.null(label.symbols)){
    diff.sig=diff[which(diff$threshold != "ns" ), ]
    len=nrow(diff.sig)
    if(len<label.max){
      label.symbols=rownames(diff.sig)
    }else{
      diff.sig=diff.sig[order(diff.sig$log2FoldChange), ]
      diff.sig= rbind(diff.sig[1:(label.max/2),], diff.sig[(len-label.max/2):len,])
      label.symbols=rownames(diff.sig)
    }
    }
    dd_text = diff[label.symbols, ]
    # print((dd_text))
    # add text
    library(ggrepel)
    p <- g1 + geom_text_repel(data=dd_text,
                       aes(x=log2FoldChange, y=-log10(padj), label=row.names(dd_text)),
                         #size=2.5, 
                         colour="black",alpha=1, max.overlaps = 50)
    return(p)
}


sc_merge_qc <- function(input_dir, type = "count"){
    if(type %!in% c("count", "multi_5_Only", "multi_5_PE")){
        cat("Only support type: count|multi_5_Only|multi_5_PE!\n")
        stop()
    }

    key_list <- list(
        count = c(
            "Estimated Number of Cells",
            "Median Genes per Cell",
            "Number of Reads",
            "Valid Barcodes",
            "Sequencing Saturation",
            "Reads Mapped to Genome",
            "Reads Mapped Confidently to Genome",
            "Fraction Reads in Cells",
            "Total Genes Detected",
            "Median UMI Counts per Cell"),
        multi_5_Only = c(40, 28, 42, 6, 4, 41, 36, 2, 45, 23, 26, 64, 72, 13, 17, 48, 56),
        multi_5_PE   = c(41, 44, 43, 6, 4, 42, 37, 2, 46, 23, 26, 65, 73, 13, 17, 49, 57)
    )

    tags = lf("x", input_dir)
    qc_table = data.frame()
    if(type == "count"){
        for(tag in tags){
            xqc_table = read.csv(paste0(input_dir, "/", tag, "/outs/metrics_summary.csv"), check.names = FALSE, header = TRUE)
            xqc_table = cbind(Sample = tag, xqc_table[, key_list[[type]]])
            qc_table = rbind(qc_table, xqc_table)
        }
    } else {
        for(tag in tags){
            xqc_table = read.csv(paste0(input_dir, "/", tag, "/outs/per_sample_outs/", tag, "/metrics_summary.csv"), check.names = FALSE)
            xqc_table = xqc_table[key_list[[type]], c(2, 5:6)]
            xqc_table = cbind(Sample = tag, xqc_table)
            qc_table = rbind(qc_table, xqc_table)
        }
    }
    return(qc_table)
}

sc_dotplot_general <- function(obj, group, file, PDF = TRUE, PNG = TRUE, species = "human", w = 10, h = 9){
    celltypes = c("T", "NK", "B", "Macro", "Mono", "DC", "Neutrophils", "Endo", "Epis", "Fibro")
    markers = unique(unlist(sapply(celltypes, selectMarkers)))
    names(markers) = NULL
    
    if(species == "mouse"){
        loadp(Hmisc)
        markers = Hmisc::capitalize(tolower(markers))
    }
    
    obj = ScaleData(obj, features = markers, verbose = FALSE)
    p <- DotPlot(obj, group.by = group, features = markers) + coord_flip()
    p <- p + setText(18)
    p <- p + theme(panel.grid = element_blank(), axis.text.x = element_text(hjust = 1,vjust = 0.5))
    p <- p + labs(x = NULL, y = NULL)
    p <- p + scale_color_gradientn(values = seq(0, 1, 0.2), colours = c('#bebcbe', '#E71908'))
    p <- p + guides(color = guide_colorbar(barwidth = unit(4.5, "cm"), title = "Average Expression"))
    p <- p + theme(legend.position = "top")
    p <- p + theme(axis.text.x = element_text(size = 15 * 0.6, hjust = 0.5))
    save.graph(p, file = file, w, h, PDF = PDF, PNG = PNG)
}

cluster_start_with_one <- function(obj, cluster_name = "seurat_clusters"){
    if(sum(obj[[cluster_name]][, 1] == 0) == 0){
        cat("Maybe you have run this function before!\n")
    } else {
        cluster = as.numeric(as.character(obj[[cluster_name]][, 1])) + 1
        cluster = factor(cluster, levels = c(1:length(unique(cluster))))
        obj[[cluster_name]] = cluster
    }
    return(obj)        
}


prepareCounts <- function(outs){
    counts = Seurat::Read10X(outs)
    features = data.table::fread(paste0(outs, "/features.tsv.gz"), data.table = FALSE, header = FALSE)
    counts = suppressWarnings(cbind(features[, -3], counts))
    colnames(counts)[1:2] = c("GeneID", "GeneName")
    write("c", counts, file = paste0(outs, "/../counts.txt"))
    write("c", counts[1:100, 1:102], file = paste0(outs, "/../counts_100_x_100.txt"))
    invisible(counts)
}

FindClusters2 <- function(obj, cluster.range = c(35, 40), by = 0.1, res = 1, verbose = FALSE){
    if(verbose)
        cat2("Find suitable resolution, start with", res, "\n")

    if(length(cluster.range) == 1){
        cluster.range = c(cluster.range, cluster.range)
    }

    sigmoid <- function(x) 1 / (1 + exp(-x))
    x = -log(10/res - 1)
    plusCounter = minusCounter = 0
    nCluster = 1
    while(nCluster < cluster.range[1] | nCluster > cluster.range[2]){
        if(verbose)
            cat2("resolution", res, "... ")
        obj <- FindClusters(obj, resolution = res, verbose = FALSE)
        nCluster = length(unique(obj$seurat_clusters))

        if(nCluster < cluster.range[1]){
            x = x + by
            plusCounter = plusCounter + 1
        } else if (nCluster > cluster.range[2]){
            x = x - by
            minusCounter = minusCounter + 1
        } else{            
            break
        }
        res = signif(sigmoid(x) * 10, 3)
        if(plusCounter & minusCounter){
            cat2("\n")
            stop("Specific cluster ranger was skipped! Try expanding the cluster range or reducing the resolution step size.")
        }

        if(verbose)
            cat2(nCluster, "clusters.\n")
    }

    obj = cluster_start_with_one(obj)
    obj@misc[["best.resolution"]] = res
    if(verbose)
        cat2("Final resolution:", res, "with", nCluster, "clusters.\n")
    return(obj)
}

numToMax <- function(x){
    x = max(x)
    ceiling(x / 10 ^ floor(log10(x) - 1) + 4) * 10 ^ floor(log10(x) - 1)
}

RemoveDoublets <- function(
    obj,
    doublet.rate,
    sample_name,
    pN = 0.25,
    pc.num = 1:30,
    use.SCT = FALSE,
    remove = TRUE
  ){

    # ref:
    # - https://www.jianshu.com/p/6770c6a05287
    # - https://www.jianshu.com/p/0cb401b1ebe6
    # - https://zhuanlan.zhihu.com/p/469335625

    ## 寻找最优pK值
    sweep.res.list <- paramSweep_v3(obj, PCs = pc.num, sct = F)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
    bcmvn <- find.pK(sweep.stats)
    pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()

    ## 排除不能检出的同源doublets，优化期望的doublets数量
    homotypic.prop <- modelHomotypic(obj$seurat_clusters)   # 最好提供celltype
    nExp_poi <- round(doublet.rate*ncol(obj)) 
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    seu.scored <- doubletFinder_v3(obj, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, 
                                   nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
    # pick out doublets
    cname <-colnames(seu.scored[[]])
    DF <- cname[grep('^DF',cname)]
    seu.scored[["doublet"]] <- as.numeric(seu.scored[[DF]]=="Doublet")
    
    table(seu.scored@meta.data$doublet)
    # remove doublets
    if(remove == "TRUE"){
        seu.removed <- subset(seu.scored, subset = doublet != 1)
    } else{
        seu.removed <- seu.scored
    }

    p1 <- DimPlot(seu.scored, group.by = DF)
    res.list <- list("plot" = p1, "obj" = seu.removed)
    return(res.list)
}

# subclusterMakers <- list(
#     "T/I/NK cell"      = c("CD3D", "CD3", "CD8A", "CD8B", "CD4", "TRGC1", "TRGC2", "TRDC", "KRT86", "KIT", 
#                            "LST1", "FGFBP2", "TYROBP", "GNLY", "FCGR3A", "KLRF1", "MKI67", "STMN1"),
#     "Epithelial cell"  = c(),
#     "B cell"           = c('IL1B', 'PTGS2', 'VEGFA', 'FCN1', 'S100A8', 'S100A9', 'S100A12', 'CD1C', 
#                            'FCER1A', 'CLEC10A', 'LILRA4', 'CLIC3', 'CXCR3', 'CCL19', 'LAMP3', 'CCR7', 
#                            'CALHM6', 'VAMP5', 'LDHA', 'CHI3L1', 'SPP1', 'FBP1', 'MMP9', 'C1QC', 'FOLR2', 
#                            'SELENOP', 'C1QB', 'TMSB4X', 'YBX1', 'DNASE1L3', 'PCLAF', 'TYMS', 'UBE2C'),
#     "Myeloid cell"           = c('CD19', 'MS4A1', 'CD79A', 'CD79B', 'FCER2', 'TCL1A', 'IGHD', 'SELL', 'NR4A1', 
#                            'NR4A2', 'CD69', 'CXCR4', 'TOMM7', 'LTB', 'CD48', 'BANK1', 'RGS13', 'LMO2', 
#                            'NEIL1', 'BCL6', 'XBP1', 'MZB1', 'IGHA1', 'IGHA2', 'IGHG1', 'IGHG2', 'STMN1', 
#                            'MKI67', 'UBE2C', 'HMGB2', 'PCLAF'),
#     "Endothelial cell" = c('CXCL12', 'SULF1', 'SEMA3G', 'GJA5', 'IGFBP3', 'PLVAP', 'PLPP3', 'COL4A1', 
#                            'COL4A2', 'COL15A1', 'IL32', 'SYNE2', 'ARGLU1', 'WSB1', 'N4BP2L2', 'ACKR1', 
#                            'SOCS3', 'GADD45B', 'HLA−DRA', 'ZFP36', 'CCL21', 'LYVE1', 'PR', 'OX1', 'TFF3', 
#                            'AKAP12', 'CPE', 'MADCAM1', 'ADGRG6', 'CLDN5'),
#     "Fibroblast"       = c('MYH11', 'ACTA2', 'MCAM', 'CAV1', 'CXCL12', 'CFD', 'MFAP5', 'IGFBP6', 'DCN', 
#                            'RSPO3', 'WNT2B', 'CCL2', 'POSTN', 'WNT5A', 'COL3A1', 'TMEM158', 'CTHRC1', 
#                            'COL6A3', 'FAP')
#                         )

###### color for use
mycolors2 = c('#E4E4E4','#FFF7EC','#FEE8C8','#FDD49E','#FDBB84','#FC8D59','#EF6548','#D7301F','#B30000','#7F0000')

col.cluster2.1 <- c('#1C79B7','#F38329','#2DA248','#DC403E','#976BA6','#8D574C','#D07DB0','#BDBF3B', '#27BDD0','#B0C7E5','#F9BA79','#A0CC8A','#F09EA1','#C7B0D1','#D4BCB8','#F7CCDD',
                   '#DBDD8D','#A3D7E5','#B04E4F','#A38A58','#ED5351','#0C8945','#0F71A7','#A82764', '#F8DAE4','#7A4B1F','#5E6BAE','#8AC997','#DAC9AC','#0F5148','#A0B5DC','#9F858E',
                   '#5C181D','#7B8380','#E8DDDA','#264220','#5AB747','#5169AE','#4B3D58','#CD428A', '#62615A','#B82129','#66762E')
col.cluster <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                 '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                 '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')

llc_cols <- c('#b3d3e8','#7dcade','#368dab', '#259b39', '#65b72f', '#a2ce89', '#e4361e', '#93646d', '#c0a26e', '#edb000', '#eb6938', '#ed996d', '#ed7372', '#e969a2', '#d3b7d0', '#825796', '#cd0081', '#e3e39a', '#7491ca', '#b85919')

#### yellow-blue color
yellow_blue_color <- c('#450458', '#48196B', '#482878', '#472F7D', '#443A83', '#3E4C8A', '#375A8C', '#34608D', '#2C718E','#287C8E', '#20928C', '#1E9A8A', '#28AE80', '#3FBC72', '#60CA60', '#90D743', '#B2DD2D', '#CDE11D', '#E7E419','#EFE51C')

stage_cols <- c('#F3BA2F', '#E5916F', '#068085', '#A72F5C', '#FBE4B8', '#E42988', '#D8ACCF', '#21863A', '#F2CF1E', '#DADADA', '#1591C0', '#844482', '#293464', '#F4B2BC', '#AFDCE0', '#897664', '#E5481A', '#9FD08C', '#748595', '#EAE443')

creatSeuratObj  <- function(countFile, metaFile){

    counts = readCounts(countFile)
    meta = readMeta(metaFile)
    obnj = CreateSeuratObject(counts = counts, meta.data = meta)
}


FeaturePlot2  <- function(obj, reduction = "umap", features, axis = TRUE, miniAxis = FALSE, legend = TRUE){
    loadp(ggplot2)
    plotL = list()
    pdata = Embeddings(obj, reduction = reduction)
    pdata = data.frame(data_x = pdata[, 1], data_y = pdata[, 2])
    for(feature in features){
        pdata$Expression = log2(obj@assays$RNA@counts[feature ,] + 1)
        p <- ggplot(pdata, aes(data_x, data_y))
        p <- p + geom_point(aes(color = Expression), size = 0.3) + setText(20)
        p <- p + scale_color_continuous(low = "lightgrey", high = "#DE1F1F")
        p <- p + labs(title = feature)
        p <- p + guides(color = guide_colorbar(barwidth = unit(0.6, "cm"), barheight = unit(3, "cm")))

        if(!axis){
            p <- p + theme_blank()
        }

        if(!legend){
            p <- p + theme(legend.position = "none")
        }

        if(miniAxis){
            x.vector = pdata[, 1]
            y.vector = pdata[, 2]
            segTable = data.frame(x = c(min(x.vector) - 1.5 - 0.03, min(x.vector) - 1.5),
                                  y = c(min(y.vector) - 1.5, min(y.vector) - 1.5),
                                  vx = c(8, 0),
                                  vy = c(0, 8))
            p <- p + geom_segment(data = segTable, mapping = aes(x = x, y = y, xend = x + vx, yend = y + vy),
                                                                 arrow = arrow(length=unit(0.25, "cm")), 
                                                                 size = 0.5)
            p <- p + annotate("text", x = min(x.vector) + 1.1, y = min(y.vector) - 2.3, label = "UMAP_1", size = 4, hjust = 0, parse = TRUE)
            p <- p + annotate("text", x = min(x.vector) - 2.5, y = min(y.vector) + 1.4, label = "UMAP_2", size = 4, hjust = 0, parse = TRUE, angle = 90)
        }
        if(reduction == "umap"){
            p <- p + labs(x = "UMAP_1", y = "UMAP_2")
        } else {
            p <- p + labs(x = "tSNE_1", y = "tSNE_2")
        }
        plotL[[feature]] <- p
    }

    loadp(gridExtra)
    res <- grid.arrange(grobs = plotL, nrow = 1) 
    return(res)
}

getMultipletRate  <- function(obj){
    # https://kb.10xgenomics.com/hc/en-us/articles/360059124751-Why-is-the-multiplet-rate-different-for-the-Next-GEM-Single-Cell-3-LT-v3-1-assay-compared-to-other-single-cell-applications-
    library_size = seq(1:10) * 1000
    multiplet_rate = c(0.8, 1.5, 2.3, 3, 3.8, 4.6, 5.3, 6.1, 6.8, 8.0) / 100
    multiplet_rate_table = data.frame(library_size,multiplet_rate)

    library_size_obj = round(dim(obj)[2] / 1000)
    library_size_obj[library_size_obj > 10] = 10
    multiplet_rate_obj = multiplet_rate_table[library_size_obj, 2]
    return(multiplet_rate_obj)
}

mergeSingleCellranger <- function(dir, min.cells = 1, min.features = 1, method = "count"){
    sample = lf("x", dir)
    dirEach = switch(method,
                count = paste0(dir, "/", sample, "/outs/filtered_feature_bc_matrix"),
                multi = paste0(dir, "/", sample, "/outs/per_sample_outs/", sample, "/count/sample_filtered_feature_bc_matrix")
                )
    # check
    notExist = dirEach[!file.exists(paste0(dirEach, "/matrix.mtx.gz"))]
    if(length(notExist) > 0){
        message("Directory below does not exist or incomplete:")
        invisible(sapply(notExist, function(x) message("\t", x)))
        stop()
    }    
    names(dirEach) = gsub("_", "-", sample)
    counts = Read10X(data.dir = dirEach)
    obj = CreateSeuratObject(counts, min.cells = min.cells, min.features = min.features)
    return(obj)
}

mkdir <- function(...){
    folders = as.character(list(...))
    for(folder in folders){
        if(!dir.exists(folder))
            dir.create(folder, recursive = TRUE)
    }
}

multiPara <- function(func, # 要执行的函数 
                    ..., # func的动态参数
                    paraL = NULL, # func的静态参数
                    bind = "c", # 结果的合并方式
                    cores = NULL # 使用的核数
                    ){

    suppressPackageStartupMessages(library(doParallel))
    if(is.null(cores)){
        cores <- detectCores(logical = FALSE)
    }

    cl <- makeCluster(cores)
    registerDoParallel(cl)

    argL <- list(...)
    result <- suppressWarnings(foreach(i = seq(length(argL[[1]])), .combine = bind) %dopar% 
        do.call(func, c(lapply(argL, `[`, i), paraL)))
    stopCluster(cl)

    invisible(result)
}



processSigMarker <- function(object, i, n = 50, write = F) {
                                                            df <- FindMarkers(object, ident.1 = i, only.pos = T)
                                                            df.sig <- subset(df, df$p_val < 0.05)
                                                            df.sig <- df.sig[with(df.sig, order(p_val)),]
                                                            print(dim(df.sig))
                                    
                                                            if (write == T)
                                                            {
                                                                out <- cbind(gene_name = rownames(df.sig), df.sig)
                                                                filname = paste0(as.character(i), '_',as.character(length(out[,1])), 'sig_markers.txt')
                                                                write.table(out, file = filname, row.names = F, col.names = T, quote = F, sep='\t')
                                                            }

                                                            topsig <- rownames(df.sig[with(df.sig, order(p_val)),][1:n,])
                                                            return(topsig)
                                                            }

QC_MT <- function(object, n = 30)
{
    object <- PercentageFeatureSet(object, pattern = "^MT-", col.name = "percent.mt")
    temp <- table(cut(object$percent.mt, breaks = seq(0, 100, 5)))
    temp <- as.data.frame(temp)
    print(ggplot(temp, aes(x = Var1, y = Freq)) + geom_bar(position = "stack", stat = "identity", width = 0.6) + gg_qc)

    print(n)
    print(dim(object)[2])
    print(dim(subset(object, percent.mt < n))[2])
    print(dim(subset(object, percent.mt < n))[2]/dim(object)[2])

    return(object)
}

QC_nFnC <- function(object, n = 0.8)
{
    object$nFeature.nCount <- object$nFeature_RNA/object$nCount_RNA
    temp <- table(cut(object$nFeature.nCount, breaks = seq(0,1,0.05)))
    temp <- as.data.frame(temp)
    print(ggplot(temp, aes(x = Var1, y = Freq)) + geom_bar(position = "stack", stat = "identity", width = 0.6) + gg_qc)

    print(n)
    print(dim(object)[2])
    print(dim(subset(object, nFeature.nCount < n))[2])
    print(dim(subset(object, nFeature.nCount < n))[2]/dim(object)[2])
    
    return(object)
}

QC_Scat <- function(object, hyl =200, hyh = 7000, vx = 4e+4)
{
    print(FeatureScatter(object, feature1 = 'nCount_RNA', feature2 = "nFeature_RNA", pt.size = 0.01) + gg_style + theme_bw() +
          geom_hline(yintercept = hyl, color = "darkred", linetype = "dashed") + 
          geom_hline(yintercept = hyh, color = "darkred", linetype = "dashed") + 
          geom_vline(xintercept = vx, color = "darkblue", linetype = "dashed"))

    print(hyl);print(hyh)
    print(dim(object)[2])
    print(dim(subset(object, nFeature_RNA > hyl & nFeature_RNA < hyh))[2])
    cat(dim(subset(object, nFeature_RNA > hyl & nFeature_RNA < hyh))[2]/dim(object)[2])

    ### nCount hist 
    temp <- table(cut(object$nCount_RNA, breaks = seq(0, 2e+5, 1e+4)))
          temp <- as.data.frame(temp)
    print(ggplot(temp, aes(x = Var1, y = Freq)) + geom_bar(position = "stack", stat = "identity", width = 0.6) + gg_qc)
    print(vx)
    print(dim(subset(object, nCount_RNA < vx))[2])
    print(dim(subset(object, nCount_RNA < vx))[2]/dim(object)[2])
    
    return(object)
}

QC_workflow <- function(object, mt_n = 30, nFnC_n = 0.8, hyl =200, hyh = 7000, vx = 4e+4, method = 'check')
{
    if (method == 'check')
    {
        object <- QC_MT(object, mt_n)
        object <- QC_nFnC(object, nFnC_n)
        object <- QC_Scat(object, hyl, hyh, vx)

        return(object)
    }
    
    if (method == 'final')
    {
        object <- subset(object, percent.mt<mt_n & nFeature_RNA > hyl & nFeature_RNA < hyh & nCount_RNA < vx & nFeature.nCount < nFnC_n)
        print(object)
        return(object)  
    }
}

readCounts <- function(file){

    library(data.table)
    counts = fread(file, check.names = FALSE)
    row.names = as.character(data.frame(counts[, 1])[, 1])
    counts[, 1] = NULL
    counts = data.frame(counts, check.names = FALSE)
    rownames(counts) = row.names
    return(counts)
}

readfile2seurat <- function(file_path) # by fox

{
    loadp(data.table, Seurat)
    data <- fread(file_path)
    rownames(data) <- data$V1
    data$V1 <- NULL
    data <- CreateSeuratObject(counts = data, min.cell = 0)
    return(data)
}

anno_cell_type <- function(obj2, anno_list, name){
    
}

selectMarkers <- function(celltype = NULL){
    markers <- list(
        B              = c("MS4A1", "MZB1", "CD79A", "SDC1"), # SDC1通常被认为是一种上皮细胞和浆细胞特异性的蛋白，但它也可以在其他类型的细胞中表达
        T              = c("CD3D", "CD3E"),
        Platelets      = c("ITGA2B", "ITGB3", "GP1BA"),
        Endo           = c("VWF", "PECAM1", "ACKR1", "CLDN5", "ENG"),
        Epis           = c("EPCAM", "CD24", "CEACAM5", "KRT8"),
        Fibro          = c("COL3A1", "DCN", "ACTA2", "MYH11"),
        Mast           = c("TPSAB1", "CPA3", "KIT"), # 归入Myeloid
        Myeloid        = c("FCGR3A", "CD14", "CD68", "S100A8", "FCN1", "TPSAB1", "CPA3", "FLT3", "LAMP3"),
        Macro          = c("CD68", "CD14", "MARCO", "MRC1", "TNF"),
        mDC            = c("CD1C", "FCER1A", "CLIC3", "CXCR3", "CCL19", "CCR7"),
        Mono           = c("S100A8", "S100A9", "S100A12", "IL1B", "PTGS2", "FCN1"),
        Neutrophils    = c("FCGR3A", "FCGR3B", "ITGAM", "CEACAM8", "FUT4", "MPO", "LYZ", "CSF3R"),
        DC             = c("FLT3", "LILRA4", "XCR1", "CLEC9A", "CD1C", "CD1E", "LAMP3", "FCER1A"), # 归入Myeloid
        Smooth_muscle  = c("ACTA2", "TAGLN", "MUSTN1", "MYH11"),
        T_I_NK         = c("CD3D", "CD3E", "CD3G", "CD2"),
        MAIT           = c("SLC4A10", "TRAV1-2", "KLRB1"),
        gdT            = c("TRDV1", "TRDV2", "TRDC", "KIR2DL4"),
        ILC            = c("KRT86", "KIT", "LST1"),
        CD8            = c("CD8A", "CD8B"),
        CD4            = c("CD4"),
        NK             = c("NCAM1", "FCGR3A", "NKG2D", "NCR1", "GZMA", "GZMB"),
        Prolif         = c("MKI67", "STMN1", "UBE2C"),
        T1_Alveolar_cell         = c("CAV1", "PDPN"),
        T2_Alveolar_cell         = c("SFTPC", "AQP5", "NKX2-1", "LYZ"),
        Club_cells         = c("SCGB1A1", "CEACAM6"),
        Ciliated_Cell         = c("FOXJ1", "DNAH5")
        )

    if(is.null(celltype)){
        cat("Availabled celltype: \n\t")
        cat(paste(names(markers)), "\n")
    } else {
        if(! celltype %in% names(markers)){
            cat("Not availabled celltype! Only support below: \n\t")
            cat(paste(names(markers)), "\n")
        } else {        
            return(markers[[celltype]])
        }
    }
}

startSC <- function(lib = FALSE){

    if(lib) loadp(Seurat, ggplot2, dbplyr, RColorBrewer)
    print(paste('Start', as.character(date()), sep = '    '))
}


test_workflow <- function(object = NULL, nPC.method = 'sigPC', N = 40, markers.1st2check = general_marker, assay = 'RNA')
                            {
                                #### help
                                if(class(object) != 'Seurat')
                                {
                                    cat("  Help for Seurat Workflow 1\n
                                           test_workflow(object, nPC.method = 'sigPC', N = 40, markers.1st2check = general_marker, assay = 'RNA')\n 
                                           nPC.method:\n 
                                           sigPC: automaticaly choosing significant PCs by ScoreJackStraw\n
                                           topPC: mannually choosing Top N PCs\n
                                           markers.1st2check: general_marker by defalt
                                           \n
                                           example: lst <- test_workflow(seurat_object)\n
                                                    seurat_object <- lst[[1]]\n
                                                    lst\n
                                           example2:lst <- test_workflow(seurat_object, 'topPC', 30)\n
                                                    seurat_object <- lst[[1]]\n
                                                    lst\n\n")

                                } 
                                #### Main
                                if(class(object) == 'Seurat')
                                {
                                    print(paste0('Started to process seurat testing workflow with ', nPC.method, ' method'))    

                                    object <- NormalizeData(object)
                                    object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
                                    object <- ScaleData(object, features = rownames(object), verbose = F)
                                    object <- RunPCA(object)

                                    #### PCA Plots
                                    DimHeatmap(object, dims = 1:15, cells = 100, balanced = TRUE)
                                    pe <- ElbowPlot(object, ndims = 50)    

                                    #### significant PCs selected by jackstraw score : nPC     
                                    object <- JackStraw(object, dims = 50)
                                    object <- ScoreJackStraw(object, dims = 1:50)
                                    pj <- JackStrawPlot(object, dims = 1:50)
                                    JackResult <- object@reductions$pca@jackstraw@overall.p.values

                                    #### nPC.method
                                    if (nPC.method == 'sigPC')
                                    {
                                        nPC <- JackResult[,'PC'][which(JackResult[,'Score'] < 0.05)]
                                    } else if (nPC.method == 'topPC')
                                    {
                                        nPC <- N
                                    }

                                    #### Run TSNE UMAP
                                    object <- RunTSNE(object, reduction = "pca", dims = nPC, seed.use = 7, verbose = F)
                                    object <- RunUMAP(object, reduction = "pca", dims = nPC, verbose = F)
                                    
                                    #### orig.ident
                                    pumap <- DimPlot(object, group.by = 'orig.ident', cols = 'grey', pt.size = 0.1) + gg_style
                                    ptsne <- DimPlot(object, reduction = 'tsne', group.by = 'orig.ident', cols = 'grey', pt.size=0.1) + gg_style

                                    #### Run SNN
                                    object <- FindNeighbors(object, reduction = "pca", dims = nPC, verbose = F)
                                    object <- FindClusters(object, resolution = seq(0.1, 1.2, by=0.1), verbose = F)
                                    
                                    if (assay == 'RNA')
                                    {
                                        p0.7umap <- DimPlot(object, group.by = 'RNA_snn_res.0.7', label=T, pt.size=0.1) + gg_style
                                        p0.7tsne <- DimPlot(object, reduction = 'tsne', group.by = 'RNA_snn_res.0.7', label=T, pt.size=0.01) + gg_style
                                    }

                                    pf <- FeaturePlot(object, reduction = 'tsne', features = markers.1st2check,  pt.size = 0.1, cols = mycolors2)
                                    
                                    ### Return Out list
                                    out_list <- list('object' = object, 'elbow.plot' = pe, 'jack.plot' = pj, 'tsne.plot'= ptsne, 'umap.plot' = pumap, 
                                                     'umap0.7' = p0.7umap, 'tsne0.7' = p0.7tsne, 'general.feature.plot' = pf, 'nPC' = nPC)
                                    
                                    print('seurat testing workflow completed')
                                    
                                    return(out_list)
                                }
                            }

subString <- function (strings, idx, sep = NULL, rev = FALSE, collapse = NULL){
    strings = as.character(strings)
    if (is.null(sep)) {
        if (rev) {
            res = as.character(sapply(strings, function(x) paste(rev(rev(strsplit(x, 
                "")[[1]])[idx]), collapse = "")))
        }
        else {
            res = as.character(sapply(strings, function(x) paste(strsplit(x, 
                "")[[1]][idx], collapse = "")))
        }
    }
    else {
        if (rev) {
            res = sapply(strsplit(strings, sep), function(x) paste(rev(rev(x)[idx]), 
                collapse = collapse))
        }
        else {
            res = sapply(strsplit(strings, sep), function(x) paste(x[idx], 
                collapse = collapse))
        }
    }
    return(res)
}

read_GPL <- function(GPL_file){
    nskip = sum(subString(readr::read_lines(GPL_file, n_max = 100), 1) == "#")
    anno_tb = data.table::fread(GPL_file, skip = nskip, data.table = FALSE)
    return(anno_tb)
}

anno_NGS <- function(expM, species = "human"){
    if(species == "human"){
        OrgDb = "org.Hs.eg.db"
    } else if(species == "mouse"){
        OrgDb = "org.Mm.eg.db"
    }

    loadp(clusterProfiler)
    library(OrgDb, character.only = TRUE)

    convert_tb = suppressWarnings(bitr(row.names(expM), fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = OrgDb))
    SYMBOL = convert_tb$SYMBOL[match(rownames(expM), convert_tb$ENSEMBL)]
    idx_NA = is.na(SYMBOL)
    expM = expM[!idx_NA, ]
    SYMBOL = SYMBOL[!idx_NA]
    expM = cbind(SYMBOL, expM)
    expM = aggregate(. ~ SYMBOL, data = expM, function(x) x[which.max(x)])
    expM = col2rownames(expM)
    return(expM)
}

anno_GEO_GPL <- function(expM, GPL, GPL_dir, method = "max", symbol_name = NULL){
  files = lf("rf", GPL_dir)  
  xfile = files[subString(basename(files), 1, "[-\\._]") == GPL][1]
  file_type = subString(basename(xfile), 1, "\\.", rev = TRUE)
  if(file_type == "txt"){
      anno_tb = read_GPL(xfile)
  } else {
    anno_tb = data.table::fread(xfile, skip = "#", data.table = FALSE)
  }

    if(is.null(symbol_name)){
        SYMBOL_col_names = c("Gene symbol", "GENE_SYMBOL", "GeneSymbol", "Symbol", "SYMBOL", "gene_assignment")
        symbol_name = SYMBOL_col_names[which(SYMBOL_col_names %in% colnames(anno_tb))]
        if(length(symbol_name) == 0){
        cat2("Can't find right 'symbol' column name. If exists exactly, please specific it to 'symbol_name'\n")
        stop()          
        }
    }      

    anno_tb = anno_tb[, c("ID", symbol_name)]  
    colnames(anno_tb) = c("ID", "SYMBOL")
    anno_tb$SYMBOL = subString(anno_tb$SYMBOL, 1, "[ /|]")
    anno_tb$SYMBOL[anno_tb$SYMBOL == "NA"] = ""
    expM_anno = anno_GEO(expM, anno_tb, method = method)
  return(expM_anno)
}

date_to_gene <- function(tdata, convert_file, reoder = TRUE, verbose = TRUE){
    convert_tb = read.table(convert_file, header = TRUE, sep = "\t")
    if(class(tdata)[1] == "character"){
        tdata = convert_tb$gene[match(tdata, convert_tb$date2)]
        cat2("Warning:", paste0(sum(is.na(tdata)), "/", length(tdata), " genes falied to convert!\n"))
    } else if(class(tdata)[1] %in% c("data.frame", "matrix")){
        idx = which(rownames(tdata) %in% convert_tb$date2)
        if(length(idx) == 0){
            cat2("Warning: all", nrow(tdata), "genes falied to convert!\n")
        }
        gene_names = convert_tb$gene[match(rownames(tdata)[idx], convert_tb$date2)]
        if(verbose) cat2(length(idx), "genes converted!\n")
        rownames(tdata)[idx] = gene_names
    }
    if(reoder) tdata = tdata[order(rownames(tdata)), ]
    return(tdata)
}

get_GEO_pheno <- function(GSE, out_dir = "00_GEO_data"){
    obj = getGEO2(GSE, out_dir)
    pheno = obj$sample_anno
    return(pheno)
}

getGEO2 <- function(GSE, out_dir = "00_GEO_data"){
    report(paste("Querying dataset:", GSE, "..."))
    if(! "GEOquery" %in% rownames(installed.packages())){
        stop()
        cat("Packages 'GEOquery' not installed! to install run 'BiocManager::install(\"GEOquery\")'")
    }

    loadp(GEOquery)
    
    mkdir(out_dir)
    options(timeout = 1000000000)
    N = Ntry = check = 10
    class(check) = "try-error"
    while(class(check) == "try-error" & Ntry > 0){
        cat2("\tTry time", 10 - Ntry + 1, "... ")
        check = try(suppressMessages(getGEO(GSE, destdir = out_dir, getGPL = FALSE)), silent = TRUE)
        Ntry = Ntry - 1
        if(class(check) == "try-error"){
            cat2("failed.\n")
        }
    }
    
    if(class(check) == "try-error"){
        message("Download failed, tried ", N, " times!")
        suppressMessages(getGEO(GSE, destdir = out_dir, getGPL = FALSE))
    } else {
        cat2("succeed.\n")
        eSet = check
    }
    
    xeset = eSet[[1]]
    GPL = xeset@annotation
    expM_raw = exprs(xeset)
    sample_anno = pData(xeset)

    # report
    cat2("- GPL:", GPL, "\n")
    cat2("- Found", nrow(sample_anno), "samples,", ncol(sample_anno), "metas.\n")
    cat2("- Writting sample_anno to", paste0(out_dir, "/", GSE, "_sample_anno.csv"), "\n")
    write.csv(sample_anno, file = paste0(out_dir, "/", GSE, "_sample_anno.csv"))
    if(class(try(suppressWarnings(library(xlsx)), silent = T)) != "try-error")
        if(R.version$os != "linux-gnu")
            table2xlsx(sample_anno, file = paste0(out_dir, "/", GSE, "_sample_anno.xlsx"))

    if(nrow(expM_raw) == 0){
        cat2("Can't found expression profile, please see\n")
        cat2(paste0("\thttps://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", GSE), "\n")
        assign(GSE, list(GPL = GPL, sample_anno = sample_anno))
    } else {
        # 自动判断是否已经log化
        qx = as.numeric(quantile(expM_raw, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
        LogC = ( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
        if(LogC){
            report("Execute log(expM + 1) ...")
            expM_raw[expM_raw < 0] = 0
            expM_raw = log2(expM_raw + 1)
        }

        # Normlization
        report("Normalize between arrays ...")
        expM_norm = limma::normalizeBetweenArrays(expM_raw)
        pdf(paste0(out_dir, "/", GSE, "_Normalization.pdf"), 8, 4)
            par(mfcol = c(1, 2))
            boxplot(expM_raw, outline = FALSE, notch = T)
            boxplot(expM_norm, outline = FALSE, notch = T)
        dev.off()

        assign(GSE, list(GPL = GPL, expM_raw = expM_raw, expM = expM_norm, sample_anno = sample_anno))
        save(list = GSE, file = paste0(out_dir, "/", GSE, ".RData"))
        
        cat2("- Successed, file save to", paste0(out_dir, "/",GSE, ".RData.\n"))        
    }
    invisible(get(GSE))
}

anno_GEO <- function(expM, anno_tb, method = "max"){

    choices = c("mean", "max", "remove", "sum")
    match.arg(method, choices)

    # check match
    n_not_match = sum(! rownames(expM) %in% anno_tb$ID)
    if(n_not_match > 0){
        cat2("- Warnnig:", paste0(n_not_match, "/", nrow(expM)), "porbes fail to matched!\n")
    } else {
        cat2("- All porbes matched!\n")
    }

    # check SYMBOL
    gene_name = anno_tb$SYMBOL[match(rownames(expM),  anno_tb$ID)]
    n_no_symbol = sum(gene_name == "", na.rm = TRUE)
    if(n_no_symbol > 0){
        cat2("- Warnning:", paste0(n_no_symbol, "/", nrow(expM)), "probes fail to annotated!\n")
    } else {
        cat2("- All porbes annotated!\n")
    }

    expM = data.frame(gene_name, expM)
    expM2 = expM[expM$gene_name != "", ]

    report(paste("Removing duplicated genes by method:", method, "..."))
    if(method == "mean"){
        expM3 = aggregate(. ~ gene_name, data = expM2, mean0)
    } else if(method == "max"){
        expM3 = aggregate(. ~ gene_name, data = expM2, function(x) x[which.max(x)])
    } else if(method == "sum"){        
        expM3 = aggregate(. ~ gene_name, data = expM2, sum)
    } else if(method == "remove"){
        expM3 = expM2[!duplicated(expM2[, 1]), ]
        expM3 = expM3[order(expM3[, 1]), ]
    }
    rownames(expM3) = expM3[, 1]
    expM3[, 1] = NULL
    return(expM3)
}

anno_GEO_local <- function(expM, GPL, anno_file, method = "mean"){

    report("Loadding annotation file ...")
    anno_obj = get(load(anno_file))

    # check GPL
    if(! GPL %in% names(anno_obj)){
        cat2("Error: Not found platform",  GPL, "in annotation file!\n")
        stop()
    }
    anno_tb = anno_obj[[GPL]]
    expM_anno = anno_GEO(expM, anno_tb)
    return(expM_anno)
}

anno_GEO_AnnoProbe <- function(expM, GPL, method = "mean"){
    # idmap包可以被 AnnoProbe 包替代
    # if(sum(c("idmap1", "idmap2", "idmap3") %in% data.frame(installed.packages())$Package) != 3){
    #   cat2("Please install annotation packages 'idmap1' 'idmap2' 'idmap3' from github.\n")
    #   stop()
    # }

    # loadp(idmap1, idmap2, idmap3)
    # anno_tb1 = try(suppressMessages(getIDs(GPL)), silent = TRUE)
    # anno_tb2 = try(suppressMessages(get_soft_IDs(GPL)), silent = TRUE)
    # anno_tb3 = try(suppressMessages(get_pipe_IDs(GPL)), silent = TRUE)
    # if(class(anno_tb1) == "data.frame"){
    #   anno_tb = anno_tb1
    # } else if(class(anno_tb2) == "data.frame"){
    #   anno_tb = anno_tb2
    # } else if(class(anno_tb3) == "data.frame"){
    #   anno_tb = anno_tb3
    # } else {
    #   cat2("Can't find annotation table, please DIY one.\n")
    #   stop()      
    # }

    if(! "AnnoProbe" %in% data.frame(installed.packages())$Package){
        cat2("Please install annotation packages 'AnnoProbe' from github.\n")
        cat2("\tdevtools::install_github(\"jmzeng1314/AnnoProbe\").\n")
        stop()
    }
    
    report("Use AnnoProbe annotation file ...")
    loadp(AnnoProbe)
    idmap <- function (gpl = "GPL570", type = "bioc", mirror = "tercent"){
        gpl = toupper(gpl)
        gpl_anno = paste(gpl, c("bioc", "soft", "pipe"), sep = "_")
        if (mirror == "tercent") {
            up = "http://49.235.27.111"
        }
        if (!checkGPL(gpl)) {
            stop("This platform is not in our list, please use our shinyAPP to custom annotate your probe sequences, or ask us to process and then update the R package!")
        }
        else {
            tryCatch(utils::data("exists_anno_list", package = "AnnoProbe"))
            gpl_anno = gpl_anno[gpl_anno %in% exists_anno_list]
            if (T) {
                tpf = paste0(paste(gpl, type, sep = "_"), ".rda")
                down = paste0("/GEOmirror/GPL/", tpf)
                download.file(paste0(up, down), tpf, mode = "wb", quiet = TRUE)
                load(tpf)
                return(get(paste(gpl, type, sep = "_")))
            }
            else {
                stop("We have that platform, but just offer other type of annotaion.")
            }
        }
    }

    anno_tb = try(suppressMessages(idmap(GPL, type = "pipe")), silent = TRUE)
    if(class(anno_tb) == "try-error"){
        cat2("- Can't find annotation table in AnnoProbe.\n")
        stop()
    }
    colnames(anno_tb) = c("ID", "SYMBOL")
    anno_tb = anno_tb[!duplicated(anno_tb$ID), ]
    anno_tb$SYMBOL = subString(anno_tb$SYMBOL, 1, "[ /|]")
    anno_tb$SYMBOL[anno_tb$SYMBOL == "NA"] = ""

    expM_anno = anno_GEO(expM, anno_tb)
    return(expM_anno)
}

anno_GEO_online <- function(expM, GPL, out_dir = "annot_dir", symbol_name = NULL){
    mkdir(out_dir)
    GPL_num = gsub("GPL", "", GPL)
    master_string = ifelse(nchar(GPL_num) < 4, "", subString(GPL_num, 1:(nchar(GPL_num) - 3)))
    master_dir = paste0("GPL", master_string, "nnn")
    url_prefix = "https://ftp.ncbi.nlm.nih.gov/geo/platforms"
    soft_url = paste0(url_prefix, "/", master_dir, "/", GPL, "/annot/", GPL, ".annot.gz")    
    soft_file = paste0(out_dir, "/", GPL, ".annot.gz")
    report("Download online annotation file ...")
    options(timeout = 10000)
    dl_res = try(download.file(soft_url, destfile = soft_file, quiet = TRUE), silent = TRUE)
    if(class(dl_res) == "try-error"){
        cat2("Can't find online annotation table. Please download annotation table manually:\n")
        cat2(paste0("\thttps://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", GPL, "\n"))
        stop()
    }
    anno_tb = suppressWarnings(data.table::fread(soft_file, data.table = FALSE))
    
    if(is.null(symbol_name)){
        SYMBOL_col_names = c("Gene symbol", "GENE_SYMBOL", "GeneSymbol", "Symbol", "SYMBOL", "gene_assignment")
        symbol_name = SYMBOL_col_names[which(SYMBOL_col_names %in% colnames(anno_tb))]
        if(length(symbol_name) == 0){
            cat2("Can't find right 'symbol' column name. If exists exactly, please specific it to 'symbol_name'\n")
            stop()          
        }
    }      
    anno_tb = anno_tb[, c("ID", symbol_name)]  
    colnames(anno_tb) = c("ID", "SYMBOL")
    anno_tb$SYMBOL = subString(anno_tb$SYMBOL, 1, "[ /|]")
    anno_tb$SYMBOL[anno_tb$SYMBOL == "NA"] = ""
    expM_anno = anno_GEO(expM, anno_tb)
    return(expM_anno)
}

cat2 <- function(...){
    cat(...)
    flush.console()
}

intersect2 <- function(...) Reduce(intersect, list(...))

run_ML_Boruta <- function(expM, Group, out_name = "Boruta", lambda = "min", seed = 1234, maxRuns = 100, pValue = 0.01){
    if(! class(expM)[1] %in% c("data.frame", "matrix")){
        stop("'expM' must be a data.frame or matrix!\n")
    }

    loadp(Boruta)
    
    if(class(Group) != "numeric"){
        stop("Group must be numeric!\n")
    }
    
    if(ncol(expM) != length(Group)){
        stop("Length of expM and Group must be same!\n")
    }
    
    x = t(expM)
    y = factor(Group)

    set.seed(seed)
    data_Boruta = data.frame(cbind(x, group = y))
    boruta = suppressMessages(Boruta(group ~ ., data = data_Boruta, doTrace = 2, maxRuns = maxRuns, pValue = pValue))

    lz = lapply(1:ncol(boruta$ImpHistory), function(i)
        boruta$ImpHistory[is.finite(boruta$ImpHistory[,i]),i])
    names(lz) = colnames(boruta$ImpHistory)
    Labels = sort(sapply(lz,median))

    pdf(paste0(out_name, "_01_Importance.pdf"), 6, 5)
        plot(boruta, xlab = "", xaxt = "n")
        axis(side = 1,las=2,labels = names(Labels), at = 1:ncol(boruta$ImpHistory), cex.axis = 0.7)
    dev.off()


    pdf(paste0(out_name, "_02_Iteration_progress.pdf"), 6, 5)
        plotImpHistory(boruta)
    dev.off()

    final.boruta = TentativeRoughFix(boruta)
    names(Labels) = gsub("\\.", "-", names(Labels))
    write.csv(Labels, file = paste0(out_name, "_03_Labels.csv"))

    Boruta_gene = getSelectedAttributes(final.boruta, withTentative = F)
    Boruta_gene = gsub("\\.", "-", Boruta_gene)
    write.csv(Boruta_gene, file = paste0(out_name, "_04_result.csv"))
    cat2("Get", length(Boruta_gene), "genes.\n")

    boruta.df = attStats(final.boruta)
    rownames(boruta.df) = gsub("\\.", "-", rownames(boruta.df))

    result_Boruta = print(boruta.df)
    write.csv(result_Boruta, file = paste0(out_name, "_05_Boruta_result_all.csv"))

    write.line(Boruta_gene, file = paste0(out_name, "_06_Boruta_genes.txt"))
    return(Boruta_gene)
}

calculate_percentile_point <- function(A, B, percentile) {
    distance <- abs(B - A)
    percentile_distance <- distance / percentile
    percentile_point <- B - percentile_distance
    return(percentile_point)
}

run_ML_Lasso <- function(expM, Group, out_name = "Lasso", lambda = "min", seed = 1234, nfolds = 10, legend = TRUE){

    loadp(glmnet, ggplot2)

    if(! class(expM)[1] %in% c("data.frame", "matrix")){
        stop("'expM' must be a data.frame or matrix!\n")
    }

    if(class(Group) != "numeric"){
        stop("Group must be numeric!\n")
    }
    
    if(ncol(expM) != length(Group)){
        stop("Length of expM and Group must be same!\n")
    }
    
  
    set.seed(1234)
    x = t(expM)
    y = factor(Group)
    fit = glmnet(x, y, family = "binomial", nlambda = 100, alpha = 1) 
    cvfit = cv.glmnet(x, y, family = "binomial", maxit = 8000, nfolds = nfolds)
    best_lamda = signif(log(cvfit$lambda.min), 6)

    pdata = broom::tidy(fit)
    pdata = pdata[pdata$term != "(Intercept)", ]
    if(length(unique(pdata$term)) > 10) legend =  FALSE
    text_y = calculate_percentile_point(min(pdata$estimate), max(pdata$estimate), 4)
    p <- ggplot(pdata, aes(log(lambda), estimate, group = term, color = term))
    p <- p + geom_line(size = 0.5)
    p <- p + geom_vline(xintercept = best_lamda, color = "grey60", alpha = 0.8, linetype = 2)
    p <- p + geom_hline(yintercept = 0)
    p <- p + ylab("Coefficients")
    p <- p + scale_x_continuous(expand = c(0.02, 0.02))
    p <- p + scale_color_manual(name = "Genes", values = col.cluster2.1)
    p <- p + setText(20, theme = "bw")
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    if(legend){
        p <- p + theme(legend.position = "top")
    } else {
        p <- p + theme(legend.position = "none")
    }
    p <- p + annotate("text", x = best_lamda, y = text_y, label = paste("Optimal Lambda =", signif(cvfit$lambda.min, 4)))
    print(p)    
    save.graph(p, file = paste0(out_name, "_01_lambda"), 7, 4.5)

    pdata_cv = broom::tidy(cvfit)
    text_y = calculate_percentile_point(min(pdata_cv$estimate), max(pdata_cv$estimate), 4)
    p <- ggplot()
    p <- p + geom_point(data = pdata_cv, aes(log(lambda), estimate))
    p <- p + geom_errorbar(data = pdata_cv, aes(x = log(lambda), ymin = conf.low, ymax = conf.high))
    p <- p + ylab("Coefficients")
    p <- p + scale_x_continuous(expand = c(0.02, 0.02))
    p <- p + setText(20, theme = "bw")
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + annotate("text", x = best_lamda, y = text_y, label = paste("Optimal Lambda =", signif(cvfit$lambda.min, 4)))
    p <- p + geom_vline(xintercept = best_lamda, color = "grey60", alpha = 0.8, linetype = 2)
    print(p)    
    save.graph(p, file = paste0(out_name, "_02_cvfit"), 7, 4.5)

    coef_min <- coef(cvfit$glmnet.fit, s = cvfit$lambda.min, exact = F)
    genes_lasso_min = rownames(coef_min)[coef_min[, 1] != 0][-1]
    cat2("Get", length(genes_lasso_min), "genes with lamda 'min'.\n")

    coef_1se <- coef(cvfit$glmnet.fit, s = cvfit$lambda.1se, exact = F)
    genes_lasso_1se = rownames(coef_1se)[coef_1se[, 1] != 0][-1]
    cat2("Get", length(genes_lasso_1se), "genes with lamda '1se'.\n")
    genes_lasso_min

    if(lambda == "min"){
        cat2("Finally, use default 'min' lamda as input of coef\n")
        cat2("\tOr, you can set 'lamda' as '1se' and rerun.\n")
        final_genes = genes_lasso_min
    } else {
        cat2("Finally, use '1se' lamda as input of coef\n")
        final_genes = genes_lasso_1se
    }
    write.line(final_genes, file = paste0(out_name, "_03_LASSO_genes.txt"))
    return(final_genes)        
}


venn <- function (..., out_name = NULL, PPT = NULL, label = NULL, 
    label.cex = 8, text.cex = 8, fill.alpha = 1, 
    fill.col = c("#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0"), 
    stroke.alpha = 1, stroke.size = 1, mar = 1, 
    show.elements = FALSE, label.sep = "\n", auto.scale = FALSE, 
    x.lim = NULL, y.lim = NULL, add_border = FALSE, hjust = NULL, stroke.linetype = "longdash",
    col.stroke = "black", show.percentage = FALSE, w = 7, h = 5) 
{
    data_list = list(...)
    if(is.null(label)){
        names(data_list) = getName(...)
    } else {
        names(data_list) = label
    }

    Len = length(data_list)
    vennM = matrix("", Len, Len, dimnames = list(getName(...), 
        getName(...)))
    for (i in 1:Len) {
        for (j in i:Len) {
            vennM[i, j] = length(intersect(data_list[[i]], data_list[[j]]))
        }
    }
    resM = apply(vennM, 1, as.numeric)
    dimnames(resM) = dimnames(vennM)
    if (Len > 4) {
        message("More than 4 objects, will not plot venn.")
        return(resM)
    }

    loadp(ggvenn)
    p <- suppressWarnings(
        ggvenn(data = data_list, 
               show_percentage = show.percentage,
               show_elements   = show.elements, 
               label_sep       = label.sep, 
               auto_scale      = auto.scale,
               set_name_size   = label.cex,
               text_size       = text.cex, 
               stroke_alpha    = stroke.alpha,
               stroke_size     = stroke.size, 
               fill_alpha      = fill.alpha, 
               fill_color      = fill.col,
               stroke_linetype = stroke.linetype
        ))
    
    if(is.null(x.lim)){
        if(length(data_list) == 4){
            nchar_left = nchar(names(data_list))[1]
            nchar_right = nchar(names(data_list))[4]
            x.lim = c(-1 - 0.25 * nchar_left, 1 + 0.25 * nchar_right)
        } else {
            x.lim = c(-1.8, 1.8)        
        }
    }

    if(is.null(y.lim)){
        if(length(data_list) == 2){
            y.lim = c(-1.2, 1.5)
        } else if(length(data_list) == 3){
            y.lim = c(-2.1, 2)        
        } else {
            y.lim = c(-1.8, 1.5)
        }
    }
   
    if(add_border) p <- p + theme_bw()
    if(!is.null(hjust)) p$layers[[3]]$data$hjust = hjust
    p <- p + scale_x_discrete(expand = c(0.05, 0.05))
    p <- p + xlim(x.lim) + ylim(y.lim)

    if(!is.null(out_name)) save.graph(p, file = out_name, w, h)


    # if(length(data_list) == 2){
    #     h = h * 0.75
    # }

    if (!is.null(PPT)) {
        require(export)
        graph2ppt(p, file = paste0(out_name, ".pptx"), width = w, height = h)
    }
    return(p)
}

run_PCA <- function(expM, Group, out_name = NULL, text.size = 20, pt.size = 1, 
    ellipse = TRUE, w = 5, h = 5, legend.ncol = NULL, color = c("firebrick3", "steelblue")
    ){
    expM = expM[apply(expM, 1, sd) != 0, ]
    pca = prcomp(t(expM), scale = TRUE)
    pca.data = data.frame(pca$x)

    loadp(ggplot2)
    p <- ggplot(pca.data, aes(x = PC1, y = PC2, color = Group))
    p <- p + geom_point(size = pt.size)
    p <- p + geom_hline(yintercept = 0, lty = 4, col = "black", lwd = 0.6)
    p <- p + geom_vline(xintercept = 0, lty = 4, col = "black", lwd = 0.6)
    p <- p + guides(color = guide_legend(override.aes = list(size = 2)))
    p <- p + setText(text.size, theme = "bw") + theme(legend.position = "top")
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + scale_color_manual(values = c("firebrick3", "steelblue"))
    if(ellipse) p <- p + stat_ellipse(aes(x = PC1, y = PC2), linetype = 2, size = 0.5, level = 0.95)
    if(!is.null(legend.ncol))
        p <- p + guides(color = guide_legend(ncol = legend.ncol))
    if(!is.null(out_name)) save.graph(p, file = out_name, w, h)
    return(p)
}

run_DEG_limma <- function(expM, Group, Case, Control, out_name = "DEG", log2FC = 1, pvalue = 0.05, padj = TRUE){

    loadp(limma)
    if(ncol(expM) != length(Group)) stop("Length not matched between Matrix and Group!")
    if(! Case %in% Group) stop(paste0("Not found Group name: ", Case, "."))
    if(! Control %in% Group) stop(paste0("Not found Group name: ", Control, "."))

    # subset
    expM = expM[apply(expM, 1, sd) != 0, ]
    idx = Group %in% c(Case, Control)
    expM = expM[, idx]
    Group = Group[idx]

    design <- model.matrix(~ 0 + factor(Group))
    colnames(design) = levels(factor(Group))
    rownames(design) = colnames(expM)
 
    code = paste0("makeContrasts(", Case, "-", Control, ", levels = design)")
    contrast.matrix = eval(parse(text = code)) 

    fit = lmFit(expM, design)
    fit2 = contrasts.fit(fit, contrast.matrix)
    fit2 = eBayes(fit2)
    DEG_tb = topTable(fit2, coef = 1, adjust = 'fdr', number = Inf)

    DEG = rep("Not_Change", nrow(DEG_tb))
    if(padj){
        DEG[DEG_tb$adj.P.Val < pvalue & DEG_tb$logFC > log2FC]  = "Up"
        DEG[DEG_tb$adj.P.Val < pvalue & DEG_tb$logFC < -log2FC] = "Down"
    } else {
        DEG[DEG_tb$P.Value < pvalue & DEG_tb$logFC > log2FC]  = "Up"
        DEG[DEG_tb$P.Value < pvalue & DEG_tb$logFC < -log2FC] = "Down"
    }
    DEG_tb$DEG = DEG
    DEG_tb = DEG_tb[order(DEG_tb$logFC, decreasing = TRUE), ]
    # DEG_tb = sortBy(DEG_tb, by = c(1, 7), decreasing  = c(TRUE, TRUE))

    # report
    cat2(paste("Number of Up regulated genes:", length(DEG[DEG == "Up"])), "\n")
    cat2(paste("Number of Down regulated genes:", length(DEG[DEG == "Down"])), "\n")
    cat2(paste("Number of Not Change regulated genes:", length(DEG[DEG == "Not_Change"])), "\n")

    if(!is.null(out_name)){
        save(DEG_tb, file = paste0(out_name, ".RData"))
        write.csv(DEG_tb, file = paste0(out_name, ".csv"))
    }        
    return(DEG_tb)
}

plot_volcano_limma <- function(
    DEG_tb, log2FC = 1, pvalue = 0.05, title = "", xtitle = "log2(Fold Change)", alpha = 0.5,
    padj = TRUE, text.cex = 20, label.cex = 2, pt.size = 0.8, label = TRUE, label.top = 10,
    w = 8, h = 7, out_name = NULL, max.overlaps = Inf, method = "pvalue", caption = FALSE ,
    xlims = NULL, ylims = NULL
    ){
    # DEG_tb = na.omit(DEG_tb)
    if(is.null(xlims)){
        x_min = floor(min(DEG_tb$logFC))
        x_max = ceiling(max(DEG_tb$logFC))
    } else {
        x_min = xlims[1]
        x_max = xlims[2]
    }

    loadp(ggplot2, ggrepel)
    if(padj){
        p <- ggplot(DEG_tb, aes(x = logFC, y = -log10(adj.P.Val), fill = DEG, color = DEG))
        if(is.null(ylims)) ylims = ceiling(max(-log10(DEG_tb$adj.P.Val)))
        caption_text = paste0("adj.P.Val < ", pvalue, " & |log2FC| > ", log2FC, "\n")
        p <- p + labs(title = title, x = xtitle, y = "-log10(FDR)")
    } else {
        p <- ggplot(DEG_tb, aes(x = logFC, y = -log10(P.Value), fill = DEG, color = DEG))
        if(is.null(ylims)) ylims = ceiling(max(-log10(DEG_tb$P.Value)))
        caption_text = paste0("P.Value < ", pvalue, " & |log2FC| > ", log2FC, "\n")
        p <- p + labs(title = title, x = xtitle, y = bquote(.("-log10(") * italic("P") * .("-value)")))
    }

    if(caption){
        p <- p + labs(subtitle = caption_text)
        p <- p + theme(plot.subtitle = element_text(size = text.cex * 0.3, hjust = 1, color = "black", face = "italic"))
    } 
    
    p <- p + geom_point(size = pt.size, alpha = alpha, shape = 21, stroke = 1)
    p <- p + xlim(x_min, x_max) + ylim(0, ylims)
    p <- p + setText(text.cex, theme = "bw", y.title.vjust = 2)
    p <- p + theme(plot.title = element_text(hjust = 0.5))
    p <- p + theme(legend.position = "top")

    N_UP = format(sum(DEG_tb$DEG == "Up"), big.mark = ",")
    N_NC = format(sum(DEG_tb$DEG == "Not_Change"), big.mark = ",")
    N_DN = format(sum(DEG_tb$DEG == "Down"), big.mark = ",")
    p <- p + scale_fill_manual(values = c(Up = "red", Down = "dodgerblue4", Not_Change = "grey"),
                                name = "", 
                                breaks = c("Down", "Not_Change", "Up"),
                                labels = paste0(c("Down", "Not Change", "Up"), ": ", c(N_DN, N_NC, N_UP)))
    p <- p + scale_color_manual(values = c(Up = "red", Down = "dodgerblue4", Not_Change = "grey"),
                                name = "", 
                                breaks = c("Down", "Not_Change", "Up"),
                                labels = paste0(c("Down", "Not Change", "Up"), ": ", c(N_DN, N_NC, N_UP)))
    p <- p + theme(legend.text = element_text(size = text.cex * 0.6))
    p <- p + guides(color = guide_legend(override.aes = list(size = text.cex * 0.1)))
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + geom_vline(xintercept = c(log2FC, -log2FC), lty = 4, col = "black", lwd = 0.6)
    p <- p + geom_hline(yintercept = -log10(pvalue), lty = 4, col = "black", lwd = 0.6)
    if(!is.null(xlims)) p <- p + xlim(xlims[1], xlims[2])
    if(!is.null(ylims)) p <- p + ylim(0, ylims)

    if(class(label) == "character"){
        label2 = intersect(rownames(DEG_tb), label)
        if(length(label2) == 0){
            cat2("Not found all specific labels!\n")
        } else {
            if(length(label2) != length(label)){
                cat2("Label(s) ", paste(setdiff(label, rownames(DEG_tb)), collapse = ","), " not found!\n")
            }
            DEG_tb$label = ifelse(rownames(DEG_tb) %in% label2, rownames(DEG_tb), "")
        }
    } else {
        if(label){
            if(method == "logFC"){
                DEG_tb = DEG_tb[order(DEG_tb$logFC, decreasing = TRUE), ]
                DEG_tb_UP = DEG_tb[DEG_tb$DEG == "Up", ]
                UP_top10 = head(rownames(DEG_tb_UP), label.top)
                DEG_tb_DN = DEG_tb[DEG_tb$DEG == "Down", ]
                DN_top10 = tail(rownames(DEG_tb_DN), label.top)
            } else if(method == "pvalue"){
                if(padj){
                    DEG_tb = DEG_tb[order(DEG_tb$adj.P.Val), ]
                    DEG_tb_UP = DEG_tb[DEG_tb$DEG == "Up", ]
                    UP_top10 = head(rownames(DEG_tb_UP), label.top)
                    DEG_tb_DN = DEG_tb[DEG_tb$DEG == "Down", ]
                    DN_top10 = head(rownames(DEG_tb_DN), label.top)            
                } else {                
                    DEG_tb = DEG_tb[order(DEG_tb$P.Value), ]
                    DEG_tb_UP = DEG_tb[DEG_tb$DEG == "Up", ]
                    UP_top10 = head(rownames(DEG_tb_UP), label.top)
                    DEG_tb_DN = DEG_tb[DEG_tb$DEG == "Down", ]
                    DN_top10 = head(rownames(DEG_tb_DN), label.top)            
                }

            }
            DEG_tb$label = ifelse(rownames(DEG_tb) %in% c(UP_top10, DN_top10), rownames(DEG_tb), "")

        }
    }

    if("label" %in% colnames(DEG_tb)){
        if(padj){
            p <- p + geom_label_repel(data = DEG_tb, aes(logFC, -log10(adj.P.Val), 
                label = label), size = label.cex, box.padding = unit(0.5, "lines"), fill = "white",
                point.padding = unit(0.5, "lines"), segment.color = "grey", 
                min.segment.length = 0,
                force= 2,
                # segment.linetype = 3,
                force_pull= 2,
                segment.color = 'black',
                segment.alpha = 0.5,
                # nudge_x= 3 - up$log2FoldChange, 
                show.legend = FALSE, max.overlaps = max.overlaps)
        } else {
            p <- p + geom_label_repel(data = DEG_tb, aes(logFC, -log10(P.Value), 
                label = label), size = label.cex, box.padding = unit(0.5, "lines"), fill = "white",
                point.padding = unit(0.5, "lines"), segment.color = "grey", 
                min.segment.length = 0,
                force= 2,
                # segment.linetype = 3,
                force_pull= 2,
                segment.color = 'black',
                segment.alpha = 0.5,
                # nudge_x= 3 - up$log2FoldChange, 
                show.legend = FALSE, max.overlaps = max.overlaps)
        }            
    }

    if(!is.null(out_name)) save.graph(p, file = out_name, w, h)
    return(p)
}

run_DEG_deseq2 <- function(expM, Group, Case, Control, out_name = NULL, log2FC = 1, pvalue = 0.05, padj = TRUE){

    loadp(DESeq2)
    if(ncol(expM) != length(Group)) stop("Length not matched between Matrix and Group!")
    if(! Case %in% Group) stop(paste0("Not found Group name: ", Case, "."))
    if(! Control %in% Group) stop(paste0("Not found Group name: ", Control, "."))

    # subset
    expM = expM[apply(expM, 1, sd) != 0, ]
    idx = Group %in% c(Case, Control)
    expM = expM[, idx]
    Group = Group[idx]

    design = factor(Group, levels = levels(factor(Group)))
    dds = DESeqDataSetFromMatrix(expM, DataFrame(design), design = ~ design)
    dds <- DESeq(dds, fitType = "local")
    DEG_tb = as.data.frame(results(dds))
    colnames(DEG_tb)[c(2, 5, 6)] = c("logFC", "P.Value", "adj.P.Val")

    DEG = rep("Not_Change", nrow(DEG_tb))
    if(padj){
        DEG[DEG_tb$adj.P.Val < pvalue & DEG_tb$logFC > log2FC]  = "Up"
        DEG[DEG_tb$adj.P.Val < pvalue & DEG_tb$logFC < -log2FC] = "Down"
    } else {
        DEG[DEG_tb$P.Value < pvalue & DEG_tb$logFC > log2FC]  = "Up"
        DEG[DEG_tb$P.Value < pvalue & DEG_tb$logFC < -log2FC] = "Down"
    }
    DEG_tb$DEG = DEG

    # report
    cat2(paste("Number of Up regulated genes:", length(DEG[DEG == "Up"])), "\n")
    cat2(paste("Number of Down regulated genes:", length(DEG[DEG == "Down"])), "\n")
    cat2(paste("Number of Not Change regulated genes:", length(DEG[DEG == "Not_Change"])), "\n")

    if(!is.null(out_name)){
        save(DEG_tb, file = paste0(out_name, ".RData"))
        write.csv(DEG_tb, file = paste0(out_name, ".csv"))
    }        
    return(DEG_tb)
}





bioForest <- function(coxFile=null,forestFile=null,forestCol=null){
  
  rt <- read.table(coxFile, header=T, sep="\t", check.names=F, row.names=1)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  height=nrow(rt)/12+5
  pdf(file=forestFile, width = 7,height = height)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
  
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="black",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol[1], forestCol[2])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.6)
  axis(1)
  dev.off()
}

format_size <- function(size){
    units = c("B", "KB", "MB", "GB", "TB", "PB", "EB")
    counter = 1
    while(size >= 1024 & counter < length(units)){
        size = size / 1024
        counter = counter + 1
    }
    size = paste0(format(trunc(size * 100, 2) / 100, big.mark = ","), " ", units[counter])
    return(size)
}

show_mem <- function(min_size = 2^20, igc = TRUE){
    mems = sort(sapply(mget(ls(envir = parent.frame()), envir = parent.frame()), pryr::object_size), decreasing = TRUE)
    mems_format = sapply(mems, format_size)
    mems_all = format_size(sum(mems))
    res = data.frame(object = names(mems_format), memory = mems_format, memory_Byte = mems, row.names = NULL)
    cat("Total memory:", mems_all, "\n")
    print(res[res$memory_Byte >= min_size, -3])
    if(igc) invisible(gc())
}

report <- function (info, type = "I", print = TRUE, file = NULL, append = TRUE) {
    type_mapping = c(I = "INFO", W = "WARN", E = "ERROR")
    if (!type %in% names(type_mapping)) {
        stop(paste0("Unsupported typs: ", type, ". Available types: I(Info)|W(Warning)|E(Error)"))
    }
    type = type_mapping[type]
    if (print) {
        cat2(paste0(type, " [", subString(as.character(Sys.time()), 1, "\\."), "] ", 
            info, "\n"))
    }
    if (!is.null(file)) {
        cat2(paste0(type, " [", subString(as.character(Sys.time()), 1, "\\."), "] ", 
            info, "\n"), file = file, append = append)
    }
}

opencwd <- function() shell.exec(getwd())

openpath <- function(path) shell.exec(path)

openfile <- function(file) shell.exec(file)

openurl <- function(url) shell.exec(url)

pheatmap_DEG <- function(expM, DEG_tb, group, ntop = 20, 
                         out_name = "Heatmap_DEG", w = 6, h = 6,
                         cluster_cols = FALSE, cluster_rows = FALSE, 
                         show_rownames = TRUE, show_colnames = FALSE){
    loadp(pheatmap)
    expM = expM[rownames(DEG_tb), ]
    pdata = expM[DEG_tb$DEG != "Not_Change", ]

    DEG_tb_UP = DEG_tb[DEG_tb$DEG == "Up", ]
    gene_UP_top = rownames(DEG_tb_UP)[order(DEG_tb_UP$adj.P.Val)[1:(ifelse(nrow(DEG_tb_UP) > ntop, ntop, nrow(DEG_tb_UP)))]]

    DEG_tb_DN = DEG_tb[DEG_tb$DEG == "Down", ]
    gene_DN_top = rownames(DEG_tb_DN)[order(DEG_tb_DN$adj.P.Val)[1:(ifelse(nrow(DEG_tb_DN) > ntop, ntop, nrow(DEG_tb_DN)))]]

    group1 = unique(group)[1]
    group2 = unique(group)[2]

    SID1 = colnames(expM)[group == group1]
    SID2 = colnames(expM)[group == group2]

    pdata = pdata[c(gene_UP_top, gene_DN_top), c(SID1, SID2)]
    annotation_col  = data.frame(Group = rep2(c(group1, group2), c(length(SID1), length(SID2))), row.names = c(SID1, SID2))
    annotation_row  = data.frame(DEG = rep2(c("Up", "Down"), c(length(gene_UP_top), length(gene_DN_top))), row.names = c(gene_UP_top, gene_DN_top))

    annotation_colors = list(Group = col.cluster2.1[2:3], DEG = c("Up" = col.cluster2.1[4], "Down" = col.cluster2.1[1]))
    names(annotation_colors[["Group"]]) = c(group1, group2)

    p <- pheatmap(pdata, 
             show_rownames = show_rownames, 
             show_colnames = show_colnames,
             scale = "row",
             cluster_cols = cluster_cols, 
             cluster_rows = cluster_rows,
             annotation_col  = annotation_col, 
             annotation_row  = annotation_row, 
             annotation_colors = annotation_colors,
             border_color = NA,
             gaps_row = c(length(gene_UP_top)), 
             gaps_col = c(length(SID1)),
             color = colorRampPalette(c("blue", "white", "red"))(50),
             annotation_names_col = FALSE,
             annotation_names_row = FALSE
    )

    p <- ggplotify::as.ggplot(p) + theme(plot.margin = margin(l = 5, r = 10))
    save.graph(p, file = out_name, w, h)
}




barplot2 <- function (data, out = NULL, PDF = TRUE, PNG = TRUE, w = 6, h = 6, 
    n = 30, title = "", padjust = TRUE, addLine = FALSE, fill = "#0089CB") 
{
    if (is.null(out)) 
        stop("Please specific parameter 'out' to save pdf file.\n")
    if (n > nrow(data)) 
        n = nrow(data)
    data = data.frame(data)[1:n, ]
    if (padjust) {
        data = data[order(data$p.adjust), ]
        data$Description = factor(data$Description, levels = rev(data$Description))
        p <- ggplot(data) + setTheme() + setText(20)
        p <- p + geom_bar(aes(x = Description, y = -log10(p.adjust)), 
            stat = "identity", width = 0.8, fill = fill, alpha = 1)
    }
    else {
        data = data[order(data$pvalue), ]
        data$Description = factor(data$Description, levels = rev(data$Description))
        p <- ggplot(data) + setTheme() + setText(20)
        p <- p + geom_bar(aes(x = Description, y = -log10(pvalue)), 
            stat = "identity", width = 0.8, fill = fill, alpha = 1)
    }
    p <- p + labs(x = "", title = title) + coord_flip() + theme(plot.title = element_text(hjust = 0))
    p <- p + theme(axis.text.y = element_text(hjust = 1)) + theme(axis.text.y = element_text(color = "black"))
    if (addLine) {
        p <- p + geom_line(aes(x = Description, y = Count, group = 1), 
            size = 0.5, color = "#BC0909")
        p <- p + geom_point(aes(x = Description, y = Count), 
            size = 1)
    }
    print(p)
    save.graph(p, file = out, w, h, PDF = PDF, PNG = PNG)
}

dotplot2 <- function (data, w = 6, h = 6, n = 20, colorHigh = "#DD1C77", x = "GeneRatioValue",
    colorLow = "#3182bd", sortBy = "p.adjust", title = "", padjust = TRUE, 
    outName = NULL, PDF = FALSE, PNG = TRUE, pathways = NULL, add_pathway = NULL, ...) 
{
    loadp(ggplot2)
    if (is.null(outName)) 
        stop("Please specific parameter 'outName' to save pdf file.\n")
    if (n > nrow(data)) 
        n = nrow(data)
    if(!is.null(pathways)){
        data = data.frame(data)[data$Description %in% pathways, ]
    } else {
        data = data.frame(data)[1:n, ]
        if(!is.null(add_pathway)) data_key = data[data$Description %in% add_pathway, ]
        if(!is.null(add_pathway)) data = rbind(data_key, data)
    }

    data = data[!duplicated(data$Description), ]
    data = data[order(data$Count, decreasing = TRUE), ]
    if (sortBy == "Count") {
        data$Description = factor(data$Description, levels = rev(data$Description))
    }
    else if (sortBy == "p.adjust") {
        if (padjust) {
            data = data[order(data$p.adjust), ]
            data$Description = factor(data$Description, levels = rev(data$Description))
        }
        else {
            data = data[order(data$pvalue), ]
            data$Description = factor(data$Description, levels = rev(data$Description))
        }
    }
    data$GeneRatioValue = as.numeric(subString(data$GeneRatio, 
        1, "/"))/as.numeric(subString(data$GeneRatio, 2, "/"))
    data = data[order(data$GeneRatioValue, decreasing = TRUE), ]

    if(x == "GeneRatioValue"){
        p <- ggplot(data, aes(x = GeneRatioValue, y = Description)) + theme_bw()
        p <- p + labs(x = "GeneRatio", y = "", title = title) + theme(plot.title = element_text(hjust = 0.5))
    } else if(x == "Count"){
        p <- ggplot(data, aes(x = Count, y = Description)) + theme_bw()
        p <- p + labs(x = "Count", y = "", title = title) + theme(plot.title = element_text(hjust = 0.5))
    }

    if (padjust) {
        p <- p + geom_point(aes(color = p.adjust, size = Count))
        p <- p + scale_colour_gradient(low = colorHigh, high = colorLow, 
            limit = c(min(data$p.adjust), max(data$p.adjust)), 
            guide = guide_colourbar(reverse = TRUE))
    }
    else {
        p <- p + geom_point(aes(color = pvalue, size = Count))
        p <- p + scale_colour_gradient(low = colorHigh, high = colorLow, 
            limit = c(min(data$pvalue), max(data$pvalue)), guide = guide_colourbar(reverse = TRUE))
    }
    p <- p + theme(axis.text.y = element_text(color = "black"))
    if (!is.null(outName)) {
        if (is.null(w)) {
            max_nchar = max(nchar(as.character(data$Description)))
            w = max_nchar * 0.13 + 2.2
        }
        save.graph(p, file = outName, w, h, PDF = PDF, PNG = PNG, ...)
    }
    invisible(p)
}

col.cluster2.1 <- c('#1C79B7','#F38329','#2DA248','#DC403E','#976BA6','#8D574C','#D07DB0','#BDBF3B', '#27BDD0','#B0C7E5','#F9BA79','#A0CC8A','#F09EA1','#C7B0D1','#D4BCB8','#F7CCDD',
                    '#DBDD8D','#A3D7E5','#B04E4F','#A38A58','#ED5351','#0C8945','#0F71A7','#A82764', '#F8DAE4','#7A4B1F','#5E6BAE','#8AC997','#DAC9AC','#0F5148','#A0B5DC','#9F858E',
                    '#5C181D','#7B8380','#E8DDDA','#264220','#5AB747','#5169AE','#4B3D58','#CD428A', '#62615A','#B82129','#66762E')

set_image <- function(w, h) options(repr.plot.width = w, repr.plot.height = h)

cat2 <- function(...){
    cat(...)
    flush.console()
}

plotVolcano <- function (diffM, DIFFcutoff = 2, FDRcutoff = 0.05, title = "Volcano Plot", 
    xtitle = "log2(FoldChange)", xlim = 5, ylim = 10, padjust = FALSE, 
    size = 0.5, pos = 8.5, label = TRUE) 
{
    library("ggplot2")
    if (padjust) {
        p <- ggplot(diffM, aes(x = log2FC, y = -log10(padj), 
            color = DEG)) + xlim(-xlim, xlim) + ylim(0, ylim)
        p <- p + theme_bw() + labs(title = title, x = xtitle, 
            y = "-log10(FDR)")
    }
    else {
        p <- ggplot(diffM, aes(x = log2FC, y = -log10(pvalue), 
            color = DEG)) + xlim(-xlim, xlim) + ylim(0, ylim)
        p <- p + theme_bw() + labs(title = title, x = xtitle, 
            y = "-log10(P-value)")
    }
    p <- p + geom_point(size = size) + setText(20) + theme(plot.title = element_text(hjust = 0.5))
    p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + scale_color_manual(values = c(UP = "red", DN = "dodgerblue4", 
        NC = "grey"))
    p <- p + theme(legend.position = "none")
    p <- p + geom_vline(xintercept = c(log2(DIFFcutoff), -log2(DIFFcutoff)), 
        lty = 4, col = "black", lwd = 0.6)
    p <- p + geom_hline(yintercept = -log10(FDRcutoff), lty = 4, 
        col = "black", lwd = 0.6)
    N_DN = sum(diffM$DEG == "DN")
    p <- p + annotate("text", x = -xlim * 2/5, y = pos, label = paste("DN:", 
        N_DN), size = 4, hjust = 1, color = "dodgerblue4", parse = TRUE)
    N_UP = sum(diffM$DEG == "UP")
    p <- p + annotate("text", x = xlim * 2/5, y = pos, label = paste("UP:", 
        N_UP), size = 4, hjust = 0, color = "red", parse = TRUE)
    if (label) {
        if(padjust){
            p <- p + geom_text_repel(data = diffM, aes(log2FC, -log10(padj), 
                label = lable), size = 3, box.padding = unit(0.5, 
                "lines"), point.padding = unit(0.8, "lines"), segment.color = "black", 
                show.legend = FALSE, max.overlaps = 50)
        } else {
            p <- p + geom_text_repel(data = diffM, aes(log2FC, -log10(pvalue), 
                label = lable), size = 3, box.padding = unit(0.5, 
                "lines"), point.padding = unit(0.8, "lines"), segment.color = "black", 
                show.legend = FALSE, max.overlaps = 50)
        }

    }
    return(p)
}

dotplotGO <- function (data, size = 28, padj = TRUE, outName = NULL, w = NULL, 
    h = 8, PDF = TRUE, PNG = TRUE, pathways = NULL) 
{
    if (padj) {
        data = data[data$p.adjust < 0.05, ]
        data = data[order(data$p.adjust), ]
    }
    else {
        data = data[data$pvalue < 0.05, ]
        data = data[order(data$pvalue), ]
    }
    if(!is.null(pathways)) data_key = data[data$Description %in% pathways, ]
    data_BP = data[data$ONTOLOGY == "BP", ]
    data_BP = data_BP[1:ifelse(nrow(data_BP) < 10, nrow(data_BP), 10), ]
    data_CC = data[data$ONTOLOGY == "CC", ]
    data_CC = data_CC[1:ifelse(nrow(data_CC) < 10, nrow(data_CC), 10), ]
    data_MF = data[data$ONTOLOGY == "MF", ]
    data_MF = data_MF[1:ifelse(nrow(data_MF) < 10, nrow(data_MF), 10), ]
    data = na.omit(rbind(data_BP, data_CC, data_MF))
    if(!is.null(pathways)) data = rbind(data_key, data)

    data$Description = factor(data$Description, levels = rev(data$Description))
    data$GeneRatioValue = as.numeric(subString(data$GeneRatio, 
        1, "/"))/as.numeric(subString(data$GeneRatio, 2, "/"))
    data = data[order(data$GeneRatioValue, decreasing = TRUE), 
        ]
    loadp(ggplot2, ggthemes)
    p <- ggplot(data, aes(x = GeneRatioValue, y = Description)) + 
        setText(25, theme = "bw")
    if (padj) {
        p <- p + geom_point(aes(color = p.adjust, size = Count))
        p <- p + scale_colour_gradient(low = "#DD1C77", high = "#3182bd", 
            limit = c(min(data$p.adjust), max(data$p.adjust)), 
            guide = guide_colourbar(reverse = TRUE))
    }
    else {
        p <- p + geom_point(aes(color = pvalue, size = Count))
        p <- p + scale_colour_gradient(low = "#DD1C77", high = "#3182bd", 
            limit = c(min(data$pvalue), max(data$pvalue)), guide = guide_colourbar(reverse = TRUE))
    }
    p <- p + labs(x = "Gene Ratio", y = "", title = "")
    p <- p + theme(axis.text.y = element_text(color = "black")) + 
        theme(axis.text.y = element_text(hjust = 1))
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    p <- p + facet_grid(ONTOLOGY ~ ., scale = "free")
    p <- p + theme(strip.text.y = element_text(size = size * 
        0.5, angle = 0))
    p <- p + theme(strip.background.y = element_rect(fill = "white"))
    if (!is.null(outName)) {
        if (is.null(w)) {
            max_nchar = max(nchar(as.character(data$Description)))
            w = max_nchar * 0.12 + 1.8
        }
        save.graph(p, file = outName, w, h, PDF = PDF, PNG = PNG)
    }
    invisible(p)
}

col2rownames <- function(counts, idx = 1, remove = TRUE, first = FALSE, sep = "_"){
    row_names = as.character(counts[, idx])
    if(sum(duplicated(row_names)) > 0){
        warning("Warnning: Duplicated rownames found!, add suffix.")
        row_names = string_redup(row_names, first = first, sep = sep)
    }
    rownames(counts) = row_names
    if(remove){
        counts[, idx] = NULL
    }
    return(counts)
}

string_redup <- function(input_strings, first = TRUE, sep = "_") {
    suffix = ave(input_strings, input_strings, FUN = function(x) if(length(x) > 1) seq(x) else 1)
    output_strings = paste0(input_strings, sep, suffix)
    if(!first) output_strings[match(unique(input_strings), input_strings)] = unique(input_strings)
    return(output_strings)
}


extract_TCGA_sample_type <- function(barcodes, return = "long"){
    sample_type_code = substr(barcodes, 14, 15)
    comp_list = list(
        "01" = "Primary Solid Tumor; TP",
        "02" = "Recurrent Solid Tumor  ; TR",
        "03" = "Primary Blood Derived Cancer - Peripheral Blood; TB",
        "04" = "Recurrent Blood Derived Cancer - Bone Marrow   ; TRBM",
        "05" = "Additional - New Primary   ; TAP",
        "06" = "Metastatic ; TM",
        "07" = "Additional Metastatic  ; TAM",
        "08" = "Human Tumor Original Cells ; THOC",
        "09" = "Primary Blood Derived Cancer - Bone Marrow ; TBM",
        "10" = "Blood Derived Normal   ; NB",
        "11" = "Solid Tissue Normal; NT",
        "12" = "Buccal Cell Normal ; NBC",
        "13" = "EBV Immortalized Normal; NEBV",
        "14" = "Bone Marrow Normal ; NBM",
        "15" = "sample type 15 ; 15SH",
        "16" = "sample type 16 ; 16SH",
        "20" = "Control Analyte; CELLC",
        "40" = "Recurrent Blood Derived Cancer - Peripheral Blood  ; TRB",
        "50" = "Cell Lines ; CELL",
        "60" = "Primary Xenograft Tissue   ; XP",
        "61" = "Cell Line Derived Xenograft Tissue ; XCL",
        "99" = "sample type 99 ; 99SH"
        )
    sample_type = unlist(comp_list[sample_type_code])
    sample_type_long = subString(sample_type, 1, "; ")
    sample_type_short = subString(sample_type, 2, "; ")

    if(return == "long"){
        return(sample_type_long)
    } else {
        return(sample_type_short)
    }

}

download_TCGA <- function(project, data_type, save.filename = NULL, directory = "./", dedup = "mean"){

    if(!suppressMessages(suppressWarnings(require(TCGAbiolinks)))){
        cat2("Please install package 'TCGAbiolinks' first:\n\tBiocManager::install(\"TCGAbiolinks\")\n")
        stop()
    }

    if(!suppressMessages(suppressWarnings(require(SummarizedExperiment)))){
        cat2("Please install package 'SummarizedExperiment' first:\n\tBiocManager::install(\"SummarizedExperiment\")\n")
        stop()
    }

    if(!suppressMessages(suppressWarnings(require(bioinfobaidu.utils)))){
        cat2("Please install package 'bioinfobaidu.utils' first:\n\tremotes::install_github(\"JiahaoWongg/bioinfobaidu.utils\")\n")
        stop()
    }

    if(data_type == "SNP"){
        if(!suppressMessages(suppressWarnings(require(maftools)))){
            cat2("Please install package 'maftools' first:\n\tBiocManager::install(\"maftools\")\n")
            stop()
        }
    }

    if(data_type == "450K"){
        if(!suppressMessages(suppressWarnings(require(sesameData)))){
            cat2("Please install package 'sesameData' first:\n\tBiocManager::install(\"sesameData\")\n")
            cat2("and install package 'sesame':\n\tBiocManager::install(\"sesame\")\n")
            stop()
        }
    }

    out_obj_name = paste0(gsub("-", "_", project), "_", data_type)

    cat2("Step-1: Query data ...\n")
    if(data_type == "mRNA"){
        query <- suppressMessages(
            GDCquery(
                project = project,
                data.category = "Transcriptome Profiling",
                data.type = "Gene Expression Quantification",
                workflow.type = "STAR - Counts"
            )
        )
    } else if (data_type == "miRNA"){
        query <- suppressMessages(
            GDCquery(
                project = project, 
                data.category = "Transcriptome Profiling",,
                data.type = "miRNA Expression Quantification"#,
                # access = "open"
            )
        )
    } else if (data_type == "CNV"){
        query <- suppressMessages(
            GDCquery(
                project = project,
                data.category = "Copy Number Variation",
                data.type = "Masked Copy Number Segment",              
                access = "open"
            )
        )
    } else if (data_type == "450K"){
        query <- suppressMessages(
            GDCquery(
                project = project, 
                data.category = "DNA Methylation", 
                data.type = "Methylation Beta Value", # Masked Intensities
                platform = "Illumina Human Methylation 450"
            )
        )
    } else if (data_type == "SNP"){
        query <- suppressMessages(
            GDCquery(
                project = project, 
                data.category = "Simple Nucleotide Variation",
                data.type = "Masked Somatic Mutation",
                access = "open"
            )
        )
    } else if (data_type == "protein"){
        query <- suppressMessages(
            GDCquery(
                project = project,
                data.category = "Proteome Profiling",
                data.type = "Protein Expression Quantification"
            )
        )
    }
    search_tb = query$results[[1]]
    cat2("\t- Found", nrow(search_tb), "samples\n")
    cat2("\t- Sample type(s):", paste0(names(table(search_tb$sample_type)), ": ", table(search_tb$sample_type), collapse = ", "), "\n")

    cat2("Step-2: Download data ...\n")
    suppressMessages(GDCdownload(query, files.per.chunk = 100, directory = directory))

    cat2("Step-3: Preapare data ...\n")
    if(is.null(save.filename)){
        save.filename = paste0(project, "_", data_type, ".RData")
    }
    
    data = invisible(suppressMessages(
        GDCprepare(query, save = TRUE, directory = directory, save.filename = save.filename)
    ))

    cat2("Step-4: Extract data ...\n")
    if(data_type == "mRNA"){
        gene_anno <- rowData(data) 
        sample_anno <- colData(data)
        sample_anno = data.frame(apply(sample_anno, 2, as.character))

        expM_count <- assay(data, "unstranded")
        expM_FPKM_raw <- log2(assay(data, "tpm_unstrand") + 1)
        expM_TPM_raw <- log2(assay(data, "fpkm_unstrand") + 1)

        if(dedup == "mean"){
            cat2("\t- Gene deduplicate using method:", dedup, "... ")
            expM_FPKM = col2rownames(aggregate(expM_FPKM_raw, by = list(factor(gene_anno$gene_name)), mean0), remove = TRUE)
            expM_TPM  = col2rownames(aggregate(expM_TPM_raw, by = list(factor(gene_anno$gene_name)), mean0), remove = TRUE)
        } else if (dedup == "unique"){
            cat2("\t- Gene deduplicate using method:", dedup, "... ")
            expM_FPKM = expM_FPKM_raw[!duplicated(gene_anno$gene_name), ]
            expM_TPM  = expM_TPM_raw[!duplicated(gene_anno$gene_name), ]
            rownames(expM_FPKM) = unique(gene_anno$gene_name)
            rownames(expM_TPM)  = unique(gene_anno$gene_name)
        }
        assign(out_obj_name, list2(data, gene_anno, sample_anno, expM_count, expM_FPKM_raw, expM_TPM_raw, expM_FPKM, expM_TPM))
    } else if (data_type == "miRNA"){
        expM_count = data[, seq(2, ncol(data), 3)]
        expM_RPM  = log2(data[, seq(3, ncol(data), 3)] + 1)
        dimnames(expM_count) = dimnames(expM_RPM) = list(data[, 1], subString(colnames(expM_count), 1, "_", rev = TRUE))
        assign(out_obj_name, list2(data, expM_count, expM_RPM))
    } else if (data_type == "CNV"){
        assign(out_obj_name, list2(data))
    } else if (data_type == "450K"){
        # https://blog.csdn.net/Ayue0616/article/details/127824296
        betaM = assay(data)
        probe_anno <- rowData(data) 
        sample_anno <- colData(data)
        sample_anno = data.frame(apply(sample_anno, 2, as.character))
        assign(out_obj_name, list2(data, betaM, probe_anno, sample_anno))
    } else if (data_type == "protein"){
        expM = col2rownames(data.frame(data, check.names = FALSE), idx = 5)
        expM = expM[, -c(1:5)]
        assign(out_obj_name, list2(data, expM))
    } else if (data_type == "SNP"){
        maf <- read.maf(data, verbose = FALSE)
        pdf(gsub(".RData", "_plotmafSummary.pdf", save.filename))
            plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
        dev.off()
        assign(out_obj_name, list2(data, maf))
    }

    res = get(out_obj_name)
    processed_out = paste0(dirname(save.filename), "/", gsub(".RData", "_processed.RData", basename(save.filename)))
    save(res, file = processed_out)

    cat2("\tSucceed!\n")
    invisible(res)
}

run_estimate <- function(expM, group = NULL, group_name = "Group", prefix = "./ESTIMATE_", Rplot = FALSE){
    if(!is.null(group)){
        if(length(group) != ncol(expM)){
            stop()
            cat(paste0("Length of 'group'(", length(group), ") do not equal to ncol of 'expM'(", ncol(expM), ")\n"))
        }
    }

    if(!suppressMessages(suppressWarnings(require(estimate)))){
        stop()
        cat0("Please install package 'estimate' first:\n\tinstall.packages(\"estimate\", repos=\"http://R-Forge.R-project.org\")")
    }

    if(!suppressMessages(suppressWarnings(require(GS)))){
        stop()
        cat0("Please install package 'GS' first:\n\tBiocManager::install(\"GS\")")
    }

    if(!suppressMessages(suppressWarnings(require(ggthemes)))){
        stop()
        cat0("Please install package 'ggthemes' first:\n\tinstall.packages(\"ggplot2\")")
    }

    cat2("- Step-1: Prepare input file ...\n")
    input_file = paste0(prefix, "_input.txt")
    write.table(expM, file = input_file, sep = "\t", quote = FALSE)

    cat2("- Step-2: Filter common genes ...\n")
    filter_file = paste0(prefix, "_filter.gct")
    filterCommonGenes(input_file, output.f = filter_file, id = "GeneSymbol")

    cat2("- Step-3: Estimate ...\n")
    estimate_file = paste0(prefix, "_score.txt")
    estimateScore(input.ds = filter_file, output.ds = estimate_file, platform = "affymetrix")
    est_res = read.table(estimate_file, skip = 2, sep = "\t", header = TRUE, check.names = FALSE)
    rownames(est_res) = est_res$NAME
    est_res[, c(1, 2)] = NULL
    colnames(est_res) = colnames(expM2)
    write("rc", est_res, estimate_file)

    rdata_file = paste0(prefix, "_estimate_score.RData")
    save(est_res, file = rdata_file)

    if(!is.null(group)){
        cat2("- Step-4: Plot ...\n")
        pdata = as.matrix(est_res)
        pdata = t(scale(t(pdata)))
        loadp(tinyarray, ggplot2)

        plot_file = paste0(prefix, "_Boxplot")
        p <- draw_boxplot(pdata, factor(group), color = c("#e5171a", "#1d4a9b"), sort = FALSE,
                     xlab = "", 
                     ylab = "Relative Enrichment Score", 
                     grouplab = group_name)
        p <- p + setText(20, x.text.angle = 35) + theme(legend.position = "top")
        save.graph(p, file = plot_file, 8, 5)

        if(Rplot){
            options(repr.plot.width = 8, repr.plot.height = 5)
            print(p)
        }
    }
    invisible(est_res)
}

run_immune_infiltration <- function(expM, group = NULL, group_name = "Group", gene_list, prefix = "./Immune_infiltration", Rplot = FALSE){
    if(!suppressMessages(suppressWarnings(require(GSVA)))){
        stop()
        cat0("Please install package 'GSVA' first:\n\tBiocManager::install(\"GSVA\")")
    }
    gsva_res <- suppressWarnings(gsva(ssgseaParam(as.matrix(expM), gene_list, normalize = TRUE)))
    save(gsva_res, file = paste0(prefix, "_28_cells.RData"))

    if(!is.null(group)){
        pdata = as.matrix(gsva_res)
        pdata = t(scale(t(pdata)))
        p <- draw_boxplot(pdata, factor(group), color = rev(c("#1d4a9b","#e5171a")), sort = FALSE,
                          xlab = "", 
                          ylab = "Relative Enrichment Score", 
                          grouplab = group_name)
        p <- p + setText(20, x.text.angle = 45) + theme(legend.position = "top")
        save.graph(p, file = paste0(prefix, "_28_cells"), 8*1.5, 4*1.5)

        if(Rplot){
            options(repr.plot.width = 8*1.5, repr.plot.height = 5*1.5)
            print(p)
        }
    }
    invisible(gsva_res)
}

run_check_point <- function(expM, group = NULL, group_name = "Group", gene_list = NULL, prefix = "./Check_point_", Rplot = FALSE){

    if(!suppressMessages(suppressWarnings(require(tinyarray)))){
        stop()
        cat0("Please install package 'tinyarray' first:\n\tinstall.packages(\"tinyarray\")")
    }

    if(is.null(gene_list)){
        gene_list = c("IDO1", "CD274", "HAVCR2", "PDCD1", "CTLA4", "LAG3", "PDCD1LG2")
    }

    gene_list_shared = intersect(gene_list, rownames(expM2))
    if(length(setdiff(gene_list, rownames(expM2))) > 0){
        cat2("Remove genes not found in expM:", paste(setdiff(gene_list, rownames(expM2)), collapse = ", "), "\n")
    }
    pdata = expM[gene_list, ]

    if(!is.null(group)){
        pdata = as.matrix(pdata)
        pdata = t(scale(t(pdata)))
        p <- draw_boxplot(pdata, factor(group), color = rev(c("#1d4a9b","#e5171a")), sort = FALSE,
                          xlab = "", 
                          ylab = "Relative Gene Expression", 
                          grouplab = group_name)
        p <- p + setText(20, x.text.angle = 45) + theme(legend.position = "top")
        save.graph(p, file = paste0(prefix, "_Boxplot"), 8, 5)

        if(Rplot){
            options(repr.plot.width = 8, repr.plot.height = 5)
            print(p)
        }
    }
    invisible(pdata)
}

run_core_pathway <- function(expM, group = NULL, group_name = "Group", gene_list, prefix = "./Core_pathway", Rplot = FALSE){

    if(!suppressMessages(suppressWarnings(require(GSVA)))){
        stop()
        cat0("Please install package 'GSVA' first:\n\tBiocManager::install(\"GSVA\")")
    }
    gsva_res <- suppressWarnings(gsva(ssgseaParam(as.matrix(expM), gene_list, normalize = TRUE)))
    save(gsva_res, file = paste0(prefix, "_core_biological_pathway.RData"))

    if(!is.null(group)){
        pdata = as.matrix(gsva_res)
        pdata = t(scale(t(pdata)))
        p <- draw_boxplot(pdata, factor(group), color = rev(c("#1d4a9b","#e5171a")), sort = FALSE,
                          xlab = "", 
                          ylab = "Relative Enrichment Score", 
                          grouplab = group_name)
        p <- p + setText(20, x.text.angle = 45) + theme(legend.position = "top")
        save.graph(p, file = paste0(prefix, "_core_biological_pathway"), 8*1.5, 4*1.5)

        if(Rplot){
            options(repr.plot.width = 8*1.5, repr.plot.height = 5*1.5)
            print(p)
        }
    }
    invisible(gsva_res)
}







