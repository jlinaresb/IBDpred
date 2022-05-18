cohort = 'morgan' #gevers #morgan

if (cohort == 'morgan') {
  valRes = readRDS('~/projects/Metagenomics/validation/validation_IBD_morgan_uc_k40')
  vi = readRDS('~/projects/Metagenomics/plots/validation/VarImp_rf_k40_morgan.rds')
  title = 'Ulcerative Colitis (Morgan et al.)'
} else if (cohort == 'gevers'){
  valRes = readRDS('~/projects/Metagenomics/validation/validation_IBD_gevers_cd_k40')
  vi = readRDS('~/projects/Metagenomics/plots/validation/VarImp_rf_k40_gevers.rds')
  title = 'Crohn´s Disease (Gevers et al.)'
}

# Run validation.r and prediction.r scripts with specific arguments!!!!!


# HeatMap
# ==========================
require(ComplexHeatmap)

head(valRes)
table(valRes$target)

res = pred$data
res = res[order(res$truth),]
head(res)
pats = rownames(res)

valRes = valRes[order(valRes$target), ]


library(circlize)
col_fun = colorRamp2(c(0, 4, 8), c("white", "gray", "red"))
col_fun(seq(0, 9))


ha = HeatmapAnnotation(Disease_status = res$truth,
                       Prediction = res$response,
                       show_legend = T,
                       annotation_name_gp = gpar(fontsize = 15),
                       col = list(Disease_status = c("NO" = "darkgoldenrod2", "YES" = "tomato1"),
                                  Prediction = c("NO" = "darkgoldenrod2", "YES" = "tomato1")))

mat = subset(valRes, select = -c(target))
feats = vi$data[order(vi$data$Importance, decreasing = T),]$features
ComplexHeatmap::Heatmap(t(mat),
                        name = 'counts',
                        show_row_names = T,
                        show_column_names = F,
                        cluster_rows = F,
                        cluster_columns = F,
                        show_row_dend = F,
                        # column_order = pats,
                        row_order = feats,
                        col = col_fun,
                        row_names_gp = gpar(fontsize = 15),
                        top_annotation = ha)
