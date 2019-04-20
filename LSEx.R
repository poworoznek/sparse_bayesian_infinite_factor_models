
gex = read.table("~/downloads/GTEx_Analysis_V6_eQTLInputFiles_geneLevelNormalizedExpressionMatrices/Pancreas_Analysis.expr.txt",
                 header = TRUE, stringsAsFactors = FALSE)
rownames(gex) = gex$Id
gex = t(as.matrix(gex[,-1]))
sgex = scale(gex)

library(GGally)
ggpairs(as.data.frame(sgex)[, sample.int(ncol(sgex),1) + 1:10])
anyNA(sgex)

source("dlfact.R")
source("mcrotfact.R")
source("clustalignplus.R")
source("permsignfact.R")
source("permuter.R")
source("signer.R")

saveRDS(gex, "gex.Rds")

subgex = gex[, 1:500]

k = 5

try = dlfact(gex, k, nrun = 1500, burn = 500)

lambda_sample = try$Lambda

sample_mean = reduce(lambda_sample, `+`)/length(lambda_sample)
rotated = mcrotfact(lambda_sample, method = "varimax", file = FALSE)
aligned = clustalignplus(rotated$samples, itermax = 1000)

label = "Sample Mean"
SampleMean = cbind(melt(sample_mean), label)
label = "Rotated Sample Mean"
RotatedMean = cbind(melt(rotated$mean), label)
label = "Aligned Sample Mean"
ProcessMean = cbind(melt(Reduce("+", aligned)/length(aligned)), label)

ggdf = rbind(SampleMean, RotatedMean, ProcessMean)

ggplot(ggdf, aes(x = Var2, y = Var1)) + 
  facet_grid(cols = vars(label)) +
  geom_tile(aes(fill=value), colour="grey20") + 
  scale_fill_gradient2(low = "#3d52bf", high = "#33961b", mid = "white") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.title = element_text(),
        plot.title = element_text(hjust = 0.5)) + 
  labs(fill = " ")

perm = sample.int(k)

switched = lapply(1:length(lambda_sample), function(ind){
  if(ind < length(lambda_sample)/2) return(lambda_sample[[ind]])
  return(lambda_sample[[ind]][, perm])
})
