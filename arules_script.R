#!/usr/bin/env Rscript
#load dataset and subscript by not fitted indexes
#install.packages("arules")
library(jsonlite)
library(stringr)
library(arules)

args <- commandArgs(trailingOnly = TRUE)
df_name <- args[1]

# df_name <- "DLBC"
# df_name <- "LAML"
df <- read.csv(sprintf("%s/no_fit.csv", df_name), header = FALSE)
no_fit_idxs <- as.vector(t(df))

# str(no_fit_idxs)
load(paste(df_name,'.mlSOM.input.Rdata',sep=""))
no_fit_df <- list(expr=data.list$indata.list$expr[no_fit_idxs,],
                  prom=data.list$indata.list$prom[no_fit_idxs,],
                  cnv=data.list$indata.list$cnv[no_fit_idxs,],
                  snv=data.list$indata.list$snv[no_fit_idxs,])
# str(no_fit_df)
# View(no_fit_df$expr)

# sample_count <- ncol(no_fit_df$expr)
gene_count <- nrow(no_fit_df$expr)

#~~~~~~~~~main~~~~~~~~~~~~~~

#arules for all genes
k <- 0
# k_all <- 0
all_gene_groups <- list()
genes_no_fit <- c()
for (gene_idx in 1:gene_count) {#gene_count
  df <- data.frame(expr=no_fit_df$expr[gene_idx,],
                   prom=no_fit_df$prom[gene_idx,],
                   cnv=no_fit_df$cnv[gene_idx,],
                   snv=no_fit_df$snv[gene_idx,])
  
  disc_df <- discretizeDF(df, methods = list(
    expr = list(method = "fixed", breaks = c(-Inf, mean(df$expr), Inf), labels = c("low", "high")), #, labels = c("low", "high")
    prom = list(method = "fixed", breaks = c(-Inf, mean(df$prom), Inf), labels = c("low", "high")),
    cnv = list(method = "fixed", breaks = c(-Inf, mean(df$cnv), Inf), labels = c("low", "high")),
    snv = list(method = "fixed", breaks = c(-Inf, 0, Inf), labels = c("low", "high"))
  ))
  
  # disc_df <- discretizeDF(df, default = list(method = "frequency", breaks = 2))
  rules <- apriori(disc_df, control = list(verbose = FALSE)) #, parameter = list(support = 0.1, confidence = 0.7)
  rules_sub <- subset(rules, subset = rhs %pin% "expr" & lift > 1) #
  nonr_rules <- rules_sub[!is.redundant(rules_sub)]
  sign_rules <- nonr_rules[is.significant(nonr_rules, disc_df, 
                                         method = "fisher",
                                         alpha = 0.05,
                                         adjust = "bonferroni")]
  # k_all = k_all + length(sign_rules)
  if (length(sign_rules) > 0){
    # print(gene_idx)
    k = k + 1
    res_rule <- sort(sign_rules, by=c("lift"))[1]
    # inspect(res_rule)
    lhs <- as(lhs(res_rule), "list")[[1]]
    rhs <- as(rhs(res_rule), "list")[[1]]
    expr <- str_extract(rhs, '\\b\\w+$')
    # print(expr)
    prom <- 0
    cnv <- 0
    snv <- 0
    for (item in lhs) {
      param_key <- str_extract(item, '\\b\\w+')
      param_value <- str_extract(item, '\\b\\w+$')
      # print(param_key)
      # print(param_value)
      switch(param_key,
             "prom" = {
               prom <- param_value
             },
             "cnv" = {
               cnv <- param_value
             },
             "snv" = {
               snv <- param_value
             }
      )
    }
    p <- ifelse(prom==expr, "+", ifelse(prom!=0, "-", 0))
    c <- ifelse(cnv==expr, "+", ifelse(cnv!=0, "-", 0))
    s <- ifelse(snv==expr, "+", ifelse(snv!=0, "-", 0))
    key <- paste(p, c, s, sep ="")
    
    value <- rownames(no_fit_df$expr)[gene_idx]
    # print(value)
    
    if (key %in% names(all_gene_groups)) {
      all_gene_groups[[key]] <- c(all_gene_groups[[key]], value)
    } else {
      all_gene_groups[[key]] <- value
    }

    # if (length(sign_rules) > 1){
    #   print(gene_idx)
    #   inspect(sign_rules)
    # }
  } 
  else { #not fitted genes
    genes_no_fit <- c(genes_no_fit, rownames(no_fit_df$expr)[gene_idx])
  }
  
}
print(all_gene_groups)
k
# k_all
write_json(all_gene_groups, sprintf("%s/gene_groups_arules.txt", df_name))
write.table(genes_no_fit, sprintf("%s/no_fit_arules.csv", df_name), sep = ",", row.names = FALSE, col.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
