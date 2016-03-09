######<------ Make supplemental figures/tables ------>######


###<- GTEx GO terms

# These GO terms were annotated using CPDB (http://cpdb.molgen.mpg.de/) -
# ConsensusPathDB-human, which integrates interaction networks in Homo sapiens
# including a variety of pathways. For our puroposes, we performed Gene 
# Ontology analysis using the option under over-representation analysis in 
# CPDB. Finally, this tool was favored over the other tools for the accessible
# output format, all categories are combined into one table including 
# genes overlapped in each category.

# The analysis was done using the top 100 driving genes in the clusters.

go_lists <- lapply(1:6, function(i) {
        go_tops <- 
                read.table(paste0("rdas/deng-clust", i, ".tab"),
                           sep = "\t",
                           header = TRUE)
        go_tops <- go_tops[go_tops$q.value < .01, ]
        #    as.character(go_top$term_name)
        with(go_tops,
             data.frame(term_goid, term_name,
                        term_category, 
                        annotated = effective_size,
                        significant = members_input_overlap_geneids) )
})
names(go_lists) <- paste0("clust_", c(1:6))

sapply(go_lists, NROW)

library(xtable)
xtable(go_lists[[1]])
xtable(go_lists[[2]])
xtable(go_lists[[3]])
xtable(go_lists[[4]])
xtable(go_lists[[5]])
xtable(go_lists[[6]])




###<- Deng et al GO terms

# The analysis was done using the top 100 driving genes in the clusters.

go_lists <- lapply(1:6, function(i) {
    go_tops <- 
        read.table(paste0("rdas/deng-clust", i, ".tab"),
                     sep = "\t",
                     header = TRUE)
    go_tops <- go_tops[go_tops$q.value < .01, ]
#    as.character(go_top$term_name)
    with(go_tops,
         data.frame(term_goid, term_name,
                    term_category, 
                    annotated = effective_size,
                    significant = members_input_overlap_geneids) )
})
names(go_lists) <- paste0("clust_", c(1:6))

sapply(go_lists, NROW)

library(xtable)
xtable(go_lists[[1]])
xtable(go_lists[[2]])
xtable(go_lists[[3]])
xtable(go_lists[[4]])
xtable(go_lists[[5]])
xtable(go_lists[[6]])
