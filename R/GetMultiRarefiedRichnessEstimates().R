#' Get replicate richness estiamtes by rarefying multiple times.
#'
#' This function produces multi-rarefied richness estiamtes from a hyloseq object.
#'
#' @param ps A phyloseq object as input.
#' @param RarefyDepth The deth to rarefy to. All samples shallower than this depth are discarded.
#' @param NumReplicates The number of rarefied replicates to produce
#' @return A tibble with mutliple richness metrics for all samples deeper than the rarefy depth. Each has all metadata from the phyloseq object along with the random seed used to generate it.
#' @examples
#' GetMultiRarefiedRichnessEstimates(ps=ps, RarefyDepth=2000)
#' @export

GetMultiRarefiedRichnessEstimates<-function(ps, RarefyDepth=2000, NumReplicates=100){ # default 10,000 min sample size

    if (all(RarefyDepth>(phyloseq::sample_sums(ps)))) {stop("RarefyDepth deeper than deepest sample")}
    ps_multiplerarefieds_list<-list()

    ps<-phyloseq::prune_samples(phyloseq::sample_sums(ps)>RarefyDepth, ps)

    for (i in 1:NumReplicates) {

            richness_rarefied<-ps %>%
                phyloseq::rarefy_even_depth(sample.size = RarefyDepth, rngseed = i,
                                    replace = TRUE, trimOTUs = TRUE, verbose = TRUE) %>%
                phyloseq::estimate_richness( measures=c("Observed", "Shannon", "Chao1"))

            meta<- ps%>%
                phyloseq::sample_data(.) %>%
                tibble::as_tibble(.) %>%
                dplyr::mutate(., "sample_names" = ps %>% phyloseq::sample_names(.) %>% make.names )

            combined_richness <- meta %>%
                dplyr::left_join(richness_rarefied%>% tibble::rownames_to_column(.),
                            by = c("sample_names"="rowname")) %>%
                dplyr::mutate_if( is.character,forcats::as_factor)
            
            ps_multiplerarefieds_list[[i]]<-combined_richness
        }

    ps_multiplerarefieds<-dplyr::bind_rows(ps_multiplerarefieds_list, .id = "RandomSeed")

    return(ps_multiplerarefieds)

}