    # data dir (misha db, redb, external files)
    sch_data_dir <<- "/net/mraid14/export/data/db/tgdb/schic2_mm9/"
    
    # misha db dir
    sch_groot <<- sprintf("%s/trackdb", sch_data_dir)

    # redb dir
    sch_redb_dir <<- sprintf("%s/redb", sch_data_dir)

    # external files dir
    sch_extfiles_dir <<- sprintf("%s/rawdata", sch_data_dir)

    #track base name
    sch_track_base <<- "scell.nextera.1CDU|scell.nextera.1CDX|scell.nextera.1CDES"

    # batch file name
    sch_batch_fn <<- "hyb_300317_2i_is_es_batch.txt"

    # pool_tn
    pool_tn <<- "scell.nextera.pool_good_hyb_2i_all_es"
 
    # insulation parameters
    ins_scale <<- 3e5
    ins_tn <<- sprintf("%s_ins_%ds", pool_tn, ins_scale)
    
    # base dir
    sch_base_dir <<- c("/net/mraid14/export/tgdata/db/tgdb/mm9/rawdata/scell_hic/cells_hyb_apr_2016/processed",
        "/net/mraid14/export/tgdata/db/tgdb/mm9/rawdata/scell_hic/cells_hyb_feb_2016/processed", 
	"/net/mraid14/export/tgdata/db/tgdb/mm9/rawdata/scell_hic/cells_hyb_idx_sort/processed")

    #table directory
    sch_table_dir <<- "/net/mraid14/export/data/users/lubling/datasets/scell/results/esh/hyb_mm9_2i_is/tables"
    sch_fig_dir <<- "/net/mraid14/export/data/users/lubling/datasets/scell/results/esh/hyb_mm9_2i_is/figs"
    sch_rdata_dir <<- "/net/mraid14/export/data/users/lubling/datasets/scell/results/esh/hyb_mm9_2i_is/rdata"

    #max_chromosome coverage bias (for filtering)
    #sch_karyo_z_thresh <<- 100
    sch_karyo_log2_enr_thresh <<- 1

    #max self level
    sch_max_non_digest <<- 0.55

    #downsamp param
    sch_min_cov <<- 20000

    #max contacts
    sch_max_cov <<- 700000

    #max trans contacts
    sch_max_trans <<- 0.15

    # mark as diploid
    sch_haploid <<- F

    #universe name
    sch_cells_group <<- "hyb_2i_idx_sort_es"

    # grouping parameters
    sch_phasing_criteria <<- "newer"

    sch_near_dists <<- 2**c(14.5, 21)
    sch_mitotic_dists <<- 2**c(21, 23.5)

    sch_near_s_start <<- 0.63
    sch_near_mid_s   <<- 0.785
    sch_mitotic_min_post_m <<- 0.3
    sch_near_max_post_m <<- 0.5
    sch_pre_m_slope <<- -1.8
    sch_pre_m_intercept <<- 1

    # a/b clustering ofn prefix
    # sch_ab_cluster_fn_pref <<- sprintf("%s/pool_domains_cluster_trans_area_kmeans_k2_d-10.000000_breaks.around.median_q.disp0.010000_2cols_min-5", sch_table_dir)
    sch_ab_cluster_fn_pref <<- sprintf("%s/paper/fig_s_domains_trans_cluster_area", sch_fig_dir)

    # epigen track names
    sch_rna_tn_rep1 <<- "rna.129_Cast.ES.pf_rep1"
    sch_rna_tn_rep2 <<- "rna.129_Cast.ES.pf_rep2"
    sch_ctcf_tns <<- c("Encode.esb4.ctcf.rep1", "Encode.esb4.ctcf.rep2")

    # dixon boundaries file name
    dixon_borders_fn <<- sprintf("%s/dixon2012/dixon2012_mESC_boundaries.txt", sch_extfiles_dir)
