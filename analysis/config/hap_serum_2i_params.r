    # UPDATE REQUIRED. data dir (were you unpacked the genomic database archive. Contains the misha db, redb, external files)
    sch_data_dir <<- "/net/mraid14/export/data/db/tgdb/schic2_mm9/"
    
    # misha db dir
    sch_groot <<- sprintf("%s/trackdb", sch_data_dir)

    # redb dir
    sch_redb_dir <<- sprintf("%s/redb", sch_data_dir)

    # external files dir
    sch_extfiles_dir <<- sprintf("%s/rawdata", sch_data_dir)

    #track base name
    sch_track_base <<- "scell.nextera.N[SX]T_[0-9]"

    # batch file name
    sch_batch_fn <<- "hap_serum_2i_es_batch.txt"

    # pool_tn
    pool_tn <<- "scell.nextera.pool_good_hap_2i_serum_es"

    # ins tn
    ins_scale <<- 3e5
    ins_tn <<- sprintf("%s_ins_%ds", pool_tn, ins_scale)

    # UPDATE REQUIRED. Directory (or several ones) holding the output of the sequence processing pipeline, or where you unpacked the supplied contact maps. 
    sch_base_dir <<- "/net/mraid14/export/data/users/lubling/datasets/scell"

    # UPDATE REQUIRED. analysis output directory
    sch_work_dir <<- "/net/mraid14/export/data/users/lubling/datasets/scell/results/esh/haploid"

    # table directory
    sch_table_dir <<- sprintf("%s/tables", sch_work_dir)

    # figures directory
    sch_fig_dir <<- sprintf("%s/figs", sch_work_dir)

    # cached data directory
    sch_rdata_dir <<- sprintf("%s/rdata", sch_work_dir)

    sch_karyo_log2_enr_thresh <<- 1

    # universe name
    sch_cells_group <<- "hap_serum_2i"

    # grouping parameters
    sch_phasing_criteria <<- "newer"
    sch_near_dists <<- 2**c(14.5, 21)
    sch_mitotic_dists <<- 2**c(21, 23.5)

    sch_near_s_start <<- 0.611
    sch_near_mid_s   <<- 0.77
    sch_mitotic_min_post_m <<- 0.3
    sch_near_max_post_m <<- 0.42
    sch_pre_m_slope <<- -1.8
    sch_pre_m_intercept <<- 1

    # a/b clustering ofn prefix
    sch_ab_cluster_fn_pref <<- sprintf("%s/paper/fig_s_domains_trans_cluster_ws", sch_fig_dir)

    # epigen track names
    sch_rna_tn_rep1 <<- "rna.129_hap.ES.pf_rep1"
    sch_rna_tn_rep2 <<- "rna.129_hap.ES.pf_rep2"

    sch_ctcf_tns <<- c("Encode.esb4.ctcf.rep1", "Encode.esb4.ctcf.rep2")
