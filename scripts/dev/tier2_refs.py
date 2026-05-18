"""Tier-2 canonical citations (downstream biology modules).

Verbatim from docs/SCIENTIFIC_AUDIT_2026-05-15.md Tier-2 section.
Consumed by scripts/dev/insert_module_references.py.
"""
from __future__ import annotations


TIER_2_REFS: dict[str, dict[str, dict[str, str]]] = {
    "rna_velocity.py": {
        "Bergen_scVelo_2020": {
            "title": "Generalizing RNA velocity to transient cell states through dynamical modeling",
            "authors": "Bergen et al.",
            "journal": "Nature Biotechnology",
            "year": "2020",
            "doi": "10.1038/s41587-020-0591-3",
            "description": "scVelo stochastic + dynamical models. Library used here.",
        },
        "LaManno_velocity_2018": {
            "title": "RNA velocity of single cells",
            "authors": "La Manno et al.",
            "journal": "Nature",
            "year": "2018",
            "doi": "10.1038/s41586-018-0414-6",
            "description": "Original spliced/unspliced RNA velocity formulation.",
        },
    },
    "trajectory.py": {
        "Wolf_PAGA_2019": {
            "title": "PAGA: graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells",
            "authors": "Wolf et al.",
            "journal": "Genome Biology",
            "year": "2019",
            "doi": "10.1186/s13059-019-1663-x",
            "description": "PAGA topology-preserving abstraction. sc.tl.paga used here.",
        },
        "Haghverdi_DPT_2016": {
            "title": "Diffusion pseudotime robustly reconstructs lineage branching",
            "authors": "Haghverdi, Buttner, Wolf et al.",
            "journal": "Nature Methods",
            "year": "2016",
            "doi": "10.1038/nmeth.3971",
            "description": "Diffusion pseudotime (sc.tl.dpt) used for ordering cells along trajectories.",
        },
    },
    "cell_cycle.py": {
        "Tirosh_cellcycle_2016": {
            "title": "Single-cell RNA-seq supports a developmental hierarchy in human oligodendroglioma",
            "authors": "Tirosh et al.",
            "journal": "Nature",
            "year": "2016",
            "doi": "10.1038/nature20123",
            "description": "Source of canonical S and G2/M cell-cycle gene sets; sc.tl.score_genes_cell_cycle implements this.",
        },
    },
    "cell_communication.py": {
        "Dimitrov_LIANA_2022": {
            "title": "Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data",
            "authors": "Dimitrov et al.",
            "journal": "Nature Communications",
            "year": "2022",
            "doi": "10.1038/s41467-022-30755-0",
            "description": "LIANA benchmark + framework. Used here via liana-py.",
        },
        "Efremova_CellPhoneDB_2020": {
            "title": "CellPhoneDB: inferring cell-cell communication from combined expression of multi-subunit ligand-receptor complexes",
            "authors": "Efremova et al.",
            "journal": "Nature Protocols",
            "year": "2020",
            "doi": "10.1038/s41596-020-0292-x",
            "description": "Source ligand-receptor resource consumed by LIANA.",
        },
    },
    "gene_regulatory_network.py": {
        "BadiaiMompel_decoupler_2022": {
            "title": "decoupleR: ensemble of computational methods to infer biological activities from omics data",
            "authors": "Badia-i-Mompel et al.",
            "journal": "Bioinformatics Advances",
            "year": "2022",
            "doi": "10.1093/bioadv/vbac016",
            "description": "decoupler ensemble framework used for TF activity scoring.",
        },
        "GarciaAlonso_DoRothEA_2019": {
            "title": "Benchmark and integration of resources for the estimation of human transcription factor activities",
            "authors": "Garcia-Alonso et al.",
            "journal": "Genome Research",
            "year": "2019",
            "doi": "10.1101/gr.240663.118",
            "description": "DoRothEA TF-target regulons consumed by decoupler.",
        },
    },
    "pathway_analysis.py": {
        "Subramanian_GSEA_2005": {
            "title": "Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles",
            "authors": "Subramanian et al.",
            "journal": "PNAS",
            "year": "2005",
            "doi": "10.1073/pnas.0506580102",
            "description": "GSEA methodology. gseapy is the Python port used here.",
        },
        "Schubert_PROGENy_2018": {
            "title": "Perturbation-response genes reveal signaling footprints in cancer gene expression",
            "authors": "Schubert et al.",
            "journal": "Nature Communications",
            "year": "2018",
            "doi": "10.1038/s41467-017-02391-6",
            "description": "PROGENy signaling-pathway responsive gene resource consumed via decoupler.",
        },
    },
    "cnv_inference.py": {
        "Patel_inferCNV_2014": {
            "title": "Single-cell RNA-seq highlights intratumoral heterogeneity in primary glioblastoma",
            "authors": "Patel et al.",
            "journal": "Science",
            "year": "2014",
            "doi": "10.1126/science.1254257",
            "description": "Chromosome-window CNV inference from scRNA expression — methodology this module re-implements via scipy.ndimage.uniform_filter1d.",
        },
        "inferCNV_Broad": {
            "title": "inferCNV — open-source CNV inference from scRNA-seq",
            "authors": "Tickle, T.L. et al. (Broad Institute)",
            "journal": "Open-source software",
            "year": "2019",
            "doi": "https://github.com/broadinstitute/inferCNV",
            "description": "R-side inferCNV is the canonical implementation; Python port here matches the windowed-mean methodology.",
        },
    },
    "composition.py": {
        "Buttner_scCODA_2021": {
            "title": "scCODA is a Bayesian model for compositional single-cell data analysis",
            "authors": "Buttner, Ostner, et al.",
            "journal": "Nature Communications",
            "year": "2021",
            "doi": "10.1038/s41467-021-27150-6",
            "description": "Bayesian compositional analysis; preferred path via pertpy. Chi-square fallback used when scCODA unavailable.",
        },
    },
    "gene_signature_scoring.py": {
        "Tirosh_scoring_2016": {
            "title": "Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq",
            "authors": "Tirosh et al.",
            "journal": "Science",
            "year": "2016",
            "doi": "10.1126/science.aad0501",
            "description": "Signature score methodology (mean expression minus matched control bin). sc.tl.score_genes implements this.",
        },
    },
    "immune_phenotyping.py": {
        "Tirosh_scoring_2016": {
            "title": "Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq",
            "authors": "Tirosh et al.",
            "journal": "Science",
            "year": "2016",
            "doi": "10.1126/science.aad0501",
            "description": "Same scoring methodology scoped to immune cell-type marker sets.",
        },
        "Zilionis_immune_atlas_2019": {
            "title": "Single-cell transcriptomics of human and mouse lung cancers reveals conserved myeloid populations",
            "authors": "Zilionis et al.",
            "journal": "Immunity",
            "year": "2019",
            "doi": "10.1016/j.immuni.2019.03.009",
            "description": "Canonical lung-cancer immune sub-typing markers underlying many phenotype calls.",
        },
    },
    "metacell.py": {
        "Persad_SEACells_2023": {
            "title": "SEACells infers transcriptional and epigenomic cellular states from single-cell genomics data",
            "authors": "Persad et al.",
            "journal": "Nature Biotechnology",
            "year": "2023",
            "doi": "10.1038/s41587-023-01716-9",
            "description": "SEACells metacell aggregation; primary backend.",
        },
    },
    "tumor_microenvironment.py": {
        "Aran_TME_2019": {
            "title": "Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage",
            "authors": "Aran et al.",
            "journal": "Nature Immunology",
            "year": "2019",
            "doi": "10.1038/s41590-018-0276-y",
            "description": "Reference-based lung TME phenotype assignments.",
        },
        "Tirosh_scoring_2016": {
            "title": "Dissecting the multicellular ecosystem of metastatic melanoma by single-cell RNA-seq",
            "authors": "Tirosh et al.",
            "journal": "Science",
            "year": "2016",
            "doi": "10.1126/science.aad0501",
            "description": "Score-gene-set framework used for TME marker scoring.",
        },
    },
    "pseudobulk_de.py": {
        "Squair_pseudobulk_2021": {
            "title": "Confronting false discoveries in single-cell differential expression",
            "authors": "Squair et al.",
            "journal": "Nature Communications",
            "year": "2021",
            "doi": "10.1038/s41467-021-25960-2",
            "description": "Benchmark demonstrating pseudobulk DE outperforms per-cell DE — methodology applied here.",
        },
    },
    "validate_cbioportal.py": {
        "Cerami_cBioPortal_2012": {
            "title": "The cBio Cancer Genomics Portal: An Open Platform for Exploring Multidimensional Cancer Genomics Data",
            "authors": "Cerami et al.",
            "journal": "Cancer Discovery",
            "year": "2012",
            "doi": "10.1158/2159-8290.CD-12-0095",
            "description": "External cancer genomics validation portal queried by this module.",
        },
    },
    "pseudo_velocity.py": {
        "project_local_proxy": {
            "title": "Project-local KNN-based pseudo-velocity proxy",
            "authors": "singlecell_factory contributors",
            "journal": "Internal documentation",
            "year": "2025",
            "doi": "PMID: 33288903",
            "description": "Lightweight transitional-probability proxy via PCA-space KNN flow when spliced/unspliced counts are unavailable. NOT a canonical RNA velocity method — see rna_velocity.py for that. PMID points to a transcription-velocity review for context only.",
        },
    },
    "evolution.py": {
        "Andor_clonal_evolution_2016": {
            "title": "Pan-cancer analysis of the extent and consequences of intratumor heterogeneity",
            "authors": "Andor et al.",
            "journal": "Nature Medicine",
            "year": "2016",
            "doi": "10.1038/nm.3984",
            "description": "Clonal evolution framework for tumor heterogeneity quantified here from CNV-derived signals.",
        },
    },
    "cell_fate.py": {
        "Lange_CellRank_2022": {
            "title": "CellRank for directed single-cell fate mapping",
            "authors": "Lange et al.",
            "journal": "Nature Methods",
            "year": "2022",
            "doi": "10.1038/s41592-021-01346-6",
            "description": "Directed transition-probability based fate mapping — methodology this module's Pearson sampling approximates.",
        },
    },
}
