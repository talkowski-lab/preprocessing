# preprocessing
Repository for QC and misc. workflows that are usually run following annotation. Also serves as a centralized location for "helper" WDLs and tasks that aren't workflows themselves but are still useful across workflows.

Each ```v0.x``` is associated with a published Dockstore workflow. 

**```ancestry-inference-hail-v01.wdl```**: https://dockstore.org/workflows/github.com/talkowski-lab/preprocessing/ancestry-inference-hail-v01
- Infer ancestry for a single cohort.

**```ancestry-inference-hail-mt-v01.wdl```**: https://dockstore.org/workflows/github.com/talkowski-lab/preprocessing/ancestry-inference-hail-mt-v01
- Infer ancestry for a single cohort, given an input that is in Hail's MatrixTable format.

**```ancestry-inference-hail-cohort-set-v01.wdl```**: https://dockstore.org/workflows/github.com/talkowski-lab/preprocessing/ancestry-inference-hail-cohort-set-v01
- Infer ancestry across multiple cohorts.

**```relatedness-hail-v01.wdl```**: https://dockstore.org/workflows/github.com/talkowski-lab/preprocessing/relatedness-hail-v01
- Infer relatedness and impute sex for a single cohort.

**```relatedness-hail-subset-samples-v01.wdl```**: https://dockstore.org/workflows/github.com/talkowski-lab/preprocessing/relatedness-hail-subset-samples-v01
- Infer relatedness and impute sex for a single cohort.
- Usually for larger cohorts (>10k samples—TODO: get a better estimate for this number??).

**```relatedness-hail-cohort-set-v01.wdl```**: https://dockstore.org/workflows/github.com/talkowski-lab/preprocessing/relatedness-hail-cohort-set-v01
- Infer relatedness and impute sex across multiple cohorts.
