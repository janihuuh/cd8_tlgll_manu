
## Activate virtual env

## Activate virtual env
me=$(whoami)
# source /Users/$me/Dropbox/lag3/applications/cellphone_db_venv/bin/activate
source /Users/janihuuh/Dropbox/AML_TIM3/applications/cpdb-venv/bin/activate
cd /Users/$me/Dropbox/lgll_sc/results/manuscript/cellphone/

## Use cellphonedb

cellphonedb method statistical_analysis --iterations=1000 --threads=28 \
    --counts-data hgnc_symbol \
    --project-name tlgll_v_non_leukemic \
    tlgll_v_non_leukemic_meta.txt tlgll_v_non_leukemic_counts.txt

cellphonedb method statistical_analysis --iterations=1000 --threads=28 \
    --counts-data hgnc_symbol \
    --project-name hyperexpanded_v_non_leukemic \
    hyperexpanded_v_non_leukemic_meta.txt hyperexpanded_v_non_leukemic_counts.txt
