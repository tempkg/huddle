import os

OUTDIR = 'out_pyscenic'
os.makedirs(OUTDIR, exist_ok=True)


rule all:
    input:
        f"{OUTDIR}/auc_mtx.csv"


rule grn:
    input:
        expr = config['exp_mtx'],
        tf_list = config['tf_list']
    output:
        adj = f"{OUTDIR}/expr_mat.adjacencies.tsv"
    params:
        seed = config['params']['grn']['seed'],
        workers = config['params']['grn']['num_workers']
    shell:
        """
        {config[tools][pyscenic]} grn \
            -m grnboost2 \
            --seed {params.seed} \
            --num_workers {params.workers} \
            -o {output.adj} \
            {input.expr} \
            {input.tf_list}
        """


rule ctx:
    input:
        adj = f"{OUTDIR}/expr_mat.adjacencies.tsv",
        dbs = config['databases'],
        anno = config['annotation_file'],
        expr = config['exp_mtx']
    output:
        regulons = f"{OUTDIR}/regulons.csv"
    params:
        workers = config['params']['ctx']['num_workers']
    shell:
        """
        {config[tools][pyscenic]} ctx \
            {input.adj} \
            {input.dbs} \
            --annotations_fname {input.anno} \
            --expression_mtx_fname {input.expr} \
            --mode "custom_multiprocessing" \
            --output {output.regulons} \
            --num_workers {params.workers}
        """


rule aucell:
    input:
        expr = config['exp_mtx'],
        regulons = f"{OUTDIR}/regulons.csv"
    output:
        auc = f"{OUTDIR}/auc_mtx.csv"
    params:
        workers = config['params']['aucell']['num_workers']
    shell:
        """
        {config[tools][pyscenic]} aucell \
            {input.expr} \
            {input.regulons} \
            -o {output.auc} \
            --num_workers {params.workers}
        """



