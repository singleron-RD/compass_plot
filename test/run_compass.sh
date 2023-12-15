source activate compass

compass --data th17_gene_expression_matrix.tsv \
    --num-processes 20 \
    --species homo_sapiens \
    --output-dir compass_outdir

