Here, we use compass to characterize cellular metabolic states based on single-cell RNA-Seq. Downstream visualization was conducted based on the author's analysis examples.

## Conda env
compass

## input data
The input gene expression matrix can be a tab-delimited text file (tsv) containing gene expression estimates (CPM, TPM, or similar scaled units) with one row per gene, one column per sample. Here, the author's data (rename as th17_geneuexpression_matrix.tsv) is used as the test data.

## Usage
1)Here is a simple example of running Compass. For more parameters, please refer to https://yoseflab.github.io/Compass/user-guide.html.
```
source activate compass

compass --data th17_gene_expression_matrix.tsv \
    --num-processes 20 \
    --species homo_sapiens \
    --output-dir compass_outdir
```

2)Visualize the results of Compass

step 1: run compass_plot.py
```
source activate compass

python ../script/compass_plot.py \
    --reaction_penalties compass_outdir/reactions.tsv \
    --metadata ../data/reaction_metadata.csv \
    --group_info group_info.tsv \
    --subsystem Glycolysis/gluconeogenesis \
    --outdir compass_plot_outdir \
    --name Glycolysis_gluconeogenesis


help info:
----
--reaction_penalties    reactions.tsv from the result of compass
--metadata              reaction metadata provided by the author
--group_info            a tsv contains group info. The first column corresponding to barcodes, The second column corresponding to group
--subsystem             subsystem used to plot, such as Glycolysis/gluconeogenesis. More subsystem can be found in the reaction metadata provided by the author
--plot_type             plot type, anno or no(t plot anno). Default: no
--outdir                output dir
--name                  name, just like subsystem, special symbols need to be removed
```

step 2: run compass_plot_update.R
```
source activate r4.1_env

Rscript ../script/compass_plot_update.R \
    --wilcox_res compass_plot_outdir/wilcox_res.tsv \
    --subsystem Glycolysis/gluconeogenesis \
    --xlab_name KD_vs_WT \
    --annotation_text_size 2 \
    --outdir compass_plot_outdir \
    --name Glycolysis_gluconeogenesis


help info:
----
--wilcox_res            wilcox result from running compass_plot.py
--subsystem             subsystem, need to be consistent with the subsystem parameter used for running compass_plot.py  
--xlab_name             xlab name of compass_plot.pdf from compass_plot.py
--annotation_text_size  annotation text size, default: 3
--outdir                output dir
--name                  name, need to be consistent with the name parameter used for running compass_plot.py
```

## results
1)wilcox_res.tsv:                                Comparison results between two groups using Wilcoxon rank sum test

2)Glycolysis_gluconeogenesis_update.pdf:         Scatter plot displays the comparison results between two groups of subsystems

3)Glycolysis_gluconeogenesis_anno_update.pdf:    Scatter plot displays the comparison results between two groups of subsystems, with metabolic process names attached
