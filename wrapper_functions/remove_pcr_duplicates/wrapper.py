_author__ = "Kyle Gontjes"
__copyright__ = "Copyright 2024, Kyle Gontjes"
__email__ = "kgontjes@umich.edu"
__license__ = "Not Determined"

import os
from snakemake.shell import shell

# Run picard
shell("picard MarkDuplicates -REMOVE_DUPLICATES true -INPUT {snakemake.input.sorted_bam} -OUTPUT {snakemake.output.bam_duplicates_removed} -METRICS_FILE {snakemake.output.picard_metrics} -CREATE_INDEX false -VALIDATION_STRINGENCY LENIENT &> {snakemake.log.picard}")
 
# Samtools sort command
shell("samtools sort {snakemake.output.bam_duplicates_removed} -m 500M -@ 0 -o {snakemake.output.final_bam} -T {snakemake.params.samtools_temp} &> {snakemake.log.post_picard_sort}") 