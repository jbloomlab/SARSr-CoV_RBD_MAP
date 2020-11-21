"""``snakemake`` file that runs entire analysis."""

# Imports ---------------------------------------------------------------------
import glob
import itertools
import os.path
import textwrap
import urllib.request

import Bio.SeqIO

import dms_variants.codonvarianttable
import dms_variants.illuminabarcodeparser

import pandas as pd

# Configuration  --------------------------------------------------------------
configfile: 'config.yaml'

# run "quick" rules locally:
localrules: make_rulegraph,
            make_summary

# Functions -------------------------------------------------------------------
def nb_markdown(nb):
    """Return path to Markdown results of notebook `nb`."""
    return os.path.join(config['summary_dir'],
                        os.path.basename(os.path.splitext(nb)[0]) + '.md')

# Information on samples and barcode runs -------------------------------------
barcode_runs = pd.read_csv(config['barcode_runs'])

# combination of the *library* and *sample* columns should be unique.
assert len(barcode_runs.groupby(['library', 'sample'])) == len(barcode_runs)

# *sample* should be the hyphen separated concatenation of
# *experiment*, *antibody*, *concentration*, and *selection*.
sample_vs_expect = (
    barcode_runs
    .assign(expect=lambda x: x[['experiment', 'antibody', 'concentration',
                                'selection']]
                             .apply(lambda r: '-'.join(r.values.astype(str)),
                                    axis=1),
            equal=lambda x: x['sample'] == x['expect'],
            )
    )
assert sample_vs_expect['equal'].all(), sample_vs_expect.query('equal != True')

# barcode runs with R1 files expanded by glob
barcode_runs_expandR1 = (
    barcode_runs
    .assign(R1=lambda x: x['R1'].str.split('; ').map(
                    lambda y: list(itertools.chain(*map(glob.glob, y)))),
            n_R1=lambda x: x['R1'].map(len),
            sample_lib=lambda x: x['sample'] + '_' + x['library'],
            )
    )
assert barcode_runs_expandR1['sample_lib'].nunique() == len(barcode_runs_expandR1)
if any(barcode_runs_expandR1['n_R1'] < 1):
    raise ValueError(f"no R1 for {barcode_runs_expandR1.query('n_R1 < 1')}")

# Rules -----------------------------------------------------------------------

# this is the target rule (in place of `all`) since it first rule listed
rule make_summary:
    """Create Markdown summary of analysis."""
    input:
        rulegraph=os.path.join(config['summary_dir'], 'rulegraph.svg'),
        variant_counts=config['variant_counts'],
        count_variants=nb_markdown('count_variants.ipynb'),
        compute_barcode_escape='results/summary/compute_barcode_escape.md',
        homolog_escape='results/summary/homolog_escape.md',
    output:
        summary = os.path.join(config['summary_dir'], 'summary.md')
    run:
        def path(f):
            """Get path relative to `summary_dir`."""
            return os.path.relpath(f, config['summary_dir'])
        with open(output.summary, 'w') as f:
            f.write(textwrap.dedent(f"""
            # Summary

            Analysis run by [Snakefile]({path(workflow.snakefile)})
            using [this config file]({path(workflow.configfiles[0])}).
            See the [README in the top directory]({path('README.md')})
            for details.

            Here is the rule graph of the computational workflow:
            ![{path(input.rulegraph)}]({path(input.rulegraph)})

            Here is the Markdown output of each notebook in the workflow:

            1. [Count variants by barcode]({path(input.count_variants)}).
               Creates a [variant counts file]({path(input.variant_counts)})
               giving counts of each barcoded variant in each sample.
               
            2. [Compute escape fractions]({path(input.compute_barcode_escape)}) for individual barcodes.
            
            3. [Determine homolog escape fraction]({path(input.homolog_escape)}) averaged across all barcodes. Generates summary plots including heatmaps.
               
            """
            ).strip())

rule make_rulegraph:
    # error message, but works: https://github.com/sequana/sequana/issues/115
    input:
        workflow.snakefile
    output:
        os.path.join(config['summary_dir'], 'rulegraph.svg')
    shell:
        "snakemake --forceall --rulegraph | dot -Tsvg > {output}"

rule homolog_escape:
    input:
        config['escape_fracs_barcodes']
    output:
        config['escape_fracs_homologs'],
        md='results/summary/homolog_escape.md',
        md_files=directory('results/summary/homolog_escape_files')
    envmodules:
        'R/3.6.2-foss-2019b'
    params:
        nb='homolog_escape.Rmd',
        md='homolog_escape.md',
        md_files='homolog_escape_files'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        mv {params.md_files} {output.md_files}
        """

rule compute_barcode_escape:
    input:
        config['variant_counts'],
        config['barcode_variant_table'],
        config['barcode_expression']
    output:
        config['escape_fracs_barcodes'],
        md='results/summary/compute_barcode_escape.md'
    envmodules:
        'R/3.6.2-foss-2019b'
    params:
        nb='compute_barcode_escape.Rmd',
        md='compute_barcode_escape.md'
    shell:
        """
        R -e \"rmarkdown::render(input=\'{params.nb}\')\";
        mv {params.md} {output.md};
        """

rule count_variants:
    """Count codon variants from Illumina barcode runs."""
    input:
        config['barcode_variant_table'],
        config['barcode_runs']
    output:
        config['variant_counts'],
        nb_markdown=nb_markdown('count_variants.ipynb')
    params:
        nb='count_variants.ipynb'
    shell:
        "python scripts/run_nb.py {params.nb} {output.nb_markdown}"
