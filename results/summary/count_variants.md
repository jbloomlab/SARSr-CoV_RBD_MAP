# Count variants in each sample
This Python Jupyter notebook counts occurrences of each barcode in each sample from Illumina barcode sequencing, and adds these counts to the variant table.

## Set up analysis
### Import Python modules.
Use [plotnine](https://plotnine.readthedocs.io/en/stable/) for ggplot2-like plotting.

The analysis relies heavily on the Bloom lab's [dms_variants](https://jbloomlab.github.io/dms_variants) package:


```python
import itertools
import multiprocessing
import multiprocessing.pool
import os
import warnings

import alignparse
import alignparse.targets

import dms_variants.codonvarianttable
from dms_variants.constants import CBPALETTE
import dms_variants.illuminabarcodeparser
import dms_variants.utils
import dms_variants.plotnine_themes

from IPython.display import display, HTML

import pandas as pd

from plotnine import *

import yaml
```

Set [plotnine](https://plotnine.readthedocs.io/en/stable/) theme to the gray-grid one defined in `dms_variants`:


```python
theme_set(dms_variants.plotnine_themes.theme_graygrid())
```

Versions of key software:


```python
print(f"Using alignparse version {alignparse.__version__}")
print(f"Using dms_variants version {dms_variants.__version__}")
```

    Using alignparse version 0.1.6
    Using dms_variants version 0.8.5


Ignore warnings that clutter output:


```python
warnings.simplefilter('ignore')
```

### Parameters for notebook
Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Make output directory if needed:


```python
os.makedirs(config['counts_dir'], exist_ok=True)
```

## Input variant table
Initialize the input table from the `process_ccs` analysis notebook:


```python
variants = pd.read_csv(config['barcode_variant_table'])

display(HTML(variants.head().to_html(index=False)))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>target</th>
      <th>library</th>
      <th>barcode</th>
      <th>substitutions</th>
      <th>variant_call_support</th>
      <th>number_of_indels</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>AncSARS1a_tree1</td>
      <td>lib1</td>
      <td>AAAAAAAAAGTGAAAG</td>
      <td>NaN</td>
      <td>3</td>
      <td>0</td>
    </tr>
    <tr>
      <td>ZXC21</td>
      <td>lib1</td>
      <td>AAAAAAAAAGTTACTA</td>
      <td>NaN</td>
      <td>4</td>
      <td>0</td>
    </tr>
    <tr>
      <td>HeB2013</td>
      <td>lib1</td>
      <td>AAAAAAAAATGAGGAC</td>
      <td>NaN</td>
      <td>3</td>
      <td>0</td>
    </tr>
    <tr>
      <td>AncAsia_tree2</td>
      <td>lib1</td>
      <td>AAAAAAAAATGCCATG</td>
      <td>NaN</td>
      <td>3</td>
      <td>0</td>
    </tr>
    <tr>
      <td>AncSARS-CoV-1_alt</td>
      <td>lib1</td>
      <td>AAAAAAAACACTTAGA</td>
      <td>NaN</td>
      <td>1</td>
      <td>0</td>
    </tr>
  </tbody>
</table>


Pull out a target sequence for matching to the barcode and flanking sequence regions.


```python
# get wildtype gene sequence for primary target
targets = alignparse.targets.Targets(seqsfile=config['amplicons'],
                                     feature_parse_specs=config['feature_parse_specs'])
```

## Setup to parse barcodes
Read data frame with list of all barcode runs.
Note how multiple R1 files are delimited by `; ` and are split out separately:


```python
print(f"Reading list of barcode runs from {config['barcode_runs']}")

barcode_runs = (pd.read_csv(config['barcode_runs'])
                .assign(R1=lambda x: x['R1'].str.split('; '))
                )
      
display(HTML(barcode_runs.to_html(index=False)))
```

    Reading list of barcode runs from data/barcode_runs.csv



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>date</th>
      <th>experiment</th>
      <th>library</th>
      <th>antibody</th>
      <th>concentration</th>
      <th>concentration_units</th>
      <th>group</th>
      <th>selection</th>
      <th>sample</th>
      <th>frac_escape</th>
      <th>cells_sorted</th>
      <th>R1</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>201106</td>
      <td>expt_68-73</td>
      <td>lib1</td>
      <td>none</td>
      <td>0</td>
      <td>ng_per_mL</td>
      <td>Vir</td>
      <td>reference</td>
      <td>expt_68-73-none-0-reference</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>[/shared/ngs/illumina/tstarr/201123_D00300_1119_BHJKJTBCX3/Unaligned/Project_tstarr/expt68_73_hom_ref_S119_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/201116_D00300_1113_BHJKKFBCX3/Unaligned/Project_tstarr/expt68_73_hom_ref_S23_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>201106</td>
      <td>expt_68</td>
      <td>lib1</td>
      <td>S309</td>
      <td>421</td>
      <td>ng_per_mL</td>
      <td>Vir</td>
      <td>escape</td>
      <td>expt_68-S309-421-escape</td>
      <td>0.490</td>
      <td>858224.0</td>
      <td>[/shared/ngs/illumina/tstarr/201123_D00300_1119_BHJKJTBCX3/Unaligned/Project_tstarr/expt70_lib1_Abneg_S126_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/201116_D00300_1113_BHJKKFBCX3/Unaligned/Project_tstarr/expt70_lib1_Abneg_S30_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>201106</td>
      <td>expt_69</td>
      <td>lib1</td>
      <td>S2E12</td>
      <td>56</td>
      <td>ng_per_mL</td>
      <td>Vir</td>
      <td>escape</td>
      <td>expt_69-S2E12-56-escape</td>
      <td>0.890</td>
      <td>1754530.0</td>
      <td>[/shared/ngs/illumina/tstarr/201123_D00300_1119_BHJKJTBCX3/Unaligned/Project_tstarr/expt71_lib1_Abneg_S129_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/201116_D00300_1113_BHJKKFBCX3/Unaligned/Project_tstarr/expt71_lib1_Abneg_S33_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>201106</td>
      <td>expt_70</td>
      <td>lib1</td>
      <td>S2X35</td>
      <td>70</td>
      <td>ng_per_mL</td>
      <td>Vir</td>
      <td>escape</td>
      <td>expt_70-S2X35-70-escape</td>
      <td>0.444</td>
      <td>864623.0</td>
      <td>[/shared/ngs/illumina/tstarr/201123_D00300_1119_BHJKJTBCX3/Unaligned/Project_tstarr/expt72_lib1_Abneg_S132_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/201116_D00300_1113_BHJKKFBCX3/Unaligned/Project_tstarr/expt72_lib1_Abneg_S36_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>201106</td>
      <td>expt_72</td>
      <td>lib1</td>
      <td>S2X58</td>
      <td>18</td>
      <td>ng_per_mL</td>
      <td>Vir</td>
      <td>escape</td>
      <td>expt_72-S2X58-18-escape</td>
      <td>0.925</td>
      <td>1770393.0</td>
      <td>[/shared/ngs/illumina/tstarr/201123_D00300_1119_BHJKJTBCX3/Unaligned/Project_tstarr/expt72_homolog_Abneg_S134_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/201116_D00300_1113_BHJKKFBCX3/Unaligned/Project_tstarr/expt72_homolog_Abneg_S38_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>201106</td>
      <td>expt_73</td>
      <td>lib1</td>
      <td>S304</td>
      <td>46</td>
      <td>ng_per_mL</td>
      <td>Vir</td>
      <td>escape</td>
      <td>expt_73-S304-46-escape</td>
      <td>0.412</td>
      <td>817857.0</td>
      <td>[/shared/ngs/illumina/tstarr/201123_D00300_1119_BHJKJTBCX3/Unaligned/Project_tstarr/expt73_homolog_Abneg_S137_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/201116_D00300_1113_BHJKKFBCX3/Unaligned/Project_tstarr/expt73_homolog_Abneg_S41_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>201109</td>
      <td>expt_74-78</td>
      <td>lib1</td>
      <td>none</td>
      <td>0</td>
      <td>ng_per_mL</td>
      <td>Vir</td>
      <td>reference</td>
      <td>expt_74-78-none-0-reference</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>[/shared/ngs/illumina/tstarr/201123_D00300_1119_BHJKJTBCX3/Unaligned/Project_tstarr/expt74_78_hom_ref_S161_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/201116_D00300_1113_BHJKKFBCX3/Unaligned/Project_tstarr/expt74_78_hom_ref_S65_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>201109</td>
      <td>expt_74</td>
      <td>lib1</td>
      <td>S2H58</td>
      <td>46</td>
      <td>ng_per_mL</td>
      <td>Vir</td>
      <td>escape</td>
      <td>expt_74-S2H58-46-escape</td>
      <td>0.893</td>
      <td>1434377.0</td>
      <td>[/shared/ngs/illumina/tstarr/201123_D00300_1119_BHJKJTBCX3/Unaligned/Project_tstarr/expt74_homolog_Abneg_S140_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/201116_D00300_1113_BHJKKFBCX3/Unaligned/Project_tstarr/expt74_homolog_Abneg_S44_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>201109</td>
      <td>expt_75</td>
      <td>lib1</td>
      <td>S2D106</td>
      <td>68</td>
      <td>ng_per_mL</td>
      <td>Vir</td>
      <td>escape</td>
      <td>expt_75-S2D106-68-escape</td>
      <td>0.890</td>
      <td>1494533.0</td>
      <td>[/shared/ngs/illumina/tstarr/201123_D00300_1119_BHJKJTBCX3/Unaligned/Project_tstarr/expt75_homolog_Abneg_S143_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/201116_D00300_1113_BHJKKFBCX3/Unaligned/Project_tstarr/expt75_homolog_Abneg_S47_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>201109</td>
      <td>expt_77</td>
      <td>lib1</td>
      <td>S2X16</td>
      <td>54</td>
      <td>ng_per_mL</td>
      <td>Vir</td>
      <td>escape</td>
      <td>expt_77-S2X16-54-escape</td>
      <td>0.891</td>
      <td>1247076.0</td>
      <td>[/shared/ngs/illumina/tstarr/201123_D00300_1119_BHJKJTBCX3/Unaligned/Project_tstarr/expt77_homolog_Abneg_S149_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/201116_D00300_1113_BHJKKFBCX3/Unaligned/Project_tstarr/expt77_homolog_Abneg_S53_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>201109</td>
      <td>expt_78</td>
      <td>lib1</td>
      <td>S2X227</td>
      <td>138</td>
      <td>ng_per_mL</td>
      <td>Vir</td>
      <td>escape</td>
      <td>expt_78-S2X227-138-escape</td>
      <td>0.933</td>
      <td>1307137.0</td>
      <td>[/shared/ngs/illumina/tstarr/201123_D00300_1119_BHJKJTBCX3/Unaligned/Project_tstarr/expt78_homolog_Abneg_S152_R1_001.fastq.gz, /shared/ngs/illumina/tstarr/201116_D00300_1113_BHJKKFBCX3/Unaligned/Project_tstarr/expt78_homolog_Abneg_S56_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>211113</td>
      <td>expt_79-81</td>
      <td>lib1</td>
      <td>none</td>
      <td>0</td>
      <td>ng_per_mL</td>
      <td>Vir</td>
      <td>reference</td>
      <td>expt_79-81-none-0-reference</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>[/shared/ngs/illumina/tstarr/210115_D00300_1154_BHK3HCBCX3/Unaligned/Project_tstarr/expt_79-81_hom_ref_S85_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>211113</td>
      <td>expt_79</td>
      <td>lib1</td>
      <td>S2H97</td>
      <td>58</td>
      <td>ng_per_mL</td>
      <td>Vir</td>
      <td>escape</td>
      <td>expt_79-S2H97-58-escape</td>
      <td>0.056</td>
      <td>107488.0</td>
      <td>[/shared/ngs/illumina/tstarr/210115_D00300_1154_BHK3HCBCX3/Unaligned/Project_tstarr/expt79_hom_Abneg_S88_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>211113</td>
      <td>expt_79strin</td>
      <td>lib1</td>
      <td>S2H97_stringent</td>
      <td>58</td>
      <td>ng_per_mL</td>
      <td>Vir</td>
      <td>escape</td>
      <td>expt_79strin-S2H97_stringent-58-escape</td>
      <td>0.007</td>
      <td>9624.0</td>
      <td>[/shared/ngs/illumina/tstarr/210115_D00300_1154_BHK3HCBCX3/Unaligned/Project_tstarr/expt79_hom_string_S89_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>211113</td>
      <td>expt_80</td>
      <td>lib1</td>
      <td>S2H13</td>
      <td>56</td>
      <td>ng_per_mL</td>
      <td>Vir</td>
      <td>escape</td>
      <td>expt_80-S2H13-56-escape</td>
      <td>0.878</td>
      <td>1364359.0</td>
      <td>[/shared/ngs/illumina/tstarr/210115_D00300_1154_BHK3HCBCX3/Unaligned/Project_tstarr/expt80_hom_Abneg_S92_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>211113</td>
      <td>expt_80strin</td>
      <td>lib1</td>
      <td>S2H13_stringent</td>
      <td>56</td>
      <td>ng_per_mL</td>
      <td>Vir</td>
      <td>escape</td>
      <td>expt_80strin-S2H13_stringent-56-escape</td>
      <td>0.843</td>
      <td>721340.0</td>
      <td>[/shared/ngs/illumina/tstarr/210115_D00300_1154_BHK3HCBCX3/Unaligned/Project_tstarr/expt80_hom_string_S93_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>211113</td>
      <td>expt_81</td>
      <td>lib1</td>
      <td>S2H14</td>
      <td>105</td>
      <td>ng_per_mL</td>
      <td>Vir</td>
      <td>escape</td>
      <td>expt_81-S2H14-105-escape</td>
      <td>0.986</td>
      <td>1874644.0</td>
      <td>[/shared/ngs/illumina/tstarr/210115_D00300_1154_BHK3HCBCX3/Unaligned/Project_tstarr/expt81_hom_Abneg_S96_R1_001.fastq.gz]</td>
    </tr>
    <tr>
      <td>211113</td>
      <td>expt_81strin</td>
      <td>lib1</td>
      <td>S2H14_stringent</td>
      <td>105</td>
      <td>ng_per_mL</td>
      <td>Vir</td>
      <td>escape</td>
      <td>expt_81strin-S2H14_stringent-105-escape</td>
      <td>0.955</td>
      <td>959111.0</td>
      <td>[/shared/ngs/illumina/tstarr/210115_D00300_1154_BHK3HCBCX3/Unaligned/Project_tstarr/expt81_hom_string_S97_R1_001.fastq.gz]</td>
    </tr>
  </tbody>
</table>


Make sure library / sample combinations are unique:


```python
assert len(barcode_runs) == len(barcode_runs.groupby(['library', 'sample']))
```

Make sure the the libraries for which we have barcode runs are all in our variant table:


```python
unknown_libs = set(barcode_runs['library']) - set(variants['library'])
if unknown_libs:
    raise ValueError(f"Libraries with barcode runs not in variant table: {unknown_libs}")
```

Now we initialize an [IlluminaBarcodeParser](https://jbloomlab.github.io/dms_variants/dms_variants.illuminabarcodeparser.html#dms_variants.illuminabarcodeparser.IlluminaBarcodeParser) for each library.

First, get the length of the barcode from the alignment target after making sure the same length for all targets:


```python
bclen = len(targets.targets[0].get_feature('barcode').seq)

assert (bclen == len(target.get_feature('barcode').seq) for target in targets.targets)

print(f"Barcodes of length {bclen}")
```

    Barcodes of length 16


The other barcode parsing params come from the config file:


```python
parser_params = config['illumina_barcode_parser_params']

display(HTML(
    pd.Series(parser_params, name='value')
    .rename_axis(index='parameter')
    .reset_index()
    .to_html(index=False)
    ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>parameter</th>
      <th>value</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>upstream</td>
      <td>GGCCGC</td>
    </tr>
    <tr>
      <td>downstream</td>
      <td></td>
    </tr>
    <tr>
      <td>minq</td>
      <td>20</td>
    </tr>
    <tr>
      <td>upstream_mismatch</td>
      <td>1</td>
    </tr>
    <tr>
      <td>downstream_mismatch</td>
      <td>0</td>
    </tr>
  </tbody>
</table>


The parser needs to know the set of valid barcodes, which are stored in the variant table and are different for each library.
So we create a different parser for each library using these valid barcode sets:


```python
# create dict keyed by library, value is parser for library
parsers = {lib: dms_variants.illuminabarcodeparser.IlluminaBarcodeParser(
                    bclen=bclen,
                    valid_barcodes=variants.loc[variants['library']==lib]['barcode'],
                    **parser_params)
           for lib in set(variants['library'])}

print('Number of valid barcodes searched for by each parser:')
display(HTML(
    pd.DataFrame([(lib, len(p.valid_barcodes)) for lib, p in parsers.items()],
                 columns=['library', 'number of valid barcodes'])
    .to_html(index=False)
    ))
```

    Number of valid barcodes searched for by each parser:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>library</th>
      <th>number of valid barcodes</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>lib2</td>
      <td>23259</td>
    </tr>
    <tr>
      <td>lib1</td>
      <td>25896</td>
    </tr>
  </tbody>
</table>


## Parse barcodes
We now parse the barcodes.
Since this will take a while, we utilize multiple CPUs via the Python [multiprocessing](https://docs.python.org/3.6/library/multiprocessing.html) module.
First, determine how many CPUs to use.
We use the minimum of the user-specified number hardcoded below and the number actually available.
(If you are running *interactively* on the Hutch cluster, you may need to reduce the number below in order to avoid an error as there is an enforced CPU limit on the home `rhino` nodes):


```python
ncpus = min(config['max_cpus'], multiprocessing.cpu_count())
print(f"Using {ncpus} CPUs")
```

    Using 16 CPUs


Parse the barcodes in parallel via a [multiprocessing.Pool](https://docs.python.org/3.6/library/multiprocessing.html#multiprocessing.pool.Pool) using all the available CPUs to get a list of the data frames with barcode counts / fates for each sample:


```python
def process_func(parser, r1files, library, sample):
    """Convenience function to be starmapped to multiprocessing pool."""
    return parser.parse(r1files, add_cols={'library': library, 'sample': sample})

# parallel computation of list of data frames
with multiprocessing.pool.Pool(processes=ncpus) as pool:
    bclist = pool.starmap(
                process_func,
                [(parsers[run.library], run.R1, run.library, run.sample)
                  for run in barcode_runs.itertuples()],
                )
```

Now concatenate the list into data frames of barcode counts and barcode fates:


```python
counts = pd.concat([samplecounts for samplecounts, _ in bclist],
                   sort=False,
                   ignore_index=True)

print('First few lines of counts data frame:')
display(HTML(counts.head().to_html(index=False)))

fates = pd.concat([samplefates for _, samplefates in bclist],
                  sort=False,
                  ignore_index=True)

print('First few lines of fates data frame:')
display(HTML(fates.head().to_html(index=False)))
```

    First few lines of counts data frame:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>barcode</th>
      <th>count</th>
      <th>library</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>AACGAAACCTCTGTCA</td>
      <td>2301</td>
      <td>lib1</td>
      <td>expt_68-73-none-0-reference</td>
    </tr>
    <tr>
      <td>ATTGGCTTACTAAATA</td>
      <td>2073</td>
      <td>lib1</td>
      <td>expt_68-73-none-0-reference</td>
    </tr>
    <tr>
      <td>CGATAATATAGACAAG</td>
      <td>2051</td>
      <td>lib1</td>
      <td>expt_68-73-none-0-reference</td>
    </tr>
    <tr>
      <td>CGCACAAACTAGTAGT</td>
      <td>1855</td>
      <td>lib1</td>
      <td>expt_68-73-none-0-reference</td>
    </tr>
    <tr>
      <td>AATGGAATATTCACAT</td>
      <td>1825</td>
      <td>lib1</td>
      <td>expt_68-73-none-0-reference</td>
    </tr>
  </tbody>
</table>


    First few lines of fates data frame:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>fate</th>
      <th>count</th>
      <th>library</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>valid barcode</td>
      <td>6480493</td>
      <td>lib1</td>
      <td>expt_68-73-none-0-reference</td>
    </tr>
    <tr>
      <td>invalid barcode</td>
      <td>2455034</td>
      <td>lib1</td>
      <td>expt_68-73-none-0-reference</td>
    </tr>
    <tr>
      <td>low quality barcode</td>
      <td>479114</td>
      <td>lib1</td>
      <td>expt_68-73-none-0-reference</td>
    </tr>
    <tr>
      <td>failed chastity filter</td>
      <td>213176</td>
      <td>lib1</td>
      <td>expt_68-73-none-0-reference</td>
    </tr>
    <tr>
      <td>unparseable barcode</td>
      <td>103814</td>
      <td>lib1</td>
      <td>expt_68-73-none-0-reference</td>
    </tr>
  </tbody>
</table>


## Examine fates of parsed barcodes
First, we'll analyze the "fates" of the parsed barcodes.
These fates represent what happened to each Illumina read we parsed:
 - Did the barcode read fail the Illumina chastity filter?
 - Was the barcode *unparseable* (i.e., the read didn't appear to be a valid barcode based on flanking regions)?
 - Was the barcode sequence too *low quality* based on the Illumina quality scores?
 - Was the barcode parseable but *invalid* (i.e., not in our list of variant-associated barcodes in the codon variant table)?
 - Was the barcode *valid*, and so will be added to variant counts.
 
First, we just write a CSV file with all the barcode fates:


```python
fatesfile = os.path.join(config['counts_dir'], 'barcode_fates.csv')
print(f"Writing barcode fates to {fatesfile}")
fates.to_csv(fatesfile, index=False)
```

    Writing barcode fates to results/counts/barcode_fates.csv


Next, we tabulate the barcode fates in wide format:


```python
display(HTML(fates
             .pivot_table(columns='fate',
                          values='count',
                          index=['library', 'sample'])
             .to_html()
             ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>fate</th>
      <th>failed chastity filter</th>
      <th>invalid barcode</th>
      <th>low quality barcode</th>
      <th>unparseable barcode</th>
      <th>valid barcode</th>
    </tr>
    <tr>
      <th>library</th>
      <th>sample</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="18" valign="top">lib1</th>
      <th>expt_68-73-none-0-reference</th>
      <td>213176</td>
      <td>2455034</td>
      <td>479114</td>
      <td>103814</td>
      <td>6480493</td>
    </tr>
    <tr>
      <th>expt_68-S309-421-escape</th>
      <td>56713</td>
      <td>520979</td>
      <td>123135</td>
      <td>42311</td>
      <td>1729027</td>
    </tr>
    <tr>
      <th>expt_69-S2E12-56-escape</th>
      <td>138430</td>
      <td>1345786</td>
      <td>314925</td>
      <td>67355</td>
      <td>4467484</td>
    </tr>
    <tr>
      <th>expt_70-S2X35-70-escape</th>
      <td>45302</td>
      <td>433349</td>
      <td>105003</td>
      <td>25151</td>
      <td>1459938</td>
    </tr>
    <tr>
      <th>expt_72-S2X58-18-escape</th>
      <td>131619</td>
      <td>1269771</td>
      <td>304638</td>
      <td>66349</td>
      <td>4206175</td>
    </tr>
    <tr>
      <th>expt_73-S304-46-escape</th>
      <td>56090</td>
      <td>558183</td>
      <td>130969</td>
      <td>30764</td>
      <td>1797869</td>
    </tr>
    <tr>
      <th>expt_74-78-none-0-reference</th>
      <td>279726</td>
      <td>3229902</td>
      <td>640755</td>
      <td>155651</td>
      <td>8478606</td>
    </tr>
    <tr>
      <th>expt_74-S2H58-46-escape</th>
      <td>101327</td>
      <td>946316</td>
      <td>233246</td>
      <td>52119</td>
      <td>3205580</td>
    </tr>
    <tr>
      <th>expt_75-S2D106-68-escape</th>
      <td>96732</td>
      <td>929023</td>
      <td>220951</td>
      <td>52851</td>
      <td>3166600</td>
    </tr>
    <tr>
      <th>expt_77-S2X16-54-escape</th>
      <td>99263</td>
      <td>956574</td>
      <td>233093</td>
      <td>49898</td>
      <td>3244146</td>
    </tr>
    <tr>
      <th>expt_78-S2X227-138-escape</th>
      <td>89899</td>
      <td>868823</td>
      <td>203773</td>
      <td>57331</td>
      <td>2839240</td>
    </tr>
    <tr>
      <th>expt_79-81-none-0-reference</th>
      <td>444944</td>
      <td>3239545</td>
      <td>546890</td>
      <td>209112</td>
      <td>8373676</td>
    </tr>
    <tr>
      <th>expt_79-S2H97-58-escape</th>
      <td>12319</td>
      <td>104863</td>
      <td>14768</td>
      <td>7395</td>
      <td>212145</td>
    </tr>
    <tr>
      <th>expt_79strin-S2H97_stringent-58-escape</th>
      <td>1848</td>
      <td>13069</td>
      <td>1197</td>
      <td>4526</td>
      <td>12221</td>
    </tr>
    <tr>
      <th>expt_80-S2H13-56-escape</th>
      <td>152085</td>
      <td>867629</td>
      <td>183161</td>
      <td>63659</td>
      <td>2917946</td>
    </tr>
    <tr>
      <th>expt_80strin-S2H13_stringent-56-escape</th>
      <td>75709</td>
      <td>470358</td>
      <td>91825</td>
      <td>33636</td>
      <td>1572057</td>
    </tr>
    <tr>
      <th>expt_81-S2H14-105-escape</th>
      <td>221026</td>
      <td>1277460</td>
      <td>254846</td>
      <td>100212</td>
      <td>4377215</td>
    </tr>
    <tr>
      <th>expt_81strin-S2H14_stringent-105-escape</th>
      <td>115497</td>
      <td>683625</td>
      <td>136974</td>
      <td>48677</td>
      <td>2308460</td>
    </tr>
  </tbody>
</table>


Now we plot the barcode-read fates for each library / sample, showing the bars for valid barcodes in orange and the others in gray.
We see that the largest fraction of barcode reads correspond to valid barcodes, and most of the others are invalid barcodes (probably because the map to variants that aren't present in our variant table since we didn't associate all variants with barcodes). The exception to this is lib2 Titeseq_03_bin3; the PCR for this sample in the original sequencing run failed, so we followed it up with a single MiSeq lane. We did not filter out the PhiX reads from this data before parsing, so these PhiX reads will deflate the fraction of valid barcode reads as expected, but does not indicate any problems.


```python
barcode_fate_plot = (
    ggplot(
        fates
        .assign(sample=lambda x: pd.Categorical(x['sample'],
                                                x['sample'].unique(),
                                                ordered=True),
                fate=lambda x: pd.Categorical(x['fate'],
                                              x['fate'].unique(),
                                              ordered=True),
                is_valid=lambda x: x['fate'] == 'valid barcode'
                ), 
        aes('fate', 'count', fill='is_valid')) +
    geom_bar(stat='identity') +
    facet_grid('sample ~ library') +
    facet_grid('sample ~ library') +
    scale_fill_manual(CBPALETTE, guide=False) +
    theme(figure_size=(1.4 * (1 + fates['library'].nunique()),
                       1.7 * (1.2 + fates['sample'].nunique())),
          axis_text_x=element_text(angle=90),
          panel_grid_major_x=element_blank()
          ) +
    scale_y_continuous(labels=dms_variants.utils.latex_sci_not,
                       name='number of reads')
    )

_ = barcode_fate_plot.draw()
```


    
![png](count_variants_files/count_variants_40_0.png)
    


## Output csv of barcode counts


```python
print(f"Writing variant counts to {config['variant_counts']}")
counts.to_csv(config['variant_counts'], index=False)
```

    Writing variant counts to results/counts/variant_counts.csv

