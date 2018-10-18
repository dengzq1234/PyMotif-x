
# PyMotif-x

Biological sequence motif discovery in Python

# Introdcution

PyMotif-x is a useable implementation the motif-x tool in the Python programming language. Motif-x is a software tool designed to extract overrepresented patterns from any sequence data set. The algorithm is an iterative strategy which builds successive motifs through comparison to a dynamic statistical background. For more information please refer to the motif-x resouce. The current implementation only supports sequences with a fixed length (i.e pre-aligned) and have a fixed central residue. For example, post-translational modification sites.

# Installation

PyMotif-x can be directly installed from github 


```python
git clone https://github.com/dengzq1234/PyMotif-x.git
```

# Usage

The script contains the function motifx is the main function. Before you start a run, you will need a foreground and background dataset of sequences. The sample dataset are extracted from another re-implementation, rmotifx. If you are a R user, please refer to rmotifx [resource](https://github.com/omarwagih/rmotifx) and its research.

The package has attached sample data for test run:
fg-data-test.txt
bg-data-test.txt 

Here, the foreground data represents phosphorylation binding sites of Casein Kinase 2. The negative data represents 10,000 random serine-centered sites.

To start the progrom, run the following in bash command:


```bash
python motifx-functions.py <foreground_filename> <background_filename> <central_residue> <minimum_ocurrences> <p_value_cutoff> <verbose> <perl_minimum_activation>
```



running on example dataset:



```bash
python motifx-functions.py fg-data-test.txt bg-data-test.txt S 20 1e-6 False False
```

The resutls returned should have the following formats:

```bash
   motif_score      motif_regex  fg_matches  fg_size  bg_matches  bg_size  fold_increase
0    31.909180  .......SD.E....          57      399          23     6039      37.509317
1    26.679736  .......S..EE...          37      342          37     6016      17.590643
2    31.909180  .......SD.D....          39      305          12     5979      63.710656
3    23.062438  .......SE.E....          24      266          32     5967      16.824248
4    15.954590  .......S..E....          56      242         325     5935       4.225811
5    24.168438  .......SE.D....          21      186          26     5610      24.361042
6    10.915342  .......S..D....          30      165         233     5584       4.357394
7     9.371511  .......SD......          25      135         224     5351       4.423776
8     7.270142  .......S.E.....          25      110         342     5127       3.407097
```

For detailed explanations of all parameters and output, check out the documentation by typing man motifx. You can also refer to the original motif-x [resource](http://motif-x.med.harvard.edu/). If you want to check its usage in R, please check: 

Wagih O, Sugiyama N, Ishihama Y, Beltrao P. (2015) Uncovering phosphorylation-based specificities through functional interaction networks (2015). Mol. Cell. Proteomics PUBMED

# Further challenges

* Allow motif discovery in non-centered k-mers
* Support more background proteome data
* Include statistic estimation

# Feedback

If you have any feedback or suggestions, please left me a message at dengziqi1234(at)gmail.com) or open an issue on github.
