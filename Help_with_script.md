## Instruction for using the script for assessing differential gene expression (DEG) adjusted for multiple comparisons

Here is get_result_table function that require 3 arguments: first_cell_type_expressions_path, second_cell_type_expressions_path, save_results_table.
If you want to analyze **DEG without multiple comparisons**, you have to give only these 3 arguments to function. 
In such case differences between mean gene expressions between cell lines will be calculated using Z-test and Ci-test.

If it is neseccary for you to include multiple comparisons to your work, you have to give 1 or 2 arguments more to this function: 

1. **multiply_comp_method** equal one of next methods (in quotation marks) (necessary to give):

bonferroni : one-step correction

sidak : one-step correction

holm-sidak : step down method using Sidak adjustments

holm : step-down method using Bonferroni adjustments

simes-hochberg : step-up method (independent)

hommel : closed method based on Simes tests (non-negative)

fdr_bh : Benjamini/Hochberg (non-negative)

fdr_by : Benjamini/Yekutieli (negative)

fdr_tsbh : two stage fdr correction (non-negative)

fdr_tsbky : two stage fdr correction (non-negative)

2. **adj_alpha** - FWER, family-wise error rate (optional); default value = 0.05

### Output will be a table, where the first column - gene names, then:

1. **DEG without multiple comparisons**: Ci-test results, Z-test p-values, Z-test results, difference in mean expressions of 2 genes.
2. **DEG with multiple comparisons**: Z-test p-values, multiply adjusted results,, adjusted p-values, difference in mean expressions of 2 genes.
