
# SRT_imputation

SRT_imputation offers the codes to reproduce the results of the publication. Also the gene-wise and spot(or cell)-wise evaluations which were done in the paper can be carried out for the users` own datasets.


## For running the reproduction

work in the ‘reproduction’ folder
```bash
from run_reprod import runQQ

runQQ()	#this step will take a while (~14 mins)

```

## For running the test and own datasets
**Input**
-	org_data: preprocessed (optionally) dataset in the dataframe structure (genes are at the columns)
-	imp_data: preprocessed (optionally) dataset in the dataframe structure (genes are at the columns)
-	g:  list of genes to use in evaluations
-	dName: name for dataset (str)
-	mName (optional): name for method (str)
-	pl_g (optional): the selected (at most 4) genes to plot (list)

**Output**
-	./QQeval/evaluations : folder to store the generated evaluations files for both gene-wise and spot(or cell)-wise
-	./ QQeval/figure : folder to store the generated figures

## For running the evaluations for test datasets
```bash
import QQevaluation as qq
org_df, imp_df, g = qq.test()
dName = 'myTest'
qq.QQevaluation(org_df, imp_df, g, dName, mName=None, pl_g =None)
```
## For running the evaluations for own datasets
```bash
import QQevaluation as qq
qq.QQevaluation(org_data, imp_data, g, dName, mName=None, pl_g=None)
```
## Lisans

[MIT](https://choosealicense.com/licenses/mit/)

  