# Project 2

## Important Information

1. Due electronically by **5:00pm, Tuesday, October 16, 2020**. As stated in the Tranche 2 Overview, 10% will be taken off with each successive late day. No submissions will be accepted after 3 days have passed.
2. Programming projects must be completed individually. You may discuss algorithms with others, but the coding should be done alone. You must explicitly name everyone with whom you discussed this project in the header comments of your code and in the collaboration attestation in the project quiz. Students must abide by the terms of the [Stanford Honor Code.Links to an external site.](https://communitystandards.stanford.edu/student-conduct-process/honor-code-and-fundamental-standard)
3. Remember to consult [Piazza (Links to an external site.)](https://piazza.com/class/j66r25o2adf4vp?), as many common questions will be asked and answered there.
4. Prior to beginning the assignment please read the [Code Policy](https://canvas.stanford.edu/courses/122868/pages/code-policy).
5. **Do not use scikit-learn**. A2 introduced this package to get you familiar with an off-the-shelf machine learning toolkit. However, this project is meant to be instructive for implementation and gaining intuition for how such methods and frameworks work.
6. Pay careful attention to instructions for naming and formatting output files. You will not receive full credit if you do not follow these instructions. We've made this autograder more flexible in how you implement the solutions, but you need to follow our directions below on the necessary class methods.
7. You will need to write additional code to answer some of the quiz questions. This code does not need to be submitted, and will not be graded. The code you submit should not perform these additional tasks. Source code must run **exactly** as specified below. Do not output any additional text to the command line.
8. **We do not provide starter code** or test cases for you to use in debugging this project. We provided these for Project 1 because that project is very difficult to test and debug *de novo*. If you understand these algorithms, you will be able to create your own test cases to check the correctness of your code. 

 

## Files to Download

You can download all files required for the project [here](https://canvas.stanford.edu/courses/122868/files/6438148/download?wrap=1). The three files provided are described in the Instructions section.

## Introduction

Each cell in the human body has the same genetic code, yet perform a diverse array of different specialized functions. From the central dogma, we know that the genetic code (DNA) is used to transcribe (or “express”) genes to RNA transcripts (messenger RNA or mRNA) that carry information about which proteins are made in a cell. Most mRNA then goes on to the ribosome to be translated into proteins. But given the availability of the mRNA as intermediate step, gene expression assays can be used to understand cell and/or tissue function, by assessing which genes are being transcribed in the cell at a given time. Gene expression can also provide insight into which genes may be involved in disease processes. In this assignment, you will analyze gene expression data collected from ectopic endometrium from healthy donors and patients with endometriosis using two computational methods: classification via K-Nearest Neighbors (KNN) and Gene Set Enrichment Analysis (GSEA).

The dataset is available on Gene Expression Omnibus, accession number [GSE25628 (Links to an external site.)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25628). (Take the time read the description at the link). In this analysis, we will compare expression from healthy controls and endometriosis patients. Endometriosis is characterized by presence of endometrial tissue outside the uterine cavity and is frequently associated with pain and infertility in women. We downloaded the expression data from this study and preprocessed it by filtering for lowly expressed genes, mapping probes to gene names, and collapsing across probes that map to the same gene. This was done using the [GEOparse (Links to an external site.)](https://geoparse.readthedocs.io/en/latest/) package.

You will use this dataset to first assess the ability of gene expression values to differentiate the healthy donor endometrium from the pathological endometrial tissue. You will implement the KNN algorithm (outlined in class) to classify each sample as disease or normal.

In the second part, you will implement GSEA (as discussed in lecture) to identify gene sets that are significantly enriched in pathological samples to get an idea of potentially activated pathways involved in the disease. We will be using pathways from the KEGG database. This is just one example of a set of pathways, there are many other pathway resources that may provide additional information and consider different sets of genes (go to [MSigDB (Links to an external site.)](http://software.broadinstitute.org/gsea/msigdb/collections.jsp#C2) if you are interested). The results from this type of analysis could be helpful downstream in understanding which genes play a role in endometriosis development and treatment.

Crispi, S., Piccolo, M. T., D'avino, A., Donizetti, A., Viceconte, R., Spyrou, M., ... & Signorile, P. G. (2013). Transcriptional profiling of endometriosis tissues identifies genes related to organogenesis defects. *Journal of cellular physiology*, 228(9), 1927-1934.

Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., ... & Mesirov, J. P. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *Proceedings of the National Academy of Sciences*, 102(43), 15545-15550.

**Learning goals**

- Understand how to implement and evaluate a classification algorithm (KNN).
- Implement the GSEA algorithm, and understand the concepts of random walks and bootstrapping.
- Understand how different (hyper)parameter choices impact algorithm performance.
- Learn about the potential and the limitations of biological interpretation of the results.

##  

## Instructions

Any data included in the output formatting examples are made up and should not be used to check code accuracy. Please pay attention to formatting requirements. Once again, **do not use built-in functions from sklearn for this project – this will be checked by our style graders and will result in loss of credit**. You can use the additional packages numpy, pandas and scipy, if desired, but this is not required.

####  

#### **Part 1: KNN**

1. Read in the expression file, `GSE25628_filtered_expression.txt`, and sample file, `GSE25628_samples.txt`, available on Canvas. The file names for this input should be command line arguments.

   1. - Run your script as:

        ```
        python3 knn.py expfile sampfile
        ```

      - The expression file is tab-delimited and contains a row for each gene and a column for each sample. We recommend checking out the python library [pandas (Links to an external site.)](https://pandas.pydata.org/) for functions to read in and analyze data matrices, if interested.

      - The samples file is tab-delimited and contains a row for each sample. The second column is the sample label, with a 0 indicating control and a 1 indicating patient.

2. Write a class called KNN to implement the algorithm to classify each sample as healthy or patient using K-Nearest Neighbors. We will use leave-one-out cross-validation (LOOCV) to assess performance, so for each sample, consider all other possible samples as neighbors based on Euclidean distance, and then compute accuracy metrics after assigning each individual based on its top $k$ neighbors, specifying the fraction needed for a positive classification, $fn$. Your code should work for any value of $k<$ the number of samples and $fn$ in $[0,1]$. Assign each sample as a patient (1) if **greater than** $fn \times k$ nearest neighbors are patients and healthy (0) otherwise. Your class can have as many functions as you want, but the following <u>functionality must be present</u> for grading purposes:

1. 1. `KNN.load_data`(`expfile`, `sampfile`) should take the paths to the expression and sample file, read them in, and store within your KNN class.
   2. `KNN.get_assignments`($k$, $fn$) should return the class assignments for all samples for given values of $k$ and $fn$ as a list of integer 0s and 1s.
   3. `KNN.calc_metrics`($k$, $fn$) should return a list of float values [sensitivity,specificity] of a KNN classifier using the given values of $k$ and $fn$. This method should run get_assignments at some point to initiate assignments used for performance metrics. 

####  

#### **Part 2. GSEA**

1. Read in the file containing KEGG gene set assignments, available on Canvas. This file name should be a command line argument.

   - - Run your script as:

       ```
       python3 gsea.py expfile sampfile keggfile
       ```

     - The KEGG file is tab-delimited and contains a row for each gene set with the name and associated genes.

     - The first column is the gene set name, the second is a link to additional info (can be ignored), and then each subsequent column contains the gene names included in that set. The number of columns with data per row will vary depending on the size of the gene set.

2. Create a class, named GSEA, to implement a given gene set enrichment analysis. Initialize the class with the expression data, sample assignments and list of gene sets. GSEA should be able to perform the following:

   1. 1. Rank genes by their log fold-change between healthy donors and endometriosis patients by taking the log2 of the ratio of mean endometriosis patient expression to mean healthy expression per gene. As the values have already been log normalized, subtracting the values will get you this ratio (you do not need to re-log-normalize the values): $\log{FC_{gene}} = \log(Patient_{gene}) -\log(Control_{gene})$
      2. Calculate the enrichment score (ES) of a gene set in the ranked gene list. You can exclude KEGG genes not in the expression data. This should implement the Brownian bridge discussed in class. As a reminder, the steps to calculate the ES are below, and this link may be helpful, Enrichment score steps](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3205944/#SD1).
         - (a) Calculate the appropriate up or down score for this particular gene set, controlling for the number of genes.
         - (b) Move down the ranked list of genes, adding to the running sum the value from (a) if the gene is in the set of interest, or subtracting from the running sum if not.
         - (c) Find the supremum of the brownian bridge as the ES. We are interested in gene sets that are enriched in up-regulated genes here, so do not take the absolute value in this case.
      3. *For your debugging only*: output the enrichment score (rounded to 2 decimal places) of each gene set, sorted with the highest score first, in a tab-delimited file named `kegg_enrichment_scores.txt` as follows (again, no header): '`KEGG_CITRATE_CYCLE_TCA_CYCLE 4.20`'
      4. Create a background distribution of enrichment scores by permuting the sample labels and re-calculating the ranked gene list 100 times and calculate a p-value for each gene set by counting the number of times in the permuted iterations a gene set had an equal or higher score than its actual score over the total number of iterations.
      5. Correct p-values for number of tests (equal to the number of gene sets) via Bonferroni

3. You can include as many functions or additional classes as you would like within to accomplish the above, but the <u>following functionality must be present</u> for grading purposes:

   1. 1. `GSEA.load_data`(`expfile`, `sampfile`, `genesets`) should take the file paths to the expression, sample and gene set data, read them in and store within the GSEA instance.
      2. `GSEA.get_gene_rank_order()` should return a list of all genes (as strings) ranked by their logFC between patient and control, with the gene with the highest logFC ordered first.
      3. `GSEA.get_enrichment_score(geneset)` should return the enrichment score, a float correct to two decimal places, for a given gene set, such as ‘KEGG_CITRATE_CYCLE_TCA_CYCLE’. This method should run get_gene_rank_order at some point to initiate enrichment calculations.
      4. `GSEA.get_sig_sets(p)` should return the list of significant gene sets (as strings), at a corrected threshold of p, by name. If no gene sets are significant, return an empty list. This method should run get_gene_rank_order and/or get_enrichment_score at some point to initiate enrichment calculations and then identify significant gene sets.

You will submit two scripts, which should be able to run as:

```
python3 knn.py expinputfile sampleinputfile
python3 gsea.py expinputfile sampleinputfile kegginputfile
```

This should produce the output as described above. No additional submission files are necessary. You can include additional functions to those described above, but make sure all the functions described are present, as the autograder will use these functions to assign points, as in the below rubric. All code should be well-documented, refer to the [style guidelines](https://canvas.stanford.edu/courses/122868/pages/code-policy).

Also note that we don't require any output written to files or on the command line. We are testing the functionality of the specified class methods only. This is in hopes of avoiding output file writing discrepancies. For debugging, we suggest using a main function to call your object's methods.

 

**Rubric**

We will use the class functions described above to implement the following tests for grading. For KNN, we will test with set values of ![LaTeX: k](https://canvas.stanford.edu/equation_images/k)k and ![LaTeX: fn](https://canvas.stanford.edu/equation_images/fn)f n in the ranges described above. You can earn partial credit if you pass a subset of the tests in a given section. Please pay attention to any rounding or output formatting instructions above.

Code total - 60 points

1. Correct implementation and evaluation of KNN - 25 points

Scoring tests

​		1. All required functions are present for the KNN class – 5 points

​		2. Correct classifications for several parameter values – 10 points

​		3. Correct sensitivity calculation for several parameter values – 5 points

​		4. Correct specificity calculation for several parameter values – 5 points

2. Correct implementation of GSEA - 35 points

Scoring tests:

​		1. All required functions are present for the GSEA class – 5 points

​		2. Differential expression ranking correct – 5 points

​		3. Enrichment scores correct – 15 points

​		4. Number of significant sets in correct range – 5 points

​		5. Significance calculation runs in under 3 minutes – 5 points

3. Style - 10 points

------

## Quiz

You will run your code to answer the questions in the quiz on Canvas. In order to answer some questions you may have to write additional functions. Any code you write for the sole purpose of answering the quiz questions should not be submitted.

------

## Submission Instructions

1. Take the quiz on Canvas.
2. Submit your scripts knn.py and gsea.py via Gradescope. Remember, we expect well-commented code.

**Please be sure you have followed all of the project directives as a portion of your grade is based on your compliance with these directives, and we use an autograder that is expecting you to follow the above instructions. You may not receive credit otherwise.**

**==PSA: the autograder is not a valid debugging strategy==**! Please write your own tests for sanity checking. This autograder takes longer to run this time around, so be cognizant with your time. Furthermore, all **functions are tested independently**, so please call previously-defined "setup" functions and utilize their returned values in downstream functions (see specific function specs above).