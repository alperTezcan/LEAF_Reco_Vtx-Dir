# Direction Reconsturction with LEAF

The Low Energy Algoritm Framework algorithm is an alternative and simple LE fitter than can be used for HyperK and SuperK. For further information: https://github.com/hyperk/LEAF/tree/master

In scope of my Bachelor's thesis, I demonstrated a method for direction reconstruction and completed certain performance tests. I presented the idea of the algorithm and the performance test results in the form of a report.

- **Analysis/**: Some analysis tools for generating the plots or getting positional information of events.
- **LEAF/**: Codes for the algoritm
- **ProducePDF_revamped/**: The improved version of the PDF generatring code for LEAF.
- **Report.pdf**: The thesis report

# How to:
1. Source RunAtStart.sh after you updated your ROOT directory.
2. Use the script ./SetupDataModel.sh to define the DataModel (if you have hk-AstroAnalysis, you should setup the global variable)
2. Enter the leaf/ repository and make clean;make
3. Enter the example repository and make clean;make
4. One example of how to run the code is set in example: test_example.sh
6. You can use shell scripts in shell/ in order to run the fitter or launch on batch.

# Useful scripts in ./macros and ./shell
You can compile with GNUMake like following in ./macros:
```
$ make ProducePDF
```
## Making tuning file (e.g. ./inputs/timePDF_Directionality_DRnew.root)
1. Produce plots by AnalyzeWSHierarchy: reads out WCSim output and makes plots.
2. Produce time PDF (and angular PDF) by ProducePDF: uses plots made by AnalyzeWSHierarchy and generate PDFs for LEAF.

```
$ AnalyzeWSHierarchy -f wcrim_hybrid.root -o plots.root
$ ProducePDF -f plots.root -o PDF.root
```

## Making generic plots
- LEAFOutputAnalysisHybrid_leafclass: read LEAF output to produce generic plots. If one uses the master branch for LEAF, please use LEAFOutputAnalysisHybrid_master