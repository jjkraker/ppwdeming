# ppwdeming 2.0.0

* Sizable update, refining functions to operate more efficiently.
* PWD_RL is no longer called by other functions, but it is retained for individual 
  application. 
* PWD_get_gh is the major update, in which the original search algorithm is 
  replaced by an improved algorithm, which we refer to as the unified approach, 
  focusing on the ratio rho=sigma/kappa.
  This results in a many-fold decrease in computing speed.
* The interconnected functions now also incorporate an optional starting 
  point for all the parameters in the search algorithm:
  PWD_get_gh, PWD_inference, and PWD_outlier include options to initialize any 
  of rho, alpha, beta, or the vector mu.
* PWD_known is also updated to reflect the full likelihood, leading to an improved
  algorithm using the optim function.
* Various housekeeping and grammatical edits were also integrated:
* Adjusted various functions to better handle missing values.
* Expand the reference to our paper by adding its title wherever we just have the DOI.
* Under PWD_get_gh Details, removed "reduces the number of optimized parameters..." 
* Under PWD_inference Description, deleted "Currently".
* Under WD_Linnet Details, added "since the underlying precision profile models are not the same."
* Removed double-spacing in message-outputs (artifact of retaining "\n" from cat, no longer needed)
* Inside PWD_outlier function, removed extraneous text output and double-spacing.
* Under PWD_outlier values, edited the description of `scalr`.
* Clarification in details for PW_outlier.
* Clarification in the description of PWD_resi.
* Streamlined the package Description.
* Standardized the argument definitions using the parameter notation (versus words) and 
  ended each with a period (as is standard) instead of a comma.
* Revised the Details section for PWD_get_gh to reflect the fewer number of parameters.
* Updated the Reference for the primary reference paper with the new doi and accepted status.
* The main change was updates to the PWD_resi function: removed Pvals, and edited 
  documentation to match, reflecting removal from the submitted paper.
* Minor documentation edits, standardizing the annotation. 
* Adjusted output name to be L instead of like for -2 log Likelihood.
* Added details for default initial values of parameters.
* Set printem defaults to reduce automatic output for PWD_outlier.
* Minor adjustments to documentation language.  
* Added Argument notes and Details for PWD_get_gh documentation,
  defining rho as the ratio sigma/kappa.  
* Fixed output formatting for `PWD_outlier`.

# ppwdeming 1.0.6

* revised individual function documentation to reference arXiv document by
arXiv DOI in format <doi:10.48550/arXiv.YYMM.NNNNN>.
* adjusted one other doi reference correspondingly.

# ppwdeming 1.0.5

* updated algorithm in function `PWD_RL` to better handle rare non-convergence.

# ppwdeming 1.0.4

* revised primary reference in DESCRIPTION file to include author names and to
  have consistent format.
* For runnable (but >10s running) code in examples, replaced \dontrun with 
  \donttest.  This was applied to functions: PWD_inference, PWD_outlier.
* Removed the default console output by: 
  (1) setting the printem to default to FALSE; and
  (2) replacing any print()/cat() with message().
  This was applied for: PWD_get_gh, PWD_inference, PWD_known, PWD_outlier, 
    PWD_resi, and WD_Linnet. 
    Relevant items are elements on the list-value of the function.
* Replaced print() with message() and stop() for possible error in WD_General.

# ppwdeming 1.0.3

* placed examples directly into function-building files;
* used \dontrun{} to keep runtime under 10s.

# ppwdeming 1.0.2

* Corrected Github address.

# ppwdeming 1.0.1

* Updated arXiv URL and Description to more closely meet auto-check.

# ppwdeming 1.0.0

* Initial CRAN submission.

# ppwdeming 0.99.0

* Minor updates to documentation.

# ppwdeming 0.0.0.9012

* Most recently updated simplified examples
