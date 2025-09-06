# ppwdeming 1.0.6.1

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
