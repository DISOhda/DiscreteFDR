# DiscreteFDR 2.0.1

-   Introduction of `mode` parameter for `hist()` function to adapt construction
    of histograms in case of conditional p-value selection.
-   Remove `amnesia` dataset (moved to `DiscreteDatasets` package).
-   Function `match.pvals()` is no longer exported.
-   Performance improvement for step-up procedures, especially for large numbers
    of tests.


# DiscreteFDR 2.0.0

-   New features:
    -   `discrete.BH()`, `DBH()`, `ADBH()` and `DBR()` are now generic
        functions. The previously existing functionality is implemented in
        `*.default` methods.
    -   `discrete.BH()`, `DBH()`, `ADBH()` and `DBR()` got 
        `*.DiscreteTestResults` methods for processing `DiscreteTestResults` R6
        class objects from package `*.DiscreteTests` directly, so they can be
        used within pipes.
    -   For consistency of new generics and methods, the first parameter
        `raw.pvalues` needed to be renamed to `test.results`.
    -   New parameter `threshold` for `discrete.BH()`, `DBH()`, `ADBH()` and
        `DBR()`. This enables selection of p-values which are smaller than or
        equal to a certain value. **Note**: observed p-values and their
        supports are then re-scaled, as the p-value distributions are now
        becoming conditional distributions. If no selection is performed (i.e.
        `threshold = 1`), `print()`, `summary()` and `plot()` outputs are as
        before. Otherwise, the now respect the re-scaled conditional
        distributions. Additionally, the `DiscreteFDR` S3 class output objects
        of the functions `discrete.BH()`, `DBH()`, `ADBH()` and `DBR()` now
        include a list `Select` with values and information regarding selection.
    -   New parameter `pCDFlist.indices` for `discrete.BH()`, `DBH()`, `ADBH()`
        and `DBR()`, which must have the same length as `pCDFlist` and may
        help increasing performance considerably. As `pCDFlist` can now include
        only unique supports, `pCDFlist.indices` must indicate the index of the
        p-values to which a given support belongs. If `pCDFlist` has the same
        length as `test.results`, it can be omitted (by setting it to `NULL`,
        the default). If users prefer using `DiscreteTestResults` objects, they
        do not have to take care of this, as unique supports and indices are
        automatically extracted from these objects.
-   New functions `generate.pvalues()` and `direct.discrete.BH()` as more
    flexible replacements for `fisher.pvalues.support()` and `fast.discrete()`.
-   Step function evaluation in C++ code has been replaced by closely optimized
    inline functions which offer performance gains of 10-50%.
    

# DiscreteFDR 1.3.7

-   Introduction of `lifecycle` mechanisms.
-   Marked `fast.Discrete()`, `fisher.pvalues.support()`, `match.pvals()`,
    `kernel_*()` and `amnesia` dataset as deprecated.
-   Various documentation updates.
-   Removal of links to `discreteMTP` packages, since it was removed from CRAN.


# DiscreteFDR 1.3.6

-   Fixed a problem with `fisher.pvalues.support` that could cause p-values to 
be wrong or NA (Thanks to Iqraa Meah).
-   Added GitHub.


# DiscreteFDR 1.3.5

-   Fixed a problem with `fisher.pvalues.support` that could cause an infinite
loop when using `alternative = two.sided` (Thanks to Lukas Jansen).
-   Changed version scheme from `x.y-z` to `x.y.z`


# DiscreteFDR 1.3-4

-   Added a `NEWS.md` file to track changes to the package.
-   Corrected a bug in `plot.DiscreteFDR` function that produced a false legend.
-   Added plausibility checks of arguments to `discrete.BH` and `DBR` functions.
