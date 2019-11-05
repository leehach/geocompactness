# geocompactness

Python code for calculating shape compactness of geographic features.

Currently implements:

* Polsby-Popper (1991) - "The Third Criterion: Compactness as a Procedural Safeguard against Partisan Gerrymandering"
* Schwartzberg (1966) - "Reapportionment, gerrymanders, and the notion of compactness"
* Reock (1961) - "A Note: Measuring Compactness as a Requirement of Legislative Apportionment"
* Convex Hull Ratio
* Moment of Inertia Shape Index (Li et al. 2013 - "An efficient measure of compactness for two-dimensional shapes and its application in regionalization problems")

There is an extensive literature on geographic compactness measures, in general and in its specific application to redistricting. Useful surveys include:

* MacEachren (1985) - "Compactness of Geographic Shape: Comparison and Evaluation of Measures"
* Young (1988) - "Measuring the compactness of legislative districts"
* Niemi et al. (1990) - "Measuring compactness and the role of a compactness standard in a test for partisan and racial gerrymandering"
* Horn et al. (1993) - "Practical application of district compactness"

This literature is somewhat fractured, and the same measures or very similar measures have often been rediscovered, and may be referred to by different names in different disciplines. The surveys above represent work by geographers (MacEachran 1985; Horn et al. 1993), legal scholars (Young 1988), and political scientists (Niemi et al. 1990).
