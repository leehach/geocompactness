# geocompactness

Python code for calculating shape compactness of geographic features.

Recommended usage:

```python
import compactness_measures as cm
```

Currently implements:

* `polsby-popper(geo)`, an isoperimetric quotient based on Polsby-Popper (1991) - "The Third Criterion: Compactness as a Procedural Safeguard against Partisan Gerrymandering".
* `schwartzberg(geo, inverse = True)`, an isoperimetric quotient based on Schwartzberg (1966) - "Reapportionment, gerrymanders, and the notion of compactness". Note that the original Schwartzberg paper proposed a measure that varied from 1 (most compact) to infinity (least compact). Today most analysts use the inverse which varies from 0 (least compact) to 1 (most compact), to make it comparable with other 0 to 1 measures. The function returns the inverse measure by default.
* `reock(geo)`, ratio of shape area to the area of the minumum bounding circle based on Reock (1961) - "A Note: Measuring Compactness as a Requirement of Legislative Apportionment". 
* `c_hull_ratio`, ratio of shape area to the area of the convex hull.
* `moment_of_inertia`, ratio of the MI (moment of inertia about the centroid) of a circle of equal area, to the MI of the shape. Moment of Inertia Shape Index (Li et al. 2013 - "An efficient measure of compactness for two-dimensional shapes and its application in regionalization problems")

There is an extensive literature on geographic compactness measures, in general and in its specific application to redistricting. Useful surveys include:

* MacEachren (1985) - "Compactness of Geographic Shape: Comparison and Evaluation of Measures"
* Young (1988) - "Measuring the compactness of legislative districts"
* Niemi et al. (1990) - "Measuring compactness and the role of a compactness standard in a test for partisan and racial gerrymandering"
* Horn et al. (1993) - "Practical application of district compactness"

This literature is somewhat fractured, and the same measures or very similar measures have often been rediscovered, and may be referred to by different names in different disciplines. The surveys above represent work by geographers (MacEachran 1985; Horn et al. 1993), legal scholars (Young 1988), and political scientists (Niemi et al. 1990).
