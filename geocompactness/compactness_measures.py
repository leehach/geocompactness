"""
compactness_measures: Functions to calculate compactness measures, and
components of compactness measures, in Euclidean space. Recommended usage

    import compactness_measures as cm

"""

import geopandas as gpd
import pandas as pd
from math import pi, sqrt
from geocompactness.smallest_enclosing_circle import make_circle

def _discrete_perimeter(geo, geo_cell):
    """Not implemented"""
    
    return None

def _continuous_perimeter(geo):
    """returns geo.length"""
    
    return geo.length

def _gross_perimeter(geo, geo_cell):
    """Not implemented"""
    
    return None

def _discrete_area(geo, geo_cell):
    """Not implemented"""
    
    return None

def _continuous_area(geo):
    """returns geo.area"""
    
    return geo.area

def perimeter(geo, geo_cell = None, gross = False):
    """
    Return perimeters of geometries in GeoSeries as Series of floats.
    
    Keyword arguments:
        geo -- GeoSeries or GeoDataFrame
        geo_cell -- GeoSeries or GeoDataFrame representing units used to build
            geo (the "container"); does not have to nest cleanly
        gross -- Calculate perimeter of geo by measuringstraight line distance 
            between nodes of geo_cell on boundary of geo (see Schwartzberg 1966)
            
    This function calculates continuous or discrete perimeter. 
    
    Continuous (Euclidean) perimeter is calculated if only geo argument is 
    provided. Currently this function just returns GeoSeries.length. 
    Future improvements could include:
        
        * Checking for lat-long coordinate system and performing geodetic
        measurement
        * Determining appropriate local CRS (most likely a State Plane or UTM
        zone) and performing calculation in that CRS.
        
    NOT YET OPERATIONALIZED: Discrete perimeter is calculated if a second
    geographic argument is provided that represents the "cells" or "building 
    blocks" of the first, larger geography.
    """

    if geo_cell is None:
        # Continuous perimeter
        return _continuous_perimeter(geo)
    elif gross:
        return _gross_perimeter(geo, geo_cell)
    else:
        return _discrete_perimeter(geo, geo_cell)

def area(geo, geo_cell = None, convex_hull = False):
    """
    Return areas of geometries in GeoSeries as Series of floats.
    
    Keyword arguments:
        geo -- GeoSeries or GeoDataFrame
        geo_cell -- GeoSeries or GeoDataFrame representing units used to build
            geo (the "container"); does not have to nest cleanly
        convex_hull -- Calculate area of convex hull of geo
        
    This function calculates continuous or area. 
    
    Continuous (Euclidean) area is calculated if only geo argument is 
    provided. Currently this function just returns GeoSeries.area. 
    Future improvements could include:
        
        * Checking for lat-long coordinate system and performing geodetic
        measurement
        * Determining appropriate local CRS (most likely a State Plane or UTM
        zone) and performing calculation in that CRS.
        
    NOT YET OPERATIONALIZED: Discrete area is calculated if a second
    geographic argument is provided that represents the "cells" or "building 
    blocks" of the first, larger geography.
    """

    if geo_cell is None:
        # Continuous area
        if convex_hull:
            return _continuous_area(geo.convex_hull)
        else:        
            return _continuous_area(geo)
    else:
        return _discrete_area(geo)
    
def polsby_popper(geo, geo_cell = None):
    """
    Returns Polsby-Popper (1991) compactness of geo as float
    
    Keyword arguments:
        geo -- GeoSeries or GeoDataFrame
        geo_cell -- GeoSeries or GeoDataFrame representing units used to build
            geo (the "container"); does not have to nest cleanly
    """
    
    return 4 * pi * area(geo, geo_cell) / (perimeter(geo, geo_cell) ** 2)

def schwartzberg(geo, inverse = True, geo_cell = None):
    """
    Returns Schwartzberg (1966) compactness of geo as float
    
    Keyword arguments:
        geo -- GeoSeries or GeoDataFrame
        inverse -- Boolean, return the original (1 to infinity) Schwartzberg,
            or the inverse (0 to 1); most analysts use the inverse, so this
            function defaults to inverse = True
        geo_cell -- GeoSeries or GeoDataFrame representing units used to build
            geo (the "container"); does not have to nest cleanly
    """
    
    schw = polsby_popper(geo, geo_cell) ** -0.5
    
    if inverse:
        return schw ** -1
    else:
        return schw

def c_hull_ratio(geo):
    
    return area(geo) / area(geo, convex_hull = True)

def reock(geo):
    """
    Returns Reock (1961) compactness of geo as float
    
    Keyword arguments:
        geo -- GeoSeries or GeoDataFrame
    """
    
    # make_circle returns x, y, r of circle. Use r to calculate area.
    mbc_area = geo.convex_hull.apply(lambda x: pi * make_circle(list(x.exterior.coords))[2] ** 2)
    return geo.area / mbc_area

def _moments(pts):
    
    if pts[0] != pts[-1]:
        pts = pts + pts[:1]
    x = [ c[0] for c in pts ]
    y = [ c[1] for c in pts ]
    sxx = syy = sxy = 0
    for i in range(len(pts) - 1):
        sxx += (y[i]**2 + y[i]*y[i+1] + y[i+1]**2) * (x[i]*y[i+1] - x[i+1]*y[i])
        syy += (x[i]**2 + x[i]*x[i+1] + x[i+1]**2) * (x[i]*y[i+1] - x[i+1]*y[i])
        sxy += (x[i]*y[i+1] + 2*x[i]*y[i] + 2*x[i+1]*y[i+1] + x[i+1]*y[i]) * (x[i]*y[i+1] - x[i+1]*y[i])

    a = coords_area(pts)
    print(a)
    
    cx, cy = coords_centroid(pts)
    print(cx, cy)
    
    return (sxx/12 - a*cy**2) + (syy/12 - a*cx**2)

    # return sxx, syy, sxy

def _moment_of_inertia(pts):
    """Returns the normalized (circle = 1) moment of inertia of a shape
    
    Keyword arguments:
        pts -- Iterable of tuples containing xy coordinate pairs of points
            that define a closed linear ring (i.e. a polygon boundary)
          
    Adapted from https://leancrew.com/all-this/2018/01/python-module-for-section-properties/
    """
    
    # Not required as we will always be passing closed linear rings
    #if pts[0] != pts[-1]:
    #    pts = pts + pts[:1]
    
    x = [ c[0] for c in pts ]
    y = [ c[1] for c in pts ]

    s = 0
    sx = sy = 0
    sxx = syy = 0 #sxy = 0

    for i in range(len(pts) - 1):
        s += x[i]*y[i+1] - x[i+1]*y[i]
        sx += (x[i] + x[i+1])*(x[i]*y[i+1] - x[i+1]*y[i])
        sy += (y[i] + y[i+1])*(x[i]*y[i+1] - x[i+1]*y[i])
        sxx += (y[i]**2 + y[i]*y[i+1] + y[i+1]**2)*(x[i]*y[i+1] - x[i+1]*y[i])
        syy += (x[i]**2 + x[i]*x[i+1] + x[i+1]**2)*(x[i]*y[i+1] - x[i+1]*y[i])
        # sxy += (x[i]*y[i+1] + 2*x[i]*y[i] + 2*x[i+1]*y[i+1] + x[i+1]*y[i])*(x[i]*y[i+1] - x[i+1]*y[i])

    a = s/2
    cx, cy = sx/(6*a), sy/(6*a)
    mi = (sxx/12 - a*cy**2) + (syy/12 - a*cx**2)  

    return abs(a**2 / (2 * pi * mi))

def _discrete_moment_of_inertia(geo, geo_cell = None, wt = 1):
    """Not implemented"""

    # Copy geo_cell so that it is not modified during function    
    tmp = geo_cell[[wt, geo_cell.geometry.name]].copy()
    
    # Calculate area of each cell, used later
    tmp["dA"] = geo_cell.area
    
    # Assign container ID using representative point, since centroid of small
    # geometry may fall outside the geometry, and therefore outside the container
    # for small geometries near edge of container
    tmp.geometry = geo_cell.representative_point()
    tmp.crs = geo_cell.crs # For some reason representative_point() destroys CRS, reassign it
    tmp = gpd.sjoin(tmp, geo[[geo.geometry.name]]) # After sjoin, order of GeoDataFrame changes. Index is preserved.

    # Replace representative point with centroid, which is needed for MI
    tmp.geometry = geo_cell.centroid

    # Calculate population-weighted centroid of each container geo
    # index_right refers to index of each feature in container
    cx = tmp.groupby("index_right").apply(lambda row: (row.geometry.x * row[wt]).sum() / row[wt].sum()).rename("cx")
    cy = tmp.groupby("index_right").apply(lambda row: (row.geometry.y * row[wt]).sum() / row[wt].sum()).rename("cy")
    tmp = tmp.merge(cx, on = "index_right").merge(cy, on = "index_right")
    
    # Calculate moment of inertia of each container geo
    tmp["centroid_distance"] = ((tmp.geometry.x - tmp["cx"])**2 + (tmp.geometry.y - tmp["cy"])**2).apply(sqrt)
    mi = tmp.groupby("index_right").apply(lambda row: (row.centroid_distance**2 * row[wt]).sum()).rename("mi")
    
    # Calculate radius of circle with same area as geo
    circle_mi = tmp.groupby("index_right").apply(lambda row: row[wt].sum() * row["dA"].sum() / (2 * pi))
    
    # Calculate compactness as ratio of MI of circle with same area and equal mass distribution to MI of container
    c_mi = circle_mi / mi
    c_mi.rename("c_mi", inplace = True)
      
    return c_mi
    
def moment_of_inertia(geo, geo_cell = None, wt = 1):
    """
    Returns moment of intertia shape index (MacEachren 1985; Li, et al. 2013) of geo as float
    
    Keyword arguments:
        geo -- GeoSeries or GeoDataFrame
        geo_cell -- GeoSeries or GeoDataFrame representing units used to build
            geo (the "container"); does not have to nest cleanly
        wt -- Int, float, or str: The weight of the unit (e.g. population
            count) in the moment of intertia calculation; this is either a
            number applied as the weight of every unit, or the name of a 
            column in geo_cell with weights; if you want the units weighted
            by area, do not supply values for geo_cell or wt
            
    The moment of intertia of a shape changes as the shape is translated in
    coordinate space. The moment of inertia about the centroid is the 
    minimum moment of intertia for a given shape. The MI of a circle of the
    same area as the shape is divided by the MI of the shape to convert to a 
    "shape index". As with other compactness measures, the final value varies
    from 0 (least compact) to 1 (most compact, a circle).
    """
    
    if geo_cell is None:
        moments = []    
        for geom in geo.geometry:
            
            moments.append(_moment_of_inertia(geom.exterior.coords))
    
        return pd.Series(moments)
    else:
        
        # wt is int or float, repeat value in "wt" column
        if type(wt) in [int, float, tuple, list, pd.Series]:
            geo_cell["wt"] = wt
            wt = "wt"

        return _discrete_moment_of_inertia(geo, geo_cell, wt)
    
        
