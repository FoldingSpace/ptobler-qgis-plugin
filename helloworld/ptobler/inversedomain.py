import numpy
import scipy.spatial

import geopandas
from shapely.geometry import Polygon
from pysal.lib.cg.alpha_shapes import alpha_shape_auto

import numba


def get_forward_projection(df):
    """
    Returns a dictionary of (lon, lat) -> (x, y)

    Dataframe is assumed to have x. y, lon, lat columns

    Arguments
    ---------
            df: a pandas.DataFrame with x, y, lon, lat columns
    Returns
    -------
            a dictionary with keys (lon, lat) values (x, y)
    """
    return dict(zip(zip(df.lon, df.lat), zip(df.x, df.y)))


def get_back_projection(df):
    """
    Returns a dictionary of (x,y) -> (lon, lat)

    Dataframe is assumed to have x. y, lon, lat columns

    Arguments
    ---------
            df: a pandas.DataFrame with x, y, lon, lat columns
    Returns
    -------
            a dictionary with keys (x, y) values (lon, lat)
    """
    return dict(zip(zip(df.x, df.y), zip(df.lon, df.lat)))


def get_convex_hull(df, xvar='x', yvar='y'):
    """
    Returns convex hull of specified attributes in df as a list of (x,y) tuples

    Note: use scipy.spatial.ConvexHull because shapely convex_hull appears unreliable
    see https://github.com/Toblerity/Shapely/issues/541 due to an upstream issue
    in GEOS

    Arguments
    ---------
            df: a pandas.DataFrame that has columns named in xvar, yvar
            xvar: name of the column containing eastings (default 'x')
            yvar: name of the column containing northings (default 'y')
    Returns
    -------
            a list of (x, y) tuples
    """
    points = numpy.array(df[[xvar, yvar]])
    hull = scipy.spatial.ConvexHull(points)
    return [tuple(points[i]) for i in hull.vertices]


def get_alpha_shape(df, xvar='x', yvar='y'):
    """
    Returns optimized alpha shape of specified attributes in df as a list of (x,y) tuples

    Arguments
    ---------
            df: a pandas.DataFrame that has columns named in xvar, yvar
            xvar: name of the column containing eastings (default 'x')
            yvar: name of the column containing northings (default 'y')
    Returns
    -------
            a list of (x, y) tuples
    """
    points = numpy.array([(xy[0], xy[1]) for xy in zip(df[xvar], df[yvar])])
    return [xy for xy in alpha_shape_auto(points).boundary.coords]


def back_project(points, Q):
    """
    Returns a list of tuples of points projected by lookup in dictionary Q

    Arguments
    ---------
            points: a list of (x, y) tuples
            Q: a dictionary with keys including the (x, y) tuples and values
                their corresponding (lon, lat) tuples
    Returns
    -------
            a list of (lon, lat) tuples
    """
    return [Q[xy] for xy in points]


def estimate_inverse_cuts(df, alpha_shape=False, method='CG', as_list=False):
    """
    Returns an estimate of the cuts or holes in the domain

    Currently only estimation by computational geometry (convex hull or alpha shape)
    is supported

    Arguments
    ---------
            df: pandas.DataFrame with x, y, lon, lat columns
            alpha_shape: bool, if True returns an alphashape, if False a convex hull
                (default True)
            method: a string indicate the method to use. 'CG' is the only supported
                method (computational geometry) at present (default 'CG')
            as_list: bool if True returns list of 2-tuples, False returns a Polygon
                (default False)
    Returns
    -------
            a geopandas.GeoDataFrame or list of 2-tuples
    """
    if method == 'CG':
        if alpha_shape:
            poly = back_project(get_alpha_shape(df), get_back_projection(df))
        else:
            poly = back_project(get_convex_hull(df), get_back_projection(df))
        if as_list:
            return poly
        else:
            poly = Polygon(poly)
            gs = geopandas.GeoSeries([poly], crs='+proj=longlat +R=6371007 +no_defs')
            return geopandas.GeoDataFrame(geometry=gs)

    return None


def estimate_inverse_domain(df, alpha_shape=False, method='CG'):
    """
    Returns an estimate of the domain of the projection

    Currently only estimation by computational geometry is supported

    Arguments
    ---------
            df: pandas.DataFrame with x, y, lon, lat columns
            alpha_shape: bool if True estimates holes in the domain as an
                alpha shape, if False uses convex hull (default False)
            method: string indicating method. Currently only computational
                geometry is supported (default 'CG')
    Returns
            a geopandas.GeoDataFrame (with holes)
    """
    if method == 'CG':
        cuts = estimate_inverse_cuts(df, alpha_shape, as_list=True)
        exterior = get_convex_hull(df, xvar='lon', yvar='lat')
        poly = Polygon(exterior, [cuts])
        gs = geopandas.GeoSeries([poly], crs='+proj=longlat +R=6371007 +no_defs')
        return geopandas.GeoDataFrame(geometry=gs)

    return None
