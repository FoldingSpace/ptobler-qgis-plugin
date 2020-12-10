import scipy
import scipy.interpolate
import scipy.spatial
import pandas
import geopandas
import shapely
import shapely.geometry
import numpy

import ptobler.inversedomain as ptobid
import ptobler.pointdata as ptobpd


SPHERE = '+proj=longlat +R=6371007 +no_defs'
PROJECTIONS = {
    'briesemeister': '+proj=ob_tran +o_proj=hammer +o_lat_p=45 +o_lon_p=-10 +lon_0=0 +R=6371007 +no_defs',
    'mollweide': '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +units=m +R=6371007 +no_defs',
}


# def empirical_projection_calculate_domain(df):
#     """
#     Returns a tuple of
#
#     Arguments
#     ---------
#             df: a DataFrame with lon, lat
#     Returns
#     -------
#
#     """
#     projpoints_longlat = geopandas.GeoDataFrame(df, geometry= geopandas.points_from_xy(df.lon, df.lat))
#     domain_longlat = geopandas.GeoSeries(
#         shapely.geometry.MultiPoint(geopandas.points_from_xy(df.lon, df.lat)).convex_hull)
#     return projpoints_longlat, domain_longlat


def empirically_project(proj_df, points_to_interpolate, discontinuities=None, estimate_cut=True,
                        method='naive', wrap=None, verbose=False):
    """
    The top level projection function - delegates the work to specific functions for different methods

    Only 'linear' and 'naive' methods currently implemented

    Arguments
    ---------
            proj_df: a dataframe with the lon, lat, x, y columns
            points_to_interpolate: a DataFrame with lon, lat atrributes of locations where
                    estimated x, y is required
            discontinuities: a geopandas.GeoDataFrame Polygon dataset where interpolation is
                    unsafe. If None (the default) uses ptobler.inversedomain.estimate_inverse_cuts
            estimate_cut: if True and no discontinuities provided one will be estimated, if
                    False no cut will be applied, default True
            method: string indicating method to use. Currently supports
                    'naive' (the default) for linear interpolation without discontinuities;
                    'linear' for linear interpolation with discontinuities;
            wrap: float amount in degrees to wrap the empirical projection at poles
                    and dateline (default None)
    Returns
    -------
            a DataFrame with lon, lat, x, y columns
    """
    if method == 'linear':
        print('Linear interpolation starting...')
        if estimate_cut and discontinuities is None:
            print('Estimating cut regions...')
            discontinuities = ptobid.estimate_inverse_cuts(proj_df)
        proj_df = add_wrap(proj_df, wrap=wrap)
        return linear_interpolation(proj_df[['lon', 'lat']],
                                    proj_df[['x', 'y']],
                                    discontinuities,
                                    points_to_interpolate,
                                    verbose=verbose)

    elif method == 'naive':
        return naive_interpolation(proj_df, points_to_interpolate)

    else:
        return None


def naive_interpolation(proj_df, points_to_interpolate):
    """
    Performs a naive linear interpolation (no allowance for interruptions)

    Arguments
    ---------
            proj_df: dataframe with lon, lat, x, y
            points_to_interpolate: a dataframe with lon, lat attributes where an estimated
                x, y is required
    Returns
    -------
            the points to interpolated DataFrame populated with interpolated x, y values
    """
    proj_x = scipy.interpolate.griddata(
        proj_df[['lon', 'lat']],  # points
        proj_df['x'],  # values
        points_to_interpolate[['lon', 'lat']],  # (grid_x, grid_y)
        method='linear')
    proj_y = scipy.interpolate.griddata(
        proj_df[['lon', 'lat']],  # points
        proj_df['y'],  # values
        points_to_interpolate[['lon', 'lat']],  # (grid_x, grid_y)
        method='linear')
    points_to_interpolate['x'] = proj_x
    points_to_interpolate['y'] = proj_y
    if len(numpy.bincount(numpy.isnan(proj_x))) > 1:
        print("Warning: Projections returned NaN.")
    return points_to_interpolate


def add_wrap(df, wrap=None):
    if wrap is None:
        return df
    # # north
    # subset = df[df.lat > (90 - wrap)][:]
    # subset.lat = 180 - subset.lat
    # subset.lon = subset.lon % 360 - 180
    # if 'ID' in subset.columns:
    #     subset.ID += subset.ID + subset.shape[0]
    # df = pandas.concat([df, subset], sort=True, ignore_index=False)
    #
    # # south
    # subset = df[df.lat < (-90 + wrap)][:]
    # subset.lat = -180 - subset.lat
    # subset.lon = subset.lon % 360 - 180
    # if 'ID' in subset.columns:
    #     subset.ID += subset.ID + subset.shape[0]
    # df = pandas.concat([df, subset], sort=True, ignore_index=False)

    # east
    subset = df[df.lon > (180 - wrap)][:]
    subset.lon = subset.lon - 360
    if 'ID' in subset.columns:
        subset.ID += subset.ID + subset.shape[0]
    df = pandas.concat([df, subset], sort=True, ignore_index=False)

    # west
    subset = df[df.lon < (-180 + wrap)][:]
    subset.lon = subset.lon + 360
    if 'ID' in subset.columns:
        subset.ID += subset.ID + subset.shape[0]
        
    return pandas.concat([df, subset], sort=True, ignore_index=False)


def get_barycentric_coordinates(transform_matrix, point_coords):
    """
    NO CLUE WHAT THIS DOES!!!

    Arguments
    ---------
            transform_matrix: ??
            point_coords: the coordinates to be transformed
    Returns
    -------
            numpy array of ??
    """
    b = transform_matrix[:2].dot(point_coords - transform_matrix[2])
    return numpy.append(b, 1-b.sum(axis=0))


def linear_interpolation(control_pt_ll, control_pt_xy,
                         interp_discontinuities=None, pts_to_interpolate_ll=None, verbose=False):
    """
    Performs linear interpolation in a triangulation of the supplied control points

    Arguments
    ---------
            control_pt_ll: pandas.DataFrame with lon, lat coords for which x, y values are known
            control_pt_xy: x, y coordinates in the projected space of the lon, lats
            interp_discontinuities: a Polygon geopandas.GeoDataFrame delineating discontinuity regions
                in the projection. Defaults to an empty geopandas.GeoDataFrame if None is supplied
            pts_to_interpolate_ll: lon, lat coordinates at which projected coordinates are required. Defaults
                to a coarse lon lat grid at 15 degree spacing if None is supplied
    :return:
    """
    print('Performing interpolation...')
    # if interp_discontinuities is None:
    #     interp_discontinuities = geopandas.GeoDataFrame(geometry=[], crs=SPHERE)

    if pts_to_interpolate_ll is None:
        pts_to_interpolate_ll = ptobpd.make_lattice(lon_num=25, lat_num=13)

    pts_to_interp = pts_to_interpolate_ll[['lon', 'lat']].to_numpy(copy=True)

    points_array = control_pt_ll[['lon', 'lat']].to_numpy(copy=True)
    tri = scipy.spatial.Delaunay(points_array)

    delaunay_triangles_gdf = geopandas.GeoDataFrame(
        geometry=[shapely.geometry.Polygon(points_array[s]) for s in tri.simplices],
        crs=SPHERE)

    if interp_discontinuities is not None:
        simplexes_intersecting_discontinuities_gdf = geopandas.sjoin(
            delaunay_triangles_gdf,
            interp_discontinuities,
            op='intersects')

        simplexes_intersecting_discontinuities_gdf.rename(
            columns={'index_right': 'intersecting_domain_edge_segment'},
            inplace=True)

        bad_simplexes = simplexes_intersecting_discontinuities_gdf['geometry'] #.drop_duplicates()
        bad_simplex_indices = bad_simplexes.index

    triangles_found = tri.find_simplex(pts_to_interp)

    results_array = []
    bad_points = []
    for p_num in range(pts_to_interp.shape[0]):
        if (triangles_found[p_num].size == 1) and (triangles_found[p_num] < 0):
            results_array.append([numpy.NaN, numpy.NaN])
            if verbose:
                print("Triangle not found for: ", pts_to_interp[p_num])
        else:
            # for one to many transformations, this needs to be modified... perhaps it is an array?
            curr_triangle = triangles_found[p_num]
            p = pts_to_interp[p_num]

            if interp_discontinuities is not None and curr_triangle in bad_simplex_indices:
                bad_points.append(p.tolist())
                results_array.append([numpy.NaN, numpy.NaN])
            else:
                bc = get_barycentric_coordinates(tri.transform[curr_triangle], p)
                curr_triangle_values = control_pt_xy.iloc[tri.simplices[curr_triangle]][{'x', 'y'}]
                interpolated_x = curr_triangle_values['x'].to_numpy(copy=True).dot(bc)
                interpolated_y = curr_triangle_values['y'].to_numpy(copy=True).dot(bc)
                results_array.append([interpolated_x, interpolated_y])

    if len(bad_points) > 0:
        if verbose:
            print("Unable to interpolate points because in bad simplex: ", bad_points)
        else:
            print("Unable to interpolate: %d" % (len(bad_points)))

    pts_to_interpolate_ll['x'] = [xy[0] for xy in results_array]
    pts_to_interpolate_ll['y'] = [xy[1] for xy in results_array]
    return pts_to_interpolate_ll


def get_xy_from_pts(pts):
    """
    Utility function, returns x, y values as two lists from a geopandas.GeoSeries of Points

    Arguments
    ---------
            pts: a Point Geoseries
    Returns
    -------
            a 2-tuple of lists of the x and y coordinates
    """
    return [p.coords[0][0] for p in pts], [p.coords[0][1] for p in pts]


def make_points(pts):
    """
    Returns a geopandas.GeoSeries of Points from a list of 2-tuples of x, y coordinates

    Arguments
    ---------
            pts: a list of 2-tuples of (x, y) coordinates
    Returns
    -------
            a Point GeoSeries
    """
    return geopandas.GeoSeries(shapely.geometry.Point(p[0], p[1]) for p in pts)


def get_proj4string(shortform=None):
    """

    :param shortform:
    :return:
    """
    if shortform is None:
        return PROJECTIONS.keys()
    else:
        return PROJECTIONS[shortform]


def proj4_project(f, folder, proj, write=False):
    """
    Arguments
    ---------
            proj: a tuple (shortform name, proj4 string)
    """
    basename = f.rsplit('/', 1)[-1].split('.')[0]
    pts = pandas.read_csv(f)
    pts = geopandas.GeoDataFrame(pts, geometry=make_points(zip(pts.lon, pts.lat)))
    pts.crs = SPHERE
    pts = pts.to_crs(proj[1])
    pts['x'], pts['y'] = get_xy_from_pts(pts.geometry.envelope)
    ptsdf = pts.drop(['geometry'], axis=1)
    if write:
        fname = folder + '/' + basename + '-p4-' + proj[0] + '.csv'
        ptsdf.to_csv(fname, index=False)
        return fname
    else:
        return ptsdf

if __name__ == '__main__':

    ## Code attempting to test the intersection of the delaunay triangulation
    ## with a cuts geodataframe causing apparent problems
    emp_proj_df = geopandas.read_file(
        '../temp_data/distancebearing/empirical_projections/dgg-7292-no-offsets-handmade-distancebearing.csv'
    )
    to_interp = geopandas.read_file(
        '../temp_data/unprojected_data/world_better.csv'
    )

    clon = -119.70
    clat = 34.416667
    cut_polys = geopandas.GeoSeries([shapely.geometry.Point(clon, clat),
                                     shapely.geometry.Point(clon + 360, -clat)]).buffer(1).unary_union

    cuts = geopandas.GeoDataFrame(geometry=geopandas.GeoSeries([cut_polys]))
    cuts.crs = SPHERE

    linear_interpolation(emp_proj_df[['lon', 'lat']],
                         emp_proj_df[['x', 'y']],
                         interp_discontinuities=cuts,
                         pts_to_interpolate_ll=to_interp[['lon', 'lat']])




def make_1x3_gdf(gdf):
    """
    Makes a 1 x 3 duplicate of the input (assumed whole world, lon, lat) input data

    The coordinates extended beyond -180/+180 lon and -90/+90 lat are transformed
    appropriately such that they are 'mirrored' in the dateline and/or poles. While
    lon, lat values are changed the x, y coordinates accompany the transformed lon,
    lat values unchanged.

    Arguments
    ---------
            gdf: a geopandas.GeoDataFrame to extend, with lon, lat coordinates
    Returns
    -------
            an extended geopandas.GeoDataFrame
    """
    gdf_1x3 = pandas.concat([
        geopandas.GeoDataFrame(gdf.copy(), geometry=geopandas.GeoSeries.translate(gdf.copy(), xoff=-360.0)),
        gdf.copy(),
        geopandas.GeoDataFrame(gdf.copy(),geometry=geopandas.GeoSeries.translate(gdf.copy(), xoff=360.0))
    ], ignore_index=True)

    # gdf_3x3 = pandas.concat(
    #     [
    #         geopandas.GeoDataFrame(
    #             gdf_row.copy(),
    #             geometry=geopandas.GeoSeries.translate(
    #                 geopandas.GeoSeries.scale(
    #                     gdf_row.copy(),
    #                     xfact=1,yfact=-1,origin=(0,0)),
    #                 xoff=180.0,
    #                 yoff=-180.0)
    #         ),
    #         gdf_row.copy(),
    #         geopandas.GeoDataFrame(
    #             gdf_row.copy(),
    #             geometry=geopandas.GeoSeries.translate(
    #                 geopandas.GeoSeries.scale(
    #                     gdf_row.copy(),
    #                     xfact=1,yfact=-1,origin=(0,0)),
    #                 xoff=180.0,
    #                 yoff=180.0)
    #         )
    #     ],
    #     ignore_index=True
    # )
    gdf_1x3.crs = gdf.crs
    # gdf_1x3['lon'] = gdf_1x3.geometry.apply(lambda x: x.coords[0][0])
    # gdf_1x3['lat'] = gdf_1x3.geometry.apply(lambda x: x.coords[0][1])
    return gdf_1x3


def make_1x3_clipped_gdf(gdf, lon_min=-200, lon_max=200):
    """
    Clips the supplied geopandas.GeoDataFrame to the supplied extents

    Arguments
    ---------
            gdf: the geopandas.GeoDataFrame to clip
            lon_min: minimum longitude (default -200)
            lon_max: maximum longitude (default 200)
    Returns
    -------
            a clipped geopandas.GeoDataFrame
    """
    bbox = shapely.geometry.Polygon([(lon_min, -90), (lon_min, 90), (lon_max, 90), (lon_max, -90)])
    bbox = geopandas.GeoDataFrame(geometry=geopandas.GeoSeries([bbox]))
    bbox.crs = gdf.crs
    return geopandas.overlay(make_1x3_gdf(gdf), bbox, how='intersection')
    # return gdf[
    #     (gdf['lon'] > lon_min) &
    #     (gdf['lon'] < lon_max)
    # ]


