"""
Module to manage making point datasets for ptobler PROJECTIONS work
"""

import random
import math
import copy
import os

import numpy as np
import pandas


def random_point():
    """
    Returns a lon lat point evenly distributed on a sphere

    Returns
    -------
            a 2-tuple (lon, lat) in degrees
    """
    lon = 2 * math.pi * random.random() - math.pi
    lat = math.acos(2 * random.random() - 1) - math.pi/2
    return math.degrees(lon), math.degrees(lat)


def random_points(npoints=1000):
    """
    Returns a list of lon lat points evenly distributed on a sphere

    Arguments
    ---------
            npoints: integer number of points required (default 1000)
    Returns
    -------
            a list of 2-tuples (lon, lat) in degrees
    """
    return [random_point() for _ in range(npoints)]


def random_points_df(npoints=1000):
    """
    Returns a pandas.DataFrame of lon lat points evenly distributed on a sphere

    Arguments
    ---------
            npoints: integer number of points required (default 1000)
    Returns
    -------
            a pandas.DataFrame with npoints (lon, lat) columns
    """
    pts = random_points(npoints)
    points = pandas.DataFrame({'lon': [p[0] for p in pts], 'lat': [p[1] for p in pts]})
    return points


def tag_with_ID_and_dir(df, direction='.'):
    """
    Part of making offsets to existing points in a pandas.DataFrame. Adds dir and ID attributes to a pandas.DataFrame.

    Arguments
    ---------
            df: a pandas.DataFrame
            direction: a string for the dir attribute (default '.')
    Returns
    -------
            The tagged pandas.DataFrame
    """
    df['dir'] = direction
    df['ID'] = list(range(df.shape[0]))
    return df


def add_offset_points(df, offset=(0, 0), direction='.'):
    """
    Returns an offset copy of the supplied pandas.DataFrame

    The copy is displaced by the specified (x, y) offset. Offset points retain
    any attributes of the original, but the 'dir' attribute is set to the
    supplied direction value.

    Arguments
    ---------
            df: a pandas.DataFrame
            offset: a 2-tuple of the x, y translation (default (0, 0))
            direction: a string for the dir attribute, e.g.if offset is (0,1)
                'N' would be an appropriate value (default '.')
    Returns
    -------
           The tagged pandas.DataFrame
    """
    df2 = copy.copy(df)
    df2.lon += offset[0]
    df2.lat += offset[1]
    df2.dir = direction
    return df2


def make_and_combine_offsets(df, d=0.001):
    """
    Returns a new pandas.DataFrame with duplicate 4x offset points

    The returned DataFrame contains the original points along with four
    N/E/S/W offsets by a distance d. Generally the points are expected to be in
    lon lat format, but any units are fine. New points will have the same ID as
    their parent points, and the dir tag will be a character from NESW as appropriate.

    Arguments
    ---------
            df: a pandas.DataFrame, which should have ID and dir attributes
            d: the distance offset to apply in each of the four directions (default 0.001)
    Returns
    -------
            A pandas.DataFrame containing the original points and four
                offset copies
    """
    # make a bunch of offset tuples left a bit, right a bit, down a bit, up a bit
    offsets = [(-d, 0), (d, 0), (0, -d), (0, d)]
    offset_copies = [add_offset_points(df, offset=o, direction=s) for o, s in zip(offsets, ['W', 'E', 'S', 'N'])]
    df2 = copy.copy(df)
    # merge in the offset copies
    for o in offset_copies:
        df2 = pandas.concat([df2, o], sort=True, ignore_index=True)
    return df2


def add_points(df=None, lons=(-180, -180, 180, 180), lats=(-90, 90, 90, -90)):
    """
    Adds the requested set of points to the supplied pandas.DataFrame (which can be None)

    Added points are inserted at the beginning of the DataFrame for more convenient
    inspection, but their ID numbers follow those in the supplied DataFrame. If df is
    None, only the specified points will be in the returned DataFrame.

    Arguments
    ---------
            df: a DataFrame to add points to
            lons: list of the longitudes of points to add (default [-180, -180, 180, 180])
            lats: list of the latitudes of points to add (default [-90, 90, 90, -90])
    Returns
            a pandas.DataFrame with ID, dir, lon, lat
    """
    pts = pandas.DataFrame({'lon': lons, 'lat': lats})
    if df is None:
        df = pts
        df = tag_with_ID_and_dir(df)
    else:
        if 'ID' in df.columns:
            n1 = max(df.ID)
            n2 = pts.shape[0]
            pts['ID'] = [x for x in range(n1 + 1, n1 + n2 + 1)]
            pts['dir'] = '.'
            df = pandas.concat([pts, df], sort=True, ignore_index=True)
        else:
            df = pandas.concat([pts, df], sort=True, ignore_index=True)
            df = tag_with_ID_and_dir(df)
    return df[['ID', 'dir', 'lon', 'lat']]


def make_random_points(df=None, num_points=1000):
    """
    Returns a point dataset as a pandas.DataFrame

    Arguments
    ---------
            df: a base dataset to augment, default None, when a random point
                dataset will be generated
            num_points: number of random points to generate if needed (default 1000)
    Returns
    -------
            a pandas.DataFrame, with ID, dir, lon, lat attributes
    """
    if df is None:
        df = random_points_df(npoints=num_points)
    else:
        df = pandas.concat([df, random_points_df(npoints=num_points)], ignore_index=True, sort=True)

    # tag the input data with ID and shift attributes
    df = tag_with_ID_and_dir(df)
    return df[['ID', 'dir', 'lon', 'lat']]


def make_lattice(df=None, lon_min=-180, lon_max=180, lon_num=25, lat_min=-90, lat_max=90, lat_num=13):
    """
    Makes a regularly spaced lattice of lon lat coordinates as a pandas.DataFrame

    Returned DataFrame has lon lat attributes, and optionally includes points
    from a provided GeoDataFrame of points (this is redundant, and it is
    recommended to use x = make_lattice, followed by add_points(x ...) to do this instead

    Arguments
    ---------
            df: a pandas.DataFrame of lon lat Points to add lattice
                points to (default None)
            lon_min: minimum longitude (inclusive) (default -180)
            lon_max: maximum longitude (inclusive) (default 180)
            lon_num: number of meridians (default 25)
            lat_min: minimum latitude (inclusive) (default -90)
            lat_max: maximum longitude (inclusive) (default 90)
            lat_num: number of parallels (default 13)
    Returns
    -------
            a pandas.DataFrame with ID dir lon lat
    """
    mypoints = []
    for lon in np.linspace(lon_min, lon_max, int(lon_num)):
        for lat in np.linspace(lat_min, lat_max, int(lat_num)):
            mypoints.append([lon, lat])

    gratpoints_df = pandas.DataFrame({'lon': [ll[0] for ll in mypoints],
                                      'lat': [ll[1] for ll in mypoints]})
    if df is not None:
        gratpoints_df = pandas.concat([df, gratpoints_df], sort=True, ignore_index=True)

    return tag_with_ID_and_dir(gratpoints_df)


def get_file_name(basename=None, grid_reqd=False, npoints=1000,
                  offsets_reqd=False, bbox_reqd=False, idx=0):
    """
    Returns a filesname suitable for subsequent processing

    Generally called by the make_point_datasets function

    Arguments
    ---------
            basename: default None
            grid_reqd: default False
            npoints: default 1000
            offsets_reqd: default False
            bbox_reqd: defaule False
            idx: default 0
    Returns
    -------
            filename string <rnd|lll>-n-<no-offsets|offsets><|bb>-idx
            examples:
                lll-329-offsets-bb-00.csv
                rnd-1000-no-offsets-00.csv
                ...
    """
    if basename is not None:
        filename = basename + '-'
    else:
        filename = "lll-" if grid_reqd else "rnd-"
    filename += str(npoints)
    filename += '-offsets' if offsets_reqd else "-no-offsets"
    filename += "-bb" if bbox_reqd else ""
    return filename + '-' + ('0' + str(idx))[-2:]


def make_point_datasets(folder, basename=None, nfiles=1, num_points=1000,
                        grid_reqd=False, grid_spec=(-180, 180, 25, -90, 90, 13),
                        offsets_reqd=False, offset_d=0.001,
                        bbox_reqd=False, bbbounds=(-180, -90, -180, 90, 180, 90, 180, -90),
                        write=True):
    """
    Makes point datasets according to the requested options

    Setting write False will return an example, setting it True will write the requested
    files and return a list of the filenames written (including the folder)

    Arguments
    ---------
            folder: (Required) folder name for the output files
            basename: a basename if desired, default None will produce auto-generated name(s)
            nfiles: number of files to write to this specification, indexed 00, 01, 02... etc
            num_points: number of points in the random version (default 1000)
            grid_reqd: True will produce a regular lattice of points according to grid_spec
            grid_spec: a list (lon_min, lon_max, lon_n, lat_min, lat_max, lat_n) with min and max
                in each direction, and tne number of meridians/parallels required. Default of
                (-180, 180, 25, -90, 90, 13) is a 15 degree spaced lattice
            offsets_reqd: True will produce the original points and 4 offset copies, each with
                the same ID as the original points, and a code NESW indicating the direction
                of offset (default False)
            offset_d: distance (in angular degrees) to offset (default 0.001)
            bbox_reqd: True will add points according to bbbounds, these will be offset to, if
                offset is requested (default False)
            bbbounds: points to add as a bounding box as (lon1, lat1, lon2, lat2, ...) Default
                (-180, -90, -180, 90, 180, 90, 180, -90). Note that any number of points can
                be 'injected' in this way, for example, the outline of a particular place might
                be requested.
            write: if True, files are written to folder, and a list of their names is returned,
                if False no file is written and an example is returned for inspection
    Returns
    -------
            if write True then a list of written filenames, if False an example pandas.DataFrame
    """
    if not os.path.exists(folder):
        os.mkdir(folder)

    nfiles = nfiles if write else 1
    filenames = []
    for i in range(nfiles):
        if bbox_reqd:
            lons = [bbbounds[i] for i in range(0, len(bbbounds), 2)]
            lats = [bbbounds[i] for i in range(1, len(bbbounds), 2)]
            gdf = add_points(None, lons, lats)
        else:
            gdf = None

        if grid_reqd:
            gdf = make_lattice(gdf, *grid_spec)
        else:
            gdf = make_random_points(gdf, num_points=num_points)

        if offsets_reqd:
            gdf = make_and_combine_offsets(gdf, d=offset_d)

        n = max(gdf.ID) + 1
        fname = folder + '/' + get_file_name(basename, grid_reqd, n, offsets_reqd, bbox_reqd, i) + '.csv'
        if write:
            gdf.to_csv(fname, index=False)
            filenames.append(fname)
    if write:
        return filenames
    else:
        return gdf


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Make random lon-lat points data')

    parser.add_argument('-n', type=int,
                        default=500, dest='NUM_POINTS',
                        help='the number of points required (default 500)')
    parser.add_argument('--nfiles', type=int,
                        default=1, dest='NUM_FILES',
                        help='the number of files to make (default 1)')
    parser.add_argument('--dest', type=str,
                        default='temp_data/input', dest='OUTPUT_FOLDER',
                        help='the output folder, which must already exist (default temp_data/input')
    parser.add_argument('--fname', type=str,
                        default='', dest='BASENAME',
                        help="basename for files, if not set defaults to 'rnd' or 'lll'")
    parser.add_argument('-o', dest='OFFSETS_REQUIRED', action='store_true',
                        help='offset points required, if dist other than 0.001 is required specify with --dist')
    parser.add_argument('--dist', type=float, dest='OFFSET_DIST',
                        default=0.001,
                        help='offset distance in degrees (default 0.001)')
    parser.add_argument('--bb', dest='BBOX_REQUIRED', action='store_true',
                        help="indicates that four 'corners' of the domain are required (defaults to world)")
    parser.add_argument('--bbox', type=float, nargs=8, dest='BBOX',
                        metavar=('lon1', 'lat1', 'lon2', 'lat2', 'lon3', 'lat3', 'lon4', 'lat4'),
                        default=[-180, -90, -180, 90, 180, 90, 180, -90],
                        help='list of 8 values specifying 4 lon lat pairs')
    parser.add_argument('--ll', dest='REG_GRID_REQUIRED', action='store_true',
                        help='indicates that a regular grid lattice is required')
    parser.add_argument('--ll-spacing', type=float, nargs=6,
                        metavar=('lon_min', 'lon_max', 'num_meridians', 'lat_min', 'lat_max', 'num_parallels'),
                        dest='LATTICE_SPEC',
                        default=[-180, 180, 25, -90, 90, 13],
                        help='list of 8 values specifying 4 lon lat pairs')

    args = parser.parse_args()

    make_point_datasets(folder=args.OUTPUT_FOLDER,
                        basename=args.BASENAME,
                        nfiles=args.NUM_FILES,
                        num_points=args.NUM_POINTS,
                        grid_reqd=args.REG_GRID_REQUIRED,
                        grid_spec=args.LATTICE_SPEC,
                        offsets_reqd=args.OFFSETS_REQUIRED,
                        offset_d=args.OFFSET_DIST,
                        bbox_reqd=args.BBOX_REQUIRED,
                        bbbounds=args.BBOX)
