import numpy as np
import copy

import pandas


def get_empty_analysis_dataframe():
    return pandas.DataFrame(columns=[
        'full_filename',
        'is_empirical_proj',
        'data_filename', 'proj_filename',
        'data_num_points', 'proj_num_points',
        'data_spatial_distribution', 'proj_spatial_distribution',
        'proj_name',
        'h_mean', 'h_stddev',
        'k_mean', 'k_stddev',
        's_mean', 's_stddev',
        'omega_mean', 'omega_stddev'
    ], data=[])


def add_projected_coord_columns(pts):
    """
    Reshapes projected dataframe with offset points for analysis

    Mainly for internal use. Unpacks the x, y associated with the various
        E, W, N, S, offsets and makes xW, yW, xE, yE, etc columns that will
        used in the analysis

    Arguments
    ---------
            pts: a dataframe with lon, lat, x, y, ID and dir attributes, as
                follows:
                    ID` = unique ID for each 'parent' point, shared with
                            its 4 shifted 'child' points
                    dir` = a code where: . indicates the source point; and
                            W, E, S, N indicates displacement from parent
                            with same ID
                    lon, lat = the longitude and latitude values
                    x, y = the projected coordinates
    Returns
    -------
            reshaped dataframe with columns for lon, lat, x, y at each of the
                    offset directions (16 additional columns in all)
    """
    results = copy.copy(pts)
    results = results[results.dir == '.']
    for dirn in ['W', 'E', 'S', 'N']:
        displaced_coord = pts[pts.dir == dirn][['ID', 'lon', 'lat', 'x', 'y']]
        results = results.merge(displaced_coord, on='ID', suffixes=('', dirn))
    return results


def handle_discontinuities(row, verbosity=0, suppress=0):
    """
    Detects unusually rapid change in the partial differences of the supplied data

    Attempts to check to see if there is an interruption in the map projection
    that runs *through* a set of offsetted points.

    Assume that deviation of a factor of 100 in either of the following pairs of distances (absolute values)
    means an interruption:
    - between x-xW and xE-x, or
    - between yN-y and y-yS

    If an interruption is detected, turn the distortion measures for the row into NaNs and flag the row.

    Arguments
    ---------
            row: a single row of data from the analysis table
            verbosity: 0 will print minimal information for user, >0 will show too much information
            suppress: suppress any results where the error_flag < 1/suppress or > suppress, default 0
                    causes no suppression of results
    Returns
    -------
            a row with suspect results tagged as above
    """
    ZERO_TOL = 10 ** -7

    if (abs(np.float64(row.x) - row.xW) < ZERO_TOL) and (abs(row.xE - row.x) < ZERO_TOL):
        dxRatioEW = 1.0
    else:
        dxRatioEW = abs(row.xE - row.x) / abs(np.float64(row.x) - row.xW)

    if (abs(np.float64(row.x) - row.xS) < ZERO_TOL) and (abs(row.xN - row.x) < ZERO_TOL):
        dxRatioNS = 1.0
    else:
        dxRatioNS = abs(row.xN - row.x) / abs(np.float64(row.x) - row.xS)

    if (abs(np.float64(row.y) - row.yW) < ZERO_TOL) and (abs(row.yE - row.y) < ZERO_TOL):
        dyRatioEW = 1.0
    else:
        dyRatioEW = abs(row.yE - row.y) / abs(np.float64(row.y) - row.yW)

    if (abs(np.float64(row.y) - row.yS) < ZERO_TOL) and (abs(row.yN - row.y) < ZERO_TOL):
        dyRatioNS = 1.0
    else:
        dyRatioNS = abs(row.yN - row.y) / abs(np.float64(row.y) - row.yS)

    offset_ratios = np.array([dxRatioEW, dxRatioNS, dyRatioEW, dyRatioNS])
    row['error_flag'] = np.max(offset_ratios)  # 'projection_interruption_likely'

    if not suppress == 0:
        if any(np.isnan(offset_ratios)) or any(offset_ratios <= (1/ suppress)) or any(offset_ratios >= suppress):
            row['cos_lat'] = float('NaN')
            row['sin_theta'] = float('NaN')
            row['a_dash'] = float('NaN')
            row['b_dash'] = float('NaN')
            row['h'] = float('NaN')
            row['k'] = float('NaN')
            row['s'] = float('NaN')
            row['omega'] = float('NaN')
            if verbosity > 0:
                print("Discountinuity around: ", int(row.lon), ",", int(row.lat))
    return row


def perform_analysis_and_add_columns(results, drop_debugging_info=False, verbosity=0, suppress=0):
    """
    Performs analysis to calculate the Tissot distortion metrics

    Also carries out the discontinuity detection using handle_discontinuities

    Arguments
    ---------
            results: the results dataframe being built
            drop_debugging_info: bool, if True throwaway detailed data when
                    done, if False retain information (default False)
            verbosity: 0: minimal information printed to console >0: too much information
            suppress: threshold for handling discontinuities (default 0 will
                    cause no suppression)
    Returns
    -------
            analysis results DataFrame
    """
    print('Basic analysis...')
    results['error_flag'] = ""

    # Partial derivatives of x, y against lon, lat
    dx_dlon = (results.xE - results.xW) / np.radians(results.lonE - results.lonW)
    dx_dlat = (results.xN - results.xS) / np.radians(results.latN - results.latS)
    dy_dlon = (results.yE - results.yW) / np.radians(results.lonE - results.lonW)
    dy_dlat = (results.yN - results.yS) / np.radians(results.latN - results.latS)

    if not drop_debugging_info:
        results['dx_dlon'] = dx_dlon
        results['dx_dlat'] = dx_dlat
        results['dy_dlon'] = dy_dlon
        results['dy_dlat'] = dy_dlat

    # Earth radius per spherical EPSG 4047
    R = 6371007
    # latitudinal cosine correction
    cos_lat = np.cos(np.radians(results.lat))

    if not drop_debugging_info:
        results['cos_lat'] = cos_lat

    # scale factors in latitudinal and longitudinal directions
    h = np.sqrt(np.square(dx_dlat) + np.square(dy_dlat)) / R
    k = np.sqrt(np.square(dx_dlon) + np.square(dy_dlon)) / (R * cos_lat)

    # area scale factor, s
    sin_theta = (dy_dlat * dx_dlon - dx_dlat * dy_dlon) / (R * R * h * k * cos_lat)
    s = h * k * sin_theta

    if not drop_debugging_info:
        results['sin_theta'] = sin_theta

    # maximum angular distortion, omega
    a_dash = np.sqrt(np.square(h) + np.square(k) + 2 * s)
    b_dash = np.sqrt(np.square(h) + np.square(k) - 2 * s)
    omega = np.degrees(2 * np.arcsin(b_dash / a_dash))

    if not drop_debugging_info:
        results['a_dash'] = a_dash
        results['b_dash'] = b_dash

    # add the results to the data table
    results['h'] = h
    results['k'] = k
    results['s'] = s
    results['omega'] = omega
    results['a'] = (a_dash + b_dash) / 2
    results['b'] = (a_dash - b_dash) / 2

    # This needs to be after all the calculations are done ->
    # Check for discontinuities in the projection.
    print('Discontinuity checking...')
    results = results.apply(handle_discontinuities, axis=1, verbosity=verbosity, suppress=suppress)

    # if we're not debugging, drop the columns we no longer need
    if drop_debugging_info:
        results = results[['ID', 'lon', 'lat', 'x', 'y', 'h', 'k', 's', 'omega', 'error_flag']]

    return results


def summary_analysis_of_files(filelist, folder, results_df=None, verbosity=0, suppress=0):
    """
    Analyses distortion for files in the list, writing results to specified folder

    Optionally adds to a summary dataframe, if provided, otherwise returns summary of
    this set of analyses

    Arguments
    ---------
            filelist: list of files to analyse - must include offsets and lon,
                lat, ID, dir, x, y atttributes
            folder: folder to write results to
            results_df: a dataframe to append summary results to
            verbosity: 0: prints minimal information to console >0 too much information!
            suppress: handling of diagnosed discontinuities default 0 will cause no suppression
    :return:
            the summary results dataframe
    """
    if results_df is None:
        results_df = get_empty_analysis_dataframe()
    for f in filelist:
        source = pandas.read_csv(f)
        res = add_projected_coord_columns(source)
        print('Running analysis on :', f)
        res = perform_analysis_and_add_columns(res, drop_debugging_info=False,
                                               verbosity=verbosity, suppress=suppress)
        basename = f.rsplit('/', 1)[-1].split('.')[0]
        res.to_csv(folder + '/' + basename + '-analysis.csv', index=False)

        # Fill out new row in the results_df dataframe
        new_row = results_df.shape[0] + 1
        filename = f.rsplit('/', 1)[-1]
        results_df.at[new_row, 'full_filename'] = filename

        if '__using__' in filename:
            results_df.at[new_row, 'is_empirical_proj'] = 'True'
            results_df.at[new_row, 'data_filename'] = filename.split('__')[0]
            results_df.at[new_row, 'proj_filename'] = filename.split('__')[-1].split('.')[0]
            results_df.at[new_row, 'proj_num_points'] = int(results_df.at[new_row, 'proj_filename'].split('-')[1])
            results_df.at[new_row, 'proj_spatial_distribution'] = results_df.at[new_row, 'proj_filename'].split('-')[0]
            results_df.at[new_row, 'proj_name'] = results_df.at[new_row, 'proj_filename'].split('-')[-1]
        else:
            results_df.at[new_row, 'is_empirical_proj'] = 'False'
            results_df.at[new_row, 'data_filename'] = ("-".join(filename.split('.')[0].split('-')[:-2]))
            results_df.at[new_row, 'proj_filename'] = ""
            results_df.at[new_row, 'proj_num_points'] = np.NaN
            results_df.at[new_row, 'proj_spatial_distribution'] = ""
            results_df.at[new_row, 'proj_name'] = filename.split('.')[0].split('-')[-1]

        results_df.at[new_row, 'data_num_points'] = int(results_df.at[new_row, 'data_filename'].split('-')[1])
        results_df.at[new_row, 'data_spatial_distribution'] = results_df.at[new_row, 'data_filename'].split('-')[0]

        results_df.at[new_row, 'h_mean'] = res['h'].mean(skipna=True)
        results_df.at[new_row, 'h_stddev'] = res['h'].std(skipna=True)
        results_df.at[new_row, 'k_mean'] = res['k'].mean(skipna=True)
        results_df.at[new_row, 'k_stddev'] = res['k'].std(skipna=True)
        results_df.at[new_row, 'a_mean'] = res['a'].mean(skipna=True)
        results_df.at[new_row, 'a_stddev'] = res['a'].std(skipna=True)
        results_df.at[new_row, 'b_mean'] = res['b'].mean(skipna=True)
        results_df.at[new_row, 'b_stddev'] = res['b'].std(skipna=True)
        results_df.at[new_row, 's_mean'] = res['s'].mean(skipna=True)
        results_df.at[new_row, 's_stddev'] = res['s'].std(skipna=True)
        results_df.at[new_row, 'omega_mean'] = res['omega'].mean(skipna=True)
        results_df.at[new_row, 'omega_stddev'] = res['omega'].std(skipna=True)
    return results_df
