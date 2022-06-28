import xarray as xr
import numpy as np


def output_refinery_old(data, X, Y, x, y):
    """
    The function selects the trajectories of the particles that are inside a
    grid cell at the start of the simulation and drops the particles from the
    xarray dataset.

    Parameters
    ----------
    data: xarray dataset
        Specifically built for parcels outputs with variables 'trajectory',
        'lat', and 'lon'
    X, Y: Either a meshgrid 2D or a 1D array
        Contians the coordinates of the grid.
    x, y: float
        Position
    Returns
    -------
    data_relevant: xarray dataset
        xrray with containing the particles that are in a grid cell at the
        start of the simulation.
    """

    if len(X.shape) == 2:
        x_range = X[0]
        y_range = Y[:, 0]

    elif len(X.shape) == 1:
        x_range = X
        y_range = Y

    lons = data['lon'][:, 0]
    lats = data['lat'][:, 0]

    lon_idx = np.digitize(lons, x_range)
    lat_idx = np.digitize(lats, y_range)

    xloc_idx = np.digitize(x, x_range)
    yloc_idx = np.digitize(y, y_range)

    intersection = (lon_idx == xloc_idx) & (lat_idx == yloc_idx)
    index = np.where(intersection)[0]

    data_relevant = data.where(data['trajectory'].isin(index), drop=True)

    return data_relevant

def output_refinery(data, X, Y, x, y):
    """
    The function selects the trajectories of the particles that are inside a grid cell
    at the start of the simulation and drops the particles from the xarray dataset. 
    
    Parameters
    ----------
    data: xarray dataset
        Specifically built for parcels outputs with variables 'trajectory', 'lat', and 'lon'
    X, Y: Either a meshgrid 2D or a 1D array
        Contians the coordinates of the grid.
    x, y: float
        Position 
    Returns
    -------
    data_relevant: xarray dataset
        xrray with containing the particles that are in a grid cell at the start of the simulation.
    """
    
    if len(X.shape) == 2:
        x_range = X[0]
        y_range = Y[:, 0]

    elif len(X.shape) == 1:
        x_range = X
        y_range = Y

    xloc_idx = np.digitize(x, x_range)
    yloc_idx = np.digitize(y, y_range)

    return data.where((data['lon'][:,0] >= x_range[xloc_idx-1]) & (data['lon'][:,0] < x_range[xloc_idx]) &
                       (data['lat'][:,0] >= y_range[yloc_idx-1]) & (data['lat'][:,0] < y_range[yloc_idx]),drop=True)


def output_refinery_nearby(data, x, y, dlon=.05, dlat=.05):
    """
    The function selects the trajectories of the particles that are within dlon/dlat distance (square)
    at the start of the simulation and drops the particles from the xarray dataset. 
    
    Parameters
    ----------
    data: xarray dataset
        Specifically built for parcels outputs with variables 'trajectory', 'lat', and 'lon'
    X, Y: Either a meshgrid 2D or a 1D array
        Contians the coordinates of the grid.
    x, y: float
        Position 
    Returns
    -------
    data_relevant: xarray dataset
        xrray with containing the particles that are in a grid cell at the start of the simulation.
    """
    

    return data.where((data['lon'][:,0] >= x-dlon) & (data['lon'][:,0] < x+dlon) &
                       (data['lat'][:,0] >= y-dlat) & (data['lat'][:,0] < y+dlat),drop=True)