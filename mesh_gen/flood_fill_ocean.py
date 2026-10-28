import numpy as np
from collections import deque

def flood_fill_ocean(mask, i0, j0, blocked_edges=None):
    """
    Identify all wet grid cells connected to a starting point.

    Parameters
    ----------
    mask : 2D numpy array
        NEMO mask: 1 (wet), 0 (land).
        Array dimensions are (t, z, j, i).
    i0, j0 : int
        Starting grid indices, using Python's 0-based indexing.
    blocked_edges : iterable of pairs, optional
        Each pair defines two adjacent grid cells across which
        the flood fill must not propagate.
        Example:
            [((j1, i1), (j2, i2)), ((j3, i3), (j4, i4))]
        Use these to block connections across the Dardanelles
        and Bosphorus without removing their wet cells.

    Returns
    -------
    region_2d : numpy.ndarray of bool, shape (y, x)
        True for horizontal grid cells belonging to the region.

    region_4d : numpy.ndarray of bool, shape (1, z, y, x)
        True for wet 3D cells belonging to the region.
        Land cells and cells outside the region are False.
    """

    _, nz, ny, nx = mask.shape

    if not (0 <= i0 < nx and 0 <= j0 < ny):
        raise ValueError("Starting point is outside the grid.")

    if not mask[0, 0, j0, i0]:
        raise ValueError("Starting point is not a wet T point.")

    # A horizontal cell is wet if it is wet at any depth.
    wet_2d = np.any(mask[0] != 0, axis=0)

    # --------------------------------------------------
    # 3. Encode barriers between adjacent cells
    # --------------------------------------------------
    # Each edge is represented by its two (j, i) endpoints.
    # This prevents crossing a strait without masking
    # either cell bordering the strait.

    blocked = set()

    if blocked_edges is not None:
        for p1, p2 in blocked_edges:
            j1, i1 = p1
            j2, i2 = p2

            if not (0 <= i1 < nx and 0 <= j1 < ny and 0 <= i2 < nx and 0 <= j2 < ny):
                raise ValueError(f"Barrier edge {p1, p2} is outside the grid.")

            # The cells must share a horizontal or vertical edge.
            #if abs(i1 - i2) + abs(j1 - j2) != 1:
            #    raise ValueError(f"Barrier cells {p1} and {p2} must be orthogonal neighbours.")

            blocked.add(frozenset((p1, p2)))

    # --------------------------------------------------
    # 4. Flood fill using four-neighbour connectivity
    # --------------------------------------------------
    region_2D = np.zeros((ny, nx), dtype=bool)
    queue = deque([(j0, i0)])
    region_2D[j0, i0] = True

    # Four-neighbour connectivity: east, west, north, south
    #neighbours = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    # Eight-neighbour connectivity:
    neighbours = [(1, 0), (-1, 0), (0, 1), (0, -1),
                  (1, 1), (-1, 1), (1, -1), (-1, -1)]

    while queue:
        j, i = queue.popleft()

        for dj, di in neighbours:
            ni = i + di
            nj = j + dj

            # Do not cross the boundaries of the array
            if not (0 <= ni < nx and 0 <= nj < ny):
                continue

            # Do not cross blocked connections.
            edge = frozenset(((j, i), (nj, ni)))
            if edge in blocked:
                continue

            # Only visit wet cells that have not been visited
            if wet_2d[nj, ni] and not region_2D[nj, ni]:
                region_2D[nj, ni] = True
                queue.append((nj, ni))
    # return 4D mask
    region_4D = mask[0, :, :, :] * region_2D[None, :, :]
    region_4D = region_4D[np.newaxis, :, :, :]
    return region_2D, region_4D


if __name__=="__main__":
    import xarray as xr
    maskfile = '/leonardo_work/RIVRD_proc/CMCC/NEMO_RUN_TEST_1979_1988/med24v1_mesh_mask_BGC.nc'
    M = xr.open_dataset(maskfile)
    #mask = M['tmask'].values[0,0,:,:]
    mask = M['tmask'].values[:]
    # a point in the Marmara Sea
    i0 = 980
    j0 = 320
    # fill in all of marmara and black sea, leave 1pt at dardanelles for boundary
    blocked_edges = [((299, 943), (299, 944))]
    M.close()
    region_2D, region_4D = flood_fill_ocean(mask, i0, j0, blocked_edges=blocked_edges)
