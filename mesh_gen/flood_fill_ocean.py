import numpy as np
from collections import deque


def flood_fill_nemo(tmask, umask, vmask, fmask, i0, j0, blocked_edges=None):
    """
    Remove a region from a NEMO C-grid using face-connected flood fill.

    Parameters
    ----------
    tmask, umask, vmask, fmask : numpy.ndarray
        Four-dimensional masks with dimensions (time, depth, y, x).
        Wet points are 1; land points are 0.

    i0, j0 : int
        Zero-based coordinates of a seed T point in the region
        to be removed.

    blocked_edges : iterable, optional
        Additional barriers between adjacent T points.

        Example:
            [((299, 943), (299, 944))]

        Each edge is represented by two (j, i) coordinates.
        Both cells remain unchanged, but the flood fill cannot
        cross the specified edge.

    Returns
    -------
    tmask_new, umask_new, vmask_new, fmask_new : numpy.ndarray
        Updated masks with the same shapes as the input arrays.

    region_2d : numpy.ndarray of bool
        Two-dimensional surface region reached by the flood fill.

    Notes
    -----
    1. The flood fill operates on surface T points.
    2. U and V masks determine whether adjacent T points are connected.
    3. T, U, V and F masks are modified through their horizontal
       surface footprints, using the original masks at every depth.
    4. Array boundaries are closed; periodic wrapping is not assumed.
    """

    # ---------------------------------------------------------
    # 1. Validate inputs
    # ---------------------------------------------------------

    for name, mask in (("tmask", tmask), ("umask", umask), ("vmask", vmask), ("fmask", fmask)):
        if not isinstance(mask, np.ndarray):
            raise TypeError(f"{name} must be a NumPy array.")

        if mask.ndim != 4:
            raise ValueError(f"{name} must have dimensions (time, depth, y, x).")

    if tmask.shape != umask.shape or tmask.shape != vmask.shape:
        raise ValueError("tmask, umask and vmask must have identical shapes.")

    nt, nz, ny, nx = tmask.shape

    if not (0 <= i0 < nx and 0 <= j0 < ny):
        raise ValueError("Seed point lies outside the grid.")

    # ---------------------------------------------------------
    # 2. Extract the surface masks
    # ---------------------------------------------------------

    # Use the first time record and the surface level.
    # The original arrays are retained for vertical propagation.

    t2d = tmask[0, 0, :, :] != 0
    u2d = umask[0, 0, :, :] != 0
    v2d = vmask[0, 0, :, :] != 0

    if not t2d[j0, i0]:
        raise ValueError("Seed point is not a wet surface T point.")

    # ---------------------------------------------------------
    # 3. Encode optional blocked edges
    # ---------------------------------------------------------

    blocked = set()

    if blocked_edges is not None:

        for p1, p2 in blocked_edges:

            j1, i1 = p1
            j2, i2 = p2

            for j, i in (p1, p2):
                if not (0 <= j < ny and 0 <= i < nx):
                    raise ValueError(f"Barrier edge {(p1, p2)} lies outside the grid.")

            if abs(j1 - j2) + abs(i1 - i2) != 1:
                raise ValueError(f"Barrier cells {p1} and {p2} must be orthogonal neighbours.")

            blocked.add(frozenset((p1, p2)))

    # ---------------------------------------------------------
    # 4. Define whether two adjacent T points are connected
    # ---------------------------------------------------------

    def can_cross(j, i, nj, ni):

        # Both T points must be wet.
        if not t2d[nj, ni]:
            return False

        if not t2d[j, i]:
            return False

        # Explicitly blocked edges are never crossed.
        edge = frozenset(((j, i), (nj, ni)))

        if edge in blocked:
            return False

        dj = nj - j
        di = ni - i

        # East: U[j, i] is the face between T[j, i] and T[j, i+1].
        if dj == 0 and di == 1:
            return u2d[j, i]

        # West: U[j, i-1] is the face between T[j, i] and T[j, i-1].
        if dj == 0 and di == -1:
            return u2d[j, i - 1]

        # North: V[j, i] is the face between T[j, i] and T[j+1, i].
        if dj == 1 and di == 0:
            return v2d[j, i]

        # South: V[j-1, i] is the face between T[j, i] and T[j-1, i].
        if dj == -1 and di == 0:
            return v2d[j - 1, i]

        # No diagonal movement.
        return False

    # ---------------------------------------------------------
    # 5. Flood fill using actual U/V face connectivity
    # ---------------------------------------------------------

    region_2d = np.zeros((ny, nx), dtype=bool)

    queue = deque([(j0, i0)])
    region_2d[j0, i0] = True

    # East, west, north, south.
    neighbours = [(0, 1), (0, -1), (1, 0), (-1, 0)]

    while queue:

        j, i = queue.popleft()

        for dj, di in neighbours:

            nj = j + dj
            ni = i + di

            # Do not wrap across array boundaries.
            if not (0 <= nj < ny and 0 <= ni < nx):
                continue

            if region_2d[nj, ni]:
                continue

            if can_cross(j, i, nj, ni):

                region_2d[nj, ni] = True
                queue.append((nj, ni))

    # ---------------------------------------------------------
    # 6. Construct the horizontal footprints to remove
    # ---------------------------------------------------------

    # T points: remove all T points reached by the flood fill.
    region_t = region_2d.copy()

    # U points:
    # U[j, i] lies between T[j, i] and T[j, i+1].
    #
    # Mask U if either adjacent T point belongs to the
    # region being removed.

    region_u = region_2d.copy()

    region_u[:, :-1] |= region_2d[:, 1:]

    # V points:
    # V[j, i] lies between T[j, i] and T[j+1, i].
    #
    # Mask V if either adjacent T point belongs to the region.

    region_v = region_2d.copy()

    region_v[:-1, :] |= region_2d[1:, :]

    # F points:
    # F[j, i] is assumed to be the northeast corner of T[j, i].
    #
    # Mask F if any of its four surrounding T points belongs
    # to the region being removed.

    region_f = region_2d.copy()

    region_f[:-1, :] |= region_2d[1:, :]
    region_f[:, :-1] |= region_2d[:, 1:]
    region_f[:-1, :-1] |= region_2d[1:, 1:]

    # ---------------------------------------------------------
    # 7. Apply the horizontal changes at every depth and time
    # ---------------------------------------------------------

    # Make copies so that the input arrays remain unchanged.

    tmask_new = tmask.copy()
    umask_new = umask.copy()
    vmask_new = vmask.copy()

    # Broadcast the 2D region through time and depth.
    tmask_new *= (~region_t)[None, None, :, :]
    umask_new *= (~region_u)[None, None, :, :]
    vmask_new *= (~region_v)[None, None, :, :]

    fmask_new = fmask.copy()
    fmask_new *= (~region_f)[None, None, :, :]

    return (tmask_new, umask_new, vmask_new, fmask_new)


if __name__=="__main__":
    import xarray as xr
    # load mask
    maskfile = '/leonardo_work/RIVRD_proc/CMCC/NEMO_RUN_TEST_1979_1988/med24v1_mesh_mask_BGC.nc'
    M = xr.open_dataset(maskfile)
    tmask = M['tmask'].values[:]
    umask = M['umask'].values[:]
    vmask = M['vmask'].values[:]
    fmask = M['fmask'].values[:]
    M.close()
    # a point in the Marmara Sea
    i0 = 980
    j0 = 320
    # close the dardanelles!
    blocked_edges = [((299, 943), (299, 944))]
    # fill in all of marmara and black sea, leave 1pt at dardanelles for boundary
    tmask_new, umask_new, vmask_new, fmask_new = flood_fill_nemo(tmask, umask, vmask, fmask, i0, j0, blocked_edges=blocked_edges)


