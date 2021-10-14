"""
Utility functions
"""
import copy

import astropy.units as u
import sunpy.image.coalignment
import sunpy.map
import skimage.feature

__all__ = ['fix_eis_pointing']


def fix_eis_pointing(
        eis_maps,
        aia_map,
        return_correlation=False,
        align_with_center=False,
        reference_index=0,
    ):
    """
    Cross-correlate AIA 171 with EIS Fe XII 195.119 to correct pointing keywords in EIS
    header
    
    Parameters
    ----------
    eis_maps: `dict`
        Each key should correspond to the line ID of the particular spectral
        line and should at least include the line ID used for coalignment
        with AIA 171
    aia_map: `~sunpy.map.Map`
        The AIA map to use for coalignment. If None, the AIA 171
        observation closest to the midpoint of the raster will be downloaded.
    return_correlation: `bool`, optional
        If true, return the cross-correlation between the EIS and AIA maps
    reference_index: `int`
        Index for map in ``eis_maps`` to use to cross-correlate with ``aia_map``
    """
    # Find alignment
    ref_coord, aia_map_resampled = align_eis_with_aia(eis_maps[reference_index], aia_map, return_aia_map=True)
    # Fix EIS maps
    eis_maps_fixed = []
    for m in eis_maps:
        # Set the pointing appropriately. Note that FITS WCS is 1-based
        meta = copy.deepcopy(m.meta)
        # TODO: confirm that this reference coord corresponds to the center of the pixel
        # as opposed to the bottom left (0.5,0.5) or top right (1.5,1.5) corner
        meta['crpix1'] = 1
        meta['crpix2'] = 1
        meta['crval1'] = ref_coord.Tx.to('arcsec').value
        meta['crval2'] = ref_coord.Ty.to('arcsec').value
        eis_maps_fixed.append(m._new_instance(m.data, meta))
    # Optionally calculate the cross-correlation
    if return_correlation:
        corr = skimage.feature.match_template(aia_map_resampled.data, eis_maps[reference_index].data, pad_input=True)
        corr_map = sunpy.map.Map(corr, aia_map_resampled.wcs)
        if align_with_center:
            cor_max = corr_map.pixel_to_world(*np.unravel_index(corr.argmax(), corr.shape)[::-1] * u.pixel)
            for m in eis_maps_fixed:
                del m.meta['crpix1']
                del m.meta['crpix2']
                m.meta['crval1'] = cor_max.Tx.to('arcsec').value
                m.meta['crval2'] = cor_max.Ty.to('arcsec').value
        return eis_maps_fixed, corr_map
    else:
        return eis_maps_fixed


def align_eis_with_aia(m_eis, m_aia, return_aia_map=False):
    """
    Calculate relative shifts to align EIS image with corresponding AIA image
    
    For a given EIS raster and a corresponding AIA image with a similar obstime,
    adjust the EIS pointing by performing a cross-correlation between the two images
    and finding the reference coordinate that 
    
    Parameters
    ----------
    m_eis : `~sunpy.map.Map`
    m_aia : `~sunpy.map.Map`
    return_aia_map : `bool`
        If True, return the AIA map resampled at the EIS resolution. This is
        useful for calculating the cross-correlation between the two maps
        later on (if needed).
    
    Returns
    -------
    reference_coord : `~astropy.coordinate.SkyCoord`
        World coordinate corresponding to the center of the lower left corner
        pixel in the EIS image. The corresponding FITS WCS pixel coordinate
        is (1,1).
    m_aia_r : `~sunpy.map.Map`
        The resampled AIA map. This is used to compute the cross-correlation
        between the the EIS and AIA maps.
    """
    n_x = (m_aia.scale.axis1 * m_aia.dimensions.x) / m_eis.scale.axis1
    n_y = (m_aia.scale.axis2 * m_aia.dimensions.y) / m_eis.scale.axis2
    m_aia_r = m_aia.resample(u.Quantity([n_x, n_y]))
    # Cross-correlate 
    # The resulting "shifts" can be interpreted as the location of the bottom left pixel of the
    # EIS raster in the pixel coordinates of the AIA image. Note that these are in array index
    # coordinates.
    yshift, xshift = sunpy.image.coalignment.calculate_shift(m_aia_r.data, m_eis.data)
    # Find the corresponding coordinate at this pixel position in the resampled AIA map
    reference_coord = m_aia_r.pixel_to_world(xshift, yshift)

    if return_aia_map:
        return reference_coord, m_aia_r
    else:
        return reference_coord