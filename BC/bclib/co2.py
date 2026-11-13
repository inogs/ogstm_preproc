import numpy as np
import matplotlib.pyplot as plt
import logging
#from bitsea.commons.grid import RegularGrid
import xarray as xr
import xesmf as xe
from scipy import interpolate
from bitsea.commons.mask import Mask
from bitsea.commons import netcdf4
from bitsea.commons import genUserDateList as DL
from pathlib import Path

class co2atm():
    """
        Class co2atm
        author : Giorgio Bolzon
    """

    def __init__(self, conf):
        """
        Initialize the co2atm class by
        - reading the CO2 data from the specified CAMS reanalysis file
        - averaging and converting it from mass mixing ratio to ppm
        - generating a map of CO2 anomalies with respect to Lampedusa, and
        - preparing it for interpolation.
          """
        logging.info("CO2 enabled")
        self.config = conf
        self.path = conf.file_co2
        self.BASE = 366.92  # ppm (mean value for Y=2000)

        with xr.open_dataset(self.path) as ds:
            latitude = ds.latitude.values[-1::-1]
            longitude = ds.longitude.values
            #grid = RegularGrid(lat=latitude, lon=longitude)
            # we are not interested in the time variability, so we take the mean over time
            co2_MMR = ds.co2.mean(dim="valid_time").load()
        # we convert from mass mixing ratio to ppm
        self.co2_ppm = 28.9644 / 44.0095 * 1.e+6 * co2_MMR
        lat0, lon0 = 35.55, 12.65 # Lampedusa
        Lampedusa=self.co2_ppm.sel(latitude=lat0,longitude=lon0,method="nearest")
        MapCO2_ZEROL = self.co2_ppm - Lampedusa
        # avoiding coastal values larger than 10 ppm
        self.MapCO2_ZEROL = MapCO2_ZEROL.where(MapCO2_ZEROL < 10, other=10)

        self.X, self.Y = np.meshgrid(longitude,latitude)
        nPoints = self.X.size
        Xpoints = np.zeros((nPoints,2),float)
        Xpoints[:,0] = self.X.ravel()
        Xpoints[:,1] = self.Y.ravel()
        self.Xpoints = Xpoints

    def regridder_2D(self):
        conf = self.config
        with xr.open_dataset(conf.file_mask) as M:
            dst_grid = xr.Dataset({
                        "lon": (("y", "x"), M['glamt'].values[:].squeeze()),
                        "lat": (("y", "x"), M['gphit'].values[:].squeeze()),
                        "mask": (("y", "x"), M['tmask'].values[0,0,:,:].astype(np.int32)),
                        })
        #
        src_grid = xr.Dataset({"lon": (("y", "x"), self.X), 
                               "lat": (("y", "x"), self.Y)})
        #
        regridder = xe.Regridder(
                src_grid,
                dst_grid,
                method = "nearest_s2d",
                periodic = False,
                unmapped_to_nan=True,
                filename="weights_co2_rgrd.nc",
                reuse_weights=False,
                )
        return regridder


    def generate(self, TheMask: Mask, plot=False):
        """ Generate CO2 files for ogstm for the whole period of simulation, with monthly variability.
        Generates also co2_timeseries.txt and co2_timeseries.png (if plot=True).
        The CO2 files are written in the output directory specified in the configuration,
        with names like CO2_yyyymm15-00:00:00.nc.
        """
        outdir= Path(self.config.dir_out)
        Monthly_variability = np.array([
        3.1121, 3.6775, 4.2504, 3.7641,
        1.8368, -0.9948, -4.4677, -6.0282,
        -4.8610, -2.3004, 0.5210, 2.3058
        ], dtype=float)

        landmask = ~ TheMask.mask_at_level(0)

        # Find the model grid index nearest to Lampedusa (lat=35.55, lon=12.65)
        lat0, lon0 = 35.55, 12.65
        dist2 = (TheMask.xlevels - lon0)**2 + (TheMask.ylevels - lat0)**2
        lamp_idx = np.unravel_index(np.argmin(dist2), dist2.shape)

        starttime=f"{self.config.simulation_start_time-1}1201-00:00:00"
        endtime=f"{self.config.simulation_end_time+1}0101-00:00:00"
        Monthly = DL.getTimeList(starttime, endtime, months=1)
        Timeseries = np.zeros(len(Monthly), dtype=float)
        Timeseries_Lampedusa = np.zeros(len(Monthly), dtype=float)

        regridder = self.regridder_2D()
        for im, m in enumerate(Monthly):
            fileOUT = outdir / f"CO2_{m.strftime('%Y%m')}15-00:00:00.nc"
            print(fileOUT)
            S_m = Monthly_variability[ m.month - 1 ]
            time_co2_value = self.BASE + 2.378 * (m.year - 2000) + S_m
            Timeseries[im] = time_co2_value

            # now, adding space variability
            M2d = self.MapCO2_ZEROL + time_co2_value
            co2=M2d.values[0,-1::-1,:]
            #f = interpolate.LinearNDInterpolator(self.Xpoints, co2.ravel())
            #co2_24  = f(TheMask.xlevels, TheMask.ylevels)
            co2_24 = regridder(co2)
            Timeseries_Lampedusa[im] = co2_24[lamp_idx]
            #co2_24[landmask] = 1.e+20
            netcdf4.write_2d_file(co2_24,'CO2',fileOUT,TheMask)
        if plot:
            plt.figure(figsize=(10,5))
            plt.plot(Monthly, Timeseries, marker='o')
            plt.title('CO2 Timeseries')
            plt.xlabel('Date')
            plt.ylabel('CO2 (ppm)')
            plt.grid(True)
            plt.tight_layout()
            plt.savefig(outdir / 'co2_timeseries.png', dpi=150)


        # Save Monthly and Timeseries to CSV file with interpolated value at Lampedusa
        with open(outdir / 'co2_timeseries.txt', 'w') as f:
            f.write("Date,CO2_scalar_ppm,CO2_Lampedusa_ppm\n")
            for m, t, tl in zip(Monthly, Timeseries, Timeseries_Lampedusa):
                f.write(f"{m.strftime('%Y%m%d')},{t:.3f},{tl:.3f}\n")



        logging.info("CO2 files written - saving co2_timeseries.txt")

if __name__ == "__main__":
    # to make it import config
    import sys
    import os
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
    import config as conf
    CO2 = co2atm(conf)
    CO2.generate(Mask.from_file(conf.file_mask))




