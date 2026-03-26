import numpy as np
import logging
from bitsea.commons.mask import Mask
#from bitsea.commons.grid import RegularGrid
import xarray as xr
from scipy import interpolate
from bitsea.commons import netcdf4
from bitsea.commons import genUserDateList as DL
from pathlib import Path
import matplotlib.pyplot as plt

class co2atm():
    """
        Class co2atm
        author : Giorgio Bolzon
    """

    def __init__(self, conf):
        logging.info("CO2 enabled")
        self.config = conf
        self.path = "co2_monthly_MED_surfatm.nc"
        self.BASE = 366.92  # ppm (mean value for Y=2000)

        with xr.open_dataset(self.path) as ds:
            latitude = ds.latitude.values[-1::-1]
            longitude = ds.longitude.values
            #grid = RegularGrid(lat=latitude, lon=longitude)
            # we are not interested in the time variability, so we take the mean over time
            co2_MMR = ds.co2.mean(dim="valid_time")
        # we convert from mass mixing ratio to ppm
        self.co2_ppm = 28.9644 / 44.0095 * 1.e+6 * co2_MMR
        lat0, lon0 = 35.55, 12.65 # Lampedusa
        Lampedusa=self.co2_ppm.sel(latitude=lat0,longitude=lon0,method="nearest")
        MapCO2_ZEROL = self.co2_ppm - Lampedusa
        # avoiding coastal values larger than 10 ppm
        self.MapCO2_ZEROL = MapCO2_ZEROL.where(MapCO2_ZEROL < 10, other=10)

        X,Y = np.meshgrid(longitude,latitude)
        nPoints = X.size
        Xpoints = np.zeros((nPoints,2),float)
        Xpoints[:,0] = X.ravel()
        Xpoints[:,1] = Y.ravel()
        self.Xpoints = Xpoints

    def generate(self, TheMask: Mask, plot=False):
        outdir= Path(self.config.dir_out)
        Monthly_variability = np.array([
        3.1121, 3.6775, 4.2504, 3.7641,
        1.8368, -0.9948, -4.4677, -6.0282,
        -4.8610, -2.3004, 0.5210, 2.3058
        ], dtype=float)

        landmask = ~ TheMask.mask_at_level(0)

        starttime=f"{self.config.simulation_start_time-1}0101-00:00:00"
        endtime=f"{self.config.simulation_end_time+1}0101-00:00:00"
        Monthly = DL.getTimeList(starttime, endtime, months=1)
        Timeseries = np.zeros(len(Monthly), dtype=float)
        for im, m in enumerate(Monthly):
            fileOUT = outdir / f"CO2_{m.strftime('%Y%m')}15-00:00:00.nc"
            print(fileOUT)
            S_m = Monthly_variability[m.month-1]
            Total_variability = self.BASE + 2.378 * (m.year - 2000) + S_m
            Timeseries[im] = Total_variability
            M2d = self.MapCO2_ZEROL + Total_variability
            co2=M2d.values[0,-1::-1,:]
            f = interpolate.LinearNDInterpolator(self.Xpoints, co2.ravel())
            co2_24  = f(TheMask.xlevels, TheMask.ylevels)
            co2_24[landmask] = 1.e+20
            netcdf4.write_2d_file(co2_24,'CO2',fileOUT,TheMask)
        if plot:
            plt.figure(figsize=(10,5))
            plt.plot(Monthly, Timeseries, marker='o')
            plt.title('CO2 Timeseries')
            plt.xlabel('Date')
            plt.ylabel('CO2 (ppm)')
            plt.grid(True)
            plt.tight_layout()
            plt.savefig('co2_timeseries.png', dpi=150)


        # Save Monthly and Timeseries to text file
        with open('co2_timeseries.txt', 'w') as f:
            f.write("# Date CO2_ppm\n")
            for m, t in zip(Monthly, Timeseries):
                f.write(f"{m.strftime('%Y%m%d')} {t:.3f}\n")



        logging.info("CO2 files written - saving co2_timeseries.txt")

if __name__ == "__main__":
    import config as conf
    CO2 = co2atm(conf)
    CO2.generate(Mask.from_file(conf.file_mask))




