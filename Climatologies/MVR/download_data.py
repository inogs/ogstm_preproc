import copernicusmarine

# no more ppn

copernicusmarine.subset(
  dataset_id="cmems_mod_med_bgc_my_4.2km-climatology_P1M-m",
  variables=["chl_avg", "chl_std", "no3_avg", "no3_std", "o2_avg", "o2_std"],
  minimum_longitude=-5.541666507720947,
  maximum_longitude=36.29166793823242,
  minimum_latitude=30.1875,
  maximum_latitude=45.97916793823242,
  start_datetime="1999-01-01T00:00:00",
  end_datetime="1999-12-01T00:00:00",
  minimum_depth=1.0182366371154785,
  maximum_depth=4152.89599609375,
)
