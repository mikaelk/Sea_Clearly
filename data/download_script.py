import pandas as pd
import os

time_array = pd.date_range('2017-01-01 00:00:00','2020-01-01 00:00:00')
str_cmd = 'python -m motuclient --motu https://my.cmems-du.eu/motu-web/Motu --service-id MEDSEA_MULTIYEAR_PHY_006_004-TDS --product-id med-cmcc-cur-rean-d --longitude-min -6 --longitude-max 36.2917 --latitude-min 30.1875 --latitude-max 45.9792 --date-min "%s 12:00:00" --date-max "%s 12:00:00"  --depth-min 1.0182 --depth-max 1.0183  --variable uo --variable vo --out-dir /home/jovyan/Data/Input/CMEMS_MED/ --out-name %s_MEDSEA_MULTIYEAR_PHY.nc --user username --pwd password'

for time_ in time_array[::-1]:

    str_date = str(time_.date())
    cmd = str_cmd % (str_date,str_date,str_date)

    os.system(cmd)
