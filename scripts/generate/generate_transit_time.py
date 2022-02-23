import ephem
import datetime

meerkat = ephem.Observer()
meerkat.lat, meerkat.long, meerkat.elevation = "-30:16:41.00", "11:49:36.0", 1054
meerkat.epoch = ephem.J2000
meerkat.date=datetime.datetime.utcnow()
meerkat_field=ephem.readdb("SCP,f|J,00:00:00.00,00:00:00:00,1054.0")
meerkat_field.compute(meerkat)
print("MEERKAT:", meerkat.next_transit(meerkat_field))

vla=ephem.Observer()
vla.date=datetime.datetime.utcnow()
vla.lat, vla.long, vla.elevation = "-30:16:41.00", "11:49:36.0", 1054
vla.epoch = ephem.J2000
vla.date=datetime.datetime.utcnow()
vla_field=ephem.readdb("SCP,f|J,00:00:00.00000,00:00:00:00000,1054.0")
vla_field.compute(vla)
print("VLA:", vla.next_transit(vla_field))
