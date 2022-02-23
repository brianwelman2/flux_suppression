import os, sys
os.chdir("/home/brian/Degrees/masters/code/python/projects/flux_suppression")

import numpy as np
import Tigger
from pyrap.tables import table
from numpy.random import uniform
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt


SPI = 0


def get_field_center(ms):
	with table(ms + "::FIELD") as t:
		phase_centre = t.getcol("PHASE_DIR")[0,0,:]
	
	return phase_centre[0], phase_centre[1] #ra0,dec0 in radians


def get_ms_freq0(ms):
	with table(ms + "::SPECTRAL_WINDOW") as ft:
		chan_freq = ft.getcol("CHAN_FREQ").flatten()
	
	return chan_freq[0]


def generate_pos(fov=None, num_sources=None, order="uniform", inradius=0.15):
	
	if order=="uniform":
		x = np.random.uniform(low=-1*np.abs(fov), high=np.abs(fov), size=num_sources) * (np.pi/180)
		y = np.random.uniform(low=-1*np.abs(fov), high=np.abs(fov), size=num_sources) * (np.pi/180)

	elif order=="circular":
		theta = np.random.uniform(low=0, high=2*np.pi, size=num_sources)
		r = np.random.uniform(low=0, high=np.abs(fov), size=num_sources)
		x = r*np.cos(theta) * (np.pi/180)
		y = r*np.sin(theta) * (np.pi/180)

	elif order=="grid":
		xx = np.linspace(-1*np.abs(fov), np.abs(fov), num_sources//10) * (np.pi/180)
		yy = np.linspace(-1*np.abs(fov), np.abs(fov), num_sources//10) * (np.pi/180)
		x1, y1 = np.meshgrid(xx, yy)
		x = x1.flatten()
		y = y1.flatten()
		arr = np.random.choice(len(x), num_sources, replace=False)
		np.random.shuffle(arr)
		x, y = x[arr], y[arr]

	elif order=="annulus":
		theta = np.random.uniform(low=0, high=2*np.pi, size=num_sources)
		R1, R2 = np.abs(fov), inradius
		r1, r2 = 1.0, R2/R1
		u = np.random.random(size=num_sources) + np.random.random(size=num_sources)
		r = np.array([2-i if i > 1 else i for i in u])
		r = np.array([r2 + i*((R1 - R2)/R2) if i < r2 else i for i in r])
		x = r*np.cos(theta) * (np.pi/180)
		y = r*np.sin(theta) * (np.pi/180)

	else:
		raise RuntimeError(f"Unknown order type '{order}'")

	return x,y


def lm2radec(ms, l, m):#l and m in radians
	rad2deg = lambda val: val * 180.0/np.pi
	ra0,dec0 = get_field_center(ms) # phase centre in radians
	rho = np.sqrt(l**2 + m**2)
	if rho==0:
		ra = ra0
		dec = dec0
	else:
		cc = np.arcsin(rho)
		ra = ra0 - np.arctan2(l*np.sin(cc), rho*np.cos(dec0)*np.cos(cc) - m*np.sin(dec0)*np.sin(cc))
		dec = np.arcsin(np.cos(cc)*np.sin(dec0) + m*np.sin(cc)*np.cos(dec0)/rho)
	return rad2deg(ra), rad2deg(dec)


def radec2lm(ms, ra, dec):
	ra0,dec0 = get_field_center(ms) # phase centre in radians
	ra = np.deg2rad(ra)
	dec = np.deg2rad(dec)

	delta_ra = ra - ra0
	l = np.cos(dec)*np.sin(delta_ra)
	m = np.sin(dec)*np.cos(dec0) - np.cos(dec)*np.sin(dec0)*np.cos(delta_ra)

	return l,m


# def lsmconvert(skymodel):
# 	"""
# 	Converting textfile to lsm.html format
# 	"""
# 	os.system(f"tigger-convert {skymodel} -f")
# 	return "%s.lsm.html"%(skymodel[:-4])


def meqskymodel(ms, point_sources, outfile, DDtags=False, emaj_s=0, emin_s=0, pa_d=0, lm_convert=True, append=False):

	if append:
		str_out = "\n"
		filemode = "a"
	else:
		if DDtags:
			str_out = "#format: name ra_d dec_d  emaj_s emin_s pa_d i q u v spi freq0 tags... \n"
		else:
			str_out = "#format: name ra_d dec_d  emaj_s emin_s pa_d i q u v spi freq0 \n"
		
		filemode = "w"
		
	freq0 = get_ms_freq0(ms)

	for i in range(len(point_sources)):
		amp, l, m = point_sources[i, 0:3]
		try:
			emaj_s, emin_s, pa_d = point_sources[i, 3:6]
		except:
			emaj_s, emin_s, pa_d = 0, 0, 0
		#m = np.absolute(m) if m < 0 else -1*m # changing the signs since meqtrees has its own coordinate system
		if lm_convert:
			ra_d, dec_d = lm2radec(ms, l, m)
		else:
			ra_d, dec_d = l, m 
		name = "A" + str(i)
		if DDtags:
			str_out += "%s %.12g %.12g %.12g %.12g %.5f %.5g 0 0 0 %.5f %f dE\n" % (name, ra_d, dec_d, emaj_s, emin_s, pa_d, amp, SPI, freq0)
		else:
			str_out += "%s %.12g %.12g %.12g %.12g %.5f %.5g 0 0 0 %.5f %f\n" % (name, ra_d, dec_d, emaj_s, emin_s, pa_d, amp, SPI, freq0)

	with open(outfile, filemode) as skymodel:
		skymodel.write(str_out)
	
	# lsm = lsmconvert(outfile)
	
	# model = Tigger.load(lsm)
	# for src in model.sources:
	# 	src.setAttribute("cluster_flux",src.flux.I)
	# 	src.setAttribute("cluster_size",1)
	# 	src.setAttribute("cluster",src.name)
	# 	src.setAttribute("cluster_lead",True)
	# model.save(lsm)



def make_skymodel(msname, nsources, peakflux, modelfrac=[30, 50, 70, 100], pareto_a=3, fov=0.5, prefix="skymodels"):
	"""
	Generate skymodels with random positions and flux following a pareto distribution
	"""

	point_sources = np.zeros((nsources, 3))

	aflux = -np.sort(-np.random.pareto(pareto_a, nsources))
	scale = aflux.max()/peakflux
	flux = aflux/scale

	cumsum = np.cumsum(flux)
	percent = cumsum*100/np.sum(flux)

	if nsources == 2:
		point_sources[0, 1], point_sources[0, 2] = 0.0, 0.0
		point_sources[1, 1], point_sources[1, 2] = generate_pos(fov=fov, num_sources=nsources - 1, order="uniform")
	else:
		point_sources[:, 1], point_sources[:, 2] = generate_pos(fov=fov, num_sources=nsources, order="uniform")

	point_sources[:, 0] = flux
	meqskymodel("ms/" + msname, point_sources, prefix + f"/model-fp-100-src-{nsources}.txt")

	for frac in modelfrac:
		k = np.argmin(np.abs(percent - frac))
		model_sources = point_sources[0:k]
		meqskymodel("ms/" + msname, model_sources, prefix + f"/model-fp-{frac}-src-{nsources}.txt")


if __name__ == "__main__":
	make_skymodel("kat7.ms", 2, 2.0, modelfrac=[], pareto_a=3, fov=0.2, prefix="skymodels/kat7/")
	make_skymodel("vlab.ms", 2, 2.0, modelfrac=[], pareto_a=3, fov=0.2, prefix="skymodels/vlab/")
	make_skymodel("vlab.ms", 100, 2.0, modelfrac=[30, 50, 70], pareto_a=3, fov=0.5, prefix="skymodels/vlab/")
	make_skymodel("meerkat.ms", 2, 2.0, modelfrac=[], pareto_a=3, fov=0.2, prefix="skymodels/meerkat/")
	make_skymodel("meerkat.ms", 100, 2.0, modelfrac=[30, 50, 70], pareto_a=3, fov=0.5, prefix="skymodels/meerkat/")