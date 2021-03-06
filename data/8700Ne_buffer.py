#!/usr/bin/env python

from pylab import *
from scipy import *
import scipy.odr as odr

filename = '8700Ne_buffer.CSV'
data = loadtxt(filename, skiprows=1, delimiter=',', usecols=[0,1])

freq = data[:,0] / 1000 * 100.0 
fluo = data[:,1] 
fluo = fluo / max(fluo)

dd = zip(freq, fluo)
freq, fluo = zip(*sorted(dd))

# Gaussian
fitfunc = lambda p, x : p[2]*exp(-(x - p[0])**2 / (2*p[1]**2)) + p[3]
## Lorentzian
#fitfunc = lambda p, x : p[2] * (p[1] / ((x - p[0])**2 + p[1]**2))/pi + p[3]
pin = [0.1, 2, -1, 1]

fitdata = odr.RealData(freq, fluo)
model = odr.Model(fitfunc)
fit = odr.ODR(fitdata, model, pin)
fit.set_job(fit_type=2)
output = fit.run()
pout = output.beta

b = pout[0]
c = pout[1]
sd_c = output.sd_beta[1]
fwhmscale = 2*sqrt(2*log(2))
fwhm = c * fwhmscale
fwhmsd = sd_c * fwhmscale
halfmax = fitfunc(pout, b+fwhm/2.0)

print "Fit contrast: %.1f%%" %(abs(pout[3]*pout[2])*100)
print "FWHM: %.1f(%.0f) Hz" %(fwhm, fwhmsd*10)


fitfreq = linspace(freq[0]-0.3,freq[-1]+0.3)
fitfluo = fitfunc(pout, fitfreq)



fig_width_pt = 510.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
params = {'backend': 'ps',
          'axes.labelsize': 20,
          'text.fontsize': 10,
          'legend.fontsize': 10,
          'xtick.labelsize': 18,
          'ytick.labelsize': 18,
          'text.usetex': True,
          'figure.figsize': fig_size}
rcParams.update(params)

axes([0.125,0.2,0.95-0.125,0.95-0.2])

figure(1)
plot(freq, fluo,'k.', linestyle='None', markersize=10)
plot(fitfreq, fitfluo,'k-')
plot([b-fwhm/2.0, b+fwhm/2.0],[halfmax, halfmax], 'k--')
text(b, halfmax+0.06, "FWHM = %.1f(%.0f) Hz" %(fwhm, fwhmsd*10), horizontalalignment='center', verticalalignment='center', fontsize=15)
xlim([fitfreq[0], fitfreq[-1]])
ylim([min(fluo)-0.05, max(fluo)+0.05])
xlabel('$100f_{rep} - f_{clock}$ (Hz)')
ylabel('Fluorescence (a.u.)')


savefig("8700Ne_buffer.eps")
show()
