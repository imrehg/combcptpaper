#!/usr/bin/env python

#from numpy import *
from pylab import *
from scipy import *
from scipy.optimize import *

filename = '8700Ne_buffer.CSV'
data = loadtxt(filename, skiprows=1, delimiter=',', usecols=[0,1])

freq = data[:,0] / 1000 * 100.0 
fluo = data[:,1] 
fluo = fluo / max(fluo) * 0.05 + 0.95

dd = zip(freq, fluo)
freq, fluo = zip(*sorted(dd))

err = array([0.005]*len(fluo))

fitfunc = lambda p, x : p[2]*exp(-(x - p[0])**2 / (2*p[1]**2)) + p[3]
#errfunc = lambda p, x, y, err : (y - fitfunc(p, x)) / err
errfunc = lambda p, x, y, err : (y - fitfunc(p, x)) / err


#pin = [0.1, 2, -5000, 5000]
pin = [0.1, 2, -1, 1]
pout, cov, info, msg, success = leastsq(errfunc, pin[:], args=(freq,fluo, err),full_output=1)


print pout
#print data
#print cov
cerr = sqrt( cov[1][1] ) * 2*log(2)*sqrt(2)
print cerr
#print sqrt(cov[0][0])
#popt, pcov = curve_fit(fitfunc, freq, fluo)
#print popt, pcov

#print pout


b = pout[0]
c = pout[1]
x1 = -sqrt(-2*c**2*log(0.5))+b
x2 = sqrt(-2*c**2*log(0.5))+b
halfmax = pout[3]+pout[2]/2
print "Fit contrast: %.1f%%" %(abs(pout[3]*pout[2])*100)
print "FWHM: %.1f Hz" %(x2-x1)


fitfreq = linspace(freq[0]-0.5,freq[-1]+0.5)
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
errorbar(freq, fluo, err, fmt='k.', linestyle='None', markersize=10)
plot(fitfreq, fitfluo,'k-')
plot([x1, x2],[halfmax, halfmax], 'k--')
text(b, halfmax+0.005, "FWHM = %.2f(%d) Hz" %(x2-x1, int(cerr*100)), horizontalalignment='center', verticalalignment='center', fontsize=15)
xlim([fitfreq[0], fitfreq[-1]])
ylim([min(fluo)-0.005, max(fluo)+0.005])
#plot([freq[0], freq[-1]],[pout[3],pout[3]])
xlabel('$100f_{rep} - f_{clock}$ (Hz)')
ylabel('Fluorescence (a.u.)')


#figure(2)
#plot(freq,errfunc(pout, freq, fluo),'.')
savefig("8700Ne_buffer.eps")
show()
