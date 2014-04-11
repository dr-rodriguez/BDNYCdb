#!/usr/bin/python
# Utilities
import warnings, glob, os, re, xlrd, cPickle, astropy.units as q, astropy.constants as ac, numpy as np, matplotlib.pyplot as plt, astropy.coordinates as apc, astrotools as a
from random import random
from heapq import nsmallest, nlargest
from scipy.interpolate import Rbf
from itertools import chain
from pysynphot import observation
from pysynphot import spectrum
warnings.simplefilter('ignore')
path = '/Users/Joe/Documents/Python/'

def app2abs(magnitude, distance):
  if isinstance(magnitude,tuple) and isinstance(distance,tuple):
    (m, sig_m), (d, sig_d) = magnitude, distance
    sig_M = np.sqrt(sig_m**2 + 25*(sig_d/d).value**2)
    M = (m-5*np.log10(d/(10*q.pc)))
    return (M, sig_M)
  else: return (magnitude-5*np.log10(distance/(10*q.pc)))  
  
def blackbody(lam, T, Flam=False, radius=1, dist=10, emitted=False):
  '''
  Given a wavelength array [um] and temperature [K], returns an array of Planck function values in [erg s-1 cm-2 A-1]
  '''
  lam, T = lam.to(q.cm), T*q.K
  I = np.pi*(2*ac.h*ac.c**2 / (lam**(4 if Flam else 5) * (np.exp((ac.h*ac.c / (lam*ac.k_B*T)).decompose()) - 1))).to(q.erg/q.s/q.cm**2/(1 if Flam else q.AA))
  return I if emitted else I*((radius*ac.R_jup)**2/(dist*q.pc)**2).decompose()

def ChiSquare(a, b, unc=None, array=False, Gtest=False, norm=True, log=True):
  a, b = [np.array(map(float,i.value)) if hasattr(i,'_unit') else np.array(map(float,i)) for i in [a,b]]
  c, variance = np.array(map(float,unc.value)) if hasattr(unc, '_unit') else np.array(map(float,np.ones(len(a)))), np.std(b)**4 # Since the standard deviation is the root of the variance
  X2 = np.array([(j*np.log(j/i)/k)**2/i for i,j,k in zip(a,b,c)]) if Gtest else np.array([((i-j)/k)**2/variance for i,j,k in zip(a,b,c)])    
  if norm: X2 = abs(np.log10(X2/min(X2))/np.log10(min(X2))) if log else X2/max(X2)
  return X2 if array else sum(X2) 
  
def contour_plot(x, y, z, best=False, figsize=(8,8), levels=20, cmap=plt.cm.jet):
  from scipy.interpolate import Rbf
  from itertools import chain
  xi, yi = np.meshgrid(np.linspace(min(x), max(x), 500), np.linspace(min(y), max(y), 25))
  rbf = Rbf(x, y, z, function='linear')
  zi = rbf(xi, yi)
  plt.figure(figsize=figsize), plt.contourf(xi, yi, zi, levels, cmap=cmap), plt.colorbar(), plt.xlim(min(x),max(x)), plt.ylim(min(y),max(y)), plt.xlabel('Teff'), plt.ylabel('log(g)')
  if best:
    coords = min(zip(*[list(chain.from_iterable(zi)),list(chain.from_iterable(xi)),list(chain.from_iterable(yi))]))[1:]
    plt.title('Teff = {}, log(g) = {}'.format(*coords)), plt.plot(*coords, c='white', marker='x')
  
def deg2sxg(ra='', dec=''):
  RA, DEC = '', ''
  if ra: RA = str(apc.angles.Angle(ra,'degree').format(unit='hour', sep=' ')) 
  if dec: DEC = str(apc.angles.Angle(dec,'degree').format(unit='degree', sep=' ')) 
  return (RA, DEC) if ra and dec else RA or DEC 
  
def dict2txt(DICT, writefile, column1='-', delim='\t', digits=None, order='', append=False):
  '''
  Given a nested dictionary *DICT*, writes a .txt file with keys as columns. 
  '''
  import csv
  D = DICT.copy()
  with open( writefile, 'a+' if append else 'w' ) as f:
    writer, w = csv.writer(f, delimiter=delim), []
    for k in D.keys():
      w.append(k)
      for i in D[k].keys():
        if digits: D[k][i] = '-' if not D[k][i] else '{:.{}f}'.format(D[k][i],digits) if isinstance(D[k][i],(float,int)) else '{:.{}f}'.format(float(D[k][i]),digits) if D[k][i].replace('.','').replace('-','').isdigit() else str(D[k][i])
        else: D[k][i] = '-' if not D[k][i] else '{}'.format(D[k][i]) if isinstance(D[k][i],(float,int)) else '{}'.format(float(D[k][i])) if D[k][i].replace('.','').replace('-','').isdigit() else str(D[k][i])
        w.append(i), w.append(str(D[k][i]))
    width = len(max(map(str,w), key=len))
    head = ['{!s:{}}'.format(column1,width)]
    headorder = order or sorted(D[D.keys()[0]].keys())
    for i in headorder: head.append('{!s:{}}'.format(i,width))
    if delim == ',': head = [i.replace(' ','') for i in head]
    writer.writerow(head)
    for i in sorted(D.keys()):
      order = order or sorted(D[i].keys())
      row = ['{!s:{}}'.format(i,width)]
      for k in order:
        if k not in D[i].keys(): D[i][k] = '' if delim==',' else '-'
        row.append('{!s:{}}'.format(D[i][k],width))
      if delim == ',': row = [i.replace(' ','') for i in row]
      writer.writerow(row)
      
def distance(coord1, coord2):
  '''
  Given n-dimensional coordinates of two points, returns the distance between them
  '''
  return np.sqrt(sum([abs(i-j)**2 for i,j in zip(coord1,coord2)]))
  
def find(filename, tree):
  '''                                                                               
  For given filename and directory tree, returns the path to the file. 
  For only file extension given as filename, returns list of paths to all files with that extnsion in that directory tree.  

  *filename*
    Filename or file extension to search for (e.g. 'my_file.txt' or just '.txt')
  *tree*
    Directory tree base to start the walk (e.g. '/Users/Joe/Documents/')
  '''
  import os
  result = []

  for root, dirs, files in os.walk(tree):
    if filename.startswith('.'):
      for f in files:
        if f.endswith(filename):
          result.append(os.path.join(root, f))
    else:  
      if filename in files:
        result.append(os.path.join(root, filename))

  return result

def get_filters(filter_directories=['{}Filters/{}/'.format(path,i) for i in ['2MASS','SDSS','WISE','IRAC','HST','Bessel','MKO']], systems=['2MASS','SDSS','WISE','IRAC','HST','Bessel','MKO']):
  '''
  Grabs all the .txt spectral response curves and returns a dictionary of wavelength array [um], filter response [unitless], effective, min and max wavelengths [um], and zeropoint [erg s-1 cm-2 A-1]. 
  '''
  files = glob.glob(filter_directories+'*.txt') if isinstance(filter_directories, basestring) else [j for k in [glob.glob(i+'*.txt') for i in filter_directories] for j in k]

  if len(files) == 0: print 'No filters in', filter_directories
  else:
    filters = {}
    for filepath in files:
      filter_name = os.path.splitext(os.path.basename(filepath))[0]
      RSR_x, RSR_y = [np.array(map(float,i)) for i in zip(*txt2dict(filepath,to_list=True,skip=['#']))]
      RSR_x, RSR_y = (RSR_x*(q.um if min(RSR_x)<100 else q.AA)).to(q.um), RSR_y*q.um/q.um
      Filt = a.filter_info(filter_name)
      filters[filter_name] = {'wav':RSR_x, 'rsr':RSR_y, 'system':Filt['system'], 'eff':Filt['eff']*q.um, 'min':Filt['min']*q.um, 'max':Filt['max']*q.um, 'ext':Filt['ext'], 'ABtoVega':Filt['ABtoVega'], 'zp':Filt['zp']*q.erg/q.s/q.cm**2/q.AA, 'zp_photon':Filt['zp_photon']/q.s/q.cm**2/q.AA }

    for i in filters.keys():
      if filters[i]['system'] not in systems: filters.pop(i)    
    return filters

def goodness(spec1, spec2, array=False, exclude=[], filt_dict=None, weighting=True, verbose=False):
  if isinstance(spec1,dict) and isinstance(spec2,dict) and filt_dict:
    bands, w1, f1, e1, f2, e2, weight, bnds = [i for i in filt_dict.keys() if all([i in spec1.keys(),i in spec2.keys()]) and i not in exclude], [], [], [], [], [], [], []
    for eff,b in sorted([(filt_dict[i]['eff'],i) for i in bands]):
      if spec1[b] and spec1[b+'_unc'] and spec2[b]: bnds.append(b), w1.append(eff), f1.append(spec1[b]), e1.append(spec1[b+'_unc']), f2.append(spec2[b]), e2.append(spec2[b+'_unc'] if b+'_unc' in spec2.keys() else 0*spec2[b].unit), weight.append((filt_dict[b]['max']-filt_dict[b]['min']) if weighting else 1)
    bands, w1, f1, e1, f2, e2, weight = map(np.array, [bnds, w, f, sig, F, Sig, weight])
    if verbose: printer(['Band','W_spec1','F_spec1','E_spec1','F_spec2','E_spec2','Weight','g-factor'],zip(*[bnds, w1, f1, e1, f2, e2, weight, weight*(f1-f2*(sum(weight*f1*f2/(e1**2 + e2**2))/sum(weight*f2**2/(e1**2 + e2**2))))**2/(e1**2 + e2**2)]))
  else:
    spec1, spec2 = [[i.value if hasattr(i,'unit') else i for i in j] for j in [spec1,spec2]]
    if exclude: spec1 = [i[idx_exclude(spec1[0],exclude)] for i in spec1]
    # (w1, f1), e1, f2, e2 = spec1[:2], spec1[2] if len(spec1)==3 else 1, np.interp(spec1[0], spec2[0], spec2[1], left=0, right=0), np.interp(spec1[0], spec2[0], spec2[2], left=0, right=0) if len(spec2)==3 else 0
    (w1, f1), e1, f2, e2 = spec1[:2], spec1[2] if len(spec1)==3 else 1, rebin_spec(spec2[0], spec2[1], spec1[0]), rebin_spec(spec2[0], spec2[2], spec1[0]) if len(spec2)==3 else 0
    weight = np.array([(w1[n+1]-w1[n])/2. if n==0 else (w1[n]-w1[n-1])/2. if n==len(w1)-1 else (w1[n+1]-w1[n-1])/2. for n,i in enumerate(w1)])
    if exclude: weight[weight>np.std(weight)] = 0
  C = sum(weight*f1*f2/(e1**2 + e2**2))/sum(weight*f2**2/(e1**2 + e2**2))
  G = weight*(f1-f2*C)**2/(e1**2 + e2**2)
  if verbose: plt.loglog(spec2[0], spec2[1]*C, 'r', label='spec2', alpha=0.6), plt.loglog(w1, f1, 'k', label='spec1', alpha=0.6), plt.loglog(w1, f2*C, 'b', label='spec2 binned', alpha=0.6), plt.grid(True), plt.legend(loc=0)
  return [G if array else sum(G), C]

def idx_include(x, include):
  try: return np.where(np.array(map(bool,map(sum, zip(*[np.logical_and(x>i[0],x<i[1]) for i in include])))))[0]
  except TypeError:
    try: return np.where(np.array(map(bool,map(sum, zip(*[np.logical_and(x>i[0],x<i[1]) for i in [include]])))))[0] 
    except TypeError: return range(len(x))

def idx_exclude(x, exclude):
  try: return np.where(~np.array(map(bool,map(sum, zip(*[np.logical_and(x>i[0],x<i[1]) for i in exclude])))))[0]
  except TypeError: 
    try: return np.where(~np.array(map(bool,map(sum, zip(*[np.logical_and(x>i[0],x<i[1]) for i in exclude])))))[0]
    except TypeError: return range(len(x))

def mag2flux(band, mag, unc=None, Flam=False, photon=False):
  '''
  For given band and magnitude, returns the flux value in [ergs][s-1][cm-2][cm-1]
  Note: Must be multiplied by wavelength in [cm] to be converted to [ergs][s-1][cm-2], not done here! 
  mag = -2.5*log10(F/zp)  =>  flux = zp*10**(-mag/2.5)
  '''
  filt = a.filter_info(band) 
  zp = filt['zp_photon' if photon else 'zp']*(1 if photon else q.erg)/q.s/q.cm**2/q.AA
  F = (zp*(filt['eff']*q.um if Flam else 1)*10**(-mag/2.5)).to((1 if photon else q.erg)/q.s/q.cm**2/(1 if Flam else q.AA))
  E = F - (zp*(filt['eff']*q.um if Flam else 1)*10**(-(mag+unc)/2.5)).to((1 if photon else q.erg)/q.s/q.cm**2/(1 if Flam else q.AA)) if unc else 1
  return [F,E] if unc else F 
  
def marginalized_distribution(data, figure='', xunits='', yunits='', xy='', color='b'):
  if figure: fig, ax1, ax2, ax3 = figure
  else:
    fig = plt.figure()
    ax1 = plt.subplot2grid((4,4), (1,0), colspan=3, rowspan=3)
    ax2, ax3 = plt.subplot2grid((4,4), (0,0), colspan=3, sharex=ax1), plt.subplot2grid((4,4), (1,3), rowspan=3, sharey=ax1)
  if xy: ax1.plot(*xy, c=color, marker='s', markersize=12)
  ax2.hist(data[0], histtype='stepfilled', align='left', alpha=0.5, color=color), ax3.hist(data[1], histtype='stepfilled', orientation='horizontal', align='left', alpha=0.5, color=color), ax1.scatter(*data, color=color), ax1.grid(True), ax2.grid(True), ax3.grid(True), ax2.xaxis.tick_top(), ax3.yaxis.tick_right()#, ax1.set_xlim(400,3000), ax1.set_ylim(2.5,6.0)
  fig.subplots_adjust(hspace=0, wspace=0), fig.canvas.draw()
  return [fig, ax1, ax2, ax3]

def modelFit(SED, phot_dict, spec_dict, dist='', filt_dict=None, exclude=[], plot=False, Rlim=(0,100), title='', weighting=True):
  '''
  For given *spectrum* [W,F,E] or dictionary of photometry, returns the best fit synthetic spectrum by varying surface gravity and effective temperature.
  '''
  fit_list, unfit_list, phot_fit, dof = [], [], isinstance(SED,dict), 0
  model_dict = phot_dict if phot_fit else spec_dict
  if phot_fit:
    for b in SED.keys():
      if 'unc' not in b:
        dof += 1
        if not SED[b] or not SED[b+'_unc']:
          SED.pop(b), SED.pop(b+'_unc')
          dof -= 1
    if all([b in SED for b in ['W3','W3_unc','W4','W4_unc']]):
      # If infrared excess, drop W4 since the models won't fit
      if SED['W4']>SED['W3']:
        SED.pop('W4'), SED.pop('W4_unc')
        dof -= 1
  for k in model_dict.keys():
    # W, F = spec_dict[k]['wavelength'], spec_dict[k]['flux']
    # good, const = goodness(SED, phot_dict[k], filt_dict=filt_dict, weighting=weighting) if phot_fit else goodness(SED, [W,F], exclude=exclude, weighting=weighting)
    # R = (dist*np.sqrt(float(const))/ac.R_jup).decompose().value
    # if R>Rlim[0] and R<Rlim[1]: fit_list.append((abs(good), k, float(const), phot_dict[k] if phot_fit else [W,F]))
    # else: unfit_list.append((abs(good), k, float(const), phot_dict[k] if phot_fit else [W,F]))
    
    W, F = spec_dict[k]['wavelength'], spec_dict[k]['flux']
    good, const = goodness(SED, phot_dict[k], filt_dict=filt_dict, weighting=weighting) if phot_fit else goodness(SED, [W,F], exclude=exclude, weighting=weighting)
    R = (dist*np.sqrt(float(const))/ac.R_jup).decompose().value
    if R>Rlim[0] and R<Rlim[1]: fit_list.append((abs(good), k, float(const), phot_dict[k] if phot_fit else [W,F]))
    else: unfit_list.append((abs(good), k, float(const), phot_dict[k] if phot_fit else [W,F]))

  top5 = nsmallest(5,fit_list)
  printer(['Goodness','Parameters','Radius' if dist else '(R/d)**2'], [[i[0], i[1], (dist*np.sqrt(i[2])/ac.R_jup).decompose()] for i in top5] if dist else top5)
  
  if plot:
    from itertools import groupby
    fig, to_print = plt.figure(), []
    for (ni,li) in [('fl',fit_list),('ul',unfit_list)]:
      for key,group in [[k,list(grp)] for k,grp in groupby(sorted(li, key=lambda x: x[1].split()[1]), lambda x: x[1].split()[1])]:
        g, p, c, d = zip(*sorted(group, key=lambda x: int(x[1].split()[0])))
        # plt.plot([int(t.split()[0]) for t in p], g, 'o' if ni=='fl' else 'x', ls='-' if ni=='fl' else 'none', color=plt.cm.jet((5.6-float(key))/2.4,1), label=key)

        plt.plot([int(t.split()[0]) for t in p], g, 'o' if ni=='fl' else 'x', ls='-' if ni=='fl' else 'none', color=plt.cm.jet((5.6-float(key))/2.4,1), label=key)

        # to_print.append([[i,j] for i,j in zip(p,g)])
    # if dof: 
      # import scipy.stats as st
      # plt.fill_between(plt.xlim(), [plt.ylim()[0],plt.ylim()[0]], [top5[0][0]+st.chi2.ppf(st.chi2.cdf(1,1),dof),top5[0][0]+st.chi2.ppf(st.chi2.cdf(1,1),dof)], color='k', alpha=0.1)
    # printer(['Params','G-factor'],[val for subl in to_print for val in subl], to_txt='/Users/Joe/Desktop/{} - {} fit.txt'.format(plot,'phot' if phot_fit else 'spec'))
    plt.legend(loc=0, ncol=2), plt.grid(True), plt.ylabel('Goodness of Fit'), plt.xlabel('Teff'), plt.suptitle(plot)
      
  P, C, D = min(fit_list)[1:]   
  res = np.average(np.diff(SED[0]))
  synW, synF = spec_dict[P]['wavelength'], spec_dict[P]['flux']
  # synF = rebin_spec(synW, synF, np.concatenate([  synW[0][np.where(synW[0]<SED[0][0])] , SED[0], synW[0][np.where(synW[0]>SED[0][-1])] ]))
  synF = rebin_spec(synW, synF, SED[0])
  return [[SED[0], synF*C], P, C]

def modelInterp(params, model_dict, filt_dict=None, plot=False):
  '''
  Returns the interpolated model atmosphere spectrum (if model_dict==spec_dict) or photometry (if model_dict==phot_dict and filt_dict provided)
  '''
  t, g = int(params.split()[0]), params.split()[1]
  p1, p2 = sorted(zip(*nsmallest(2,[[abs(int(k.split()[0])-t),k] for k in model_dict.keys() if g in k]))[1])
  t1, t2 = int(p1.split()[0]), int(p2.split()[0])
  if filt_dict:
    D = {}
    for i in filt_dict.keys():
      try: 
        D[i] = model_dict[p2][i]+(model_dict[p1][i]-model_dict[p2][i])*(t**4-t2**4)/(t1**4-t2**4)
        if plot: plt.loglog(filt_dict[i]['eff'], model_dict[p1][i], 'bo', label=p1 if i=='J' else None, alpha=0.7), plt.loglog(filt_dict[i]['eff'], D[i], 'ro', label=params if i=='J' else None, alpha=0.7), plt.loglog(filt_dict[i]['eff'], model_dict[p2][i], 'go', label=p2 if i=='J' else None, alpha=0.7)
      except KeyError: pass
    if plot: plt.legend(loc=0), plt.grid(True)
    return D
  else:
    w1, f1, w2, f2 = model_dict[p1]['wavelength'], model_dict[p1]['flux'], model_dict[p2]['wavelength'], model_dict[p2]['flux']
    if len(f1)!=len(f2): f2 = np.interp(w1, w2, f2, left=0, right=0)
    F = f2+(f1-f2)*(t**4-t2**4)/(t1**4-t2**4)
    if plot: plt.loglog(w1, f1, '-b', label=p1, alpha=0.7), plt.loglog(w1, F, '-r', label=params, alpha=0.7), plt.loglog(w1, f2, '-g', label=p2, alpha=0.7), plt.legend(loc=0), plt.grid(True)
    return [w1,F]
  
def modelReplace(spectrum, model, replace=[], Flam=False, tails=False, plot=False):
  '''
  Returns the given *spectrum* with the tuple ranges in *replace* replaced by the given *model*.
  '''
  if tails: replace += [(0.1,spectrum[0][0].value),(spectrum[0][-1].value,100)]
  fW, fF, fE = [i[idx_include(model[0].value,replace)] for i in [model[0].value, model[1].value, model[1].value]]
  W, F, E = [i[idx_exclude(spectrum[0].value,replace)] for i in spectrum]
  W, F, E = map(np.array,zip(*sorted(zip(*[np.concatenate(i) for i in [[W.value,fW],[F.value,fF],[E.value,fE]]]), key=lambda x: x[0])))
  if plot: plt.figure(), plt.loglog(*model, color='k', alpha=0.3), plt.loglog(*spectrum[:2], color='b'), plt.loglog(W, F, color='k', ls='--'), plt.legend(loc=0)
  return [i*j for i,j in zip([W,F,E],[k.unit for k in spectrum])]
  
def norm_spec(spectrum, template, exclude=[], include=[]):
  '''
  Returns *spectrum* with [W,F] or [W,F,E] normalized to *template* [W,F] or [W,F,E].
  Wavelength range tuples provided in *exclude* argument are ignored during normalization, i.e. exclude=[(0.65,0.72),(0.92,0.97)].
  '''                                                          
  S, T = scrub([i.value if hasattr(i,'_unit') else i for i in spectrum]), scrub([i.value if hasattr(i,'_unit') else i for i in template])
  S0, T0 = [i[idx_include(S[0],[(T[0][0],T[0][-1])])] for i in S], [i[idx_include(T[0],[(S[0][0],S[0][-1])])] for i in T]
  if exclude: S0, T0 = [[i[idx_exclude(j[0],exclude)] for i in j] for j in [S0,T0]]
  if include: S0, T0 = [[i[idx_include(j[0],include)] for i in j] for j in [S0,T0]]
  try: norm = np.trapz(T0[1], x=T0[0])/np.trapz(np.interp(T0[0],*S0[:2]), x=T0[0])
  except ValueError: norm = 1            
  S[1] = S[1]*norm                                                                              
  try: S[2] = S[2]*norm                                                        
  except IndexError: pass
  return S

def normalize(spectra, template, composite=True, plot=False, SNR=100, exclude=[], trim=[], replace=[], D_Flam=None):
  '''
  Normalizes a list of *spectra* with [W,F,E] or [W,F] to a *template* spectrum.
  Returns one normalized, composite spectrum if *composite*, else returns the list of *spectra* normalized to the *template*.
  '''    
  if not template: 
    spectra = [unc(i, SNR=SNR) for i in sorted(spectra, key=lambda x: x[1][-1])]
    template = spectra.pop()
                                                                                    
  if trim:
    all_spec = [template]+spectra
    for n,x1,x2 in trim: all_spec[n] = [i[idx_exclude(all_spec[n][0],[(x1,x2)])] for i in all_spec[n]]
    template, spectra = all_spec[0], all_spec[1:]
  
  (W, F, E), normalized = template, []
  for S in spectra: normalized.append(norm_spec([i.value if hasattr(i,'_unit') else i for i in unc(S, SNR=SNR)], [W,F,E], exclude=exclude+replace))
  if plot: plt.loglog(W, F, alpha=0.5), plt.fill_between(W, F-E, F+E, alpha=0.1)
    
  if composite:
    for n,(w,f,e) in enumerate(normalized):
      tries = 0
      while tries<3:
        IDX, idx = np.where(np.logical_and(W<w[-1],W>w[0]))[0], np.where(np.logical_and(w>W[0],w<W[-1]))[0]
        if not any(IDX):
          normalized.pop(n), normalized.append([w,f,e])
          tries += 1
        else:
          (W0, F0, E0), (w0, f0, e0) = [i[IDX] for i in [W,F,E]], [i[idx] for i in [w,f,e]]
          f0, e0 = np.interp(W0, w0, f0), np.interp(W0, w0, e0)
          if exclude:
            Eidx = idx_include(W0,exclude)
            keep, E0[Eidx] = E0[Eidx], 1E-30
          f_mean = np.array([np.average([fl,FL], weights=[1/er,1/ER]) for fl,er,FL,ER in zip(f0,e0,F0,E0)])
          if exclude: E0[Eidx] = keep
          e_mean = np.sqrt(e0**2 + E0**2)
          spec1, spec2 = min([W,F,E], [w,f,e], key=lambda x: x[0][0]), max([W,F,E], [w,f,e], key=lambda x: x[0][-1])
          spec1, spec2 = [i[np.where(spec1[0]<W0[0])[0]] for i in spec1], [i[np.where(spec2[0]>W0[-1])[0]] for i in spec2]
          W, F, E = [np.concatenate([i,j,k]) for i,j,k in zip(spec1,[W0,f_mean,e_mean],spec2)]
          tries = 5

    if replace: W, F, E = modelReplace([W,F,E], replace=replace, D_Flam=D_Flam)

  if plot:
    for w,f,e in normalized: plt.loglog(w, f, alpha=0.5), plt.fill_between(w, f-e, f+e, alpha=0.2)
    if composite: plt.loglog(W, F, '--', c='k', lw=1), plt.fill_between(W, F-E, F+E, color='k', alpha=0.2)
    plt.yscale('log', nonposy='clip')

  if not composite: normalized.insert(0, unc(scrub(template), SNR=SNR))
  else: normalized = [[W,F,E]]
  return normalized[0][:len(template)] if composite else [i[:len(template)] for i in normalized]

def pi2pc(parallax): 
  if isinstance(parallax,(tuple,list)):
    pi, sig_pi = parallax[0]*q.arcsec/1000., parallax[1]*q.arcsec/1000.
    d, sig_d = (1*q.pc*q.arcsec)/pi, sig_pi*q.pc*q.arcsec/pi**2
    return (d, sig_d)
  else: return (1*q.pc*q.arcsec)/(parallax*q.arcsec/1000.)

def printer(labels, values, format='', truncate=150, to_txt=None, highlight=[]):
  '''
  Prints a nice table of *values* with *labels* with auto widths else maximum width if *same* else *col_len* if specified. 
  '''
  def red(t): print "\033[01;31m{0}\033[00m".format(t),
  print '\r'
  values = [["-" if not i else "{:.6g}".format(i) if isinstance(i,(float,int)) else i if isinstance(i,(str,unicode)) else "{:.6g} {}".format(i.value,i.unit) if hasattr(i,'unit') else i for i in j] for j in values]
  auto, txtFile = [max([len(i) for i in j])+2 for j in zip(labels,*values)], open(to_txt, 'a') if to_txt else None
  lengths = format if isinstance(format,list) else [min(truncate,i) for i in auto]
  col_len = [max(auto) for i in lengths] if format=='max' else [150/len(labels) for i in lengths] if format=='fill' else lengths
  for l,m in zip(labels,col_len):
    print str(l)[:truncate].ljust(m),
    if to_txt: txtFile.write(str(l)[:truncate].replace(' ','').ljust(m))
  for row_num,v in enumerate(values):
    print '\n',
    if to_txt: txtFile.write('\n') 
    for col_num,(k,j) in enumerate(zip(v,col_len)):
      if (row_num,col_num) in highlight: red(str(k)[:truncate].ljust(j))
      else: print str(k)[:truncate].ljust(j),
      if to_txt: txtFile.write(str(k)[:truncate].replace(' ','').ljust(j))
  print '\n'

def rebin_spec(wave, specin, wavnew, waveunits='um'):
    spec = spectrum.ArraySourceSpectrum(wave=wave, flux=specin)
    f = np.ones(len(wave))
    filt = spectrum.ArraySpectralElement(wave, f, waveunits=waveunits)
    obs = observation.Observation(spec, filt, binset=wavnew, force='taper')
    return obs.binflux

def rgb_image(images, save=''):
  '''
  Saves an RGB false color image at *save* made from a stack of three *images*
  From the APLpy (Apple Pie) module (http://aplpy.readthedocs.org/en/latest/howto_rgb.html)
  '''
  import aplpy
  aplpy.make_rgb_image(images,save)
  
def separation(ra1, dec1, ra2, dec2):
  '''
  Given coordinates *ra1*, *dec1*, *ra2*, *dec2* of two objects, returns the angular separation in arcseconds.
  '''
  if isinstance(ra1,str): ra1 = float(ra1) if ra1.isdigit() else sxg2deg(ra=ra1)
  if isinstance(dec1,str): dec1 = float(dec1) if dec1.isdigit() else sxg2deg(dec=dec1)
  if isinstance(ra2,str): ra2 = float(ra2) if ra2.isdigit() else sxg2deg(ra=ra2)
  if isinstance(dec2,str): dec2 = float(dec2) if dec2.isdigit() else sxg2deg(dec=dec2) 

  try: return (float(apc.angles.AngularSeparation(ra1, dec1, ra2, dec2, q.degree).format(decimal=True,unit='degree'))*q.degree).to(q.arcsec).value
  except TypeError: return None

def sameName(name1, name2, chars=4):
  '''
  Boolean: Given names of two objects, checks that they have a certain number of name characters in common 
  Note: discounts '2' in object names with '2MASS' or '2m'

  *chars*
   Number of consecutive characters to match, 4 by default. (e.g '2m0355' and '2MASSJ0355+1234' have '0355' in common with chars=4)
  '''
  import re
  def clean(name):
    for i in ['2MASS', '2mass', '2M', '2m']:
      try: name = name.replace(i,'')
      except AttributeError: pass
    return name

  n1, n2 = re.sub('\D','',clean(str(name1))), re.sub('\D','',clean(str(name2)))  
  return True if re.sub('\D','',clean(str(name1)))[:chars] == re.sub('\D','',clean(str(name2)))[:chars] else False

def scrub(data):
  '''
  For input data [w,f,e] or [w,f] returns the list with NaN, negative, and zero flux (and corresponsing wavelengths and errors) removed. 
  '''
  data = [i for i in data if isinstance(i,(np.ndarray,q.Quantity))]
  return [i[np.where((data[1].value>0) & (~np.isnan(data[1].value)))] if hasattr(i,'_unit') else i[np.where((data[1]>0) & (~np.isnan(data[1])))] for i in data]

def smooth(x,beta):
  """
  Smooths a spectrum *x* using a Kaiser-Bessel smoothing window of narrowness *beta* (~1 => very smooth, ~100 => not smooth) 
  """
  window_len = 11
  s = np.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]]
  w = np.kaiser(window_len,beta)
  y = np.convolve(w/w.sum(), s, mode='valid')
  return y[5:len(y)-5]*(x.unit if hasattr(x, 'unit') else 1)
  
def str2Q(x,target=''):
  '''
  Given a string of units unconnected to a number, returns the units as a quantity to be multiplied with the number. 
  Inverse units must be represented by a forward-slash prefix or negative power suffix, e.g. inverse square seconds may be "/s2" or "s-2" 

  *x*
    The units as a string, e.g. str2Q('W/m2/um') => np.array(1.0) * W/(m**2*um)
  *target*
    The target units as a string if rescaling is necessary, e.g. str2Q('Wm-2um-1',target='erg/s/cm2/cm') => np.array(10000000.0) * erg/(cm**3*s)
  '''
  if x:       
    def Q(IN):
      OUT = 1
      text = ['Jy', 'erg', '/s', 's-1', 's', '/um', 'um-1', 'um', '/cm2', 'cm-2', 'cm2', '/cm', 'cm-1', 'cm', '/A', 'A-1', 'A', 'W', '/m2', 'm-2', 'm2', '/m', 'm-1', 'm', '/Hz', 'Hz-1']
      vals = [q.Jy, q.erg, q.s**-1, q.s**-1, q.s, q.um**-1, q.um**-1, q.um, q.cm**-2, q.cm**-2, q.cm**2, q.cm**-1, q.cm**-1, q.cm, q.AA**-1, q.AA**-1, q.AA, q.W, q.m**-2, q.m**-2, q.m**2, q.m**-1, q.m**-1, q.m, q.Hz**-1, q.Hz**-1]
      for t,v in zip(text,vals):
        if t in IN:
          OUT = OUT*v
          IN = IN.replace(t,'')
      return OUT

    unit = Q(x)
    if target:
      z = str(Q(target)).split()[-1]
      try:
        unit = unit.to(z)
      except ValueError:
        print "{} could not be rescaled to {}".format(unit,z)

    return unit 
  else:
    return 1 
      
def squaredError(a, b, c):
  '''
  Computes the squared error of two arrays. Pass to scipy.optimize.fmin() to find least square or use scipy.optimize.leastsq()
  '''
  a -= b
  a *= a 
  c = np.array([1 if np.isnan(e) else e for e in c])
  return sum(a/c)

def sxg2deg(ra='', dec=''):
  RA, DEC = '', ''
  if ra: RA = float(apc.angles.Angle(ra, unit='hour').format(decimal=True, precision=8))
  if dec: DEC = float(apc.angles.Angle(dec, unit='degree').format(decimal=True, precision=8))
  return (RA, DEC) if ra and dec else RA or DEC

def tails(spectrum, model, plot=False):
  '''
  Appends the Wein and Rayleigh-Jeans tails of the *model* to the given *spectrum*
  '''
  start, end = np.where(model[0]<spectrum[0][0])[0], np.where(model[0]>spectrum[0][-1])[0]
  final = [np.concatenate(i) for i in [[model[0][start],spectrum[0],model[0][end]], [model[1][start],spectrum[1],model[1][end]], [np.zeros(len(start)),spectrum[2],np.zeros(len(end))]]]  
  if plot: plt.loglog(*spectrum[:2]), plt.loglog(model[0],model[1]), plt.loglog(*final[:2], color='k', ls='--')
  return final

def txt2dict(txtfile, delim='', skip=[], ignore=[], to_list=False, all_str=False, obj_col=0, key_row=0, start=1):
  '''
  For given *txtfile* returns a parent dictionary with keys from *obj_col* and child dictionaries with keys from *key_row*, delimited by *delim* character.
  Characters to *ignore*, entire lines to *skip*, and the data row to *start* at can all be specified.
  Floats and integers are returned as numbers unless *all_str* is set True.
  '''
  def replace_all(text, dic):
    for i in dic: text = text.replace(i,' ')
    return text
    
  txt = open(txtfile)
  d = filter(None,[[j.strip() for j in replace_all(i,ignore).split(delim or None)] for i in txt if not any([i.startswith(c) for c in skip])])
  txt.close()
  
  for i in d: i.insert(0,i.pop(obj_col))
  keys = d[key_row][1:]

  if all_str: return d if to_list else {row[0]:{k:str(v).replace('\"','').replace('\'','') for k,v in zip(keys,row[1:])} for row in d[start:]}
  else: return d if to_list else {row[0]:{k:float(v) if v.replace('-','').replace('.','').isdigit() and '-' not in v[1:] else str(v).replace('\"','').replace('\'','') if isinstance(v,unicode) else True if v=='True' else False if v=='False' else v.replace('\"','').replace('\'','') for k,v in zip(keys,row[1:])} for row in d[start:]}

def try_except(success, failure, exceptions):
  '''
  Replaces the multi-line try/except statement with a function
  '''
  try:
    return success() if callable(success) else success
  except exceptions or Exception:
    return failure() if callable(failure) else failure      

def unc(spectrum, SNR=100):
  '''
  Removes NaNs negatives and zeroes from *spectrum* arrays of form [W,F] or [W,F,E].
  Generates E at signal to noise *SNR* for [W,F] and replaces NaNs with the same for [W,F,E]. 
  '''
  S = scrub(spectrum)
  # if len(S)==3: S[2][np.where(np.isnan(S[2]))] = S[1][np.where(np.isnan(S[2]))]/SNR
  if len(S)==3:
    try: S[2] = np.array([i/SNR if np.isnan(j) else j for i,j in zip(S[1],S[2])], dtype='float32')
    except TypeError: S[2] = np.array(S[1]/SNR)
  elif len(S)==2: S.append(np.array(S[1]/SNR))
  return S

def xl2dict(filepath, sheet=1, obj_col=0, key_row=0, start=1, manual_keys=''):
  workbook = xlrd.open_workbook(filepath)
  column_names = manual_keys or [str(i) for i in workbook.sheet_by_index(sheet).row_values(key_row)]
  objects = workbook.sheet_by_index(sheet).col_values(obj_col)[start:]
  if manual_keys: values = [workbook.sheet_by_index(sheet).col_values(n)[start:] for n in range(len(manual_keys))]
  else: values = [workbook.sheet_by_index(sheet).col_values(c)[start:] for c in [column_names.index(i) for i in column_names]]
  return {str(obj): {str(cn):str(val.encode('utf-8')) if isinstance(val,unicode) else val for cn,val in zip(column_names,value)} for obj,value in zip(objects,zip(*values))}
