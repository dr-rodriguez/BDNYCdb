#!/usr/bin/python
# Utilities
<<<<<<< HEAD
import warnings, re, cPickle, numpy as np, matplotlib.pyplot as plt, astropy as ap
=======
import warnings, re, cPickle, astropy.units as q, numpy as np, matplotlib.pyplot as plt
>>>>>>> 7fd44e7c8b020a7c710a3cf3fca8ee420757e6cf
from random import random
warnings.simplefilter('ignore')
path = '/Users/Joe/Documents/Python/'

def app2abs(mag, dist, app=True): return (mag-(1 if app else -1)*5*np.log10(dist/(10*q.pc))).value.item()                                                      

def ChiSquare(a, b, unc=None, array=False, Gtest=False):
  a, b = [np.array(map(float,i.value)) if hasattr(i,'_unit') else np.array(map(float,i)) for i in [a,b]]
  c, variance = np.array(map(float,unc.value)) if hasattr(unc, '_unit') else np.array(map(float,np.ones(len(a)))), np.std(b)**4 # Since the standard deviation is the root of the variance
  X2 = np.array([(j*np.log(j/i)/k)**2/i for i,j,k in zip(a,b,c)]) if Gtest else np.array([((i-j)/k)**2/variance for i,j,k in zip(a,b,c)])    
  return X2 if array else sum(X2)
  
def dict2txt(DICT, writefile, column1='-', delim='\t', digits=5, order=''):
  '''
  Given a nested dictionary *DICT*, writes a .txt file with keys as columns. 
  '''
  import csv
  D = DICT.copy()
  with open( writefile, 'w' ) as f:
    writer, w = csv.writer(f, delimiter=delim), []
    for k in D.keys():
      w.append(k)
      for i in D[k].keys():
        D[k][i] = '-' if not D[k][i] else round(D[k][i],digits) if isinstance(D[k][i],(float,int)) else str(D[k][i]) if isinstance(D[k][i],unicode) else D[k][i]
        w.append(i), w.append(str(D[k][i]))
    width = len(max(w, key=len))
    head = ['{!s:{}}'.format(column1,width)]
    headorder = order or D[D.keys()[0]].keys()
    for i in headorder: head.append('{!s:{}}'.format(i,width))
    writer.writerow(head)
    for i in D.keys():
      order = order or sorted(D[i].keys())
      row = ['{!s:{}}'.format(i,width)]
      for k in order:
        if k not in D[i].keys(): D[i][k] = '-'
        row.append('{!s:{}}'.format(D[i][k],width))
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

def goodness(spectrum, model, radius=1, dist=1):
  (w, f, sig), (W, F) = spectrum, model
  return sum((f-(F*((radius/dist)**2)))/sig)

def idx_include(x, include):
  try: return np.where(np.array(map(bool,map(sum, zip(*[np.logical_and(x>i[0],x<i[1]) for i in include])))))[0]
  except TypeError: return range(len(x))

def idx_exclude(x, exclude):
  try: return np.where(~np.array(map(bool,map(sum, zip(*[np.logical_and(x>i[0],x<i[1]) for i in exclude])))))[0]
  except TypeError: return range(len(x))

def mag2flux(band, mag, unc=None, Flam=False, photon=False):
  '''
  For given band and magnitude, returns the flux value in [ergs][s-1][cm-2][cm-1] or [Jy] if Jy=True
  Note: Must be multiplied by wavelength in [cm] to be converted to [ergs][s-1][cm-2], not done here! 
  mag = -2.5*log10(F/zp)  =>  flux = zp*10**(-mag/2.5)
  '''
  filt = a.filter_info(band) 
  zp = filt['zp_photon' if photon else 'zp']*(1 if photon else q.erg)/q.s/q.cm**2/q.Angstrom
  F = (zp*(filt['eff']*q.um if Flam else 1)*10**(-mag/2.5)).to((1 if photon else q.erg)/q.s/q.cm**2/(1 if Flam else q.Angstrom))
  E = F - (zp*(filt['eff']*q.um if Flam else 1)*10**(-(mag+unc)/2.5)).to((1 if photon else q.erg)/q.s/q.cm**2/(1 if Flam else q.Angstrom)) if unc else 1
  return [F,E] if unc else F

def modelFit(spectrum, exclude=[], Flam=True, SNR=50, D_Flam=None):
  '''
  For given *spectrum* [W,F,E] returns the best fit synthetic spectrum by varying surface gravity and effective temperature.
  '''
  S, spec_list = D_Flam or cPickle.load(open('/Users/Joe/Documents/Python/Pickles/synSpec_Flam_3000.p',"rb")), []
  
  if isinstance(spectrum,dict):
    from syn_phot.syn_phot import get_filters, color_table
    filters, CT = get_filters(), color_table(photon=True)
    
    for p in CT:
      wavs, modelMags, specMags, bands = [], [], [], []
      for b in [b for b in list(set(spectrum.keys()).intersection(CT[p].keys())) if all([spectrum[b],CT[p][b]])]:
        wavs.append(filters[b]['eff']), modelMags.append(CT[p][b]), specMags.append(spectrum[b]), bands.append(b)
      (wavs, modelMags, specMags, bands), norm = map(np.array,zip(*sorted(zip(wavs,modelMags,specMags,bands)))), spectrum['J']/CT[p]['J']
      spec_list.append((ChiSquare(modelMags*norm, specMags), p, norm))

  else:    
    spec = [i.value if hasattr(i,'_unit') else i for i in unc(spectrum, SNR=SNR)]
    wave, flux, error = [i[idx_exclude(spec[0],exclude)] for i in spec] if exclude else spec

    for k in S.keys():
      model = np.interp(wave, S[k]['W'], S[k]['F'], left=0, right=0)
      norm = np.trapz(flux)/np.trapz(model)
      spec_list.append((ChiSquare(model*norm, flux, unc=error), k, norm))

  from heapq import nsmallest
  X2, params, norm = min(spec_list)
  printer(['Chi-square','Parameters'], nsmallest(5,spec_list))
  return [S[params]['W'], (S[params]['F']*norm*(S[params]['W'] if Flam else 1)).to(q.erg/q.s/q.cm**2/(1 if Flam else q.Angstrom)), params]

def norm_spec(spectrum, template, exclude=[]):
  '''
  Returns *spectrum* with [W,F] or [W,F,E] normalized to *template* [W,F] or [W,F,E].
  Wavelength range tuples provided in *exclude* argument are ignored during normalization, i.e. exclude=[(0.65,0.72),(0.92,0.97)].
  '''                                                          
  S, T = scrub(spectrum), scrub(template)
  S0, T0 = [i[idx_include(S[0],[(T[0][0],T[0][-1])])] for i in S], [i[idx_include(T[0],[(S[0][0],S[0][-1])])] for i in T]
  if exclude: S0, T0 = [[i[idx_exclude(j[0],exclude)] for i in j] for j in [S0,T0]]
  norm = np.trapz(T0[1], x=T0[0])/np.trapz(np.interp(T0[0],*S0[:2])*(S[1].unit if hasattr(S[1],'_unit') else 1), x=T0[0])                 
  S[1] = S[1]*norm                                                                              
  try: S[2] = S[2]*norm                                                        
  except IndexError: pass
  return S

def normalize(spectra, template, composite=True, plot=False, SNR=100, exclude=[], trim=[], modelReplace=[], D_Flam=None):
  '''
  Normalizes a list of *spectra* with [W,F,E] or [W,F] to a *template* spectrum.
  Returns one normalized, composite spectrum if *composite*, else returns the list of *spectra* normalized to the *template*.
  '''    
  if not template: 
    spectra = sorted(spectra, key=lambda x: x[1][-1])
    template = spectra.pop()
                                                                                    
  if trim:
    all_spec = [template]+spectra
    for n,x1,x2 in trim: all_spec[n] = [i[idx_exclude(all_spec[n][0],[(x1,x2)])] for i in all_spec[n]]
    template, spectra = all_spec[0], all_spec[1:]
  
  (W, F, E), normalized = [i.value if hasattr(i,'_unit') else i for i in unc(template, SNR=SNR)], []
  for S in spectra: normalized.append(norm_spec([i.value if hasattr(i,'_unit') else i for i in unc(S, SNR=SNR)], [W,F,E], exclude=exclude+modelReplace))
  if plot: plt.loglog(W, F, alpha=0.5), plt.fill_between(W, F-E, F+E, alpha=0.1)
    
  if composite:
    for w,f,e in normalized:
      IDX, idx = np.where(np.logical_and(W<w[-1],W>w[0]))[0], np.where(np.logical_and(w>W[0],w<W[-1]))[0]
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

    if modelReplace:
      (fitW, fitF, fitP) = modelFit([W,F,E], exclude=exclude+modelReplace, Flam=False, SNR=SNR, D_Flam=D_Flam or cPickle.load(open('/Users/Joe/Documents/Python/Pickles/synSpec_Flam_3000.p',"rb")))
      if plot: plt.loglog(fitW[::20], fitF[::20], alpha=0.3)
      fitW, fitF = [i[idx_include(W,modelReplace)] for i in [W,np.interp(W, fitW, fitF, left=0, right=0)]]
      W, F = [i[idx_exclude(W,modelReplace)] for i in [W,F]]
      W, F, E = map(np.array,zip(*sorted(zip(*[np.concatenate(i) for i in [[W,fitW],[F,fitF],[E]]]), key=lambda x: x[0])))

  if plot:
    for w,f,e in normalized: plt.loglog(w, f, alpha=0.5), plt.fill_between(w, f-e, f+e, alpha=0.2)
    if composite: plt.loglog(W, F, '--', c='k', lw=1), plt.fill_between(W, F-E, F+E, color='k', alpha=0.2)
    plt.yscale('log', nonposy='clip')

  if not composite: normalized.insert(0, unc(scrub(template), SNR=SNR))
  else: normalized = [[W,F,E]]
  return normalized[0][:len(template)] if composite else [i[:len(template)] for i in normalized]

def pi2pc(pi): return (1*q.pc*q.arcsec)/(pi*q.arcsec/1000.)

def printer(labels, values, format='', to_txt=None):
  '''
  Prints a nice table of *values* with *labels* with auto widths else maximum width if *same* else *col_len* if specified. 
  '''
  print '\r'
  values = [["None" if not i else "{:.10g}".format(i) if isinstance(i,(float,int)) else i if isinstance(i,(str,unicode)) else "{:.10g} {}".format( float(i.value if hasattr(i,'_value') else i), str(i.unit.to_string() if hasattr(i,'_unit') else '').split()[1]) for i in j] for j in values]
  auto, txtFile = [max([len(i) for i in j])+2 for j in zip(labels,*values)], open(to_txt, 'a') if to_txt else None
  lengths = format if isinstance(format,list) else auto
  col_len = [max(auto) for i in lengths] if format=='max' else [150/len(labels) for i in lengths] if format=='fill' else lengths
  for l,m in zip(labels,col_len):
    print str(l).ljust(m),
    if to_txt: txtFile.write(str(l).replace(' ','').ljust(m))
  for v in values:
    print '\n',
    if to_txt: txtFile.write('\n') 
    for k,j in zip(v,col_len):
      print str(k).ljust(j),
      if to_txt: txtFile.write(str(k).replace(' ','').ljust(j))
  print '\n'

def rgb_image(images, save=''):
  '''
  Saves an RGB false color image at *save* made from a stack of three *images*
  From the APLpy (Apple Pie) module (http://aplpy.readthedocs.org/en/latest/howto_rgb.html)
  '''
  import aplpy
  aplpy.make_rgb_image(images,save)
  
def sameCoords(ra1, dec1, ra2, dec2, deg=0.001389):
  '''
  Boolean: Given coordinates of two objects, checks that they are within a certain number of degrees of each other (i.e. the same object)

  *deg*
   Number of degrees of tolerance (0.01667 degrees = 1', 0.001389 degrees = 5" by default)
  '''
  if isinstance(ra1,str): ra1 = float(ra1) if ra1.isdigit() else a.HMS2deg(ra=ra1)
  if isinstance(dec1,str): dec1 = float(dec1) if dec1.isdigit() else a.HMS2deg(dec=dec1)
  if isinstance(ra2,str): ra2 = float(ra2) if ra2.isdigit() else a.HMS2deg(ra=ra2)
  if isinstance(dec2,str): dec2 = float(dec2) if dec2.isdigit() else a.HMS2deg(dec=dec2) 

  try: return True if (ra2-deg <= ra1 <= ra2+deg) and (dec2-deg <= dec1 <= dec2+deg) else False
  except TypeError: return False

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
  return [i[np.where((data[1]>0) & (~np.isnan(data[1])))] for i in data]

def smooth(x,beta):
  """
  Smooths a spectrum *x* using a Kaiser-Bessel smoothing window of narrowness *beta* (~1 => very smooth, ~100 => not smooth) 
  """
  window_len = 11
  s = np.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]]
  w = np.kaiser(window_len,beta)
  y = np.convolve(w/w.sum(), s, mode='valid')
  return y[5:len(y)-5]

def specInterp(teff, logg, D, bin=25, Flam=False):
  greater = min([d for d in D.keys() if (D[d]['Teff'] > teff) and (D[d]['logg'] == logg)])
  lesser = max([d for d in D.keys() if (D[d]['Teff'] < teff) and (D[d]['logg'] == logg)])
  Tg, Tl, W, Wl, Fg, Fl = D[greater]['Teff'], D[lesser]['Teff'], D[greater]['W'], D[lesser]['W'], D[greater]['F'], D[lesser]['F']
  if len(W) != len(Wl): Fl = np.interp(W,Wl,Fl)*q.erg/q.s/q.cm**2 if Flam else np.interp(W,Wl,Fl)*q.erg/q.s/q.cm**3
  F = Fl + (Fg-Fl)*(float(teff)**4-float(Tl)**4)/(float(Tg)**4-float(Tl)**4)
  
  # print lesser, greater 
  # plt.plot(W[::bin],Fg[::bin],color='r',label=greater), plt.plot(W[::bin],Fl[::bin],color='b',label=lesser)
  # plt.plot(W[::bin],F[::bin],ls='--',color='k',lw=2,label='Weighted Avg')
  # plt.xscale('log'), plt.yscale('log'), plt.grid(True), plt.legend(loc=0) 
  
  return (W[::bin], F[::bin])

def str2Q(x,target=''):
  '''
  Given a string of units unconnected to a number, returns the units as a quantity to be multiplied with the number. 
  Inverse units must be represented by a forward-slash prefix or negative power suffix, e.g. inverse square seconds may be "/s2" or "s-2" 

  *u*
    The units as a string, e.g. str2Q('W/m2/um') => np.array(1.0) * W/(m**2*um)
  *target*
    The target units as a string if rescaling is necessary, e.g. str2Q('Wm-2um-1',target='erg/s/cm2/cm') => np.array(10000000.0) * erg/(cm**3*s)
  '''
  if x:       
    def Q(IN):
      OUT = 1
      text = ['erg', '/s', 's-1', 's', '/um', 'um-1', 'um', '/cm2', 'cm-2', 'cm2', '/cm', 'cm-1', 'cm', \
              '/A', 'A-1', 'A', 'W', '/m2', 'm-2', 'm2', '/m', 'm-1', 'm', '/Hz', 'Hz-1']
      vals = [q.erg, q.s**-1, q.s**-1, q.s, q.um**-1, q.um**-1, q.um, q.cm**-2, q.cm**-2, q.cm**2, q.cm**-1, q.cm**-1, q.cm, 
              q.Angstrom**-1, q.Angstrom**-1, q.Angstrom, q.W, q.m**-2, q.m**-2, q.m**2, q.m**-1, q.m**-1, q.m, q.Hz**-1, q.Hz**-1]
      for t,v in zip(text,vals):
        if t in IN:
          OUT = OUT*v
          IN = IN.replace(t,'')
      return OUT

    unit = Q(x)
    if target:
      q = str(Q(target)).split()[-1]
      try:
        unit = unit.to(q)
      except ValueError:
        print "{} could not be rescaled to {}".format(unit,q)

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
  from astropy import coordinates as apc
  RA, DEC = '', ''
  if ra: RA = float(apc.angles.Angle(ra, unit='hour').format(decimal=True, precision=8))
  if dec: DEC = float(apc.angles.Angle(dec, unit='degree').format(decimal=True, precision=8))
  return (RA, DEC) if ra and dec else RA or DEC

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
  if len(S)==3: S[2] = np.array([(i/SNR) if np.isnan(j) else j for i,j in zip(*S[1:])], dtype='float32')
  elif len(S)==2: S.append(np.array([(i/SNR) for i in S[1]], dtype='float32'))
  return S

def xl2dict(filepath, sheet=1):
  import re, xlrd
  workbook, D = xlrd.open_workbook(filepath), {}
  column_names = [str(i) for i in workbook.sheet_by_index(sheet).row_values(0)]
  col_idx = [column_names.index(i) for i in column_names]
  values = [workbook.sheet_by_index(sheet).col_values(c)[1:] for c in col_idx]
  for value in zip(*values): D[value[0]] = {cn:val for cn,val in zip(column_names[1:],value[1:])}
  return D