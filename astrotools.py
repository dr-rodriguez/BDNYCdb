# Slimmed down version of astrotools necessary to interact with the SQL database
import astropy.io.ascii as ad, matplotlib.pyplot as plt, numpy as np, astropy.io.fits as pf

def filter_info(band):
  '''
  (By Joe Filippazzo)
   
  Effective, min, and max wavelengths in [um] and zeropoint in [erg s-1 cm-2 A-1] and [photon s-1 cm-2 A-1] for SDSS, Bessel, 2MASS, IRAC and WISE filters IN THE VEGA SYSTEM. Values from SVO filter profile service.
  
  *band*
      Name of filter band (e.g. 'J' from 2MASS, 'W1' from WISE, etc.) or list of filter systems (e.g. ['SDSS','2MASS','WISE'])
  '''
  Filters = { "FUV":     { 'eff': 0.154226, 'min': 0.134032, 'max': 0.180643, 'zp': 6.486734e-09, 'zp_photon': 5.035932e+02, 'toVega':0,      'ext': 2.62,   'system': 'GALEX' },
              "NUV":     { 'eff': 0.227437, 'min': 0.169252, 'max': 0.300667, 'zp': 4.511628e-09, 'zp_photon': 5.165788e+02, 'toVega':0,      'ext': 2.94,   'system': 'GALEX' },

              "U":       { 'eff': 0.357065, 'min': 0.303125, 'max': 0.417368, 'zp': 3.656264e-09, 'zp_photon': 6.576522e+02, 'toVega':0.0915, 'ext': 1.56,   'system': 'Bessel' }, 
              "B":       { 'eff': 0.437812, 'min': 0.363333, 'max': 0.549706, 'zp': 6.286883e-09, 'zp_photon': 1.385995e+03, 'toVega':0.0069, 'ext': 1.31,   'system': 'Bessel' }, 
              "V":       { 'eff': 0.544579, 'min': 0.473333, 'max': 0.687500, 'zp': 3.571744e-09, 'zp_photon': 9.837109e+02, 'toVega':0,      'ext': 1.02,   'system': 'Bessel' },
              "R":       { 'eff': 0.641420, 'min': 0.550435, 'max': 0.883333, 'zp': 2.157178e-09, 'zp_photon': 6.971704e+02, 'toVega':0.0018, 'ext': 0.83,   'system': 'Bessel' },
              "I":       { 'eff': 0.797880, 'min': 0.704167, 'max': 0.916667, 'zp': 1.132454e-09, 'zp_photon': 4.549636e+02, 'toVega':-0.0014,'ext': 0.61,   'system': 'Bessel' },

              "u":       { 'eff': 0.3543,   'min': 0.304828, 'max': 0.402823, 'zp': 3.617963e-09, 'zp_photon': 6.546739e+02, 'toVega':0.91,   'ext': 1.58,   'system': 'SDSS' },  # AB to Vega transformations from Blanton et al. (2007)
              "g":       { 'eff': 0.4770,   'min': 0.378254, 'max': 0.554926, 'zp': 5.491077e-09, 'zp_photon': 1.282871e+03, 'toVega':-0.08,  'ext': 1.23,   'system': 'SDSS' },  # AB to Vega transformations from Blanton et al. (2007)
              "r":       { 'eff': 0.6231,   'min': 0.541534, 'max': 0.698914, 'zp': 2.528924e-09, 'zp_photon': 7.794385e+02, 'toVega':0.16,   'ext': 0.89,   'system': 'SDSS' },  # AB to Vega transformations from Blanton et al. (2007)
              "i":       { 'eff': 0.7625,   'min': 0.668947, 'max': 0.838945, 'zp': 1.409436e-09, 'zp_photon': 5.278550e+02, 'toVega':0.37,   'ext': 0.68,   'system': 'SDSS' },  # AB to Vega transformations from Blanton et al. (2007)
              "z":       { 'eff': 0.9134,   'min': 0.796044, 'max': 1.083325, 'zp': 8.501067e-10, 'zp_photon': 3.807540e+02, 'toVega':0.54,   'ext': 0.52,   'system': 'SDSS' },  # AB to Vega transformations from Blanton et al. (2007)
              
              "DES_u":   { 'eff': 0.3543,   'min': 0.304828, 'max': 0.402823, 'zp': 5.360165e-09, 'zp_photon': 1.038526e+03, 'toVega':0,   'ext': 0,   'system': 'DES' },
              "DES_g":   { 'eff': 0.4770,   'min': 0.378254, 'max': 0.554926, 'zp': 5.215897e-09, 'zp_photon': 1.243521e+03, 'toVega':0,   'ext': 0,   'system': 'DES' },
              "DES_r":   { 'eff': 0.6231,   'min': 0.541534, 'max': 0.698914, 'zp': 2.265389e-09, 'zp_photon': 7.234969e+02, 'toVega':0,   'ext': 0,   'system': 'DES' },
              "DES_i":   { 'eff': 0.7625,   'min': 0.668947, 'max': 0.838945, 'zp': 1.235064e-09, 'zp_photon': 4.820083e+02, 'toVega':0,   'ext': 0,   'system': 'DES' },
              "DES_z":   { 'eff': 0.9134,   'min': 0.796044, 'max': 1.083325, 'zp': 8.072777e-10, 'zp_photon': 3.712548e+02, 'toVega':0,   'ext': 0,   'system': 'DES' },
              "DES_Y":   { 'eff': 1.0289,   'min': 0.304828, 'max': 0.402823, 'zp': 6.596909e-10, 'zp_photon': 3.280450e+02, 'toVega':0,   'ext': 0,   'system': 'DES' },
              
              "J":       { 'eff': 1.2350,   'min': 1.080647, 'max': 1.406797, 'zp': 3.129e-10,    'zp_photon': 1.943482e+02, 'toVega':0,      'ext': 0.0166, 'system': '2MASS' }, # ZP from Cohen et al. (2003)
              "H":       { 'eff': 1.6620,   'min': 1.478738, 'max': 1.823102, 'zp': 1.133e-10,    'zp_photon': 9.437966e+01, 'toVega':0,      'ext': 0.0146, 'system': '2MASS' }, # ZP from Cohen et al. (2003)
              "Ks":      { 'eff': 2.1590,   'min': 1.954369, 'max': 2.355240, 'zp': 4.283e-11,    'zp_photon': 4.664740e+01, 'toVega':0,      'ext': 0.0710, 'system': '2MASS' }, # ZP from Cohen et al. (2003)

              "MKO_Y":   { 'eff': 1.02894,  'min': 0.9635,   'max': 1.1025,   'zp': 5.869238e-10, 'zp_photon': 3.033632e+02, 'toVega':0,      'ext': 0.41,   'system': 'MKO' },
              "MKO_J":   { 'eff': 1.250,    'min': 1.148995, 'max': 1.348332, 'zp': 3.01e-10,     'zp_photon': 1.899569e+02, 'toVega':0,      'ext': 0.30,   'system': 'MKO' },   # eff and ZP from Tokunaga & Vacca (2005)
              "MKO_H":   { 'eff': 1.644,    'min': 1.450318, 'max': 1.808855, 'zp': 1.18e-10,     'zp_photon': 9.761983e+01, 'toVega':0,      'ext': 0.20,   'system': 'MKO' },   # eff and ZP from Tokunaga & Vacca (2005)
              "MKO_K":   { 'eff': 2.198,    'min': 1.986393, 'max': 2.397097, 'zp': 4.00e-11,     'zp_photon': 4.488476e+01, 'toVega':0,      'ext': 0.12,   'system': 'MKO' },   # eff and ZP from Tokunaga & Vacca (2005)
              "MKO_L":   { 'eff': 3.754,    'min': 3.326622, 'max': 4.207764, 'zp': 5.31e-12,     'zp_photon': 1.016455e+01, 'toVega':0,      'ext': 0.06,   'system': 'MKO' },   # eff and ZP from Tokunaga & Vacca (2005)
              "MKO_M":   { 'eff': 4.702,    'min': 4.496502, 'max': 4.865044, 'zp': 2.22e-12,     'zp_photon': 5.305197e+00, 'toVega':0,      'ext': 0.05,   'system': 'MKO' },   # eff and ZP from Tokunaga & Vacca (2005)

              "DENIS_I": { 'eff': 0.78621,  'min': 0.7007,   'max': 0.9140,   'zp': 1.182102e-09, 'zp_photon': 4.681495e+02, 'toVega':0,      'ext': 0.63,   'system': 'DENIS' },
              "DENIS_J": { 'eff': 1.22106,  'min': 1.0508,   'max': 1.3980,   'zp': 3.190256e-10, 'zp_photon': 1.961698e+02, 'toVega':0,      'ext': 0.31,   'system': 'DENIS' },
              "DENIS_Ks":{ 'eff': 2.14650,  'min': 1.9474,   'max': 2.3979,   'zp': 4.341393e-11, 'zp_photon': 4.691482e+01, 'toVega':0,      'ext': 0.13,   'system': 'DENIS' },              

              "W1":      { 'eff': 3.3526,   'min': 2.754097, 'max': 3.872388, 'zp': 8.1787e-12,   'zp_photon': 1.375073e+01, 'toVega':0,      'ext': 0.07,   'system': 'WISE' }, # eff and ZP from Jarrett et al. (2011)
              "W2":      { 'eff': 4.6028,   'min': 3.963326, 'max': 5.341360, 'zp': 2.4150e-12,   'zp_photon': 5.586982e+00, 'toVega':0,      'ext': 0.05,   'system': 'WISE' }, # eff and ZP from Jarrett et al. (2011)
              "W3":      { 'eff': 11.5608,  'min': 7.443044, 'max': 17.26134, 'zp': 6.5151e-14,   'zp_photon': 3.567555e-01, 'toVega':0,      'ext': 0.06,   'system': 'WISE' }, # eff and ZP from Jarrett et al. (2011)
              "W4":      { 'eff': 22.0883,  'min': 19.52008, 'max': 27.91072, 'zp': 5.0901e-15,   'zp_photon': 5.510352e-02, 'toVega':0,      'ext': 0.02,   'system': 'WISE' }, # eff and ZP from Jarrett et al. (2011)

              "[3.6]":   { 'eff': 3.507511, 'min': 3.129624, 'max': 3.961436, 'zp': 6.755364e-12, 'zp_photon': 1.192810e+01, 'toVega':0,      'ext': 0.07,   'system': 'IRAC' },
              "[4.5]":   { 'eff': 4.436578, 'min': 3.917328, 'max': 5.056057, 'zp': 2.726866e-12, 'zp_photon': 6.090264e+00, 'toVega':0,      'ext': 0.05,   'system': 'IRAC' },
              "[5.8]":   { 'eff': 5.628102, 'min': 4.898277, 'max': 6.508894, 'zp': 1.077512e-12, 'zp_photon': 3.052866e+00, 'toVega':0,      'ext': 0.04,   'system': 'IRAC' },
              "[8]":     { 'eff': 7.589159, 'min': 6.299378, 'max': 9.587595, 'zp': 3.227052e-13, 'zp_photon': 1.232887e+00, 'toVega':0,      'ext': 0.03,   'system': 'IRAC' },              
              "[24]":    { 'eff': 23.20960, 'min': 19.88899, 'max': 30.93838, 'zp': 3.935507e-15, 'zp_photon': 4.598249e-02, 'toVega':0,      'ext': 0.02,   'system': 'MIPS' },

              "G":       { 'eff': 0.60,     'min': 0.321,    'max': 1.103,    'zp': 2.862966e-09, 'zp_photon': 8.053711e+02, 'toVega':0,      'ext': 0,      'system': 'GAIA'},
              "BP":      { 'eff': 0.55,     'min': 0.321,    'max': 0.680,    'zp': 4.298062e-09, 'zp_photon': 1.067265e+03, 'toVega':0,      'ext': 0,      'system': 'GAIA'},
              "RP":      { 'eff': 0.75,     'min': 0.627,    'max': 1.103,    'zp': 1.294828e-09, 'zp_photon': 4.948727e+02, 'toVega':0,      'ext': 0,      'system': 'GAIA'},

              "F090M":   { 'eff': 0.897360, 'min': 0.784317, 'max': 1.013298, 'zp': 8.395228e-10, 'zp_photon': 3.792477e+02, 'toVega':0,      'ext': 0.51,   'system': 'HST' },
              "F110W":   { 'eff': 1.059175, 'min': 0.782629, 'max': 1.432821, 'zp': 4.726040e-10, 'zp_photon': 2.519911e+02, 'toVega':0,      'ext': 0.39,   'system': 'HST' },
              "F140W":   { 'eff': 1.364531, 'min': 1.185379, 'max': 1.612909, 'zp': 2.133088e-10, 'zp_photon': 1.465263e+02, 'toVega':0,      'ext': 0.26,   'system': 'HST' },
              "F164N":   { 'eff': 1.646180, 'min': 1.629711, 'max': 1.663056, 'zp': 1.109648e-10, 'zp_photon': 9.195720e+01, 'toVega':0,      'ext': 0.19,   'system': 'HST' },
              "F170M":   { 'eff': 1.699943, 'min': 1.579941, 'max': 1.837134, 'zp': 1.015711e-10, 'zp_photon': 8.692163e+01, 'toVega':0,      'ext': 0.18,   'system': 'HST' },
              "F190N":   { 'eff': 1.898486, 'min': 1.880845, 'max': 1.917673, 'zp': 6.957714e-11, 'zp_photon': 6.649628e+01, 'toVega':0,      'ext': 0.15,   'system': 'HST' },
              "F215N":   { 'eff': 2.148530, 'min': 2.128579, 'max': 2.168078, 'zp': 4.355167e-11, 'zp_photon': 4.710529e+01, 'toVega':0,      'ext': 0.13,   'system': 'HST' },
              "F336W":   { 'eff': 0.332930, 'min': 0.295648, 'max': 0.379031, 'zp': 3.251259e-09, 'zp_photon': 5.486427e+02, 'toVega':0,      'ext': 1.70,   'system': 'HST' },
              "F390N":   { 'eff': 0.388799, 'min': 0.384000, 'max': 0.393600, 'zp': 5.673647e-09, 'zp_photon': 1.143901e+03, 'toVega':0,      'ext': 1.48,   'system': 'HST' },
              "F475W":   { 'eff': 0.470819, 'min': 0.386334, 'max': 0.556272, 'zp': 5.331041e-09, 'zp_photon': 1.260475e+03, 'toVega':0,      'ext': 1.21,   'system': 'HST' },
              "F555W":   { 'eff': 0.533091, 'min': 0.458402, 'max': 0.620850, 'zp': 4.062007e-09, 'zp_photon': 1.061011e+03, 'toVega':0,      'ext': 1.05,   'system': 'HST' },
              "F625W":   { 'eff': 0.626619, 'min': 0.544589, 'max': 0.709961, 'zp': 2.478260e-09, 'zp_photon': 7.679998e+02, 'toVega':0,      'ext': 0.68,   'system': 'HST' },
              "F656N":   { 'eff': 0.656368, 'min': 0.653838, 'max': 0.658740, 'zp': 1.434529e-09, 'zp_photon': 4.737886e+02, 'toVega':0,      'ext': 0.81,   'system': 'HST' },
              "F673N":   { 'eff': 0.673224, 'min': 0.667780, 'max': 0.678367, 'zp': 1.908442e-09, 'zp_photon': 6.499706e+02, 'toVega':0,      'ext': 0.78,   'system': 'HST' },
              "F775W":   { 'eff': 0.765263, 'min': 0.680365, 'max': 0.863185, 'zp': 1.323662e-09, 'zp_photon': 5.055354e+02, 'toVega':0,      'ext': 0.65,   'system': 'HST' },
              "F850LP":  { 'eff': 0.963736, 'min': 0.832000, 'max': 1.100000, 'zp': 8.069014e-10, 'zp_photon': 3.706372e+02, 'toVega':0,      'ext': 0.46,   'system': 'HST' }}    
  
  if isinstance(band,list):
    for i in Filters.keys(): 
      if Filters[i]['system'] not in band: Filters.pop(i)
    return Filters
  elif isinstance(band,str):
    return Filters[band]

def read_spec(specFiles, errors=True, atomicron=False, negtonan=False, plot=False, linear=False, wlog=False, verbose=True):
    '''
    (by Alejandro N |uacute| |ntilde| ez, Jocelyn Ferrara)
    
    Read spectral data from fits or ascii files. It returns a list of numpy arrays with wavelength in position 0, flux in position 1 and error values (if requested) in position 2. More than one file name can be provided simultaneously.
    
    **Limitations**: Due to a lack of set framework for ascii file headers, this function assumes ascii files to have wavelength in column 1, flux in column 2, and (optional) error in column 3. Ascii spectra are assumed to be linear, so the kwarg *linear* is disabled for ascii files. Fits files that have multiple spectral orders will not be interpreted correctly with this function.
    
    *specFiles*
      String with fits file name (with full path); it can also be a python list of file names.
    *errors*
      Boolean, whether to return error values for the flux data; return nans if unavailable.
    *atomicron*
      Boolean, if wavelength units are in Angstrom, whether to convert them to microns.
    *negtonan*
      Boolean, whether to set negative flux values equal to zero.
    *plot*
      Boolean, whether to plot the spectral data, including error bars when available.
    *linear*
      Boolean, whether to return spectrum only if it is linear. If it cannot verify linearity, it will assume linearity.
    *verbose*
      Boolean, whether to print warning messages.
    '''
    
    # 1. Convert specFiles into a list type if it is only one file name
    if isinstance(specFiles, str):
        specFiles = [specFiles,]
    
    try:
        specFiles[0]
    except TypeError:
        print 'File name(s) in invalid format.'
        return
    
    # 2. Initialize array to store spectra
    specData = [None] * len(specFiles)
    
    # 3. Loop through each file name:
    for spFileIdx,spFile in enumerate(specFiles):
        
        # 3.1 Determine the type of file it is
        isFits = False
        ext = spFile[-4:].lower()
        if ext == 'fits' or ext == '.fit':
            isFits = True
            
        # 3.2. Get data from file
        if isFits:
            try:
                fitsData, fitsHeader = pf.getdata(spFile, header=True)
            except IOError:
                print 'Could not open ' + str(spFile) + '.'
                continue
        # Assume ascii file otherwise (isFits = False)
        else:
            try:
                aData = ad.open(spFile)
                specData[spFileIdx] = [aData[0].tonumpy(), aData[1].tonumpy()]
                if len(aData) >= 3 and errors:
                    specData[spFileIdx].append(aData[2].tonumpy())
                # # Check (when header available) whether data is linear.
                # if aData.header:
                #     lindex = str(aData.header).upper().find('LINEAR')
                #     if lindex == -1:
                #         isLinear = False
                #     else:
                #         isLinear = True
                #     if linear and not isLinear:
                #         if verbose:
                #             print 'Data in ' + spFile + ' is not linear.'
                #         return
            except IOError:
                print 'Could not open ' + str(spFile) + '.'
                continue
        
        # 3.3. Check if data in fits file is linear
        if isFits:
            KEY_TYPE = ['CTYPE1']
            setType  = set(KEY_TYPE).intersection(set(fitsHeader.keys()))
            if len(setType) == 0:
                if verbose:
                    print 'Data in ' + spFile + ' assumed to be linear.'
                isLinear = True
            else:
                valType = fitsHeader[setType.pop()]
                if valType.strip().upper() == 'LINEAR':
                    isLinear = True
                else:
                    isLinear = False
            if linear and not isLinear:
                if verbose:
                    print 'Data in ' + spFile + ' is not linear.'
                return
        
        # 3.4. Get wl, flux & error data from fits file
        #      (returns wl in pos. 0, flux in pos. 1, error values in pos. 2)
        if isFits:
            specData[spFileIdx] = __get_spec(fitsData, fitsHeader, spFile, errors, \
                                             verb=verbose)
            if specData[spFileIdx] is None:
                continue
        
            # Generate wl axis when needed
            if specData[spFileIdx][0] is None:
                specData[spFileIdx][0] = __create_waxis(fitsHeader, len(specData[spFileIdx][1]), spFile, wlog=wlog, verb=verbose)
            # If no wl axis generated, then clear out all retrieved data for object
            if specData[spFileIdx][0] is None:
                specData[spFileIdx] = None
                continue
        
        # 3.5. Convert units in wl-axis from Angstrom into microns if desired
        if atomicron:
            if specData[spFileIdx][0][-1] > 8000:
                specData[spFileIdx][0] = specData[spFileIdx][0] / 10000
        
        # 3.6. Set negative flux values equal to zero (next step sets them to nans)
        if negtonan:
            negIdx = np.where(specData[spFileIdx][1] < 0)
            if len(negIdx[0]) > 0:
                specData[spFileIdx][1][negIdx] = 0
                if verbose:
                    print '%i negative data points found in %s.' \
                            % (len(negIdx[0]), spFile)
        
        # 3.7. Set zero flux values as nans (do this always)
        zeros = np.where(specData[spFileIdx][1] == 0)
        if len(zeros[0]) > 0:
            specData[spFileIdx][1][zeros] = np.nan
        
    
    # 4. Plot the spectra if desired
    if plot:
        plot_spec(specData, ploterrors=True)
    
    # 5. Clear up memory
    fitsHeader = ''
    fitsData   = ''
    
    return specData

def specType(SpT):
  '''
  (By Joe Filippazzo)
  
  Converts between float and letter M, L, T and Y spectral types (e.g. 14.5 => 'L4.5' and 'T3' => 23).
  
  *SpT*
    Float spectral type between 0.0 and 39.9 or letter/number spectral type between M0.0 and Y9.9
  '''
  if isinstance(SpT,str) and SpT[0] in ['M','L','T','Y'] and float(SpT[1:]) < 10:
    try:
      return [l+float(SpT[1:]) for m,l in zip(['M','L','T','Y'],[0,10,20,30]) if m == SpT[0]][0]
    except ValueError:
      print "Spectral type must be a float between 0 and 40 or a string of class M, L, T or Y."
  elif isinstance(SpT,float) or isinstance(SpT,int) and 0.0 <= SpT < 40.0:
    return '{}{}'.format('MLTY'[int(SpT//10)], (int(SpT) if int(SpT)==SpT and isinstance(SpT, float) else SpT) % 10)
  else:
    return SpT
 
def __create_waxis(fitsHeader, lenData, fileName, verb=True, wlog=False):
    # Define key names in
    KEY_MIN  = ['COEFF0','CRVAL1']         # Min wl
    KEY_DELT = ['COEFF1','CDELT1','CD1_1'] # Delta of wl
    KEY_OFF  = ['LTV1']                    # Offset in wl to subsection start
    
    # Find key names for minimum wl, delta, and wl offset in fits header
    setMin  = set(KEY_MIN).intersection(set(fitsHeader.keys()))
    setDelt = set(KEY_DELT).intersection(set(fitsHeader.keys()))
    setOff  = set(KEY_OFF).intersection(set(fitsHeader.keys()))
    
    # Get the values for minimum wl, delta, and wl offset, and generate axis
    if len(setMin) >= 1 and len (setDelt) >= 1:
        nameMin = setMin.pop()
        valMin  = fitsHeader[nameMin]
        
        nameDelt = setDelt.pop()
        valDelt  = fitsHeader[nameDelt]
        
        if len(setOff) == 0:
            valOff = 0
        else:
            nameOff = setOff.pop()
            valOff  = fitsHeader[nameOff]
        
        # generate wl axis
        if nameMin == 'COEFF0' or wlog==True:
            # SDSS fits files
            wAxis = 10 ** (np.arange(lenData) * valDelt + valMin)
        else:
            wAxis = (np.arange(lenData) * valDelt) + valMin - (valOff * valDelt)
        
    else:
        wAxis = None
        if verb:
            print 'Could not re-create wavelength axis for ' + fileName + '.'
    
    return wAxis


def __get_spec(fitsData, fitsHeader, fileName, errorVals, verb=True):
    if errorVals:
        validData = [None] * 3
    else:
        validData = [None] * 2
    
    # Identify number of data sets in fits file
    dimNum = len(fitsData)
    
    # Identify data sets in fits file
    fluxIdx  = None
    waveIdx  = None
    sigmaIdx = None
    
    if dimNum == 1:
        fluxIdx = 0
    elif dimNum == 2:
        if len(fitsData[0]) == 1:
            sampleData = fitsData[0][0][20]
        else:
            sampleData = fitsData[0][20]
        if sampleData < 0.0001:
            # 0-flux, 1-unknown
            fluxIdx  = 0
        else:
            waveIdx = 0
            fluxIdx = 1
    elif dimNum == 3:
        waveIdx  = 0
        fluxIdx  = 1
        sigmaIdx = 2
    elif dimNum == 4:
    # 0-flux clean, 1-flux raw, 2-background, 3-sigma clean
        fluxIdx  = 0
        sigmaIdx = 3
    elif dimNum == 5:
    # 0-flux, 1-continuum substracted flux, 2-sigma, 3-mask array, 4-unknown
        fluxIdx  = 0
        sigmaIdx = 2
    elif dimNum > 10:
    # Implies that only one data set in fits file: flux
        fluxIdx = -1
        if np.isscalar(fitsData[0]):
            fluxIdx = -1
        elif len(fitsData[0]) == 2:
        # Data comes in a xxxx by 2 matrix (ascii origin)
            tmpWave = []
            tmpFlux = []
            for pair in fitsData:
                tmpWave.append(pair[0])
                tmpFlux.append(pair[1])
            fitsData = [tmpWave,tmpFlux]
            fitsData = np.array(fitsData)
            
            waveIdx = 0
            fluxIdx = 1
        else:
        # Indicates that data is structured in an unrecognized way
            fluxIdx = None
    else:
        fluxIdx = None
        
    # Fetch wave data set from fits file
    if fluxIdx is None:
    # No interpretation known for fits file data sets
        validData = None
        if verb:
            print 'Unable to interpret data in ' + fileName + '.'
        return validData
    else:
        if waveIdx is not None:
            if len(fitsData[waveIdx]) == 1:
            # Data set may be a 1-item list
                validData[0] = fitsData[waveIdx][0]
            else:
                validData[0] = fitsData[waveIdx]
    
    # Fetch flux data set from fits file
    if fluxIdx == -1:
        validData[1] = fitsData
    else:
        if len(fitsData[fluxIdx]) == 1:
            validData[1] = fitsData[fluxIdx][0]
        else:
            validData[1] = fitsData[fluxIdx]
    
    # Fetch sigma data set from fits file, if requested
    if errorVals:
        if sigmaIdx is None:
            validData[2] = np.array([np.nan] * len(validData[1]))
        else:
            if len(fitsData[sigmaIdx]) == 1:
                validData[2] = fitsData[sigmaIdx][0]
            else:
                validData[2] = fitsData[sigmaIdx]
        
        # If all sigma values have the same value, replace them with nans
        if validData[2][10] == validData[2][11] == validData[2][12]:
            validData[2] = np.array([np.nan] * len(validData[1]))
    
    return validData