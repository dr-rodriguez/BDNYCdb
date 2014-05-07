# Slimmed down version of astrotools necessary to interact with the SQL database
import astropy.io.ascii as ad, matplotlib.pyplot as plt, numpy as np, astropy.io.fits as pf

def filter_info(band):
  '''
  (By Joe Filippazzo)
   
  Effective, min, and max wavelengths in [um] and zeropoint in [erg s-1 cm-2 A-1] and [photon s-1 cm-2 A-1] for SDSS, Bessel, 2MASS, IRAC and WISE filters. Values from SVO filter profile service.
  
  *band*
      Name of filter band (e.g. 'J' from 2MASS, 'W1' from WISE, etc.) or list of filter systems (e.g. ['SDSS','2MASS','WISE'])
  '''
  Filters = { "FUV":    { 'eff': 0.154226, 'min': 0.134032, 'max': 0.180643, 'zp': 6.486815e-09, 'zp_photon': 5.0360e+02, 'ABtoVega':0, 'ext': 2.62, 'system': 'GALEX' },
              "NUV":    { 'eff': 0.227437, 'min': 0.169252, 'max': 0.300667, 'zp': 4.511632e-09, 'zp_photon': 5.1658e+02, 'ABtoVega':0, 'ext': 2.94, 'system': 'GALEX' },
              "U":      { 'eff': 0.357065, 'min': 0.303125, 'max': 0.417368, 'zp': 3.698499e-09, 'zp_photon': 6.5743e+02, 'ABtoVega':0, 'ext': 1.56, 'system': 'Bessel' }, 
              "B":      { 'eff': 0.437812, 'min': 0.363333, 'max': 0.549706, 'zp': 6.180252e-09, 'zp_photon': 1.3786e+03, 'ABtoVega':0, 'ext': 1.31, 'system': 'Bessel' }, 
              "V":      { 'eff': 0.544579, 'min': 0.473333, 'max': 0.687500, 'zp': 3.521649e-09, 'zp_photon': 9.8165e+02, 'ABtoVega':0, 'ext': 1.02, 'system': 'Bessel' },
              "R":      { 'eff': 0.641420, 'min': 0.550435, 'max': 0.883333, 'zp': 1.882309e-09, 'zp_photon': 6.3446e+02, 'ABtoVega':0, 'ext': 0.83, 'system': 'Bessel' },
              "I":      { 'eff': 0.797880, 'min': 0.704167, 'max': 0.916667, 'zp': 1.132658e-09, 'zp_photon': 4.5495e+02, 'ABtoVega':0, 'ext': 0.61, 'system': 'Bessel' },
              "u":      { 'eff': 0.3543,   'min': 0.304828, 'max': 0.402823, 'zp': 3.652243e-09, 'zp_photon': 6.5486e+02, 'ABtoVega':0.91, 'ext': 1.58, 'system': 'SDSS' },
              "g":      { 'eff': 0.4770,   'min': 0.378254, 'max': 0.554926, 'zp': 5.400197e-09, 'zp_photon': 1.2828e+03, 'ABtoVega':-0.08, 'ext': 1.23, 'system': 'SDSS' },
              "r":      { 'eff': 0.6231,   'min': 0.541534, 'max': 0.698914, 'zp': 2.503279e-09, 'zp_photon': 7.7944e+02, 'ABtoVega':0.16, 'ext': 0.89, 'system': 'SDSS' }, 
              "i":      { 'eff': 0.7625,   'min': 0.668947, 'max': 0.838945, 'zp': 1.398122e-09, 'zp_photon': 5.2785e+02, 'ABtoVega':0.37, 'ext': 0.68, 'system': 'SDSS' }, 
              "z":      { 'eff': 0.9134,   'min': 0.796044, 'max': 1.083325, 'zp': 8.440047e-10, 'zp_photon': 3.8076e+02, 'ABtoVega':0.54, 'ext': 0.52, 'system': 'SDSS' }, 
              "J":      { 'eff': 1.2350,   'min': 1.080647, 'max': 1.406797, 'zp': 3.110467e-10, 'zp_photon': 1.9433e+02, 'ABtoVega':0, 'ext': 0.0166, 'system': '2MASS' },
              "H":      { 'eff': 1.6620,   'min': 1.478738, 'max': 1.823102, 'zp': 1.135347e-10, 'zp_photon': 9.4381e+01, 'ABtoVega':0, 'ext': 0.0146, 'system': '2MASS' },
              "Ks":     { 'eff': 2.1590,   'min': 1.954369, 'max': 2.355240, 'zp': 4.278708e-11, 'zp_photon': 4.6646e+01, 'ABtoVega':0, 'ext': 0.0710, 'system': '2MASS' },
              "MKO_J":  { 'eff': 1.241608, 'min': 1.148995, 'max': 1.348332, 'zp': 3.039310e-10, 'zp_photon': 1.8997e+02, 'ABtoVega':0, 'ext': 0.30, 'system': 'MKO' },
              "MKO_H":  { 'eff': 1.615118, 'min': 1.450318, 'max': 1.808855, 'zp': 1.200648e-10, 'zp_photon': 9.7621e+01, 'ABtoVega':0, 'ext': 0.20, 'system': 'MKO' },
              "MKO_K":  { 'eff': 2.181858, 'min': 1.986393, 'max': 2.397097, 'zp': 4.086224e-11, 'zp_photon': 4.4883e+01, 'ABtoVega':0, 'ext': 0.12, 'system': 'MKO' },
              "MKO_L'": { 'eff': 3.732206, 'min': 3.326622, 'max': 4.207764, 'zp': 5.432256e-12, 'zp_photon': 1.0164e+01, 'ABtoVega':0, 'ext': 0.06, 'system': 'MKO' },
              "MKO_M'": { 'eff': 4.664426, 'min': 4.496502, 'max': 4.865044, 'zp': 2.251348e-12, 'zp_photon': 5.3049e+00, 'ABtoVega':0, 'ext': 0.05, 'system': 'MKO' },
              "DENIS_I":  { 'eff': 0.78621,   'min': 0.7007, 'max': 0.9140, 'zp': 1.182035e-09, 'zp_photon': 4.6817e+02, 'ABtoVega':0, 'ext': 0.63, 'system': 'DENIS' },
              "DENIS_J":  { 'eff': 1.22106,   'min': 1.0508, 'max': 1.3980, 'zp': 3.190860e-10, 'zp_photon': 1.9619e+02, 'ABtoVega':0, 'ext': 0.31, 'system': 'DENIS' },
              "DENIS_Ks": { 'eff': 2.14650,   'min': 1.9474, 'max': 2.3979, 'zp': 4.341472e-11, 'zp_photon': 4.6916e+01, 'ABtoVega':0, 'ext': 0.13, 'system': 'DENIS' },              
              "W1":     { 'eff': 3.4,      'min': 2.754097, 'max': 3.872388, 'zp': 8.033283e-12, 'zp_photon': 1.3751e+01, 'ABtoVega':0, 'ext': 0.07, 'system': 'WISE' },
              "W2":     { 'eff': 4.6,      'min': 3.963326, 'max': 5.341360, 'zp': 2.385688e-12, 'zp_photon': 5.5870e+00, 'ABtoVega':0, 'ext': 0.05, 'system': 'WISE' },
              "W3":     { 'eff': 12,       'min': 7.443044, 'max': 17.26134, 'zp': 5.533295e-14, 'zp_photon': 3.5676e-01, 'ABtoVega':0, 'ext': 0.06, 'system': 'WISE' },
              "W4":     { 'eff': 22,       'min': 19.52008, 'max': 27.91072, 'zp': 4.891490e-15, 'zp_photon': 5.5099e-02, 'ABtoVega':0, 'ext': 0.02, 'system': 'WISE' },
              "[3.6]":  { 'eff': 3.507511, 'min': 3.129624, 'max': 3.961436, 'zp': 6.660957e-12, 'zp_photon': 1.1928e+01, 'ABtoVega':0, 'ext': 0.07, 'system': 'IRAC' },
              "[4.5]":  { 'eff': 4.436578, 'min': 3.917328, 'max': 5.056057, 'zp': 2.685475e-12, 'zp_photon': 6.0902e+00, 'ABtoVega':0, 'ext': 0.05, 'system': 'IRAC' },
              "[5.8]":  { 'eff': 5.628102, 'min': 4.898277, 'max': 6.508894, 'zp': 1.056775e-12, 'zp_photon': 3.0529e+00, 'ABtoVega':0, 'ext': 0.04, 'system': 'IRAC' },
              "[8]":    { 'eff': 7.589159, 'min': 6.299378, 'max': 9.587595, 'zp': 3.089364e-13, 'zp_photon': 1.2329e+00, 'ABtoVega':0, 'ext': 0.03, 'system': 'IRAC' },              
              "F336W":  { 'eff': 0.332930, 'min': 0.295648, 'max': 0.379031, 'zp': 3.251255e-09, 'zp_photon': 5.4864e+02, 'ABtoVega':0, 'ext': 1.70, 'system': 'HST' },
              "F390N":  { 'eff': 0.388799, 'min': 0.384000, 'max': 0.393600, 'zp': 5.673651e-09, 'zp_photon': 1.1439e+03, 'ABtoVega':0, 'ext': 1.48, 'system': 'HST' },
              "F475W":  { 'eff': 0.470819, 'min': 0.386334, 'max': 0.556272, 'zp': 5.331064e-09, 'zp_photon': 1.2605e+03, 'ABtoVega':0, 'ext': 1.21, 'system': 'HST' },
              "F555W":  { 'eff': 0.533091, 'min': 0.458402, 'max': 0.620850, 'zp': 4.062008e-09, 'zp_photon': 1.0610e+03, 'ABtoVega':0, 'ext': 1.05, 'system': 'HST' },
              "F625W":  { 'eff': 0.626619, 'min': 0.544589, 'max': 0.709961, 'zp': 2.478257e-09, 'zp_photon': 7.6800e+02, 'ABtoVega':0, 'ext': 0.68, 'system': 'HST' },
              "F656N":  { 'eff': 0.656368, 'min': 0.653838, 'max': 0.658740, 'zp': 1.427989e-09, 'zp_photon': 4.7157e+02, 'ABtoVega':0, 'ext': 0.81, 'system': 'HST' },
              "F673N":  { 'eff': 0.673224, 'min': 0.667780, 'max': 0.678367, 'zp': 1.908457e-09, 'zp_photon': 6.4997e+02, 'ABtoVega':0, 'ext': 0.78, 'system': 'HST' },
              "F775W":  { 'eff': 0.765263, 'min': 0.680365, 'max': 0.863185, 'zp': 1.323661e-09, 'zp_photon': 5.0553e+02, 'ABtoVega':0, 'ext': 0.65, 'system': 'HST' },
              "F850LP": { 'eff': 0.963736, 'min': 0.832000, 'max': 1.100000, 'zp': 8.068953e-10, 'zp_photon': 3.7063e+02, 'ABtoVega':0, 'ext': 0.46, 'system': 'HST' }}    
  
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
  
  Converts between float and letter/number M, L, T and Y spectral types (e.g. 14.5 => 'L4.5' and 'T3' => 23).
  
  *SpT*
    Float spectral type between 0.0 and 39.9 or letter/number spectral type between M0.0 and Y9.9
  '''
  if isinstance(SpT,str) and SpT[0] in ['M','L','T','Y'] and float(SpT[1:]) < 10:
    try:
      return [l+float(SpT[1:]) for m,l in zip(['M','L','T','Y'],[0,10,20,30]) if m == SpT[0]][0]
    except ValueError:
      print "Spectral type must be a float between 0 and 40 or a string of class M, L, T or Y."
  elif isinstance(SpT,float) or isinstance(SpT,int) and 0.0 <= SpT < 40.0:
    return '{}{}'.format('MLTY'[int(SpT//10)], SpT % 10.)
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