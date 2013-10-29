#!/usr/bin/python
# BDNYC database
import io, xlrd, cPickle, sqlite3 as sql, pyfits as pf, numpy as np, matplotlib.pyplot as plt
# from SEDplots import SEDalt as SED
from astrotools import utilities as u, astrotools as a
path, dbpath = '/Users/Joe/Documents/Python/', '/Users/Joe/Dropbox/BDNYCdb/BDNYC.db'

# ==============================================================================================================================================
# ================================= Data Tables ================================================================================================
# ==============================================================================================================================================

def load_all():
  # Before loading the database, read it to memory, generate the new one and the cross-match for missing objects?
  # I.e. how do I ensure that records loaded manually don't get erased when the database is regenerated?
  load_telescopes(), load_systems(), load_instruments(), load_spectral_types(), load_publications(), load_sources(), load_photometry(), load_spectra(), load_parallaxes(), load_proper_motions(), load_radial_velocities()
    
def load_parallaxes():
  '''
  Get parallaxes in milliarcseconds from parallax data spreadsheet
  '''
  (con, cur), phot, parallaxes = get_db(dbpath, modify=True), SED.parallax_data(phot=True), []
  cur.execute("DROP TABLE IF EXISTS parallaxes"), cur.execute("CREATE TABLE parallaxes (id INTEGER PRIMARY KEY, source_id INTEGER, parallax REAL, parallax_unc REAL, publication_id INTEGER)")    
  workbook, parallaxes = xlrd.open_workbook('/Users/Joe/Documents/Python/Photometry/Parallax_Data.xlsx'), []
  ras, decs, (pi, pi_unc, pi_ref) = [float(i) for i in workbook.sheet_by_index(0).col_values(7)[3:]], [float(i) for i in workbook.sheet_by_index(0).col_values(8)[3:]], [workbook.sheet_by_index(0).col_values(c)[3:] for c in [9,10,15]]
  PI = {obj:{'ra':ras[idx], 'dec':decs[idx], 'pi':pi[idx] if isinstance(pi[idx],float) else None, 'pi_unc':pi_unc[idx] if isinstance(pi_unc[idx],float) else None, 'pi_ref':None if pi_ref[idx]=='null' else str(pi_ref[idx])} for idx,obj in enumerate([str(i) for i in workbook.sheet_by_index(0).col_values(0)[3:]])}

  for ID,RA,DEC in cur.execute("SELECT id,ra,dec FROM sources").fetchall():
    for obj in PI.keys():
      st = PI[obj]
      if u.sameCoords(RA, DEC, st['ra'], st['dec']):
        if st['pi']: parallaxes.append((None, ID, st['pi'], st['pi_unc'], st['pi_ref']))

  cur.executemany("INSERT INTO parallaxes VALUES (?, ?, ?, ?, ?)", parallaxes), con.commit(), con.close()
  
def load_photometry():
  '''  
  Get photometry by RA and DEC from 2MASS, WISE, SDSS, UKIDSS, IRAC surveys (5" search cone)
  '''
  con, cur = get_db(dbpath, modify=True)
  cur.execute("DROP TABLE IF EXISTS photometry")
  cur.execute("CREATE TABLE photometry (id INTEGER PRIMARY KEY, source_id INTEGER, band TEXT, magnitude REAL, magnitude_unc REAL, system INTEGER, telescope_id INTEGER, instrument_id INTEGER, publication_id INTEGER, comments TEXT)")

  # SDSS Photometry
  SS, photometry = u.txt2dict('/Users/Joe/Documents/Python/Photometry/Surveys/SDSS Photometry.txt'), []
  for unum in SS.keys():
    if SS[unum]['type'] == 'STAR':
      for band in ['modelMag_u','modelMag_g','modelMag_r','modelMag_i','modelMag_z']:
        try: photometry.append((None, unum2id(unum), band[-1], SS[unum][band], SS[unum][band[:8]+'Err'+band[8:]], 1, 1, None, None, None))        
        except KeyError: pass

  # WISE/2MASS Photometry
  WM = u.txt2dict('/Users/Joe/Documents/Python/Photometry/Surveys/WISE2MASS Photometry.txt', ignore=['\\'], obj_col=3, start=3)
  lookup = {'j_m_2mass':'J', 'j_m':'J', 'h_m_2mass':'H', 'h_m':'H', 'k_m_2mass':'Ks', 'k_m':'Ks', 'w1mpro':'W1', 'w2mpro':'W2', 'w3mpro':'W3', 'w4mpro':'W4'}
  for unum in WM.keys():
    for band in lookup.keys():
      idx = 2 if 'mass' in band else 3
      try: photometry.append((None, unum2id(unum), lookup[band], WM[unum][band], WM[unum][band[:idx]+'sig'+band[idx:]], idx, idx, None, None, None))
      except KeyError: pass
      try: cur.execute("UPDATE sources SET designation=? WHERE unum=?", (WM[unum]['designation'],unum))
      except IndexError: pass
      try: cur.execute("UPDATE sources SET shortname=? WHERE unum=?", (shortname(WM[unum]['designation']),unum))
      except IndexError: pass

  # # Just 2MASS Photometry
  # M = u.txt2dict('/Users/Joe/Documents/Python/Photometry/Surveys/2MASS Photometry.txt', ignore=['\\'], obj_col=3, start=3)
  # for unum in M.keys():
  #   # if J,H,Ks not in photometry...
  #   # Fix comprehension to get uncertainty
  #   try:
  #     for band in lookup.keys():
  #       photometry.append((None, unum2id(unum), lookup[band], M[unum][band], M[unum][band[:idx]+'sig'+band[idx:]], 2, None, None, None))
  #   except KeyError: pass

  cur.executemany("INSERT INTO photometry VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", photometry)
  con.commit(), con.close()
  
def load_proper_motions():
  (con, cur), phot, proper_motions = get_db(dbpath, modify=True), SED.parallax_data(phot=True), []
  cur.execute("DROP TABLE IF EXISTS proper_motions"), cur.execute("CREATE TABLE proper_motions (id INTEGER PRIMARY KEY, source_id INTEGER, proper_motion_ra REAL, proper_motion_ra_unc REAL, proper_motion_dec REAL, proper_motion_dec_unc REAL, spectrum_id INTEGER, publication_id INTEGER)")
  objects = cur.execute("SELECT id,ra,dec,designation,unum FROM sources")
  for num,ra,dec,designation,unum in objects:
    for p in phot.keys():
      if u.sameCoords(ra, dec, phot[p]['RA'], phot[p]['Dec'], deg=0.001389): # 5" tolerance
        proper_motions.append((None, num, phot[p]['m'], phot[p]['sm'], phot[p]['m'], phot[p]['sm'], None, phot[p]['Ref']))
  cur.executemany("INSERT INTO proper_motions VALUES (?, ?, ?, ?, ?, ?, ?, ?)", proper_motions), con.commit(), con.close()

# def load_radial_velocities():
#   con, cur = get_db(dbpath, modify=True)
#   cur.execute("DROP TABLE IF EXISTS radial_velocities"), cur.execute("CREATE TABLE radial_velocities (id INTEGER PRIMARY KEY, source_id INTEGER, radial_velocity REAL, radial_velocity_unc REAL, spectrum_id INTEGER, publication_id INTEGER)")
#   radial_velocities = []
#   cur.executemany("INSERT INTO radial_velocities VALUES (?, ?, ?, ?, ?, ?)", radial_velocities), con.commit(), con.close()
      
def load_sources():
  ''' 
  Get sources from the original database unums
  '''
  (con, cur), old_db = get_db(dbpath, modify=True), a.browse_db()
  cur.execute("DROP TABLE IF EXISTS sources"), cur.execute("CREATE TABLE sources (id INTEGER PRIMARY KEY, ra REAL, dec REAL, designation TEXT, publication_id INTEGER, comments TEXT, unum TEXT, shortname TEXT, names TEXT)")
  sources = [(None, old_db[unum]['RA'], old_db[unum]['DEC'], None, None, None, unum, None, old_db[unum]['name'].replace(' ','') if old_db[unum]['name'] else None) for unum in old_db.keys()]
  cur.executemany("INSERT INTO sources VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?)", sources), con.commit(), con.close()

def load_spectra():
  '''
  Get spectra from the 1624 FITS files
  '''
  (con, cur), FITS, spectra = get_db(dbpath, modify=True), fitsFiles(), []
  cur.execute("DROP TABLE IF EXISTS spectra"), cur.execute("CREATE TABLE spectra (id INTEGER PRIMARY KEY, source_id INTEGER, wavelength ARRAY, wavelength_units TEXT, flux ARRAY, flux_units TEXT, unc ARRAY, snr ARRAY, wavelength_order INTEGER, publication_id INTEGER, obs_date TEXT, instrument TEXT, telescope TEXT, airmass REAL, filename TEXT, comment TEXT)")

  for filename in FITS.keys():
    f = FITS[filename]
    try: source_id = unum2id(FITS[filename]['unum'])
    except KeyError: source_id = None
    try: wav, flx, err = u.unc(a.read_spec(f['path'], errors=True, atomicron=True, negtonan=True, verbose=False)[0])
    except (TypeError,IndexError):
      print filename
      wav, flx, err = None, None, None
    try: snr = flx/err if any(flx) and any(err) else None
    except (TypeError,IndexError): snr = None
    spectra.append((None, source_id, wav, f['xunits'], flx, f['yunits'], err, snr, None, None, f['date'], f['instr'], f['scope'], f['airmass'], filename, f['name']))

  cur.executemany("INSERT INTO spectra VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", spectra), con.commit(), con.close()

def load_systems():
  con, cur = get_db(dbpath, modify=True)
  cur.execute("DROP TABLE IF EXISTS systems"), cur.execute("CREATE TABLE systems (id INTEGER PRIMARY KEY, name TEXT)")
  systems = [(None,system) for system in ['SDSS','2MASS','WISE','AB','ST','Vega','MKO','CIT','Johnson','Cousins','Landolt','Stromgren']]
  cur.executemany("INSERT INTO systems VALUES (?, ?)", systems), con.commit(), con.close()

def load_telescopes():
  con, cur = get_db(dbpath, modify=True)
  cur.execute("DROP TABLE IF EXISTS telescopes"), cur.execute("CREATE TABLE telescopes (id INTEGER PRIMARY KEY, name TEXT)")
  telescopes = [(None,telescope) for telescope in ['SDSS','2MASS','WISE','UKIDSS','Pan-STARRS','HST','Spitzer','NASA IRTF','Keck I','Keck II','KP 4m','KP 2.1m','KP Bok','MMT','CTIO 1.5m','CTIO 4m','Gemini-North','Gemini-South','ESO VLT U2','ARC 3.5m','Subaru']]
  cur.executemany("INSERT INTO telescopes VALUES (?, ?)", telescopes), con.commit(), con.close()
  
def load_instruments():
  con, cur = get_db(dbpath, modify=True)
  cur.execute("DROP TABLE IF EXISTS instruments"), cur.execute("CREATE TABLE instruments (id INTEGER PRIMARY KEY, name TEXT)")
  instruments = [(None,instrument) for instrument in ['R-C Spec','IRAC','GMOS-N','GMOS-S','FORS1','LRIS','SPeX, IRTF Spectrograph','LDSS3-Two','FOCAS']]
  cur.executemany("INSERT INTO instruments VALUES (?, ?)", instruments), con.commit(), con.close()

def load_publications():  
  # Google Doc of shortname, ADSURL, and description? Easy to import. Example of description?
  con, cur = get_db(dbpath, modify=True)
  cur.execute("DROP TABLE IF EXISTS publications"), cur.execute("CREATE TABLE publications (id INTEGER PRIMARY KEY, shortname TEXT, ADSURL TEXT, description TEXT)")
  publications = []
  cur.executemany("INSERT INTO publications VALUES (?, ?, ?, ?)", publications), con.commit(), con.close()

def load_spectral_types():
  # 12.5 or 'L2.5' format? What about Greek letter gravity indicators, 'peculiar', 'blue/red', etc?
  con, cur = get_db(dbpath, modify=True)
  cur.execute("DROP TABLE IF EXISTS spectral_types"), cur.execute("CREATE TABLE spectral_types (id INTEGER PRIMARY KEY, source_id INTEGER, spectral_type REAL, publication_id INTEGER, regime TEXT, adopted BOOLEAN)")
  workbook, spectral_types = xlrd.open_workbook('/Users/Joe/Documents/Python/Photometry/Parallax_Data.xlsx'), []
  ras, decs, (OPT, OPTref, IR, IRref) = [float(i) for i in workbook.sheet_by_index(0).col_values(7)[3:]], [float(i) for i in workbook.sheet_by_index(0).col_values(8)[3:]], [workbook.sheet_by_index(0).col_values(c)[3:] for c in [82,83,85,86]]
  SpT = {obj:{'ra':ras[idx], 'dec':decs[idx], 'OPT':OPT[idx] if isinstance(OPT[idx],float) else None, 'OPTref':None if OPTref[idx]=='null' else str(OPTref[idx]), 'IR':IR[idx] if isinstance(IR[idx],float) else None, 'IRref':None if IRref[idx]=='null' else str(IRref[idx])} for idx,obj in enumerate([str(i) for i in workbook.sheet_by_index(0).col_values(0)[3:]])}
  
  # Match objects by RA and DEC and load spectral types with source_ids
  for ID,RA,DEC in cur.execute("SELECT id,ra,dec FROM sources").fetchall():
    for obj in SpT.keys():
      st = SpT[obj]
      if u.sameCoords(RA, DEC, st['ra'], st['dec']):
        if st['OPT']: spectral_types.append((None, ID, st['OPT'], st['OPTref'], 'optical', False))
        if st['IR']: spectral_types.append((None, ID, st['IR'], st['IRref'], 'IR', False))
  
  cur.executemany("INSERT INTO spectral_types VALUES (?, ?, ?, ?, ?, ?)", spectral_types), con.commit(), con.close()

# ==============================================================================================================================================
# ================================= Synthetic Spectra Database (separate from BDNYC.db database) ===============================================
# ==============================================================================================================================================

def load_synthetic_spectra():
  import os, glob
  from syn_phot import syn_phot as s
  (con, cur), files, bt_settl = get_db('/Users/Joe/Documents/Python/Models/model_atmospheres.db', modify=True), glob.glob('/Users/Joe/Documents/Python/Models/BT-Settl_M-0.0_a+0.0/*.spec.7'), []
  for f in files:
    obj = s.read_btsettl(f, Flam=False, radius=1, dist=10)
    bt_settl.append((None, obj['Teff'], obj['logg'], obj['W'].magnitude, obj['F'].magnitude, obj['B'].magnitude))
    print "{} {}".format(obj['Teff'], obj['logg'])
  cur.execute("DROP TABLE IF EXISTS bt_settl"), cur.execute("CREATE TABLE bt_settl (id INTEGER PRIMARY KEY, teff INTEGER, logg REAL, wavelength ARRAY, flux ARRAY, blackbody ARRAY)")    
  cur.executemany("INSERT INTO bt_settl VALUES (?, ?, ?, ?, ?, ?)", bt_settl), con.commit(), con.close()

# ==============================================================================================================================================
# ================================= BDdb.py functions ==========================================================================================
# ==============================================================================================================================================

def adapt_array(arr):
  out = io.BytesIO()
  np.save(out, arr)
  out.seek(0)
  # http://stackoverflow.com/a/3425465/190597 (R. Hill)
  return buffer(out.read())

def convert_array(text):
  out = io.BytesIO(text)
  out.seek(0)
  return np.load(out)

# Converts np.array to TEXT when inserting
sql.register_adapter(np.ndarray, adapt_array)
# Converts TEXT to np.array when selecting
sql.register_converter("ARRAY", convert_array)

def fitsFiles(pickle=False, directories=['Spectra/nir_spectra/','Spectra/optical_spectra/','Spectra/IRTF_Library/']):
  if pickle:  
    import glob, os
    
    # Find every .fits file in a directory tree
    FILES, DICT, old_db = [item for sublist in [glob.glob(path+folder+'*.fits') for folder in directories] for item in sublist], {}, a.browse_db()

    # Open Kelle's unum/filename list
    unumFile, UandF = open(path+'Spectra/unum_lookup.txt'), {}
    for i in unumFile: UandF[i.split()[1]] = i.split()[0]                                       
    unumFile.close() 

    for f in FILES:
      filename, header, unum = os.path.basename(f), pf.getheader(f), None

      # RA and DEC
      try:
        if isinstance(header['RA'],float): ra = header['RA']
        elif isinstance(header['RA'],str): ra = float(header['RA']) if header['RA'].isdigit() else a.HMS2deg(ra=header['RA'].replace(':',' ')) 
      except: ra = ''
      try:
        if isinstance(header['DEC'],float): dec = header['DEC']
        elif isinstance(header['DEC'],str): dec = float(header['DEC']) if header['DEC'].isdigit() else a.HMS2deg(dec=header['DEC'].replace(':',' '))
      except: dec = ''
      
      # Check for unum in Kelle's list, then in filename, then by RA and Dec
      if filename in UandF: unum = UandF[filename]
      elif 'U' in filename.upper():
        idx = filename.upper().index('U')
        num = filename[idx+1:idx+6]
        if num.isdigit() and len(num)==5: unum = filename[idx:idx+6].upper()
      elif ra and dec:
        for U in old_db.keys():
          if u.sameCoords(old_db[U]['RA'], old_db[U]['DEC'], ra, dec): unum = U

      # x- and y-units
      try:
        xunits = header['XUNITS'] 
        if 'microns' in xunits or 'Microns' in xunits or 'um' in xunits: xunits = 'um'
      except KeyError:
        try:
           if header['BUNIT']: xunits = 'um'
        except KeyError: xunits = ''
      try: yunits = header['YUNITS'].replace(' ','')
      except KeyError:
        try: yunits = header['BUNIT'].replace(' ','')
        except KeyError: yunits = ''
      if 'erg' in yunits and 'A' in yunits: yunits = 'ergs-1cm-2A-1'
      elif 'erg' in yunits and 'um' in yunits: yunits = 'ergs-1cm-2um-1'
      elif 'W' in yunits and 'um' in yunits: yunits = 'Wm-2um-1'
      elif 'W' in yunits and 'A' in yunits: yunits = 'Wm-2A-1'

      # Date, object name, telescope and instrument
      try: date = header['DATE_OBS']
      except KeyError:
        try: date = header['DATE-OBS']
        except KeyError: date = ''
      try: obj = header['OBJECT']
      except KeyError: obj = ''
      try: scope = header['TELESCOP']
      except KeyError: scope = ''
      try: instr = header['INSTRUME']
      except KeyError: instr = ''
      try: airmass = header['AIRMASS']
      except Exception: airmass = ''
      try: name = old_db[unum]['name']
      except KeyError: name = None

      # Return dictionary of all objects by filename
      DICT[filename] = { 'xunits': xunits, 'yunits': yunits, 'header': header, 'RA': ra, 'DEC': dec, 'date': date, \
                             'airmass':airmass, 'instr': instr, 'scope': scope, 'path': f, 'unum':unum, 'name':name} 

    cPickle.dump(DICT, open(path+'Pickles/fitsFiles.p',"wb"))
    print "{} files added to {}!".format(len(DICT), 'fitsFiles.p')
  else:
    return cPickle.load(open(path+'Pickles/fitsFiles.p','rb'))

def get_db(path, rows=False, modify=False):
  con = sql.connect(path, isolation_level=None, detect_types=sql.PARSE_DECLTYPES)
  con.text_factory = str
  if rows: con.row_factory = sql.Row
  return [con, con.cursor()] if modify else con.cursor()

def spec_arrays(data, rows=False, SNR=False):
  if rows: return [np.array(map(float,i.split())) for i in [data['wavelength'], data['flux'], data['unc']]+([data['SNR']] if SNR else [])]
  else: return [np.array(map(float,i.split())) for i in list(data)]

def inventory(ID=''):
  db = get_db(dbpath)
  if ID: D = db.execute("SELECT id, unum, shortname, ra, dec,  \
                        (SELECT COUNT(*) FROM spectra WHERE spectra.source_id=sources.id), \
                        (SELECT COUNT(*) FROM photometry WHERE photometry.source_id=sources.id), \
                        (SELECT parallax from parallaxes WHERE parallaxes.source_id=sources.id) \
                        FROM sources WHERE id=?", [ID]).fetchall()
  else: D = db.execute("SELECT id, unum, shortname, ra, dec, \
                       (SELECT COUNT(*) FROM spectra WHERE spectra.source_id=sources.id), \
                       (SELECT COUNT(*) FROM photometry WHERE photometry.source_id=sources.id), \
                       (SELECT parallax from parallaxes WHERE parallaxes.source_id=sources.id) FROM sources").fetchall()
  u.printer(['id','unum','shortname','ra','dec','Spec Count','Phot Count','Parallax'], D)

def unum2id(unum): return {str(k):j for j,k in [tuple(i) for i in get_db(dbpath).execute("SELECT id,unum FROM sources")]}[unum] 

def shortname(desig): return desig[desig.index('J')+1:desig.index('J')+5] + desig[desig.index('-' if '-' in desig else '+'):desig.index('-' if '-' in desig else '+')+5]