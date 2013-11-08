#!/usr/bin/python
# BDNYC database
import io, os, glob, xlrd, cPickle, BDdb, sqlite3 as sql, pyfits as pf, numpy as np, matplotlib.pyplot as plt, utilities as u
from astrotools import astrotools as a
path = '/Users/Joe/Documents/Python/'
db = BDdb.get_db('/Users/Joe/Dropbox/BDNYCdb/BDNYC.db')

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
  phot, parallaxes = SED.parallax_data(phot=True), []
  db.query.execute("DROP TABLE IF EXISTS parallaxes"), db.query.execute("CREATE TABLE parallaxes (id INTEGER PRIMARY KEY, source_id INTEGER, parallax REAL, parallax_unc REAL, publication_id INTEGER)")    
  workbook, parallaxes = xlrd.open_workbook('/Users/Joe/Documents/Python/Photometry/Parallax_Data.xlsx'), []
  ras, decs, (pi, pi_unc, pi_ref) = [float(i) for i in workbook.sheet_by_index(0).col_values(7)[3:]], [float(i) for i in workbook.sheet_by_index(0).col_values(8)[3:]], [workbook.sheet_by_index(0).col_values(c)[3:] for c in [9,10,15]]
  PI = {obj:{'ra':ras[idx], 'dec':decs[idx], 'pi':pi[idx] if isinstance(pi[idx],float) else None, 'pi_unc':pi_unc[idx] if isinstance(pi_unc[idx],float) else None, 'pi_ref':None if pi_ref[idx]=='null' else str(pi_ref[idx])} for idx,obj in enumerate([str(i) for i in workbook.sheet_by_index(0).col_values(0)[3:]])}

  for ID,RA,DEC in db.query.execute("SELECT id,ra,dec FROM sources").fetchall():
    for obj in PI.keys():
      st = PI[obj]
      if u.sameCoords(RA, DEC, st['ra'], st['dec']):
        if st['pi']: parallaxes.append((None, ID, st['pi'], st['pi_unc'], st['pi_ref']))

  db.query.executemany("INSERT INTO parallaxes VALUES (?, ?, ?, ?, ?)", parallaxes), db.modify.commit()
  
def load_photometry():
  '''  
  Get photometry by RA and DEC from 2MASS, WISE, SDSS, UKIDSS, IRAC surveys (5" search cone)
  '''
  db.query.execute("DROP TABLE IF EXISTS photometry")
  db.query.execute("CREATE TABLE photometry (id INTEGER PRIMARY KEY, source_id INTEGER, band TEXT, magnitude REAL, magnitude_unc REAL, system INTEGER, telescope_id INTEGER, instrument_id INTEGER, publication_id INTEGER, comments TEXT)")

  # SDSS Photometry
  SS, photometry = u.txt2dict('/Users/Joe/Documents/Python/Photometry/Surveys/SDSS Photometry.txt', all_str=True), []
  for unum in SS.keys():
    if SS[unum]['type'] == 'STAR':
      for band in ['modelMag_u','modelMag_g','modelMag_r','modelMag_i','modelMag_z']:
        try: photometry.append((None, unum2id(unum), band[-1], SS[unum][band], SS[unum][band[:8]+'Err'+band[8:]], 1, 1, None, None, None))        
        except KeyError: pass

  # WISE/2MASS Photometry
  WM = u.txt2dict('/Users/Joe/Documents/Python/Photometry/Surveys/WISE2MASS Photometry.txt', skip=['\\'], ignore=['|'], obj_col=3, start=3)
  lookup = {'j_m_2mass':'J', 'j_m':'J', 'h_m_2mass':'H', 'h_m':'H', 'k_m_2mass':'Ks', 'k_m':'Ks', 'w1mpro':'W1', 'w2mpro':'W2', 'w3mpro':'W3', 'w4mpro':'W4'}
  for unum in WM.keys():
    for band in lookup.keys():
      idx = 2 if 'mass' in band else 3
      try: photometry.append((None, unum2id(unum), lookup[band], WM[unum][band], WM[unum][band[:idx]+'sig'+band[idx:]], idx, idx, None, None, None))
      except KeyError: pass
      try: db.query.execute("UPDATE sources SET designation=? WHERE unum=?", (WM[unum]['designation'],unum))
      except IndexError: pass
      try: db.query.execute("UPDATE sources SET shortname=? WHERE unum=?", (shortname(WM[unum]['designation']),unum))
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

  db.query.executemany("INSERT INTO photometry VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", photometry), db.modify.commit()
  
def load_proper_motions():
  phot, proper_motions = SED.parallax_data(phot=True), []
  db.query.execute("DROP TABLE IF EXISTS proper_motions"), db.query.execute("CREATE TABLE proper_motions (id INTEGER PRIMARY KEY, source_id INTEGER, proper_motion_ra REAL, proper_motion_ra_unc REAL, proper_motion_dec REAL, proper_motion_dec_unc REAL, spectrum_id INTEGER, publication_id INTEGER)")
  objects = db.query.execute("SELECT id,ra,dec,designation,unum FROM sources")
  for num,ra,dec,designation,unum in objects:
    for p in phot.keys():
      if u.sameCoords(ra, dec, phot[p]['RA'], phot[p]['Dec'], deg=0.001389): # 5" tolerance
        proper_motions.append((None, num, phot[p]['m'], phot[p]['sm'], phot[p]['m'], phot[p]['sm'], None, phot[p]['Ref']))
  db.query.executemany("INSERT INTO proper_motions VALUES (?, ?, ?, ?, ?, ?, ?, ?)", proper_motions), db.modify.commit()

# def load_radial_velocities():
#   con, cur = get_db(dbpath, modify=True)
#   db.query.execute("DROP TABLE IF EXISTS radial_velocities"), db.query.execute("CREATE TABLE radial_velocities (id INTEGER PRIMARY KEY, source_id INTEGER, radial_velocity REAL, radial_velocity_unc REAL, spectrum_id INTEGER, publication_id INTEGER)")
#   radial_velocities = []
#   db.query.executemany("INSERT INTO radial_velocities VALUES (?, ?, ?, ?, ?, ?)", radial_velocities), db.modify.commit()
      
def load_sources():
  ''' 
  Get sources from the original database unums
  '''
  from Old_BDNYC_Database import BDNYC

  # Initialize the database
  f = open('/Users/Joe/Documents/Python/Modules/Old_BDNYC_Database/BDNYCData.txt','rb')
  bdnyc = cPickle.load(f)
  f.close()
  old_db = bdnyc.browse()

  db.query.execute("DROP TABLE IF EXISTS sources"), db.query.execute("CREATE TABLE sources (id INTEGER PRIMARY KEY, ra REAL, dec REAL, designation TEXT, publication_id INTEGER, comments TEXT, unum TEXT, shortname TEXT, names TEXT)")
  for unum in old_db.keys():
    try: db.query.execute("INSERT INTO sources VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?)", (None, u.sxg2deg(ra=old_db[unum]['RA']), u.sxg2deg(dec=old_db[unum]['DEC']), None, None, None, unum, None, old_db[unum]['name'] or None)), db.modify.commit()
    except ValueError: pass

def check(filename, head=False):
  paths = u.find(filename,'/Users/Joe/Documents/Python/Spectra')
  try:
    for f in paths:
      print f
      if head: print pf.getheader(f)
      wav, flx, err = u.unc(a.read_spec(f, errors=True, atomicron=True, negtonan=True, verbose=False)[0])
      plt.plot(wav,flx), plt.fill_between(wav, flx-err, flx+err), plt.title(filename)
  except KeyError:
    print filename    

def load_spectra(directories=['Spectra/nir_spectra/','Spectra/optical_spectra/','Spectra/IRTF_Library/']):
  '''
  Get spectra from the 1624 FITS files
  '''
  db.query.execute("DROP TABLE IF EXISTS spectra"), db.query.execute("CREATE TABLE spectra (id INTEGER PRIMARY KEY, source_id INTEGER, wavelength ARRAY, wavelength_units TEXT, flux ARRAY, flux_units TEXT, unc ARRAY, snr ARRAY, wavelength_order INTEGER, regime TEXT, publication_id INTEGER, obs_date TEXT, instrument_id INTEGER, telescope_id INTEGER, airmass REAL, filename TEXT, comment TEXT, header HEADER)")

  # Add FITS files to database
  for f in [item for sublist in [glob.glob(path+folder+'*.fits') for folder in directories] for item in sublist]:
    try: db.add_fits(f)
    except: print "Couldn't add file {}".format(f)
  
  # for f in [item for sublist in [glob.glob('/Users/Joe/Documents/Python/Spectra/ascii_files/') for folder in directories] for item in sublist]

def load_systems():
  db.query.execute("DROP TABLE IF EXISTS systems"), db.query.execute("CREATE TABLE systems (id INTEGER PRIMARY KEY, name TEXT)")
  systems = [(None,system) for system in ['SDSS','2MASS','WISE','AB','ST','Vega','MKO','CIT','Johnson','Cousins','Landolt','Stromgren']]
  db.query.executemany("INSERT INTO systems VALUES (?, ?)", systems), db.modify.commit()

def load_telescopes():
  db.query.execute("DROP TABLE IF EXISTS telescopes"), db.query.execute("CREATE TABLE telescopes (id INTEGER PRIMARY KEY, name TEXT)")
  telescopes = [(None,telescope) for telescope in ['SDSS','2MASS','WISE','UKIDSS','Pan-STARRS','HST','Spitzer','NASA IRTF','Keck I','Keck II','KP 4m','KP 2.1m','KP Bok','MMT','CTIO 1.5m','CTIO 4m','Gemini-North','Gemini-South','ESO VLT U2','ARC 3.5m','Subaru']]
  db.query.executemany("INSERT INTO telescopes VALUES (?, ?)", telescopes), db.modify.commit()
  
def load_instruments():
  db.query.execute("DROP TABLE IF EXISTS instruments"), db.query.execute("CREATE TABLE instruments (id INTEGER PRIMARY KEY, name TEXT)")
  instruments = [(None,instrument) for instrument in ['R-C Spec','IRAC','GMOS-N','GMOS-S','FORS1','LRIS','SPeX, IRTF Spectrograph','LDSS3-Two','FOCAS']]
  db.query.executemany("INSERT INTO instruments VALUES (?, ?)", instruments), db.modify.commit()

def load_publications():  
  # Google Doc of shortname, ADSURL, and description? Easy to import. Example of description?
  db.query.execute("DROP TABLE IF EXISTS publications"), db.query.execute("CREATE TABLE publications (id INTEGER PRIMARY KEY, shortname TEXT, ADSURL TEXT, description TEXT)")
  publications = []
  db.query.executemany("INSERT INTO publications VALUES (?, ?, ?, ?)", publications), db.modify.commit()

def load_spectral_types():
  # 12.5 or 'L2.5' format? What about Greek letter gravity indicators, 'peculiar', 'blue/red', etc?
  db.query.execute("DROP TABLE IF EXISTS spectral_types"), db.query.execute("CREATE TABLE spectral_types (id INTEGER PRIMARY KEY, source_id INTEGER, spectral_type REAL, publication_id INTEGER, regime TEXT, adopted BOOLEAN)")
  workbook, spectral_types = xlrd.open_workbook('/Users/Joe/Documents/Python/Photometry/Parallax_Data.xlsx'), []
  ras, decs, (OPT, OPTref, IR, IRref) = [float(i) for i in workbook.sheet_by_index(0).col_values(7)[3:]], [float(i) for i in workbook.sheet_by_index(0).col_values(8)[3:]], [workbook.sheet_by_index(0).col_values(c)[3:] for c in [82,83,85,86]]
  SpT = {obj:{'ra':ras[idx], 'dec':decs[idx], 'OPT':OPT[idx] if isinstance(OPT[idx],float) else None, 'OPTref':None if OPTref[idx]=='null' else str(OPTref[idx]), 'IR':IR[idx] if isinstance(IR[idx],float) else None, 'IRref':None if IRref[idx]=='null' else str(IRref[idx])} for idx,obj in enumerate([str(i) for i in workbook.sheet_by_index(0).col_values(0)[3:]])}
  
  # Match objects by RA and DEC and load spectral types with source_ids
  for ID,RA,DEC in db.query.execute("SELECT id,ra,dec FROM sources").fetchall():
    for obj in SpT.keys():
      st = SpT[obj]
      if u.sameCoords(RA, DEC, st['ra'], st['dec']):
        if st['OPT']: spectral_types.append((None, ID, st['OPT'], st['OPTref'], 'optical', False))
        if st['IR']: spectral_types.append((None, ID, st['IR'], st['IRref'], 'IR', False))
  
  db.query.executemany("INSERT INTO spectral_types VALUES (?, ?, ?, ?, ?, ?)", spectral_types), db.modify.commit()

# ==============================================================================================================================================
# ================================= Adapters and converters for special database data types ====================================================
# ==============================================================================================================================================

def adapt_array(arr):
  out = io.BytesIO()
  np.save(out, arr), out.seek(0)
  return buffer(out.read())
  
def convert_array(text):
  out = io.BytesIO(text)
  out.seek(0)
  return np.load(out)

def adapt_header(header):
  return repr([repr(r) for r in header.ascardlist()])

def convert_header(text):
  hdu = pf.PrimaryHDU()
  for i in [eval(l) for l in eval(text)]:
    if i[0]: hdu.header.append(i)
  return hdu.header

sql.register_adapter(np.ndarray, adapt_array)
sql.register_adapter(pf.header.Header, adapt_header)
sql.register_converter("ARRAY", convert_array)
sql.register_converter("HEADER", convert_header)

# ==============================================================================================================================================
# ================================= Synthetic Spectra Database (separate from BDNYC.db database) ===============================================
# ==============================================================================================================================================

def load_synthetic_spectra():
  import os, glob
  from syn_phot import syn_phot as s
  syn, files, bt_settl = BDdb.get_db('/Users/Joe/Documents/Python/Models/model_atmospheres.db'), glob.glob('/Users/Joe/Documents/Python/Models/BT-Settl_M-0.0_a+0.0/*.spec.7'), []
  for f in files:
    obj = s.read_btsettl(f, Flam=False, radius=1, dist=10)
    bt_settl.append((None, obj['Teff'], obj['logg'], obj['W'].magnitude, obj['F'].magnitude, obj['B'].magnitude))
    print "{} {}".format(obj['Teff'], obj['logg'])
  syn.query.execute("DROP TABLE IF EXISTS bt_settl"), syn.query.execute("CREATE TABLE bt_settl (id INTEGER PRIMARY KEY, teff INTEGER, logg REAL, wavelength ARRAY, flux ARRAY, blackbody ARRAY)")    
  syn.query.executemany("INSERT INTO bt_settl VALUES (?, ?, ?, ?, ?, ?)", bt_settl), syn.modify.commit(), syn.modify.close()
