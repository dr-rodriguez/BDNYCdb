#!/usr/bin/python
# BDNYC database
import io, os, glob, xlrd, cPickle, BDdb, sqlite3 as sql, pyfits as pf, numpy as np, matplotlib.pyplot as plt, utilities as u
from astrotools import astrotools as a
path = '/Users/Joe/Documents/Python/'
db = BDdb.get_db('/Users/Joe/Dropbox/BDNYCdb/BDNYC.db')

# ==============================================================================================================================================
# ================================= Fetch Data =================================================================================================
# ==============================================================================================================================================

def parallax_data():
  MK = ['name','Discovery_ref','SpT','flag','flag_ref','HST/AO_inst','HST/AO_ref','ra','dec','pi','pi_sig','mu','mu_sig','PA','PA_sig','pi_ref','alpha_J2000','delta_J2000','Epoch_UT','Epoch_JD','2mass_J','2mass_J_sig','2mass_J_ref','2mass_H','2mass_H_sig','2mass_H_ref','2mass_Ks','2mass_Ks_sig','2mass_Ks_ref','MKO_Y','MKO_Y_sig','MKO_Y_ref','MKO_J','MKO_J_sig','MKO_J_ref','MKO_H','MKO_H_sig','MKO_H_ref','MKO_K','MKO_K_sig','MKO_K_ref','MKO_L','MKO_L_sig','MKO_L_ref','MKO_M','MKO_M_sig','MKO_M_ref','SDSS_u','SDSS_u_sig','SDSS_g','SDSS_g_sig','SDSS_r','SDSS_r_sig','SDSS_i','SDSS_i_sig','SDSS_z','SDSS_z_sig','SDSS_ref','WISE_W1','WISE_W1_sig','WISE_W2','WISE_W2_sig','WISE_W3','WISE_W3_sig','WISE_W4','WISE_W4_sig','WISE_nb','WISE_na','WISE_cc','WISE_ext','WISE_var','WISE_qual','IRAC_CH1','IRAC_CH1_sig','IRAC_CH2','IRAC_CH2_sig','IRAC_CH3','IRAC_CH3_sig','IRAC_CH4','IRAC_CH4_sig','IRAC_ref','OPT_SpT_str','OPT_SpT_num','OPT_SpT_ref','IR_SpT_str','IR_SpT_num','IR_SpT_ref','OPT/IR_SpT','J-2massJ','H-2massH','K-2massK','Y-2massJ','CH4s-H','syn_color_ref','comp_sep_arcsec','comp_mu_alpha','comp_mu_alpha_sig','comp_mu_delta','comp_mu_delta_sig','comp_ref','V_tan','Name_LaTeX','ID']
  D, db = u.xl2dict(path+'Photometry/Parallax_Data.xlsx', sheet=0, start=3, key_row=1, manual_keys=MK), BDdb.get_db('/Users/Joe/Dropbox/BDNYCdb/BDNYC.db')
  sources = db.query.execute("SELECT id, unum, ra, dec FROM sources WHERE ra!='' AND dec!=''").fetchall()
  for d in D.keys():
    ra, dec, name, matches = D[d]['ra'], D[d]['dec'], D[d]['name'], []
    for s in sources:
      if u.separation(ra, dec, s[2], s[3]) < 5: matches.append(s)
    if matches: D[d]['source_id'], D[d]['unum'] = sorted(matches, key=lambda x: np.sqrt(x[2]**2 + x[3]**2))[0][:2]
    else: D[d]['source_id'], D[d]['unum'] = '', ''
  return D

# ==============================================================================================================================================
# ================================= Populate db Tables =========================================================================================
# ==============================================================================================================================================

def load_all():
  # Before loading the database, read it to memory, generate the new one and the cross-match for missing objects?
  # I.e. how do I ensure that records loaded manually don't get erased when the database is regenerated?
  load_telescopes(), load_systems(), load_instruments(), load_spectral_types(), load_publications(), load_sources(), load_photometry(), load_spectra(), load_parallaxes(), load_proper_motions(), load_radial_velocities()
    
def load_parallaxes():
  '''
  Get parallaxes in milliarcseconds from parallax data spreadsheet
  '''
  PD, parallaxes = parallax_data(), []
  db.query.execute("DROP TABLE IF EXISTS parallaxes"), db.query.execute("CREATE TABLE parallaxes (id INTEGER PRIMARY KEY, source_id INTEGER, parallax REAL, parallax_unc REAL, publication_id INTEGER)")    

  for obj in PD.keys(): parallaxes.append((None, PD[obj]['source_id'], PD[obj]['pi'], PD[obj]['pi_sig'], PD[obj]['pi_ref']))

  db.query.executemany("INSERT INTO parallaxes VALUES (?, ?, ?, ?, ?)", parallaxes), db.modify.commit()
  
def load_photometry():
  '''  
  Get photometry by RA and DEC from 2MASS, WISE, SDSS, UKIDSS, IRAC surveys (5" search cone)
  '''
  db.query.execute("DROP TABLE IF EXISTS photometry")
  db.query.execute("CREATE TABLE photometry (id INTEGER PRIMARY KEY, source_id INTEGER, band TEXT, magnitude REAL, magnitude_unc REAL, system INTEGER, telescope_id INTEGER, instrument_id INTEGER, publication_id INTEGER, comments TEXT)")

  # SDSS Photometry
  SS, photometry = u.txt2dict(path+'Photometry/Surveys/SDSS Photometry_2.txt', all_str=True, delim=','), []
  for ID in SS.keys():
    for band in ['modelMag_u','modelMag_g','modelMag_r','modelMag_i','modelMag_z']:
      try: photometry.append((None, ID, band[-1], SS[ID][band], SS[ID][band[:8]+'Err'+band[8:]], 1, 1, None, None, None))        
      except KeyError: pass

  # WISE/2MASS Photometry
  WM = u.txt2dict(path+'Photometry/Surveys/WISE2MASS Photometry_2.txt', skip=['\\'], ignore=['|'], obj_col=3, start=3)
  lookup = {'j_m_2mass':'J', 'j_m':'J', 'h_m_2mass':'H', 'h_m':'H', 'k_m_2mass':'Ks', 'k_m':'Ks', 'w1mpro':'W1', 'w2mpro':'W2', 'w3mpro':'W3', 'w4mpro':'W4'}
  for ID in WM.keys():
    for band in lookup.keys():
      inst, idx = 2 if 'mass' in band else 3, 3 if 'mass' in band else 2
      try: photometry.append((None, ID, lookup[band], WM[ID][band], WM[ID][band[:idx]+'sig'+band[idx:]], inst, inst, None, None, None))
      except KeyError: pass
      try: db.query.execute("UPDATE sources SET designation=? WHERE id=?", (WM[ID]['designation'],ID))
      except IndexError: pass
      try: db.query.execute("UPDATE sources SET shortname=? WHERE id=?", (BDdb.shortname(WM[ID]['designation']),ID))
      except IndexError: pass

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
  f = open(path+'Modules/Old_BDNYC_Database/BDNYCData.txt','rb')
  bdnyc = cPickle.load(f)
  f.close()
  old_db = bdnyc.browse()

  db.query.execute("DROP TABLE IF EXISTS sources"), db.query.execute("CREATE TABLE sources (id INTEGER PRIMARY KEY, ra REAL, dec REAL, designation TEXT, publication_id INTEGER, comments TEXT, unum TEXT, shortname TEXT, names TEXT)")
  for unum in old_db.keys():
    try: db.query.execute("INSERT INTO sources VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?)", (None, u.sxg2deg(ra=old_db[unum]['RA']), u.sxg2deg(dec=old_db[unum]['DEC']), None, None, None, unum, None, old_db[unum]['name'] or None)), db.modify.commit()
    except ValueError: pass

def check(filename, head=False):
  paths = u.find(filename,path+'Spectra')
  try:
    for f in paths:
      print f
      wav, flx, err = u.unc(a.read_spec(f, errors=True, atomicron=True, negtonan=True, verbose=False)[0])
      plt.plot(wav,flx), plt.fill_between(wav, flx-err, flx+err), plt.title(filename)
      if head: return pf.getheader(f)
  except KeyError:
    print filename    

def load_spectra(directories=['Spectra/nir_spectra/','Spectra/optical_spectra/','Spectra/IRTF_Library/']):
  '''
  Get spectra from the 1624 FITS files and Emily's ascii spectra
  '''
  db.query.execute("DROP TABLE IF EXISTS spectra"), db.query.execute("CREATE TABLE spectra (id INTEGER PRIMARY KEY, source_id INTEGER, wavelength ARRAY, wavelength_units TEXT, flux ARRAY, flux_units TEXT, unc ARRAY, snr ARRAY, wavelength_order INTEGER, regime TEXT, publication_id INTEGER, obs_date TEXT, instrument_id INTEGER, telescope_id INTEGER, airmass REAL, filename TEXT, comment TEXT, header HEADER)")

  # Add FITS files to database
  for f in [item for sublist in [glob.glob(path+folder+'*.fits') for folder in directories] for item in sublist]:
    try: db.add_fits(f)
    except: print "Couldn't add file {}".format(f)
  
  # Add high-res ascii files to database
  for f in [[''.join([root,d]) for d in dirs] for root,dirs,files in os.walk(path+'Spectra/ascii_files/')][0]:
    try:
      ascii = glob.glob(f+'/*ascii_hc').pop()
      try: snr = glob.glob(f+'/snr.dat').pop()
      except: snr = ''
      db.add_ascii(ascii, snrPath=snr)
    except: print "Couldn't add files from {}".format(f) 
  
  # Add med-res ascii files to database
  for f in glob.glob(path+'Spectra/ascii_files/medres/*.dat'):
    try: db.add_ascii(f)
    except: print "Couldn't add file {}".format(f)

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
  PD, spec_types = parallax_data(), []
  db.query.execute("DROP TABLE IF EXISTS spectral_types"), db.query.execute("CREATE TABLE spectral_types (id INTEGER PRIMARY KEY, source_id INTEGER, spectral_type REAL, publication_id INTEGER, regime TEXT, adopted BOOLEAN)")
  for obj in PD.keys():
    if PD[obj]['OPT_SpT_num']!='NaN' or PD[obj]['OPT_SpT_num']!='null': spec_types.append((None, PD[obj]['source_id'], PD[obj]['OPT_SpT_num'], PD[obj]['OPT_SpT_ref'], 'OPT', None))
    if PD[obj]['IR_SpT_num']!='NaN' or PD[obj]['IR_SpT_num']!='null': spec_types.append((None, PD[obj]['source_id'], PD[obj]['IR_SpT_num'], PD[obj]['IR_SpT_ref'], 'IR', None))
  
  db.query.executemany("INSERT INTO spectral_types VALUES (?, ?, ?, ?, ?, ?)", spec_types), db.modify.commit()

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
  import os, glob, quantities as q
  from syn_phot import syn_phot as s
  syn, files, bt_settl = BDdb.get_db(path+'Models/model_atmospheres.db'), glob.glob(path+'Models/BT-Settl_M-0.0_a+0.0/*.spec.7'), []
  for f in files:
    obj = s.read_btsettl(f, Flam=False, radius=1, dist=10)
    bt_settl.append((None, obj['Teff'], obj['logg'], obj['W'].magnitude, obj['F'].magnitude, ((obj['W']*obj['F']).rescale(q.erg/q.s/q.cm**2)).magnitude, obj['B'].magnitude))
    print "{} {}".format(obj['Teff'], obj['logg'])
  syn.query.execute("DROP TABLE IF EXISTS bt_settl"), syn.query.execute("CREATE TABLE bt_settl (id INTEGER PRIMARY KEY, teff INTEGER, logg REAL, W ARRAY, F ARRAY, lamF ARRAY, blackbody ARRAY)")    
  syn.query.executemany("INSERT INTO bt_settl VALUES (?, ?, ?, ?, ?, ?, ?)", bt_settl), syn.modify.commit(), syn.modify.close()
