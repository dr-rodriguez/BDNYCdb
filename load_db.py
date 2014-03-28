#!/usr/bin/python
# BDNYC database
import io, os, glob, xlrd, cPickle, BDdb, sqlite3 as sql, astropy.io.fits as pf, astropy.units as q, astropy.constants as ac, numpy as np, matplotlib.pyplot as plt, utilities as u, astrotools as a
path = '/Users/Joe/Documents/Python/'

# ==============================================================================================================================================
# ==============================================================================================================================================
# ==============================================================================================================================================
# ==============================================================================================================================================
# ==============================================================================================================================================
# ==============================================================================================================================================
# ==============================================================================================================================================
# ==============================================================================================================================================
# ==============================================================================================================================================
# ==============================================================================================================================================
# ==============================================================================================================================================
# ==============================================================================================================================================
# ==============================================================================================================================================
# ==============================================================================================================================================
# ============================================= ADMIN ONLY! PLEASE DO NOT EDIT THIS FILE! ======================================================
# ==============================================================================================================================================
# =============================== load_db is only to be used to generate the database from scratch. ============================================
# ==============================================================================================================================================
# =================================== Running these scripts will delete the current database! ==================================================
# ==============================================================================================================================================
# ==============================================================================================================================================
# ==============================================================================================================================================
# ==============================================================================================================================================
# ==============================================================================================================================================
# ==============================================================================================================================================
# ==============================================================================================================================================
# ==============================================================================================================================================
# ==============================================================================================================================================
# ==============================================================================================================================================
# ==============================================================================================================================================
# ==============================================================================================================================================
# ==============================================================================================================================================
# ==============================================================================================================================================


# ==============================================================================================================================================
# ================================= Fetch Data =================================================================================================
# ==============================================================================================================================================

def parallax_data():
  MK = ['name','Discovery_ref','SpT','flag','flag_ref','HST/AO_inst','HST/AO_ref','ra','dec','pi','pi_sig','mu','mu_sig','PA','PA_sig','pi_ref','alpha_J2000','delta_J2000','Epoch_UT','Epoch_JD','2mass_J','2mass_J_sig','2mass_J_ref','2mass_H','2mass_H_sig','2mass_H_ref','2mass_Ks','2mass_Ks_sig','2mass_Ks_ref','MKO_Y','MKO_Y_sig','MKO_Y_ref','MKO_J','MKO_J_sig','MKO_J_ref','MKO_H','MKO_H_sig','MKO_H_ref','MKO_K','MKO_K_sig','MKO_K_ref','MKO_L','MKO_L_sig','MKO_L_ref','MKO_M','MKO_M_sig','MKO_M_ref','SDSS_u','SDSS_u_sig','SDSS_g','SDSS_g_sig','SDSS_r','SDSS_r_sig','SDSS_i','SDSS_i_sig','SDSS_z','SDSS_z_sig','SDSS_ref','WISE_W1','WISE_W1_sig','WISE_W2','WISE_W2_sig','WISE_W3','WISE_W3_sig','WISE_W4','WISE_W4_sig','WISE_nb','WISE_na','WISE_cc','WISE_ext','WISE_var','WISE_qual','IRAC_CH1','IRAC_CH1_sig','IRAC_CH2','IRAC_CH2_sig','IRAC_CH3','IRAC_CH3_sig','IRAC_CH4','IRAC_CH4_sig','IRAC_ref','OPT_SpT_str','OPT_SpT_num','OPT_SpT_ref','IR_SpT_str','IR_SpT_num','IR_SpT_ref','OPT/IR_SpT','J-2massJ','H-2massH','K-2massK','Y-2massJ','CH4s-H','syn_color_ref','comp_sep_arcsec','comp_mu_alpha','comp_mu_alpha_sig','comp_mu_delta','comp_mu_delta_sig','comp_ref','V_tan','Name_LaTeX','ID']
  D, db, lookup = u.xl2dict(path+'Photometry/Parallax_Data.xlsx', sheet=0, start=3, key_row=1, manual_keys=MK), BDdb.get_db('/Users/Joe/Dropbox/BDNYCdb/BDNYC.db'), u.txt2dict(path+'Photometry/parallax_data_lookup.txt', to_list=True)
  id_unum = db.query.execute("SELECT id, unum FROM sources").fetchall()
  for d in D.keys():
    D[d]['source_id'] = D[d]['unum'] = ''
    
    # First check the lookup list
    for U,obj in lookup:
      if d==obj: 
        D[d]['unum'] = U
        try: D[d]['source_id'] = {str(k):j for j,k in [tuple(i) for i in id_unum]}[D[d]['unum']]
        except KeyError: pass
    
    # If no match, try to match by RA and Dec
    if not D[d]['source_id'] and not D[d]['unum']:
      ra, dec, name, matches = D[d]['alpha_J2000'], D[d]['delta_J2000'], D[d]['name'], []
      for s in db.query.execute("SELECT id, unum, ra, dec FROM sources WHERE ra AND dec").fetchall():
        if u.separation(ra, dec, s[2], s[3]) < 5: matches.append(s)
      if matches: D[d]['source_id'], D[d]['unum'] = sorted(matches, key=lambda x: u.separation(ra, dec, x[2], x[3]))[0][:2]
  
  print '{} / 259 matches.'.format(len([d for d in D if D[d]['unum']]))
  return D, [d for d in D if not D[d]['unum']]

def coordFix(txt='/Users/Joe/Desktop/matches.txt'):
  import astropy.units as q, astropy.coordinates.angles as A
  P, results = parallax_data(), []
  for p in P.keys():
    if P[p]['unum']:
      try:
        match = db.dict.execute("SELECT * FROM sources WHERE unum=?", [P[p]['unum']]).fetchone()
        results.append([P[p]['unum'], P[p]['name'], match['names'] or match['designation'], P[p]['alpha_J2000'], match['ra'], abs((P[p]['alpha_J2000']-match['ra'])*q.degree).to(q.arcsec), P[p]['delta_J2000'], match['dec'], abs((P[p]['delta_J2000']-match['dec'])*q.degree).to(q.arcsec), (np.sqrt((P[p]['alpha_J2000']-match['ra'])**2 + (P[p]['delta_J2000']-match['dec'])**2)*q.degree).to(q.arcsec), u.separation(P[p]['alpha_J2000'], P[p]['delta_J2000'], match['ra'], match['dec'])]) 
      except IndexError: pass
  u.printer(['unum','xl Name','db Name','xl RA (deg)','db RA (deg)','RA diff (")','xl DEC (deg)','db DEC (deg)','DEC diff (")','Pathagoras (")', 'Great Circle (")'], results, to_txt=txt)
  below, above = [i for i in zip(*results)[10] if i<=5], [i for i in zip(*results)[10] if i>5]
  plt.hist(sorted(below),bins=5), plt.hist(sorted(above),bins=45,color='r'), plt.grid(True), plt.xlabel('Separation'), plt.ylabel('Count')
  print "{} matched, {} unmatched => {} percent under 5 arcsec".format(len(below),len(above),len(below)*1./len(results))

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
  lookup = {'j_m_2mass':'J', 'j_m':'J', '2mass_J':'J', 'h_m_2mass':'H', 'h_m':'H', '2mass_H':'H', 'k_m_2mass':'Ks', 'k_m':'Ks', '2mass_Ks':'Ks', 'w1mpro':'W1', 'WISE_W1':'W1', 'w2mpro':'W2', 'WISE_W2':'W2', 'w3mpro':'W3', 'WISE_W3':'W3', 'w4mpro':'W4', 'WISE_W4':'W4', 'IRAC_CH1':'[3.6]', 'IRAC_CH2':'[4.5]', 'IRAC_CH3':'[5.8]', 'IRAC_CH4':'[8]', 'MKO_Y':'Y', 'MKO_J':'J', 'MKO_H':'H', 'MKO_K':'K', 'MKO_L':'L', 'MKO_M':'M'}

  # SDSS Photometry
  SS, photometry = u.txt2dict(path+'Photometry/Surveys/SDSS Photometry_2.txt', all_str=True, delim=','), []
  for ID in SS.keys():
    for band in ['modelMag_u','modelMag_g','modelMag_r','modelMag_i','modelMag_z']:
      try: photometry.append((None, ID, band[-1], SS[ID][band], SS[ID][band[:8]+'Err'+band[8:]], 1, 1, None, None, None))        
      except KeyError: pass

  # WISE/2MASS Photometry
  WM = u.txt2dict(path+'Photometry/Surveys/WISE2MASS Photometry_2.txt', skip=['\\'], ignore=['|'], obj_col=3, start=3)
  for ID in WM.keys():
    for band in lookup.keys():
      try:
        if WM[ID][band] and WM[ID][band]!='null' and WM[ID][band]!='NaN':
          inst, idx = 2 if 'mass' in band else 3, 3 if 'mass' in band else 2
          try: photometry.append((None, ID, lookup[band], WM[ID][band], WM[ID][band[:idx]+'sig'+band[idx:]], inst, inst, None, None, None))
          except KeyError: pass
          try: db.query.execute("UPDATE sources SET designation=? WHERE id=?", (WM[ID]['designation'],ID))
          except IndexError: pass
          try: db.query.execute("UPDATE sources SET shortname=? WHERE id=?", (BDdb.shortname(WM[ID]['designation']),ID))
          except IndexError: pass
      except KeyError: pass

  P = parallax_data()
  import random
  K = [k for k in P[random.choice(P.keys())].keys() if any([i in k for i in ['MKO','SDSS','WISE','IRAC','2mass']]) and '-' not in k]
  for p in P:
    if P[p]['source_id']:
      for k in K:
        if P[p][k] and P[p][k]!='null' and P[p][k]!='NaN' and 'sig' not in k and 'ref' not in k:
          try: 
            system = db.query.execute("SELECT id FROM systems WHERE name LIKE '%{}%'".format(k.split('_')[0])).fetchone()[0]
            scope = 1 if 'SDSS' in k else 2 if '2mass' in k else 3 if 'WISE' in k else 6 if 'IRAC' in k else 9 if 'MKO' in k else None
            photometry.append((None, P[p]['source_id'], lookup[k], P[p][k], P[p][k+'_sig'], system, scope, None, P[p][k+'_ref'], 'From Parallax_Data.xlsx'))
          except: pass
  
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
  db.query.execute("DROP TABLE IF EXISTS spectral_types"), db.query.execute("CREATE TABLE spectral_types (id INTEGER PRIMARY KEY, source_id INTEGER, spectral_type REAL, gravity TEXT, publication_id INTEGER, regime TEXT, adopted BOOLEAN)")
  for obj in PD.keys():
    if PD[obj]['OPT_SpT_num']!='NaN' and PD[obj]['OPT_SpT_num']!='null': spec_types.append((None, PD[obj]['source_id'], PD[obj]['OPT_SpT_num'], PD[obj]['OPT_SpT_ref'], 'OPT', None))
    if PD[obj]['IR_SpT_num']!='NaN' and PD[obj]['IR_SpT_num']!='null': spec_types.append((None, PD[obj]['source_id'], PD[obj]['IR_SpT_num'], PD[obj]['IR_SpT_ref'], 'IR', None))
  
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

def load_log():
  db.query.execute("DROP TABLE IF EXISTS changelog")
  db.query.execute("CREATE TABLE changelog (id INTEGER PRIMARY KEY, event_date DATE, event_type TEXT, before_event TEXT, after_event TEXT)")
  for table in [str(i[0]) for i in db.query.execute("select name from sqlite_master where type = 'table'").fetchall()]:
    for field in [str(i[1]) for i in db.query.execute("PRAGMA table_info({})".format(table)).fetchall()]:
      insert = "CREATE TRIGGER db_insert_{0}_{1} AFTER INSERT ON {0} BEGIN INSERT INTO changelog (event_date, event_type, before_event, after_event) VALUES (DATETIME('NOW'), 'INSERT', None, NEW.{1}); END".format(table,field)
      update = "CREATE TRIGGER db_update_{0}_{1} AFTER UPDATE ON {0} BEGIN INSERT INTO changelog (event_date, event_type, before_event, after_event) VALUES (DATETIME('NOW'), 'UPDATE', OLD.{1}, NEW.{1}); END".format(table,field)
      delete = "CREATE TRIGGER db_delete_{0}_{1} DELETE ON {0} BEGIN INSERT INTO changelog (event_date, event_type, before_event, after_event) VALUES (DATETIME('NOW'), 'DELETE', OLD.{1}, None); END".format(table,field)
      db.query.execute("DROP TRIGGER IF EXISTS {0}.db_insert_{0}_{1}".format(table,field)), db.query.execute("DROP TRIGGER IF EXISTS {0}.db_update_{0}_{1}".format(table,field)), db.query.execute("DROP TRIGGER IF EXISTS {0}.db_delete_{0}_{1}".format(table,field))
      db.query.execute(insert), db.query.execute(update), db.query.execute(delete)
  db.modify.commit()

# ==============================================================================================================================================
# ================================= Synthetic Spectra Database (separate from BDNYC.db database) ===============================================
# ==============================================================================================================================================

def load_marley_spectra():
  syn, files, models = BDdb.get_db('/Users/Joe/Dropbox/BDNYCdb/model_atmospheres.db'), glob.glob(path+'Models/Marley Models/sp_t*'), []
  syn.query.execute("DROP TABLE IF EXISTS marley_saumon"), syn.query.execute("CREATE TABLE marley_saumon (id INTEGER PRIMARY KEY, teff INTEGER, logg REAL, f_sed INTEGER, k_zz REAL, wavelength ARRAY, flux ARRAY)")    

  def read_model(filepath):
    fn = os.path.splitext(os.path.basename(filepath))[0][4:]
    t, f, g, kzz = int(fn.split('g')[0]), 2 if 'f' in fn else None, round(np.log10(int(fn.split('g')[1].split('nc')[0].split('f')[0])*100),1), float(fn.split('kzz')[1]) if 'kzz' in fn else None
    W, F = zip(*u.txt2dict(filepath, to_list=True, start=2))
    W = np.array(map(float, W), dtype='float32')*q.um
    return [W[::-1], (np.array(map(float, F), dtype='float32')*q.erg/q.s/q.cm**2/q.Hz*(ac.c/W**2)).to(q.erg/q.s/q.cm**2/q.AA)[::-1], t, f, g, kzz]
  
  for f in files:
    W, F, t, f, g, kzz = read_model(f)
    syn.query.execute("INSERT INTO marley_saumon VALUES (?, ?, ?, ?, ?, ?, ?)", (None, t, g, f, kzz, W.value, F.value)), syn.modify.commit()
    print "{} {} {} {}".format(t, g, f, kzz)
  syn.modify.close()

def load_bt_settl_spectra(year=2013):
  syn, files, bt_settl = BDdb.get_db('/Users/Joe/Dropbox/BDNYCdb/model_atmospheres.db'), glob.glob(path+'Models/BT-Settl_M-0.0_a+0.0_{}/*.spec.7'.format(year)), []
  syn.query.execute("DROP TABLE IF EXISTS bt_settl_{}".format(year)), syn.query.execute("CREATE TABLE bt_settl_{} (id INTEGER PRIMARY KEY, teff INTEGER, logg REAL, wavelength ARRAY, flux ARRAY)".format(year))    
  
  def read_btsettl(filepath):
    obj, Widx, (T, G) = {}, [], map(float, os.path.splitext(os.path.basename(filepath))[0].replace('lte','').split('-')[:2])
    T, data = int(T*100), open(filepath, 'r')
    lines = [i.split()[:3] for i in data]
    data.close()
    W, F = [[i[idx] for i in lines[::20]] for idx in [0,1]]
    W = (np.array([float(w) for w in W])*q.AA).to(q.um)
    for n,w in enumerate(W):
      if (w.value <= 0.2) or (w.value >= 40.0): Widx.append(n)                                                     
    W, F = np.delete(W,Widx)*q.um, np.delete(F,Widx)
    return [T, G, W, ((10**np.array([float(f.replace('D','E')) for f in F]))*q.erg/q.s/q.cm**3).to(q.erg/q.s/q.cm**2/q.AA)]
    
  for f in files:
    obj = read_btsettl(f)
    try:
      syn.query.execute("INSERT INTO bt_settl_{} VALUES (?, ?, ?, ?, ?)".format(year), (None, obj[0], obj[1], obj[2].value, obj[3].value)), syn.modify.commit()
      print "{} {}".format(obj[0], obj[1])
    except IOError: print "Failed: {} {}".format(obj[0], obj[1])
  syn.modify.close()    