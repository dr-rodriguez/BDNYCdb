#!/usr/bin/python
# BDNYC database
import io, os, sqlite3 as sql, numpy as np, matplotlib.pyplot as plt, pyfits as pf, utilities as u, astrotools.astrotools as a

class get_db:
  def __init__(self, dbpath):
    con = sql.connect(dbpath, isolation_level=None, detect_types=sql.PARSE_DECLTYPES)
    con.text_factory = str
    self.modify = con
    self.query = con.cursor()
    self.dict = con.cursor()
    self.dict.row_factory = sql.Row
    
  def add_data(self, CSV, table):
    '''
    Adds data in *CSV* file to the specified database *table*. Note column names (row 1 of CSV file) must match table fields to insert though order doesn't matter.
    '''
    data, fields, insert = u.txt2dict(CSV, all_str=True, delim=',', to_list=True), zip(*self.query.execute("PRAGMA table_info({})".format(table)).fetchall())[1], []
    columns, query = data.pop(0), "INSERT INTO {} VALUES({})".format(table, ','.join('?'*len(fields)))
    for row in data:
      values = [None for i in fields]
      for field in fields: values[fields.index(field)] = row[columns.index(field)] if field in columns and field!='id' else None
      insert.append(tuple(values))
    u.printer(fields,insert), self.query.executemany(query, insert), self.modify.commit()
    print "{} records added to the {} table.".format(len(data),table.upper())
    
  def add_ascii(self, asciiPath, snrPath=''):
    filename = os.path.basename(asciiPath)
    (name, wav_order, date), data = filename.replace('_',' ').replace('.',' ').split()[:3], zip(*u.txt2dict(asciiPath, to_list=True))
    wav, flx = [np.array(i, dtype='float32') for i in data]
    try:
      snr = np.array(zip(*u.txt2dict(snrPath, to_list=True, start=1))[1][1:], dtype='float32')
      err = flx/snr
    except: snr = err = ''
    xunits, yunits = 'um', 'normalized'
    regime = pub_id = instr = scope = airmass = header = source_id = ''
    
    q = "SELECT id, unum FROM sources WHERE names LIKE '%{0}%' OR shortname LIKE '{0}'".format(name)
    result = self.query.execute(q).fetchall()
    if result:
      if len(result)>1: print result
      else: source_id, unum = result[0]
    
    try: self.query.execute("INSERT INTO spectra VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (None, source_id, wav, xunits, flx, yunits, err, snr, wav_order, regime, pub_id, date, instr, scope, airmass, filename, name, header)), self.modify.commit()
    except: print "Couldn't add file {} to database.".format(filename)

  def add_fits(self, fitsPath, verbose=False):
    '''
    Checks the header of the *fitsFile* and inserts it into the database, with source_id if possible.
    '''
    filename, header, lookup, unum = os.path.basename(fitsPath), pf.getheader(fitsPath), u.txt2dict('/Users/Joe/Documents/Python/Spectra/unum_lookup.txt'), None
    
    # RA and DEC (This is rejecting too many files! Why same RA and DEC for different files, e.g. see Terminal?)
    try: RA, DEC = header['RA'], header['DEC']
    except KeyError: RA, DEC = '', ''
    try:
      ra = RA if isinstance(RA,float) else float(RA) if RA.replace('+','').replace('-','').replace('.','').isdigit() else u.sxg2deg(ra=RA)
      dec = DEC if isinstance(DEC,float) else float(DEC) if DEC.replace('+','').replace('-','').replace('.','').isdigit() else u.sxg2deg(dec=DEC)
      # print "{} {} => {} {}".format(RA,DEC,ra,dec)
    except ValueError:
      print "{}: {} {}".format(filename,RA,DEC)
      ra, dec = '', ''

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
    try:
      n = header['TELESCOP'].lower()
      scope = 5 if 'hst' in n else 6 if 'spitzer' in n else 7 if 'irtf' in n else 9 if 'keck' in n and 'ii' in n else 8 if 'keck' in n and 'i' in n else 10 if 'kp' in n and '4' in n else 11 if 'kp' in n and '2' in n else 12 if 'bok' in n else 13 if 'mmt' in n else 14 if 'ctio' in n and '1' in n else 15 if 'ctio' in n and '4' in n else 16 if 'gemini' in n and 'north' in n else 17 if 'gemini' in n and 'south' in n else 18 if 'vlt' in n else 19 if '3.5m' in n else 20 if 'subaru' in n else 0
    except KeyError: scope = ''
    try: 
      i = header['INSTRUME'].lower()
      instr = 1 if 'r-c spec' in i or 'test' in i or 'nod' in i else 2 if 'gmos-n' in i else 3 if 'gmos-s' in i else 4 if 'fors' in i else 5 if 'lris' in i else 6 if 'spex' in i else 7 if 'ldss3' in i else 8 if 'focas' in i else 0
    except KeyError: instr = ''
    try: airmass = header['AIRMASS']
    except: airmass = 0
    try: name = old_db[unum]['name']
    except: name = ''
    pub_id, wav_order = '', ''
    
    unum = has_unum(filename) or has_unum(name) or has_unum(obj)
    for U in lookup:
      if filename == lookup[U]['filename']: unum = U
    
    if unum:
      try: source_id = {str(k):j for j,k in [tuple(i) for i in self.query.execute("SELECT id,unum FROM sources")]}[unum]
      except KeyError: source_id = None
    else: source_id = None
    
    q = "SELECT id, unum FROM sources WHERE names LIKE '%{0}%' OR shortname LIKE '{0}'".format(filename.replace('.fits',''))
    result = self.query.execute(q).fetchall()
    if result:
      if len(result)>1: print result
      else: source_id, unum = result[0]
    
    # try:
    wav, flx, err = u.unc(a.read_spec(fitsPath, errors=True, atomicron=True, negtonan=True, verbose=False)[0])
    regime = 'OPT' if wav[0]<0.8 and wav[-1]<1.2 else 'NIR' if wav[0]<1.2 and wav[-1]>2 else 'MIR' if wav[-1]>3 else None     
    try: snr = flx/err if any(flx) and any(err) else None
    except (TypeError,IndexError): snr = None
    self.query.execute("INSERT INTO spectra VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (None, source_id, wav, xunits, flx, yunits, err, snr, wav_order, regime, pub_id, date, instr, scope, airmass, filename, name, header)), self.modify.commit()
    if verbose: u.printer(['filename','source_id', 'xunits', 'yunits', 'regime', 'date', 'instr', 'scope', 'airmass', 'name'],[[filename, source_id, xunits, yunits, regime, date, instr, scope, airmass, name]])
    # except: 
      # print [filename, source_id, xunits, yunits, date, instr, scope, airmass, name]
  
  def clean_up(self, table):
    '''
    Removes exact duplicates from the specified *table* keeping the record with the lowest id.
    '''
    query = "DELETE FROM {0} WHERE id NOT IN (SELECT min(id) FROM {0} GROUP BY {1})".format(table,', '.join(zip(*self.query.execute("PRAGMA table_info({})".format(table)).fetchall())[1][1:]))  
    self.query.execute(query), self.modify.commit()
  
  def inventory(self, ID='', with_pi=False, plot=True):
    '''
    Prints a summary of all objects in the database at *dbpath*. If *ID* prints only that object's summary and plots if *plot*
    '''
    if ID: D = self.query.execute("SELECT id, unum, ra, dec, (SELECT filename FROM spectra WHERE spectra.source_id=sources.id AND spectra.regime='OPT'), (SELECT filename FROM spectra WHERE spectra.source_id=sources.id AND spectra.regime='NIR' OR spectra.regime='MIR'), (SELECT COUNT(*) FROM photometry WHERE photometry.source_id=sources.id), (SELECT parallax from parallaxes WHERE parallaxes.source_id=sources.id), (SELECT spectral_type FROM spectral_types WHERE spectral_types.source_id=sources.id AND regime='optical'), (SELECT spectral_type FROM spectral_types WHERE spectral_types.source_id=sources.id AND regime='IR') FROM sources WHERE id=?", [ID]).fetchall()
    else:
      if with_pi: D = self.query.execute("SELECT s.id, s.unum, s.ra, s.dec, (SELECT filename FROM spectra WHERE spectra.source_id=s.id AND spectra.regime='OPT'), (SELECT filename FROM spectra WHERE spectra.source_id=s.id AND spectra.regime='NIR' OR spectra.regime='MIR'), (SELECT COUNT(*) FROM photometry WHERE photometry.source_id=s.id), (SELECT parallax from parallaxes WHERE parallaxes.source_id=s.id), (SELECT spectral_type FROM spectral_types WHERE spectral_types.source_id=s.id AND regime='optical'), (SELECT spectral_type FROM spectral_types WHERE spectral_types.source_id=s.id AND regime='IR') FROM sources AS s JOIN parallaxes AS p ON p.source_id=s.id WHERE p.parallax!='None'").fetchall()
      else: D = self.query.execute("SELECT id, unum, ra, dec, (SELECT filename FROM spectra WHERE spectra.source_id=sources.id AND spectra.regime='OPT'), (SELECT filename FROM spectra WHERE spectra.source_id=sources.id AND spectra.regime='NIR' OR spectra.regime='MIR'), (SELECT COUNT(*) FROM photometry WHERE photometry.source_id=sources.id), (SELECT parallax from parallaxes WHERE parallaxes.source_id=sources.id), (SELECT spectral_type FROM spectral_types WHERE spectral_types.source_id=sources.id AND regime='optical'), (SELECT spectral_type FROM spectral_types WHERE spectral_types.source_id=sources.id AND regime='IR') FROM sources").fetchall()
    if D:
      u.printer(['id','unum','ra','dec','Optical','NIR','Phot Count','Parallax','OPT','IR'], D)
      if ID and plot:
        for i in self.query.execute("SELECT * FROM spectra WHERE source_id=?", [ID]).fetchall(): plt.figure(), plt.loglog(i[2], i[4], c='b'), plt.fill_between(i[2], i[4]-i[6], i[4]+i[6], color='k', alpha=0.2), plt.xlim(0.4,3.0), plt.grid(True), plt.yscale('log',nonposy='clip'), plt.figtext(0.15,0.88,'{}\n{}\n{}\n{}'.format(i[14],i[12],i[11],i[10]),verticalalignment='top')
    else: print "No sources found."

# ==============================================================================================================================================
# ================================= Adapters and converters for special data types =============================================================
# ==============================================================================================================================================
  
def adapt_array(arr):
  out = io.BytesIO()
  np.save(out, arr), out.seek(0)
  return buffer(out.read())

def adapt_header(header): return repr([repr(r) for r in header.ascardlist()])

def convert_array(text):
  out = io.BytesIO(text)
  out.seek(0)
  return np.load(out)

def convert_header(text):
  hdu = pf.PrimaryHDU()
  for i in [eval(l) for l in eval(text)]:
    if i[0]: hdu.header.append(i)
  return hdu.header

sql.register_adapter(np.ndarray, adapt_array), sql.register_adapter(pf.header.Header, adapt_header)
sql.register_converter("ARRAY", convert_array), sql.register_converter("HEADER", convert_header)

# ==============================================================================================================================================
# ================================= Little helper functions ====================================================================================
# ==============================================================================================================================================

def has_unum(text):
  if 'U' in text.upper():
    idx = text.upper().index('U')
    num = text[idx+1:idx+6]
    return text[idx:idx+6].upper() if num.isdigit() and len(num)==5 else None
  else: return None    

def shortname(desig): return desig[desig.index('J')+1:desig.index('J')+5] + desig[desig.index('-' if '-' in desig else '+'):desig.index('-' if '-' in desig else '+')+5]
