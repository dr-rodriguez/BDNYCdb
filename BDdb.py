#!/usr/bin/python
# BDNYC database
import io, os, sqlite3 as sql, numpy as np, matplotlib.pyplot as plt, astropy.io.fits as pf, utilities as u, astrotools as a
path = '/Users/Joe/Documents/Python/'

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
    Adds data in *CSV* file to the specified database *table*. Note column names (row 1 of CSV file) must match table fields to insert though order and completeness don't matter.
    '''
    data, fields, insert = u.txt2dict(CSV, all_str=True, delim=',', to_list=True), zip(*self.query.execute("PRAGMA table_info({})".format(table)).fetchall())[1], []
    columns, query = data.pop(0), "INSERT INTO {} VALUES({})".format(table, ','.join('?'*len(fields)))
    for row in data:
      values = [None for i in fields]
      for field in fields: values[fields.index(field)] = row[columns.index(field)] if field in columns and field!='id' else None
      insert.append(tuple(values))
    u.printer(fields,insert), self.query.executemany(query, insert), self.modify.commit()
    print "{} records added to the {} table.".format(len(data),table.upper())
   
  def add_ascii(self, asciiPath, header_chars=['#'], skip=[], start=0, source_id='', wavelength_unit='', flux_unit='', publication_id='', obs_date='', wavelength_order='', instrument_id='', telescope_id='', airmass=0, comment=''): 
    filename, data = os.path.basename(asciiPath), zip(*u.txt2dict(asciiPath, to_list=True, skip=header_chars+skip))
    wavelength, flux, unc = [np.array(i, dtype='float32') for i in data][start:start+3]
    snr, regime = flux/unc, 'OPT' if wavelength[0]<0.8 and wavelength[-1]<1.2 else 'NIR' if wavelength[0]<1.2 and wavelength[-1]>2 else 'MIR' if wavelength[-1]>3 else None
    try:
      hdu, h = pf.PrimaryHDU(), [i for i in pf.open(asciiPath) if any([i.startswith(char) for char in header_chars])]
      hdu.header.append(('COMMENT',' '.join([' '.join(i) for n,i in enumerate(h) if any([h[n][0].startswith(char) for char in header_chars])]),''))
      header = hdu.header
    except: header = ''
    try:
      self.query.execute("INSERT INTO spectra VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (None, source_id, wavelength, wavelength_unit, flux, flux_unit, unc, snr, wavelength_order, regime, publication_id, obs_date, instrument_id, telescope_id, airmass, filename, comment, header)), self.modify.commit()
      u.printer(['source_id','wavelength_unit','flux_unit','regime','publication_id','obs_date', 'instrument_id', 'telescope_id', 'airmass', 'filename', 'comment'],[[source_id, wavelength_unit, flux_unit, regime, publication_id, obs_date, instrument_id, telescope_id, airmass, filename, comment]])
    except: print "Couldn't add spectrum to database."
    
  def add_BDSS_ascii(self, asciiPath, snrPath=''):
    '''
    Given the path to an ascii data file *asciiPath* with filename format NAME_ORDER.DATE inserts spectrum into the database. If *snrPath* is provided, generates uncertainty and inserts both arrays.
    '''
    filename, data, lookup = os.path.basename(asciiPath), zip(*u.txt2dict(asciiPath, to_list=True, skip=['#'])), u.txt2dict(path+'Spectra/unum_lookup.txt', to_list=True)
    try: (name, wav_order, date) = filename.replace('.dat','').replace('_',' ').replace('.',' ').split()[:3]
    except ValueError: name, wav_order, date = filename.replace('_',' ').replace('.',' ').split()[0], '', ''
    wav, flx = [np.array(i, dtype='float32') for i in data][:2]
    try:
      snr = np.array(zip(*u.txt2dict(snrPath, to_list=True, start=1))[1][1:], dtype='float32')
      err = flx/snr
    except: snr = err = ''
    xunits, yunits, instr, scope, h = 'um', 'normalized', 9, 9, u.txt2dict(asciiPath, to_list=True)
    pub_id = airmass = ''
    hdu, regime = pf.PrimaryHDU(), 'OPT' if wav[0]<0.8 and wav[-1]<1.2 else 'NIR' if wav[0]<1.2 and wav[-1]>2 else 'MIR' if wav[-1]>3 else None
    hdu.header.append(('COMMENT',' '.join([' '.join(i) for n,i in enumerate(h) if h[n][0].startswith('#')]),''))
    
    unum = has_unum(filename) or has_unum(name)
    for U,fn in lookup:
      if filename == fn: unum = U
    
    if unum:
      try: source_id = {str(k):j for j,k in [tuple(i) for i in self.query.execute("SELECT id, unum FROM sources")]}[unum]
      except KeyError: pass
    else:
      q = "SELECT id, unum FROM sources WHERE names LIKE '%{0}%' OR shortname LIKE '%{0}%' OR designation LIKE '%{0}%'".format(filename.replace('_ascii_hc','').replace('.dat',''))
      result = self.query.execute(q).fetchall()
      if result:
        if len(result)>1: print name, shortname, unum, result
        else: source_id, unum = result[0]
      else: source_id = unum = ''
    
    try: self.query.execute("INSERT INTO spectra VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (None, source_id, wav, xunits, flx, yunits, err, snr, wav_order, regime, pub_id, date, instr, scope, airmass, filename, name, hdu.header)), self.modify.commit()
    except: print "Couldn't add file {} to database.".format(filename)

  def add_fits(self, fitsPath, verbose=False):
    '''
    Checks the header of the *fitsFile* and inserts it into the database, with source_id if possible.
    '''
    filename, header, lookup, unum = os.path.basename(fitsPath), pf.getheader(fitsPath), u.txt2dict(path+'Spectra/unum_lookup.txt', to_list=True), None
    
    # RA and DEC (This is rejecting too many files! Why same RA and DEC for different files, e.g. see Terminal?)
    try: RA, DEC = header['RA'], header['DEC']
    except KeyError: RA, DEC = '', ''
    try:
      ra = RA if isinstance(RA,float) else float(RA) if RA.replace('+','').replace('-','').replace('.','').isdigit() else u.sxg2deg(ra=RA)
      dec = DEC if isinstance(DEC,float) else float(DEC) if DEC.replace('+','').replace('-','').replace('.','').isdigit() else u.sxg2deg(dec=DEC)
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
      instr = 1 if 'r-c spec' in i or 'test' in i or 'nod' in i else 2 if 'gmos-n' in i else 3 if 'gmos-s' in i else 4 if 'fors' in i else 5 if 'lris' in i else 6 if 'spex' in i else 7 if 'ldss3' in i else 8 if 'focas' in i else 9 if 'nirspec' in i else 0
    except KeyError: instr = ''
    try: airmass = header['AIRMASS']
    except: airmass = 0
    try: name = old_db[unum]['name']
    except: name = ''
    pub_id, wav_order = '', ''
    
    unum = has_unum(filename) or has_unum(name) or has_unum(obj)
    for U,fn in lookup:
      if filename == fn: unum = U
    
    if unum:
      try: source_id = {str(k):j for j,k in [tuple(i) for i in self.query.execute("SELECT id,unum FROM sources")]}[unum]
      except KeyError: 
        try:
          self.query.execute("INSERT INTO sources VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?)", (None, ra, dec, None, None, None, unum, None, None)), self.modify.commit()
          q = "SELECT * FROM sources WHERE unum='{}'".format(unum)
          source_id = self.dict.execute(q).fetchone()['id']
        except: source_id = None
    else:
      q = "SELECT id, unum FROM sources WHERE names LIKE '%{0}%' OR shortname LIKE '%{0}%' OR designation LIKE '%{0}%'".format(filename.replace('.fits',''))
      result = self.query.execute(q).fetchall()
      if result:
        if len(result)>1: print name, shortname, unum, result
        else: source_id, unum = result[0]
      else: source_id = unum = ''

    try:
      wav, flx, err = u.unc(a.read_spec(fitsPath, errors=True, atomicron=True, negtonan=True, verbose=False)[0])
      regime = 'OPT' if wav[0]<0.8 and wav[-1]<1.2 else 'NIR' if wav[0]<1.2 and wav[-1]>2 else 'MIR' if wav[-1]>3 else None     
      try: snr = flx/err if any(flx) and any(err) else None
      except (TypeError,IndexError): snr = None
      try: self.query.execute("INSERT INTO spectra VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (None, source_id, wav, xunits, flx, yunits, err, snr, wav_order, regime, pub_id, date, instr, scope, airmass, filename, name, header)), self.modify.commit()
      except: self.query.execute("INSERT INTO spectra VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (None, source_id, wav, xunits, flx, yunits, err, snr, wav_order, regime, pub_id, date, instr, scope, airmass, filename, name, None)), self.modify.commit()
      if verbose: u.printer(['filename','source_id', 'xunits', 'yunits', 'regime', 'date', 'instr', 'scope', 'airmass', 'name'],[[filename, source_id, xunits, yunits, regime, date, instr, scope, airmass, name]])
    except: 
      print [filename, source_id, xunits, yunits, date, instr, scope, airmass, name]
  
  def clean_up(self, table):
    '''
    Removes exact duplicates from the specified *table* keeping the record with the lowest id.
    '''
    query = "DELETE FROM {0} WHERE id NOT IN (SELECT min(id) FROM {0} GROUP BY {1})".format(table,', '.join(zip(*self.query.execute("PRAGMA table_info({})".format(table)).fetchall())[1][1:]))  
    self.query.execute(query), self.modify.commit()
  
  def inventory(self, ID='', with_pi=False, SED=False, plot=True, data=False):
    '''
    Prints a summary of all objects in the database at *dbpath*. If *ID* prints only that object's summary and plots if *plot*
    '''
    q = "SELECT id, unum, designation, ra, dec, (SELECT COUNT(*) FROM spectra WHERE spectra.source_id=sources.id AND spectra.regime='OPT'), (SELECT COUNT(*) FROM spectra WHERE spectra.source_id=sources.id AND spectra.regime='NIR'), (SELECT COUNT(*) FROM photometry WHERE photometry.source_id=sources.id), (SELECT parallax from parallaxes WHERE parallaxes.source_id=sources.id), (SELECT parallax_unc from parallaxes WHERE parallaxes.source_id=sources.id), (SELECT spectral_type FROM spectral_types WHERE spectral_types.source_id=sources.id AND regime='OPT'), (SELECT spectral_type FROM spectral_types WHERE spectral_types.source_id=sources.id AND regime='IR') FROM sources"
    if ID: q += ' WHERE id={}'.format(ID)
    elif with_pi and not ID: q = "SELECT s.id, s.unum, s.designation, s.ra, s.dec, (SELECT COUNT(*) FROM spectra WHERE spectra.source_id=s.id AND spectra.regime='OPT'), (SELECT COUNT(*) FROM spectra WHERE spectra.source_id=s.id AND spectra.regime='NIR'), (SELECT COUNT(*) FROM photometry WHERE photometry.source_id=s.id), (SELECT parallax from parallaxes WHERE parallaxes.source_id=s.id), (SELECT parallax_unc from parallaxes WHERE parallaxes.source_id=sources.id), (SELECT spectral_type FROM spectral_types WHERE spectral_types.source_id=s.id AND regime='OPT'), (SELECT spectral_type FROM spectral_types WHERE spectral_types.source_id=s.id AND regime='IR') FROM sources AS s JOIN parallaxes AS p ON p.source_id=s.id WHERE p.parallax!='None'"
    elif SED and not ID: q = "SELECT s.id, s.unum, s.designation, s.ra, s.dec, (SELECT COUNT(*) FROM spectra WHERE spectra.source_id=s.id AND spectra.regime='OPT'), (SELECT COUNT(*) FROM spectra WHERE spectra.source_id=s.id AND spectra.regime='NIR'), (SELECT COUNT(*) FROM photometry WHERE photometry.source_id=s.id), (SELECT parallax from parallaxes WHERE parallaxes.source_id=s.id), (SELECT parallax_unc from parallaxes WHERE parallaxes.source_id=sources.id), (SELECT spectral_type FROM spectral_types WHERE spectral_types.source_id=s.id AND regime='OPT'), (SELECT spectral_type FROM spectral_types WHERE spectral_types.source_id=s.id AND regime='IR') FROM sources AS s JOIN parallaxes AS p ON p.source_id=s.id WHERE p.parallax!='None' AND (SELECT COUNT(*) FROM spectra WHERE spectra.source_id=s.id AND spectra.regime='NIR')>0 AND (SELECT COUNT(*) FROM photometry WHERE photometry.source_id=s.id) >0"
    D = self.query.execute(q).fetchall()
    if D:
      u.printer(['id','unum','name','ra','dec','Optical','NIR','Phot Count','Pi','Pi_unc','OPT','IR'], D)
      if ID and plot:
        for i in self.dict.execute("SELECT s.wavelength,s.flux,s.unc,s.filename,t.name,i.name,s.obs_date FROM spectra AS s JOIN telescopes AS t JOIN instruments as i ON s.instrument_id=i.id AND s.telescope_id=t.id WHERE s.source_id=?", [ID]).fetchall(): 
          plt.figure(), plt.loglog(i[0], i[1], c='b'), plt.xlim(0.4,3.0), plt.grid(True), plt.yscale('log', nonposy='clip')
          try: 
            Y = plt.ylim()
            plt.fill_between(i[0], i[1]-i[2], i[1]+i[2], alpha=0.2, color='b'), plt.ylim(Y)
          except: pass
          try: plt.figtext(0.15,0.88, '{}\n{}\n{}\n{}'.format(i[3].replace('_',' '),i[4],i[5],i[6]), verticalalignment='top')
          except: pass
      if data: return D
    else: print "No sources found."
    
  def identify(self, search):
    try: q = "SELECT id,ra,dec,designation,unum,shortname,names FROM sources WHERE ra BETWEEN "+str(search[0]-0.01667)+" AND "+str(search[0]+0.01667)+" AND dec BETWEEN "+str(search[1]-0.01667)+" AND "+str(search[1]+0.01667)
    except TypeError: q = "SELECT id,ra,dec,designation,unum,shortname,names FROM sources WHERE names like '%"+search+"%' or designation like '%"+search+"%'"
    u.printer(['id','ra','dec','designation','unum','short','names'], self.query.execute(q).fetchall())

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

def associate(unum, filename, lookupPath=''):  
  if len(unum) == 6 and unum.startswith('U') and unum[1:].isdigit():
    lookup = lookupPath or path+'Spectra/unum_lookup.txt'
    with open(lookup, "a") as unumFile:
      unumFile.write('\n'+unum+'\t'+filename) 
    print "File {} linked to object {}.".format(filename,unum)
  else: print "Could not update list with unum '{}' and filename '{}'.".format(unum,filename)

def has_unum(text):
  if 'U' in text.upper():
    idx = text.upper().index('U')
    num = text[idx+1:idx+6]
    return text[idx:idx+6].upper() if num.isdigit() and len(num)==5 else None
  else: return None    

def shortname(desig): return desig[desig.index('J')+1:desig.index('J')+5] + desig[desig.index('-' if '-' in desig else '+'):desig.index('-' if '-' in desig else '+')+5]
