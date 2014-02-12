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
   
  def add_ascii(self, asciiPath, header_chars=['#'], skip=[], start=0, source_id='', wavelength_units='', flux_units='', publication_id='', obs_date='', wavelength_order='', instrument_id='', telescope_id='', airmass=0, comment=''): 
    filename, data = os.path.basename(asciiPath), zip(*u.txt2dict(asciiPath, to_list=True, skip=header_chars+skip))
    wavelength, flux, unc = [np.array(i, dtype='float32') for i in data][start:start+3]
    snr, regime = flux/unc, 'OPT' if wavelength[0]<0.8 and wavelength[-1]<1.2 else 'NIR' if wavelength[0]<1.2 and wavelength[-1]>2 else 'MIR' if wavelength[-1]>3 else None
    try:
      hdu, h = pf.PrimaryHDU(), [i for i in pf.open(asciiPath) if any([i.startswith(char) for char in header_chars])]
      hdu.header.append(('COMMENT',' '.join([' '.join(i) for n,i in enumerate(h) if any([h[n][0].startswith(char) for char in header_chars])]),''))
      header = hdu.header
    except: header = ''
    try:
      self.query.execute("INSERT INTO spectra VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (None, source_id or unum2id, wavelength, wavelength_units, flux, flux_units, unc, snr, wavelength_order, regime, publication_id, obs_date, instrument_id, telescope_id, airmass, filename, comment, header)), self.modify.commit()
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

  def add_fits(self, fitsPath, source_id='', wavelength_units='', flux_units='', publication_id='', obs_date='', wavelength_order='', instrument_id='', telescope_id='', airmass=0, comment=''):
    '''
    Checks the header of the *fitsFile* and inserts it into the database, with source_id if possible.
    '''
    filename, header = os.path.basename(fitsPath), pf.getheader(fitsPath)

    # x- and y-units
    if not wavelength_units:
      try:
        wavelength_units = header['XUNITS'] 
        if 'microns' in wavelength_units or 'Microns' in wavelength_units or 'um' in wavelength_units: wavelength_units = 'um'
      except KeyError:
        try:
           if header['BUNIT']: wavelength_units = 'um'
        except KeyError: wavelength_units = ''
    if not flux_units:
      try: flux_units = header['YUNITS'].replace(' ','')
      except KeyError:
        try: flux_units = header['BUNIT'].replace(' ','')
        except KeyError: flux_units = ''
    if 'erg' in flux_units and 'A' in flux_units: flux_units = 'ergs-1cm-2A-1' if 'erg' in flux_units and 'A' in flux_units else 'ergs-1cm-2um-1' if 'erg' in flux_units and 'um' in flux_units else 'Wm-2um-1' if 'W' in flux_units and 'um' in flux_units else 'Wm-2A-1' if 'W' in flux_units and 'A' in flux_units else ''

    # Date, object name, telescope and instrument
    if not obs_date:
      try: obs_date = header['DATE_OBS']
      except KeyError:
        try: obs_date = header['DATE-OBS']
        except KeyError: date = ''
    if not telescope_id:
      try:
        n = header['TELESCOP'].lower()
        telescope_id = 5 if 'hst' in n else 6 if 'spitzer' in n else 7 if 'irtf' in n else 9 if 'keck' in n and 'ii' in n else 8 if 'keck' in n and 'i' in n else 10 if 'kp' in n and '4' in n else 11 if 'kp' in n and '2' in n else 12 if 'bok' in n else 13 if 'mmt' in n else 14 if 'ctio' in n and '1' in n else 15 if 'ctio' in n and '4' in n else 16 if 'gemini' in n and 'north' in n else 17 if 'gemini' in n and 'south' in n else 18 if 'vlt' in n else 19 if '3.5m' in n else 20 if 'subaru' in n else 0
      except KeyError: telescope_id = ''
    if not instrument_id:
      try: 
        i = header['INSTRUME'].lower()
        instrument_id = 1 if 'r-c spec' in i or 'test' in i or 'nod' in i else 2 if 'gmos-n' in i else 3 if 'gmos-s' in i else 4 if 'fors' in i else 5 if 'lris' in i else 6 if 'spex' in i else 7 if 'ldss3' in i else 8 if 'focas' in i else 9 if 'nirspec' in i else 0
      except KeyError: instrument_id = ''
    try: airmass = header['AIRMASS']
    except: airmass = 0
    
    try:
      wav, flx, err = u.unc(a.read_spec(fitsPath, errors=True, atomicron=True, negtonan=True, verbose=False)[0])
      regime = 'OPT' if wav[0]<0.8 and wav[-1]<1.2 else 'NIR' if wav[0]<1.2 and wav[-1]>2 else 'MIR' if wav[-1]>3 else None     
      try: snr = flx/err if any(flx) and any(err) else None
      except (TypeError,IndexError): snr = None
      try: self.query.execute("INSERT INTO spectra VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (None, source_id, wav, wavelength_units, flx, flux_units, err, snr, wavelength_order, regime, publication_id, obs_date, instrument_id, telescope_id, airmass, filename, comment, header)), self.modify.commit()
      except: self.query.execute("INSERT INTO spectra VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (None, source_id, wav, wavelength_units, flx, flux_units, err, snr, wavelength_order, regime, publication_id, obs_date, instrument_id, telescope_id, airmass, filename, comment, None)), self.modify.commit()
      u.printer(['filename','source_id', 'xunits', 'yunits', 'regime', 'date', 'instr', 'scope', 'airmass', 'name'],[[filename, source_id, wavelength_units, flux_units, regime, obs_date, instrument_id, telescope_id, airmass, comment]])
    except: 
      print [filename, source_id, wavelength_units, flux_units, obs_date, instrument_id, telescope_id, airmass, comment]
  
  def clean_up(self, table):
    '''
    Removes exact duplicates from the specified *table* keeping the record with the lowest id.
    '''
    # First pass to delete only exact duplicates, i.e every column value is identical.
    columns = zip(*self.query.execute("PRAGMA table_info(spectra)").fetchall())[1]
    query = "DELETE FROM {0} WHERE id NOT IN (SELECT min(id) FROM {0} GROUP BY {1})".format(table,', '.join(columns))  
    self.query.execute(query), self.modify.commit()
    
    # Second pass to delete duplicate spectra if the flux arrays are identical, choosing the record with the greatest number of column values
    if table=='spectra':
      # score = "SELECT (CASE WHEN {}".format(' IS NOT null THEN 1 ELSE 0 END) + (CASE WHEN '.join(columns))+" IS NOT null THEN 1 ELSE 0 END),id FROM spectra"      
      self.query.execute("DELETE FROM spectra WHERE id NOT IN (SELECT max(id) FROM spectra GROUP BY flux)"), self.modify.commit()      
      
  def inventory(self, ID='', unum='', with_pi=False, SED=False, plot=False, data=False):
    '''
    Prints a summary of all objects in the database. Input string or list of strings in *ID* or *unum* for specific objects.
    '''
    q = "SELECT sources.id, sources.unum, sources.designation, sources.ra, sources.dec, (SELECT COUNT(*) FROM spectra WHERE spectra.source_id=sources.id), (SELECT COUNT(*) FROM spectra WHERE spectra.source_id=sources.id AND spectra.regime='OPT'), (SELECT COUNT(*) FROM spectra WHERE spectra.source_id=sources.id AND spectra.regime='NIR'), (SELECT COUNT(*) FROM photometry WHERE photometry.source_id=sources.id), (SELECT parallax from parallaxes WHERE parallaxes.source_id=sources.id), (SELECT parallax_unc from parallaxes WHERE parallaxes.source_id=sources.id), (SELECT spectral_type FROM spectral_types WHERE spectral_types.source_id=sources.id AND regime='OPT'), (SELECT spectral_type FROM spectral_types WHERE spectral_types.source_id=sources.id AND regime='IR'), (SELECT gravity from spectral_types WHERE spectral_types.source_id=sources.id) FROM sources"
    if ID or unum:
      if unum: 
        unums = unum if isinstance(unum,list) else [unum]
        IDS = zip(*self.dict.execute("SELECT * FROM sources WHERE unum IN ({})".format("'"+"','".join(unums)+"'")).fetchall())[0]
      else: IDS = ID if isinstance(ID,list) else [ID]  
      q += ' WHERE id IN ({})'.format("'"+"','".join(map(str,IDS))+"'")
    elif with_pi: q += " JOIN parallaxes ON parallaxes.source_id=sources.id WHERE parallaxes.parallax!='None'"
    elif SED: q += " JOIN parallaxes ON parallaxes.source_id=sources.id WHERE parallaxes.parallax!='None' AND (SELECT COUNT(*) FROM spectra WHERE spectra.source_id=sources.id AND spectra.regime='NIR')>0 AND (SELECT COUNT(*) FROM photometry WHERE photometry.source_id=sources.id) >0" 
    try:
      D = self.query.execute(q).fetchall()
      if D:
        u.printer(['id','unum','name','ra','dec','Spec Count','Optical','NIR','Phot Count','Pi','Pi_unc','OPT','IR','grav'], D)
        if plot:
          for I in ID:
            spectra = self.dict.execute("SELECT s.wavelength,s.flux,s.unc,s.filename,t.name,i.name,s.obs_date,s.source_id FROM spectra AS s JOIN telescopes AS t JOIN instruments as i ON s.instrument_id=i.id AND s.telescope_id=t.id WHERE s.source_id={}".format(I)).fetchall() or self.dict.execute("SELECT wavelength,flux,filename FROM spectra WHERE source_id={}".format(I)).fetchall()
            for i in spectra: plt.figure(), plt.rc('text', usetex=False, fontsize=12), plt.loglog(i[0], i[1], c='b'), plt.grid(True), plt.yscale('log', nonposy='clip'), plt.title('source_id = {}'.format(i[7])), plt.figtext(0.15,0.88, '{}\n{}\n{}\n{}'.format(i[3].replace('_',' '),i[4],i[5],i[6]), verticalalignment='top')
        if data: return D
      else: print "No sources found{}.".format(' with id '+str(ID) if ID else '')
    except: pass
    
  def identify(self, search):
    try: q = "SELECT id,ra,dec,designation,unum,shortname,names FROM sources WHERE ra BETWEEN "+str(search[0]-0.01667)+" AND "+str(search[0]+0.01667)+" AND dec BETWEEN "+str(search[1]-0.01667)+" AND "+str(search[1]+0.01667)
    except TypeError: q = "SELECT id,ra,dec,designation,unum,shortname,names FROM sources WHERE names like '%"+search+"%' or designation like '%"+search+"%'"
    u.printer(['id','ra','dec','designation','unum','short','names'], self.query.execute(q).fetchall())

  def merge(self, conflicted, master):
    pass
    # Compare conflicted and master for unique data in conflicted.
    # Print all new/modified records for visual inspection and store in change_log table.
    # Pull out all unique conflicted records and append them to their respective tables with new ids.
    # Save database master and move conflicted copy to Resolved Version Histories directory.
    
  def output_spectrum(self, spectrum_id, fmt='ascii', filename=None):
    '''
    Prints a file of the spectrum with id *spectrum_id* to an *ascii* or *fits* file.
    '''
    data = self.dict.execute("SELECT * FROM spectra WHERE id={}".format(spectrum_id)).fetchone()
    if data:
      import csv
      fn, header = '/Users/Joe/Desktop/{}.txt'.format(filename or data['filename']), repr([repr(r) for r in data['header'].ascardlist()])
      if header:
        with open(fn, 'w' ) as f:
          keys, vals, coms = zip(*[eval(l) for l in eval(header)])
          klen, vlen, clen = [len(max(['1' if isinstance(j,bool) and j else '0' if isinstance(j,bool) else str(j) for j in i], key=len)) for i in [keys,vals,coms]]
          for k,v,c in zip(keys,vals,coms): csv.writer(f, delimiter='\t').writerow(['# {!s:{}}'.format(k,klen)+'= {!s:{}}'.format(v,vlen)+' / {!s:{}}'.format(c,clen)])
          csv.writer(f, delimiter='\t').writerow([' '])
      u.dict2txt({str(w):{'flux [{}]'.format(data['flux_units']):str(f), 'unc [{}]'.format(data['flux_units']):str(e)} for w,f,e in zip(data['wavelength'],data['flux'],data['unc'])}, fn, column1='# wavelength [{}]'.format(data['wavelength_units']), append=True)
    else: print "No spectrum found with id {}".format(spectrum_id)

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
