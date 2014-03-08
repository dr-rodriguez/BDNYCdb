#!/usr/bin/python
# BDNYC database
import io, os, itertools, sqlite3 as sql, numpy as np, matplotlib.pyplot as plt, astropy.io.fits as pf, utilities as u, astrotools as a

class get_db:
  def __init__(self, dbpath):
    if os.path.isfile(dbpath):
      con = sql.connect(dbpath, isolation_level=None, detect_types=sql.PARSE_DECLTYPES)
      con.text_factory = str
      self.modify = con
      self.query = con.cursor()
      self.dict = con.cursor()
      self.dict.row_factory = sql.Row
    else: print "Sorry, no such file '{}'".format(dbpath)
        
  def add_data(self, CSV, table):
    '''
    Adds data in *CSV* file to the specified database *table*. Note column names (row 1 of CSV file) must match table fields to insert though order and completeness don't matter.
    '''    
    data, (fields, types) = u.txt2dict(CSV, all_str=True, delim=',', to_list=True), zip(*self.query.execute("PRAGMA table_info({})".format(table)).fetchall())[1:3]
    columns, insert, update = data.pop(0), [], []
    for row in data:
      values = [None for i in fields]
      for field in fields: values[fields.index(field)] = row[columns.index(field)] if field in columns and row[columns.index(field)] else None
      update.append(tuple(values)) if values[fields.index('id')] else insert.append(tuple(values))
    if insert:
      u.printer(fields, insert, truncate=15)
      self.query.executemany("INSERT INTO {} VALUES({})".format(table, ','.join('?'*len(fields))), insert), self.modify.commit()
      print "{} new records added to the {} table.".format(len(insert),table.upper())
    if update:
      for item in update:
        record = self.query.execute("SELECT * FROM {} WHERE id={}".format(table, item[fields.index('id')])).fetchone()
        if any([str(i)!=str(r) and t!='ARRAY' and r and not i for i,r,t in zip(item,record,types)]):
          u.printer(['']+list([d for d,t in zip(fields,types) if t!='ARRAY']), [['CURRENT']+list([d for d,t in zip(record,types) if t!='ARRAY']),['NEW_REC']+list([d for d,t in zip(item,types) if t!='ARRAY'])], truncate=15)
          replace = raw_input("Replace [r], complete [c], or keep [Press *Enter*] current record? : ")
          if replace.lower()=='r': self.query.execute("UPDATE {} SET {} WHERE id={}".format(table, ", ".join(["=".join([i,"'"+str(j)+"'" if j and t=='TEXT' else j if j and t in ['INTEGER','NUMERIC'] else 'NULL']) for i,j,t in zip(fields,item,types) if t!='ARRAY']), item[fields.index('id')])), self.modify.commit()
          elif replace.lower()=='c' and any([i and not r for i,r in zip(item,record)]): self.query.execute("UPDATE {} SET {} WHERE id={}".format(table, ", ".join(["=".join([i,"'"+str(j)+"'" if t=='TEXT' else str(j) if t in ['INTEGER','NUMERIC'] else 'NULL']) for i,j,t,r in zip(fields,item,types,record) if t!='ARRAY' and j and not r]), item[fields.index('id')])), self.modify.commit()
          else: pass
    
    self.clean_up(table)
   
  def add_ascii(self, asciiPath, snrPath='', header_chars=['#'], skip=[], start=0, source_id='', unum='', wavelength_units='', flux_units='', publication_id='', obs_date='', wavelength_order='', instrument_id='', telescope_id='', airmass=0, comment=''): 
    filename, data = os.path.basename(asciiPath), zip(*u.txt2dict(asciiPath, to_list=True, skip=header_chars+skip))
    wavelength, flux = [np.array(i, dtype='float32') for i in data][start:start+2]
    try:
      snr = np.array(zip(*u.txt2dict(snrPath, to_list=True, start=1))[1][1:], dtype='float32') if snrPath else ''
      unc = flx/snr if snrPath else np.array(data[2], dtype='float32')
    except: snr = unc = ''
    regime = 'OPT' if wavelength[0]<0.8 and wavelength[-1]<1.2 else 'NIR' if wavelength[0]<1.2 and wavelength[-1]>2 else 'MIR' if wavelength[-1]>3 else None
    try:
      hdu, h = pf.PrimaryHDU(), [i for i in pf.open(asciiPath) if any([i.startswith(char) for char in header_chars])]
      hdu.header.append(('COMMENT',' '.join([' '.join(i) for n,i in enumerate(h) if any([h[n][0].startswith(char) for char in header_chars])]),''))
      header = hdu.header
    except: header = ''
    if not source_id and unum: source_id = self.dict.execute("SELECT * FROM sources WHERE unum='{}'".format(unum)).fetchone()['id']
    try:
      self.query.execute("INSERT INTO spectra VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (None, source_id, wavelength, wavelength_units, flux, flux_units, unc, snr, wavelength_order, regime, publication_id, obs_date, instrument_id, telescope_id, airmass, filename, comment, header)), self.modify.commit()
      u.printer(['source_id','wavelength_unit','flux_unit','regime','publication_id','obs_date', 'instrument_id', 'telescope_id', 'airmass', 'filename', 'comment'],[[source_id, wavelength_units, flux_units, regime, publication_id, obs_date, instrument_id, telescope_id, airmass, filename, comment]])
      self.clean_up('spectra')
    except: 
    	print "Couldn't add spectrum to database."

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
      try: source_id = {str(k):j for j,k in [tuple(i) for i in self.query.execute("SELECT id, unum FROM sources").fetchall()]}[unum]
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
        except KeyError:
          try: obs_date = header['DATE']
          except KeyError: obs_date = ''
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
      self.clean_up('spectra')
    except: print "Couldn't add fits file {}".format(fitsPath); print [filename, source_id, wavelength_units, flux_units, obs_date, instrument_id, telescope_id, airmass, comment]
  
  def clean_up(self, table):
    '''
    Removes exact duplicates from the specified *table* keeping the record with the lowest id. Then find duplicates and prompt for conflict resolution.
    '''
    columns = zip(*self.query.execute("PRAGMA table_info({})".format(table)).fetchall())[1]
    self.query.execute("DELETE FROM {0} WHERE {1}".format(table, column[1]+' IS NULL' if len(columns)==2 else (' IS NULL AND '.join(columns[1:])+' IS NULL')))
    self.query.execute("DELETE FROM {0} WHERE id NOT IN (SELECT min(id) FROM {0} GROUP BY {1})".format(table,', '.join(columns[1:]))), self.modify.commit()
    
    if table!='spectra': # Until I fix this SQL query!
      # Print out similar records and prompt to delete possible duplicates.
      dups = list(set([item for sublist in [[i[:len(i)/2],i[len(i)/2:]] for i in self.query.execute("SELECT t1.*, t2.* FROM {0} t1 JOIN {0} t2 ON {1} WHERE ({2}) BETWEEN 2 AND {3} AND t1.{4}!=t2.{4}".format(table, ' OR '.join(['t1.{0}=t2.{0}'.format(c) for c in columns[1:]]), '+'.join(["(CASE WHEN t1.{0}=t2.{0} AND t1.{0} IS NOT NULL AND t1.{0} NOT IN ('None','null','NaN','') THEN 1 ELSE 0 END)".format(c) for c in columns[1:]]),len(columns[1:]),columns[0])).fetchall()] for item in sublist]))
      if dups:
        dups = [[repr(i) if isinstance(i, np.ndarray) else i for i in d] for d in dups]
        print '\nPossible duplicates in {} table '.format(table.upper())+'='*(125-len(table))
        u.printer(columns, dups, truncate=15)
        if table=='sources': print "Resolve any issues for the SOURCES table using the SQLite Database Browser."
        else:
          delete = raw_input("Record ids to delete? [Press *Enter* for none] : ")
          if delete:
            if all([i.isdigit() for i in delete]):
              delete = map(int, delete.split(','))
              for ID in delete:
                if ID in map(int, [i[0] for i in dups]):
                  sure = raw_input('Are you sure you want to delete record {}? [y/n] : '.format(ID))
                  if sure.lower()=='y':                  
                    try: self.query.execute("DELETE FROM {} WHERE {}={}".format(table,columns[0],ID)), self.modify.commit()
                    except IOError: print 'Could not delete record {}'.format(ID)
                else: print "Id {} not in list of duplicates, bro.".format(ID)
            else: print 'Just comma separated integers, please! No records deleted.'
  
  def identify(self, search):
    try: q = "SELECT id,ra,dec,designation,unum,shortname,names FROM sources WHERE ra BETWEEN "+str(search[0]-0.01667)+" AND "+str(search[0]+0.01667)+" AND dec BETWEEN "+str(search[1]-0.01667)+" AND "+str(search[1]+0.01667)
    except TypeError: q = "SELECT id,ra,dec,designation,unum,shortname,names FROM sources WHERE names like '%"+search+"%' or designation like '%"+search+"%'"
    u.printer(['id','ra','dec','designation','unum','short','names'], self.query.execute(q).fetchall()) 
      
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
          for I in IDS:
            spectra = self.dict.execute("SELECT s.wavelength,s.flux,s.unc,s.filename,t.name,i.name,s.obs_date,s.source_id FROM spectra AS s JOIN telescopes AS t JOIN instruments as i ON s.instrument_id=i.id AND s.telescope_id=t.id WHERE s.source_id={}".format(I)).fetchall() or self.dict.execute("SELECT wavelength,flux,filename FROM spectra WHERE source_id={}".format(I)).fetchall()
            for i in spectra: plt.figure(), plt.rc('text', usetex=False, fontsize=12), plt.loglog(i[0], i[1], c='b'), plt.grid(True), plt.yscale('log', nonposy='clip'), plt.title('source_id = {}'.format(i[7])), plt.figtext(0.15,0.88, '{}\n{}\n{}\n{}'.format(i[3].replace('_',' '),i[4],i[5],i[6]), verticalalignment='top')
        if data: return D
      else: print "No sources found{}.".format(' with id '+str(ID) if ID else '')
    except: pass

  def merge(self, conflicted):
    if os.path.isfile(conflicted):
      con, master = get_db(conflicted), self.query.execute("PRAGMA database_list").fetchall()[0][2]
      con.query.execute("ATTACH DATABASE '{}' AS m".format(master)), self.query.execute("ATTACH DATABASE '{}' AS c".format(conflicted)), con.query.execute("ATTACH DATABASE '{}' AS c".format(conflicted)), self.query.execute("ATTACH DATABASE '{}' AS m".format(master))
            
      for table in zip(*self.query.execute("SELECT * FROM sqlite_master WHERE type='table'").fetchall())[1]:
        columns, records = zip(*self.query.execute("PRAGMA table_info({})".format(table)).fetchall())[1], self.query.execute("SELECT Count(*) FROM {}".format(table)).fetchone()[0]
        data = map(list, con.query.execute("SELECT * FROM (SELECT 1 AS db, {0} FROM m.{1} UNION ALL SELECT 2 AS db, {0} FROM c.{1}) GROUP BY {0} HAVING COUNT(*)=1 AND db=2".format(','.join(columns[1:]),table)).fetchall())

        if data:
          print "Merging table {}".format(table.upper())
          for i in data: i[0] = None
          self.query.executemany("INSERT INTO {} VALUES({})".format(table, ','.join(['?' for c in columns])), data), self.modify.commit()
          print "{} records added to {}:".format(len(data), master)
          u.printer(columns, [[repr(i) for i in d] for d in data], truncate=15)
          self.clean_up(table)
        else: print "{} table identical.".format(table.upper())
      
      con.query.execute("DETACH DATABASE c"), self.query.execute("DETACH DATABASE c"), con.query.execute("DETACH DATABASE m"), self.query.execute("DETACH DATABASE m"), con.modify.close()
    else: print "File '{}' not found!".format(conflicted)
    
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
