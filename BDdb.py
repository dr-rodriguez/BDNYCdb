#!/usr/bin/python
# BDNYC database
import io, os, itertools, warnings, sqlite3 as sql, numpy as np, matplotlib.pyplot as plt, astropy.io.fits as pf, utilities as u, astrotools as a
warnings.simplefilter('ignore')

class get_db:
  def __init__(self, dbpath):
    '''
    Initialize the database at **dbpath**.  
    '''
    if os.path.isfile(dbpath):
      con = sql.connect(dbpath, isolation_level=None, detect_types=sql.PARSE_DECLTYPES)
      con.text_factory = str
      self.modify = con
      self.query = con.cursor()
      self.dict = con.cursor()
      self.dict.row_factory = sql.Row
    else: print "Sorry, no such file '{}'".format(dbpath)
        
  def add_data(self, CSV, table, multiband=False):
    '''
    Adds data in **CSV** file to the specified database **table**. Note column names (row 1 of CSV file) must match table fields to insert, however order and completeness don't matter.
    '''   
    data, (columns, types), insert, update = u.txt2dict(CSV, all_str=True, delim=',', to_list=False if multiband and table=='photometry' else True, start=0), zip(*self.query.execute("PRAGMA table_info({})".format(table)).fetchall())[1:3], [], []
    if multiband and table=='photometry':
      bands = u.get_filters().keys()
      for row in data.keys():
        for i in bands:
          if all([k in data[row].keys() for k in [i,i+'_unc']]):
            values = [data[row].get(col,None) for col in columns]
            values[columns.index('magnitude')], values[columns.index('magnitude_unc')], values[columns.index('band')], values[columns.index('source_id')] = data[row].get(i), data[row].get(i+'_unc'), i, row
            if values[1].isdigit(): update.append(tuple(values)) if values[columns.index('id')] else insert.append(tuple(values))
    else:
      data_columns = data.pop(0)
      for row in data:
        values = [None for i in columns]
        values[0] = sorted(list(set(range(1,self.query.execute("SELECT max(id) FROM {}".format(table)).fetchone()[0]+2))-set(zip(*self.query.execute("SELECT id FROM {}".format(table)).fetchall())[0])))[0]
        for col in columns: values[columns.index(col)] = row[data_columns.index(col)] if col in data_columns and row[data_columns.index(col)] else None
        update.append(tuple(values)) if values[columns.index('id')] else insert.append(tuple(values))
    if insert:
      u.printer(columns, insert, truncate=30, empties=True), self.query.executemany("INSERT INTO {} VALUES({})".format(table, ','.join('?'*len(columns))), insert), self.modify.commit()
      print "{} new records added to the {} table.".format(len(insert),table.upper())
    if update:
      u.printer(['Command','Result'],[['-'*30,'-'*100],['[column name]','Display full record entry for that column without taking action'],['k','Keeps both records and assigns second one new id if necessary'],['r','Replaces all columns of first record with second record values'],['r [column name] [column name]...','Replaces specified columns of first record with second record values'],['c','Complete empty columns of first record with second record values where possible'],['[Enter]','Keep first record and delete second'],['abort','Abort merge of current table, undo all changes, and proceed to next table']], title=' ')
      for item in update:
        record = self.query.execute("SELECT * FROM {} WHERE id={}".format(table, item[columns.index('id')])).fetchone()
        if record: compare_records(self, table, columns, record, item)
    self.clean_up(table)
   
  def add_ascii(self, asciiPath, source_id, snrPath='', header_chars=['#'], skip=[], start=0, wavelength_units='', flux_units='', publication_id='', obs_date='', wavelength_order='', instrument_id='', telescope_id='', mode_id='', airmass=0, comment=''): 
    '''
    Adds an ascii spectrum to the *spectra* table given an **asciiPath**. Any *spectra* table columns besides *wavelength*, *flux*, *unc*, *snr* arrays can be specified as arguments.
    '''
    filename, data = os.path.basename(asciiPath), zip(*u.txt2dict(asciiPath, to_list=True, skip=header_chars+skip))
    wavelength, flux = [np.array(i, dtype='float32') for i in data][start:start+2]
    try:
      snr = np.array(zip(*u.txt2dict(snrPath, to_list=True, start=1))[-1][1:], dtype='float32') if snrPath else ''
      unc = flx/snr if snrPath else np.array(data[2], dtype='float32')
    except: snr = unc = ''
    wave_const = 10000 if wavelength[0]<100 else 1
    regime = 'OPT' if wavelength[0]<0.8*wave_const and wavelength[-1]<1.2*wave_const else 'NIR' if wavelength[0]<1.2*wave_const and wavelength[-1]>2*wave_const else 'MIR' if wavelength[-1]>3*wave_const else None
    try:
      h = [[i.strip().replace('# ','').replace('\n','') for i in j.replace('=',' /').replace('COMMENT','COMMENT / ').split(' /')] for j in open(asciiPath) if any([j.startswith(char) for char in header_chars])]
      for n,i in enumerate(h): 
        if len(i)==1: h.pop(n)
        elif len(i)==2: h[n].append('')
        elif len(i)>=4: h[n] = [h[n][0],h[n][1],' '.join(h[n][2:])] 
      hdu = pf.PrimaryHDU()
      for i in h:
        try: hdu.header.append(tuple(i))
        except ValueError: pass
      header = hdu.header
    except: header = None
    # except IOError: header = pf.PrimaryHDU().header
    try:
      # Get the lowest spec_id from the spectra table
      spec_id = sorted(list(set(range(1,self.query.execute("SELECT max(id) FROM spectra").fetchone()[0]+2))-set(zip(*self.query.execute("SELECT id FROM spectra").fetchall())[0])))[0]
      try: self.query.execute("INSERT INTO spectra VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (spec_id, source_id, wavelength, wavelength_units, flux, flux_units, unc, snr, wavelength_order, regime, publication_id, obs_date, instrument_id, telescope_id, mode_id, airmass, filename, comment, header)), self.modify.commit()
      except: self.query.execute("INSERT INTO spectra VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (spec_id, source_id, wavelength, wavelength_units, flux, flux_units, unc, snr, wavelength_order, regime, publication_id, obs_date, instrument_id, telescope_id, mode_id, airmass, filename, comment, None)), self.modify.commit()
      u.printer(['spec_id','source_id','wavelength_unit','flux_unit','regime','publication_id','obs_date', 'instrument_id', 'telescope_id', 'mode_id', 'airmass', 'filename', 'comment'],[[spec_id,source_id, wavelength_units, flux_units, regime, publication_id, obs_date, instrument_id, telescope_id, mode_id, airmass, filename, comment]], empties=True)
      # self.clean_up('spectra')
    except IOError: 
      print "Couldn't add spectrum to database."

  def add_fits(self, fitsPath, source_id, unc_fitsPath='', wavelength_units='', flux_units='', publication_id='', obs_date='', wavelength_order='', regime='', instrument_id='', telescope_id='', mode_id='', airmass=0, comment='', wlog=False, SDSS=False):
    '''
    Checks the header of the **fitsFile** and inserts the data with **source_id**.
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
        n = header['TELESCOP'].lower() if isinstance(header['TELESCOP'],str) else ''
        telescope_id = 5 if 'hst' in n else 6 if 'spitzer' in n else 7 if 'irtf' in n else 9 if 'keck' in n and 'ii' in n else 8 if 'keck' in n and 'i' in n else 10 if 'kp' in n and '4' in n else 11 if 'kp' in n and '2' in n else 12 if 'bok' in n else 13 if 'mmt' in n else 14 if 'ctio' in n and '1' in n else 15 if 'ctio' in n and '4' in n else 16 if 'gemini' in n and 'north' in n else 17 if 'gemini' in n and 'south' in n else 18 if 'vlt' in n else 19 if '3.5m' in n else 20 if 'subaru' in n else 21 if ('mag' in n and 'ii' in n) or ('clay' in n) else 22 if ('mag' in n and 'i' in n) or ('baade' in n) else None
      except KeyError: telescope_id = ''
    if not instrument_id:
      try: 
        i = header['INSTRUME'].lower()
        instrument_id = 1 if 'r-c spec' in i or 'test' in i or 'nod' in i else 2 if 'gmos-n' in i else 3 if 'gmos-s' in i else 4 if 'fors' in i else 5 if 'lris' in i else 6 if 'spex' in i else 7 if 'ldss3' in i else 8 if 'focas' in i else 9 if 'nirspec' in i else 0
      except KeyError: instrument_id = ''
    try: airmass = header['AIRMASS']
    except: airmass = 0
    
    try:
      if SDSS:
        data = pf.open(fitsPath, memmap=True)[1].data
        flx, wav, err = map(np.array,zip(*data)[:3])
        flx, wav, err, wavelength_units, flux_units = flx*10**-17, 10**wav, np.sqrt(1/err)*10**-17, 'A', 'ergs-1cm-2A-1'
      else:
        data = a.read_spec(fitsPath, errors=True, atomicron=True, negtonan=True, verbose=False, wlog=wlog)[0]
        wav, flx = data[:2]
        try: err = a.read_spec(unc_fitsPath, errors=True, atomicron=True, negtonan=True, verbose=False, wlog=wlog)[0][1] if unc_fitsPath else data[2]
        except: err = ''
      try: snr = flx/err if any(flx) and any(err) else None
      except (TypeError,IndexError): snr = None

      if not regime:
        if wav[0]<500 or wavelength_units=='um': regime = 'OPT' if wav[0]<0.8 and wav[-1]<1.2 else 'NIR' if wav[0]<1.2 and wav[-1]>2 else 'MIR' if wav[-1]>2.5 else None     
        else: regime = 'OPT' if wav[0]<8000 and wav[-1]<12000 else 'NIR' if wav[0]<12000 and wav[-1]>20000 else 'MIR' if wav[-1]>25000 else None     

      spec_id = sorted(list(set(range(1,self.query.execute("SELECT max(id) FROM spectra").fetchone()[0]+2))-set(zip(*self.query.execute("SELECT id FROM spectra").fetchall())[0])))[0]
      self.query.execute("INSERT INTO spectra VALUES(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (spec_id, source_id, wav, wavelength_units, flx, flux_units, err, snr, wavelength_order, regime, publication_id, obs_date, instrument_id, telescope_id, mode_id, airmass, filename, comment, header)), self.modify.commit()
      u.printer(['spec_id','source_id','wavelength_unit','flux_unit','regime','publication_id','obs_date', 'instrument_id', 'telescope_id', 'mode_id', 'airmass', 'filename', 'comment'],[[spec_id,source_id, wavelength_units, flux_units, regime, publication_id, obs_date, instrument_id, telescope_id, mode_id, airmass, filename, comment]], empties=True)
      # self.clean_up('spectra')
    except KeyError: print "Couldn't add fits file {}".format(fitsPath); print [filename, source_id, wavelength_units, flux_units, obs_date, instrument_id, telescope_id, mode_id, airmass, comment]
  
  def clean_up(self, table):
    '''
    Removes exact duplicates, blank records or data without a *source_id* from the specified **table**. Then finds possible duplicates and prompts for conflict resolution.
    '''
    (columns, types), dup, ignore, I = zip(*self.query.execute("PRAGMA table_info({})".format(table)).fetchall())[1:3], 1, [], ''
    
    # Delete blank records, exact duplicates, or data without a source_id
    self.query.execute("DELETE FROM {0} WHERE ({1})".format(table, columns[1]+' IS NULL' if len(columns)==2 else (' IS NULL AND '.join(columns[1:])+' IS NULL')))
    self.query.execute("DELETE FROM {0} WHERE id NOT IN (SELECT min(id) FROM {0} GROUP BY {1})".format(table,', '.join(columns[1:]))), self.modify.commit()
    if 'source_id' in columns: self.query.execute("DELETE FROM {0} WHERE source_id IS NULL OR source_id IN ('null','None','')".format(table))
    
    if table in ['sources','spectra','photometry','spectral_types','radial_velocities','parallaxes','proper_motions']:
      # Define columns to check for uniqueness
      primary, secondary, ignore = ['flux','magnitude','parallax','spectral_type','proper_motion_ra','proper_motion_dec','radial_velocity'], ['band','regime'], []
      
      # Find non-unique records and run BDdb.compare_records()
      while dup:
        # Check SOURCES table for records with similar RA and Dec but different ids that are not companions or system components.
        if table=='sources': dup = self.query.execute(  "SELECT t1.*, t2.* FROM sources t1 JOIN sources t2 WHERE t1.id!=t2.id \
                                                        AND (t1.companions IS NULL OR (t1.companions IS NOT NULL AND t1.companions NOT LIKE '%' || CAST(t2.id AS TEXT) || '%')) AND (t2.companions IS NULL OR (t2.companions IS NOT NULL AND t2.companions NOT LIKE '%' || CAST(t1.id AS TEXT) || '%')) \
                                                        AND (t1.components IS NULL OR (t1.components IS NOT NULL AND t1.components NOT LIKE '%' || CAST(t2.id AS TEXT) || '%')) AND (t2.components IS NULL OR (t2.components IS NOT NULL AND t2.components NOT LIKE '%' || CAST(t1.id AS TEXT) || '%')) \
                                                        AND (t1.ra BETWEEN t2.ra-0.0007 AND t2.ra+0.0007) AND (t1.dec BETWEEN t2.dec-0.00077 AND t2.dec+0.0007)\
                                                        {}".format(' AND '+' AND '.join(['(t1.id NOT IN ({0}) AND t2.id NOT IN ({0}))'.format(','.join(map(str,i))) for i in ignore]) if ignore else '')).fetchone() 

        # Check all other tables for records with identical primary and secondary column values but different ids.        
        else: dup = self.query.execute("SELECT t1.*, t2.* FROM {0} AS t1 JOIN {0} AS t2 ON t1.source_id=t2.source_id WHERE t1.id!=t2.id AND t1.{1}=t2.{1}{2}{3}".format(table, [c for c in columns if c in primary].pop(), ' AND t1.{0}=t2.{0}'.format([c for c in columns if c in secondary].pop()) if [c for c in columns if c in secondary] else '', ' AND '+' AND '.join(['(t1.id NOT IN ({0}) AND t2.id NOT IN ({0}))'.format(','.join(map(str,i))) for i in ignore]) if ignore else '')).fetchone()        
        
        # Compare potential duplicates and prompt user for action on each
        if dup and dup[:len(dup)/2][0]!=dup[len(dup)/2:][0]:
          I = compare_records(self, table, columns, dup[:len(dup)/2], dup[len(dup)/2:], delete=True)
          if isinstance(I,list): ignore.append(I)
          elif I=='undo': pass # Add this functionality!
          elif I=='abort': break
          else: pass
    
    # Finish or abort table merge
    if I=='abort': 
      print '\nAborted merge of {} table. Undoing all changes.\n'.format(table.upper())
      return 'abort'
    else: print 'Finished clean up on {} table.'.format(table.upper())

  def edit_columns(self, table, columns, types):
    '''
    Rearrange, add or delete columns from database **table** with desired ordered list of **columns** and corresponding data **types**.
    '''
    types[0] = 'INTEGER PRIMARY KEY'
    self.query.execute("ALTER TABLE {0} RENAME TO TempOldTable".format(table)), self.query.execute("CREATE TABLE {0} ({1})".format(table, ', '.join(['{} {}'.format(c,t) for c,t in zip(columns,types)]))), self.query.execute("INSERT INTO {0} ({1}) SELECT {1} FROM TempOldTable".format(table, ','.join([c for c in list(zip(*self.query.execute("PRAGMA table_info(TempOldTable)").fetchall())[1]) if c in columns]))), self.query.execute("DROP TABLE TempOldTable")
    
  def header(self, spectrum_id_or_path):
    '''
    Prints the header information for the given **spectrum_id_or_path**.
    '''
    if isinstance(spectrum_id_or_path,int):
      try: 
        H = self.dict.execute("SELECT * FROM spectra WHERE id={}".format(spectrum_id_or_path)).fetchone()['header']
        if H: return H
        else: print 'No header for spectrum {}'.format(spectrum_id_or_path)
      except TypeError: print 'No spectrum with id {}'.format(spectrum_id_or_path)
    elif os.path.isfile(spectrum_id_or_path):
      if spectrum_id_or_path.endswith('.fits'):
        return pf.getheader(spectrum_id_or_path)
      else:
        txt, H = open(spectrum_id_or_path), []
        for i in txt: 
          if i.startswith('#'): H.append(i)
        txt.close()
        print ''.join(H) if H else 'No header for spectrum {}'.format(spectrum_id_or_path)
    else: print 'No such file {}'.format(spectrum_id_or_path)

  def identify(self, search):
    '''
    For **search** input of (ra,dec) decimal degree tuple, i.e. '(12.3456,-65.4321)', returns all sources within 1 arcminute.
    For **search** input of text string, i.e. 'vb10', returns all sources with case-insensitive partial text matches in *names* or *designation* columns.
    '''
    try: q = "SELECT id,ra,dec,designation,unum,shortname,names FROM sources WHERE ra BETWEEN "+str(search[0]-0.01667)+" AND "+str(search[0]+0.01667)+" AND dec BETWEEN "+str(search[1]-0.01667)+" AND "+str(search[1]+0.01667)
    except TypeError: q = "SELECT id,ra,dec,designation,unum,shortname,components,companions,names FROM sources WHERE REPLACE(names,' ','') like '%"+search.replace(' ','')+"%' or designation like '%"+search+"%' or unum like '%"+search+"%' or shortname like '%"+search+"%'"
    results = self.query.execute(q).fetchall()
    if results: 
      if len(results)==1: self.inventory(int(results[0][0]))
      else: u.printer(['id','ra','dec','designation','unum','short','components','companions','names'], results, truncate=40, empties=True)
    else: print "No objects found by {}".format(search)
      
  def inventory(self, ID, verbose=True, plot=False, data=False):
    '''
    Prints a summary of all objects in the database. Input string or list of strings in **ID** or **unum** for specific objects.
    '''
    if ID:
      q = "SELECT sources.id, sources.unum, sources.designation, sources.ra, sources.dec, (SELECT COUNT(*) FROM spectra WHERE spectra.source_id=sources.id), (SELECT COUNT(*) FROM spectra WHERE spectra.source_id=sources.id AND spectra.regime='OPT'), (SELECT COUNT(*) FROM spectra WHERE spectra.source_id=sources.id AND spectra.regime='NIR'), (SELECT COUNT(*) FROM photometry WHERE photometry.source_id=sources.id), (SELECT parallax from parallaxes WHERE parallaxes.source_id=sources.id), (SELECT parallax_unc from parallaxes WHERE parallaxes.source_id=sources.id), (SELECT spectral_type FROM spectral_types WHERE spectral_types.source_id=sources.id AND regime='OPT'), (SELECT spectral_type FROM spectral_types WHERE spectral_types.source_id=sources.id AND regime='IR'), (SELECT gravity from spectral_types WHERE spectral_types.source_id=sources.id) FROM sources"
      IDS = ID if isinstance(ID,list) else [ID]
      
      if len(IDS)==1 and verbose:
        for table in ['sources']+[t for t in zip(*self.query.execute("SELECT * FROM sqlite_master WHERE type='table'").fetchall())[1] if t!='sources']:
          columns, types = map(list,zip(*self.query.execute("PRAGMA table_info({})".format(table)).fetchall())[1:3])
          if 'source_id' in columns or table=='sources':
            data = map(list, self.query.execute("SELECT * FROM {} WHERE {}".format(table,'id={}'.format(ID) if table=='sources' else 'source_id={}'.format(ID))).fetchall())
            for d in data+[columns,types]:
              if table!='sources': d.pop(1)
              if table=='spectra': d.pop(1), d.pop(2), d.pop(3), d.pop(3), d.pop(-4)
            if data: u.printer([c.replace('wavelength_units','W').replace('flux_units','F').replace('comment','com').replace('header','head').replace('wavelength_order','ord').replace('lication_id','').replace('rument_id','').replace('escope_id','').replace('mode_id','mode') for c in columns], [['Yes' if t=='HEADER' or c=='comment' and v else str(v)[2:8] if t=='ARRAY' and v is not '' else v for c,t,v in zip(columns,types,d)] for d in data] if table=='spectra' else data, truncate=20 if table=='sources' else 50, title='SOURCE: '+data[0][-3] if table=='sources' else table.upper(), empties=True)
      else:
        try:
          D = self.query.execute(q+' WHERE id IN ({})'.format("'"+"','".join(map(str,IDS))+"'")).fetchall()
          if D:
            u.printer(['id','unum','name','ra','dec','Spec Count','Optical','NIR','Phot Count','Pi','Pi_unc','OPT','IR','grav'], D, empties=True)
            if data: return D
          else: print "No sources found{}.".format(' with id '+str(ID) if ID else '')
        except IndexError: pass
      if plot:
        for I in IDS:
          for i in self.dict.execute("SELECT * FROM spectra WHERE source_id={}".format(I)).fetchall(): self.plot_spectrum(i['id'])

  def merge(self, conflicted, tables=[], diff_only=True):
    '''
    Merges specific **tables** or all tables of **conflicted** databse into the master database.
    '''
    if os.path.isfile(conflicted):
      # Load and attach master and conflicted databases
      con, master, reassign = get_db(conflicted), self.query.execute("PRAGMA database_list").fetchall()[0][2], {}
      con.query.execute("ATTACH DATABASE '{}' AS m".format(master)), self.query.execute("ATTACH DATABASE '{}' AS c".format(conflicted)), con.query.execute("ATTACH DATABASE '{}' AS c".format(conflicted)), self.query.execute("ATTACH DATABASE '{}' AS m".format(master))
      
      # Drop any backup tables from failed merges
      for table in tables: self.query.execute("DROP TABLE IF EXISTS Backup_{0}".format(table))
      
      # Gather user data to add to CHANGELOG table
      import socket, datetime
      user, machine_name, date, modified_tables = raw_input('Please enter your name : '), socket.gethostname(), datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), []
      
      # Print instructions for user
      u.printer(['Command','Result'],[['-'*30,'-'*100],['[column name]','Display full record entry for that column without taking action'],['k','Keeps both records and assigns second one new id if necessary'],['r','Replaces all columns of first record with second record values'],['r [column name] [column name]...','Replaces specified columns of first record with second record values'],['c','Complete empty columns of first record with second record values where possible'],['[Enter]','Keep first record and delete second'],['abort','Abort merge of current table, undo all changes, and proceed to next table']], title=' ')
      
      # Merge table by table, starting with SOURCES
      tables = tables or ['sources']+[t for t in zip(*self.query.execute("SELECT * FROM sqlite_master WHERE name NOT LIKE '%Backup%' AND type='table'{}".format(" AND name IN ({})".format("'"+"','".join(tables)+"'") if tables else '')).fetchall())[1] if t!='sources']
      for table in tables:
        # Get column names and data types from master table and column names from conflicted table
        (columns, types), conflicted_cols = zip(*self.query.execute("PRAGMA table_info({})".format(table)).fetchall())[1:3], zip(*con.query.execute("PRAGMA table_info({})".format(table)).fetchall())[1]
                
        if any([i not in columns for i in conflicted_cols]):
          # Abort table merge if conflicted has new columns not present in master. New columns must be added to the master database first via db.edit_columns().
          print "\nMerge of {0} table aborted since conflicted copy has columns {1} not present in master.\nAdd new columns to master with BDdb.edit_columns() and try again.\n".format(table.upper(),[i for i in conflicted_cols if i not in columns])
        
        else:
          # Add new columns from master table to conflicted table if necessary
          if any([i not in conflicted_cols for i in columns]): con.query.execute("DROP TABLE IF EXISTS Conflicted_{0}".format(table)), con.query.execute("ALTER TABLE {0} RENAME TO Conflicted_{0}".format(table)), con.query.execute("CREATE TABLE {0} ({1})".format(table, ', '.join(['{} {}'.format(c,t) for c,t in zip(columns,types)]))), con.query.execute("INSERT INTO {0} ({1}) SELECT {1} FROM Conflicted_{0}".format(table, ','.join(conflicted_cols))), con.query.execute("DROP TABLE Conflicted_{0}".format(table)), con.modify.commit()
        
          # Pull unique records from conflicted table
          data = map(list, con.query.execute("SELECT * FROM (SELECT 1 AS db, {0} FROM m.{2} UNION ALL SELECT 2 AS db, {0} FROM c.{2}) GROUP BY {1} HAVING COUNT(*)=1 AND db=2".format(','.join(columns),','.join(columns[1:]),table)).fetchall())

          if data:
            if diff_only:
              u.printer(columns, [[repr(i) for i in d] for d in data], truncate=15 if table=='spectra' else 20)
            else:
              # Make temporary table copy so changes can be undone at any time
              self.query.execute("DROP TABLE IF EXISTS Backup_{0}".format(table)), self.query.execute("ALTER TABLE {0} RENAME TO Backup_{0}".format(table)), self.query.execute("CREATE TABLE {0} ({1})".format(table, ', '.join(['{} {}'.format(c,t) for c,t in zip(columns,types)]))), self.query.execute("INSERT INTO {0} ({1}) SELECT {1} FROM Backup_{0}".format(table, ','.join(columns)))

              # Create a dictionary of any reassigned ids from merged SOURCES tables and replace applicable source_ids in other tables.
              print "\nMerging {} tables.\n".format(table.upper())
              try: count = self.query.execute("SELECT MAX(id) FROM {}".format(table)).fetchone()[0]+1
              except TypeError: count = 1
              for n,i in enumerate([d[1:] for d in data]):
                if table=='sources': reassign[i[0]] = count
                elif 'source_id' in columns and i[1] in reassign.keys(): i[1] = reassign[i[1]]
                else: pass
                i[0] = count
                data[n] = i
                count += 1
            
              # Insert unique conflicted records into master and run BDdb.clean_up()
              self.query.executemany("INSERT INTO {} VALUES({})".format(table, ','.join(['?' for c in columns])), data), self.modify.commit()
              print "{} records added to {} table at '{}':".format(len(data), table, master)
              u.printer(columns, [[repr(i) for i in d] for d in data], truncate=15 if table=='spectra' else 20)
              abort = self.clean_up(table)
          
              # Undo all changes to table if merge is aborted. Otherwise, push table changes to master.
              if abort: self.query.execute("DROP TABLE {0}".format(table)), self.query.execute("ALTER TABLE Backup_{0} RENAME TO {0}".format(table)), self.modify.commit()
              else: self.query.execute("DROP TABLE Backup_{0}".format(table)), self.modify.commit(), modified_tables.append(table.upper())
          
          else: print "{} tables identical.".format(table.upper())
      
      # Add data to CHANGELOG table
      if not diff_only:
        user_description = raw_input('\nPlease describe the changes made in this merge : ')
        self.query.execute("INSERT INTO changelog VALUES(?, ?, ?, ?, ?, ?, ?)", (None, date, user, machine_name, ', '.join(modified_tables), user_description, os.path.basename(conflicted))), self.modify.commit()
      
      # Finish up and detach
      print "\nMerge complete!" if not diff_only else "\nDiff complete. No changes made to either database."
      con.query.execute("DETACH DATABASE c"), self.query.execute("DETACH DATABASE c"), con.query.execute("DETACH DATABASE m"), self.query.execute("DETACH DATABASE m"), con.modify.close()
    else: print "File '{}' not found!".format(conflicted)
    
  def output_spectrum(self, spectrum_id, filepath):
    '''
    Prints a file of the spectrum with id **spectrum_id** to an ascii file with specified **filepath**.
    '''
    data = self.dict.execute("SELECT * FROM spectra WHERE id={}".format(spectrum_id)).fetchone()
    if data:
      import csv
      fn, header = '{}{}.txt'.format(filepath, data['filename'] or spectrum_id), repr([repr(r) for r in data['header'].ascardlist()]) if data['header'] else None
      if data['header']:
        with open(fn, 'w' ) as f:
          keys, vals, coms = zip(*[eval(l) for l in eval(header)])
          klen, vlen, clen = [len(max(['1' if isinstance(j,bool) and j else '0' if isinstance(j,bool) else str(j) for j in i], key=len)) for i in [keys,vals,coms]]
          for k,v,c in zip(keys,vals,coms): csv.writer(f, delimiter='\t').writerow(['# {!s:{}}'.format(k,klen)+'= {!s:{}}'.format(v,vlen)+' / {!s:{}}'.format(c,clen)])
          csv.writer(f, delimiter='\t').writerow([' '])
      u.dict2txt({str(w):{'flux [{}]'.format(data['flux_units']):str(f), 'unc [{}]'.format(data['flux_units']):str(e)} for w,f,e in zip(data['wavelength'],data['flux'],data['unc'])}, fn, column1='# wavelength [{}]'.format(data['wavelength_units']), append=True)
    else: print "No spectrum found with id {}".format(spectrum_id)
  
  def plot_spectrum(self, ID):
    '''
    Plots spectrum **ID** from SPECTRA table.
    '''
    i = self.dict.execute("SELECT * FROM spectra WHERE id={}".format(ID)).fetchone()
    if i:
      try:
        plt.figure(), plt.rc('text', usetex=False), plt.loglog(i['wavelength'], i['flux'], c='b', label='spec_id: {}'.format(i['id'])), plt.grid(True), plt.yscale('log', nonposy='clip'), plt.title('source_id = {}'.format(i['source_id'])), plt.figtext(0.15,0.88, '{}\n{}\n{}\n{}'.format(i['filename'],self.query.execute("SELECT name FROM telescopes WHERE id={}".format(i['telescope_id'])).fetchone()[0] if i['telescope_id'] else '',self.query.execute("SELECT name FROM instruments WHERE id={}".format(i['instrument_id'])).fetchone()[0] if i['instrument_id'] else '',i['obs_date']), verticalalignment='top'), plt.xlabel('[{}]'.format(i['wavelength_units'])), plt.ylabel('[{}]'.format(i['flux_units'])), plt.legend(loc=8, frameon=False)
        X, Y = plt.xlim(), plt.ylim()
        try: plt.fill_between(i['wavelength'], i['flux']-i['unc'], i['flux']+i['unc'], color='b', alpha=0.3), plt.xlim(X), plt.ylim(Y)
        except TypeError: print 'No uncertainty array for spectrum {}'.format(ID)
      except IOError: print "Couldn't print spectrum {}".format(ID)
    else: print "No spectrum {} in the SPECTRA table.".format(ID)

  def lookup(self, table, ids=None, concatenate='', delim='/'):
    '''
    Quickly look up records from the specified *table* and list *ids* to limit results. Specify column values to *concatenate* into a string.
    '''
    if type(ids)==str and table.lower()=='publications' and ',' not in ids:
      return self.query.execute("SELECT * FROM publications WHERE shortname LIKE '%{0}%' OR description LIKE '%{0}%'".format(ids)).fetchall()
    else:
      if type(ids)==int: ids = [ids]
      if type(ids)==str and ',' in ids: ids = ids.split(',')
      if type(ids)==list and not ids: ids = ''
      results = self.query.execute("SELECT {} FROM {} WHERE id IN ({})".format(concatenate or '*',table, ','.join(map(str,ids)))).fetchall() if ids else self.query.execute("SELECT * FROM {}".format(table)).fetchall()
      return '' if not results else delim.join(map(str,zip(*results)[0])) if concatenate else results

# ==============================================================================================================================================
# ================================= Adapters and converters for special data types =============================================================
# ==============================================================================================================================================
  
def adapt_array(arr):
  '''
  Turns an ARRAY string stored in the database back into a Numpy array.
  '''
  out = io.BytesIO()
  np.save(out, arr), out.seek(0)
  return buffer(out.read())

def adapt_header(header): 
  '''
  Turns a HEADER string stored in the database back into a FITS header.
  '''
  return repr([repr(r) for r in header.ascardlist() if type(r[1]) in [int,float,str]]) if hasattr(header,'cards') else None

def convert_array(text):
  '''
  Converts a Numpy array into an ARRAY string to put into the database.
  '''
  out = io.BytesIO(text)
  out.seek(0)
  return np.load(out)

def convert_header(text):
  '''
  Converts a FITS feader into a HEADER string to put into the database.
  '''
  if text:
    hdu = pf.PrimaryHDU()
    for i in [eval(l) for l in eval(text)]:
      if i[0]: 
        try: hdu.header.append(i)
        except: pass
    return hdu.header
  else: return None

sql.register_adapter(np.ndarray, adapt_array), sql.register_adapter(pf.header.Header, adapt_header)
sql.register_converter("ARRAY", convert_array), sql.register_converter("HEADER", convert_header)

# ==============================================================================================================================================
# ================================= Little helper functions ====================================================================================
# ==============================================================================================================================================

def compare_records(db, table, columns, old, new, options=['r','c','k','sql'], delete=False):
  pold, pnew = ['{:.3g}...{:.3g}'.format(i[0],i[-1]) if isinstance(i, np.ndarray) else i for i in old], ['{:.3g}...{:.3g}'.format(i[0],i[-1]) if isinstance(i, np.ndarray) else i for i in new]
  u.printer(columns, [pold,pnew], truncate=20, empties=True, highlight=list(itertools.chain.from_iterable([[(i,n+1) for n,(ov,nv) in enumerate(zip(pold[1:],pnew[1:])) if ov!=nv] for i in range(2)])))
  replace = raw_input("Keep both records [k]? Or replace [r], complete [c], or keep only [Press *Enter*] record {}? (Type column name to inspect or 'help' for options): ".format(old[0]))

  while replace.lower() in columns or replace.lower()=='help':
    if replace.lower() in columns: u.printer(['id',replace.lower()], [[i for idx,i in enumerate(pold) if idx in [0,columns.index(replace.lower())]],[i for idx,i in enumerate(pnew) if idx in [0,columns.index(replace.lower())]]], empties=True)    
    elif replace.lower()=='help': u.printer(['Command','Result'],[['-'*30,'-'*100],['[column name]','Display full record entry for that column without taking action'],['k','Keep both records and assign second one new id if necessary'],['r','Replace all columns of first record with second record values'],['r [column name] [column name]...','Replace specified columns of first record with second record values'],['c','Complete empty columns of first record with second record values where possible'],['[Enter]','Keep first record and delete second'],['abort','Abort merge of current table, undo all changes, and proceed to next table']], title=' ')
    replace = raw_input("Keep both records [k]? Or replace [r], complete [c], or keep only [Press *Enter*] record {}? (Type column name to inspect or 'help' for options): ".format(old[0]))

  if replace and (all([i in list(columns)+options for i in replace.lower().split()]) or replace.lower().startswith('sql')):
    if replace.lower().startswith('r') and 'r' in options:
      if replace.lower()=='r':
        sure = raw_input('Are you sure you want to replace record {} with record {}? [y/n] : '.format(old[0],new[0]))
        if sure.lower()=='y':
          empty_cols, new_vals = zip(*[['{}=?'.format(e),n] for e,n in zip(columns[1:],new[1:])])
          if delete: db.query.execute("DELETE FROM {} WHERE id={}".format(table, new[0]))
          db.query.execute("UPDATE {} SET {} WHERE id={}".format(table, ','.join(empty_cols), old[0]), tuple(new_vals)), db.modify.commit()
      elif all([i in list(columns)+options for i in replace.lower().split()]):
        empty_cols, new_vals = zip(*[['{}=?'.format(e),n] for e,n in zip(columns[1:],new[1:]) if e in replace])
        if empty_cols:
          if delete: db.query.execute("DELETE FROM {} WHERE id={}".format(table, new[0]))
          db.query.execute("UPDATE {} SET {} WHERE id={}".format(table, ','.join(empty_cols), old[0]), tuple(new_vals)), db.modify.commit()
    elif replace.lower()=='c' and 'c' in options:
      try:
        empty_cols, new_vals = zip(*[['{}=?'.format(e),n] for e,o,n in zip(columns[1:],old[1:],new[1:]) if repr(o).lower() in ['','none','null'] and repr(n).lower() not in ['','none','null']])
        if delete: db.query.execute("DELETE FROM {} WHERE id={}".format(table, new[0]))
        db.query.execute("UPDATE {} SET {} WHERE id={}".format(table, ','.join(empty_cols), old[0]), tuple(new_vals)), db.modify.commit()
      except:
        if delete: db.query.execute("DELETE FROM {} WHERE id={}".format(table, new[0])), db.modify.commit()
        else: pass
    elif replace.lower()=='k' and 'k' in options: return [old[0],new[0]]
    elif replace.lower().startswith('sql ') and 'sql' in options: 
      try: db.query.execute(replace[4:]), db.modify.commit()
      except: pass 
  elif replace.lower()=='abort': return 'abort'
  elif not replace:
    if delete: db.query.execute("DELETE FROM {} WHERE id={}".format(table, new[0])), db.modify.commit()
    else: pass

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
