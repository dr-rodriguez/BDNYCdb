#!/usr/bin/python
# BDNYC database
import io, sqlite3 as sql, numpy as np, matplotlib.pyplot as plt, pyfits as pf, utilities as u

class get_db:
  def __init__(self, dbpath, rows=False):
    con = sql.connect(dbpath, isolation_level=None, detect_types=sql.PARSE_DECLTYPES)
    con.text_factory = str
    if rows: con.row_factory = sql.Row
    self.modify = con
    self.query = con.cursor()
    
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
  
  def clean_up(self, table):
    '''
    Removes exact duplicates from the specified *table* keeping the record with the lowest id.
    '''
    query = "DELETE FROM {0} WHERE id NOT IN (SELECT min(id) FROM {0} GROUP BY {1})".format(table,', '.join(zip(*self.query.execute("PRAGMA table_info({})".format(table)).fetchall())[1][1:]))  
    self.query.execute(query), self.modify.commit()
  
  def inventory(self, ID='', with_pi=True, plot=True):
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
# ================================= Adapters and converters for special data types ====================================================
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