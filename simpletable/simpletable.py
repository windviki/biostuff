"""
SimpleTable: simple wrapper around `pytables`_ hdf5
------------------------------------------------------------------------------

.. _`pytables`: http://pytables.org

This module removes some of the boiler-plate code required to use the excellent `pytables`_
module to save and access structured data.

Example Usage::

  >>> from simpletable import SimpleTable
  >>> import tables

define a table as a subclass of simple table.

  >>> RGB = tables.Enum(list('RGB'))
  >>> class ATable(SimpleTable):
  ...     x = tables.Float32Col()
  ...     y = tables.Float32Col()
  ...     name = tables.StringCol(16)
  ...     color = tables.EnumCol(RGB, 'R', 'uint8')

instantiate with: args: filename, tablename

  >>> tbl = ATable('test_docs.h5', 'atable1')

insert as with pytables:

  >>> row = tbl.row
  >>> for i in range(50):
  ...    row['x'], row['y'] = i, i * 10
  ...    row['name'] = "name_%i" % i
  ...    # NOTE how we have to manually translate the enum column.
  ...    row['color'] = RGB['G']
  ...    row.append()
  >>> tbl.flush()

can have the enum cols automatically translated using `insert`

  >>> data = {'x': 1000, 'y': 2000, 'color': 'G', 'name': 'flintstone'}
  >>> tbl.insert(data, row)
  >>> row.append()
  >>> tbl.flush()

there is also `insert_many()` method with takes an iterable
of dicts with keys matching the colunns (x, y, name) in this
case.

query the data (query() alias of tables' readWhere()
note that pytables sends back the data with enum cols as they were
and does nothing to translate them to their original values.

  >>> tbl.query('(x > 4) & (y < 70)') #doctest: +NORMALIZE_WHITESPACE
  array([(1, 'name_5', 5.0, 50.0), (1, 'name_6', 6.0, 60.0)],
         dtype=[('color', '|u1'), ('name', '|S16'), ('x', '<f4'), ('y', '<f4')])

get translated enumcols in an iterator with the .q  method.

  >>> r = tbl.q('x == 1000') # doctest: +NORMALIZE_WHITESPACE
  >>> r # doctest: +ELLIPSIS
  <generator ...>

  >>> list(r)
  [{'color': 'G', 'x': 1000.0, 'name': 'flintstone', 'y': 2000.0}]

or use the `translate_enum` method

  >>> for row_with_enum in tbl.query('(x > 4) & (y < 70)'):
  ...     tbl.translate_enum(row_with_enum)
  {'color': 'G', 'x': 5.0, 'name': 'name_5', 'y': 50.0}
  {'color': 'G', 'x': 6.0, 'name': 'name_6', 'y': 60.0}

Note that using `q` or `translate_enum` will affect performance.
"""

import tables
_filter = tables.Filters(complib="lzo", complevel=1, shuffle=True)

class SimpleTable(tables.Table):
    def __init__(self, file_name, table_name, description=None, 
                 group_name='default', mode='a', title="", filters=_filter, 
                 expectedrows=512000):

        f = tables.openFile(file_name, mode)
        self.uservars = None

        if group_name is None: group_name = 'default'
        parentNode = f._getOrCreatePath('/' + group_name, True)

        if table_name in parentNode: # existing table
            description = None
        elif description is None: # pull the description from the attrs
            description = dict(self._get_description())

        tables.Table.__init__(self, parentNode, table_name, 
                       description=description, title=title,
                       filters=filters,
                       expectedrows=expectedrows,
                       _log=False)
        self._c_classId = self.__class__.__name__
        self.enums = dict([(k, self.getEnum(k)) for k, v in self.coldescrs.iteritems() if isinstance(v, tables.EnumCol)])

    def _get_description(self):
        # pull the description from the attrs
        for attr_name in dir(self):
            if attr_name[0] == '_': continue
            try:
                attr = getattr(self, attr_name)
            except:
                continue
            if isinstance(attr, tables.Atom):
                yield attr_name, attr

    def insert_many(self, data_generator, attr=False):
        row = self.row
        cols = self.colnames
        if not attr:
            for d in data_generator:
                for c in cols:
                    row[c] = d[c]
                row.append()
        else:
            for d in data_generator:
                for c in cols:
                    row[c] = getattr(d, c)
                row.append()
        self.flush()
    
    def insert(self, data, row, attr=False):
        for col in self.colnames:
            if col in self.enums:
                row[col] = self.enums[col][data[col]]
            else:
                row[col] = data[col]
        
    query = tables.Table.readWhere

    def translate_enum(self, row):
        d = {}
        for col in self.colnames:
            d[col] = row[col]
        for col in self.enums:
            d[col] = self.enums[col](row[col])
        return d

    def q(self, *args):
        qq = self.query(*args)
        for row in qq:
            yield self.translate_enum(row)


# convience sublcass that i use a lot.
class BlastTable(SimpleTable):
      query      = tables.StringCol(5)
      subject    = tables.StringCol(5)

      pctid      = tables.Float32Col()
      hitlen     = tables.UInt16Col()
      nmismatch  = tables.UInt16Col()
      ngaps      = tables.UInt16Col()

      qstart     = tables.UInt32Col()
      qstop      = tables.UInt32Col()
      sstart     = tables.UInt32Col()
      sstop      = tables.UInt32Col()

      evalue     = tables.Float64Col()
      score      = tables.Float32Col()


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    import os
    #os.unlink('test_docs.h5')
