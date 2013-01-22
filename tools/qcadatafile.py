import numpy
import re
import json
import copy
import datetime
from os import walk as oswalk
import fnmatch
import sys

hasMongo = False
try:
    from pymongo import MongoClient
    hasMongo = True
except ImportError:
    pass

# TODO: change date, so it conforms to RFC 2822 (email), e.g. 'Thu, 28 Jun 2001
#       14:17:15 +0000', or use '2013-01-11T19:28:20Z'; and use UTC

class QcaDatafileError (Exception):
    def __init__ (self, value):
        self.value = value
    def __str__ (self):
        return repr(self.value)

class QcaDatafile:
    EMPTY_DOCUMENT = {
            'info': {
                'program': None, 
                'version': None, 
                'date': None,
                'dataFormat': None
            },
            'params': None,
            'data': None
        }

    def __init__ (self, filename):
        self._filename = filename
        self._indices = []
        self.parse(filename)
    
    def parse (self, filename):
        f = open(filename)
        empty = 0
        fpos = 0
        fpos_end = 0
        ln = 1
        ln_end = 0
        ids = self._indices
        # start first index and first block
        ids.append([[(ln, fpos)]])
        
        for l in f:
            if l.strip(' \t') == '\n':
                if empty == 0:
                    ln_end = ln
                    fpos_end = fpos
                empty += 1
            
            else:
                cindex = ids[-1]
                cblock = cindex[-1]
                if empty > 0:
                    if (len(cblock) < 2):
                        # no data => begin of data = end of block
                        cblock.append((ln_end, fpos_end))
                    # end of current block
                    cblock.append((ln_end, fpos_end))
                    if empty == 1:
                        # start new block
                        cindex.append([(ln, fpos)])
                    elif empty >= 2:
                        # start new index and new block
                        ids.append([[(ln, fpos)]])
                    if l[0] != '#':
                        # no header for new block
                        ids[-1][-1].append((ln, fpos))
                elif empty == 0 and l[0] != '#' and len(cblock) < 2:
                    # end of header and start of data
                    cblock.append((ln, fpos))
                empty = 0
            
            ln += 1
            fpos += len(l)

        # close last block
        cindex = ids[-1]
        cblock = cindex[-1]
        if (len(cblock) < 2):
            cblock.append(cblock[0])
        cblock.append((ln, fpos))
        f.close()

    def info (self):
        ids = self._indices
        for i, index in enumerate(ids):
            print "Index", i
            print "   " + str(index[-1][-1][0] - index[0][0][0]) + " lines"
            print
            for j, block in enumerate(index):
                print "   Block", j
                print "      " + str(block[-1][0] - block[0][0]) + " lines"
                print "      Begin Header", block[0]
                print "      Begin Data", block[1]
                print "      End  ", block[2]
            print

    def getAsString (self, index, block=-1):
        start = 0
        end = 0
        if block == -1:
            start = self._indices[index][0][0][1]
            end = self._indices[index][-1][-1][1]
        else:
            start = self._indices[index][block][0][1]
            end = self._indices[index][block][-1][1]
        return self._readfile(start, end)

    def getAsMatrix (self, index, block=-1):
        ls = self.getAsString(index, block)
        return numpy.loadtxt(self._lines(ls))

    def getHeader (self, index):
        # read header from file
        start = self._indices[index][0][0][1]
        end = self._indices[index][0][1][1]
        header = self._readfile(start, end)
        # remove leading pound character
        header = re.sub(r'^# ?', '', header, 0, re.MULTILINE)
        # parse header
        version = None
        date = None
        params = None
        tableHeaders = []
        m = re.search(r'^program version: (.+)$', header, re.MULTILINE)
        if m: version = m.group(1)
        m = re.search(r'^date: (.+)$', header, re.MULTILINE)
        if m: date = m.group(1)
        m = re.search(r'^{.*^}', header, re.MULTILINE | re.DOTALL)
        if m: 
            paramsJson = m.group(0)
            params = json.loads(paramsJson)
            params = self._updateLeaves(params, self._convert)
        m = re.search(r'\n?.+$', header)
        if m:
            lastLine = m.group(0)
            tableHeaders = lastLine.split()
        return {
                'version': version,
                'date': date,
                'params': params,
                'tableHeaders': tableHeaders
               };

    def getBody (self, index):
        # read body from file
        # we ignore the blocks and just read the data of all blocks as one body
        start = self._indices[index][0][1][1]
        end = self._indices[index][-1][2][1]
        body = self._readfile(start, end)
        # remove comments, leading white spaces, and empty lines
        body = re.sub(r'^#.+$\n?', '', body, 0, re.MULTILINE)
        body = re.sub(r'^\s*', '', body, 0, re.MULTILINE)
        body = re.sub(r'\n\n', '\n', body, 0, re.MULTILINE)
        # parse body
        rows = []
        for l in body.splitlines():
            cols = [float(n) for n in l.split()]
            rows.append(cols)
        return rows

    def get (self, index):
        h = self.getHeader(index)
        b = self.getBody(index)
        if len(b) < 1:
            raise QcaDatafileError('Body is empty')
        if len(h['tableHeaders']) != len(b[0]):
            raise QcaDatafileError('Number of table headers does not match number ' +
                    'of table columns (Index {:d})'.format(index))
        o = copy.deepcopy(self.EMPTY_DOCUMENT)
        o['info']['version'] = h['version']
        o['info']['date'] = h['date']
        o['params'] = h['params']
        o['data'] = {'headers': h['tableHeaders'], 'rows': b}
        if not (h['version'] is None and 
                h['date'] is None and 
                h['params'] is None):
            o['info']['program'] = 'runQca'
            o['info']['dataFormat'] = 'QcaDataFileV1'
        return o

    def getAll (self):
        os = []
        l = copy.deepcopy(self.EMPTY_DOCUMENT)
        for i, index in enumerate(self._indices):
            o = self.get(i)
            if (o['info']['date'] is None and 
                o['info']['version'] is None and
                o['params'] is None):
                    o['info'] = copy.deepcopy(l['info'])
                    o['params'] = copy.deepcopy(l['params'])
            else:
                l['info'] = copy.deepcopy(o['info'])
                l['params'] = copy.deepcopy(o['params'])
            os.append(o)
        return os

    def toJson (self, prettyPrint=True):
        if prettyPrint:
            return json.dumps(self.getAll(), indent=2, separators=(',', ': '))
        else:
            return json.dumps(self.getAll())

    def _readfile (self, start, end):
        f = open(self._filename)
        f.seek(start)
        ls = ""
        while True:
            l = f.readline()
            if l != "" and f.tell() <= end:
                ls += l
            else:
                break
        f.close()
        return ls
    
    def _lines (self, lsString):
        ls = lsString.splitlines()
        for l in ls:
            yield l

    def _convert (self, s):
        try: 
            return int(s)
        except ValueError:
            try:
                return float(s)
            except ValueError:
                return s

    def _updateLeaves (self, d, f):
        if isinstance(d, dict):
            n = {}
            for k, v in d.iteritems():
                n[k] = self._updateLeaves(v, f)
            return n
        elif isinstance(d, list):
            n = []
            for v in d:
                n.append(self._updateLeaves(v, f))
            return n
        else:
            return f(d)

class PrintJson:
    def __init__ (self, prettyPrint=True):
        self._prettyPrint = prettyPrint
        self._count = 0

    def __call__ (self, o):
        if self._count > 0:
            print(',')
        if self._prettyPrint:
            print(json.dumps(o, indent=2, separators=(',', ': ')))
        else:
            print(json.dumps(o))
        self._count += 1

class ObjectIntoMongo:
    def __init__ (self):
        c = MongoClient()
        self._db = c.physics

    def __call__ (self, o):
        # Convert date string to a Datetime object
        # We ignore the timezone information, e.g. -0600
        if o['info']['date'] is not None:
            dateString = o['info']['date']
            date = datetime.datetime.strptime(dateString[:-6], "%a %b %d %H:%M:%S %Y")
            o['info']['date'] = date
        self._db.qca.insert(o)

def processFiles (directory, function):
    for d, ds, fs in oswalk(directory):
        for f in fnmatch.filter(fs, '*.dat'):
            sys.stderr.write('Processing ' + d + '/' + f)
            try:
                df = QcaDatafile(d + '/' + f)
                os = df.getAll()
                for o in os:
                    function(o)
                sys.stderr.write(' [Ok]\n')
            except QcaDatafileError as e:
                sys.stderr.write(' [Failed]\n')
                sys.stderr.write('    ' + str(e) + '\n')

def filesToJson (directory, prettyPrint=True):
    print('[')
    processFiles(directory, PrintJson(prettyPrint))
    print(']')

def filesToMongo (directory):
    processFiles(directory, ObjectIntoMongo())

if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == 'toJson':
        d = sys.argv[2]
        filesToJson(d)
    elif len(sys.argv) == 3 and sys.argv[1] == 'toMongo':
        if not hasMongo:
            sys.stderr.write('Sorry, PyMongo is not installed.\n')
            sys.exit(-1)
        d = sys.argv[2]
        filesToMongo(d)
    else:
        sys.stderr.write('Usage: python ' + sys.argv[0] + ' toJson|toMongo\n')