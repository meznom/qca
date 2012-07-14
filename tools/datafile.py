import numpy

class Datafile:
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

    def get (self, index, block=-1):
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
        ls = self.get(index, block)
        return numpy.loadtxt(self._lines(ls))

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
