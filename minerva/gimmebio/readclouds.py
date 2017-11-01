from .fastx import *

################################################################################

class GemcodeMismatchException( Exception):
    pass

class NoGemcodeException( Exception):
    pass

class ReadCloudException(Exception):
    pass

class ChromiumReadPair( ReadPair):

    __slots__ = ('barcode',)

    def __init__(self, r1, r2):
        super( ChromiumReadPair, self).__init__(r1,r2)
        self.barcode = None
        for tag in r1.tags:
            if 'BX:' in tag:
                self.barcode = tag
                break
        if self.barcode is None:
            raise NoGemcodeException()
        for tag in r2.tags:
            if 'BX:' in tag:
                if tag != self.barcode:
                    raise GemcodeMismatchException()
    @classmethod
    def fromReadPair(cls, readPair):
        return ChromiumReadPair( readPair.r1, readPair.r2)


class ReadCloud:

    __slots__ = ('barcode', 'readPairs')
    
    def __init__(self, barcode, readPairs=[]):
        self.barcode = barcode
        self.readPairs = []
        
        for rP in readPairs:
            self.addPair(rP)
            
    def addPair(self, cRP):
        if type(cRP) != ChromiumReadPair:
            raise TypeError()
        if self.barcode is not None:
            if self.barcode != cRP.barcode:
                raise GemcodeMismatchException()
        else:
            self.barcode = cRP.barcode
        self.readPairs.append(cRP)

    def getSeqs(self):
        out = []
        for rP in self.readPairs:
            out.append(rP.r1.seq)
            out.append(rP.r2.seq)
        return out
        
    def __iter__(self):
        return iter(self.readPairs)

    def __str__(self):
        out = ''
        for rp in self.readPairs:
            out += str(rp) + '\n'
        return out

    def __len__(self):
        return len(self.readPairs)
    
################################################################################

def iterReadClouds( filelike):
    rC = ReadCloud(None)
    for rP in iterFastq( filelike, interleaved=True):
        try:
            cRP = ChromiumReadPair.fromReadPair(rP)
        except  NoGemcodeException as nge:
            continue

        try:
            rC.addPair(cRP)
        except GemcodeMismatchException as gme:
            yield rC
            rC = ReadCloud( cRP.barcode)
            rC.addPair(cRP)
    yield rC
    
