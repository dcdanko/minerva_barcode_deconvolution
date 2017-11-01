import sys


class ProgressBar:

    def __init__(self, total, filelike=sys.stderr, length=100, writeEvery=10):
        self.total = total
        self.events = 0.0
        self.filelike = filelike
        self.length = length
        self.writeEvery = writeEvery
        

        
    def write(self):
        prefix = '{:,} / {:,} '
        prefLength = len(prefix.format(self.total, self.total))
        prefix = prefix.format(self.events, self.total)
        prefix += ' '*(prefLength - len(prefix))
        
        lengthLeft = self.length - prefLength
        lengthLeft -= 2
        p = self.events / self.total
        filled = int(p * lengthLeft)
        bar = '|' + '#'*filled + '-'*(lengthLeft-filled) + '|'
        out = '\r' + prefix + bar
        self.filelike.write(out)

    def finish(self):
        self.write()
        self.filelike.write('\n')

    def increment(self):
        self.events += 1
        if (self.events % self.writeEvery) == 0:
            self.write()
        if self.events == self.total:
            self.finish()
        
