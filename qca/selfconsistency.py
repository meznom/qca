from . import __version__
from collections import OrderedDict

class SelfConsistency(object):
    def __init__(self, f=None):
        self.N_max = 1000
        self.epsilon = 1E-8
        self.f = f
        self.args_start = None
        self.alpha = 0.1
        
        self.i = 0
        
        self.program = 'SelfConsistency'
        self.version = __version__
        self.results = OrderedDict()
    
    def init(self):
        pass
    
    def run(self, args_start_=None):
        if args_start_ is not None:
            self.args_start = args_start_
        
        # self consistency loop
        args = self.args_start
        args_new = args
        args_old = args
        for self.i in range(1,self.N_max+1):
            args_new = self.f(args)
            args_old = args
            args = [o + self.alpha*(n-o) for n,o in zip(args_new,args_old)]
            cs = [(abs(n-o) < self.epsilon) for n,o in zip(args_new,args_old)]
            if cs.count(False) == 0:
                break
        
        self.f(args_new)
        self.results['args'] = args_new
        self.results['N_iteration'] = self.i
        return (args,self.i)

    def __getstate__(self):
        i = OrderedDict()
        i['parameters'] = OrderedDict()
        i['parameters']['alpha'] = self.alpha
        i['parameters']['epsilon'] = self.epsilon
        i['parameters']['N_max'] = self.N_max
        i['parameters']['args_start'] = self.args_start
        # If f does not support our serialization protocol, just use a string
        # representation
        if hasattr(self.f, '__getstate__'):
            i['parameters']['f'] = self.f
        else:
            i['parameters']['f'] = str(self.f)
        i['results'] = self.results
        return i

    def __setstate__(self,i):
        self.__init__()
        self.alpha = i['parameters']['alpha']
        self.epsilon = i['parameters']['epsilon']
        self.N_max = i['parameters']['N_max']
        self.args_start = i['parameters']['args_start']
        self.f = i['parameters']['f']
        self.results = i['results']
