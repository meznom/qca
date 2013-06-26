import _qca
from collections import OrderedDict

class Layout(_qca.Layout):
    def __getstate__(self):
        i = OrderedDict()
        i['r_sites'] = self.r_sites
        i['r_charges'] = self.r_charges
        i['charges'] = self.charges
        i['epc'] = self.epc
        return i
