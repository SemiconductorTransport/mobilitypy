from .database import database
import numpy as np

## ==============================================================================
class _AlloyParams:
    def __init__(self, compositions=None, binaries=['AlN', 'GaN'], alloy='AlGaN'):
        self.comps_ = compositions
        self.bins_ = binaries
        self.alloy_ = alloy

    def _get_ternary_params(self):
        assert len(self.bins_) == 2, 'Provide two binary compounds'
        bin_1_params_db = database.get(self.bins_[0])
        bin_2_params_db = database.get(self.bins_[1])
        alloy_params_db = database.get(self.alloy_)
        self.alloy_params_ = {}
        for key, bowing in alloy_params_db.items():
            self.alloy_params_[key] = self.comps_ * bin_1_params_db.get(key) +\
            (1-self.comps_) * bin_2_params_db.get(key) - bowing*self.comps_*(1-self.comps_)
        #print(self.alloy_params_)
            
    def _get_alloy_params(self, system='ternary'):
        if self.comps_ is None:
            self.comps_ = np.linspace(0.01, 0.99, 101)
        elif isinstance(self.comps_, float) or isinstance(self.comps_, int):
            self.comps_ = np.array([self.comps_])
        if system == 'ternary':
            self._get_ternary_params()