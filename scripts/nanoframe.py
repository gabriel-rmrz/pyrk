import awkward as awk
from coffea.analysis_objects import JaggedCandidateArray
import uproot
import uproot_methods
from pdb import set_trace
from fnmatch import fnmatch

class NanoFrame():
    'Simple class that provides a lazy interface with the NanoAODs'
    def __init__(self, *infiles, branches = []):
        if all(isinstance(i, dict) for i in infiles):
            self.tts = infiles
            self.keys_ = set(self.tts[0].keys())
            self.dict_like_ = True
        elif all(not isinstance(i, dict) for i in infiles):
            self.ufs = map(uproot.open, infiles)
            self.tts = [i['Events'] for i in self.ufs]
            self.keys_ = set([i.decode() for i in self.tts[0].keys()])
            self.dict_like_ = False
        else:
            raise RuntimeError('Cannot mix files and dicts!')

        if branches:
            self.keys_ = set(
                i for i in self.keys_ 
                if any(fnmatch(i, branch) for branch in branches)
                )
        self.cache_ = set()
        self.used_branches_ = set()
        self.table_ = awk.Table()

    def array(self, key):
        self.used_branches_.add(key)
        return awk.concatenate([
            i.array(key) if not self.dict_like_ else i[key]
            for i in self.tts
        ])

    def __getitem__(self, key):
        if key in self.cache_:
            return self.table_[key]
        elif key in self.keys_:
            ret = self.array(key)
            self.table_[key] = ret
            self.cache_.add(key)
            return self.table_[key]
        else:
            branch = key + '_'
            subset = [k for k in self.keys_ if k.startswith(branch)]
            info = {i.replace(branch, '') : self.array(i) for i in subset}
            counter = 'n' + key
            counts = 0
            if counter in self.keys_:
                counts = self.array(counter)
                for name, branch in info.items():
                    if not (branch.counts == counts).all():
                        raise ValueError(f'Key {name} does not have the right shape')
                info = {i : j.content for i, j in info.items()}
            #check that everything is there to make a p4
            if all(i in info for i in ['pt', 'eta', 'phi', 'mass']): # FIXME! wrong logic
                ret = JaggedCandidateArray.candidatesfromcounts(
                    counts,
                    **info
                )
            else:
                ret = awk.Table(**info)
                if all(i in info for i in ['pt', 'phi']): 
                    ret['p4'] = uproot_methods.TLorentzVectorArray.from_ptetaphi(
                        ret['pt'], 0, ret['phi'], 0
                    )
            
            self.table_[key] = ret
            self.cache_.add(key)
            return self.table_[key]

    def __setitem__(self, key, val):
        self.table_[key] = val

    @property
    def used_branches(self):
        return sorted(list(self.used_branches_))

    @property
    def objects(self):
        'branch groups available'
        return list(
            set(
                i.split('_')[0] for i in self.keys_ 
                if '_' in i
            )
        )
        
    @property
    def keys(self):
        'Keys available for loading and usage'
        return list(self.keys_)
    
    @property
    def columns(self):
        'columns already loaded'
        return self.table_.columns
