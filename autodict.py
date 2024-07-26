class autodict(dict):
    def __init__(self, *args, **kwargs):
        super(autodict, self).__init__(*args, **kwargs)

    def __getitem__(self, key):
        val = dict.__getitem__(self, key)
        return val

    def __setitem__(self, key, val):
        pass
