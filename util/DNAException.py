
class ConfigFileException(Exception):
    def __init__(self, info):
        self.info = info

    def __str__(self):
        return self.info
    

class EncodeException(Exception):
    def __init__(self, info):
        self.info = info

    def __str__(self):
        return self.info
    
class DecodeException(Exception):
    def __init__(self, info):
        self.info = info

    def __str__(self):
        return self.info