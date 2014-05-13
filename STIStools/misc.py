import os

def IRAF2filepath(IRAFstring, oref='.'):
    '''Change format of  path to reference file from file header

    Parameters
    ----------
    IRAFstring : string
        value of keyword in fits header
    oref : string
        file path to reference files (in IRAF that is the environment
        variable ``oref``)

    Returns
    -------
    filepath : string
        Path and file name
    '''
    return os.path.join(*IRAFstring.replace('oref', oref).split('$'))
