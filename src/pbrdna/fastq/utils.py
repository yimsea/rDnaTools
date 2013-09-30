import numpy as np

from pbcore.io.FastqIO import FastqRecord

def meanP( record ):
    try:
        assert isinstance(record, FastqRecord)
    except:
        raise TypeError("Record is not a FastqRecord!")
    pValues = [10**-(i/10) for i in np.float32(record.quality)]
    return sum(pValues) / len(pValues)

def meanPQv( record ):
    try:
        assert isinstance(record, FastqRecord)
    except:
        raise TypeError("Record is not a FastqRecord!")
    return pValueToQv( meanP(record) )

def pValueToQv( pValue ):
    return -10 * np.log10( pValue )
