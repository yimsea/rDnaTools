import os, sys, logging
from collections import namedtuple

BlasrM1 = namedtuple('BlasrM1', ['qname', 'tname', 'qstrand', 'tstrand',
                                 'score', 'pctsimilarity', 
                                 'tstart', 'tend', 'tlength',
                                 'qstart', 'qend', 'qlength',
                                 'ncells'])

log = logging.getLogger()

def file_exists( filename ):
    return os.path.exists( filename ) and os.path.getsize( filename ) > 0

def all_files_exist( filenames ):
    return all( map(file_exists, filenames) )

def is_executable( filepath ):
    if filepath is None:
        return False
    return os.path.isfile( filepath ) and os.access( filepath, os.X_OK )

def get_zmw( read ):
    parts = read.split('/')
    return '/'.join(parts[0:2])

def split_root_from_ext( input_file ):
    root, ext = os.path.splitext( input_file )
    if ext == '.h5':
        root, ext = os.path.splitext( root )
        return (root, ext+'.h5')
    return (root, ext)

def predict_output( input_file, output_type ):
    root, ext = os.path.splitext( input_file )
    return '{0}.{1}'.format(root, output_type)

def return_empty():
    return []

def create_directory( dirname ):
    try:
        os.mkdir( dirname )
    except OSError:
        pass
    if not os.path.isdir( dirname ):
        msg = 'Could not create directory "%s"!' % dirname
        log.error( msg )
        raise OSError( msg )
    return dirname

def which(program):
    """
    Find and return path to local executables  
    """
    fpath, fname = os.path.split(program)
    if fpath:
        if is_executable(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exeFile = os.path.join(path, program)
            if is_executable(exeFile):
                return exeFile
    return None

def validate_input( filename, allowed_suffixes ):
    if isinstance( allowed_suffixes, str ):
        allowed_suffixes = [allowed_suffixes]
    elif isinstance( allowed_suffixes, list ):
        pass
    else:
        msg = 'Allowed suffixes must be either String or List!'
        log.error( msg )
        raise TypeError( msg )
    # First we check whether the file has a valid suffix
    if not any( [filename.endswith(suffix) for suffix in allowed_suffixes] ):
        msg = '"%s" does not have an allowed suffix!' % filename
        log.error( msg )
        raise ValueError( msg )
    # Next we check whether the input file exists where specified
    if not file_exists( filename ):
        msg = 'File %s does not exist!' % filename
        log.error( msg )
        raise OSError( msg )
    # Finally we return the absolute path to the file
    return os.path.abspath( filename )

def validate_output( filename ):
    if filename in [sys.stdout, sys.stderr]:
        return filename
    return os.path.abspath( filename )

def validate_executable( executable ):
    """
    Return the path to an executable if it is valid, otherwise error
    """
    path = which( executable )
    if path is None:
        msg = '"%s" is not a valid executable!' % executable
        log.error( msg )
        raise ValueError( msg )
    return path

def validate_int( variable, value, minimum=None, maximum=None ):
    """
    Validate an integer by confirming it's type and range
    """
    if not isinstance( value, int ):
        msg = '%s is not an integer!' % variable
        log.error( msg )
        raise TypeError( msg )
    if minimum and value < minimum:
        msg = '%s is below the minimum value (%s)!' % (variable, minimum)
        log.error( msg )
        raise ValueError( msg )
    if maximum and value > maximum:
        msg = '%s exceeds the maximum value (%s)!' % (variable, maximum)
        log.error( msg )
        raise ValueError( msg )

def validate_float( variable, value, minimum=None, maximum=None ):
    """
    Validate a float by confirming it's type and range
    """
    if not isinstance( value, float ):
        msg = '%s is not a floating point!' % variable
        log.error( msg )
        raise TypeError( msg )
    if minimum and value < minimum:
        msg = '%s is below the minimum value (%s)!' % (variable, minimum)
        log.error( msg )
        raise ValueError( msg )
    if maximum and value > maximum:
        msg = '%s exceeds the maximum value (%s)!' % (variable, maximum)
        log.error( msg )
        raise ValueError( msg )
