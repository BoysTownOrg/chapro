import ctypes
import os
import sys

def get_lib_ext():
    if sys.platform == 'darwin':
        return '.dylib'
    elif sys.platform == 'win32':
        return '.dll'
    return '.so'

# Assuming build dir is at the root
current_dir = os.path.dirname(os.path.abspath(__file__))
lib_path = os.path.join(current_dir, '..', 'build', f'libchapro{get_lib_ext()}')

try:
    libchapro = ctypes.CDLL(lib_path)
except OSError as e:
    print(f"Error loading {lib_path}. Did you compile CHAPRO using CMake first?", file=sys.stderr)
    raise e

# Setup cha_version
libchapro.cha_version.argtypes = []
libchapro.cha_version.restype = ctypes.c_char_p

def version():
    """Returns the version of CHAPRO."""
    return libchapro.cha_version().decode('utf-8')

if __name__ == '__main__':
    print(f"Loaded CHAPRO library version: {version()}")
