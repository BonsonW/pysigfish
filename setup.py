try:
    from setuptools import setup, Extension
    from setuptools.command.install import install

except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

import sys

include_dirs = []
# import numpy as np

# creat dummy closures for metadata reading on first parse of setup.py
# that way it picks up the requirements and installs them, then can use them
# for the install.
try:
    import numpy as np
    include_dirs = ['include/', 'sigfish/include/', 'sigfish/slow5lib/include/', 'sigfish/slow5lib/include/', 'sigfish/slow5lib/thirdparty/streamvbyte/include/', np.get_include()]
except ImportError:
    include_dirs = ['include/', 'sigfish/include/', 'sigfish/slow5lib/include/', 'sigfish/slow5lib/include/', 'sigfish/slow5lib/thirdparty/streamvbyte/include/']
    def np(*args, ** kwargs ):
        import numpy as np
        return np(*args, ** kwargs)

try:
    from Cython.Build import build_ext
except ImportError:
    def build_ext(*args, ** kwargs ):
        from Cython.Build import build_ext
        return build_ext(*args, ** kwargs)

# from Cython.Build import build_ext

#adapted from https://github.com/lh3/minimap2/blob/master/setup.py

sources=['python/pysigfish.pyx',
         'sigfish/src/cdtw.c',
         'sigfish/src/events.c',
         'sigfish/src/genref.c',
         'sigfish/src/jnn.c',
         'sigfish/src/rjnn.c',
         'sigfish/src/misc.c',
         'sigfish/src/model.c',
         'sigfish/src/sigfish.c',
         'sigfish/src/thread.c',
         'sigfish/slow5lib/src/slow5.c',
         'sigfish/slow5lib/src/slow5_press.c',
         'sigfish/slow5lib/src/slow5_misc.c',
         'sigfish/slow5lib/src/slow5_idx.c',
         'sigfish/slow5lib/thirdparty/streamvbyte/src/streamvbyte_zigzag.c',
         'sigfish/slow5lib/thirdparty/streamvbyte/src/streamvbyte_decode.c',
         'sigfish/slow5lib/thirdparty/streamvbyte/src/streamvbyte_encode.c']
depends=['python/pysigfish.pxd',
         'python/pysigfish.h',
         'sigfish/src/cdtw.h',
         'sigfish/src/error.h',
         'sigfish/src/jnn.h',
         'sigfish/src/rjnn.h',
         'sigfish/src/khash.h',
         'sigfish/src/kseq.h',
         'sigfish/src/ksort.h',
         'sigfish/src/misc.h',
         'sigfish/src/model.h',
         'sigfish/src/ref.h',
         'sigfish/src/stat.h',
         'sigfish/slow5lib/include/slow5/slow5.h',
         'sigfish/slow5lib/include/slow5/slow5_defs.h',
         'sigfish/slow5lib/include/slow5/slow5_error.h',
         'sigfish/slow5lib/include/slow5/slow5_press.h',
         'sigfish/slow5lib/include/slow5/klib/khash.h',
         'sigfish/slow5lib/include/slow5/klib/kvec.h',
         'sigfish/slow5lib/src/slow5_extra.h',
         'sigfish/slow5lib/src/slow5_idx.h',
         'sigfish/slow5lib/src/slow5_misc.h',
         'sigfish/slow5lib/src/klib/ksort.h',
         'thirdparty/streamvbyte/include/streamvbyte.h',
         'thirdparty/streamvbyte/include/streamvbyte_zigzag.h']
extra_compile_args = ['-g', '-Wall', '-O2', '-std=c99']
# extra_compile_args = []
# os.environ["CFLAGS"] = '-g -Wall -O2 -std=c99'

# arch=platform.machine()
# if arch in ["aarch64", "arm64"]:
#     extra_compile_args.append('-D__ARM_NEON__')
# elif arch in ["aarch64"]:
# 	extra_compile_args.append('-mfpu=neon')
# elif arch in ["x86_64"]:
#     extra_compile_args.extend(['-mssse3'])   # WARNING: ancient x86_64 CPUs don't have SSSE3

libraries = ['m', 'z']
library_dirs = ['.']

extensions = [Extension('pysigfish',
                  sources = sources,
                  depends = depends,
                  extra_compile_args = extra_compile_args,
                  libraries = libraries,
                  include_dirs = include_dirs,
                  library_dirs = library_dirs,
                  language = 'c' )]

def readme():
	with open('docs/pysigfish_api/pysigfish.md') as f:
		return f.read()


setup(
    name = 'pysigfish',
    version='0.0.1',
    url = 'https://github.com/Psy-Fer/pysigfish',
    description='pysigfish python bindings',
    long_description=readme(),
    long_description_content_type='text/markdown',
    author='James Ferguson, Hasindu Gamaarachchi',
    author_email='j.ferguson@garvan.org.au',
    maintainer='James Ferguson',
    maintainer_email='j.ferguson@garvan.org.au',
    license = 'MIT',
    keywords = ['nanopore', 'signal', 'read_until'],
    ext_modules=extensions,
    cmdclass= {'build_ext': build_ext},
    classifiers = [
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: C',
        'Programming Language :: Cython',
        'Programming Language :: Python :: 3',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    python_requires='>=3.8.16',
    install_requires=["numpy"],
    setup_requires=["Cython", "numpy"]
    )
