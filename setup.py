#!/usr/bin/env python

from os.path import join

def configuration(parent_package='',top_path=None):

    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)

    if fortran:
        config.add_subpackage('odespy')
    config.get_version(join('odespy', 'version.py'))

    return config

if __name__ == '__main__':
    import sys
    if '--no-fortran' in sys.argv:
        sys.argv.remove('--no-fortran') # interfers with distutils cmlargs
        fortran = False
    else:
        fortran = True

    from numpy.distutils.core import setup
    #setup(**configuration(top_path='').todict())
    setup(
        name='odespy',
        url='...',
        download_url='...',
        license='GPL',
        author='Liwei Wang and Hans Petter Langtangen',
        author_email='hpl@simula.no',
        configuration=configuration)

