#!/usr/bin/env python

from os.path import join

def configuration(parent_package='',top_path=None):

    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)

    config.add_subpackage('odespy')
    config.get_version(join('odespy', 'version.py'))

    return config

if __name__ == '__main__':
    from odespy.version import full_version
    import sys
    fortran = True
    if '--no-fortran' in sys.argv:
        sys.argv.remove('--no-fortran') # intefers with distutils sys.argv use
        fortran = False

    name = 'odespy'
    download_url = 'https://github.com/hplgit/odespy'
    url = 'http://hplgit.github.io/odespy/doc/web/index.html'
    author = 'Liwei Wang and Hans Petter Langtangen'
    author_email = 'hpl@simula.no'
    license = 'GPL'
    version = full_version

    if fortran:
        from numpy.distutils.core import setup
        #setup(**configuration(top_path='').todict())
        setup(
            name=name,
            version=version,
            url=url,
            download_url=download_url,
            license=license,
            author=author,
            author_email=author_email,
            configuration=configuration)
    else:
        # Run plain distutils
        from distutils.core import setup
        setup(
            name=name,
            version=version,
            url=url,
            download_url=download_url,
            license=license,
            author=author,
            author_email=author_email,
            description='',
            packages=['odespy'],
            )
