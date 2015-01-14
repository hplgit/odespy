#!/usr/bin/env python

fortran = True  # compile Fortran libraries?

from os.path import join
import sys

def configuration(parent_package='', top_path=None):

    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    config = Configuration('odespy', parent_package, top_path)

    if fortran:
        config.add_library('_odepack',
                           sources=[join('odepack','solve_odepack.f'),
                                    join('odepack','opkd*.f')])

        config.add_library('_rkc',
                           sources=[join('rkc','solve_rkc.f'),
                                    join('rkc','rkc.f')])
        config.add_library('_rkf45',
                           sources=[join('rkf45','advance_rkf45.f'),
                                    join('rkf45','rkf45.f'),
                                    join('rkf45','rkf45_associate.f')])

        config.add_library('_radau5',
                           sources=[join('radau5','advance_radau5.f'),
                                    join('radau5','radau5.f'),
                                    join('radau5','radaua.f'),])


        # Extensions:
        config.add_extension('_odepack',
                             sources=[join('odepack','odepack.pyf')],
                             libraries=['_odepack'])
        config.add_extension('_rkc',
                             sources=[join('rkc','rkc.pyf')],
                             libraries=['_rkc'])
        config.add_extension('_rkf45',
                             sources=[join('rkf45','rkf45.pyf')],
                             libraries=['_rkf45'])
        config.add_extension('_radau5',
                             sources=[join('radau5','radau5.pyf')],
                             libraries=['_radau5'])

    config.add_data_dir('tests')
    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())

