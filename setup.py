from setuptools import setup, find_packages

setup(name='hm',
      version=0.1,
      description='Classes and methods for hydrological model development',
      url='https://github.com/simonmoulds/hm',
      author='Simon Moulds',
      author_email='sim.moulds@gmail.com',
      license='GPL',
      package_dir= {'hm' : 'hm', 'hm.pcraster' : 'hm/pcraster'},
      packages=['hm', 'hm.pcraster'],
      # packages=find_packages(),
      zip_safe=False
)
      
