from setuptools import setup, find_packages

setup(name='biosapi',
      version='0.1.0',
      description='BIOS API',
      url='https://github.com/Fxe/biosapi',
      author='Filipe Liu',
      author_email='fliu@anl.gov',
      license='MIT',
      packages=find_packages(),
      install_requires=[
          "requests >= 2.22.0"
      ],
      zip_safe=False)