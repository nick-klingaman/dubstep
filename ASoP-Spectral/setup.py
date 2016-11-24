from setuptools import setup


def readme():
    with open('README.rst') as f:
        return f.read()


setup(name='ASoP1_spectral',
      version='1.0',
      description='Analysing Scales of Precipitation: Spectral',
      long_description=readme(),
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: Apache 2.0 License',
        'Programming Language :: Python :: 2.7',
        'Topic :: Text Processing :: Linguistic',
      ],
      keywords='ASoP',
      author='Gill Martin',
      author_email='gill.martin@metoffice.gov.uk',
      license='Apache2.0',
      packages=['ASoP1_spectral'],
      test_suite='nose.collector',
      tests_require=['nose', 'nose-cover3'],
      include_package_data=True,
      zip_safe=False)

