from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='seqcluster',
      version='0.99',
      description='Small RNA-seq analysis',
      long_description=readme(),
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Topic :: Text Processing :: Linguistic',
      ],
      keywords='RNA-seq miRNA snoRNA NGS',
      url='http://github.com/lpantano/seqcluster',
      author='Lorena Pantano',
      author_email='lpantano@iscb.org',
      license='MIT',
      packages=['seqcluster'],
      install_requires=[
          'pybedtools',
          'pysam'
      ],
      test_suite='nose',
      entry_points={
          'console_scripts': ['seqcluster=seqcluster.command_line:main'],
      },
      include_package_data=True,
      zip_safe=False)