"""small RNA-seq pipeline"""

from setuptools import setup, find_packages

def readme():
    with open('README.rst') as f:
        return f.read()

with open("requirements.txt", "r") as f:
        install_requires = [x.strip() for x in f.readlines() if not x.startswith("#")]


setup(name='seqcluster',
      version='1.1.13',
      description='Small RNA-seq pipeline',
      long_description=readme(),
      classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
      ],
      keywords='RNA-seq miRNA snoRNA NGS',
      url='http://github.com/lpantano/seqcluster',
      author='Lorena Pantano',
      author_email='lpantano@iscb.org',
      license='MIT',
      packages=find_packages(),
      test_suite='nose',
      entry_points={
          'console_scripts': ['seqcluster=seqcluster.command_line:main', 'seqcluster_install=seqcluster.install:main'],
      },
      package_data={'seqclsuter.templates': ['*']},
      include_package_data=True,
      zip_safe=False)
