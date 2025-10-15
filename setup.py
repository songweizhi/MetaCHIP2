import os
from setuptools import setup


def version():

    setup_dir = os.path.dirname(os.path.realpath(__file__))
    version_file = open(os.path.join(setup_dir, 'MetaCHIP2', 'VERSION'))

    return version_file.readline().strip()


__long_description__ = '''

MetaCHIP2: community-level HGT identification pipeline

Dr Shan Zhang
Dr Weizhi Song (songwz03@gmail.com)
Department of Ocean Science (OCES), 
Hong Kong University of Science and Technology (HKUST), 
Hong Kong
'''


setup(name="MetaCHIP2",
      version=version(),
      long_description=__long_description__,
      license="GPL3+",
      author="Shan Zhang, Torsten Thomas and Weizhi Song",
      author_email="songwz03@gmail.com",
      keywords="Bioinformatics Metagenomics HGT_detection",
      description="HGT detection pipeline",
      url="https://github.com/songweizhi/MetaCHIP2",
      packages=['MetaCHIP2'],
      package_data={'': ['*.r', '*.R', '*.py', 'VERSION', 'Ranger-DTL.mac', 'Ranger-DTL.linux']},
      include_package_data=True,
      install_requires=['biopython', 'matplotlib', 'numpy', 'scipy', 'reportlab', 'ete3', 'DendroPy', 'pycirclize', 'pandas', 'pytz'],
      scripts=['bin/MetaCHIP2'])

