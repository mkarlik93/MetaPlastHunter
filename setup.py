from setuptools import setup, find_packages

if __name__ == '__main__':
    setup(
        name='MetaPlastHunter',
        version='0.0.1',
        author='Michal Karlicki',
        packages=find_packages(exclude=['tests']),
        include_package_data=True,
        url='https://github.com/mkarlik93/MetaPlastHunter',
        license='LICENSE.txt',
        description='Efficent pipeline for fast searching of chloroplast reads in huge metagenomic datasets',
        long_description=open('README.md').read(),
#        entry_points={'console_scripts': ['HiCEnterprise=HiCEnterprise.__main__:main']},
        keywords=['bioinformatics', 'metagenomes', 'chloroplast', 'classification'],
        classifiers=[
            'Development Status :: 2 - Pre-Alpha'
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU License',
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
            'Topic :: Scientific/Engineering :: Bio-Informatics'
        ],
        install_requires=[
            "numpy",
            "scipy",
            "ete3",
            'pandas',
            "matplotlib"
        ],
        setup_requires=['pytest-runner<=3.0.1'],
        tests_require=['pytest'],

    )
