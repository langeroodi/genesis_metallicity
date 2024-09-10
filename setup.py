from setuptools import setup, find_packages

setup(
    name='genesis_metallicity',
    version='0.1',
    packages=find_packages(),
    description='non-parametric electron temperature and gas-phase metallicity estimation',
    author='Danial Langeroodi',
    author_email='dlangaroudi@gmail.com',
    url='will_be_updated',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent'],
    include_package_data=True,
    package_data={'genesis_metallicity': ['data/*.pkl']},
    python_requires='>=3.6',
)
