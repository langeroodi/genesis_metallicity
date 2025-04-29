from setuptools import setup, find_packages

with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

with open('requirements.txt') as f:
    required = f.read().splitlines()

VERSION = '1.1.0'
def write_version_py(filename='genesis_metallicity/version.py'):
    with open(filename, 'w') as f:
        f.write(f"__version__ = '{VERSION}'\n")
write_version_py()

setup(
    name='genesis_metallicity',
    version=VERSION,
    packages=find_packages(),
    description='non-parametric gas-phase metallicity and electron temperature estimation',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Danial Langeroodi',
    author_email='dlangaroudi@gmail.com',
    url='https://github.com/langeroodi/genesis_metallicity',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent'],
    include_package_data=True,
    package_data={'genesis_metallicity': ['data/*.pkl']},
    python_requires='>=3.6',
    install_requires=required,
    license='MIT',
    license_files=('LICENSE',),
)
