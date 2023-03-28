import codecs
import os.path
from setuptools import setup
import subprocess

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            vers = line.split(delim)[1]
            if vers.endswith('-dev'):
                vers += '0+' + _get_repo_hash()
            return vers
    else:
        raise RuntimeError("Unable to find version string.")

def _get_repo_hash():
    """
    Talk to git and find out the tag/hash of our latest commit
    Add '_uc' to the end if the repo has uncommited changes
    """
    try:
        p = subprocess.Popen("git rev-parse --short HEAD".split(' '),
                            stdout=subprocess.PIPE)
    except EnvironmentError:
        print("Couldn't run git to get a version number for setup.py")
        return "detached"

    ver = p.communicate()[0].decode().strip()
    try:
        p = subprocess.Popen("git diff --exit-code".split(' '),
                            stdout=subprocess.PIPE)
        p.communicate()
        if p.returncode:
            ver += '_uc'
    except EnvironmentError:
        pass

    return ver


setup(
    name='laytr',
    version=get_version('laytr/__init__.py'),
    author="ACEnglish",
    author_email="acenglish@gmail.com",
    url="https://github.com/ACEnglish/laytr",
    packages=['laytr'],
    license='MIT',
    description="Library for variant benchmarking stratification",
    long_description=open('README.md', encoding='UTF-8').read(),
    long_description_content_type='text/markdown',
    entry_points={
      'console_scripts': [
         'laytr = laytr.__main__:main'
      ]
    },
    install_requires=[
        "truvari >= 4.0",
        "matplotlib >= 3.7",
        "minisom >= 2.3",
        "seaborn >= 0.12",
        "plotly >= 5.13",
        "scipy >= 1.10",
        "statsmodels >= 0.13.5",
    ],
)
