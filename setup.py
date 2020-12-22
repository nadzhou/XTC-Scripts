from setuptools import setup 

with open("requirements.txt", "r") as file: 
    requires = file.read()

requires = requires.split("\n")

print(requires)

setup(
    name='XTC script',
    version='1.0',
    scripts=['xtc_read.py'],
    install_requires=requires,
)