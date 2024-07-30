## python3 setup.py bdist_wheel
## python3 -m twine upload dist/* 
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
     name='mobilitypy',  
     author="Badal Mondal",
     author_email="badalmondal.chembgc@gmail.com",
     maintainer="Badal Mondal",
     maintainer_email="badalmondal.chembgc@gmail.com",
     description="mobilitypy: Python package for mobility calculations in semiconductors",
     long_description=long_description,
     long_description_content_type="text/markdown",
     install_requires=['numpy', 'scipy', 'matplotlib' ,'pandas'],
     url="https://github.com/SemiconductorTransport/mobilitypy",
     packages=setuptools.find_packages(),
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
         "Operating System :: OS Independent",
     ],
 )
