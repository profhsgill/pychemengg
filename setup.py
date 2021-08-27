import setuptools
import os

with open("Readme.md", "r") as fh:
    long_description = fh.read()


# Read version number from Version.txt
# __location__ = os.path.realpath(
#                     os.path.join(os.getcwd(), os.path.dirname(__file__)))
# path_to_version_file =  os.path.dirname(os.path.dirname(__location__))
# path_to_version_file = os.path.dirname(os.path.dirname(os.getcwd()))
version_file = "pychemengg/Version.txt"
quotes = ('"', "'")
with open(version_file,"r") as file:
    for line in file:
        if line.startswith('__version__'):
            _, _, current_version = line.split()
            for quote in quotes:
                if current_version.startswith(quote):
                    current_version = current_version.replace(quote,'')            
            break


setuptools.setup(
    name="pychemengg", 
    version=current_version,
    author="Prof. Harvinder Singh Gill",
    author_email="profhsgill@gmail.com",
    description="A framework for problem solving and critical thinking in chemical engineering.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/profhsgill",
    keywords=["chemical engineering", "chemical-engineering", 
              "heat transfer","heat-transfer", 
              "material balances","material-balances",
              "chemical reaction engineering",
              "thermodynamics", "reactions",
              "reactors", "reaction-engineering",
              "process control", "fluid flow"],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Environment :: Web Environment",
        "Intended Audience :: Education",
        "Intended Audience :: End Users/Desktop",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Topic :: Education",
        "Topic :: Scientific/Engineering",],
    install_requires=["numpy", "scipy"],
    python_requires='>=3',
    include_package_data=True,
    zip_safe=False)