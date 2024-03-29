import setuptools

with open("README.md", "r") as fh:
  long_description = fh.read()

setuptools.setup(
  name="VASPsrc",
  version="1.0.0",
  author="GAO Yifan",
  author_email="gawcista@gmail.com",
  description="Python3 scripts for processing VASP output data",
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://github.com/gawcista/vaspsrc",
  packages=setuptools.find_packages(),
  classifiers=[
  "Programming Language :: Python :: 3",
  "License :: OSI Approved :: MIT License",
  "Operating System :: OS Independent",
  ],
)