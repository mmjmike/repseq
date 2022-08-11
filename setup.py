from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()


setup(
    name="repseq",
    url="https://github.com/mmjmike/repseq",
    author="Mikhail Myshkin",
    author_email="mikhail.myshkin@phystech.edu",
    version="0.0.1",
    description="Analyze TCR/BCR repertoires",
    py_modules=["cluster", "stats", "process", "overlap"],
    package_dir={"": "src"},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        ],
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires = [
            "networkx",
        ],

)
