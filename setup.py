import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="dig",
    version="0.0.1",
    author="Wenyin Wei, Shaocheng Liu, Liang Liao",
    author_email="wenyin.wei.ww@gmail.com",
    description="Dig, short for diagnostic, is used in EAST fusion machine to help process diagnostic data. For the time being, dig supports some general functionality like drawing the overview plot of experiment parameters.  ",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/WenyinWei/dig",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
    install_requires=[
        "numpy", "scipy",
        "matplotlib", "plotly",
        "Deprecated", "appdata" # appdata for "mhd_settings.ini" location
    ]
)