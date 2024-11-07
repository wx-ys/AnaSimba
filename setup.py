import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


with open('requirements.txt') as f:
    install_requires = f.read().splitlines()

setuptools.setup(
    name = "AnaSimba",
    version = "1.0.0",
    author = "Shuai Lu",
    author_email = "lushuai@stu.xmu.edu.cn",
    description = "Simba simulation data analysis",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://github.com/wx-ys/AnaSimba",
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    license = 'MIT',
    packages = setuptools.find_packages(where='src'),
    package_dir = {'': 'src'},
    python_requires = ">=3.8",
    install_requires = install_requires,
    include_package_data = True,
)