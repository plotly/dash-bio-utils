import io
from setuptools import setup

setup(
    name='dash_bio_utils',
    version='1.0.0',
    url='http://github.com/plotly/dash-bio-utils',
    author='The Plotly Team',
    author_email='dashbio@plot.ly',
    packages=['dash_bio_utils'],
    include_package_data=True,
    description='Simple parsing tools that supplement dash-bio.',
    long_description=io.open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    install_requires=[],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"
    ]
)
