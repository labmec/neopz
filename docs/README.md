# Documentation

With these source files you can build the documentation yourself to read it locally.

## Dependencies for generating Doxygen output

[Doxygen](https://www.doxygen.nl/index.html) is the only dependency.

## Dependencies for generating documentation using Sphinx

Apart from [Sphinx](https://www.sphinx-doc.org/en/master/), you also will need [Breathe](https://breathe.readthedocs.io/) for extending Sphinx to be able to use Doxygen's output.

Using `pip` (or `pip3`, depending on your Python setup) you can install the additional dependencies with:

```sh
pip3 install sphinx=4.2.0
pip3 install sphinx-rtd-theme
pip3 install breathe==4.31.0
pip3 install sphinx-sitemap
pip3 install sphinxcontrib-bibtex
```
## Generating the Sphinx documentation

While the Doxygen's documentation is far more complete, it can also be quite overwhelming, and some revision is needed. The documentation built with Sphinx is meant to be more concise and just enough to get you familiar with, and curious about, the NeoPZ library.

You can generate the documentation by running the following targets

```sh
make run-doxygen
make run-sphinx
```

And the result can be seen at `build_dir/docs/sphinx/index.html`.