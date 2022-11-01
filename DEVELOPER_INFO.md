
# Some notes for developers and testers


## Create a virtual environment for FPEX0

Create a virtual environment in commonly used .venv folder,  
displaying the (FPEX0) prompt when activated.

```
python -m venv .venv --prompt FPEX0
```

Activate the venv as follows:
```
source .venv/bin/activate
```

Deactivate by typing
```
deactivate
```



## Install FPEX0 from Test-PyPI:

The sympy package is not found on the Test-PyPI,  
so we have to tell pip to also use the official PyPI:

```
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple fpex0
```



## Workflows 
For collaborators the repository includes GitHub Actions workflows to simplify package management.

### Build Test Publish to Test-PyPI
This workflow builds, tests and finally publishes the package to the [test package index](https://test.pypi.org). 
<br> It is triggered on every tag being pushed that begins with "tv", or by manual activation (at the
Actions board on the repository).

### Build Publish to PyPI
If no problems occured, you can use this workflow to publish the package to the actual
[Python Package Index](https://pypi.org). 
<br> It is triggered on every tag being pushed that begins with "rv", or by manual activation.

### Creating version tags automatically
To create and push a correct version-tag scripts/deploy.sh can be used on your local machine. Make sure 
that your local repository is up to date before every usage as the script reads the version-info 
locally, but the published package is built using the version-info on the remote. 
Run `bash deploy.sh -t` for a test tag and `bash deploy.sh -r` for a release tag. <br>
For more detailed info run `bash deploy.sh --help`.