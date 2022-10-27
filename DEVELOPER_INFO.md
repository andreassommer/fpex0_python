
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



## Workflow for testing and publishing to Test-PyPI

For collaborators the repository includes a GitHub Actions workflow that can builds, tests and finally 
publishes the package to test.pypi.org. The workflow is triggered when pushing tags that begin with "v".
<br>
To create and push a correct version-tag scripts/deploy.sh can be used on your local machine. Make sure 
that your local repository is up to date before every usage as the script reads the version-info 
locally, but the published package is built using the version-info on the remote. 
Run `deploy.sh --help` for help.