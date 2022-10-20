
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
