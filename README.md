# myPRL

A simple Pressure Ruby Luminescence spectra fitting 
program for pressure determination in high pressure experiments.

The program uses the International Practical Pressure Scale (IPPS) "Ruby2020" proposed by the AIRAPT task group (see [Shen *et al.* High Pressure Research, **40**:3, 299-314](https://doi.org/10.1080/08957959.2020.1791107) for details).

Do not hesitate to contact me if you want more information, to report a bug, or to propose improvements.

## Use (GNU/Linux)

### 1) Get a copy of the code:

`$ git clone https://github.com/alexisforestier/myPRL.git`

`$ cd myPRL`

### 2) (Optional) Set up a virtual environment in the code folder:

`$ python3 -m venv .venv`

`$ source .venv/bin/activate`

### 3) Install the required dependencies:

`$ python3 -m pip install -r requirements.txt `

or manually install the required non-native python packages: *numpy*, *pandas*, *matplotlib*, and *scipy*.

### 4) Run as a script:

`$ python3 myPRL.py`

If you are using a virtual environment use

`$ deactivate`

to quit it.

### 5) (Optional) Create an executable using the following command

`$ python3 -m PyInstaller myPRL.spec`

Preferably use a virtual environment if you plan to create an executable as it will only include the required dependencies. A **dist** folder will be created containing the executable.

## Use (Windows)

### 1) Download the source code as a .zip file and extract it where you want

### 2) (Optional) Set up a virtual environment in the code folder:

Open the windows command prompt, use `cd` to navigate to the code folder you extracted before and type:

`> python3 -m venv myvenv`

`> .\myvenv\Scripts\activate.bat`

### 3) Install the required dependencies

In the windows command prompt in the program directory, type:

`> python3 -m pip install -r requirements.txt`

or manually install the required non-native python packages: *numpy*, *pandas*, *matplotlib*, and *scipy*.

### 4) Run the code as a script

`> python3 myPRL.py` 

or run `myPRL.py` using your favorite python interpreter.

### 5) (Optional) Create an executable using the following command

`> python3 -m PyInstaller myPRL.spec`

Preferably use a virtual environment if you plan to create an executable as it will only include the required dependencies. A **dist** folder will be created containing the `.exe` file. You may have to run the command prompt as admin. 