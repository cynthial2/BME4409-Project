# BME4409-Project

This is the source code for our BME4409 project. This is a modified version of the code found in the paper our project is based on that simulates a 24-hour period for a T1D patient undergoing physical exercise, found here: https://gitlab.com/csb.ethz/t1d-exercise-model

How to run the code:

1. Download all files as a zip folder to preserve directory structure
2. Create the environment needed to run the code. This can be done using the provided environment.yml file or by using these the original requirements provided

Method 1: Using the environment.yml file to create the environment in the command line
- Navigate to the directory the environment.yml is located in
- Type: conda env create -f environment.yml
- After the environment is finished being created, it should be named: t1d-fullday-modified
- Type: conda activate t1d-fullday-modified
- This should activate the environment, and you should be able to run the code (provided the python files are located in the same directory as the env file)
- Type: python fullday_simulation.py
- This should run the code and you should see the menu

Method 2: Create your own environment using the original requirements provided
- Python (3.7)
- Python packages: numpy (1.19.2), scipy (1.5.0), pandas (1.3.2), matplotlib (3.2.2)

3. Choose a stimulation to test. There are 4 options (baseline, insulin injections, epinephrine injection, periodic carbohydrate consumption during marathon)
4. The generated graphs should be located in the graphs folder
