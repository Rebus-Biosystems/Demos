# Demos
A repository to host a collection python-based demos 

# Installation
  - Install Python 3.8 or later on your machine
  - Install VSCode or Similar IDE
  - In a terminal run:
    - python -m pip install --user --upgrade pip 
    - python -m pip install --user virtualenv 
    - python -m venv env_demo
    - (Activate the environment)
      - Windows:  .\env_demo\Scripts\activate
      - Mac: env_demo/Source/Scripts/Activate 
    - pip3 install "napari[all]"
    - pip3 install -r requirements.txt 

# Running the demo
  - Run visualize_all_genes.py 
  - Provide path to the Table folder via command line 
  - example: visualize_all_genes.py' '--rootdir' 'J:\Data\HiFi\2021-08-13_Amsbio mouse brain_Esper HiFi_Imaging buffer_30DM genes\Results\Tables' 


