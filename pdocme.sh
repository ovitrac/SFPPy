# Generate documentation of main classes
# INRAE\olivier.vitrac@agroparistech.fr

# revision history 2022-02-14 $


export PYTHONPATH='/home/olivi/natacha/python':'/home/olivi/anaconda3/lib/python3.9':'/home/olivi/anaconda3/lib/python3.9/lib-dynload':'/home/olivi/anaconda3/lib/python3.9/site-packages':'/home/olivi/anaconda3/lib/python3.9/site-packages/locket-0.2.1-py3.9.egg':'/home/olivi/anaconda3/lib/python3.9/site-packages/IPython/extensions':'/home/olivi/.ipython'

pdoc -f --html -o ./html/ ./patankar/layer.py
pdoc -f --html -o ./html/ ./patankar/senspatankar.py
pdoc -f --html -o ./html/ ./patankar/food.py
pdoc -f --html -o ./html/ ./patankar/__init__.py
