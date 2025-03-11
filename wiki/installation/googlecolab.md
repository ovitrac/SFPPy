# Using SFPPy on Google Colab without Installation 🍏⏩🍎

This tutorial explains how to use `SFPPy` in a [Google Colab](https://colab.research.google.com/) environment without a formal installation. It covers:

1. Cloning the `SFPPy` repository (first time) or updating it (subsequent times).
2. Listing the directory contents.
3. Running an example script.
4. Adding the `SFPPy` folder to Python's `sys.path` to allow direct imports.
5. Run you own code

---

## 1. Clone or Update the Repository

Before running `SFPPy`, you need to ensure the repository is present and up-to-date. The following code checks whether the `SFPPy` folder exists in your current working directory. If it does not exist, it performs a `git clone`. If it does, it updates the repository using `git pull`.

```python
import os

if not os.path.exists("SFPPy"):
    # First time: clone the repository
    !git clone https://github.com/ovitrac/SFPPy.git
else:
    # Subsequent times: update the repository (if needed)
    !git -C SFPPy pull
```



---

## 2. List the Directory Contents

After cloning or updating, verify that the `SFPPy` folder is present by listing the current directory contents:

```bash
%ls
```

---

## 3. Run an Example Script

Test the installation by running an example script provided by SFPPy. Assuming `example1.py` is in the current working directory (or adjust the path as needed):

```bash
!python example1.py
```

---

## 4. Add the SFPPy Folder to `sys.path`

To ensure you can import modules from SFPPy without a traditional installation, add the folder to Python's `sys.path` if it is not already included. The code below retrieves the current working directory and appends it to `sys.path` if necessary:

```python
import sys
import os

SFPPyfolder = os.getcwd()  # Get the current working directory
print(f"Current path: {SFPPyfolder}")

if SFPPyfolder not in sys.path:
    sys.path.append(SFPPyfolder)
    print(f"Added {SFPPyfolder} to sys.path")
else:
    print(f"{SFPPyfolder} already in sys.path")
```





---

## 5. Execute Your SFPPy Code Directly in a Colab Cell

Once you have set up SFPPy, you can simply write or paste any SFPPy-related code into a cell and run it. Below is an example that checks if you can import substances from SFPPy.

### Example: Importing Substances

In this example, we import the `migrant` function from `patankar.loadpubchem` and use it to process a list of substances. You can continue your script as usual after this test.

```python
# %% Place your code below
from patankar.loadpubchem import migrant
list_of_migrants = ["toluene", "limonene", "BHT", "Irganox 1076"]
m = [migrant(s) for s in list_of_migrants]
print(m)
```

> ⚠️ *Note: All imports of SFPPy modules must include the `patankar` prefix.*
---

## Conclusion

By following these steps, you can work with SFPPy in Google Colab without performing a traditional installation. Clone the repository the first time, update it in subsequent sessions, list the contents to confirm its presence, run an example to verify functionality, and finally, ensure the folder is in your Python path for seamless imports





***

<div style="border: 2px solid #4CAF50; border-radius: 8px; padding: 10px; background: linear-gradient(to right, #4CAF50, #FF4D4D); color: white; text-align: center; font-weight: bold;">
  <span style="font-size: 20px;">🍏⏩🍎 <strong>SFPPy for Food Contact Compliance and Risk Assessment</strong></span><br>
  Contact <a href="mailto:olivier.vitrac@gmail.com" style="color: #fff; text-decoration: underline;">Olivier Vitrac</a> for questions |
  <a href="https://github.com/ovitrac/SFPPy" style="color: #fff; text-decoration: underline;">Website</a> |
  <a href="https://ovitrac.github.io/SFPPy/" style="color: #fff; text-decoration: underline;">Documentation</a>
</div>
