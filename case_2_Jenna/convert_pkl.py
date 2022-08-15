import os 
import dill 
import pickle 

def convert(old_pkl):
    """
    Convert a Python 2 pickle to Python 3
    """
    # Make a name for the new pickle
    new_pkl = os.path.splitext(os.path.basename(old_pkl))[0]+"_p3.pkl"

    # Convert Python 2 "ObjectType" to Python 3 object
    dill._dill._reverse_typemap["ObjectType"] = object

    # Open the pickle using latin1 encoding
    with open(old_pkl, "rb") as f:
        loaded = pickle.load(f, encoding="bytes")

    # Re-save as Python 3 pickle
    with open(new_pkl, "wb") as outfile:
        pickle.dump(loaded, outfile)

if __name__ == "__main__":

    convert("case_2/encoded_rxn_dict_with_cond_p2.pkl")