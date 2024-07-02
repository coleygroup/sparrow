from abc import ABC, abstractmethod
from typing import Dict
import subprocess
import pandas as pd
from pathlib import Path

class RxnClass(ABC):
    def __init__(self):
        self.status_log = None

    @abstractmethod
    def get_rxn_class(rxn: str) -> str:
        """Maps single reaction to its class, if it has one. If not, returns string 0."""

    @abstractmethod
    def get_rxn_classes(rxn_file: Path) -> Dict:
        """Converts a list of reactions into a dictionary of keys, corresponding to reaction classes, which map to lists of reactions"""

class NameRxnClass(RxnClass):
    def __init__(self, 
                 dir: str = None,
                 tmp: str = './tmp/'):
        self.dir = dir
        self.tmp = Path(tmp)
        self.tmp.mkdir(parents=True, exist_ok=True)
        self.no_class_num = 0
        super().__init__()

    def get_rxn_class(self, rxn, attempt=0):
        if attempt > 5:
            print(f"Classifying {rxn} failed")
            self.no_class_num += 1
            return f"Unclassified_{self.no_class_num}"
        if rxn[0] != '>':
            tmp_in = self.tmp / "rxn.smi"
            tmp_out = self.tmp / "rxn.smi.out"
            with open(tmp_in, "w") as f: 
                f.write(rxn)
            # input = open("test.smi", 'w')
            # input.write(rxn)

            cmd_str = " ".join([self.dir + "/namerxn", "-completer", "-addrxnname", "-osmi", str(tmp_in), "-o", str(tmp_out)])
            subprocess.Popen(cmd_str, shell=True)
            
            class_num = ""
            with open(tmp_out, "r") as f:
                output = f.readlines()
            # output = open("test.smi.out", 'r')
            
            try: 
                for line in output:
                    class_num = line.strip().split()[1]
            except:
                return self.get_rxn_class(rxn, attempt + 1)
            
            if class_num != '' and class_num != None:
                return class_num
            
        self.no_class_num += 1
        return f"Unclassified_{self.no_class_num}"

    def get_rxn_classes(self, rxn_file):
        tmp_out = self.tmp / "rxn.smi.out"
        cmd_str = " ".join(["../" + self.dir + "/namerxn", "-completer", "-addrxnname", "-osmi", str(rxn_file), "-o", str(tmp_out)])
        subprocess.Popen(cmd_str, shell=True)

        classes = {}
        with open(tmp_out, "r") as f: 
            output = f.readlines()

        for line in output: 
            class_num = ""
            line = line.strip().split()
            if len(line) >= 2:
                class_num = line[1]
            if class_num == "" or class_num == None:
                self.no_class_num += 1
                class_num = "Unclassified by NameRxn" + str(self.no_class_num)

            if class_num not in classes.keys():
                classes[class_num] = []
            if line != []:
                classes[class_num].append(line[0])
        return classes
    
class LookupClass(RxnClass):
    def __init__(self, 
                 csv_path: str = None):
        # TODO: test the csv -> dict conversion
        df = pd.read_csv(csv_path)
        self.dict = dict(zip(df.iloc[:, 0], df.iloc[:, 1]))
        super().__init__()

    def get_rxn_class(self, rxn):
        return self.dict.get(rxn, "None")
    
    def get_rxn_classes(self, rxns):
        classes = {}
        for rxn in rxns:
            class_name = self.get_rxn_class(rxn)
            if class_name not in classes.keys():
                classes[class_name] = []
            classes[class_name].append(rxn)
        return classes