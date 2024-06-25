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
                 dir: str = None):
        self.dir = dir
        self.no_class_num = 0
        super().__init__()

    def get_rxn_class(self, rxn, attempt=0):
        if attempt > 3:
            print("Classifying " + rxn + " failed")
            self.no_class_num += 1
            return "Unclassified by NameRxn" + str(self.no_class_num)
        if rxn[0] != '>':
            input = open("test.smi", 'w')
            input.write(rxn)

            cmd_str = " ".join(["../" + self.dir + "/namerxn", "-completer", "-addrxnname", "-osmi", "test.smi", "-o", "test.smi.out"])
            subprocess.Popen(cmd_str, shell=True)
            output = open("test.smi.out", 'r')
            class_num = ""
            try: 
                for line in output:
                    class_num = line.strip().split()[1]
                    # print(class_num)
            except:
                return self.get_rxn_class(rxn, attempt + 1)
            
            if class_num != '' and class_num != None:
                return class_num
            
        self.no_class_num += 1
        return "Unclassified by NameRxn" + str(self.no_class_num)

    def get_rxn_classes(self, rxn_file):
        cmd_str = " ".join(["../" + self.dir + "/namerxn", "-completer", "-addrxnname", "-osmi", str(rxn_file), "-o", "test.smi.out"])
        subprocess.Popen(cmd_str, shell=True)

        classes = {}
        output = open("test.smi.out", 'r')

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
        # test the csv -> dict conversion
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