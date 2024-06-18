from abc import ABC, abstractmethod
from typing import List, Tuple, Dict, Union
import subprocess
import os
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
            # test with different namerxn string -> plus make the input string the requirement for using name rxn

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
        # convert the csv_path to the dict of rxn to class then save as the self.dict here
        self.dict = {}
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
    
# if __name__ == '__main__':
#     classifier = NameRxnClass('HazELNut')
#     classifier.get_rxn_class("CC(C)(C)OC(=O)N1C[C@@H](O)C[C@H]1C(=O)O.O=C(O)c1ccccc1-c1ccccc1>>CC(C)(C)OC(=O)N1C[C@@H](OC(=O)c2ccccc2-c2ccccc2)C[C@H]1C(=O)O")
#     print("Finished!")
